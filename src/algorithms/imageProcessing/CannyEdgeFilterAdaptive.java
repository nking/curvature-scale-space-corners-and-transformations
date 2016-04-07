package algorithms.imageProcessing;

import algorithms.compGeometry.HoughTransform;
import algorithms.imageProcessing.util.PairIntWithIndex0;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import java.util.logging.Logger;
import java.util.HashSet;
import algorithms.util.PairInt;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;

/**
 * The CannyEdge filter is an algorithm to operate on an image to
 * replace objects with their edges.  
 * 
 * The class began by following the general advice given in
 * "Performance Analysis of Adaptive Canny Edge Detector 
 * Using Bilateral Filter" by Rashmi, Kumar, Jaiswal, and Saxena, 
 * but made modifications afterwards.
 * Their paper has the following qualities: 
 *<pre>
 * -- instead of a Gaussian filter, uses a bilateral filter
 * -- Otsu's method is used for adaptive threshold levels
 * -- they note that their algorithm performs better on high quality images but 
 *    does not find edges well in low signal-to-noise images.
 * </pre>
 * This class uses 2 one dimensional binomial filters for smoothing.
 * 
 * <pre>
 * Usage:
 * Note, by default, a histogram equalization is not performed.
 * By default, the number of neighbor histogram levels is 1 in the 
 * adaptive thresholding.
 * To use the filter in adaptive mode, use filter.overrideDefaultNumberOfLevels
 * and a number 16 of higher is recommended.
 * 
 * To see a difference in the adaptive approach, run this class on the test
 * image susan-in.gif using filter.overrideDefaultNumberOfLevels(16) 
 * 
 * To adjust the filter to remove lower intensity edges, use
 * filter.setOtsuScaleFactor.  The default factor is 0.45 
 * (determined w/ a checkerboard image).
 * For the Lena test image for example, one might prefer only the brightest
 * edges, so use a higher setting than the default.
 * 
 * The "line drawing mode" is present as a convenience method to process images
 * that are line drawings or shapes filled with solid colors, for the case when 
 * the user needs the theta image for example.  If only the edges are needed, 
 * it would be faster for the user to use 
 *     ImageProcessor.applyAdaptiveMeanThresholding(mask, 1);
 *     followed by
 *     ZhangSuenLineThinner.applyFilter(filterProducts.getGradientXY());
 * Note that line drawing mode starts with the 2 layer filter here, though.
 * </pre>
 * Not ready for use yet...
 * 
 * @author nichole
 */
public class CannyEdgeFilterAdaptive {
              
    /** the factor that the low threshold is below the high threshold in the 
    2 layer filter.
    */
    protected float factorBelowHighThreshold = 2.f;
           
    private EdgeFilterProducts filterProducts = null;
    
    private boolean performHistEq = false;
    
    private boolean performNonMaxSuppr = true;
    
    private boolean debug = false;
    
    private boolean restoreJunctions = false;
    
    private boolean useAlternateNMS = false;
    
    /**
     * the sigma from the blur combined with the sigma present in the gradient
     * are present in this variable by the end of processing.
     * The variable is used to interpret resolution of theta angles, for example.
     * The sigmas are combined using: sigma^2 = sigma_1^2 + sigma_2^2.
     * The FWHM of a gaussian is approx 2.35 * sigma.
     * (HWZI is about 5 * sigma, by the way).
     * So, for the default use of the filter, a sigma of 1 combined with sqrt(1)/2
     * results in a minimum resolution of 3 pixels, hence about 19 degrees.
     */
    private double approxProcessedSigma = 0;
    
    /**
     * default number of levels used in the histogram used to find the
     * high threshold in the 2 - layer thresholding.  it's the number
     *  of levels given to Otsu's multi-level thresholding.
     * for a value of 1, the histogram is performed over the entire image.
     * for a value of 2, histograms are calculated for 2 sets of calculations
     * over the image, one in which the average intensity of neighbors is
     * in bin1 and the other in bin2, etc... the histograms are indexed
     * by the average values of the neighbors.
     */
    private int numberOfLevelsForHistogram = 1;
    
    private boolean useAdaptive2Layer = true;

    private boolean lineDrawingMode = false;
        
    private float otsuScaleFactor = 0.75f;//0.65f;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public CannyEdgeFilterAdaptive() {        
    }
    
    public void setToDebug() {
        debug = true;
    }
    
    /**
     * to enable more complete contours, use this to restore pixels that were
     * removed during non-maximum suppression that disconnected edges and
     * have values above the low threshold of the 2 layer adaptive filter.
     */
    public void setToRestoreJunctions() {
        restoreJunctions = true;
    }
    
    public void setToUseLineDrawingMode() {
        lineDrawingMode = true;      
    }
    
    /**
     * applies a histogram equalization before any processing in order to rescale
     * the data to use the entire range of 0 to 255.
     */
    public void setToPerformHistogramEqualization() {
        performHistEq = true;
    }
    
    public void setToUseAlternateNonMaximumSuppression() {
        useAlternateNMS = true;
    }
    
    /**
     * by default this is 0.45.
     */
    public void setOtsuScaleFactor(float factor) {
        otsuScaleFactor = factor;
    }
    
    public void setToNotUseNonMaximumSuppression() {
        performNonMaxSuppr = false;
    }
    
    /**
     * override the default number of levels used in the histogram used to find the
     * high threshold in the 2 - layer thresholding, where the default is 1.  
     * it's the number of levels given to Otsu's multi-level thresholding.
     * for a value of 1, the histogram is performed over the entire image.
     * for a value of 2, histograms are calculated for 2 sets of calculations
     * over the image, one in which the average intensity of neighbors is
     * in bin1 and the other in bin2, etc... the histograms are indexed
     * by the average values of the neighbors.
     */
    public void overrideDefaultNumberOfLevels(int nLevels) {
        this.numberOfLevelsForHistogram = nLevels;
    }
    
    public void setToUseSingleThresholdIn2LayerFilter() {
        useAdaptive2Layer = false;
    }
    
    /**
     * override the default factor of low threshold below high threshold, which
     * is 2.
     * @param factor 
     */
    public void override2LayerFactorBelowHighThreshold(float factor) {
        factorBelowHighThreshold = factor;
    }
    
    public void applyFilter(final GreyscaleImage input) {
        
        if (input.getWidth() < 3 || input.getHeight() < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }
        
        if (performHistEq) {
            HistogramEqualization hEq = new HistogramEqualization(input);
            hEq.applyFilter();
        }
        
        if (lineDrawingMode) {
            otsuScaleFactor = 1.0f;
            apply2LayerFilter(input, new HashSet<PairIntWithIndex0>(), input);
            if (debug) {
                MiscDebug.writeImage(input, "_after_2_layer_");
            }
            GreyscaleImage mask = input.copyImage();
            for (int i = 0; i < input.getNPixels(); ++i) {
                int v = mask.getValue(i);
                mask.setValue(i, 255 - v);
            }                
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.applyAdaptiveMeanThresholding(mask, 1);
            for (int i = 0; i < mask.getNPixels(); ++i) {
                int v = mask.getValue(i);
                if (v == 255) {
                    input.setValue(i, 0);
                }
                mask.setValue(i, 255 - v);
            }
            if (debug) {
                MiscDebug.writeImage(input, "_adapt_masked_");
            }
            
            filterProducts = createDiffOfGaussians(input);
            approxProcessedSigma = Math.sqrt(
                approxProcessedSigma*approxProcessedSigma + (0.678*0.678));
            
            if (debug) {
                MiscDebug.writeImage(filterProducts.getGradientXY(), "_gXY_");
                MiscDebug.writeImage(filterProducts.getTheta(), "_theta_");
            }
            
            filterProducts.getGradientXY().resetTo(mask);
            
            ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
            lt.applyFilter(filterProducts.getGradientXY());
            
            //exploreHoughTransforms(filterProducts.getGradientXY());
            
            input.resetTo(filterProducts.getGradientXY());
            
            return;
        }
        
        // (1) smooth image using separable binomial filters
        SIGMA sigma = SIGMA.ONE;
        if (sigma.equals(SIGMA.ONE)) {
            ATrousWaveletTransform at = new ATrousWaveletTransform();
            GreyscaleImage smoothed = at.smoothToSigmaOne(input);
            input.resetTo(smoothed);
            approxProcessedSigma = 1;
        } else if (sigma.equals(SIGMA.ZEROPOINTSEVENONE)) {
            ATrousWaveletTransform at = new ATrousWaveletTransform();
            GreyscaleImage smoothed = at.smoothToSigmaZeroPointSevenOne(input);
            input.resetTo(smoothed);
            approxProcessedSigma = Math.sqrt(2.)/2.;
        } else {
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.blur(input, sigma, 0, 255);
            approxProcessedSigma = SIGMA.getValue(sigma);
        }
        
        //(2) create gradient
        // uses a binomial filter for a first derivative gradient, sobel.
        filterProducts = createGradient(input);
        
        GreyscaleImage gradientCopyBeforeThinning = filterProducts.getGradientXY().copyImage();

        approxProcessedSigma = Math.sqrt(
            approxProcessedSigma*approxProcessedSigma + (1./4.));
        
        //MiscDebug.writeImage(filterProducts.getGradientXY(), "_pre_nms_");
        
        Set<PairIntWithIndex0> removedDisconnecting = new HashSet<PairIntWithIndex0>();
        
        //(3) non-maximum suppression
        if (performNonMaxSuppr) {
            applyNonMaximumSuppression(filterProducts, removedDisconnecting);                 
        }
                
        //(4) adaptive 2 layer filter                        
        apply2LayerFilter(filterProducts.getGradientXY(), removedDisconnecting,
            gradientCopyBeforeThinning);
        
        if (restoreJunctions) {
            int minResolution = (int)Math.ceil(2.35 * approxProcessedSigma);
            int minResolvableAngle = (int)Math.ceil(Math.atan2(1, minResolution) * 180./Math.PI);
            MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
            curveHelper.additionalThinning45DegreeEdges2(
                filterProducts.getTheta(), filterProducts.getGradientXY(), 
                minResolvableAngle);
        }
        
        // is this necessary?
        for (int i = 0; i < filterProducts.getGradientXY().getNPixels(); ++i) {
            int v = filterProducts.getGradientXY().getValue(i);
            if (v < 0) {
                filterProducts.getGradientXY().setValue(i, 0);
            }
        }
        
        //exploreHoughTransforms(filterProducts.getGradientXY());
                       
        input.resetTo(filterProducts.getGradientXY());
    }
   
    /*
    Tracing edges through the image and hysteresis thresholding

       It uses 2 "lower" threshold limits, an upper and lower edge limit.

        points below low threshhold are rejected.

        points above the highest threshhold are "sure-edge" pixels.

        for points between highest and low threshold: 
            if connected to "sure-edge"
                pixels, are part of an edge, else discarded.

        Canny recommends:
            high limit
           -----------   should be 2 to 3 based on SNR
            low limit
    
       Note that pixels that are attached to the image boundaries should not
       be nulled in the 2nd stage of the layer if that would disconnect a line.
    
       The threshold is determined adaptively using Otsu's multilevel threshold
        calculation.
    */
    protected void apply2LayerFilter(final GreyscaleImage gradientXY,
        Set<PairIntWithIndex0> removedDisconnecting,
        GreyscaleImage gradientCopyBeforeThinning) {
        
        int n = gradientXY.getNPixels();
        int w = gradientXY.getWidth();
        int h = gradientXY.getHeight();
        
        if (w < 3 || h < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }
    
        int nLevels = numberOfLevelsForHistogram;
        
        OtsuThresholding ot = new OtsuThresholding();
        
        int[] pixelThresholds = null;
        if (useAdaptive2Layer && (nLevels > 1)) {
            pixelThresholds = ot.calculateMultiLevelThreshold256(gradientXY, nLevels);
        }
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        /*
        the paper suggests extending the association of a bright "sure edge"
        pixel to a radius of 2 neighbors (here neighbors have intensity > t_low).
        
        The traversal used in other methods is to visit each pixel and if
           its intensity > thigh, add the pixel and then visit it's neighbors
           to see if their intensities are > tLow.
           (there are redundant visits, but the results are idempotent).
        
        The extended association is better handled from the perspective of
        only adding the pixel currently visited using rules:
            if v < t_low
                continue
            else if v > t_high
                set value in output
            else
                if any neighbor has v > t_high
                    set value of current pixel in output
            else 
                if any neighbor has value > t_low
                   search the 5x5 region and if any neighbor has v > t_high,
                      set value in current pixel in output
        
        Note: 
           If using nLevels > 1, should consider the different results
        of using the threshold determined by the current pixel versus the
        thresholds extracted from the perspective of the "sure edge" pixels
        only, and then a pattern of searching around it.
        TODO: write a second method for the later to try that too...
        */
                
        float tHigh = 0;
        float tLow = 0;
        if (!useAdaptive2Layer || (nLevels == 1)) {
            tHigh = otsuScaleFactor * ot.calculateBinaryThreshold256(gradientXY);
            tLow = tHigh/factorBelowHighThreshold;
        }
        
        GreyscaleImage img2 = gradientXY.createWithDimensions();
        
        for (int i = 0; i < gradientXY.getNPixels(); ++i) {
            
            if (useAdaptive2Layer && (nLevels > 1)) {
                
                tHigh = otsuScaleFactor * pixelThresholds[i];
            
                tLow = tHigh/factorBelowHighThreshold;
            }
            
            int v = gradientXY.getValue(i);
            
            if (v < tLow) {
                continue;
            } else if (v > tHigh) {
                img2.setValue(i, v);
                continue;
            }
            
            int x = gradientXY.getCol(i);
            int y = gradientXY.getRow(i);
            
            boolean foundHigh = false;
            boolean foundMid = false;
            
            for (int k = 0; k < dxs.length; ++k) {                
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                    continue;
                }
                int v2 = gradientXY.getValue(x2, y2);
                if (v2 > tHigh) {
                    foundHigh = true;
                    break;
                } else if (v2 > tLow) {
                    foundMid = true;
                }
            }
            if (foundHigh) {
                img2.setValue(i, v);
                continue;
            }
            if (!foundMid) {
                continue;
            }
            // search the 5 by 5 region for a "sure edge" pixel
            for (int dx = -2; dx <= 2; ++dx) {
                int x2 = x + dx;
                if ((x2 < 0) || (x2 > (w - 1))) {
                    continue;
                }
                for (int dy = -2; dy <= 2; ++dy) {
                    int y2 = y + dy;
                    if ((y2 < 0) || (y2 > (h - 1))) {
                        continue;
                    }
                    if (x2 == x && y2 == y) {
                        continue;
                    }
                    int v2 = gradientXY.getValue(x2, y2);
                    if (v2 > tHigh) {
                        img2.setValue(i, v);
                        foundHigh = true;
                        break;
                    }
                }
                if (foundHigh) {
                    break;
                }
            }
        }
        
        gradientXY.resetTo(img2);
        
        if (useAlternateNMS) {
            MiscDebug.writeImage(gradientXY, "_before_corrections");
            applyPostLineThinningCorrections(gradientXY, 
                gradientCopyBeforeThinning);
            MiscDebug.writeImage(gradientXY, "_before_restore_junctions");
        }
        
        if (restoreJunctions) {
        
            // while have information about the thresholds, will use them to
            // make decisions about restoring pixels that disconnected edges.
            PairInt[][] neighborCoordOffsets
                = AbstractLineThinner.createCoordinatePointsForEightNeighbors(
                    0, 0);

            for (PairIntWithIndex0 p : removedDisconnecting) {
                                
                if (ImageSegmentation.doesDisconnect(gradientXY,
                    neighborCoordOffsets, p.getX(), p.getY())) {

                    int x = p.getX();
                    int y = p.getY();
                    int i = gradientXY.getIndex(x, y);

                    if (useAdaptive2Layer && (nLevels > 1)) {
                        tHigh = otsuScaleFactor * pixelThresholds[i];
                        tLow = tHigh/factorBelowHighThreshold;
                    }

                    if (p.getPixIndex() > tLow) {
                        if (isAdjacentToAHorizOrVertLine(gradientXY, x, y, 3)) {
                            gradientXY.setValue(x, y, 255);
                        }
                    }
                }
            }
            MiscDebug.writeImage(gradientXY, "_after_restore_junctions");
            if (useAlternateNMS) {
                applyPostLineThinningCorrections(gradientXY, 
                    gradientCopyBeforeThinning);
                MiscDebug.writeImage(gradientXY, "_after_corrections_2");
            }
        }
    }
    
    /**
     * convolve the image with a Sobel X 1D kernel which is the same as a 
     * Gaussian first derivative with sigma = sqrt(1)/2.
     * FWHM=2.355*sigma
     * 
     * @param input
     * @return 
     */
    private GreyscaleImage getGradientX1D(final GreyscaleImage input) {
                
        return getGradient1D(input, true);
    }
    
    /**
     * convolve the image with a Sobel Y 1D kernel which is the same as a 
     * Gaussian first derivative with sigma = sqrt(1)/2.
     * 
     * @param input
     * @return 
     */
    private GreyscaleImage getGradientY1D(final GreyscaleImage input) {
                
        return getGradient1D(input, false);
    }
    
    /**
     * convolve the image with a Sobel 1D kernel which is the same as a 
     * Gaussian first derivative with sigma = sqrt(1)/2.
     * 
     * @param input
     * @return 
     */
    private GreyscaleImage getGradient1D(final GreyscaleImage input, 
        boolean calculateForX) {
        
        log.fine("getGradientID calculateForX=" + calculateForX);
                
        // 0.5f, -0.0f, -0.5f
        float[] kernel = Gaussian1DFirstDeriv.getKernel(
            //SIGMA.ZEROPOINTSEVENONE);
            SIGMA.ZEROPOINTFIVE);
                
        GreyscaleImage output = input.copyToSignedImage();
        
        apply1DKernelToImage(output, kernel, calculateForX);
        
        return output;
    }
    
    private void apply1DKernelToImage(final GreyscaleImage input, 
        float[] kernel, boolean calculateForX) {
        
        log.fine("apply1DKernelToImage calculateForX=" + calculateForX);
        
        GreyscaleImage output = input.copyImage();
      
        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();
        
        for (int i = 0; i < input.getWidth(); i++) {
           
            for (int j = 0; j < input.getHeight(); j++) {
                
                double conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, calculateForX);
                
                int g = (int) conv;
                
                // because the values may be combined with other images or
                // involved in other operations such as adding in quadrature,
                // don't limit the range to be between 0 and 255
                output.setValue(i, j, g);
            }
        }
        
        input.resetTo(output);
    }
    
    /**
     * construct the gradient in X, gradient in Y, their combined average and
     * theta image
     * using two 1-D passes of a Sobel 1D kernel which is the same as a 
     * Gaussian first derivative with sigma = sqrt(1)/2 where FWHM=2.355*sigma.
     * The theta image has range 0 t 360.
     * 
     * @param img
     * @return 
     */
    protected EdgeFilterProducts createGradient(final GreyscaleImage img) {
        
        GreyscaleImage gX, gY, g, theta;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        gX = getGradientX1D(img);

        gY = getGradientY1D(img);

        g = imageProcessor.combineConvolvedImages(gX, gY);
        
        // the theta is in range 0 to 360
        theta = imageProcessor.computeTheta360_0(gX, gY);
        
        if (debug) {
            // when replace the aspect library, put these renders in the 
            //   equivalent replacement
            long ts = MiscDebug.getCurrentTimeFormatted();
            MiscDebug.writeImage(theta, "_theta_" + ts);
            MiscDebug.writeImage(gX, "_gX_" + ts);
            MiscDebug.writeImage(gY, "_gY_" + ts);
            /*
            int x = 37; int y = 163;
            System.out.println("(" + x + ", " + y + ") math.atan2(" + gY.getValue(x, y)
                + "," + gX.getValue(x, y) + ")*180./math.pi=" + 
                theta.getValue(x, y));*/
        }
        
        EdgeFilterProducts efp = new EdgeFilterProducts();
        efp.setGradientX(gX);
        efp.setGradientY(gY);
        efp.setGradientXY(g);
        efp.setTheta(theta);
        
        return efp;
    }

    /**
     * note, uses small gaussian, so should only be used when very little or
     * no smoothing has occurred (line drawings, for example might be best
     * handled with this and no smoothing).
     * @param img
     * @return 
     */
    protected EdgeFilterProducts createDiffOfGaussians(final GreyscaleImage img) {
        
        GreyscaleImage gX, gY, g, theta;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        gX = createGradientFromDiffOfGauss(img, true);

        gY = createGradientFromDiffOfGauss(img, false);
            
        g = imageProcessor.combineConvolvedImages(gX, gY);
        
        // the theta is in range 0 to 360
        theta = imageProcessor.computeTheta360_0(gX, gY);
        
        EdgeFilterProducts efp = new EdgeFilterProducts();
        efp.setGradientX(gX);
        efp.setGradientY(gY);
        efp.setGradientXY(g);
        efp.setTheta(theta);
        
        return efp;
    }
    
    private GreyscaleImage createGradientFromDiffOfGauss(
        final GreyscaleImage img, boolean calculateForX) {
        /*
        what looks really good is computationally too long:
            g1 is convolved by a kernel > 1000*g0 kernel which uses sigma=1
        */
        
        /* 
        for line drawings, use sigma1 = 0.42466090014400953f and sigma2 = 2*sigma1
        */
        
        float resultOne = 0.42466090014400953f;
        
        float sigma = 1.f * resultOne;
        
        float[] kernel = Gaussian1D.getKernel(sigma);

        float[] kernel2 = Gaussian1D.getKernel(sigma * 1.6f);
            
        GreyscaleImage g0 = img.copyImage();
        
        GreyscaleImage g1 = img.copyImage();
          
        apply1DKernelToImage(g0, kernel, calculateForX);
                
        apply1DKernelToImage(g1, kernel2, calculateForX);
       
        ImageProcessor ImageProcessor = new ImageProcessor();
        
        GreyscaleImage g = ImageProcessor.subtractImages(g1, g0);
        
        return g;
    }

    private void applyAlternateNonMaximumSuppression(EdgeFilterProducts filterProducts,
        Set<PairIntWithIndex0> disconnectingRemovals) {
        
         NonMaximumSuppression nms = new NonMaximumSuppression();
            
         //TODO: radius can be adjusted from 1 to higher
         nms.<PairIntWithIndex0>nonmaxsup(filterProducts.getGradientXY(),
            filterProducts.getTheta(), 1.5, false, disconnectingRemovals);
         
    }
      
    private void applyNonMaximumSuppression(EdgeFilterProducts filterProducts,
        Set<PairIntWithIndex0> disconnectingRemovals) {
        
        if (useAlternateNMS) {
            
            applyAlternateNonMaximumSuppression(filterProducts, disconnectingRemovals);
            
            return;
        }
        
        
        int minResolution = (int)Math.ceil(2.35 * approxProcessedSigma);
        int minResolvableAngle = (int)Math.ceil(Math.atan2(1, minResolution) * 180./Math.PI);
        
        GreyscaleImage theta = filterProducts.getTheta();
        
        GreyscaleImage img = filterProducts.getGradientXY();
                   
        int n = img.getNPixels();
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        PairInt[][] neighborCoordOffsets
            = AbstractLineThinner.createCoordinatePointsForEightNeighbors(
                0, 0);
        
        /*
         *           Y
         *          90
         *     135   |    +45
         *           |
         *   180------------ 0   X
         *           |
         *    225    |   315
         *          270
        */
        for (int x = 1; x < (w - 1); ++x) {
            for (int y = 1; y < (h - 1); ++y) {
            
                int pixIdx = theta.getInternalIndex(x, y);
           
                int t = theta.getValue(pixIdx);

                int compDx1, compDy1, compDx2, compDy2;
                
                // angles within minResolvableAngle of horizontal or vertical
                // are orthogonal to the pattern present for other pixels
                
                if ((Math.abs(t - 0) < minResolvableAngle) || 
                    (Math.abs(180 - t) < minResolvableAngle) || 
                    (Math.abs(360 - t) < minResolvableAngle)) {
                    // measured 0 or 180, but should be 90 or 270
                    compDx1 = 1;
                    compDy1 = 0;
                    compDx2 = -1;
                    compDy2 = 0;
                } else if ((Math.abs(90 - t) < minResolvableAngle) ||
                    (Math.abs(270 - t) < minResolvableAngle) ) {
                    // measured 90 or 270, but should be 0 or 180
                    compDx1 = 0;
                    compDy1 = 1;
                    compDx2 = 0;
                    compDy2 = -1;
                } else if ((Math.abs(t - 90) < 22.5) || 
                    (Math.abs(270 - t) < 22.5)) {
                    // is within 45 degree cone centered at 90 or 270
                    compDx1 = 0;
                    compDy1 = 1;
                    compDx2 = 0;
                    compDy2 = -1;
                } else if ((Math.abs(t - 0) < 22.5) ||
                    (Math.abs(180 - t) < 22.5) | (Math.abs(360 - t) < 22.5)) {
                    // is within 45 degree cone centered at 0 or 180
                    compDx1 = -1;
                    compDy1 = 0;
                    compDx2 = 1;
                    compDy2 = 0;
                } else if ((Math.abs(t - 45) < 22.5) || 
                    (Math.abs(225 - t) < 22.5)) {
                    // is within 45 degree cone centered at 45 or 225
                    compDx1 = 1; 
                    compDy1 = 1;
                    compDx2 = -1;
                    compDy2 = -1;
                } else if ((Math.abs(t - 135) < 22.5) || 
                    (Math.abs(315 - t) < 22.5)) {
                    // is within 45 degree cone centered at 135 or 315
                    compDx1 = -1; 
                    compDy1 = 1;
                    compDx2 = 1;
                    compDy2 = -1;
                } else {
                    throw new IllegalStateException("Error in algorithm bounds: " 
                        +  " t=" + t);
                }

                int v = img.getValue(x, y);
            
                int v1 = img.getValue(x + compDx1, y + compDy1);
                int v2 = img.getValue(x + compDx2, y + compDy2);
            
                if ((v < v1) || (v < v2)) {
                    
                    if (restoreJunctions) {
                        boolean disconnects = ImageSegmentation.doesDisconnect(img,
                            neighborCoordOffsets, x, y);

                        if (disconnects) {
                            PairIntWithIndex0 p = new PairIntWithIndex0(x, y, v);
                            disconnectingRemovals.add(p);
                        }
                    }
                    
                    img.setValue(x, y, 0);
                }
            }
        }
        
        // because minResolvableAngle for default settings is smaller than 45/2., 
        // the resolution correction is usually not necessary, but does no
        // harm to include it 
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        curveHelper.additionalThinning45DegreeEdges2(theta, img, minResolvableAngle);
    }
    
    private void exploreHoughTransforms(GreyscaleImage gradientXY) {
        
        int n = gradientXY.getNPixels();
        int w = gradientXY.getWidth();
        int h = gradientXY.getHeight();
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        for (int i = 0; i < n; ++i) {
            if (gradientXY.getValue(i) > 0) {
                points.add(new PairInt(gradientXY.getCol(i), gradientXY.getRow(i)));
            }
        }
        
        HoughTransform ht = new HoughTransform();
        Map<Set<PairInt>, PairInt> lines = ht.findContiguousLines(points, 3);
        
        int count = 0;
        Image tmp = gradientXY.copyToColorGreyscale();
        for (int i = 0; i < tmp.getNPixels(); ++i) {
            if (tmp.getR(i) > 0) {
                tmp.setRGB(i, 0, 0, 0);
            } else {
                tmp.setRGB(i, 255, 255, 255);
            }
        }
        for (Entry<Set<PairInt>, PairInt> entry : lines.entrySet()) {
            Set<PairInt> group = entry.getKey();
            int[] clr = ImageIOHelper.getNextRGB(count);
            for (PairInt p : group) {
                tmp.setRGB(p.getX(), p.getY(), clr[0], clr[1], clr[2]);
            }
            count++;
        }
        MiscDebug.writeImage(tmp, "_lines_");
    }

    private void correctThetaBelowResolution(GreyscaleImage theta, 
        GreyscaleImage gXYMask, double theApproxProcessedSigma) {
        
        int minResolution = (int)Math.ceil(2.35 * theApproxProcessedSigma);
        int minResolvableAngle = (int)Math.ceil(Math.atan2(1, minResolution) * 180./Math.PI);
        
        int n = theta.getNPixels();
        int w = theta.getWidth();
        int h = theta.getHeight();
        
        int mr = (int)Math.ceil(this.approxProcessedSigma * 2.35);
        Image tmp = new Image(theta.getWidth(), theta.getHeight());
        for (int i = 0; i < n; ++i) {
            
            if (gXYMask.getValue(i) == 0) {
                continue;
            }
            
            int t = theta.getValue(i);

            if (
                (Math.abs(t - 0) < minResolvableAngle) || 
                (Math.abs(180 - t) < minResolvableAngle) ||
                (Math.abs(360 - t) < minResolvableAngle) ||
                (Math.abs(90 - t) < minResolvableAngle) ||
                (Math.abs(270 - t) < minResolvableAngle) ) {
                
                if ((Math.abs(filterProducts.getGradientX().getValue(i)) <= mr) || 
                    (Math.abs(filterProducts.getGradientY().getValue(i)) <= mr)) {
                    
                    tmp.setRGB(i, 255, 0, 0);
                    
                } else {
                    
                    tmp.setRGB(i, 0, 255, 0);
                }
            }
        }
        
        MiscDebug.writeImage(tmp, "_gXY_BELOW_RESOLUTION_");
        
        for (int i = 0; i < n; ++i) {
            
            if (gXYMask.getValue(i) == 0) {
                continue;
            }
            
            int t = theta.getValue(i);
            
            if (
                (Math.abs(t - 0) < minResolvableAngle) || 
                (Math.abs(180 - t) < minResolvableAngle) ||
                (Math.abs(360 - t) < minResolvableAngle) ||
                (Math.abs(90 - t) < minResolvableAngle) ||
                (Math.abs(270 - t) < minResolvableAngle) ) {
                
                int t1 = t - 90;
                int t2 = t + 90;
                if (t1 > minResolvableAngle) {
                    theta.setValue(i, t1);
                } else if (t2 < (360 - minResolvableAngle)) {
                    theta.setValue(i, t2);
                } else {
                    throw new IllegalStateException("Error in algorithm.  theta="
                        + t);
                }
            }            
        }
    }

    private boolean isAdjacentToAHorizOrVertLine(GreyscaleImage gradientXY, 
        int x, int y, int minLineSize) {
        
        /*
           1  1  1
           1  0  1 
           1  1  1
        
        if find a neighbor, check it for being part of _ or | of size minLineSize        
        */
        int w = gradientXY.getWidth();
        int h = gradientXY.getHeight();
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        for (int k = 0; k < dxs.length; ++k) {
            int x2 = x + dxs[k];
            int y2 = y + dys[k];
            if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                continue;
            }
            if (gradientXY.getValue(x2, y2) == 0) {
                continue;
            }
            
            /*
            if (x2, y2) is part of a vertical line, return true
            
              (x2, y2) 
                 *   
                 *    
            */
            int c = 0;
            for (int dy2 = ((-minLineSize) + 1); dy2 < 0; ++dy2) {
                int y3 = y2 + dy2;
                if ((y3 < 0) || (y3 > (h - 1)) || (gradientXY.getValue(x2, y3) == 0)) {
                    continue;
                }
                c++;
            }
            if (c == (minLineSize - 1)) {
                return true;
            }
            /*            
                   *     
                (x2, y2)
                   *    
            */
            c = 0;
            for (int dy2 = (-minLineSize/2); dy2 <= (minLineSize/2); ++dy2) {
                int y3 = y2 + dy2;
                if ((y3 < 0) || (y3 > (h - 1)) || (gradientXY.getValue(x2, y3) == 0)
                    || (y3 == y2)) {
                    continue;
                }
                c++;
            }
            if (c == (minLineSize - 1)) {
                return true;
            }
            /*            
                   *
                   *
               (x2, y2)
            */
            c = 0;
            for (int dy2 = 1; dy2 < minLineSize; ++dy2) {
                int y3 = y2 + dy2;
                if ((y3 < 0) || (y3 > (h - 1)) || (gradientXY.getValue(x2, y3) == 0)) {
                    continue;
                }
                c++;
            }
            if (c == (minLineSize - 1)) {
                return true;
            }
            
            /*
            if (x2, y2) is part of a horizontal line, return true
            
              *   *   (x2, y2) 
            */
            c = 0;
            for (int dx2 = ((-minLineSize) + 1); dx2 < 0; ++dx2) {
                int x3 = x2 + dx2;
                if ((x3 < 0) || (x3 > (w - 1)) || (gradientXY.getValue(x3, y2) == 0)) {
                    continue;
                }            
                c++;
            }
            if (c == (minLineSize - 1)) {
                return true;
            }
            
            /*            
                  *  (x2, y2)   *
            */
            c = 0;
            for (int dx2 = (-minLineSize/2); dx2 <= (minLineSize/2); ++dx2) {
                int x3 = x2 + dx2;
                if ((x3 < 0) || (x3 > (w - 1)) || (gradientXY.getValue(x3, y2) == 0)
                    || (x3 == x2)) {
                    continue;
                }
                c++;
            }
            if (c == (minLineSize - 1)) {
                return true;
            }
            
            /*            
                  (x2, y2)   *    *
            */
            c = 0;
            for (int dx2 = 1; dx2 < minLineSize; ++dx2) {
                int x3 = x2 + dx2;
                if ((x3 < 0) || (x3 > (w - 1)) || (gradientXY.getValue(x3, y2) == 0)) {
                    continue;
                }
                c++;
            }
            if (c == (minLineSize - 1)) {
                return true;
            }
        }
        
        return false;
    }

    public GreyscaleImage createGradientPSFForTesting() {
                
        float[] k = Gaussian1D.getBinomialKernelSigmaZeroPointFive();
        
        GreyscaleImage gX = new GreyscaleImage(3, 3);
        gX.setValue(1, 1, 127);        
        apply1DKernelToImage(gX, k, true);
        
        GreyscaleImage gY = new GreyscaleImage(3, 3);
        gY.setValue(1, 1, 127);
        apply1DKernelToImage(gY, k, false);
        
        ImageProcessor imageProcessor = new ImageProcessor();
        GreyscaleImage g = imageProcessor.combineConvolvedImages(gX, gY);
        
        return g;
    }

    /**
     * @return the filterProducts
     */
    public EdgeFilterProducts getFilterProducts() {
        return filterProducts;
    }

    public void setSetters(CannyEdgeFilterSettings settings) {
        
        if (settings.getNormalizeByHistogram()) {
            setToPerformHistogramEqualization();
        }
        
        if (settings.getUseLineDrawingMode()) {
            setToUseLineDrawingMode();
        }
    }

    private void applyPostLineThinningCorrections(GreyscaleImage gradientXY,
        GreyscaleImage valuesBeforeThinning) {

        Set<PairInt> correctedPoints = new HashSet<PairInt>();
        
        for (int i = 0; i < gradientXY.getWidth(); ++i) {
            for (int j = 0; j < gradientXY.getHeight(); ++j) {
                if (gradientXY.getValue(i, j) > 0) {
                    correctedPoints.add(new PairInt(i, j));
                }
            }
        }
        
        int n0 = gradientXY.getWidth();
        int n1 = gradientXY.getHeight();
        
        PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
        pltc.correctForHolePattern100(correctedPoints, n0, n1);
        pltc.correctForLineHatHoriz(correctedPoints, n0, n1);
        pltc.correctForLineHatVert(correctedPoints, n0, n1);
        
        GreyscaleImage out = gradientXY.createWithDimensions();
        for (PairInt p : correctedPoints) {
                                    
            int v = valuesBeforeThinning.getValue(p);
            
            out.setValue(p.getX(), p.getY(), v);
        }
        
        gradientXY.resetTo(out);
    }
    
}
