package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import java.util.logging.Logger;
import java.util.HashSet;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

/**
 * NOTE: need to edit these comments.  for this color version of canny,
 * using the "L" and "C" images of LCH CIELUV color space.
 * May also add use of "H"...still experimenting.
 * 
 * 
 * 
 * The CannyEdge filter is an algorithm to operate on an image to
 * replace objects with their edges.  
 * 
 * The class began by following the general advice given in
 * "Performance Analysis of Adaptive Canny Edge Detector 
 * Using Bilateral Filter" by Rashmi, Kumar, Jaiswal, and Saxena, 
 * but made modifications afterwards and added a "C" gradient of colorspace
 * LCH to the greyscale gradient.
 * 
 * Their paper has the following qualities: 
 *<pre>
 * -- instead of a Gaussian filter, uses a bilateral filter
 * -- an adaptive threshold algorithm is used in the 2 layer filter
</pre>
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
 * image susan-in.gif using filter.overrideToUseAdaptiveThreshold() 
 * 
 * To adjust the filter to remove lower intensity edges, use
 * filter.setOtsuScaleFactor.  The default factor is 0.45 
 * (determined w/ a checkerboard image).
 * For the Lena test image for example, one might prefer only the brightest
 * edges, so use a higher setting than the default.
 * 
 * Not ready for use yet...
 * 
 * @author nichole
 */
public class CannyEdgeColorAdaptive {
              
    /** the factor that the low threshold is below the high threshold in the 
    2 layer filter.
    */
    protected float factorBelowHighThreshold = 2.f;
           
    private EdgeFilterProducts filterProducts = null;
    
    private EdgeFilterProducts lFilterProducts = null;
    
    private EdgeFilterProducts cFilterProducts = null;
    
    /**
     * by default, adds the gradients of L and C from colorspace LCH, scaled
     * to their own color space vector maxima then the maximum possible in their
     * color vectors is mapped to 255.
     * If "scaleGrdients" is set to true, the individual sobel gradients are
     * scaled so that their maximum values present in the image are transformed
     * to values 255.
     */
    private boolean scaleGradients = false;
        
    private boolean performNonMaxSuppr = true;
    
    private boolean debug = false;
    
    private boolean restoreJunctions = true;
            
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
    
    private boolean useAdaptiveThreshold = false;
    
    private boolean useAdaptive2Layer = true;
        
    private float otsuScaleFactor = 0.75f;//0.65f;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected boolean useLineThinner = true;
    
    public CannyEdgeColorAdaptive() {        
    }
    
    public void setToDebug() {
        debug = true;
    }
    
    /**
     * by default, adds the gradients of L and C from colorspace LCH, scaled
     * to their own color space vector maxima then the maximum possible in their
     * color vectors is mapped to 255.
     * If "scaleGrdients" is set to true, the individual sobel gradients are
     * scaled so that their maximum values present in the image are transformed
     * to values 255.
     */
    public void overrideToScaleGradients() {
        scaleGradients = true;
    }
    
    public void overrideToNotUseLineThinner() {
        useLineThinner = false;
    }
    
    /**
     * to enable more complete contours, use this to restore pixels that were
     * removed during non-maximum suppression that disconnected edges and
     * have values above the low threshold of the 2 layer adaptive filter.
     */
    public void setToNotRestoreJunctions() {
        restoreJunctions = false;
    }
    
    /**
     * by default this is 0.45.
     * @param factor
     */
    public void setOtsuScaleFactor(float factor) {
        otsuScaleFactor = factor;
    }
    
    public void setToNotUseNonMaximumSuppression() {
        performNonMaxSuppr = false;
    }
    
    /**
     * set this to use the adaptive threshold in the 2 layer
     * filter.  it adjusts the threshold by regions of size
     * 15.  Note that if the image has alot of noise, this
     * will include alot of noise in the result.
     */
    public void overrideToUseAdaptiveThreshold() {
        useAdaptiveThreshold = true;
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
    
    /**
     * apply the filter.  note that unlike the other canny filters in this
     * project, the input is not modified.
     * @param input 
     */
    public void applyFilter(Image input) {
        
        if (input.getWidth() < 3 || input.getHeight() < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }
            
        // (1) smooth image using separable binomial filters
        SIGMA sigma = null;//SIGMA.ZEROPOINTSEVENONE;
        if (sigma == null) {
            //no smoothing
            approxProcessedSigma = 0.4;
        } else if (sigma.equals(SIGMA.ONE)) {
            ATrousWaveletTransform at = new ATrousWaveletTransform();
            Image smoothed = at.smoothToSigmaOne(input);
            input = smoothed;
            approxProcessedSigma = 1;
        } else if (sigma.equals(SIGMA.ZEROPOINTSEVENONE)) {
            ATrousWaveletTransform at = new ATrousWaveletTransform();
            Image smoothed = at.smoothToSigmaZeroPointSevenOne(input);
            input = smoothed;
            approxProcessedSigma = Math.sqrt(2.)/2.;
        } else {
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.blur(input, sigma, 0, 255);
            approxProcessedSigma = SIGMA.getValue(sigma);
        }
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        GreyscaleImage[] lch = imageProcessor.createLCHForLUV(input);
                
        //(2) create gradient
        // uses a binomial filter for a first derivative gradient, sobel.
        EdgeFilterProducts[] eps = createGradients(lch);
        //EdgeFilterProducts[] eps = createDiffOfGaussians(lch);
        lFilterProducts = eps[0];
        cFilterProducts = eps[1];
        
        GreyscaleImage[] gXYsOrig = new GreyscaleImage[]{
            lFilterProducts.getGradientXY().copyImage(),
            cFilterProducts.getGradientXY().copyImage()
        };

        approxProcessedSigma = Math.sqrt(
            approxProcessedSigma*approxProcessedSigma + (1./4.));
                
        List<Set<PairInt>> rmvdDisconnecting = new ArrayList<Set<PairInt>>();
        rmvdDisconnecting.add(new HashSet<PairInt>());
        rmvdDisconnecting.add(new HashSet<PairInt>());
        
        //(3) non-maximum suppression
        if (performNonMaxSuppr) {
            for (int i = 0; i < 2; ++i) {
                applyNonMaximumSuppression(eps[i], 
                    rmvdDisconnecting.get(i));
            }
        }
        
        if (debug) {
            for (int i = 0; i < 2; ++i) {
                MiscDebug.writeImage(eps[i].getGradientXY(),
                    "_after_nms_" + i);
            }
        }
        
        //(4) adaptive 2 layer filter
        this.filterProducts = apply2LayerFilter(eps, rmvdDisconnecting, gXYsOrig);
                
        if (restoreJunctions) {
            
            //TODO: consdier whether should apply this
            //   to the L and C grdients separately
            
            int minResolution = (int)Math.ceil(2.35 * approxProcessedSigma);
            int minResolvableAngle = (int)Math.ceil(
                Math.atan2(1, minResolution) * 180./Math.PI);
            if (minResolvableAngle < 0) {
                minResolvableAngle *= -1;
            }
            
            MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
            curveHelper.additionalThinning45DegreeEdges2(
                filterProducts.getTheta(), filterProducts.getGradientXY(), 
                minResolvableAngle);
        }
                                           
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
    protected EdgeFilterProducts apply2LayerFilter(EdgeFilterProducts[] lcFilterProducts, 
        List<Set<PairInt>> lcRemovedDisconnecting,
        GreyscaleImage[] gXYsOrig) {
                
        int n = lcFilterProducts[0].getGradientXY().getNPixels();
        int w = lcFilterProducts[0].getGradientXY().getWidth();
        int h = lcFilterProducts[0].getGradientXY().getHeight();
        
        if (w < 3 || h < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }
           
        ImageProcessor imageProcessor = new ImageProcessor();

        // 1st dimension is 0 or 1 and is L or C data
        double[][][] lcThreshImg = null;
       
        if (useAdaptive2Layer && useAdaptiveThreshold) {
            
            lcThreshImg = new double[2][][];
            
            AdaptiveThresholding th = new AdaptiveThresholding();
            
            for (int i = 0; i < 2; ++i) {
            
                lcThreshImg[i] = th.createAdaptiveThresholdImage(
                    imageProcessor.copy(
                        lcFilterProducts[i].getGradientXY()), 15, 0.2);
                        
                if (debug) {//DEBUG
                    double[][] imgCp = imageProcessor.copy(
                        lcFilterProducts[i].getGradientXY());
                    for (int ii = 0; ii < w; ++ii) {
                        for (int jj = 0; jj < h; ++jj) {
                            double t = lcThreshImg[i][ii][jj];
                            if (imgCp[ii][jj] > t) {
                                imgCp[ii][jj] = 255.;
                            } else {
                                imgCp[ii][jj] = 0;
                            }
                        }
                    }
                    MiscDebug.writeImage(imgCp, "img_" + i + "_thresholded_.png");
                    MiscDebug.writeImage(lcThreshImg[i], "img_" + i 
                        + "_adaptive_threshold_.png");
                }
            }
        }
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        /*
        the paper suggests extending the association of a bright "sure edge"
        pixel to a radius of 2 neighbors (here neighbors have intensity > t_low).
        
        The traversal used in other methods is to visit each pixel and if
           its intensity > tHigh, add the pixel and then visit it's neighbors
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
               
        float[] tHs = new float[2];
        float[] tLs = new float[2];
        if (!useAdaptive2Layer || !useAdaptiveThreshold) {
            OtsuThresholding ot = new OtsuThresholding();
            for (int i = 0; i < 2; ++i) {
                tHs[i] = otsuScaleFactor * ot.calculateBinaryThreshold256(
                    lcFilterProducts[i].getGradientXY());
                tLs[i] = tHs[i]/factorBelowHighThreshold;
            }
        }
       
        GreyscaleImage combXY = lcFilterProducts[0]
            .getGradientX().createWithDimensions();
                
        // store pixels w/ v > tHigh
        // and store any of it's neighbors w/ v > tLow
        
        GreyscaleImage[] gs = new GreyscaleImage[]{
            lcFilterProducts[0].getGradientXY(),
            lcFilterProducts[1].getGradientXY()};

        double[] tHs0 = new double[2];
        double[] tLs0 = new double[2];
        
        if (lcThreshImg == null) {
            for (int i = 0; i < 2; ++i) {
                tHs0[i] = tHs[i];
                tLs0[i] = tLs[i];
            }
        }
        
        double[] vs = new double[2];
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {

                for (int i = 0; i < 2; ++i) {
                    vs[i] = gs[i].getValue(x, y);
                }

                if (lcThreshImg != null) {
                    for (int i = 0; i < 2; ++i) {
                        tHs0[i] = lcThreshImg[i][x][y];
                        tLs0[i] = tHs0[i]/2;
                    }
                }

                // both "L" and "C" gradient pixels are below high threshold
                if (vs[0] < tHs0[0] && vs[1] < tHs0[1]) {
                    continue;
                }

                combXY.setValue(x, y, 255);

                // store any adjacent w/ v > tLow
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }
                    for (int i = 0; i < 2; ++i) {
                        if (gs[i].getValue(x2, y2) > tLs0[i]) {
                            combXY.setValue(x2, y2, 255);
                            break;
                        } 
                    }
                }
            }
        }
        
        applyPostLineThinningCorrections(combXY);
        
        for (int i = 0; i < 2; ++i) {            
            gs[i].resetTo(combXY);
        }

        if (debug) {
            MiscDebug.writeImage(combXY, "_after_linethinning_1_");
        }
        
        if (restoreJunctions) {
        
            // while have information about the thresholds, will use them to
            // make decisions about restoring pixels that disconnected edges.
            PairInt[][] neighborCoordOffsets
                = AbstractLineThinner.createCoordinatePointsForEightNeighbors(
                    0, 0);

            for (int i0 = 0; i0 < 2; ++i0) {
            
                Set<PairInt> removedDisconnecting = lcRemovedDisconnecting.get(i0);
                
                for (PairInt p : removedDisconnecting) {

                    if (ImageSegmentation.doesDisconnect(combXY,
                        neighborCoordOffsets, p.getX(), p.getY())) {

                        int x = p.getX();
                        int y = p.getY();

                        if (lcThreshImg != null) {
                            for (int ii = 0; ii < 2; ++ii) {
                                tHs0[ii] = lcThreshImg[ii][x][y];
                                tLs0[ii] = tHs0[ii] / 2;
                            }
                        }

                        int v = gs[i0].getValue(p);

                        if (v > tHs0[i0]) {
                            if (isAdjacentToAHorizOrVertLine(combXY, x, y, 3)) {
                                combXY.setValue(x, y, 255);
                            }
                        }
                    }
                }
            }
     
            if (debug) {
                MiscDebug.writeImage(combXY, "_after_restore_junctions_");
            }
            
            applyPostLineThinningCorrections(combXY);
            
            if (debug) {
                MiscDebug.writeImage(combXY, "_after_linethinning_2_");
            }
        }
        
        GreyscaleImage combX = lcFilterProducts[0].getGradientX()
            .copyImage();
        GreyscaleImage combY = lcFilterProducts[0].getGradientY()
            .copyImage();

        // need to make a combined gradient in X and Y, then calc theta
        //   image, based upon combXY results.
        int nBins = 40;
        int binSz = (int)Math.ceil(512./(double)nBins);
        for (int i = 0; i < 2; ++i) {
            // plotting histograms to see how to make
            int[] hX = new int[nBins];
            GreyscaleImage gX = lcFilterProducts[i].getGradientX();
            for (int pixIdx = 0; pixIdx < n; ++pixIdx) {
                int v = gX.getValue(pixIdx);
                int bin = (v - -255)/binSz;
                hX[bin]++;
            }
            int[] hY = new int[nBins];
            GreyscaleImage gY = lcFilterProducts[i].getGradientY();
            for (int pixIdx = 0; pixIdx < n; ++pixIdx) {
                int v = gY.getValue(pixIdx);
                int bin = (v - -255)/binSz;
                hY[bin]++;
            }
            System.out.println(i + ": \n  hX=" 
                + Arrays.toString(hX) + " \n  hY=" + Arrays.toString(hY)
                + "\n  " + Arrays.toString(tHs));
            
            //Need to look at this in more detail to decide whether the
            //  addition of 2 canny edges of lch[0] and lch[1] added
            //     together are better than this pixel by pixel
            //     2 layer filter approach.
            //     the approach in this class currently is not as good for
            //     low contrast images.
        }
        
        // usually, the color gradient seems to be smaller values,
        // so until see the histograms and impl scaling logic,
        // will use tHs
        float cFactor = tHs[0]/tHs[1];
        
        for (int pixIdx = 0; pixIdx < n; ++pixIdx) {
            int v1 = combX.getValue(pixIdx);
            int v2 = Math.round(cFactor * lcFilterProducts[1]
                .getGradientX().getValue(pixIdx));
            
            // instead of averaging, will take the strongest
            if (v2 > v1) {
                if (v2 > 255) {
                    combX.setValue(pixIdx, 255);
                } else {
                    combX.setValue(pixIdx, v2);
                }
            }
            
            // --- gradient Y ---
            
            v1 = combY.getValue(pixIdx);
            v2 = Math.round(cFactor * lcFilterProducts[1]
                .getGradientY().getValue(pixIdx));
            
            // instead of averaging, will take the strongest
            if (v2 > v1) {
                if (v2 > 255) {
                    combY.setValue(pixIdx, 255);
                } else {
                    combY.setValue(pixIdx, v2);
                }
            }
        }

        GreyscaleImage combTheta = 
            imageProcessor.computeTheta180(combX, combY);
        
        EdgeFilterProducts efp = new EdgeFilterProducts();
        efp.setGradientX(combX);
        efp.setGradientY(combY);
        efp.setGradientXY(combXY);
        efp.setTheta(combTheta);
        
        return efp;
    }
    
    /**
     * convolve the image with a Sobel X 1D kernel which is the same as a 
     * Gaussian first derivative with sigma = sqrt(1)/2.
     * FWHM=2.355*sigma
     * 
     * @param input
     * @return 
     */
    private GreyscaleImage createGradientX1D(final GreyscaleImage input) {
                
        return createGradient1D(input, true);
    }
    
    /**
     * convolve the image with a Sobel Y 1D kernel which is the same as a 
     * Gaussian first derivative with sigma = sqrt(1)/2.
     * 
     * @param input
     * @return 
     */
    private GreyscaleImage getGradientY1D(final GreyscaleImage input) {
                
        return createGradient1D(input, false);
    }
    
    /**
     * convolve the image with a Sobel 1D kernel which is the same as a 
     * Gaussian first derivative with sigma = sqrt(1)/2.
     * 
     * @param input
     * @return 
     */
    private GreyscaleImage createGradient1D(final GreyscaleImage input, 
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
     * The theta image has range 0 to 180.
     * 
     * @param lch
     * @return 
     */
    protected EdgeFilterProducts[] createGradients(final GreyscaleImage[] lch) {
    
        GreyscaleImage gX, gY, g, theta;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        EdgeFilterProducts[] eps = new EdgeFilterProducts[2];
        
        for (int i = 0; i < 2; ++i) {
        
            gX = createGradientX1D(lch[i]);

            gY = getGradientY1D(lch[i]);
        
            g = imageProcessor.combineConvolvedImages(gX, gY);
        
            // the theta is in range 0 to 180
            theta = imageProcessor.computeTheta180(gX, gY);
        
            if (debug) {
                // when replace the aspect library, put these renders in the 
                //   equivalent replacement
                long ts = MiscDebug.getCurrentTimeFormatted();
                MiscDebug.writeImage(theta, "_theta_" + i + "_" + ts);
                MiscDebug.writeImage(gX, "_gX_" + i + "_" + ts);
                MiscDebug.writeImage(gY, "_gY_" + i + "_" + ts);
            }
        
            EdgeFilterProducts efp = new EdgeFilterProducts();
            efp.setGradientX(gX);
            efp.setGradientY(gY);
            efp.setGradientXY(g);
            efp.setTheta(theta);
            
            eps[i] = efp;
        }
        
        return eps;
    }

    /**
     * note, uses small gaussian, so should only be used when very little or
     * no smoothing has occurred (line drawings, for example might be best
     * handled with this and no smoothing).
     * @param lch
     * @return 
     */
    protected EdgeFilterProducts[] createDiffOfGaussians(final 
        GreyscaleImage[] lch) {
        
        EdgeFilterProducts[] eps = new EdgeFilterProducts[2];
        
        GreyscaleImage gX, gY, g, theta;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        for (int i = 0; i < 2; ++i) {
            
            gX = createGradientFromDiffOfGauss(lch[i], true);

            gY = createGradientFromDiffOfGauss(lch[i], false);

            g = imageProcessor.combineConvolvedImages(gX, gY);

            // the theta is in range 0 to 180
            theta = imageProcessor.computeTheta180(gX, gY);

            EdgeFilterProducts efp = new EdgeFilterProducts();
            efp.setGradientX(gX);
            efp.setGradientY(gY);
            efp.setGradientXY(g);
            efp.setTheta(theta);
            
            eps[i] = efp;
        }
        
        return eps;
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

    private void applyNonMaximumSuppression(EdgeFilterProducts filterProducts,
        Set<PairInt> disconnectingRemovals) {
        
         NonMaximumSuppression nms = new NonMaximumSuppression();
         
         //TODO: radius can be adjusted from 1 to higher
         nms.nonmaxsup(filterProducts.getGradientXY(),
            filterProducts.getTheta(), 1.2, disconnectingRemovals);
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

    /**
     * get the filter products for gradient and orientation.
     * note that the orientation image has values between 0 and 180.
     * @return the filterProducts
     */
    public EdgeFilterProducts getFilterProducts() {
        return filterProducts;
    }

    private void applyPostLineThinningCorrections(GreyscaleImage gradientXY) {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.applyThinning(gradientXY, false);
        
        if (useLineThinner) {
            ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
            lt.applyFilter(gradientXY);
        }        
    }

    private int countAboveZero(GreyscaleImage img) {
        
        int n = 0;
        for (int i = 0; i < img.getNPixels(); ++i) {
            if (img.getValue(i) > 0) {
                n++;
            }
        }
        
        return n;
    }

}
