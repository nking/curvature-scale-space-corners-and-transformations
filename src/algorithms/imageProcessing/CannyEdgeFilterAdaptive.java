package algorithms.imageProcessing;

import algorithms.misc.Misc;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * The CannyEdge filter is an algorithm to operate on an image to
 * replace objects with their edges.  
 * 
 * The class began by following the general advice given in
 * "Performance Analysis of Adaptive Canny Edge Detector 
 * Using Bilateral Filter" by Rashmi, Kumar, Jaiswal, and Saxena, but made some
 * modifications.
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
 * By default, the number of levels is 1.
 * To use the filter adaptive mode, use filter.overrideDefaultNumberOfLevels
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
 * </pre>
 * Not ready for use yet... may improve the line thinning...
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

    /*
    for gradient, default is sobel kernel, but can override to use a small
    difference of gaussian.
    */
    private boolean useDiffOfGauss = false;
        
    private float otsuScaleFactor = 0.45f;//0.4f;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public CannyEdgeFilterAdaptive() {        
    }
    
    public void setToUseDiffOfGauss() {
        
        //TODO: consider dilate before erosion
        
        this.useDiffOfGauss = true;        
    }
    
    /**
     * applies a histogram equalization before any processing in order to rescale
     * the data to use the entire range of 0 to 255.
     */
    public void setToPerformHistogramEqualization() {
        performHistEq = true;
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
        
        // (1) smooth image using separable binomial filters
        SIGMA sigma = SIGMA.ONE;
        if (sigma.equals(SIGMA.ONE)) {
            ATrousWaveletTransform at = new ATrousWaveletTransform();
            GreyscaleImage smoothed = at.smoothToSigmaOne(input);
            input.resetTo(smoothed);
        } else if (sigma.equals(SIGMA.ZEROPOINTSEVENONE)) {
            ATrousWaveletTransform at = new ATrousWaveletTransform();
            GreyscaleImage smoothed = at.smoothToSigmaZeroPointSevenOne(input);
            input.resetTo(smoothed);
        } else {
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.blur(input, sigma, 0, 255);
        }
        
        //(2) create gradient
        EdgeFilterProducts filterProducts;
        if (useDiffOfGauss) {
            filterProducts = createDiffOfGaussians(input);
        } else {
            filterProducts = createGradient(input);
        }
        
        //(3) non-maximum suppression
        if (performNonMaxSuppr) {
            applyNonMaximumSuppression(filterProducts);
        }
        
        //(4) adaptive 2 layer filter                        
        apply2LayerFilter(filterProducts);

        input.resetTo(filterProducts.getGradientXY());
        
        // is this necessary?
        for (int i = 0; i < input.getNPixels(); ++i) {
            int v = input.getValue(i);
            if (v < 0) {
                input.setValue(i, 0);
            }
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
    protected void apply2LayerFilter(final EdgeFilterProducts edgeProducts) {
        
        GreyscaleImage input = edgeProducts.getGradientXY();
        int n = input.getNPixels();
        int w = input.getWidth();
        int h = input.getHeight();
        
        if (w < 3 || h < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }
        
        int nLevels = numberOfLevelsForHistogram;
        
        OtsuThresholding ot = new OtsuThresholding();
        
        int[] pixelThresholds = null;
        if (useAdaptive2Layer && (nLevels > 1)) {
            pixelThresholds = ot.calculateMultiLevelThreshold256(input, nLevels);
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
            tHigh = otsuScaleFactor * ot.calculateBinaryThreshold256(input);
            tLow = tHigh/factorBelowHighThreshold;
        }
        
        GreyscaleImage img2 = input.createWithDimensions();
        
        for (int i = 0; i < input.getNPixels(); ++i) {
            
            if (useAdaptive2Layer && (nLevels > 1)) {
                
                tHigh = otsuScaleFactor * pixelThresholds[i];
            
                tLow = tHigh/factorBelowHighThreshold;
            }
            
            int v = input.getValue(i);
            
            if (v < tLow) {
                continue;
            } else if (v > tHigh) {
                img2.setValue(i, v);
                continue;
            }
            
            int x = input.getCol(i);
            int y = input.getRow(i);
            
            boolean foundHigh = false;
            boolean foundMid = false;
            
            for (int k = 0; k < dxs.length; ++k) {                
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                    continue;
                }
                int v2 = input.getValue(x2, y2);
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
                    int v2 = input.getValue(x2, y2);
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

        input.resetTo(img2);
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
        
        EdgeFilterProducts efp = new EdgeFilterProducts();
        efp.setGradientX(gX);
        efp.setGradientY(gY);
        efp.setGradientXY(g);
        efp.setTheta(theta);
        
        return efp;
    }

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

    private void applyNonMaximumSuppression(EdgeFilterProducts filterProducts) {
        
        /*
        TODO:
        in order to reduce effects of double lines due to blurring before and 
        during edge creation,
        may want to copy the gradient image, then apply nms to the gradient,
        then correct for double lines and then replace lost junctions using
        the copied gradient image for a reference.
        */
        
        GreyscaleImage theta = filterProducts.getTheta();
        
        GreyscaleImage img = filterProducts.getGradientXY();
        
        int n = img.getNPixels();
        
        int w = img.getWidth();
        int h = img.getHeight();
                
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
        for (int i = 0; i < n; ++i) {
            
            int x = theta.getCol(i);
            int y = theta.getRow(i);
            if ((x == 0) || (y == 0) || (x == (w - 1)) || (y == (h - 1))) {
                continue;
            }
            
            int t = theta.getValue(i);
         
            int compDx1, compDy1, compDx2, compDy2;

            if ((t < 23) || (t > 337) || ((t > 157) && (t < 203))) {
                // in range of +0 or +180
                compDx1 = 1;
                compDy1 = 0;
                compDx2 = -1;
                compDy2 = 0;
            } else if (((t > 22) && (t < 68)) || ((t > 202) && (t < 248))) {
                // in range of +45 or +225
                compDx1 = 1;
                compDy1 = 1;
                compDx2 = -1;
                compDy2 = -1;
            } else if (((t > 67) && (t < 113)) || ((t > 247) && (t < 293))) {
                // in range of +90 or +270
                compDx1 = 1;
                compDy1 = 1;
                compDx2 = -1;
                compDy2 = -1;
            } else { //if (((t > 112) && (t < 158)) || ((t > 292) && (t < 338))) {
                // in range of +135 or +315
                compDx1 = -1;
                compDy1 = 1;
                compDx2 = 1;
                compDy2 = -1;
            }
            
            int v = img.getValue(x, y);
            
            int v1 = img.getValue(x + compDx1, y + compDy1);
            int v2 = img.getValue(x + compDx2, y + compDy2);
            
            if ((v < v1) || (v < v2)) {
                img.setValue(x, y, 0);
            }
        }
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        curveHelper.additionalThinning45DegreeEdges2(theta, img);
    }
    
}
