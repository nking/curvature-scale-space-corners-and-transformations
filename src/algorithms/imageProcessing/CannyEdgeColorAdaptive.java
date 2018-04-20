package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import java.util.Set;
import algorithms.misc.Misc;
import java.util.HashSet;
import java.util.logging.Logger;

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
 * 
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
            
        ImageProcessor imageProcessor = new ImageProcessor();
        
        GreyscaleImage[] lch = imageProcessor.createLCHForLUV(input);
        
        CannyEdgeFilterAdaptive2 cannyL = new CannyEdgeFilterAdaptive2();
        CannyEdgeFilterAdaptive2 cannyC = new CannyEdgeFilterAdaptive2();
        if (!useLineThinner) {
            cannyL.overrideToNotUseLineThinner();
            cannyC.overrideToNotUseLineThinner();
        }
        if (useAdaptiveThreshold) {
            cannyL.overrideToUseAdaptiveThreshold();
            cannyC.overrideToUseAdaptiveThreshold();
        }
        
        cannyL.override2LayerFactorBelowHighThreshold(factorBelowHighThreshold);
        cannyC.override2LayerFactorBelowHighThreshold(factorBelowHighThreshold);
        
        if (!performNonMaxSuppr) {
            cannyL.setToNotUseNonMaximumSuppression();
            cannyC.setToNotUseNonMaximumSuppression();
        }
        
        if (!restoreJunctions) {
            cannyL.setToNotRestoreJunctions();
            cannyC.setToNotRestoreJunctions();
        }
        
        if (!useAdaptive2Layer) {
            cannyL.setToUseSingleThresholdIn2LayerFilter();
            cannyC.setToUseSingleThresholdIn2LayerFilter();
        }
        
        if (!useLineThinner) {
            cannyL.overrideToNotUseLineThinner();
            cannyC.overrideToNotUseLineThinner();
        }
        
        if (debug) {
            cannyL.setToDebug();
            cannyC.setToDebug();
        }
        
        cannyL.setOtsuScaleFactor(otsuScaleFactor);
        cannyC.setOtsuScaleFactor(otsuScaleFactor);
        //cannyC.setOtsuScaleFactor(1.0f);  // less sensitive and less noisey
        
        cannyL.applyFilter(lch[0]);
        cannyC.applyFilter(lch[1]);
        
        EdgeFilterProducts edgeProductsL = cannyL.getFilterProducts();
        EdgeFilterProducts edgeProductsC = cannyC.getFilterProducts();
        
        // DEBUG: temporary look at recalculating the L thresholds
        //        to filter out scaled C values to reduce noise.
        //        assuming not adaptive for now.
        int tLowL = edgeProductsL.getGradientXY().min();
        
        float lFactor = 255.f/(float)edgeProductsL.getGradientXY().max();
        float cFactor = 255.f/(float)edgeProductsC.getGradientXY().max();
        
        GreyscaleImage combXY = edgeProductsL.getGradientXY()
            .createWithDimensions();
        
        GreyscaleImage combX = edgeProductsL.getGradientX()
            .createWithDimensions();
        
        GreyscaleImage combY = edgeProductsL.getGradientY()
            .createWithDimensions();
        
        int n = combXY.getNPixels();
        
        int v0, v1, v, vx, vy;
        for (int i = 0; i < n; ++i) {
            v0 = edgeProductsL.getGradientXY().getValue(i);
            v1 = edgeProductsC.getGradientXY().getValue(i);
            
            v0 = Math.round(v0 * lFactor);
            if (v0 > 255) {
                v0 = 255;
            }
            v1 = Math.round(v1 * cFactor);
            if (v1 > 255) {
                v1 = 255;
            }
            
            if (cFactor > 1) {
                if (v1 < tLowL) {
                    v1 = 0;
                }
            }
                        
            // choosing the largest of both instead of avg
            if (v0 > v1) {
                v = v0;
                vx = edgeProductsL.getGradientX().getValue(i);
                vy = edgeProductsL.getGradientY().getValue(i);
            } else {
                v = v1;
                vx = edgeProductsC.getGradientX().getValue(i);
                vy = edgeProductsC.getGradientY().getValue(i);
            }
            combXY.setValue(i, v);
            
            combX.setValue(i, vx);
            combY.setValue(i, vy);
        }
        
        GreyscaleImage combTheta = 
            imageProcessor.computeTheta180(combX, combY);
        
        EdgeFilterProducts efp = new EdgeFilterProducts();
        efp.setGradientX(combX);
        efp.setGradientY(combY);
        efp.setGradientXY(combXY);
        efp.setTheta(combTheta);
         
        this.filterProducts = efp;
    }
   
    /**
     * get the filter products for gradient and orientation.
     * note that the orientation image has values between 0 and 180.
     * @return the filterProducts
     */
    public EdgeFilterProducts getFilterProducts() {
        return filterProducts;
    }

    private class CannyEdgeFilterAdaptive2 {

        private float factorBelowHighThreshold = 2.f;

        private EdgeFilterProducts filterProducts = null;

        private boolean performNonMaxSuppr = true;

        private boolean debug = false;
        
        private boolean restoreJunctions = true;
 
        private boolean useHigherThresholdIfNeeded = false;

        private double approxProcessedSigma = 0;

        private boolean useAdaptiveThreshold = false;

        private boolean useAdaptive2Layer = true;

        private boolean lineDrawingMode = false;

        private float otsuScaleFactor = 0.75f;//0.65f;

        private Logger log = Logger.getLogger(this.getClass().getName());

        private boolean useLineThinner = true;
        
        public CannyEdgeFilterAdaptive2() {
        }

        public void setToDebug() {
            debug = true;
        }

        public void overrideToNotUseLineThinner() {
            useLineThinner = false;
        }

        public void setToNotRestoreJunctions() {
            restoreJunctions = false;
        }

        public void setToUseLineDrawingMode() {
            lineDrawingMode = true;
        }
        
        public void setOtsuScaleFactor(float factor) {
            otsuScaleFactor = factor;
        }

        public void setToNotUseNonMaximumSuppression() {
            performNonMaxSuppr = false;
        }

        public void setToUseHigherThresholdIfNeeded() {
            useHigherThresholdIfNeeded = true;
        }

        public void overrideToUseAdaptiveThreshold() {
            useAdaptiveThreshold = true;
        }
        
        public void setToUseSingleThresholdIn2LayerFilter() {
            useAdaptive2Layer = false;
        }

        /**
         * override the default factor of low threshold below high threshold,
         * which is 2.
         *
         * @param factor
         */
        public void override2LayerFactorBelowHighThreshold(float factor) {
            factorBelowHighThreshold = factor;
        }

        public void applyFilter(final GreyscaleImage input) {

            if (input.getWidth() < 3 || input.getHeight() < 3) {
                throw new IllegalArgumentException("images should be >= 3x3 in size");
            }

            if (lineDrawingMode) {
                useAdaptive2Layer = true;
                useAdaptiveThreshold = false;
                //apply2LayerFilter(input, new HashSet<PairInt>(), input);
                if (debug) {
                    GreyscaleImage imgcp = input.copyImage();
                    for (int i = 0; i < imgcp.getNPixels(); ++i) {
                        int v = imgcp.getValue(i);
                        if (v > 0) {
                            imgcp.setValue(i, 255);
                        }
                    }
                    MiscDebug.writeImage(imgcp, "_after_2_layer_");
                }

                ImageProcessor imageProcessor = new ImageProcessor();

                GreyscaleImage in = imageProcessor.closing(input);
                imageProcessor.blur(in, SIGMA.ZEROPOINTFIVE);

                // convert all > 0 to 1's for 'closing' operation
                for (int pixIdx = 0; pixIdx < in.getNPixels(); ++pixIdx) {
                    if (in.getValue(pixIdx) > 0) {
                        in.setValue(pixIdx, 1);
                    }
                }

                in = imageProcessor.closing(in);
                
                {//DEBUG
                    long ts = MiscDebug.getCurrentTimeFormatted();
                    GreyscaleImage tmp = in.copyImage();
                    tmp.multiply(255);
                    MiscDebug.writeImage(tmp, "_line_" + ts);
                }

                int nnzs = imageProcessor.countNonZeroes(in);
                while (true) {
                    int nnzs2 = imageProcessor.applyThinning2(in);
                    if (nnzs == nnzs2) {
                        break;
                    }
                    nnzs = nnzs2;
                }

                {//DEBUG
                    long ts = MiscDebug.getCurrentTimeFormatted();
                    GreyscaleImage tmp = in.copyImage();
                    tmp.multiply(255);
                    MiscDebug.writeImage(tmp, "_line_thinned_" + ts);
                }
                
                // either convert all > 0 back to original values, or to a large
                //   enough value to make a strong gradient
                for (int pixIdx = 0; pixIdx < in.getNPixels(); ++pixIdx) {
                    if (in.getValue(pixIdx) > 0) {
                        in.setValue(pixIdx, 200);
                    }
                }

                filterProducts = createGradientForLineMode(in);

                //TODO: update these or remove them...left from
                // a previous version of the code
                approxProcessedSigma = Math.sqrt(
                    approxProcessedSigma * approxProcessedSigma + (0.678 * 0.678));

                input.resetTo(filterProducts.getGradientXY());

                return;
            }

            // (1) smooth image using separable binomial filters
            SIGMA sigma = SIGMA.ONE;
            if (sigma.equals(SIGMA.ONE)) {
                ATrousWaveletTransform at = new ATrousWaveletTransform();
                GreyscaleImage smoothed = at.smoothToLevel1B3Spline(input);
                input.resetTo(smoothed);
                approxProcessedSigma = 1;
            } else if (sigma.equals(SIGMA.ZEROPOINTSEVENONE)) {
                ATrousWaveletTransform at = new ATrousWaveletTransform();
                GreyscaleImage smoothed = at.smoothFirstLevelTriangle(input);
                input.resetTo(smoothed);
                approxProcessedSigma = Math.sqrt(2.) / 2.;
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
                approxProcessedSigma * approxProcessedSigma + (1. / 4.));

            Set<PairInt> removedDisconnecting = new HashSet<PairInt>();

            if (debug) {
                MiscDebug.writeImage(filterProducts.getGradientXY(), "_before_nms_");
            }
            
            //(3) non-maximum suppression
            if (performNonMaxSuppr) {
                applyNonMaximumSuppression(filterProducts, removedDisconnecting);
            }

            if (debug) {
                MiscDebug.writeImage(filterProducts.getGradientXY(), "_after_nms_");
            }

            //(4) adaptive 2 layer filter
            apply2LayerFilter(filterProducts.getGradientXY(), removedDisconnecting,
                gradientCopyBeforeThinning);
            
            if (restoreJunctions) {
                int minResolution = (int) Math.ceil(2.35 * approxProcessedSigma);
                int minResolvableAngle = (int) Math.ceil(
                    Math.atan2(1, minResolution) * 180. / Math.PI);
                if (minResolvableAngle < 0) {
                    minResolvableAngle *= -1;
                }

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

            input.resetTo(filterProducts.getGradientXY());
        }
        
        protected void apply2LayerFilter(final GreyscaleImage gradientXY,
            Set<PairInt> removedDisconnecting,
            GreyscaleImage gradientCopyBeforeThinning) {

            int n = gradientXY.getNPixels();
            int w = gradientXY.getWidth();
            int h = gradientXY.getHeight();

            if (w < 3 || h < 3) {
                throw new IllegalArgumentException("images should be >= 3x3 in size");
            }

            ImageProcessor imageProcessor = new ImageProcessor();

            double[][] threshImg = null;
            
            if (useAdaptive2Layer && useAdaptiveThreshold) {
                AdaptiveThresholding th = new AdaptiveThresholding();
                threshImg = th.createAdaptiveThresholdImage(
                    imageProcessor.copy(gradientXY), 15, 0.2);
                if (debug) {//DEBUG
                    double[][] imgCp = imageProcessor.copy(gradientXY);
                    for (int i = 0; i < w; ++i) {
                        for (int j = 0; j < h; ++j) {
                            double t = threshImg[i][j];
                            if (imgCp[i][j] > t) {
                                imgCp[i][j] = 255.;
                            } else {
                                imgCp[i][j] = 0;
                            }
                        }
                    }
                    MiscDebug.writeImage(imgCp, "img_a_thresholded_.png");
                    MiscDebug.writeImage(threshImg, "img_a_adaptive_threshold_.png");
                }
            }

            int[] dxs = Misc.dx8;
            int[] dys = Misc.dy8;
            
            float tHigh = 0;
            float tLow = 0;
            if (!useAdaptive2Layer || !useAdaptiveThreshold) {
                OtsuThresholding ot = new OtsuThresholding();

                double[][] g = new double[w][];
                for (int i = 0; i < w; ++i) {
                    g[i] = new double[h];
                    for (int j = 0; j < h; ++j) {
                        g[i][j] = gradientXY.getValue(i, j);
                    }
                }
                int nBins = 256 / 5;
                float t = (float) ot.calculateBinaryThreshold2D(g, nBins);

                tHigh = otsuScaleFactor * t;
                tLow = tHigh / factorBelowHighThreshold;
            }

            GreyscaleImage img2 = gradientXY.createWithDimensions();

            // store pixels w/ v > tHigh
            // and store any of it's neighbors w/ v > tLow
            
            for (int x = 0; x < w; ++x) {
                for (int y = 0; y < h; ++y) {

                    double v = gradientXY.getValue(x, y);

                    double tHigh0, tLow0;
                    if (threshImg != null) {
                        tHigh0 = threshImg[x][y];
                        tLow0 = tHigh0 / 2;
                    } else {
                        tHigh0 = tHigh;
                        tLow0 = tLow;
                    }

                    if (v < tHigh0) {
                        continue;
                    }

                    img2.setValue(x, y, 255);
                    
                    // store any adjacent w/ v > tLow
                    for (int k = 0; k < dxs.length; ++k) {
                        int x2 = x + dxs[k];
                        int y2 = y + dys[k];
                        if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                            continue;
                        }
                        double v2 = gradientXY.getValue(x2, y2);
                        if (v2 > tLow0) {
                            img2.setValue(x2, y2, 255);
                        }
                    }
                }
            }

            if (debug) {
                MiscDebug.writeImage(gradientXY, "_before_2layer_");
                MiscDebug.writeImage(img2, "_in_2layer_");
            }

            gradientXY.resetTo(img2);

            if (debug) {
                MiscDebug.writeImage(gradientXY, "_after_linethinning_1_");
            }
             
            if (restoreJunctions) {

                // while have information about the thresholds, will use them to
                // make decisions about restoring pixels that disconnected edges.
                PairInt[][] neighborCoordOffsets
                    = AbstractLineThinner.createCoordinatePointsForEightNeighbors(
                        0, 0);

                for (PairInt p : removedDisconnecting) {

                    if (ImageSegmentation.doesDisconnect(gradientXY,
                        neighborCoordOffsets, p.getX(), p.getY())) {

                        int x = p.getX();
                        int y = p.getY();
                        int i = gradientXY.getIndex(x, y);

                        double tHigh0, tLow0;
                        if (threshImg != null) {
                            tHigh0 = threshImg[x][y];
                            tLow0 = tHigh0 / 2;
                        } else {
                            tHigh0 = tHigh;
                            tLow0 = tLow;
                        }

                        int v = gradientCopyBeforeThinning.getValue(p);

                        if (v > tLow0) {
                            if (isAdjacentToAHorizOrVertLine(gradientXY, x, y, 3)) {
                                gradientXY.setValue(x, y, 255);
                            }
                        }
                    }
                }

                if (debug) {
                    MiscDebug.writeImage(gradientXY, "_after_restore_junctions_");
                }

                if (debug) {
                    MiscDebug.writeImage(gradientXY, "_after_linethinning_2_");
                }
            }
        }

        private GreyscaleImage getGradientX1D(final GreyscaleImage input) {

            return getGradient1D(input, true);
        }

        private GreyscaleImage getGradientY1D(final GreyscaleImage input) {

            return getGradient1D(input, false);
        }
        
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
        
        protected EdgeFilterProducts createGradient(final GreyscaleImage img) {

            GreyscaleImage g, theta;

            ImageProcessor imageProcessor = new ImageProcessor();

            GreyscaleImage[] gXgY = imageProcessor.createCentralDifferenceGradients(img);

            g = imageProcessor.combineConvolvedImages(gXgY[0], gXgY[1]);

            // the theta is in range 0 to 180
            theta = imageProcessor.computeTheta180(gXgY[0], gXgY[1]);

            if (debug) {
                // when replace the aspect library, put these renders in the
                //   equivalent replacement
                long ts = MiscDebug.getCurrentTimeFormatted();
                MiscDebug.writeImage(theta, "_theta_" + ts);
                MiscDebug.writeImage(gXgY[0], "_gX_" + ts);
                MiscDebug.writeImage(gXgY[1], "_gY_" + ts);
                /*
            int x = 37; int y = 163;
            System.out.println("(" + x + ", " + y + ") math.atan2(" + gY.getValue(x, y)
                + "," + gX.getValue(x, y) + ")*180./math.pi=" +
                theta.getValue(x, y));*/
            }

            EdgeFilterProducts efp = new EdgeFilterProducts();
            efp.setGradientX(gXgY[0]);
            efp.setGradientY(gXgY[1]);
            efp.setGradientXY(g);
            efp.setTheta(theta);

            return efp;
        }
        
        protected EdgeFilterProducts createGradientForLineMode(final GreyscaleImage img) {

            GreyscaleImage g, theta;

            ImageProcessor imageProcessor = new ImageProcessor();

            GreyscaleImage[] gXgY = imageProcessor.createCentralDifferenceGradients(img);

            g = imageProcessor.combineConvolvedImages(gXgY[0], gXgY[1]);

            // the theta is in range 0 to 180
            theta = imageProcessor.computeTheta180(gXgY[0], gXgY[1]);

            if (debug) {
                // when replace the aspect library, put these renders in the
                //   equivalent replacement
                long ts = MiscDebug.getCurrentTimeFormatted();
                GreyscaleImage tmp = g.copyImage();
                MiscDebug.writeImage(tmp, "_gXY_" + ts);
                MiscDebug.writeImage(theta, "_theta_" + ts);
                MiscDebug.writeImage(gXgY[0], "_gX_" + ts);
                MiscDebug.writeImage(gXgY[1], "_gY_" + ts);
                /*
            int x = 37; int y = 163;
            System.out.println("(" + x + ", " + y + ") math.atan2(" + gY.getValue(x, y)
                + "," + gX.getValue(x, y) + ")*180./math.pi=" +
                theta.getValue(x, y));*/
            }

            EdgeFilterProducts efp = new EdgeFilterProducts();
            efp.setGradientX(gXgY[0]);
            efp.setGradientY(gXgY[1]);
            efp.setGradientXY(g);
            efp.setTheta(theta);

            return efp;
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
                for (int dy2 = (-minLineSize / 2); dy2 <= (minLineSize / 2); ++dy2) {
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
                for (int dx2 = (-minLineSize / 2); dx2 <= (minLineSize / 2); ++dx2) {
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
        
        public EdgeFilterProducts getFilterProducts() {
            return filterProducts;
        }

        private void applyPostLineThinningCorrections(GreyscaleImage gradientXY) {

            ImageProcessor imageProcessor = new ImageProcessor();
            GreyscaleImage tmp = gradientXY.copyImage();
            for (int i = 0; i < tmp.getNPixels(); ++i) {
                if (tmp.getValue(i) > 0) {
                    tmp.setValue(i, 1);
                }
            }
            imageProcessor.applyThinning2(tmp);

            for (int i = 0; i < tmp.getNPixels(); ++i) {
                if (tmp.getValue(i) == 0) {
                    gradientXY.setValue(i, 0);
                }
            }

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
 }
