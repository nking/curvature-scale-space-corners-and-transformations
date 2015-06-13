package algorithms.imageProcessing;

import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.util.Errors;
import algorithms.misc.MiscMath;
import java.util.logging.Logger;

/**
 * The CannyEdge filter is an algorithm to operate on an image to
 * replace objects with their edges.
 * 
 * The implementation here is slightly different/
 * 
 * The images are first normalized by histogram equalization.
 * 
 * Then the x and y gradients of the image are found.
 * 
 * The combination of the x and y gradients are used as the input for the
 * rest of the algorithm.
 * 
 * "hysteresis thresholding" is performed as a 2 layer filter to remove the 
 * lowest intensity pixels in an image and then to remove pixels not connected 
 * to the highest intensity pixels.
 * 
 * A line thinner is used here instead of "non-maximal" suppression to
 * thin the edges.
 
   An additional stage of line thinning is applied to edges that are at a
   45 degree angle as suggested by the ECSS paper.
 
 * @author nichole
 */
public class CannyEdgeFilter {
        
    /**
     * represents the original image histogram if histogram equalization was
     * not performed, else, represents the image histogram just after
     * histogram equalization
     */
    private HistogramHolder imgHistogram = null;
        
    protected boolean histogramEqualizationCheckFinished = false;
        
    protected boolean doNotNormalizeByHistogram = false;
    
    protected float highThreshold = 5.0f;
    
    public static float defaultLowThreshold = 
        LowIntensityRemovalFilter.defaultLowThreshFactor;
    
    public static float defaultOutdoorLowThreshold = 3.5f;
    
    protected float lowThreshold = defaultLowThreshold;
        
    private boolean useLineDrawingMode = false;
    
    private GreyscaleImage gXY = null;
    
    private GreyscaleImage gX = null;
    
    private GreyscaleImage gY = null;
    
    private GreyscaleImage gTheta = null;
        
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    // ErosionFilter
    private Class<? extends ILineThinner> lineThinnerClass = ZhangSuenLineThinner.class;
    
    protected boolean useOutdoorMode = false;
    
    protected int[] shrinkToSize = null;
    
    public CannyEdgeFilter() {        
    }
 
    /**
     * apply histogram normalization before processing.  For some images, this
     * will increase the contrast of fainter features.
     */
    public void doNotPerformHistogramEqualization() {
        this.doNotNormalizeByHistogram = true;
    }
    
    public void overrideHighThreshold(float thresh) {
        highThreshold = thresh;
    }
    
    public void overrideLowThreshold(float thresh) {
        lowThreshold = thresh;
    }
    
    /**
    Line Drawing Mode uses the  difference of Gaussians instead of 1st deriv Gaussian 
    to try to keep the lines thinner from the start.

        (1) convolve the image with a Gaussian sigma=1 kernel as 1D 
            operations in X and Y separately.
        (2) convolve those with a Gaussian sigma=1 kernel again and
            that results in an image convolved in total by sqrt(2).
        (3) The X gradient is the second convolution - the first, so
            that is effectively a factor of sqrt(2) = 1.414 which is
            close to the recommended 1.6 to mimic a Laplacian filter.
            The same is done for Y.
    */
    public void useLineDrawingMode() {
        
        //TODO: consider dilate before erosion
        
        this.useLineDrawingMode = true;
        
        highThreshold = 2.0f;
    }
       
    public void setLineThinnerClass(Class<? extends ILineThinner> cls) {
        lineThinnerClass = cls;
    }
    
    public void useOutdoorMode() {
        
        useOutdoorMode = true;
        
        highThreshold = 50.0f;
        
        lowThreshold = defaultOutdoorLowThreshold;
    }
    
    public void setFilterImageTrim(int xOffset, int yOffset, int width, int
        height) {
        shrinkToSize = new int[]{xOffset, yOffset, width, height};
    }
    
    public void setSetters(CannyEdgeFilterSettings settings) {
        
        if (settings == null) {
            return;
        }
        
        if (settings.getUseOutdoorMode()) {
            useOutdoorMode = true;
        }
        
        if (settings.getDoNotNormalizeByHistogram()) {
            doNotNormalizeByHistogram = true;
        }
        
        if (settings.getUseLineDrawingMode()) {
            useLineDrawingMode = true;
        }
        
        if (settings.getOverrideHighThreshold()) {
            highThreshold = settings.getHighThreshold();
        }
        
        if (settings.getShrinkToSize() != null) {
            shrinkToSize = settings.getShrinkToSize();
        }
    }
    
    public void applyFilter(final GreyscaleImage input) {
        
        if (input.getWidth() < 3 || input.getHeight() < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        if (shrinkToSize != null) {
            
            imageProcessor.shrinkImage(input, shrinkToSize);
            
        } else {
            
            int[] offsetsXY = imageProcessor.shrinkImageToFirstNonZeros(input);
            
            if (offsetsXY[0] > 0 || offsetsXY[1] > 0) {
                shrinkToSize = new int[]{offsetsXY[0], offsetsXY[1],
                    input.getWidth(), input.getHeight()};
            }
        }
        
        if (useLineDrawingMode) {
            // line drawings should be thinned first
            applyLineThinnerFilter(input);
        }
             
        if (useOutdoorMode) {
            imageProcessor.blur(input, 2.0f); //3.0
        }

        applyHistogramEqualization(input);
        
        //[gX, gY, gXY, theta
        GreyscaleImage[] gradientProducts = createGradientProducts(input);
        
        gX = gradientProducts[0].copyImage();
        
        gY = gradientProducts[1].copyImage();
        
        gXY = gradientProducts[2].copyImage();
        
        gTheta = gradientProducts[3].copyImage();
                
        input.resetTo(gradientProducts[2]);
           
        if (!useLineDrawingMode) {
            apply2LayerFilter(input);
        }
       
        applyLineThinnerFilter(input);
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        curveHelper.additionalThinning45DegreeEdges(gradientProducts[3], input);
    }
    
    /**
     * after applyFilter has already been used, this method can be used after
     * changing the settings such as overrideLowThreshold to re-process the
     * image products with the changed settings.  Note that gradientXYImage
     * is modified in the process, but gradientThetaImage and imgHistogram
     * are not.
     * 
     * @param gradientXYImage
     * @param gradientThetaImage
     * @param imgHistogram 
     */
    public void reApply2LayerFilter(final GreyscaleImage gradientXYImage, 
        final GreyscaleImage gradientThetaImage, HistogramHolder imgHistogram) {
        
        if (gradientXYImage.getWidth() < 3 || gradientXYImage.getHeight() < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }
        
        if (!useLineDrawingMode) {
            apply2LayerFilter(gradientXYImage);
        }
        
        applyLineThinnerFilter(gradientXYImage);
        
        if (lineThinnerClass.equals(ErosionFilter.class)) {
            
            MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
            curveHelper.additionalThinning45DegreeEdges(gradientThetaImage.copyImage(), 
                gradientXYImage);
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
    
    */
    protected void apply2LayerFilter(final GreyscaleImage input) {
        
        if (input.getWidth() < 3 || input.getHeight() < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }
        
        LowIntensityRemovalFilter filter2 = new LowIntensityRemovalFilter();
        
        if (useOutdoorMode) {
            // expecting values 3.5 to 5.0
            filter2.overrideLowThresholdFactor(lowThreshold);
        }
                
        ImageStatistics stats = filter2.removeLowIntensityPixels(input,
            imgHistogram);
        
        double lowThresh = stats.getLowThresholdApplied();
                
        log.info("2-layer filter: low thresh=" + lowThresh);
        
        double threshold2 = lowThresh * highThreshold;
        
        // count number of pixels between lowThresh and threshold2 and
        // above threshold2.  the later should help scale highThreshold
        // factor from 3 to 5 when needed.
        int n0 = ImageStatisticsHelper.countPixels(input, (int)lowThresh, 
            (int)threshold2);
        
        int n1 = ImageStatisticsHelper.countPixels(input,
            (int)threshold2, 255);
        
        double r = (double)n1/(double)n0;
        
        /*
        NOTE: if wanted the final result to include slightly lower intensity
        pixels w/o changing the noise level, these parameters seem to give the
        next best results:
            gaussian sigma = SIGMA.ZEROPOINTSEVONEONE in getGradient1D
            lowThreshFactor = 2.0 in LowIntensityRemovalFilter
            threshold2 *= 0.4 here
        */

        log.info("threshold2=" + threshold2 + " n0=" + n0 + " n1=" + n1 + 
            " n1/n0=" + r);
        
        GreyscaleImage img2 = input.createWithDimensions();
        
        // find points that are "sure-edge" points, above threshold2
        for (int i = 0; i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {
                int pixG = input.getValue(i, j);
                if (pixG >= threshold2) {
                    img2.setValue(i, j, pixG);
                }
            }
        }
        
        int w = img2.getWidth();
        int h = img2.getHeight();
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        
        //TODO:  define "connected" as connected only to "sure-edge" points or
        //       connected to any point in the image in progress?
        // for now, choosing the first and keeping them separate.
        
        GreyscaleImage img3 = input.createWithDimensions();
        
        for (int i = 0; i < input.getWidth(); i++) {
            
            for (int j = 0; j < input.getHeight(); j++) {
                
                int pixG = input.getValue(i, j);
                
                if (img2.getValue(i, j) > 0) {                    
                
                    img3.setValue(i, j, pixG);
                
                } else {
                                        
                    for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                        
                        int x = i + dxs[nIdx];
                        int y = j + dys[nIdx];
                        if ((x<0) || (y<0) || (x>(w-1)) || (y>(h-1))) {
                            continue;
                        }
                        
                        int v2 = img2.getValue(x, y);
                        
                        // to add smaller lower intensity details, can use weak association too:
                        /*int v3 = img3.getValue(x, y);
                        if (v3 > 0) {
                            img3.setValue(i, j, pixG);
                            break;
                        }*/
                        
                        if (v2 > 0) {
                            img3.setValue(i, j, pixG);
                            break;
                        }
                    }
                }
            }
        }

        input.resetTo(img3);
    }

    private GreyscaleImage createGradientXFromDiffOfGauss(
        final GreyscaleImage img) {
        
        return createGradientFromDiffOfGauss(img, true);
    }

    private GreyscaleImage createGradientYFromDiffOfGauss(
        final GreyscaleImage img) {
        
        return createGradientFromDiffOfGauss(img, false);
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
        
        apply1DKernelToImage(g0, kernel, calculateForX);
        
        GreyscaleImage g1 = img.copyImage();
        
        apply1DKernelToImage(g1, kernel2, calculateForX);
       
        ImageProcessor ImageProcessor = new ImageProcessor();
        
        GreyscaleImage g = ImageProcessor.subtractImages(g1, g0);
        
        return g;
    }
    
    private void applyGaussian(final GreyscaleImage img, float sigma) {
       
        float[] kernel = Gaussian1D.getKernel(sigma);
                
        apply1DKernelToImage(img, kernel, true);
                
        apply1DKernelToImage(img, kernel, false);
       
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
                
        GreyscaleImage output = input.copyImage();
        
        apply1DKernelToImage(output, kernel, calculateForX);
        
        return output;
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
    
    void applyLineThinnerFilter(final GreyscaleImage img) {
                
        ILineThinner lineThinner;
        
        try {
            
            lineThinner = lineThinnerClass.newInstance();
            
            if (gXY != null) {
                if (shrinkToSize != null) {
                    GreyscaleImage gXY2 = gXY.copyImage();
                    ImageProcessor imageProcessor = new ImageProcessor();
                    imageProcessor.shrinkImage(gXY2, shrinkToSize);
                    lineThinner.setEdgeGuideImage(gXY2);
                } else {
                    lineThinner.setEdgeGuideImage(gXY);
                }
            }
            
            lineThinner.applyFilter(img);
            
        } catch (InstantiationException ex) {
            throw new IllegalStateException(
                "could not instantiate line thinner: " +
                lineThinnerClass.getSimpleName());
        } catch (IllegalAccessException ex) {
            throw new IllegalStateException(
                "could not instantiate line thinner: " +
                lineThinnerClass.getSimpleName());
        }        
    }

    private void removeOnePixelSpanningBorders(final GreyscaleImage img) {
        
        // remove 1 pixel edges on borders that extend entire length
        boolean foundBoundaryLine = true;
        
        for (int col = 0; col < img.getWidth(); col++) {
            if (img.getValue(col, 0) == 0) {
                foundBoundaryLine = false;
                break;
            }
        }
        if (foundBoundaryLine) {
            for (int col = 0; col < img.getWidth(); col++) {
                img.setValue(col, 0, 0);
            }
        }
        
        for (int col = 0; col < img.getWidth(); col++) {
            if (img.getValue(col, 1) == 0) {
                foundBoundaryLine = false;
                break;
            }
        }
        if (foundBoundaryLine) {
            for (int col = 0; col < img.getWidth(); col++) {
                img.setValue(col, 1, 0);
            }
        }
        
        foundBoundaryLine = true;
        for (int col = 0; col < img.getWidth(); col++) {
            if (img.getValue(col, img.getHeight() - 1) == 0) {
                foundBoundaryLine = false;
                break;
            }
        }
        if (foundBoundaryLine) {
            for (int col = 0; col < img.getWidth(); col++) {
                img.setValue(col, img.getHeight() - 1, 0);
            }
        }
        
        foundBoundaryLine = true;
        for (int col = 0; col < img.getWidth(); col++) {
            if (img.getValue(col, img.getHeight() - 2) == 0) {
                foundBoundaryLine = false;
                break;
            }
        }
        if (foundBoundaryLine) {
            for (int col = 0; col < img.getWidth(); col++) {
                img.setValue(col, img.getHeight() - 2, 0);
            }
        }
        
        foundBoundaryLine = true;
        for (int row = 1; row < (img.getHeight() - 1); row++) {
            if (img.getValue(0, row) == 0) {
                foundBoundaryLine = false;
                break;
            }
        }
        if (foundBoundaryLine) {
            for (int row = 0; row < img.getHeight(); row++) {
                img.setValue(0, row, 0);
            }
        }
        foundBoundaryLine = true;
        for (int row = 1; row < (img.getHeight() - 1); row++) {
            if (img.getValue(1, row) == 0) {
                foundBoundaryLine = false;
                break;
            }
        }
        if (foundBoundaryLine) {
            for (int row = 0; row < img.getHeight(); row++) {
                img.setValue(1, row, 0);
            }
        }
        
        foundBoundaryLine = true;
        for (int row = 1; row < (img.getHeight() - 1); row++) {
            if (img.getValue(img.getWidth() - 1, row) == 0) {
                foundBoundaryLine = false;
                break;
            }
        }
        if (foundBoundaryLine) {
            for (int row = 0; row < img.getHeight(); row++) {
                img.setValue(img.getWidth() - 1, row, 0);
            }
        }
        
        foundBoundaryLine = true;
        for (int row = 1; row < (img.getHeight() - 1); row++) {
            if (img.getValue(img.getWidth() - 2, row) == 0) {
                foundBoundaryLine = false;
                break;
            }
        }
        if (foundBoundaryLine) {
            for (int row = 0; row < img.getHeight(); row++) {
                img.setValue(img.getWidth() - 2, row, 0);
            }
        }
        
    }
    
    /**
     * construct the gradient in X, gradient in Y and the theta image from the 
     * given img and return the results as new GreyscaleImage[]{gX, gY, gXY, theta}.
     * 
     * @param img
     * @return 
     */
    protected GreyscaleImage[] createGradientProducts(final GreyscaleImage img) {
        
        GreyscaleImage gX, gY, g, theta;
        
        ImageProcessor ImageProcessor = new ImageProcessor();
        
        if (useLineDrawingMode) {
            
            gX = createGradientXFromDiffOfGauss(img);

            gY = createGradientYFromDiffOfGauss(img);

            removeOnePixelSpanningBorders(gX);
            
            removeOnePixelSpanningBorders(gY);
            
            theta = ImageProcessor.computeTheta(gX, gY);

            g = ImageProcessor.combineConvolvedImages(gX, gY);
            
        } else {
        
            gX = getGradientX1D(img);
            
            gY = getGradientY1D(img);
            
            removeOnePixelSpanningBorders(gX);
            
            removeOnePixelSpanningBorders(gY);
            
            theta = ImageProcessor.computeTheta(gX, gY);
            
            g = ImageProcessor.combineConvolvedImages(gX, gY);
       
        }
        
        return new GreyscaleImage[]{gX, gY, g, theta};
    }

    protected void applyHistogramEqualization(final GreyscaleImage input) {
         
        if (!doNotNormalizeByHistogram && !histogramEqualizationCheckFinished) {
            
            boolean doEqualization = histogramEqualizationIsNeeded(input);
            
            histogramEqualizationCheckFinished = true;
            
            log.fine("doEqualization=" + doEqualization);
            
            if (!doEqualization) {
                return;
            }
            
            HistogramEqualization hEq = new HistogramEqualization(input);
            hEq.applyFilter();
                                    
            // redo the histogram for future use
            HistogramHolder hist = createImageHistogram(input);
            
            imgHistogram = hist;            
        }
    }
    
    private HistogramHolder createImageHistogram(final GreyscaleImage input) {
        
        float[] pixValues = new float[input.getNPixels()];
        for (int i = 0; i < input.getNPixels(); i++) {
            pixValues[i] = input.getValues()[i];
        }

        float[] simulatedErrors = Errors.populateYErrorsBySqrt(pixValues);

        HistogramHolder hist = Histogram.createSimpleHistogram(
            0.0f, 256.f, 10, pixValues, simulatedErrors);

        return hist;
    }

    /**
     * determine if the range of image values should be stretched to fit
     * the possible range 0-255 and return true if so.  This has the side
     * effect of setting the instance variable imgHistogram to the generated
     * histogram needed to decide whether normalization is needed.
     * 
     * @param input
     * @return 
     */
    private boolean histogramEqualizationIsNeeded(final GreyscaleImage input) {
        
        HistogramHolder hist = createImageHistogram(input);
        
        this.imgHistogram = hist;
        
        float max = MiscMath.findMax(input.getValues());

        float min = MiscMath.findMin(input.getValues());
      
        if (max < (255.f * 0.8f)) {
            return true;
        }
        if (min > (255.f * 0.2f)) {
            return true;
        }

        double sum = 0;
        for (int i = 0; i < input.getNPixels(); i++) {
            sum += input.getValues()[i];
        }
        double avg = sum/input.getNPixels();
        
        int yMaxIdx = MiscMath.findYMaxIndex(hist.getYHist());

        float mode = hist.getXHist()[yMaxIdx];

        float ySumBelowMid = 0;
        float ySumAboveMid = 0;
        float yTot = 0;
        int len = hist.getXHist().length;
        for (int i = 0; i < len; i++) {
            float y = hist.getYHist()[i];
            if (hist.getXHist()[i] < 127) {
                ySumBelowMid += y;
            } else {
                ySumAboveMid += y;
            }
            yTot += y;
        }
                   
        if (min == 0 && max == 255) {
            // might be line drawings
            float yFrac1 = ySumAboveMid/yTot;
            float yFrac0 = ySumBelowMid/yTot;
            if (Math.abs(yFrac1 - yFrac0) > 0.25) {
                return true;
            }
            return false;
        }        
        
        // NOTE: this lower limit to min might need to be adjusted
        if (((min/255.f) < 0.1) && (mode < 127.)) {
            return false;
        }
        
        if ((mode/127.) < 0.9f) {
            //sometimes the mode is dominated by empty pixels, so check
            // the mean too.
            if ((avg/127.) < 0.9f) {
                                
                return true;
            }
        }
        
        return false;
    }

    public GreyscaleImage getGradientX() {
        return gX;
    }
    
    public GreyscaleImage getGradientY() {
        return gY;
    }
    
    public GreyscaleImage getGradientXY() {
        return gXY;
    }
    
    public GreyscaleImage getTheta() {
        return gTheta;
    }

    /**
     * @return the imgHistogram
     */
    public HistogramHolder getImgHistogram() {
        return imgHistogram;
    }
}
