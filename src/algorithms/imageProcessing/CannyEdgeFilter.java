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
    protected HistogramHolder imgHistogram = null;
        
    protected boolean histogramEqualizationCheckFinished = false;
    
    protected boolean histogramEqualizationWasPerformed = false;
    
    protected boolean doNotNormalizeByHistogram = false;
    
    protected float highThreshold = 5.0f;
    
    private boolean useLineDrawingMode = false;
    
    private GreyscaleImage gXY = null;
    
    private GreyscaleImage thetaXY = null;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    private Class<? extends ILineThinner> lineThinnerClass = ErosionFilter.class;
    
    protected boolean useOutdoorMode = false;
    
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
        
        this.useLineDrawingMode = true;
        
        highThreshold = 2.0f;
    }
       
    public void setLineThinnerClass(Class<? extends ILineThinner> cls) {
        lineThinnerClass = cls;
    }
    
    public void useOutdoorMode() {
        
        useOutdoorMode = true;
        
        highThreshold = 50.0f;
    }
    
    public void applyFilter(final GreyscaleImage input) {
        
        if (input.getWidth() < 3 || input.getHeight() < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }
                
        ImageProcesser imageProcesser = new ImageProcesser();
        
        imageProcesser.shrinkImageToFirstNonZeros(input);
             
        if (useOutdoorMode) {
            imageProcesser.blur(input, 2.0f); //3.0
        }

        applyHistogramEqualization(input);
        
        //[gX, gY, gXY, theta
        GreyscaleImage[] gradientProducts = createGradientProducts(input);
        
        gXY = gradientProducts[2].copyImage();
        
        thetaXY = gradientProducts[3].copyImage();
        
        input.resetTo(gradientProducts[2]);
           
        if (!useLineDrawingMode) {
            apply2LayerFilter(input);
        }
        
        applyLineThinnerFilter(input);
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        curveHelper.additionalThinning45DegreeEdges(gradientProducts[3], input);
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
       (TODO: the same correction might need to be applied to the low
       threshold intensity filter too.)
    
    */
    protected void apply2LayerFilter(final GreyscaleImage input) {
        
        if (input.getWidth() < 3 || input.getHeight() < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }
                  
        LowIntensityRemovalFilter filter2 = new LowIntensityRemovalFilter();
        
        if (useOutdoorMode) {
            filter2.overrideLowThresholdFactor(3.5f);
        }
                
        ImageStatistics stats = filter2.removeLowIntensityPixels(input,
            imgHistogram);
        
        double lowThresh = stats.getLowThresholdApplied();
                
        double threshold2 = lowThresh * highThreshold;
        
        // count number of pixels between lowThresh and threshold2 and
        // above threshold2.  the later should help scale highThreshold
        // factor from 3 to 5 when needed.
        int n0 = ImageStatisticsHelper.countPixels(input, (int)lowThresh, 
            (int)threshold2);
        
        int n1 = ImageStatisticsHelper.countPixels(input,
            (int)threshold2, 255);
        
        double r = (double)n1/(double)n0;
        
        log.fine("threshold2=" + threshold2 + " n0=" + n0 + " n1=" + n1 + 
            " n1/n0=" + r);
        
        GreyscaleImage img2 = new GreyscaleImage(input.getWidth(), 
            input.getHeight());
        
        // find points that are "sure-edge" points
        for (int i = 0; i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {
                int pixG = input.getValue(i, j);
                if (pixG >= threshold2) {
                    img2.setValue(i, j, pixG);
                }
            }
        }
        
        //TODO:  define "connected" as connected only to "sure-edge" points or
        //       connected to any point in the image in progress?
        // for now, choosing the first and keeping them separate.
        
        GreyscaleImage img3 = new GreyscaleImage(input.getWidth(), 
            input.getHeight());
        
        for (int i = 0; i < input.getWidth(); i++) {
            
            for (int j = 0; j < input.getHeight(); j++) {
                
                int pixG = input.getValue(i, j);
                
                if (img2.getValue(i, j) > 0) {                    
                
                    img3.setValue(i, j, pixG);
                
                } else if (pixG <= threshold2) {
                    
                    boolean found = true;
                    
                    for (int ii = -1; ii < 2; ii++) {
                        if (((i + ii) < 0) || ((i + ii) > 
                            (img3.getWidth() - 1))) {
                            continue;
                        }
                        for (int jj = -1; jj < 2; jj++) {
                            if (((j + jj) < 0) || ((j + jj) > 
                                (img3.getHeight() - 1))) {
                                continue;
                            }
                            if (img2.getValue(i + ii, j + jj) > 0) {
                                img3.setValue(i, j, pixG);
                                found = true;
                                break;
                            }
                        }
                        if (found) {
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
       
        ImageProcesser imageProcesser = new ImageProcesser();
        
        GreyscaleImage g = imageProcesser.subtractImages(g1, g0);
        
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
    private GreyscaleImage getGradientXID(final GreyscaleImage input) {
                
        return getGradientID(input, true);
    }
    
    /**
     * convolve the image with a Sobel Y 1D kernel which is the same as a 
     * Gaussian first derivative with sigma = sqrt(1)/2.
     * 
     * @param input
     * @return 
     */
    private GreyscaleImage getGradientYID(final GreyscaleImage input) {
                
        return getGradientID(input, false);
    }
    
    /**
     * convolve the image with a Sobel 1D kernel which is the same as a 
     * Gaussian first derivative with sigma = sqrt(1)/2.
     * 
     * @param input
     * @return 
     */
    private GreyscaleImage getGradientID(final GreyscaleImage input, 
        boolean calculateForX) {
        
        log.fine("getGradientID calculateForX=" + calculateForX);
                
        // 0.5f, -0.0f, -0.5f
        float[] kernel = Gaussian1DFirstDeriv.getKernel(
            SIGMA.ZEROPOINTFIVE);
                
        GreyscaleImage output = input.copyImage();
        
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
    
    private void applyLineThinnerFilter(final GreyscaleImage img) {
                
        ILineThinner lineThinner;
        
        try {
            
            lineThinner = lineThinnerClass.newInstance();
            
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
     * given img and return the results as new GreyscalImage[]{gX, gY, gXY, theta}.
     * 
     * @param img
     * @return 
     */
    protected GreyscaleImage[] createGradientProducts(final GreyscaleImage img) {
        
        GreyscaleImage gX, gY, g, theta;
        
        ImageProcesser imageProcesser = new ImageProcesser();
        
        if (useLineDrawingMode) {
            
            gX = createGradientXFromDiffOfGauss(img);

            gY = createGradientYFromDiffOfGauss(img);

            removeOnePixelSpanningBorders(gX);
            
            removeOnePixelSpanningBorders(gY);
            
            theta = imageProcesser.computeTheta(gX, gY);

            g = imageProcesser.combineConvolvedImages(gX, gY);
            
        } else {
        
            gX = getGradientXID(img);
            
            gY = getGradientYID(img);
            
            removeOnePixelSpanningBorders(gX);
            
            removeOnePixelSpanningBorders(gY);
            
            theta = imageProcesser.computeTheta(gX, gY);
            
            g = imageProcesser.combineConvolvedImages(gX, gY);
       
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
            
            histogramEqualizationWasPerformed = true;
                        
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

    public GreyscaleImage getGradientXY() {
        return gXY;
    }
    
    public GreyscaleImage getThetaXY() {
        return thetaXY;
    }
}
