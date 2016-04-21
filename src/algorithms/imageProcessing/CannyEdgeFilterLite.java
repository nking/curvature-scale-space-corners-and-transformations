package algorithms.imageProcessing;

import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.Misc;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

/**
   this version of the Canny edge detector has only a few of the
   characteristics of the standard canny edge detector.
   -- creates a gradient image
   -- applies a 2 layer hysteresis filter to remove low intensity
      pixels.
   -- applies a line thinner instead of the non-maximal suppression
      algorithm 

   no orientation image is available as a result.
  
 * @author nichole
 */
public class CannyEdgeFilterLite {
              
    protected float highThreshold = 2.0f;
    
    public static float defaultLowThreshold = 
        LowIntensityRemovalFilter.defaultLowThreshFactor;
        
    protected float lowThreshold = defaultLowThreshold;
    
    protected double lowThresholdApplied2Layer = Float.MIN_VALUE;
            
    protected GreyscaleImage gradientXY = null;
    
    private boolean useDiffOfGauss = false;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    // ErosionFilter
    private Class<? extends ILineThinner> lineThinnerClass = ZhangSuenLineThinner.class;
    
    protected boolean useOutdoorMode = false;
    
    protected boolean overrideToUseSobel = false;
    
    protected int[] shrinkToSize = null;
    
    public CannyEdgeFilterLite() {        
    }
    
    public void setToUseSobel() {
        this.overrideToUseSobel = true;
    }
       
    public void setLineThinnerClass(Class<? extends ILineThinner> cls) {
        lineThinnerClass = cls;
    }
    
    public void setToUseDiffOfGauss() {
        
        //TODO: consider dilate before erosion
        
        this.useDiffOfGauss = true;        
    }
    
    public void overrideHighThreshold(float thresh) {
        highThreshold = thresh;
    }
    
    public void overrideLowThreshold(float thresh) {
        lowThreshold = thresh;
    }
    
    public void applyFilter(final GreyscaleImage input) {
        
        if (input.getWidth() < 3 || input.getHeight() < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }
        
        if (overrideToUseSobel) {
            gradientXY = input.copyImage();
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.applySobelKernel(gradientXY);
        } else if (useDiffOfGauss) {
            gradientXY = createDiffOfGaussians(input);
        } else {
            gradientXY = createGradient(input);
        }
                        
        input.resetTo(gradientXY);
                        
        apply2LayerFilter(input);

        for (int i = 0; i < input.getNPixels(); ++i) {
            int v = input.getValue(i);
            if (v < 0) {
                input.setValue(i, 0);
            }
        }
    }
    
    /**
     * apply the canny edge filter only to the region with bounding box
     * lower left corner (xLL, yLL) to upper right corner (xUR, yUR) where
     * xUR > yLL and yUR > yLL and only to pixels which hold a sentinel value
     * in the output (to avoid repeating a calculation).
     * It's the user's responsibility to make sure that all points in the region
     * are present in the input.
     * @param input values from image for the bounded region (xLL, yLL) to (xUR, yUR).
     * @param output image of edges calculated for the given bounding box points.
     * @param xLL lower left corner x coordinate of bounding box.  xLL is less than xUR.
     * @param yLL lower left corner y coordinate of bounding box.  yLL is less than yUR.
     * @param xUR upper right corner x coordinate of bounding box.  xUR is larger than xLL.
     * @param yUR upper right corner y coordinate of bounding box.  yUR is larger than yLL.
     */
    public void applyFilterToRegion(final Map<PairInt, Integer> input,
        final GreyscaleImage output, int xLL, int yLL, int xUR, int yUR) {
        
        if (!(xLL < xUR)) {
            throw new IllegalArgumentException("xLL must be less than xUR");
        }
        if (!(yLL < yUR)) {
            throw new IllegalArgumentException("yLL must be less than yUR");
        }
        
        Map<PairInt, Integer> gradientValues = new HashMap<PairInt, Integer>();
        if (overrideToUseSobel) {
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.applySobelKernel(input,
                xLL, yLL, xUR, yUR, gradientValues);
        } else if (useDiffOfGauss) {
            createDiffOfGaussians(input, xLL, yLL, xUR, yUR, gradientValues);
        } else {
            createGradient(input, xLL, yLL, xUR, yUR, gradientValues);
        }
                     
        apply2LayerFilter(gradientValues, xLL, yLL, xUR, yUR,
            output.getWidth(), output.getHeight());
        
        for (int x = xLL; x <= xUR; ++x) {
            for (int y = yLL; y < yUR; ++y) {
                PairInt p = new PairInt(x, y);
                Integer v = gradientValues.get(p);
                if (v != null) {
                    output.setValue(p.getX(), p.getY(), v.intValue());
                }
            }
        }
        
        gradientXY = output;
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
        
        HistogramHolder imgHistogram = createImageHistogram(input);
        
        LowIntensityRemovalFilter filter2 = new LowIntensityRemovalFilter();
        
        if (lowThreshold != defaultLowThreshold) {
            filter2.overrideLowThresholdFactor(lowThreshold);
        }
        
        ImageStatistics stats = filter2.removeLowIntensityPixels(input,
            imgHistogram);
        
        lowThresholdApplied2Layer = stats.getLowThresholdApplied();
                
        log.info("2-layer filter: low thresh=" + lowThresholdApplied2Layer);
        
        double threshold2 = lowThresholdApplied2Layer * highThreshold;
        
        // count number of pixels between lowThresh and threshold2 and
        // above threshold2.  the later should help scale highThreshold
        // factor from 3 to 5 when needed.
        int n0 = ImageStatisticsHelper.countPixels(input, (int)lowThresholdApplied2Layer, 
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
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
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
    
    /**
     * note that settings for lowThreshold, lowThresholdApplied2Layer, and
     * highThreshold are not used because otsu's threshold algorithm is used
     * instead.
     * 
     */
    protected void apply2LayerFilter(Map<PairInt, Integer> gradientValues,
        int xLL, int yLL, int xUR, int yUR, int width, int height) {
        
        float otsuScaleFactor = 0.75f;
        float factorBelowHighThreshold = 2.f;
        
        OtsuThresholding ot = new OtsuThresholding();
        float tHigh = tHigh = otsuScaleFactor * ot.calculateBinaryThreshold256(
            gradientValues);
        float tLow = tHigh/factorBelowHighThreshold;
        
        Map<PairInt, Integer> output = new HashMap<PairInt, Integer>();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        for (int x = xLL; x <= xUR; ++x) {
            for (int y = yLL; y <= yUR; ++y) {
                
                PairInt p = new PairInt(x, y);
        
                Integer v = gradientValues.get(p);
                
                if (v < tLow) {
                    continue;
                } else if (v > tHigh) {
                    output.put(p, v);
                    continue;
                }
            
                boolean foundHigh = false;
                boolean foundMid = false;

                for (int k = 0; k < dxs.length; ++k) {                
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    PairInt p2 = new PairInt(x2, y2);
                    Integer v2 = gradientValues.get(p2);
                    if (v2 == null) {
                        continue;
                    }
                    if (v2.intValue() > tHigh) {
                        foundHigh = true;
                        break;
                    } else if (v2.intValue() > tLow) {
                        foundMid = true;
                    }
                }
                if (foundHigh) {
                    output.put(p, v);
                    continue;
                }
                if (!foundMid) {
                    continue;
                }
                // search the 5 by 5 region for a "sure edge" pixel
                for (int dx = -2; dx <= 2; ++dx) {
                    int x2 = x + dx;
                    for (int dy = -2; dy <= 2; ++dy) {
                        int y2 = y + dy;
                        if (x2 == x && y2 == y) {
                            continue;
                        }
                        PairInt p2 = new PairInt(x2, y2);
                        Integer v2 = gradientValues.get(p2);
                        if (v2 == null) {
                            continue;
                        }
                        if (v2.intValue() > tHigh) {
                            output.put(p, v);
                            foundHigh = true;
                            break;
                        }
                    }
                    if (foundHigh) {
                        break;
                    }
                }
            }
        }
                
        gradientValues.clear();
        
        gradientValues.putAll(output);
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
    
    void applyLineThinnerFilter(final GreyscaleImage img) {
                
        ILineThinner lineThinner;
        
        try {
            
            lineThinner = lineThinnerClass.newInstance();
            
            lineThinner.setEdgeGuideImage(gradientXY);
                        
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
     * construct the gradient in X, gradient in Y and combine them
     * 
     * @param img
     * @return 
     */
    protected GreyscaleImage createGradient(final GreyscaleImage img) {
        
        GreyscaleImage gX, gY, g;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        gX = getGradientX1D(img);

        gY = getGradientY1D(img);

        removeOnePixelSpanningBorders(gX);

        removeOnePixelSpanningBorders(gY);

        g = imageProcessor.combineConvolvedImages(gX, gY);
        
        return g;
    }

    protected void createGradient(final Map<PairInt, 
        Integer> input, int xLL, int yLL, int xUR, int yUR,
        Map<PairInt, Integer> outputGradientValues) {
            
        ImageProcessor imageProcessor = new ImageProcessor();
        
        float[] kernel = Gaussian1DFirstDeriv.getKernel(
            SIGMA.ZEROPOINTFIVE);
        
        imageProcessor.applyKernel(input, xLL, yLL, xUR, yUR, kernel,
            outputGradientValues);
    }
    
    protected GreyscaleImage createDiffOfGaussians(final GreyscaleImage img) {
        
        GreyscaleImage gX, gY, g;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        gX = createGradientFromDiffOfGauss(img, true);

        gY = createGradientFromDiffOfGauss(img, false);

        removeOnePixelSpanningBorders(gX);
            
        removeOnePixelSpanningBorders(gY);
            
        g = imageProcessor.combineConvolvedImages(gX, gY);
        
        return g;
    }
    
    protected void createDiffOfGaussians(final Map<PairInt, Integer> input,
        int xLL, int yLL, int xUR, int yUR,
        Map<PairInt, Integer> outputGradientValue
        ) {
                
        Map<PairInt, Integer> gX = createGradientFromDiffOfGauss(input,
            xLL, yLL, xUR, yUR, true);

        Map<PairInt, Integer> gY = createGradientFromDiffOfGauss(input,
            xLL, yLL, xUR, yUR, false);
        
        for (int x = xLL; x <= xUR; ++x) {
            for (int y = yLL; y < yUR; ++y) {
                
                PairInt p = new PairInt(x, y);
                
                int vX = gX.get(p).intValue();

                int vY = gY.get(p).intValue();
                
                int v = (int) Math.round(Math.sqrt(vX * vX + vY * vY));

                outputGradientValue.put(p, Integer.valueOf(v));
            }
        }        
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
    
    private Map<PairInt, Integer> createGradientFromDiffOfGauss(
        final Map<PairInt, Integer> pointValues,
        int xLL, int yLL, int xUR, int yUR, 
        boolean calculateForX) {
        
        /*
        what looks really good is computationally too long:
            g1 is convolved by a kernel > 1000*g0 kernel which uses sigma=1
        */
        
        /* 
        for line drawings, use sigma1 = 0.42466090014400953f and sigma2 = 2*sigma1
        */
        
        float resultOne = 0.42466090014400953f;
        
        float sigma = 1.f * resultOne;
        
        float[] kernel0 = Gaussian1D.getKernel(sigma);

        float[] kernel2 = Gaussian1D.getKernel(sigma * 1.6f);
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Map<PairInt, Integer> c0 = imageProcessor.applyKernel(pointValues,
            xLL, yLL, xUR, yUR,
            kernel0, calculateForX);
                    
        Map<PairInt, Integer> c2 = imageProcessor.applyKernel(pointValues,
            xLL, yLL, xUR, yUR, 
            kernel2, calculateForX);
        
        // subtract the two
        Map<PairInt, Integer> output = new HashMap<PairInt, Integer>();
        
        for (int xp = xLL; xp <= xUR; ++xp) {
            for (int yp = yLL; yp <= yUR; ++yp) {
                
                PairInt p = new PairInt(xp, yp);
                
                int v0 = c0.get(p).intValue();
                int v2 = c2.get(p).intValue();
                output.put(p, v2 - v0);
            }
        }
            
        return output;
    }
    
    private HistogramHolder createImageHistogram(final GreyscaleImage input) {
        
        float[] pixValues = new float[input.getNPixels()];
        for (int i = 0; i < input.getNPixels(); i++) {
            pixValues[i] = input.getValue(i);
        }

        float[] simulatedErrors = Errors.populateYErrorsBySqrt(pixValues);

        HistogramHolder hist = Histogram.createSimpleHistogram(
            0.0f, 256.f, 10, pixValues, simulatedErrors);

        return hist;
    }

}
