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
        
        CannyEdgeFilterAdaptive cannyL = new CannyEdgeFilterAdaptive();
        CannyEdgeFilterAdaptive cannyC = new CannyEdgeFilterAdaptive();
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
        cannyC.setOtsuScaleFactor(1.0f);
        
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

}
