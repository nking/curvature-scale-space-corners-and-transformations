package algorithms.imageProcessing;

import java.util.logging.Logger;

/**
  calculate the transformation function for an image to perform histogram
  equalization on it, that is change the range of intensities to fill the
  possible range.
  
 * @author nichole
 */
public class HistogramEqualization {
    
    protected GreyscaleImage img = null;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    private final long[] aHist = new long[256];
    
    private final long[] aHistC = new long[256];
    
    protected long aHistCMin = Long.MAX_VALUE;
    
    /**
     * operate on the image to perform histogram normalization of it.
     * 
     * @param input 
     */
    public HistogramEqualization(GreyscaleImage input) {
        
        this.img = input;
    }
    
    /**
     * apply the transformation function to the image to stretch the intensities
     * to use the full range for the color band 0 to 255.  Note that a slight 
     * difference
     * has been applied with standard equation.  The standard scales all
     * values > 0 to between 0 and 255.  This scales all values > 0 to
     * between 1 and 255 so a count remains significant and a 0 remains a zero.
     */
    public void applyFilter() {
        
        calculateHistogram();
        
        calculateCumulativeHistogram();
        
        applyTransformationFunction();
    }
   
    protected void calculateHistogram() {
        
        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
                int a = img.getValue(i, j);
                aHist[a]++;
            }
        }
    }
    
    protected void calculateCumulativeHistogram() {
        
        aHistC[0] = aHist[0];
        
        for (int i = 1; i < aHist.length; i++) {
            aHistC[i] = aHistC[i - 1] + aHist[i];
            
            if ((aHistC[i] > 0) && (aHistC[i] < aHistCMin)) {
                aHistCMin = aHistC[i];
            }
        }
        
    }
    
    /**
     * apply the transformation function.  Note that a slight difference
     * has been applied with standard equation.  The standard scales all
     * values > 0 to between 0 and 255.  This scales all values > 0 to
     * between 1 and 255 so a count remains significant and a 0 remains a zero.
     */
    protected void applyTransformationFunction() {
        /*
                      ( cdf(v) - cdf_min            )
          h(v) = round( ----------------- * (L - 1) )
                      ( (M * N) - cdf_min           )
        */      
        
        // normalize values to be between 0 and 255... seems like it should
        // be 1 to 255 so a pixel w/ signal retains a signal.
        long sz = img.getWidth() * img.getHeight();
        
        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
                
                int a = img.getValue(i, j);
                
                long aC = aHistC[a];
                
                /*
                the implied max should be explicit...sz pixels has something 
                to do with that
                
                sz - min = (x - min)*255
                
                sz - min = (x - min)*(255 - 1)
                
                */
                
                int aT = (a == 0) ? 0 : 
                    Math.round((aC - aHistCMin)*254.f/(sz - aHistCMin)) + 1;
               
                img.setValue(i, j, aT);
            }
        }        
    }
    
    /**
     * return the values for the histogram
     * @return the rHist
     */
    public long[] getHist() {
        return aHist;
    }

    /**
     * return the values for the cumulative histogram
     * @return the rHistC
     */
    public long[] getHistC() {
        return aHistC;
    }

    public long getHistCMin() {
        return aHistCMin;
    }
}
