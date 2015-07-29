package algorithms.imageProcessing;

import java.util.logging.Logger;

/**
  calculate the transformation function for an image to perform histogram
  equalization on it, that is change the range of intensities to fill the
  possible range.
  
 * @author nichole
 */
public class HistogramEqualizationForColor {
    
    protected Image img = null;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    private final long[] rHist = new long[256];
    
    private final long[] rHistC = new long[256];
    
    protected long rHistCMin = Long.MAX_VALUE;
    
    private final long[] gHist = new long[256];
    
    private final long[] gHistC = new long[256];
    
    protected long gHistCMin = Long.MAX_VALUE;
    
    private final long[] bHist = new long[256];
    
    private final long[] bHistC = new long[256];
    
    protected long bHistCMin = Long.MAX_VALUE;
    
    private boolean finished = false;
    
    /**
     * operate on the image to perform histogram normalization of it.
     * 
     * @param input 
     */
    public HistogramEqualizationForColor(Image input) {
        
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
        
        if (finished) {
            return;
        }
        
        calculateHistogram();
        
        calculateCumulativeHistogram();
        
        applyTransformationFunction();
        
        finished = true;
    }
   
    protected void calculateHistogram() {
        
        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
                int r = img.getR(i, j);
                rHist[r]++;
                int g = img.getG(i, j);
                gHist[g]++;
                int b = img.getB(i, j);
                bHist[b]++;
            }
        }
    }
    
    protected void calculateCumulativeHistogram() {
        
        rHistC[0] = rHist[0];        
        for (int i = 1; i < rHist.length; i++) {
            rHistC[i] = rHistC[i - 1] + rHist[i];
            if ((rHistC[i] > 0) && (rHistC[i] < rHistCMin)) {
                rHistCMin = rHistC[i];
            }
        }
        gHistC[0] = gHist[0];        
        for (int i = 1; i < gHist.length; i++) {
            gHistC[i] = gHistC[i - 1] + gHist[i];
            if ((gHistC[i] > 0) && (gHistC[i] < gHistCMin)) {
                gHistCMin = gHistC[i];
            }
        }
        bHistC[0] = bHist[0];        
        for (int i = 1; i < bHist.length; i++) {
            bHistC[i] = bHistC[i - 1] + bHist[i];
            if ((bHistC[i] > 0) && (bHistC[i] < bHistCMin)) {
                bHistCMin = bHistC[i];
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
                
                int r = img.getR(i, j);
                long rC = rHistC[r];
                
                int g = img.getG(i, j);
                long gC = gHistC[g];
                
                int b = img.getB(i, j);
                long bC = bHistC[b];
                
                /*
                the implied max should be explicit...
                
                sz - min = (x - min)*255
                
                sz - min = (x - min)*(255 - 1)
                
                */
                
                int rT = (r == 0) ? 0 : 
                    Math.round((rC - rHistCMin)*254.f/(sz - rHistCMin)) + 1;
                
                int gT = (g == 0) ? 0 : 
                    Math.round((gC - gHistCMin)*254.f/(sz - gHistCMin)) + 1;
                
                int bT = (b == 0) ? 0 : 
                    Math.round((bC - bHistCMin)*254.f/(sz - bHistCMin)) + 1;
               
                img.setRGB(i, j, rT, gT, bT);
            }
        }        
    }
    
    /**
     * return the values for the histogram
     * @return the rHist
     */
    public long[] getRHist() {
        return rHist;
    }

    /**
     * return the values for the cumulative histogram
     * @return the rHistC
     */
    public long[] getRHistC() {
        return rHistC;
    }

    public long getRHistCMin() {
        return rHistCMin;
    }
}
