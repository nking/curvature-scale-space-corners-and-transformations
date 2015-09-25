package algorithms.compGeometry.clustering.distanceTransform;

import algorithms.compGeometry.clustering.distanceTransform.util.MiscMath;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.util.Errors;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class CriticalDensitySolver {
    
    private boolean debug = false;
    
    public CriticalDensitySolver() {
    }
    
    public void setToDebug() {
        debug = true;
    }
    
    /**
     * using histograms of inverse sqrt of distance transform, find the center 
     * of the first peak and return it, else return 0
     * (0 as a critical density should result in an infinite critical separation
     * so no clusters).
     * @param distTrans
     * @param nPoints
     * @param width
     * @param height
     * @return 
     */
    float findCriticalDensity(int[][] distTrans, int nPoints, int width, 
        int height) {
        
        int w = distTrans.length;
        int h = distTrans[0].length;
        
        float[] values = new float[w*h];
        int count2 = 0;
        for (int i0 = 0; i0 < w; ++i0) {
            for (int j0 = 0; j0 < h; ++j0) {
                int v = distTrans[i0][j0];
                values[count2] = (float)(1./Math.sqrt(v));
                count2++;
            }
        }
        
        float critDens = findCriticalDensity(values);
       
        return critDens;
    }
    
    /**
     * using histograms, find the center of the first peak and return it, else
     * return 0;
     * @param values
     * @return 
     */
    protected float findCriticalDensity(float[] values) {
         
        float[] vErrors = Errors.populateYErrorsBySqrt(values);

        List<HistogramHolder> histList = new ArrayList<HistogramHolder>();

        float xl = MiscMath.findMax(values);//4.0f;
        int nb = 40;
        
        /*
        steps xl down as needed and changes nb to 20 when small
        */
        
        int hc = 0;
        
        int maxHC = 10;
        
        while (hc < maxHC) {
            
            HistogramHolder hist = Histogram.createSimpleHistogram(
                0, xl, nb, values, vErrors);
        
            if (debug) {
                try {
                    hist.plotHistogram("clstr", "_cluster_" + hc);
                } catch (IOException ex) {
                    Logger.getLogger(CriticalDensitySolver.class.getName()).
                        log(Level.SEVERE, null, ex);
                }
            }
            
            if (hist == null || hist.getXHist() == null || hist.getXHist().length == 0) {
                break;
            }
            
            histList.add(hist);
            
            int len = hist.getXHist().length;
            
            int yMaxIdx = MiscMath.findYMaxIndex(hist.getYHist());
            
            int yMax = hist.getYHist()[yMaxIdx];
            int yLast = hist.getYHist()[len - 1];
            
            int yLimit = yMax/10;
            
            int yLimitIdx = -1;
            
            if (yLast < yLimit) {
                // shorten xmax and try again
                for (int i = (len - 1); i > -1; --i) {
                    int y = hist.getYHist()[i];
                    if (y >= yLimit) {
                        yLimitIdx = i;
                    }
                }
                if (yLimitIdx > -1) {
                    xl = hist.getXHist()[yLimitIdx];
                } else {
                    xl = hist.getYHist()[1];
                }
            } else {
                break;
            }
            
            hc++;
        }
        
        if (histList.isEmpty()) {
            // should default null result be density such that there are no clusters (~0)
            // or every point is its own cluster (infinity)
            return 0;
        }
        
        // find area of first peak, start search after first bin which
        // sometimes has a delta function.
        //TODO: this may need to be revised
        int firstNonZeroIdx = -1;
        int firstZeroAfterPeakIdx = -1;
        
        HistogramHolder hist = histList.get(histList.size() - 1);
        int len = hist.getXHist().length;
        
        for (int i = 1; i < len; ++i) {
            
            int y = hist.getYHist()[i];
            
            if (firstNonZeroIdx == -1) {
                if (y > 0) {
                    firstNonZeroIdx = i;
                }
            } else {
                // the start of the first peak has been found so look for the end
                if (y == 0) {
                    firstZeroAfterPeakIdx = i;
                    break;
                }
            }
        }
        
        // find weighted x from firstNonZeroIdx to firstZeroAfterPeakIdx
        float yPeakSum = 0;
        for (int i = firstNonZeroIdx; i < firstZeroAfterPeakIdx; ++i) {
            yPeakSum += hist.getYHist()[i];
        }
        
        float weightedX = 0;
        for (int i = firstNonZeroIdx; i < firstZeroAfterPeakIdx; ++i) {
            float w = hist.getYHistFloat()[i]/yPeakSum;
            weightedX += (w * hist.getXHist()[i]);
        }
        
        return weightedX;
    }
}
