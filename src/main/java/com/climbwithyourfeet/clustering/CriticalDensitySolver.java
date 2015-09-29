package com.climbwithyourfeet.clustering;

import com.climbwithyourfeet.clustering.util.Histogram;
import com.climbwithyourfeet.clustering.util.HistogramHolder;
import com.climbwithyourfeet.clustering.util.MiscMath;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class CriticalDensitySolver {
        
    private boolean debug = false;
    
    /**
     *
     */
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     *
     */
    public CriticalDensitySolver() {
    }
    
    /**
     *
     */
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
        
        float critDens = findCriticalDensity(distTrans);

        return critDens;
        
    }
    
    /**
     * using histograms of 1/sqrt(distanceTransform[i][j]), find the center of 
     * the first peak and return it, else
     * return 0;
     * @param values
     * @return 
     */
    protected float findCriticalDensity(float[] values) {
         
        float[] vErrors = Histogram.populateYErrorsBySqrt(values);

        List<HistogramHolder> histList = new ArrayList<HistogramHolder>();

        float xl = MiscMath.findMax(values);
        int nb = 40;
        
        /*
        steps xl down as needed and changes nb to 20 when small
        */
        
        int hc = 0;
        
        int maxHC = 5;
        
        boolean breakOnNext = false;
        
        while (hc < maxHC) {
                        
            HistogramHolder hist = Histogram.createSimpleHistogram(
                0, xl, nb, values, vErrors);
        
            if (debug) {
                String outFileSuffix = "_cluster_" + hc;
                hist.plotHistogram("clstr", outFileSuffix);
            }
            
            if (hist == null || hist.getXHist() == null || hist.getXHist().length == 0) {
                break;
            }
            
            histList.add(hist);
            
            if (breakOnNext) {
                break;
            }
            
            int len = hist.getXHist().length;
            
            double areaH0 = MiscMath.calculateArea(hist, 0, (len/2) - 1);
            double areaH1 = MiscMath.calculateArea(hist, (len/2), len - 1);
                        
            if ((areaH1/areaH0) > 0.75) {
                // decrease the number of bins
                if (nb <= 10) {
                    break;
                }
                nb /= 2;
                hc++;
                continue;
            }
            
            int yMaxIdx = MiscMath.findYMaxIndex(hist.getYHist());
                        
            int yMax = hist.getYHist()[yMaxIdx];
            int yLast = hist.getYHist()[len - 1];
            
            int yLimit = yMax/10;
                        
            if (hc > 0) {
                yLimit = yMax/15;
            }
            
            int yLimitIdx = -1;
                        
            if (yLast < yLimit) {
                // shorten xmax and try again
                for (int i = (len - 1); i > -1; --i) {
                    int y = hist.getYHist()[i];
                    if (y >= yLimit) {
                        yLimitIdx = i;
                        break;
                    }
                }
                if (nb == 40) {
                    nb = 20;
                }
                float halfBin =  0.5f*(hist.getXHist()[1] -  hist.getXHist()[0]);
                if (yLimitIdx > -1) {
                    // a work around to next getting stuck at same size:
                    if ((len - yLimitIdx) < 5) {
                        yLimitIdx = len - 5;
                    }
                    float tmp = hist.getXHist()[yLimitIdx] - halfBin;
                    if ((xl/tmp) > 15) {
                        // extreme zoom-in of high peak near idx=0
                        if (yLimitIdx == 0) {
                            xl = hist.getXHist()[1] - halfBin;
                        } else {
                            xl = tmp;
                        }
                        breakOnNext = true;
                    } else if (tmp < halfBin) {
                        // extreme zoom-in of high peak near idx=0
                        xl = halfBin;
                        breakOnNext = true;
                    } else if (tmp < xl) {
                        xl = tmp;
                    }
                } else {
                    xl = hist.getXHist()[1] - halfBin;
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
        
        HistogramHolder hist = histList.get(histList.size() - 1);
        int len = hist.getXHist().length;
        
        // find area of first peak, start search after first bin which
        // sometimes has a delta function.
        //TODO: this may need to be revised
        int firstNonZeroIdx = -1;
        int firstZeroAfterPeakIdx = len - 1;
        
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
                
        if (firstNonZeroIdx == -1) {
            return 0;
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
                
        //wanting an answer that is a little higher 
        // than the weghted center but still within the bounds of the peak.
        /*int nh = (firstZeroAfterPeakIdx - firstNonZeroIdx)/2;
        double areaH0 = MiscMath.calculateArea(hist, firstNonZeroIdx, nh);
        double areaH1 = MiscMath.calculateArea(hist, nh + 1, firstZeroAfterPeakIdx);
        */
        float frac9 = (0.9f * hist.getXHist()[firstZeroAfterPeakIdx] )
            + (0.1f * hist.getXHist()[firstNonZeroIdx]);
        
        if (weightedX < frac9) {
            weightedX = 0.5f * (weightedX + frac9);
        }
        
        return weightedX;
    }

    private float findCriticalDensity(int[][] distTrans) {
        
        int w = distTrans.length;
        int h = distTrans[0].length;
        
        float[] values = new float[w * h];
        int count2 = 0;
        for (int i0 = 0; i0 < w; ++i0) {
            for (int j0 = 0; j0 < h; ++j0) {
                int v = distTrans[i0][j0];
                values[count2] = (float) (1. / Math.sqrt(v));
                count2++;
            }
        }

        return findCriticalDensity(values);
    }

}
