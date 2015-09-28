package com.climbwithyourfeet.clustering;

import com.climbwithyourfeet.clustering.util.Histogram;
import com.climbwithyourfeet.clustering.util.HistogramHolder;
import com.climbwithyourfeet.clustering.util.MiscMath;
import com.climbwithyourfeet.clustering.util.PairInt;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
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
        
        TreeMap<Integer, Integer> freqMap = new TreeMap<Integer, Integer>();
        int minValue = Integer.MAX_VALUE;
        int maxValue = Integer.MIN_VALUE;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int v = distTrans[i][j];
                if (v < minValue) {
                    minValue = v;
                }
                if (v > maxValue) {
                    maxValue = v;
                }
                Integer key = Integer.valueOf(v);
                Integer c = freqMap.get(key);
                if (c == null) {
                    freqMap.put(key, Integer.valueOf(1));
                } else {
                    freqMap.put(key, Integer.valueOf(c.intValue() + 1));
                }
            }
        }
        Logger log = Logger.getLogger(DTClusterFinder.class.getName());
        log.info("dt minValue=" + minValue + " maxValue=" + maxValue);

        if (freqMap.size() == 1) {
            float dens = (float)(1./Math.sqrt(maxValue));
            return dens;
        }
        
        /*
        wanting to use freq map when the distance between clusters is large
        resulting in missing the separations of the largest voids.
        
        This may need to be revised, but it looks like histograms should be
        used when the freq map increases from last bin to first,
        else the freq map should be used, specifically the last bin of the
        frequency map.
        
        */
        
        boolean increasing = true;
        int lastCount = 0;
        for (Entry<Integer, Integer> entry : freqMap.entrySet()) {
            int count = entry.getValue().intValue();
            if (count < lastCount) {
                increasing = false;
                break;
            }
            lastCount = count;
            log.info("value=" + entry.getKey() + " count=" + entry.getValue());
        }
        
        //TODO: need to decide the comparison value
        int lastBinCount = freqMap.lastEntry().getValue().intValue();
        float lastBinFrac = (float)lastBinCount/(float)nPoints;
        
        log.info("lastBinFrac=" + lastBinFrac);
        
        if (increasing || (lastBinCount < 5) || (lastBinFrac < 0.5)) {
            
            float critDens = findCriticalDensity(distTrans);

            return critDens;
        }
        
        float critDens = findCriticalDensity(distTrans, freqMap);
        
        return critDens;
    }
    
    /**
     * using distance transform to find the critical density for clustering
     * return 0;
     * @return 
     */
     protected float findCriticalDensity(int[][] distTrans, 
         TreeMap<Integer, Integer> freqMap) {
         
         /*
         TODO:
             For datasets in which the clusters are much further away from one 
             another than the largest value in the frequency map,
             those furthest value points should be collected.
             We want the representative distances between those furthest value
             points (they uniquely span the distances between the clusters, or
             at least the closest among them as groups does).
             Need to group those furthest value points by the assoc radius 
             derived so far:
                 clustered within the distance 2./(2.5/math.sqrt(maxValue)).
             Then the distances between those largest value point groups should 
             be calculated (their coords are the centroids of the groups. so
             they are sometimes the centers of the cluster they surround).
             Then the final critical density would be 1./maxValueSepDistance.
         */
        int w = distTrans.length;
        int h = distTrans[0].length;

        Integer key = freqMap.lastKey();
        int maxValue = key.intValue();        
        
        // --- iterating until have more than 1 x coordinate to span a void ---
        Set<PairInt> maxValuePoints = new HashSet<PairInt>();
        Set<Integer> xs = new HashSet<Integer>();
        int nIter = 0;
        int nMax = freqMap.size();
        while (((nIter == 0) || (xs.size() == 1)) && (nIter < nMax)) {
            
            for (int i = 0; i < w; ++i) {
                for (int j = 0; j < h; ++j) {
                    int v = distTrans[i][j];
                    if (v == key.intValue()) {
                        maxValuePoints.add(new PairInt(i, j));
                        xs.add(Integer.valueOf(i));
                    }
                }
            }
            key = freqMap.lowerKey(key);
            nIter++;
        }
        
        DTGroupFinder fndr = new DTGroupFinder();
        float tmpCritDensity = (float) (1. / Math.sqrt(maxValue));
        fndr.calculateGroups(tmpCritDensity, maxValuePoints);
        
        int nGroups = fndr.getNumberOfGroups();
        
        List<Set<PairInt>> groups = new ArrayList<Set<PairInt>>();
        List<Float> centroidsX = new ArrayList<Float>();
        List<Float> centroidsY = new ArrayList<Float>();
        for (int i = 0; i < nGroups; ++i) {
            Set<PairInt> group = fndr.getGroup(i);
            groups.add(group);
            double[] xyCen = MiscMath.calculateXYCentroids(group);
            centroidsX.add(Float.valueOf((float) xyCen[0]));
            centroidsY.add(Float.valueOf((float) xyCen[1]));
        }
        
        /*
         can use closest pair algorithm on centroidsX, centroidsY
         and then critDensity = 1./closestDistance
         */
        ClosestPair closestPair = new ClosestPair();
        ClosestPair.ClosestPairFloat result = closestPair.findClosestPair(
            centroidsX, centroidsY);
        
        if (result.point0 == null) {
            return findCriticalDensity(distTrans);
        }
        
        float critDens = 1.f / result.separation;
        
        // density is due to 2 points:
        critDens /= 2.f;
        
        Logger.getLogger(this.getClass().getName()).info(
            "freqMap critDens=" + critDens);

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
                float half =  0.5f*(hist.getXHist()[1] -  hist.getXHist()[0]);
                if (yLimitIdx > -1) {
                    // a work around to next getting stuck at same size:
                    if ((len - yLimitIdx) < 5) {
                        yLimitIdx = len - 5;
                    }
                    float tmp = hist.getXHist()[yLimitIdx] - half;
                    if ((xl/tmp) > 15) {
                        // extreme zoom-in of high peak near idx=0
                        xl = tmp;
                        breakOnNext = true;
                    } else if (tmp < half) {
                        // extreme zoom-in of high peak near idx=0
                        xl = half;
                        breakOnNext = true;
                    } else if (tmp < xl) {
                        xl = tmp;
                    }
                } else {
                    xl = hist.getYHist()[1] - half;
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
