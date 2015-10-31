package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * A quick algorithm for matching points around ideal closed curves in
 * two images.  The algorithm uses point ordering around the curve and
 * needs the same number of corners in both curves.
 * 
 * The runtime complexity is O(N_corners).
 * 
 * It can be used and followed by ClosedCurveCornerMatcher if it fails.
 * 
 * @author nichole
 */
public class ClosedCurveCornerMatcher0 {
    
    private final Logger log = Logger.getLogger(this.getClass().getName());

    private boolean solverHasFinished = false;

    protected final int dither = 3;//6
    
    //protected final int degreeIntervals = 20;

    //protected final int rotationTolerance = 20;

    private TransformationPair4 solutionTransformationPair = null;
    
    public ClosedCurveCornerMatcher0() {
    }

    /**
     * find solution using corners paired by curve order.
     * The algorithm requires c1.size() == c2.size().
     * with runtime O(n).
     * 
     */
    public boolean matchCorners(final IntensityFeatures features1,
        final IntensityFeatures features2, final List<CornerRegion> corners1, 
        final List<CornerRegion> corners2, boolean cornersAreAlreadyCCW,
        GreyscaleImage img1, GreyscaleImage img2) {
        
        if (solverHasFinished) {
            throw new IllegalStateException(
            "matchContours cannot be invoked more than once");
        }

        List<CornerRegion> c1 = new ArrayList<CornerRegion>(corners1.size());
        List<CornerRegion> c2 = new ArrayList<CornerRegion>(corners2.size());        
        if (c1.size() != c2.size()) {
            throw new IllegalArgumentException("c1 and c2 must be same size");
        }
        c1.addAll(corners1);
        c2.addAll(corners2);
        if (!cornersAreAlreadyCCW) {
            sortCornersToCCW(c2);
            sortCornersToCCW(c2);
        }
        
        int deltaIdx = 0;
        
        FeatureMatcher featureMatcher = new FeatureMatcher();
        
        List<TransformationPair4> trList = new ArrayList<TransformationPair4>();        
        
        while (deltaIdx < c1.size()) {
            
            TransformationPair4 transformationPair 
                = new TransformationPair4(0, deltaIdx, c1.size());
            
            // store matched and cost
            for (int i = 0; i < c1.size(); ++i) {
                
                int j = i + deltaIdx;
                if (j > (c2.size() - 1)) {
                    j = j - c2.size();
                }
                
                CornerRegion region1 = c1.get(i);
                CornerRegion region2 = c2.get(j);
                
                FeatureComparisonStat compStat = null;
                try {
                    compStat = featureMatcher.ditherAndRotateForBestLocation(
                        features1, features2, region1, region2, dither, img1, 
                        img2);
                } catch (CornerRegion.CornerRegionDegneracyException ex) {
                    log.fine(ex.getMessage());
                }
               
                if (compStat != null) {
                    if (compStat.getSumIntensitySqDiff() < compStat.getImg2PointIntensityErr()) {
                        transformationPair.addMatched(region1, region2, compStat);
                    }
                }                
            }
            
            // need at least 3 point for a cost evaluation
            if (transformationPair.getMatchedCornerRegions1().size() > 2) {
                trList.add(transformationPair);
            }
            
            deltaIdx++;
        }
        
        // chose best solution if there are more than 1
        int nMaxMatched = Integer.MIN_VALUE;
        double minCost = Double.MAX_VALUE;
        int minCostIdx = -1;
        for (int i = 0; i < trList.size(); ++i) {
            
            TransformationPair4 transP = trList.get(i);
            
            List<Integer> removedIndexes = 
                FeatureMatcher.removeDiscrepantThetaDiff(
                    transP.getMatchedCompStats());
            
            if (!removedIndexes.isEmpty()) {
                for (int ii = removedIndexes.size() - 1; ii > -1; --ii) {
                    int idx = removedIndexes.get(ii);
                    transP.getMatchedCornerRegions1().remove(idx);
                    transP.getMatchedCornerRegions2().remove(idx);
                }
            }
            
            if (transP.getMatchedCompStats().size() < 3) {
                continue;
            }
            
            double cost = calculateCombinedIntensityStat(transP.getMatchedCompStats());
            transP.setCost(cost);
            
            if (((transP.getMatchedCornerRegions1().size() >= nMaxMatched) &&
                (cost < minCost)) ||
                ((cost == minCost) 
                && (transP.getMatchedCornerRegions1().size() > nMaxMatched))) {
                
                minCost = cost;
                minCostIdx = i;
                nMaxMatched = transP.getMatchedCornerRegions1().size();
            }
        }
        
        solverHasFinished = true;
        
        if (minCostIdx > -1) {
            
            solutionTransformationPair = trList.get(minCostIdx);
            
            return true;
        }
        
        return false;
    }

    protected double calculateCombinedIntensityStat(List<FeatureComparisonStat> 
        compStats) {
        
        double sum = 0;
        
        for (FeatureComparisonStat compStat : compStats) {
            sum += compStat.getSumIntensitySqDiff();
        }
        
        sum /= (double) compStats.size();
        
        return sum;
    }

    /**
     * get the curvature scale space images shift between the
     * first set of contours and the second set.
     * @return
     */
    public TransformationPair4 getTransformationPair() {
        return solutionTransformationPair;
    }
    
    private void sortCornersToCCW(List<CornerRegion> corners) {
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        PairIntArray cornerXY = new PairIntArray();
        for (int ii = 0; ii < corners.size(); ++ii) {
            CornerRegion cr = corners.get(ii);
            cornerXY.add(cr.getX()[cr.getKMaxIdx()], cr.getY()[cr.getKMaxIdx()]);
        }
        boolean isCW = curveHelper.curveIsOrderedClockwise(cornerXY);
        if (isCW) {
            int n = corners.size();
            if (n > 1) {
                int end = n >> 1;
                // 0 1 2 3 4
                for (int ii = 0; ii < end; ii++) {
                    int idx2 = n - ii - 1;
                    CornerRegion swap = corners.get(ii);
                    corners.set(ii, corners.get(idx2));
                    corners.set(idx2, swap);
                }
            }
        }
    }
    
}
