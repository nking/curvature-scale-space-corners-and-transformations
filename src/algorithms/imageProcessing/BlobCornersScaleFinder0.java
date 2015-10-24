package algorithms.imageProcessing;

import algorithms.compGeometry.clustering.FixedDistanceGroupFinder;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * class to attempt to find Euclidean scale transformation between image1
 * and image2 using corners of a blob and an assumption of
 * ordered matching between them.
 * The method is fast, but cannot be used to compare the same object boundaries
 * if the boundaries have different numbers of corners.
 * 
 * @author nichole
 */
public class BlobCornersScaleFinder0 extends AbstractBlobScaleFinder {

    protected void filterToSameNumberIfPossible(
        PairIntArray curve1, List<CornerRegion> regions1,
        PairIntArray curve2, List<CornerRegion> regions2) {
        
        int dist = 2;
        
        if (regions1.size() != regions2.size()) {
            filterCorners(curve1, regions1, dist);
            filterCorners(curve2, regions2, dist);
        }
    }
    
    protected void filterCorners(PairIntArray curve, 
        List<CornerRegion> regions, int dist) {
        
        /*
        if there are more than 1 corner within dist of 2 or so of on another,
        remove all except strongest corner.
        */
        List<Set<Integer>> closeCornerIndexes = findCloseCorners(dist, regions);
        
        if (closeCornerIndexes.isEmpty()) {
            return;
        }
        
        List<Integer> remove = new ArrayList<Integer>();
        for (Set<Integer> set : closeCornerIndexes) {
            float maxK = Float.MIN_VALUE;
            Integer maxKIndex = null;
            for (Integer index : set) {
                CornerRegion cr = regions.get(index.intValue());
                float k = cr.getK()[cr.getKMaxIdx()];
                if (k > maxK) {
                    maxK = k;
                    maxKIndex = index;
                }
            }
            for (Integer index : set) {
                if (!index.equals(maxKIndex)) {
                    remove.add(index);
                }
            }
        }
        
        if (remove.size() > 1) {
            Collections.sort(remove);
        }
        for (int i = (remove.size() - 1); i > -1; --i) {
            regions.remove(remove.get(i).intValue());
        }
    }
    
    private static List<Set<Integer>> findCloseCorners(int tolD, 
        List<CornerRegion> regions) {
       
        float[] x = new float[regions.size()];
        float[] y = new float[regions.size()];
        for (int i = 0; i < regions.size(); ++i) {
            CornerRegion cr = regions.get(i);            
            x[i] = cr.getX()[cr.getKMaxIdx()];
            y[i] = cr.getY()[cr.getKMaxIdx()];
        }
        
        List<Set<Integer>> close = new ArrayList<Set<Integer>>();
        
        FixedDistanceGroupFinder groupFinder = new FixedDistanceGroupFinder(x, y);

        groupFinder.findGroupsOfPoints(tolD);
        
        int nGroups = groupFinder.getNumberOfGroups();
        
        for (int i = 0; i < nGroups; ++i) {
            Set<Integer> group = groupFinder.getGroupIndexes(i);
            if (group.size() > 1) {
                close.add(group);
            }
        }
        
        return close;
    }
  
    public TransformationParameters solveForScale(
        BlobCornerHelper img1Helper, IntensityFeatures features1,
        SegmentationType type1, boolean useBinned1,
        BlobCornerHelper img2Helper, IntensityFeatures features2,
        SegmentationType type2, boolean useBinned2,
        float[] outputScaleRotTransXYStDev) {

        GreyscaleImage img1 = img1Helper.imgHelper.getGreyscaleImage(useBinned1);
            
        GreyscaleImage img2 = img2Helper.imgHelper.getGreyscaleImage(useBinned2);

        List<List<CornerRegion>> corners1List = 
            img1Helper.getPerimeterCorners(type1, useBinned1);
        
        List<List<CornerRegion>> corners2List = 
            img2Helper.getPerimeterCorners(type2, useBinned2);
        
        List<Set<PairInt>> blobs1 = img1Helper.imgHelper.getBlobs(type1, useBinned1);
        List<Set<PairInt>> blobs2 = img2Helper.imgHelper.getBlobs(type2, useBinned2);
        List<PairIntArray> perimeters1 = img1Helper.imgHelper.getBlobPerimeters(type1, useBinned1);
        List<PairIntArray> perimeters2 = img2Helper.imgHelper.getBlobPerimeters(type2, useBinned2);
        
        int n1 = corners1List.size();
        int n2 = corners2List.size();
        
/*
MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
  
int idx1 = 14;
Integer index1 = Integer.valueOf(idx1);
Set<PairInt> blob1 = blobs1.get(idx1);
double[] xyCen1 = curveHelper.calculateXYCentroids(blob1);
List<CornerRegion> cr1 = new ArrayList<CornerRegion>(
    corners1List.get(Integer.valueOf(index1)));                
        
int idx2 = 13;
Integer index2 = Integer.valueOf(idx2);
Set<PairInt> blob2 = blobs2.get(idx2);
double[] xyCen2 = curveHelper.calculateXYCentroids(blob2);
List<CornerRegion> cr2 = new ArrayList<CornerRegion>(
    corners2List.get(Integer.valueOf(index2))); 

log.info("index1=" + index1.toString() + " index2=" + index2.toString()
+ " xyCen1=" + Arrays.toString(xyCen1) + " xyCen2=" + Arrays.toString(xyCen2));

try {
ImageExt img1C = img1.createColorGreyscaleExt();
ImageExt img2C = img2.createColorGreyscaleExt();
//ImageIOHelper.addCurveToImage(perimeter1, img1C, 0, 0, 0, 255);
ImageIOHelper.addAlternatingColorCornerRegionsToImage(cr1, img1C, 0, 0, 1);
           
//ImageIOHelper.addCurveToImage(perimeter2, img2C, 1, 0, 0, 255);
ImageIOHelper.addAlternatingColorCornerRegionsToImage(cr2, img2C, 0, 0, 1);

String bin = ResourceFinder.findDirectory("bin");
ImageIOHelper.writeOutputImage(bin + "/debug_corners_1.png", img1C);
ImageIOHelper.writeOutputImage(bin + "/debug_corners_2.png", img2C);

ImageIOHelper.writeLabeledCornerRegions(cr1, 0, 0, "/debug_labeled_corners_1.png");
ImageIOHelper.writeLabeledCornerRegions(cr2, 0, 0, "/debug_labeled_corners_2.png");

int z = 1;//1,3 in contours2 should be averaged?
} catch(IOException e) {
}
*/
        Map<Integer, TransformationPair4> trMap = new HashMap<Integer, TransformationPair4>();
        
        for (int i = 0; i < n1; ++i) {
            
            List<CornerRegion> corners1 = corners1List.get(i);
            int nc1 = corners1.size();
            if (nc1 < 3) {
                continue;
            }
            
            TransformationPair4 bestMatches = null;
            
            for (int j = 0; j < n2; ++j) {
                
                List<CornerRegion> corners2 = corners2List.get(j);
                int nc2 = corners2.size();
                if (nc2 < 3) {
                    continue;
                }
                
                int minN = Math.min(nc1, nc2);
                int diffNC = Math.abs(nc1 - nc2);
                
                if (((minN > 10) && (diffNC > 3)) || ((minN <= 10) && (diffNC > 2))){
                    continue;
                }
                
                List<CornerRegion> c1;
                List<CornerRegion> c2;
                
                if (nc1 == nc2) {
                    c1 = corners1;
                    c2 = corners2;
                } else {
                    c1 = new ArrayList<CornerRegion>(corners1);
                    c2 = new ArrayList<CornerRegion>(corners2);                  
                    filterToSameNumberIfPossible(perimeters1.get(i), c1,
                        perimeters2.get(j), c2);
                    
                    if (c1.size() != c2.size()) {
                        continue;
                    }
                }
                
                ClosedCurveCornerMatcher0 matcher 
                    = new ClosedCurveCornerMatcher0(features1, features2, c1, 
                        c2, true);
                
                boolean solved = matcher.matchCorners();
                
                if (solved) {
                    
                    if (bestMatches == null) {
                        
                        bestMatches = matcher.getTransformationPair();
                        bestMatches.setCornerListIndex1(i);
                        bestMatches.setCornerListIndex2(j);
                        
                    } else {
                        
                        TransformationPair4 tr = matcher.getTransformationPair();
                        tr.setCornerListIndex1(i);
                        tr.setCornerListIndex2(j);
                        
                        int bN = bestMatches.getMatchedCornerRegions1().size();
                        int cN = tr.getMatchedCornerRegions1().size();
                        
                        if (((cN >= bN) && (tr.getCost() < bestMatches.getCost()))
                            || ((tr.getCost() == bestMatches.getCost()) && (cN > bN))) {
                            
                            bestMatches = tr;
                        }
                    }
                }
            }
            
            if (bestMatches != null) {
                trMap.put(Integer.valueOf(i), bestMatches);
            }
        }
              
        if (trMap.isEmpty()) {
            return null;
        }
        
        int binFactor1 = img1Helper.imgHelper.getBinFactor(useBinned1);
        int binFactor2 = img2Helper.imgHelper.getBinFactor(useBinned2);
        
        // find best (lowest cost) and combine others with it if similar
        double minCost = Double.MAX_VALUE;
        Integer minCostKey = null;
        
        for (Entry<Integer, TransformationPair4> entry : trMap.entrySet()) {
            TransformationPair4 tr = entry.getValue();
            if (tr.getCost() < minCost) {
                minCost = tr.getCost();
                minCostKey = entry.getKey();
            }
        }
        
        TransformationPair4 minCostTR = trMap.remove(minCostKey);
                
        TransformationParameters minCostParams = calculateTransformation(
            binFactor1, binFactor2, minCostTR.getMatchedCompStats(), 
            outputScaleRotTransXYStDev);
        
        List<FeatureComparisonStat> combine = new ArrayList<FeatureComparisonStat>();
        
        for (Entry<Integer, TransformationPair4> entry : trMap.entrySet()) {
            
            TransformationPair4 tr = entry.getValue();
            
            if (rotationIsConsistent(minCostParams, tr.getMatchedCompStats(), 20)) {
                
                float[] outputScaleRotTransXYStDev00 = new float[4];
                TransformationParameters params = calculateTransformation(
                    binFactor1, binFactor2, tr.getMatchedCompStats(), 
                    outputScaleRotTransXYStDev00);
                
                if (params != null) {
                    if (areSimilar(minCostParams, params, 25)) {
                        combine.addAll(tr.getMatchedCompStats());
                    }
                }
            }
        }
        
        if (combine.isEmpty()) {
            return minCostParams;
        }
        
        combine.addAll(minCostTR.getMatchedCompStats());
        
        TransformationParameters combinedParams = calculateTransformation(
            binFactor1, binFactor2, combine, outputScaleRotTransXYStDev);
        
        return combinedParams;
    }

}
