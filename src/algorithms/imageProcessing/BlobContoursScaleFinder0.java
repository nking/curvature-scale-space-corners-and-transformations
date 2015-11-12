package algorithms.imageProcessing;

import algorithms.imageProcessing.util.MiscStats;
import algorithms.misc.MiscDebug;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * class to attempt to find Euclidean scale transformation between image1
 * and image2 using curvature scale space contours and an assumption of
 * ordered matching between the peaks.
 * The method is fast, but needs blobs with concave and convex boundaries
 * to produce strong inflection points.
 * 
 * @author nichole
 */
public class BlobContoursScaleFinder0 extends AbstractBlobScaleFinder {

    protected void filterToSameNumberIfPossible(
        PairIntArray curve1, List<BlobPerimeterRegion> regions1,
        PairIntArray curve2, List<BlobPerimeterRegion> regions2) {
        
        //TODO: this may change to improve results for different resolution
        // data for example... may need to have increasing tolerances
        // that are different for contours1 than contours2
        
        int tolT = 3; 
        float tolD = 4.0f;
        
        List<Integer> rmIndexes = BlobsAndContours.combineOverlappingPeaks0(
            tolT, tolD, curve1, regions1);
        
        rmIndexes = BlobsAndContours.combineOverlappingPeaks0(tolT, tolD, curve2, 
            regions2);
       
    }
    
    protected List<CurvatureScaleSpaceContour> copyContours(
        List<CurvatureScaleSpaceContour> contours) {
        
        List<CurvatureScaleSpaceContour> c = new ArrayList<CurvatureScaleSpaceContour>();
        
        for (int i = 0; i < contours.size(); ++i) {
            c.add(contours.get(i).copy());
        }
        
        return c;
    }
    
    public MatchingSolution solveForScale(
        BlobContourHelper img1Helper, IntensityFeatures features1,
        SegmentationType type1, boolean useBinned1,
        BlobContourHelper img2Helper, IntensityFeatures features2,
        SegmentationType type2, boolean useBinned2) {

        List<List<BlobPerimeterRegion>> bpr1List = 
            img1Helper.getPerimeterRegions(type1, useBinned1);
        
        List<List<BlobPerimeterRegion>> bpr2List = 
            img2Helper.getPerimeterRegions(type2, useBinned2);
        
        List<PairIntArray> perimeters1 = img1Helper.imgHelper.getBlobPerimeters(type1, useBinned1);
        List<PairIntArray> perimeters2 = img2Helper.imgHelper.getBlobPerimeters(type2, useBinned2);
        
        GreyscaleImage img1 = img1Helper.imgHelper.getGreyscaleImage(useBinned1);
        GreyscaleImage img2 = img2Helper.imgHelper.getGreyscaleImage(useBinned2);
        
        int binFactor1 = img1Helper.imgHelper.getBinFactor(useBinned1);
        int binFactor2 = img2Helper.imgHelper.getBinFactor(useBinned2);
        
        if (debug) {
            Collection<BlobPerimeterRegion> set1 = new ArrayList<BlobPerimeterRegion>();
            Collection<BlobPerimeterRegion> set2 = new ArrayList<BlobPerimeterRegion>();
            for (int i = 0; i < bpr1List.size(); ++i) {
                set1.addAll(bpr1List.get(i));
            }
            for (int i = 0; i < bpr2List.size(); ++i) {
                set2.addAll(bpr2List.get(i));
            }
            try {
                MiscDebug.writeImage(set1, img1.copyToColorGreyscale(), "contour_regions_1");
                MiscDebug.writeImage(set2, img2.copyToColorGreyscale(), "contour_regions_2");
            } catch (IOException ex) {
                Logger.getLogger(BlobContoursScaleFinder0.class.getName()).log(Level.SEVERE, null, ex);
            }            
        }
        
        int n1 = bpr1List.size();
        int n2 = bpr2List.size();
        
/*        
MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
float[] xPoints1 = new float[perimeters1.size()];
float[] yPoints1 = new float[perimeters1.size()];
double[][] xy1 = new double[perimeters1.size()][2];
for (int i = 0; i < perimeters1.size(); ++i) {
xy1[i] = curveHelper.calculateXYCentroids(perimeters1.get(i));
xPoints1[i] = (float)xy1[i][0];
yPoints1[i] = (float)xy1[i][1];
}
float[] xPoints2 = new float[perimeters2.size()];
float[] yPoints2 = new float[perimeters2.size()];
double[][] xy2 = new double[perimeters2.size()][2];
for (int i = 0; i < perimeters2.size(); ++i) {
xy2[i] = curveHelper.calculateXYCentroids(perimeters2.get(i));
xPoints2[i] = (float)xy2[i][0];
yPoints2[i] = (float)xy2[i][1];
}
ScatterPointPlotterPNG plotter = new ScatterPointPlotterPNG();
plotter.plotLabeledPoints(0, img1.getWidth(), 0, img1.getHeight(), xPoints1, yPoints1,
"img1", "X", "Y");
ScatterPointPlotterPNG plotter2 = new ScatterPointPlotterPNG();
plotter2.plotLabeledPoints(0, img2.getWidth(), 0, img2.getHeight(), xPoints2, yPoints2,
"img2", "X", "Y");
try {
    plotter.writeToFile("img1_labelled.png");
    plotter2.writeToFile("img2_labelled.png");
} catch (IOException ex) {
    Logger.getLogger(BlobContoursScaleFinder0.class.getName()).log(Level.SEVERE, null, ex);
}
*/
        Map<Integer, TransformationPair3> trMap = new HashMap<Integer, TransformationPair3>();
        
        Map<Integer, TransformationParameters> pMap = new HashMap<Integer,
            TransformationParameters>();
        
        for (int i = 0; i < n1; ++i) {
            
            List<BlobPerimeterRegion> bpr1 = bpr1List.get(i);
            int nc1 = bpr1.size();
            if (nc1 < 3) {
                continue;
            }
            
            TransformationPair3 bestMatches = null;
            
            for (int j = 0; j < n2; ++j) {
                
                List<BlobPerimeterRegion> bpr2 = bpr2List.get(j);
                int nc2 = bpr2.size();
                if (nc2 < 3) {
                    continue;
                }
                
                int minN = Math.min(nc1, nc2);
                int diffNC = Math.abs(nc1 - nc2);
                
                if (((minN > 10) && (diffNC > 3)) || ((minN <= 10) && (diffNC > 2))){
                    continue;
                }
                
                List<BlobPerimeterRegion> c1p;
                List<BlobPerimeterRegion> c2p;
                
                if (nc1 == nc2) {
                    c1p = bpr1;
                    c2p = bpr2;
                } else {
                    c1p = copyRegions(bpr1);
                    c2p = copyRegions(bpr2);                  
                    filterToSameNumberIfPossible(perimeters1.get(i), c1p,
                        perimeters2.get(j), c2p);
                    
                    if (c1p.size() != c2p.size()) {
                        continue;
                    }
                }
            
                //matching algorithm here for c1p and c2p
                ClosedCurveContourMatcher0 matcher 
                    = new ClosedCurveContourMatcher0();
                
                boolean solved = matcher.matchCorners(features1, features2, 
                    c1p, c2p, img1, img2);
                
                if (solved) {

                    assert(matcher.getTransformationPair3() != null);

                    if (bestMatches == null) {
                        
                        bestMatches = matcher.getTransformationPair3();
                        bestMatches.setContourIndex1(i);
                        bestMatches.setContourIndex2(j);
                        
                    } else {
                        
                        TransformationPair3 tr = matcher.getTransformationPair3();
                        tr.setContourIndex1(i);
                        tr.setContourIndex2(j);
                        
                        int bN = bestMatches.getMatchedContourRegions1().size();
                        int cN = tr.getMatchedContourRegions1().size();
                        
                        if (((cN >= bN) && (tr.getCost() < bestMatches.getCost()))
                            || ((tr.getCost() == bestMatches.getCost()) && (cN > bN))) {
                            
                            bestMatches = tr;
                        }
                    }
                }
            }
            
            if (bestMatches != null) {
                
                trMap.put(Integer.valueOf(i), bestMatches);
                
                TransformationParameters params = calculateTransformation(
                    binFactor1, binFactor2, bestMatches.getMatchedCompStats(), 
                    new float[4]);
                
                if (params == null) {
                    continue;
                }
                
                if ((bestMatches.getMatchedContourRegions1().size() > 4) && 
                    (i < (n1 - 1))) {
                    if (MiscStats.standardDeviationsAreSmall(params)) {
                        
                        double c = calculateCombinedIntensityStat(
                            bestMatches.getMatchedCompStats());
                        log.info("MATCHED EARLY: combined compStat=" + c);
                        
                        MatchingSolution soln = new MatchingSolution(params,
                            bestMatches.getMatchedCompStats());
                        
                        return soln;
                    }
                }
                
                pMap.put(Integer.valueOf(i), params);
            }
        }
        
        if (trMap.isEmpty()) {
            return null;
        }
        
//TODO: consider changing to the biparite matching then most frequent params filter
        
        // find best (lowest cost) and combine others with it if similar
        double minCost = Double.MAX_VALUE;
        Integer minCostKey = null;
        
        for (Entry<Integer, TransformationPair3> entry : trMap.entrySet()) {
            TransformationPair3 tr = entry.getValue();
            if (tr.getCost() < minCost) {
                minCost = tr.getCost();
                minCostKey = entry.getKey();
            }
        }
        
        assert(minCostKey != null);
        
        TransformationPair3 minCostTR = trMap.remove(minCostKey);
                
        TransformationParameters minCostParams = calculateTransformation(
            binFactor1, binFactor2, minCostTR.getMatchedCompStats(), 
            new float[4]);
        
        List<FeatureComparisonStat> combine = new ArrayList<FeatureComparisonStat>();
    
        if (minCostParams != null) {
            for (Entry<Integer, TransformationPair3> entry : trMap.entrySet()) {

                TransformationPair3 tr = entry.getValue();

                if (rotationIsConsistent(minCostParams, tr.getMatchedCompStats(), 20)) {

                    TransformationParameters params = pMap.get(entry.getKey());

                    if (params != null) {
                        if (areSimilar(minCostParams, params, 25)) {
                            combine.addAll(tr.getMatchedCompStats());
                        }
                    }
                }
            }
        }
        
        if (combine.isEmpty()) {

            if (minCostParams == null) {
                return null;
            }
            
            MatchingSolution soln = new MatchingSolution(minCostParams,
                minCostTR.getMatchedCompStats());
            
            return soln;
        }
        
        combine.addAll(minCostTR.getMatchedCompStats());
        
        TransformationParameters combinedParams = calculateTransformation(
            binFactor1, binFactor2, combine, new float[4]);
        
        if (combinedParams == null) {
            return null;
        }
        
        MatchingSolution soln = new MatchingSolution(combinedParams,
            combine);
        
        return soln;
    }

    private List<BlobPerimeterRegion> copyRegions(List<BlobPerimeterRegion> 
        contourPoints) {
        
        List<BlobPerimeterRegion> copy = new ArrayList<BlobPerimeterRegion>();
        
        for (int i = 0; i < contourPoints.size(); ++i) {
            copy.add(contourPoints.get(i).copy());
        }
        
        return copy;
    }

}
