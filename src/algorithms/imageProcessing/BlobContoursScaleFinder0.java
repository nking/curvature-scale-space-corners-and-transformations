package algorithms.imageProcessing;

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
 * and image2 using curvature scale space contours and an assumption of
 * ordered matching between the peaks.
 * The method is fast, but needs blobs with concave and convex boundaries
 * to produce strong inflection points.
 * 
 * @author nichole
 */
public class BlobContoursScaleFinder0 extends AbstractBlobScaleFinder {

    protected void filterToSameNumberIfPossible(
        PairIntArray curve1, List<CurvatureScaleSpaceContour> contours1, 
        List<BlobPerimeterRegion> regions1,
        PairIntArray curve2, List<CurvatureScaleSpaceContour> contours2,
        List<BlobPerimeterRegion> regions2) {
        
        //TODO: this may change to improve results for different resolution
        // data for example... may need to have increasing tolerances
        // that are different for contours1 than contours2
        
        // the curves have already had tolT=0.04 and tolD=6 applied to them
        float tolT = 0.055f; 
        float tolD = 9.5f;
        
        List<Integer> rmIndexes = BlobsAndContours.combineOverlappingPeaks(
            tolT, tolD, curve1, contours1);
        
        if (rmIndexes != null && !rmIndexes.isEmpty()) {
            for (int i = (rmIndexes.size() - 1); i > -1; --i) {
                regions1.remove(rmIndexes.get(i).intValue());
            }
        }
        
        rmIndexes = BlobsAndContours.combineOverlappingPeaks(tolT, tolD, curve2, 
            contours2);
       
        if (rmIndexes != null && !rmIndexes.isEmpty()) {
            for (int i = (rmIndexes.size() - 1); i > -1; --i) {
                regions2.remove(rmIndexes.get(i).intValue());
            }
        }
    }
    
    protected List<CurvatureScaleSpaceContour> copyContours(
        List<CurvatureScaleSpaceContour> contours) {
        
        List<CurvatureScaleSpaceContour> c = new ArrayList<CurvatureScaleSpaceContour>();
        
        for (int i = 0; i < contours.size(); ++i) {
            c.add(contours.get(i).copy());
        }
        
        return c;
    }
    
    public TransformationParameters solveForScale(
        BlobContourHelper img1Helper, IntensityFeatures features1,
        SegmentationType type1, boolean useBinned1,
        BlobContourHelper img2Helper, IntensityFeatures features2,
        SegmentationType type2, boolean useBinned2,
        float[] outputScaleRotTransXYStDev) {

        List<List<CurvatureScaleSpaceContour>> contours1List = 
            img1Helper.getPerimeterContours(type1, useBinned1);
        
        List<List<CurvatureScaleSpaceContour>> contours2List = 
            img2Helper.getPerimeterContours(type2, useBinned2);
        
        List<Set<PairInt>> blobs1 = img1Helper.imgHelper.getBlobs(type1, useBinned1);
        List<Set<PairInt>> blobs2 = img2Helper.imgHelper.getBlobs(type2, useBinned2);
        List<PairIntArray> perimeters1 = img1Helper.imgHelper.getBlobPerimeters(type1, useBinned1);
        List<PairIntArray> perimeters2 = img2Helper.imgHelper.getBlobPerimeters(type2, useBinned2);
        
        GreyscaleImage img1 = img1Helper.imgHelper.getGreyscaleImage(useBinned1);
        GreyscaleImage img2 = img2Helper.imgHelper.getGreyscaleImage(useBinned2);
        
        Map<Integer, List<BlobPerimeterRegion>> contours1PointMaps = 
            new HashMap<Integer, List<BlobPerimeterRegion>>();

        Map<Integer, List<BlobPerimeterRegion>> contours2PointMaps = 
            new HashMap<Integer, List<BlobPerimeterRegion>>();
         
        int n1 = contours1List.size();
        int n2 = contours2List.size();
        
        Map<Integer, TransformationPair3> trMap = new HashMap<Integer, TransformationPair3>();
        
        for (int i = 0; i < n1; ++i) {
            
            List<CurvatureScaleSpaceContour> contour1 = contours1List.get(i);
            int nc1 = contour1.size();
            if (nc1 < 3) {
                continue;
            }
            
            List<BlobPerimeterRegion> contour1Points = contours1PointMaps.get(Integer.valueOf(i));
            if (contour1Points == null) {
                contour1Points = extractBlobPerimeterRegions(i, contour1,
                    perimeters1.get(i), blobs1.get(i));
                contours1PointMaps.put(Integer.valueOf(i), contour1Points);
            }
            
            TransformationPair3 bestMatches = null;
            
            for (int j = 0; j < n2; ++j) {
                
                List<CurvatureScaleSpaceContour> contour2 = contours2List.get(j);
                int nc2 = contour2.size();
                if (nc2 < 3) {
                    continue;
                }
                
                int minN = Math.min(nc1, nc2);
                int diffNC = Math.abs(nc1 - nc2);
                
                if (((minN > 10) && (diffNC > 3)) || ((minN <= 10) && (diffNC > 2))){
                    continue;
                }
                
                List<BlobPerimeterRegion> contour2Points = contours2PointMaps.get(Integer.valueOf(j));
                if (contour2Points == null) {
                    contour2Points = extractBlobPerimeterRegions(j, contour2,
                        perimeters2.get(j), blobs2.get(j));
                    contours2PointMaps.put(Integer.valueOf(j), contour2Points);
                }
            
                List<CurvatureScaleSpaceContour> c1;
                List<CurvatureScaleSpaceContour> c2;
                List<BlobPerimeterRegion> c1p;
                List<BlobPerimeterRegion> c2p;
                
                if (nc1 == nc2) {
                    c1 = contour1;
                    c2 = contour2;
                    c1p = contour1Points;
                    c2p = contour2Points;
                } else {
                    c1 = copyContours(contour1);
                    c2 = copyContours(contour2);
                    c1p = copyRegions(contour1Points);
                    c2p = copyRegions(contour2Points);                  
                    filterToSameNumberIfPossible(perimeters1.get(i), c1, c1p,
                        perimeters2.get(j), c2, c2p);
                    
                    if (c1.size() != c2.size()) {
                        continue;
                    }
                }
            
                //matching algorithm here for c1p and c2p
                ClosedCurveContourMatcher0 matcher 
                    = new ClosedCurveContourMatcher0();
                
                boolean solved = matcher.matchCorners(features1, features2, 
                    c1p, c2p, img1, img2);
                
                if (solved) {
                    
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
            }
        }
        
 /*
 MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
  
int idx1 = 0;
Integer index1 = Integer.valueOf(idx1);
Set<PairInt> blob1 = blobs1.get(idx1);
double[] xyCen1 = curveHelper.calculateXYCentroids(blob1);
List<BlobPerimeterRegion> cr1 = new ArrayList<BlobPerimeterRegion>(
    contours1PointMaps.get(Integer.valueOf(index1)));                
ClosedCurveContourMatcher0.sortRegion(cr1);
        
int idx2 = 5;
Integer index2 = Integer.valueOf(idx2);
Set<PairInt> blob2 = blobs2.get(idx2);
double[] xyCen2 = curveHelper.calculateXYCentroids(blob2);
List<BlobPerimeterRegion> cr2 = new ArrayList<BlobPerimeterRegion>(
    contours2PointMaps.get(Integer.valueOf(index2)));                
ClosedCurveContourMatcher0.sortRegion(cr2);

log.info("index1=" + index1.toString() + " index2=" + index2.toString()
+ " xyCen1=" + Arrays.toString(xyCen1) + " xyCen2=" + Arrays.toString(xyCen2));

try {
ImageExt img1C = img1.createColorGreyscaleExt();
ImageExt img2C = img2.createColorGreyscaleExt();
//ImageIOHelper.addCurveToImage(perimeter1, img1C, 0, 0, 0, 255);
ImageIOHelper.addAlternatingColorRegionsToImage(cr1, img1C, 0, 0, 1);
           
//ImageIOHelper.addCurveToImage(perimeter2, img2C, 1, 0, 0, 255);
ImageIOHelper.addAlternatingColorRegionsToImage(cr2, img2C, 0, 0, 1);

String bin = ResourceFinder.findDirectory("bin");
ImageIOHelper.writeOutputImage(bin + "/debug_contours_1.png", img1C);
ImageIOHelper.writeOutputImage(bin + "/debug_contours_2.png", img2C);

ImageIOHelper.writeLabeledRegions(cr1, 0, 0, "/debug_labeled_contours_1.png");
ImageIOHelper.writeLabeledRegions(cr2, 0, 0, "/debug_labeled_contours_2.png");

int z = 1;//1,3 in contours2 should be averaged?
} catch(IOException e) {
} */       
        if (trMap.isEmpty()) {
            return null;
        }
        
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
        
        int binFactor1 = img1Helper.imgHelper.getBinFactor(useBinned1);
        int binFactor2 = img2Helper.imgHelper.getBinFactor(useBinned2);
        
        TransformationPair3 minCostTR = trMap.remove(minCostKey);
                
        TransformationParameters minCostParams = calculateTransformation(
            binFactor1, binFactor2, minCostTR.getMatchedCompStats(), 
            outputScaleRotTransXYStDev);
        
        List<FeatureComparisonStat> combine = new ArrayList<FeatureComparisonStat>();
        
        for (Entry<Integer, TransformationPair3> entry : trMap.entrySet()) {
            
            TransformationPair3 tr = entry.getValue();
            
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

    private List<BlobPerimeterRegion> copyRegions(List<BlobPerimeterRegion> 
        contourPoints) {
        
        List<BlobPerimeterRegion> copy = new ArrayList<BlobPerimeterRegion>();
        
        for (int i = 0; i < contourPoints.size(); ++i) {
            copy.add(contourPoints.get(i).copy());
        }
        
        return copy;
    }

}
