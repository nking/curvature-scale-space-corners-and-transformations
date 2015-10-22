package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairFloat;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

/**
 *
 * @author nichole
 */
public class BlobContoursScaleFinder0 extends AbstractBlobScaleFinder {

    protected List<List<CurvatureScaleSpaceContour>> filteredToSameNumber(
        List<CurvatureScaleSpaceContour> contour1,
        List<CurvatureScaleSpaceContour> contour2) {
        
        /*
        NOTE: this processing can happen at a higher level to cache
        for all contours.
        
        averaging very close points by scale-free length (and x,y):
        -- average to one point if in t-space their locations are within 
           delta t=0.025 of each other (and close in x,y space, dist ~5 or 6)
        */
        
        List<List<CurvatureScaleSpaceContour>> filtered = 
            new ArrayList<List<CurvatureScaleSpaceContour>>(2);
        
        throw new UnsupportedOperationException("not yet implemented");
        //return filtered;
    }
    
    public TransformationParameters solveForScale(
        BlobContourHelper img1Helper, IntensityFeatures features1,
        SegmentationType type1, boolean useBinned1,
        BlobContourHelper img2Helper, IntensityFeatures features2,
        SegmentationType type2, boolean useBinned2,
        float[] outputScaleRotTransXYStDev) {

        GreyscaleImage img1 = img1Helper.imgHelper.getGreyscaleImage(useBinned1);
            
        GreyscaleImage img2 = img2Helper.imgHelper.getGreyscaleImage(useBinned2);

        List<List<CurvatureScaleSpaceContour>> contours1List = 
            img1Helper.getPerimeterContours(type1, useBinned1);
        
        List<List<CurvatureScaleSpaceContour>> contours2List = 
            img2Helper.getPerimeterContours(type2, useBinned2);
        
        List<Set<PairInt>> blobs1 = img1Helper.imgHelper.getBlobs(type1, useBinned1);
        List<Set<PairInt>> blobs2 = img2Helper.imgHelper.getBlobs(type2, useBinned2);
        List<PairIntArray> perimeters1 = img1Helper.imgHelper.getBlobPerimeters(type1, useBinned1);
        List<PairIntArray> perimeters2 = img2Helper.imgHelper.getBlobPerimeters(type2, useBinned2);

       
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        // plot the contours
        // or specific contours to look at how best to filter to common
        // number of members
        
        //exploring first before further specification
        
        int idx1 = 0;
        Integer index1 = Integer.valueOf(idx1);
        PairIntArray perimeter1 = perimeters1.get(idx1);
        Set<PairInt> blob1 = blobs1.get(idx1);
        List<CurvatureScaleSpaceContour> contours1 = contours1List.get(idx1);
        double[] xyCen1 = curveHelper.calculateXYCentroids(perimeter1);
        
        int idx2 = 5;
        Integer index2 = Integer.valueOf(idx2);
        PairIntArray perimeter2 = perimeters2.get(idx2);
        Set<PairInt> blob2 = blobs2.get(idx2);
        List<CurvatureScaleSpaceContour> contours2 = contours2List.get(idx2);
        double[] xyCen2 = curveHelper.calculateXYCentroids(perimeter2);
        
log.info("index1=" + index1.toString() + " index2=" + index2.toString()
+ " xyCen1=" + Arrays.toString(xyCen1) + " xyCen2=" + Arrays.toString(xyCen2));

log.info(
String.format("[%d](%d,%d) [%d](%d,%d)  nCurvePoints=%d, %d", 
idx1, (int)Math.round(xyCen1[0]), (int)Math.round(xyCen1[1]),
idx2, (int)Math.round(xyCen2[0]), (int)Math.round(xyCen2[1]),
perimeter1.getN(), perimeter2.getN()));               

try {
ImageExt img1C = img1.createColorGreyscaleExt();
ImageExt img2C = img2.createColorGreyscaleExt();
//ImageIOHelper.addCurveToImage(perimeter1, img1C, 0, 0, 0, 255);
ImageIOHelper.addAlternatingColorContoursToImage(contours1, img1C, 0, 0, 1);
           
//ImageIOHelper.addCurveToImage(perimeter2, img2C, 1, 0, 0, 255);
ImageIOHelper.addAlternatingColorContoursToImage(contours2, img2C, 0, 0, 1);

String bin = ResourceFinder.findDirectory("bin");
ImageIOHelper.writeOutputImage(bin + "/debug_contours_1.png", img1C);
ImageIOHelper.writeOutputImage(bin + "/debug_contours_2.png", img2C);

ImageIOHelper.writeLabeledContours(perimeter1, contours1, 0, 0, 
    "/debug_labeled_contours_1.png");
ImageIOHelper.writeLabeledContours(perimeter2, contours2, 0, 0,
    "/debug_labeled_contours_2.png");

int z = 1;
} catch(IOException e) {
}

        /* this just needs to use the O(N) pattern of matching
         0 with 0, 0 with 1, etc
        and use the MatchedPointsTransformationCalculator method for
        matched points to remove outliers make pairs of the matched points 
        and solve the transformation
        */

        throw new UnsupportedOperationException("not yet implemented");

    }

}
