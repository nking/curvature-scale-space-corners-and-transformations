package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class BlobContoursScaleFinder0 extends AbstractBlobScaleFinder {

    protected void filterToSameNumberIfPossible(
        PairIntArray curve1, List<CurvatureScaleSpaceContour> contours1, 
        PairIntArray curve2, List<CurvatureScaleSpaceContour> contours2) {
        
        //TODO: this may change to improve results for different resolution
        // data for example... may need to have increasing tolerances
        // that are different for contours1 than contours2
        
        // the curves have already had tolT=0.04 and tolD=6 applied to them
        float tolT = 0.055f; 
        float tolD = 9.5f;
        
        BlobsAndContours.combineOverlappingPeaks(tolT, tolD, curve1, contours1);
        
        BlobsAndContours.combineOverlappingPeaks(tolT, tolD, curve2, contours2);
       
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
        
        int idx2 = 3;
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

int z = 1;//1,3 in contours2 should be averaged?
} catch(IOException e) {
}
       
        int n1 = contours1List.size();
        int n2 = contours2List.size();
        
        for (int i = 0; i < n1; ++i) {
            
            List<CurvatureScaleSpaceContour> contour1 = contours1List.get(i);
            int nc1 = contour1.size();
            if (nc1 < 3) {
                continue;
            }
            
            for (int j = 0; j < n2; ++j) {
                
                List<CurvatureScaleSpaceContour> contour2 = contours2List.get(j);
                int nc2 = contour2.size();
                if (nc2 < 3) {
                    continue;
                }
                
                if (Math.abs(nc1 - nc2) > 2) {
                    continue;
                }
                List<CurvatureScaleSpaceContour> c1;
                List<CurvatureScaleSpaceContour> c2;
                if (nc1 == nc2) {
                    c1 = contour1;
                    c2 = contour2;
                } else {
                    c1 = copyContours(contour1);
                    c2 = copyContours(contour2);
                    filterToSameNumberIfPossible(perimeters1.get(i), c1, 
                        perimeters2.get(j), c2);
                    
                    if (c1.size() != c2.size()) {
                        continue;
                    }
                }
                
                //matching algorithm here for c1 and 
            }
        }

        throw new UnsupportedOperationException("not yet implemented");

    }

}
