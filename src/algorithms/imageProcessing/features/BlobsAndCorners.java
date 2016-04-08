package algorithms.imageProcessing.features;

import algorithms.imageProcessing.*;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.scaleSpace.CSSCornerMaker;
import algorithms.imageProcessing.scaleSpace.ScaleSpaceCurve;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.MiscDebug;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * holds the blobs, their ordered perimeters, and the scale space image
 * contours of the perimeters
 *
 * @author nichole
 */
public class BlobsAndCorners  {
        
    /**
     * get a list of corner regions for the perimeters extracted from the
     * blobs.  (Note that the default construction of perimeters already
     * orders them in counter clockwise manner.  That's necessary for
     * consistent use of the corner region orientation.)
     * @param blobPerimeterHelper
     * @param type
     * @param useBinnedImage
     * @param factorIncreaseForCurvatureMinimum
     * @return 
     */
    public static List<List<CornerRegion>> populateCorners(
        final BlobPerimeterCornerHelper blobPerimeterHelper,
        SegmentationType type, boolean useBinnedImage,
        final float factorIncreaseForCurvatureMinimum) {

        List<PairIntArray> perimeterLists = blobPerimeterHelper.getBlobPerimeters(
            type, useBinnedImage);
        
        PairIntArray allCorners = new PairIntArray();
        
        int width = useBinnedImage ? 
            blobPerimeterHelper.getGreyscaleImageBinned().getWidth() :
            blobPerimeterHelper.getGreyscaleImage().getWidth();
        
        int height = useBinnedImage ? 
            blobPerimeterHelper.getGreyscaleImageBinned().getHeight() :
            blobPerimeterHelper.getGreyscaleImage().getHeight();
        
        List<List<CornerRegion>> cornerRegionLists = populateCorners(
            perimeterLists, useBinnedImage,
            factorIncreaseForCurvatureMinimum,
            width, height);
        
        if (blobPerimeterHelper.isInDebugMode()) {
            MiscDebug.writeEdgesAndCorners(perimeterLists, cornerRegionLists,
                1, blobPerimeterHelper.getGreyscaleImage(useBinnedImage),
                "blob_corners_" + blobPerimeterHelper.getDebugTag() + "_" 
                + MiscDebug.getCurrentTimeFormatted());
        }
        
        return cornerRegionLists;
    }
    
    /**
     * get a list of corner regions for the perimeters extracted from the
     * blobs.  (Note that the default construction of perimeters already
     * orders them in counter clockwise manner.  That's necessary for
     * consistent use of the corner region orientation.)
     * @param perimeterLists
     * @param useBinnedImage
     * @param factorIncreaseForCurvatureMinimum
     * @param width
     * @param height
     * @return 
     */
    public static List<List<CornerRegion>> populateCorners(
        List<PairIntArray> perimeterLists, boolean useBinnedImage,
        final float factorIncreaseForCurvatureMinimum,
        int width, int height) {

        List<List<CornerRegion>> cornerRegionLists = 
            new ArrayList<List<CornerRegion>>();
        
        PairIntArray allCorners = new PairIntArray();
        
        Map<Integer, List<CornerRegion> > indexCornerRegionMap = 
            findCornersInScaleSpaceMaps(perimeterLists, allCorners,
            factorIncreaseForCurvatureMinimum,
            width, height);
        
        int nTotCR = 0;
        
        for (int i = 0; i < perimeterLists.size(); ++i) {
            List<CornerRegion> list = indexCornerRegionMap.get(Integer.valueOf(i));
            if (list == null) {
                // have to keep ordered, parallel indexes
                cornerRegionLists.add(new ArrayList<CornerRegion>());
            } else {
                
                // detailed check for ordering the perimeters.
                
                // TODO: consider corrections to this for perimeter also at higher level
                // then it's not needed here
                
                PairIntArray curve = perimeterLists.get(i);
                
                // cannot order if no edge points were present for the curveIndexes
                if (curve.getN() > 0) {
                    
                    orderCornersCCW(curve, list);
                
                    nTotCR += list.size();
                }
                
                cornerRegionLists.add(list);
            }
        }
        
        assert(perimeterLists.size() == cornerRegionLists.size());
        
        return cornerRegionLists;
    }
    
    public static void orderCornersCCW(PairIntArray curve, List<CornerRegion>
        list) {
        
        // cannot order corners if curveIndexes are missing from CornerRegions
        if (curve.getN() == 0) {
            return;
        }
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
                
        PairIntArray tmpCornerChk = new PairIntArray();

        for (int ii = 0; ii < list.size(); ++ii) {
            CornerRegion cr = list.get(ii);
            if (ii == 0 && cr.getIndexWithinCurve() != 0) {
                // adding test point for clockwise calculation
                int idx = cr.getIndexWithinCurve()/2;
                tmpCornerChk.add(curve.getX(idx), curve.getY(idx));
            }
            tmpCornerChk.add(cr.getX()[cr.getKMaxIdx()], cr.getY()[cr.getKMaxIdx()]);
        }

        if (list.get(list.size() - 1).getIndexWithinCurve() < (curve.getN() - 1)) {
            // adding test point for clockwise calculation
            int idx = curve.getN() - 1;
            tmpCornerChk.add(curve.getX(idx), curve.getY(idx));
        }
        boolean isCW = curveHelper.curveIsOrderedClockwise(tmpCornerChk);
        if (isCW) {
            int n = list.size();
            if (n > 1) {
                int end = n >> 1;
                for (int ii = 0; ii < end; ii++) {
                    int idx2 = n - ii - 1;
                    CornerRegion swap = list.get(ii);
                    list.set(ii, list.get(idx2));
                    list.set(idx2, swap);
                }
            }
        }
    }

    /**
     * Find corners in the edges by creating scale space maps for each edge that
     * range from a maximum value determined in getMaxSIGMAForECSS() and
     * determine the corners in the highest sigma of those maps and refine the
     * corner locations in the smaller sigma maps.  The corners are found using
     * the curvature minima and maxima points in the curve.  A lower threshold
     * is determined and used during the maxima finding and minima comparison.
     * Each corner candidate is larger than one of the adjacent minima by
     * a factor such as 2 or 3.
     * The results are set in the instance variable corners.
     * The returned variable is the set of scale space maps which might be
     * useful for other purposes, but are no longer needed for the corner
     * determination.
     *
     * @param theEdges
     * @param outputCorners
     * @param factorIncreaseForCurvatureMinimum
     * @param imageWidth
     * @param imageHeight
     * @return scale space maps for each edge
     */
    protected static Map<Integer, List<CornerRegion> >
        findCornersInScaleSpaceMaps(final List<PairIntArray> theEdges, 
        final PairIntArray outputCorners,
        final float factorIncreaseForCurvatureMinimum, int imageWidth,
        int imageHeight) {

        CSSCornerMaker cornerMaker = new CSSCornerMaker(imageWidth, imageHeight);
        
        cornerMaker.increaseFactorForCurvatureMinimum(
            factorIncreaseForCurvatureMinimum);
        
        // an empty map to allow re-use of method which uses junctions if given
        Map<Integer, Set<PairIntWithIndex>> noJunctions = 
            new HashMap<Integer, Set<PairIntWithIndex>>();
        
        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > sMap =
            cornerMaker.findCornersInScaleSpaceMaps(theEdges, noJunctions,
                outputCorners); 
        
        return cornerMaker.getEdgeCornerRegionMap();        
    }

}
