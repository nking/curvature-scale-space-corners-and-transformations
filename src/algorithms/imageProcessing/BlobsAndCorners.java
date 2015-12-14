package algorithms.imageProcessing;

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
     * @param outdoorMode
     * @param enableJaggedLineCorrections
     * @param factorIncreaseForCurvatureMinimum
     * @param correctLineArtifacts
     * @return 
     */
    public static List<List<CornerRegion>> populateCorners(
        final BlobPerimeterHelper blobPerimeterHelper,
        SegmentationType type, boolean useBinnedImage,
        final boolean outdoorMode,
        final boolean enableJaggedLineCorrections,
        final float factorIncreaseForCurvatureMinimum,
        final boolean correctLineArtifacts) {

        List<List<CornerRegion>> cornerRegionLists = 
            new ArrayList<List<CornerRegion>>();
        
        List<PairIntArray> perimeterLists = blobPerimeterHelper.getBlobPerimeters(
            type, useBinnedImage);
        
        PairIntArray allCorners = new PairIntArray();
        
        int width = useBinnedImage ? 
            blobPerimeterHelper.getGreyscaleImageBinned().getWidth() :
            blobPerimeterHelper.getGreyscaleImage().getWidth();
        
        int height = useBinnedImage ? 
            blobPerimeterHelper.getGreyscaleImageBinned().getHeight() :
            blobPerimeterHelper.getGreyscaleImage().getHeight();
        
        Map<Integer, List<CornerRegion> > indexCornerRegionMap = 
            findCornersInScaleSpaceMaps(perimeterLists, outdoorMode, allCorners,
            enableJaggedLineCorrections, factorIncreaseForCurvatureMinimum,
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
                
                orderCornersCCW(curve, list);
                
                nTotCR += list.size();
                
                cornerRegionLists.add(list);
            }
        }
        
        if (correctLineArtifacts) {
            int binFactor = blobPerimeterHelper.getBinFactor();
            int thetaTol = 1;
            int radiusTol = 7;
            
           //use hough transform for lines to remove corners from line artifacts
           CornerCorrector.removeCornersFromLineArtifacts(perimeterLists, 
               cornerRegionLists, thetaTol, radiusTol, width, height);
        }
    
        if (blobPerimeterHelper.isInDebugMode()) {
            MiscDebug.writeEdgesAndCorners(perimeterLists, cornerRegionLists,
                1, blobPerimeterHelper.getGreyscaleImage(useBinnedImage),
                "blob_corners_" + blobPerimeterHelper.getDebugTag() + "_" 
                + MiscDebug.getCurrentTimeFormatted());
        }

        assert(perimeterLists.size() == cornerRegionLists.size());
        
        return cornerRegionLists;
    }
    
    public static void orderCornersCCW(PairIntArray curve, List<CornerRegion>
        list) {
        
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
     * @param doUseOutdoorMode
     * @param outputCorners
     * @param enableJaggedLineCorrections
     * @param factorIncreaseForCurvatureMinimum
     * @param imageWidth
     * @param imageHeight
     * @return scale space maps for each edge
     */
    protected static Map<Integer, List<CornerRegion> >
        findCornersInScaleSpaceMaps(final List<PairIntArray> theEdges, 
        final boolean doUseOutdoorMode, final PairIntArray outputCorners,
        final boolean enableJaggedLineCorrections,
        final float factorIncreaseForCurvatureMinimum, int imageWidth,
        int imageHeight) {

        CSSCornerMaker cornerMaker = new CSSCornerMaker(imageWidth, imageHeight);
        
        if (enableJaggedLineCorrections) {
            cornerMaker.enableJaggedLineCorrections();
        } else {
            cornerMaker.disableJaggedLineCorrections();
        }
       
        cornerMaker.increaseFactorForCurvatureMinimum(
            factorIncreaseForCurvatureMinimum);
        
        // an empty map to allow re-use of method which uses junctions if given
        Map<Integer, Set<PairIntWithIndex>> noJunctions = 
            new HashMap<Integer, Set<PairIntWithIndex>>();
        
        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > sMap =
            cornerMaker.findCornersInScaleSpaceMaps(theEdges, noJunctions,
                doUseOutdoorMode, outputCorners); 
        
        return cornerMaker.getEdgeCornerRegionMap();        
    }

}
