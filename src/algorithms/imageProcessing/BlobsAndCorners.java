package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * holds the blobs, their ordered perimeters, and the scale space image
 * contours of the perimeters
 *
 * @author nichole
 */
public class BlobsAndCorners  {
    
    public static List<List<CornerRegion>> populateCorners(
        final BlobPerimeterHelper blobPerimeterHelper,
        SegmentationType type, boolean useBinnedImage,
        final boolean outdoorMode,
        final boolean enableJaggedLineCorrections,
        final float factorIncreaseForCurvatureMinimum) {

        List<List<CornerRegion>> cornerLists = 
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
                
        for (int i = 0; i < perimeterLists.size(); ++i) {
            List<CornerRegion> list = indexCornerRegionMap.get(Integer.valueOf(i));
            if (list == null) {
                // have to keep ordered, parallel indexes
                cornerLists.add(new ArrayList<CornerRegion>());
            } else {
                
                // detailed check for ordering the perimeters.
                // TODO: consider corrections to this for perimeter also at higher level
                // then it's not needed here
                MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
                PairIntArray cornerXY = new PairIntArray();
                for (int ii = 0; ii < list.size(); ++ii) {
                    CornerRegion cr = list.get(ii);
                    cornerXY.add(cr.getX()[cr.getKMaxIdx()], cr.getY()[cr.getKMaxIdx()]);
                }
                boolean isCW = curveHelper.curveIsOrderedClockwise(cornerXY);
                if (isCW) {
                    int n = cornerXY.getN();
                    if (n > 1) {
                        int end = n >> 1;
                        // 0 1 2 3 4
                        for (int ii = 0; ii < end; ii++) {
                            int idx2 = n - ii - 1;
                            CornerRegion swap = list.get(ii);
                            list.set(ii, list.get(idx2));
                            list.set(idx2, swap);
                        }
                    }
                }
                
                cornerLists.add(list);
            }
        }
        
        if (blobPerimeterHelper.isInDebugMode()) {
            MiscDebug.writeEdgesAndCorners(perimeterLists, cornerLists,
                1, blobPerimeterHelper.getGreyscaleImage(useBinnedImage),
                "blob_corners_" + blobPerimeterHelper.getDebugTag() + "_" 
                + MiscDebug.getCurrentTimeFormatted());
        }

        assert(perimeterLists.size() == cornerLists.size());
        
        return cornerLists;
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
        
        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > sMap =
            cornerMaker.findCornersInScaleSpaceMaps(theEdges, 
                doUseOutdoorMode, outputCorners); 
        
        return cornerMaker.getEdgeCornerRegionMap();        
    }

}
