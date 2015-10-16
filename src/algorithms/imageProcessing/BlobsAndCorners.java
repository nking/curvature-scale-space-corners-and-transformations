package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
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

    protected final BlobsAndPerimeters blobsAndPerimeters;
    
    private final List<List<CornerRegion>> cornerLists = 
        new ArrayList<List<CornerRegion>>();

    protected boolean enableJaggedLineCorrections = true;
    
    protected float factorIncreaseForCurvatureMinimum = 1.f;
            
    private boolean outdoorMode = false;
    
    /**
     * constructor, does all the work of extracting blobs, perimeter, and
     * scale space image corners.
     *
     * @param imgGrey
     * @param imgSeg
     * @param smallestGroupLimit
     * @param largestGroupLimit
     * @param type
     * @param segmentedToLineDrawing true if image contains mostly white
     * pixels and black lines for object contours (this is the result of
     * segmentation using adaptive mean thresholding, for example)
     */
    public BlobsAndCorners(GreyscaleImage imgGrey, GreyscaleImage imgSeg,
        int smallestGroupLimit, int largestGroupLimit, 
        SegmentationType type, boolean segmentedToLineDrawing) {
        
        blobsAndPerimeters = new BlobsAndPerimeters(imgGrey, imgSeg, 
            smallestGroupLimit, largestGroupLimit, type, 
            segmentedToLineDrawing);
        
        populateCorners();
    }
    
    /**
     * constructor for using in debug mode, does all the work of extracting
     * blobs, perimeter, and scale space image corners.
     *
     * @param imgGrey
     * @param imgSeg
     * @param smallestGroupLimit
     * @param largestGroupLimit
     * segmentedToLineDrawing true if image contains mostly white
     * pixels and black lines for object contours (this is the result of
     * segmentation using adaptive mean thresholding, for example)
     * @param type
     * @param segmentedToLineDrawing
     * @param debugTag
     */
    public BlobsAndCorners(GreyscaleImage imgGrey, GreyscaleImage imgSeg,
        int smallestGroupLimit, int largestGroupLimit, SegmentationType type,
        boolean segmentedToLineDrawing, String debugTag) {

        blobsAndPerimeters = new BlobsAndPerimeters(imgGrey, imgSeg, 
            smallestGroupLimit, largestGroupLimit, type, 
            segmentedToLineDrawing, debugTag);
        
        populateCorners();
    }

    protected void populateCorners() {

        cornerLists.clear();

        List<PairIntArray> perimeterLists = blobsAndPerimeters.getBlobOrderedPerimeters();
        
        PairIntArray allCorners = new PairIntArray();
        
        Map<Integer, List<CornerRegion> > indexCornerRegionMap = 
            findCornersInScaleSpaceMaps(perimeterLists, outdoorMode, allCorners);
        
        List<Integer> remove = new ArrayList<Integer>();
        
        for (int i = 0; i < perimeterLists.size(); ++i) {
            List<CornerRegion> list = indexCornerRegionMap.get(Integer.valueOf(i));
            if (list == null) {
                remove.add(i);
            } else {
                cornerLists.add(list);
            }
        }
        
        for (int i = (remove.size() - 1); i > -1; --i) {
            perimeterLists.remove(remove.get(i).intValue());
        }

        assert(perimeterLists.size() == cornerLists.size());
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
     * @return scale space maps for each edge
     */
    protected Map<Integer, List<CornerRegion> >
    findCornersInScaleSpaceMaps(final List<PairIntArray> theEdges, 
        final boolean doUseOutdoorMode, final PairIntArray outputCorners) {

        CSSCornerMaker cornerMaker = new CSSCornerMaker();
        
        if (enableJaggedLineCorrections) {
            cornerMaker.enableJaggedLineCorrections();
        }
       
        cornerMaker.increaseFactorForCurvatureMinimum(
            factorIncreaseForCurvatureMinimum);
        
        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > sMap =
            cornerMaker.findCornersInScaleSpaceMaps(theEdges, 
                doUseOutdoorMode, outputCorners); 
        
        return cornerMaker.getEdgeCornerRegionMap();        
    }

    /**
     * @return the corners
     */
    public List<List<CornerRegion>> getCorners() {
        return cornerLists;
    }

    public List<Set<PairInt>> getBlobs() {
        return blobsAndPerimeters.getBlobs();
    }

    public List<PairIntArray> getBlobOrderedPerimeters() {
        return blobsAndPerimeters.getBlobOrderedPerimeters();
    }
    
    public void enableJaggedLineCorrections() {
        enableJaggedLineCorrections = true;
    }
    public void disableJaggedLineCorrections() {
        enableJaggedLineCorrections = false;
    }
    
    public void increaseFactorForCurvatureMinimum(float factor) {
        factorIncreaseForCurvatureMinimum = factor;
    }
    
    public void resetFactorForCurvatureMinimum() {
        factorIncreaseForCurvatureMinimum = 1.f;
    }
}
