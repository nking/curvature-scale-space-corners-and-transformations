package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * holds the blobs, their ordered perimeters, and the scale space image
 * contours of the perimeters
 *
 * @author nichole
 */
public class BlobsAndContours  {

    protected final BlobsAndPerimeters blobsAndPerimeters;
    
    private final List<List<CurvatureScaleSpaceContour>> contours = 
        new ArrayList<List<CurvatureScaleSpaceContour>>();

    /**
     * constructor, does all the work of extracting blobs, perimeter, and
     * scale space image contours.
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
    public BlobsAndContours(GreyscaleImage imgGrey, GreyscaleImage imgSeg,
        int smallestGroupLimit, int largestGroupLimit, 
        SegmentationType type, boolean segmentedToLineDrawing) {
        
        blobsAndPerimeters = new BlobsAndPerimeters(imgGrey, imgSeg, 
            smallestGroupLimit, largestGroupLimit, type, 
            segmentedToLineDrawing);
        
        populateContours();
    }
    
    /**
     * constructor for using in debug mode, does all the work of extracting
     * blobs, perimeter, and scale space image contours.
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
    public BlobsAndContours(GreyscaleImage imgGrey, GreyscaleImage imgSeg,
        int smallestGroupLimit, int largestGroupLimit, 
        SegmentationType type,
        boolean segmentedToLineDrawing, String debugTag) {

        blobsAndPerimeters = new BlobsAndPerimeters(imgGrey, imgSeg, 
            smallestGroupLimit, largestGroupLimit, type, 
            segmentedToLineDrawing, debugTag);
        
        populateContours();
    }

    protected void populateContours() {

        contours.clear();

        List<ScaleSpaceCurveImage> scaleSpaceImages
            = new ArrayList<ScaleSpaceCurveImage>();

        boolean setToExtractWeakCurvesTooIfNeeded = false;
        
        if (blobsAndPerimeters.getSegmentationType().equals(
            SegmentationType.BINARY)) {
            // might need to restrict this to the binned images
            setToExtractWeakCurvesTooIfNeeded = true;
        }

        boolean allContoursZero = true;

        List<PairIntArray> blobOrderedPerimeters = 
            blobsAndPerimeters.getBlobOrderedPerimeters();
        
        for (int edgeIndex = 0; edgeIndex < blobOrderedPerimeters.size(); ++edgeIndex) {

            PairIntArray edge = blobOrderedPerimeters.get(edgeIndex);
    
MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
double[] xycen = curveHelper.calculateXYCentroids(edge);

            ScaleSpaceCurveImage sscImg =
                CurvatureScaleSpaceInflectionSingleEdgeMapper.createScaleSpaceImage(
                    edge, edgeIndex);

            scaleSpaceImages.add(sscImg);

            List<CurvatureScaleSpaceContour> c =
                CurvatureScaleSpaceInflectionSingleEdgeMapper.populateContours(
                sscImg, edgeIndex, setToExtractWeakCurvesTooIfNeeded, edge);
            
            /*if (c.isEmpty() && (somewhat large)) {
                c = CurvatureScaleSpaceInflectionSingleEdgeMapper.populateContours(
                    sscImg, edgeIndex, true, edge);
            }*/

            //TODO: consider extracting smaller contours for each edge if needed
            // instead of only when there are no strong contours
            
            contours.add(c);

            if (!c.isEmpty()) {
                allContoursZero = false;
            }
        }

        //TODO: consider whether want to include weak contours
        if (allContoursZero) {

            setToExtractWeakCurvesTooIfNeeded = true;

            contours.clear();

            for (int edgeIndex = 0; edgeIndex < blobOrderedPerimeters.size(); ++edgeIndex) {

                PairIntArray edge = blobOrderedPerimeters.get(edgeIndex);
                
                ScaleSpaceCurveImage sscImg = scaleSpaceImages.get(edgeIndex);

                List<CurvatureScaleSpaceContour> c =
                    CurvatureScaleSpaceInflectionSingleEdgeMapper.populateContours(
                    sscImg, edgeIndex, setToExtractWeakCurvesTooIfNeeded, edge);

                contours.add(c);
            }
        }

        assert(blobOrderedPerimeters.size() == contours.size());
    }

    /**
     * @return the contours
     */
    public List<List<CurvatureScaleSpaceContour>> getContours() {
        return contours;
    }

    public List<Set<PairInt>> getBlobs() {
        return blobsAndPerimeters.getBlobs();
    }

    public List<PairIntArray> getBlobOrderedPerimeters() {
        return blobsAndPerimeters.getBlobOrderedPerimeters();
    }
}
