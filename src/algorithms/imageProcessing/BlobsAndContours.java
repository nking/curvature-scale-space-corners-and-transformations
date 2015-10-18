package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.List;

/**
 * holds the blobs, their ordered perimeters, and the scale space image
 * contours of the perimeters
 *
 * @author nichole
 */
public class BlobsAndContours  {
    
    public static List<List<CurvatureScaleSpaceContour>> populateContours(
        final BlobPerimeterHelper blobPerimeterHelper,
        SegmentationType type, boolean useBinnedImage) {

        List<List<CurvatureScaleSpaceContour>> contours 
            = new ArrayList<List<CurvatureScaleSpaceContour>>();

        List<ScaleSpaceCurveImage> scaleSpaceImages
            = new ArrayList<ScaleSpaceCurveImage>();

        boolean setToExtractWeakCurvesTooIfNeeded = false;
        
        if (type.equals(SegmentationType.BINARY)) {
            // might need to restrict this to the binned images
            setToExtractWeakCurvesTooIfNeeded = true;
        }

        boolean allContoursZero = true;

        List<PairIntArray> blobOrderedPerimeters = 
            blobPerimeterHelper.getBlobPerimeters(type, useBinnedImage);
        
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
        
        return contours;
    }

}
