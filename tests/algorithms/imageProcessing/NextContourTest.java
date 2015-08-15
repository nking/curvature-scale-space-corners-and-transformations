package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class NextContourTest {
    
    public NextContourTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testInternalDataStructures0() throws Exception {
        
        List<CurvatureScaleSpaceContour> contours = 
            getImageContours("closed_curve.png");
        
        List<CurvatureScaleSpaceContour> alreadyVisited = 
            new ArrayList<CurvatureScaleSpaceContour>();
        
        NextContour nextContour = new NextContour(contours, alreadyVisited);
        
        assertTrue(nextContour.origContours.size() == contours.size());
        assertTrue(nextContour.remainingContours.size() 
            == contours.size());
        
        //assert that internal data structure is sorted by descending sigma
        float lastSigma = Float.MAX_VALUE;
        Iterator<CurvatureScaleSpaceContour> iter = nextContour.origContours.iterator();
        while (iter.hasNext()) {
            
            CurvatureScaleSpaceContour contour = iter.next();
            
            float sigma = contour.getPeakSigma();
            
            assertTrue(sigma <= lastSigma);
            
            lastSigma = sigma;
        }
        
        //assert that internal data structure is sorted by descending sigma
        lastSigma = Float.MAX_VALUE;
        for (CurvatureScaleSpaceContour contour : nextContour.remainingContours) {
                        
            float sigma = contour.getPeakSigma();
            
            assertTrue(sigma <= lastSigma);
            
            lastSigma = sigma;
        }
      
        
        // ======= using the find methods has side effect of removing the
        //         contour from the lookup datastructures
        CurvatureScaleSpaceContour contour = 
            nextContour.findTallestContourWithinScaleSpace();
                
        assertNotNull(contour);
        
        assertTrue(contour.getPeakSigma() == contours.get(0).getPeakSigma());
        assertTrue(contour.getEdgeNumber() == contours.get(0).getEdgeNumber());
                        
        // assert internal look-up structures don't have this found contour now        
        assertTrue(nextContour.remainingContours.size() == (contours.size() - 1));
        assertFalse(nextContour.remainingContours.contains(contours.get(0)));
        
        contour = 
            nextContour.findTheNextSmallestUnvisitedSibling(contours.get(1));
        
        assertNotNull(contour);
        
        assertTrue(contour.getPeakSigma() == contours.get(2).getPeakSigma());
        assertTrue(contour.getEdgeNumber() == contours.get(2).getEdgeNumber());
        
        contour = nextContour.findTheNextSmallestUnvisitedSibling(contour);
        assertNull(contour);
        
        // assert internal look-up structures don't have this found contour now        
        assertTrue(nextContour.remainingContours.size() == (contours.size() - 2));
        assertFalse(nextContour.remainingContours.contains(contours.get(2)));
        
        // ====== get the last contour remaining in lookups
        contour = nextContour.findTallestContourWithinScaleSpace();
        
        assertNotNull(contour);
        
        assertTrue(contour.getPeakSigma() == contours.get(1).getPeakSigma());
        assertTrue(contour.getEdgeNumber() == contours.get(1).getEdgeNumber());
        
        // assert internal look-up structures don't have this found contour now
       assertTrue(nextContour.remainingContours.isEmpty());
        
    }
    
    protected List<CurvatureScaleSpaceContour> getImageContours(String fileName)
        throws Exception {

        List<CurvatureScaleSpaceContour> contours = 
            new ArrayList<CurvatureScaleSpaceContour>();

        String filePath = ResourceFinder.findFileInTestResources(fileName);

        ImageExt img = ImageIOHelper.readImageExt(filePath);

        CurvatureScaleSpaceImageMaker instance
            = new CurvatureScaleSpaceImageMaker(img);

        instance.initialize();
        
        instance.useLineDrawingMode();
        
        List<PairIntArray> curves = instance.getClosedCurves();

        for (int i = 0; i < curves.size(); i++) {

            PairIntArray curve = curves.get(i);

            Map<Float, ScaleSpaceCurve> scaleSpaceMap
                = instance.createScaleSpaceMetricsForEdge2(curve);

            ScaleSpaceCurveImage scaleSpaceImage
                = instance.convertScaleSpaceMapToSparseImage(scaleSpaceMap, i,
                  curve.getN());
                 
            ContourFinder contourFinder = new ContourFinder();

            List<CurvatureScaleSpaceContour> result = 
                contourFinder.findContours(scaleSpaceImage, i);
        
            contours.addAll(result);
        }
        
        return contours;
    }
}
