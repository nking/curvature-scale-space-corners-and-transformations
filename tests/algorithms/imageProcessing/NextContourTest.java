package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
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
        
        NextContour nextContour = new NextContour(contours, 
            alreadyVisited);
        
        assertTrue(nextContour.unvisitedContours.size() == contours.size());
                
        //assert that internal data structure is sorted by descending sigma
        float lastSigma = Float.MAX_VALUE;
        
        Iterator<CurvatureScaleSpaceContour> iter = 
            nextContour.unvisitedContours.iterator();
        
        while (iter.hasNext()) {
            
            CurvatureScaleSpaceContour contour = iter.next();
            
            float sigma = contour.getPeakSigma();
            
            assertTrue(sigma <= lastSigma);
            
            lastSigma = sigma;
        }
        
        assertTrue(nextContour.getMatchedContours1().isEmpty());
        assertTrue(nextContour.getMatchedContours2().isEmpty());
        
        //TODO: add an assert that t is increasing for same sigma
        
        // ======= using the find methods has side effect of removing the
        //         contour from the lookup datastructures
        CurvatureScaleSpaceContour contour = 
            nextContour.findTallestContourWithinAScaleSpace();
                
        assertNotNull(contour);
        
        assertTrue(contour.getPeakSigma() == contours.get(0).getPeakSigma());
        assertTrue(contour.getEdgeNumber() == contours.get(0).getEdgeNumber());
                
        // assert internal look-up structures don't have this found contour now        
        assertFalse(nextContour.unvisitedContours.contains(contour));
        
        nextContour.addMatchedContours(contour, new CurvatureScaleSpaceContour(
            12345, 0.12345f));
        
        assertTrue(nextContour.getMatchedContours1().size() == 1);
        assertTrue(nextContour.getMatchedContours2().size() == 1);
        assertTrue(nextContour.getMatchedContours1().get(0).equals(contour));
        assertTrue(nextContour.getMatchedContours2().get(0).getPeakSigma() 
            == 12345);
                
        // ======= using the find methods has side effect of removing the
        //         contour from the lookup datastructures
        contour = 
            nextContour.findTheNextSmallestUnvisitedSibling(contour);
        
        assertNotNull(contour);
        
        assertTrue(contour.getPeakSigma() == contours.get(2).getPeakSigma());
        assertTrue(contour.getEdgeNumber() == contours.get(2).getEdgeNumber());
        
        // assert internal look-up structures don't have this found contour now        
        assertFalse(nextContour.unvisitedContours.contains(contour));
        
        nextContour.addMatchedContours(contour, new CurvatureScaleSpaceContour(
            2345, 0.2345f));
        
        assertTrue(nextContour.getMatchedContours1().size() == 2);
        assertTrue(nextContour.getMatchedContours2().size() == 2);
        assertTrue(nextContour.getMatchedContours1().get(1).equals(contour));
        assertTrue(nextContour.getMatchedContours2().get(1).getPeakSigma() 
            == 2345);
        
        // ====== get the last contour remaining in lookups
        contour = 
            nextContour.findTallestContourWithinAScaleSpace();
        
        assertNotNull(contour);
        
        assertTrue(contour.getPeakSigma() == contours.get(1).getPeakSigma());
        assertTrue(contour.getEdgeNumber() == contours.get(1).getEdgeNumber());
        
        // assert internal look-up structures don't have this found contour now
        assertFalse(nextContour.unvisitedContours.contains(contour));
        
        nextContour.addMatchedContours(contour, new CurvatureScaleSpaceContour(
            345, 0.345f));
        
        assertTrue(nextContour.getMatchedContours1().size() == 3);
        assertTrue(nextContour.getMatchedContours2().size() == 3);
        assertTrue(nextContour.getMatchedContours1().get(2).equals(contour));
        assertTrue(nextContour.getMatchedContours2().get(2).getPeakSigma() 
            == 345);
    }
    
    protected List<CurvatureScaleSpaceContour> getImageContours(String fileName)
        throws Exception {

        List<CurvatureScaleSpaceContour> contours = 
            new ArrayList<CurvatureScaleSpaceContour>();

        String filePath = ResourceFinder.findFileInTestResources(fileName);

        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath);

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
                = instance.convertScaleSpaceMapToSparseImage(scaleSpaceMap, i);
                 
            ContourFinder contourFinder = new ContourFinder();

            List<CurvatureScaleSpaceContour> result = 
                contourFinder.findContours(scaleSpaceImage, i);
        
            contours.addAll(result);
        }
        
        return contours;
    }
}
