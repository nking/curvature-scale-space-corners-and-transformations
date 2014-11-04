package algorithms.imageProcessing;

import algorithms.ResourceFinder;
import java.util.ArrayList;
import java.util.HashMap;
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
        
        Map<Integer, List<Integer> > curveIndexToContour =  
            new HashMap<Integer, List<Integer> >();
        
        for (int i = 0; i < contours.size(); i++) {
            
            CurvatureScaleSpaceContour contour = contours.get(i);
            
            Integer curveIdx = Integer.valueOf(contour.getEdgeNumber());
            
            List<Integer> indexes = curveIndexToContour.get(curveIdx);
            if (indexes == null) {
                indexes = new ArrayList<Integer>();
                curveIndexToContour.put(curveIdx, indexes);
            }
            
            indexes.add(Integer.valueOf(i));            
        }
        
        List<CurvatureScaleSpaceContour> alreadyVisited = 
            new ArrayList<CurvatureScaleSpaceContour>();
        
        NextContour nextContour = new NextContour(contours, true, 
            curveIndexToContour, alreadyVisited);
        
        assertTrue(nextContour.origContours.size() == contours.size());
        assertTrue(nextContour.curveIndexToOrigContours.size() 
            == curveIndexToContour.size());
        
        assertTrue(nextContour.contourIndex.size() == contours.size());
        
        assertTrue(nextContour.curveList.length == curveIndexToContour.size());
        
        //assert that internal data structure is sorted by descending sigma
        float lastSigma = Float.MAX_VALUE;
        for (int i = 0; i < nextContour.origContours.size(); i++) {
            
            CurvatureScaleSpaceContour contour = 
                nextContour.origContours.get(i);
            
            float sigma = contour.getPeakSigma();
            
            assertTrue(sigma <= lastSigma);
            
            lastSigma = sigma;
        }
        
        
        //assert that internal data structure is sorted by descending sigma
        lastSigma = Float.MAX_VALUE;
        for (PairInt ci : nextContour.contourIndex) {
                        
            int curveIdx = ci.getX();
            
            int contourIdx = ci.getY();
            
            CurvatureScaleSpaceContour contour = 
                nextContour.origContours.get(contourIdx);
            
            assertTrue(contour.getEdgeNumber() == curveIdx);
            
            float sigma = contour.getPeakSigma();
            
            assertTrue(sigma <= lastSigma);
            
            assertTrue(nextContour.curveList[curveIdx].length > 0);
            boolean found = false;
            for (int j = 0; j < nextContour.curveList[curveIdx].length; j++) {
                PairInt ci2 = nextContour.curveList[curveIdx][j];
                if ((ci2.getX() == ci.getX()) && (ci2.getY() == ci.getY())) {
                    found = true;
                    break;
                }
            }
            assertTrue(found);
            
            lastSigma = sigma;
        }
        
        int nTot = 0;
        for (int i = 0; i < nextContour.curveList.length; i++) {
            
            PairInt[] indexes = nextContour.curveList[i];
            
            for (int j = 0; j < indexes.length; j++) {
                PairInt ci = indexes[j];
                assertTrue(nextContour.contourIndex.contains(ci));
                int idx = ci.getY();
                assertTrue(idx > -1 && idx < nextContour.origContours.size());
            }
            
            nTot += indexes.length;
        }
        
        assertTrue(nTot == contours.size());
        
        
        // ======= using the find methods has side effect of removing the
        //         contour from the lookup datastructures
        CurvatureScaleSpaceContour contour = 
            nextContour.findTallestContourWithinAScaleSpace(0);
                
        assertNotNull(contour);
        
        assertTrue(contour.getPeakSigma() == contours.get(0).getPeakSigma());
        assertTrue(contour.getEdgeNumber() == contours.get(0).getEdgeNumber());
                
        // assert internal look-up structures don't have this found contour now        
        assertTrue(nextContour.contourIndex.size() == (contours.size() - 1));
        assertTrue(nextContour.curveList.length == curveIndexToContour.size());
        assertTrue(nextContour.curveList[0].length == 3);
        assertNull(nextContour.curveList[0][0]);
        
        
        // ======= using the find methods has side effect of removing the
        //         contour from the lookup datastructures
        PairInt target = nextContour.curveList[0][1];
        contour = 
            nextContour.findTheNextSmallestUnvisitedSibling(target);
        
        assertNotNull(contour);
        
        assertTrue(contour.getPeakSigma() == contours.get(2).getPeakSigma());
        assertTrue(contour.getEdgeNumber() == contours.get(2).getEdgeNumber());
        
        // assert internal look-up structures don't have this found contour now        
        assertTrue(nextContour.contourIndex.size() == 1);
        assertTrue(nextContour.curveList.length == curveIndexToContour.size());
        assertTrue(nextContour.curveList[0].length == 3);
        assertNull(nextContour.curveList[0][0]);
        assertNull(nextContour.curveList[0][2]);
        
        // ====== get the last contour remaining in lookups
        contour = 
            nextContour.findTallestContourWithinAScaleSpace(0);
        
        assertNotNull(contour);
        
        assertTrue(contour.getPeakSigma() == contours.get(1).getPeakSigma());
        assertTrue(contour.getEdgeNumber() == contours.get(1).getEdgeNumber());
        
        // assert internal look-up structures don't have this found contour now
        assertTrue(nextContour.contourIndex.isEmpty());
        assertTrue(nextContour.curveList.length == curveIndexToContour.size());
        assertNull(nextContour.curveList[0]);
        
    }
    
    @Test
    public void testInternalDataStructures1() throws Exception {
        
        List<CurvatureScaleSpaceContour> contours = 
            getImageContours("closed_curve.png");
        
        Map<Integer, List<Integer> > curveIndexToContour =  
            new HashMap<Integer, List<Integer> >();
        
        for (int i = 0; i < contours.size(); i++) {
            
            CurvatureScaleSpaceContour contour = contours.get(i);
            
            Integer curveIdx = Integer.valueOf(contour.getEdgeNumber());
            
            List<Integer> indexes = curveIndexToContour.get(curveIdx);
            if (indexes == null) {
                indexes = new ArrayList<Integer>();
                curveIndexToContour.put(curveIdx, indexes);
            }
            
            indexes.add(Integer.valueOf(i));            
        }
        
        List<CurvatureScaleSpaceContour> alreadyVisited = 
            new ArrayList<CurvatureScaleSpaceContour>();
        alreadyVisited.add(contours.get(0));
        
        NextContour nextContour = new NextContour(contours, true,
            curveIndexToContour, alreadyVisited);
        
        assertTrue(nextContour.origContours.size() == contours.size());
        assertTrue(nextContour.curveIndexToOrigContours.size() 
            == curveIndexToContour.size());
        
        assertTrue(nextContour.contourIndex.size() == (contours.size() - 1));
        
        assertTrue(nextContour.curveList.length == curveIndexToContour.size());
        
        //assert that internal data structure is sorted by descending sigma
        float lastSigma = Float.MAX_VALUE;
        for (int i = 0; i < nextContour.origContours.size(); i++) {
            
            CurvatureScaleSpaceContour contour = 
                nextContour.origContours.get(i);
            
            float sigma = contour.getPeakSigma();
            
            assertTrue(sigma <= lastSigma);
            
            lastSigma = sigma;
        }
        
        
        //assert that internal data structure is sorted by descending sigma
        lastSigma = Float.MAX_VALUE;
        for (PairInt ci : nextContour.contourIndex) {
                        
            int curveIdx = ci.getX();
            
            int contourIdx = ci.getY();
            
            CurvatureScaleSpaceContour contour = 
                nextContour.origContours.get(contourIdx);
            
            assertTrue(contour.getEdgeNumber() == curveIdx);
            
            float sigma = contour.getPeakSigma();
            
            assertTrue(sigma <= lastSigma);
            
            assertTrue(nextContour.curveList[curveIdx].length > 0);
            boolean found = false;
            for (int j = 0; j < nextContour.curveList[curveIdx].length; j++) {
                PairInt ci2 = nextContour.curveList[curveIdx][j];
                if ((ci2.getX() == ci.getX()) && (ci2.getY() == ci.getY())) {
                    found = true;
                    break;
                }
            }
            assertTrue(found);
            
            lastSigma = sigma;
        }
        
        int nTot = 0;
        for (int i = 0; i < nextContour.curveList.length; i++) {
            
            PairInt[] indexes = nextContour.curveList[i];
            
            for (int j = 0; j < indexes.length; j++) {
                PairInt ci = indexes[j];
                assertTrue(nextContour.contourIndex.contains(ci));
                int idx = ci.getY();
                assertTrue(idx > -1 && idx < nextContour.origContours.size());
            }
            
            nTot += indexes.length;
        }
        
        assertTrue(nTot == (contours.size() - 1));
        
    }

    protected List<CurvatureScaleSpaceContour> getImageContours(String fileName)
        throws Exception {

        List<CurvatureScaleSpaceContour> contours = 
            new ArrayList<CurvatureScaleSpaceContour>();

        String filePath = ResourceFinder.findFileInTestResources(fileName);

        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath);

        CurvatureScaleSpaceImageMaker instance
            = new CurvatureScaleSpaceImageMaker(img, true);

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
