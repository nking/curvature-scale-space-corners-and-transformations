package algorithms.imageProcessing;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class CurvatureScaleSpaceContourTest {
    
    public CurvatureScaleSpaceContourTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testContour() {
        
        float sigma = 32;
        float t = 0.5f;
        
        int x = 123;
        int y = 1290;
        
        CurvatureScaleSpaceImagePoint point0 = new CurvatureScaleSpaceImagePoint(
            sigma, t, x, y);
        
        assertTrue(point0.getSigma() == sigma);
        assertTrue(point0.getScaleFreeLength() == t);
        assertTrue(point0.getXCoord() == x);
        assertTrue(point0.getYCoord() == y);
        
        CurvatureScaleSpaceImagePoint point1 = new CurvatureScaleSpaceImagePoint(
            sigma, t + 0.25f, x + 20, y - 100);
        
        CurvatureScaleSpaceImagePoint[] points = new 
            CurvatureScaleSpaceImagePoint[]{point0, point1};
            
        CurvatureScaleSpaceContour contour = new CurvatureScaleSpaceContour
            (sigma, t);
        
        contour.setPeakDetails(points);
        
        assertTrue(contour.getPeakSigma() == sigma);
        assertTrue(contour.getPeakScaleFreeLength() == t);
        assertTrue(contour.getPeakDetails().length == 2);
        
        assertTrue(contour.toString().contains("sigma"));
    }

    public static void main(String[] args) {
        
        try {
            
            CurvatureScaleSpaceContourTest test = 
                new CurvatureScaleSpaceContourTest();
            
            test.testContour();
                        
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
