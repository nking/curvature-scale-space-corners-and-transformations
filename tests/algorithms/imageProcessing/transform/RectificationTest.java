package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Rectification.RectifiedPoints;
import algorithms.matrix.MatrixUtil;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertNotNull;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class RectificationTest extends TestCase {
    
    public RectificationTest() {
    }

    /**
     * Test of epipolar method, of class Rectification.
     */
    public void testEpipolar() throws Exception {
        System.out.println("testEpipolar");
        double[][] k1 = Zhang98Data.getIntrinsicCameraMatrix();
        double[][] k2 = MatrixUtil.copy(k1);
        //x1, x2 size is 3 X 256
        double[][] x1 = Zhang98Data.getObservedFeaturesInImage(1);
        double[][] x2 = Zhang98Data.getObservedFeaturesInImage(5);
                
        RectifiedPoints result = Rectification.epipolar(k1, k2, x1, x2);
        
        assertNotNull(result);
    }

    /**
     * Test of rectify method, of class Rectification.
     */
    @Test
    public void testRectify() throws Exception {
        System.out.println("rectify");
        double[][] x1 = null;
        double[][] x2 = null;
        double oX = 0.0;
        double oY = 0.0;
        Rectification.RectifiedImage expResult = null;
        Rectification.RectifiedImage result = Rectification.rectify(x1, x2, oX, oY);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of hWarp method, of class Rectification.
     */
    @Test
    public void testHWarp() throws Exception {
        System.out.println("hWarp");
        double[][] img = null;
        double[][] h = null;
        Rectification.RectifiedImage expResult = null;
        Rectification.RectifiedImage result = Rectification.hWarp(img, h);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
    
}
