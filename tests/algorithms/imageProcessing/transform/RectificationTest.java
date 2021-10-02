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
    
}
