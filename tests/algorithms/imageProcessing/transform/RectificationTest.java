package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Rectification.RectifiedPoints;
import algorithms.matrix.MatrixUtil;
import java.security.NoSuchAlgorithmException;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertNotNull;
import no.uib.cipr.matrix.NotConvergedException;
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
                
        RectifiedPoints result;
        
        result = Rectification.epipolar(k1, k2, x1, x2);
        assertNotNull(result);
        
        double oX = k1[0][2]; 
        double oY = k1[1][2];
        result = Rectification.rectify(x1, x2, oX, oY);
        
        int i;
        double[][] x1R= result.getX1();
        double[][] x2R= result.getX2();
        int n = x1R[0].length;
        System.out.println("rectified");
        for (i = 0; i < n; ++i) {
            System.out.printf("%d) (%.1f, %.1f, %.1f)  (%.1f, %.1f, %.1f)\n",
                i, x1R[0][i], x1R[1][i], x1R[2][i],
                x2R[0][i], x2R[1][i], x2R[2][i]);
        }
      
    }
    
}
