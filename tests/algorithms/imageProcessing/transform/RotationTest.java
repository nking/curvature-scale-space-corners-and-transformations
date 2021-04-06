
package algorithms.imageProcessing.transform;

import algorithms.util.FormatArray;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class RotationTest extends TestCase {
    
    public RotationTest() {
    }

    public void testCreateRodriguesFormulaRotationMatrix() {
        
        //from http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example.html
        
        double[] v = new double[]{-1.451113, -1.827059, -0.179105};
        double[][] r = Rotation.createRodriguesFormulaRotationMatrix(v);
        
        double[][] expected = new double[3][3];
        expected[0] = new double[]{-0.43583, 0.875946, -0.480436};
        expected[1] = new double[]{0.765974, 0.338032, 0.546825};
        expected[2] = new double[]{0.641392, -0.344170, -0.685684};
        
        System.out.printf("r=\n%s\n", FormatArray.toString(r,"%.3e"));
        
        assertEquals(expected.length, r.length);
        double tol = 1e-4;
        double diff;
        for (int i = 0; i < r.length; ++i) {
            for (int j = 0; j < r[i].length; ++j) {
                diff = Math.abs(expected[i][j] - r[i][j]);
                assertTrue(diff < tol);
            }
        }
    }
    
}
