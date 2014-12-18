package algorithms.imageProcessing.util;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;
import org.ejml.simple.*;
    
/**
 *
 * @author nichole
 */
public class MatrixUtilTest {
    
    public MatrixUtilTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    @Test
    public void testDot() throws Exception {
       
        double[][] m1 = new double[2][3];
        m1[0] = new double[]{0, 1, 0};
        m1[1] = new double[]{1000, 100, 10};
        
        double[][] m2 = new double[3][2];
        m2[0] = new double[]{2, 1};
        m2[1] = new double[]{3, 0};
        m2[2] = new double[]{4, 0};
         
        /*
        0     1     0     2  1
        1000  100  10     3  0
                          4  0
        
        0*2    + 1*3   + 0*0     0*1    +  1*0   +  0*0
        1000*2 + 100*3 + 10*4    1000*1 +  100*0 + 10*0
        */
       
        double[][] m = MatrixUtil.dot(new SimpleMatrix(m1), 
            new SimpleMatrix(m2));
        
        assertTrue(m.length == 2);
        assertTrue(m[0].length == 2);
        assertTrue(m[1].length == 2);
        
        assertTrue(m[0][0] == 3);
        assertTrue(m[1][0] == 2340);
        assertTrue(m[0][1] == 0);
        assertTrue(m[1][1] == 1000);
    }
    
    @Test
    public void testMultiply() throws Exception {

        /*
         * | 1 2 |  | 4 |
         * | 2 3 |  | 3 |
         * | 3 4 |
         *
         *   = | (4 + 6)   |
         *     | (8 + 9)   |
         *     | (12 + 12) |
        */
        double[][] a = new double[2][3];
        a[0] = new double[]{1, 2, 3};
        a[1] = new double[]{2, 3, 4};

        double[] b = new double[]{4, 3};

        double[] m = MatrixUtil.multiply(a, b);

        assertTrue( m[0] ==  10 );
        assertTrue( m[1] ==  17 );
        assertTrue( m[2] ==  24 );
    }

}
