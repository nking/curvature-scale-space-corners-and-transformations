package algorithms.imageProcessing.util;

import java.util.Arrays;
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

        double[][] a = new double[2][3];
        a[0] = new double[]{1, 2, 3};
        a[1] = new double[]{2, 3, 4};

        double[] b = new double[]{4, 3};

        double[] m = MatrixUtil.multiply(a, b);

        assertTrue( m[0] ==  10 );
        assertTrue( m[1] ==  17 );
        assertTrue( m[2] ==  24 );
     
        double[][] c = new double[3][];
        for (int i = 0; i < 3; i++) {
            c[i] = new double[3];
            for (int j = 0; j < 3; j++) {
                c[i][j] = 1;
            }
        }
        c[1][2] = 0;
        
           
        /*
        a b c      p0 p1 p2
        d e f      p3 p4 p5
                   p6 p7 p8        
        a*p0+... a*p a*p
        d*p0+... d*p d*p
        
        1 2 3    1 1 1
        2 3 4    1 1 0
                 1 1 1
        
        (1*1 + 2*1 + 3*1)  (1*1 + 2*1 + 3*1)  (1*1 + 0 + 3*1)
        (2*1 + 3*1 + 4*1)  (2*1 + 3*1 + 4*1)  (2*1 + 0 + 4*1)
        
        6  6 4
        9  9 6
        */
        
        double[][] d = MatrixUtil.multiply(a, c);

        assertTrue(d[0][0] == 6);
        assertTrue(d[0][1] == 6);
        assertTrue(d[0][2] == 4);
        assertTrue(d[1][0] == 9);
        assertTrue(d[1][1] == 9);
        assertTrue(d[1][2] == 6);
        
        /*
        example:  m is 1 2 3
                       4 5 6
                       7 8 9
        
                  n is 100  101  1
                       200  201  1
        
        multiply m by transpose of n:
        
        1 2 3     100  200
        4 5 6     101  201
        7 8 9      1   1
        
        (1*100 + 2*101 + 3*1)   (1*200 + 2*201 + 3*1)     305   605
        (4*100 + 5*101 + 6*1)   (4*200 + 5*201 + 6*1)  =  911  1811
        (7*100 + 8*101 + 9*1)   (7*200 + 8*201 + 9*1)    1517  3017
        */
        double[][] aa = new double[3][];
        aa[0] = new double[]{1, 2, 3};
        aa[1] = new double[]{4, 5, 6};
        aa[2] = new double[]{7, 8, 9};
        
        double[][] bb = new double[2][];
        bb[0] = new double[]{100, 101, 1};
        bb[1] = new double[]{200, 201, 1};
        
        double[][] cc = MatrixUtil.multiplyByTranspose(aa, bb);
        
        assertTrue(cc[0][0] == 305);
        assertTrue(cc[0][1] == 605);
        assertTrue(cc[1][0] == 911);
        assertTrue(cc[1][1] == 1811);
        assertTrue(cc[2][0] == 1517);
        assertTrue(cc[2][1] == 3017);
    }

    @Test
    public void testAdd() throws Exception {

        double[] a = new double[]{1, 2, 3, 4};
        double[] b = new double[]{100, 100, 100, 100};

        double[] expected = new double[]{101, 102, 103, 104};
        
        double[] c = MatrixUtil.add(a, b);
        
        assertTrue(Arrays.equals(expected, c));
    }
    
    @Test
    public void testAdd2() throws Exception {

        float[] a = new float[]{1, 2, 3, 4};
        float[] b = new float[]{100, 100, 100, 100};

        float[] expected = new float[]{101, 102, 103, 104};
        
        float[] c = MatrixUtil.add(a, b);
        
        assertTrue(Arrays.equals(expected, c));
    }
    
    @Test
    public void testSubtract() throws Exception {

        float[] a = new float[]{100, 100, 100, 100};
        float[] b = new float[]{1, 2, 3, 4};

        float[] expected = new float[]{99, 98, 97, 96};
        
        float[] c = MatrixUtil.subtract(a, b);
        
        assertTrue(Arrays.equals(expected, c));
    }
}
