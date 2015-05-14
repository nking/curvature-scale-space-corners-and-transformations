package algorithms.imageProcessing;

import algorithms.misc.Complex;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class FFTTest extends TestCase {

    public FFTTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    public void testBitReverseCopy() throws Exception {
        
        // n = 7
        double[] a = new double[]{0, 1, 2, 3, 4, 5, 6};
        
        //int nBits = MiscMath.numberOfBits(a.length - 1);
        
        // 0 000  -->  000  0
        // 1 001       100  4
        // 2 010       010  2
        // 3 011       110  6
        // 4 100       001  1
        // 5 101       101  5
        // 6 110       011  3
        
        FFT fft = new FFT();
        
        double[] r = fft.bitReverseCopy(a);
        assertTrue(a[0] == r[0]);
        assertTrue(a[1] == r[4]);
        assertTrue(a[2] == r[2]);
        assertTrue(a[3] == r[6]);
        assertTrue(a[4] == r[1]);
        assertTrue(a[5] == r[5]);
        assertTrue(a[6] == r[3]);
        
        Complex[] a1 = new Complex[a.length];
        a1[0] = new Complex(0, 0);
        a1[1] = new Complex(1, 0);
        a1[2] = new Complex(2, 0);
        a1[3] = new Complex(3, 0);
        a1[4] = new Complex(4, 0);
        a1[5] = new Complex(5, 0);
        a1[6] = new Complex(6, 0);
        Complex[] r1 = fft.bitReverseCopy(a1);
        assertTrue(a1[0].re() == r1[0].re());
        assertTrue(a1[1].re() == r1[4].re());
        assertTrue(a1[2].re() == r1[2].re());
        assertTrue(a1[3].re() == r1[6].re());
        assertTrue(a1[4].re() == r1[1].re());
        assertTrue(a1[5].re() == r1[5].re());
        assertTrue(a1[6].re() == r1[3].re());
    }
    
    public void testFFT() {

        double[] x = new double[]{
            -0.03480425839330703, 0.07910192950176387,
             0.7233322451735928, 0.1659819820667019};

        FFT fft = new FFT();
        
        double[] expectedReal = new double[]{
             0.466806, -0.379068, 0.221722, -0.379068};

        Complex[] x4 = new Complex[4];
        x4[0] = new Complex(x[0], 0);
        x4[1] = new Complex(x[1], 0);
        x4[2] = new Complex(x[2], 0);
        x4[3] = new Complex(x[3], 0);
        Complex[] y4 = fft.fft(x4);
        for (int i = 0; i < y4.length; i++) {
            double v1Diff = y4[i].re() - expectedReal[i];
            assertTrue( Math.abs( v1Diff ) < 0.01 );
        }
        Complex[] y4Inv = fft.fft(y4, false);
        for (int i = 0; i < y4.length; i++) {
            double v1Diff = y4Inv[i].re() - x[i];
            assertTrue( Math.abs( v1Diff ) < 0.01 );
        }
        
        
        //-----------
        x = new double[]{0, 1, 1, 0};

        int[] x2 = new int[]{0, 1, 1, 0};
        double[] y = fft.fft(x2);

        assertTrue(y.length == 4);
        
        
        //-----------
        x = new double[]{1, 1, 1, 1};
        x2 = new int[]{1, 1, 1, 1};
        
        y = fft.fft(x2);

        assertTrue(y.length == 4);
        
        //-----------
        x = new double[]{0, 0, 0, 1, 1, 0, 0, 0};
        x2 = new int[]{0, 0, 0, 1, 1, 0, 0, 0};
        y = fft.fft(x2);
        
        x = new double[]{1, 0, 0, 1, 1, 0, 0, 1};
        x2 = new int[]{1, 0, 0, 1, 1, 0, 0, 1};
        y = fft.fft(x2);
        
    }
}
