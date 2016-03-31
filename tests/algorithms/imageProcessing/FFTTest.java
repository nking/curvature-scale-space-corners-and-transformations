package algorithms.imageProcessing;

import algorithms.misc.Complex;
import algorithms.misc.MiscMath;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class FFTTest extends TestCase {

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
    
    public void testFFT_delta() {
        
        FFT fft = new FFT();
        
        int n = 16;
        Complex[] a = new Complex[n];
        
        for (int i = 0; i < n; ++i) {
            a[i] = new Complex(16, 0);
        }
        
        Complex[] aTr = fft.fft(a, true);

        // should be a delta function
        assertTrue(aTr.length == n);
        assertEquals(64.0, aTr[0].re());
        assertEquals(0.0, aTr[0].im());

        for (int i = 1; i < n; ++i) {  
            assertEquals(0.0, aTr[i].re());
            assertEquals(0.0, aTr[i].im());
        }
                       
        fft = new FFT();
        aTr = fft.fft(aTr, true);

        for (int i = 0; i < n; ++i) {  
            assertTrue( Math.abs(a[i].re() - aTr[i].re()) < 0.01);
            assertEquals(0.0, aTr[i].im());
        }
        
    }

    public void testFFT_triangular() throws Exception {

        boolean makePlot = false;
        
        FFT fft = new FFT();
        
        algorithms.util.PolygonAndPointPlotter plotter = null;
        
        if (makePlot) {
            plotter = new algorithms.util.PolygonAndPointPlotter();
        }

        /*
        triangular function
        base is width w and height 2w        
        */
        
        int n = 1 << 2;
        int w = 1;
        Complex[] a = new Complex[n];
        
        a[0] = new Complex(1, 0);
        a[1] = new Complex(2*w, 0);
        a[2] = new Complex(1, 0);
        a[3] = new Complex(1, 0);
        
        if (makePlot) {
            float[] x = new float[a.length];
            float[] y = new float[x.length];
            for (int ii = 0; ii < a.length; ++ii) {
                x[ii] = ii;
                y[ii] = (float)a[ii].re();
            }
            float[] xPolygon = x;
            float[] yPolygon = y;
            float minY = MiscMath.findMin(y);
            float maxY = MiscMath.findMax(y);
            plotter.addPlot(-1, x[y.length - 1] + 1, minY, maxY, x, y, xPolygon,
                yPolygon, "a");
        }
        
        Complex[] aTr = fft.fft(a, true);
        
        if (makePlot) {
            float[] x = new float[a.length];
            float[] y = new float[x.length];
            for (int ii = 0; ii < aTr.length; ++ii) {
                x[ii] = ii;
                y[ii] = (float) aTr[ii].re();
            }
            float[] xPolygon = x;
            float[] yPolygon = y;
            float minY = MiscMath.findMin(y);
            float maxY = MiscMath.findMax(y);
            plotter.addPlot(-1, x[y.length - 1] + 1, minY, maxY, x, y, xPolygon,
                yPolygon, "FFT of a");

            plotter.writeFile();
        }

        assertTrue( Math.abs(aTr[0].re() - 2.5) < 0.001);
        double v = Math.abs(aTr[1].re());
        assertTrue( v < 0.01);
        assertTrue( aTr[2].re() < aTr[1].re());   
        
        // the FFT of triangle is the sync squared function 
        fft = new FFT();
        aTr = fft.fft(aTr, false);
        
        for (int i = 0; i < n; ++i) {  
            assertTrue( Math.abs(a[i].re() - aTr[i].re()) < 0.01);
            assertEquals(0.0, aTr[i].im());
        }
    }
    
    public void testFFT_signum() throws Exception {

        boolean makePlot = false;
        
        FFT fft = new FFT();
        
        algorithms.util.PolygonAndPointPlotter plotter = null;
        
        if (makePlot) {
            plotter = new algorithms.util.PolygonAndPointPlotter();
        }

        /*
        signum function
        */
        int n = 1 << 2;
        Complex[] a = new Complex[n];
        
        a[0] = new Complex(1, 0);
        a[1] = new Complex(1, 0);
        a[2] = new Complex(-1, 0);
        a[3] = new Complex(-1, 0);
        
        if (makePlot) {
            float[] x = new float[a.length];
            float[] y = new float[x.length];
            for (int ii = 0; ii < a.length; ++ii) {
                x[ii] = ii;
                y[ii] = (float)a[ii].re();
            }
            float[] xPolygon = x;
            float[] yPolygon = y;
            float minY = MiscMath.findMin(y);
            float maxY = MiscMath.findMax(y);
            plotter.addPlot(-1, x[y.length - 1] + 1, minY, maxY, x, y, xPolygon,
                yPolygon, "a");
        }
        
        Complex[] aTr = fft.fft(a, true);
        
        if (makePlot) {
            float[] x = new float[a.length];
            float[] y = new float[x.length];
            for (int ii = 0; ii < aTr.length; ++ii) {
                x[ii] = ii;
                y[ii] = (float) aTr[ii].re();
            }
            float[] xPolygon = x;
            float[] yPolygon = y;
            float minY = MiscMath.findMin(y);
            float maxY = MiscMath.findMax(y);
            plotter.addPlot(-1, x[y.length - 1] + 1, minY, maxY, x, y, xPolygon,
                yPolygon, "FFT of a");

            plotter.writeFile();
        }

        assertTrue( Math.abs(aTr[0].re() - 0) < 0.001);
        assertTrue( Math.abs(aTr[1].re() - 1) < 0.001);
        assertTrue( Math.abs(aTr[2].re() - 0) < 0.001);
        assertTrue( Math.abs(aTr[3].re() - 1) < 0.001);
        
        // the FFT of triangle is the sync squared function 
        fft = new FFT();
        aTr = fft.fft(aTr, false);
        
        for (int i = 0; i < n; ++i) {  
            assertTrue( Math.abs(a[i].re() - aTr[i].re()) < 0.01);
            assertTrue( Math.abs(a[i].im() - aTr[i].im()) < 0.01);
        }
    }
    
}
