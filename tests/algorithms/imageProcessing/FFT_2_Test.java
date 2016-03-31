package algorithms.imageProcessing;

import algorithms.misc.Complex;
import junit.framework.TestCase;
import thirdparty.ca.uol.aig.fftpack.Complex1D;
import thirdparty.ca.uol.aig.fftpack.ComplexDoubleFFT;

/**
 *
 * @author nichole
 */
public class FFT_2_Test extends TestCase {

    public void testFFT_compare() {
                
        int n = 16;
        Complex[] a = new Complex[n];
        Complex1D x2 = new Complex1D();
        
        for (int k = 0; k < 3; ++k) {
            
            if (k == 0) {
                // delta function
                n = 16;
                a = new Complex[n];
                x2 = new Complex1D();
                x2.x = new double[n];
                x2.y = new double[n];
                for (int i = 0; i < n; ++i) {
                    a[i] = new Complex(16, 0);
                    x2.x[i] = 16;
                    x2.y[i] = 0;
                }
            } else if (k == 1) {
                // triangle function
                n = 1 << 2;
                int w = 1;
                a = new Complex[n];
                a[0] = new Complex(1, 0);
                a[1] = new Complex(2 * w, 0);
                a[2] = new Complex(1, 0);
                a[3] = new Complex(1, 0);
                x2 = new Complex1D();
                x2.x = new double[n];
                x2.y = new double[n];
                x2.x[0] = 1;
                x2.y[0] = 0;
                x2.x[1] = 2 * w;
                x2.y[1] = 0;
                x2.x[2] = 1;
                x2.y[2] = 0;
                x2.x[3] = 1;
                x2.y[3] = 0;
            } else if (k == 2) {
                //signum function
                n = 1 << 2;
                a = new Complex[n];
                a[0] = new Complex(1, 0);
                a[1] = new Complex(1, 0);
                a[2] = new Complex(-1, 0);
                a[3] = new Complex(-1, 0);
                x2 = new Complex1D();
                x2.x = new double[n];
                x2.y = new double[n];
                x2.x[0] = 1;
                x2.y[0] = 0;
                x2.x[1] = 1;
                x2.y[1] = 0;
                x2.x[2] = -1;
                x2.y[2] = 0;
                x2.x[3] = -1;
                x2.y[3] = 0;
            }
        
            FFT fft = new FFT();
            fft.setToNotNormalize();
            ComplexDoubleFFT fft2 = new ComplexDoubleFFT(n);
                
            Complex[] aTr = fft.fft(a, true);
            fft2.ft(x2);
            for (int i = 0; i < n; ++i) {
                double r1 = aTr[i].re();
                double r2 = x2.x[i];

                double im1 = aTr[i].im();
                double im2 = x2.y[i];

                assertTrue(Math.abs(r1 - r2) < 0.01);
                assertTrue(Math.abs(im1 - im2) < 0.01);
            }

            double norm = 1./(double)n;

            // --- apply inverse FFT -----
            aTr = fft.fft(aTr, false);
            fft2.bt(x2);
            for (int i = 0; i < n; ++i) {

                double r0 = a[i].re();

                double r1 = aTr[i].re();
                double r2 = x2.x[i];

                double im1 = aTr[i].im();
                double im2 = x2.y[i];

                assertTrue(Math.abs(r1 - r2) < 0.01);
                assertTrue(Math.abs(im1 - im2) < 0.01);
                assertTrue(Math.abs(r0 - (r1*norm)) < 0.01);
            }
        }
    }
    
}
