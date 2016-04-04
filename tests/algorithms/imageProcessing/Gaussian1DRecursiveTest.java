package algorithms.imageProcessing;

import algorithms.util.PolygonAndPointPlotter;
import algorithms.misc.MiscMath;
import java.io.IOException;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class Gaussian1DRecursiveTest extends TestCase {
    
    public Gaussian1DRecursiveTest() {
    }

    public void testCompareRecursiveConvolveGaussian1D() throws IOException {
        
        boolean doPlot = false;
                              
        float[] gSIGMATWO = Gaussian1D.getKernel(SIGMA.TWO);
        float[] gSIGMAFOUR = Gaussian1D.getKernel(SIGMA.FOUR);
        
        int n1 = 2 * gSIGMAFOUR.length - 1;
        int n1Alt = 2 * gSIGMATWO.length - 1;
        if (n1Alt > n1) {
            n1 = n1Alt;
        }
        int ns1 = (n1 - gSIGMAFOUR.length)/2;
        
        int n0 = n1;
        int ns0 = (n1 - gSIGMATWO.length)/2;
        
        float[] x0 = new float[n0];
        float[] y0 = new float[n0];
        
        float[] x1 = new float[n1];
        float[] y1 = new float[n1];
        
        float norm1 = 1.f/MiscMath.findMax(gSIGMAFOUR);
        
        int idx = 0;
        for (int i = 0; i < n0; ++i) {
            x0[i] = i;
            if (i < ns0 || (idx > (gSIGMATWO.length - 1))) {
                y0[i] = 0;
            } else {
                y0[i] = gSIGMATWO[idx];
                idx++;
            }
        }
        
        idx = 0;
        for (int i = 0; i < n1; ++i) {
            x1[i] = i;
            if (i < ns1 || (idx > (gSIGMAFOUR.length - 1))) {
                y1[i] = 0;
            } else {
                y1[i] = norm1 * gSIGMAFOUR[idx];
                idx++;
            }
        }
        
        float fwhm0 = GaussianHelperForTests.measureFWHM(x0, y0);
               
        PolygonAndPointPlotter plotter = null;
        
        if (doPlot) {
            
            plotter = new PolygonAndPointPlotter();

            float xMax, yMax, xMin, yMin;

            xMin = MiscMath.findMin(x0);
            xMax = MiscMath.findMax(x0);
            yMin = MiscMath.findMin(y0);
            yMax = MiscMath.findMax(y0);
            plotter.addPlot(xMin, xMax, yMin, yMax, x0, y0, x0, y0, " SIGMA=2");
        }
        
        // ---- convolve y0 with itself, widening sigm by factor sqrt(2) each time -------
        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();
        float[] kernel = Arrays.copyOf(y0, y0.length);
        for (int i = 0; i < n0; ++i) {
            float convolve = kernel1DHelper.convolvePointWithKernel(y0, i, 
                kernel);
            y0[i] = convolve;
        }
        kernel = Arrays.copyOf(y0, y0.length);
        for (int i = 0; i < n0; ++i) {
            float convolve = kernel1DHelper.convolvePointWithKernel(y0, i, 
                kernel);
            y0[i] = convolve;
        }

        // normalize so that peak is at 1
        float norm0 = 1.f/MiscMath.findMax(y0);
        for (int i = 0; i < n0; ++i) {
            y0[i] *= norm0;
        }

        float fwhm2 = GaussianHelperForTests.measureFWHM(x1, y1);
        float fwhm1 = GaussianHelperForTests.measureFWHM(x0, y0);

        if (doPlot) {
            float xMax, yMax, xMin, yMin;
            
            // compare y0 with y1
            xMin = MiscMath.findMin(x0);
            xMax = MiscMath.findMax(x0);
            yMin = MiscMath.findMin(y0);
            yMax = MiscMath.findMax(y0);
            plotter.addPlot(xMin, xMax, yMin, yMax, x0, y0, x0, y0, 
                " SIGMA=2 convolved with itself twice");

            xMin = MiscMath.findMin(x1);
            xMax = MiscMath.findMax(x1);
            yMin = MiscMath.findMin(y1);
            yMax = 1.2f*Math.max(MiscMath.findMax(y1), MiscMath.findMax(y0));
            plotter.addPlot(xMin, xMax, yMin, yMax, x1, y1, x0, y0, 
                " SIGMA=4 w/ overplot of conv");

            String filePath = plotter.writeFile();
        }
        
        /*
        variance^2 = variance_a^2 + variance_b^2 + 2*covariance(a, b)
        
        for gaussian, variance = sigma^2
        */
        
        assertEquals(fwhm2, fwhm1);            
    }

    public void testCompareRecursiveConvolveFirstDerivGaussian1D() throws IOException {

        //NOTE: this shows the kernel is not recursive
        
        boolean doPlot = false;
        
        /*
        temporarily flipping the negative y's to positive for easier viewing
        */
                              
        float[] gSIGMATWO = Gaussian1DFirstDeriv.getKernel(SIGMA.TWO);
        float[] gSIGMAFOUR = Gaussian1DFirstDeriv.getKernel(SIGMA.FOUR);
        
        // reverse the polarity of the y's
        int idx0First = -1;
        int idx0Last = -1;
        for (int i = 0; i < gSIGMATWO.length; ++i) {
            if (gSIGMATWO[i] < 0) {
                if (idx0First == -1) {
                    idx0First = i;
                } else {
                    idx0Last = i;
                }
                gSIGMATWO[i] *= -1;
            }
        }
        int idx1First = -1;
        int idx1Last = -1;
        for (int i = 0; i < gSIGMAFOUR.length; ++i) {
            if (gSIGMAFOUR[i] < 0) {
                if (idx1First == -1) {
                    idx1First = i;
                } else {
                    idx1Last = i;
                }
                gSIGMAFOUR[i] *= -1;
            }
        }
        
        int n1 = 200;//2 * gSIGMAFOUR.length - 1;
        
        int ns1 = (n1 - gSIGMAFOUR.length)/2;
        
        int n0 = n1;
        int ns0 = (n1 - gSIGMATWO.length)/2;
        
        float[] x0 = new float[n0];
        float[] y0 = new float[n0];
        
        float[] x1 = new float[n1];
        float[] y1 = new float[n1];
        
        float norm1 = 1.f/MiscMath.findMax(gSIGMAFOUR);
        
        int idx = 0;
        for (int i = 0; i < n0; ++i) {
            x0[i] = i;
            if (i < ns0 || (idx > (gSIGMATWO.length - 1))) {
                y0[i] = 0;
            } else {
                y0[i] = gSIGMATWO[idx];
                idx++;
            }
        }
        
        idx = 0;
        for (int i = 0; i < n1; ++i) {
            x1[i] = i;
            if (i < ns1 || (idx > (gSIGMAFOUR.length - 1))) {
                y1[i] = 0;
            } else {
                y1[i] = norm1 * gSIGMAFOUR[idx];
                idx++;
            }
        }
        
        float fwhm0 = GaussianHelperForTests.measureFWHM(x0, y0);
               
        PolygonAndPointPlotter plotter = null;
        
        if (doPlot) {
            
            plotter = new PolygonAndPointPlotter();

            float xMax, yMax, xMin, yMin;

            xMin = MiscMath.findMin(x0);
            xMax = MiscMath.findMax(x0);
            yMin = MiscMath.findMin(y0);
            yMax = MiscMath.findMax(y0);
            plotter.addPlot(xMin, xMax, yMin, yMax, x0, y0, x0, y0, " SIGMA=2");
        }
        
        // ---- convolve y0 with sigma one, twice to result in conv with sigma*2 -------
        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();
        int cTimes = 4;
        for (int j = 0; j < cTimes; ++j) {
            float[] kernel = Arrays.copyOf(y0, y0.length);
            for (int i = 0; i < n0; ++i) {
                float convolve = kernel1DHelper.convolvePointWithKernel(y0, i, 
                    kernel);
                y0[i] = convolve;
            }
        }
        
        float fwhm2 = GaussianHelperForTests.measureFWHM(x1, y1);
        float fwhm1 = GaussianHelperForTests.measureFWHM(x0, y0);

        if (doPlot) {
            float xMax, yMax, xMin, yMin;
            
            xMin = MiscMath.findMin(x0);
            xMax = MiscMath.findMax(x0);
            yMin = MiscMath.findMin(y0);
            yMax = MiscMath.findMax(y0);
            plotter.addPlot(xMin, xMax, yMin, yMax, x0, y0, x0, y0, 
                " SIGMA=2 convolved with itself " + cTimes + " times");

            xMin = MiscMath.findMin(x1);
            xMax = MiscMath.findMax(x1);
            yMin = MiscMath.findMin(y1);
            yMax = MiscMath.findMax(y1);
            plotter.addPlot(xMin, xMax, yMin, yMax, x1, y1, x1, y1, 
                " SIGMA=4 w/ overplot of conv");

            String filePath = plotter.writeFile();
        }
        
        System.out.println("fwhm0=" + fwhm0 + " fwhm1=" + fwhm1 + " fwhm2=" + fwhm2);
        
        // can see it's not a recursive function...
        
    }
}
