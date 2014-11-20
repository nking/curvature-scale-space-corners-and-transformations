package algorithms.imageProcessing;

import algorithms.util.PolygonAndPointPlotter;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.util.Arrays;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * class to assert correct use of a kernel with an offset from 0
 * 
 * @author nichole
 */
public class GaussianTest {
    
    public GaussianTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    @Test
    public void testGaussian1DConvolve0() throws IOException {
        
        // examine the error in convolution.
        // convolving a dirac delta function with a gaussian that should
        //   result in a multiplication by '1'.
        // 0.005 of the area in the resulting profile is displaced into the
        // wings of the profile, so the error is pretty small for the intended
        // use.
                                      
        float sigma = 0.42466090014400953f;
        
        float mu0 = 0.0f;
        float[] gSIGMA = Gaussian1D.getKernel(sigma);
     
        float mu1 = 400.0f;
        
        // 0 1 2 3 4
        int n = 21;
        int h = n >> 1;
        PairIntArray a0 = new PairIntArray(n);
        for (int i = 0; i < n; i++) {
            int xc = i - h;
            if (i == h) {
                a0.add(xc, 1);
            } else {
                a0.add(xc, 0);
            }
        }
        
        Kernel1DHelper k1dh = new Kernel1DHelper();
        
        // keeping the 1D convolutions separate
        PairIntArray convolvedA0X = new PairIntArray(n);
        
        // convolving the X then the Y to make a 2-D set.
        PairIntArray convolvedA0XY = new PairIntArray(n);
        
        for (int i = 0; i < n; i++) {
                        
            int x = (int)k1dh.convolvePointWithKernel(a0, i, gSIGMA, true);
                        
            //int y = (int)k1dh.convolvePointWithKernel(a0, i, gSIGMA, false);
                        
            convolvedA0X.add(x, a0.getY(i));
        }
        
        for (int i = 0; i < n; i++) {
                        
            int y = (int)k1dh.convolvePointWithKernel(convolvedA0X, i, gSIGMA, 
                false);
                        
            convolvedA0XY.add(convolvedA0X.getX(i), y);
        }
        
        // make the same a0 test vector, but displaced by 400 in x
        PairIntArray a1 = a0.copy();
        
        for (int i = 0; i < a1.getN(); i++) {
            a1.set(i, (int)(i + mu1), a1.getY(i));
        }
        
        // keeping the 1D convolutions separate
        PairIntArray convolvedA1X = new PairIntArray(n);
        PairIntArray convolvedA1Y = new PairIntArray(n);
        
        // convolving the X then the Y to make a 2-D set.
        PairIntArray convolvedA1XY = new PairIntArray(n);
        
        for (int i = 0; i < a1.getN(); i++) {
                        
            int x = (int)k1dh.convolvePointWithKernel(a1, i, gSIGMA, true);
                        
            int y = (int)k1dh.convolvePointWithKernel(a1, i, gSIGMA, false);
            
            convolvedA1X.add(x, a1.getY(i));
            
            convolvedA1Y.add(a1.getX(i), y);
        }
        
        for (int i = 0; i < a1.getN(); i++) {
                        
            int y = (int)k1dh.convolvePointWithKernel(convolvedA1Y, i, gSIGMA, 
                false);
                        
            convolvedA1XY.add(convolvedA1X.getX(i), y);
        }
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        
        // plot the unconvolved functions
        // plot the unconvolved with convolved mixed (a0x with convolvedA0Y)...
        
        plotter.addPlot(MiscMath.findMin(a0.getX()), MiscMath.findMax(a0.getX()),
            MiscMath.findMin(a0.getY()), MiscMath.findMax(a0.getY()), 
            a0.getX(), a0.getY(), a0.getX(), a0.getY(), "a");
        
        plotter.addPlot(
            MiscMath.findMin(convolvedA0X.getX()), 
            MiscMath.findMax(convolvedA0X.getX()), 
            MiscMath.findMin(convolvedA0X.getY()), 
            MiscMath.findMax(convolvedA0X.getY()),
            convolvedA0X.getX(), convolvedA0X.getY(), 
            convolvedA0X.getX(), convolvedA0X.getY(),
            "a conv w/ 1D X gaussian(mu=" + mu0 + ", sigma=" + sigma 
            + ")");
        
        plotter.addPlot(
            MiscMath.findMin(convolvedA0XY.getX()), 
            MiscMath.findMax(convolvedA0XY.getX()), 
            MiscMath.findMin(convolvedA0XY.getY()), 
            MiscMath.findMax(convolvedA0XY.getY()),
            convolvedA0XY.getX(), convolvedA0XY.getY(), 
            convolvedA0XY.getX(), convolvedA0XY.getY(), 
            "a conv w/ 2D gaussian(mu=" + mu0 + ", sigma=" + sigma
            + ")");
        
        plotter.addPlot(
            MiscMath.findMin(a1.getX()),
            MiscMath.findMax(a1.getX()),
            MiscMath.findMin(a1.getY()),
            MiscMath.findMax(a1.getY()),
            a1.getX(), a1.getY(), a1.getX(), a1.getY(),
            "a offset by " + mu1);
        
        plotter.addPlot(
            mu1 - (convolvedA1X.getN() >> 1),
            //MiscMath.findMin(convolvedA1X.getX()),
            MiscMath.findMax(convolvedA1X.getX()), 
            MiscMath.findMin(convolvedA1X.getY()), 
            MiscMath.findMax(convolvedA1X.getY()),
            convolvedA1X.getX(), convolvedA1X.getY(), 
            convolvedA1X.getX(), convolvedA1X.getY(),
            "a offset by " + mu1 + 
            "conv w/ 1D X gaussian(mu=" + mu1 + ", sigma=" + sigma 
            + ")");
        
        plotter.addPlot(
            MiscMath.findMin(convolvedA1Y.getX()),
            MiscMath.findMax(convolvedA1Y.getX()), 
            MiscMath.findMin(convolvedA1Y.getY()), 
            MiscMath.findMax(convolvedA1Y.getY()),
            convolvedA1Y.getX(), convolvedA1Y.getY(),
            convolvedA1Y.getX(), convolvedA1Y.getY(),
            "a offset by " + mu1 + 
            "conv w/ 1D Y gaussian(mu=" + mu1 + ", sigma=" + sigma 
            + ")");
        
        plotter.addPlot(
            MiscMath.findMin(convolvedA1XY.getX()),
            MiscMath.findMax(convolvedA1XY.getX()), 
            MiscMath.findMin(convolvedA1XY.getY()), 
            MiscMath.findMax(convolvedA1XY.getY()),
            convolvedA1XY.getX(), convolvedA1XY.getY(),
            convolvedA1XY.getX(), convolvedA1XY.getY(),
            "a conv w/ 2D gaussian(mu=" + mu1 + ", sigma=" + sigma
            + ")");
        
        String filePath = plotter.writeFile();
        
        System.out.println(filePath);
        
        double area0 = GaussianHelperForTests.areaUnderTheCurve(
            a0.getX(), a0.getY());
        double area1 = GaussianHelperForTests.areaUnderTheCurve(
            convolvedA1XY.getX(), convolvedA1XY.getY());
        
        double ratio = (area0 < area1) ? area0/area1 : area1/area0;

        System.out.println("area under original curve=" + area0 
            + " area under convolved curve=" + area1 + " ratio=" + ratio);
        
        // the error is smaller that that due to rounding to pixels (integers)
        assertTrue(ratio > 0.94);
           
        float eps = 2.f*(convolvedA1X.getX(1) - convolvedA1X.getX(0));
        
        float expectedFWHM = GaussianHelperForTests.expectedFWHM(sigma);
        
        float measuredFWHMLarge = GaussianHelperForTests.measureFWHM(
            convolvedA1XY.getX(), convolvedA1XY.getY());
        
        float measuredFWHMSmall = 
            GaussianHelperForTests.measureFWHMRoundToSmaller(
                convolvedA1XY.getX(), convolvedA1XY.getY());
        
        float measuredFWHM = (measuredFWHMLarge + measuredFWHMSmall)/2.f;
         
        System.out.println("Gaussian w/ sigma=" + sigma + " expected FWHM=" + 
            expectedFWHM + " measured=" + measuredFWHM);
        
        assertTrue(Math.abs(expectedFWHM - measuredFWHM) < eps);
        
    }
    
    @Test
    public void testGaussian1DConvolve1() throws IOException {
        
        // examine convolution w/ first deriv of Gaussian.
        // convolving a dirac delta function with a gaussian that should
        //   result in a multiplication by '1'.
                                 
        float sigma = 0.42466090014400953f;
        
        float mu0 = 0.0f;
        float[] gSIGMA = Gaussian1DFirstDeriv.getKernel(sigma);
     
        float mu1 = 400.0f;
        
        // 0 1 2 3 4
        int n = 21;
        int h = n >> 1;
        PairIntArray a0 = new PairIntArray(n);
        for (int i = 0; i < n; i++) {
            int xc = i - h;
            if (i == h) {
                a0.add(xc, 1);
            } else {
                a0.add(xc, 0);
            }
        }
         
        // make the same a0 test vector, but displaced by 400 in x
        PairIntArray a1 = a0.copy();
        
        for (int i = 0; i < a1.getN(); i++) {
            a1.set(i, (int)(i + mu1), a1.getY(i));
        }
                
        Kernel1DHelper k1dh = new Kernel1DHelper();
        
        // keeping the 1D convolutions separate
        PairIntArray convolvedA0X = new PairIntArray(n);
        
        // convolving the X then the Y to make a 2-D set.
        PairIntArray convolvedA0XY = new PairIntArray(n);
        
        for (int i = 0; i < n; i++) {
            int x = (int)k1dh.convolvePointWithKernel(a0, i, gSIGMA, true);                                                
            convolvedA0X.add(x, a0.getY(i));
        }
        
        for (int i = 0; i < n; i++) {
            int y = (int)k1dh.convolvePointWithKernel(convolvedA0X, i, gSIGMA, 
                false);                        
            convolvedA0XY.add(convolvedA0X.getX(i), y);
        }
       
        // keeping the 1D convolutions separate
        PairIntArray convolvedA1X = new PairIntArray(n);
        PairIntArray convolvedA1Y = new PairIntArray(n);
        
        // convolving the X then the Y to make a 2-D set.
        PairIntArray convolvedA1XY = new PairIntArray(n);
        
        for (int i = 0; i < a1.getN(); i++) {
            int x = (int)k1dh.convolvePointWithKernel(a1, i, gSIGMA, true);                        
            int y = (int)k1dh.convolvePointWithKernel(a1, i, gSIGMA, false);
            convolvedA1X.add(x, a1.getY(i));            
            convolvedA1Y.add(a1.getX(i), y);
        }
        
        for (int i = 0; i < a1.getN(); i++) {
            int y = (int)k1dh.convolvePointWithKernel(convolvedA1X, i, gSIGMA, 
                false);                        
            convolvedA1XY.add(convolvedA1X.getX(i), y);
        }
       
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        
        // plot the unconvolved functions
        // plot the unconvolved with convolved mixed (a0x with convolvedA0Y)...
        
        plotter.addPlot(MiscMath.findMin(a0.getX()), MiscMath.findMax(a0.getX()),
            MiscMath.findMin(a0.getY()), MiscMath.findMax(a0.getY()), 
            a0.getX(), a0.getY(), a0.getX(), a0.getY(), "a");
        
        plotter.addPlot(
            MiscMath.findMin(convolvedA0X.getX()), 
            MiscMath.findMax(convolvedA0X.getX()), 
            MiscMath.findMin(convolvedA0X.getY()), 
            MiscMath.findMax(convolvedA0X.getY()),
            convolvedA0X.getX(), convolvedA0X.getY(), 
            convolvedA0X.getX(), convolvedA0X.getY(),
            "a conv w/ 1D X gaussianFirstDeriv(mu=" + mu0 + ", sigma=" + sigma 
            + ")");
        
        plotter.addPlot(
            MiscMath.findMin(convolvedA0XY.getX()), 
            MiscMath.findMax(convolvedA0XY.getX()), 
            MiscMath.findMin(convolvedA0XY.getY()), 
            MiscMath.findMax(convolvedA0XY.getY()),
            convolvedA0XY.getX(), convolvedA0XY.getY(), 
            convolvedA0XY.getX(), convolvedA0XY.getY(), 
            "a conv w/ 2D gaussianFirstDeriv(mu=" + mu0 + ", sigma=" + sigma
            + ")");
        
        plotter.addPlot(
            MiscMath.findMin(a1.getX()),
            MiscMath.findMax(a1.getX()),
            MiscMath.findMin(a1.getY()),
            MiscMath.findMax(a1.getY()),
            a1.getX(), a1.getY(), a1.getX(), a1.getY(),
            "a offset by " + mu1);
        
        plotter.addPlot(
            MiscMath.findMin(convolvedA1X.getX()), 
            MiscMath.findMax(convolvedA1X.getX()), 
            MiscMath.findMin(convolvedA1X.getY()), 
            MiscMath.findMax(convolvedA1X.getY()),
            convolvedA1X.getX(), convolvedA1X.getY(), 
            convolvedA1X.getX(), convolvedA1X.getY(),
            "a offset by " + mu1 + 
            "conv w/ 1D X gaussianFirstDeriv(mu=" + mu1 + ", sigma=" + sigma 
            + ")");
        
        plotter.addPlot(
            MiscMath.findMin(convolvedA1Y.getX()),
            MiscMath.findMax(convolvedA1Y.getX()), 
            MiscMath.findMin(convolvedA1Y.getY()), 
            MiscMath.findMax(convolvedA1Y.getY()),
            convolvedA1Y.getX(), convolvedA1Y.getY(),
            convolvedA1Y.getX(), convolvedA1Y.getY(),
            "a offset by " + mu1 + 
            "conv w/ 1D Y gaussianFirstDeriv(mu=" + mu1 + ", sigma=" + sigma 
            + ")");
        
        plotter.addPlot(
            MiscMath.findMin(convolvedA1XY.getX()),
            MiscMath.findMax(convolvedA1XY.getX()), 
            MiscMath.findMin(convolvedA1XY.getY()), 
            MiscMath.findMax(convolvedA1XY.getY()),
            convolvedA1XY.getX(), convolvedA1XY.getY(),
            convolvedA1XY.getX(), convolvedA1XY.getY(),
            "a conv w/ 2D gaussianFirstDeriv(mu=" + mu1 + ", sigma=" + sigma
            + ")");
        
        String filePath = plotter.writeFile();
        
        System.out.println(filePath);
        
        double area0 = GaussianHelperForTests.areaUnderTheCurve(
            a0.getX(), a0.getY());
        double area1 = GaussianHelperForTests.areaUnderTheCurve(
            convolvedA1XY.getX(), convolvedA1XY.getY());
        
        System.out.println("a1 x = " + Arrays.toString(a1.getX()));
        System.out.println("a1 y = " + Arrays.toString(a1.getY()));
        
        System.out.println("convolvedA0X x = " + Arrays.toString(
            convolvedA0X.getX()));
        System.out.println("convolvedA1X y = " + Arrays.toString(
            convolvedA0X.getY()));
        
        System.out.println("convolvedA1X x = " + Arrays.toString(
            convolvedA1X.getX()));
        System.out.println("convolvedA1X y = " + Arrays.toString(
            convolvedA1X.getY()));
        System.out.println("convolvedA1Y x = " + Arrays.toString(
            convolvedA1Y.getX()));
        System.out.println("convolvedA1Y y = " + Arrays.toString(
            convolvedA1Y.getY()));
        
        System.out.println("convolvedA1XY x = " + Arrays.toString(
            convolvedA1XY.getX()));
        System.out.println("convolvedA1XY y = " + Arrays.toString(
            convolvedA1XY.getY()));
        
        System.out.println("area under original curve=" + area0 
            + " area under convolved curve=" + area1);
        
        // the error is smaller that that due to rounding to pixels (integers)
        //assertTrue(Math.abs(area1 - 0) < 0.01);
    }
    
    @Test
    public void testGaussian1DConvolve2() throws IOException {
        
        float sigma = 0.42466090014400953f;
        
        float mu0 = 0.0f;
        float[] gSIGMA = Gaussian1DSecondDeriv.getKernel(sigma);
     
        float mu1 = 400.0f;
        
        // 0 1 2 3 4
        int n = 21;
        int h = n >> 1;
        PairIntArray a0 = new PairIntArray(n);
        for (int i = 0; i < n; i++) {
            int xc = i - h;
            if (i == h) {
                a0.add(xc, 1);
            } else {
                a0.add(xc, 0);
            }
        }
        
        Kernel1DHelper k1dh = new Kernel1DHelper();
        
        // keeping the 1D convolutions separate
        PairIntArray convolvedA0X = new PairIntArray(n);
        
        // convolving the X then the Y to make a 2-D set.
        PairIntArray convolvedA0XY = new PairIntArray(n);
        
        for (int i = 0; i < n; i++) {
                        
            int x = (int)k1dh.convolvePointWithKernel(a0, i, gSIGMA, true);
                        
            //int y = (int)k1dh.convolvePointWithKernel(a0, i, gSIGMA, false);
                        
            convolvedA0X.add(x, a0.getY(i));
        }
        
        for (int i = 0; i < n; i++) {
                        
            int y = (int)k1dh.convolvePointWithKernel(convolvedA0X, i, gSIGMA, 
                false);
                        
            convolvedA0XY.add(convolvedA0X.getX(i), y);
        }
        
        // make the same a0 test vector, but displaced by 400 in x
        PairIntArray a1 = a0.copy();
        
        for (int i = 0; i < a1.getN(); i++) {
            a1.set(i, (int)(i + mu1), a1.getY(i));
        }
        
        // keeping the 1D convolutions separate
        PairIntArray convolvedA1X = new PairIntArray(n);
        PairIntArray convolvedA1Y = new PairIntArray(n);
        
        // convolving the X then the Y to make a 2-D set.
        PairIntArray convolvedA1XY = new PairIntArray(n);
        
        for (int i = 0; i < a1.getN(); i++) {
                        
            int x = (int)k1dh.convolvePointWithKernel(a1, i, gSIGMA, true);
                        
            int y = (int)k1dh.convolvePointWithKernel(a1, i, gSIGMA, false);
            
            convolvedA1X.add(x, a1.getY(i));
            
            convolvedA1Y.add(a1.getX(i), y);
        }
        
        for (int i = 0; i < a1.getN(); i++) {
                        
            int y = (int)k1dh.convolvePointWithKernel(convolvedA1X, i, gSIGMA, 
                false);
                        
            convolvedA1X.add(convolvedA1X.getX(i), y);
        }
       
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        
        // plot the unconvolved functions
        // plot the unconvolved with convolved mixed (a0x with convolvedA0Y)...
        
        plotter.addPlot(MiscMath.findMin(a0.getX()), MiscMath.findMax(a0.getX()),
            MiscMath.findMin(a0.getY()), MiscMath.findMax(a0.getY()), 
            a0.getX(), a0.getY(), a0.getX(), a0.getY(), "a");
        
        plotter.addPlot(
            MiscMath.findMin(convolvedA0X.getX()), 
            MiscMath.findMax(convolvedA0X.getX()), 
            MiscMath.findMin(convolvedA0X.getY()), 
            MiscMath.findMax(convolvedA0X.getY()),
            convolvedA0X.getX(), convolvedA0X.getY(), 
            convolvedA0X.getX(), convolvedA0X.getY(),
            "a conv w/ 1D X gaussianSecondDeriv(mu=" + mu0 + ", sigma=" + sigma 
            + ")");
        
        plotter.addPlot(
            MiscMath.findMin(convolvedA0XY.getX()), 
            MiscMath.findMax(convolvedA0XY.getX()), 
            MiscMath.findMin(convolvedA0XY.getY()), 
            MiscMath.findMax(convolvedA0XY.getY()),
            convolvedA0XY.getX(), convolvedA0XY.getY(), 
            convolvedA0XY.getX(), convolvedA0XY.getY(), 
            "a conv w/ 2D gaussianSecondDeriv(mu=" + mu0 + ", sigma=" + sigma
            + ")");
        
        plotter.addPlot(
            MiscMath.findMin(a1.getX()),
            MiscMath.findMax(a1.getX()),
            MiscMath.findMin(a1.getY()),
            MiscMath.findMax(a1.getY()),
            a1.getX(), a1.getY(), a1.getX(), a1.getY(),
            "a offset by " + mu1);
        
        plotter.addPlot(
            MiscMath.findMin(convolvedA1X.getX()), 
            MiscMath.findMax(convolvedA1X.getX()), 
            MiscMath.findMin(convolvedA1X.getY()), 
            MiscMath.findMax(convolvedA1X.getY()),
            convolvedA1X.getX(), convolvedA1X.getY(), 
            convolvedA1X.getX(), convolvedA1X.getY(),
            "a offset by " + mu1 + 
            "conv w/ 1D X gaussianSecondDeriv(mu=" + mu1 + ", sigma=" + sigma 
            + ")");
        
        plotter.addPlot(
            MiscMath.findMin(convolvedA1Y.getX()),
            MiscMath.findMax(convolvedA1Y.getX()), 
            MiscMath.findMin(convolvedA1Y.getY()), 
            MiscMath.findMax(convolvedA1Y.getY()),
            convolvedA1Y.getX(), convolvedA1Y.getY(),
            convolvedA1Y.getX(), convolvedA1Y.getY(),
            "a offset by " + mu1 + 
            "conv w/ 1D Y gaussianSecondDeriv(mu=" + mu1 + ", sigma=" + sigma 
            + ")");
        
        plotter.addPlot(
            MiscMath.findMin(convolvedA1XY.getX()),
            MiscMath.findMax(convolvedA1XY.getX()), 
            MiscMath.findMin(convolvedA1XY.getY()), 
            MiscMath.findMax(convolvedA1XY.getY()),
            convolvedA1XY.getX(), convolvedA1XY.getY(),
            convolvedA1XY.getX(), convolvedA1XY.getY(),
            "a conv w/ 2D gaussianSecondDeriv(mu=" + mu1 + ", sigma=" + sigma
            + ")");
        
        String filePath = plotter.writeFile();
        
        System.out.println(filePath);
        
        double area0 = GaussianHelperForTests.areaUnderTheCurve(
            a0.getX(), a0.getY());
        double area1 = GaussianHelperForTests.areaUnderTheCurve(
            convolvedA1XY.getX(), convolvedA1XY.getY());
        
        System.out.println("area under original curve=" + area0 
            + " area under convolved curve=" + area1);
        
        // the error is smaller that that due to rounding to pixels (integers)
        //assertTrue(Math.abs(area1 - 1) < 0.01);
    }
    
}
