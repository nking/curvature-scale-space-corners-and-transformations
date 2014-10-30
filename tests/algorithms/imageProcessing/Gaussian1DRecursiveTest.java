package algorithms.imageProcessing;

import algorithms.PolygonAndPointPlotter;
import algorithms.misc.MiscMath;
import java.io.IOException;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class Gaussian1DRecursiveTest {
    
    public Gaussian1DRecursiveTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testCompareRecursiveConvolveGaussian1D() throws IOException {
        
        System.out.println("testCompareRecursiveConvolveGaussian1D");
        
        /*
        a is a dirac delta function w/ peak at value 1
        
        a convolved with gauss(sigma=2) = a smoothed by sigma=2
        
        (a convolved with gauss(sigma=2)) convolved w/ gauss(sigma=2) = a 
           smoothed by sigma=2*sqrt(2)
       */
        
        // kernel for 2*sqrt(2) is size 13
        int n = 15;
        int h = n >> 1;
        PairIntArray a = new PairIntArray(n);
        for (int i = 0; i < n; i++) {
            int xc = i - h;
            if (i == h) {
                a.add(xc, 255);
            } else {
                a.add(xc, 0);
            }
        }

        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();
                      
        float[] gSIGMATWO = Gaussian1D.getKernel(SIGMA.TWO);
        
        // ===== convolve 'a' with a gaussian kernel of sigma=2 ========
        // ===== as a 1-D convolution with the y-axis
        // 1-D convolution of 'a' X-axis w/ a sigma=2 1D Gaussian kernel
        PairIntArray ac2 = new PairIntArray(n);
        for (int i = 0; i < n; i++) {            
            double convY = kernel1DHelper.convolvePointWithKernel(a, i, 
                gSIGMATWO, false);
            ac2.add(a.getX(i), (int)convY);
        }
        
        // ====== convolving 'ac2' w/ a sigma=2 kernel should produce same
        //        results as 'a' convolved by a sigma=2*sqrt(2).
        PairIntArray ac2c2 = new PairIntArray(n);
        for (int i = 0; i < n; i++) {            
            double convY = kernel1DHelper.convolvePointWithKernel(ac2, i, 
                gSIGMATWO, false);
            ac2c2.add(ac2.getX(i), (int)convY);
        }
       
        // ==== convolve 'a' by a sigma=2*sqrt(2) kernel
        float[] gSIGMATWOSQRT2 = Gaussian1D.getKernel(SIGMA.TWOSQRT2);
        PairIntArray ac2sqrt2 = new PairIntArray(n);
        for (int i = 0; i < n; i++) {            
            double convY = kernel1DHelper.convolvePointWithKernel(a, i, 
                gSIGMATWOSQRT2, false);
            ac2sqrt2.add(a.getX(i), (int)convY);
        }
        
        // compare  ac2c2  to  ac2sqrt2
      
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        
        float xMax, yMax, xMin, yMin;
        
        xMin = MiscMath.findMin(ac2.getX());
        xMax = MiscMath.findMax(ac2.getX());
        yMin = MiscMath.findMin(ac2.getY());
        yMax = MiscMath.findMax(ac2.getY());
        plotter.addPlot(xMin, xMax, yMin, yMax, 
            ac2.getX(), ac2.getY(), ac2.getX(), ac2.getY(),
            "a conv w/ gaussian1d(sigma=2)");
        String filePath = plotter.writeFile();
        
        xMin = MiscMath.findMin(ac2c2.getX());
        xMax = MiscMath.findMax(ac2c2.getX());
        yMin = MiscMath.findMin(ac2c2.getY());
        yMax = MiscMath.findMax(ac2c2.getY());
        plotter.addPlot(xMin, xMax, yMin, yMax, 
            ac2c2.getX(), ac2c2.getY(), ac2c2.getX(), ac2c2.getY(),
            "a conv w/ g(2) conv w/ g(2)");
        plotter.writeFile();
        
        xMin = MiscMath.findMin(ac2sqrt2.getX());
        xMax = MiscMath.findMax(ac2sqrt2.getX());
        yMin = MiscMath.findMin(ac2sqrt2.getY());
        yMax = MiscMath.findMax(ac2sqrt2.getY());
        plotter.addPlot(xMin, xMax, yMin, yMax, 
            ac2sqrt2.getX(), ac2sqrt2.getY(), 
            ac2sqrt2.getX(), ac2sqrt2.getY(),
            "a conv w/ g(2*sqrt(2))");
        plotter.writeFile();
    
         System.out.println(filePath);
         
        // assert that peaks are the same and FWHM are the same
        float yMax0 = Float.MIN_VALUE;
        float yMax1 = Float.MIN_VALUE;
        int yMax0Idx = -1;
        int yMax1Idx = -1;
        for (int i = 0; i < ac2c2.getN(); i++) {
            if (ac2c2.getY(i) > yMax0) {
                yMax0 = ac2c2.getY(i);
                yMax0Idx = i;
            }
        }
        for (int i = 0; i < ac2sqrt2.getN(); i++) {
            if (ac2sqrt2.getY(i) > yMax1) {
                yMax1 = ac2sqrt2.getY(i);
                yMax1Idx = i;
            }
        }
        assertTrue(yMax0Idx == yMax1Idx);
        assertTrue(Math.abs(yMax0 - yMax1) < 0.01);
        
        float FWHM0 = GaussianHelperForTests.measureFWHM(ac2c2.getX(), 
            ac2c2.getY());
        
        float FWHM1 = GaussianHelperForTests.measureFWHM(ac2sqrt2.getX(), 
            ac2sqrt2.getY());
                
        float expectedFWHM = (float)(SIGMA.getValue(SIGMA.TWOSQRT2) 
            * 2.f * Math.sqrt(2*Math.log(2)));
       
        float eps = 2*(ac2c2.getX(4) - ac2c2.getX(3));
        
        assertTrue(Math.abs(FWHM0 - expectedFWHM) < eps);
        assertTrue(Math.abs(FWHM1 - expectedFWHM) < eps);       
                        
    }

}
