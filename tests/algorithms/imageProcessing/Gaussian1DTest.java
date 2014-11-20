package algorithms.imageProcessing;

import algorithms.util.PolygonAndPointPlotter;
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
public class Gaussian1DTest {
    
    public Gaussian1DTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testEstimateHWZI() {
                        
        SIGMA sigma = SIGMA.TWO;
        
        float s = SIGMA.getValue(sigma);
        assertTrue(s == 2);
        
        float t = 0.0F;
        
        float fractionMax = 0.5F;
                
        float expResult = (float)(s * Math.sqrt(2.f * Math.log(2)) + t);
        
        float hwhmHalf = Gaussian1D.estimateHWZI(sigma, fractionMax);
        
        assertTrue(expResult == hwhmHalf);
        
        float result = Gaussian1D.estimateHWZI(SIGMA.ONE, 0.01f);
        
        assertTrue(Math.abs(result - 3) < 0.2);
        assertTrue(result > 2.6);
        
    }

    @Test
    public void testGetNormalization() {
                                
        float[] g = Gaussian1D.getKernel(SIGMA.ONE);
        
        float sum = 0;
        for (float gg : g) {
            sum += gg;
        }
        
        assertTrue(Math.round(sum) == 1);
    }
  
    @Test
    public void testCheckAllBinomialKernels() throws IOException {
                        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        
        for (SIGMA sigma : SIGMA.values()) {
            
            float[] gb = Gaussian1D.getBinomialKernel(sigma);
            if (gb == null) {
                continue;
            }
            float mu = 0;
            int halfWidthInPixels = gb.length >> 1; 
            float[] xb = new float[gb.length];
            for (int i = 0; i < xb.length; i++) {
                xb[i] = i + (mu - halfWidthInPixels);
            }
            
            float[] g = Gaussian1D.getKernel(SIGMA.getValue(sigma));
            halfWidthInPixels = g.length >> 1; 
            float[] x = new float[g.length];
            for (int i = 0; i < x.length; i++) {
                x[i] = i + (mu - halfWidthInPixels);
            }
            
            System.out.println("sigma=" + sigma);
            
            float yMin = MiscMath.findMin(gb);
            if (yMin > 0) {
                yMin = 0;
            } else if (yMin < 0) {
                yMin *= 0.9f;
            }
            float yMin2 = MiscMath.findMin(g);
            if (yMin2 < yMin) {
                yMin = yMin2;
            }
            
            float yMax = 1.1f*MiscMath.findMax(gb);
            float yMax2 = MiscMath.findMax(g);
            if (yMax2 > yMax) {
                yMax = yMax2;
            }
            float xMax = MiscMath.findMax(xb);
            if (xMax < 1) {
                xMax = 1;
            }
            float xMax2 = MiscMath.findMax(x);
            if (xMax2 > xMax) {
                xMax = xMax2;
            }
            float xMin = MiscMath.findMin(xb);
            if (xMin > 1) {
                xMin = 0;
            }
            float xMin2 = MiscMath.findMin(x);
            if (xMin2 < xMin) {
                xMin = xMin2;
            }
            
            double area = GaussianHelperForTests.areaUnderTheCurve(xb, gb);
            float FWHM = GaussianHelperForTests.measureFWHM(xb, gb);
            float expectedFWHM = GaussianHelperForTests.expectedFWHM(
                SIGMA.getValue(sigma));
            float eps;
            if (xb.length > 3) {
                eps = 2.f*((xb[3] - xb[2]) + (xb[2] - xb[1]) + (xb[1] - xb[0]))/3.f;
            } else if (xb.length == 2) {
                eps = 2.f*(xb[1] - xb[0]);
            } else {
                eps = 1;
            }
            
            plotter.addPlot(xMin, xMax, yMin, yMax, xb, gb, x, g,
                "gauss(sigma=" + sigma.toString() + " area=" + area +")");
            plotter.writeFile();
            
            assertTrue(Math.abs(area - 1) < 0.1);
            if (FWHM > -1) {
                float diff = Math.abs(FWHM - expectedFWHM);
                assertTrue(diff < eps);
            }
        }
        int z = 1;
        for (SIGMA sigma : SIGMA.values()) {
            
            float[] gb = Gaussian1DFirstDeriv.getBinomialKernel(sigma);
            if (gb == null) {
                continue;
            }
            float mu = 0;
            int halfWidthInPixels = gb.length >> 1; 
            float[] xb = new float[gb.length];
            for (int i = 0; i < xb.length; i++) {
                xb[i] = i + (mu - halfWidthInPixels);
            }
            
            float[] g = Gaussian1DFirstDeriv.getKernel(SIGMA.getValue(sigma));
            halfWidthInPixels = g.length >> 1; 
            float[] x = new float[g.length];
            for (int i = 0; i < x.length; i++) {
                x[i] = i + (mu - halfWidthInPixels);
            }
            
            float yMin = MiscMath.findMin(gb);
            if (yMin > 0) {
                yMin = 0;
            } else if (yMin < 0) {
                yMin *= 0.9f;
            }
            float yMin2 = MiscMath.findMin(g);
            if (yMin2 < yMin) {
                yMin = yMin2;
            }
            
            float yMax = 1.1f*MiscMath.findMax(gb);
            float yMax2 = MiscMath.findMax(g);
            if (yMax2 > yMax) {
                yMax = yMax2;
            }
            float xMax = MiscMath.findMax(xb);
            if (xMax < 1) {
                xMax = 1;
            }
            float xMax2 = MiscMath.findMax(x);
            if (xMax2 > xMax) {
                xMax = xMax2;
            }
            float xMin = MiscMath.findMin(xb);
            if (xMin > 1) {
                xMin = 0;
            }
            float xMin2 = MiscMath.findMin(x);
            if (xMin2 < xMin) {
                xMin = xMin2;
            }
            
            double area = GaussianHelperForTests.areaUnderTheCurve(xb, gb);
            
            plotter.addPlot(xMin, xMax, yMin, yMax, xb, gb, x, g,
                "gaussFirstDeriv(sigma=" + sigma.toString() + " area=" + area 
                + ")");
            plotter.writeFile();
            
            z = 1;
            //assertTrue(Math.abs(area) < 0.1);
        }
        z = 1;
        for (SIGMA sigma : SIGMA.values()) {
            
            float[] gb = Gaussian1DSecondDeriv.getBinomialKernel(sigma);
            if (gb == null) {
                continue;
            }
            float mu = 0;
            int halfWidthInPixels = gb.length >> 1; 
            float[] xb = new float[gb.length];
            for (int i = 0; i < xb.length; i++) {
                xb[i] = i + (mu - halfWidthInPixels);
            }
            
            float[] g = Gaussian1DSecondDeriv.getKernel(SIGMA.getValue(sigma));
            halfWidthInPixels = g.length >> 1; 
            float[] x = new float[g.length];
            for (int i = 0; i < x.length; i++) {
                x[i] = i + (mu - halfWidthInPixels);
            }
            
            float yMin = MiscMath.findMin(gb);
            if (yMin > 0) {
                yMin = 0;
            } else if (yMin < 0) {
                yMin *= 0.9f;
            }
            float yMin2 = MiscMath.findMin(g);
            if (yMin2 < yMin) {
                yMin = yMin2;
            }
            
            float yMax = 1.1f*MiscMath.findMax(gb);
            float yMax2 = MiscMath.findMax(g);
            if (yMax2 > yMax) {
                yMax = yMax2;
            }
            float xMax = MiscMath.findMax(xb);
            if (xMax < 1) {
                xMax = 1;
            }
            float xMax2 = MiscMath.findMax(x);
            if (xMax2 > xMax) {
                xMax = xMax2;
            }
            float xMin = MiscMath.findMin(xb);
            if (xMin > 1) {
                xMin = 0;
            }
            float xMin2 = MiscMath.findMin(x);
            if (xMin2 < xMin) {
                xMin = xMin2;
            }
            
            double area = GaussianHelperForTests.areaUnderTheCurve(xb, gb);
            
            plotter.addPlot(xMin, xMax, yMin, yMax, xb, gb, x, g,
                "gaussSecondDeriv(sigma=" + sigma.toString() + " area=" 
                    + area +")");
            plotter.writeFile();
            
            assertTrue(Math.abs(area) < 0.1);
        }
        
        z = 1;
    }
    
    @Test
    public void testCheckAllKernels() throws IOException {
        
        System.out.println("testCheckAllKernels");
                
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        
        for (SIGMA sigma : SIGMA.values()) {
            
            float[] gb = Gaussian1D.getKernel(sigma);
            if (gb == null) {
                continue;
            }
            float mu = 0;
            int halfWidthInPixels = gb.length >> 1; 
            float[] xb = new float[gb.length];
            for (int i = 0; i < xb.length; i++) {
                xb[i] = i + (mu - halfWidthInPixels);
            }
         
            System.out.println("sigma=" + sigma);
            
            float yMin = MiscMath.findMin(gb);
            if (yMin > 0) {
                yMin = 0;
            } else if (yMin < 0) {
                yMin *= 1.1f;
            }
            float yMax = 1.1f*MiscMath.findMax(gb);            
            float xMax = MiscMath.findMax(xb);
            if (xMax < 1) {
                xMax = 1;
            }
            float xMin = MiscMath.findMin(xb);
            if (xMin > 1) {
                xMin = 0;
            }
            
            double area = GaussianHelperForTests.areaUnderTheCurve(xb, gb);
            float FWHM = GaussianHelperForTests.measureFWHM(xb, gb);
            float expectedFWHM = GaussianHelperForTests.expectedFWHM(
                SIGMA.getValue(sigma));
            float eps;
            if (xb.length > 3) {
                eps = 2.f * ((xb[3] - xb[2]) + (xb[2] - xb[1]) + (xb[1] - xb[0]))/3.f;
            } else if (xb.length == 2) {
                eps = 2.f * (xb[1] - xb[0]);
            } else {
                eps = 1;
            }
            
            plotter.addPlot(xMin, xMax, yMin, yMax, null, null, xb, gb,
                "gauss(sigma=" + sigma.toString() + " area=" + area +")");
            plotter.writeFile();
            
            assertTrue(Math.abs(area - 1) < 0.1);
            if (FWHM > -1) {
                float diff = Math.abs(FWHM - expectedFWHM);
                assertTrue(diff < eps);
            }
        }
        int z = 1;
        for (SIGMA sigma : SIGMA.values()) {
            
            float[] gb = Gaussian1DFirstDeriv.getKernel(sigma);
            if (gb == null) {
                continue;
            }
            float mu = 0;
            int halfWidthInPixels = gb.length >> 1; 
            float[] xb = new float[gb.length];
            for (int i = 0; i < xb.length; i++) {
                xb[i] = i + (mu - halfWidthInPixels);
            }
            
            float yMin = MiscMath.findMin(gb);
            if (yMin > 0) {
                yMin = 0;
            } else if (yMin < 0) {
                yMin *= 1.1f;
            }
            float yMax = 1.1f*MiscMath.findMax(gb);            
            float xMax = MiscMath.findMax(xb);
            if (xMax < 1) {
                xMax = 1;
            }
            float xMin = MiscMath.findMin(xb);
            if (xMin > 1) {
                xMin = 0;
            }
            
            double area = GaussianHelperForTests.areaUnderTheCurve(xb, gb);
            
            plotter.addPlot(xMin, xMax, yMin, yMax, null, null, xb, gb,
                "gaussFirstDeriv(sigma=" + sigma.toString() + " area=" + area 
                + ")");
            plotter.writeFile();
            
            assertTrue(Math.abs(area) < 0.1);
        }
        z = 1;
        for (SIGMA sigma : SIGMA.values()) {
            
            float[] gb = Gaussian1DSecondDeriv.getKernel(sigma);
            if (gb == null) {
                continue;
            }
            float mu = 0;
            int halfWidthInPixels = gb.length >> 1; 
            float[] xb = new float[gb.length];
            for (int i = 0; i < xb.length; i++) {
                xb[i] = i + (mu - halfWidthInPixels);
            }
            
            float yMin = MiscMath.findMin(gb);
            if (yMin > 0) {
                yMin = 0;
            } else if (yMin < 0) {
                yMin *= 1.1f;
            }
            float yMax = 1.1f*MiscMath.findMax(gb);            
            float xMax = MiscMath.findMax(xb);
            if (xMax < 1) {
                xMax = 1;
            }
            float xMin = MiscMath.findMin(xb);
            if (xMin > 1) {
                xMin = 0;
            }
            
            double area = GaussianHelperForTests.areaUnderTheCurve(xb, gb);
if (sigma.ordinal() == SIGMA.TWO.ordinal()) {
    for (int i = 0; i < gb.length; i++) {
        float gs = gb[i];
        if (Math.abs(gs) < 1E-3) {
            System.out.print(", (float)" + gs);
        } else {
            System.out.print(", " + gs + "f");
        }
    }
    System.out.println("");
}          
            plotter.addPlot(xMin, xMax, yMin, yMax, null, null, xb, gb,
                "gaussSecondDeriv(sigma=" + sigma.toString() + " area=" 
                    + area +")");
            plotter.writeFile();
            
            assertTrue(Math.abs(area) < 0.1);
        }
        
        z = 1;
    }
    
    @Test
    public void testGetKernels() throws IOException {
        
        System.out.println("testGetKernels");
                
        SIGMA sigma = SIGMA.ONE;

        while ((sigma != null) && (sigma.compareTo(
            SIGMA.ONEHUNDREDANDTWENTYEIGHTSQRT2) < 1)) {
            
            float[] g = Gaussian1D.getKernel(sigma);
            
            float[] x = new float[g.length];
            
            int h = g.length >> 1;
            
            int yMaxIdx = -1;
            float yMax = Float.MIN_VALUE;
            float sum = 0;
            for (int i = 0; i < g.length; i++) {
                int xc = i - h;
                x[i] = xc;
                sum += g[i];
                if (g[i] > yMax) {
                    yMax = g[i];
                    yMaxIdx = i;
                }
            }
            assertTrue(Math.abs(sum - 1.0f) < 0.1);
            
            int yHalfMax0Idx = -1;
            int yHalfMax1Idx = -1;
            for (int i = 0; i < yMaxIdx; i++) {
                if (g[i] < (yMax/2.)) {
                    yHalfMax0Idx = i;
                }
            }
            for (int i = (g.length - 1); i > yMaxIdx; i--) {
                if (g[i] < (yMax/2.)) {
                    yHalfMax1Idx = i;
                }
            }
            
            // error of up to 1 pixel rounding on each side
            float eps = 2.f*((x[3] - x[2]) + (x[2] - x[1]) + (x[1] - x[0]))/3.f;
            float FWHM = x[yHalfMax1Idx] - x[yHalfMax0Idx];
            
            //exp(-(x_0 - mu)^2/(2*o~^2)) = 1/2
            //-(x_0 - mu)^2/(2*o~^2) = -ln(2)
            //(x_0 - mu)^2 = ln(2) * (2*o~^2)
            //(x_0 - mu) = sqrt(ln(2) * 2*o~^2)
            // x_0 = mu + o~*sqrt(ln(2) * 2)
            float expectedFWHM = (float)(SIGMA.getValue(sigma) * 2.f * 
                Math.sqrt(2*Math.log(2)));
            
            System.out.println("sigma=" + SIGMA.getValue(sigma) +
                " FWHM=" + FWHM + " expected=" + expectedFWHM + " eps=" + eps);
            
            assertTrue(Math.abs(expectedFWHM - FWHM) < eps);
            
            sigma = SIGMA.increment(sigma);            
        }
        
    }
    
    @Test
    public void testGetKernel() throws IOException {
        
        System.out.println("testGetKernel");
        
        float t = 0;
        
        SIGMA sigma = SIGMA.FOUR;
        //float s = SIGMA.getValue(SIGMA.ZEROPOINTSEVENONE);//0.42466090014400953f;
        
        boolean usePreprepared = false;
        
        float[] g;
        if (usePreprepared) {
            g = Gaussian1D.getKernel(sigma);
        } else {
            g = Gaussian1D.getKernel(sigma);
        }
        
        float[] x = new float[g.length];
        
        int h = g.length >> 1;
               
        float sum = 0;
        for (int i = 0; i < g.length; i++) {
            int xc = i - h;
            x[i] = xc;
            sum += g[i];
        }
                
        float yMax = MiscMath.findMax(g);
        float xMax = MiscMath.findMax(x);
        
        PolygonAndPointPlotter plotter = 
            new PolygonAndPointPlotter(-1*xMax - 1, xMax + 1, 0, yMax);
        
        double area = GaussianHelperForTests.areaUnderTheCurve(x, g);
        
        plotter.addPlot(x, g, x, g, "gaussian 1-d for sigma="  + sigma + " sum="
            + sum + " area=" + area);
        
        String filePath = plotter.writeFile();
        
        System.out.println(filePath);
        
        assertTrue(Math.abs(sum - 1) < 0.1);
        
        assertTrue(Math.abs(area - 1) < 0.1);
        
        //=======================================
        
        int nKernelPoints = g.length;
        
        if ((nKernelPoints & 1) == 0) {
            nKernelPoints--;
        }
        
        float[] gfd;
        if (usePreprepared) {
            gfd = Gaussian1DFirstDeriv.getKernel(sigma);
        } else {
            gfd = Gaussian1DFirstDeriv.getKernel(sigma, t, nKernelPoints);
        }
        
        x = new float[gfd.length];
                
        h = gfd.length >> 1;
        sum = 0;
        for (int i = 0; i < gfd.length; i++) {           
            x[i] = (i - h);
            sum += gfd[i];
        }
                
        yMax = MiscMath.findMax(gfd);
        float yMin = MiscMath.findMin(gfd);
        float xMin = MiscMath.findMin(x);
        xMax = MiscMath.findMax(x);
        
        area = GaussianHelperForTests.areaUnderTheCurve(x, gfd);
        
        plotter = 
            new PolygonAndPointPlotter(xMin - 1, xMax + 1, yMin, yMax);
        
        plotter.addPlot(x, gfd, x, gfd, 
            "gaussian 1-d first derivative for sigma=" + sigma + " sum=" + sum
             + " area=" + area);
        
        filePath = plotter.writeFile2();
        
        System.out.println(filePath);
        
        assertTrue(Math.abs(sum) < 0.1);
        
        assertTrue(Math.abs(area) < 0.1);
        
        //=========================================
                
        float[] gsd;
        if (usePreprepared) {
            gsd = Gaussian1DSecondDeriv.getKernel(sigma);
        } else {
            gsd = Gaussian1DSecondDeriv.getKernel(sigma, t, nKernelPoints);
        }
        
        x = new float[gsd.length];
        
        h = gsd.length >> 1;
        sum = 0;
        for (int i = 0; i < gfd.length; i++) {             
            x[i] = (i - h);
            sum += gsd[i];
        }
        
        yMax = MiscMath.findMax(gsd);
        yMin = MiscMath.findMin(gsd);
        xMin = MiscMath.findMin(x);
        xMax = MiscMath.findMax(x);
        
        area = GaussianHelperForTests.areaUnderTheCurve(x, gsd);
        
        plotter = 
            new PolygonAndPointPlotter(xMin - 1, xMax + 1, yMin, yMax);
        
        plotter.addPlot(x, gsd, x, gsd, 
            "gaussian 1-d second derivative for sigma=" + sigma + " sum=" + sum
            + " area=" + area);
        
        filePath = plotter.writeFile3();
        
        System.out.println(filePath);
                
        assertTrue(Math.abs(sum) < 0.1);
        
        assertTrue(Math.abs(area) < 0.1);
                
    }
    
    @Test
    public void testGetKernel2() throws IOException {
        
        System.out.println("testGetKernel2");
                
        SIGMA sigma = SIGMA.ONE;
                        
        float[] g = Gaussian1D.getKernel(sigma);
        
        float[] x = new float[g.length];
        
        int h = g.length >> 1;
               
        float sum = 0;
        for (int i = 0; i < g.length; i++) {
            int xc = i - h;
            x[i] = xc;
            sum += g[i];
        }
                
        float yMax = MiscMath.findMax(g);
        float xMax = MiscMath.findMax(x);
        
        PolygonAndPointPlotter plotter = 
            new PolygonAndPointPlotter(-1*xMax - 1, xMax + 1, 0, yMax);
        
        plotter.addPlot(x, g, x, g, "gaussian 1-d for sigma="  + sigma + " sum=" 
            + sum);
        
        String filePath = plotter.writeFile();
        
        System.out.println(filePath);
        
        assertTrue(Math.abs(sum - 1) < 0.1);
        
        int yMaxIdx = MiscMath.findYMaxIndex(g);
        assertTrue(Math.abs(g[yMaxIdx] - 0.4) < 0.05f);
        assertTrue(x[yMaxIdx] == 0);
        
        //=======================================
       
        float[] gfd = Gaussian1DFirstDeriv.getKernel(sigma);
        
        x = new float[gfd.length];
                
        h = gfd.length >> 1;
        sum = 0;
        for (int i = 0; i < gfd.length; i++) {           
            x[i] = (i - h);
            sum += gfd[i];
        }
        
        assertTrue(Math.abs(sum) < 0.1);
        
        yMax = MiscMath.findMax(gfd);
        float yMin = MiscMath.findMin(gfd);
        float xMin = MiscMath.findMin(x);
        xMax = MiscMath.findMax(x);
        
        plotter = 
            new PolygonAndPointPlotter(xMin - 1, xMax + 1, yMin, yMax);
        
        plotter.addPlot(x, gfd, x, gfd, 
            "gaussian 1-d first derivative for sigma=" + sigma.toString());
        
        filePath = plotter.writeFile2();
        
        System.out.println(filePath);
        
        assertTrue(x[h] == 0);
        assertTrue(gfd[h] == 0);
        assertTrue(Math.abs(sum) < 0.1);
        yMaxIdx = MiscMath.findYMaxIndex(gfd);
        assertTrue(Math.abs(gfd[yMaxIdx] - 0.24) < 0.05f);
        assertTrue(x[yMaxIdx] == -1);
        int yMinIdx = MiscMath.findYMinIndex(gfd);
        assertTrue(Math.abs(gfd[yMinIdx] - -0.24) < 0.05f);
        assertTrue(x[yMinIdx] == 1);
        
        //=========================================
                
        float[] gsd = Gaussian1DSecondDeriv.getKernel(sigma);
        
        x = new float[gsd.length];
        
        h = gsd.length >> 1;
        sum = 0;
        for (int i = 0; i < gfd.length; i++) {             
            x[i] = (i - h);
            sum += gsd[i];
        }
        
        yMax = MiscMath.findMax(gsd);
        yMin = MiscMath.findMin(gsd);
        xMin = MiscMath.findMin(x);
        xMax = MiscMath.findMax(x);
        
        plotter = 
            new PolygonAndPointPlotter(xMin - 1, xMax + 1, yMin, yMax);
        
        plotter.addPlot(x, gsd, x, gsd, 
            "gaussian 1-d second derivative for sigma=" + sigma);
        
        filePath = plotter.writeFile3();
        
        System.out.println(filePath);
                
        assertTrue(Math.abs(sum) < 0.1);
        
        assertTrue(x[h] == 0);
        assertTrue(gfd[h] == 0);
        assertTrue((int)sum == 0);
        yMaxIdx = MiscMath.findYMaxIndex(gfd);
        //assertTrue(Math.abs(gfd[yMaxIdx] - 0.16) < 0.05f);
        //assertTrue(x[yMaxIdx] == -2);*/
    }
    
    @Test
    public void testGetHalfKernelUsingBinomialFilter() throws IOException {
        
        System.out.println("testGetHalfKernelUsingBinomialFilter");
                
        double[] halfKernel = Gaussian1D.getHalfKernelUsingBinomialFilterSigmaOne();
        
        int n = halfKernel.length*2 - 1;
        int h = halfKernel.length - 1;
        int hm1 = halfKernel.length - 1;
        
        float[] x = new float[n];
        float[] g = new float[n];
        
        float sum = 0;
        for (int i = 0; i < halfKernel.length; i++) {
            
            double k = halfKernel[i]; // 0 1 2 3 4 5
            
            int x0 = i - h;
            x[i] = x0;
            g[i] = (float)k;
            sum += g[i];
            
            if (i < hm1) {
                int idx = n - i - 1;
                int x1 = (idx - h);
                x[idx] = x1;
                g[idx] = (float)k;
                sum += g[idx];
            }
        }
                
        float yMax = MiscMath.findMax(g);
        
        PolygonAndPointPlotter plotter = 
            new PolygonAndPointPlotter(x[0] - 1, x[x.length - 1] + 1, 0, yMax);
        
        plotter.addPlot(x, g, x, g, "gaussian 1-d for sigma="  + SIGMA.TWO + " sum=" 
            + sum);
        
        String filePath = plotter.writeFile();
    }
  
}
