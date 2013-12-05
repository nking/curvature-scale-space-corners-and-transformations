package algorithms.curves;

import algorithms.misc.MiscMath;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.TestCase;

public class DerivGEVTest extends TestCase {

    public void estDerivWRTX() throws Exception {
        
        float k = 1.80f;
        float sigma = 0.85f;
        float mu = 0.441f;
        float yConst = 1;
        
        float[] xp = new float[30];

        for (int i = 0; i < xp.length; i++) {
            xp[i] = (float)i/xp.length;
        }
       
        GeneralizedExtremeValue gev = new GeneralizedExtremeValue(new float[0], 
            new float[0], new float[0], new float[0]);
        float[] yGEV = gev.generateCurve(xp, k, sigma, mu);
                    
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0.f, 1.0f, 0f, 1.3f);
        plotter.addPlot(xp, yGEV, null, null, null, null, "");
        String filePath = plotter.writeFile();
       
        
        // generate a curve, and look at several points to see where
        //   the slope is > 0,        == 0 ,    and < 0
        //                x=0.02    x=0.0333    x=0.7
        
        /*
        for (int i = 0; i < xp.length; i++) {
            System.out.println( xp[i] + ":" + DerivGEV.derivWRTX(yConst, mu, k, sigma, xp[i]));
        }
        */
        
        double d = DerivGEV.derivWRTX(yConst, mu, k, sigma, xp[0]);
        assertTrue( d > 0);
        
        // top of curve where slope should be near zero
        d = DerivGEV.derivWRTX(yConst, mu, k, sigma, (float)(0.04278)); // d=-0.000589
        assertTrue( Math.abs(d) < 0.01);
        
        // slope should be near -1 at x near 0.4
        d = DerivGEV.derivWRTX(yConst, mu, k, sigma, 0.4f);
        assertTrue( Math.abs(d - -1) < 0.1);
        
        d = DerivGEV.derivWRTX(yConst, mu, k, sigma, 0.7f);
        assertTrue( d < 0);
    }
    
    public void testDerivWRTK() throws Exception {
        
        // k is the shape parameter
        
        float k = 1.80f;
        float sigma = 0.85f;
        float mu = 0.441f;
        float yConst = 1;
        float xPoint = 0.4f; // point at which slope is approx -1 for dy/dx
        
        float[] xp = new float[30];

        for (int i = 0; i < xp.length; i++) {
            xp[i] = (float)i/xp.length;
        }
       
        float factor = 0.0001f;
        
        GeneralizedExtremeValue gev = new GeneralizedExtremeValue(new float[0], 
            new float[0], new float[0], new float[0]);
        float[] yGEV0 = gev.generateCurve(xp, k - (k*factor), sigma, mu);
        float[] yGEV1 = gev.generateCurve(xp, k, sigma, mu);
        float[] yGEV2 = gev.generateCurve(xp, k + (k*factor), sigma, mu);
                    
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0.f, 1.0f, 0f, 1.3f);
        
        for (int i = 0; i < yGEV1.length; i++) {
            Double deriv = DerivGEV.derivWRTK(yConst, mu, k, sigma, xp[i]);
            if (deriv != null) {
                System.out.println( String.format(" df(k)/dk= %8.4f  x=%.3f", deriv.doubleValue(), xp[i]));
            }
        }
        
        
        plotter.addPlot(xp, yGEV0, null, null, null, null, "k-deltaK");
        plotter.addPlot(xp, yGEV1, null, null, null, null, "k=" + k);
        plotter.addPlot(xp, yGEV2, null, null, null, null, "k+deltaK");
        String filePath = plotter.writeFile();
       
        double last = Integer.MIN_VALUE;
        for (int i = 1; i < 10; i++) {
            double delta = (k/10.) * i;
            float d = (float) (k + delta);
                        
            // looks like the deriv increases with increasing shape k for this region of the curve
            
            Double deriv = DerivGEV.derivWRTK(yConst, mu, d, sigma, xPoint);
            
            if (deriv != null) {
                System.out.println("k=" + d + "  derivWRTK=" + deriv);
            }
            
            assertNotNull(deriv);
            
            assertTrue(last < deriv.doubleValue());
            
            Double deriv2 = DerivGEV.estimateDerivUsingDeltaK(mu, k, sigma, xPoint);
            
            System.out.println("  compare derivWRTK to  estimateDerivUsingDeltaK " + deriv + " : " + deriv2);
            
            last = deriv.doubleValue();           
        }
        
        
        //TODO:  test for other curves and portion
        
    }

    public void estDerivWRTSigma() throws Exception {
        
        System.out.println("testDerivWRTSigma");
        
        // sigma is the scale factor
        
        float k = 1.80f;
        float sigma = 0.85f;
        float mu = 0.441f;
        float yConst = 1;
        float xPoint = 0.4f; // point at which slope is approx -1 for dy/dx
        
        float[] xp = new float[30];

        for (int i = 0; i < xp.length; i++) {
            xp[i] = (float)i/xp.length;
        }
        
        float factor = 0.0001f;
        
        float[] yGEV0 = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, (sigma - (sigma*factor)), mu);
        float[] yGEV1 = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, sigma, mu);
        float[] yGEV2 = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, (sigma + (sigma*factor)), mu);
                    
        float[] yGEV = GeneralizedExtremeValue.genCurve(xp, k, sigma, mu);
        int idx = MiscMath.findYMaxIndex(yGEV);
        idx = 12;
        float yFactor = yGEV[idx];
        float xForYFactor = xp[idx];
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0.f, 1.0f, 0f, 1.3f);
        plotter.addPlot(xp, yGEV0, null, null, null, null, "sigma-delta");
        plotter.addPlot(xp, yGEV1, null, null, null, null, "sigma=" + sigma);
        plotter.addPlot(xp, yGEV2, null, null, null, null, "sigma+delta");
        String filePath = plotter.writeFile();
       
        double last = Integer.MIN_VALUE;
        for (int i = 1; i < 10; i++) {
            double delta = (sigma*factor) * i;
            float d = (float) (sigma + delta);
                        
            // looks like the deriv increases with increasing scale sigma for this region of the curve
            
            Double deriv = DerivGEV.derivWRTSigma(yConst, mu, k, d, xPoint);
            if (deriv != null) {
                System.out.println("sigma=" + d + "  derivWRTSigma=" + deriv);
            }
            
            assertNotNull(deriv);
            
            assertTrue(last < deriv.doubleValue());
            
            Double deriv2 = DerivGEV.estimateDerivUsingDeltaSigma(mu, k, sigma, xPoint);
            
            System.out.println("  compare derivWRTSigma to  estimateDerivUsingDeltaSigma " + deriv + " : " + deriv2);
            
            last = deriv.doubleValue();           
        }
        System.out.println("==> yFactor=" + yFactor + " delta=" + sigma*factor);
        
    }

    public void estDerivWRTMu() throws Exception {
        
        System.out.println("testDerivWRTMu");
        
        // mu is the location variable
        
        float k = 1.80f;
        float sigma = 0.85f;
        float mu = 0.441f;
        float yConst = 1;
        float xPoint = 0.4f; // point at which slope is approx -1 for dy/dx
        
        float[] xp = new float[30];

        for (int i = 0; i < xp.length; i++) {
            xp[i] = (float)i/xp.length;
        }
       
        GeneralizedExtremeValue gev = new GeneralizedExtremeValue(new float[0], 
            new float[0], new float[0], new float[0]);
                    
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0.f, 1.0f, 0f, 1.3f);
               
        float factor = 0.0001f;
        
        double last = Integer.MIN_VALUE;
        for (int i = 1; i < 10; i++) {
            double delta = (mu*factor) * i;
            float d = (float) (mu + delta);
                                    
            Double deriv = DerivGEV.derivWRTMu(yConst, d, k, sigma, 0.7f);
            if (deriv != null) {
                System.out.println("mu=" + d + "  derivWRTMu=" + deriv);
            }
            float[] yGEV1 = gev.generateCurve(xp, k, sigma, d);
            plotter.addPlot(xp, yGEV1, null, null, null, null, "mu=" + d);
            
            assertNotNull(deriv);
            
            assertTrue(last < deriv.doubleValue());
            

            Double deriv2 = DerivGEV.estimateDerivUsingDeltaMu(mu, k, sigma, xPoint);
            
            System.out.println("  compare derivWRTMu to  estimateDerivUsingDeltaMu " + deriv + " : " + deriv2);
            
            last = deriv.doubleValue();       
        }
        String filePath = plotter.writeFile();
    }

}
