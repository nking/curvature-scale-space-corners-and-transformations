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
        
        for (int i = 0; i < xp.length; i++) {
            System.out.println( xp[i] + ":" + DerivGEV.derivWRTX(yConst, mu, k, sigma, xp[i]));
        }
        
        double d = DerivGEV.derivWRTX(yConst, mu, k, sigma, xp[0]);
        assertTrue( d > 0);
        
        d = DerivGEV.derivWRTX(yConst, mu, k, sigma, (float)(0.0428));
        assertTrue( Math.abs(d) < 0.1);
        
        d = DerivGEV.derivWRTX(yConst, mu, k, sigma, 0.7f);
        assertTrue( d < 0);
    }
    
    public void testDerivWRTK() throws Exception {
        
        // k is the shape parameter
        
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
        
        for (int i = 0; i < xp.length; i++) {
            System.out.println( xp[i] + ":" + DerivGEV.derivWRTK(yConst, mu, k, sigma, xp[i]));
        }
       
        double d = DerivGEV.derivWRTK(yConst, mu, k, sigma, xp[0]);
        assertTrue( d > 0);
       /* 
        d = DerivGEV.derivWRTK(yConst, mu, k, sigma, (float)(0.0428));
        assertTrue( Math.abs(d) < 0.1);
        
        d = DerivGEV.derivWRTK(yConst, mu, k, sigma, 0.7f);
        assertTrue( d < 0);*/
    }

    public void estDerivWRTSigma() throws Exception {
        
        float k = 0.000198f;    // 1E-5 to 1E-3
        float sigma = 0.115f;   // 0.025 to 0.5
        float mu = 0.17647058f; // normalized to 0:1 scale

        float[] xp = new float[20];
        for (int i = 0; i < xp.length; i++) {
            xp[i] = i;
        }
        GeneralizedExtremeValue gev = new GeneralizedExtremeValue(new float[0], 
            new float[0], new float[0], new float[0]);
        float[] yGEV = gev.generateCurve(xp, k, sigma, mu);
            
        float yConst = 1;
        
        double d = DerivGEV.derivWRTSigma(yConst, mu, k, sigma, xp[1]);
        
        assertTrue(d > 0);
    }

    public void estDerivWRTMu() throws Exception {
        
        float k = 0.000198f;    // 1E-5 to 1E-3
        float sigma = 0.115f;   // 0.025 to 0.5
        float mu = 0.17647058f; // normalized to 0:1 scale

        float[] xp = new float[20];
        for (int i = 0; i < xp.length; i++) {
            xp[i] = i;
        }
        GeneralizedExtremeValue gev = new GeneralizedExtremeValue(new float[0], 
            new float[0], new float[0], new float[0]);
        float[] yGEV = gev.generateCurve(xp, k, sigma, mu);
            
        float yConst = 1;
        
        double d = DerivGEV.derivWRTMu(yConst, mu, k, sigma, xp[1]);
        
        assertTrue(d > 0);
    }

}
