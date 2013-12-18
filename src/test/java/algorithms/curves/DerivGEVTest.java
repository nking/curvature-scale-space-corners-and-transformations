package algorithms.curves;

import java.security.SecureRandom;
import java.util.Arrays;
import java.util.logging.Logger;

import algorithms.misc.MiscMath;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.TestCase;

public class DerivGEVTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public void testDerivWRTX() throws Exception {
        
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
            log.fine( xp[i] + ":" + DerivGEV.derivWRTX(yConst, mu, k, sigma, xp[i]));
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
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed( seed);
        log.info("seed=" + seed);
        
        for (int i = 0; i < 100; i++) {
            float x = sr.nextFloat()*sr.nextInt(120000);
            Double deriv = DerivGEV.derivWRTX(yConst, mu, k, sigma, x);
            log.fine( x + ":" + deriv);
        }
        
        /*
        float s = 0.85f;
        s = 1.9f;
        
        k = 1.0f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, sigma, mu);                    
        plotter.addPlot(xp, yGEV, null, null, null, null, "k=1");
        filePath = plotter.writeFile();
        
        k = 1.5f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, sigma, mu);                    
        plotter.addPlot(xp, yGEV, null, null, null, null, "k=1.5");
        filePath = plotter.writeFile();
        
        k = 2.0f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, sigma, mu);                    
        plotter.addPlot(xp, yGEV, null, null, null, null, "k=2");
        filePath = plotter.writeFile();
        
        k = 2.5f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, sigma, mu);                    
        plotter.addPlot(xp, yGEV, null, null, null, null, "k=2.5");
        filePath = plotter.writeFile();
        
        k = 3.0f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, sigma, mu);                    
        plotter.addPlot(xp, yGEV, null, null, null, null, "k=3");
        filePath = plotter.writeFile();
        
        k = 3.5f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, sigma, mu);                    
        plotter.addPlot(xp, yGEV, null, null, null, null, "k=3.5");
        filePath = plotter.writeFile();
        
        k = 4.0f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, sigma, mu);                    
        plotter.addPlot(xp, yGEV, null, null, null, null, "k=4");
        filePath = plotter.writeFile();
        
        k = 4.5f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, sigma, mu);                    
        plotter.addPlot(xp, yGEV, null, null, null, null, "");
        filePath = plotter.writeFile();
        */
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
        
        /*
        // ====== exploring known GEV curve and expected residual suggested for a change in k ====
        float minChiSqSum = Float.MAX_VALUE;
        int minChiSqSumIdx = 0;
        float[] yg0 = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, sigma, mu);
        float[] yg0e = new float[yg0.length];
        Arrays.fill(yg0e, 0.03f);
        float s2  = 0.85f; //1.07f;
        float mu2 = 0.441f;//0.59f;
        float k2  = 1.80f - 1.7f; //2.75f;
        float avgResidK = 0.f;
        int yMaxModelIdx = MiscMath.findYMaxIndex(GeneralizedExtremeValue.generateNormalizedCurve(xp, k2, s2, mu2));
        log.fine("yMaxModelIdx=" + yMaxModelIdx);
        double[] suggestedK = new double[xp.length];
        for (int i = 0; i < xp.length; i++) {            
            Double d = DerivGEV.derivWRTK(yConst, mu2, k2, s2, xp[i]);
            double deltaK = 0.0001f;
            Double d2 = DerivGEV.derivWRTK(yConst, mu2, (float)(k2 - deltaK), s2, xp[i]);
            Double dd = (d2-d)/deltaK;            
            double preconditionedResidual = d2/dd;
            float chiSqSum = DerivGEV.chiSqSum((float)(k2 + preconditionedResidual), s2, mu2, xp, yg0, yg0e);
            float chiSqSum2 = DerivGEV.chiSqSum((float)(k2 - preconditionedResidual), s2, mu2, xp, yg0, yg0e);
            log.fine( String.format("x[%d]=%4.3f  (d/dk=%4.5f, d2/dkdk=%4.5f) ==> (+%4.4f  chiSqSum=%4.4f) (-%4.4f  chiSqSum=%4.4f)", 
                i, xp[i], d, dd, preconditionedResidual, chiSqSum, preconditionedResidual, chiSqSum2));
            avgResidK += preconditionedResidual;
            suggestedK[i] = (chiSqSum < chiSqSum2) ? preconditionedResidual : -1.f*preconditionedResidual;
            if (suggestedK[i] < minChiSqSum) {
                minChiSqSum = (float) suggestedK[i];
                minChiSqSumIdx = i;
            }
        }
        avgResidK /= xp.length;
        float stdDevSuggestedK = 0;
        for (int i = 0; i < xp.length; i++) {
            stdDevSuggestedK += Math.pow((suggestedK[i] - avgResidK), 2);
        }
        stdDevSuggestedK = (float) (Math.sqrt(stdDevSuggestedK/(xp.length - 1.0f)));//N-1 because had to calculate mean from the data
        float avgKWithoutOutliers = 0;
        int countWithoutOutliers = 0;
        for (int i = 0; i < xp.length; i++) {
            float diff = (float)suggestedK[i] - avgResidK;
            if (Math.abs(diff - avgResidK) < 3.*stdDevSuggestedK) {
                avgKWithoutOutliers += suggestedK[i];
                countWithoutOutliers++;
            } else {
                log.fine("  exclude " + suggestedK[i]);
            }
        }
        avgKWithoutOutliers /= countWithoutOutliers;
        log.fine("===> avg d/dk/d2/dkdk = " + avgKWithoutOutliers 
            + "  d/dk/d2/dkdk at minChiSqSum = " + suggestedK[minChiSqSumIdx]
            + "  d/dk/d2/dkdk at yMax=" + suggestedK[yMaxModelIdx]);
        // results show it's still the best approach to take the minChiSqSum solution from the min change in d/dk
        //    the problem with this approach is we are actually needing to know what deltaK to test chisqsum for
        //    and the curve gives us d/dk modified by d2/dkdk as the suggested step in either direction that
        //    will create a change in the gev.
        //    the correct answer would test both + modified d/dk and - d/dk
        //  is the idx at position +1 past model y max always the best answer for d/dk?
        // ===> looks like a good pattern would be to get the modified d/dk at position yMaxIdx + 1
        //      and test whether adding that or subtracting it improves the chi sq sum and take the best answer.
        */
        
        float factor = 0.0001f;
        
        GeneralizedExtremeValue gev = new GeneralizedExtremeValue(new float[0], 
            new float[0], new float[0], new float[0]);
        float[] yGEV0 = gev.generateCurve(xp, k - (k*factor), sigma, mu);
        float[] yGEV1 = gev.generateCurve(xp, k, sigma, mu);
        float[] yGEV2 = gev.generateCurve(xp, k + (k*factor), sigma, mu);
                    
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0.f, 1.0f, 0f, 1.3f);
     
        
        plotter.addPlot(xp, yGEV0, null, null, null, null, "k-deltaK");
        plotter.addPlot(xp, yGEV1, null, null, null, null, "k=" + k);
        plotter.addPlot(xp, yGEV2, null, null, null, null, "k+deltaK");
        String filePath = plotter.writeFile();
       
        double last = Integer.MIN_VALUE;
        for (int i = 1; i < 10; i++) {
            double delta = (k*factor) * i;
            float d = (float) (k + delta);
                        
            // looks like the deriv increases with increasing shape k for this region of the curve
            
            Double deriv = DerivGEV.derivWRTK(yConst, mu, d, sigma, xPoint);
            
            if (deriv != null) {
                log.fine("k=" + d + "  derivWRTK=" + deriv);
            }
            
            assertNotNull(deriv);
            
            // these are all nearly the same value since x is the same and delta k is small
            //assertTrue(last < deriv.doubleValue());
            
            Double deriv2 = DerivGEV.estimateDerivUsingDeltaK(mu, k, sigma, xPoint);
            
            log.fine("  compare derivWRTK to  estimateDerivUsingDeltaK " + deriv + " : " + deriv2);
            
            last = deriv.doubleValue();           
        }
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed( seed);
        log.info("seed=" + seed);
        
        for (int i = 0; i < 100000; i++) {
            float x = sr.nextFloat()*sr.nextInt(120000);
            while (x == 0) {
                x = sr.nextFloat()*sr.nextInt(120000);
            }
            Double deriv = DerivGEV.derivWRTK(yConst, mu, k, sigma, x);
            assertNotNull(deriv);
            Double deriv2 = DerivGEV.estimateDerivUsingDeltaK(mu, k, sigma, x);
            assertNotNull(deriv2);
                        
            float m = sr.nextFloat()*sr.nextInt((int)Math.ceil(x));
            float s = sr.nextFloat()*sr.nextInt(10);
            while (s == 0) {
                s = sr.nextFloat()*sr.nextInt(10);
            }
            float k_ = sr.nextFloat()*sr.nextInt(10);
            while (k_ == 0) {
                k_ = sr.nextFloat()*sr.nextInt(10);
            }
            deriv = DerivGEV.derivWRTK(yConst, m, k_, s, x);
            assertNotNull(deriv);

            //log.info( x + ":" + deriv);
        }
        
    }

    public void testDerivWRTSigma() throws Exception {
        
        log.info("testDerivWRTSigma");
        
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
        
        /*
        // ====== exploring known GEV curve and expected residual suggested for a change in k ====
        float minChiSqSum = Float.MAX_VALUE;
        int minChiSqSumIdx = 0;
        float[] yg0 = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, sigma, mu);
        float[] yg0e = new float[yg0.length];
        Arrays.fill(yg0e, 0.03f);
        float s2  = sigma - 0.4f;
        float mu2 = mu;
        float k2  = k;
        float avgResid = 0.f;
        int yMaxModelIdx = MiscMath.findYMaxIndex(GeneralizedExtremeValue.generateNormalizedCurve(xp, k2, s2, mu2));
        log.fine("yMaxModelIdx=" + yMaxModelIdx);
        double[] suggested = new double[xp.length];
        for (int i = 0; i < xp.length; i++) {            
            Double d = DerivGEV.derivWRTSigma(yConst, mu2, k2, s2, xp[i]);
            
            // modification by 2nd derivs for sigma is more complex
            // (∂^2f/∂sigma∂sigma)
            double delta = 0.0001f*s2;
            Double d2 = DerivGEV.derivWRTSigma(yConst, mu2, k2, (float)(s2 + delta), xp[i]);
            Double dd = (d2-d)/delta;
            
            double preconditionedResidual = DerivGEV.calculatePreconditionerModifiedResidualSigma(yConst, mu2, k2, sigma, xp[i]);
            
            float chiSqSum = DerivGEV.chiSqSum(k2, (float)(s2 + preconditionedResidual), mu2, xp, yg0, yg0e);
            float chiSqSum2 = DerivGEV.chiSqSum(k2, (float)(s2 - preconditionedResidual), mu2, xp, yg0, yg0e);
            log.info( String.format("x[%d]=%4.3f  (d/ds=%4.5f, d2/dsds=%4.5f) ==> (+%4.4f  chiSqSum=%4.4f) (-%4.4f  chiSqSum=%4.4f)", 
                i, xp[i], d, dd, preconditionedResidual, chiSqSum, preconditionedResidual, chiSqSum2));
            avgResid += preconditionedResidual;
            suggested[i] = (chiSqSum < chiSqSum2) ? preconditionedResidual : -1.f*preconditionedResidual;
            if (suggested[i] < minChiSqSum) {
                minChiSqSum = (float) suggested[i];
                minChiSqSumIdx = i;
            }
        }
        avgResid /= xp.length;
        float stdDevSuggested = 0;
        for (int i = 0; i < xp.length; i++) {
            stdDevSuggested += Math.pow((suggested[i] - avgResid), 2);
        }
        stdDevSuggested = (float) (Math.sqrt(stdDevSuggested/(xp.length - 1.0f)));//N-1 because had to calculate mean from the data
        float avgKWithoutOutliers = 0;
        int countWithoutOutliers = 0;
        for (int i = 0; i < xp.length; i++) {
            float diff = (float)suggested[i] - avgResid;
            if (Math.abs(diff - avgResid) < 3.*stdDevSuggested) {
                avgKWithoutOutliers += suggested[i];
                countWithoutOutliers++;
            } else {
                log.fine("  exclude " + suggested[i]);
            }
        }
        avgKWithoutOutliers /= countWithoutOutliers;
        log.fine("===> avg d/ds/d2/dsds = " + avgKWithoutOutliers 
            + "  d/ds/d2/dsds at minChiSqSum = " + suggested[minChiSqSumIdx]
            + "  d/ds/d2/dsds at yMax=" + suggested[yMaxModelIdx]);
        // ===> looks like a good pattern would be to get the modified d/dsigma at position yMaxIdx + 1
        //      and test whether adding that or subtracting it improves the chi sq sum and take the best answer.
        // it's quick and gives an answer in the right direction.
        */
        
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
                log.fine("sigma=" + d + "  derivWRTSigma=" + deriv);
            }
            
            assertNotNull(deriv);
            
            assertTrue(last < deriv.doubleValue());
            
            Double deriv2 = DerivGEV.estimateDerivUsingDeltaSigma(mu, k, sigma, xPoint);
            
            log.fine("  compare derivWRTSigma to  estimateDerivUsingDeltaSigma " + deriv + " : " + deriv2);
            
            last = deriv.doubleValue();           
        }
        log.fine("==> yFactor=" + yFactor + " delta=" + sigma*factor);
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed( seed);
        log.info("seed=" + seed);
        
        for (int i = 0; i < 100000; i++) {
            float x = sr.nextFloat()*sr.nextInt(120000);
            while (x == 0) {
                x = sr.nextFloat()*sr.nextInt(120000);
            }
            Double deriv = DerivGEV.derivWRTSigma(yConst, mu, k, sigma, x);
            assertNotNull(deriv);
            Double deriv2 = DerivGEV.estimateDerivUsingDeltaSigma(mu, k, sigma, x);
            assertNotNull(deriv2);
            //log.fine( x + ":" + deriv);
            
            float m = sr.nextFloat()*sr.nextInt((int)Math.ceil(x));
            float s = sr.nextFloat()*sr.nextInt(10);
            while (s == 0) {
                s = sr.nextFloat()*sr.nextInt(10);
            }
            float k_ = sr.nextFloat()*sr.nextInt(10);
            while (k_ == 0) {
                k_ = sr.nextFloat()*sr.nextInt(10);
            }
            deriv = DerivGEV.derivWRTSigma(yConst, m, k_, s, x);
            assertNotNull(deriv);
        }
    }

    public void testDerivWRTMu() throws Exception {
        
        log.info("testDerivWRTMu");
        
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
        /*
        // ====== exploring known GEV curve and expected residual suggested for a change in k ====
        float minChiSqSum = Float.MAX_VALUE;
        int minChiSqSumIdx = 0;
        float[] yg0 = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, sigma, mu);
        float[] yg0e = new float[yg0.length];
        Arrays.fill(yg0e, 0.03f);
        float s2  = sigma;
        float mu2 = mu  - 0.4f;
        float k2  = k;
        float avgResid = 0.f;
        int yMaxModelIdx = MiscMath.findYMaxIndex(GeneralizedExtremeValue.generateNormalizedCurve(xp, k2, s2, mu2));
        log.fine("yMaxModelIdx=" + yMaxModelIdx);
        double[] suggested = new double[xp.length];
        for (int i = 0; i < xp.length; i++) {            
            Double d = DerivGEV.derivWRTSigma(yConst, mu2, k2, s2, xp[i]);
            
            // modification by 2nd derivs for mu is more complex
            // (∂^2f/∂mu∂mu)
            double delta = 0.0001f*mu2;
            Double d2 = DerivGEV.derivWRTMu(yConst, (float)(mu2 + delta), k2, s2, xp[i]);
            Double dd = (d2-d)/delta;
            
            double preconditionedResidual = DerivGEV.calculatePreconditionerModifiedResidualMu(yConst, mu2, k2, sigma, xp[i]);
            
            float chiSqSum = DerivGEV.chiSqSum(k2, s2, (float)(mu2 + preconditionedResidual), xp, yg0, yg0e);
            float chiSqSum2 = DerivGEV.chiSqSum(k2, s2, (float)(mu2 - preconditionedResidual), xp, yg0, yg0e);
            log.fine( String.format("x[%d]=%4.3f  (d/dm=%4.5f, d2/dmdm=%4.5f) ==> (+%4.4f  chiSqSum=%4.4f) (-%4.4f  chiSqSum=%4.4f)", 
                i, xp[i], d, dd, preconditionedResidual, chiSqSum, preconditionedResidual, chiSqSum2));
            avgResid += preconditionedResidual;
            suggested[i] = (chiSqSum < chiSqSum2) ? preconditionedResidual : -1.f*preconditionedResidual;
            if (suggested[i] < minChiSqSum) {
                minChiSqSum = (float) suggested[i];
                minChiSqSumIdx = i;
            }
        }
        avgResid /= xp.length;
        float stdDevSuggested = 0;
        for (int i = 0; i < xp.length; i++) {
            stdDevSuggested += Math.pow((suggested[i] - avgResid), 2);
        }
        stdDevSuggested = (float) (Math.sqrt(stdDevSuggested/(xp.length - 1.0f)));//N-1 because had to calculate mean from the data
        float avgKWithoutOutliers = 0;
        int countWithoutOutliers = 0;
        for (int i = 0; i < xp.length; i++) {
            float diff = (float)suggested[i] - avgResid;
            if (Math.abs(diff - avgResid) < 3.*stdDevSuggested) {
                avgKWithoutOutliers += suggested[i];
                countWithoutOutliers++;
            } else {
                log.fine("  exclude " + suggested[i]);
            }
        }
        avgKWithoutOutliers /= countWithoutOutliers;
        log.fine("===> avg d/dm/d2/dmdm = " + avgKWithoutOutliers 
            + "  d/dm/d2/dmdm at minChiSqSum = " + suggested[minChiSqSumIdx]
            + "  d/dm/d2/dmdm at yMax=" + suggested[yMaxModelIdx]);
        // ===> looks like a good pattern would be to get the modified d/dsigma at position yMaxIdx + 1
        //      and test whether adding that or subtracting it improves the chi sq sum and take the best answer.
        // it's quick and gives an answer in the right direction.
       */
        float factor = 0.0001f;
        
        GeneralizedExtremeValue gev = new GeneralizedExtremeValue(new float[0], 
            new float[0], new float[0], new float[0]);
                    
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0.f, 1.0f, 0f, 1.3f);        
        
        double last = Integer.MIN_VALUE;
        for (int i = 1; i < 10; i++) {
            double delta = (mu*factor) * i;
            float d = (float) (mu + delta);
                                    
            Double deriv = DerivGEV.derivWRTMu(yConst, d, k, sigma, xPoint);
            if (deriv != null) {
                log.info("mu=" + d + "  derivWRTMu=" + deriv);
            }
            float[] yGEV1 = gev.generateCurve(xp, k, sigma, d);
            plotter.addPlot(xp, yGEV1, null, null, null, null, "mu=" + d);
            
            assertNotNull(deriv);
            
            //assertTrue(last < deriv.doubleValue());
            

            Double deriv2 = DerivGEV.estimateDerivUsingDeltaMu(mu, k, sigma, xPoint);
            
            log.fine("  compare derivWRTMu to  estimateDerivUsingDeltaMu " + deriv + " : " + deriv2);
            
            last = deriv.doubleValue();       
        }
        String filePath = plotter.writeFile();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed( seed);
        log.info("seed=" + seed);
        
        for (int i = 0; i < 100000; i++) {
            float x = sr.nextFloat()*sr.nextInt(120000);
            while (x == 0) {
                x = sr.nextFloat()*sr.nextInt(120000);
            }
            Double deriv = DerivGEV.derivWRTMu(yConst, mu, k, sigma, x);
            assertNotNull(deriv);
            Double deriv2 = DerivGEV.estimateDerivUsingDeltaMu(mu, k, sigma, x);
            assertNotNull(deriv2);
            //log.fine( x + ":" + deriv);
            
            float m = sr.nextFloat()*sr.nextInt((int)Math.ceil(x));
            float s = sr.nextFloat()*sr.nextInt(10);
            while (s == 0) {
                s = sr.nextFloat()*sr.nextInt(10);
            }
            float k_ = sr.nextFloat()*sr.nextInt(10);
            while (k_ == 0) {
                k_ = sr.nextFloat()*sr.nextInt(10);
            }
            deriv = DerivGEV.derivWRTMu(yConst, m, k_, s, x);
            assertNotNull(deriv);
        }
    }

}
