package algorithms.imageProcessing;

import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 * @author nichole
 */
public class CurvatureTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public CurvatureTest(String testName) {
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
    
    public void testComputeCurvature_circle() throws Exception {
                
        // test that the curvature for every point on a circle is ~1/r

        /*
        create test circles: 
           -- with enough points so that the largest sigma kernels still 
              completely cover them.
           -- of a radius and intervals that result in a measurable size
              near or above the FWHM of a gaussian of size sigma
              (presumably, features smaller than that are diluted by the
              convolution, so are not as good for testing for the expected
              values of curvature).
        
        Gradient constructed from kernel of size sigma should be sensitive to
        scale sizes roughly larger than the FWHM of a gaussian with given sigma.
        
        half width of fraction max 
            = sigma * Math.sqrt(-2. * Math.log(fractionMax))
        HWHM = sigma*1.177 
        FWHM = 2.3548 * sigma
        
        So features that vary over scales < 2.3548 * sigma (diluted by the
            gaussian convolution) are less easy to detect.
        
        To test curvature for a circle:  we measure 1/radius.
        For each point in the circle, the ds along the circle must be measurable
        for the sigma.
        
        math.cos theta ranging from 0 to 1
           r * [0 - 1] >= 2.3548 * sigma
        
           if testing up to sigma=4, then a ds of 10 is good.
        
           ds = 10;  
           2 * pi * r = nPoints * ds
           r = (nPoints * ds)/(2 * pi)
        
        */
        
        double dTheta = 10.0;
        int n = (int)(360.f/dTheta);
        float xc0 = 20.0f;
        float yc0 = 20.0f;
        float xc1 = 400.0f;
        float yc1 = 400.0f;
        float r = (float)(((float)n) * 10.f/(2. * Math.PI));//10.0f;
        
        float expectedCurvature = (1.f/r);
        
        log.info("r=" + r + " 1/r=" + expectedCurvature);
        
        PairIntArray xy0 = new PairIntArray(n);
        PairIntArray xy1 = new PairIntArray(n);
        
        // need each pixel to be changing, so a dTheta of right size is needed
        // so dTheta >= 10 is a good choice
        for (int i = 0; i < n; i++) {
            /*
            (x-xc)^2 + (y-yc)^2 = r
            x = xc + r*cos(theta)
            y = yc + r*sin(theta)
            */
            double thetaRadians = (i*dTheta)/180.;
            
            int x0 = (int) (xc0 + r * Math.cos(thetaRadians));
            int y0 = (int) (yc0 + r * Math.sin(thetaRadians));
            xy0.add(x0, y0);
                        
            int x1 = (int) (xc1 + r * Math.cos(thetaRadians));
            int y1 = (int) (yc1 + r * Math.sin(thetaRadians));
            xy1.add(x1, y1);
        }
                                
        ScaleSpaceCurvature scaleSpaceHelper = new ScaleSpaceCurvature();
        
        for (SIGMA sigma : SIGMA.values()) {
            
            if (sigma.ordinal() == SIGMA.FOURSQRT2.ordinal()) {
                break;
            }
            
            /* uncomment to print only the binomial kernels
            if (null == Gaussian1DFirstDeriv.getBinomialKernel(sigma)) {
                continue;
            }*/        
            
            log.info("sigma=" + sigma.toString() + ")");
        
            ScaleSpaceCurve curve0 = scaleSpaceHelper.computeCurvature(xy0, 
                sigma, SIGMA.getValue(sigma));

            ScaleSpaceCurve curve1 = scaleSpaceHelper.computeCurvature(xy1,
                sigma, SIGMA.getValue(sigma));

            /*
            expecting curvature is 1/r ~ 0.1

            binomial kernel was important in getting the correct answer.
            */
            log.info("CIRCLE RESULTS: (sigma=" + sigma.toString() + ")");
            int h = n >> 1;
            for (int i = 0; i < n; i++) {

                log.info(i + ") " + curve0.getK(i) + " : " + curve1.getK(i));

                if (i > h) {
                    if (i < (n - h - 1)) {
                        assertTrue(Math.abs(curve0.getK(i) - expectedCurvature) 
                            < (0.75*expectedCurvature));
                        assertTrue(Math.abs(curve1.getK(i) - expectedCurvature) 
                            < (0.75*expectedCurvature));
                    }
                }
            }
        }
        
        log.info("ComputeCurvature_line");
       
        // test that the curvature for every point on a line is ~0
        
        // 1 horizontal line, vertical line, and diagonal line
        //n = 18;
        PairIntArray xyH = new PairIntArray(n);
        PairIntArray xyV = new PairIntArray(n);
        PairIntArray xyD = new PairIntArray(n);
        
        for (int i = 0; i < n; i++) {
            
            xyH.add(i, 100);
            
            xyV.add(10, i);
            
            xyD.add(i, i);
        }
                                        
        for (SIGMA sigma : SIGMA.values()) {
            
            if (sigma.ordinal() == SIGMA.FOUR.ordinal()) {
                break;
            }
            
            /*if (null == Gaussian1DFirstDeriv.getBinomialKernel(sigma)) {
                continue;
            }*/   
            
            log.info("sigma=" + sigma.toString() + ")");
        
            // x first deriv =  -1, x second deriv = 0, 
            // y first deriv = 0, y second deriv = 0,
            // k = 0
            ScaleSpaceCurve curveH = scaleSpaceHelper.computeCurvature(xyH, 
                sigma, SIGMA.getValue(sigma));

            // x first deriv =  0, x second deriv = 0, 
            // y first deriv = 1, y second deriv = 0,
            // k = 0
            ScaleSpaceCurve curveV = scaleSpaceHelper.computeCurvature(xyV, 
                sigma, SIGMA.getValue(sigma));

            // x first deriv =  1, x second deriv = 0, 
            // y first deriv = -1, y second deriv = 0,
            // k = 0.1
            ScaleSpaceCurve curveD = scaleSpaceHelper.computeCurvature(xyD, 
                sigma, SIGMA.getValue(sigma));

            log.info("LINE RESULTS: (sigma=" + sigma.toString() 
                + ")");
            
            int h = n >> 1;
            
            for (int i = 0; i < n; i++) {
                log.info(i + ") " + curveH.getK(i) + " : " + curveV.getK(i) 
                    + " : " + curveD.getK(i));
                if (i > h) {
                    if (i < (n - h - 1)) {
                        assertTrue(Math.abs(curveH.getK(i)) < 0.001);
                        assertTrue(Math.abs(curveV.getK(i)) < 0.001);
                        assertTrue(Math.abs(curveD.getK(i)) < 0.001);
                    }
                }
            }
        }
    }
    
    public static void main(String[] args) {
        
        try {
            CurvatureTest test = new CurvatureTest("CurvatureTest");
            test.testComputeCurvature_circle();
        } catch(Exception e) {
            System.out.println(e.getMessage());
        }
    }
}
