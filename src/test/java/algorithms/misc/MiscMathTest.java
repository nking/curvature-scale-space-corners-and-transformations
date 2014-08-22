package algorithms.misc;

import algorithms.util.ArrayPair;
import java.security.SecureRandom;
import java.util.logging.Logger;

import algorithms.util.PolygonAndPointPlotter;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MiscMathTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    /**
     * Test of findPowerOf10 method, of class MiscMath.
     */
    public void testFindPowerOf10() {

        assertTrue(MiscMath.findPowerOf10(0.1f) == -1);

        assertTrue(MiscMath.findPowerOf10(0.11f) == -1);

        assertTrue(MiscMath.findPowerOf10(0.55f) == -1);

        int p0;
        int p1;
        p0 = MiscMath.findPowerOf10(0.099999999f);
        p1 = MiscMath.findPowerOf10_2(0.099999999f);
        //assertTrue(MiscMath.findPowerOf10(0.099999999f) == -2); <=== precision error before argument reaches method.

        p0 = MiscMath.findPowerOf10(0.01f);
        p1 = MiscMath.findPowerOf10_2(0.01f);
        assertTrue(p0 == -2);

        assertTrue(MiscMath.findPowerOf10(0.011f) == -2);

        assertTrue(MiscMath.findPowerOf10(0.001f) == -3);

        //assertTrue(MiscMath.findPowerOf10(0.0099999999f) == -3); <=== precision error before argument reaches method.

        assertTrue(MiscMath.findPowerOf10(0.f) == 0);

        assertTrue(MiscMath.findPowerOf10(1.f) == 0);

        assertTrue(MiscMath.findPowerOf10(10.f) == 1);

        assertTrue(MiscMath.findPowerOf10(11.f) == 1);

        assertTrue(MiscMath.findPowerOf10(100.f) == 2);

        assertTrue(MiscMath.findPowerOf10(-3.1f) == 0);

        assertTrue(MiscMath.findPowerOf10(-31.1f) == 1);

        assertTrue(MiscMath.findPowerOf10(-310.1f) == 2);

        assertTrue(MiscMath.findPowerOf10(-0.1f) == -1);
    }

    /**
     * Test of roundDownByLargestPower method, of class MiscMath.
     */
    public void testRoundDownByLargestPower() {

        assertTrue(MiscMath.roundDownByLargestPower(0.1f) == 0.1f);

        assertTrue(MiscMath.roundDownByLargestPower(0.11f) == 0.1f);

        assertTrue(MiscMath.roundDownByLargestPower(3.1f) == 3.0f);

        assertTrue(MiscMath.roundDownByLargestPower(31.1f) == 30.0f);

        assertTrue(MiscMath.roundDownByLargestPower(31.0f) == 30.0f);

        assertTrue(MiscMath.roundDownByLargestPower(-3.1f) == -4.0f);

        assertTrue(MiscMath.roundDownByLargestPower(-31.1f) == -40.0f);

        assertTrue(MiscMath.roundDownByLargestPower(-31.0f) == -40.0f);
    }

    /**
     * Test of roundUpByExponent method, of class MiscMath.
     */
    public void testRoundUpByLargestPower() {

        assertTrue(MiscMath.roundUpByLargestPower(31.1f) == 40.0f);

        assertTrue(MiscMath.roundUpByLargestPower(0.11f) == 0.2f);

        assertTrue(MiscMath.roundUpByLargestPower(-0.011f) == -0.02f);

        assertTrue(MiscMath.roundUpByLargestPower(-0.11f) == -0.2f);

        assertTrue(MiscMath.roundUpByLargestPower(310.1f) == 400.0f);

        assertTrue(MiscMath.roundUpByLargestPower(-3.1f) == -4.0f);

        assertTrue(MiscMath.roundUpByLargestPower(3.1f) == 4.0f);

        assertTrue(MiscMath.roundUpByLargestPower(10.0f) == 10.0f);

        assertTrue(MiscMath.roundUpByLargestPower(5.0f) == 5.0f);
    }
    
    public void testTaylor() throws Exception {
        double s;
        int i;
        double onePlusX = 1 + -0.001;
        i = 4;
        s = MiscMath.taylor(onePlusX, i);
        assertTrue(Math.abs(s) < 0.01);
        
        i = 11;
        s = MiscMath.taylor(onePlusX, i);
        assertTrue(Math.abs(s) < 0.01);
        
        i = 16;
        s = MiscMath.taylor(onePlusX, i);
        assertTrue(Math.abs(s) < 0.01);
        
        
        onePlusX = 1 + -1.001;
        i = 4;
        s = MiscMath.taylor(onePlusX, i);       // -1.836
        assertTrue(Math.abs(s - -2.09) < 0.1);
        
        i = 11;
        s = MiscMath.taylor(onePlusX, i);
        assertTrue(Math.abs(s - -3.03) < 0.1);
        
        i = 16;
        s = MiscMath.taylor(onePlusX, i);
        assertTrue(Math.abs(s - -3.397) < 0.1);
        
        
        s = MiscMath.taylor(2.00001, 4);
        s = MiscMath.taylor(2.00001, 11);
        s = MiscMath.taylor(2.00001, 16);
        
        onePlusX = 1 + -1.00001;
        i = 4;
        s = MiscMath.taylor(onePlusX, i);       // -1.836
        assertTrue(Math.abs(s - -2.083) < 0.1);
        
        i = 11;
        s = MiscMath.taylor(onePlusX, i);
        assertTrue(Math.abs(s - -3.02) < 0.1);
        
        i = 16;
        s = MiscMath.taylor(onePlusX, i);
        assertTrue(Math.abs(s - -3.38) < 0.1);
    }
    
    public void testPoissonRandom() throws Exception {
                
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        
        long seed = System.currentTimeMillis();
        //seed = 1386750505246l;

        sr.setSeed(seed);
        log.info("SEED=" + seed);
        
        int nIter = 100;
        
        int lambda = 2;
        
        float[] values = new float[nIter];
        float[] valueErrors = new float[nIter];
        
        for (int i = 0; i < nIter; i++) {
        
            int rInt = MiscMath.poissonRandom(sr, lambda);
        
            values[i] = rInt;
            valueErrors[i] = 0.03f*rInt;            
        }
        
        int nBins = (int) MiscMath.findMax(values);
        
        HistogramHolder hist = Histogram.createSimpleHistogram(nBins, values, valueErrors);
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(hist.getXHist(), hist.getYHistFloat(), null, null, "");
        plotter.writeFile();
    }
    
    public void test0() throws Exception {
        // constructor isn't used, this is to complete the coverage.
        MiscMath mm = new MiscMath();
    }

    /**
     * Test of findMax method, of class MiscMath.
     */
    public void testFindMax_floatArr() {

    }

    /**
     * Test of findMax method, of class MiscMath.
     */
    public void testFindMax_intArr() {

    }

    /**
     * Test of findMin method, of class MiscMath.
     */
    public void testFindMin() {

    }

    /**
     * Test of findYMinIndex method, of class MiscMath.
     */
    public void testFindYMinIndex() {

    }

    /**
     * Test of findYMaxIndex method, of class MiscMath.
     */
    public void testFindYMaxIndex_floatArr() {

    }

    /**
     * Test of findYMaxIndex method, of class MiscMath.
     */
    public void testFindYMaxIndex_floatArr_int() {

    }
}
