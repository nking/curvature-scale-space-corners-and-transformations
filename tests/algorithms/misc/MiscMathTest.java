package algorithms.misc;

import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.security.SecureRandom;
import java.util.logging.Logger;
import algorithms.util.PolygonAndPointPlotter;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
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
    
    public void testGetAvgAndStDev() throws Exception {
        
        int[] x = new int[]{2, 2, 2, 4, 4, 4};
        
        double[] avgAndStDev = MiscMath.getAvgAndStDev(x, x.length);
        assertEquals(avgAndStDev[0], 3.);
        assertTrue(Math.abs(avgAndStDev[1] - 1.1) < 0.01);
        
        avgAndStDev = MiscMath.getAvgAndStDev(x);
        assertEquals(avgAndStDev[0], 3.);
        assertTrue(Math.abs(avgAndStDev[1] - 1.1) < 0.01);
    }
    
    public void testGet20NeighborOffsets() throws Exception {
                
        PairIntArray offsets = MiscMath.get20NeighborOffsets();
        
        assertTrue(offsets.getN() == 20);
        
        Set<PairInt> expected = new HashSet<PairInt>();
        for (int dx = -2; dx <= 2; ++dx) {
            for (int dy = -2; dy <= 2; ++dy) {
                if (dx == 0 && dy == 0) {
                    continue;
                }
                expected.add(new PairInt(dx, dy));
            }
        }
        //remove the 4 corners        
        expected.remove(new PairInt(-2, -2));
        expected.remove(new PairInt(-2, 2));
        expected.remove(new PairInt(2, -2));
        expected.remove(new PairInt(2, 2));
        
        assertTrue(expected.size() == 20);
        
        for (int i = 0; i < offsets.getN(); ++i) {
            PairInt p = new PairInt(offsets.getX(i), offsets.getY(i));
            assertTrue(expected.remove(p));
        }
        
        assertTrue(expected.isEmpty());
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
        plotter.addPlot(hist.getXHist(), hist.getYHistFloat(),
            Errors.populateYErrorsBySqrt(hist.getXHist()),
            Errors.populateYErrorsBySqrt(hist.getYHistFloat()), "");
        plotter.writeFile();
    }

    public void test0() throws Exception {
        // constructor isn't used, this is to complete the coverage.
        MiscMath mm = new MiscMath();
    }

    public void testComputeNDivNMinusK() throws Exception {

        // n! / (n - k)!
        long result = MiscMath.computeNDivNMinusK(2, 1);
        assertTrue(result == 2);

        result = MiscMath.computeNDivNMinusK(4, 3);
        assertTrue(result == 4*3*2);

        result = MiscMath.computeNDivNMinusK(4, 2);
        assertTrue(result == 4*3);

        result = MiscMath.computeNDivNMinusK(6, 3);
        assertTrue(result == 6*5*4);
    }

    public void testFactorial() throws Exception {

        long result = MiscMath.factorial(2);
        assertTrue(result == 2);

        result = MiscMath.factorial(4);
        assertTrue(result == 4*3*2);

        result = MiscMath.factorial(5);
        assertTrue(result == 5*4*3*2);

        result = MiscMath.factorial(1);
        assertTrue(result == 1);

        result = MiscMath.factorial(0);
        assertTrue(result == 0);
    }
    
    public void testFactorialBigInteger() throws Exception {

        BigInteger expected = new BigInteger("87178291200");
        
        BigInteger result = MiscMath.factorialBigInteger(14);

        assertTrue(result.compareTo(expected) == 0);
    }
    
    public void testcomputeNDivKTimesNMinusK() throws Exception {
        
        int n, k;
        long nComb;
        
        n = 6;
        k = 2;
        nComb = MiscMath.computeNDivKTimesNMinusK(n, k);
        assertEquals(15, nComb);
        
        n = 2;
        k = 1;
        nComb = MiscMath.computeNDivKTimesNMinusK(n, k);
        assertEquals(2, nComb);
        
        n = (int)Math.sqrt(Integer.MAX_VALUE);
        k = 2;
        nComb = MiscMath.computeNDivKTimesNMinusK(n, k);
        assertEquals(((n*(n-1))/2), nComb);
        
        n = (int)Math.sqrt(Integer.MAX_VALUE);
        k = 2;
        nComb = MiscMath.computeNDivKTimesNMinusKBigInteger(n, k);
        assertEquals(((n*(n-1))/2), nComb);
        
    }
    
    /*
    protected static long computeNDivKTimesNMinusKBigInteger(int n, int k) {
    */

    public void testReadSetBits() throws Exception {

        List<Integer> setBits = new ArrayList<Integer>();

        Long bitstring = Long.valueOf(3);

        MiscMath.readSetBits(bitstring, setBits);

        assertTrue(setBits.size() == 2);
        assertTrue(setBits.get(0).intValue() == 0);
        assertTrue(setBits.get(1).intValue() == 1);

        //5	       101        ==> '0' '2'
        setBits.clear();
        bitstring = Long.valueOf(5);

        MiscMath.readSetBits(bitstring, setBits);

        assertTrue(setBits.size() == 2);
        assertTrue(setBits.get(0).intValue() == 0);
        assertTrue(setBits.get(1).intValue() == 2);

        //6	       110        ==> '1' '2'
        setBits.clear();
        bitstring = Long.valueOf(6);

        MiscMath.readSetBits(bitstring, setBits);

        assertTrue(setBits.size() == 2);
        assertTrue(setBits.get(0).intValue() == 1);
        assertTrue(setBits.get(1).intValue() == 2);

        //15	      1111        ==> '0' '1' '2' '3'
        setBits.clear();
        bitstring = Long.valueOf(15);

        MiscMath.readSetBits(bitstring, setBits);

        assertTrue(setBits.size() == 4);
        assertTrue(setBits.get(0).intValue() == 0);
        assertTrue(setBits.get(1).intValue() == 1);
        assertTrue(setBits.get(2).intValue() == 2);
        assertTrue(setBits.get(3).intValue() == 3);
    }

    public void testCubicRoots() throws Exception {

        double a0 = 1;
        double a1 = 6;
        double a2 = -4;
        double a3 = -24;

        double[] roots = MiscMath.solveCubicRoots(a0, a1, a2, a3);
        assertTrue(roots.length == 3);
        assertTrue(Math.abs(roots[0] - 2) < 0.1);
        assertTrue(Math.abs(roots[1] - -6) < 0.1);
        assertTrue(Math.abs(roots[2] - -2) < 0.1);

        a0 = 2;
        a1 = 3;
        a2 = -11;
        a3 = -6;
        roots = MiscMath.solveCubicRoots(a0, a1, a2, a3);
        assertTrue(roots.length == 3);
        assertTrue(Math.abs(roots[0] - 2) < 0.1);
        assertTrue(Math.abs(roots[1] - -3) < 0.1);
        assertTrue(Math.abs(roots[2] - -0.5) < 0.1);

        a0 = 1;
        a1 = -7;
        a2 = 4;
        a3 = 12;
        roots = MiscMath.solveCubicRoots(a0, a1, a2, a3);
        assertTrue(roots.length == 3);
        assertTrue(Math.abs(roots[0] - 6) < 0.1);
        assertTrue(Math.abs(roots[1] - -1) < 0.1);
        assertTrue(Math.abs(roots[2] - 2) < 0.1);
    }

    public void testFindStrongestPeakIndexes() throws Exception {

        float[] xHist = new float[]{
            -279.96667f, -239.90002f, -199.83334f, -159.76666f,
            -119.700005f, -79.63334f, -39.566673f, 0.5000076f,
            40.566658f, 80.63331f, 120.69999f, 160.76666f,
            200.83331f, 240.89996f, 280.96667f};

        int[] yHist = new int[] {0, 24, 141, 259, 482, 894, 1110, 1239,
            890, 608, 302, 151, 56, 0, 0};

        HistogramHolder h = new HistogramHolder();
        h.setYHist(yHist);
        h.setXHist(xHist);

        List<Integer> peakIndexes = MiscMath.findStrongPeakIndexes(h, 0.1f);

        assertNotNull(peakIndexes);
        assertTrue(peakIndexes.size() == 1);
        assertTrue(peakIndexes.get(0).intValue() == 7);

        List<Integer> sortedPeakIndexes =
            MiscMath.findStrongPeakIndexesDescSort(h, 0.09f);
        assertNotNull(sortedPeakIndexes);
        assertTrue(sortedPeakIndexes.size() == 1);
        assertTrue(sortedPeakIndexes.get(0).intValue() == 7);

        // this should have 3 peaks
        xHist = new float[]{
            -466.63333f, -399.9f, -333.16666f, -266.43335f, -199.70001f,
            -132.96667f, -66.23337f, 0.49996567f, 67.23331f, 133.96664f,
            200.69998f, 267.43332f, 334.1666f, 400.89993f, 467.63327f
        };

        yHist = new int[] {
            115, 265, 385, 560, 484, 581, 748, 809,
            661, 375, 371, 379, 304, 109, 10
        };

        h = new HistogramHolder();
        h.setYHist(yHist);
        h.setXHist(xHist);

        peakIndexes = MiscMath.findStrongPeakIndexes(h, 0.09f);

        assertNotNull(peakIndexes);
        assertTrue(peakIndexes.size() == 2);
        assertTrue(peakIndexes.get(0).intValue() == 3);
        assertTrue(peakIndexes.get(1).intValue() == 7);

        sortedPeakIndexes = MiscMath.findStrongPeakIndexesDescSort(h, 0.09f);

        assertNotNull(sortedPeakIndexes);
        assertTrue(sortedPeakIndexes.size() == 2);
        assertTrue(sortedPeakIndexes.get(0).intValue() == 7);
        assertTrue(sortedPeakIndexes.get(1).intValue() == 3);
    }

    public void testChooseRandomly() throws Exception {

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");

        long seed = System.currentTimeMillis();
        //seed = 1386750505246l;

        sr.setSeed(seed);
        log.info("SEED=" + seed);

        int k = 7;

        int[] selected = new int[k];

        //----
        int nMax = 7;

        MiscMath.chooseRandomly(sr, selected, nMax);

        Set<Integer> set = new HashSet<Integer>();
        for (int i = 0; i < selected.length; i++) {
            int s = selected[i];
            Integer sel = Integer.valueOf(s);
            assertFalse(set.contains(sel));
            set.add(sel);
        }

        //------
        nMax = 2*k;
        MiscMath.chooseRandomly(sr, selected, nMax);

        set = new HashSet<Integer>();
        for (int i = 0; i < selected.length; i++) {
            int s = selected[i];
            Integer sel = Integer.valueOf(s);
            assertFalse(set.contains(sel));
            set.add(sel);
        }

        //------
        nMax = sr.nextInt(500);
        while (nMax < 3*k) {
            nMax = sr.nextInt(500);
        }

        MiscMath.chooseRandomly(sr, selected, nMax);

        set = new HashSet<Integer>();
        for (int i = 0; i < selected.length; i++) {
            int s = selected[i];
            Integer sel = Integer.valueOf(s);
            assertFalse(set.contains(sel));
            set.add(sel);
        }
    }

    public void testNumberOfBits() throws Exception {

        for (int i = -256; i >= 0; i--) {
            String bitstring = Integer.toBinaryString(i);
            int expected = bitstring.length();
            int result = MiscMath.numberOfBits(i);
            assertTrue(expected == result);
        }

        for (int i = 0; i <= 256; i++) {
            String bitstring = Integer.toBinaryString(i);
            int expected = bitstring.length();
            int result = MiscMath.numberOfBits(i);
            assertTrue(expected == result);
        }
    }

    public void testBitReverse() throws Exception {

        int max = 1 << 4;
        int nBits = MiscMath.numberOfBits(max);

        for (int i = 0; i <= max; i++) {

            String bitstring = Integer.toBinaryString(i);
            while (bitstring.length() < nBits) {
                bitstring = "0" + bitstring;
            }
            char[] revBitstring = bitstring.toCharArray();
            for (int ii = 0; ii < (bitstring.length()/2); ii++) {
                int idx2 = bitstring.length() - ii - 1;
                char swap = revBitstring[ii];
                revBitstring[ii] = revBitstring[idx2];
                revBitstring[idx2] = swap;
            }

            String rBitstring = String.valueOf(revBitstring);

            int rev = Integer.parseInt(rBitstring, 2);

            int result = MiscMath.bitReverse(i, nBits);

            assertTrue(rev == result);
        }

    }

    public void testIsAPowerOf2() throws Exception {

        assertTrue(MiscMath.isAPowerOf2(2));

        assertTrue(MiscMath.isAPowerOf2(1<<2));

        assertTrue(MiscMath.isAPowerOf2(1<<7));

        assertFalse(MiscMath.isAPowerOf2(71));

        assertFalse(MiscMath.isAPowerOf2(17));
    }

    public void testCalculatePolarTheta() throws Exception {
        //public static double calculatePolarTheta(float x, float y) {

        /*
         * <pre>    90(=pi/2)
         *            |
         *            |
         *   180 ----------- 0
         *   (=pi)    |
         *            |
         *          270(=3pi/2)
         * </pre>
        */

        double eps = 0.001;

        float x, y;
        double theta, expectedTheta;

        x = 1;
        y = 0;
        expectedTheta = 0;
        theta = MiscMath.calculatePolarTheta(x, y);
        assertTrue(Math.abs(expectedTheta - theta) < eps);

        x = 0;
        y = 1;
        expectedTheta = 90. * (Math.PI/180.);
        theta = MiscMath.calculatePolarTheta(x, y);
        assertTrue(Math.abs(expectedTheta - theta) < eps);

        x = -1;
        y = 0;
        expectedTheta =  180. * (Math.PI/180.);
        theta = MiscMath.calculatePolarTheta(x, y);
        assertTrue(Math.abs(expectedTheta - theta) < eps);

        x = 0;
        y = -1;
        expectedTheta =  270. * (Math.PI/180.);
        theta = MiscMath.calculatePolarTheta(x, y);
        assertTrue(Math.abs(expectedTheta - theta) < eps);

        double[] d = new double[] {
            30. * (Math.PI/180.),
            60. * (Math.PI/180.),
            120. * (Math.PI/180.),
            150. * (Math.PI/180.),
            210. * (Math.PI/180.),
            240. * (Math.PI/180.),
            300. * (Math.PI/180.),
            330. * (Math.PI/180.)
        };

        for (int i = 3; i < d.length; i++) {
            x = (float)Math.cos(d[i]);
            y = (float)Math.sin(d[i]);
            expectedTheta = d[i];
            theta = MiscMath.calculatePolarTheta(x, y);
            assertTrue(Math.abs(expectedTheta - theta) < eps);
        }

    }

    /*
    public static int[] writeDegreeIntervals(int rotStart, int rotStop,
        int rotDelta) {
    */
    public void testWriteDegreeIntervals() throws Exception {

        float eps = 0.001f;
        float startRot = 0;
        float stopRot = 359;
        float deltaRot = 10;
        float[] expectedRot = new float[]{
            0, 10, 20, 30, 40, 50, 60, 70, 80, 90,
            100, 110, 120, 130, 140, 150, 160, 170, 180, 190,
            200, 210, 220, 230, 240, 250, 260, 270, 280, 290,
            300, 310, 320, 330, 340, 350
        };
        float[] result = MiscMath.writeDegreeIntervals(startRot, stopRot, deltaRot);
        for (int i = 0; i < result.length; ++i) {
            float r = result[i];
            float e = expectedRot[i];
            assertTrue(Math.abs(r - e) < eps);
        }

        startRot = 350;
        stopRot = 5;
        deltaRot = 1;
        expectedRot = new float[]{350, 351, 352, 353, 354, 355, 356, 357, 358,
            359, 0, 1, 2, 3, 4, 5};
        result = MiscMath.writeDegreeIntervals(startRot, stopRot, deltaRot);
        for (int i = 0; i < result.length; ++i) {
            float r = result[i];
            float e = expectedRot[i];
            assertTrue(Math.abs(r - e) < eps);
        }
    }
    
    public void testWriteToBigEndianBytes() {
                
        for (long i = 0; i < 62; ++i) {
            
            long v = 1 << i;
            
            byte[] bytes = MiscMath.writeToBigEndianBytes(v);
            
            BigInteger b = new BigInteger(bytes);
            
            long r = b.longValueExact();
            
            assertTrue(r == v);
        }
    }
    
    public void testLexicographicallyOrderedPoints() throws Exception {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        log.info("SEED=" + seed);
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        int limit = 1 << 10;
        
        for (int i = 0; i < 1000; ++i) {
            int x = sr.nextInt(limit);
            int y = sr.nextInt(limit);
            
            points.add(new PairInt(x, y));
            
            // occassionally, add a few with same x but different y
            if (((i % 100) == 0) && sr.nextBoolean()) {
                int prevY = y;
                while (y == prevY) {
                    y = sr.nextInt(limit);
                    points.add(new PairInt(x, y));
                }
            }
        }
                
        Set<PairInt> points0 = new HashSet<PairInt>(points);
        Set<PairInt> points1 = new HashSet<PairInt>(points);
        
        LinkedHashSet<PairInt> ordered0 = MiscMath.lexicographicallyOrderPointsByScan(
            points0, limit, limit);
        
        LinkedHashSet<PairInt> ordered1 = MiscMath.lexicographicallyOrderPointsBySort(
            points1);
        
        int prevX = Integer.MIN_VALUE;
        int prevY = Integer.MIN_VALUE;
        for (PairInt p : ordered0) {
            int x = p.getX();
            int y = p.getY();
            assertTrue((prevX < x) || ((prevX == x) && (prevY < y)));
            prevX = x;
            prevY = y;
        }
        
        prevX = Integer.MIN_VALUE;
        prevY = Integer.MIN_VALUE;
        for (PairInt p : ordered1) {
            int x = p.getX();
            int y = p.getY();
            assertTrue((prevX < x) || ((prevX == x) && (prevY < y)));
            prevX = x;
            prevY = y;
        }
    }
    
    public void testRescale() throws Exception {
        /*
        public static int[] rescale(double[] a, int vi, int vf) {
        */
        
        double[] a = new double[]{2, 3, 4, 5};
        int[] scaled = MiscMath.rescale(a, 0, 10);
        
        assertEquals(4, scaled.length);
        
        assertEquals(0, scaled[0]);
        assertEquals(3, scaled[1]);
        assertEquals(7, scaled[2]);
        assertEquals(10, scaled[3]);
        
    }
}
