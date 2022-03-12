package algorithms.compGeometry.convexHull;

import algorithms.compGeometry.convexHull.GrahamScanLong.CH;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.util.PairIntArray;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.logging.Logger;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.framework.TestCase;

public class GrahamScanLongTest extends TestCase {

    protected static Logger log = Logger.getLogger(GrahamScanLongTest.class.getName());
    
    public void testScan() throws Exception {

        /*            2,6*
         *
         *
         *     0,2*   2,2 3,2
         *                      7,1*
         *            2,0*
         *         
         */
        int n = 6;
        long[] x = new long[]{0, 2, 7, 2, 2, 3};
        long[] y = new long[]{2, 2, 1, 6, 0, 2};
        
        CH ch = GrahamScanLong.computeHull(x, y);
        
        //System.out.printf("ch=%s\n", ch.toString());
        
        assertEquals(5, ch.getXH().length);

        // clockwise order
        long[] expectedxx = new long[]{2, 7, 2, 0, 2};
        long[] expectedyy = new long[]{0, 1, 6, 2, 0};

        for (int i = 0; i < expectedxx.length; i++) {
            assertEquals(expectedxx[i], ch.getXH()[i]);
            assertEquals(expectedyy[i], ch.getYH()[i]);
        }

    }
    
    public void testScanExceptions() throws Exception {

        boolean threwException = false;
        
        try {
            CH ch = GrahamScanLong.computeHull(null, null);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        long[] x = new long[]{0, 2};
        long[] y = new long[]{2, 2};
        
        try {
            CH ch = GrahamScanLong.computeHull(x, y);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        x = new long[]{0, 2, 2};
        y = new long[]{2, 2, 2};
        
        try {
            CH ch = GrahamScanLong.computeHull(x, y);
        } catch (GrahamScanTooFewPointsException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
    }

    public void testCalculateConvexHull6() throws Exception {

        int ntries = 1;
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.nanoTime();
        System.out.println("seed=" + seed);
        sr.setSeed(seed);

        long xMin = 10;
        long yMin = 10;
        long xMax = 1000;
        long yMax = 1000;
            
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xMin, xMax, yMin, yMax);
        
        for (int i = 0; i < ntries; i++) {
            int n = 1000;
            long[] x = new long[n];
            long[] y = new long[n];
            
            float[] xF = new float[n];
            float[] yF = new float[n];

            float maxRadius = 200;

            for (int ii = 0; ii < n; ii++) {

                float radius = maxRadius * sr.nextFloat();
                double angle = 360. * sr.nextDouble();

                float[] xy = calculateXAndYFromXcYcAndRadius(600.f, 400.f, radius, angle);

                xF[ii] = xy[0];
                yF[ii] = xy[1];
                x[ii] = (long)(Math.round(xy[0]));
                y[ii] = (long)(Math.round(xy[1]));
            }
            
            CH ch = GrahamScanLong.computeHull(x, y);
            
            int nH = ch.getXH().length;
                        
            PairIntArray curve = new PairIntArray(nH - 1);
            float[] xHull = new float[nH];
            float[] yHull = new float[nH];
            for (int ii = 0; ii < nH; ++ii) {
                xHull[ii] = ch.getXH()[ii];
                yHull[ii] = ch.getYH()[ii];
                if (ii < (nH - 1)) {
                    curve.add((int)xHull[ii], (int)yHull[ii]);
                }
            }
            plotter.addPlot(xF, yF, xHull, yHull, "gs");
            plotter.writeFile();
            
            MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
            boolean isCCW = !curveHelper.curveIsOrderedClockwise(curve);
            assertTrue(isCCW);
        }

        for (int i = 0; i < ntries; i++) {
            int n = 1000;
            PairIntArray xy = new PairIntArray(n);
            
            int maxRadius = 200;

            createRandomPointsAroundCenter(maxRadius, n, 600, 400, xy);

        }
        
    }
    
    protected void createRandomPointsAroundCenter(int maxRadius,
        int numberOfPoints, int xc, int yc, PairIntArray xyPoints) throws 
        NoSuchAlgorithmException {

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(System.nanoTime());

        for (int i = 0; i < numberOfPoints; i++) {

            float radius = maxRadius * sr.nextFloat();
            double angle = 360. * sr.nextDouble();

            float[] xy = calculateXAndYFromXcYcAndRadius((float)xc, (float)yc, 
                radius, angle);

            xyPoints.add(Math.round(xy[0]), Math.round(xy[1]));
        }
    }

    static float[] calculateXAndYFromXcYcAndRadiusCCW(float xc, float yc, float radius, double angleInDegreesFromYEQ0XGT0) {

        return calculateXAndYFromXcYcAndRadius(xc, yc, radius, 360 - angleInDegreesFromYEQ0XGT0);
    }
    
    public static float[] calculateXAndYFromXcYcAndRadius(float xc, float yc, float radius, double angleInDegreesFromYEQ0XGT0) {

        double dx = radius * Math.cos(angleInDegreesFromYEQ0XGT0 * (Math.PI/180.f));
        double dy = radius * Math.sin(angleInDegreesFromYEQ0XGT0 * (Math.PI/180.f));

        float x = (float) (xc + dx);
        float y = (float) (yc - dy);
        return new float[]{x, y};
    }

    /**
     * Test suite
     *
     * @return static Test
     */
    public static Test suite() {
        log.fine("Creating a TestSuite for GrahamScan");
        return new TestSuite(GrahamScanLongTest.class);
    }

    /**
     * Set up a Junit test runner
     *
     * @param args Not used.
     */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(suite());
    }
}
