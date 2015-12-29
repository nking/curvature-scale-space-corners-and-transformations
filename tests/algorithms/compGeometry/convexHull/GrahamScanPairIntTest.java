package algorithms.compGeometry.convexHull;

import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.logging.Logger;
import algorithms.util.PolygonAndPointPlotter;
import java.util.List;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.framework.TestCase;

public class GrahamScanPairIntTest extends TestCase {

    protected static Logger log = Logger.getLogger(GrahamScanPairIntTest.class.getName());
    
    public void testScan() throws Exception {

        /*            2,6
         *
         *
         *     0,2   2,2 3,2
         *                      7,1
         *            2,0
         *    7
         *    6   <>
         *    5   
         *    4       
         *    3       
         *    <> <> <>        
         *    1             <>
         *    0 1<> 3 4 5 6 7        
         */
        PairInt[] points = new PairInt[6];
        points[0] = new PairInt(0, 2);
        points[1] = new PairInt(2, 2);
        points[2] = new PairInt(7, 1);
        points[3] = new PairInt(2, 6);
        points[4] = new PairInt(2, 0);
        points[5] = new PairInt(3, 2);

        GrahamScanPairInt<PairInt> scan = new GrahamScanPairInt<PairInt>();
        scan.computeHull(points);
        List<PairInt> hull = scan.getHull();

        assertTrue(hull.size() == 5);

        // clockwise order
        float[] expectedxx = new float[]{0, 2, 7, 2, 0};
        float[] expectedyy = new float[]{2, 6, 1, 0, 2};

        for (int i = 0; i < expectedxx.length; i++) {
            assertTrue(Math.abs(expectedxx[i] - hull.get(i).getX()) < 0.01);
            assertTrue(Math.abs(expectedyy[i] - hull.get(i).getY()) < 0.01);
        }

        PairIntArray curve = new PairIntArray(hull.size());
        for (int i = 0; i < hull.size(); ++i) {
            int x = hull.get(i).getX();
            int y = hull.get(i).getY();
            curve.add(x, y);
        }
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        boolean isCW = curveHelper.curveIsOrderedClockwise(curve);
        assertTrue(isCW);
        
        points = new PairInt[5];
        points[0] = new PairInt(0, 2);
        points[1] = new PairInt(2, 6);
        points[2] = new PairInt(7, 1);
        points[3] = new PairInt(2, 0);
        points[4] = new PairInt(0, 2);
        curve = new PairIntArray(points.length);
        for (int i = 0; i < points.length; ++i) {
            int x = points[i].getX();
            int y = points[i].getY();
            curve.add(x, y);
        }
        isCW = curveHelper.curveIsOrderedClockwise2(curve);
        assertTrue(isCW);
    }
    
    public void testScanExceptions() throws Exception {

        boolean threwException = false;
        
        try {
            GrahamScanPairInt<PairInt> scan = new GrahamScanPairInt<PairInt>();
            scan.computeHull(null);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        PairInt[] points = new PairInt[2];
        points[0] = new PairInt(0, 2);
        points[1] = new PairInt(2, 2);
        
        try {
            GrahamScanPairInt<PairInt> scan = new GrahamScanPairInt<PairInt>();
            scan.computeHull(points);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        points = new PairInt[3];
        points[0] = new PairInt(0, 2);
        points[1] = new PairInt(2, 2);
        points[2] = new PairInt(2, 2);
        
        try {
            GrahamScanPairInt<PairInt> scan = new GrahamScanPairInt<PairInt>();
            scan.computeHull(points);
        } catch (GrahamScanTooFewPointsException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        try {
            GrahamScanPairInt<PairInt> scan = new GrahamScanPairInt<PairInt>();
            scan.populateHull();
        } catch (GrahamScanTooFewPointsException e) {
            threwException = true;
        }
        assertTrue(threwException);
    }

    public void testCalculateConvexHull6() throws Exception {

        int ntries = 1;

        float xMin = 10;
        float yMin = 10;
        float xMax = 1000;
        float yMax = 1000;
            
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xMin, xMax, yMin, yMax);
        
        for (int i = 0; i < ntries; i++) {
            int n = 1000;
            float[] x = new float[n];
            float[] y = new float[n];

            float maxRadius = 200;

            createRandomPointsAroundCenter(maxRadius, n, 600.f, 400.f, x, y, 0);

            GrahamScanPairInt<PairInt> scan = new GrahamScanPairInt<PairInt>();
            
            PairInt[] points = new PairInt[n];
            for (int ii = 0; ii < n; ++ii) {
                points[ii] = new PairInt(Math.round(x[ii]), Math.round(y[ii]));
            }
            
            scan.computeHull(points);
            
            List<PairInt> hull = scan.getHull();
            
            PairIntArray curve = new PairIntArray(hull.size() - 1);
            float[] xHull = new float[hull.size()];
            float[] yHull = new float[hull.size()];
            for (int ii = 0; ii < hull.size(); ++ii) {
                xHull[ii] = hull.get(ii).getX();
                yHull[ii] = hull.get(ii).getY();
                if (ii < (hull.size() - 1)) {
                    curve.add(hull.get(ii).getX(), hull.get(ii).getY());
                }
            }
            
            plotter.addPlot(x, y, xHull, yHull, "gs");
            plotter.writeFile();
            
            MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
            boolean isCW = curveHelper.curveIsOrderedClockwise(curve);
            assertTrue(isCW);
        }

        for (int i = 0; i < ntries; i++) {
            int n = 1000;
            PairIntArray xy = new PairIntArray(n);
            
            int maxRadius = 200;

            createRandomPointsAroundCenter(maxRadius, n, 600, 400, xy);

        }
        
    }

    protected void createRandomPointsAroundCenter(float maxRadius,
        int numberOfPoints, float xc, float yc, double[] x, double[] y, 
        int xyStartOffset) throws NoSuchAlgorithmException {

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(System.nanoTime());

        for (int i = 0; i < numberOfPoints; i++) {

            float radius = maxRadius * sr.nextFloat();
            double angle = 360. * sr.nextDouble();

            float[] xy = calculateXAndYFromXcYcAndRadius(xc, yc, radius, angle);

            x[xyStartOffset + i] = xy[0];
            y[xyStartOffset + i] = xy[1];
        }
    }

    protected void createRandomPointsAroundCenter(float maxRadius,
        int numberOfPoints, float xc, float yc, float[] x, float[] y, 
        int xyStartOffset) throws NoSuchAlgorithmException {

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(System.nanoTime());

        for (int i = 0; i < numberOfPoints; i++) {

            float radius = maxRadius * sr.nextFloat();
            double angle = 360. * sr.nextDouble();

            float[] xy = calculateXAndYFromXcYcAndRadius(xc, yc, radius, angle);

            x[xyStartOffset + i] = xy[0];
            y[xyStartOffset + i] = xy[1];
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
        return new TestSuite(GrahamScanPairIntTest.class);
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
