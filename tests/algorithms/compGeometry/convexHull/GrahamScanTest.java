package algorithms.compGeometry.convexHull;

import algorithms.util.PairIntArray;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.logging.Logger;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.framework.TestCase;

public class GrahamScanTest extends TestCase {

    protected static Logger log = Logger.getLogger(GrahamScanTest.class.getName());
    
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
        float[] x = new float[]{0, 2, 7, 2, 2, 3};
        float[] y = new float[]{2, 2, 1, 6, 0, 2};

        GrahamScan scan = new GrahamScan();
        scan.computeHull(x, y);
        float[] chx = scan.getXHull();
        float[] chy = scan.getYHull();

        assertTrue(chx.length == 5);
        assertTrue(chy.length == 5);

        // clockwise order
        float[] expectedxx = new float[]{0, 2, 7, 2, 0};
        float[] expectedyy = new float[]{2, 6, 1, 0, 2};

        for (int i = 0; i < expectedxx.length; i++) {
            assertTrue(Math.abs(expectedxx[i] - chx[i]) < 0.01);
            assertTrue(Math.abs(expectedyy[i] - chy[i]) < 0.01);
        }

        //--------------------
        x = new float[]{0, 2, 7, 2, 2, 3};
        y = new float[]{2, 2, 1, 6, 0, 2};
        
        PairIntArray xy = new PairIntArray();
        for (int i = 0; i < x.length; i++) {
            xy.add(Math.round(x[i]), Math.round(y[i]));
        }
        
    }
    
    public void testScanExceptions() throws Exception {

        boolean threwException = false;
        float[] x = new float[]{0, 2, 7, 2, 2, 3};
        float[] y = new float[]{2, 2, 1, 6, 0};
        
        try {
            GrahamScan scan = new GrahamScan();
            scan.computeHull(null, y);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        try {
            GrahamScan scan = new GrahamScan();
            scan.computeHull(x, null);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);

        
        try {
            GrahamScan scan = new GrahamScan();
            scan.computeHull(x, y);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        x = new float[]{0, 2};
        y = new float[]{2, 2};
        try {
            GrahamScan scan = new GrahamScan();
            scan.computeHull(x, y);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        
        float[] xx = new float[]{0, 2, 2};
        float[] yy = new float[]{2, 2, 2};
        try {
            GrahamScan scan = new GrahamScan();
            scan.computeHull(xx, yy);
        } catch (GrahamScanTooFewPointsException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        try {
            GrahamScan scan = new GrahamScan();
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

            GrahamScan scan = new GrahamScan();
            scan.computeHull(x, y);
            float[] xHull = scan.xHull;
            float[] yHull = scan.yHull;
            
            plotter.addPlot(x, y, xHull, yHull, "gs");
            plotter.writeFile();
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
        return new TestSuite(GrahamScanTest.class);
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
