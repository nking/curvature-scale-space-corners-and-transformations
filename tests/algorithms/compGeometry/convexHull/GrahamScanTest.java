package algorithms.compGeometry.convexHull;

import algorithms.util.PairIntArray;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.logging.Logger;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.framework.TestCase;

/**
 * adapted from 
 * https://code.google.com/p/two-point-correlation/source/browse/src/main/java/algorithms/compGeometry/convexHull/GrahamScan.java
 * under MIT License (MIT), Nichole King 2013
 */
public class GrahamScanTest extends TestCase {

    protected static Logger log = Logger.getLogger(GrahamScanTest.class.getName());
    
    /**
     * @see TestCase#setUp()
     */
    protected void setUp() throws Exception {
        super.setUp();
    }

    /**
     * @see TestCase#tearDown()
     */
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void testScan() throws Exception {

        /*            2,6
         *
         *
         *     0,2    2,2 3,2
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
        scan.computeHull(x,y);
        float[] chx = scan.xHull;
        float[] chy = scan.yHull;

        assertTrue(chx.length == 5);
        assertTrue(chy.length == 5);

        // reversed because the cvhull storage is a stack:
        float[] expectedxx = new float[]{2, 7, 2, 0, 2};
        float[] expectedyy = new float[]{0, 1, 6, 2, 0};

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
        
        GrahamScanInt scanInt = new GrahamScanInt();
        scanInt.computeHull(xy);
        int[] chxint = scanInt.xHull;
        int[] chyint = scanInt.yHull;

        assertTrue(chxint.length == 5);
        assertTrue(chyint.length == 5);
        
        for (int i = 0; i < expectedxx.length; i++) {
            assertTrue(Math.abs(expectedxx[i] - chxint[i]) < 0.01);
            assertTrue(Math.abs(expectedyy[i] - chyint[i]) < 0.01);
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

            GrahamScanInt scanInt = new GrahamScanInt();
            scanInt.computeHull(xy);
            int[] xHull = scanInt.xHull;
            int[] yHull = scanInt.yHull;

            plotter.addPlot(xy.getX(), xy.getY(), xHull, yHull, "gs int");
            plotter.writeFile();
        }
        
    }

    public void testCalculateConvexHull7() throws Exception {

        float[] x = new float[]{2503.197f, 2364.4219f, 2562.4644f, 2562.4644f, 2335.7095f};
        float[] y = new float[]{638.1283f,  783.9287f,  799.73816f, 799.73816f, 1011.5769f};
        float xMin = 0;
        float yMin = 0;
        float xMax = 3000;
        float yMax = 3000;

        GrahamScan scan = new GrahamScan();
        scan.computeHull(x, y);
        float[] xHull = scan.xHull;
        float[] yHull = scan.yHull;

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xMin, xMax, yMin, yMax);
        plotter.addPlot(x, y, xHull, yHull, null);
        plotter.writeFile();


        x = new float[]{2288.6953f, 2260.7102f, 2261.9114f, 2105.9382f,
            2370.0544f, 2252.9143f, 2226.3518f, 2326.1133f, 2416.4487f,
            2311.942f, 2313.0254f, 2302.8582f, 2424.5393f, 2345.3005f,
            2368.0608f, 2316.284f, 2286.185f, 2309.2407f, 2072.8574f,
            2185.007f, 2190.7693f, 2190.7693f, 2520.3206f, 2511.2515f,
            2400.5505f, 2313.9111f, 2305.864f, 2458.3481f, 2172.3674f,
            2172.3674f, 2300.0789f, 2122.8271f, 2385.6226f, 2385.6226f};
        y = new float[]{1377.0491f, 1425.5775f, 1429.8291f, 1484.899f,
            1491.344f, 1502.3098f, 1513.2986f, 1550.7242f, 1565.6873f,
            1585.4646f, 1588.9802f, 1592.6627f, 1597.1661f, 1597.2826f,
            1604.8295f, 1606.774f, 1613.5851f, 1637.6198f, 1537.3373f,
            1655.0707f, 1668.9308f, 1668.9308f, 1655.198f, 1666.6757f,
            1692.4037f, 1715.7584f, 1727.0739f, 1732.0771f, 1754.0664f,
            1754.0664f, 1767.475f, 1760.0884f, 1823.989f, 1823.989f};

        scan = new GrahamScan();
        scan.computeHull(x, y);
        xHull = scan.xHull;
        yHull = scan.yHull;

        plotter.addPlot(x, y, xHull, yHull, null);
        plotter.writeFile();


        x = new float[]{1327.5411f, 1351.577f, 1340.271f, 1429.236f, 1129.4193f, 1265.5983f, 1323.0608f, 1518.2924f, 1169.6117f, 1115.9628f, 1356.29f, 1199.6755f, 1280.1757f, 1210.8116f, 1393.1355f, 1309.0094f, 1141.6178f, 1264.3749f, 1293.9006f, 1289.5648f, 1306.3184f, 1254.3723f, 1290.6552f, 1320.1703f, 1272.5048f, 1341.5015f, 1272.0065f, 1272.0065f, 1553.6871f, 1337.17f, 1157.7848f, 1079.6444f, 1211.8549f, 1285.4514f, 1285.4514f, 1288.3295f, 1347.0162f, 1377.6075f, 1449.616f, 1528.4755f, 1505.8534f, 1505.8534f, 1303.468f, 1303.468f};
        y = new float[]{364.9084f, 364.9535f, 378.23294f, 428.76364f, 468.40115f, 483.19467f, 483.31006f, 507.48907f, 507.96622f, 508.03027f, 517.006f, 526.424f, 527.9386f, 530.8042f, 534.5966f, 553.1233f, 553.1453f, 556.70966f, 577.2935f, 577.3543f, 587.17694f, 592.3263f, 598.7822f, 608.82764f, 610.1282f, 613.0668f, 616.05597f, 616.05597f, 607.1783f, 662.53656f, 643.2337f, 644.67f, 655.7444f, 671.61334f, 671.61334f, 693.8654f, 699.3863f, 720.7201f, 705.19806f, 681.2691f, 697.5069f, 697.5069f, 784.3944f, 784.3944f};
        scan = new GrahamScan();
        scan.computeHull(x, y);
        xHull = scan.xHull;
        yHull = scan.yHull;

        plotter.addPlot(x, y, xHull, yHull, null);
        plotter.writeFile();


        x = new float[]{317.7767f, 374.70038f, 193.74489f, 379.76953f, 233.85587f, 380.68527f, 267.06616f, 257.47916f, 364.01483f, 441.82767f, 374.0536f, 406.5328f, 319.8891f, 259.20575f, 419.7319f, 404.36636f, 354.01508f, 328.1624f, 421.01f, 270.25766f, 344.21442f, 146.69083f, 258.98032f, 351.951f, 467.74384f, 286.0171f, 257.51688f, 315.08923f, 195.02649f, 319.97906f, 412.0982f, 335.4689f, 348.5402f, 339.97723f, 376.19168f, 333.0743f, 519.27f, 502.12686f, 470.9786f, 394.68344f, 122.9748f, 130.33989f, 320.04355f, 217.01782f, 217.01782f, 310.62234f, 310.62234f, 309.89655f, 309.89655f, 197.49863f, 241.50131f, 328.87814f, 549.2164f, 484.4831f, 484.4831f, 271.2044f, 307.63132f, 93.54473f, 194.36359f, 373.30283f, 373.30283f, 452.29074f, 452.29074f, 363.5025f, 351.63565f, 284.12476f, 284.12476f};
        y = new float[]{1287.3486f, 1295.4368f, 1341.8411f, 1353.923f, 1357.3695f, 1362.5659f, 1364.9083f, 1374.9125f, 1400.0546f, 1408.0156f, 1410.1508f, 1412.9609f, 1415.606f, 1445.7646f, 1447.352f, 1449.7291f, 1449.7787f, 1452.8091f, 1461.6896f, 1466.5848f, 1474.3411f, 1475.45f, 1478.2322f, 1490.5894f, 1490.644f, 1495.9879f, 1496.1282f, 1519.1797f, 1521.6104f, 1524.1167f, 1530.9297f, 1537.4369f, 1539.5703f, 1541.8643f, 1547.5017f, 1547.8535f, 1483.9647f, 1516.1166f, 1533.2858f, 1550.3828f, 1550.9169f, 1554.0568f, 1559.364f, 1563.0688f, 1563.0688f, 1599.6744f, 1599.6744f, 1619.5978f, 1619.5978f, 1624.3027f, 1640.6165f, 1641.0421f, 1600.8939f, 1661.6283f, 1661.6283f, 1669.3585f, 1673.9454f, 1624.8038f, 1701.1073f, 1703.9318f, 1703.9318f, 1716.5507f, 1716.5507f, 1742.5619f, 1756.4325f, 1751.3546f, 1751.3546f};
        scan = new GrahamScan();
        scan.computeHull(x, y);
        xHull = scan.xHull;
        yHull = scan.yHull;

        plotter.addPlot(x, y, xHull, yHull, null);
        plotter.writeFile();


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

    private void createRandomCoordinates(double[] x, double[] y, double xMin, double xMax, double yMin, double yMax) 
        throws NoSuchAlgorithmException {

        SecureRandom r = SecureRandom.getInstance("SHA1PRNG");
        r.setSeed(System.nanoTime());

        int n = x.length;

        double xrange = xMax - xMin;
        double yrange = yMax - yMin;

        for (int i = 0; i < n; i++) {
            double r1 = r.nextDouble();
            double r2 = r.nextDouble();
            x[i] = (xrange * r1) + xMin;
            y[i] = (yrange * r2) + yMin;
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
