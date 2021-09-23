package algorithms.util;

import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.misc.MiscMath;
import algorithms.statistics.CDFStandardNormal;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class LinearRegressionTest extends TestCase {
    
    public LinearRegressionTest() {
        super(LinearRegressionTest.class.getName());
    }
    
    long seed;
    SecureRandom sr = null;
    
    public void setUp() throws NoSuchAlgorithmException {
        
        seed = System.currentTimeMillis();
        System.out.println("SEED=" + seed);
        
        sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(seed);
    }
    
    public void tearDown() {
    }

    public void test0() throws Exception {
       
        // 2 parallel diagonal lines
        
        PairIntArray dxdy = new PairIntArray();
        double slope = 2.0;
        double yIntercept1 = (110 - slope*10);
        double yIntercept2 = (120 - slope*10);
        for (int x = 10; x < 100; x++) {
            /*
            (y1 - y0)/(x1 - x0) = slope
            y1 - y0 = slope*(x1 - x0);
            y1 = y0 + slope*(x1 - x0);
            y1 = (y0 - slope*x0) + slope*x1
            */
            double dy = slope*(double)x;
            int y = (int)Math.round(yIntercept1 + dy);
            
            dxdy.add(x, y);
            y = (int)Math.round(yIntercept2 + dy);
            dxdy.add(x, y);
        }
        
        LinearRegression instance = new LinearRegression();
        instance.plotTheLinearRegression(dxdy.getX(), dxdy.getY());
        
        float[] yInterceptAndSlope = 
            instance.calculateTheilSenEstimatorParams(dxdy.getX(), dxdy.getY());
        
        double expectedYIntercept = (yIntercept1 + yIntercept2)/2.;
        double expectedSlope = slope;
        
        assertTrue(Math.abs(yInterceptAndSlope[0] - expectedYIntercept) < 1);
        
        assertTrue(Math.abs(yInterceptAndSlope[1] - expectedSlope) < 0.1);
        
    }
    
    public void test00() throws Exception {
       
        // small random deviations from a straight line
                         
        float slope = 2.0f;
        float yIntercept1 = (110.f - slope*10f);
        
        float s = 0.25f;
        PairFloatArray dxdy = generateLine(0, yIntercept1, slope, s, 10, 100);
        
        
        LinearRegression instance = new LinearRegression();
        instance.plotTheLinearRegression(dxdy.getX(), dxdy.getY());
        
        float[] yInterceptAndSlope = 
            instance.calculateTheilSenEstimatorParams(dxdy.getX(), dxdy.getY());
        
        double expectedYIntercept = yIntercept1;
        double expectedSlope = slope;
        
        //assertTrue(Math.abs(yInterceptAndSlope[0] - expectedYIntercept) < Math.sqrt(d));
        
        assertTrue(Math.abs(yInterceptAndSlope[1] - expectedSlope) < 0.2);
        
    }
    
    public void test1() {
       
        PairIntArray dxdy = new PairIntArray();
        readThielSenTestData(dxdy);
        
        LinearRegression instance = new LinearRegression();
        instance.plotTheLinearRegression(dxdy.getX(), dxdy.getY());
        
        float[] yInterceptAndSlope = 
            instance.calculateTheilSenEstimatorParams(dxdy.getX(), dxdy.getY());
        
        if (true) {
            return;
        }
        assertTrue(Math.abs(yInterceptAndSlope[0] - 2900.f) < 50.);
        
        assertTrue(Math.abs(yInterceptAndSlope[1] - 1.0f) < 0.03);
        
        /*
        calculate the distance between the points and the regression line
           as an average and a standard deviation.
        */
        float lineX0 = MiscMath.findMin(dxdy.getX());
        float lineX1 = MiscMath.findMax(dxdy.getX());
        float lineY0 = yInterceptAndSlope[0] + yInterceptAndSlope[1] * lineX0;
        float lineY1 = yInterceptAndSlope[0] + yInterceptAndSlope[1] * lineX1;
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        double[] dist = new double[dxdy.getN()];
        for (int i = 0; i < dist.length; i++) {
            dist[i] = curveHelper.distanceFromPointToALine(lineX0, lineY0, 
                lineX1, lineY1, dxdy.getX(i), dxdy.getY(i));
        }
        
        double[] avgAndStDev = MiscMath.getAvgAndStDev(dist);
        
        // remove ~ 26 outliers out of 103 points
        int z = 1;
    }
    
    public void test0_float() throws Exception {
       
        // 2 parallel diagonal lines
        
        PairFloatArray dxdy = new PairFloatArray();
        float slope = 2.0f;
        float yIntercept1 = (110 - slope*10.f);
        float yIntercept2 = (120 - slope*10.f);
        
        for (int x = 10; x < 100; x++) {
            /*
            (y1 - y0)/(x1 - x0) = slope
            y1 - y0 = slope*(x1 - x0);
            y1 = y0 + slope*(x1 - x0);
            y1 = (y0 - slope*x0) + slope*x1
            */
            float dy = slope*(float)x;
            float y = yIntercept1 + dy;
            
            dxdy.add(x, y);
            y = yIntercept2 + dy;
            dxdy.add(x, y);
        }
        
        LinearRegression instance = new LinearRegression();
        instance.plotTheLinearRegression(dxdy.getX(), dxdy.getY());
        
        float[] yInterceptAndSlope = 
            instance.calculateTheilSenEstimatorParams(dxdy.getX(), dxdy.getY());
        
        double expectedYIntercept = (yIntercept1 + yIntercept2)/2.;
        double expectedSlope = slope;
        
        assertTrue(Math.abs(yInterceptAndSlope[0] - expectedYIntercept) < 1);
        
        assertTrue(Math.abs(yInterceptAndSlope[1] - expectedSlope) < 0.1);
        
    }
    
    private PairFloatArray generateLine(float x0, float y0, float slope, double sigma,
        float xStart, float xEnd) {
        /*
           |      slope=2
           |      .
           |     .
           |____._____     <-- x-intercept at x=10, y=0
               .
              .
             .
            . <--- y-intercept at x0=0, y0=-20
        
        slope=(y-y0)/(x-x0) => (x-x0)*slope=(y-y0) => y = y0 + (x-x0)*slope
        */
            
        double n01, xn, yn;
        double xi, yi;
        PairFloatArray dxdy = new PairFloatArray();
        for (float x = xStart; x < xEnd; x++) {
            yi = y0 + (x-x0)*slope;
            xi = x;
            if (sr.nextBoolean()) {
                //random sample from unit norm gaussian.
                n01 = CDFStandardNormal.approxInverseShort(sr.nextDouble());
                // and N(0,1) ~ (N(mean,sigma) - mean)/sigma
                // N(mean,sigma) = x + sigma*n01
                xi += (n01*sigma);
                n01 = CDFStandardNormal.approxInverseShort(sr.nextDouble());
                yi += (n01*sigma);
            }
            dxdy.add((float)xi, (float)yi);
        }
        return dxdy;
    }
    
    public void test00_float() throws Exception {
       
        // small random deviations from a straight line
        
        /*
           |      slope=2
           |      .
           |     .
           |____._____     <-- x-intercept at x=10, y=0
               .
              .
             .
            . <--- y-intercept at x0=0, y0=-20
        
        slope=(y-y0)/(x-x0) => (x-x0)*slope=(y-y0) => y = y0 + (x-x0)*slope
        */
                
        float slope = 2.0f;
        float yIntercept1 = -20f;
        float s = 0.25f;
        PairFloatArray dxdy = generateLine(0, yIntercept1, slope, s, 10, 100);
        
        LinearRegression instance = new LinearRegression();
        instance.plotTheLinearRegression(dxdy.getX(), dxdy.getY());
        
        float[] yInterceptAndSlope = 
            instance.calculateTheilSenEstimatorParams(dxdy.getX(), dxdy.getY());
        
        double expectedYIntercept = yIntercept1;
        double expectedSlope = slope;
        
        assertTrue(Math.abs(yInterceptAndSlope[0] 
            - expectedYIntercept) < Math.sqrt(s));
        
        assertTrue(Math.abs(yInterceptAndSlope[1] - expectedSlope) < 0.2);
        
    }
    
    public void test1_float() {
       
        PairFloatArray dxdy = new PairFloatArray();
        readThielSenTestData(dxdy);
        
        LinearRegression instance = new LinearRegression();
        instance.plotTheLinearRegression(dxdy.getX(), dxdy.getY());
        
        float[] yInterceptAndSlope = 
            instance.calculateTheilSenEstimatorParams(dxdy.getX(), dxdy.getY());
        
        if (true) {
            return;
        }
        assertTrue(Math.abs(yInterceptAndSlope[0] - 2900.f) < 50.);
        
        assertTrue(Math.abs(yInterceptAndSlope[1] - 1.0f) < 0.03);
        
        /*
        calculate the distance between the points and the regression line
           as an average and a standard deviation.
        */
        float lineX0 = MiscMath.findMin(dxdy.getX());
        float lineX1 = MiscMath.findMax(dxdy.getX());
        float lineY0 = yInterceptAndSlope[0] + yInterceptAndSlope[1] * lineX0;
        float lineY1 = yInterceptAndSlope[0] + yInterceptAndSlope[1] * lineX1;
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        double[] dist = new double[dxdy.getN()];
        for (int i = 0; i < dist.length; i++) {
            dist[i] = curveHelper.distanceFromPointToALine(lineX0, lineY0, 
                lineX1, lineY1, dxdy.getX(i), dxdy.getY(i));
        }
        
        double[] avgAndStDev = MiscMath.getAvgAndStDev(dist);
        
        // remove ~ 26 outliers out of 103 points
        int z = 1;
    }
    
    
    /*
    using data points from the svg available at:
    http://commons.wikimedia.org/wiki/File:Thiel-Sen_estimator.svg
    */
    private void readThielSenTestData(PairIntArray xy1) {
        
        xy1.add((int) Math.round(28.347), (int) Math.round(3050.818));
        xy1.add((int) Math.round(56.693), (int) Math.round(2940.701));
        xy1.add((int) Math.round(85.04), (int) Math.round(2942.112));
        xy1.add((int) Math.round(113.386), (int) Math.round(2922.829));
        xy1.add((int) Math.round(141.732), (int) Math.round(2817.593));
        xy1.add((int) Math.round(170.079), (int) Math.round(2862.549));
        xy1.add((int) Math.round(198.425), (int) Math.round(2783.021));
        xy1.add((int) Math.round(226.772), (int) Math.round(2829.397));
        xy1.add((int) Math.round(255.118), (int) Math.round(2710.108));
        xy1.add((int) Math.round(283.465), (int) Math.round(2756.589));
        xy1.add((int) Math.round(311.811), (int) Math.round(2716.387));
        xy1.add((int) Math.round(340.158), (int) Math.round(2735.223));
        xy1.add((int) Math.round(368.504), (int) Math.round(2439.297));
        xy1.add((int) Math.round(396.851), (int) Math.round(2474.149));
        xy1.add((int) Math.round(425.197), (int) Math.round(2390.326));
        xy1.add((int) Math.round(453.543), (int) Math.round(2410.0));
        xy1.add((int) Math.round(481.89), (int) Math.round(2571.953));
        xy1.add((int) Math.round(510.236), (int) Math.round(2320.793));
        xy1.add((int) Math.round(538.583), (int) Math.round(2334.67));
        xy1.add((int) Math.round(566.929), (int) Math.round(2366.697));
        xy1.add((int) Math.round(595.276), (int) Math.round(2252.961));
        xy1.add((int) Math.round(623.622), (int) Math.round(2392.574));
        xy1.add((int) Math.round(651.969), (int) Math.round(2368.743));
        xy1.add((int) Math.round(680.315), (int) Math.round(2244.386));
        xy1.add((int) Math.round(708.662), (int) Math.round(2206.536));
        xy1.add((int) Math.round(737.008), (int) Math.round(2073.293));
        xy1.add((int) Math.round(765.354), (int) Math.round(2270.166));
        xy1.add((int) Math.round(793.701), (int) Math.round(2052.446));
        xy1.add((int) Math.round(822.047), (int) Math.round(2195.055));
        xy1.add((int) Math.round(850.394), (int) Math.round(2187.155));
        xy1.add((int) Math.round(878.74), (int) Math.round(2189.327));
        xy1.add((int) Math.round(907.087), (int) Math.round(2152.691));
        xy1.add((int) Math.round(935.433), (int) Math.round(1885.666));
        xy1.add((int) Math.round(963.78), (int) Math.round(1899.318));
        xy1.add((int) Math.round(992.126), (int) Math.round(1897.532));
        xy1.add((int) Math.round(1020.473), (int) Math.round(1925.292));
        xy1.add((int) Math.round(1048.819), (int) Math.round(1855.75));
        xy1.add((int) Math.round(1077.166), (int) Math.round(1831.943));
        xy1.add((int) Math.round(1105.512), (int) Math.round(1706.527));
        xy1.add((int) Math.round(1133.858), (int) Math.round(1730.246));
        xy1.add((int) Math.round(1162.205), (int) Math.round(1762.077));
        xy1.add((int) Math.round(1190.551), (int) Math.round(1844.738));
        xy1.add((int) Math.round(1218.898), (int) Math.round(1797.771));
        xy1.add((int) Math.round(1247.244), (int) Math.round(1753.082));
        xy1.add((int) Math.round(1275.591), (int) Math.round(1659.165));
        xy1.add((int) Math.round(1303.937), (int) Math.round(2902.67));
        xy1.add((int) Math.round(1332.284), (int) Math.round(1697.364));
        xy1.add((int) Math.round(1360.63), (int) Math.round(1495.082));
        xy1.add((int) Math.round(1388.977), (int) Math.round(1437.04));
        xy1.add((int) Math.round(1417.323), (int) Math.round(1440.659));
        xy1.add((int) Math.round(1445.669), (int) Math.round(1518.253));
        xy1.add((int) Math.round(1474.016), (int) Math.round(1359.016));
        xy1.add((int) Math.round(1502.363), (int) Math.round(1947.112));
        xy1.add((int) Math.round(1530.709), (int) Math.round(1500.887));
        xy1.add((int) Math.round(1559.055), (int) Math.round(1376.082));
        xy1.add((int) Math.round(1587.402), (int) Math.round(1377.504));
        xy1.add((int) Math.round(1615.748), (int) Math.round(1374.184));
        xy1.add((int) Math.round(1644.094), (int) Math.round(1567.445));
        xy1.add((int) Math.round(1672.441), (int) Math.round(1233.486));
        xy1.add((int) Math.round(1700.787), (int) Math.round(1342.985));
        xy1.add((int) Math.round(1729.135), (int) Math.round(1259.289));
        xy1.add((int) Math.round(1757.48), (int) Math.round(1171.99));
        xy1.add((int) Math.round(1785.826), (int) Math.round(1106.017));
        xy1.add((int) Math.round(1814.174), (int) Math.round(1194.191));
        xy1.add((int) Math.round(1842.52), (int) Math.round(1136.358));
        xy1.add((int) Math.round(1870.867), (int) Math.round(1048.459));
        xy1.add((int) Math.round(1899.213), (int) Math.round(1163.572));
        xy1.add((int) Math.round(1927.559), (int) Math.round(953.636));
        xy1.add((int) Math.round(1955.906), (int) Math.round(733.188));
        xy1.add((int) Math.round(1984.252), (int) Math.round(983.902));
        xy1.add((int) Math.round(2012.598), (int) Math.round(994.518));
        xy1.add((int) Math.round(2040.945), (int) Math.round(904.925));
        xy1.add((int) Math.round(2069.291), (int) Math.round(2826.269));
        xy1.add((int) Math.round(2097.639), (int) Math.round(222.962));
        xy1.add((int) Math.round(2125.984), (int) Math.round(778.727));
        xy1.add((int) Math.round(2154.33), (int) Math.round(919.042));
        xy1.add((int) Math.round(2182.678), (int) Math.round(846.093));
        xy1.add((int) Math.round(2211.023), (int) Math.round(2381.066));
        xy1.add((int) Math.round(2239.371), (int) Math.round(280.228));
        xy1.add((int) Math.round(2267.717), (int) Math.round(574.438));
        xy1.add((int) Math.round(2296.062), (int) Math.round(2033.601));
        xy1.add((int) Math.round(2324.41), (int) Math.round(746.172));
        xy1.add((int) Math.round(2352.756), (int) Math.round(2672.973));
        xy1.add((int) Math.round(2381.104), (int) Math.round(635.111));
        xy1.add((int) Math.round(2409.449), (int) Math.round(2871.514));
        xy1.add((int) Math.round(2437.795), (int) Math.round(1411.262));
        xy1.add((int) Math.round(2466.143), (int) Math.round(519.237));
        xy1.add((int) Math.round(2494.488), (int) Math.round(2792.811));
        xy1.add((int) Math.round(2522.834), (int) Math.round(2478.352));
        xy1.add((int) Math.round(2551.182), (int) Math.round(1248.189));
        xy1.add((int) Math.round(2579.527), (int) Math.round(57.313));
        xy1.add((int) Math.round(2607.875), (int) Math.round(411.902));
        xy1.add((int) Math.round(2636.221), (int) Math.round(1371.529));
        xy1.add((int) Math.round(2664.566), (int) Math.round(992.871));
        xy1.add((int) Math.round(2692.914), (int) Math.round(1382.086));
        xy1.add((int) Math.round(2721.26), (int) Math.round(2469.339));
        xy1.add((int) Math.round(2749.607), (int) Math.round(1194.956));
        xy1.add((int) Math.round(2777.953), (int) Math.round(2433.983));
        xy1.add((int) Math.round(2806.299), (int) Math.round(2840.234));
        xy1.add((int) Math.round(2834.646), (int) Math.round(377.596));
        xy1.add((int) Math.round(2862.992), (int) Math.round(528.729));
        xy1.add((int) Math.round(2891.338), (int) Math.round(2862.222));
        xy1.add((int) Math.round(2919.686), (int) Math.round(2778.425));
    }
   
    /*
    using data points from the svg available at:
    http://commons.wikimedia.org/wiki/File:Thiel-Sen_estimator.svg
    */
    private void readThielSenTestData(PairFloatArray xy1) {

        xy1.add(28.347f, 3050.818f);
        xy1.add(56.693f, 2940.701f);
        xy1.add(85.04f, 2942.112f);
        xy1.add(113.386f, 2922.829f);
        xy1.add(141.732f, 2817.593f);
        xy1.add(170.079f, 2862.549f);
        xy1.add(198.425f, 2783.021f);
        xy1.add(226.772f, 2829.397f);
        xy1.add(255.118f, 2710.108f);
        xy1.add(283.465f, 2756.589f);
        xy1.add(311.811f, 2716.387f);
        xy1.add(340.158f, 2735.223f);
        xy1.add(368.504f, 2439.297f);
        xy1.add(396.851f, 2474.149f);
        xy1.add(425.197f, 2390.326f);
        xy1.add(453.543f, 2410.0f);
        xy1.add(481.89f, 2571.953f);
        xy1.add(510.236f, 2320.793f);
        xy1.add(538.583f, 2334.67f);
        xy1.add(566.929f, 2366.697f);
        xy1.add(595.276f, 2252.961f);
        xy1.add(623.622f, 2392.574f);
        xy1.add(651.969f, 2368.743f);
        xy1.add(680.315f, 2244.386f);
        xy1.add(708.662f, 2206.536f);
        xy1.add(737.008f, 2073.293f);
        xy1.add(765.354f, 2270.166f);
        xy1.add(793.701f, 2052.446f);
        xy1.add(822.047f, 2195.055f);
        xy1.add(850.394f, 2187.155f);
        xy1.add(878.74f, 2189.327f);
        xy1.add(907.087f, 2152.691f);
        xy1.add(935.433f, 1885.666f);
        xy1.add(963.78f, 1899.318f);
        xy1.add(992.126f, 1897.532f);
        xy1.add(1020.473f, 1925.292f);
        xy1.add(1048.819f, 1855.75f);
        xy1.add(1077.166f, 1831.943f);
        xy1.add(1105.512f, 1706.527f);
        xy1.add(1133.858f, 1730.246f);
        xy1.add(1162.205f, 1762.077f);
        xy1.add(1190.551f, 1844.738f);
        xy1.add(1218.898f, 1797.771f);
        xy1.add(1247.244f, 1753.082f);
        xy1.add(1275.591f, 1659.165f);
        xy1.add(1303.937f, 2902.67f);
        xy1.add(1332.284f, 1697.364f);
        xy1.add(1360.63f, 1495.082f);
        xy1.add(1388.977f, 1437.04f);
        xy1.add(1417.323f, 1440.659f);
        xy1.add(1445.669f, 1518.253f);
        xy1.add(1474.016f, 1359.016f);
        xy1.add(1502.363f, 1947.112f);
        xy1.add(1530.709f, 1500.887f);
        xy1.add(1559.055f, 1376.082f);
        xy1.add(1587.402f, 1377.504f);
        xy1.add(1615.748f, 1374.184f);
        xy1.add(1644.094f, 1567.445f);
        xy1.add(1672.441f, 1233.486f);
        xy1.add(1700.787f, 1342.985f);
        xy1.add(1729.135f, 1259.289f);
        xy1.add(1757.48f, 1171.99f);
        xy1.add(1785.826f, 1106.017f);
        xy1.add(1814.174f, 1194.191f);
        xy1.add(1842.52f, 1136.358f);
        xy1.add(1870.867f, 1048.459f);
        xy1.add(1899.213f, 1163.572f);
        xy1.add(1927.559f, 953.636f);
        xy1.add(1955.906f, 733.188f);
        xy1.add(1984.252f, 983.902f);
        xy1.add(2012.598f, 994.518f);
        xy1.add(2040.945f, 904.925f);
        xy1.add(2069.291f, 2826.269f);
        xy1.add(2097.639f, 222.962f);
        xy1.add(2125.984f, 778.727f);
        xy1.add(2154.33f, 919.042f);
        xy1.add(2182.678f, 846.093f);
        xy1.add(2211.023f, 2381.066f);
        xy1.add(2239.371f, 280.228f);
        xy1.add(2267.717f, 574.438f);
        xy1.add(2296.062f, 2033.601f);
        xy1.add(2324.41f, 746.172f);
        xy1.add(2352.756f, 2672.973f);
        xy1.add(2381.104f, 635.111f);
        xy1.add(2409.449f, 2871.514f);
        xy1.add(2437.795f, 1411.262f);
        xy1.add(2466.143f, 519.237f);
        xy1.add(2494.488f, 2792.811f);
        xy1.add(2522.834f, 2478.352f);
        xy1.add(2551.182f, 1248.189f);
        xy1.add(2579.527f, 57.313f);
        xy1.add(2607.875f, 411.902f);
        xy1.add(2636.221f, 1371.529f);
        xy1.add(2664.566f, 992.871f);
        xy1.add(2692.914f, 1382.086f);
        xy1.add(2721.26f, 2469.339f);
        xy1.add(2749.607f, 1194.956f);
        xy1.add(2777.953f, 2433.983f);
        xy1.add(2806.299f, 2840.234f);
        xy1.add(2834.646f, 377.596f);
        xy1.add(2862.992f, 528.729f);
        xy1.add(2891.338f, 2862.222f);
        xy1.add(2919.686f, 2778.425f);
    }
}
