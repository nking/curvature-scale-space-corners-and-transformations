package algorithms.util;

import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.misc.MiscMath;
import java.security.SecureRandom;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class LinearRegressionTest {
    
    public LinearRegressionTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
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
    
    @Test
    public void test00() throws Exception {
       
        // small random deviations from a straight line
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(System.currentTimeMillis());
        
        int d = 10;
        
        PairIntArray dxdy = new PairIntArray();
        double slope = 2.0;
        double yIntercept1 = (110 - slope*10);
        for (int x = 10; x < 100; x++) {
            /*
            (y1 - y0)/(x1 - x0) = slope
            y1 - y0 = slope*(x1 - x0);
            y1 = y0 + slope*(x1 - x0);
            y1 = (y0 - slope*x0) + slope*x1
            */
            double dy = slope*(double)x;
            int y = (int)Math.round(yIntercept1 + dy);
            
            if (sr.nextBoolean()) {
                int xr = sr.nextInt(d);
                int yr = sr.nextInt(d);
                if (sr.nextBoolean()) {
                    xr *= -1;
                }
                if (sr.nextBoolean()) {
                    yr *= -1;
                }
                x += xr;
                y += yr;
            }
            
            dxdy.add(x, y);
        }
        
        LinearRegression instance = new LinearRegression();
        instance.plotTheLinearRegression(dxdy.getX(), dxdy.getY());
        
        float[] yInterceptAndSlope = 
            instance.calculateTheilSenEstimatorParams(dxdy.getX(), dxdy.getY());
        
        double expectedYIntercept = yIntercept1;
        double expectedSlope = slope;
        
        assertTrue(Math.abs(yInterceptAndSlope[0] - expectedYIntercept) < Math.sqrt(d));
        
        assertTrue(Math.abs(yInterceptAndSlope[1] - expectedSlope) < 0.2);
        
    }
    
    @Test
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
    
}
