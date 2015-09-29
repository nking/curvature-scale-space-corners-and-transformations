package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.util.PolygonAndPointPlotter;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.logging.Logger;

/**
 * 
 * @author nichole
 */
public class ConvexHullTest extends BaseTwoPointTest {

    boolean debug = true;

    /**
     *
     */
    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    /**
     *
     * @throws Exception
     */
    public void test0() throws Exception {

        log.info("test0()");

        // set of points with small subtended angles between them to assert
        // that the hull does not create slivers or spikes when the angle
        // tolerance for equality is too large of a number
        
        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(
            xmin, xmax, ymin, ymax);
        
        int nIter = 1;
        
        for (int i = 0; i < nIter; i++) {

            float[] xg = new float[]{55.830814f, 60.528652f, 61.03494f, 
                62.117302f, 62.33406f, 64.04489f, 64.83354f, 65.33296f, 
                66.928055f, 67.72369f, 69.087234f, 70.32792f, 70.43887f, 
                70.48301f, 71.225716f, 71.89872f, 72.91253f, 73.84599f, 
                67.854416f, 74.58776f, 76.92973f, 77.28116f, 77.92835f, 
                78.33889f, 74.51303f, 74.6322f, 76.50901f, 80.28216f, 
                86.431435f, 85.49746f, 91.74197f, 93.50279f, 71.42536f, 
                62.20092f, 64.59583f, 65.21205f, 65.21965f, 66.838844f, 
                57.57771f, 61.85385f, 59.28523f, 60.204033f, 67.35809f};
            
            float[] yg = new float[]{248.67566f, 248.40112f, 244.66013f, 
                256.1919f, 256.76428f, 239.83064f, 249.93317f, 258.84848f, 
                247.46191f, 253.68329f, 246.47108f, 253.5513f, 251.97011f, 
                254.19048f, 255.46196f, 258.77994f, 257.9417f, 244.37459f, 
                239.33615f, 241.79202f, 251.23628f, 242.14833f, 246.61655f, 
                252.84854f, 259.9717f, 261.69385f, 255.72478f, 259.24185f, 
                258.03833f, 263.31223f, 263.60165f, 260.29755f, 262.81522f, 
                262.80103f, 268.09784f, 262.6013f, 264.82138f, 262.13403f, 
                266.63126f, 271.03107f, 233.90057f, 237.91388f, 232.07506f};
            
            AxisIndexer indexer = new AxisIndexer();
            indexer.sortAndIndexX(xg, yg, xg.length);
            
            GrahamScan gs = new GrahamScan();
            
            gs.computeHull(indexer.getX(), indexer.getY());
            
            plotter.addPlot(indexer.getX(), indexer.getY(), gs.getXHull(), 
                gs.getYHull(), "");
            
            plotter.writeFile();
            
            assertNotNull(gs.getXHull());
            assertNotNull(gs.getYHull());
            assertTrue(gs.getXHull().length == 9);
            assertTrue(gs.getYHull().length == 9);
        }
        
        log.info("\n start computing stats for all sets");
    }
    
    /**
     *
     * @throws Exception
     */
    public void test1() throws Exception {

        log.info("test1()");
        
        // testing that random distributions don't produce exceptions
        
        int niter = 100;
        
        int xMin = 100;
        int xMax = 1000;

        int yMin = 100;
        int yMax = 1000;
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(
            xMin, xMax, yMin, yMax);

        RunStats[] stats = runAlgorithms(niter,
            xMin, xMax, yMin, yMax, plotter);

        StringBuffer summary = new StringBuffer();

        for (RunStats stat : stats) {

            stat.calculateFirstMomentRuntimeStats();

            summary.append(stat.getSummary()).append("\n");
        }

        log.info(summary.toString());
        
        String plotFilePath = plotter.writeFile();
        
        log.info("plot file = " + plotFilePath);
    }
    
    /**
     *
     * @param nIterations
     * @param xMin
     * @param xMax
     * @param yMin
     * @param yMax
     * @param plotter
     * @return
     * @throws NoSuchAlgorithmException
     * @throws IOException
     */
    public RunStats[] runAlgorithms(int nIterations,
        float xMin, float xMax, float yMin, float yMax, 
        PolygonAndPointPlotter plotter) throws 
        NoSuchAlgorithmException, IOException, GrahamScanTooFewPointsException {

        RunStats grahamStats = new RunStats(nIterations, "GrahamScan");

        for (int i = 0; i < nIterations; i++) {

            AxisIndexer indexer = createIndexerWithRandomPoints(
                xMin, xMax, yMin, yMax);
                        
            runGrahamScan(grahamStats, indexer.getX(), indexer.getY(), plotter);
        }

        return new RunStats[]{grahamStats};
    }

    /**
     *
     * @param runStats
     * @param x
     * @param y
     * @param plotter
     * @throws GrahamScanTooFewPointsException
     */
    public void runGrahamScan(RunStats runStats, float[] x, float[] y,
        PolygonAndPointPlotter plotter) 
        throws GrahamScanTooFewPointsException {

        runStats.recordSystemStats();

        long startTime = System.nanoTime();

        GrahamScan scan = new GrahamScan();

        scan.computeHull(x, y);

        double runtime = (System.nanoTime() - startTime)/1000000000.;

        runStats.addRuntime(runtime);
        
        assertNotNull(scan.getXHull());
        assertNotNull(scan.getYHull());
        assertTrue(scan.getXHull().length > 0);
        assertTrue(scan.getYHull().length > 0);
        
        plotter.addPlot(x, y, scan.getXHull(), scan.getYHull(), "");
        
    }
}
