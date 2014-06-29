package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.clustering.twopointcorrelation.AxisIndexer;
import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator;
import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION;
import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.curves.GEVYFit;
import algorithms.misc.HistogramHolder;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.logging.Logger;

/**
 * 
 * @author nichole
 */
public class ConvexHullTest extends BaseTwoPointTest {

    boolean debug = true;

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public void testConvexHull() throws Exception {

        log.info("testConvexHull()");

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        //SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        //srr.setSeed(System.currentTimeMillis());
        //long seed = srr.nextLong();

        long seed = System.currentTimeMillis();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
     
        seed = 1387019322723l;
        sr.setSeed(seed);
        log.info("SEED=" + seed);
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);
        
        int nIter = 1;
        
        for (int i = 0; i < nIter; i++) {

            int[] clusterNumbers = new int[]{300, 100, 30};
            
            int nBackgroundPoints = 100;
                            
            CLUSTER_SEPARATION clusterSeparation = CLUSTER_SEPARATION.LARGE;
                            
            AxisIndexer indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                clusterNumbers, nBackgroundPoints, clusterSeparation);
            
            /*
            int count = 0;
            for (int ii = 0; ii < indexer.getNumberOfPoints(); ii++) {
                if (indexer.getY()[ii] > 90 && indexer.getX()[ii] < 100) {
                    count++;
                }
            }
            float[] xt = new float[count];
            float[] yt = new float[count];
            float[] xte = new float[count];
            float[] yte = new float[count];
            count = 0;
            for (int ii = 0; ii < indexer.getNumberOfPoints(); ii++) {
                if (indexer.getY()[ii] > 90 && indexer.getX()[ii] < 100) {
                    xt[count] = indexer.getX()[ii];
                    yt[count] = indexer.getY()[ii];
                    xte[count] = indexer.getXErrors()[ii];
                    yte[count] = indexer.getYErrors()[ii];
                    count++;
                }
            }
            AxisIndexer tmp = new AxisIndexer();
            tmp.sortAndIndexXThenY(xt, yt, xte, yte, count);
            indexer = tmp;
            */
            
            float[] xg = new float[]{55.830814f, 60.528652f, 61.03494f, 62.117302f, 62.33406f, 64.04489f, 64.83354f, 65.33296f, 66.928055f, 67.72369f, 69.087234f, 70.32792f, 70.43887f, 70.48301f, 71.225716f, 71.89872f, 72.91253f, 73.84599f, 67.854416f, 74.58776f, 76.92973f, 77.28116f, 77.92835f, 78.33889f, 74.51303f, 74.6322f, 76.50901f, 80.28216f, 86.431435f, 85.49746f, 91.74197f, 93.50279f, 71.42536f, 62.20092f, 64.59583f, 65.21205f, 65.21965f, 66.838844f, 57.57771f, 61.85385f, 59.28523f, 60.204033f, 67.35809f};
            
            float[] yg = new float[]{248.67566f, 248.40112f, 244.66013f, 256.1919f, 256.76428f, 239.83064f, 249.93317f, 258.84848f, 247.46191f, 253.68329f, 246.47108f, 253.5513f, 251.97011f, 254.19048f, 255.46196f, 258.77994f, 257.9417f, 244.37459f, 239.33615f, 241.79202f, 251.23628f, 242.14833f, 246.61655f, 252.84854f, 259.9717f, 261.69385f, 255.72478f, 259.24185f, 258.03833f, 263.31223f, 263.60165f, 260.29755f, 262.81522f, 262.80103f, 268.09784f, 262.6013f, 264.82138f, 262.13403f, 266.63126f, 271.03107f, 233.90057f, 237.91388f, 232.07506f};
            
            indexer = new AxisIndexer();
            indexer.sortAndIndexX(xg, yg, xg.length);
            
            GrahamScan gs = new GrahamScan();
            
            gs.computeHull(indexer.getX(), indexer.getY());
            
            plotter.addPlot(indexer.getX(), indexer.getY(), gs.getXHull(), gs.getYHull(), "");
            plotter.writeFile();
            
        }
        
        log.info("\n start computing stats for all sets");


        log.info("SEED=" + seed);
    }
}
