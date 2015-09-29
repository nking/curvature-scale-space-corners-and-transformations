package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION;
import algorithms.curves.GEVYFit;
import java.security.SecureRandom;
import java.util.logging.Logger;

/**
 unit test specifically for examining the distributions of 4 extreme examples of 
 points:
    -- a single group with no background points and the clustering is not a
       function of distance from centers of clusters
    -- a single group w/ no background points and the points are radially 
       clustered toward center of cluster
    -- a few groups of varying size with no background points and no radial 
       clustering
    -- background points only without groups
 
 The goals are:
    improve the histograms so that the histograms are always well formed and 
    if possible, distinguishable from one another.
 
    provide regression tests for the goal to easily verify that future 
    improvements in the project form histograms correctly.
     
    compare the GEV fitting algorithm results (non-quadratic conjugate gradient
    solver versus downhill simplex method).
    
 * @author nichole
 */
public class DistributionsTest extends BaseTwoPointTest {

    boolean debug = true;

    boolean writeToTmpData = false;

    /**
     *
     */
    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    /**
     *
     * @throws Exception
     */
    public void test_Find_Clusters_Stats() throws Exception {

        log.info("test_Find_Clusters_Stats()");
        
        // randomly generates point distributions to assert that exceptions
        // aren't thrown and creates plot to compare visually the fits
        // between the 2 fitting algorithms.

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;
        
        TwoPointCorrelationPlotter tpcPlotter = new TwoPointCorrelationPlotter(
            xmin, xmax, ymin, ymax);
        TwoPointCorrelationPlotter tpcPlotter3 = new TwoPointCorrelationPlotter(
            xmin, xmax, ymin, ymax);
        
        long seed = System.currentTimeMillis();
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        //seed = 1404074772766l;
        sr.setSeed(seed);
        log.info("SEED=" + seed);

        int nSwitches = 6;

        int nIterPerBackground = 1;

        int m = nIterPerBackground*nSwitches;

        AxisIndexer indexer = null;

        int count = 0;

        for (int i = 0; i < nSwitches; i++) {
            
            for (int ii = 0; ii < nIterPerBackground; ii++) {

                float xdiff = xmax - xmin;
                float ydiff = ymax - ymin;
                
                switch(i) {
                    case 0: {
                        // create a single group of uniformly distributed points 
                        // that does not fill the boundaries and does not have 
                        // any background points outside of group.
                        indexer = createIndexerWithRandomPoints(sr, 
                            (float)(xmin + 0.25*xdiff), (float)(xmax - 0.25*xdiff), 
                            (float)(ymin + 0.25*ydiff), (float)(ymax - 0.25*ydiff),
                            new int[]{300}, 0, CLUSTER_SEPARATION.MODERATE);
                        break;
                    }
                    case 1: {
                        // create a single group of radially distributed points 
                        // that does not fill the boundaries and does not have 
                        // any background points outside of group.
                        float maximumRadius = 0.18f * xdiff;
                        indexer = createIndexerWithRandomPointsAroundCenterWithDSquared(
                            sr, 300, xmin, xmax, ymin, ymax, maximumRadius);
                        break;
                    }
                    case 2: {
                        // 2 uniformly filled groups of random points
                        indexer = createIndexerWithRandomPoints(sr,
                            xmin, xmax, ymin, ymax,
                            3, 200, 300, 0.01f);
                        
                        break;
                    }
                    case 3: {
                        // background points only, no groups
                        //createIndexerWithRandomPoints(SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
                        //   int[] nClusters, int nBackgroundPoints, CLUSTER_SEPARATION clusterSeparation) {
                        indexer = createIndexerWithRandomPoints(sr, 
                            xmin, xmax, ymin, ymax,
                            new int[]{0}, 1000, CLUSTER_SEPARATION.MODERATE);
                        break;
                    }
                    case 4: {
                        // background with 6 groups in it
                        int nClusters = 6;
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, 
                            ymin, ymax, nClusters, 200, 300, 1f);
                        break;
                    }
                    case 5: {
                        // 7 groups of uniformly filled random points and a
                        // relatively high density background
                        int nClusters = 7;
                        indexer = createIndexerWithRandomPoints(sr,
                            xmin, xmax, ymin, ymax, nClusters, 200, 300, 3.f);
                        
                        break;
                    }
                                        
                    default:
                        break;
                }

                log.info(" " + i + ":" + ii + " (" + indexer.nXY + " points) ... ");

                if (writeToTmpData) {
                    writeIndexerToTmpData(indexer, count);
                }

                for (int iii = 0; iii <= 1; iii++) {
                    
                    TwoPointCorrelation twoPtC = new TwoPointCorrelation(indexer);
                    
                    if (iii == 1) {
                        twoPtC.setUseDownhillSimplexHistogramFitting();
                    }

                    twoPtC.logPerformanceMetrics();

                    //twoPtC.setThresholdFactorToThree();

                    //twoPtC.setBackground(0.2f, 0.1f);
                    twoPtC.calculateBackground();

                    twoPtC.findClusters();

                    String plotLabel = "";

                    TwoPointVoidStats stats = (TwoPointVoidStats)twoPtC.backgroundStats;
                    if (stats != null) {          
                        GEVYFit bestFit = stats.bestFit;
                        if (bestFit != null) {
                            plotLabel = String.format(
                                "*(%d %d) best k=%.4f sigma=%.4f mu=%.4f chiSqSum=%.6f chst=%.1f",
                                i, ii, bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), 
                                bestFit.getChiSqSum(), bestFit.getChiSqStatistic());
                        }
                    }

                    if (iii == 0) {
                        // plot for non-quadratic conjugate gradient solver
                        tpcPlotter.addPlot(twoPtC, plotLabel);
                        String filePath = tpcPlotter.writeFile();
                        log.fine("filePath=" + filePath);
                    } else {
                        // plot for the downhill simplex method
                        tpcPlotter3.addPlot(twoPtC, plotLabel);
                        String filePath = tpcPlotter3.writeFile3();
                        log.fine("filePath3=" + filePath);
                    }
                }
            }
        }

        log.info("\n start computing stats for all sets");

        count = 0;

        log.info("SEED=" + seed);
    }
}
