package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION;
import algorithms.curves.GEVYFit;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.util.PolygonAndPointPlotter;
import java.security.SecureRandom;
import java.util.logging.Logger;

/**
 unit test specifically for examining the distributions of 4 extreme examples of points:
    -- a single group with no background points and the clustering is not a
       function of distance from centers of clusters
    -- a single group w/ no background points and the points are radially clustered toward center
       of cluster
    -- a few groups of varying size with no background points and no radial clustering
    -- background points only without groups
 
 The goals are:
    improve the histograms so that the histograms are always well formed and 
    if possible, distinguishable from one another.
 
    provide regression tests for the goal to easily verify that future improvements in the project
    form histograms correctly.
 
 Dependencies on project components:
    It uses findVoids() to create the array of densities for all points.
    
 * @author nichole
 */
public class DistributionsTest extends BaseTwoPointTest {

    boolean debug = true;

    boolean writeToTmpData = false;

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public void test_Find_Clusters_Stats() throws Exception {

        log.info("test_Find_Clusters_Stats()");

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        PolygonAndPointPlotter plotter2 = new PolygonAndPointPlotter();
        TwoPointCorrelationPlotter plotter3 = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);
        
        long seed = System.currentTimeMillis();
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        seed = 1386991714037l;
        sr.setSeed(seed);
        log.info("SEED=" + seed);

        int nSwitches = 5;

        int nIterPerBackground = 1;

        int m = nIterPerBackground*nSwitches;

        DoubleAxisIndexer indexer = null;

        int count = 0;
        
        int nClusters = 3;

        for (int i = 0; i < nSwitches; i++) {
            
            for (int ii = 0; ii < nIterPerBackground; ii++) {

                float xdiff = xmax - xmin;
                float ydiff = ymax - ymin;
                
                switch(i) {
                    case 0: {
                        // create a single group of uniformly distributed points that does not fill the boundaries
                        // and does not have any background points outside of group.
                        //createIndexerWithRandomPoints(SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
                        //   int[] nClusters, int nBackgroundPoints, CLUSTER_SEPARATION clusterSeparation) {
                        indexer = createIndexerWithRandomPoints(sr, 
                            (float)(xmin + 0.25*xdiff), (float)(xmax - 0.25*xdiff), 
                            (float)(ymin + 0.25*ydiff), (float)(ymax - 0.25*ydiff),
                            new int[]{300}, 0, CLUSTER_SEPARATION.MODERATE);
                        break;
                    }
                    case 1: {
                        // create a single group of radially distributed points that does not fill the boundaries
                        // and does not have any background points outside of group.
                        float maximumRadius = 0.18f * xdiff;
                        indexer = createIndexerWithRandomPointsAroundCenterWithDSquared(
                            sr, 300, xmin, xmax, ymin, ymax, maximumRadius);
                        break;
                    }
                    case 2: {
                        // 2 groups of uniformly filled random points
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
                        // background with 3 groups in it
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            nClusters, 200, 300, 1f);
                        break;
                    }
                                        
                    default:
                        break;
                }

                log.info(" " + i + ":" + ii + " (" + indexer.nXY + " points) ... ");

                if (writeToTmpData) {
                    writeIndexerToTmpData(indexer, count);
                }

                String plotLabel = "";
                plotter.addPlot(xmin, xmax, ymin, ymax, indexer.getX(), indexer.getY(), null, null, " " + i);
                plotter.writeFile();
                
                
                TwoPointCorrelation twoPtC = new TwoPointCorrelation(
                    generator.x, generator.y, generator.xErrors, generator.yErrors, generator.x.length);

                twoPtC.setDebug(true);                
              
                twoPtC.logPerformanceMetrics();
                
                //twoPtC.setBackground(0.2f, 0.1f);
                twoPtC.calculateBackground();
                
                TwoPointVoidStats stats = (TwoPointVoidStats)twoPtC.backgroundStats;
                
                if (stats != null) {
                    
                    float[] densities = stats.voidFinder.getTwoPointDensities();
                    float[] densityErrors = stats.voidFinder.getTwoPointDensityErrors();
                    
                   
                    HistogramHolder simpleHistogram = 
                        Histogram.calculateSturgesHistogramRemoveZeroTail(densities, densityErrors);
                    plotter2.addPlot(
                        simpleHistogram.getXHist(), simpleHistogram.getYHistFloat(),
                        simpleHistogram.getXErrors(), simpleHistogram.getYErrors(), null, null,
                        " " + i );
                    
                    HistogramHolder histogram = stats.statsHistogram;
        
                    GEVYFit bestFit = stats.bestFit;
                    if (bestFit != null) {
                        // label needs:  x10, peak,  mean/peak, median/mean and x80/median
                        plotLabel = String.format(
                            "*(%d %d) best k=%.4f sigma=%.4f mu=%.4f chiSqSum=%.6f chst=%.1f",
                            i, ii, bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), bestFit.getChiSqSum(), bestFit.getChiSqStatistic()
                        );
                        if (debug) {
                            log.info(plotLabel + " findVoid sampling=" + stats.getSampling().name());
                        }
                        plotter2.addPlot(histogram.getXHist(), histogram.getYHistFloat(),
                            histogram.getXErrors(), histogram.getYErrors(), 
                            bestFit.getOriginalScaleX(), bestFit.getOriginalScaleYFit(),
                            plotLabel);
                    } else {
                        plotter2.addPlot(histogram.getXHist(), histogram.getYHistFloat(),
                            histogram.getXErrors(), histogram.getYErrors(), null, null,
                            "*" + i );
                    }
                    
                    String filePath = plotter2.writeFile2();
                    System.out.println("filePath=" + filePath);
                   
                }
                
                twoPtC.findClusters();
                
                twoPtC.calculateHullsOfClusters();

                plotter3.addPlot(twoPtC, plotLabel);
                String filePath = plotter3.writeFile3();
                System.out.println("filePath=" + filePath);
            }
        }

        log.info("\n start computing stats for all sets");

        count = 0;

        log.info("SEED=" + seed);
    }
}
