package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION;
import algorithms.curves.GEVYFit;
import algorithms.misc.HistogramHolder;
import algorithms.util.ResourceFinder;
import java.security.SecureRandom;
import java.util.logging.Logger;

/**
 * class to test the TwoPointCorrelation class on larger datasets with and without clusters.
 * 
 * @author nichole
 */
public class FindClusters2Test extends BaseTwoPointTest {

    boolean debug = true;

    boolean writeToTmpData = false;

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public void test_Find_Clusters_Stats() throws Exception {

        log.info("test_Find_Clusters_Stats()");

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);

        //SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        //srr.setSeed(System.currentTimeMillis());
        //long seed = srr.nextLong();

        long seed = System.currentTimeMillis();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
     
        //long seed = 1386750505246l;

        sr.setSeed(seed);
        log.info("SEED=" + seed);

        // a long running test to calculate and print the stats of fits
        //  for sparse, moderate, and densely populated backgrounds,
        //  all with the same number of clusters and cluster points, though
        //  randomly distributed.

        int nSwitches = 2;

        int nIterPerBackground = 1;

        int m = nIterPerBackground*nSwitches;

        DoubleAxisIndexer indexer = null;

        int count = 0;
        
        int nClusters = 3;

        for (int i = 0; i < nSwitches; i++) {

            for (int ii = 0; ii < nIterPerBackground; ii++) {                
                                
                switch(i) {
                    case 0:                        
                        // essentially white noise
                        
                        int minimumNumberOfPointsPerCluster = 30;
                        int maximumNumberOfPointsPerCluster = 60;
                        float backgroundPointFractionToClusters = 100.0f;
                        
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            nClusters, minimumNumberOfPointsPerCluster, maximumNumberOfPointsPerCluster, 
                            backgroundPointFractionToClusters);
                        
                        break;
                    case 1:
                        
                        int[] clusterNumbers = new int[]{1000, 300, 100};
                        int nBackgroundPoints = 10000;
                        CLUSTER_SEPARATION clusterSeparation = CLUSTER_SEPARATION.LARGE;
                        
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            clusterNumbers, nBackgroundPoints, clusterSeparation);
                            
                        break;
                   
                    default:
                        break;
                }

                indexer.sortAndIndexXThenY(generator.x, generator.y,
                    generator.xErrors, generator.yErrors, generator.x.length);

                log.info(" " + count + " (" + indexer.nXY + " points) ... ");


                if (writeToTmpData) {
                    // write to tmpdata if need to use in tests improve fits, histogram etc
                    String str = String.valueOf(count);
                    while (str.length() < 3) {
                        str = "0" + str;
                    }
                    String fileNamePostfix = "_clusters_" + str + ".dat";
                    String fileName = CreateClusterDataTest.indexerFileNamePrefix + fileNamePostfix;
                    String filePath = ResourceFinder.getAFilePathInTmpData(fileName);
                    CreateClusterDataTest.writeIndexer(filePath, indexer);
                }

                log.info(" " + count + " (" + indexer.nXY + " points) ... ");

                TwoPointCorrelation twoPtC = new TwoPointCorrelation(
                    generator.x, generator.y,
                    generator.xErrors, generator.yErrors, generator.x.length);

                //twoPtC.setDebug(true);
                
//twoPtC.setUseDownhillSimplexHistogramFitting();
              
                twoPtC.logPerformanceMetrics();
                twoPtC.calculateBackground();
                twoPtC.findClusters();

                TwoPointVoidStats stats = (TwoPointVoidStats)twoPtC.backgroundStats;
                HistogramHolder histogram = stats.statsHistogram;

                String plotLabel = null;

                GEVYFit bestFit = stats.bestFit;
                if (bestFit != null) {
                    
                    // label needs:  x10, peak,  mean/peak, median/mean and x80/median
                    plotLabel = String.format(
                        "  (%d %d) best k=%.4f sigma=%.4f mu=%.4f chiSqSum=%.6f chst=%.1f",
                        i, ii, bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), bestFit.getChiSqSum(), bestFit.getChiSqStatistic()
                    );
                    if (debug) {
                        log.info(plotLabel + " findVoid sampling=" + stats.getSampling().name());
                    }
                }
                
                if (false) { // for print out to improve fit using NonQuadraticConjugateGradientSolverTest
                    if (i == 0 && ii == 0) {
                        StringBuilder xsb = new StringBuilder();
                        StringBuilder ysb = new StringBuilder();
                        StringBuilder xesb = new StringBuilder();
                        StringBuilder yesb = new StringBuilder();
    
                        for (int z = 0; z < histogram.getYHist().length; z++) {
                            if (z > 0) {
                                xsb.append("f, ");
                                ysb.append("f, ");
                                xesb.append("f, ");
                                yesb.append("f, ");
                            }
                            xsb.append(histogram.getXHist()[z]);
                            ysb.append(histogram.getYHist()[z]);
                            xesb.append(histogram.getXErrors()[z]);
                            yesb.append(histogram.getYErrors()[z]);
                        }
                        System.out.println("float[] x = new float[]{"  + xsb.append("f").toString() + "};");
                        System.out.println("float[] y = new float[]{"  + ysb.append("f").toString() + "};");
                        System.out.println("float[] xe = new float[]{" + xesb.append("f").toString() + "};");
                        System.out.println("float[] ye = new float[]{" + yesb.append("f").toString() + "};");
                        int z = 1;
                    }
                }
                
                twoPtC.calculateHullsOfClusters();

                plotter.addPlot(twoPtC, plotLabel);
                //plotter.addPlotWithoutHull(twoPtC, plotLabel);
                plotter.writeFile();

                if (i == 0) {
                    //assertTrue(twoPtC.getNumberOfGroups() >= nClusters);
                    if (twoPtC.getNumberOfGroups() > 0) {
                        log.severe("Note:  for seed=" + seed + " and i=" + i + ", ii=" + ii 
                            + " solution should not find clusters: " + nClusters);
                    }
                }

                count++;
            }
        }

        log.info("\n start computing stats for all sets");

        count = 0;

        log.info("SEED=" + seed);
    }
}
