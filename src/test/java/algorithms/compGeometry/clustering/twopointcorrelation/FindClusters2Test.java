package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION;
import algorithms.curves.GEVYFit;
import algorithms.misc.HistogramHolder;
import algorithms.util.ResourceFinder;
import java.security.SecureRandom;
import java.util.Arrays;
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

        long seed = System.currentTimeMillis();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
     
        seed = 1386750505246l;

        sr.setSeed(seed);
        log.info("SEED=" + seed);

        // a long running test to calculate and print the stats of fits
        //  for sparse, moderate, and densely populated backgrounds,
        //  all with the same number of clusters and cluster points, though
        //  randomly distributed.

        int nSwitches = 3;

        int nIterPerBackground = 1;

        int m = nIterPerBackground * nSwitches;

        DoubleAxisIndexer indexer = null;

        int count = 0;
        
        int nClusters = 3;
        
        /* creating 3 cases of moderate background density:
         *    case 0:    90*100 = 9000 background points
         *    case 1:  the case 0 9000 background points + 3 clusters that have a range of points from 30 to 60 with large separation
         *    case 2:  the case 0 9000 background points + 3 clusters that have 1000, 300, 100 points           with large separation
         */
        
        int numberOfBackgroundPoints = 9000;
        
        CLUSTER_SEPARATION clusterSeparation = CLUSTER_SEPARATION.LARGE;
        
        for (int ii = 0; ii < nIterPerBackground; ii++) { 
            
            float[] xb = new float[numberOfBackgroundPoints];
            float[] yb = new float[numberOfBackgroundPoints];

            int xyStartOffset = 0;
            
            createRandomPointsInRectangle(sr, numberOfBackgroundPoints,
                xmin, xmax, ymin, ymax, xb, yb, xyStartOffset);
            
            float[] xbe = new float[numberOfBackgroundPoints];
            float[] ybe = new float[numberOfBackgroundPoints];
            for (int i = 0; i < numberOfBackgroundPoints; i++) {
                // simulate x error as a percent error of 0.03 for each bin
                xbe[i] = xb[i] * 0.03f;
                ybe[i] = (float) (Math.sqrt(yb[i]));
            }
            
            
            for (int i = 0; i < nSwitches; i++) {              
                                
                switch(i) {
                    
                    case 0: {
                        indexer = new DoubleAxisIndexer();
                        indexer.sortAndIndexXThenY(xb, yb, xbe, ybe, xbe.length);
                        break;
                    }
                    case 1: {
                        int[] clusterNumbers = new int[]{60, 45, 30};
                        
                        int tot = numberOfBackgroundPoints + 60 + 45 + 30;
                        
                        xb = Arrays.copyOf(xb, tot);
                        yb = Arrays.copyOf(yb, tot);
                        
                        float[] xbc = new float[clusterNumbers.length];
                        float[] ybc = new float[clusterNumbers.length];
                        
                        generator.createRandomClusters(sr, xmin, xmax, ymin, ymax,
                            clusterNumbers, clusterSeparation, xb, yb,xbc, ybc, numberOfBackgroundPoints);
                        
                        xbe = Arrays.copyOf(xbe, tot);
                        ybe = Arrays.copyOf(ybe, tot);
                        for (int j = numberOfBackgroundPoints; j < tot; j++) {
                            // simulate x error as a percent error of 0.03 for each bin
                            xbe[j] = xb[j] * 0.03f;
                            ybe[j] = (float) (Math.sqrt(yb[j]));
                        }
                        
                        indexer = new DoubleAxisIndexer();
                        indexer.sortAndIndexXThenY(xb, yb, xbe, ybe, xbe.length);
                        
                        break;
                    }
                    case 2: {
                        int[] clusterNumbers = new int[]{1000, 300, 100};
                        
                        int tot = numberOfBackgroundPoints + 1000 + 300 + 100;
                        
                        xb = Arrays.copyOf(xb, tot);
                        yb = Arrays.copyOf(yb, tot);
                        
                        float[] xbc = new float[clusterNumbers.length];
                        float[] ybc = new float[clusterNumbers.length];
                        
                        generator.createRandomClusters(sr, xmin, xmax, ymin, ymax,
                            clusterNumbers, clusterSeparation, xb, yb,xbc, ybc, numberOfBackgroundPoints);
                        
                        xbe = Arrays.copyOf(xbe, tot);
                        ybe = Arrays.copyOf(ybe, tot);
                        for (int j = numberOfBackgroundPoints; j < tot; j++) {
                            // simulate x error as a percent error of 0.03 for each bin
                            xbe[j] = xb[j] * 0.03f;
                            ybe[j] = (float) (Math.sqrt(yb[j]));
                        }
                        
                        indexer = new DoubleAxisIndexer();
                        indexer.sortAndIndexXThenY(xb, yb, xbe, ybe, xbe.length);
                        
                        break;
                    }

                    default:
                        break;
                }
               
                /* to zoom in to confirm density estimate visually:
                DoubleAxisIndexer tmp = new DoubleAxisIndexer();
                tmp.sortAndIndexXThenY(
                    Arrays.copyOf(indexer.getX(), 100), Arrays.copyOf(indexer.getY(), 100),
                    Arrays.copyOf(indexer.getXErrors(), 100), Arrays.copyOf(indexer.getYErrors(), 100),
                    100);
                indexer = tmp;
                */
                
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

                TwoPointCorrelation twoPtC = new TwoPointCorrelation(indexer);

                twoPtC.setDebug(true);
                
//twoPtC.setUseDownhillSimplexHistogramFitting();
              
                twoPtC.logPerformanceMetrics();
                //twoPtC.setBackground(0.5f, 0.02f);
                //twoPtC.setBackground(0.28f, 0.02f);
                twoPtC.calculateBackground();
                twoPtC.findClusters();

                String plotLabel = "";
                
                if (twoPtC.backgroundStats != null) {
                    TwoPointVoidStats stats = (TwoPointVoidStats)twoPtC.backgroundStats;
                    HistogramHolder histogram = stats.statsHistogram;
        
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
                        if (i == 1 && ii == 0) {
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
                }
                
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
