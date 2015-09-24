package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION;
import algorithms.curves.GEVYFit;
import algorithms.imageProcessing.DistanceTransform;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.security.SecureRandom;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

/**
 * @author nichole
 */
public class FindClustersTest extends BaseTwoPointTest {

    boolean debug = false;

    boolean writeToTmpData = false;

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public void testFindClustersStats() throws Exception {

        log.info("testFindClustersStats()");

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

        //seed = 1387775326745l;

        //sr.setSeed(-2384802679227907254l);
        log.info("SEED=" + seed);

        // a long running test to calculate and print the stats of fits
        //  for sparse, moderate, and densely populated backgrounds,
        //  all with the same number of clusters and cluster points, though
        //  randomly distributed.

        int nSwitches = 3;

        int nIterPerBackground = 3;

        AxisIndexer indexer = null;

        int count = 0;

        for (int ii = 0; ii < nIterPerBackground; ii++) {
            for (int i = 0; i < nSwitches; i++) {
                try {
                    switch(i) {
                        case 0:
                            //~100
                            indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                                //10, 100, 110, 100.0f);
                                3, 33, 33, 0.1f);
                            break;
                        case 1:
                            //~1000
                            indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                                3, 33, 33, 10f);
                            break;
                        case 2:
                            // 100*100
                            indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                                3, 30, 60, 100.0f);
                            break;
                        case 3: {
                            int[] clusterNumbers = new int[]{1000, 300, 100};
                            int nBackgroundPoints = 10000;
                            CLUSTER_SEPARATION clusterSeparation = CLUSTER_SEPARATION.LARGE;

                            indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                                clusterNumbers, nBackgroundPoints, clusterSeparation);
                            break;
                        }
                        default: {
                            int[] clusterNumbers = new int[]{1000, 300, 100};
                            int nBackgroundPoints = 100000;
                            CLUSTER_SEPARATION clusterSeparation = CLUSTER_SEPARATION.LARGE;

                            indexer = generator.createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                                clusterNumbers, nBackgroundPoints, clusterSeparation);
                            break;
                        }
                    }
                    //TODO: consider reverting this code... and adding the ability to plot to the 
                    // performance metrics
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

                    long t0 = System.currentTimeMillis();


                    Set<PairInt> points = new HashSet<PairInt>();
                    for (int k = 0; k < indexer.getNXY(); ++k) {
                        PairInt p = new PairInt(Math.round(indexer.getX()[k]), 
                            Math.round(indexer.getY()[k]));
                        points.add(p);
                    }

                    TwoPointCorrelation twoPtC = new TwoPointCorrelation(indexer);

                    //twoPtC.setDebug(true);
                    
                    if (i == 4) {
                        //twoPtC.setAllowRefinement();
                    }
             
                    //twoPtC.logPerformanceMetrics();
                    //twoPtC.useFindMethodForDataWithoutBackgroundPoints();
                    twoPtC.calculateBackground();
                    //twoPtC.setBackground(0.03f, 0.01f);
                    
                    twoPtC.findClusters();

                    long t1 = (System.currentTimeMillis() - t0)/1000;
        
                    log.info(i + ")  ====> total RT(sec) = " + t1 + " nPoints=" + indexer.getNumberOfPoints());
        
                    
                    String plotLabel = "";

                    TwoPointVoidStats stats = (TwoPointVoidStats)twoPtC.backgroundStats;
                    
                    if (twoPtC.backgroundStats != null && 
                        (twoPtC.backgroundStats instanceof TwoPointVoidStats)) {
                        
                        HistogramHolder histogram = stats.statsHistogram;

                        GEVYFit bestFit = stats.bestFit;
                        if (bestFit != null) {
                            // label needs:  x10, peak,  mean/peak, median/mean and x80/median
                            plotLabel = String.format(
                                "(%d %d) k=%.4f s=%.4f m=%.4f chSq=%.6f chst=%.1f",
                                i, ii, bestFit.getK(), bestFit.getSigma(), 
                                bestFit.getMu(), bestFit.getChiSqSum(), 
                                bestFit.getChiSqStatistic()
                            );
                            if (debug) {
                                log.info(plotLabel + " findVoid sampling=" 
                                + stats.getSampling().name());
                            }
                        }
                    }

                    plotter.addPlot(twoPtC, plotLabel);
                    plotter.writeFile();                    
                    
                    // ----- a look at the distance transform ----
                    int[] minMaxXY = MiscMath.findMinMaxXY(points);
                    int w = minMaxXY[1] + 1;
                    int h = minMaxXY[3] + 1;
                    DistanceTransform dtr = new DistanceTransform();
                    int[][] dt = dtr.applyMeijsterEtAl(points, w, h);
                    float[] values = new float[dt.length*dt[0].length];
                    int count2 = 0;
                    for (int i0 = 0; i0 < dt.length; ++i0) {
                        for (int j0 = 0; j0 < dt[0].length; ++j0) {
                            int v = dt[i0][j0];
                            values[count2] = (float)(1./Math.sqrt(v));
                            count2++;
                        }
                    }
                    float xl = 1.5f;
                    if (i == 0) {
                        xl = 0.05f;
                    }
                    int nb = 20;
                    float[] vErrors = Errors.populateYErrorsBySqrt(values);
                    HistogramHolder hist = Histogram.createSimpleHistogram(
                        0, xl, nb, values, vErrors);
                    hist.plotHistogram("clstr", "_cluster_" + i + "_" + ii);

                } catch(Throwable e) {
                    log.severe(e.getMessage());
                }

                count++;
            }
        }

        log.info("\n start computing stats for all sets");

        count = 0;

        log.info("SEED=" + seed);
    }
}
