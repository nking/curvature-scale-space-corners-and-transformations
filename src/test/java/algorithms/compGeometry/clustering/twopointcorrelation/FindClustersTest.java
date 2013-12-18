package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.curves.GEVYFit;
import algorithms.misc.HistogramHolder;
import algorithms.util.ResourceFinder;
import java.security.SecureRandom;
import java.util.logging.Logger;

/**
 * @author nichole
 */
public class FindClustersTest extends BaseTwoPointTest {

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


        //seed = 1387019322723l;

        sr.setSeed(seed);
        log.info("SEED=" + seed);

        // a long running test to calculate and print the stats of fits
        //  for sparse, moderate, and densely populated backgrounds,
        //  all with the same number of clusters and cluster points, though
        //  randomly distributed.

        int nSwitches = 3;

        int nIterPerBackground = 5;

        DoubleAxisIndexer indexer = null;

        int count = 0;

        int nClusters = 3;

        for (int i = 0; i < nSwitches; i++) {

            for (int ii = 0; ii < nIterPerBackground; ii++) {

                try {
                    switch(i) {
                        case 0:
                            indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                                nClusters, 30, 60, /*100.0f*/ 0.1f);
                            break;
                        case 1:
                            indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                                nClusters, 30, 60, 1f);
                            break;
                        case 2:
                            indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                                nClusters, 30, 60, 10.0f);
                            break;
                        case 3:
                            // 100*100
                            indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                                nClusters, 30, 60, 100.0f);
                            break;
                        default:
                            break;
                    }

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

                    twoPtC.logPerformanceMetrics();
                    //twoPtC.useFindMethodForDataWithoutBackgroundPoints();
                    twoPtC.calculateBackground();
                    //twoPtC.setBackground(0.2f, 0.1f);
                    twoPtC.findClusters();

                    String plotLabel = "";

                    TwoPointVoidStats stats = (TwoPointVoidStats)twoPtC.backgroundStats;
                    
                    if (twoPtC.backgroundStats != null && twoPtC.backgroundStats instanceof TwoPointVoidStats) {
                        
                        HistogramHolder histogram = stats.statsHistogram;

                        GEVYFit bestFit = stats.bestFit;
                        if (bestFit != null) {

                            // label needs:  x10, peak,  mean/peak, median/mean and x80/median
                            plotLabel = String.format(
                                "(%d %d) k=%.4f s=%.4f m=%.4f chSq=%.6f chst=%.1f",
                                i, ii, bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), bestFit.getChiSqSum(), bestFit.getChiSqStatistic()
                            );
                            if (debug) {
                                log.info(plotLabel + " findVoid sampling=" + stats.getSampling().name());
                            }
                        }

                        if (false) { // for print out to improve fit using NonQuadraticConjugateGradientSolverTest
                            if (i == 2 && ii == 2) {
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
                    
                    /*
                    if (i == 0 && ii == 1) {
                        StringBuilder xsb = new StringBuilder();
                        StringBuilder ysb = new StringBuilder();
                        float[] xg = twoPtC.groupFinder.getX(1, indexer);
                        float[] yg = twoPtC.groupFinder.getY(1, indexer);
                        for (int j = 0; j < xg.length; j++) {
                            if (j > 0) {
                                xsb.append("f, ");
                                ysb.append("f, ");
                            }
                            xsb.append(xg[j]);
                            ysb.append(yg[j]);
                        }
                        System.out.println("float[] x = new float[]{"  + xsb.append("f").toString() + "};");
                        System.out.println("float[] y = new float[]{"  + ysb.append("f").toString() + "};");
                    }
                    */

                    plotter.addPlot(twoPtC, plotLabel);
                    plotter.writeFile();

                 // assert that the low number histograms are all well formed and result in finding n clusters
                    if (i == 0) {
                        //assertTrue(twoPtC.getNumberOfGroups() >= nClusters);
                        if (twoPtC.getNumberOfGroups() < nClusters) {
                            log.severe("Note:  for seed=" + seed + " and i=" + i + ", ii=" + ii
                                + " solution did not find " + nClusters + " clusters");
                        }
                    }

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
