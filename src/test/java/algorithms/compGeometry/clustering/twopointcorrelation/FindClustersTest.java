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

        SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        srr.setSeed(System.currentTimeMillis());
        long seed = srr.nextLong();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(seed);
        //sr.setSeed(7202122544352439191l);

        log.info("SEED=" + seed);

        // a long running test to calculcate and print the stats of fits
        //  for sparse, moderate, and densely populated backgrounds,
        //  all with the same number of clusters and cluster points, though
        //  randomly distributed.

        int nSwitches = 3;

        int nIterPerBackground = 3;

        int m = nIterPerBackground*nSwitches;

        DoubleAxisIndexer indexer = null;

        int count = 0;

        for (int i = 0; i < nSwitches; i++) {

            for (int ii = 0; ii < nIterPerBackground; ii++) {

                switch(i) {
                    case 0:
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 30, 60, /*100.0f*/ 0.1f);
                        break;
                    case 1:
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 30, 60, 1f);
                        break;
                    case 2:
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 30, 60, 10.0f);
                        break;
                    case 3:
                        // 100*100
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 30, 60, 100.0f);
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

                System.out.println(" " + count + " (" + indexer.nXY + " points) ... ");

                TwoPointCorrelation twoPtC = new TwoPointCorrelation(
                    generator.x, generator.y,
                    generator.xErrors, generator.yErrors, generator.x.length);

                twoPtC.setDebug(true);

                twoPtC.logPerformanceMetrics();
                twoPtC.calculateBackground();
                twoPtC.findClusters();

                TwoPointVoidStats stats = (TwoPointVoidStats)twoPtC.backgroundStats;
                HistogramHolder histogram = stats.statsHistogram;

                String plotLabel = null;

                System.out.print(" storing statistics ");
                GEVYFit bestFit = stats.bestFit;
                if (bestFit != null) {
                    float mean = bestFit.getXMean();
                    float median = bestFit.getXMedian();
                    float peak = bestFit.getXPeak();
                    float x05 = bestFit.getX05Percent();
                    float x10 = bestFit.getX10Percent();
                    float x80 = bestFit.getX80Percent();
                    float x95 = bestFit.getX95Percent();

                    // === stats for plot labels =====
                    float meanDivPeak = mean/peak;
                    float medianDivMean = median/mean;
                    float x80DivMedian = x80/median;
                    // label needs:  x10, peak,  mean/peak, median/mean and x80/median
                    plotLabel = String.format(
                        "  (%d %d) x10=%.4f peak=%.4f av/peak=%.2f med/av=%.2f chst=%.1f",
                        i, ii, x10, peak, meanDivPeak, medianDivMean, bestFit.getChiSqStatistic()
                    );
                    if (debug) {
                        System.out.println(plotLabel + " findVoid sampling=" + stats.getSampling().name());
                    }
                }

                twoPtC.calculateHullsOfClusters();

                plotter.addPlot(twoPtC, plotLabel);
                plotter.writeFile();

                if (debug) {

                    if (twoPtC.backgroundStats instanceof TwoPointVoidStats) {

                        HistogramHolder hist = ((TwoPointVoidStats)twoPtC.backgroundStats).statsHistogram;

                        log.info("\n   (" + i + " " + ii + ")");
                        StringBuffer s0 = new StringBuffer("  x = new float[]{");
                        StringBuffer s1 = new StringBuffer("  y = new float[]{");
                        for (int iii = 0; iii < hist.getXHist().length; iii++) {
                            if (iii > 0) {
                                s0.append(", ");
                                s1.append(", ");
                            }
                            s0.append(hist.getXHist()[iii]).append("f");
                            s1.append(hist.getYHistFloat()[iii]).append("f");
                        }
                        s0.append("};");
                        s1.append("};");
                        log.info(s0.toString());
                        log.info(s1.toString());
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
