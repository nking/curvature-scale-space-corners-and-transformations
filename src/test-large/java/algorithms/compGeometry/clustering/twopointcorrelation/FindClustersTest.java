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

    boolean debug = false;

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

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);

        SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        srr.setSeed(System.currentTimeMillis());
        long seed = srr.nextLong();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        //sr.setSeed(seed);
        sr.setSeed(7202122544352439191l);

        log.info("SEED=" + seed);

        // a long running test to calculcate and print the stats of fits
        //  for sparse, moderate, and densely populated backgrounds,
        //  all with the same number of clusters and cluster points, though
        //  randomly distributed.

        int nSwitches = 4;

        int nIterPerBackground = 10;

        int m = nIterPerBackground*nSwitches;

        float[] means = new float[m];
        float[] medians = new float[m];
        float[] peaks = new float[m];
        float[] x05s = new float[m];
        float[] x10s = new float[m];
        float[] x80s = new float[m];
        float[] x95s = new float[m];
        float[] chiSqStats = new float[m];

        AxisIndexer indexer = null;

        int count = 0;

        for (int i = 0; i < nSwitches; i++) {

            /*if (i != 3) {
                continue;
            }*/
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
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                             3, 150, 300, 20.0f);
                        break;
                    default:
                        // 100*100
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 30, 60, 100.0f);
                        //SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
                        //int numberOfClusters, int minimumNumberOfPointsPerCluster, int maximumNumberOfPointsPerCluster,
                        //float backgroundPointFractionToClusters
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

                twoPtC.setDebug(true);

                //twoPtC.logPerformanceMetrics();
                //twoPtC.calculateBackground();
                twoPtC.findClusters();

                TwoPointVoidStats stats = (TwoPointVoidStats)twoPtC.backgroundStats;
                HistogramHolder histogram = stats.statsHistogram;

                String plotLabel = null;

                GEVYFit bestFit = stats.bestFit;
                if (bestFit != null) {
                    float mean = bestFit.getXMean();
                    float median = bestFit.getXMedian();
                    float peak = bestFit.getXPeak();
                    float x05 = bestFit.getX05Percent();
                    float x10 = bestFit.getX10Percent();
                    float x80 = bestFit.getX80Percent();
                    float x95 = bestFit.getX95Percent();

                    means[count] = mean;
                    medians[count] = median;
                    peaks[count] = peak;
                    x05s[count] = x05;
                    x10s[count] = x10;
                    x80s[count] = x80;
                    x95s[count] = x95;

                    // === stats for plot labels =====
                    float meanDivPeak = mean/peak;
                    float medianDivMean = median/mean;
                    float x80DivMedian = x80/median;
                    chiSqStats[count] = bestFit.getChiSqStatistic();
                    // label needs:  x10, peak,  mean/peak, median/mean and x80/median
                    plotLabel = String.format(
                        "  (%d %d) x10=%.4f peak=%.4f av/peak=%.2f med/av=%.2f chst=%.1f",
                        i, ii, x10, peak, meanDivPeak, medianDivMean, chiSqStats[count]
                    );
                    if (debug) {
                        log.info(plotLabel + " findVoid sampling=" + stats.getSampling().name());
                    }
                }

                
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
        float[] meanDivPeaks = new float[nSwitches];
        float[] medianDivMeans = new float[nSwitches];
        float[] x80DivMedians = new float[nSwitches];
        float[] x80DivMeans = new float[nSwitches];

        float[] meanDivPeaksSD = new float[nSwitches];
        float[] medianDivMeansSD = new float[nSwitches];
        float[] x80DivMediansSD = new float[nSwitches];
        float[] x80DivMeansSD = new float[nSwitches];

        for (int i = 0; i < nSwitches; i++) {

            float meanDivPeakSum = 0;
            float medianDivMeanSum = 0;
            float x80DivMedianSum = 0;
            float x80DivMeanSum = 0;

            log.info("[" + i + "]");

            for (int j = 0; j < nIterPerBackground; j++) {

                int n = i*nIterPerBackground + j;

                float peak = peaks[n];
                float mean = means[n];
                float median = medians[n];

                float x05 = x05s[n];
                float x10 = x10s[n];
                float x80 = x80s[n];
                float x95 = x95s[n];

                float meanDivPeak = mean/peak;
                float medianDivMean = median/mean;
                float x80DivMedian = x80/median;
                float x80DivMean = x80/mean;

                String line = String.format(
                    "   (%d) peak=%.4f mean=%.4f median=%.4f x05=%.4f x10=%.4f x80=%.4f x95=%.4f mean/peak=%.2f median/mean=%.2f x80/median=%.2f x80/mean=%.2f chist=%.1f",
                    j, peak, mean, median, x05, x10, x80, x95, meanDivPeak, medianDivMean, x80DivMedian, x80DivMean, chiSqStats[n]
                );
                //log.info(line);
                log.info(line);

                meanDivPeakSum += meanDivPeak;
                medianDivMeanSum += medianDivMean;
                x80DivMedianSum += x80DivMedian;
                x80DivMeanSum += x80DivMean;
            }
            meanDivPeakSum /= nIterPerBackground;
            medianDivMeanSum /= nIterPerBackground;
            x80DivMedianSum /= nIterPerBackground;
            x80DivMeanSum /= nIterPerBackground;

            meanDivPeaks[i] = meanDivPeakSum;
            medianDivMeans[i] = medianDivMeanSum;
            x80DivMedians[i] = x80DivMedianSum;
            x80DivMeans[i] = x80DivMeanSum;

            double meanDivPeakSumSD = 0;
            double medianDivMeanSumSD = 0;
            double x80DivMedianSumSD = 0;
            double x80DivMeanSumSD = 0;

            for (int j = 0; j < nIterPerBackground; j++) {

                float peak = peaks[i*nIterPerBackground + j];

                float mean = means[i*nIterPerBackground + j];
                float median = medians[i*nIterPerBackground + j];

                float x05 = x05s[i*nIterPerBackground + j];
                float x10 = x10s[i*nIterPerBackground + j];
                float x80 = x80s[i*nIterPerBackground + j];
                float x95 = x95s[i*nIterPerBackground + j];

                double meanDivPeak = mean/peak;
                double medianDivMean = median/mean;
                double x80DivMedian = x80/median;
                double x80DivMean = x80/mean;

                meanDivPeakSumSD +=  Math.pow((meanDivPeak - meanDivPeakSum), 2);
                medianDivMeanSumSD += Math.pow((medianDivMean - medianDivMeanSum), 2);
                x80DivMedianSumSD += Math.pow((x80DivMedian - x80DivMedianSum), 2);
                x80DivMeanSumSD += Math.pow((x80DivMean -x80DivMeanSum), 2);
            }

            meanDivPeakSumSD = (float) Math.sqrt(meanDivPeakSumSD/(nIterPerBackground - 1.0f));
            medianDivMeanSumSD = (float) Math.sqrt(medianDivMeanSumSD/(nIterPerBackground - 1.0f));
            x80DivMedianSumSD = (float) Math.sqrt(x80DivMedianSumSD/(nIterPerBackground - 1.0f));
            x80DivMeanSumSD = (float) Math.sqrt(x80DivMeanSumSD/(nIterPerBackground - 1.0f));

            meanDivPeaksSD[i] = (float)meanDivPeakSumSD;
            medianDivMeansSD[i] = (float)medianDivMeanSumSD;
            x80DivMediansSD[i] = (float)x80DivMedianSumSD;
            x80DivMeansSD[i] = (float)x80DivMeanSumSD;
        }

        log.info("Final stats:");
        for (int i = 0; i < nSwitches; i++) {

            String line = String.format(
                "    mean/peak=%.2f +- %.4f   median/mean=%.2f +- %.4f    x80/median=%.2f +- %.4f    x80/median=%.2f +- %.4f",
                meanDivPeaks[i], meanDivPeaksSD[i],
                medianDivMeans[i], medianDivMeansSD[i],
                x80DivMedians[i], x80DivMediansSD[i],
                x80DivMeans[i], x80DivMeansSD[i]
            );
            //log.info(line);
            log.info(line);
        }

        log.info("SEED=" + seed);
    }

}
