package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.curves.GEVYFit;
import algorithms.util.PolygonAndPointPlotter;
import java.security.SecureRandom;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class FittingStatsTest extends BaseTwoPointTest {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());
    
    public void test_Fit_Stats() throws Exception {

        log.info("test_Fit_Stats()");

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        PolygonAndPointPlotter plotter2 = new PolygonAndPointPlotter();

        SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        srr.setSeed(System.currentTimeMillis());
        long seed = srr.nextLong();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(seed);
        //sr.setSeed(6236290033146721436l);

        // a long running test to calculcate and print the stats of fits
        //  for sparse, moderate, and densely populated backgrounds,
        //  all with the same number of clusters and cluster points, though
        //  randomly distributed.

        // ===== see tests in tests_larger for larger number of iterations ====
        int nIterPerBackground = 0;

        int m = nIterPerBackground*3;

        float[] means = new float[m];
        float[] medians = new float[m];
        float[] peaks = new float[m];
        float[] x10s = new float[m];
        float[] x80s = new float[m];
        float[] x95s = new float[m];

        DoubleAxisIndexer indexer = null;

        int count = 0;

        for (int i = 0; i < 3; i++) {

            for (int ii = 0; ii < nIterPerBackground; ii++) {

                switch(i) {
                    case 0:
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 30, 60, 0.1f);
                        break;
                    case 1:
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 30, 60, 1f);
                        break;
                    case 2:
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 30, 60, 10.0f);
                        break;
                    default:
                        break;
                }

                indexer.sortAndIndexXThenY(x, y, xErrors, yErrors, x.length);

                log.info(" " + count + "(" + indexer.nXY + " points) ... ");

                TwoPointVoidStats stats = new TwoPointVoidStats(indexer);
                stats.setDebug(false);
                stats.setUseCompleteSampling(false);
                stats.calc();

                stats.plotFit(plotter2);

                System.out.print(" storing statistics ");

                GEVYFit bestFit = stats.bestFit;

                if (bestFit != null) {
                    float[] xf = bestFit.getX();
                    float mean = xf[ bestFit.getxMeanIndex()];
                    float median = xf[ bestFit.getXMedianIndex()];
                    float peak = xf[ bestFit.getXPeakIndex()];
                    float x10 = xf[ bestFit.getX10PercentIndex()];
                    float x80 = xf[ bestFit.getX80PercentIndex()];
                    float x95 = xf[ bestFit.getX95PercentIndex()];

                    means[count] = mean;
                    medians[count] = median;
                    peaks[count] = peak;
                    x10s[count] = x10;
                    x80s[count] = x80;
                    x95s[count] = x95;
                }

                count++;
            }
        }

        log.info("\n start computing stats for all sets");

        count = 0;
        float[] meanDivPeaks = new float[3];
        float[] medianDivMeans = new float[3];
        float[] x80DivMedians = new float[3];
        float[] x80DivMeans = new float[3];

        float[] meanDivPeaksSD = new float[3];
        float[] medianDivMeansSD = new float[3];
        float[] x80DivMediansSD = new float[3];
        float[] x80DivMeansSD = new float[3];

        for (int i = 0; i < 3; i++) {

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

                float x10 = x10s[n];
                float x80 = x80s[n];
                float x95 = x95s[n];

                float meanDivPeak = mean/peak;
                float medianDivMean = median/mean;
                float x80DivMedian = x80/median;
                float x80DivMean = x80/mean;

                String line = String.format(
                    "   (%d) peak=%.4f mean=%.4f median=%.4f x10=%.4f x80=%.4f x95=%.4f mean/peak=%.2f median/mean=%.2f x80/median=%.2f x80/mean=%.2f",
                    j, peak, mean, median, x10, x80, x95, meanDivPeak, medianDivMean, x80DivMedian, x80DivMean
                );
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
        for (int i = 0; i < 3; i++) {

            String line = String.format(
                "    mean/peak=%.2f +- %.4f   median/mean=%.2f +- %.4f    x80/median=%.2f +- %.4f    x80/median=%.2f +- %.4f",
                meanDivPeaks[i], meanDivPeaksSD[i],
                medianDivMeans[i], medianDivMeansSD[i],
                x80DivMedians[i], x80DivMediansSD[i],
                x80DivMeans[i], x80DivMeansSD[i]
            );

            log.info(line);
        }
    }

}
