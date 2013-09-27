package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.misc.HistogramHolder;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;
import java.security.SecureRandom;
import java.util.logging.Logger;

/**
  estimates performance metrics for TwoPointCorrelation as
  the dependencies upon N of the amount of memory needed and the runtime.

  @author nichole
 */
public class PerformanceMetricsRunner {

    protected static void printSystemStats() {

        Logger log = Logger.getLogger(PerformanceMetricsRunner.class.getName());

        Runtime rt = Runtime.getRuntime();
        int numProcessors = rt.availableProcessors();
        long freeMemory = rt.freeMemory();
        long totalMemory = rt.totalMemory();
        long maxMemory = rt.maxMemory();

        String str = String.format(
            "Total amount of memory for use=%12d, free memory=%12d, constrained to max=%12d",
            totalMemory, freeMemory, maxMemory);

        log.info(str);

        rt = null;

        OperatingSystemMXBean bean = ManagementFactory.getOperatingSystemMXBean();

        str = String.format("arch=%s, os =%s %s, num processors=%d, current ave system load=%.3f",
            bean.getArch(), bean.getName(), bean.getVersion(), bean.getAvailableProcessors(),
            bean.getSystemLoadAverage());

        log.info(str);

        bean = null;
    }

    protected static void run(DoubleAxisIndexer indexer) throws IOException, TwoPointVoidStatsException {

        Logger log = Logger.getLogger(PerformanceMetricsRunner.class.getName());

        log.info(" (" + indexer.nXY + " points) ... ");


        TwoPointCorrelation clusterFinder = new TwoPointCorrelation(
            indexer.getX(), indexer.getY(),
            indexer.getXErrors(), indexer.getYErrors(), indexer.x.length);

        clusterFinder.setDebug(false);
        clusterFinder.logPerformanceMetrics();

        clusterFinder.calculateBackground();

        clusterFinder.findClusters();

        clusterFinder.calculateHullsOfClusters();
    }

    public static void main(String[] args) throws Exception {

        // generator background and stars with nPoints = 100, 1000, 10000, 100000 pausing in between

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        int nIter = 3;

        for (int ii = 0; ii < nIter; ii++) {

            for (int i = 0; i < 3; i++) {

                printSystemStats();

                RandomClusterAndBackgroundGenerator generator = new RandomClusterAndBackgroundGenerator();

                SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
                srr.setSeed( System.currentTimeMillis() );
                long seed = srr.nextLong();

                SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
                sr.setSeed( seed );
                //sr.setSeed(-2384802679227907254l);

                DoubleAxisIndexer indexer = null;

                int count = 0;

                switch(i) {
                    case 0:
                        //~100
                        indexer = generator.createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 33, 33, 0.1f);
                        break;
                    case 1:
                        //~1000
                        indexer = generator.createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 33, 33, 10f);
                        break;
                    case 2:
                        // 100*100
                        indexer = generator.createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 30, 60, 100.0f);
                        break;
                    case 3:
                        // 100*1000
                        indexer = generator.createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 30, 60, 1000.0f);
                        break;
                    default:
                        break;
                }

                indexer.sortAndIndexXThenY(generator.x, generator.y,
                    generator.xErrors, generator.yErrors, generator.x.length);

                generator = null;

                run(indexer);
            }
        }
    }
}
