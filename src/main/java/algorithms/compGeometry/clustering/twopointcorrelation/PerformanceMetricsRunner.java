package algorithms.compGeometry.clustering.twopointcorrelation;

import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;
import java.security.SecureRandom;
import java.util.logging.Logger;

import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION;

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

        long t0 = System.currentTimeMillis();

        TwoPointCorrelation clusterFinder = new TwoPointCorrelation(
            indexer.getX(), indexer.getY(),
            indexer.getXErrors(), indexer.getYErrors(), indexer.getX().length);

        //clusterFinder.setUseDownhillSimplexHistogramFitting();

        clusterFinder.setDebug(false);
        clusterFinder.logPerformanceMetrics();

        clusterFinder.calculateBackground();

        clusterFinder.findClusters();
        
        long t1 = (System.currentTimeMillis() - t0)/1000;
        
        log.info("  ====> total RT(sec) = " + t1);
    }

    public static void main(String[] args) throws Exception {

        // generator background and stars with nPoints = 100, 1000, 10000, 100000 pausing in between

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        int nIter = 3;

        for (int ii = 0; ii < nIter; ii++) {

            for (int i = 0; i < 5; i++) {

                printSystemStats();

                RandomClusterAndBackgroundGenerator generator = new RandomClusterAndBackgroundGenerator();

                SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
                srr.setSeed( System.currentTimeMillis() );
                long seed = srr.nextLong();

                SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
                //sr.setSeed( seed );
                sr.setSeed(-2384802679227907254l);

                DoubleAxisIndexer indexer = null;

                int count = 0;

                switch(i) {
                    case 0:
                        //~100
                        indexer = generator.createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            //10, 100, 110, 100.0f);
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
                    case 3: {
                        int[] clusterNumbers = new int[]{1000, 300, 100};
                        int nBackgroundPoints = 10000;
                        CLUSTER_SEPARATION clusterSeparation = CLUSTER_SEPARATION.LARGE;
                        
                        indexer = generator.createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            clusterNumbers, nBackgroundPoints, clusterSeparation);
                        break;
                    }
                    default: {
                        int[] clusterNumbers = new int[]{2000, 3000, 100};
                        int nBackgroundPoints = 100000;
                        CLUSTER_SEPARATION clusterSeparation = CLUSTER_SEPARATION.LARGE;
                        
                        indexer = generator.createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            clusterNumbers, nBackgroundPoints, clusterSeparation);
                        break;
                    }
                }

                indexer.sortAndIndexXThenY(generator.x, generator.y,
                    generator.xErrors, generator.yErrors, generator.x.length);

                generator = null;

                run(indexer);
            }
        }
    }
}
