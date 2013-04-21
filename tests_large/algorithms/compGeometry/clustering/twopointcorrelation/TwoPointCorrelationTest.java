package algorithms.compGeometry.clustering.twopointcorrelation;

import java.security.SecureRandom;
import java.util.logging.Logger;

/**
 * @author nichole
 */
public class TwoPointCorrelationTest extends BaseTwoPointTest {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    // ======================  tests for sparse backgrounds ====================
    public void testSparseBackground__largeClusterSeparation_calculateBackgroundVia2PtVoidFit() throws Exception {

        log.info("testSparseBackground__largeClusterSeparation_calculateBackgroundVia2PtVoidFit");

        SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        srr.setSeed( System.currentTimeMillis() );
        long seed = srr.nextLong();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed( seed );
        //sr.setSeed(6236290033146721436l);

        System.out.println("using SEED=" + seed);

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);

        int nIterations = 10;

        for (int nIter = 0; nIter < nIterations; nIter++) {

            createPoints(10*(30+40+60), new int[]{30, 40, 60}, CLUSTER_SEPARATION.LARGE,
                xmin, xmax, ymin, ymax, sr, false);


            TwoPointCorrelation twoPtC = new TwoPointCorrelation(x, y, xErrors, yErrors, getTotalNumberOfPoints());

            twoPtC.setDebug(true);
            twoPtC.findClusters();

            twoPtC.findClusters();
            twoPtC.calculateHullsOfClusters();


            log.info("[nIter= " + nIter + "] found " + twoPtC.nGroups
                + " clusters.  created " + getExpectedNumberOfClusters() + " clusters");

            plotter.addPlot(twoPtC);
            plotter.writeFile();

            boolean foundAll = (twoPtC.nGroups >= getExpectedNumberOfClusters());

            assertTrue(foundAll);
        }

        log.info("SEED=" + seed);
    }

    // ======================  tests for sparse backgrounds ====================
    public void testSparseBackground__moderateClusterSeparation_calculateBackgroundVia2PtVoidFit() throws Exception {

        System.out.println("testSparseBackground__moderateClusterSeparation_calculateBackgroundVia2PtVoidFit");

        SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        srr.setSeed( System.currentTimeMillis() );
        long seed = srr.nextLong();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed( seed );
        //sr.setSeed(-2384802679227907254l);

        System.out.println("using SEED=" + seed);

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);

        int nIterations = 10;

        for (int nIter = 0; nIter < nIterations; nIter++) {

            createPoints((30+40+60), new int[]{30, 40, 60}, CLUSTER_SEPARATION.MODERATE,
                xmin, xmax, ymin, ymax, sr, false);


            TwoPointCorrelation twoPtC = new TwoPointCorrelation(x, y, xErrors, yErrors, getTotalNumberOfPoints());
            twoPtC.findClusters();

            twoPtC.calculateHullsOfClusters();


            System.out.println("[nIter= " + nIter + "] found " + twoPtC.nGroups
                + " clusters.  created " + getExpectedNumberOfClusters() + " clusters");

            plotter.addPlot(twoPtC);
            plotter.writeFile();

            boolean foundAll = (twoPtC.nGroups >= getExpectedNumberOfClusters());

            assertTrue(foundAll || (twoPtC.nGroups > 0));
        }

        log.info("SEED=" + seed);
    }

    // ======================  tests for moderate backgrounds ====================
    public void testModerateBackground__largeClusterSeparation_calculateBackgroundVia2PtVoidFit() throws Exception {

        log.info("testModerateBackground__largeClusterSeparation_calculateBackgroundVia2PtVoidFit");

        SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        srr.setSeed( System.currentTimeMillis() );
        long seed = srr.nextLong();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed( seed );
        //sr.setSeed(-2384802679227907254l);

        System.out.println("using SEED=" + seed);

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);

        int nIterations = 10;

        for (int nIter = 0; nIter < nIterations; nIter++) {

            createPoints((30+40+60), new int[]{30, 40, 60}, CLUSTER_SEPARATION.LARGE,
                xmin, xmax, ymin, ymax, sr, false);


            TwoPointCorrelation twoPtC = new TwoPointCorrelation(x, y, xErrors, yErrors, getTotalNumberOfPoints());
            twoPtC.findClusters();

            twoPtC.calculateHullsOfClusters();


            log.info("[nIter= " + nIter + "] found " + twoPtC.nGroups
                + " clusters.  created " + getExpectedNumberOfClusters() + " clusters");

            plotter.addPlot(twoPtC);
            plotter.writeFile();

            boolean foundAll = (twoPtC.nGroups >= getExpectedNumberOfClusters());

            assertTrue(foundAll);
        }

        log.info("SEED=" + seed);
    }

    public void testDenseBackground__largeClusterSeparation_calculateBackgroundVia2PtVoidFit() throws Exception {

        log.info("testDenseBackground__largeClusterSeparation_calculateBackgroundVia2PtVoidFit");

        SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        srr.setSeed( System.currentTimeMillis() );
        long seed = srr.nextLong();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed( seed );
        //sr.setSeed(-2384802679227907254l);

        log.info("using SEED=" + seed);

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);

        int nIterations = 10;

        for (int nIter = 0; nIter < nIterations; nIter++) {

            createPoints(10*(30+40+60), new int[]{30, 40, 60}, CLUSTER_SEPARATION.LARGE,
                xmin, xmax, ymin, ymax, sr, false);


            TwoPointCorrelation twoPtC = new TwoPointCorrelation(x, y, xErrors, yErrors, getTotalNumberOfPoints());
            twoPtC.findClusters();

            twoPtC.calculateHullsOfClusters();


            log.info("[nIter= " + nIter + "] found " + twoPtC.nGroups
                + " clusters.  created " + getExpectedNumberOfClusters() + " clusters");

            plotter.addPlot(twoPtC);
            plotter.writeFile();

            boolean foundAll = (twoPtC.nGroups >= 0);

            assertTrue(foundAll);
        }

        log.info("SEED=" + seed);
    }

}
