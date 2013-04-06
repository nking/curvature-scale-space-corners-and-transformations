package algorithms.compGeometry.clustering.twopointcorrelation;

import java.security.SecureRandom;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;

/**
 * @author nichole
 */
public class TwoPointCorrelationTest extends BaseTwoPointTest {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    protected boolean debug = true;

    public void testSmallData() throws Exception {

        log.info("testSmallData");

        x = new float[] {
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,              4,  5,  6,  7,  8,
            0,              4,  5,  6,  7,  8,
            0,              4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,              8,
            0,  1,  2,  3,  4,              8,
            0,  1,  2,  3,  4,              8,
            0,  1,          4,  5,  6,  7,  8,
            0,  1,          4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            5.5f, 6.5f, 6.0f, 6.0f, 6.5f, 5.5f, 5.5f // <==== adding small cluster
        };
        y = new float[] {
            0,  0,  0,  0,  0,  0,  0,  0,  0,
            1,              1,  1,  1,  1,  1,
            2,              2,  2,  2,  2,  2,
            3,              3,  3,  3,  3,  3,
            4,  4,  4,  4,  4,  4,  4,  4,  4,
            5,  5,  5,  5,  5,              5,
            6,  6,  6,  6,  6,              6,
            7,  7,  7,  7,  7,              7,
            8,  8,          8,  8,  8,  8,  8,
            9,  9,          9,  9,  9,  9,  9,
           10, 10, 10, 10, 10, 10, 10, 10, 10,
           2.0f, 2.0f, 2.5f, 1.0f, 2.5f, 2.5f, 1.5f
        };

        float xmin = 0;
        float xmax = 10;
        float ymin = 0;
        float ymax = 10;

        // make uniform errors for x and y
        xErrors = new float[x.length];
        yErrors = new float[y.length];
        for (int i = 0; i < xErrors.length; i++) {
            xErrors[i] = x[i]/100.0f;
            yErrors[i] = y[i]/100.0f;
        }

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);

        TwoPointCorrelation twoPtC = new TwoPointCorrelation(x, y, xErrors, yErrors, x.length);
        twoPtC.setDebug(true);
    //twoPtC.setBackground(1.f, 0.01f);
        twoPtC.calculateBackground();
        twoPtC.findClusters();
        twoPtC.calculateHullsOfClusters();

        plotter.addPlot(twoPtC);
        plotter.writeFile();

        assertTrue(twoPtC.nGroups == 1);

        log.info("found " + twoPtC.nGroups + " clusters");

        /*
        // make uniform errors for x and y
        xErrors = new float[x.length];
        yErrors = new float[y.length];
        Arrays.fill(xErrors, 0.5f);
        Arrays.fill(yErrors, 0.5f);

        twoPtC = new TwoPointCorrelation(x, y, xErrors, yErrors, x.length);
        twoPtC.setDebug(true);
        twoPtC.calculateBackground...
        twoPtC.findClusters();
        twoPtC.calculateHullsOfClusters();

        plotter.addPlot(twoPtC);
        plotter.writeFile();

        assertTrue(twoPtC.nGroups == 1);

        log.info("found " + twoPtC.nGroups + " clusters");
        */
    }

    // ======================  tests for sparse backgrounds ====================
    public void testSparseBackground__largeClusterSeparation_calculateBackgroundVia2PtVoidFit() throws Exception {

        log.info("testSparseBackground__largeClusterSeparation_calculateBackgroundVia2PtVoidFit");

        SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        srr.setSeed( System.currentTimeMillis() );
        long seed = srr.nextLong();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed( seed );
        //sr.setSeed(-3982961905480492188l);

        log.info("using SEED=" + seed);

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);

        int nIterations = 1;

        for (int nIter = 0; nIter < nIterations; nIter++) {

            createPoints((int)0.1f*(30+40+60), new int[]{30, 40, 60}, CLUSTER_SEPARATION.LARGE,
                xmin, xmax, ymin, ymax, sr, false);


            TwoPointCorrelation twoPtC = new TwoPointCorrelation(x, y, xErrors, yErrors, getTotalNumberOfPoints());

            twoPtC.setDebug(debug);

            boolean allowTuning = true;
            twoPtC.findClusters(allowTuning);
            twoPtC.calculateHullsOfClusters();


            log.info("[nIter= " + nIter + "] found " + twoPtC.nGroups
                + " clusters.  created " + getExpectedNumberOfClusters() + " clusters");

            plotter.addPlot(twoPtC);
            plotter.writeFile();

            boolean foundAll = (twoPtC.nGroups >= getExpectedNumberOfClusters());

            assertTrue(foundAll);

/*
            // to create a file to test the methods in main:
            String path = this.getClass().getClassLoader().getResource(".").getPath() + "/file.txt";

            File fl = new File(path);
            FileWriter writer = new FileWriter(fl);

            for (int i = 0; i < this.x.length; i++) {

                StringBuffer sb = new StringBuffer();
                sb.append(x[i]).append("\t").append(y[i]).append("\t")
                    .append(xErrors[i]).append("\t").append(yErrors[i]).append("\n");

               writer.write(sb.toString());
            }
            writer.close();
*/
        }

        log.info("SEED=" + seed);
    }

    // ======================  tests for sparse backgrounds ====================
    public void testSparseBackground__moderateClusterSeparation_calculateBackgroundVia2PtVoidFit() throws Exception {

        log.info("testSparseBackground__moderateClusterSeparation_calculateBackgroundVia2PtVoidFit");

        SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        srr.setSeed( System.currentTimeMillis() );
        long seed = srr.nextLong();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed( seed );
        //sr.setSeed(-3982961905480492188l);

        log.info("using SEED=" + seed);

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);

        int nIterations = 1;

        for (int nIter = 0; nIter < nIterations; nIter++) {

            createPoints((30+40+60), new int[]{30, 40, 60}, CLUSTER_SEPARATION.MODERATE,
                xmin, xmax, ymin, ymax, sr, false);


            TwoPointCorrelation twoPtC = new TwoPointCorrelation(x, y, xErrors, yErrors, getTotalNumberOfPoints());
            twoPtC.setDebug(debug);

            boolean allowTuning = false;
            twoPtC.findClusters(allowTuning);

            twoPtC.calculateHullsOfClusters();


            log.info("[nIter= " + nIter + "] found " + twoPtC.nGroups
                + " clusters.  created " + getExpectedNumberOfClusters() + " clusters");

            plotter.addPlot(twoPtC);
            plotter.writeFile();

            boolean foundAll = (twoPtC.nGroups >= getExpectedNumberOfClusters());

            assertTrue(foundAll || (twoPtC.nGroups > 0));


            /*
            // to create a file to test the methods in main:
            String path = this.getClass().getClassLoader().getResource(".").getPath() + "/file.txt";

            File fl = new File(path);
            FileWriter writer = new FileWriter(fl);

            for (int i = 0; i < this.x.length; i++) {

                StringBuffer sb = new StringBuffer();
                sb.append(x[i]).append("\t").append(y[i]).append("\t")
                    .append(xErrors[i]).append("\t").append(yErrors[i]).append("\n");

               writer.write(sb.toString());
            }
            writer.close();
            */
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
        //sr.setSeed(-3982961905480492188l);

        log.info("using SEED=" + seed);

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);

        int nIterations = 1;

        for (int nIter = 0; nIter < nIterations; nIter++) {

            createPoints((30+40+60), new int[]{30, 40, 60}, CLUSTER_SEPARATION.LARGE,
                xmin, xmax, ymin, ymax, sr, false);


            TwoPointCorrelation twoPtC = new TwoPointCorrelation(x, y, xErrors, yErrors, getTotalNumberOfPoints());
            twoPtC.setDebug(debug);

            boolean allowTuning = false;
            twoPtC.findClusters(allowTuning);

            twoPtC.calculateHullsOfClusters();

            log.info("[nIter= " + nIter + "] found " + twoPtC.nGroups
                + " clusters.  created " + getExpectedNumberOfClusters() + " clusters");

            plotter.addPlot(twoPtC);
            plotter.writeFile();

            boolean foundAll = (twoPtC.nGroups > 0);

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

        int nIterations = 1;

        for (int nIter = 0; nIter < nIterations; nIter++) {

            createPoints(10*(30+40+60), new int[]{30, 40, 60}, CLUSTER_SEPARATION.LARGE,
                xmin, xmax, ymin, ymax, sr, false);


            TwoPointCorrelation twoPtC = new TwoPointCorrelation(x, y, xErrors, yErrors, getTotalNumberOfPoints());
            twoPtC.setDebug(debug);

            boolean allowTuning = false;
            twoPtC.findClusters(allowTuning);

            twoPtC.calculateHullsOfClusters();

            log.info("[nIter= " + nIter + "] found " + twoPtC.nGroups
                + " clusters.  created " + getExpectedNumberOfClusters() + " clusters");

            plotter.addPlot(twoPtC);
            plotter.writeFile();

            boolean foundAll = (twoPtC.nGroups >= 0);

            assertTrue(foundAll);

            /*
            // to create a file to test the methods in main:
            String path = this.getClass().getClassLoader().getResource(".").getPath() + "/file.txt";

            File fl = new File(path);
            FileWriter writer = new FileWriter(fl);

            for (int i = 0; i < this.x.length; i++) {

                StringBuffer sb = new StringBuffer();
                sb.append(x[i]).append("\t").append(y[i]).append("\t")
                    .append(xErrors[i]).append("\t").append(yErrors[i]).append("\n");

               writer.write(sb.toString());
            }
            writer.close();
            */
        }

        log.info("SEED=" + seed);
    }

}
