package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.util.PolygonAndPointPlotter;
import java.security.SecureRandom;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class FittingStatsTest extends BaseTwoPointTest {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public void testFitStats() throws Exception {

        log.info("testFitStats()");

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
        int nIterPerBackground = 1;

        int m = nIterPerBackground*3;
        
        AxisIndexer indexer = null;

        int count = 0;

        for (int i = 0; i < 4; i++) {

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
                    case 3: {
                        createPoints((30+40+60), new int[]{30, 40, 60},
                            RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION.LARGE,
                            xmin, xmax, ymin, ymax, sr, false);
                        indexer = new AxisIndexer();
                        indexer.sortAndIndexX(generator.x, generator.y,
                            generator.xErrors, generator.yErrors, generator.x.length);
                        break;
                    }
                    default:
                        break;
                }

                indexer.sortAndIndexX(generator.x, generator.y,
                    generator.xErrors, generator.yErrors, generator.x.length);

                log.info(" " + count + "(" + indexer.nXY + " points) ... ");

                TwoPointVoidStats stats = new TwoPointVoidStats(indexer);
                stats.setDebug(false);
                stats.calc();

                stats.plotFit(plotter2);

                count++;
            }
        }

        log.info("\n start computing stats for all sets");        
    }

}
