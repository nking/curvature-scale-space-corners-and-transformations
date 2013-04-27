package algorithms.compGeometry.clustering.twopointcorrelation;

import java.security.SecureRandom;
import java.util.logging.Logger;
import static junit.framework.Assert.*;

/**
 *
 * @author nichole
 */
public class TwoPointVoidStats3Test extends BaseTwoPointTest {

    boolean debug = true;

    protected Logger log = Logger.getLogger(this.getClass().getName());

    public void testFindVoidsRoughRangeSearch() throws Exception {

        log.info("testFindVoidsRoughRangeSearch");

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

        System.out.println("*1*");
        
        createPoints((30+40+60), new int[]{30, 40, 60},
            RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION.LARGE,
            xmin, xmax, ymin, ymax, sr, false);

        System.out.println("*2*");
        
        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(generator.x, generator.y,
            generator.xErrors, generator.yErrors, generator.x.length);

        System.out.println("*3*");
        
        TwoPointVoidStats stats = new TwoPointVoidStats(indexer);
        stats.setDebug(false);
        
        stats.setUseSemiCompleteRangeSampling();

        System.out.println("*4*");
        
        stats.calc();

        System.out.println("*5*");
        
        assertNotNull(stats.getBestFit());
    }
}
