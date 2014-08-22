package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.util.ArrayPair;

import java.util.logging.Logger;
import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;
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

    public void testFindVoidsSemiCompleteSearch() throws Exception {

        log.info("testFindVoidsSemiCompleteSearch");

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

        log.fine("*1*");
        
        createPoints((30+40+60), new int[]{30, 40, 60},
            RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION.LARGE,
            xmin, xmax, ymin, ymax, sr, false);

        log.fine("*2*");
        
        AxisIndexer indexer = new AxisIndexer();
        indexer.sortAndIndexX(generator.x, generator.y,
            generator.xErrors, generator.yErrors, generator.x.length);

        log.fine("*3*");
        
        TwoPointVoidStats stats = new TwoPointVoidStats(indexer);
        stats.setDebug(false);
        //stats.setStandardDeviationFactor(2.5f);
        
        log.fine("*4*");
        
        stats.calc();

        log.fine("*5*");
        
        assertNotNull(stats.getBestFit());
    }
}
