package algorithms.curves;

import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 * a tool to find common curves among the large range of possible curves
 * in the GEV function which depends upon k, sigma, and (x-mu).
 * 
 * @author nichole
 */
public class GEVSimilarityToolTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected boolean debug = true;

    protected boolean enable = true;
    
    
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    public void test0() throws Exception {

        log.info("test0()");

        if (!enable) {
            return;
        }
        
        int persistFileNum = 2;
        
        boolean usePersisted = false;
        boolean persist = true;
              
        if (usePersisted && persist) {
            System.err.println("Cannot have usePersisted=true and persist=true");
            return;
        }
        
        GEVSimilarityTool tool = new GEVSimilarityTool();

        if (usePersisted) {
            tool.readPersisted(persistFileNum);
        } else {
            tool.calculateCurveDiffs();
        }
        
        if (persist) {
            tool.persist(persistFileNum);
        }
        
        tool.plotResults();
    }
    
}
