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

        log.info("generate tmpdata2/similar_curve_parameters.txt");

        if (!enable) {
            return;
        }
        
        GEVSimilarityTool tool = new GEVSimilarityTool();
        
        tool.resetForNWithinTen(20.0f);

        tool.calculateCurveDiffs();
        
        tool.plotResults();
    }
    
}
