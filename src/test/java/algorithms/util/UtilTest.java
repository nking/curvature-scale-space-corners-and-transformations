package algorithms.util;


import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class UtilTest extends TestCase {
    
    /**
     *
     * @throws Exception
     */
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    /**
     *
     * @throws Exception
     */
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    /**
     *
     * @throws Exception
     */
    public void test() throws Exception {
        
        long availMem = Util.getAvailableHeapMemory();
        assertTrue(availMem > 0);
        
        // for coverage:
        Util u = new Util();
    }
}
