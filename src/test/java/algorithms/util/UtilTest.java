package algorithms.util;


import junit.framework.TestCase;

public class UtilTest extends TestCase {
    
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    public void test() throws Exception {
        
        long availMem = Util.getAvailableHeapMemory();
        assertTrue(availMem > 0);
        
        // for coverage:
        Util u = new Util();
    }
}
