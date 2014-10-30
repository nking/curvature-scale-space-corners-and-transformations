package algorithms.imageProcessing;

import junit.framework.TestCase;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class StackTest extends TestCase {

    public StackTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void testPushPop() {

        System.out.println("testPushPop");


        Stack<Integer> stack = new Stack<Integer>();
        
        assertTrue(stack.isEmpty());
        
        int n = 100;
        for (int i = 0; i < n; i++) {
            stack.push(Integer.valueOf(i));
        }
        
        assertFalse(stack.isEmpty());
        
        for (int i = (n - 1); i > -1; i--) {
            
            Object peek = stack.peek();
            
            Object obj = stack.pop();
            
            assertTrue(((Integer)obj).intValue() == i);
            
            assertEquals(peek, obj);
        }
        
        assertTrue(stack.isEmpty());
    }

    public static void main(String[] args) {
        try {
            
            StackTest test = new StackTest("StackTest");
            
            test.testPushPop();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
