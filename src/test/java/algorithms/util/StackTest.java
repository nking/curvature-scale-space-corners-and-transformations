package algorithms.util;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class StackTest extends TestCase {

    /**
     *
     * @param testName
     */
    public StackTest(String testName) {
        super(testName);
    }

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
     */
    public void testInsertAndPop() {

        Stack stack = new Stack();
        
        int[] tasks = new int[5];

        for (int i = 0; i < 5; i++) {

            tasks[i] = i + 10;
            
            stack.insert(tasks[i]);
        }

        int count = tasks.length;

        while (!stack.isEmpty()) {

            int task = stack.pop().getKey();

            assertTrue(tasks[count - 1] == task);

            count--;
        }
        
        
        boolean threwException = false;
        try {
            stack.insert(-1);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        stack = new Stack();
        assertNull(stack.pop());
    }

}
