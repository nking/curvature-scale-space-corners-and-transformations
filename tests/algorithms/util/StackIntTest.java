package algorithms.util;

import junit.framework.TestCase;

/**
  adapted from 
  http://code.google.com/p/two-point-correlation/source/browse/src/test/java/algorithms/util/StackTest.java
  under MIT License (MIT), Nichole King 2013 
 * @author nichole
 */
public class StackIntTest extends TestCase {

    public StackIntTest(String testName) {
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

    public void testInsertAndPop() {

        StackInt stack = new StackInt();
        
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
        
        stack = new StackInt();
        assertNull(stack.pop());
    }

}
