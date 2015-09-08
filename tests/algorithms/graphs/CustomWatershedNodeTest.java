package algorithms.graphs;

import algorithms.util.PairInt;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class CustomWatershedNodeTest extends TestCase {
    
    public CustomWatershedNodeTest() {
    }

    public void testInsertOutgoing() {
        
        PairInt nodeLocation = new PairInt(0, 0);
        
        CustomWatershedNode node = new CustomWatershedNode(nodeLocation, 2);
        
        node.insertOutgoing(new PairInt(0, 1));
        
        node.insertOutgoing(new PairInt(1, 1));
        
        assertTrue(node.getConnectedNumber() == 2);
        
        assertTrue(node.get(0).equals(new PairInt(0, 1)));
        
        assertTrue(node.get(1).equals(new PairInt(1, 1)));
        
        node.reset(0, new PairInt(1, 0));
        
        assertTrue(node.get(0).equals(new PairInt(1, 0)));
        
        node.setToResolved(new PairInt(1, 1));
        
        assertTrue(node.isResolved());
        
        assertTrue(node.getResolved().equals(new PairInt(1, 1)));
        
        boolean caughtException = false;
        try {
            node.insertOutgoing(new PairInt(0, 10));
        } catch (Throwable t) {
            caughtException = true;
        }
        assertTrue(caughtException);
    }

}
