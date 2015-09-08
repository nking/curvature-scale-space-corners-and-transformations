package algorithms.graphs;

import algorithms.util.PairInt;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class CustomWatershedDAGTest extends TestCase {
    
    public CustomWatershedDAGTest() {
    }

    public void testInsert() {
       
        CustomWatershedDAG dag = new CustomWatershedDAG();
        
        int[] diffInt = new int[]{1, 2, 10};
        PairInt[] points = new PairInt[]{new PairInt(2, 1), new PairInt(1, 2),
            new PairInt(2, 3)
        };
        
        PairInt key = new PairInt(2, 2);
        
        dag.orderAndInsert(key, diffInt, points, 3);
        
        assertFalse(dag.isResolved(key));

        int nConnected = dag.getConnectedNumber(key);
        assertTrue(nConnected == 3);
        
        assertTrue(dag.contains(key));
        
        assertFalse(dag.contains(new PairInt(100000, 200000)));
        
        PairInt p0 = dag.getConnectedNode(key, 0);
        assertTrue(p0.equals(new PairInt(2, 3)));
        
        PairInt p1 = dag.getConnectedNode(key, 1);
        assertTrue(p1.equals(new PairInt(1, 2)));
        
        PairInt p2 = dag.getConnectedNode(key, 2);
        assertTrue(p2.equals(new PairInt(2, 1)));
    }

}
