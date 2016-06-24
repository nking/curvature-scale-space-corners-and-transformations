package algorithms.search;

import algorithms.util.PairFloat;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class KDTreeFloatTest extends TestCase {
    
    public KDTreeFloatTest() {
    }
 
    public void test() {
        
		float[] x = new float[] {5, 5, 9, 3, 4, 1, 7, 2};
    	float[] y = new float[] {5, 5, 6, 6, 9, 1, 9, 9};
    	
    	int lastUsableIndex = KDTreeFloat.reduceToUniqueWithMoveUp(x, y);
    	
		assertEquals(6, lastUsableIndex);
		
		float[] expectedx = new float[] {5, 9, 3, 4, 1, 7, 2};
    	float[] expectedy = new float[] {5, 6, 6, 9, 1, 9, 9};
    	
    	for (int i = 0; i < expectedx.length; i++) {
    		assertEquals(expectedx[i], x[i]);
    		assertEquals(expectedy[i], y[i]);
    	}
		
	}
    
    public void test2() {
        
		float[] x = new float[] {5, 5, 9, 3, 4, 1, 7, 2};
    	float[] y = new float[] {5, 5, 6, 6, 9, 1, 9, 9};
    	
        KDTreeFloat kdTree = new KDTreeFloat(x, y, false);
    	        
        Set<PairFloat> nodes = 
            kdTree.findNearestNeighbor(6f, 8f);
    	assertNotNull(nodes);
        assertEquals(1, nodes.size());
        PairFloat node = nodes.iterator().next();
    	assertEquals(7f, node.getX());
    	assertEquals(9f, node.getY());
	}
    
    public void test3() {
        
		float[] x = new float[] {5, 5, 9, 13, 14, 21, 27, 32};
    	float[] y = new float[] {5, 5, 6,  6,  9,  1,  9,  9};
    	
        KDTreeFloat kdTree = new KDTreeFloat(x, y, true);
    	
        //kdTree.printTree();
        
        Set<PairFloat> nodes = 
            kdTree.findNearestNeighbor(12f, 7f);
    	assertNotNull(nodes);
        assertEquals(1, nodes.size());
        PairFloat node = nodes.iterator().next();
    	assertEquals(13f, node.getX());
    	assertEquals(6f, node.getY());
	}
    
    public void test4() {
        
		float[] x = new float[] {9, 3, 4, 1, 7, 2};
    	float[] y = new float[] {6, 6, 9, 1, 9, 9};
    	
        KDTreeFloat kdTree = new KDTreeFloat(x, y, true);
    	
        //kdTree.printTree();
       
        Set<PairFloat> nodes = 
            kdTree.findNearestNeighbor(5f, 7f);
    	assertNotNull(nodes);
        assertEquals(2, nodes.size());
        
        Set<PairFloat> expected = new HashSet<PairFloat>();
        expected.add(new PairFloat(3, 6));
        expected.add(new PairFloat(4, 9));
        
        for (PairFloat p : nodes) {
            assertTrue(expected.remove(p));
        }
        assertTrue(expected.isEmpty());
        
	}
}
