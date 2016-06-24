package algorithms.search;

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
    	        
        KDTreeNodeFloat node = 
            kdTree.findNearestNeighbor(6f, 8f);
    	assertNotNull(node);
    	assertEquals(7f, node.x);
    	assertEquals(9f, node.y);
	}
    
    public void test3() {
        
		float[] x = new float[] {5, 5, 9, 13, 14, 21, 27, 32};
    	float[] y = new float[] {5, 5, 6,  6,  9,  1,  9,  9};
    	
        KDTreeFloat kdTree = new KDTreeFloat(x, y, true);
    	
        kdTree.printTree();
        
        KDTreeNodeFloat node = kdTree.findNearestNeighbor(12f, 7f);
    	assertNotNull(node);
    	assertEquals(13f, node.x);
    	assertEquals(6f, node.y);
	}
}
