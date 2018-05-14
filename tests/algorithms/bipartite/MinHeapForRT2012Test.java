package algorithms.bipartite;

import algorithms.imageProcessing.HeapNode;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MinHeapForRT2012Test extends TestCase {
    
    public MinHeapForRT2012Test() {
    }
    
    /**
     * Test of printLastKnownMinMax method, of class MinHeapForRT2012.
     */
    public void test0() throws NoSuchAlgorithmException {
        
        SecureRandom random = SecureRandom.getInstanceStrong();
        long seed = System.currentTimeMillis();
        System.out.println("SEED=" + seed);
        random.setSeed(seed);
        
        int nTests = 1;
        int maxValue = 150000;
        int k = 0;
        boolean typeFibonacci = false;
        boolean typeYFT = false;
        
        while ((!typeFibonacci || !typeYFT) && (k < 4)) {
            
            int nT = nTests;
            int maxV = maxValue * (k + 1);
            int approxN = maxV * 10;
            int maxNumberOfBits = (int)Math.ceil(Math.log(maxV)/Math.log(2));
            
            for (int t = 0; t < nT; ++t) {
                
                MinHeapForRT2012 heap = new MinHeapForRT2012(maxV, approxN,
                    maxNumberOfBits);
                
                System.out.println(heap.toString());
                if (heap.toString().contains("ibonacci")) {
                    typeFibonacci = true;
                } else {
                    typeYFT = true;
                }
                System.out.println("maxV=" + maxV + " approxN=" + approxN + 
                    " maxNBits=" + maxNumberOfBits);
                
                for (int i = 0; i < approxN; ++i) {
                    int key = random.nextInt(maxV);
                    HeapNode node = new HeapNode(key);
                    node.setData(Integer.valueOf(key));
                    heap.insert(node);
                }
                
                assertEquals(approxN, (int)heap.getNumberOfNodes());
                
                int lastKey = Integer.MIN_VALUE;
                for (int i = 0; i < approxN; ++i) {
                    HeapNode node = heap.extractMin();
                    int key = (int)node.getKey();
                    assertTrue(key >= lastKey);
                    assertEquals(key, ((Integer)node.getData()).intValue());
                    lastKey = key;
                }
            }
            k++;
        }
    }
    
}
