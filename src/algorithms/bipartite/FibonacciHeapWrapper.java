package algorithms.bipartite;

import algorithms.bipartite.MinCostUnbalancedAssignment.PathNode;
import algorithms.imageProcessing.Heap;
import algorithms.imageProcessing.HeapNode;

/**
 *
 * @author nichole
 */
public class FibonacciHeapWrapper {

    private final int nApprox;
    
    private final Heap[] heaps;
    
    private long lastKnownMinKey = 0;
    private int lastKnownMinKeyIdx = 0;
    private long lastKnownMaxKey = -1;
    
    private final int binSz;
    
    private int n = 0;
    
    public FibonacciHeapWrapper(int nEstimate, int maxC) {
        
        int pow2 = (int)Math.ceil(Math.log(nEstimate)/Math.log(2));
    
        nApprox = nEstimate;
        
        if (pow2 > 6 && (maxC > 1)) {
            
            // make multiple maps separated by sequential partitions of values
            int nBins = pow2;
            
            if (maxC < pow2) {
                //pow2 range is 5 to 31
                nBins = maxC;
            }
            
            binSz = (int)Math.ceil((float)maxC/(float)nBins);
            
            heaps = new Heap[nBins];
            
        } else {
            heaps = new Heap[1];
            binSz = Integer.MAX_VALUE;
        }
    }
    
    public void insert(PathNode node) {
        
        int key = (int)node.getKey();
        
        int binIdx;
        if (heaps.length == 1) {
            binIdx = 0;
        } else {
            binIdx = key/binSz;
        }
        
        if (heaps[binIdx] == null) {
            heaps[binIdx] = new Heap();
        }
        
        heaps[binIdx].insert(node);
        
        n++;
        
        if (key < lastKnownMinKey) {
            lastKnownMinKey = key;
            lastKnownMinKeyIdx = binIdx;
        }
        if (key > lastKnownMaxKey) {
            lastKnownMaxKey = key;
        }
    }
    
    public void decreaseKey(PathNode node, long key2) {
        
        int key = (int)node.getKey();
        
        int binIdx;
        if (heaps.length == 1) {
            binIdx = 0;
        } else {
            binIdx = key/binSz;
        }
        
        heaps[binIdx].remove(node);
                
        node.setKey(key2);
        
        if (heaps.length > 1) {
            binIdx = (int)key2/binSz;
        }
     
        if (heaps[binIdx] == null) {
            heaps[binIdx] = new Heap();
        }
        
        heaps[binIdx].insert(node);
        
        if (key < lastKnownMinKey) {
            lastKnownMinKey = key;
            lastKnownMinKeyIdx = binIdx;
        }
    }
    
    public PathNode extractMin() {
        
        if (n == 0) {
            return null;
        }
     
        for (int i = (int)lastKnownMinKey; i < heaps.length; ++i) {
            
            if ((heaps[i] == null) || (heaps[i].getNumberOfNodes() == 0)) {
                continue;
            }
            
            HeapNode node = heaps[i].extractMin();
            lastKnownMinKey = i;
            n--;
            
            return (PathNode)node;
        }
        
        return null;
    }
    
    public long getNumberOfNodes() {
        return n;
    }
}
