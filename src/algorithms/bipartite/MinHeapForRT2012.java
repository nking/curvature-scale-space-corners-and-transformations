package algorithms.bipartite;

import algorithms.YFastTrie;
import algorithms.imageProcessing.HeapNode;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.util.logging.Logger;

/**
 * A min heap for the MinCostUnbalancedAssignment.
 * It uses the "YFastTrie min priority queue algorithm" 
 * pattern by default if the VM has enough memory, else
 * uses a Fibonacci Heap.
 * All operations for the "YFastTrie" are constant time.
 * 
 * The Fibonacci Heap has O(1) operations excepting
 * extractMin which is O(lg_2(N_nodes)).
 * 
 * @author nichole
 */
public class MinHeapForRT2012 {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    // 1 = Fibonacci, 2 = YFastTrie
    private final int algorithm;
        
    private int lastKnownMinKey0 = 0;
    
    // for use in tuning the capacity
    private int lastKnownMaxKey0 = 0;
    
    private final FibonacciHeapWrapper heap1;

    private final YFastTrieWrapper heap2;
    
    void printLastKnownMinMax() {
        log.fine("min=" + lastKnownMinKey0
            + " max=" + lastKnownMaxKey0);
    }
    
    /**
     * 
     * @param capacity estimate of maximum value to store.
     * @param approxN approximate number of nodes expected
     * to be in the heap as a rough maximum at a given time.
     * (it's used to help determine which algorithm to use
     * internally).
     * 
     * IA min heap for the MinCostUnbalancedAssignment.
     * It uses the "YFastTrie min priority queue algorithm" 
     * pattern by default if the VM has enough memory, else
     * uses a Fibonacci Heap.
     * All operations for the "YFastTrie" are constant time.
     * 
     * The Fibonacci Heap has O(1) operations excepting
     * extractMin which is O(lg_2(N_nodes)).
     * @param maxNumberOfBits
     * 
     */
    public MinHeapForRT2012(int capacity, int approxN, int maxNumberOfBits) {

        //use yfasttrie if theres enough memory        
        long totalMemory = Runtime.getRuntime().totalMemory();
        MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
        long heapUsage = mbean.getHeapMemoryUsage().getUsed();
        long avail = totalMemory - heapUsage;

        long[] yftEstimate = YFastTrie.estimateSizeOnHeap(capacity, 
                maxNumberOfBits);
        
        log.fine("avail=" + avail + " yftEst=" + yftEstimate[1] + " < " +
            (yftEstimate[1] < avail));
        
        if (yftEstimate[1] < avail) {
            // wanting the base of the prefix tree to be filled
            // to improve performance.   for larger N
            
            algorithm = 2;
            
            heap2 = new YFastTrieWrapper(capacity);
            
            heap1 = null;
            
        } else {
            
            algorithm = 1;
        
            heap1 = new FibonacciHeapWrapper(approxN, capacity);
        
            heap2 = null;
        }
    }
    
    public void insert(HeapNode node) {
        
        if (node.getKey() < 0) {
            throw new IllegalArgumentException("key must be >= 0");
        }
        
        switch(algorithm) {
            case 1:
                insert1(node);
                break;
            default:
                insert2(node);
                break;
        }
    }
    
    public HeapNode extractMin() {
        
        switch(algorithm) {
            case 1:
                return extractMin1();
            default:
                return extractMin2();
        }
        
    }
    
    private HeapNode extractMin1() {
        
        HeapNode node = heap1.extractMin();
        if (node != null) {
            return node;
        } else {
            return null;
        }
    }
    
    private HeapNode extractMin2() {        
        return heap2.extractMin();
    }
    
    private void insert1(HeapNode node) {
         
        int key = (int)node.getKey();
        
        heap1.insert(node);
        
        log.fine("insert into minHeap at key =" + key);        
    }
    
    private void insert2(HeapNode node) {
         
        int key = (int)node.getKey();
        
        heap2.insert(node);
        
        log.fine("insert into minHeap at key =" + key);        
    }
    
    public void decreaseKey(HeapNode node, long key2) {
    
        switch(algorithm) {
            case 1:
                decreaseKey1(node, key2);
                break;
            default:
                decreaseKey2(node, key2);
                break;
        }
    }
     
    private void decreaseKey1(HeapNode node, long key2) {

        log.fine("decreaseKey in fibHeap from key=" + 
            node.getKey() + " to key=" + key2);
        
        heap1.decreaseKey(node, key2);
    }
    
    private void decreaseKey2(HeapNode node, long key2) {

        log.fine("decreaseKey in fibHeap from key=" + 
            node.getKey() + " to key=" + key2);
        
        heap2.decreaseKey(node, key2);
    }
    
    public long getNumberOfNodes() {
        
        switch(algorithm) {
            case 1:
                return heap1.getNumberOfNodes();
            default:
                return heap2.getNumberOfNodes();
        }
    }

    @Override
    public String toString() {
    
        StringBuilder sb = new StringBuilder();
        sb.append("min heap type = ");
        switch(algorithm) {
            case 1:
                sb.append("Fibonacci Heap"); break;
            default:
                sb.append("YFastTrie min priority queue"); break;
        }
        sb.append(". size=").append(getNumberOfNodes());
        
        return sb.toString();
    }
    
    
}
