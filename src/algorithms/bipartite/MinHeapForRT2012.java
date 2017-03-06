package algorithms.bipartite;

import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.HeapNode;
import java.util.logging.Logger;

/**
 * a min heap for the MinCostUnbalancedAssignment.
 * It uses the "Dial algorithm" pattern by default
 * (which is similar to counting sort in part),
 * but if the requested capacity is larger than
 * 46300, a Fibonacci Heap is used internally instead.
 * All operations for the "Dial" mode are essentially O(1),
 * else, the extractMin operation is O(lg_2(N_nodes))
 * for the internal Fibonacci heap.
 * 
 * NOTE: in the future, would like to include MLB for
 * another option that improves the extractMin runtime
 * complexity.
 * 
 * @author nichole
 */
public class MinHeapForRT2012 {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    // 0 = use Dial, 1 = use Fibonacci, 2 = use XFastTrie
    private final int algorithm;
    
    private final DoubleLinkedCircularList[] heap0;
    
    private int lastKnownMinKey0 = 0;
    
    // for use in tuning the capacity
    private int lastKnownMaxKey0 = 0;
    
    private long n0 = 0;
    
    private final FibonacciHeapWrapper heap1;

    private final XFastTrieWrapper heap2;
    
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
     * If capacity is less than 46300 (might lower this as memory will begin
     * to affect performance at high capacity),
     * then a minimum priority monotonic bucket queue which uses the node keys as priorities
     * is used.  creation of the structure has runtime complexity O(capacity), 
     * but thereafter, all operations are O(1).
     * If capacity is higher than 46300, then the approxN is used to 
     * decide between an XFastTrie and a FibonacciHeap.
     * 
     * TODO: implement a multi-level bucket such as Cherkassky, Goldberg,
     * and Silverstein 1999 or Raman 1996 which has a faster extractMin
     * than a Fibonacci Heap.
     */
    public MinHeapForRT2012(int capacity, int approxN) {
                
        if (capacity < 46300) {
            //TODO:  need MLB for this...
        
            // the 1 level Dial algorithm has O(1) inserts and
            //    constant time extractMin.
            algorithm = 0;
        
            heap0 = new DoubleLinkedCircularList[capacity];
            
            heap1 = null;
            
            heap2 = null;

            log.fine("useDial=true cap=" + capacity + " approxN=" + approxN);
            
        } else if (approxN > ((1<<(6-1))-1)) {
            // wanting the base of the prefix tree to be filled
            // to improve performance.   for larger N and
            // this range of maxC the XFastTrie has better extractMin.
            // could reduce the conditional to '5' instead of '6' 
            
            algorithm = 2;
            
            heap2 = new XFastTrieWrapper(capacity);
            
            heap0 = null;
            
            heap1 = null;
            
        } else {
            
            algorithm = 1;
        
            heap0 = null;
            
            heap1 = new FibonacciHeapWrapper(approxN, capacity);
        
            heap2 = null;
        }
    }
    
    public void insert(HeapNode node) {
        
        if (node.getKey() < 0) {
            throw new IllegalArgumentException("key must be >= 0");
        }
        
        switch(algorithm) {
            case 0:
                insert0(node);
                break;
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
            case 0:
                return extractMin0();
            case 1:
                return extractMin1();
            default:
                return extractMin2();
        }
        
    }
    
    private HeapNode extractMin0() {
    
        for (int i = lastKnownMinKey0; i < heap0.length; ++i) {
            DoubleLinkedCircularList bucket = heap0[i];
            if (bucket != null && (bucket.getNumberOfNodes() > 0)) {
                HeapNode node = bucket.getSentinel().getLeft();
                bucket.remove(node);
                lastKnownMinKey0 = i;
                n0--;
                return (HeapNode)node;
            }
        }
        
        return null;
    }
    
    private HeapNode extractMin1() {
        
        HeapNode node = heap1.extractMin();
        if (node != null) {
            return (HeapNode)node;
        } else {
            return null;
        }
    }
    
    private HeapNode extractMin2() {        
        return heap2.extractMin();
    }
    
    private void insert0(HeapNode node) {
         
        int key = (int)node.getKey();
        
        DoubleLinkedCircularList bucket = heap0[key];
        if (bucket == null) {
            bucket = new DoubleLinkedCircularList();
            heap0[key] = bucket;
        }
        
        log.fine("insert into minHeap at key =" + node.toString());
       
        bucket.insert(node);
        n0++;
        
        if (key < lastKnownMinKey0) {
            lastKnownMinKey0 = key;
        }
        if (key > lastKnownMaxKey0) {
            lastKnownMaxKey0 = key;
        }
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
            case 0:
                decreaseKey0(node, key2);
                break;
            case 1:
                decreaseKey1(node, key2);
                break;
            default:
                decreaseKey2(node, key2);
                break;
        }
    }
     
    private void decreaseKey0(HeapNode node, long key2) {

        log.fine("decreaseKey in minHeap from key=" + 
            node.getKey() + " to key=" + key2);
        
        int prevKey = (int)node.getKey();
        heap0[prevKey].remove(node);
        
        node.setKey(key2);
        
        insert0(node);
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
            case 0:
                return n0;
            case 1:
                return heap1.getNumberOfNodes();
            default:
                return heap2.getNumberOfNodes();
        }
    }
}
