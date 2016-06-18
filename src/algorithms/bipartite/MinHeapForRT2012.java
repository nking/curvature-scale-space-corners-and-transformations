package algorithms.bipartite;

import algorithms.bipartite.MinCostUnbalancedAssignment.PathNode;
import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.Heap;
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
 *  NOTE: in the future, will replace the fibonacci heap here
    with MLB.
 * @author nichole
 */
public class MinHeapForRT2012 {

    
    /*    
    notes for the bipartite weighted min-cost matching:
    for best insert and extractMin until have MLB implemented:
    if maxC < approx 4000,
       the dial algorithm has insert O(1) and extractMin O(1)
    else 
   
     * For large N, the XFastTrie by itself (plus an external
     * hashmao to store HeapNodes) is a better result because
     * the term O(w-l) will be smaller.
     * 
     * For mid to small N, the fibonacci heap has better insert and
     * extractMin performance (O(1) and O(log_2(N)), respectively). 
    */
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    private final boolean useDial;
    
    private final DoubleLinkedCircularList[] heap0;
    
    private int lastKnownMinKey0 = 0;
    
    // for use in tuning the capacity
    private int lastKnownMaxKey0 = 0;
    
    private long n0 = 0;
    
    private final Heap heap1;
    
    void printLastKnownMinMax() {
        log.info("min=" + lastKnownMinKey0
            + " max=" + lastKnownMaxKey0);
    }
    
    public MinHeapForRT2012(int capacity) {
        
        // using OO composition instead of specialization
        
        if (capacity < 46300) {
            
            useDial = true;
        
            heap0 = new DoubleLinkedCircularList[capacity];
            
            heap1 = null;

            log.fine("useDial=" + useDial + " cap=" + capacity);
             
        } else {
            
            useDial = false;
        
            heap0 = null;
            
            heap1 = new Heap();
        }
    }
    
    public void insert(PathNode node) {
        
        if (useDial) {
            insert0(node);
        } else {
            insert1(node);
        }
    }
    
    public PathNode extractMin() {
        
        if (useDial) {
            return extractMin0();
        } else {
            return extractMin1();
        }
    }
    
    private PathNode extractMin0() {
    
        for (int i = lastKnownMinKey0; i < heap0.length; ++i) {
            DoubleLinkedCircularList bucket = heap0[i];
            if (bucket != null && (bucket.getNumberOfNodes() > 0)) {
                HeapNode node = bucket.getSentinel().getLeft();
                bucket.remove(node);
                lastKnownMinKey0 = i;
                n0--;
                return (PathNode)node;
            }
        }
        
        return null;
    }
    
    private PathNode extractMin1() {
        
        HeapNode node = heap1.extractMin();
        if (node != null) {
            return (PathNode)node;
        } else {
            return null;
        }
    }
    
    private void insert0(PathNode node) {
         
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
    
    private void insert1(PathNode node) {
         
        int key = (int)node.getKey();
        
        heap1.insert(node);
        
        log.fine("insert into minHeap at key =" + key);        
    }
    
    public void decreaseKey(PathNode node, long key2) {
    
        if (useDial) {
            decreaseKey0(node, key2);
        } else {
            decreaseKey1(node, key2);
        }
    }
     
    private void decreaseKey0(PathNode node, long key2) {

        log.fine("decreaseKey in minHeap from key=" + 
            node.getKey() + " to key=" + key2);
        
        int prevKey = (int)node.getKey();
        heap0[prevKey].remove(node);
        
        node.setKey(key2);
        
        insert0(node);
    }
    
    private void decreaseKey1(PathNode node, long key2) {

        log.fine("decreaseKey in fibHeap from key=" + 
            node.getKey() + " to key=" + key2);
        
        heap1.decreaseKey(node, key2);
    }
    
    public long getNumberOfNodes() {
        if (useDial) {
            return n0;
        } else {
            return heap1.getNumberOfNodes();
        }
    }
}
