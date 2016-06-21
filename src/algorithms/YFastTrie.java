package algorithms;

import algorithms.imageProcessing.HeapNode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import thirdparty.ods.Integerizer;
import thirdparty.ods.XFastTrie;
import thirdparty.ods.XFastTrieNode;
import thirdparty.ods.SSet;

/**
 * NOTE: NOT READY FOR USE.  Until the red-black trees internally
 * are replaced by a data structure that has operations
 * faster than O(log_2(N_number_of_nodes)), the XFastTrie by itself
 * is a better choice (along with a supplemental hashmap to store
 * nodes).
 * 
 * Note, have not read the Willard paper yet, just a few online
 * lecture notes to implement this.  A couple suggest that the
 * performance of each red black tree is O(log_2(maxC)), but that
 * would require the partitions to be based on the number of points
 * in the maps and that would quickly be very many tree maps for
 * a large maxC (each with number of points being w).
 * 
 * @author nichole
 */
public class YFastTrie implements SSet<HeapNode>{

    /*    
    designing from browsing a few different lecture notes
    online. the yfast trie uses same w and maxC as
    the XFastTrie.
      - creates w red black trees to hold inserted heap nodes.
        the w trees each have range size of maxC/w and
        start from 0 extending to last one holding maxC.
      - each tree has a representative if it has any nodes
        and those are stored
        in the XFastTrie of this YFastTrie.
      - because the XFastTrie only holds w xft values,
        the space complexity is reduced.
-------------------------------
YFastTrie

   - w bits set by maximum expected value to be added.
   - one XFastTrie to hold the representives (at most w in number)
   - w red black trees to keep ordered points.
     - because some of the items added may have more than
       one with same key value, the values in the red black tree
       will be linked lists.

    NOTE: topics to consider for improvements:
          the distribution of rb trees, that is their partitions,
          could be improved dynamically.
          For example, if maxC were value 127, but the majority
          of nodes at some point in time were in bin 0 at values
          near 4, one would prefer to divide that tree 
          into more than one tree to speed up searches.
          This begins to look like a good reason to
          compare to multi-level-buckets.  The only implementation
          I could find was the Andrew Goldberg MLB offered 
          under a license that is freely available for non-commercial
          use, else need to contact for permission... not wanting
          to include mixed license restrictions for now...
          (so I didn't download and read the code.  am reading
          his 2 papers on the subject, but they depend upon other
          papers too, so gathering all the specs for the MLB
          algorithms is not complete...)
          -- one possible work around without making dynamic
          partitions in the YFastTrie would be to know or 
          estimate the population of data ahead of time and 
          then make separate YFastTrie's for manually partitioned 
          data (changing zero-points, a.k.a. bias levels as needed
          before and after use of more than one YFastTrie)
    */
    
    private int n = 0;
    
    private final int w;
    
    private final int maxC;
    
    private final int binSz;

    private final XFastTrie<XFastTrieNode<Integer>, Integer> xft;
    
    private final Map<Integer, HeapNode> xftReps = 
        new HashMap<Integer, HeapNode>();
    
    private final List<TreeMap<Integer, LinkedList<HeapNode>>> rbs;
    
    public YFastTrie(int wBits) {
        if (wBits < 32 && wBits > 1) {
            this.w = wBits;
        } else {
            throw new IllegalStateException("wBits "
                + " shoulw be greater than 1 and less than 32");
        }
        maxC = (1 << (w - 1)) - 1;
        // NOTE: if change out TreeMap to a data structure that
        //       is faster than O(w) and does not depend on N,
        //       this may change, but for now, the number of bins
        //       will be w.
        //       the TreeMap operations for evenly distributed
        //       data are currently then O(log_2(N/w)).
        binSz = (int)Math.ceil((float)maxC/(float)w);
        rbs = new ArrayList<TreeMap<Integer, LinkedList<HeapNode>>>(w);
        for (int i = 0; i < w; ++i) {
            rbs.add(new TreeMap<Integer, LinkedList<HeapNode>>());
        }
        
        XFastTrieNode<Integer> clsNode = new XFastTrieNode<Integer>();
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
        
        xft = new XFastTrie<XFastTrieNode<Integer>, Integer>(clsNode, it, w);
    }
    
    public YFastTrie() {
        
        this.w = 32;
        
        maxC = (1 << (w - 1)) - 1;
        // NOTE: if change out TreeMap to a data structure that
        //       is faster than O(w) and does not depend on N,
        //       this may change, but for now, the number of bins
        //       will be w.
        //       the TreeMap operations for evenly distributed
        //       data are currently then O(log_2(N/w)).
        binSz = (int)Math.ceil((float)maxC/(float)w);
        rbs = new ArrayList<TreeMap<Integer, LinkedList<HeapNode>>>(w);
        for (int i = 0; i < w; ++i) {
            rbs.add(new TreeMap<Integer, LinkedList<HeapNode>>());
        }
        
        XFastTrieNode<Integer> clsNode = new XFastTrieNode<Integer>();
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
        
        xft = new XFastTrie<XFastTrieNode<Integer>, Integer>(clsNode, it, w);
    }

    /**
     * runtime complexity is roughly O(log_2(N/w))
     * @param node
     * @param index 
     */
    private void addToRBTree(HeapNode node, int index) {
        
        TreeMap<Integer, LinkedList<HeapNode>> map =
            rbs.get(index);
        
        assert(map != null);
        
        Integer key = Integer.valueOf(Integer.valueOf((int)node.getKey()));
        
        // O(log_2(N/w))
        LinkedList<HeapNode> list = map.get(key);
    
        if (list == null) {
            list = new LinkedList<HeapNode>();
            // O(log_2(N/w))
            map.put(key, list);
        }
        
        list.add(node);
    }
    
    /**
     * runtime complexity is roughly O(log_2(N/w))
     * @param node
     * @param index 
     */
    private boolean deleteFromRBTree(HeapNode node, int index) {
        
        TreeMap<Integer, LinkedList<HeapNode>> map =
            rbs.get(index);
        
        assert(map != null);
        
        Integer key = Integer.valueOf((int)node.getKey());

        // O(log_2(N/w))
        LinkedList<HeapNode> list = map.get(key);
    
        if (list == null) {
            return false;
        }
        
        //O(1)
        list.remove(node);
        
        if (list.size() == 0) {
            // O(log_2(N/w))
            map.remove(key);
        }
        
        return true;
    }
    
    /**
     * runtime complexity is roughly
     *    O(log_2(w)) + O(w-l) + O(log_2(N/w))
     * where N is total number of nodes in the YFastTrie,
     * w is the maximum value possible in bit length,
     * and l is the number of levels in the prefix trie 
     * already filled with other entries.
     * 
     * NOTE:
     * For small maxC, the "Dial algorithm" has best insert and
     * extractMin runtime complexities (O(1) and O(1), respectively).
     * 
     * For large N, the XFastTrie by itself (plus an external
     * hashmao to store HeapNodes) is a better choice because
     * the term O(w-l) will be smaller.
     * 
     * For mid to small N, the fibonacci heap has better insert and
     * extractMin performance (O(1) and O(log_2(N)), respectively).
     * 
     * @param node a heap node with key >= 0 and having bit length 
     * less than or equal to w.
     * @return 
     */
    @Override
    public boolean add(HeapNode node) {

        if (node.getKey() < 0) {
            throw new IllegalArgumentException("node.key must "
                + "be greater than or equal to 0");
        } else if (node.getKey() > maxC) {
            throw new IllegalArgumentException("node.key must "
                + "be less than " + maxC);
        }
        
        int index = (int)node.getKey()/binSz;
        
        HeapNode existingRepr = xftReps.get(Integer.valueOf(index));
        
        if (existingRepr == null) {
            // insert is O(log_2(w)) + O(l-w)
            xft.add(Integer.valueOf((int)node.getKey()));
            xftReps.put(Integer.valueOf(index), node);
        } else if (node.getKey() < existingRepr.getKey()) {
            // delete is O(log_2(w)) + O(l-w)
            // insert is O(log_2(w)) + O(l-w)
            xft.remove(Integer.valueOf((int)existingRepr.getKey()));
            xft.add(Integer.valueOf((int)node.getKey()));
            xftReps.put(Integer.valueOf(index), node);
        }
        
        // O(log_2(N/w))
        addToRBTree(node, index);
        
        n++;
        
        return true;
    }

    /**
     * runtime complexity is roughly
     *    O(log_2(w)) + O(w-l) + O(log_2(N/w))
     * where N is total number of nodes in the YFastTrie,
     * w is the maximum value possible in bit length,
     * and l is the number of levels in the prefix trie 
     * already filled with other entries.
     * 
     * @param node
     * @return 
     */
    @Override
    public boolean remove(HeapNode node) {
        
        if (node.getKey() < 0) {
            throw new IllegalArgumentException("node.key must "
                + "be greater than or equal to 0");
        } else if (node.getKey() > maxC) {
            throw new IllegalArgumentException("node.key must "
                + "be less than " + maxC);
        }
        
        int index = (int)node.getKey()/binSz;
        
        // O(log_2(N/w))
        boolean removed = deleteFromRBTree(node, index);
        
        if (!removed) {
            return false;
        }
        
        TreeMap<Integer, LinkedList<HeapNode>> map =
            rbs.get(index);
        
        HeapNode existingRepr = xftReps.get(Integer.valueOf(index));
        
        if (node.getKey() == existingRepr.getKey()) {
            if (map.isEmpty()) {
                // delete is O(log_2(w)) + O(w-l)
                xft.remove(Integer.valueOf((int)existingRepr.getKey()));
                xftReps.remove(Integer.valueOf(index));
            } else {
                // O(log_2(N/w))
                Entry<Integer, LinkedList<HeapNode>> 
                    entry = map.firstEntry();                
                LinkedList<HeapNode> list = entry.getValue();
                
                HeapNode node2 = list.getFirst();
                int key2 = (int)node2.getKey();
                // delete is O(log_2(w)) + O(w-l)
                // insert is O(log_2(w)) + O(w-l)
                xft.remove(Integer.valueOf((int)existingRepr.getKey()));
                xft.add(Integer.valueOf(key2));
                xftReps.put(Integer.valueOf(index), node2); 
            }
        }
        
        n--;
        
        return true;
    }

    /**
     * runtime complexity is roughly O(log_2(N/w))
     * @param node
     * @return 
     */
    @Override
    public HeapNode find(HeapNode node) {
        
        //TODO: revist to improve runtimes
        
        if (node.getKey() < 0) {
            throw new IllegalArgumentException("node.key must "
                + "be greater than or equal to 0");
        } else if (node.getKey() > maxC) {
            throw new IllegalArgumentException("node.key must "
                + "be less than " + maxC);
        }
        
        int key = (int)node.getKey();
        
        int index = (int)node.getKey()/binSz;
                
        TreeMap<Integer, LinkedList<HeapNode>> map =
            rbs.get(index);
        
        // O(log_2(N/w))
        LinkedList<HeapNode> list = map.get(Integer.valueOf(key));
        if (list == null) {
            return null;
        }
        
        return list.getFirst();
    }

    /**
     * runtime complexity is O(log_2(w)) + O(log_2(N/w)).
     * @param node
     * @return 
     */
    @Override
    public HeapNode predecessor(HeapNode node) {
    
        if (node.getKey() < 0) {
            throw new IllegalArgumentException("node.key must "
                + "be greater than or equal to 0");
        } else if (node.getKey() > maxC) {
            throw new IllegalArgumentException("node.key must "
                + "be less than " + maxC);
        }
        
        return predecessor((int)node.getKey());
    }
    
    /**
     * runtime complexity is O(log_2(w)) + O(log_2(N/w)).
     * @param nodeKeyIdx
     * @return 
     */
    public HeapNode predecessor(int nodeKeyIdx) {

        //TODO: revisit to reduce runtime complexity
        
        if (nodeKeyIdx < 0) {
            throw new IllegalArgumentException("node.key must "
                + "be greater than or equal to 0");
        } else if (nodeKeyIdx > maxC) {
            throw new IllegalArgumentException("node.key must "
                + "be less than " + maxC);
        }
        
        Integer nodeKey = Integer.valueOf(nodeKeyIdx);
        
        int nodeIndex = nodeKeyIdx/binSz;
        
        boolean isAMinimum = xft.find(nodeKey) != null;
        
        /*
        if the node is not a minima, the answer is in
           the node's map if its size is larger > 1
        */
        if (!isAMinimum && (rbs.get(nodeIndex).size() > 1)) {
        
            TreeMap<Integer, LinkedList<HeapNode>> map =
                rbs.get(nodeIndex);
          
            // O(log_2(N/w))
            return map.lowerEntry(nodeKey).getValue().getFirst();
        }
        
        //O(log_2(w))
        Integer prev = xft.predecessor(nodeKey);
        if (prev == null) {
            return null;
        }
        
        int prev0Index = prev.intValue()/binSz;
            
        TreeMap<Integer, LinkedList<HeapNode>> map =
            rbs.get(prev0Index);
          
        // O(log_2(N/w))
        Entry<Integer, LinkedList<HeapNode>> list =
            map.lowerEntry(nodeKey);
        
        if (list == null) {
            return null;
        }
        
        return list.getValue().getFirst();
    }
    
    @Override
    public HeapNode successor(HeapNode node) {
                
        if (node.getKey() < 0) {
            throw new IllegalArgumentException("node.key must "
                + "be greater than or equal to 0");
        } else if (node.getKey() > maxC) {
            throw new IllegalArgumentException("node.key must "
                + "be less than " + maxC);
        }
        
        return successor((int)node.getKey());
    }

    /**
     * runtime complexity is roughly O(log_2(w)) + O(log_2(N/w))
     * @param nodeKeyIdx
     * @return 
     */
    public HeapNode successor(final int nodeKeyIdx) {
        
        if (nodeKeyIdx < 0) {
            throw new IllegalArgumentException("node.key must "
                + "be greater than or equal to 0");
        } else if (nodeKeyIdx > maxC) {
            throw new IllegalArgumentException("node.key must "
                + "be less than " + maxC);
        }
        
        Integer nodeKey = Integer.valueOf(nodeKeyIdx);
        
        // because the representatives are minima for their
        // trees, the successor query to xft would not necessarily
        // find the tree that node is in.
        // 
        // if node is the max within its tree,
        //    then xft.successor is the successor
        //    because its the minimim of next populated tree.
        // else
        //    the node's map size is > 1 and
        //    the successor is learned from the query on the tree map.
 
        int nodeIndex = nodeKey.intValue()/binSz;
        
        HeapNode existingRepr = xftReps.get(Integer.valueOf(nodeIndex));

        if ((existingRepr != null) &&
            (existingRepr.getKey() == nodeKey.intValue()) &&
            (rbs.get(nodeIndex).size() > 1)) {
            
            // if map size is > 1, the answer is in this map
            // they're in the same tree so the answer is too
            
            TreeMap<Integer, LinkedList<HeapNode>> nodeMap =
                rbs.get(nodeIndex);
            
            // O(log_2(N/w))
            Entry<Integer, LinkedList<HeapNode>> entry 
                = nodeMap.higherEntry(nodeKey);
            
            return entry.getValue().getFirst();
        }
        
        Integer next0 = xft.successor(nodeKey);
        
        if (next0 == null) {
            return null;
        }
        
        int next0Index = next0.intValue()/binSz;
        
        // if nodeKey is the last that could be placed
        // in it's bin, then the answer must be in the
        // successor tree.
       
        Integer prev0 = xft.predecessor(nodeKey);
        
        if ((prev0 == null) ||
            (nodeIndex == next0Index) ||
            ((nodeKeyIdx % binSz) - 1 == 0) || (prev0 == null)
            ) {
            
            // they're in the same tree so the answer is too
            TreeMap<Integer, LinkedList<HeapNode>> nodeMap =
                rbs.get(next0Index);
            
            // O(log_2(N/w))
            Entry<Integer, LinkedList<HeapNode>> entry 
                = nodeMap.higherEntry(nodeKey);
            
            if (entry == null) {
                return null;
            }
            
            return entry.getValue().getFirst();
        }
        
        int prev0Index = prev0.intValue()/binSz;
        
        TreeMap<Integer, LinkedList<HeapNode>> nodeMap =
            rbs.get(prev0Index);
            
        if (nodeMap.isEmpty()) {
            return null;
        }
        
        // O(log_2(N/w))
        Entry<Integer, LinkedList<HeapNode>> entry 
            = nodeMap.higherEntry(nodeKey);
          
        if (entry != null) {
            return entry.getValue().getFirst();
        }
           
        // else, must be in the successor tree if anywhere
        nodeMap = rbs.get(next0Index);

        // O(log_2(N/w))
        entry = nodeMap.higherEntry(nodeKey);

        if (entry == null) {
            return null;
        }

        return entry.getValue().getFirst();
    }

    /**
     * runtime complexity is O(log_2(w)) 
     * @return 
     */
    @Override
    public HeapNode minimum() {
        
        //O(log_2(w))
        Integer value = xft.minimum();
                
        if (value == null) {
            // no nodes in trie, hence none in trees?
            assert(xft.size() == 0);
            return null;
        }
        
        int index = value.intValue()/binSz;
        
        HeapNode minNode = xftReps.get(index);
        
        return minNode;
    }

    /**
     * runtime complexity is roughly O(log_2(w)) + O(log_2(N/w))
     * @return 
     */
    @Override
    public HeapNode maximum() {
        
        //O(log_2(w))
        Integer qIndex = xft.find(Integer.valueOf(maxC));
        
        if (qIndex != null) {

            // the maximum in last map can be returned.
            int nodeIndex = rbs.size() - 1;
                
            //TODO: if a treeMap.get() is faster than
            // a treeMap.lastEntry, change these:
            
            TreeMap<Integer, LinkedList<HeapNode>> nodeMap =
                rbs.get(nodeIndex);
        
            // O(log_2(N/w))
            Entry<Integer, LinkedList<HeapNode>> entry =
                nodeMap.lastEntry();
        
            return entry.getValue().getLast();
        }
        
        // else, xft.predecessor finds the max tree
        // then the last entry in it is the return
        
        //O(log_2(w))
        qIndex = xft.predecessor(Integer.valueOf(maxC));
        
        if (qIndex == null) {
            // no entries in trie, hence none in trees?
            assert(xft.size() == 0);
            return null;
        }
        
        int index = qIndex.intValue()/binSz;
        
        TreeMap<Integer, LinkedList<HeapNode>> map =
            rbs.get(index);
        
        // O(log_2(N/w))
        Entry<Integer, LinkedList<HeapNode>> entry =
            map.lastEntry();
        
        return entry.getValue().getLast();
    }
    
    /**
     * runtime complexity is roughly
     *    O(log_2(w)) + O(w-l) + O(log_2(N/w))
     * where N is total number of nodes in the YFastTrie,
     * w is the maximum value possible in bit length,
     * and l is the number of levels in the prefix trie 
     * already filled with other entries.
     * 
     * * NOTE:
     * For small maxC, the "Dial algorithm" has best insert and
     * extractMin runtime complexities (O(1) and O(1), respectively).
     * 
     * For large N, the XFastTrie by itself (plus an external
     * hashmao to store HeapNodes) is a better result because
     * the term O(w-l) will be smaller.
     * 
     * For mid to small N, the fibonacci heap has better insert and
     * extractMin performance (O(1) and O(log_2(N)), respectively).
     * 
     * @return 
     */
    public HeapNode extractMinimum() {
        
        //O(log_2(w))
        HeapNode min = minimum();

        if (min == null) {
            assert(xft.size() == 0);
            return null;
        }
        
        //TODO: could combine these two operations to reduce a query
        
        //O(log_2(w)) + O(w-l) + O(log_2(N/w))
        remove(min);
        
        return min;
    }
    
    /**
     * runtime complexity is roughly
     *    O(log_2(w)) + O(w-l) + O(log_2(N/w))
     * @return 
     */
    public HeapNode extractMaximum() {
        
        HeapNode max = maximum();

        if (max == null) {
            assert(xft.size() == 0);
            return null;
        }
        
        //TODO: could combine these two operations to reduce a query
        
        remove(max);
        
        return max;
    }
    
    @Override
    public int size() {
        return n;
    }
    
}
