package algorithms;

import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.util.Map.Entry;
import java.util.TreeMap;
import thirdparty.ods.Integerizer;
import thirdparty.ods.XFastTrie;
import thirdparty.ods.XFastTrieNode;

/**
 * NOTE: NOT READY FOR USE.  
 * 
 * from wikipedia
 *     https://en.wikipedia.org/wiki/Y-fast_trie
 *  
 * a y-fast trie is a data structure for storing 
 * integers from a bounded domain. It supports exact and predecessor 
 * or successor queries in time O(log log M), using O(n) space, 
 * where n is the number of stored values and M is the maximum 
 * value in the domain. 
 * The structure was proposed by Dan Willard in 1982[1] to decrease 
 * the O(n log M) space used by an x-fast trie.
   
   The Y-Fast trie has the ordered associative array operations + successor and
   predecessor.
   
   NOTE that the runtime complexities listed are not yet achieved.
   The current operations depend upon settings, but should be 
   constant runtime complexity .lte. O(10).

   Find(k): find the value associated with the given key.
       runtime complexity is O(log log(M))
   Successor(k): find the key/value pair with the smallest key larger than or 
       equal to the given key.
       runtime complexity is O(log log(M))
   Predecessor(k): find the key/value pair with the largest key less than or 
       equal to the given key.
       runtime complexity is O(log log(M))
   Insert(k, v): insert the given key/value pair.
       runtime complexity is O(log log(M))
   Delete(k): remove the key/value pair with the given key.
       runtime complexity is O(log log(M))
 
 * Note, have not read the Willard paper yet, just a few online
 * lecture notes to implement this.  
 * 
 * the performance of each binary search tree is at most O(log_2(maxC/w))
 * where maxC is the maximum value to be stored.
 * 
 * @author nichole
 */
public class YFastTrie {

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
    
    private int nBins;

    private final XFastTrie<XFastTrieNode<Integer>, Integer> xft;
    
    // key = bin index, value = repr value.
    // each repr value is the minimum stored in the bin.
    private final TIntIntMap xftReps = new TIntIntHashMap();
    
    // there are w items in rbs
    // each list item is a sorted binary search tree of numbers in that bin.
    //    the value is the tree holds the number of times that number is present.
    private final TIntObjectMap<TreeMap<Integer, Integer>> rbs;

    private boolean chooseByN = true;

    public YFastTrie(int wBits) {
        
        if (wBits < 31 && wBits > 1) {
            this.w = wBits;
        } else {
            throw new IllegalStateException("wBits "
                + " shoulw be greater than 1 and less than 32");
        }
        maxC = (1 << w) - 1;
        
        /*
        LG binsize = (int)Math.ceil((float)maxC/(float)w);
        MID binsize = 10
        
        w,     binsz,     n,             rt
        31,    69273666,  31,            26.   LG
        31,    5.0,       429496729.6,   2.32  MID
        
        10,    102,       10,            6.67  LG
        10,    4.0,       256.0,         2.0   MID
        
        aiming for a runtime of about O(10) or better without increasing n too
        much.
        
        alternatively, could keep n as large as possible within good
        performance range, hence the binsz will be small, hence the
        runtime will be small.
        
        n            TreeMaps
        n * binSz    objects in TreeMaps
        
        heapsize: 100's or more MB
        */
        
        if (chooseByN) {
            binSz = chooseBinSizeByN();
        } else {
            int tmpLg = (int) Math.ceil((float) maxC / (float) w);
            if ((Math.log(tmpLg) / Math.log(2)) < 10) {
                binSz = tmpLg;
            } else {
                binSz = 1024;
            }
        }
        
        nBins = (int)Math.ceil((float)maxC/(float)binSz);
                
        System.out.println("nBins=" + nBins + "  rt of ops=" +
            (Math.log(binSz)/Math.log(2)));
        
        rbs = new TIntObjectHashMap<TreeMap<Integer, Integer>>();
         
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
        
        this.w = 30;
        
        maxC = (1 << w) - 1;
        
        if (chooseByN) {
            binSz = chooseBinSizeByN();
        } else {
            int tmpLg = (int) Math.ceil((float) maxC / (float) w);
            if ((Math.log(tmpLg) / Math.log(2)) < 10) {
                binSz = tmpLg;
            } else {
                binSz = 1024;
            }
        }
        
        nBins = (int)Math.ceil((float)maxC/(float)binSz);
                
        System.out.println("nBins=" + nBins + "  rt of ops=" +
            (Math.log(binSz)/Math.log(2)));
        
        rbs = new TIntObjectHashMap<TreeMap<Integer, Integer>>();
        
        XFastTrieNode<Integer> clsNode = new XFastTrieNode<Integer>();
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
        
        xft = new XFastTrie<XFastTrieNode<Integer>, Integer>(clsNode, it, w);
    }
    
    protected TreeMap<Integer, Integer> getTreeMap(int index) {
        Integer key = Integer.valueOf(index);
        TreeMap<Integer, Integer> map = rbs.get(key);
        if (map == null) {
            map = new TreeMap<Integer, Integer>();
            rbs.put(key, map);
        }
        return map;
    }

    /**
     * 
     * @param node
     * @param index 
     */
    private void addToRBTree(int node, int index) {
        
        TreeMap<Integer, Integer> map = getTreeMap(index);
        
        assert(map != null);
        
        Integer key = Integer.valueOf(node);
        
        Integer multiplicity = map.get(key);
    
        if (multiplicity == null) {
            multiplicity = Integer.valueOf(1);
        } else {
            multiplicity = Integer.valueOf(1 + multiplicity.intValue());
        }
        
        map.put(key, multiplicity);        
    }
    
    /**
     * 
     * @param node
     * @param index 
     */
    private boolean deleteFromRBTree(int node, int index) {
                
        TreeMap<Integer, Integer> map = getTreeMap(index);
        
        assert(map != null);
        
        Integer key = Integer.valueOf(node);

        Integer multiplicity = map.get(key);
    
        if (multiplicity == null) {
            return false;
        }
        
        if (multiplicity.intValue() > 0) {
            multiplicity = Integer.valueOf(multiplicity.intValue() - 1);
            if (multiplicity.intValue() > 0) {
                map.put(key, multiplicity);
            }
        }
        if (multiplicity.intValue() == 0) {
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
     * @param node a number >= 0 and having bit length 
     * less than or equal to w.
     * @return 
     */
    public boolean add(int node) {

        if (node < 0) {
            throw new IllegalArgumentException("node must "
                + "be greater than or equal to 0");
        } else if (node > maxC) {
            throw new IllegalArgumentException("node.key must "
                + "be less than " + maxC + " node=" + node);
        }
        
        int index = node/binSz;
        
        int existingRepr = xftReps.get(index);
                
        if (!xftReps.containsKey(index)) {
            // insert is O(log_2(w)) + O(l-w)
            xft.add(Integer.valueOf(node));
            xftReps.put(index, node);
        } else if (node < existingRepr) {
            // delete is O(log_2(w)) + O(l-w)
            // insert is O(log_2(w)) + O(l-w)
            xft.remove(Integer.valueOf(existingRepr));
            xft.add(Integer.valueOf(node));
            xftReps.put(index, node);
        }
                
        addToRBTree(node, index);
        
        n++;
        
        return true;
    }

    /**
     * 
     * @param node
     * @return 
     */
    public boolean remove(int node) {
        
        if (node < 0) {
            throw new IllegalArgumentException("node must "
                + "be greater than or equal to 0");
        } else if (node > maxC) {
            throw new IllegalArgumentException("node must "
                + "be less than " + maxC);
        }
        
        int index = node/binSz;
                
        boolean removed = deleteFromRBTree(node, index);
                
        if (!removed) {
            return false;
        }
        
        if (!xftReps.containsKey(index)) {
            return false;
        }
        
        TreeMap<Integer, Integer> map = getTreeMap(index);
      
        int existingRepr = xftReps.get(index);
      
        if (map.isEmpty()) {
            // just deleted the last item so remove from rbs
            // delete is O(log_2(w)) + O(w-l)
            if (xftReps.containsKey(index)) {
                xft.remove(Integer.valueOf(existingRepr));
                xftReps.remove(index);
            }
        } else if (node == existingRepr) {
            
            //existingRepr is maintained as the minimum in the bin,
            //   so if a node w/ this value is removed and the multiplicity
            //      was 1, need to assign a new repr
            
            // O(log_2(N/w))
            Integer multiplicity = map.get(Integer.valueOf(node));
    
            if (multiplicity == null) {
                // remove the current repr and assign a new one
                // delete is O(log_2(w)) + O(w-l)
                xft.remove(Integer.valueOf(existingRepr));
                xftReps.remove(index);
            
                // O(log_2(N/w))
                Entry<Integer, Integer> minEntry = map.firstEntry(); 
                xft.add(minEntry.getKey());
                xftReps.put(index, minEntry.getKey()); 
            }            
        }
        
        n--;
        
        return true;
    }

    /**
     * @param node
     * @return 
     */
    public int find(int node) {
                
        if (node < 0) {
            throw new IllegalArgumentException("node must "
                + "be greater than or equal to 0");
        } else if (node > maxC) {
            throw new IllegalArgumentException("node must "
                + "be less than " + maxC);
        }
                
        int index = node/binSz;
                
        TreeMap<Integer, Integer> map = getTreeMap(index);
        
        Integer multiplicity = map.get(Integer.valueOf(node));
        if (multiplicity == null) {
            return -1;
        }
        
        return node;
    }

    /**
     * @param node
     * @return value preceding node, else -1 if there is not one
     */
    public int predecessor(int node) {
    
        if (node < 0) {
            throw new IllegalArgumentException("node must "
                + "be greater than or equal to 0");
        } else if (node > maxC) {
            throw new IllegalArgumentException("node must "
                + "be less than " + maxC);
        }
        
        Integer nodeKey = Integer.valueOf(node);
        
        int nodeIndex = node/binSz;
        
        boolean isAMinimum = xft.find(nodeKey) != null;
        
        /*
        if the node is not a minima, the answer is in
           the node's map if its size is larger > 1
        */
        if (!isAMinimum && (rbs.get(nodeIndex).size() > 1)) {
        
            TreeMap<Integer, Integer> map = getTreeMap(nodeIndex);
          
            Entry<Integer, Integer> pred = map.lowerEntry(nodeKey);
            if (pred == null) {
                return -1;
            }
            
            return pred.getKey().intValue();
        }
       
        // else, predeccessor is in the closest bin < nodeIndex that has
        //    items in it.
                
        Integer prev = xft.predecessor(nodeKey);
        if (prev == null) {
            return -1;
        }
        
        int prev0Index = prev.intValue()/binSz;
            
        TreeMap<Integer, Integer> map = getTreeMap(prev0Index);
        
        Entry<Integer, Integer> lastItem = map.lastEntry();
               
        if (lastItem == null) {
            return -1;
        }
        
        return lastItem.getKey();
    }
    
    public int successor(int node) {
                
        if (node < 0) {
            throw new IllegalArgumentException("node must "
                + "be greater than or equal to 0");
        } else if (node > maxC) {
            throw new IllegalArgumentException("node must "
                + "be less than " + maxC);
        }
        
        Integer nodeKey = Integer.valueOf(node);
        
        int nodeIndex = node/binSz;
        
        boolean isAMinimum = xft.find(nodeKey) != null;
        
        if (isAMinimum) {
            // if tree size > 1, the next key is the successor
            // else, the xft sucessor to nodeIndex is the successor
            
            TreeMap<Integer, Integer> nodeMap = getTreeMap(nodeIndex);
            
            if (nodeMap.size() > 1) {
                Entry<Integer, Integer> successor = nodeMap.higherEntry(nodeKey);
                assert(successor != null);
                return successor.getKey();
            }
            
            Integer successorRepr = xft.successor(nodeKey);
            if (successorRepr == null) {
                return -1;
            }
            
            // the successor representative is then the next value
            return successorRepr;
        }
        
        // else, the node is not a repr
        //   if there is a tree successor to the node, that is the successor
        //   else, the xft successor to nodeIndex is the successor
        
        TreeMap<Integer, Integer> nodeMap = getTreeMap(nodeIndex);
            
        Entry<Integer, Integer> sEntry = nodeMap.higherEntry(nodeKey);
        
        if (sEntry != null) {
            return sEntry.getKey();
        }
        
        Integer successorRepr = xft.successor(nodeKey);
        if (successorRepr == null) {
            return -1;
        }

        // the successor representative is then the next value
        return successorRepr;
    }

    /**
     * runtime complexity is O(log_2(w)) 
     * @return minimum, else -1 if empty
     */
    public int minimum() {
        
        if (xft.size() == 0) {
            return -1;
        }
        
        Integer repr = xft.minimum();
        
        assert(repr != null);
       
        return repr.intValue();
    }

    /**
     * @return maximum, else -1 if empty
     */
    public int maximum() {
        
        if (xft.size() == 0) {
            return -1;
        }
        
        Integer maxRepr = xft.maximum();
        
        assert(maxRepr != null);
        
        int index = maxRepr.intValue()/binSz;
        
        TreeMap<Integer, Integer> map = getTreeMap(index);
        
        assert(map != null);
        
        Entry<Integer, Integer> lastItem = map.lastEntry();
        
        assert(lastItem != null);
        
        return lastItem.getKey();
    }
    
    /**
     * TODO: calc runtime complexity again
     * 
     * @return minumum, else -1 if empty
     */
    public int extractMinimum() {
        
        //O(log_2(w))
        int min = minimum();

        if (min == -1) {
            assert(xft.size() == 0);
            return -1;
        }
                
        remove(min);
        
        return min;
    }
    
    /**
     * TODO: calc runtime complexity again
     * 
     * @return maximum, else -1 if empty
     */
    public int extractMaximum() {
        
        int max = maximum();

        if (max == -1) {
            assert(xft.size() == 0);
            return -1;
        }
                
        remove(max);
        
        return max;
    }
    
    public int size() {
        return n;
    }

    private int chooseBinSizeByN() {
        
        long totalMemory = Runtime.getRuntime().totalMemory();
        MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
        long heapUsage = mbean.getHeapMemoryUsage().getUsed();
        long avail = totalMemory - heapUsage;

        long n = avail/32;
        
        // n = maxC/binsz
        int bs = (int)(maxC/n);
        
        double rt = Math.log(bs)/Math.log(2);
        
        if (rt < 10) {
            if (bs > 10) {
                //this is the number of items, that is capacity, of each map
                return bs;
            }
        }
        
        // else, fall back to using the default for rt = O(10)
        
        return (int)Math.ceil((float)maxC/(float)w);
        
        /*
        LG binsize = (int)Math.ceil((float)maxC/(float)w);
        MID binsize = 10
        
        w,     binsz,     n,             rt
        31,    69273666,  31,            26.   LG
        31,    5.0,       429496729.6,   2.32  MID
        
        10,    102,       10,            6.67  LG
        10,    4.0,       256.0,         2.0   MID
        
        aiming for a runtime of about O(10) or better without increasing n too
        much.
        
        alternatively, could keep n as large as possible within good
        performance range, hence the binsz will be small, hence the
        runtime will be small.
        
        n            TreeMaps
        n * binSz    objects in TreeMaps
        
        heapsize: 100's or more MB
        */
        
    }
    
    protected int getBinSz() {
        return binSz;
    }
}
