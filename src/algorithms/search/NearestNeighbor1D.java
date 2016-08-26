package algorithms.search;

import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.HashSet;
import java.util.Set;
import thirdparty.ods.Integerizer;
import thirdparty.ods.XFastTrie;
import thirdparty.ods.XFastTrieNode;

/**
 * a nearest neighbor's algorithm using XFastTrie
 * for predecessor and successor queries
 * on spatially indexed numbers.
 * 
 * The algorithm performs better on dense data
 * (that is because the base of the prefix tree
 * is filled, leaving smaller number of nodes to
 * create in linear time).
 * The queries depend upon the maximum of x to be queried.
 * 
 * @author nichole
 */
public class NearestNeighbor1D {
    
    private final XFastTrie<XFastTrieNode<Integer>, Integer> xbt;
            
    private final int maxIndex;
    
    /**
     * 
     * @param maxX maximum x value of any data point including
     *    those to be queries
     */
    public NearestNeighbor1D(int maxX) {

        maxIndex = maxX + 1;
        
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
        
        int maxW = 1 + (int)Math.ceil(Math.log(maxIndex)/Math.log(2));
        xbt = new XFastTrie<XFastTrieNode<Integer>, Integer>(
            new XFastTrieNode<Integer>(), it, maxW);
    }
    
    public void insert(int index) {
        xbt.add(Integer.valueOf(index));
    }
    
    /**
    <pre>
      best case: O(1) if key exists in
         the trie, else 2 * O(log_2(maxW)).
            
     Note: maxW = 1 + Math.ceil(Math.log(maxX)/Math.log(2));            
     </ore>
    
     * @param x
     */
    public TIntSet findClosest(final int x) {        
        return findClosest(x, Integer.MAX_VALUE);
    }
    
    /**
    <pre>
      runtime complexity is
         best case: O(1) if key exists in
         the trie, else 2 * O(log_2(maxW)).
                     
      Note: maxW = 1 + Math.ceil(Math.log(maxX * maxY)/Math.log(2));
     </ore>
    
     * @param x
     * @param dMax
     * @return a set of values within dMax that are the 
     * closest points, else returns an empty set
     */
    public TIntSet findClosest(int x, int dMax) {
        
        if (x >= maxIndex) {
            throw new IllegalArgumentException("x cannot be larger than "
                + " maxX given in constructor, " + (maxIndex-1));
        }
         
        TIntSet closestIndexes = new TIntHashSet();
        
        Integer srch = Integer.valueOf(x);
        
        //O(1)
        Integer q = xbt.find(srch);
        if (q != null) {
            closestIndexes.add(srch.intValue());
            return closestIndexes;
        }
        
        if (dMax == 0) {
            return closestIndexes;
        }
                        
        //O(log_2(maxW))
        Integer predecessor = xbt.predecessor(srch);
        Integer successor = xbt.successor(srch);
        
        double dp2 = dist(x, predecessor);
        double ds2 = dist(x, successor);
    
        if (dp2 == ds2 && (ds2 <= dMax)) {
            closestIndexes.add(predecessor.intValue());
            closestIndexes.add(successor.intValue());
        } else if ((dp2 < ds2) && (dp2 <= dMax)) {
            closestIndexes.add(predecessor.intValue());
        } else if ((ds2 < dp2) && (ds2 <= dMax)) {
            closestIndexes.add(successor.intValue());
        }
        
        return closestIndexes;
    }

    private double dist(int x, Integer p2) {
        
        if (p2 == null) {
            return Double.MAX_VALUE;
        }
        
        int diff = p2.intValue() - x;
        if (diff < 0) {
            diff *= -1;
        }
        
        return diff;
    }

}
