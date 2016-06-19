package algorithms.bipartite;

import algorithms.bipartite.MinCostUnbalancedAssignment.PathNode;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import thirdparty.ods.Integerizer;
import thirdparty.ods.XFastTrie;
import thirdparty.ods.XFastTrieNode;

/**
 * a wrapper for the XFastTrie to provide methods of a minimum
 * heap that can handle more than one node of the same key.
 * @author nichole
 */
public class XFastTrieWrapper {

    private final int w;
    
    private final int maxC;
    
    private final XFastTrie<XFastTrieNode<Integer>, Integer> xft;
    
    private final Map<Integer, Set<PathNode>> map =
        new HashMap<Integer, Set<PathNode>>();
    
    private long lastKnownMinKey = 0;
    private long lastKnownMaxKey = -1;
        
    private int n = 0;
    
    public XFastTrieWrapper(int maxC) {
            
        w = (int)Math.ceil(Math.log(maxC)/Math.log(2)) + 1;
        
        this.maxC = maxC;
        
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
        
        XFastTrieNode<Integer> node = new XFastTrieNode<Integer>();
        
        xft = new XFastTrie<XFastTrieNode<Integer>, Integer>(node, it, w);
    }
    
    public int getW() {
        return w;
    }
    
    /**
     * runtime complexity for best case is O(1) when there
     * is already a similar key in the XFastTrie, else is
     * is O(log_2(w)) + O(w-l)
     * where w is the number of
     * bits set in the constructor
     * and l is the prefix tree already filled leading
     * up to the value x.  The performance of the XFastTrie
     * increases when more nodes are in it (can see that in the
     * l term).
     * @param node 
     */
    public void insert(PathNode node) {
        
        Integer key = Integer.valueOf((int)node.getKey());
        
        Set<PathNode> set = map.get(key);
        
        if (set == null) {
            // O(log_2(w)) + O(w-l)
            boolean added = xft.add(key);
            assert(added);
            
            set = new HashSet<PathNode>();
            map.put(key, set);
        }
        
        //O(1)
        set.add(node);
        
        n++;
        
        if (key < lastKnownMinKey) {
            lastKnownMinKey = key;
        }
        if (key > lastKnownMaxKey) {
            lastKnownMaxKey = key;
        }
    }
    
    /**
     * runtime complexity for best case is O(1) when there
     * is already a similar key in the XFastTrie, else is
     * is O(log_2(w)) + O(w-l)
     * where w is the number of
     * bits set in the constructor
     * and l is the prefix tree already filled leading
     * up to the value x.  The performance of the XFastTrie
     * increases when more nodes are in it (can see that in the
     * l term).
     * @param node 
     */
    public void decreaseKey(PathNode node, long key2) {

        Integer key = Integer.valueOf((int)node.getKey());
        
        Set<PathNode> set0 = map.get(key);
        
        assert(set0 != null);
        
        set0.remove(node);
        
        if (set0.size() == 0) {
            boolean removed = xft.remove(key);
            assert(removed);
            map.remove(key);
        }
                        
        node.setKey(key2);
        
        Integer index2 = Integer.valueOf((int)key2);
        
        Set<PathNode> set2 = map.get(index2);
         
        if (set2 == null) {
            boolean added = xft.add(index2);
            assert(added);
            
            set2 = new HashSet<PathNode>();
            map.put(index2, set2);
        }
        
        set2.add(node);
     
        if (key < lastKnownMinKey) {
            lastKnownMinKey = key;
        }
    }
    
    /**
     * runtime complexity for best case is O(1) when there
     * is already a similar key in the XFastTrie, else is
     * is O(log_2(w)) + O(w-l)
     * where w is the number of
     * bits set in the constructor
     * and l is the prefix tree already filled leading
     * up to the value x.  The performance of the XFastTrie
     * increases when more nodes are in it (can see that in the
     * l term).
     * @return node 
     */
    public PathNode extractMin() {
        
        if (n == 0) {
            return null;
        }
     
        Integer key = xft.minimum();
        
        Set<PathNode> set = map.get(key);
        
        PathNode node = set.iterator().next();
        set.remove(node);
        
        if (set.isEmpty()) {
            map.remove(key);
            xft.remove(key);
        }
       
        lastKnownMinKey = key.intValue();
        n--;
        
        return (PathNode)node;
    }
    
    public long getNumberOfNodes() {
        return n;
    }
}
