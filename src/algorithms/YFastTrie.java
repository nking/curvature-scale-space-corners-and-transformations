package algorithms;

import algorithms.imageProcessing.HeapNode;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;
import thirdparty.ods.Integerizer;
import thirdparty.ods.XFastTrie;
import thirdparty.ods.XFastTrieNode;
import thirdparty.ods.SSet;

/**
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
   - choosing representatives:
      - will choose the min within each tree (min heap is current use case)
      - this needs to be updated on insert and on delete
        because the remaining queries rely on these being minimima
   - XFastTrieNode:
      - rbIndex is needed as additional field
      - an array of size w rb key values is needed to check whether
        a representative exists quickly.
   - insert node:
        - the rb tree is found by index = key/w (no query to xft)
        - if there is no representative for the tree,
          one is added,
          else if key is smaller than existing repr, replace existing
        - node is added to rb tree at index key/w
   - delete node:
        - the rb tree is found by index = key/w (no query to xft)
        - remove key from rb tree.
        - if its linked list is empty, remove the key from rb tree.
        - if this key is the repr for the tree,
          remove the repr from xdt and array.
          - note that if this is the case, there are no more keys
            in the rb tree.
            assert that.
    - find value:
       - xft.find for exact value in xft
            and if found, return node from found rb tree index and key
       - else 
          - find the largest min value (==repr) smaller than value.
            (== xft.predecessor). 
            then with the found rbtree index, search the tree.
   - extractMin:
     - find min of the repr keys:
        xft.minimum gives rbIndex and node value.  get the node from
        the rb tree and delete it
   - extractMax:
     - find max of the repr keys:
        xft.maximum gives rbIndex and node value.  get the node from
        the rb tree and delete it
   - successor value:
       - xft.successor to find repr
          r    r    r
         |    |    |    |
         - rbIndex is the tree the successor should be
           in unless the value is the maximum of the tree
           in which case the successor is the min value of next tree,
           else is successor of value in current rb tree.
   - predecessor value:
       - xft.predecessor to find repr
          r    r    r
         |    |    |    |
         - if value is the minimum value in rb tree,
             return max value of preceding populated rb tree,
           else, return predecessor of value within the rb tree.
    - minimum value:
        - xft.minimum
    - maximum value:
        - xft.find for value maxC and if found, return node
        - else,
          xft.predecessor(maxC) find the rb tree index,
          then max of that tree is the maximum
    
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
          under a license that prohibits commercial use
          (so I didn't download and read the code.  am reading
          his 2 papers on the subject, but they depend upon other
          papers too, so gathering all the specs for his
          algorithm is not complete...)
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

    private final XFastTrie<YXFTNode, Integer> xft;
    
    private final int[] rbKeys;
    
    private final List<TreeMap<Integer, LinkedList<HeapNode>>> rbs;
    
    public YFastTrie(int wBits) {
        if (wBits < 32 && wBits > 1) {
            this.w = wBits;
        } else {
            throw new IllegalStateException("wBits "
                + " shoulw be greater than 1 and less than 32");
        }
        maxC = (1 << (w - 1)) - 1;
        rbKeys = new int[w];
        rbs = new ArrayList<TreeMap<Integer, LinkedList<HeapNode>>>(w);
        for (int i = 0; i < w; ++i) {
            rbs.add(new TreeMap<Integer, LinkedList<HeapNode>>());
        }
        
        YXFTNode clsNode = new YXFTNode();
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
        
        xft = new XFastTrie<YXFTNode, Integer>(clsNode, it, w);
    }
    
    public YFastTrie() {
        
        this.w = 32;
        
        maxC = (1 << (w - 1)) - 1;
        rbKeys = new int[w];
        rbs = new ArrayList<TreeMap<Integer, LinkedList<HeapNode>>>(w);
        for (int i = 0; i < w; ++i) {
            rbs.add(new TreeMap<Integer, LinkedList<HeapNode>>());
        }
        
        YXFTNode clsNode = new YXFTNode();
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
        
        xft = new XFastTrie<YXFTNode, Integer>(clsNode, it, w);
    }

    @Override
    public HeapNode find(HeapNode x) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean add(HeapNode x) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean remove(HeapNode x) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public HeapNode predecessor(HeapNode x) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public HeapNode successor(HeapNode x) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public HeapNode minimum() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public HeapNode maximum() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public HeapNode extractMinimum() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public HeapNode extractMaximum() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    @Override
    public int size() {
        return n;
    }
    
    public class YXFTNode extends XFastTrieNode<Integer> {
        protected int rbIndex;
        public YXFTNode() {};
        public int getRBIndex() {
            return rbIndex;
        }
        public void setRBIndex(int idx) {
            this.rbIndex = idx;
        }
    }
}
