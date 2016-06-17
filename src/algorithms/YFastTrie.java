package algorithms;

/**
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
        in the XFastTrie of this yFastTrie.
      - because the xfast trie only holds w xft values,
         the space complexity is reduced.
-------------------------------
YFastTrie
-------------------------------
-w : int // num bits
-maxC : int // max vlue corresponding to w
-binSize : int // maxC/w
-bt : XFastTrie
-rbKeys : int[w]
-rbs : List<TreeMap<Integer, LinkedList<Hea[Node>>> // size is w
--------------------------------
+insert(HeapNode) : boolean
+extractMin(HeapNode) : HeapNode
+extractMax(HeapNode) : HeapNode
+successor(int key) : HeapNpde
+predecessor(int key) : HeapNpde
+find(int key) : HeaNode
+minimum() : HeaNode
+maximum() : HeaNode
+delete(HeapNode) : boolean
+size() : int
---------------------------------
note: find, min, max, successor are O(log log maxC)
      insert, delete are O(1)
---------------------------------
    
yfast trie
   - w bits set by maximum expected value to be added.
   - one xfast trie to hold the representives (at most w in number)
   - w red black trees to keep ordered points.
     - because some of the items added may have more than
       one with same key value, the values in teh red black tree
       will be linked lists.
   - choosing representatives:
      - will choose the min within each tree (min heap is current use case)
      - this needs to be updated on insert and on delete
        because the remaining queries rely on these being minimima
   - xfasttrie node:
      - rbIndex is needed as additional field OR
        (that may require edits to sfastttrie to
         make sure data is copied when new nodes are created...removed most of those).
      - an array of size w rb key values is needed to check whether
        a representative exists quickly.
   - insert node:
        - the rb tree is found by index = maxC/w (no query to xft)
        - if there is no representative for the tree,
          one is added,
          else if key is smaller than existing repr, replace existing
        - node is added to rb tree at index key/w
   - delete node:
        - the rb tree is found by index = maxC/w (no query to xft)
        - remove key from rb tree.
        - if its linked list is empty, remove the key from rb tree.
        - if this key is the repr for the tree,
          remove the repr from xdt and array.
          - note that if this is the case, there are no more keys
            in the rb tree.
            assert that.
    - find value:
       - search for exact value in xft
            and if found, return node from found index
       - else search the repr which are all min values of their trees.
       - find the largest min value (==repr) smaller than value.
         search the repr rbIndex tree.
   - extractMin:
     - find min of the repr keys:
        xft.minimum gives rbIndex and node value.  get the node from
        the rb tree and delete it
   - successor value:
       - find successor repr in xfastrie:
          r    r    r
         |    |    |    |
         - rbIndex is the tree the successor should be
           in unless the value is the maximum of the tree
           in which case the successor is the min value of next tree,
           else is successor of value in current rb tree.
   - predecessor value:
       - find predecessor repr in xfastrie:
          r    r    r
         |    |    |    |
         - if value is the minimum value in rb tree,
             return max value of preceding populated rb trr,
           else, return predecessor of value within the rb tree.
    puased here...
    */
}
