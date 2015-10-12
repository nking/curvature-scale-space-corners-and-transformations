package algorithms.imageProcessing;

import algorithms.util.Stack;
import java.util.logging.Logger;

/**
 * Class contains a Fibonacci heap, that is, a loose collection of trees based 
 * upon binomial heaps, hence satisfying the minimum heap property:
 *     object.child.key >= object.key.
 * 
 * With a Fibonacci heap, the minimum key of the entire heap is always at the
 * top of one of the trees.
 *
 * Fibonacci heap potential = t + 2m
 *     where t = number of trees
 *           m = number of marked nodes.  (marked when node has been recently 
 *               made a child of another node and >= 1 of it's own children 
 *               have been cut.  root nodes are never marked.)
 * 
 * <pre>
 * Runtime complexity:
 * 
 *    Find-minimum is O(1) amortized time because there is always an instance
 * reference to it.
 *
 *    Insert, decrease key work in constant amortized time.
 *
 *    Delete and delete minimum work in O(log n) amortized time.
 *
 *    Extract-min and delete should be used sparingly for this structure to be
 * best utilized and are usually implemented as O(log_2 N).
 *
 * This was implemented following pseudo-code from 
 * "Introduction to Algorithms", by Cormen, Leiserson, Rivest, & Stein
 * on the Fibonacci Heap.
 * </pre>
 * 
 * @author nichole
 */
public class Heap {

    protected Logger log = Logger.getLogger(this.getClass().getName());
    
	/** circular doubly linked list of minimum nodes for their respective 
    min-heap-ordered trees */
	private DoubleLinkedCircularList rootList = new DoubleLinkedCircularList();

	/** root of tree containing a minimum key.  it's null for an empty tree */
	protected HeapNode minimumNode = null;
	
	protected long n = 0;

    /**
     * insert node into heap.  runtime is O(1).  makes no attempt to consolidate
     * tree.
     *
     * @param node
     */
    public void insert(HeapNode node) {
    	if (node.getKey() == DoubleLinkedCircularList.noValue) {
            throw new IllegalArgumentException(
                "node.key must be set before insert into heap");
        }
    	node.setNumberOfChildren(0);
        node.setParent(null);
        node.removeChildren();
        node.setLeft(node);
        node.setRight(node);
        node.setMark(false);

        // concatenate root list containing node with this.rootList
        rootList.insert(node);
            
        if ((minimumNode == null) || (node.getKey() < minimumNode.getKey())) {
            minimumNode = node;
        }
        
        n++;
    }

    public boolean isEmpty() {
        return (minimumNode == null);
    }
    
    /**
     * extract minimum from the heap.
     *
     * @return
     */
    public HeapNode extractMin() {
        
        int sentinel = DoubleLinkedCircularList.sentinelKey;
        
    	HeapNode z = minimumNode;

        if (z == null) {
            return z;
        }

        //save reference to right of minimum node
        HeapNode nextMin = z.getRight();
        
        // detach each child and add it to heap
        HeapNode x = z.getChildren().getSentinel().getRight();
          
        // for each child x of z
        while (x.getKey() != sentinel) {
            HeapNode next = x.getRight();
            x.setParent(null);
            rootList.insert(x);
            x = next;
        }

        rootList.remove(z);
    
        if (z.equals(nextMin)) {
            minimumNode = null;
        } else {
            minimumNode = nextMin;
            consolidate();
        }

        n--;
            
        return z;
    }

    void consolidate() {
        
    	// D[n] = max degree of any node = lg_2(n) = lg_2(Integer.MAX) = 31
        //int maxDegree = (int) (Math.log(this.n)/Math.log(2));
        int maxDegree = 31;

        HeapNode[] a = new HeapNode[maxDegree];

        HeapNode w = rootList.getSentinel().getRight();

        // n*m*(constants)
        while (w.getKey() != DoubleLinkedCircularList.sentinelKey) {

            HeapNode x = w;
            
            // because the x.right gets changed in link(), nab the next 
            // reference before link
            HeapNode next = w.getRight();
                                  
            int d = x.getNumberOfChildren();
            
            assert(d <= maxDegree);

            // is there another node of the same degree, that is, has the 
            // same number of children?
            while ((d < a.length) && (a[d] != null)) {
                
                HeapNode y = a[d];

                if (x.getKey() > y.getKey()) {
                    HeapNode tmp = x;
                    x = y;
                    y = tmp;                   
                }

                // link removes y (which has a key larger than x now) from 
                // rootList and adds it as a child of x
                link(y, x);
                
                a[d] = null;
                d++;
            }
            if (d < a.length) {
                a[d] = x;
            } else {
                throw new IllegalStateException("maxDegree=" + maxDegree 
                + " but d is " + d);
            }

            w = next;
        }

        minimumNode = null;
        
        // remove all from root list:
        rootList.resetSentinel();
        rootList.number = 0;
                
        for (int i = 0; i < a.length; i++) {
            if (a[i] != null) {
            
            	rootList.insert(a[i]);
                            	
                if ((minimumNode == null) || (a[i].getKey() < minimumNode.getKey()) ) {
                    minimumNode = a[i];
                }
            }
        }
    }

    void link(HeapNode y, HeapNode x) {
    	 // moves y to a child position of x
        rootList.remove(y);
        x.addChild(y);
        y.setParent(x);
        y.setMark(false);
    }

    /**
     * decrease key for node x
     *
     * runtime is O(1)
     *
     * @param x
     * @param decreaseToThisKey
     */
    public void decreaseKey(HeapNode x, long decreaseToThisKey) {
        if (decreaseToThisKey > x.getKey()) {
            throw new IllegalArgumentException(
                "key cannot be larger than x.key");
        }
        x.setKey(decreaseToThisKey);
        HeapNode y = x.getParent();
        if ((y != null) && (x.getKey() < y.getKey())) {
            cut(x, y);
            cascadingCut(y);
        }
        if (x.getKey() < minimumNode.getKey()) {
            minimumNode = x;
        }
    }

    /**
     * removes child node from tree and starts a new one with it.
     *
     * @param x
     * @param y 
     */
    protected void cut(HeapNode x, HeapNode y) {
        // remove x from child list of y and decrement y.degree
        y.removeChild(x);

        // add x to root list
        rootList.insert(x);
        x.setParent(null);
        x.setMark(false);
    }

    /**
     * c*O(1)
     *
     * @param y
     */
    protected void cascadingCut(HeapNode y) {
        HeapNode z = y.getParent();
        if (z != null) {
            if (!y.isMark()) {
                y.setMark(true);
            } else {
                cut(y, z);
                cascadingCut(z);
            }
        }
    }

    protected void remove(HeapNode x) {
        // runtime O(1)
        decreaseKey(x, DoubleLinkedCircularList.minValue);
        
        extractMin();
    }

    /**
     * a depth first search of all nodes (that is, it descends all children 
     * of a node before proceeding to the next in the current doubly linked 
     * circular list).
     * 
     * @param key
     * @return
     */
    HeapNode search(long key) {
        
        HeapNode node = rootList.getSentinel().getRight();
        
        int sentinel = DoubleLinkedCircularList.sentinelKey;
         
        Stack<HeapNode> stack = new Stack<HeapNode>();
        
        while (!stack.isEmpty() || (node.getKey() != sentinel)) {
            if (node.getKey() != sentinel) {

                stack.push(node);

                node = node.getRight();

            } else {
                
                node = stack.pop();

                //System.out.println(node.key);
                if (node.getKey() == key) {
                    return node;
                }

                node = node.getChildren().getSentinel().getRight();
            }
        }
        
        return null;
    }

    public DoubleLinkedCircularList getRootList() {
        return rootList ;
    }
    
    protected void printRootList() {
        StringBuffer sb = new StringBuffer();
        HeapNode t = this.rootList.getSentinel().getRight();
        while (t.getKey() != DoubleLinkedCircularList.sentinelKey) {
            if ((minimumNode == null) || (t.getKey() < minimumNode.getKey())) {
                minimumNode = t;
            }
            String str = String.format("%d", t.getKey());
            if (sb.length() > 0) {
                sb.append(" ");
            }
            sb.append(str);
            t = t.getRight();
        }
        sb.insert(0, "Looking for " + minimumNode + " :");
        log.info(sb.toString());
    }
}
