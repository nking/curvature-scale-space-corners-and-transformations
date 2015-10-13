package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import algorithms.util.Stack;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.CharBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import java.nio.file.Files;
import static java.nio.file.StandardOpenOption.APPEND;
import static java.nio.file.StandardOpenOption.WRITE;
import java.nio.charset.CharsetEncoder;
import java.nio.charset.CoderResult;
import java.nio.file.OpenOption;
import java.util.EnumSet;
import java.util.logging.Level;
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

	private long n = 0;

    /**
     * insert node into heap.  runtime is O(1).  makes no attempt to consolidate
     * tree.
     *
     * @param node
     */
    public void insert(HeapNode node) {
    	if (node.getKey() == DoubleLinkedCircularList.noValue) {
            throw new IllegalArgumentException(
                "node.key must be set before insert into heap." +
                " must have value != DoubleLinkedCircularList.noValue");
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

        long sentinel = DoubleLinkedCircularList.sentinelKey;

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
 if (n < 0) {
     int zz = 1;
 }
        return z;
    }

    public long getNumberOfNodes() {
        return n;
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

        // search rootList using in-order traversal

        HeapNode node = rootList.getSentinel().getRight();

        long sentinel = DoubleLinkedCircularList.sentinelKey;

        Stack<HeapNode> stack = new Stack<HeapNode>();

        while (!stack.isEmpty() || (node.getKey() != sentinel)) {
            if (node.getKey() != sentinel) {

                stack.push(node);

                node = node.getRight();

            } else {

                node = stack.pop();

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

    public void printRootList() {
        StringBuilder sb = new StringBuilder(
            String.format("(n=%d rootList.n=%d) ", n, rootList.number));
        sb.append(" minimumNode=");
        if (minimumNode != null) {
            sb.append(minimumNode);
        }
        sb.append(";  rootList=");
        HeapNode t = this.rootList.getSentinel().getRight();
        while (t.getKey() != DoubleLinkedCircularList.sentinelKey) {
            String str = String.format("%d", t.getKey());
            if (sb.length() > 0) {
                sb.append(" ");
            }
            sb.append(str);
            t = t.getRight();
        }

        log.info(sb.toString());
    }

    private BufferedWriter debugWriter = null;
    
    private BufferedWriter createWriter() throws IOException {
                
        String bin = ResourceFinder.findDirectory("bin");
        String fileName = "debug_heap_" + System.currentTimeMillis() + ".txt";
        String filePath = bin + "/" + fileName;
        File file = new File(filePath);
        
        BufferedWriter writer = Files.newBufferedWriter(file.toPath(), 
            Charset.forName("US-ASCII"));
        
        return writer;
    }

    @Override
    protected void finalize() throws Throwable {
        try {
            closeDebug();
        } finally {
            super.finalize();
        }
    }
    
    private void closeDebug() {
        if (debugWriter != null) {
            try {
                debugWriter.close();
            } catch (IOException ex) {
                log.severe(ex.getMessage());
            }
        }
    }
    
    public void printHeapToTestOut(String label) {
        
        boolean print = true;
        
        if (debugWriter == null) {
            try {
                debugWriter = createWriter();
            } catch (IOException ex) {
                closeDebug();
                print = false;
            }
        }
        
        if (print) {
            try {
                printHeap(label, debugWriter);
            } catch (IOException ex) {
                log.severe(ex.getMessage());
            }
        }
    }
    
    public void printHeap(String label, BufferedWriter writer) throws IOException {
        
        int bufferSize = 1024;//2 * 72 * 4;
                        
        if (label != null) {
            char[] c = label.toCharArray();
            writer.write(c, 0, c.length);
            writer.write("\n");
            writer.flush();
        }

        // traverse heap using in-order traversal
        
        char[] c = String.format("(n=%d rootList.n=%d) ", n, rootList.number).toCharArray();
        writer.write(c, 0, c.length);
        c = " minimumNode=".toCharArray();
        writer.write(c, 0, c.length);
        if (minimumNode != null) {
            c = minimumNode.toString().toCharArray();
            writer.write(c, 0, c.length);
        }
        c = ";  heap=\n".toCharArray();
        writer.write(c, 0, c.length);
        writer.flush();
            
        // pre-order traversal of the heap

        HeapNode node = rootList.getSentinel().getRight();

        long sentinel = DoubleLinkedCircularList.sentinelKey;

        Stack<HeapNode> stack = new Stack<HeapNode>();

        int currentLevel = -1;
        
        StringBuilder sb = new StringBuilder(bufferSize);

        int nIter = 0;
        
        while (!stack.isEmpty() || (node.getKey() != sentinel)) {
            
            nIter++;
            
            if (node.getKey() != sentinel) {

                currentLevel++;
                
                if (sb.length() > 72) {
                    sb.append("\n");
                    c = sb.toString().toCharArray();
                    writer.write(c, 0, c.length);
                    if (nIter % 100 == 0) {
                        writer.flush();
                    }
                    sb = new StringBuilder(bufferSize);
                    if (currentLevel > 0) {
                        sb.append("    ");
                    }
                }

                sb.append(" ").append("[").append(currentLevel).append("] key=");
                if (node.getKey() == Long.MAX_VALUE) {
                    sb.append("M");
                } else {
                    sb.append(node.getKey());
                }

                stack.push(node);

                node = node.getChildren().getSentinel().getRight();

            } else {

                node = stack.pop();

                boolean eol = (currentLevel == 0);
                if (!eol) {
                    eol = true;
                    int nSb = sb.length();
                    if (nSb > 1) {
                        int c0 = sb.charAt(nSb - 1);
                        int c1 = sb.charAt(nSb - 2);
                        int space = (int)' ';
                        if (c0 == space && c1 == space) {
                            eol = false;
                        }
                    }
                }
                if (!eol) {
                    if (sb.length() > 72) {
                        sb.append("\n");
                        c = sb.toString().toCharArray();
                        writer.write(c, 0, c.length);
                        if (nIter % 100 == 0) {
                            writer.flush();
                        }
                        sb = new StringBuilder("    ");
                    }
                }
                if (eol) {
                    sb.append("\n");
                    c = sb.toString().toCharArray();
                    writer.write(c, 0, c.length);
                    if ((nIter % 100) == 0) {
                        writer.flush();
                    }
                    sb = new StringBuilder();
                    if (currentLevel > 0) {
                        sb.append("    ");
                    }
                }

                currentLevel--;

                node = node.getRight();
            }
        }

        c = sb.toString().toCharArray();
        writer.write(c, 0, c.length);
        writer.flush();
    }

}
