package algorithms.imageProcessing;

import java.security.SecureRandom;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class HeapTest extends TestCase {
    
    Logger log = Logger.getLogger(this.getName());

    public HeapTest(String testName) {
        super(testName);
    }

    public void testInsert() {

        Heap h = new Heap();

        h.insert(new HeapNode(3));
        assertTrue(h.getRootList().getSentinel().getRight().getKey() == 3);
        assertTrue(h.getRootList().getSentinel().getRight().getRight().getKey() == DoubleLinkedCircularList.sentinelKey);
        assertTrue(h.getRootList().getSentinel().getRight().getLeft().getKey() == DoubleLinkedCircularList.sentinelKey);
        assertTrue(h.minimumNode.getKey() == 3);

        h.insert(new HeapNode(1));
        assertTrue(h.getRootList().getSentinel().getRight().getKey() == 1);
        assertTrue(h.getRootList().getSentinel().getRight().getRight().getKey() == 3);
        assertTrue(h.getRootList().getSentinel().getRight().getRight().getRight().getKey() == DoubleLinkedCircularList.sentinelKey);
        assertTrue(h.getRootList().getSentinel().getRight().getLeft().getKey() == DoubleLinkedCircularList.sentinelKey);
        assertTrue(h.minimumNode.getKey() == 1);

        h.insert(new HeapNode(2));
        assertTrue(h.getRootList().getSentinel().getRight().getKey() == 2);
        assertTrue(h.getRootList().getSentinel().getRight().getRight().getKey() == 1);
        assertTrue(h.getRootList().getSentinel().getRight().getRight().getRight().getKey() == 3);
        assertTrue(h.getRootList().getSentinel().getRight().getRight().getRight().getRight().getKey() == DoubleLinkedCircularList.sentinelKey);
        assertTrue(h.getRootList().getSentinel().getRight().getLeft().getKey() == DoubleLinkedCircularList.sentinelKey);
        assertTrue(h.minimumNode.getKey() == 1);
    }

    public void testInsert2() {
        
        // example from cormen et al. chap 20, Fibonacci heaps, Fig 20.2

    	HeapNode node23 = new HeapNode(23);

    	HeapNode node7 = new HeapNode(7);

    	HeapNode node3 = new HeapNode(3);


    	HeapNode node18 = new HeapNode(18);
    	HeapNode node39 = new HeapNode(39);
    	node18.addChild(node39);

    	HeapNode node38 = new HeapNode(38);
    	HeapNode node41 = new HeapNode(41);
    	node38.addChild(node41);

    	HeapNode node52 = new HeapNode(52);

    	node3.addChild(node38);
    	node3.addChild(node52);
    	node3.addChild(node18);


    	HeapNode node17 = new HeapNode(17);
    	HeapNode node30 = new HeapNode(30);
    	node17.addChild(node30);

    	HeapNode node24 = new HeapNode(24);

    	HeapNode node26 = new HeapNode(26);
    	HeapNode node35 = new HeapNode(35);
    	node26.addChild(node35);

    	HeapNode node46 = new HeapNode(46);

    	node24.addChild(node46);
    	node24.addChild(node26);

    	node18.setMark(true);
    	node39.setMark(true);
    	node26.setMark(true);

        Heap h = new Heap();
        h.insert(node24);
        h.insert(node17);
        h.insert(node3);
        h.insert(node23);
        h.insert(node7);

        assertTrue(h.minimumNode.getKey() == node3.getKey());
        assertTrue(h.getRootList().number == 5);
        assertTrue(h.getRootList().getSentinel().getRight().getKey() == 7);
        assertTrue(h.getRootList().getSentinel().getRight().getRight().getKey() == 23);
        assertTrue(h.getRootList().getSentinel().getRight().getRight().getRight().getKey() == 3);
        assertTrue(h.getRootList().getSentinel().getRight().getRight().getRight().getRight().getKey() == 17);
        assertTrue(h.getRootList().getSentinel().getRight().getRight().getRight().getRight().getRight().getKey() == 24);

        HeapNode node21 = new HeapNode(21);
        h.insert(node21);

        assertTrue(h.minimumNode.getKey() == node3.getKey());
        assertTrue(h.getRootList().number == 6);

    }

    public void testConsolidate() {
        Heap h = new Heap();
        h.insert(new HeapNode(3));
        h.insert(new HeapNode(1));
        h.insert(new HeapNode(2));

        //  [2] -- [1] -- [3]

        assertTrue(h.getRootList().getSentinel().getRight().getKey() == 2);
        assertTrue(h.getRootList().getSentinel().getRight().getRight().getKey() == 1);
        assertTrue(h.getRootList().getSentinel().getRight().getRight().getRight().getKey() == 3);

        h.consolidate();

        // [3] -- [1]
        //         |
        //        [2]
/*log.info("h.getRootList().getSentinel().right.key=" + h.getRootList().getSentinel().right.key);
log.info("h.getRootList().getSentinel().right.right.key=" + h.getRootList().getSentinel().right.right.key);
log.info("h.getRootList().getSentinel().right.getChildren().getSentinel().right.key=" + h.getRootList().getSentinel().right.getChildren().getSentinel().right.key);
log.info("h.getRootList().getSentinel().right.right.right.key=" + h.getRootList().getSentinel().right.right.right.key);
log.info("h.getRootList().getSentinel().right.left.key=" + h.getRootList().getSentinel().right.left.key);
log.info("h.minimumNode.key=" + h.minimumNode.key);*/

        assertTrue(h.getRootList().getSentinel().getRight().getKey() == 1);
        assertTrue(h.getRootList().getSentinel().getRight().getRight().getKey() == 3);
        // link places '2' into a child node of '1'
        assertTrue(h.getRootList().getSentinel().getRight().getChildren().getSentinel().getRight().getKey() == 2);
        assertTrue(h.getRootList().getSentinel().getRight().getRight().getRight().getKey() == DoubleLinkedCircularList.sentinelKey);
        assertTrue(h.getRootList().getSentinel().getRight().getLeft().getKey() == DoubleLinkedCircularList.sentinelKey);
        assertTrue(h.minimumNode.getKey() == 1);
    }

    public void testExtractMin() {
        
        // extractmin brings the children of deleted node up to the rootList
        
        Heap h = new Heap();
        h.insert(new HeapNode(3));
        h.insert(new HeapNode(1));
        h.insert(new HeapNode(2));
        
        h.printHeapToTestOut("testExtractMin");

        HeapNode min1 = h.extractMin();
        assertNotNull(min1);
        assertTrue(min1.getKey() == 1);
        //    [2]     after first extraction, we have:
        //     |
        //    [3]

        h.printHeapToTestOut("--");
        
        HeapNode min2 = h.extractMin();
        assertNotNull(min2);
        assertTrue(min2.getKey() == 2);
        
        h.printHeapToTestOut("--");
        
        assertTrue(h.minimumNode.getKey() == 3);
        
        h.printRootList();
    }
    
    public void testExtractMin2() {
        
        // example from cormen et al. chap 20, Fibonacci heaps, Fig 20.3

        HeapNode node23 = new HeapNode(23);

        HeapNode node7 = new HeapNode(7);

        HeapNode node3 = new HeapNode(3);


        HeapNode node18 = new HeapNode(18);
        HeapNode node39 = new HeapNode(39);
        //node18.addChild(node39);

        HeapNode node38 = new HeapNode(38);
        HeapNode node41 = new HeapNode(41);
        //node38.addChild(node41);

        HeapNode node52 = new HeapNode(52);

        //node3.addChild(node38);
        //node3.addChild(node52);
        //node3.addChild(node18);


        HeapNode node17 = new HeapNode(17);
        HeapNode node30 = new HeapNode(30);
        //node17.addChild(node30);

        HeapNode node24 = new HeapNode(24);

        HeapNode node26 = new HeapNode(26);
        HeapNode node35 = new HeapNode(35);
        //node26.addChild(node35);

        HeapNode node46 = new HeapNode(46);

        //node24.addChild(node46);
        //node24.addChild(node26);
        
        HeapNode node21 = new HeapNode(21);

        node18.setMark(true);
        node39.setMark(true);
        node26.setMark(true);

        Heap h = new Heap();
        h.insert(node3);
        h.insert(node17);
        h.insert(node24);
        h.insert(node23);
        h.insert(node7);
        h.insert(node21);
        
        assertTrue(h.getNumberOfNodes() == 6);
        
        
        // the insert removes the children, and unmarks the nodes so redo that for test conditions
        node3.addChild(node38);
        node3.addChild(node52);
        node3.addChild(node18);
        node17.addChild(node30);
        node24.addChild(node46);
        node24.addChild(node26);
        node18.addChild(node39);
        node38.addChild(node41);       
        node26.addChild(node35);
                
        node18.setMark(true);
        node39.setMark(true);
        node26.setMark(true);
        
        //heap.n will be wrong

        assertTrue(h.minimumNode.getKey() == node3.getKey());
        assertTrue(h.getRootList().number == 6);
        
        h.printHeapToTestOut("--");
        
        HeapNode min = h.extractMin();
        
        h.printHeapToTestOut("--");
        
        // results are different from the figure 20.3 because I use a DoubleLinkedCircularList which adds nodes to the right of the sentinel.
        
        assertTrue(min.getKey() == 3);
        
        assertTrue(h.minimumNode.getKey() == 7);
        assertTrue(h.getRootList().number == 3);
        
        assertTrue(h.minimumNode.getNumberOfChildren() == 3);
        
        assertTrue(h.minimumNode.getRight().getNumberOfChildren() == 2);
        
        assertTrue(h.minimumNode.getRight().getRight().getNumberOfChildren() == 1);
        
//         * *key=7 (children=3)
//    [junit]     key=18 (children=2)
//    [junit]         key=38 (children=1)
//    [junit]             key=41 (children=0)
//    [junit]         key=39 (children=0)
//    [junit]     key=21 (children=1)
//    [junit]         key=52 (children=0)
//    [junit]     key=23 (children=0)
//    [junit] *key=24 (children=2)
//    [junit]     key=26 (children=1)
//    [junit]         key=35 (children=0)
//    [junit]     key=46 (children=0)
//    [junit] *key=17 (children=1)
//    [junit]     key=30 (children=0)
        
    }

    public void testExtractMin4() {
        Heap h = new Heap();
        h.insert(new HeapNode(0));
        h.insert(new HeapNode(Long.MAX_VALUE));
        h.insert(new HeapNode(Long.MAX_VALUE));
        h.insert(new HeapNode(Long.MAX_VALUE));
        h.insert(new HeapNode(Long.MAX_VALUE));

        HeapNode min1 = h.extractMin();
        assertNotNull(min1);
        assertTrue(min1.getKey() == 0);

        HeapNode min2 = h.extractMin();
        assertNotNull(min2);
        assertTrue(min2.getKey() == Long.MAX_VALUE);
    }
   
    public void testDecreaseKey0() {
        
        //example from Cormen et al. chap 20, Fibonacci heaps, Fig 20.4

        HeapNode node7 = new HeapNode(7);
        HeapNode node24 = new HeapNode(24);
        HeapNode node17 = new HeapNode(17);
        HeapNode node23 = new HeapNode(23);
        HeapNode node26 = new HeapNode(26);
        HeapNode node35 = new HeapNode(35);
        HeapNode node46 = new HeapNode(46);
        HeapNode node30 = new HeapNode(30);
        HeapNode node18 = new HeapNode(18);
        HeapNode node21 = new HeapNode(21);
        HeapNode node39 = new HeapNode(39);
        HeapNode node52 = new HeapNode(52);
        HeapNode node38 = new HeapNode(38);
        HeapNode node41 = new HeapNode(41);

        Heap h = new Heap();
        h.insert(node38);
        h.insert(node18);
        h.insert(node7);
                
        assertTrue(h.getNumberOfNodes() == 3);
        
        // the insert removes the children, and unmarks the nodes so redo that for test conditions
        node7.addChild(node23);
        node7.addChild(node17);
        node7.addChild(node24);
        node24.addChild(node46);
        node24.addChild(node26);
        node26.addChild(node35);
        node17.addChild(node30);
        node18.addChild(node39);       
        node18.addChild(node21);
        node21.addChild(node52);
        node38.addChild(node41);
        
        // heap.n will now be wrong
        
        node18.setMark(true);
        node39.setMark(true);
        node26.setMark(true);

        assertTrue(h.minimumNode.getKey() == 7);
        assertTrue(h.getRootList().number == 3);
                
        h.printHeapToTestOut("--");
        
        h.decreaseKey(node46, 15);
        
        h.printHeapToTestOut("--");
        
//         *  *key=15 (children=0)
//    [junit] *key=7 (children=3)
//    [junit]     key=24 (children=1)
//    [junit]         key=26 (children=1)
//    [junit]             key=35 (children=0)
//    [junit]     key=17 (children=1)
//    [junit]         key=30 (children=0)
//    [junit]     key=23 (children=0)
//    [junit] *key=18 (children=2)
//    [junit]     key=21 (children=1)
//    [junit]         key=52 (children=0)
//    [junit]     key=39 (children=0)
//    [junit] *key=38 (children=1)
//    [junit]     key=41 (children=0)

        
        assertTrue(h.minimumNode.getKey() == 7);
        assertTrue(h.getRootList().number == 4);
        
        HeapNode t = h.search(46);
        assertNull(t);
        t = h.search(15);
        assertNotNull(t);
        assertTrue(t.getKey() == 15);
        
        h.decreaseKey(node35, 5);
        
        h.printHeapToTestOut("--");
        
        //
 //        *  *key=24 (children=0)
 //   [junit] *key=26 (children=0)
 //   [junit] *key=5 (children=0)
 //   [junit] *key=15 (children=0)
 //   [junit] *key=7 (children=2)
 //   [junit]     key=17 (children=1)
 //   [junit]         key=30 (children=0)
 //   [junit]     key=23 (children=0)
 //   [junit] *key=18 (children=2)
 //   [junit]     key=21 (children=1)
 //   [junit]         key=52 (children=0)
 //   [junit]     key=39 (children=0)
 //   [junit] *key=38 (children=1)
 //   [junit]     key=41 (children=0)
        
        assertTrue(h.minimumNode.getKey() == 5);
        assertTrue(h.getRootList().number == 7);
        
        t = h.search(35);
        assertNull(t);
        t = h.search(5);
        assertNotNull(t);
        assertTrue(t.getKey() == 5);
    }

    public void testDecreaseKey() {
        Heap h = new Heap();
        h.insert(new HeapNode(3));
        HeapNode m = new HeapNode(10);
        h.insert(m);
        h.insert(new HeapNode(2));

        h.decreaseKey(m, 1);

        HeapNode min = h.extractMin();
        assertNotNull(min);
        assertTrue(min.getKey() == 1);

    }

    public void testDelete() {
        HeapNode t;
        
        Heap h = new Heap();
        h.insert(new HeapNode(3));
        HeapNode m = new HeapNode(10);
        h.insert(m);
        h.insert(new HeapNode(2));

        t = h.search(2);
        assertNotNull(t);
        
        h.remove(m);

        t = h.search(10);
        assertNull(t);
        t = h.search(3);
        assertNotNull(t);
        t = h.search(2);
        assertNotNull(t);
        
        HeapNode min = h.extractMin();
        assertTrue(min.getKey() == 2);
        
        assertTrue(h.minimumNode.getKey() == 3);
        
        t = h.search(3);
        assertNotNull(t);
        
        t = h.search(2);
        assertNull(t);
    }

    public void testRandomInsertRemove() throws Exception {

    	int niter = 100;

	    SecureRandom r = SecureRandom.getInstance("SHA1PRNG");
	    r.setSeed(System.nanoTime());

	    Heap heap = new Heap();

	    Set<HeapNode> searchForNodes = new HashSet<HeapNode>();

	    for (int i = 0; i < niter; i++) {

	    	int num = r.nextInt(Integer.MAX_VALUE);

	    	HeapNode node = new HeapNode();
	    	node.setKey(num);

	    	heap.insert(node);

            searchForNodes.add(node);
	    }

	    HeapNode node = heap.extractMin();
	    while (node != null) {
	    	searchForNodes.remove(node);
	    	node = heap.extractMin();
	    }

	    assertTrue(searchForNodes.isEmpty());
    }

    public void testRandomInsertRemove_0() throws Exception {

    	int niter = 100;

	    SecureRandom r = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.nanoTime();
        //seed=1393747922105894000l;
	    r.setSeed(seed);
        
        System.out.println("seed=" + seed);

	    Heap heap = new Heap();

	    int[] insertedKeys = new int[niter];
        
	    for (int i = 0; i < niter; i++) {

            int num = r.nextInt(Integer.MAX_VALUE);
	    	
	    	insertedKeys[i] = num;

	    	HeapNode node = new HeapNode();
	    	node.setKey(num);

	    	heap.insert(node);
	    }

        Arrays.sort(insertedKeys);
        
        for (int i = 0; i < insertedKeys.length; i++) {

            Integer expKey = insertedKeys[i];

            HeapNode node = heap.extractMin();
            
            assertTrue(node.getKey() == expKey);
        }

        // decreaseKey
        heap = new Heap();
        for (int i = 0; i < niter; i++) {

	    	int num = r.nextInt();

	    	HeapNode node = new HeapNode();
	    	node.setKey(num);

	    	heap.insert(node);

            int decreaseToThisKey = num - 1;

            heap.decreaseKey(node, decreaseToThisKey);

            assertTrue(node.getKey() == decreaseToThisKey);
	    }
    }

}
