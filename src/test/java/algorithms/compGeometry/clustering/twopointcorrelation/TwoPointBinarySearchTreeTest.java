package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.clustering.twopointcorrelation.TwoPointBinarySearchTree.Node;
import java.util.ArrayList;
import java.util.List;
import java.util.Stack;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TwoPointBinarySearchTreeTest extends TestCase {

    public TwoPointBinarySearchTreeTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    /**
     * Test of storeIfDoesNotContain method, of class TwoPointIdentities.
     */
    public void testInsert() {
        int count = 0;
        ITwoPointIdentity identities = new TwoPointBinarySearchTree();
        for (int i = 1; i < 121; i++) {
            for (int j = (i + 1); j < 121; j++) {
                ((TwoPointBinarySearchTree)identities).insert(i, j);
                count++;
            }
        }
        boolean stored = false;
        for (int i = 1; i < 121; i++) {
            for (int j = (i + 1); j < 121; j++) {
                Node node = ((TwoPointBinarySearchTree)identities).search(i, j);
                assertTrue(node.a0 == i && node.a1 == j);
                stored = identities.storeIfDoesNotContain(i, j);
                assertFalse(stored);
            }
        }
        
        assertTrue(((TwoPointBinarySearchTree)identities).n == count);
        
        List<Node> nodes = preOrderTraversal(
            ((TwoPointBinarySearchTree)identities).root);
        
        assertTrue(nodes.size() == count);
    }

    /**
     * Test of storeIfDoesNotContain method, of class TwoPointIdentities.
     */
    public void testStoreIfDoesNotContain() {

        int pointIndex0 = 0;
        int pointIndex1 = 10;

        ITwoPointIdentity identities = new TwoPointBinarySearchTree();
        boolean stored = identities.storeIfDoesNotContain(pointIndex0, pointIndex1);
        assertTrue(stored);

        stored = identities.storeIfDoesNotContain(pointIndex0, pointIndex1);
        assertFalse(stored);

        for (int i = 1; i < 121; i++) {
            for (int j = (i + 1); j < 121; j++) {
                stored = identities.storeIfDoesNotContain(i, j);
                assertTrue(stored);
            }
        }

        for (int i = 1; i < 121; i++) {
            for (int j = (i + 1); j < 121; j++) {
                stored = identities.storeIfDoesNotContain(i, j);
                assertFalse(stored);
            }
        }
    }

    protected List<Node> preOrderTraversal(Node node) {
            
        List<Node> list = new ArrayList<Node>();

        //root, left subtree, right subtree
        Stack<Node> stack = new Stack<Node>();
        while (!stack.isEmpty() || (node != null)) {
            if (node != null) {
                list.add(node);
                stack.push(node);
                node = node.getLeft();
            } else {
                node = stack.pop();
                node = node.getRight();
            }
        }
        return list;
    }
}
