package algorithms.util;

import junit.framework.TestCase;

public class DisjointSet2Test extends TestCase {

    public void test0() throws Exception {
        
        DisjointSet2Helper disjointSetHelper = new DisjointSet2Helper();
        
        DisjointSet2Node<String> x = new DisjointSet2Node<String>(new String("c"));
        
        DisjointSet2Node<String> xTree = disjointSetHelper.makeSet(x);
        
        assertTrue(x.getParent().equals(x));
        assertTrue(x.getRank() == 0);
        assertTrue(xTree.getParent().equals(x));
        assertTrue(xTree.getRank() == 0);
        
       
        DisjointSet2Node<String> x2 = new DisjointSet2Node<String>();
        x2.setMember(new String("h"));
        
        DisjointSet2Node<String> x2Tree = disjointSetHelper.makeSet(x2);
        
        assertTrue(x2.getParent().equals(x2));
        assertTrue(x2.getRank() == 0);
        assertTrue(x2Tree.getParent().equals(x2));
        assertTrue(x2Tree.getRank() == 0);
        
        xTree = disjointSetHelper.union(xTree, x2Tree);
        // when the ranks are equal, the 1st becomes parent
        assertTrue(x.getParent().equals(x));
        assertTrue(x2.getParent().equals(x));
        assertTrue(xTree.getParent().equals(x));
        assertTrue(xTree.getRank() == 1);

        
        DisjointSet2Node<String> x3 = new DisjointSet2Node<String>();
        x3.setMember(new String("b"));
        
        DisjointSet2Node<String> x3Tree = disjointSetHelper.makeSet(x3);
        
        assertTrue(x3.getParent().equals(x3));
        assertTrue(x3.getRank() == 0);
        assertTrue(x3Tree.getParent().equals(x3));
        assertTrue(x3Tree.getRank() == 0);
        
        
        xTree = disjointSetHelper.union(xTree, x3Tree);
        assertTrue(x.getParent().equals(x));
        assertTrue(x2.getParent().equals(x));
        assertTrue(x3.getParent().equals(x));
        assertTrue(xTree.getParent().equals(x));
        // the rank doesn't increase unless they have equal ranks
        assertTrue(xTree.getRank() == 1);
        
        
        DisjointSet2Node<String> x4 = new DisjointSet2Node<String>();
        x4.setMember(new String("e"));
        
        DisjointSet2Node<String> x4Tree = disjointSetHelper.makeSet(x4);
        
         xTree = disjointSetHelper.union(xTree, x4Tree);
        assertTrue(x.getParent().equals(x));
        assertTrue(x2.getParent().equals(x));
        assertTrue(x3.getParent().equals(x));
        assertTrue(x4.getParent().equals(x));
        assertTrue(xTree.getParent().equals(x));
        assertTrue(xTree.getRank() == 1);
    }
}