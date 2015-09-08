package algorithms.disjointSets;

import java.util.HashSet;
import java.util.Set;

import junit.framework.TestCase;

public class DisjointSetTest extends TestCase {

    public void test0() throws Exception {
        
        DisjointSetHelper disjointSetHelper = new DisjointSetHelper();
        
        DisjointSetNode<String> x = new DisjointSetNode<String>();
        x.setMember(new String("c"));
        
        DisjointSet<String> xList = disjointSetHelper.makeSet(x);
        
        assertTrue(x.getRepresentative().equals(x));
        assertTrue(xList.getHead().equals(x));
        assertTrue(xList.getTail().equals(x));
        assertTrue(xList.getNumberOfNodes() == 1);
        
        
        
        DisjointSetNode<String> x2 = new DisjointSetNode<String>();
        x2.setMember(new String("h"));
        
        DisjointSet<String> x2List = disjointSetHelper.makeSet(x2);
        
        assertTrue(x2.getRepresentative().equals(x2));
        assertTrue(x2List.getHead().equals(x2));
        assertTrue(x2List.getTail().equals(x2));
        
        
        xList = disjointSetHelper.union(xList, x2List);
        assertTrue(x.getRepresentative().equals(x));
        assertTrue(x2.getRepresentative().equals(x));
        assertTrue(xList.getHead().equals(x));
        assertTrue(xList.getTail().equals(x2));
        assertTrue(xList.getNumberOfNodes() == 2);
        
        
        DisjointSetNode<String> x3 = new DisjointSetNode<String>();
        x3.setMember(new String("b"));
        
        DisjointSet<String> x3List = disjointSetHelper.makeSet(x3);
        
        assertTrue(x3.getRepresentative().equals(x3));
        assertTrue(x3List.getHead().equals(x3));
        assertTrue(x3List.getTail().equals(x3));
        
        xList = disjointSetHelper.union(xList, x3List);
        assertTrue(x.getRepresentative().equals(x));
        assertTrue(x2.getRepresentative().equals(x));
        assertTrue(x3.getRepresentative().equals(x));
        assertTrue(xList.getHead().equals(x));
        assertTrue(xList.getTail().equals(x3));
        assertTrue(xList.getNumberOfNodes() == 3);
        
        
        DisjointSetNode<String> x4 = new DisjointSetNode<String>();
        x4.setMember(new String("e"));
        
        DisjointSet<String> x4List = disjointSetHelper.makeSet(x4);
        
        assertTrue(x4.getRepresentative().equals(x4));
        assertTrue(x4List.getHead().equals(x4));
        assertTrue(x4List.getTail().equals(x4));
        
        xList = disjointSetHelper.union(xList, x4List);
        assertTrue(x.getRepresentative().equals(x));
        assertTrue(x2.getRepresentative().equals(x));
        assertTrue(x3.getRepresentative().equals(x));
        assertTrue(x4.getRepresentative().equals(x));
        assertTrue(xList.getHead().equals(x));
        assertTrue(xList.getTail().equals(x4));
        assertTrue(xList.getNumberOfNodes() == 4);
        
        
        Set<String> set = new HashSet<>();
        set.add(x.getMember());
        set.add(x2.getMember());
        set.add(x3.getMember());
        set.add(x4.getMember());
        
        DisjointSetNode<String> current = xList.getHead();
        while (current != null) {
            if (set.contains(current.getMember())) {
                set.remove(current.getMember());
            }
            current = current.getNext();
        }
        assertTrue(set.isEmpty());
    }

}
