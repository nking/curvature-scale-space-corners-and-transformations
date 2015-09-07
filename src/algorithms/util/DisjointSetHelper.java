package algorithms.util;

/**
 * a disjoint set implemented with linked lists.
 * each set is a linked list.
 * 
 * based upon pseudocode from "Introduction to Cormen et al." by Cormen et al.
 * 
 * @author nichole
 */
public class DisjointSetHelper {

    /**
     * make a set out of the given node.
     * runtime complexity is O(1).
     * 
     * @param x
     * @return
     */
    public DisjointSet makeSet(DisjointSetNode x) {
        x.representative = x;
        DisjointSet list = new DisjointSet();
        list.head = x;
        list.tail = x;
        list.numberOfNodes = 1;
        return list;
    }
    
    /**
     * find the set representative for the given node.
     * runtime complexity is O(1).
     * @param x
     * @return
     */
    public DisjointSetNode findSet(DisjointSetNode x) {
        return x.representative;
    }
    
    /**
     * append the shorter list onto the end of the longer's list.
     * runtime complexity is  O(N_shorter).
     * @param x
     * @return
     */
    public DisjointSet union(DisjointSet x, DisjointSet y) {
        
        DisjointSet longer, shorter;
       
        if (x.numberOfNodes >= y.numberOfNodes) {
            longer = x;
            shorter = y;
        } else {
            longer = y;
            shorter = x;
        }
        
        // add next references to longer
        longer.tail.next = shorter.head;
        
        DisjointSetNode latest = shorter.head;
        while (latest != null) {
            latest.representative = longer.head;
            latest = latest.next;
        }
        longer.tail = shorter.tail;
        
        longer.setNumberOfNodes(longer.getNumberOfNodes() + shorter.getNumberOfNodes());
        
        return longer;
    }
    
    public static String print(DisjointSet x) {
        
        DisjointSetNode current = x.getHead();
        
        StringBuilder sb = new StringBuilder();
        while (current != null) {
            if (sb.length() > 0) {
                sb.append(",");
            }
            sb.append(current.getMember().toString());
            current = current.getNext();
        }
        
        return sb.toString();
    }
  
}
