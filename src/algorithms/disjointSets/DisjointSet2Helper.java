package algorithms.disjointSets;

import java.util.ArrayList;
import java.util.List;

/**
 * a disjoint set implemented as a forest.
 * 
 * based upon pseudocode from "Introduction to Algorithms" by Cormen et al.
 * 
 * @author nichole
 */
public class DisjointSet2Helper {

     /**
     * make a set out of the given node.
     * runtime complexity is O(1).
     * @param x
     * @return
     */
    public <T> DisjointSet2Node<T> makeSet(DisjointSet2Node<T> x) {
        x.setParent(x);
        x.setRank(0);
        return x;
    }
    
    /**
     * <pre>
     * find the set representative for the given node.  As a side effect,
     * also updates x and all of it's ancestors with the found result so that
     * they directly point to the top most parent and subsequent lookups are
     * faster (this is path compression).
     * runtime complexity:
     *     the method uses iteration.  
     *     if we represent x_height as the number of nodes between x's tree 
     *     root and x, we have 2 times x_height iterations of statements,
     * so O(x_height).   With path compression here, the amoritzed
     * worst-case running time is Θ(α(n)) where α is the inverse Ackermann 
     * function.  
     * 
     * The inverse Ackermann function is α(n) = min{k : A_k(1) ≥ n}.
     * For most purposes, α(n) = O(1) so then the amoritized running time is O(1).
     * </pre>
     * 
     * @param x
     * @return
     */
    public <T> DisjointSet2Node<T> findSet(DisjointSet2Node<T> x) {
                
        // iterative
        if (!x.equals(x.getParent())) {
            
            List<DisjointSet2Node<T>> update = new ArrayList<DisjointSet2Node<T>>();
            
            DisjointSet2Node<T> parent = x;
            while (!parent.equals(parent.getParent())) {
                update.add(parent);
                parent = parent.getParent();
            }
            
            // update the nodes with parent
            for (DisjointSet2Node<T> node : update) {
                node.setParent(parent);
            }
  
        }
        
        return x.getParent();
    }
    
    /**
      assigns to x and y the same parent from both of their parents, choosing
      the one with largest rank.
          
       Runtime complexity is O(1).
       
     * @param x
     * @param y
     * @return the root found to be the one with equal number of nodes or more nodes
     */
    private <T> DisjointSet2Node<T> link(DisjointSet2Node<T> x, DisjointSet2Node<T> y) {
        
        if (x.equals(y)) {
            return x;
        }
        
        DisjointSet2Node<T> parent;
        if (x.getRank() >= y.getRank()) {
            parent = x;
            y.setParent(parent);
            if (x.getRank() == y.getRank()) {
                parent.setRank(parent.getRank() + 1);
            }
        } else {
            parent = y;
            x.setParent(parent);
        }
        return parent;
    }
    
    /**
     * append the shorter list onto the end of the longer list.
     * The amortized runtime complexity is O(1) due to the path compression
     * in findSet.
     * The method is also known as "union-find" because it uses findSet 
     * as the first steps before link.
     * 
     * @param x
     * @param y
     * @return
     */
    public <T> DisjointSet2Node<T> union(DisjointSet2Node<T> x, DisjointSet2Node<T> y) {
        
        x = findSet(x);
        
        y = findSet(y);
        
        DisjointSet2Node<T> parent = link(x, y);

        return parent;
    }
    
    public static <T> String print(DisjointSet2Node<T> x) {
                
        StringBuilder sb = new StringBuilder();
        
        sb.append(x.toString());
        
        DisjointSet2Node<T> p = x.getParent();
        
        while (p != null) {
            DisjointSet2Node<T> nextP = p.getParent();
            if (nextP.equals(p)) {
                break;
            }
            sb.append("parent->").append(p.toString());
            p = nextP;
        }
        
        return sb.toString();
    }
  
}
