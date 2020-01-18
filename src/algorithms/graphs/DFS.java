package algorithms.graphs;

import algorithms.MultiArrayMergeSort;
import algorithms.QuickSort;
import algorithms.util.SimpleLinkedListNode;
import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

/**
   DFS

   searches the full depth of a graph or subgraph when possible first then
      backtracks to the unexplored edges and unexplored nodes repeating until
      all nodes are visited.  unlike BFS, it may contain many predecessor trees, 
      that is a predecessor forest of nodes that are the shortest from the 
      source to each reachable node.  for this reason, DFS searches can need a 
      lot of memory.

   average runtime is approx O(|E|), worst case runtime: O(|V| + |E|)
   worst case space needed: O(|V|)

   implemented following Cormen et al. "Introduction To Algorithms"

 * @author nichole
 */
public class DFS {

    /**
     * adjacency matrix with connected i->j indicated by the index and each
     *    node in the linked list, respectively.
     * for example, adjacent to node 3 is found via directedEdges[3] as all in the linked list.
     */
    protected final SimpleLinkedListNode[] directedEdges;

    /** 
     * holds state for whether a node has been visited.  0 = not visited,
     * 1 = visiting now, 2 = was visited.
    */
    protected final int[] visited;

    /**
     * time when node is first discovered
     */
    protected final int[] td;

    /**
     * time when node's adjacency list has all been visited
     */
    protected final int[] tf;
   
    protected final int[] predecessor;

    protected int time;
    
    /**
     * @param theDirectedEdges  adjacency matrix with connected i->j indicated 
     * by the index and each node in the linked list, respectively.
     * Note that the key of each node is expected to be the same as it's index
     * in the adjacency matrix.
     * For example, adjacent to node 3 is found via directedEdges[3] as all in 
     * the linked list.
     */
    public DFS(SimpleLinkedListNode[] theDirectedEdges) {
        this.directedEdges = theDirectedEdges;
        this.visited = new int[directedEdges.length];
        this.td = new int[directedEdges.length];
        this.tf = new int[directedEdges.length];
        this.predecessor = new int[directedEdges.length];
    }

    public void walk() {

        init();

        for (int u = 0; u < directedEdges.length; u++) {
            if (visited[u] == 0) {
                visit(u);
            }
        }        
    }
    
    private void init() {
        for (int u = 0; u < directedEdges.length; u++) {
            visited[u] = 0;
            predecessor[u] = -1;
        }
        time = 0;
    }
    
    private void visit(int u) {
        visited[u] = 1;
        time++;
        //System.out.println("visit(" + u + ") to set td=" + time);
        td[u] = time;

        SimpleLinkedListNode adjacent = directedEdges[u];
        
        while (adjacent != null && adjacent.getKey() != -1) {
            int v = adjacent.getKey();
            if (visited[v] == 0) {
                predecessor[v] = u;
                visit(v);                 
            }
            adjacent = adjacent.getNext();
        }
        visited[u] = 2;
        time++;
        tf[u] = time;
        //System.out.println("   visited(" + u + ") to set tf=" + time);
    }
   
      /**
     * return the indexes in order of the starts of their traversals
     * @return 
     */
    public int[] getOrderedBeginIndexes() {
        return sortForIndexes(td);
    }
    
    private int[] sortForIndexes(int[] a) {
        assert(a.length == this.directedEdges.length);
        int[] idxs = new int[a.length];
        for (int i = 0; i < idxs.length; ++i) {
            idxs[i] = i;
        }
        QuickSort.sortBy1stArg(a, idxs);
        return idxs;
    }
    
    /**
     * return the indexes in order of the ends of their traversal
     * @return 
     */
    public int[] getOrderedEndIndexes() {
        return sortForIndexes(tf);
    }
    
    public int[] getTd() {
        return td;
    }

    public int[] getTf() {
        return tf;
    }
    
}
