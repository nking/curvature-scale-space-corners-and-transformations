package algorithms.graphs;

import algorithms.QuickSort;
import algorithms.util.SimpleLinkedListNode;
import java.util.ArrayList;
import java.util.Arrays;
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
public class DFSIterative {
    /**
     * adjacency matrix with connected i->j indicated by the index and each
     *    node in the linked list, respectively.
     * for example, adjacent to node 3 is found via directedEdges[3] as all in the linked list.
     */
    private SimpleLinkedListNode[] directedEdges;
    
    /** 
     * holds state for whether a node has been visited.  0 = not visited,
     * 1 = visiting now, 2 = was visited.
    */
    private int[] visited;

    /**
     * time when node is first discovered
     */
    private int[] td;

    /**
     * time when node's adjacency list has all been visited
     */
    private int[] tf;
   
    private int[] predecessor;

    private int time;

    public DFSIterative() {
        
    }

    /**
     * @param directedEdges  adjacency matrix with connected i->j indicated 
     * by the index and each node in the linked list, respectively.
     * Note that the key of each node is expected to be the same as it's index
     * in the adjacency matrix.
     * For example, adjacent to node 3 is found via directedEdges[3] as all in 
     * the linked list.
     */
    public void walk(SimpleLinkedListNode[] directedEdges) {
        this.directedEdges = Arrays.copyOf(directedEdges, directedEdges.length);
        visited = new int[directedEdges.length];
        td = new int[directedEdges.length];
        tf = new int[directedEdges.length];
        predecessor = new int[directedEdges.length];
        Arrays.fill(td, -1);
        Arrays.fill(tf, -1);
        Arrays.fill(predecessor, -1);
        time = 0;
        
        for (int u = 0; u < directedEdges.length; u++) {
            if (visited[u] == 0) {
                walk(u);
            }
        }
    }
    
    private void walk(int u) {
        
        Stack<Snapshot> stack = new Stack<Snapshot>();
        Snapshot current;
        
        System.out.println("*load method frame for " + u);
        
        current = new Snapshot(u);
        current.stage = 0;
        stack.push(current);
        
        while(!stack.empty()) {
            
            current = stack.pop();
            
            System.out.println(current.toString());
            
            switch(current.stage) {
                case 0: { 
                    // before recursion is invoked
                    visited[current.node] = 1;
                    time++;
                    System.out.println("  0: visiting " + current.node + " to set td=" + time);
                    td[current.node] = time;
                    
                    current.stage = 1;
                    stack.push(current);
                    
                    System.out.format("  0: push onto stack u=%d\n", current.node);
                            
                    SimpleLinkedListNode next = directedEdges[current.node];
                    
                    if (next != null && next.getKey() != -1) {
                        
                        int v = next.getKey();
                        
                        directedEdges[current.node].delete(next);
                                                      
                        if (visited[v] == 0) {
                            
                            predecessor[v] = current.node;
                            
                            Snapshot newSnapshot = new Snapshot(v);
                            newSnapshot.stage = 0;
                            stack.push(newSnapshot);

                            System.out.format("   0: and push onto stack v=%d\n", v);
                            System.out.println("   0: [v: " + newSnapshot.toString() + "]");
  
                            continue;
                        }
                    }
                    break;
                }
                case 1: {
                    System.out.println(" 1: have all child links been visited?  snap="
                       + current.toString());
                    
                    SimpleLinkedListNode next = directedEdges[current.node];
                    if (next != null && next.getKey() != -1) {
                        
                        int v = next.getKey();
                        
                        System.out.format(" 1: there is a child link %d\n", v);
                        
                        directedEdges[current.node].delete(next);
                        
                        current.stage = 1;
                        stack.push(current);

                        System.out.format("  0: push onto stack u=%d\n", current.node);
                                                      
                        if (visited[v] == 0) {
                            
                            predecessor[v] = current.node;
                            
                            Snapshot newSnapshot = new Snapshot(v);
                            newSnapshot.stage = 0;
                            stack.push(newSnapshot);

                            System.out.format("   1: and push onto stack v=%d\n", v);
                            System.out.println("   1: [v: " + newSnapshot.toString() + "]");
  
                            continue;
                        }
                        
                        continue;
                    }
                    
                    visited[current.node] = 2;
                    time++;
                    tf [current.node] = time;
                    System.out.format(" 1: end visit to %d, set tf=%d\n",
                            current.node, time);

                    break;
                }
            }
        }
        
        //---------------------------------------
        
    }
    
    private class Snapshot {
        
        /**
         * index of current snapshot within this instance's arrays.
         */
        protected final int node;
                
        protected int stage = 0;
                        
        public Snapshot(int u) {
            this.node = u;
        }
                
        public Snapshot(Snapshot s) {
            this.stage = s.stage;
            this.node = s.node;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("node=").append(Integer.toString(node))
                .append(", stage=").append(Integer.toString(stage))
                .append(", prev=").append(Integer.toString(predecessor[node]))
                .append(", visited=").append(Integer.toString(visited[node]))
            ;
            return sb.toString();
        }
        
    }
    
    /**
     * return the indexes in order of the starts of their traversals
     * @return 
     */
    public int[] getOrderedBeginIndexes() {
        if (td == null) {
            return null;
        }
        return sortForIndexes(td);
    }
    
    private int[] sortForIndexes(int[] a) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (directedEdges == null) {
            return null;
        }
        assert(a.length == directedEdges.length);
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
        if (tf == null) {
            return null;
        }
        return sortForIndexes(tf);
    }
    
    public int[] getTd() {
        if (td == null) {
            return null;
        }
        return td;
    }

    public int[] getTf() {
        if (tf == null) {
            return null;
        }
        return tf;
    }
    
}
