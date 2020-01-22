package algorithms.graphs;

import algorithms.QuickSort;
import algorithms.disjointSets.DisjointSet2Helper;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
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
   
  NOTE: have added use of maps and sets to determine the independent groups
     during the traversal.
     Also using a disjoint set to keep track of the topmost parent in a the
     independent connected set.
     Disjoint sets use path compression to update representatives quickly.
     (THIS IS NOT FINISHED YET).
    
           3
       1       2
    1: mkset, map:(1, set1)
    2: mkset, map:(2, set2)
    3: mkset, map:(3, set3)
       3's links: merge 1 into 3, rm set1, map: 1, set3, rmset1 
   
 * @author nichole
 */
public class DFSWithIndependentSets {
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
    
    /**
     * key = graph index.
     * value = top-most predecessor index.
     */
    private TIntObjectHashMap<DisjointSet2Node<Integer>> indexParentMap;
    
    /**
     * key = top-most predecessor index
     * value = sets of paths through the code present as sub-sequences in tf
     */
    private TIntObjectHashMap<TIntHashSet> parentGroupMap;
    
    /**
     * key = top-most predecessor index
     * value = set of connected nodes.  each set is independent of one another
     *    (no connecting edges between them).
     */
    private TIntObjectHashMap<TIntHashSet> independentGroupsMap;
    
    private DisjointSet2Helper disjointSetHelper;
    
    private int time;

    public DFSWithIndependentSets() {
        
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
        indexParentMap = new TIntObjectHashMap<DisjointSet2Node<Integer>>();
        parentGroupMap = new TIntObjectHashMap<TIntHashSet>();
        independentGroupsMap = new TIntObjectHashMap<TIntHashSet>();
        
        disjointSetHelper = new DisjointSet2Helper();
        
        for (int u = 0; u < directedEdges.length; u++) {
            if (visited[u] == 0) {
                walk(u);
            }
        }
    }
    
    private void walk(int u) {
        
        Stack<Snapshot> stack = new Stack<Snapshot>();
        Snapshot current;
        
        //System.out.println("*load method frame for " + u);
        
        current = new Snapshot(u);
        current.stage = 0;
        stack.push(current);
        
        while(!stack.empty()) {
            
            current = stack.pop();
            
            //System.out.println(current.toString());
            
            switch(current.stage) {
                case 0: { 
                    // before recursion is invoked
                    visited[current.node] = 1;
                    time++;
                    //System.out.println("  0: visiting " + current.node + " to set td=" + time);
                    td[current.node] = time;
                    
                    current.stage = 1;
                    stack.push(current);
                    
                    addToMap(current.node);
                    
                    //System.out.format("  0: push onto stack u=%d\n", current.node);
                            
                    SimpleLinkedListNode next = directedEdges[current.node];
                    
                    if (next != null && next.getKey() != -1) {
                        
                        int v = next.getKey();
                        
                        directedEdges[current.node].delete(next);
                                                      
                        if (visited[v] == 0) {
                            
                            predecessor[v] = current.node;
                            
                            addToMap(v);
                            
                            Snapshot newSnapshot = new Snapshot(v);
                            newSnapshot.stage = 0;
                            stack.push(newSnapshot);

                            //System.out.format("   0: and push onto stack v=%d\n", v);
                            //System.out.println("   0: [v: " + newSnapshot.toString() + "]");
  
                            continue;
                        } else if (predecessor[v] == -1) {
                            // in case the instance graph is not ordered top-down
                            predecessor[v] = current.node;
                            
                            addToMap(v);
                        }
                        //NOTE: could make another change here to add one more
                        //   else if conditional for the case of
                        //     visited[v] > 0 and predecessor[v] > -1
                        //     and have knowledge that predecessor path for
                        //     predecessor[v] is shorter than predecessor[u]
                        //     so would want to reset predecessor[v] to u in that case.
                        //     Currently making another version of DFS which adds disjoint sets,
                        //     so this might be addressable there.
                    }
                    break;
                }
                case 1: {
                    //System.out.println(" 1: have all child links been visited?  snap="
                    //   + current.toString());
                    
                    SimpleLinkedListNode next = directedEdges[current.node];
                    if (next != null && next.getKey() != -1) {
                        
                        int v = next.getKey();
                        
                        //System.out.format(" 1: there is a child link %d\n", v);
                        
                        directedEdges[current.node].delete(next);
                        
                        current.stage = 1;
                        stack.push(current);

                        //System.out.format("  0: push onto stack u=%d\n", current.node);
                                                      
                        if (visited[v] == 0) {
                            
                            predecessor[v] = current.node;
                            
                            addToMap(v);
                            
                            Snapshot newSnapshot = new Snapshot(v);
                            newSnapshot.stage = 0;
                            stack.push(newSnapshot);

                            //System.out.format("   1: and push onto stack v=%d\n", v);
                            //System.out.println("   1: [v: " + newSnapshot.toString() + "]");
  
                            continue;
                        } else if (predecessor[v] == -1) {
                            // in case the instance graph is not ordered top-down
                            predecessor[v] = current.node;
                            
                            addToMap(v);
                        }
                        
                        continue;
                    }
                    
                    visited[current.node] = 2;
                    time++;
                    tf[current.node] = time;
                    //System.out.format(" 1: end visit to %d, set tf=%d\n",
                    //    current.node, time);

                    break;
                }
            }
        }
    }

    /**
     * update indexParentMap and independentGroupsMap for the given node and
     * presence or absence of predecessor.
     * @param node 
     */
    private void addToMap(int node) {
        
        /*
                3
            1       2
        1: mkset, map:(1, set1)
        2: mkset, map:(2, set2)
        3: mkset, map:(3, set3),
            3's links: merge 1 into 3, rm set1, map: 1, set3, rmset1 
        
        Uodate indexParentMap and independentGroupsMap.
        
        key = graph index. value = top-most predecessor index.
        TIntObjectHashMap<DisjointSet2Node<Integer>> indexParentMap
     
        key = top-most predecessor index. value = set of connected nodes.  
              each set is independent of one another
              (no connecting edges between them).
        TIntObjectHashMap<TIntHashSet> independentGroupsMap;
        */
        
        DisjointSet2Node<Integer> connectedRepr = indexParentMap.get(node);
        if (connectedRepr == null) {
            connectedRepr = new DisjointSet2Node<Integer>(node);
            connectedRepr = disjointSetHelper.makeSet(connectedRepr);
        }
        
        TIntHashSet indepSet = parentGroupMap.get(node);
        if (indepSet == null) {
            indepSet = new TIntHashSet();
            indepSet.add(node);
        } else {
            parentGroupMap.remove(node);
        }
        int parentNode = node;
        if (predecessor[node] > -1) {
            DisjointSet2Node<Integer> prevRepr = indexParentMap.get(predecessor[node]);
            assert(prevRepr != null);
            DisjointSet2Node<Integer> parent = prevRepr.getParent();
            parentNode = parent.getMember();
            
            TIntHashSet indepSet2 = parentGroupMap.get(parentNode);
            assert(indepSet2 != null);
            indepSet.addAll(indepSet2);
            //independentGroupsMap.remove(parentNode);
            
            connectedRepr = disjointSetHelper.union(prevRepr, connectedRepr);
            connectedRepr.setParent(parent);
            
            indexParentMap.put(parentNode, connectedRepr);
        } 
        
        parentGroupMap.put(parentNode, indepSet); 
        indexParentMap.put(node, connectedRepr);
    }
    
    private class Snapshot {
        
        /**
         * index of current snapshot within DFSIterative instance's arrays.
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
        return td;
    }

    public int[] getTf() {
        return tf;
    }
    
    /**
     * get a map of sequential path members in tf
     * @return a map w/ keys = parent index of independent set, value = 
     *    indexes of members in a sequential path in tf
     */
    public TIntObjectHashMap<TIntHashSet> getSequentialTFGroups() {
        return parentGroupMap;
    }
    public int getParentIndexForSequentialTFGroups(int node) {
        DisjointSet2Node<Integer> parentNode = indexParentMap.get(node);
        if (parentNode == null) {
            throw new IllegalStateException("there was no entry for node=" + node);
        }
        return parentNode.getParent().getMember();
    }

    /*    
    public TIntObjectHashMap<TIntHashSet> getIndependentGroups() {
    }
    public int getParentIndexForIndependentGroup(int node) {
    }
    */
    
    /**
     * get predecessor indexes
     * @return get predecessor indexes
     */
    public int[] getPredecessorIndexes() {
        if (predecessor == null) {
            return null;
        }
        return Arrays.copyOf(predecessor, predecessor.length);
    }
    
    /**
     * creates a string of the independent sets in format:
     *    parent=(top-most index), set=(indexes in independent set)
     * @return 
     */
    public String printIndependentSets() {
        
        StringBuilder sb = new StringBuilder();
        if (parentGroupMap == null) {
            return sb.toString();
        }
        
        TIntObjectIterator<TIntHashSet> iter = this.parentGroupMap.iterator();
        
        for (int ii = this.parentGroupMap.size(); ii-- > 0;) {
            iter.advance();
            int parentNode = iter.key();
            TIntHashSet indepSet = iter.value();
            sb.append("parent=").append(Integer.toString(parentNode)).append("; set=");
            TIntIterator iter2 = indepSet.iterator();
            while (iter2.hasNext()) {
                int idx = iter2.next();
                sb.append(Integer.toString(idx)).append(", ");
            }
            sb.append("\n");
        }
        return sb.toString();
    }
}
