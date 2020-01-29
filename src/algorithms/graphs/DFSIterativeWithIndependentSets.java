package algorithms.graphs;

import algorithms.QuickSort;
import algorithms.disjointSets.DisjointSet2Helper;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;

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
   
  NOTE: have added use of disjoint sets to determine the independent sets
     during the traversal.
     Disjoint sets use path compression to update representatives quickly,
     hence, all connected components can point to the same parent.
         
 * @author nichole
 */
public class DFSIterativeWithIndependentSets {
    /**
     * adjacency matrix with connected i->j indicated by the index and each
     *    node in the linked list, respectively.
     * for example, adjacent to node 3 is found via directedEdges[3] as all in the linked list.
     */
    private SimpleLinkedListNode[] g;
    
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
     * value = disjoint set holder, parent.member of internal set is key for parentGroupMap
     */
    private TIntObjectHashMap<DisjointSetHolder> indexDJSetMap;
    
    /**
     * key = top-most predecessor index
     * value = sets of paths through the code present as sub-sequences in tf
     */
    private TIntObjectHashMap<TIntHashSet> parentGroupMap;
    
    private DisjointSet2Helper disjointSetHelper;
    
    private int time;
    
    private Logger log = Logger.getLogger(getClass().getSimpleName());
    
    private Level logLevel = Level.FINE;

    public DFSIterativeWithIndependentSets() {
        
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
        if (directedEdges == null) {
            throw new IllegalArgumentException("directedEdges cannot be null");
        }
        g = directedEdges.clone();
        for (int i = 0; i < g.length; ++i) {
            g[i] = new SimpleLinkedListNode(directedEdges[i]);
        }
        visited = new int[g.length];
        td = new int[g.length];
        tf = new int[g.length];
        predecessor = new int[g.length];
        Arrays.fill(td, -1);
        Arrays.fill(tf, -1);
        Arrays.fill(predecessor, -1);
        time = 0;
        indexDJSetMap = new TIntObjectHashMap<DisjointSetHolder>();
        parentGroupMap = new TIntObjectHashMap<TIntHashSet>();
        
        disjointSetHelper = new DisjointSet2Helper();
        
        for (int u = 0; u < g.length; u++) {
            if (visited[u] == 0) {
                walk(u);
            }
        }
        
        populateGroupMap();
    }
    
    private void walk(int u) {
        
        Stack<Snapshot> stack = new Stack<Snapshot>();
        Snapshot current;
        
        log.log(logLevel, "*load method frame for " + u);
        
        current = new Snapshot(u);
        current.stage = 0;
        stack.push(current);
        
        while(!stack.empty()) {
            
            current = stack.pop();
            
            log.log(logLevel, current.toString());
            
            switch(current.stage) {
                case 0: { 
                    // before recursion is invoked
                    visited[current.node] = 1;
                    time++;
                    log.log(logLevel, "  stage 0: visiting " + current.node + " to set td=" + time);
                    td[current.node] = time;
                    
                    current.stage = 1;
                    stack.push(current);
                    
                    addToMap(current.node, predecessor[current.node]);
                    
                    log.log(logLevel, 
                        String.format("  stage 0: push onto stack u=%d\n", current.node));
                            
                    SimpleLinkedListNode next = g[current.node];
                    
                    if (next != null && next.getKey() != -1) {
                        
                        int v = next.getKey();
                        
                        g[current.node].delete(next);
                                                      
                        if (visited[v] == 0) {
                            
                            predecessor[v] = current.node;
                            
                            addToMap(v, predecessor[v]);
                            
                            Snapshot newSnapshot = new Snapshot(v);
                            newSnapshot.stage = 0;
                            stack.push(newSnapshot);

                            log.log(logLevel, 
                                String.format("   stage 0: and push onto stack v=%d\n", v));
                            log.log(logLevel, "   stage 0: [v: " + newSnapshot.toString() + "]");
  
                            continue;
                        } else if (predecessor[v] == -1) {
                            // in case the instance graph is not ordered top-down
                            predecessor[v] = current.node;
                            
                            addToMap(v, predecessor[v]);
                        } else {
                            addToMap(v, current.node);
                        }
                    }
                    break;
                }
                case 1: {
                    log.log(logLevel, " stage 1: have all child links been visited?  snap="
                       + current.toString());
                    
                    SimpleLinkedListNode next = g[current.node];
                    if (next != null && next.getKey() != -1) {
                        
                        int v = next.getKey();
                        
                        log.log(logLevel, 
                            String.format(" stage 1: there is a child link %d\n", v));
                        
                        g[current.node].delete(next);
                        
                        current.stage = 1;
                        stack.push(current);

                        log.log(logLevel, 
                            String.format("  stage 1: push onto stack u=%d\n", current.node));
                                                      
                        if (visited[v] == 0) {
                            
                            predecessor[v] = current.node;
                            
                            addToMap(v, predecessor[v]);
                            
                            Snapshot newSnapshot = new Snapshot(v);
                            newSnapshot.stage = 0;
                            stack.push(newSnapshot);

                            log.log(logLevel, String.format(
                                "   stage 1: and push onto stack v=%d\n", v));
                            log.log(logLevel, "   stage 1: [v: " + newSnapshot.toString() + "]");
  
                            continue;
                        } else if (predecessor[v] == -1) {
                            // in case the instance graph is not ordered top-down
                            predecessor[v] = current.node;
                            
                            addToMap(v, current.node);
                        } else {
                            addToMap(v, current.node);
                        }
                        
                        continue;
                    } else {
                        addToMap(current.node, predecessor[current.node]);
                    }
                    
                    visited[current.node] = 2;
                    time++;
                    tf[current.node] = time;
                    log.log(logLevel, String.format(" stage 1: end visit to %d, set tf=%d\n",
                        current.node, time));

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
    private void addToMap(int nodeIdx, int prevIdx) {
        
        /*
        5 cases:
            prevIdx exists:
                (0) nodeIdx exists and is the same as prevIdx
                    (in this case, no updates are necessary)
                (1) nodeIdx exists.
                    (in this case, nodeIdx data was stored upon visit from path
                     of another parent node.  This additional parent node, prevIdx,
                     has existing data, so the 2 datasets need to be merged)
                (2) nodeIdx does not exist.
                    (in this case, nodeIdx data is created, added to prevIdx
                     data, and uses the existing djset to store the merged values in
                     index map.)
        
            prevIdx does not exist:
               (in this case, nodeIdx data was stored upon visit from path
                  of a parent node and this somehow got invoked during
                  the loop of the edge indexes.  no need to update anything)
               (3) nodeIdx exists.
                   (this should not happen)
               (4) nodeIdx does not exist
                   (in this case, new data is created for nodeIdx and stored in maps)
        
        key = graph index. value = holder for DJSet which has parent idx as key for parentGroupMap.
        TIntObjectHashMap<DisjointSetHolder> indexDJSetMap
     
        TIntObjectHashMap<TIntHashSet> parentGroupMap
        */
        
        if (prevIdx > -1) {
            // merge existing prev with (1) existing or (2) new node
            DisjointSetHolder prevRef = indexDJSetMap.get(prevIdx);
            assert(prevRef != null);
            DisjointSet2Node<Integer> prevDJSet = prevRef.set;
            if (prevDJSet == null) {
                throw new IllegalStateException("prevDJSet should not be null");
            }
            DisjointSet2Node<Integer> prevParent = prevDJSet.getParent();
            int prevParentIdx = prevParent.getMember();
                       
            if (indexDJSetMap.contains(nodeIdx)) {
                if (nodeIdx == prevIdx) {
                    //Case (0)
                    return;
                }
                // Case (1) merge sets
                DisjointSetHolder nodeRef = indexDJSetMap.get(nodeIdx);
                DisjointSet2Node<Integer> nodeDJSet = nodeRef.set;
                
                if (prevDJSet.equals(nodeDJSet)) {
                    // Case (?):  
                    return;
                }
                
                log.log(logLevel, "  merge: nodeIdx=" + nodeIdx + " prevIdx=" +
                    prevIdx + " predecessor[" + nodeIdx + "]=" +
                    predecessor[nodeIdx] + 
                    "\n   nodeDJSet=" + nodeDJSet.toString() +
                    "\n    prevDJSet=" + prevDJSet.toString()                        
                );
                
                prevDJSet = disjointSetHelper.unionChooseY(nodeDJSet, prevDJSet);
                
                prevRef.set = prevDJSet;
                nodeRef.set = prevDJSet;
                log.log(logLevel, "  merged: " +
                    "\n    prevDJSet=" + prevDJSet.toString()                        
                );
                
                // indexParentMap entries already existed for both nodes,
                //   so only needed to update existing values,
                //   no need for indexParentMap.put
                
                //Note that the parent set node will be pruned after all map additions
                //in O(N) by traversing the index map to keep only existing parent
                //nodes in the parentGroupMap.
                
            } else {
                // Case (2) create new data and add it to prev
                indexDJSetMap.put(nodeIdx, prevRef);
                                
                DisjointSet2Node<Integer> temp = new DisjointSet2Node<Integer>(nodeIdx);
                temp = disjointSetHelper.makeSet(temp);
                
                log.log(logLevel, "  merge: nodeIdx=" + nodeIdx + " prevIdx=" +
                    prevIdx + "    \npredecessor[" + nodeIdx + "]=" +
                    predecessor[nodeIdx] + 
                    "\n    temp=" + temp.toString() +
                    "]n    prevDJSet=" + prevDJSet.toString()                        
                );
                
                prevDJSet = disjointSetHelper.unionChooseY(temp, prevDJSet);
                
                prevRef.set = prevDJSet;
                
                log.log(logLevel, "  merged: " +
                    "\n    prevDJSet=" + prevDJSet.toString()                        
                );
                
                indexDJSetMap.put(nodeIdx, prevRef);
            }
            
        } else {
            // prev does not exist so this must be a top-most node
            if (indexDJSetMap.contains(nodeIdx)) {
                // Case (3) prevIdx does not exist but nodeIdx does
                return;
            } else {
                // Case (4) create data for nodeIdx and store it for itself and as it's own parent
                DisjointSetHolder nodeRef = new DisjointSetHolder();
                nodeRef.set = new DisjointSet2Node<Integer>(nodeIdx);
                nodeRef.set = disjointSetHelper.makeSet(nodeRef.set);
                assert(nodeRef.set.getParent().getMember().intValue() == nodeIdx);
                indexDJSetMap.put(nodeIdx, nodeRef);

                TIntHashSet nodeGroupSet = new TIntHashSet();
                nodeGroupSet.add(nodeIdx);
            }
        }
        log.log(logLevel, "   addToMap results:" + 
            "\n    prevIdx=" + prevIdx + " prevRef=" + indexDJSetMap.get(prevIdx) +
            "\n    nodeIdx=" + nodeIdx + " nodeRef=" + indexDJSetMap.get(nodeIdx));
    }

    private void populateGroupMap() {
        
        TIntObjectIterator<DisjointSetHolder> iter = indexDJSetMap.iterator();
        for (int ii = indexDJSetMap.size(); ii-- > 0;) {
            iter.advance();
            int idx = iter.key();
            int pIdx = iter.value().set.getParent().getMember();
            TIntHashSet set = parentGroupMap.get(pIdx);
            if (set == null) {
                set = new TIntHashSet();
                parentGroupMap.put(pIdx, set);
            }
            set.add(idx);
        }
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
     * class to hold a disjoint set node so that the reference stored as a value
     * in a java HashMap is updated and shared among more than one instance,
     * rather than obsolete after a union creates a new parent node for other
     * elements.
     */
    private class DisjointSetHolder {
        protected DisjointSet2Node<Integer> set = null;

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("hash: @").append(Integer.toHexString(hashCode()));
            if (set != null) {
                sb.append(" ").append(set.toString());
            }
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
        if (g == null) {
            return null;
        }
        assert(a.length == g.length);
        a = Arrays.copyOf(a, a.length);
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
    public TIntObjectHashMap<TIntHashSet> getIndependentSets() {
        return parentGroupMap;
    }
    public int getParentIndexForIndependentSets(int node) {
        if (!indexDJSetMap.contains(node)) {
            throw new IllegalStateException("there was no entry for node=" + node);
        }
        DisjointSet2Node<Integer> parentNode = indexDJSetMap.get(node).set;
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
