package algorithms.graphs;

import algorithms.QuickSort;
import algorithms.disjointSets.DisjointSet2Helper;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
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

*   NOTE: have added use of disjoint sets to determine the independent sets
     during the traversal.
     Disjoint sets use path compression to update representatives quickly,
     hence, all connected components can point to the same parent.
     * 
 * @author nichole
 */
public class DFSWithIndependentSets {

    /**
     * adjacency matrix with connected i->j indicated by the index and each
     *    node in the linked list, respectively.
     * for example, adjacent to node 3 is found via directedEdges[3] as all in the linked list.
     */
    protected SimpleLinkedListNode[] g;

    /** 
     * holds state for whether a node has been visited.  0 = not visited,
     * 1 = visiting now, 2 = was visited.
    */
    protected int[] visited;

    /**
     * time when node is first discovered
     */
    protected int[] td;

    /**
     * time when node's adjacency list has all been visited
     */
    protected int[] tf;
   
    protected int[] predecessor;

    protected int time = 0;
    
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
            
    private Logger log = Logger.getLogger(getClass().getSimpleName());
    
    private Level logLevel = Level.FINE;
    
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
        if (directedEdges == null) {
            throw new IllegalArgumentException("directedEdges cannot be null");
        }
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
                visit(u);
            } else {
                addToMap(u, predecessor[u]);
            }
        }
        
        populateGroupMap();
        
        modifySequenceList();
    }
    
    private void visit(int u) {
        log.log(logLevel, "  visit u=" + u);
        
        visited[u] = 1;
        time++;
        //System.out.println("  visiting " + u + " to set td=" + time);
        td[u] = time;
        
        addToMap(u, predecessor[u]);

        SimpleLinkedListNode next = g[u];
        
        while (next != null && next.getKey() != -1) {
            int v = next.getKey();
            log.log(logLevel, "        v=" + v);
            if (visited[v] == 0) {
                predecessor[v] = u;
                visit(v);
            } else {
                addToMap(v, u);
            }
            next = next.getNext();
        }
        //addToMap(u, predecessor[u]);
        visited[u] = 2;
        time++;
        tf[u] = time;
                
        //System.out.println("  visited " + u + ") to set tf=" + time);
    }
    
    private void modifySequenceList() {
        assert(tf != null);
        int p, t;
        
        for (int i = 1; i < tf.length; ++i) {
            p = predecessor[i];
            t = tf[i];
            //if ()
        }
    }
    
    /**
     * update indexParentMap and independentGroupsMap for the given node and
     * presence or absence of predecessor.
     * @param node 
     */
    private void addToMap(int nodeIdx, int prevIdx) {
        
        log.log(logLevel, "  addToMap: nodeIdx=" + nodeIdx + " prevIdx=" +
            prevIdx + " predecessor[" + nodeIdx + "]=" +
            predecessor[nodeIdx]);
        
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
        
        NOTE: these objects and thier names will change soon.
        when finished, will have data structure for sequential indexes in tf,
        and will have data structures holding independent sets.
        
        key = graph index. value = holder for DJSet which has parent idx as key for parentGroupMap.
        TIntObjectHashMap<DisjointSetHolder> indexDJSetMap
     
        TIntObjectHashMap<TIntHashSet> parentGroupMap
        
        key = top-most predecessor index. value = set of connected nodes.  
              each set is independent of one another
              (no connecting edges between them).
        TIntObjectHashMap<TIntHashSet> independentGroupsMap;
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
                
                if (prevRef.equals(nodeRef)) {
                    // Case (?):  
                    return;
                }
                if (prevDJSet.equals(nodeDJSet) || 
                    prevDJSet.getMember().intValue() == nodeDJSet.getMember().intValue()) {
                    // Case (?):  
                    return;
                }
                
                log.log(logLevel, "   merge: nodeIdx=" + nodeIdx + " prevIdx=" +
                    prevIdx + "  predecessor[" + nodeIdx + "]=" +
                    predecessor[nodeIdx] + 
                    "\n          nodeDJSet=" + nodeDJSet.toString() +
                    "\n          prevDJSet=" + prevDJSet.toString()                        
                );
                
                prevDJSet = disjointSetHelper.unionChooseY(nodeDJSet, prevDJSet);
                
                prevRef.set = prevDJSet;
                nodeRef.set = prevDJSet;
                log.log(logLevel, "    merged: " +
                    "\n          prevDJSet=" + prevDJSet.toString()                        
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
                
                log.log(logLevel, "    merge: nodeIdx=" + nodeIdx + " prevIdx=" +
                    prevIdx + " predecessor[" + nodeIdx + "]=" + predecessor[nodeIdx] + 
                    "\n           temp=" + temp.toString() +
                    "\n           prevDJSet=" + prevDJSet.toString()                        
                );
                
                prevDJSet = disjointSetHelper.unionChooseY(temp, prevDJSet);
                
                prevRef.set = prevDJSet;
                
                log.log(logLevel, "   merged: " +
                    "\n          prevDJSet=" + prevDJSet.toString()                        
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
            "\n         prevIdx=" + prevIdx + " prevRef=" + indexDJSetMap.get(prevIdx) +
            "\n         nodeIdx=" + nodeIdx + " nodeRef=" + indexDJSetMap.get(nodeIdx));
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
