package algorithms.bipartite;

import algorithms.bipartite.MinCostUnbalancedAssignment.Forest;
import algorithms.bipartite.MinCostUnbalancedAssignment.LeftNode;
import algorithms.bipartite.MinCostUnbalancedAssignment.RightNode;
import algorithms.bipartite.MinCostUnbalancedAssignment.PathNode;
import algorithms.bipartite.MinCostUnbalancedAssignment.SinkNode;
import algorithms.bipartite.MinCostUnbalancedAssignment.SourceNode;
import static algorithms.bipartite.MinCostUnbalancedAssignment.extractNodes;
import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.Heap;
import algorithms.imageProcessing.HeapNode;
import algorithms.mst.PrimsMST;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * NOTE:  this class is NOT READY FOR USE YET
 * 
 * @author nichole
 */
public class HopcroftKarpRT2012 {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * NOT READY FOR USE.
     runtime complexity O(m * sqrt(s)) where m is number of edges
     and s is the size of the matching whose target size
     may be less than the maximum matchable
     * @param g
     * @return 
     */
    public Map<Integer, Integer> findMaxMatching(Graph g, int s) {
    
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
                 
        ResidualDigraph rM = new ResidualDigraph(g, m);
        
        int prevMSize = m.size();
        
        int nIter = 0;
            
        // estimate lambda for the dial array length by
        // looking at the number of edges with more than
        // one connection.
        int lambda = estimateLambda(rM);
                
        while (true) {
                        
            boolean augmented = buildForestAndAugment(rM, lambda);

            if (!augmented) {
                return m;
            }
            
            /*            
            let P = {P1, P2, ...Pk} be a maximum set of vertex-disjoint
               shortest augmenting paths with respect to M
            M = the symmetric difference between M and
               (P1 union P2 union ...Pk)
            
            now applying the symmetric difference to m and m2
            */
            
            Map<Integer, Integer> m2 = rM.extractMatchings();
                        
            //announce(M is a matching)
            
            log.info("nIter=" + nIter + " m2.size=" + m2.size()
                + " m.size=" + m.size());
            /*
            // debug:
            for (Entry<Integer, Integer> entry : m2.entrySet()) {
                log.info("m2 match= " + entry.getKey() + "->" + entry.getValue());
            }
            for (Entry<Integer, Integer> entry : m.entrySet()) {
                log.info("m match= " + entry.getKey() + "->" + entry.getValue());
            }
            */
            
            if (m2.size() >= s) {
                return m2;
            }
            
            Map<Integer, Integer> tmpM = new HashMap<Integer, Integer>(m);
                    
            Set<Integer> mR = new HashSet<Integer>(m.values());
            // if m2 has a match that does not conflict with m,
            // add it to m (== vertex disjoint)
            for (Map.Entry<Integer, Integer> entry : m2.entrySet()) {
                Integer key = entry.getKey();
                Integer value = entry.getValue();
                if (!tmpM.containsKey(key) && !mR.contains(value)) {
                    tmpM.put(key, value);
                }
            }
            
            log.fine("union size=" + tmpM.size());
            
            if (tmpM.size() >= s) {
                return tmpM;
            }
            
            // remove the intersection of m and m2
            for (Map.Entry<Integer, Integer> entry : m.entrySet()) {
                Integer key = entry.getKey();
                if (m2.containsKey(key) &&
                    m2.get(key).equals(entry.getValue())) {
                    tmpM.remove(key);
                }
            }
            
            m = tmpM;
            log.fine("symmetric diff of m and m2.size=" + m.size());
            /*for (Entry<Integer, Integer> entry : m.entrySet()) {
                log.info("sym diff m match= " + entry.getKey() + "->" + entry.getValue());
            }*/
            
            if (m.size() >= s) {
                return m;
            }
            
            if (prevMSize == m.size()) {
                // user estimate of matching size might be an
                // overestimate.
                log.warning("m.size=" + m.size() + " which is"
                    + " less than the requsted size s=" +s);
                return m;
            }
            
            rM = new ResidualDigraph(g, m);
            
            assert(prevMSize < m.size());
            
            prevMSize = m.size();
            ++nIter;            
        }        
    }
    
    private int estimateLambda(ResidualDigraph rM) {
        //TODO: this may need to be revised
        int n = 0;
        for (Map.Entry<Integer, Set<Integer>> entry :
            rM.getForwardLinksRM().entrySet()) {
            if (entry.getValue().size() > 1) {
                n++;
            }
        }
        return (n + 2);
    }

    /**
     * find augmenting paths of minimal length.
     * by building the shortest path forest.
     * 
     * implementing section 3.4, pg 22 of Ramshaw and Tarjan 2012.
     * 
     * @param rM the residual digraph built from the matching 
     * graph M of graph G.
     * 
     * @return a forest of trees whose forest array indexes 
     * are the path lengths and whose items hold at the root, the
     * remaining maidens, that is unmatched left nodes (a.k.a. G's
     * X nodes).
     * @param lambda the length to use for a countins sort of
     * augmenting path lengths
     */
    protected boolean buildForestAndAugment(
        final ResidualDigraph rM, int lambda) {
      
        Forest forest = new Forest(lambda);
        
        MinHeapForRT2012 heap = new MinHeapForRT2012(4);
        
        /*
        this class is invoked as the first step in a
        hopcroft-karp matching algorithm and there are
        at first no matching nodes,
        so all X nodes are maiden nodes.
        but other uses my give the method a residual graph
        that has matched arcs in it.
        
        The code below treats each maiden X node as a single
        source shortest path problem.
        
        The storage of the found paths is a little more complex
        in that an array indexed by calculated path lengths
        is created and each item is a doubly linked list of
        paths (this is array indexing pattern of the "Dial" 
        algorithm and it's similar to part of "counting sort").
        so forest[0] holds a doubly linked list.
           each item in that doubly linked list is a maiden 
             node without a predecessor and having dist=0.        
        
        Note that book-keeping for the predecessors requires
        a copy be made upon insert of a node into the heap
        so that all possible foward links are present in 
        the heap.
        */
     
        // init all nodes to inf length
        Map<Integer, LeftNode> leftNodes = new HashMap<Integer, LeftNode>();
        Map<Integer, RightNode> rightNodes = new HashMap<Integer, RightNode>();
        for (int i = 0; i < rM.getNRight(); ++i) {
            Integer rNode = Integer.valueOf(i);
            RightNode node = new RightNode();
            node.setKey(Long.MAX_VALUE);
            node.setData(rNode);
            rightNodes.put(rNode, node);
        }
        for (int i = 0; i < rM.getNLeft(); ++i) {
            Integer lNode = Integer.valueOf(i);
            LeftNode node = new LeftNode();
            node.setKey(Long.MAX_VALUE);
            node.setData(lNode);
            leftNodes.put(lNode, node);
        }
        
        /*
        for the bfs and tracking "visited", need to track
        individually for each maiden as its own single shortest
        path.  doing so for the right nodes
        key = maiden node index (== X index)
        value = visited nodes along key path
        */
        Map<Integer, Set<Integer>> vXY = new HashMap<Integer, Set<Integer>>();
        
        Set<Integer> augmentedLeft = new HashSet<Integer>();
        Set<Integer> augmentedRight = new HashSet<Integer>();
        long prevKey = -1;
        long lastAugKey = -1;
        
        // married X nodes
        Set<Integer> matchedLeft = new HashSet<Integer>(
            rM.getBackwardLinksRM().values());    
 
        // for all maidens
        // set key to 0, then ScanAndAdd(index)
        Set<Integer> maidens = new HashSet<Integer>();
        for (int i = 0; i < rM.getNLeft(); ++i) {
            Integer lNode = Integer.valueOf(i);
            if (matchedLeft.contains(lNode)) {
                continue;
            }
            maidens.add(lNode);
            LeftNode node = leftNodes.get(lNode);
            node.setKey(0);
            vXY.put(lNode, new HashSet<Integer>());
            prevKey = scanAndAdd(heap, forest, rM, 
                rightNodes, node, prevKey,
                augmentedLeft, augmentedRight, 
                node, vXY.get(lNode));
            assert(prevKey == 0L);
        }
     
        int nRight = rM.getNRight();
        
        log.fine("done adding " + maidens.size() + 
            " maiden nodes and their links to heap");
        
        // at this point, the maidens are all in index 0 of
        // the forest trees.
        
        int nIter = 0;
        
        while (heap.getNumberOfNodes() > 0) {
                        
            // in the heap are men not in the forest who are in
            // an alternating path from a maiden.
            // the key is the length of shortest path so far
            PathNode y = heap.extractMin();
            assert(y instanceof RightNode);
            assert(y.getData() != null);
            
            log.fine("heap.size=" + heap.getNumberOfNodes());
            
            log.fine("extractMin=" + y.toString());
            
            Integer yIndex = (Integer)(y.getData());
        
            if (augmentedRight.contains(yIndex)) {
                continue;
            }
            
            assert(y.topPredecessor != null);
            
            Integer topIndex = (Integer)y.topPredecessor.getData();
        
            long currentKey = forest.add(y, prevKey);
            
            if (currentKey > prevKey) {
                log.info("augment forest[" + prevKey + "] ");
                //forest[0] are the maiden nodes alone
                if (prevKey > 0) {
                    //debug(forest);
                    augmentPath(rM, forest, prevKey, augmentedLeft,
                        augmentedRight);
                }
                forest.removeFirstItem((int)prevKey);
                lastAugKey = prevKey;
                prevKey = currentKey;
            }
            
            /*
            if y is married then
                x := wife of y;
                set l(x) := l(y) and 
                ScanAndAdd(x); 
            else
                exit(bachelor β := y reached); 
            fi;
            */
            
            // the married nodes are the keys in the backward
            // links of the residual digraph
            Integer xIndex = rM.getBackwardLinksRM().get(yIndex);
            if (xIndex != null) {
                
                LeftNode xNode = leftNodes.get(xIndex);
                
                xNode.setKey(y.getKey());
                if (xNode.pathPredecessor == null) {
                    xNode.pathPredecessor = y;
                    xNode.topPredecessor = y.topPredecessor;
                }
                 
                currentKey = scanAndAdd(heap, forest, 
                    rM, rightNodes, xNode, prevKey,
                    augmentedLeft, augmentedRight, y.topPredecessor, 
                    vXY.get(topIndex));
                
                if (currentKey > prevKey) {
                    log.info("augment forest[" + prevKey + "] ");
                    if (prevKey > 0) {
                        //debug(forest);
                        augmentPath(rM, forest, prevKey, augmentedLeft,
                            augmentedRight);
                    }
                    forest.removeFirstItem((int)prevKey);
                    lastAugKey = prevKey;
                    prevKey = currentKey;    
                }       
                
            } else {
                //exit(bachelor β := y reached);
                log.fine("bachelor y=" + y.toString());
                //break;
                
                if (maidens.contains(topIndex)) {
                    maidens.remove(topIndex);                    
                } else if (maidens.isEmpty()) {
                    log.fine("last maiden's bachelor reached");
                    //debug(forest);
                    break;
                }
            }
        
            log.fine("nIter=" + nIter);
        }
        
        if (lastAugKey < prevKey) {
            log.info("augment forest[" + prevKey + "] ");
            if (prevKey > 0) {
                //debug(forest);
                augmentPath(rM, forest, prevKey, augmentedLeft,
                    augmentedRight);
            }
            lastAugKey = prevKey;
        }
       
        log.fine("lastAugKey=" + lastAugKey);
        
        return !(augmentedLeft.isEmpty() && augmentedRight.isEmpty());
    }

    /**    
    given a maiden (left) node, search forward links to 
    right nodes in alternating paths
    to insert right nodes into the heap
    or update their heap keys for shorter paths.
    
    note that any yNodes inserted into the heap internally,
    have their keys updated.
    
     * @param heap
     * @param forest
     * @param rM
     * @param yNodes
     * @param xNode 
     */
    private long scanAndAdd(MinHeapForRT2012 heap, Forest forest,
        ResidualDigraph rM, Map<Integer, RightNode> yNodes, 
        LeftNode xNode, long prevKey,
        Set<Integer> augmentedLeft,
        Set<Integer> augmentedRight,
        LeftNode topNode, Set<Integer> visitedY) {

        Integer xIndex = (Integer)(xNode.getData());
      
        if (augmentedLeft.contains(xIndex)) {
            return prevKey;
        }
        
        long lX = xNode.getKey();
        assert(lX < Long.MAX_VALUE);
        
        Set<Integer> forwardLinks = rM.getForwardLinksRM().get(xIndex);
        
        if (forwardLinks != null) {
        
            for (Integer yIndex : forwardLinks) {

                if (visitedY.contains(yIndex) ||
                    augmentedRight.contains(yIndex)) {
                    continue;
                }
                visitedY.add(yIndex);
                
                //link length = net cost of the edge 
                //    (lp(x ⇒ y) = cp(X,Y))
                //link length of backward link Y->X, is 0, lp(y ⇒ x) := 0.

                // l(x) and l(y) are the keys in the heap node

                RightNode yNode = yNodes.get(yIndex);
                long lOld = yNode.getKey();
                assert(((Integer)yNode.getData()).intValue() ==
                    yIndex.intValue());
                
                //L := l(x) + lp(x ⇒ y) 
                //   = l(x) + cp(x, y)
                //   = l(x) + c(x, y) − pd(x) + pd(y)
                // from page 41 regarding first steps in FlowAssign:  
                // "We first perform some initialization. 
                // Ignoring the edges weights, we use Hopcroft-Karp 
                // to look for some matching of size t."
                // so, c(x,y) which is the edge weight in Graph g,
                // is 1 here by this statement.
                int cp = 1;

                long ell = lX + cp;
                if (ell < lOld) {
                    if (lOld == Long.MAX_VALUE) {
                        assert(yNode instanceof RightNode);
                        yNode = (RightNode)yNode.copy();
                        yNode.pathPredecessor = xNode.copy();
                        if (yNode.pathPredecessor.topPredecessor != null) {
                            assert (yNode.pathPredecessor.topPredecessor.getData().equals(
                                topNode.getData()));
                            yNode.topPredecessor
                                = yNode.pathPredecessor.topPredecessor;
                        } else {
                            yNode.topPredecessor = topNode;
                        }
                        yNode.setKey(ell);
                        heap.insert(yNode);
                        log.fine(String.format("HEAP insert: %s",
                            yNode.toString()));
                    } else {
                        Integer prev = 
                            (yNode.pathPredecessor == null) ?
                            null :
                            (Integer)yNode.pathPredecessor.getData();
                        Integer top = 
                            (yNode.topPredecessor == null) ?
                            null :
                            (Integer)yNode.topPredecessor.getData();
                        assert(prev != null);
                        assert(top != null);
                        if (prev.intValue() != xIndex) {
                            yNode.pathPredecessor = xNode;
                            yNode.topPredecessor = xNode.topPredecessor;
                        }
                        heap.decreaseKey(yNode, ell);
                        log.fine(String.format("HEAP decr: %s",
                            yNode.toString()));
                    }
                }
            }
        }
        
        //add x to the forest;
        return forest.add(xNode, prevKey);
    }
    
    private void augmentPath(ResidualDigraph rM, Forest forest, 
        long foresIdx, Set<Integer> augmentedLeft, 
        Set<Integer> augmentedRight) {

        /*
        will extract the paths, that is the branches in the
        tree forest[forestIdx]
        and find a maximal set of compatible paths within
        those.
        
        could use kruskal's or prims to create the set of
        edges (created from building the minimum spanning tree).
        then augment using each edge, but 
        skipping those violating vertex-disjoint.
        
        since prim's runtime complexity is better for
        |v| < |E|, will use that.
        
        */
        
        // ---- extract the forest tree as edges ---
        
        DoubleLinkedCircularList tree = forest.get((int)foresIdx);

        assert(tree != null);
        
        Set<PairInt> edges = new HashSet<PairInt>();
        long n = tree.getNumberOfNodes();
        HeapNode node = tree.getSentinel();
        log.fine("number of branches in tree=" + n);
 
        int j = 0;
        while (j < n) {
            node = node.getLeft();
            
            List<PathNode> path = MinCostUnbalancedAssignment.extractNodes((PathNode)node);
            if (path.size() < 2) {
                ++j;
                continue;
            }
            //debugPath(path);
            //discard if not vertex disjoint from augmented sets
            boolean skip = false;
            List<PairInt> tmp = new ArrayList<PairInt>();
            for (int ii = 0; ii < (path.size() - 1); ++ii) {
                PathNode node1 = path.get(ii);
                PathNode node2 = path.get(ii + 1);
                // index1 is the left index of arc
                // index2 is the right index of the arc
                Integer index1, index2;
                if (node1 instanceof LeftNode) {
                    index1 = (Integer)node1.getData();
                    index2 = (Integer)node2.getData();
                } else {
                    index1 = (Integer)node2.getData();
                    index2 = (Integer)node1.getData();
                }
                if (augmentedLeft.contains(index1) ||
                    augmentedRight.contains(index2)) {
                    skip = true;
                    break;
                }
                tmp.add(new PairInt(index1.intValue(), index2.intValue()));
            }
            if (!skip) {
                edges.addAll(tmp);
            }
            ++j;
        }
        
        if (edges.isEmpty()) {
            return;
        }

        List<PairInt> edges2;
        if (edges.size() > 2) {
            //edges2 = filterUsingMST(edges);
            edges2 = filter(edges);
        } else {
            edges2 = new ArrayList<PairInt>(edges);
        }

        //augment M along path;

        //pg 11: "But our augmenting paths will be paths 
        //in an auxiliary graph called the residual digraph"

        //"We augment along a tight augmenting path 
        //by swapping the status of the edges that 
        //underlie its links, saturating the idle edges
        //and idling the saturated ones. This process 
        //increases the size of the matching by exactly 1. 
        //It marries off the maiden and bachelor at the 
        //ends of the augmenting path, thus reducing the
        //number of places where future augmenting paths 
        //can start or end.

        /*
        separating edits into 2 lists:
           undo "saturated" then make idle saturated
        */
        List<PairInt> undoSaturated = new ArrayList<PairInt>();
        List<PairInt> makeSaturated = new ArrayList<PairInt>();
        
        for (PairInt edge : edges2) {
            Integer leftIndex = Integer.valueOf(edge.getX());
            Integer rightIndex = Integer.valueOf(edge.getY());
            Integer bLeftIndex = rM.getBackwardLinksRM().get(
                rightIndex);
            if (bLeftIndex != null && bLeftIndex.equals(leftIndex)) {
                undoSaturated.add(edge);
            } else {
                makeSaturated.add(edge);
            }
        }
                
        for (PairInt edge : undoSaturated) {

            Integer leftIndex = Integer.valueOf(edge.getX());
            Integer rightIndex = Integer.valueOf(edge.getY());

            // remove backwards link and mapping
            rM.getBackwardLinksRM().remove(rightIndex);

            Set<Integer> rIndexes = rM.getForwardLinksRM().get(leftIndex);
            // create a forward link        
            if (rIndexes == null) {
                rIndexes = new HashSet<Integer>();
                rM.getForwardLinksRM().put(leftIndex, rIndexes);
            }
            rIndexes.add(rightIndex);

            log.fine("augmented to remove :" + leftIndex + " to " +
                rightIndex);
            
            augmentedLeft.add(leftIndex);
            augmentedRight.add(rightIndex);
        }
        
        Set<Integer> tmpAR = new HashSet<Integer>();
        Set<Integer> tmpAL = new HashSet<Integer>();
        
        for (PairInt edge : makeSaturated) {

            Integer rightIndex = Integer.valueOf(edge.getY());

            if (tmpAR.contains(rightIndex)) {
                continue;
            }
            
            Integer leftIndex = Integer.valueOf(edge.getX());
            
            if (tmpAL.contains(leftIndex)) {
                continue;
            }
            
            // assert saturated mapping for right node doesn't exist
            Integer v2 = rM.getBackwardLinksRM().get(rightIndex);
            /*if (v2 != null) {
                continue;
            }*/
            assert(v2 == null);            
            
            Set<Integer> rIndexes = rM.getForwardLinksRM().get(leftIndex);
        
            boolean forwardFound = (rIndexes != null) && rIndexes.contains(rightIndex);
        
            if (forwardFound) {
                // remove existing "idle" forward link
                rIndexes.remove(rightIndex);

                // create a backward link and matched mapping
                rM.getBackwardLinksRM().put(rightIndex, leftIndex);

                log.fine("augmented to add :" + leftIndex + " to " +
                    rightIndex);
                
                tmpAL.add(leftIndex);
                tmpAR.add(rightIndex);
            }
        }
        augmentedRight.addAll(tmpAR);
        augmentedLeft.addAll(tmpAL);
       
    }
    
    private String debug(Set<Integer> aL, Set<Integer> aR) {
        StringBuilder sb = new StringBuilder();
        sb.append("aL=[");
        for (Integer index : aL) {
            sb.append(index).append(",");
        }
        sb.append("]");
        sb.append(" aR=[");
        for (Integer index : aR) {
            sb.append(index).append(",");
        }
        sb.append("]");
        return sb.toString();
    }
    
    private List<PairInt> filterUsingMST(Set<PairInt> edges) {
        
        // ---- prepare edges as a single graph and
        //      adjacency map to give to prim's
        
        // re-number the vertexes and make an adjacency map or matrix
        // populate these from Set<PairInt> edges
        int nVertexes = 0;
        Map<Integer, Integer> leftToNew = new HashMap<Integer, Integer>();
        Map<Integer, Integer> revLeftToNew = new HashMap<Integer, Integer>();
        for (PairInt edge : edges) {
            Integer index = Integer.valueOf(edge.getX());
            if (!leftToNew.containsKey(index)) {
                Integer index2 = Integer.valueOf(nVertexes);
                leftToNew.put(index, index2);
                revLeftToNew.put(index2, index);
                nVertexes++;
            }
        }
        int nL = nVertexes;        
        Map<Integer, Integer> rightToNew = new HashMap<Integer, Integer>();
        Map<Integer, Integer> revRightToNew = new HashMap<Integer, Integer>();
        
        for (PairInt edge : edges) {
            Integer index = Integer.valueOf(edge.getY());
            if (!rightToNew.containsKey(index)) {
                Integer index2 = Integer.valueOf(nVertexes);
                rightToNew.put(index, index2);
                revRightToNew.put(index2, index);
                nVertexes++;
            }
        }
        
        Map<Integer, Set<PairInt>> adjCostMap =
            new HashMap<Integer, Set<PairInt>>();
        
        for (PairInt edge : edges) {
            Integer index1 = Integer.valueOf(edge.getX());
            index1 = leftToNew.get(index1);
            Integer index2 = Integer.valueOf(edge.getY());
            index2 = rightToNew.get(index2);
            
            Set<PairInt> set2 = adjCostMap.get(index1);
            if (set2 == null) {
                set2 = new HashSet<PairInt>();
                adjCostMap.put(index1, set2);
            }
            // using a cost of 1 for all edges
            set2.add(new PairInt(index2.intValue(), 1));
            
            set2 = adjCostMap.get(index2);
            if (set2 == null) {
                set2 = new HashSet<PairInt>();
                adjCostMap.put(index2, set2);
            }
            // using a cost of 1 for all edges
            set2.add(new PairInt(index1.intValue(), 1));
        }
                
        // use prim's mst to make a maximal set of edges
        PrimsMST prims = new PrimsMST();
        prims.calculateMinimumSpanningTree(
            nVertexes, adjCostMap); 
        int[] prev = prims.getPrecessorArray();
        log.fine("predecessors=" + Arrays.toString(prev));
        
        edges.clear();
        
        for (int idx2 = 0; idx2 < prev.length; ++idx2) {
            int idx1 = prev[idx2];
            if (idx1 == -1) {
                continue;
            }
            Integer index1, index2;
            if (idx2 < nL) {
                //[left] = right
                index1 = Integer.valueOf(idx2);
                index2 = Integer.valueOf(idx1);
            } else {
                //[right] = left
                index1 = Integer.valueOf(idx1);
                index2 = Integer.valueOf(idx2);
            }
            
            index1 = revLeftToNew.get(index1);            
            index2 = revRightToNew.get(index2);
         
            log.fine(" passed filter1: " + index1 + ":" + index2);
            
            edges.add(new PairInt(index1.intValue(),
                index2.intValue()));            
        }
        
        List<PairInt> edges2 = new ArrayList<PairInt>();
        int nIter = 0;
        int nc = 0;
        while ((nIter == 0) || (nc > 0)) {
            nIter++;
            Map<Integer, Integer> lF = new HashMap<Integer, Integer>();
            Map<Integer, Integer> rF = new HashMap<Integer, Integer>();
            for (PairInt edge : edges) {
                Integer index1 = Integer.valueOf(edge.getX());
                Integer index2 = Integer.valueOf(edge.getY());
                Integer count = lF.get(index1);
                if (count == null) {
                    lF.put(index1, Integer.valueOf(1));
                } else {
                    lF.put(index1, Integer.valueOf(count.intValue() + 1));
                }
                count = rF.get(index2);
                if (count == null) {
                    rF.put(index2, Integer.valueOf(1));
                } else {
                    rF.put(index2, Integer.valueOf(count.intValue() + 1));
                }
            }
            
            // add to edges2  edges w/
            // left indexes w/ frequency=1
            // right indexes w/ frequency=1
            // and remove those from edges
            
            nc = 0;
            
            for (Map.Entry<Integer, Integer> entry : lF.entrySet()) {
                if (entry.getValue().intValue() == 1) {
                    int idx1 = entry.getKey().intValue();
                    PairInt p0 = null;
                    for (PairInt p : edges) {
                        if (p.getX() == idx1) {
                            p0 = p;
                            break;
                        }
                    }
                    assert (p0 != null);
                    edges2.add(p0);
                    edges.remove(p0);
                    nc++;
                }
            }
            
            for (Map.Entry<Integer, Integer> entry : rF.entrySet()) {
                if (entry.getValue().intValue() == 1) {
                    int idx2 = entry.getKey().intValue();
                    PairInt p0 = null;
                    for (PairInt p : edges) {
                        if (p.getY() == idx2) {
                            p0 = p;
                            break;
                        }
                    }
                    if (p0 != null) {
                        edges2.add(p0);
                        edges.remove(p0);
                        nc++;
                    }
                }
            }
            log.fine("edges.size=" + edges.size() + " edges2.size="
                + edges2.size());
        }
        
        edges2.addAll(edges);

        return edges2;
    }
    
    private List<PairInt> filter(Set<PairInt> edges) {
        
        List<PairInt> edges2 = new ArrayList<PairInt>();
        int nIter = 0;
        int nc = 0;
        while ((nIter == 0) || (nc > 0)) {
            nIter++;
            Map<Integer, Integer> lF = new HashMap<Integer, Integer>();
            Map<Integer, Integer> rF = new HashMap<Integer, Integer>();
            for (PairInt edge : edges) {
                Integer index1 = Integer.valueOf(edge.getX());
                Integer index2 = Integer.valueOf(edge.getY());
                Integer count = lF.get(index1);
                if (count == null) {
                    lF.put(index1, Integer.valueOf(1));
                } else {
                    lF.put(index1, Integer.valueOf(count.intValue() + 1));
                }
                count = rF.get(index2);
                if (count == null) {
                    rF.put(index2, Integer.valueOf(1));
                } else {
                    rF.put(index2, Integer.valueOf(count.intValue() + 1));
                }
            }
            
            // add to edges2  edges w/
            // left indexes w/ frequency=1
            // right indexes w/ frequency=1
            // and remove those from edges
            
            nc = 0;
            
            for (Map.Entry<Integer, Integer> entry : lF.entrySet()) {
                if (entry.getValue().intValue() == 1) {
                    int idx1 = entry.getKey().intValue();
                    PairInt p0 = null;
                    for (PairInt p : edges) {
                        if (p.getX() == idx1) {
                            p0 = p;
                            break;
                        }
                    }
                    assert (p0 != null);
                    edges2.add(p0);
                    edges.remove(p0);
                    nc++;
                }
            }
            
            for (Map.Entry<Integer, Integer> entry : rF.entrySet()) {
                if (entry.getValue().intValue() == 1) {
                    int idx2 = entry.getKey().intValue();
                    PairInt p0 = null;
                    for (PairInt p : edges) {
                        if (p.getY() == idx2) {
                            p0 = p;
                            break;
                        }
                    }
                    if (p0 != null) {
                        edges2.add(p0);
                        edges.remove(p0);
                        nc++;
                    }
                }
            }
            log.fine("edges.size=" + edges.size() + " edges2.size="
                + edges2.size());
        }
        
        edges2.addAll(edges);

        return edges2;
    }
    
    private void debug(Forest forest) {
        
        for (Integer key : forest.getKeys()) {

            DoubleLinkedCircularList tree = forest.get(key.intValue());
            
            long n = tree.getNumberOfNodes();
            HeapNode node = tree.getSentinel();
            int j = 0;
            while (j < n) {
                node = node.getLeft();
                
                List<PathNode> path = extractNodes((PathNode)node);
                int n2 = path.size();
                for (int ii = 0; ii < n2; ++ii) {
                    PathNode node1 = path.get(ii);
                    log.info("forest[" + key + "] tree branch[" 
                        + j + "] node[" + ii + "]=" + node1.toString());
                }                    
                j++;
            }
        }            
    }
    
}
