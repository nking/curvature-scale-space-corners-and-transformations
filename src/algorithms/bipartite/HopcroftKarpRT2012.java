package algorithms.bipartite;

import algorithms.bipartite.MinCostUnbalancedAssignment.Forest;
import algorithms.bipartite.MinCostUnbalancedAssignment.LeftNode;
import algorithms.bipartite.MinCostUnbalancedAssignment.RightNode;
import algorithms.bipartite.MinCostUnbalancedAssignment.PathNode;
import static algorithms.bipartite.MinCostUnbalancedAssignment.extractNodes;
import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.HeapNode;
import algorithms.mst.PrimsMST;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 * NOTE:  this class is NOT READY FOR USE YET.
 * It needs an implementation of "multi-level buckets"
 * to reduce the minHeap's extractMin and should possibly
 * be refactored to use primitives.
 * 
 * @author nichole
 */
public class HopcroftKarpRT2012 {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * NOT READY FOR USE.
     runtime complexity O(m * sqrt(s)) where m is number of edges
     and s is the size of the matching whose target size
     may be less than the maximum matchable.
     * It needs an implementation of "multi-level buckets"
     * to reduce the minHeap's extractMin and should possibly
     * be refactored to use primitives.
     * @param g
     * @return 
     */
    public TIntIntMap findMaxMatching(Graph g, int s) {
    
        TIntIntMap m = new TIntIntHashMap();
                 
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
            
            TIntIntMap m2 = rM.extractMatchings();
                        
            //announce(M is a matching)
            
            log.fine("nIter=" + nIter + " m2.size=" + m2.size()
                + " m.size=" + m.size());
            
            if (m2.size() >= s) {
                return m2;
            }
            
            TIntIntMap tmpM = new TIntIntHashMap(m);
            
            TIntSet mR = new TIntHashSet(m.values());
            // if m2 has a match that does not conflict with m,
            // add it to m (== vertex disjoint)

            TIntIntIterator iter = m2.iterator();
            for (int i = m2.size(); i-- > 0;) {
                iter.advance();
                int key = iter.key();
                int value = iter.value();
                if (!tmpM.containsKey(key) && !mR.contains(value)) {
                    tmpM.put(key, value);
                }
            }
            
            log.fine("union size=" + tmpM.size());
            
            if (tmpM.size() >= s) {
                return tmpM;
            }
            
            // remove the intersection of m and m2
            iter = m.iterator();
            for (int i = m.size(); i-- > 0;) {
                iter.advance();
                int key = iter.key();
                if (m2.containsKey(key) &&
                    m2.get(key)== iter.value()) {
                    tmpM.remove(key);
                }
            }
            
            m = tmpM;
            log.fine("symmetric diff of m and m2.size=" + m.size());
            
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

        TIntObjectIterator<TIntSet> iter = rM.getForwardLinksRM().
            iterator();

        for (int i = rM.getForwardLinksRM().size(); i-- > 0;) {
            iter.advance();
            TIntSet set = iter.value();
            if (set.size() > 1) {
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
     * @param lambda the length to use for a counting sort of
     * augmenting path lengths
     */
    protected boolean buildForestAndAugment(
        final ResidualDigraph rM, int lambda) {
      
        Forest forest = new Forest(lambda);
        
        MinHeapForRT2012 heap = new MinHeapForRT2012(4, 
            rM.countOfForwardBipartiteLinks());
        
        /*
        this class is invoked as the first step in a
        hopcroft-karp matching algorithm and there are
        at first no matching nodes,
        so all X nodes are maiden nodes.
        but other uses may give the method a residual graph
        that has matched arcs in it.
        
        The code below treats each maiden X node as a single
        source shortest path problem.
           
        Note that book-keeping for the predecessors requires
        a copy be made upon insert of a node into the heap
        so that all possible foward links are present in 
        the heap.
        */
     
        // init all nodes to inf length
        PathNodes pathNodes = new PathNodes(rM.getNLeft(), 
            rM.getNRight());
        
        TIntObjectMap<LeftNode> leftNodes 
            = pathNodes.getLeftNodes();
        TIntObjectMap<RightNode> rightNodes 
            = pathNodes.getRightNodes();
       
        /*
        for the bfs and tracking "visited", need to track
        individually for each maiden as its own single shortest
        path.  doing so for the right nodes
        key = maiden node index (== X index)
        value = visited nodes along key path
        */
        TIntObjectMap<TIntSet> vXY = 
            new TIntObjectHashMap<TIntSet>();
        
        TIntSet augmentedLeft = new TIntHashSet();
        TIntSet augmentedRight = new TIntHashSet();
        long prevKey = -1;
        long lastAugKey = -1;
        
        // married X nodes
        TIntSet matchedLeft = new TIntHashSet(
            rM.getBackwardLinksRM().values());    
 
        // for all maidens
        // set key to 0, then ScanAndAdd(index)
        TIntSet maidens = new TIntHashSet();
        for (int i = 0; i < rM.getNLeft(); ++i) {
            if (matchedLeft.contains(i)) {
                continue;
            }
            maidens.add(i);
            LeftNode node = leftNodes.get(i);
            node.setKey(0);
            vXY.put(i, new TIntHashSet());
            prevKey = scanAndAdd(heap, forest, rM, 
                rightNodes, node, prevKey,
                augmentedLeft, augmentedRight, 
                node, vXY.get(i));
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
            PathNode y = (PathNode) heap.extractMin();
            assert(y instanceof RightNode);
            
            log.fine("heap.size=" + heap.getNumberOfNodes());
            
            log.fine("extractMin=" + y.toString());
            
            int yIdx = y.index;
        
            if (augmentedRight.contains(yIdx)) {
                continue;
            }
            
            assert(y.topPredecessor != null);
            
            int topIdx = y.topPredecessor.index;
            
            long currentKey = forest.add(y, prevKey);
            
            if (currentKey > prevKey) {
                log.fine("augment forest[" + prevKey + "] ");
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
            if (rM.getBackwardLinksRM().containsKey(yIdx)) {
                int xIdx = rM.getBackwardLinksRM().get(yIdx);
                LeftNode xNode = leftNodes.get(xIdx);
                
                xNode.setKey(y.getKey());
                if (xNode.pathPredecessor == null) {
                    xNode.pathPredecessor = y;
                    xNode.topPredecessor = y.topPredecessor;
                }
                 
                currentKey = scanAndAdd(heap, forest, 
                    rM, rightNodes, xNode, prevKey,
                    augmentedLeft, augmentedRight, y.topPredecessor, 
                    vXY.get(topIdx));
                
                if (currentKey > prevKey) {
                    log.fine("augment forest[" + prevKey + "] ");
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
                
                if (maidens.contains(topIdx)) {
                    maidens.remove(topIdx);                    
                } else if (maidens.isEmpty()) {
                    log.fine("last maiden's bachelor reached");
                    //debug(forest);
                    break;
                }
            }
        
            log.fine("nIter=" + nIter);
        }
        
        if (lastAugKey < prevKey) {
            log.fine("augment forest[" + prevKey + "] ");
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
        ResidualDigraph rM, TIntObjectMap<RightNode> yNodes, 
        LeftNode xNode, long prevKey,
        TIntSet augmentedLeft,
        TIntSet augmentedRight,
        LeftNode topNode, TIntSet visitedY) {

        int xIdx = xNode.index;
      
        if (augmentedLeft.contains(xIdx)) {
            return prevKey;
        }
        
        long lX = xNode.getKey();
        assert(lX < Long.MAX_VALUE);
        
        TIntSet forwardLinks 
            = rM.getForwardLinksRM().get(xIdx);
        
        if (forwardLinks != null) {
            TIntIterator iter = forwardLinks.iterator();
            while (iter.hasNext()) {
                int yIdx = iter.next();
                if (visitedY.contains(yIdx) ||
                    augmentedRight.contains(yIdx)) {
                    continue;
                }
                visitedY.add(yIdx);
                
                //link length = net cost of the edge 
                //    (lp(x ⇒ y) = cp(X,Y))
                //link length of backward link Y->X, is 0, lp(y ⇒ x) := 0.

                // l(x) and l(y) are the keys in the heap node

                RightNode yNode = yNodes.get(yIdx);
                long lOld = yNode.getKey();
                assert(yNode.index == yIdx);
                
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
                            assert (yNode.pathPredecessor.topPredecessor.index ==
                                topNode.index);
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
                        int prev = 
                            (yNode.pathPredecessor == null) ?
                            null :
                            yNode.pathPredecessor.index;
                        int top = 
                            (yNode.topPredecessor == null) ?
                            null :
                            yNode.topPredecessor.index;
                        assert(prev != -1);
                        assert(top != -1);
                        if (prev != xIdx) {
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
        long foresIdx, TIntSet augmentedLeft, 
        TIntSet augmentedRight) {

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
                int idx1, idx2;
                if (node1 instanceof LeftNode) {
                    idx1 = node1.index;
                    idx2 = node2.index;
                } else {
                    idx1 = node2.index;
                    idx2 = node1.index;
                }
                if (augmentedLeft.contains(idx1) ||
                    augmentedRight.contains(idx2)) {
                    skip = true;
                    break;
                }
                tmp.add(new PairInt(idx1, idx2));
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
            int leftIdx = edge.getX();
            int rightIdx = edge.getY();
            if (rM.getBackwardLinksRM().containsKey(rightIdx)) {
                int bLeftIdx = rM.getBackwardLinksRM().get(
                    rightIdx);
                if (bLeftIdx == leftIdx) {
                    undoSaturated.add(edge);
                } else {
                    makeSaturated.add(edge);
                }
            } else {
                makeSaturated.add(edge);
            }
        }
                
        for (PairInt edge : undoSaturated) {

            int leftIdx = edge.getX();
            int rightIdx = edge.getY();

            // remove backwards link and mapping
            rM.getBackwardLinksRM().remove(rightIdx);

            TIntSet rIndexes = rM.getForwardLinksRM().get(leftIdx);
            // create a forward link        
            if (rIndexes == null) {
                rIndexes = new TIntHashSet();
                rM.getForwardLinksRM().put(leftIdx, rIndexes);
            }
            rIndexes.add(rightIdx);

            log.fine("augmented to remove :" + leftIdx 
                + " to " + rightIdx);
            
            augmentedLeft.add(leftIdx);
            augmentedRight.add(rightIdx);
        }
        
        TIntSet tmpAR = new TIntHashSet();
        TIntSet tmpAL = new TIntHashSet();
        
        for (PairInt edge : makeSaturated) {

            int rightIdx = edge.getY();

            if (tmpAR.contains(rightIdx)) {
                continue;
            }
            
            int leftIdx = edge.getX();
            
            if (tmpAL.contains(leftIdx)) {
                continue;
            }
            
            // assert saturated mapping for right node doesn't exist
            assert(rM.getBackwardLinksRM().containsKey(rightIdx));
            int v2 = rM.getBackwardLinksRM().get(rightIdx);
            
            TIntSet rIndexes = rM.getForwardLinksRM().get(leftIdx);
        
            boolean forwardFound = (rIndexes != null) && 
                rIndexes.contains(rightIdx);
        
            if (forwardFound) {
                // remove existing "idle" forward link
                rIndexes.remove(rightIdx);

                // create a backward link and matched mapping
                rM.getBackwardLinksRM().put(rightIdx, leftIdx);

                log.fine("augmented to add :" + leftIdx + " to " +
                    rightIdx);
                
                tmpAL.add(leftIdx);
                tmpAR.add(rightIdx);
            }
        }
        augmentedRight.addAll(tmpAR);
        augmentedLeft.addAll(tmpAL);
       
    }
    
    private String debug(TIntSet aL, TIntSet aR) {
        StringBuilder sb = new StringBuilder();
        sb.append("aL=[");
        TIntIterator iter = aL.iterator();
        while (iter.hasNext()) {
            int idx = iter.next();
            sb.append(Integer.toString(idx)).append(",");
        }
        sb.append("]");
        sb.append(" aR=[");
        iter = aR.iterator();
        while (iter.hasNext()) {
            int idx = iter.next();
            sb.append(Integer.valueOf(idx)).append(",");
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
        TIntIntMap leftToNew = new TIntIntHashMap();
        TIntIntMap revLeftToNew = new TIntIntHashMap();
        for (PairInt edge : edges) {
            int idx = edge.getX();
            if (!leftToNew.containsKey(idx)) {
                leftToNew.put(idx, nVertexes);
                revLeftToNew.put(nVertexes, idx);
                nVertexes++;
            }
        }
        int nL = nVertexes;        
        TIntIntMap rightToNew = new TIntIntHashMap();
        TIntIntMap revRightToNew = new TIntIntHashMap();
        
        for (PairInt edge : edges) {
            int idx = edge.getY();
            if (!rightToNew.containsKey(idx)) {
                rightToNew.put(idx, nVertexes);
                revRightToNew.put(nVertexes, idx);
                nVertexes++;
            }
        }
        
        TIntObjectMap<TIntIntMap> adjCostMap =
            new TIntObjectHashMap<TIntIntMap>();
        
        for (PairInt edge : edges) {
            int idx1 = edge.getX();
            idx1 = leftToNew.get(idx1);
            int idx2 = edge.getY();
            idx2 = rightToNew.get(idx2);
            
            TIntIntMap map1 = adjCostMap.get(idx1);
            if (map1 == null) {
                map1 = new TIntIntHashMap();
                adjCostMap.put(idx1, map1);
            }
            // using a cost of 1 for all edges
            map1.put(idx2, 1);
            
            TIntIntMap map2 = adjCostMap.get(idx2);
            if (map2 == null) {
                map2 = new TIntIntHashMap();
                adjCostMap.put(idx2, map2);
            }
            // using a cost of 1 for all edges
            map2.put(idx1, 1);
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
            int t1, t2;
            if (idx2 < nL) {
                //[left] = right
                t1 = idx2;
                t2 = idx1;
            } else {
                //[right] = left
                t1 = idx1;
                t2 = idx2;
            }
            
            t1 = revLeftToNew.get(t1);            
            t2 = revRightToNew.get(t2);
         
            log.fine(" passed filter1: " + t1 + ":" + t2);
            
            edges.add(new PairInt(t1, t2));            
        }
        
        List<PairInt> edges2 = new ArrayList<PairInt>();
        int nIter = 0;
        int nc = 0;
        while ((nIter == 0) || (nc > 0)) {
            nIter++;
            TIntIntMap lF = new TIntIntHashMap();
            TIntIntMap rF = new TIntIntHashMap();
            for (PairInt edge : edges) {
                int idx1 = edge.getX();
                int idx2 = edge.getY();
                
                if (!lF.containsKey(idx1)) {
                    lF.put(idx1, 1);
                } else {
                    int count = lF.get(idx1);
                    lF.put(idx1, count + 1);
                }
                
                if (!rF.containsKey(idx2)) {
                    rF.put(idx2, 1);
                } else {
                    int count = rF.get(idx2);
                    rF.put(idx2, count + 1);
                }
            }
            
            // add to edges2  edges w/
            // left indexes w/ frequency=1
            // right indexes w/ frequency=1
            // and remove those from edges
            
            nc = 0;
            
            TIntIntIterator iter = lF.iterator();
            for (int i = lF.size(); i-- > 0;) {
                iter.advance();                
                if (iter.value() == 1) {
                    int idx1 = iter.key();
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
            
            TIntIntIterator iter2 = rF.iterator();
            for (int i = rF.size(); i-- > 0;) {
                iter2.advance();                
                if (iter2.value() == 1) {
                    int idx2 = iter2.key();
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
            TIntIntMap lF = new TIntIntHashMap();
            TIntIntMap rF = new TIntIntHashMap();
            for (PairInt edge : edges) {
                int idx1 = edge.getX();
                int idx2 = edge.getY();
                if (lF.containsKey(idx1)) {
                    lF.put(idx1, 1);
                } else {
                    int count = lF.get(idx1);
                    lF.put(idx1, count + 1);
                }
                
                if (rF.containsKey(idx2)) {
                    rF.put(idx2, 1);
                } else {
                    int count = rF.get(idx2);
                    rF.put(idx2, count + 1);
                }
            }
            
            // add to edges2  edges w/
            // left indexes w/ frequency=1
            // right indexes w/ frequency=1
            // and remove those from edges
            
            nc = 0;
            
            TIntIntIterator iter = lF.iterator();
            for (int i = lF.size(); i-- > 0;) {
                iter.advance();
                if (iter.value() == 1) {
                    int idx1 = iter.key();
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
            
            iter = rF.iterator();
            for (int i = rF.size(); i-- > 0;) {
                iter.advance();                
                if (iter.value() == 1) {
                    int idx2 = iter.key();
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
