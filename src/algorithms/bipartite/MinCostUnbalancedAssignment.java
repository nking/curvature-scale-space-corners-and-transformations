package algorithms.bipartite;

import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.Heap;
import algorithms.imageProcessing.HeapNode;
import algorithms.util.PairInt;
import algorithms.util.TrioInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

/**
 * A solver for the min-cost, unbalanced, bipartite
 * assignment problem that also uses weight scaling
 * to solve the perfect and imperfect
 * assignment problems, but not incremental
 * with a runtime complexity of
 <pre>
 O(m * sqrt(n) * log(n * C)) 
 where m is the number of edges (a.k.a. arcs) in the graph,
 n is the maximum number of nodes in the two graphs to be matched,
 C is a constant greater than or equal to the maximum 
 edge weight and is greater than or equal to 1.

 The weight scaling is what allows the algorithm to achieve
 much better performance with respect to s
 where s is the number of matched nodes.
 
 The code below follows the paper of
 "On Minimum-Cost Assignments in Unbalanced Bipartite Graphs"
 by Ramshaw and Tarjan, 2012
 for their algorithm, FlowAssign, Refine, and other methods
 used in the paper.
 </pre>
 * @author nichole
 */
public class MinCostUnbalancedAssignment {

    private Logger log = Logger.getLogger(this.getClass().getName());

    /*
    an alternating path is in the residual digraph, and 
       alternates between forward and backward.
    augmenting paths are in the residual digraph.
       an augmenting path starts at a maiden and ends at
       a bachelor (therefore, first and last links are
       forward links).  the links from Y to X are already
       matched (married).
   - an augmenting path takes forward steps to add a new 
     edge to the matched AND takes backward steps to 
     remove an edge from the matched.
   - an augmenting path is "tight" when all of the edges
     that underlie its links are tight, that is, have 
     zero net cost.
   - c(P) is the cost of augmenting path P.
     it's the sum of costs of the arcs of the path in 
     the flow network.
   - Let x_i_j = 1   if i is assigned to j,  else =0
     to augment along path p is to replace x_i_j by 1 for
     each arc (i,j) in P directed from X to Y, and to
     replace each x_i_j by 0 for each arc (j, i) in P
     directed from Y to X.
     - suppose x' is the matching obtained from x after
       augmenting along path P.
       the cost of x' is cx' = cx + c(P)

    Definition 2-7 on page 15 defines a proper pseudoflow
    of f and cp.
    */
    
    /**
     * class specializing a fibonacci heap node to identify
     * a left node, a.k.a. X node
     */
    private class LeftNode extends PathNode {
        public LeftNode() {
            this.id = "LeftNode";
        }        
        public LeftNode copy() {
            LeftNode node = new LeftNode();
            node.setKey(getKey());
            node.setData(getData());
            node.pathPredecessor = pathPredecessor;
            return node;
        }
    }
    
    public static abstract class PathNode extends HeapNode {
        PathNode pathPredecessor = null;
        public abstract PathNode copy();
        String id = "";
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append(id).append(" key=").append(Long.toString(getKey()));
            sb.append(" index=").append(getData().toString());
            PathNode prev = pathPredecessor;
            while (prev != null) {
                sb.append(" [prev=").append(prev.toString())
                    .append("]");
                prev = prev.pathPredecessor;
            }
            return sb.toString();
        }
    }
    
    /**
     * class specializing a fibonacci heap node to identify
     * a right node, a.k.a. Y node
     */
    private class RightNode extends PathNode {
        public RightNode() {
            this.id = "RightNode";
        }
        public RightNode copy() {
            RightNode node = new RightNode();
            node.setKey(getKey());
            node.setData(getData());
            node.pathPredecessor = pathPredecessor;
            return node;
        }
    }
    
    /**
     * class specializing a fibonacci heap node to identify
     * a source
     */
    private class SourceNode extends PathNode {
        public SourceNode() {
            this.id = "SourceNode";
        }        
        public SourceNode copy() {
            SourceNode node = new SourceNode();
            node.setKey(getKey());
            node.setData(getData());
            node.pathPredecessor = pathPredecessor;
            return node;
        }
    }
    /**
     * class specializing a fibonacci heap node to identify
     * a sink node
     */
    private class SinkNode extends PathNode {
        public SinkNode() {
            this.id = "SinkNode";
        }        
        public SinkNode copy() {
            SinkNode node = new SinkNode();
            node.setKey(getKey());
            node.setData(getData());
            node.pathPredecessor = pathPredecessor;
            return node;
        }
    }
    
    private int estimateLambda(ResidualDigraph rM) {
        //TODO: this may need to be revised
        int n = 0;
        for (Entry<Integer, Set<Integer>> entry :
            rM.getForwardLinksRM().entrySet()) {
            if (entry.getValue().size() > 1) {
                n++;
            }
        }
        return (n + 2);
    }

    /*
    input to FlowAssign is a bipartite graph G with integral 
    edge weights and a target size t. 
    
    s := min(t, ν(G))
    initialization uses Hopcroft-Karp and ignores the edge weights
    to find any matching of size t.
    If t exceeds nu(G) a warning is logged and a matching
    of size s is returned.
    (a matching of size s in the graph G corresponds to an 
    integral flow f of value |f| = s in the flow network N_G)
    
    
    input to Refine is a flow f of value s and prices p
    that together make all arcs (qε)-proper, where q is an integer 
    parameter. Refine builds a new flow f′, also of value s, 
    and prices p′ that together make all arcs ε-proper. 
    
    Note:  In FlowAssign, saturated edges will be kept proper, 
       but not necessarily tight. So, for a saturated edge 
       (x, y), we will know only that cp(x, y) ≤ 0, and we 
       will need to define lp(y ⇒ x) := −cp(x, y), to keep 
       our link lengths nonnegative. That is, for an 
       alternating path A, we will define
       lp(A):= summation over arcs X->Y (cp(X,Y)) -
                  summation over arcs Y->X (cp(X,Y))
    
    FlowAssign (G, t)
      (M, s) := HopcroftKarp(G, t);
      convert M into an integral flow f on N_G with |f| = s; 
      set ε := ε and, 
      for all nodes v in N_G, 
        set pd(v) := 0; 
      while ε > ε do
        ε := ε/q;
        Refine(f,p,ε); od;
      round prices to integers that make all arcs proper; 
    
    
    Refine(f,p,ε)
      S := {the s women who are matched in f};
      D := {the s men who are matched in f};
      convert the s bipartite arcs that are saturated in 
        f to idle; 
      raise the prices p, as in Figure 7.4, to make all arcs 
        ε-proper; 
      int h := s;
      while h > 0 do
        build a shortest-path forest from the current surpluses S, 
          stopping when a current deficit in D is reached;
        raise prices at forest nodes by multiples of ε, 
          shortening the discovered augmenting path to length 0;
        find a maximal set P of length-0 augmenting paths 
          that are compatible, as defined in Section 8.3;
        augment f along each of the paths in P in turn, thereby 
          reducing |S| = |D| = h by |P|;
      od;
      
    */

    /**
     * 
     * @param g bipartite graph with integral weights.
     * @return 
     */
    public Map<Integer, Integer> flowAssign(Graph g) {

        int sz = Math.min(g.getNLeft(), g.getNRight());
        
        Map<Integer, Integer> m = hopcroftKarp(g, sz);
        
        g = g.copyToCreateSourceSink();
        
        FlowNetwork gFlow = new FlowNetwork(g, m);
       
        //TODO: estimate of eps may need to be revised
        
        int s = m.size();
         
        // q >= 2.  consider q = O(log_2(n)) pg 41 par 3.
        int q = 2;
        
        // pg 44
        // since all costs are integers, can set eps < 1/6s
        // where s is size of m
        float eps = 1.f/(6.f * (float)s);
        // looking ahead at link lengths, one needs eps <= 0.5
        // in order to distinguish between the smallest
        // difference of integer costs, which is 1.
        
        // the loop conitional requires the initial eps
        // to be larger than epsBar.
        // TODO: calculate eps_bar correctly...
        // until then, just using a counter and limit
        // of the rough expected runtime
        int rIter = (int)(Math.log(s * gFlow.getMaxC())/Math.log(q));
        
        log.info("eps=" + eps + " rIter=" + rIter + " maxC=" + 
            gFlow.getMaxC());
               
        // all nodes V in gFlow have prices = 0
        
        int nIterR = 0;
                
        while (nIterR < rIter) {
            
            // pg 44, assertions I1, I2, I3, and I4
            assert(gFlow.assertFlowValue(s));
            assert(gFlow.assertPricesAreQuantizedEps(eps));
            assert(gFlow.integralFlowIsEpsProper(eps));
            assert(gFlow.assertSaturatedBipartiteIsEpsSnug(eps));
            
            eps /= (float)q;
            
            refine(gFlow, s, eps, q);
        
            ++nIterR;
        }
        
        // assert nIter is approx log_q(s * maxC)
        
        // round prices to integers that make all arcs proper
        roundFinalPrices(gFlow);
        
        return m;
    }
    
    protected void refine(FlowNetwork gFlow, int s, float eps,
        int q) {
        
        // S = left nodes matched in gFlow
        List<Integer> surplus = new ArrayList<Integer>();
        
        // D = right nodes matched in gFlow
        List<Integer> deficit = new ArrayList<Integer>();
        
        gFlow.getMatchedLeftRight(surplus, deficit);
        
        // set the flow of saturated bipartite arcs to 0
        for (int i = 0; i < surplus.size(); ++i) {
            int idx1 = surplus.get(i);
            int idx2 = deficit.get(i);
            gFlow.getFlow().put(new PairInt(idx1, idx2), 
                Float.valueOf(0));
        }
        
        /*
        see Figure 7.4 on pg 53.
        raise prices so that every arc in 
        gFlow becomes eps-proper, for the resulting pseudoflow 
        f and for the new, smaller value of eps.
        */
        gFlow.raisePricesUntilEpsProper(eps, q);
        
        //in [0] holds the length of the terminating deficit
        // node which is also the forest index it was
        // added to.
        int[] terminatingDeficitIdx = new int[1];
        int h = 2;
        int nHIter = 0;
        while (h > 0) {
          
            log.info("nHIter=" + nHIter);
            
            ResidualDigraph2 rF = new ResidualDigraph2(gFlow);

            // build a shortest-path forest from the current surpluses S, 
            // stopping when a current deficit in D is reached;
            // (pg 55)
            DoubleLinkedCircularList[] forest = 
                buildForest2(gFlow, rF, surplus, deficit, eps,
                terminatingDeficitIdx);
         
            log.info("l(termIndex)=" + terminatingDeficitIdx[0]);
         
            debug(forest);
           
            // raise prices at forest nodes by multiples of ε, 
            // shortening the discovered augmenting path to length 0;
            // array i[v] for all nodes in V (==left nodes)
            // if in forest: i(v) := l(termDefIdx) - l(v), else i(v)=0.
            // then for all nodes v: pd'(v) = pd(v) + i(v)*eps, 
            // ==> cp'(v, w) = cp(v, w) + (i(w) - i(v))*eps
            // ==> raising the price pd(v) at *some* node v in NG 
            //     by eps 
            //     lowers by 1 the length of any link in the 
            //     residual digraph Rf that leaves v 
            //     and raises by 1 the length of any link that enters v.
            // ==> ?? ignoring backwardLinks w->v???
            
            // NOTE: w.r.t. v nodes, I'm assuming that they are only
            //       those in the forest or connected to it,
            //       because if not, all combinations of links
            //       would need to be calculated which would defeat
            //       one of the main points of buildForest2 which
            //       stops the calculations early upon first deficit node.
            
            // --- for price and path length changes:
            //
            //could traverse forest and make maps
            //    key    value
            //    v->w   list of locations as pairints of forest/branch
            //    w->v       same
            //    source->v  same
            //    w->sink    same
            //then traverse forest making price changes in FlowNetwork
            //    and making link length changes to the
            //    forest nodes by looking up where they
            //    are in the new maps.
            //    the instructions dont say, but should 
            //    any changes be applied to nodes further
            //    down a branch?
            //then from the forest, make an adjacency list
            //    (and lists to include source and sink?)
            //    of the links with lengths = 0;
            
            
            int[] incrLeft = new int[gFlow.getNLeft()];
            int[] incrRight = new int[gFlow.getNRight()];

            Map<Integer, Set<TrioInt>> leftToRightLoc
                = new HashMap<Integer, Set<TrioInt>>();
            Map<Integer, Set<TrioInt>> rightToLeftLoc
                = new HashMap<Integer, Set<TrioInt>>();
            Map<Integer, TrioInt> sourceToLeftLoc
                = new HashMap<Integer, TrioInt>();
            Map<Integer, TrioInt> leftToSourceLoc
                = new HashMap<Integer, TrioInt>();
            Map<Integer, TrioInt> rightToSinkLoc
                = new HashMap<Integer, TrioInt>();
            Map<Integer, TrioInt> sinkToRightLoc
                = new HashMap<Integer, TrioInt>();
            
            Map<Integer, Map<Integer, List<PathNode>>> forestTreeMap = 
                new HashMap<Integer, Map<Integer, List<PathNode>>>();
            
            int lt = terminatingDeficitIdx[0];
            
            // traverse trees in the forest
            for (int forestIdx = 0; forestIdx < forest.length; ++forestIdx) {
                DoubleLinkedCircularList tree = forest[forestIdx];
                if (tree == null) {
                    continue;
                }
                long n = tree.getNumberOfNodes();
                HeapNode node = tree.getSentinel();
                int treeIdx = 0;
                while (treeIdx < n) {
                    node = node.getRight(); 
                    List<PathNode> path = extractNodes(node);
                    insert(forestTreeMap, forestIdx, treeIdx, path);
                    for (int branchIdx = 0; branchIdx < (path.size() - 1); ++branchIdx) {
                        PathNode node2 = path.get(branchIdx);
                        PathNode node1 = path.get(branchIdx + 1);
                        int l1 = (int)node1.getKey();
                        int l2 = (int)node2.getKey();
                        int idx1 = ((Integer)node1.getData()).intValue();
                        int idx2 = ((Integer)node2.getData()).intValue();
                        if (node1 instanceof LeftNode) {
                            assert(!(node2 instanceof LeftNode));
                            //i(v) := l(termDefIdx) - l(v)
                            incrLeft[idx1] = lt - l1;
                            if (node2 instanceof RightNode) {
                                insert(leftToRightLoc, node1, node2,
                                    forestIdx, treeIdx, branchIdx);
                                incrRight[idx2] = lt - l2;
                            } else {
                                assert(node2 instanceof SourceNode);
                                assert(leftToSourceLoc.get(
                                    (Integer)node1.getData()) == null);
                                leftToSourceLoc.put(
                                    (Integer)node1.getData(),
                                    new TrioInt(forestIdx, treeIdx, branchIdx));
                                incrLeft[idx2] = lt - l2;
                            }
                        } else if (node1 instanceof RightNode) {
                            assert(!(node2 instanceof RightNode));
                            //i(v) := l(termDefIdx) - l(v)
                            incrRight[idx1] = lt - l1;
                            if (node2 instanceof LeftNode) {
                                insert(rightToLeftLoc, node1, node2,
                                    forestIdx, treeIdx, branchIdx);
                                incrLeft[idx2] = lt - l2;
                            } else {
                                assert(node2 instanceof SinkNode);
                                assert(rightToSinkLoc.get(
                                    (Integer)node1.getData()) == null);
                                rightToSinkLoc.put(
                                    (Integer)node1.getData(),
                                    new TrioInt(forestIdx, treeIdx, branchIdx));
                                incrRight[idx2] = lt - l2;
                            }
                        } else if (node1 instanceof SourceNode) {
                            assert(node2 instanceof LeftNode);
                            assert(sourceToLeftLoc.get(
                                (Integer)node2.getData()) == null);
                            sourceToLeftLoc.put(
                                (Integer)node2.getData(),
                                new TrioInt(forestIdx, treeIdx, branchIdx));                            
                            //i(v) := l(termDefIdx) - l(v)
                            incrLeft[idx1] = lt - l1;
                            incrLeft[idx2] = lt - l2;
                        } else {
                            //node1 is a sink node
                            assert(node2 instanceof RightNode);
                            assert(sinkToRightLoc.get(
                                (Integer)node2.getData()) == null);
                            sinkToRightLoc.put(
                                (Integer)node2.getData(),
                                new TrioInt(forestIdx, treeIdx, branchIdx));
                            //i(v) := l(termDefIdx) - l(v)
                            incrRight[idx1] = lt - l1;
                            incrRight[idx2] = lt - l2;
                        }
                    }
                    treeIdx++;
                }
            }
        
            //pd'(v) = pd(v) + i(v)*eps,
            for (int i = 0; i < gFlow.getNLeft(); ++i) {
                float incr = eps * incrLeft[i];
                if (incr > 0) {
                    gFlow.addToLeftPrice(i, incr);
                    // -1 to links leaving v
                    // +1 to links entering v
                    
                    Integer index = Integer.valueOf(i);
                    // left to right
                    Set<TrioInt> set = leftToRightLoc.get(index);
                    if (set != null) {
                        for (TrioInt loc : set) {
                            PathNode node = forestTreeMap
                                .get(Integer.valueOf(loc.getX()))
                                .get(Integer.valueOf(loc.getY()))
                                .get(loc.getZ());
                            long key = node.getKey();
                            node.setKey(key - 1);
                        }
                    }
                    if (leftToSourceLoc.containsKey(index)) {
                        TrioInt loc = leftToSourceLoc.get(index);
                        PathNode node = forestTreeMap
                            .get(Integer.valueOf(loc.getX()))
                            .get(Integer.valueOf(loc.getY()))
                            .get(loc.getZ());
                        long key = node.getKey();
                        node.setKey(key - 1);
                    }
                    
                    //links entering v
                    set = rightToLeftLoc.get(index);
                    if (set != null) {
                        for (TrioInt loc : set) {
                            PathNode node = forestTreeMap
                                .get(Integer.valueOf(loc.getX()))
                                .get(Integer.valueOf(loc.getY()))
                                .get(loc.getZ());
                            long key = node.getKey();
                            node.setKey(key + 1);
                        }
                    }
                    if (sourceToLeftLoc.containsKey(index)) {
                        TrioInt loc = sourceToLeftLoc.get(index);
                        PathNode node = forestTreeMap
                            .get(Integer.valueOf(loc.getX()))
                            .get(Integer.valueOf(loc.getY()))
                            .get(loc.getZ());
                        long key = node.getKey();
                        node.setKey(key + 1);
                    }
                }
            }
            for (int i = 0; i < gFlow.getNRight(); ++i) {
                Integer index = Integer.valueOf(i);
                float incr = eps * incrRight[i];
                if (incr > 0) {
                    gFlow.addToRightPrice(i, incr);
                    // -1 to links leaving v
                    // +1 to links entering v
                    Set<TrioInt> set = rightToLeftLoc.get(index);
                    if (set != null) {
                        for (TrioInt loc : set) {
                            PathNode node = forestTreeMap
                                .get(Integer.valueOf(loc.getX()))
                                .get(Integer.valueOf(loc.getY()))
                                .get(loc.getZ());
                            long key = node.getKey();
                            node.setKey(key - 1);
                        }
                    }
                    if (rightToSinkLoc.containsKey(index)) {
                        TrioInt loc = rightToSinkLoc.get(index);
                        PathNode node = forestTreeMap
                            .get(Integer.valueOf(loc.getX()))
                            .get(Integer.valueOf(loc.getY()))
                            .get(loc.getZ());
                        long key = node.getKey();
                        node.setKey(key - 1);
                    }
                    
                    //links entering v
                    set = leftToRightLoc.get(index);
                    if (set != null) {
                        for (TrioInt loc : set) {
                            PathNode node = forestTreeMap
                                .get(Integer.valueOf(loc.getX()))
                                .get(Integer.valueOf(loc.getY()))
                                .get(loc.getZ());
                            long key = node.getKey();
                            node.setKey(key + 1);
                        }
                    }
                    if (sinkToRightLoc.containsKey(index)) {
                        TrioInt loc = sinkToRightLoc.get(index);
                        PathNode node = forestTreeMap
                            .get(Integer.valueOf(loc.getX()))
                            .get(Integer.valueOf(loc.getY()))
                            .get(loc.getZ());
                        long key = node.getKey();
                        node.setKey(key + 1);
                    }
                }
            }
            
            //assert from pg 52 
            //       I1', I2, I3, I4, on FlowNetwork
            //       and I5 on ResidualDigraph2
            
            assert(gFlow.assertFlowValueIncludingSrcSnk(s));
            assert(gFlow.assertPricesAreQuantizedEps(eps));
            assert(gFlow.integralFlowIsEpsProper(eps));
            assert(gFlow.assertSaturatedBipartiteIsEpsSnug(eps));
            // assert I5: Rf has no cycles of length zero.
            
            log.info("in refine after prices and link length changes");
            debug(forest);
            
            //Let Rf0 denote that subgraph of the residual digraph 
            //formed by links of length zero.
           
            
            
            // --- Sect 8.3, create maximal set of compatible augmenting paths
            //
            //then apply pseudocode from pg 62 , Figure 8.2
            //   with input = length 0 adj list, the surplus
            //   list and the deficit list
            //   to create the maximal set of compatible paths.
            //
            //   state that needs to be tracked for a vertex:
            //    - visited (== marked)
            //    - identity (== Left or Right or Source or Sink)
            
            
            //augment f along each of the paths in P in turn, thereby 
            //   reducing |S| = |D| = h by |P|;
            
        if (true)
          throw new UnsupportedOperationException("not yet implemented");
            
            ++nHIter;
        }
    }
    
    /**
     * 
     * @param gFlow
     * @param surplus
     * @param deficit
     * @param eps
     * @param terminatingDeficitIdx output array of size 1
     * to hold th length of the key of the terminating deficit
     * node (which is also the index of the forest tree it 
     * was added to).
     * @return 
     */
    protected DoubleLinkedCircularList[] buildForest2(
        final FlowNetwork gFlow, ResidualDigraph2 rF,
        List<Integer> surplus,
        List<Integer> deficit, float eps,
        int[] terminatingDeficitIdx) {
        
        Set<Integer> d = new HashSet<Integer>(deficit);
    
        //TODO: revisit this.
        int lambda = 3 * Math.min(gFlow.getNLeft(), 
            gFlow.getNRight());
        
        DoubleLinkedCircularList[] forest 
            = new DoubleLinkedCircularList[lambda];

        //TODO: revisit this
        // need to use a limit which will always include
        // the shortest path length at any time.
        // they're quantized w/ eps.
        DoubleLinkedCircularList[] minHeap = new DoubleLinkedCircularList[lambda];
            
        Map<Integer, LeftNode> leftNodes = new HashMap<Integer, LeftNode>();
        Map<Integer, RightNode> rightNodes = new HashMap<Integer, RightNode>();
        for (int i = 0; i < gFlow.getNLeft(); ++i) {
            Integer index = Integer.valueOf(i);
            LeftNode node = new LeftNode();
            node.setKey(Long.MAX_VALUE);
            node.setData(index);
            leftNodes.put(index, node);
        }
        // NOTE: pseudocode does not suggest init of Right to inf
        for (int i = 0; i < gFlow.getNRight(); ++i) {
            Integer index = Integer.valueOf(i);
            RightNode node = new RightNode();
            node.setKey(Long.MAX_VALUE);
            node.setData(index);
            rightNodes.put(index, node);
        }        
        
        SourceNode sourceNode = new SourceNode();
        sourceNode.setKey(Long.MAX_VALUE);
        sourceNode.setData(Integer.valueOf(gFlow.getNLeft()));
        
        SinkNode sinkNode = new SinkNode();
        sinkNode.setKey(Long.MAX_VALUE);
        sinkNode.setData(Integer.valueOf(gFlow.getNRight()));
       
        for (Integer sigma : surplus) {
            LeftNode sNode = leftNodes.get(sigma);
            sNode.setKey(0);
            insertIntoHeap(minHeap, sNode);
        }
               
        do {
            PathNode node1 = extractMinFromHeap(minHeap);
            Integer index1 = (Integer)node1.getData();
            int idx1 = index1.intValue();
         
            log.info("extractMin=" + node1.toString());
            
            boolean node1IsLeft = (node1 instanceof LeftNode);
            boolean nodeIsSource = (node1 instanceof SourceNode);
            boolean nodeIsSink = (node1 instanceof SinkNode);
            
            //scan:
            if (nodeIsSource) {        
                // scan forward source links
                for (Integer index2 : rF.getForwardLinksSourceRM()) {
                    handleSourceForwardLink(minHeap, gFlow, 
                        (SourceNode)node1, index2, 
                        leftNodes, lambda, eps); 
                }
            } else if (nodeIsSink) {
                // scan backward sink links
                for (Integer index2 : rF.getBackwardLinksSinkRM()) {
                    handleSinkBackwardLink(minHeap, gFlow, 
                        (SinkNode)node1, index2, rightNodes, lambda, eps);
                }
            } else if (node1IsLeft) {
                // scan the bipartite arcs forward
                Set<Integer> indexes2 = rF.getForwardLinksRM().get(
                    index1);
                if (indexes2 != null) {
                    for (Integer index2 : indexes2) {
                        handleRight(minHeap, gFlow, node1, index2, 
                            rightNodes, lambda, eps);                        
                    }
                }
           
                // if there's a source link
                if (rF.getBackwardLinksSourceRM().contains(index1)) {
                    // insert a copy of the source node
                    SourceNode sNode2 = sourceNode.copy();
                    handleSourceBackwardLink(minHeap, gFlow, 
                        node1, sNode2, lambda, eps);
                }
            } else {
                // node1 is a RightNode
                Integer index2 = rF.getBackwardLinksRM().get(index1);
                if (index2 != null) {
                    handleLeft(minHeap, gFlow, node1, index2, 
                       leftNodes, lambda, eps);                     
                }
                // if there is a sink link
                if (rF.getForwardLinksSinkRM().contains(index1)) {
                    // insert a copy of the sink node
                    SinkNode sNode2 = sinkNode.copy();
                    handleSink(minHeap, gFlow, node1, sNode2, 
                        lambda, eps);
                }
            }
                
            PathNode node1Cp = node1.copy();

            log.info("add to forest key=" + node1Cp.toString());
            
            //add v to the forest;
            addToForest(forest, node1Cp);

            HeapNode rootOfTD = forest[(int)node1Cp.getKey()].getSentinel()
                .getRight();
            Integer rootIndex = (Integer)rootOfTD.getData();
            terminatingDeficitIdx[0] = (int)node1Cp.getKey();
                
            if (d.contains(index1) && !node1IsLeft) {
                break;
            }
            
        } while (true);
        
        /*
        link lengths in residual digrph 2:
           - a forward link v->w has length
               lp(v->w) = Math.ceil(cp(v,w)/eps)
               (an idle arc will have lp(v->w) >= 0)
           - a backward link w->v has length
               lp(w->v) = 1 - Math.ceil(cp(v,w)/eps)
               (value will be >= 0)
        */
        
        return forest;
    }
    
    /**
     * NOTE: method does not yet account for the finer details 
     * of matching size from pg 41, paragraph 3.
     * find a maximal matching of size s for the left
     * and right nodes in the bipartite graph g
     * @param g
     * @param s
     * @return 
     */
    protected Map<Integer, Integer> hopcroftKarp(Graph g,
        int s) {
        
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
        
        if (false) {             
            //temporarily, replacing w/ O(m * sqrt(n))
            HopcroftKarp hk = new HopcroftKarp();
            int[] matched = hk.hopcroftKarpV0(new GraphWithoutWeights(g));
            log.info("matched=" + Arrays.toString(matched));
            for (int i = 0; i < matched.length; ++i) {
                int v = matched[i];
                if (v > -1) {
                    m.put(Integer.valueOf(i), Integer.valueOf(v));
                }
            }
            return m;
        }
                
        ResidualDigraph rM = new ResidualDigraph(g, m);
                
        return hopcroftKarp(g, rM, s);
    }
    
    /**
     * NOTE: this method is not yet throughly tested
     * @param g
     * @param rM
     * @param s
     * @return 
     */
    protected Map<Integer, Integer> hopcroftKarp(Graph g, 
        ResidualDigraph rM, int s) {
       
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
        
        int prevMSize = m.size();
        
        int nIter = 0;
            
        // estimate lambda for the dial array length by
        // looking at the number of edges with more than
        // one connection.
        int lambda = estimateLambda(rM);

        while (true) {
            
            DoubleLinkedCircularList[] augmentingPaths =
                buildForest(rM, lambda);

            if (allAreEmpty(augmentingPaths)) {
                return m;
            }
                       
            Set<PairInt> augmented = new HashSet<PairInt>();

            // start at i=1, because i=0 is not alternating? 
            for (int i = 1; i < augmentingPaths.length; ++i) {

                DoubleLinkedCircularList tree = augmentingPaths[i];

                if (tree == null) {
                    continue;
                }

                //augment M along path;
                /*
                pg 11: "But our augmenting paths will be paths 
                in an auxiliary graph called the residual digraph"

                "We augment along a tight augmenting path 
                by swapping the status of the edges that 
                underlie its links, saturating the idle edges
                and idling the saturated ones. This process 
                increases the size of the matching by exactly 1. 
                It marries off the maiden and bachelor at the 
                ends of the augmenting path, thus reducing the
                number of places where future augmenting paths 
                can start or end.
                */

                // use getRight() to traverse LIFO order which
                //   makes removing redundant paths easier
                //   currently.
                
                // each item in tree is a HeapNode whose path
                //    is included as pathPredecessor of the node.

                long n = tree.getNumberOfNodes();
                HeapNode node = tree.getSentinel();
                int j = 0;
                while (j < n) {
                    
                    node = node.getRight();
                    
                    /*
                    node is a path.
                    
                    The following is a path example of length 1. 
                    XB is current HeapNode and it is a LeftNode.
                    Its pathPedecessor is YA and it's a RightNode.
                    So, for this node, we ascend until pathPredecessor is null.
                    XA
                       \
                         YA
                       /
                    XB
                    
                    */
                    List<PathNode> path = extractNodes(node);

                    for (int ii = 0; ii < (path.size() - 1); ++ii) {
                        
                        PathNode node1 = path.get(ii);
                        PathNode node2 = path.get(ii + 1);
                        
                        log.fine("forest[" + i + "] tree branch[" 
                            + j + "] node[" + ii + "]=" + node1.toString());
                        
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
                        
                        PairInt pair = new PairInt(index1.intValue(),
                            index2.intValue());
                        if (augmented.contains(pair)) {
                            break;
                        }
                        swapLinkExistence(rM, index1, index2);
                        augmented.add(pair);
                    }                    
                    j++;
                }
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
            log.info("nIter=" + nIter + " m2.size=" + m2.size());

            // debug:
            for (Entry<Integer, Integer> entry : m2.entrySet()) {
                log.info("m2 match= " + entry.getKey() + "->" + entry.getValue());
            }
            for (Entry<Integer, Integer> entry : m.entrySet()) {
                log.info("m match= " + entry.getKey() + "->" + entry.getValue());
            }
                        
            // remove the intersection of m and m2
            for (Entry<Integer, Integer> entry : m.entrySet()) {
                Integer key = entry.getKey();
                if (m2.containsKey(key) &&
                    m2.get(key).equals(entry.getValue())) {
                    m2.remove(key);
                }
            }
            Set<Integer> m2R = new HashSet<Integer>(m2.values());
            // if m has a match that does not conflict with m2,
            // add it to m2 (== vertex disjoint)
            for (Entry<Integer, Integer> entry : m.entrySet()) {
                Integer key = entry.getKey();
                Integer value = entry.getValue();
                if (!m2.containsKey(key) && !m2R.contains(value)) {
                    m2.put(key, value);
                }
            }
            m = m2;
            log.info("symmetric diff of m and m2.size=" + m.size());
            for (Entry<Integer, Integer> entry : m.entrySet()) {
                log.info("sym diff m match= " + entry.getKey() + "->" + entry.getValue());
            }
            
            if (m.size() >= s) {
                return m;
            }
            
            rM = new ResidualDigraph(g, m);
                        
            assert (prevMSize < m.size());
            prevMSize = m.size();
            ++nIter;            
        }        
    }
    
    private void roundFinalPrices(FlowNetwork gFlow) {
        //see pg 46
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    /*
    NOTES
    ------
    G has vertices and edges 
       (and edge is usually written (X,Y)),
    N_G, the flow network, has nodes and arcs 
       (and arc is written X → Y, and all arcs go forward.) 
    residual digraph has nodes and links, 
       (and link is written X ⇒ Y, and some links 
       go forward while others go backward.)
       - Given a matching M in G, the residual digraph 
         R_M has nodes and links that correspond precisely 
         to the vertices and edges of G, 
         but each saturated edge (X, Y) 
           becomes a backward-directed link Y ⇒ X 
           in the residual digraph R_M, 
         while each idle edge (X, Y) becomes a 
           forward-directed link X ⇒ Y. 
         That is, idle arcs (f=0) become left-to-right links in R_M, 
         while saturated arcs (f=1) become right-to-left links. 
         Note that the residual digraph R_M depends only on the 
         matching M, not on the prices p.    
    
    
       - length or distance in the residual digraph
         is a term usable when the link has length >= 0
         and the underlying edges are "tight" (have 0 net cost).
   **    - to find the shortest possible augmenting path,
           can use a variant of Dijkstra's algorithm.   
           Note that if X->Y is a forward link R_M, corresponding
           to an idle edge (x,y), 
              define the length of that link 
                to be the net cost of the edge (lp(x ⇒ y) = cp(X,Y)).
              define the length of any backward link y ⇒ x 
                to be zero: lp(y ⇒ x) := 0. 
              define the length of an alternating path to be 
                 the sum of the lengths of its links.
              Note that the term link-count is in contrast, the
                number of links regardless of edge cost.
    
    A pseudoflow is a flux in which the flow f(X,Y) along each 
      arc X → Y is nonnegative and satisfies the unit-capacity 
      constraint: 0 ≤ f (X, Y) ≤ 1. 
      - an arc is idle in f if it's a pseudoflow w/ f(X,Y) = 0 
      - an arc is saturated if f(X,Y) = 1, (cp(x, y) ≤ 0)
      - else arc is fractional flow with 0 < f(X,Y) < 1  
The value of a flow f, denoted |f|, is the total 
    flow out of the source, which is also the total 
    flow into the sink – and, for that matter, 
    the total flow over all of the bipartite arcs.
A flux f is integral when, for every arc X → Y, 
    the flow f(X,Y) over that arc is an integer. 
    We will typically be dealing with fluxes, 
    pseudoflows, and flows that are integral. 
    If a pseudoflow on some flow network is integral, 
    then every arc is either idle or saturated — 
    no arcs have fractional flow.
Given any flux f in the flow network N_G, 
    we define the cost of that flux to be the sum 
    of the costs of its arcs:
        c(f) = summation over arcs in N_G ( f(X,Y) * c(X,Y) )
Matchings in G are integral flows in N_G
     and number of matchings in M = |f|
     and c(M) = c(f)
     Thus, a min-cost matching of some size s corresponds to 
     a min-cost integral flow of value s.
-------    
    */
    
    /**
     * find augmenting paths of minimal length.
     * by building the shortest path forest.
     * 
     * implementing section 3.4, pg 22 of Ramshaw and Tarjan 2012.
     * 
     * TODO: edit to accept lambda, the largest link length
     * in the forest.
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
    protected DoubleLinkedCircularList[] buildForest(
        final ResidualDigraph rM, int lambda) {
        
        DoubleLinkedCircularList[] forest 
            = new DoubleLinkedCircularList[lambda];
        
        Heap heap = new Heap();
        
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
        paths.
        so forest[0] holds a doubly linked list.
           each item in that doubly linked list is a maiden 
             node without a predecessor and having dist=0.        
        */
        
        // init all nodes to inf length
        Map<Integer, LeftNode> leftNodes = new HashMap<Integer, LeftNode>();
        Map<Integer, RightNode> rightNodes = new HashMap<Integer, RightNode>();
        for (Integer rNode : rM.getRightRM()) {
            RightNode node = new RightNode();
            node.setKey(Long.MAX_VALUE);
            node.setData(rNode);
            rightNodes.put(rNode, node);
        }
        for (Integer lNode : rM.getLeftRM()) {
            LeftNode node = new LeftNode();
            node.setKey(Long.MAX_VALUE);
            node.setData(lNode);
            leftNodes.put(lNode, node);
        }
        
        // married X nodes
        Set<Integer> matchedLeft = new HashSet<Integer>(
            rM.getBackwardLinksRM().values());
        
  //TODO: Fredman and Tarjan [10]      
   
        // for all maidens
        // set key to 0, then ScanAndAdd(index)
        for (Integer lNode : rM.getLeftRM()) {
            if (matchedLeft.contains(lNode)) {
                continue;
            }
            LeftNode node = leftNodes.get(lNode);
            node.setKey(0);
            scanAndAdd(heap, forest, rM, rightNodes, node);
        }
        
        // at this point, the maidens are all in index 0 of
        // the forest trees.
        
        while (!heap.isEmpty()) {
                        
            // in the heap are men not in the forest who are in
            // an alternating path from a maiden.
            // the key is the length of shortest path so far
            HeapNode y = heap.extractMin();
            assert(y instanceof RightNode);
            assert(y.getData() != null);
        
            log.fine("heap.size=" + heap.getNumberOfNodes());
            
            addToForest(forest, y);
                    
            /*
            if y is married then
                x := wife of y;
                set l(x) := l(y) and 
                ScanAndAdd(x); 
             else
                exit(bachelor β := y reached); 
             fi;
            */
            Integer yIndex = (Integer)(y.getData());
            // the married nodes are the keys in the backward
            // links of the residual digraph
            Integer xIndex = rM.getBackwardLinksRM().get(yIndex);
            if (xIndex != null) {
                
                LeftNode xNode = leftNodes.get(xIndex);
                RightNode yNode = rightNodes.get(yIndex);
                
                xNode.setKey(yNode.getKey());
                xNode.pathPredecessor = yNode;
                scanAndAdd(heap, forest, rM, rightNodes, xNode);
                
            } else {
                // break early for finding a complete short augmenting path
                //exit(bachelor β := y reached);
                log.info("bachelor y index=" + y.getData());
                break;
            }
        }
        
        return forest;
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
    private void scanAndAdd(
        Heap heap, DoubleLinkedCircularList[] forest,
        ResidualDigraph rM, Map<Integer, RightNode> yNodes, 
        LeftNode xNode) {
        
        long lX = xNode.getKey();
        Integer x = (Integer)(xNode.getData());
        
        Set<Integer> forwardLinks = rM.getForwardLinksRM().get(x);
        for (Integer y : forwardLinks) {
            
            //link length = net cost of the edge 
            //    (lp(x ⇒ y) = cp(X,Y))
            //link length of backward link Y->X, is 0, lp(y ⇒ x) := 0.
            
            // l(x) and l(y) are the keys in the heap node
            
            RightNode yNode = yNodes.get(y);
            long lY = yNode.getKey();
            assert(((Integer)(yNode.getData())).intValue() 
                == y.intValue());
            
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
            long lOld = lY;
            if (ell < lOld) {
                lY = ell;            
                yNode.pathPredecessor = xNode;
                if (lOld == Long.MAX_VALUE) {
                    yNode.setKey(lY);
                    heap.insert(yNode);
                    log.fine(String.format("HEAP insert: %s",
                        yNode.toString()));
                } else {
                    heap.decreaseKey(yNode, lY);
                    log.fine(String.format("HEAP decr: %s",
                        yNode.toString()));
                }
            }
        }
            
        //add x to the forest;
        addToForest(forest, xNode);
    }
    
    private void addToForest(DoubleLinkedCircularList[] forest,
        HeapNode node) {
        
        /*
        Section 8.1, pg 57
        We ignore any paths we find whose lengths exceed forest.length. 
        We also maintain an integer B, which stores the 
        value l(v) for the node v that was most recently 
        added to the forest. 
        We add nodes v to the forest in nondecreasing 
        order of l(v), so B never decreases. 
        To implement insert(v,k), we add v to the list Q[k].
        */
        
        // NOTE: if the same node is inserted more than once,
        // the forest will be corrupted
        
        if (node.getKey() < forest.length) {
            int k = (int) node.getKey();
            DoubleLinkedCircularList list = forest[k];
            if (list == null) {
                list = new DoubleLinkedCircularList();
                forest[k] = list;
            }
            list.insert(node);
            
            log.fine("add to forest[" + k + "] " 
                + " node " + node.toString());
        }
    }

    private boolean allAreEmpty(DoubleLinkedCircularList[] 
        augmentingPaths) {
        
        for (int i = 0; i < augmentingPaths.length; ++i) {
            DoubleLinkedCircularList path = augmentingPaths[i];
            if (path != null) {
                return false;
            }
        }
        
        return true;
    }
    
    private void swapLinkExistence(ResidualDigraph rM, 
        Integer leftIndex, Integer rightIndex) {
        
        Set<Integer> rIndexes = rM.getForwardLinksRM().get(leftIndex);
        
        boolean forwardFound = (rIndexes != null) && rIndexes.contains(rightIndex);
        
        if (forwardFound) {
            
            // remove existing "idle" forward link
            rIndexes.remove(rightIndex);
            
            // create a backward link and matched mapping
            rM.getBackwardLinksRM().put(rightIndex, leftIndex);
            
            log.fine("augment to add :" + leftIndex + " to " +
                rightIndex);
            
            return;
        }
        
        log.fine("augment to remove :" + leftIndex + " to " +
            rightIndex);
        
        // assert that a backward link exists
        Integer v2 = rM.getBackwardLinksRM().get(rightIndex);
        assert(v2 != null);
        assert(v2.equals(leftIndex));
        
        // remove backwards link and mapping
        rM.getBackwardLinksRM().remove(rightIndex);
        
        // create a forward link        
        if (rIndexes == null) {
            rIndexes = new HashSet<Integer>();
            rM.getForwardLinksRM().put(leftIndex, rIndexes);
        }
        rIndexes.add(rightIndex);
    }
    
    List<PathNode> extractNodes(HeapNode node) {
        
        List<PathNode> nodes = new ArrayList<PathNode>();
        
        PathNode node1 = (PathNode) node;
        while (node1 != null) {
            nodes.add(node1);
            node1 = node1.pathPredecessor;
        }
                
        return nodes;
    }
    
    private void insertIntoHeap(DoubleLinkedCircularList[] minHeap, 
        PathNode node) {
         
        int key = (int)node.getKey();
        
        DoubleLinkedCircularList bucket = minHeap[key];
        if (bucket == null) {
            bucket = new DoubleLinkedCircularList();
            minHeap[key] = bucket;
        }
        
        log.info("insert into heap=" + node);
        
        bucket.insert(node);        
    }
    
    private PathNode extractMinFromHeap(
        DoubleLinkedCircularList[] minHeap) {

        for (DoubleLinkedCircularList bucket : minHeap) {
            if (bucket != null && (bucket.getNumberOfNodes() > 0)) {
                HeapNode node = bucket.getSentinel().getLeft();
                bucket.remove(node);
                return (PathNode)node;
            }
        }
        
        return null;
    }

    private void decreaseKeyInHeap(DoubleLinkedCircularList[] 
        minHeap, PathNode node2, long lTot) {

        log.info("decreaseKey=" + node2);
        
        int prevKey = (int)node2.getKey();
        minHeap[prevKey].remove(node2);
        
        node2.setKey(lTot);
        
        insertIntoHeap(minHeap, node2);
    }

    private void debug(DoubleLinkedCircularList[] forest) {
        
        for (int i = 0; i < forest.length; ++i) {

            DoubleLinkedCircularList tree = forest[i];
            if (tree == null) {
                continue;
            }
            long n = tree.getNumberOfNodes();
            HeapNode node = tree.getSentinel();
            int j = 0;
            while (j < n) {
                node = node.getRight();                    
                List<PathNode> path = extractNodes(node);
                for (int ii = 0; ii < (path.size() - 1); ++ii) {
                    PathNode node1 = path.get(ii);
                    PathNode node2 = path.get(ii + 1);
                    log.info("forest2[" + i + "] tree branch[" 
                    + j + "] node[" + ii + "]=" + node1.toString());
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
                    /*   
                    PairInt pair = new PairInt(index1.intValue(),
                        index2.intValue());
                    if (augmented.contains(pair)) {
                        break;
                    }
                    swapLinkExistence(rM, index1, index2);
                    augmented.add(pair);
                    */
                }                    
                j++;
            }
        }            
    }
    
    private void handleSourceForwardLink(
        DoubleLinkedCircularList[] minHeap,
        FlowNetwork gFlow, SourceNode node1,
        Integer index2, Map<Integer, LeftNode> leftNodes,
        int lambda, float eps) {
    
        long l1 = node1.getKey();
        int idx1 = ((Integer)node1.getData()).intValue();
        int idx2 = index2.intValue();
        
        LeftNode node2 = leftNodes.get(index2);
        float cp = gFlow.calcSourceNetCost(idx2);
        
        long lp = (long) Math.ceil(cp / eps);
        long lTot = l1 + lp;
        long lOld = node2.getKey();
        if ((lTot < lambda) && (lTot < lOld)) {
            node2.pathPredecessor = node1;
            if (lOld == Long.MAX_VALUE) {
                node2.setKey(lTot);
                insertIntoHeap(minHeap, node2);
            } else {
                decreaseKeyInHeap(minHeap, node2, lTot);
            }
        }
    }
    
    private void handleLeft(DoubleLinkedCircularList[] minHeap,
        FlowNetwork gFlow, PathNode node1,
        Integer index2, Map<Integer, LeftNode> leftNodes,
        int lambda, float eps) {
    
        long l1 = node1.getKey();
        int idx1 = ((Integer)node1.getData()).intValue();
        int idx2 = index2.intValue();
        
        LeftNode node2 = leftNodes.get(index2);
        float cp = gFlow.calcNetCost(index2.intValue(), idx1);
        long lp = 1 - (long) Math.ceil(cp / eps);
        long lTot = l1 + lp;
        long lOld = node2.getKey();
        if ((lTot < lambda) && (lTot < lOld)) {
            node2.pathPredecessor = node1;
            if (lOld == Long.MAX_VALUE) {
                node2.setKey(lTot);
                insertIntoHeap(minHeap, node2);
            } else {
                decreaseKeyInHeap(minHeap, node2, lTot);
            }
        }
    }
    
    private void handleRight(DoubleLinkedCircularList[] minHeap,
        FlowNetwork gFlow, PathNode node1,
        Integer index2, Map<Integer, RightNode> rightNodes,
        int lambda, float eps) {
        
        long l1 = node1.getKey();
        int idx1 = ((Integer)node1.getData()).intValue();
        int idx2 = index2.intValue();
        
        RightNode node2 = rightNodes.get(index2);
        float cp = gFlow.calcNetCost(idx1, idx2);
        long lp = (long) Math.ceil(cp / eps);
        long lTot = l1 + lp;
        long lOld = node2.getKey();     
        if ((lTot < lambda) && (lTot < lOld)) {
            node2.pathPredecessor = node1;
            if (lOld == Long.MAX_VALUE) {
                node2.setKey(lTot);
                insertIntoHeap(minHeap, node2);
            } else {
                decreaseKeyInHeap(minHeap, node2, lTot);
            }
        }
    }
    
    private void handleSinkBackwardLink(
        DoubleLinkedCircularList[] minHeap,
        FlowNetwork gFlow, SinkNode node1,
        Integer index2, Map<Integer, RightNode> rightNodes,
        int lambda, float eps) {
        
        long l1 = node1.getKey();
        int idx1 = ((Integer)node1.getData()).intValue();
        int idx2 = index2.intValue();
        
        RightNode node2 = rightNodes.get(index2);
        float cp = gFlow.calcSinkNetCost(idx2);
        long lp = 1 - (long) Math.ceil(cp / eps);
        long lTot = l1 + lp;
        long lOld = node2.getKey();
        if ((lTot < lambda) && (lTot < lOld)) {
            node2.pathPredecessor = node1;
            if (lOld == Long.MAX_VALUE) {
                node2.setKey(lTot);
                insertIntoHeap(minHeap, node2);
            } else {
                decreaseKeyInHeap(minHeap, node2, lTot);
            }
        }
    }
    
    private void handleSourceBackwardLink(
        DoubleLinkedCircularList[] minHeap,
        FlowNetwork gFlow, PathNode node1,
        SourceNode node2, int lambda, float eps) {
        
        long l1 = node1.getKey();
        int idx1 = ((Integer)node1.getData()).intValue();
        int idx2 = ((Integer)node2.getData()).intValue();

        float cp = gFlow.calcSourceNetCost(idx1);
        long lp = 1 - (long) Math.ceil(cp / eps);
        
        long lTot = l1 + lp;
        long lOld = node2.getKey();
        if ((lTot < lambda) && (lTot < lOld)) {
            node2.pathPredecessor = node1;
            if (lOld == Long.MAX_VALUE) {
                node2.setKey(lTot);
                insertIntoHeap(minHeap, node2);
            } else {
                decreaseKeyInHeap(minHeap, node2, lTot);
            }
        }
    }
    
    private void handleSink(DoubleLinkedCircularList[] minHeap,
        FlowNetwork gFlow, PathNode node1,
        SinkNode node2, int lambda, float eps) {
        
        long l1 = node1.getKey();
        int idx1 = ((Integer)node1.getData()).intValue();
        int idx2 = ((Integer)node2.getData()).intValue();

        float cp = gFlow.calcSinkNetCost(idx1);
        long lp = (long) Math.ceil(cp / eps);
        
        long lTot = l1 + lp;
        long lOld = node2.getKey();
        if ((lTot < lambda) && (lTot < lOld)) {
            node2.pathPredecessor = node1;
            if (lOld == Long.MAX_VALUE) {
                node2.setKey(lTot);
                insertIntoHeap(minHeap, node2);
            } else {
                decreaseKeyInHeap(minHeap, node2, lTot);
            }
        }
    }
    
    private void insert(Map<Integer, Map<Integer, List<PathNode>>> 
        forestTreeMap, int forestIdx, int treeIdx, 
        List<PathNode> path) {
        
        Integer index0 = Integer.valueOf(forestIdx);
        Integer index1 = Integer.valueOf(treeIdx);
        
        Map<Integer, List<PathNode>> map0 = 
            forestTreeMap.get(index0);
        if (map0 == null) {
            map0 = new HashMap<Integer, List<PathNode>>();
            forestTreeMap.put(index0, map0);
        }
        
        map0.put(index1, path);
        
    }
    
    private void insert(Map<Integer, Set<TrioInt>> map, 
        PathNode node1, PathNode node2, 
        int forestIdx, int treeIdx, int branchIdx) {
        
        TrioInt trioInt = new TrioInt(forestIdx, treeIdx, branchIdx);
        
        Integer key = (Integer)node1.getData();
        
        Set<TrioInt> set = map.get(key);
        if (set == null) {
            set = new HashSet<TrioInt>();
            map.put(key, set);
        }
        
        set.add(trioInt);
    }

}
