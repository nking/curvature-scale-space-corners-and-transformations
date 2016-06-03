package algorithms.bipartite;

import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.Heap;
import algorithms.imageProcessing.HeapNode;
import algorithms.misc.Misc;
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
import java.util.Stack;
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
 * 
 Note that the min-cost is the optimized goal, so if a 
 graph is given that has a possible solution with more 
 matches than the number of min-cost matches, this 
 * solver will not completely match the nodes in that
 * graph.
 * (in that case, one might want to create a pre-processing
 * stage to filter the graph...pre-matching those edges
 * that must be matched inspite of cost because no other
 * edges are connected to those nodes for example.
 * might consider the very fast state management
 * present in the conflict analysis of the fastest boolean
 * sat solvers to design a fast pre-filter).
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
    public static class LeftNode extends PathNode {
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
        // m = 0 is unmarked, m=1 is marked
        int m = 0;
        public abstract PathNode copy();
        // id is only used in toString fr debug statements
        String id = "";
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append(id).append(" key=").append(Long.toString(getKey()));
            sb.append(" index=").append(getData().toString());
            PathNode prev = pathPredecessor;
            while (prev != null) {
                sb.append("\n    [prev=").append(prev.toString())
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
    public static class RightNode extends PathNode {
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
    public static class SourceNode extends PathNode {
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
    public static class SinkNode extends PathNode {
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
        
        if (g.getSourceNode() == -1) {
            g = g.copyToCreateSourceSink();
        }
        
        {//DEBUG
            log.info("hopcroft-karp matches:");
            for (Entry<Integer, Integer> entry : m.entrySet()) {
                log.info("hk matched " + entry.getKey() + " to "
                    + entry.getValue());
            }
        }
        
        FlowNetwork gFlow = new FlowNetwork(g, m);
        //assert(gFlow.printFlowValueIncludingSrcSnk(m.size()));
        
        log.info("init FlowNetwork:");
        gFlow.printNetCosts();
        
        //TODO: estimate of eps may need to be revised
        
        int s = m.size();

        // q >= 2.  consider q = O(log_2(n)) pg 41 par 3.
        int q = 2;

        // pg 32, the weight scaling techinique
        // starts w/ eps ~ maxC and reduces to ~1/s or 1/(6*s)
        // w/ nIter ~ log_q(s*maxC)
        // pg 44
        // eps_up = q^(e_up) where eps_up is smallest power of
        //    q that exceeds maxC
        // e_up * math.log(q)
        int e_up = 1 + (int)Math.floor(Math.log(gFlow.getMaxC())/Math.log(q));
        double eps_up = Math.pow(q, e_up);
        
        int e_down = -(1 + (int)Math.floor(
            Math.log(s + 2)/Math.log(q)));
        
        double eps_down = Math.pow(q, e_down);
        
        int eps = 1 + (int)Math.floor(eps_up);
        if (eps > q) {
            // first iteration through net costs math.ceil(cp/eps)
            // will be able to distinguish between integer differences
            // in cost of 1 if the eps it receives is q-1 = 1.
            // since eps is divided by q before refine,
            // eps should be approx q here.
            // TODO: needs testing for a range of max(cost) w.r.t. q
            eps = q;
        }
        
        // expected number of iterations without a constant factor
        int rIter = (int)(Math.log(s * gFlow.getMaxC())/Math.log(q));
                  
        int nIterR = 0;
        
        log.info("eps=" + eps + " rIter=" + rIter + " maxC=" + 
            gFlow.getMaxC() + " eps_down=" + eps_down);
               
        // all nodes V in gFlow have prices = 0
        
        while (eps > eps_down) {
            
            log.info("nIterR=" + nIterR + " s=" + s);
            
            // pg 44, assertions I1, I2, I3, and I4
            assert(gFlow.assertFlowValue(s));
            assert(gFlow.assertPricesAreQuantizedEps(eps));
            // if using a smaller eps than maxC, cannot assert
            //    these on first round
            if ((nIterR > 0) || (eps > gFlow.getMaxC())) {
                assert(gFlow.integralFlowIsEpsProper(eps));
                assert(gFlow.assertSaturatedBipartiteIsEpsSnug(eps));
            }
            
            eps /= (float)q;
            
            int ext = refine(gFlow, s, eps, q);
            
            if (ext > 0) {
                m = gFlow.extractMatches();
                // if returning gFlow instead:
                roundFinalPrices(gFlow, eps_down);
                gFlow.printFlowValueIncludingSrcSnk(m.size());
                return m;
            }
            
            ++nIterR;
        }
        
        // assert nIter is approx log_q(s * maxC)
        
        // round prices to integers that make all arcs proper
        roundFinalPrices(gFlow, eps_down);
        
        m = gFlow.extractMatches();
        
        return m;
    }
    
    protected int refine(FlowNetwork gFlow, int s, float eps,
        int q) {
        
        log.info("at start of refine, s=" + s);
        
        assert(gFlow.printFlowValueIncludingSrcSnk(s));
        
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
            log.info("surplus idx=" + Integer.toString(idx1));
        }

        assert(surplus.size() == s);
        
        /*
        see Figure 7.4 on pg 53.
        raise prices so that every arc in 
        gFlow becomes eps-proper, for the resulting pseudoflow 
        f and for the new, smaller value of eps.
        */
        gFlow.raisePricesUntilEpsProper(eps, q);
        
        log.info("after raise prices:");
        gFlow.printNetCosts();
        assert(gFlow.printFlowValueIncludingSrcSnk(s));        
        
        //in [0] holds the length of the terminating deficit
        // node which is also the forest index it was
        // added to.
        int[] terminatingDeficitIdx = new int[1];
        int h = s;
        int nHIter = 0;
        while (h > 0) {
            
            log.info("nHIter=" + nHIter + " h=" + h);
            
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
            
            assert(gFlow.printFlowValueIncludingSrcSnk(s));
            
            modifyPricesAndPathLengths(gFlow, forest, 
                terminatingDeficitIdx[0], eps);
            
            log.info("after modify prices and link lengths");
            debug(forest);
            
            //assert from pg 52 
            //       I1', I2, I3, I4, on FlowNetwork
            //       and I5 on ResidualDigraph2
            
            assert(gFlow.assertFlowValueIncludingSrcSnk(s));
            assert(gFlow.assertPricesAreQuantizedEps(eps));
            assert(gFlow.integralFlowIsEpsProper(eps));
            assert(gFlow.assertSaturatedBipartiteIsEpsSnug(eps));
            // assert I5: Rf has no cycles of length zero.
            
            log.info("srching for zero length paths in this:");
            debug(forest);
            
            //Let Rf0 denote that subgraph of the residual digraph 
            //formed by links of length zero.
            List<LinkedList<PathNode>> zeroLengthLinks =
                extractZeroLengthLinks(forest);
            
            assert(gFlow.printSurplusAndDeficit());
            
            // --- Sect 8.3, create maximal set of compatible augmenting paths
            List<LinkedList<PathNode>> cPaths =
                findMaximalSetOfCompatiblePaths(
                gFlow.getNLeft(), gFlow.getNRight(),
                gFlow.getSinkNode(), gFlow.getSourceNode(),
                zeroLengthLinks,surplus, 
                new HashSet<Integer>(deficit));
            
            if (cPaths.isEmpty()) {
                log.severe("did not find an augmenting path.  h=" + h);
                
                //gFlow.printSaturatedLinks();
                
                /*
                NOTE: when this has number of saturated arcs
                less than s, where s is the number that we know
                is the largest possible matching from hopcroft-karp
                which only cares about connectedness to maximize
                the number of matches, there were higher cost 
                edges which had smaller connectivity so dropped
                out of the solution.
                because the algorithm goal is min-cost, can
                just return this smaller set of matchings than s.
                */
                
                return 1;
            }
            
            debug(cPaths);
            
            //augment f along each of the paths in P in turn, thereby 
            //   reducing |S| = |D| = h by |P|;
            augmentFlow(gFlow, cPaths, surplus, deficit);
            
            log.info("after augment flow:");
            gFlow.printSaturatedLinks();
            
            h = surplus.size();
            
            // pg 63 assert I1', I2, I34  (not I3, but I4?)
            assert(gFlow.assertFlowValueIncludingSrcSnk(s));
            assert(gFlow.assertPricesAreQuantizedEps(eps));
            assert(gFlow.integralFlowIsEpsProper(eps));
            
            ++nHIter;
        }
        
        return 0;
    }
    
    /**
     * 
     * @param gFlow
     * @param rF
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
        List<Integer> surplus, List<Integer> deficit, float eps,
        int[] terminatingDeficitIdx) {
        
        log.info("buildForest2");
        
        Set<Integer> d = new HashSet<Integer>(deficit);
    
        //TODO: revisit this. 
        int lambda = (int)Math.ceil(3 * (gFlow.getMaxC()/eps));
        if (lambda < 4) {
            lambda = 4;
        }
        
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
            if (node1 == null) {
                break;
            }
            Integer index1 = (Integer)node1.getData();
            int idx1 = index1.intValue();
         
            log.info("extractMin = " + node1.toString());
            
            boolean node1IsLeft = (node1 instanceof LeftNode);
            boolean nodeIsSource = (node1 instanceof SourceNode);
            boolean nodeIsSink = (node1 instanceof SinkNode);
            
            //scan all links leaving node1:
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
     * and right nodes in the bipartite graph g.
     * Note that if the maximum matched size does not reach
     * the requested size s, but has not grown during internal
     * iterations, the result is returned for that smaller size.
     * @param g
     * @param s
     * @return 
     */
    protected Map<Integer, Integer> hopcroftKarp(Graph g,
        int s) {
        
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
        
        if (false) { 
            //runtime complexity is O(m * sqrt(max matching size))
            //    where m is number of edges in graph
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
        
        // runtime complexity O(m * sqrt(s)) where m is number of edges
        // and s is the size of the matching whose target size
        // may be less than the maximum matchable
        ResidualDigraph rM = new ResidualDigraph(g, m);
                
        return hopcroftKarp(g, rM, s);
    }
    
    /**
     * NOT READY FOR USE
     * NOTE: this method is not yet throughly tested.
     * Note that if s is an overestimate and the algorithm
     * finds a repeated size upon internal iteration, it
     * will return that smaller size match.
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

        //NOTE: this is not correct.  does not always produce
        //  maximal matchings
        
        while (true) {
            
            DoubleLinkedCircularList[] augmentingPaths =
                buildForest(rM, lambda);

            if (allAreEmpty(augmentingPaths)) {
                return m;
            }
            
            //debug(augmentingPaths);
                       
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
            log.info("nIter=" + nIter + " m2.size=" + m2.size()
                + " m.size=" + m.size());

            // debug:
            for (Entry<Integer, Integer> entry : m2.entrySet()) {
                log.info("m2 match= " + entry.getKey() + "->" + entry.getValue());
            }
            for (Entry<Integer, Integer> entry : m.entrySet()) {
                log.info("m match= " + entry.getKey() + "->" + entry.getValue());
            }
    
            if (m2.size() >= s) {
                return m2;
            }
            
            Map<Integer, Integer> tmpM = new HashMap<Integer, Integer>(m);
                    
            Set<Integer> mR = new HashSet<Integer>(m.values());
            // if m2 has a match that does not conflict with m,
            // add it to m (== vertex disjoint)
            for (Entry<Integer, Integer> entry : m2.entrySet()) {
                Integer key = entry.getKey();
                Integer value = entry.getValue();
                if (!tmpM.containsKey(key) && !mR.contains(value)) {
                    tmpM.put(key, value);
                }
            }
            
            // remove the intersection of m and m2
            for (Entry<Integer, Integer> entry : m.entrySet()) {
                Integer key = entry.getKey();
                if (m2.containsKey(key) &&
                    m2.get(key).equals(entry.getValue())) {
                    tmpM.remove(key);
                }
            }
            
            m = tmpM;
            log.info("symmetric diff of m and m2.size=" + m.size());
            for (Entry<Integer, Integer> entry : m.entrySet()) {
                log.info("sym diff m match= " + entry.getKey() + "->" + entry.getValue());
            }
            
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
    
    private void roundFinalPrices(FlowNetwork gFlow,
        double eps_down) {
        
        //TODO: this may need revision
        
        //see pg 46
        //round prices to integers with
        // pd^~(v) = math.floor(pd(v) + k * eps_down)
        // where k is a carefully chosen integer in
        // the range (0 to 1)/eps_down
        // 
        // k is chosen from 0 to 16, excluding k=10
        //
        // each currently "improper" arc for saturated arcs
        //   is used to examine
        //   k as fractional part of pd(v) or pd(w).
        //   the bad value :
        //      the largest k that causes pd(v) to round down
        //          == same k that causes pd(w) to round up
        // the arcs that are scanned are chosen to be the source
        //   to left arcs.
        //   the scan pattern is to calculate whether each
        //      saturated arc is "proper" or "improper" for
        //      all values of k, that is integers in range
        //      [0 to 1)/eps_down are marked first as good
        //      and then as bad when a bad value of k is found.
        //
        //      then from the integers that are still marked good,
        //      choose that as a good value for k
        //      and adjust the prices using it.
        //
        //      the final resulting f will then be min-cost
        
        // definition of saturated and improper in
        // Figure 6.1, p 43
        //   improper is net cost > 0
        //   eps-proper is net cost > -eps
        
        int sourceNode = gFlow.getSourceNode();
        
        Set<Integer> ks = new HashSet<Integer>();
        int maxK = (int)Math.ceil(1./eps_down);
        for (int k = 0; k < maxK; ++k) {
            
            boolean keep = true;
        
            float factor = (float)eps_down * k;
            
            for (Integer xIndex : gFlow.getSourceForwardArcs()) {
                int xIdx = xIndex.intValue();
                float unitFlow = gFlow.getSourceToLeftFlow(xIdx);
                if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated
                    //cp = cost - pdX + pdY so incr pdX                    
                    int cost = gFlow.getSourceToLeftCost(xIdx);
                    double pdX = gFlow.getLeftPrice(sourceNode);
                    double pdY = gFlow.getLeftPrice(xIdx);
                    
                    // values causing pdX to round down
                    // while pdY rounds up
                    // are not proper
                    double pdX2 = Math.floor(pdX + factor);
                    double pdY2 = Math.floor(pdY + factor);
                        
                    if ((pdX2 < pdX) && (pdY2 > pdY)) {
                        keep = false;
                        break;
                    }
                    
                }
            }
            if (keep) {
                ks.add(Integer.valueOf(k));
            }
        }
        
        if (ks.size() == maxK) {
            //TODO: revisit this
            // no improper arcs
            return;
        }
        
        assert(!ks.isEmpty());
        
        // pick from ks to modify prices
        //pd^~(v) = math.floor(pd(v) + k * eps_down)
        
        int k = ks.iterator().next().intValue();
        gFlow.addToAllPrices((float)(k * eps_down));        
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
     * Note that have edited it to not exit at first
     * unmatched y, to let it continue to fill all augmentable
     * paths.  the runtime complexity needs to be
     * adjusted for this.
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
        
        /*
        NOTE: bookkeeping for the predecessors could be improved
        here.
        */
        
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
        log.info("done adding maiden nodes to heap");
        
        // at this point, the maidens are all in index 0 of
        // the forest trees.
       
        int nIter = 0;
                
        while (!heap.isEmpty()) {
                        
            // in the heap are men not in the forest who are in
            // an alternating path from a maiden.
            // the key is the length of shortest path so far
            PathNode y = (PathNode)heap.extractMin();
            assert(y instanceof RightNode);
            assert(y.getData() != null);
        
            log.fine("heap.size=" + heap.getNumberOfNodes());
            
            log.info("extractMin=" + y.toString());
            
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
                
                xNode.setKey(y.getKey());
                if (xNode.pathPredecessor == null) {
                    xNode.pathPredecessor = y;
                }
                                
                scanAndAdd(heap, forest, rM, rightNodes, xNode);
                                
            } else {
                //exit(bachelor β := y reached);
                log.info("bachelor y=" + y.toString());
                if (y.getKey() > 1) {
                    break;
                }
            }
            
            log.info("nIter=" + nIter);
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
    private long scanAndAdd(
        Heap heap, DoubleLinkedCircularList[] forest,
        ResidualDigraph rM, Map<Integer, RightNode> yNodes, 
        LeftNode xNode) {
        
        long lX = xNode.getKey();
        assert(lX < Long.MAX_VALUE);
        Integer x = (Integer)(xNode.getData());
        
        Set<Integer> forwardLinks = rM.getForwardLinksRM().get(x);
        for (Integer y : forwardLinks) {
            
            //link length = net cost of the edge 
            //    (lp(x ⇒ y) = cp(X,Y))
            //link length of backward link Y->X, is 0, lp(y ⇒ x) := 0.
            
            // l(x) and l(y) are the keys in the heap node
            
            RightNode yNode = yNodes.get(y);
            long lOld = yNode.getKey();
            assert(((Integer)(yNode.getData())).intValue() == y.intValue());
            
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
                yNode = yNode.copy();
                yNode.pathPredecessor = xNode.copy();
                if (lOld == Long.MAX_VALUE) {
                    yNode.setKey(ell);
                    heap.insert(yNode);
                    log.info(String.format("HEAP insert: %s",
                        yNode.toString()));
                } else {
                    heap.decreaseKey(yNode, ell);
                    log.info(String.format("HEAP decr: %s",
                        yNode.toString()));
                }
            }
        }
        
        //add x to the forest;
        addToForest(forest, xNode);
   
        return lX;
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
    
    public static List<PathNode> extractNodes(HeapNode node) {
        
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

        log.info("decreaseKey=" + node2 + " to key=" + lTot);
        
        int prevKey = (int)node2.getKey();
        minHeap[prevKey].remove(node2);
        
        node2.setKey(lTot);
        
        insertIntoHeap(minHeap, node2);
    }
    
    private void debug(List<LinkedList<PathNode>> cPaths) {

        log.info("cPaths.size=" + cPaths.size());
        
        for (int i = 0; i < cPaths.size(); ++i) {
            StringBuilder sb = new StringBuilder("i=");
            sb.append(Integer.toString(i)).append(" ");
            LinkedList<PathNode> path = cPaths.get(i);
            for (PathNode node : path) {
                sb.append(node.toString()).append(" ");
            }
            log.info(sb.toString());
        }
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
                int n2 = path.size();
                for (int ii = 0; ii < n2; ++ii) {
                    PathNode node1 = path.get(ii);
                    log.info("forest2[" + i + "] tree branch[" 
                        + j + "] node[" + ii + "]=" + node1.toString());
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
        log.info(String.format("(source,%d) cp=%.1f lp=%s lTot=%s lOld=%s",
            idx2, cp, Long.toString(lp), Long.toString(lTot), 
            Long.toString(lOld)));
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
        log.info(String.format("(%d,%d) cp=%.1f lp=%s lTot=%s lOld=%s",
            idx2, idx1, cp, Long.toString(lp), Long.toString(lTot), 
            Long.toString(lOld)));
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
        log.info(String.format("(%d,%d) cp=%.1f lp=%s lTot=%s lOld=%s",
            idx1, idx2, cp, Long.toString(lp), Long.toString(lTot), 
            Long.toString(lOld)));
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
        log.info(String.format("(%d,sink) cp=%.1f lp=%s lTot=%s lOld=%s",
            idx2, cp, Long.toString(lp), Long.toString(lTot), 
            Long.toString(lOld)));
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
        log.info(String.format("(source,%d) cp=%.1f lp=%s lTot=%s lOld=%s",
            idx1, cp, Long.toString(lp), Long.toString(lTot), 
            Long.toString(lOld)));
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
        log.info(String.format("(%d,sink) cp=%.1f lp=%s lTot=%s lOld=%s",
            idx1, cp, Long.toString(lp), Long.toString(lTot), 
            Long.toString(lOld)));
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
    
    private void pruneToLength0Links(
        Map<Integer, Map<Integer, List<PathNode>>> forestTreeMap, 
        Map<Integer, Set<TrioInt>> leftToRightLoc, 
        Map<Integer, Set<TrioInt>> rightToLeftLoc, 
        Map<Integer, TrioInt> sourceToLeftLoc, 
        Map<Integer, TrioInt> leftToSourceLoc, 
        Map<Integer, TrioInt> rightToSinkLoc, 
        Map<Integer, TrioInt> sinkToRightLoc) {

        Map<Integer, Set<Integer>> leftToRight = new
            HashMap<Integer, Set<Integer>>();
        for (Entry<Integer, Set<TrioInt>> entry :
            leftToRightLoc.entrySet()) {
            Integer leftIndex = entry.getKey();
            for (TrioInt loc : entry.getValue()) {
                PathNode node = forestTreeMap
                    .get(Integer.valueOf(loc.getX()))
                    .get(Integer.valueOf(loc.getY()))
                    .get(loc.getZ());
                if (node.getKey() == 0L) {
                    Set<Integer> indexes2 = leftToRight.get(leftIndex);
                    if (indexes2 == null) {
                        indexes2 = new HashSet<Integer>();
                        leftToRight.put(leftIndex, indexes2);
                    }
                    indexes2.add((Integer)node.getData());
                }
            }
        }
    }
    
    private void modifyPricesAndPathLengths(FlowNetwork gFlow,
        DoubleLinkedCircularList[] forest, int lt,
        float eps) {

        Map<Integer, Integer> incrLeft = 
            new HashMap<Integer, Integer>();
        Map<Integer, Integer> incrRight = 
            new HashMap<Integer, Integer>();
        
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

        Map<Integer, Map<Integer, List<PathNode>>> forestTreeMap
            = new HashMap<Integer, Map<Integer, List<PathNode>>>();

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
                Misc.<PathNode>reverse(path);
                insert(forestTreeMap, forestIdx, treeIdx, path);
                for (int branchIdx = 0; branchIdx < (path.size() - 1); ++branchIdx) {
                    PathNode node1 = path.get(branchIdx);
                    PathNode node2 = path.get(branchIdx + 1);
                    int l1 = (int) node1.getKey();
                    int l2 = (int) node2.getKey();
                    int idx1 = ((Integer) node1.getData()).intValue();
                    int idx2 = ((Integer) node2.getData()).intValue();
                    if (node1 instanceof LeftNode) {
                        assert (!(node2 instanceof LeftNode));
                        //i(v) := l(termDefIdx) - l(v)
                        if ((lt - l1) != 0) {
                            incrLeft.put(Integer.valueOf(idx1),
                                Integer.valueOf(lt - l1));
                        }
                        if (node2 instanceof RightNode) {
                            insert(leftToRightLoc, node1, node2,
                                forestIdx, treeIdx, branchIdx + 1);
                            if ((lt - l2) != 0) {
                                incrRight.put(Integer.valueOf(idx2),
                                    Integer.valueOf(lt - l2));
                            }
                        } else {
                            assert (node2 instanceof SourceNode);
                            assert (leftToSourceLoc.get(
                                (Integer) node1.getData()) == null);
                            leftToSourceLoc.put(
                                (Integer) node1.getData(),
                                new TrioInt(forestIdx, treeIdx, 
                                    branchIdx + 1));
                            if ((lt - l2) != 0) {
                                incrLeft.put(Integer.valueOf(idx2),
                                    Integer.valueOf(lt - l2));
                            }
                        }
                    } else if (node1 instanceof RightNode) {
                        assert (!(node2 instanceof RightNode));
                        //i(v) := l(termDefIdx) - l(v)
                        if ((lt - l1) != 0) {
                            incrRight.put(Integer.valueOf(idx1),
                                Integer.valueOf(lt - l1));
                        }
                        if (node2 instanceof LeftNode) {
                            insert(rightToLeftLoc, node1, node2,
                                forestIdx, treeIdx, branchIdx + 1);
                            if ((lt - l2) != 0) {
                                incrLeft.put(Integer.valueOf(idx2),
                                    Integer.valueOf(lt - l2));
                            }
                        } else {
                            assert (node2 instanceof SinkNode);
                            assert (rightToSinkLoc.get(
                                (Integer) node1.getData()) == null);
                            rightToSinkLoc.put(
                                (Integer) node1.getData(),
                                new TrioInt(forestIdx, treeIdx, 
                                    branchIdx + 1));
                            if ((lt - l2) != 0) {
                                incrRight.put(Integer.valueOf(idx2),
                                    Integer.valueOf(lt - l2));
                            }
                        }
                    } else if (node1 instanceof SourceNode) {
                        assert (node2 instanceof LeftNode);
                        assert (sourceToLeftLoc.get(
                            (Integer) node2.getData()) == null);
                        sourceToLeftLoc.put(
                            (Integer) node2.getData(),
                            new TrioInt(forestIdx, treeIdx, 
                                branchIdx + 1));
                        //i(v) := l(termDefIdx) - l(v)
                        if ((lt - l1) != 0) {
                            incrLeft.put(Integer.valueOf(idx1),
                                Integer.valueOf(lt - l1));
                        }
                        if ((lt - l2) != 0) {
                            incrLeft.put(Integer.valueOf(idx2),
                                Integer.valueOf(lt - l2));
                        }
                    } else {
                        //node1 is a sink node
                        assert (node2 instanceof RightNode);
                        assert (sinkToRightLoc.get(
                            (Integer) node2.getData()) == null);
                        sinkToRightLoc.put(
                            (Integer) node2.getData(),
                            new TrioInt(forestIdx, treeIdx, 
                                branchIdx + 1));
                        //i(v) := l(termDefIdx) - l(v)
                        if ((lt - l1) != 0) {
                            incrRight.put(Integer.valueOf(idx1),
                                Integer.valueOf(lt - l1));
                        }
                        if ((lt - l2) != 0) {
                            incrRight.put(Integer.valueOf(idx2),
                                Integer.valueOf(lt - l2));
                        }
                    }
                }
                treeIdx++;
            }
        }

        //pd'(v) = pd(v) + i(v)*eps
        for (Entry<Integer, Integer> entry : incrLeft.entrySet()) {
            Integer index = entry.getKey();
            float incr = eps * entry.getValue();
            gFlow.addToLeftPrice(index.intValue(), incr);
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
        
        for (Entry<Integer, Integer> entry : incrRight.entrySet()) {
            Integer index = entry.getKey();
            float incr = eps * entry.getValue().intValue();
            gFlow.addToRightPrice(index.intValue(), incr);
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
    
    private List<LinkedList<PathNode>>
        findMaximalSetOfCompatiblePaths(
            int nLeft, int nRight, int sourceNodeIdx,
            int sinkNodeIdx,
            List<LinkedList<PathNode>> pathLinkLists,
            List<Integer> surplus, Set<Integer> deficit) {
        
        //pseudocode from pg 62 , Figure 8.2
        //   with input = length 0 adj list, the surplus
        //   list and the deficit list
        //   to create the maximal set of compatible paths.
        //
        //   state that needs to be tracked for a vertex:
        //    - visited (== marked)
        //    - identity (== Left or Right or Source or Sink)

        Map<PathNode, Set<PathNode>> pathLinksMap =
            createNewAdjacencyMap(pathLinkLists,
            nLeft, nRight, sourceNodeIdx, sinkNodeIdx);
        
        Set<LeftNode> surp = new HashSet<LeftNode>();
        Set<RightNode> def = new HashSet<RightNode>();
        makeSurplusAndDeficitSubSets(pathLinksMap, 
            surplus, deficit, surp, def);
        
        List<LinkedList<PathNode>> augPaths 
            = new ArrayList<LinkedList<PathNode>>();
        
        Stack<PathNode> stack = new Stack<PathNode>();
        
        A: while (!surp.isEmpty()) {
            LeftNode xNode = surp.iterator().next();
            surp.remove(xNode);
            stack.add(xNode);
            B: while (!stack.isEmpty()) {
                PathNode v = stack.peek();
                log.info("v=stack.top=" + v.toString());
                boolean isLeftNode = (v instanceof LeftNode);
                boolean isRightNode = (v instanceof RightNode);
                if (isLeftNode || isRightNode) {
                    v.m = 1;
                }
                if (isRightNode && def.contains(v)) {
                    // add all of stack to P as an augmenting Path
                    LinkedList<PathNode> augPath = new
                        LinkedList<PathNode>();
                    while (!stack.isEmpty()) {
                        augPath.add(stack.pop());
                    }
                    if (!augPath.isEmpty()) {
                        augPaths.add(augPath);
                    }
                    continue A;
                }
                Set<PathNode> set = pathLinksMap.get(v);
                if (set != null) {
                    while (!set.isEmpty()) {
                        PathNode w = set.iterator().next();
                        set.remove(w);
                        log.info("stack.push w = " + w.toString());
                        if (w.m == 0) {
                            stack.push(w);
                            continue B;
                        }
                    }
                }
                // do nothing with v
                stack.pop();
            }
        }
        
        return augPaths;
        
        /*
        a stack K is needed.
        a (linked list) List L[v] is needed for each node v to 
            hold the w's for links with v->w length = 0.
            (==the rF0 adjacency lists)

        for v in NG do set v to be unmarked;
        set stack K to empty;
        Surp := S;
            A: while Surp not empty do
            x := first(Surp); Surp := rest(Surp);
            push(x, K);
            B: while K not empty do
                v := top(K);
                if v != source and v != sink then mark(v) fi;
                if v in D then
                    add K to P as an augmenting path;
                    set K to empty;
                    goto A;
                fi;
                C: while L[v] not empty do
                    w := first(L[v]); L[v] := rest(L[v]);
                    if w not marked then
                        push(w,K); goto B fi; od;
                    pop(K);
                od;
            od;
        */
    }
    
    private List<LinkedList<PathNode>> extractZeroLengthLinks(
        DoubleLinkedCircularList[] forest) {

        // each linked list is a pair of 2 path nodes whose link length = 0
        List<LinkedList<PathNode>> pathLists = new
            ArrayList<LinkedList<PathNode>>();
        
        // traverse trees in the forest
        for (int forestIdx = 1; forestIdx < forest.length; ++forestIdx) {
            DoubleLinkedCircularList tree = forest[forestIdx];
            if (tree == null) {
                continue;
            }
            long n = tree.getNumberOfNodes();
            HeapNode node = tree.getSentinel();
            int treeIdx = 0;
            while (treeIdx < n) {
                node = node.getRight();
                List<PathNode> path = MinCostUnbalancedAssignment.extractNodes(node);
                int n2 = path.size();
                Misc.<PathNode>reverse(path);
                for (int branchIdx = 0; branchIdx < (n2 - 1); ++branchIdx) {
                    PathNode node1 = path.get(branchIdx);
                    PathNode node2 = path.get(branchIdx + 1);
                    int l1 = (int) node1.getKey();
                    int l2 = (int) node2.getKey();
                    if (l2 != 0) {
                        continue;
                    }
                    LinkedList<PathNode> link =
                        new LinkedList<PathNode>();
                    link.add(node1);
                    link.add(node2);
                    pathLists.add(link);
                }
                treeIdx++;
            }
        }
        
        return pathLists;
    }

    private void augmentFlow(FlowNetwork gFlow, 
        List<LinkedList<PathNode>> cPaths,
        List<Integer> surplus, List<Integer> deficit) {

        int hBefore = surplus.size();
        assert(hBefore == deficit.size());
        
        // see pg 63, Sect 8.4
        
        // - reverse forward-backward of links along paths
        //   AND idle-saturated states of the underlying arcs
        
        for (LinkedList<PathNode> link : cPaths) {
            assert(link.size() == 2);
            PathNode node2 = link.pollFirst();
            PathNode node1 = link.pollFirst();
            assert(link.isEmpty());
            int idx1 = ((Integer)node1.getData()).intValue();
            int idx2 = ((Integer)node2.getData()).intValue();
        
            if (node1 instanceof SourceNode) {
                assert(!(node2 instanceof SourceNode));
                gFlow.augmentSourceToLeftFlowAndArc(idx2);
            } else if (node1 instanceof SinkNode) {
                assert(!(node2 instanceof SinkNode));
                gFlow.augmentRightToSinkFlowAndArc(idx2);
            } else if (node1 instanceof LeftNode) {
                assert(!(node2 instanceof LeftNode));
                if (node2 instanceof RightNode) {
                    gFlow.augmentFlowAndArc(idx1, idx2);
                } else if (node2 instanceof SourceNode) {
                    gFlow.augmentSourceToLeftFlowAndArc(idx1);
                } else if (node2 instanceof SinkNode) {
                    gFlow.augmentRightToSinkFlowAndArc(idx1);
                }
            } else {
                assert(!(node2 instanceof RightNode));
                if (node2 instanceof LeftNode) {
                    gFlow.augmentFlowAndArc(idx2, idx1);
                } else if (node2 instanceof SourceNode) {
                    gFlow.augmentSourceToLeftFlowAndArc(idx1);
                } else if (node2 instanceof SinkNode) {
                    gFlow.augmentRightToSinkFlowAndArc(idx1);
                }
            }                
        }
        
        assert(gFlow.printFlowValueIncludingSrcSnk(hBefore));

        surplus.clear();
        deficit.clear();
        //gFlow.getMatchedLeftRight(surplus, deficit);
        gFlow.getSurplusLeftIndexes(surplus);
        gFlow.getDeficitRightIndexes(deficit);
        
        log.info("nSurplus=" + surplus.size() + " nDeficit="
            + deficit.size());
        
        int hAfter = surplus.size();
        assert((hBefore - hAfter) == cPaths.size());
    }

    private Map<PathNode, Set<PathNode>> createNewAdjacencyMap(
        List<LinkedList<PathNode>> pathLinkLists,
        int nLeft, int nRight, int sinkNodeIdx, int sourceNodeIdx) {
        
        Map<Integer, LeftNode> leftNodes = 
            new HashMap<Integer, LeftNode>();
        
        Map<Integer, RightNode> rightNodes = 
            new HashMap<Integer, RightNode>();
        
        SinkNode sinkNode = new SinkNode();
        sinkNode.setData(Integer.valueOf(sinkNodeIdx));
        
        SourceNode sourceNode = new SourceNode();
        sourceNode.setData(Integer.valueOf(sourceNodeIdx));
        
        Map<PathNode, Set<PathNode>> pathLinksMap =
            new HashMap<PathNode, Set<PathNode>>();
        for (LinkedList<PathNode> link : pathLinkLists) {
            PathNode node1 = link.pollFirst();
            PathNode node2 = link.pollFirst();
            assert(link.isEmpty());
            // replace node1 and node2 
            //     with nodes above to have same identity
            //     and field m
            Integer index1 = (Integer)node1.getData();
            Integer index2 = (Integer)node2.getData();
            if (node1 instanceof SinkNode) {
                node1 = sinkNode;
            } else if (node1 instanceof SourceNode) {
                node1 = sourceNode;
            } else if (node1 instanceof LeftNode) {
                node1 = leftNodes.get(index1);
                if (node1 == null) {
                    node1 = new LeftNode();
                    node1.setData(index1);
                    node1.m = 0;
                    leftNodes.put(index1, (LeftNode)node1);
                }
            } else {
                node1 = rightNodes.get(index1);
                if (node1 == null) {
                    node1 = new RightNode();
                    node1.setData(index1);
                    node1.m = 0;
                    rightNodes.put(index1, (RightNode)node1);
                }
            }
            if (node2 instanceof SinkNode) {
                node2 = sinkNode;
            } else if (node2 instanceof SourceNode) {
                node2 = sourceNode;
            } else if (node2 instanceof LeftNode) {
                node2 = leftNodes.get(index2);
                if (node2 == null) {
                    node2 = new LeftNode();
                    node2.setData(index2);
                    node2.m = 0;
                    leftNodes.put(index2, (LeftNode)node2);
                }
            } else {
                node2 = rightNodes.get(index2);
                if (node2 == null) {
                    node2 = new RightNode();
                    node2.setData(index2);
                    node2.m = 0;
                    rightNodes.put(index2, (RightNode)node2);
                }
            }
            
            Set<PathNode> values = pathLinksMap.get(node1);
            if (values == null) {
                values = new HashSet<PathNode>();
                pathLinksMap.put(node1, values);
            }
            values.add(node2);
        }

        return pathLinksMap;
    }
    
    private void makeSurplusAndDeficitSubSets(
        Map<PathNode, Set<PathNode>> pathLinksMap, 
        List<Integer> surplus, Set<Integer> deficit, 
        Set<LeftNode> outputSurplus, 
        Set<RightNode> outputDeficit) {

        Set<Integer> sSet = new HashSet<Integer>(surplus);
        Set<Integer> dSet = new HashSet<Integer>(deficit);
        
        for (Entry<PathNode, Set<PathNode>> entry : pathLinksMap.entrySet()) {
            PathNode node1 = entry.getKey();
            Integer index1 = (Integer)node1.getData();
            if (node1 instanceof LeftNode && sSet.contains(index1)) {
                outputSurplus.add((LeftNode)node1);
            } else if (node1 instanceof RightNode && dSet.contains(index1)) {
                outputDeficit.add((RightNode)node1);
            }
            for (PathNode node2 : entry.getValue()) {
                Integer index2 = (Integer)node2.getData();
                if (node2 instanceof LeftNode && sSet.contains(index2)) {
                    outputSurplus.add((LeftNode)node2);
                } else if (node2 instanceof RightNode && dSet.contains(index2)) {
                    outputDeficit.add((RightNode)node2);
                }
            }
        }
    
    }

}
