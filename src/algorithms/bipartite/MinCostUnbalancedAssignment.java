package algorithms.bipartite;

import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.Heap;
import algorithms.imageProcessing.HeapNode;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
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
        
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("LeftNode key=").append(Long.toString(getKey()));
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
    
    public static class PathNode extends HeapNode {
        PathNode pathPredecessor = null;
    }
    
    /**
     * class specializing a fibonacci heap node to identify
     * a right node, a.k.a. Y node
     */
    private class RightNode extends PathNode {

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("RightNode key=").append(Long.toString(getKey()));
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
     * @param g bipartite graph with integral weights
     * @return 
     */
    public Map<Integer, Integer> flowAssign(Graph g) {

        //TODO: add size restriction t
        Map<Integer, Integer> m = hopcroftKarp(g);
        
        FlowNetwork gFlow = new FlowNetwork(g, m);
        
        int s = m.size();
         
        // q >= 2
        int q = 2;
        
        // since all costs are integers, can set eps < 1/6s
        // where s is size of m
        // NOTE: compare this to the maxC in gFlow.  pg 32 suggests eps=maxC
        float eps = 1.f/(6.f * (float)s);
        
        float eBar = 1.f + (float)Math.floor(Math.log(
            gFlow.getMaxC())/Math.log(q));
        //e_bar = 1 + log_q(C)
        //eps_bar = q^(e_bar)
        float epsBar = (float)Math.pow(eBar, q);
        log.info("eps=" + eps + " epsBar=" + epsBar + " maxC=" + 
            gFlow.getMaxC());
        
        assert(eps > gFlow.getMaxC());
       
        // all nodes V in gFlow have prices = 0
        
        while (eps > epsBar) {
            
            // pg 44, assertions I1, I2, I3, and I4
            assert(gFlow.assertFlowValue(s));
            assert(gFlow.assertPricesAreQuantizedEps(eps));
            assert(gFlow.integralFlowIsEpsProper(eps));
            assert(gFlow.assertSaturatedBipartiteIsEpsSnug(eps));
            
            eps /= (float)q;
            
            refine(gFlow, s, eps, q);
        }
        
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
         
        //in [0] holds the index for the deficit node which
        // terminates the surplus loop,
        //in [1] holds the surplus vertex index at the root 
        // of the tree that the terminating deficit node was 
        // added to.
        int[] terminatingDeficitIdx = new int[1];
        int h = 2;
        while (h > 0) {
          
            ResidualDigraph2 rF = new ResidualDigraph2(gFlow);

            // build a shortest-path forest from the current surpluses S, 
            // stopping when a current deficit in D is reached;
            // (pg 55)
            DoubleLinkedCircularList[] forest = 
                buildForest2(gFlow, rF, surplus, deficit, eps,
                terminatingDeficitIdx);
                        
            //raise prices at forest nodes by multiples of ε, 
            //shortening the discovered augmenting path to length 0;
            raisePrices2(gFlow, rF, forest, terminatingDeficitIdx,
                eps);
            
            /*
            assert from pg 52 
                   I1', I2, I3, I4, on FlowNetwork
                   and I5 on ResidualDigraph2
            */
            
            /*
            find a maximal set P of length-0 augmenting paths 
               that are compatible, as defined in Section 8.3;
            augment f along each of the paths in P in turn, thereby 
               reducing |S| = |D| = h by |P|;
            */
        }
       
        throw new UnsupportedOperationException("not yet implemented");
    /*    
        lengths of links:
       - quantization into units of eps is used for net cost:
         for a forward link: lp(v->w) = Math.ceil(cp(v, w)/eps)
         for a backward link: lp(w->v) = 1 - Math.ceil(cp(v, w)/eps)
         both results are >= 0        
    */
        /*
        a surplus of f is a node (not the sink node), that has
           entering flow > exiting flow.
        a deficit of f is a node (not the source node), that has
           exiting flow > entering flow.
        
        For a woman x in X, 
            let the left stub to x be the pseudoflow that 
            saturates the left-dummy arc |- -> x, 
            but leaves all other arcs idle. 
        For a man y in Y, 
            the right stub from 
            y saturates only the right-dummy arc
            y -> -|. 
        Any pseudoflow f that arises in Refine is the sum 
           of some flow, some left-stubs, and some right-stubs. 
           The flow component, which we denote f, encodes 
           the partial matching that Refine has constructed so far, 
           during this scaling phase. 
           - We initialize fˆ to zero, so this matching starts 
             out empty. The left-stubs remember those women who 
             were matched at the end of the previous phase and 
             who have not yet been either matched or replaced 
             during this phase.
             - Those women are the surpluses of the pseudoflow f, 
             and they constitute the set S. 
             - The right-stubs remember the previously matched 
               men in a similar way. Those men are the deficits of f, 
             and they constitute D.
             -- then |fˆ| = s - h, where h = |S| == |D|
                (if h = 0. fˆ = f so is an integral flow)
        */
    }
    
    /**
     * 
     * @param gFlow
     * @param surplus
     * @param deficit
     * @param eps
     * @param terminatingDeficitIdx output array of size 2
     * to hold first the index for the deficit node which
     * terminates the surplus loop, and then second to hold
     * the surplus vertex index at the root of the tree that
     * the terminating deficit node was added to.
     * @return 
     */
    protected DoubleLinkedCircularList[] buildForest2(
        final FlowNetwork gFlow, ResidualDigraph2 rF,
        List<Integer> surplus,
        List<Integer> deficit, float eps,
        int[] terminatingDeficitIdx) {
        
        if (true) {
            throw new IllegalArgumentException("not yet implemented");
        }
        
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
    
//these lengths should be in residual graph2
        
        Map<Integer, LeftNode> leftNodes = new HashMap<Integer, LeftNode>();
        Map<Integer, RightNode> rightNodes = new HashMap<Integer, RightNode>();
        for (int i = 0; i < gFlow.getNRight(); ++i) {
            Integer index = Integer.valueOf(i);
            RightNode node = new RightNode();
            node.setKey(Long.MAX_VALUE);
            node.setData(index);
            rightNodes.put(index, node);
        }
        for (int i = 0; i < gFlow.getNLeft(); ++i) {
            Integer index = Integer.valueOf(i);
            LeftNode node = new LeftNode();
            node.setKey(Long.MAX_VALUE);
            node.setData(index);
            leftNodes.put(index, node);
        }
                  
        for (Integer sigma : surplus) {
    
      //build a new node instead of get from list?
          
            LeftNode sNode = leftNodes.get(sigma);
            sNode.setKey(0);
            
            //insert(sigma, 0);
            insertIntoHeap(minHeap, sNode);
               
            PathNode node1 = extractMinFromHeap(minHeap);
            Integer index1 = (Integer)node1.getData();
            int idx1 = index1.intValue();
            final long l1 = node1.getKey();
            
            boolean node1IsLeft = (node1 instanceof LeftNode);

  //extracted node should not need to be copied
            
            PathNode node1Cp;
            if (node1IsLeft) {
                node1Cp = new LeftNode();
            } else {
                node1Cp = new RightNode();
            }
            node1Cp.setKey(node1.getKey());
            node1Cp.setData(node1.getData());
            node1Cp.pathPredecessor = node1.pathPredecessor;
            
            do {
                //scan:
                if (node1IsLeft) {
                    Set<Integer> indexes2 = rF.getForwardLinksRM().get(index1);
                    if (indexes2 != null) {
                        for (Integer index2 : indexes2) {
                            
                            RightNode node2 = rightNodes.get(index2);
                            float cp = gFlow.calcNetCost(idx1, index2.intValue());
                            long lp = (long)Math.ceil(cp/eps);
                            long lTot = l1 + lp;
                            long lOld = node2.getKey();
                            if ((lTot < lambda) && (lTot < lOld)) {
                                node2.pathPredecessor = node1Cp;
                                if (lOld == Long.MAX_VALUE) {
                                    node2.setKey(lTot);
                                    insertIntoHeap(minHeap, node2);
                                } else {
                                    decreaseKeyInHeap(minHeap, node2, lTot);
                                }
                            }
                            /*
                             L := l(v) + lp(v -> w); 
                             L_old := l(w);
                             if L <= lambda and L < L_old 
                                then set l(w) := L;
                                if Lold = inf
                                   then insert(w, L)
                                else 
                                   decrease-key(w, L);
                             */
                        }
                    }
                } else {
                    Integer index2 = rF.getBackwardLinksRM().get(index1); 
                    if (index2 != null) {
                        LeftNode node2 = leftNodes.get(index2);
                        float cp = gFlow.calcNetCost(index2.intValue(), idx1);
                        long lp = 1 - (long)Math.ceil(cp/eps);
                        long lTot = l1 + lp;
                        long lOld = node2.getKey();
                        if ((lTot < lambda) && (lTot < lOld)) {
                            node2.pathPredecessor = node1Cp; 
                            if (lOld == Long.MAX_VALUE) {
                                node2.setKey(lTot);
                                insertIntoHeap(minHeap, node2);
                            } else {
                                decreaseKeyInHeap(minHeap, node2, lTot);
                            }
                        }
                    }
                }
                //add v to the forest;
                addToForest(forest, node1Cp);
                
                terminatingDeficitIdx[0] = idx1;
                terminatingDeficitIdx[1] = 
                    ((Integer)(forest[(int)node1Cp.getKey()].getSentinel()
                    .getLeft().getData())).intValue();
            } while (!d.contains(index1));
        }
        
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
    
    //TODO: add t as limit for size
    protected Map<Integer, Integer> hopcroftKarp(Graph g) {
        
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
        
        if (true) {             
            //temporarily, replacing w/ O(m * sqrt(n))
            HopcroftKarp hk = new HopcroftKarp();
            int[] matched = hk.hopcroftKarpV0(createUnweightedGraph(g));
            log.info("matched=" + Arrays.toString(matched));
            for (int i = 0; i < matched.length; ++i) {
                int v = matched[i];
                if (v > -1) {
                    m.put(Integer.valueOf(i), Integer.valueOf(v));
                }
            }
            return m;
        }
                
        ResidualDigraph rM = createResidualGraph(g, m);
        
        return hopcroftKarp(g, rM);
    }
    
    protected Map<Integer, Integer> hopcroftKarp(Graph g, 
        ResidualDigraph rM) {
       
        if (true) {
            throw new UnsupportedOperationException("not yet implemented");
        }
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
        
        while (true) {
        
            DoubleLinkedCircularList[] augmentingPaths =
                buildForest(rM);

            if (allAreEmpty(augmentingPaths)) {
                return m;
            }

            // start at i=1, because i=0 is not alternating? 
            for (int i = 1; i < augmentingPaths.length; ++i) {

                DoubleLinkedCircularList tree = augmentingPaths[i];

                if (tree == null) {
                    continue;
                }

                //NOTE: not sure the logic is correct here.

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

                // HeapNode getLeft() traverses the tree nodes in FIFO order.
                
                // each item in tree is a HeapNode whose path
                //    is included as pathPredecessor of the node.
                
                long n = tree.getNumberOfNodes();
                HeapNode node = tree.getSentinel();
                int nCurrent = 0;
                for (int j = 0; j < n; ++j) {
                    
                    node = node.getLeft();
                    
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
                        
                        log.info("forest[" + i + "] tree branch[" 
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
                        
                        swapLinkExistence(rM, m, index1, index2);                        
                    }
                }
               
                //announce(M is a matching)
                log.info("m.size=" + m.size());
            }
        } 
        //return m;
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
     */
    protected DoubleLinkedCircularList[] buildForest(
        final ResidualDigraph rM) {
        
        //TODO: revisit this.
        int lambda = 2 * Math.min(rM.getLeftRM().size(), rM.getRightRM().size());
        
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
             node without a predecessor and havinf dist=0.
        
        
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
        
  //TODO: see Fredman and Tarjan [10]      
   
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
        
            log.info("heap.size=" + heap.getNumberOfNodes());
            
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
                
                // not necessary to update original key,
                //  but it does show presence in larger tree
                //  when node.key > forest key
                xNode.setKey(yNode.getKey());
                
                // make a copy in case it's already in forest
                LeftNode xNode2 = new LeftNode();
                xNode2.setKey(yNode.getKey());
                xNode2.setData(xNode.getData());
                xNode2.pathPredecessor = yNode;
                // any rightNodes inserted into heap get decreased keys
                scanAndAdd(heap, forest, rM, rightNodes, xNode2);
                
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
                    log.info(String.format("HEAP insert: %s",
                        yNode.toString()));
                } else {
                    heap.decreaseKey(yNode, lY);
                    log.info(String.format("HEAP decr: %s",
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
            
            String str = 
                node.getClass().getSimpleName().contains("Left") ?
                " Left" : " Right";
            log.info("add to forest[" + k + "] " + str 
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
    
    protected ResidualDigraph createResidualGraph(Graph g, 
        Map<Integer, Integer> m) {
                
        ResidualDigraph rM = new ResidualDigraph();
        
        for (int i = 0; i < g.getNLeft(); ++i) {
            rM.getLeftRM().add(Integer.valueOf(i));
        }
        
        for (int i = 0; i < g.getNRight(); ++i) {
            rM.getRightRM().add(Integer.valueOf(i));
        }
        
        for (Entry<PairInt, Integer> entry : g.getEdgeWeights().entrySet()) {
            
            PairInt p = entry.getKey();
            
            // presumably, edge weight of 0 is "not connected"
            if (entry.getValue().intValue() == 0) {
                continue;
            }
            
            Integer x = Integer.valueOf(p.getX());
            Integer y = Integer.valueOf(p.getY());
            
            Set<Integer> ys = rM.getForwardLinksRM().get(x);
            if (ys == null) {
                ys = new HashSet<Integer>();
                rM.getForwardLinksRM().put(x, ys);
            }
            
            ys.add(y);
        }
        
        for (Entry<Integer, Integer> entry : m.entrySet()) {
            Integer x = entry.getKey();
            Integer y = entry.getValue();
            
            if (rM.getForwardLinksRM().containsKey(x)) {
                rM.getForwardLinksRM().get(x).remove(y);
            }
            
            rM.getBackwardLinksRM().put(y, x);
        }

        return rM;
    }
        
    private GraphWithoutWeights createUnweightedGraph(Graph g) {

        GraphWithoutWeights g2 = new GraphWithoutWeights(
            g.getNLeft(), g.getNRight());
        Map<Integer, Set<Integer>> adjMap = g2.getAdjacencyMap();
        
        for (Entry<PairInt, Integer> entry : 
            g.getEdgeWeights().entrySet()) {
            
            Integer index1 = entry.getKey().getX();
            Integer index2 = entry.getKey().getY();
            
            Set<Integer> indexes2 = adjMap.get(index1);
            if (indexes2 == null) {
                indexes2 = new HashSet<Integer>();
                adjMap.put(index1, indexes2);
            }
            indexes2.add(index2);
        }
        
        return g2;
    }
    
    private void swapLinkExistence(ResidualDigraph rM, 
        Map<Integer, Integer> m, Integer leftIndex, 
        Integer rightIndex) {
        
        Set<Integer> rIndexes = rM.getForwardLinksRM().get(leftIndex);
        
        boolean forwardFound = (rIndexes != null) && rIndexes.contains(rightIndex);
        
        if (forwardFound) {
            
            // remove existing "idle" forward link
            rIndexes.remove(rightIndex);
            
            // create a backward link and matched mapping
            rM.getBackwardLinksRM().put(rightIndex, leftIndex);
            m.put(leftIndex, rightIndex);
            
            return;
        }
        
        // assert that a backward link exists
        Integer v2 = rM.getBackwardLinksRM().get(rightIndex);
        assert(v2 != null && v2.equals(leftIndex));
        
        // remove backwards link and mapping
        rM.getBackwardLinksRM().remove(rightIndex);
        m.remove(leftIndex, rightIndex);
        
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
        
        // NOTE: a buicket in the heap may already have a
        // node in it, so linked lists are used.
        // a doubly linked list is chosen so that the
        // list can be read as FIFO
        
        int key = (int)node.getKey();
        
        DoubleLinkedCircularList bucket = minHeap[key];
        if (bucket == null) {
            bucket = new DoubleLinkedCircularList();
            minHeap[key] = bucket;
        }
        
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

        int prevKey = (int)node2.getKey();
        minHeap[prevKey].remove(node2);
        
        node2.setKey(lTot);
        
        insertIntoHeap(minHeap, node2);
    }

    private void raisePrices2(FlowNetwork gFlow, ResidualDigraph2 rF,
        DoubleLinkedCircularList[] forest, 
        int[] terminatingDeficitIdx, float eps) {
        
         /*
        raise prices at forest nodes by multiples of ε, 
            shortening the discovered augmenting path to length 0;
            - In a round of price increases in the main loop of 
              Refine, raising the price pd(v) at some node v in NG 
              by eps lowers by 1 the length of any link in the 
              residual digraph Rf that leaves v and raises by 1 
              the length of any link that enters v.
            - A link v -> w leaving v is either forward or backward.
              If it is forward, the arc underlying it is the idle 
              arc v -> w. Raising the price at v by eps lowers the 
              net cost of this arc by eps, which, 
              by equation (7-3), lowers the length of the link by 1.
              If the link v -> w is backward, the arc underlying it 
              is the saturated arc w -> v. 
              Raising the price at v by eps raises the net cost of 
              this arc by eps, which, by equation (7-4), also lowers 
              the length of the link by 1.
            -- raise prices: 
              For each node v in the shortest-path forest, 
              we set the new (dispose) price pd'(v) by
              pd'(v) := p(d) + (l(termDefIdx) - l(v))*eps
              where termDefIdx is the deficit whose discovery 
              halted the growth of the shortest-path forest. 
              (reminder that path lengths are multiples of eps).
                 - Let sigma be the surplus at the root of the tree 
              that termDefIdx joins. 
            index (tree) that termDefIdx was last inserted into.  
              We have pd'(sigma) = pd(sigma) + l(termDefIdx)*eps, 
              but pd'(termDefIdx) = pd(termDefIdx). 
              So Prop 7-5 tells us that our price increases 
              shorten the path from sigma to termDefIdx by 
                 l(termDefIdx) length units. 
              Since l(termDefIdx) was the length of that path before 
              our price increases, its length after the increases 
              will be zero. 
              - If our invariants are preserved, 
              all of the links along that path must end up of 
              length 0, meaning that all of the underlying arcs 
              are eps-tight. As for our invariants, it’s clear 
              that I1' and I2 continue to hold. 
              But establishing the other three invariants takes 
              more work.
               For each node v in the network NG,
               define i(v) to be the multiple of eps by which 
               we raise the price pd(v). 
               -- So, for nodes v in the forest: 
                   i(v) := l(termDefIdx) - l(v), while, 
               -- for nodes v not in the forest:
                   i(v) := 0. 
               ** We then have the repricing formula 
                   pd'(v) = pd(v) + i(v)*eps, for all nodes v. 
               
               Using this repricing formula, we can express 
               the impact of our repricings on an arc
               v->w without needing to know whether or not the 
               nodes v and w lie in the forest. We have
                   cp(v, w) = c(v, w) - pd(v) + pd(w)
                   cp'(v, w) = c(v, w) - pd'(v) + pd'(w)
                 ==> cp'(v, w) = cp(v, w) + (i(w) - i(v))*eps
            -- after price incr,
                   assert from pg 52 
                   I1', I2, I3, I4, on FlowNetwork
                   and I5 on ResidualDigraph2
            -- 
        */
        
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
