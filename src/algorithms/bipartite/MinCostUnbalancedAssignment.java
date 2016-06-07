package algorithms.bipartite;

import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.HeapNode;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.TrioInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
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
 * A solver for the min-cost, unbalanced, weighted bipartite
 * assignment problem that uses weight scaling
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
 by Ramshaw and Tarjan, 2012.

 <e>NOTE: for the best performance, the user should make sure
 that the maximum cost in their graph is less than roughly
 46340 because that is the limit wherein the an internal
 data structure switches to using an O(lg2(N_nodes))
 * extract min operation rather than the approximiately
 * O(1) for smaller maximum costs.
 </b>
* 
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
 * 
 * A strength of the algorithm to note is that it does not
 * artificially double the graph in order to handle
 * unequally sized left and right sets (it does not use
 * Bipartite double cover).
 </pre>
 * @author nichole
 */
public class MinCostUnbalancedAssignment {

    private Logger log = Logger.getLogger(this.getClass().getName());

    private FlowNetwork finalFN = null;

    /**
     * class specializing a fibonacci heap node to identify
     * a left node, a.k.a. X node
     */
    static class LeftNode extends PathNode {
        public LeftNode() {
            this.id = "LeftNode";
        }
        @Override
        public PathNode construct() {
            return new LeftNode();
        }
    }
    
    static abstract class PathNode extends HeapNode {
        PathNode pathPredecessor = null;
        LeftNode topPredecessor = null;
        // m = 0 is unmarked, m=1 is marked
        int m = 0;
        public abstract PathNode construct();
        // id is only used in toString fr debug statements
        String id = "";        
        public PathNode copy() {
            PathNode node = construct();
            node.setKey(getKey());
            node.setData(getData());
            if (pathPredecessor != null) {
                PathNode p = pathPredecessor;
                PathNode pNode = node;
                while (p != null) {
                    pNode.pathPredecessor = p.copy();
                    pNode = pNode.pathPredecessor;
                    p = p.pathPredecessor;
                }
            }
            if (topPredecessor != null) {
                node.topPredecessor = topPredecessor;
            }
            return node;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append(id).append(" key=").append(Long.toString(getKey()));
            sb.append(" index=").append(getData().toString());
            PathNode prev = pathPredecessor;
            while (prev != null) {
                String str = String.format("\n   [prev=%s %d %s]", 
                    prev.id, (int)prev.getKey(), prev.getData());
                sb.append(str);
                prev = prev.pathPredecessor;
            }
            if (topPredecessor != null) {
                sb.append(" *top=").append(topPredecessor.toString());
            }
            return sb.toString();
        }
    }
    
    /**
     * class specializing a fibonacci heap node to identify
     * a right node, a.k.a. Y node
     */
    static class RightNode extends PathNode {
        public RightNode() {
            this.id = "RightNode";
        }
        @Override
        public PathNode construct() {
            return new RightNode();
        }
    }
    
    /**
     * class specializing a fibonacci heap node to identify
     * a source
     */
    static class SourceNode extends PathNode {
        public SourceNode() {
            this.id = "SourceNode";
        }        
        @Override
        public PathNode construct() {
            return new SourceNode();
        }
    }
    
    /**
     * class specializing a fibonacci heap node to identify
     * a sink node
     */
    static class SinkNode extends PathNode {
        public SinkNode() {
            this.id = "SinkNode";
        }        
        @Override
        public PathNode construct() {
            return new SinkNode();
        }
    }
    
    /*
    input to FlowAssign is a bipartite graph G with integer 
    edge weights greater than zero and a target size t. 
    
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
     * match the left and right vertices in graph g by
     * minimum cost assignment and return the mappings.
     * 
     * @param g bipartite graph with integer weights having
     * values greater than zero.
     * @return map of indexes of left nodes matched to indexes of
     * right nodes
     */
    public Map<Integer, Integer> flowAssign(Graph g) {

        validateGraph(g);
        
        if (g.getNLeft() == 1 && g.getNRight() == 1) {
            return singleNodesSolution(g);
        }
        
        // a first guess at the maximal matching size
        int sz = Math.min(g.getNLeft(), g.getNRight());

long t0 = System.currentTimeMillis();        
        // hopcroft-karp produces a maximal matching of nodes
        // without using edge weights, just uses connectivity
        Map<Integer, Integer> m = hopcroftKarp(g, sz);
long t1 = System.currentTimeMillis();
long tSec = (t1 - t0)/1000;
System.out.println(tSec + " sec for hopcroftkarp");
        
        // if the code was given a graph created without 
        // using source and sink, need to transform that here.
        if (g.getSourceNode() == -1) {
            g = g.copyToCreateSourceSink();
        }
        
        /*
        {//DEBUG
            log.info("hopcroft-karp matches:");
            for (Entry<Integer, Integer> entry : m.entrySet()) {
                log.info("hk matched " + entry.getKey() + " to "
                    + entry.getValue());
            }
        }*/
        
        FlowNetwork gFlow = new FlowNetwork(g, m);
        //assert(gFlow.printFlowValueIncludingSrcSnk(m.size()));
        
        log.info("init FlowNetwork:");
        //gFlow.printNetCosts();
        
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
        
        float eps = 1.f + (int)Math.floor(eps_up);
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
        
        // S = left nodes matched in gFlow
        List<Integer> surplus = new ArrayList<Integer>();
        
        // D = right nodes matched in gFlow
        List<Integer> deficit = new ArrayList<Integer>();
        
        gFlow.getMatchedLeftRight(surplus, deficit);
        
        while ((eps > eps_down) && (nIterR < 2*rIter)) {
            
            log.info("nIterR=" + nIterR + " s=" + s + " eps=" + eps);
            
            // pg 44, assertions I1, I2, I3, and I4
            assert(gFlow.assertFlowValue(s));
            assert(gFlow.assertPricesAreQuantizedEps(eps));
            // if using a smaller eps than maxC, cannot assert
            //    these on first round
            if ((nIterR > 0) || (eps > gFlow.getMaxC())) {
                assert(gFlow.integralFlowIsEpsProper(eps));
                assert(gFlow.assertSaturatedBipartiteIsEpsSnug(eps));
            }
            
            eps /= ((float)q);
            
            int ext = refine(gFlow, s, eps, q, surplus, deficit);
            
            if (ext > 0) {
                m = gFlow.extractMatches();
                roundFinalPrices(gFlow, eps_down);
                //gFlow.printFlowValueIncludingSrcSnk(m.size());
                finalFN = gFlow;
                return m;
            }
            
            ++nIterR;
        }
        
        // assert nIter is approx log_q(s * maxC)

t0 = System.currentTimeMillis(); 
        // round prices to integers that make all arcs proper
        roundFinalPrices(gFlow, eps_down);
t1 = System.currentTimeMillis();  
tSec = (t1 - t0)/1000;
System.out.println(tSec + " sec for roundFinalPrices");
        
        m = gFlow.extractMatches();
        
        finalFN = gFlow;
        
        return m;
    }
    
    protected int refine(FlowNetwork gFlow, int s, float eps,
        int q, List<Integer> surplus, List<Integer> deficit) {
        
        log.fine("at start of refine, s=" + s + " eps=" + eps
            + " q=" + q);
        
        //assert(gFlow.printFlowValueIncludingSrcSnk(s));
        
        // set the flow of saturated bipartite arcs to 0
        for (int i = 0; i < surplus.size(); ++i) {
            int idx1 = surplus.get(i);
            int idx2 = deficit.get(i);
            gFlow.getFlow().put(new PairInt(idx1, idx2), 
                Float.valueOf(0));
            log.fine("surplus idx=" + Integer.toString(idx1));
        }
        
        /*
        see Figure 7.4 on pg 53.
        raise prices so that every arc in 
        gFlow becomes eps-proper, for the resulting pseudoflow 
        f and for the new, smaller value of eps.
        */
long t00 = System.currentTimeMillis();
        gFlow.raisePricesUntilEpsProper(eps, q);
long t11 = System.currentTimeMillis();
long tSec = (t11 - t00)/1000;
System.out.println(tSec + " sec for raisePricesUntilEpsProper");

        log.fine("after raise prices:");
        //gFlow.printNetCosts();
        //assert(gFlow.printFlowValueIncludingSrcSnk(s));        
        
        //in [0] holds the length of the terminating deficit
        // node which is also the forest index it was
        // added to.
        int[] terminatingDeficitIdx = new int[1];
        int h = s;
        int nHIter = 0;
t00 = System.currentTimeMillis();
        while (h > 0) {
            
            log.info("nHIter=" + nHIter + " h=" + h);
            
            ResidualDigraph2 rF = new ResidualDigraph2(gFlow);
long t0 = System.currentTimeMillis();
            // build a shortest-path forest from the current surpluses S, 
            // stopping when a current deficit in D is reached;
            // (pg 55)

            Forest forest = 
                buildForest2(gFlow, rF, surplus, deficit, eps,
                terminatingDeficitIdx);
long t1 = System.currentTimeMillis();
tSec = (t1 - t0)/10;
System.out.println(tSec + " 100ths of sec for buildForest2");

            //debug(forest);
           
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
            
            //assert(gFlow.printFlowValueIncludingSrcSnk(s));
         
            log.fine("l(termIndex)=" + terminatingDeficitIdx[0]);
            
t0 = System.currentTimeMillis();

            modifyPricesAndPathLengths(gFlow, rF, forest, 
                terminatingDeficitIdx[0], eps);

t1 = System.currentTimeMillis();
tSec = (t1 - t0)/10;
System.out.println(tSec + " 100ths of sec for modifyPricesAndPathLengths");

            //log.info("after modify prices and link lengths");
            //debug(forest);
            
            //assert from pg 52 
            //       I1', I2, I3, I4, on FlowNetwork
            //       and I5 on ResidualDigraph2
            
            assert(gFlow.assertFlowValueIncludingSrcSnk(s));
            assert(gFlow.assertPricesAreQuantizedEps(eps));
            assert(gFlow.integralFlowIsEpsProper(eps));
            assert(gFlow.assertSaturatedBipartiteIsEpsSnug(eps));
            // assert I5: Rf has no cycles of length zero.
            
            //Let Rf0 denote that subgraph of the residual digraph 
            //formed by links of length zero.
t0 = System.currentTimeMillis();
            List<LinkedList<PathNode>> zeroLengthLinks =
                extractZeroLengthLinks(forest);
t1 = System.currentTimeMillis();  
tSec = (t1 - t0)/10;
System.out.println(tSec 
+ " 100ths of sec for extractZeroLengthLinks");

            log.fine("srching for zero length paths:");
            //debug(forest);
            //debug(zeroLengthLinks);
            
            //assert(gFlow.printSurplusAndDeficit());
            
t0 = System.currentTimeMillis();              
            // --- Sect 8.3, create maximal set of compatible augmenting paths
            List<LinkedList<PathNode>> cPaths =
                findMaximalSetOfCompatiblePaths(
                gFlow.getNLeft(), gFlow.getNRight(),
                gFlow.getSourceNode(), gFlow.getSinkNode(),
                zeroLengthLinks,surplus, 
                new HashSet<Integer>(deficit));
t1 = System.currentTimeMillis();  
tSec = (t1 - t0)/10;
System.out.println(tSec 
+ " 100ths of sec for findMaximalSetOfCompatiblePaths");

            log.fine("cPaths.size=" + cPaths.size());
            if (cPaths.isEmpty()) {
                
                log.warning("did not find an augmenting path.  h=" + h);
                
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
            
            //debug(cPaths);
  
t0 = System.currentTimeMillis();              
            //augment f along each of the paths in P in turn, thereby 
            //   reducing |S| = |D| = h by |P|;
            //NOTE that surplus and deficit are modified and updated
            // within augmentFlow
            augmentFlow(gFlow, cPaths, surplus, deficit);
t1 = System.currentTimeMillis();  
tSec = (t1 - t0)/10;
System.out.println(tSec + " 100ths of sec for augmentFlow");

            //log.info("after augment flow:");
            //gFlow.printSaturatedLinks();

//surplus.clear();deficit.clear();
//gFlow.getMatchedLeftRight(surplus, deficit);

            h = surplus.size();
            
            // pg 63 assert I1', I2, I34  (not I3, but I4?)
            assert(gFlow.assertFlowValueIncludingSrcSnk(s));
            assert(gFlow.assertPricesAreQuantizedEps(eps));
            assert(gFlow.integralFlowIsEpsProper(eps));
            
            ++nHIter;
        }
t11 = System.currentTimeMillis();
tSec = (t11 - t00)/1000;
System.out.println(tSec + " sec for h block, nHIter=" + nHIter);

        return 0;
    }
    
    public static class Forest {
    
        // replaces the array DoubleLinkedCircularList[]
        // with a sparsely populated structure and
        // an ordered list of keys
        
        //NOTE that the use of this structure is monotonically
        // increasing,
        // that is, the keys for inserts are always
        // the same or increasing, so one only needs to compare
        // with last item in orderedKeys.
     
        private final int upperKeyLimit;
        
        // only allowing a remove operation for top of list
        private LinkedList<Integer> orderedKeys = new LinkedList<Integer>();
        
        private Map<Integer, DoubleLinkedCircularList> forest =
            new HashMap<Integer, DoubleLinkedCircularList>();
       
        public Forest(int upperKeyLimit) {
            this.upperKeyLimit = upperKeyLimit;
        }
        
        public DoubleLinkedCircularList get(int key) {
            return forest.get(Integer.valueOf(key));
        }
        
        public void removeFirstItem(int key) {
           Integer firstKey = orderedKeys.getFirst();
           if (firstKey == null) {
               return;
           }
           if (firstKey.intValue() == key) {
               orderedKeys.removeFirst();
               forest.remove(firstKey);
           }
        }
        
        public long add(PathNode node, long lastKey) {
            
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

            if (node.getKey() < upperKeyLimit) {
                
                int k = (int) node.getKey();
                Integer key = Integer.valueOf(k);
                                
                DoubleLinkedCircularList list = forest.get(key);
                
                if (list == null) {
                    list = new DoubleLinkedCircularList();
                    orderedKeys.add(key);
                    forest.put(key, list);
                }

                //NOTE: since some of the nodes are still possibly
                // nodes still present in the fibonacci heap,
                // one should copy nodes here to help not corrupt the
                // heap nodes.
                node = node.copy();

                list.insert(node);

                lastKey = node.getKey();
            }

            return lastKey;
        }
        
        public List<Integer> getKeys() {
            return new ArrayList<Integer>(orderedKeys);
        }
    }
    
    /**
     * 
     * @param gFlow
     * @param rF
     * @param surplus
     * @param deficit
     * @param eps
     * @param terminatingDeficitIdx output array of size 1
     * to hold the length of the key of the terminating deficit
     * node (which is also the index of the forest tree it 
     * was added to).
     * @return 
     */
    protected Forest buildForest2(
        final FlowNetwork gFlow, ResidualDigraph2 rF,
        List<Integer> surplus, List<Integer> deficit, float eps,
        int[] terminatingDeficitIdx) {
        
        log.fine("buildForest2");
        
long t0 = System.currentTimeMillis();

        Set<Integer> d = new HashSet<Integer>(deficit);
    
        //TODO: revisit this. 
        int lambda = (int)Math.ceil(3 * (gFlow.getMaxC()/eps));
        if (lambda < 4) {
            lambda = 4;
        }
        log.fine("buildForest2 forest length lambda is set to " + 
            lambda);
        
        long lastKey = -1;
        
        // sparsely populated holder for the 
        // DoubleLinkedCircularList trees
        Forest forest = new Forest(lambda);

        //TODO: revisit this
        // need to use a limit which will always include
        // the shortest path length at any time.
        // they're quantized w/ eps.
        MinHeapForRT2012 minHeap = new MinHeapForRT2012(lambda);
         
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
            minHeap.insert(sNode);
        }
        
long t1 = System.currentTimeMillis();
long tSec = (t1 - t0)/10;
System.out.println(tSec + " 100ths of sec for "
    + "init in buildForest2 heap.size=" +
    minHeap.getNumberOfNodes());

        do {
            PathNode node1 = minHeap.extractMin();
            if (node1 == null) {
                terminatingDeficitIdx[0] = (int)lastKey;
                break;
            }
            Integer index1 = (Integer)node1.getData();
            int idx1 = index1.intValue();
         
            log.fine("extractMin = " + node1.toString());
            
            boolean node1IsLeft = (node1 instanceof LeftNode);
            boolean nodeIsSource = (node1 instanceof SourceNode);
            boolean nodeIsSink = (node1 instanceof SinkNode);
            
            //scan all links leaving node1:
            if (nodeIsSource) {        
                // scan forward source links
                for (Integer index2 : rF.getForwardLinksSourceRM()) {
                    handlePlusLink(minHeap, 
                        node1, leftNodes.get(index2),
                        gFlow.calcSourceNetCost(index2.intValue()),
                        lambda, eps); 
                }
            } else if (nodeIsSink) {
                // scan backward sink links
                for (Integer index2 : rF.getBackwardLinksSinkRM()) {
                    handleMinusLink(minHeap, node1, 
                        rightNodes.get(index2), 
                        gFlow.calcSinkNetCost(index2.intValue()), 
                        lambda, eps);
                }
            } else if (node1IsLeft) {
                // scan the bipartite arcs forward
                Set<Integer> indexes2 = rF.getForwardLinksRM().get(
                    index1);
                if (indexes2 != null) {
                    for (Integer index2 : indexes2) {
                        handlePlusLink(minHeap, node1, 
                            rightNodes.get(index2), 
                            gFlow.calcNetCost(idx1, index2.intValue()),
                            lambda, eps);
                    }
                }
           
                // if there's a source link
                if (rF.getBackwardLinksSourceRM().contains(index1)) {
                    // insert a copy of the source node
                    PathNode sNode2 = sourceNode.copy();
                    handleMinusLink(minHeap, node1, sNode2, 
                        gFlow.calcSourceNetCost(idx1),
                        lambda, eps);
                }
            } else {
                // node1 is a RightNode
                Integer index2 = rF.getBackwardLinksRM().get(index1);
                if (index2 != null) {
                    handleMinusLink(minHeap, node1, 
                       leftNodes.get(index2), 
                       gFlow.calcNetCost(index2.intValue(), idx1),
                       lambda, eps);                     
                }
                // if there is a sink link
                if (rF.getForwardLinksSinkRM().contains(index1)) {
                    // insert a copy of the sink node
                    PathNode sNode2 = sinkNode.copy();
                    handlePlusLink(minHeap, node1, sNode2, 
                        gFlow.calcSinkNetCost(idx1),
                        lambda, eps);
                }
            }
            
            log.fine("add to forest key=" + node1.toString());
            
            //add v to the forest;
            lastKey = forest.add(node1, lastKey);

            if (d.contains(index1) && !node1IsLeft) {
                terminatingDeficitIdx[0] = (int)lastKey;
                if (lastKey != (int)node1.getKey()) {
                    throw new IllegalArgumentException(
                    "forest was not made large enough to hold last key");
                }
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
                
        if (true) { 
            //runtime complexity O(m * sqrt(n))
            HopcroftKarp hk = new HopcroftKarp();
            int[] matched = hk.hopcroftKarpV0(new GraphWithoutWeights(g));
            log.fine("matched=" + Arrays.toString(matched));
            Map<Integer, Integer> m = new HashMap<Integer, Integer>();
            for (int i = 0; i < matched.length; ++i) {
                int v = matched[i];
                if (v > -1) {
                    m.put(Integer.valueOf(i), Integer.valueOf(v));
                }
            }
            return m;
        }
       
        HopcroftKarpRT2012 hk = new HopcroftKarpRT2012();
        Map<Integer, Integer> m = hk.findMaxMatching(g, s);
       
        return m;
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

        Integer v2 = rM.getBackwardLinksRM().get(rightIndex);
        
        if (forwardFound) {
            
            if (v2 != null) {
                return;
            }

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
    
    public static List<PathNode> extractNodes(PathNode node) {
        
        List<PathNode> nodes = new ArrayList<PathNode>();
        
        PathNode node1 = node;
        while (node1 != null) {
            nodes.add(node1);
            node1 = node1.pathPredecessor;
        }
                
        return nodes;
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
                    log.info("forest2[" + key + "] tree branch[" 
                        + j + "] node[" + ii + "]=" + node1.toString());
                }                    
                j++;
            }
        }            
    }
    
    private void handlePlusLink(MinHeapForRT2012 minHeap, 
        PathNode node1, PathNode node2, float cp, int lambda, 
        float eps) {
    
        long l1 = node1.getKey();
        
        long lp = (long) Math.ceil(cp / eps);
        long lTot = l1 + lp;
        long lOld = node2.getKey();
        if ((lTot < lambda) && (lTot < lOld)) {
            node2.pathPredecessor = node1;
            if (lOld == Long.MAX_VALUE) {
                node2.setKey(lTot);
                minHeap.insert(node2);
            } else {
                minHeap.decreaseKey(node2, lTot);
            }
        }
    }
    
    private void handleMinusLink(MinHeapForRT2012 minHeap,
        PathNode node1, PathNode node2, float cp,
        int lambda, float eps) {
    
        long l1 = node1.getKey();
        
        long lp = 1 - (long) Math.ceil(cp / eps);
        long lTot = l1 + lp;
        long lOld = node2.getKey();
        if ((lTot < lambda) && (lTot < lOld)) {
            node2.pathPredecessor = node1;
            if (lOld == Long.MAX_VALUE) {
                node2.setKey(lTot);
                minHeap.insert(node2);
            } else {
                minHeap.decreaseKey(node2, lTot);
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
     
    private List<List<PathNode>> extractPathNodes(Forest forest,
        int forestIdx) {
        
        List<List<PathNode>> pathLists = new ArrayList<List<PathNode>>();
                    
        DoubleLinkedCircularList tree = forest.get(forestIdx);
        
        if (tree == null) {
            return null;
        }
        
        long n = tree.getNumberOfNodes();
        HeapNode node = tree.getSentinel();
        int branchIdx = 0;
        while (branchIdx < n) {
            node = node.getRight();
            List<PathNode> path = extractNodes((PathNode)node);
            Misc.<PathNode>reverse(path);
            pathLists.add(path);
            ++branchIdx;
        }
        
        return pathLists;
    }
    
    private void modifyPricesAndPathLengths(
        FlowNetwork gFlow, ResidualDigraph2 rF,
        Forest forest, int lt, float eps) {

        debug(forest);
        
        // extract the paths for forest[lt] and modify them
        //  while applying implied price changes to gFlow.
                
        /*
        sect 8.2, pg 57
       
        float cp = cost - pdX + pdY;
        normally: link lengths in residual digrph 2:
           - a forward link v->w has length
               lp(v->w) = Math.ceil(cp(v,w)/eps)
           - a backward link w->v has length
               lp(w->v) = 1 - Math.ceil(cp(v,w)/eps)
         
        for each left node, that is v, node in forest,
            p0d(v) := p(d) + eps*(l(term) - l(v))
        
          a node inserted into the forest at branch=1 or higher
          has a key and a predecessoer which is a maiden node
          or that predecessor has a predecessor which is, etc.
          the inserted node has a key greater than or equal to
          its predecessors and the inserted node has key l(term)
          if it was in the last branch added to.
          so changes to link lengths by + or - l(term) must 
          be - for those to enable them to be link=0 links
          and hence candidates for augmenting paths.
        
          more specifically,
             after the first round of buildForest2, the links
             inserted into the forest at forest[1] are
             "idle" links from surplus to deficit nodes.
             The inserted Right node w/ key=1 has a 
             predecessor Left node w/ key=0.
                Lft(0) -> Rgt(1) (is "idle" in this example)
             
                the Rgt key length which was assigned using
                   leftKey + Math.ceil(cp(left, right) / eps)
                   where cp = costXY - pdX + pdY
                -> reduce Rgt(1) key by l(term) - l(x), that is
                   by 1 - 0
                   -> by first increasing the price pd(X) (that is the
                      left node price)
                      by (l(term) - l(X))*eps
                      (that reduces the netcost along XY)
        
            or p0d(v) := p(d) + (l(term) - l(x))*eps
            pdX += (l(term) - l(x))*eps which is pdX += eps here
            then pdY += (l(term) - l(y))*eps doesn't change
            recalculated cost is lower by +eps
            so l(y) is lowered by 1 (from (l(term) - l(x)))

            that link is then handled...

            to address "saturated links" and multiple links in a path:

            N2     N1
            L <---- R
            cp = gFlow.calcNetCost(index2.intValue(), idx1)
            where cp = costXY - pdX + pdY
            lp = 1 - math.ceil(cp/eps)
            key of node2 = key of node1 + lp
            l(term) = node2.key
            pd_N2 += 0.  pd_N1 += (l(term) - l(N1))*eps
            cp is increased by (l(term) - l(N1))*eps
                so the term math.ceil(cp/eps) is + (l(term) - l(N1))
            so lp is -((l(term) - l(N1)))
            reducing key of node2 by (l(term) - l(N1))
        
            to continue that if there is a previous node to N1:
            N2     N1
            L <---- R  <---- prevL
        
                    then for N1, an idle link into it from node prevL.
                      we already have pd_N1 += (l(term) - l(N1))*eps
                      pd_prevL += (l(term) - l(prevL))*eps
                      cp change is -(l(term) - l(prevL))*eps
                           +(l(term) - l(N1))*eps
                           = (l(prevL)-l(N1))*eps
                      then the N1.key changes by (l(prevL)-l(N1))
                      which should bring it to zero.
        
            Looks like the price changes to the FlowNetwork gFlow
            should be applied using a single path at a time
            (so that price changes of nodes within that path are
            only applied once for that path).
        
            For another path with the same node in it,
            one can see that cannot re-fetch the prices from
            gFlow to re-calculate the netcost, but must instead
            only use the path link lengths to derive relative
            changes to apply to the gFlow prices.  gFlow is write only
            for these price changes.
        
        
             One thing that is not obvious is that if 
             there was a tree of paths at forest[1] and
             also at forest[2] and then l(term) is the key
             at root of forest[2] which is '2',
             paths from forest[1] should be discarded?
             since l(term) - l(forest[1].root.key) will be > 0
                 couldn't find that stated in the paper yet, so
             the changes applied to gFlow even from forest
             trees with smaller keys might still be necessary even
             though the resulting links will not be 0,
             but those may have redundant links to highest key tree.
             need to work further through examples...
             
        */
        
        List<List<PathNode>> ltPathLinks = null;
        
        if (forest.getKeys().size() > 2) {
            debug(forest);
            // if lt is > forestKey > 0, are there redundant paths
            //   present in different index trees?
            int z = 1;
        }
        
        List<Integer> forestIndexes = forest.getKeys();
        assert(forestIndexes.get(forestIndexes.size() - 1).intValue()
            == lt);
        
        // apply changes to the lt tree, then the loer index trees
        for (int idx = (forestIndexes.size() - 1); idx >= 0; --idx) {
            
            int forestIdx = forestIndexes.get(idx).intValue();
            
            List<List<PathNode>> pathLists = extractPathNodes(forest, 
                forestIdx);
            
            if (pathLists == null) {
                continue;
            }
            
            if (forestIdx == lt) {
                ltPathLinks = pathLists;
            }
        
            for (int pathIdx = 0; pathIdx < pathLists.size(); ++pathIdx) {
            
                List<PathNode> path = pathLists.get(pathIdx);
            
                Map<Integer, Integer> leftPriceIncreases =
                    new HashMap<Integer, Integer>();
            
                Map<Integer, Integer> rightPriceIncreases =
                    new HashMap<Integer, Integer>();
            
                for (int i = 0; i < (path.size() - 1); ++i) {
                    PathNode node1 = path.get(i);
                    PathNode node2 = path.get(i + 1);
                    int l1 = (int) node1.getKey();
                    int l2 = (int) node2.getKey();
                    int index1 = (Integer) node1.getData();
                    int index2 = (Integer) node2.getData();
                    
                    boolean node1IsLeft = (node1 instanceof LeftNode);
                    boolean node1IsSource = (node1 instanceof SourceNode);
                    boolean node1IsSink = (node1 instanceof SinkNode);
                    boolean node2IsLeft = (node2 instanceof LeftNode);
                    boolean node2IsSource = (node2 instanceof SourceNode);
                    boolean node2IsSink = (node2 instanceof SinkNode);
                    
                    // determine if node1 -> node2 is "idle" or "saturated"

                    if (node1IsSource) {
                        assert (node2IsLeft);
                        if (rF.getBackwardLinksSourceRM().contains(index2)) {
                            // saturated  node1(source) <--- node2(Left)

                        } else {
                            assert (rF.getForwardLinksSourceRM().contains(index2));
                            // idle  node1(source) --> node2(left)

                        }
                    } else if (node1IsSink) {
                        assert(node2 instanceof RightNode);
                        if (rF.getBackwardLinksSinkRM().contains(index2)) {
                            // "saturated" node2(right) --> node1(sink)
                            
                        } else {
                            assert(rF.getForwardLinksSinkRM().contains(index2));
                            // "idle" node2(right) --> node1(sink)
                        
                        }
                    } else if (node1IsLeft) {
                        assert(!(node2 instanceof LeftNode));
                        if (node2IsSource) {
                            if (rF.getBackwardLinksSourceRM().contains(index1)) {
                                // saturated  node2(source) <--- node1(Left)
                                
                            } else {
                                assert(rF.getForwardLinksSourceRM().contains(index1));
                                // idle  node2(source) --> node1(left)
                            
                            }
                        } else {
                            // else node2 is Right
                            if (rF.getBackwardLinksRM().containsKey(index2)
                                && rF.getBackwardLinksRM().get(index2).equals(index1)) {
                                // saturated link node1(left) <-- node2(right)
                            
                            } else {
                                assert(rF.getForwardLinksRM().get(index1)
                                    .contains(index2));
                                // idle link node1(left) --> node2(right)
                                
                            }
                        }
                    } else {
                        // else node1 is Right
                        assert(!(node2 instanceof RightNode));
                        if (node2IsSink) {
                            if (rF.getBackwardLinksSinkRM().contains(index1)) {
                                // "saturated" node1(right) <-- node2(sink)

                            } else {
                                assert (rF.getForwardLinksSinkRM().contains(index1));
                                // "idle" node1(right) --> node2(sink)
                                
                            }
                            
                        } else {
                            assert(node2IsLeft);
                            if (rF.getBackwardLinksRM().containsKey(index1)
                                && rF.getBackwardLinksRM().get(index1)
                                .equals(index2)) {
                                // "saturated" node2(left) <-- node1(right)
                                
                            } else {
                                assert (rF.getForwardLinksRM().get(index2) != null
                                && rF.getForwardLinksRM().get(index2)
                                .contains(index1));
                                // "idle" node2(left) --> node1(right) 
                                
                            }
                        }
                    }    
                }
            }
        }
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    private List<LinkedList<PathNode>>
        findMaximalSetOfCompatiblePaths(
            int nLeft, int nRight, int sourceNodeIdx,
            int sinkNodeIdx,
            List<LinkedList<PathNode>> pathLinkLists,
            Collection<Integer> surplus, Set<Integer> deficit) {
        
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
                log.fine("v=stack.top=" + v.toString());
                boolean isLeftNode = (v instanceof LeftNode);
                boolean isRightNode = (v instanceof RightNode);
                if (isLeftNode || isRightNode) {
                    v.m = 1;
                }
                if (isRightNode && def.contains((RightNode)v)) {
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
                        log.fine("stack.push w = " + w.toString());
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
        Forest forest) {

        // each linked list is a pair of 2 path nodes whose link length = 0
        List<LinkedList<PathNode>> pathLists = new
            ArrayList<LinkedList<PathNode>>();
        
        // traverse trees in the forest
        for (Integer forestKey : forest.getKeys()) {
            
            DoubleLinkedCircularList tree = forest.get(forestKey);
            
            long n = tree.getNumberOfNodes();
            HeapNode node = tree.getSentinel();
            int treeIdx = 0;
            while (treeIdx < n) {
                node = node.getRight();
                List<PathNode> path = extractNodes((PathNode)node);
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
        
        //assert(gFlow.printFlowValueIncludingSrcSnk(hBefore));

        surplus.clear();
        deficit.clear();
        gFlow.getSurplusLeftIndexes(surplus);
        gFlow.getDeficitRightIndexes(deficit);
        
        log.fine("nSurplus=" + surplus.size() + " nDeficit="
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
    
    public FlowNetwork getFinalFlowNetwork() {
        return finalFN;
    }
    
    private void makeSurplusAndDeficitSubSets(
        Map<PathNode, Set<PathNode>> pathLinksMap, 
        Collection<Integer> surplus, Set<Integer> deficit, 
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
    
    private void debugPath(List<PathNode> path) {
        StringBuilder sb = new StringBuilder("path=");
        for (PathNode node : path) {
            String str = String.format("%s %s", 
                node.id, node.getData().toString());
            sb.append(str).append(", ");
        }
        log.info(sb.toString());
    }
    
    private void validateGraph(Graph g) {

        StringBuilder errors = new StringBuilder();
        if (g.getNLeft() == 0) {
            errors.append("graph cannot have 0 left nodes.");
        } 
        if (g.getNRight() == 0) {
            errors.append(" graph cannot have 0 rght nodes.");
        }
        for (Entry<PairInt, Integer> entry : g.getEdgeWeights().entrySet()) {
            if (entry.getValue().intValue() < 1) {
                errors.append(" all edge weights must be > 0.");
                break;
            }
        }
        
        if (g.getEdgeWeights().isEmpty()) {
            errors.append(" g must have at least one edge.");
        }
        
        if (errors.length() > 0) {
            throw new IllegalArgumentException(errors.toString());
        }
    }

    private Map<Integer, Integer> singleNodesSolution(Graph g) {
    
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
        
        if (g.getEdgeWeights().containsKey(new PairInt(0, 0))) {
            m.put(Integer.valueOf(0), Integer.valueOf(0));
        }
        
        throw new IllegalStateException("graph only had 1 left"
            + " and right node, but no edge between them.");
    }
}
