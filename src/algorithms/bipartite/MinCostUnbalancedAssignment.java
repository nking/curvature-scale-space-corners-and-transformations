package algorithms.bipartite;

import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.HeapNode;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.QuadInt;
import algorithms.util.TrioInt;
import gnu.trove.iterator.TIntFloatIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.map.TIntFloatMap;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntFloatHashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
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
import java.util.logging.Level;
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
(HPL-2012-40)
 
* 
 Note that the min-cost is the optimized goal, so if a 
 graph is given that has a possible solution with more 
 matches than the number of min-cost matches, this 
 * solver will not completely match the nodes in that
 * graph
 * (it is not a maximum cardinality, min-cost matching).
 * (To preprocess the data for larger cardinality, one 
 * might want ro consider pre-matching those edges
 * that must be matched inspite of cost because no other
 * edges are connected to certain nodes else would be unmatched.
 * might consider the very fast state management
 * present in the conflict analysis of the fastest boolean
 * sat solvers to design a fast pre-filter).
 * 
 * A strength of the algorithm to note is that it does not
 * artificially double the graph in order to handle
 * unequally sized left and right sets (that is, it does not use
 * Bipartite double cover).
 </pre>
 * 
 * Note, the class is not thread safe.  For example, the instance
 * of NetworkFlow passed to methods should only be used by
 * that method and thread at the given time.
 * 
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
        public void resetExceptData() {
            m = 0;
            pathPredecessor = null;
            topPredecessor = null;
            setKey(DoubleLinkedCircularList.noValue);
            setRight(null);
            setLeft(null);
            setParent(null);
            setChildren(null);
            setMark(false);
            setNumberOfChildren(0);
        }
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
    public TIntIntMap flowAssign(Graph g) {

        validateGraph(g);
                
        if (g.getNLeft() == 1 && g.getNRight() == 1) {
            return singleNodesSolution(g);
        }
        
        // a first guess at the maximal matching size
        int sz = Math.min(g.getNLeft(), g.getNRight());

long t0 = System.currentTimeMillis();

        // hopcroft-karp produces a maximal matching of nodes
        // without using edge weights, just uses connectivity
        TIntIntMap m = hopcroftKarp(g, sz);

long t1 = System.currentTimeMillis();
long tSec = (t1 - t0)/1000;
System.out.println(tSec + " sec for hopcroftkarp");

        // if the code was given a graph created without 
        // using source and sink, need to transform that here.
        if (g.getSourceNode() == -1) {
            g = g.copyToCreateSourceSink();
        }
        
        GraphUtil.condenseEdgeWeights(g);
        
        FlowNetwork gFlow = new FlowNetwork(g, m);
        //assert(gFlow.printFlowValueIncludingSrcSnk(m.size()));
        
        log.info("init FlowNetwork:");
        //gFlow.printNetCosts();
        
        //TODO: estimate of eps may need to be revised
        /*
        the discussion of setting eps seems to be for a max-cost problem.
        
        for the numerical resolution that enables the smallest cost key
        to have math.ceil(cost/eps) = 1, eps at that tme must be
        equal to or less than the minimum cost key.  
        since eps is divided by q before first use,
        eps = q * minCost.
        The problem is then that lambda as the key range and hence size
        of minHeap in buildForest2 would need to be larger by a factor
        of eps_from_max_cost/eps_from_min_cost.  Instead of doing that,
        could continue to use the eps as defined by max cost, but
        further sort the items within a bucket by the eps derived from
        the minimum cost.
           SO, need a double key minHeap insert, and a extractMin by order
           from a bucket.
           (the eps_from_max_cost is the eps suggested by the authors,
           and that would be the main key + 1 to use to locate the
           minHeap bucket.  Then order within the bucket needs to
           be determined by the eps_from_min_cost.
           SO, need to think of how to scale the 2nd key so that can use
           the "coutning sort" O(1) performance.
           MLB implementations exist...
        
        another caveat is that the top=down assertions of "eprs-proper" from large eps to 
        smaller would possibly fail if eps were set with the minimum cost of
        the flow.
        */
        
        int s = m.size();

        // q >= 2.  consider q = O(log_2(n)) pg 41 par 3.
        // the 2nd, shorter paper states:
        // q > 1 ,  "q=8 or q=16" might be good choices
        int q = 8;//1 + (int)Math.floor(Math.log(s)/Math.log(2));
        
        // pg 32, the weight scaling techinique
        // starts w/ eps ~ maxC and reduces to ~1/s or 1/(6*s)
        // w/ nIter ~ log_q(s*maxC)
        // pg 44
        // eps_up = q^(e_up) where eps_up is smallest power of
        //    q that exceeds maxC
        // e_up * math.log(q)
                
         // expected number of iterations without a constant factor
        
        int e_up_large = 1 + (int)Math.floor(Math.log(gFlow.getMaxC())/Math.log(q));
        double eps_up_large = Math.pow(q, e_up_large - 1);
        float epsLarge = 1.f + (int)Math.floor(eps_up_large);
        
        int e_up = 1 + (int)Math.floor(Math.log(gFlow.getMinC())/Math.log(q));
        double eps_up = Math.pow(q, e_up - 1);
        int rIter = (int)(Math.log(s * gFlow.getMinC())/Math.log(q));
        if (rIter == 0) {
            rIter = 1;
        }
        
        int e_down = -(1 + (int)Math.floor( Math.log(s + 2)/Math.log(q)));
        double eps_down = Math.pow(q, e_down);

        float eps = 1.f + (int)Math.floor(eps_up);
                
        int nIterR = 0;
        
        log.info("eps=" + eps + " epsLarge=" + epsLarge
            + " rIter=" + rIter 
            + " minC=" + gFlow.getMinC() + " maxC=" + gFlow.getMaxC() 
            + " eps_down=" + eps_down);

        // all nodes V in gFlow have prices = 0
        
t0 = System.currentTimeMillis();
 
        gFlow.createPathNodes();
 
t1 = System.currentTimeMillis();
tSec = (t1 - t0);
System.out.println(tSec + " msec for createPathNodes");
       
        
        while ((epsLarge > eps_down) && (nIterR < 2*rIter)) {
            
            log.info("nIterR=" + nIterR + " s=" + s + " eps=" + eps
                + " epsLarge=" + epsLarge);
            
            // pg 44, assertions I1, I2, I3, and I4
            assert(gFlow.assertFlowValue(s));
            assert(gFlow.assertPricesAreQuantizedEps(eps));
            // if using a smaller eps than maxC, cannot assert
            //    these on first round
            if ((nIterR > 0) || (eps > gFlow.getMaxC())) {
                //assert(gFlow.integralFlowIsEpsProper(eps));
                assert(gFlow.integralBipartiteFlowIsEpsProper(eps));
                assert(gFlow.assertSaturatedBipartiteIsEpsSnug(eps));
            }
            
            eps /= ((float) q);
            epsLarge /= ((float)q);

            int ext = refine(gFlow, s, eps, epsLarge, q);

            if (ext > 0) {
long t2 = System.currentTimeMillis();
                m = gFlow.extractMatches();
long t3 = System.currentTimeMillis();
                roundFinalPrices(gFlow, eps_down);
long t4 = System.currentTimeMillis();
System.out.println((t4 - t3) + " for mmsec for roundFinalPrices");
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
        float epsLarge, int q) {
        
        log.fine("at start of refine, s=" + s + " eps=" + eps
            + " q=" + q);
        
        //assert(gFlow.printFlowValueIncludingSrcSnk(s));
        
        // S = left nodes matched in gFlow
        TIntSet surplus = new TIntHashSet();
        
        // D = right nodes matched in gFlow
        TIntSet deficit = new TIntHashSet();
       
        gFlow.zeroTheMatchedBipartiteFlow(surplus, deficit);

        // assert I4 and I5
        assert(gFlow.assertSaturatedBipartiteIsEpsSnug(epsLarge));
        
long t00 = System.currentTimeMillis();

        /*
        see Figure 7.4 on pg 53.
        raise prices so that every arc in 
        gFlow becomes eps-proper, for the resulting pseudoflow 
        f and for the new, smaller value of eps.
        this makes I3 true.
        */
        gFlow.raisePricesUntilEpsProper(epsLarge, q);

long t11 = System.currentTimeMillis();
long tSec = (t11 - t00)/1000;
System.out.println(tSec + " sec for raisePricesUntilEpsProper");

        //log.info("after raise prices, w/ eps=" + eps);
        //gFlow.printNetCosts();
        //assert(gFlow.printFlowValueIncludingSrcSnk(s));        
        
        //in [0] holds the length of the terminating deficit
        //    node which is also the forest index it was
        //    added to.
        int h = s;
        int nHIter = 0;
        
t00 = System.currentTimeMillis();

        // this should only execute sqrt(s) times
 
        while (h > 0) {
            
            log.fine("nHIter=" + nHIter + " h=" + h + " eps=" + eps
                + " epsLarge-" + epsLarge);
            
            ResidualDigraph2 rF = new ResidualDigraph2(gFlow);
long t0 = System.currentTimeMillis();
            // build a shortest-path forest from the current surpluses S, 
            // stopping when a current deficit in D is reached;
            // (pg 55)

            // pg 57, estimate lambda as size for dial arrays
            int lambda = (4*q + 4)*s/h;
            
            // map of surplus indexes and the forest index where
            //   their paths are stored
            TIntIntMap spIndexes = new TIntIntHashMap();
            
            float maxDivMin = gFlow.getMaxC()/gFlow.getMinC();
            lambda = 1 + (int)Math.floor(maxDivMin * lambda);
       
            Forest forest = 
                buildForest2(gFlow, rF, surplus, deficit, eps,
                spIndexes, lambda);
long t1 = System.currentTimeMillis();
tSec = (t1 - t0);
System.out.println(tSec + " msec for buildForest2");
           
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
            
            //assert(gFlow.printFlowValueIncludingSrcSnk(s));
                     
t0 = System.currentTimeMillis();

            log.fine("before modify prices:");
            //debug(forest);
            //gFlow.printNetCosts();
                
            List<PathsAndPrices> zeroLengthPaths = 
                modifyPathLengths(gFlow, rF, forest, 
                spIndexes, eps);
           
t1 = System.currentTimeMillis();
tSec = (t1 - t0);
System.out.println(tSec + " msec for modifyPricesAndPathLengths");

            if (zeroLengthPaths.isEmpty()) {
                log.warning("extractedPaths is empty.  h=" + h);                               
                // h > 0 so the matching is unbalanced and one-sided
                // perfect or is imperfect
                return 1;
            }
            
            /*
            NOTE: 
            paper has already modified prices by this point and then 
            asserts from pg 52:
            I1', I2, I3, I4, on FlowNetwork
            
            This method is instead, caching the price changes and 
            applying them only for the paths which are found as the
            maximal set of length 0 disjoimt vertex paths
            (then will apply assertions).
            */
                      
t0 = System.currentTimeMillis();              
            // --- Sect 8.3, create maximal set of compatible augmenting paths
            
            List<LinkedList<PathNode>> cPaths =
                findMaximalSetOfCompatiblePaths(gFlow,
                zeroLengthPaths, surplus, deficit);

t1 = System.currentTimeMillis();  
tSec = (t1 - t0);
System.out.println(tSec 
+ " msec for findMaximalSetOfCompatiblePaths");

t0 = System.currentTimeMillis();

            raisePricesForMaximalSet(gFlow, cPaths, zeroLengthPaths);

t1 = System.currentTimeMillis();  
tSec = (t1 - t0);
System.out.println(tSec 
+ " msec for raisePricesForMaximalSet");
            
            // postponed assertions:
            // assert from pg 52 
            //       I1', I2, I3, I4, on FlowNetwork
            //       and I5 on ResidualDigraph2
            assert(gFlow.assertFlowValueIncludingSrcSnk(s));
            assert(gFlow.assertPricesAreQuantizedEps(eps));
            //assert(gFlow.integralFlowIsEpsProper(epsLarge));
            //assert(gFlow.assertSaturatedBipartiteIsEpsSnug(eps));
            //assert I5: Rf has no cycles of length zero
            
            log.fine("cPaths.size=" + cPaths.size());
            if (cPaths.isEmpty()) {
                log.warning("did not find an augmenting path.  h=" + h);                               
                // h > 0 so the matching is unbalanced and one-sided
                // perfect or is imperfect
                return 1;
            }
        
t0 = System.currentTimeMillis();              
            //augment f along each of the paths in P in turn, thereby 
            //   reducing |S| = |D| = h by |P|;
            //NOTE that surplus and deficit are modified and updated
            // within augmentFlow
            augmentFlow(gFlow, cPaths);
            
t1 = System.currentTimeMillis();  
tSec = (t1 - t0);
System.out.println(tSec + " msec for augmentFlow");

            surplus.clear();
            deficit.clear();
            gFlow.getSurplusLeftIndexes(surplus);
            gFlow.getDeficitRightIndexes(deficit);

            log.fine("nSurplus=" + surplus.size() + " nDeficit="
                + deficit.size());

            log.fine("after augmentation (eps=" + eps + " epsLarge=" + epsLarge);
            //debug(forest);
            //gFlow.printNetCosts();
            
            // pg 63 assert I1', I2, I3
            assert(gFlow.assertFlowValueIncludingSrcSnk(s));
            assert(gFlow.assertPricesAreQuantizedEps(eps));
            //assert(gFlow.integralFlowIsEpsProper(epsLarge));  
            assert(gFlow.integralBipartiteFlowIsEpsProper(eps));  
            assert(gFlow.assertSaturatedBipartiteIsEpsSnug(eps));
                        
            //log.info("after augment flow:");
            //gFlow.printSaturatedLinks();

            h = surplus.size();

            // NOTE: revisit this.  adding an early exit if all of the
            // integral flow is "proper"
            if (h == 0) {
                boolean isProper = gFlow.integralFlowIsProper();
                if (isProper) {
                    return 2;
                }
            }
            
            gFlow.resetPathNodes(Long.MAX_VALUE);
            
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
        
        private TIntObjectMap<DoubleLinkedCircularList> forest =
            new TIntObjectHashMap<DoubleLinkedCircularList>();
       
        public Forest(int upperKeyLimit) {
            this.upperKeyLimit = upperKeyLimit;
        }
        
        public DoubleLinkedCircularList get(int key) {
            return forest.get(key);
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
        
        /**
         * node and it's predecessors are copied and node is inserted
         * into the forest.
         * @param node
         * @param lastKey
         * @return 
         */
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
                                
                DoubleLinkedCircularList list = forest.get(k);
                
                if (list == null) {
                    list = new DoubleLinkedCircularList();
                    orderedKeys.add(k);
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
     * runtime complexity is 
     * 
     * O(surplus * 
     *   (minHeap.extractMin +
     *   lg2(avg_n_connections_out_of_surplus_node * minHeap.insert (or decrKey))))
     * 
     * The minHeap extractMin could be improved with multi-level buckets.
     * 
     * @param gFlow
     * @param rF
     * @param surplus
     * @param deficit
     * @param eps
     * @terminatingKeys map of keu=surplus index, value = terminating 
     * deficit forest index
     * @return 
     */
    protected Forest buildForest2(
        final FlowNetwork gFlow, ResidualDigraph2 rF,
        TIntSet surplus, TIntSet deficit, float eps,
        TIntIntMap terminatingKeys, int lambda) {
        
        log.fine("buildForest2");

        log.fine("nSurplus=" + surplus.size() + " nDeficit="
            + deficit.size());
        
        terminatingKeys.clear();
        
long t0 = System.currentTimeMillis();

        /*
        if a surplus node has been inserted into the forest, any further
        paths with higher end path keys and with same surplus node as a 
        predecessort will be skipped because it's a longer path.
        TODO: there is overlapping logic in surplusForestKeys
        and terminatingKeys so this could be made more concise later.
        */
        TIntIntMap surplusForestKeys = new TIntIntHashMap();
        
        if (lambda < 4) {
            lambda = 4;
        }
        log.fine("buildForest2 forest length lambda is set to " + 
            lambda + " surplus.size=" + surplus.size());
        
        long lastKey = -1;
        
        // sparsely populated holder for the 
        // DoubleLinkedCircularList trees
        Forest forest = new Forest(lambda);
     
        MinHeapForRT2012 minHeap = new MinHeapForRT2012(lambda,
            rF.countOfForwardBipartiteLinks());
           
        PathNodes pathNodes = gFlow.getPathNodes();
     
        TIntObjectMap<LeftNode> leftNodes = pathNodes.getLeftNodes();
        
        TIntObjectMap<RightNode> rightNodes = pathNodes.getRightNodes();
        
        SourceNode sourceNode = pathNodes.getSourceNode();
        
        SinkNode sinkNode = pathNodes.getSinkNode();
        
long tIns = 0;
        TIntIterator iter = surplus.iterator();
        while (iter.hasNext()) {
            int sigma = iter.next();        
            LeftNode sNode = leftNodes.get(sigma);
            sNode.setKey(0);
 long ta0 = System.currentTimeMillis();
            minHeap.insert(sNode);
  tIns += (System.currentTimeMillis() - ta0);
        }
       
long t1 = System.currentTimeMillis();
System.out.println(tIns + " for msec minHeap.insert in init");   
long tSec = (t1 - t0);
System.out.println(tSec + " msec for init in buildForest2 heap.size=" +
minHeap.getNumberOfNodes());

        // V * minHeap,extractMin +
        //     V * n_edges_per_V * minHeap.insert (or decr Key)
        do {
            PathNode node1 = minHeap.extractMin();
            if (node1 == null) {
                //TODO: revisit this for whether additional handling of
                //   lastKey is needed
                break;
            }
            Integer index1 = (Integer)node1.getData();
                                
            if (terminatingKeys.containsKey(index1.intValue())) {
                continue;
            } else if (
                node1.pathPredecessor != null &&
                node1.pathPredecessor instanceof RightNode &&
                terminatingKeys.containsKey(
                (Integer)node1.pathPredecessor.getData()
                )) {
                continue;
            } else if (node1.pathPredecessor != null) {
                if (node1.pathPredecessor instanceof LeftNode &&
                surplusForestKeys.containsKey(
                (Integer)node1.pathPredecessor.getData()
                )) {
                    continue;
                } else if (node1.pathPredecessor.pathPredecessor 
                instanceof LeftNode && surplusForestKeys.containsKey(
                (Integer)node1.pathPredecessor.pathPredecessor.getData()
                )) {
                    continue;
                }
            }

            log.fine("extractMin = " + node1.toString());
        
            if (node1.getKey() > lambda) {
                // since this path did not reach a deficit
                // should not store the last key?
                //terminatingKeys.put(index1, Integer.valueOf((int) lastKey));
                break;
            }
            
            /*
            TODO: for first invocation within refine and first invocation
            of refine, can exclude the source and sink arc inserts here.
            */
            
            int idx1 = index1.intValue();
             
            boolean node1IsLeft = (node1 instanceof LeftNode);
            boolean nodeIsSource = (node1 instanceof SourceNode);
            boolean nodeIsSink = (node1 instanceof SinkNode);
            
            //scan all links leaving node1:
            if (nodeIsSource) {
                // scan forward source links
                TIntIterator iter2 = rF.getForwardLinksSourceRM().iterator();
                while (iter2.hasNext()) {
                    int idx2 = iter2.next();
                    if (surplusForestKeys.containsKey(idx2)) {
                        continue;
                    }
                    handlePlusLink(minHeap, node1, 
                        leftNodes.get(idx2),
                        gFlow.calcSourceNetCost(idx2), lambda, eps); 
                }
            } else if (nodeIsSink) {
                // scan backward sink links
                TIntIterator iter2 = rF.getBackwardLinksSinkRM().iterator();
                while (iter2.hasNext()) {
                    int idx2 = iter2.next();                
                    handleMinusLink(minHeap, node1, 
                        rightNodes.get(idx2), 
                        gFlow.calcSinkNetCost(idx2), 
                        lambda, eps);
                }
            } else if (node1IsLeft) {
                // scan the bipartite arcs forward
                TIntSet indexes2 = rF.getForwardLinksRM().get(idx1);
                if (indexes2 != null) {
                    TIntIterator iter2 = indexes2.iterator();
                    while (iter2.hasNext()) {
                        int idx2 = iter2.next();                    
                        float cp = gFlow.calcNetCost(idx1, idx2);
                        
                        handlePlusLink(minHeap, node1, rightNodes.get(idx2), 
                            cp, lambda, eps);
                    }
                }
           
                // if there's a source link
                if (rF.getBackwardLinksSourceRM().contains(idx1)) {
                    // insert a copy of the source node
                    PathNode sNode2 = sourceNode.copy();
                    handleMinusLink(minHeap, node1, sNode2, 
                        gFlow.calcSourceNetCost(idx1),
                        lambda, eps);
                }
            } else {
                // node1 is a RightNode
                if (rF.getBackwardLinksRM().containsKey(idx1)) {
                    int idx2 = rF.getBackwardLinksRM().get(idx1);
                    handleMinusLink(minHeap, node1, 
                       leftNodes.get(idx2), 
                       gFlow.calcNetCost(idx2, idx1),
                       lambda, eps);                     
                }
                // if there is a sink link
                if (rF.getForwardLinksSinkRM().contains(idx1)) {
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

            if (lastKey > 0) {
                // if a single shortest path root node is being inserted:
                if ((node1.pathPredecessor != null) &&
                    (node1.pathPredecessor instanceof LeftNode) && 
                    surplus.contains(
                        (Integer)node1.pathPredecessor.getData())) {
                    surplusForestKeys.put(
                        (Integer)node1.pathPredecessor.getData(), 
                        (int)lastKey);
                }
            }
            
            // store each shortest oath's terminating forest index
            if (deficit.contains(idx1) && !node1IsLeft) {
                
                terminatingKeys.put(idx1, (int)lastKey);
                
                log.fine("terminatingKeys.size=" + terminatingKeys.size()
                   + " deficit.size=" + deficit.size());
                
                if (terminatingKeys.size() == deficit.size()) {
                    break;
                }                
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
    protected TIntIntMap hopcroftKarp(Graph g, int s) {
                
        if (true) { 
            //runtime complexity O(m * sqrt(n))
            HopcroftKarp hk = new HopcroftKarp();
            int[] matched = hk.hopcroftKarpV0(new GraphWithoutWeights(g));
            log.fine("matched=" + Arrays.toString(matched));
            TIntIntMap m = new TIntIntHashMap();
            for (int i = 0; i < matched.length; ++i) {
                int v = matched[i];
                if (v > -1) {
                    m.put(i, v);
                }
            }
            return m;
        }
       
        HopcroftKarpRT2012 hk = new HopcroftKarpRT2012();
        TIntIntMap m = hk.findMaxMatching(g, s);
       
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
        
        TIntSet ks = new TIntHashSet();
        int maxK = (int)Math.ceil(1./eps_down);
        for (int k = 0; k < maxK; ++k) {
            
            boolean keep = true;
        
            float factor = (float)eps_down * k;
            
            TIntIterator iter = gFlow.getSourceForwardArcs().iterator();
            while (iter.hasNext()) {
                int xIdx = iter.next();
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
                ks.add(k);
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
        
        int kV = ks.iterator().next();
        gFlow.addToAllPrices((float)(kV * eps_down));        
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
    
    /**
     * NOTE any paths in the tree with only a single link are discarded
     * @param forest
     * @param forestIdx
     * @return 
     */
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
            
            if (path.size() > 1) {
                Misc.<PathNode>reverse(path);
                pathLists.add(path);
            }
            ++branchIdx;
        }
        
        return pathLists;
    }
    
    private class PathsAndPrices {
        List<PathNode> path;
        TIntFloatMap leftPriceIncreases;
        TIntFloatMap rightPriceIncreases;
    }

    private List<PathsAndPrices> modifyPathLengths(
        FlowNetwork gFlow, ResidualDigraph2 rF,
        Forest forest, TIntIntMap spIndexes, float eps) {
 
        /*
        NOTE:
        In this method, the links are reduced to 0 along terminating
        paths and thw implied price changes are cached until later.
        When the maximal set of vertex disjoint paths are chosen,
        will apply their cached price change to the FlowNetwork gFlow.
        
        Need to return these data in a list:
            List<PathNode> paths            
            Map<Integer, Float> leftPriceIncreases
            Map<Integer, Float> rightPriceIncreases
        */
        
        List<PathsAndPrices> output = new ArrayList<PathsAndPrices>();
        
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
          has a key and a predecessor which is a maiden node
          or that predecessor has a predecessor which is, etc.
             for example Right key=1 has predecessor Left key=0.
             the Right node is extracted from the forest and the
             remaining path is derived from the recurence of
             predecessor nodes.  Those are reversed.
             The resulting path is then 
                 (Left w/ key=0 then Right w/ key=1.
        
          more specifically,
             after the first round of buildForest2, the links
             inserted into the forest at forest[1] are
             "idle" links from surplus to deficit nodes.
             The inserted Right node w/ key=1 has a 
             predecessor Left node w/ key=0 and is extracted and reversed:
                Lft(0) -> Rgt(1) (is "idle" in this example)  
        */
        
        List<Integer> forestIndexes = forest.getKeys();
        int n0 = forestIndexes.size();
        
        if (n0 == 0 || spIndexes.isEmpty()) {
            return output;
        }

        TIntSet spForestIndexes = new TIntHashSet(spIndexes.values());
        TIntIterator iter = spForestIndexes.iterator();
        while (iter.hasNext()) {
            int ltForestKey = iter.next();
            long lt = forest.get(ltForestKey).getSentinel()
                .getRight().getKey();
           
            List<List<PathNode>> pathLists = extractPathNodes(forest, 
                ltForestKey);
            
            assert(pathLists != null);
                        
            for (int pathIdx = 0; pathIdx < pathLists.size(); ++pathIdx) {
            
                List<PathNode> path = pathLists.get(pathIdx);

                log.fine("path " + pathIdx + " (lt=" + lt + ")");
                if (log.isLoggable(Level.FINE)) {                
                    debugPath(path);
                }
                
                // these get stored for every path, but not applied to
                // the gFlow at this time
                TIntFloatMap leftPriceIncreases = new TIntFloatHashMap();
            
                TIntFloatMap rightPriceIncreases = new TIntFloatHashMap();
            
                for (int i = 0; i < (path.size() - 1); ++i) {
                    PathNode node1 = path.get(i);
                    PathNode node2 = path.get(i + 1);
                    int l1 = (int) node1.getKey();
                    int l2 = (int) node2.getKey();
                    int idx1 = ((Integer) node1.getData()).intValue();
                    int idx2 = ((Integer) node2.getData()).intValue();
                    
                    boolean node1IsRight = (node1 instanceof RightNode);
                    boolean node1IsLeft = (node1 instanceof LeftNode);
                    boolean node1IsSource = (node1 instanceof SourceNode);
                    boolean node1IsSink = (node1 instanceof SinkNode);
                    boolean node2IsRight = (node2 instanceof RightNode);
                    boolean node2IsLeft = (node2 instanceof LeftNode);
                    boolean node2IsSource = (node2 instanceof SourceNode);
                    boolean node2IsSink = (node2 instanceof SinkNode);
                    long delta2 = (lt - l2);
                    long delta1 = (lt - l1);
                    
                    // determine if node1 -> node2 is "idle" or "saturated"

                    if (node1IsRight) {
                        assert(!node2IsRight);
                        if (node2IsLeft) {
                            if (!leftPriceIncreases.containsKey(idx2) && delta2 > 0) {
                                float p2 = (delta2 * eps);
                                leftPriceIncreases.put(idx2, p2);
                            }
                            if (!rightPriceIncreases.containsKey(idx1) && delta1 > 0) {
                                float p1 = (delta1 * eps);
                                rightPriceIncreases.put(idx1, p1);
                            }
                            
                            if (rF.getBackwardLinksRM().containsKey(idx1)
                                && rF.getBackwardLinksRM().get(idx1)
                                == idx2) {

                                // "saturated" node2(left) <-- node1(right)
                                l2 -= (lt - l1);
                                
                                if (l2 < 0) {
                                    //X=N2  Y=N1
                                    log.severe("forest at time of error:");
                                    debug(forest);
                                    throw new IllegalStateException(
                                        "for saturated link X=" + node2.getData() + " to Y="
                                        + node1.getData() + " results in l(X)=" + l2);
                                }
                                
                                node2.setKey(l2);
                            } else {
                                assert (rF.getForwardLinksRM().get(idx2) != null
                                && rF.getForwardLinksRM().get(idx2)
                                .contains(idx1));

                                // "idle" node2(left) --> node1(right) 
                                l1 -= (lt - l2);
                                
                                if (l1 < 0) {
                                    //X=N2 Y=N1
                                    log.severe("forest at time of error:");
                                    debug(forest);
                                    throw new IllegalStateException(
                                        "for idle link X=" + node2.getData() + " to Y="
                                        + node1.getData() + " results in l(Y)=" + l1);
                                }
                                
                                node1.setKey(l1);
                            }
                        } else if (node2IsSink) {
                            // node2 is sink so add to existing price changes?
                            
                            if (!rightPriceIncreases.containsKey(idx1) && delta1 > 0) {
                                float p1 = (delta1 * eps);
                                rightPriceIncreases.put(idx1, p1);
                            }
                            if (rF.getBackwardLinksSinkRM().contains(idx1)) {
                                // "saturated" node1(right) <-- node2(sink)
                                l1 -= (lt - l2);
                                
                                if (l1 < 0) {
                                    // Y=N1  Y=N2=SINK
                                    log.severe("forest at time of error:");
                                    debug(forest);
                                    throw new IllegalStateException(
                                        "for saturated link Y=" + node1.getData() + " to sink="
                                        + node2.getData() + " results in l(Y)=" + l1);
                                }
                                
                                node1.setKey(l1);
                            } else {
                                assert (rF.getForwardLinksSinkRM().contains(idx1));
                                // "idle" node1(right) --> node2(sink)
                                l2 -= (lt - l1);
                                
                                if (l2 < 0) {
                                    // Y=N1  Y=N2=SINK
                                    log.severe("forest at time of error:");
                                    debug(forest);
                                    throw new IllegalStateException(
                                        "for idle link Y=" + node1.getData() + " to sink="
                                        + node2.getData() + " results in l(sink)=" + l2);
                                }
                                
                                node2.setKey(l2);
                            }
                        }
                        // end of node1 is right node
                    } else if (node1IsLeft) {
                        assert(!node2IsLeft);
                        if (node2IsRight) {   
    
                            if (!rightPriceIncreases.containsKey(idx2) && delta2 > 0) {
                                float p2 = (delta2 * eps);
                                rightPriceIncreases.put(idx2, p2);
                            }
                            if (!leftPriceIncreases.containsKey(idx1) && delta1 > 0) {
                                float p1 = (delta1 * eps);
                                leftPriceIncreases.put(idx1, p1);
                            }

                            if (rF.getBackwardLinksRM().containsKey(idx2)
                                && rF.getBackwardLinksRM().get(idx2) == idx1) {

                                // saturated link node1(left) <-- node2(right)
                                l1 -= (lt - l2);
                                
                                if (l1 < 0) {
                                    //X=N1  Y=N2
                                    log.info("forest at time of error:");
                                    debug(forest);
                                    throw new IllegalStateException(
                                        "for saturated link X=" + node1.getData() + " to Y="
                                        + node2.getData() + " results in l(X)=" + l1);
                                }
                                
                                node1.setKey(l1);
                            } else {
                                assert(rF.getForwardLinksRM().get(idx1)
                                    .contains(idx2));
                                // idle link node1(left) --> node2(right)
                                l2 -= (lt - l1);
                                
                                if (l2 < 0) {
                                    //X=N1  Y=N2
                                    log.info("forest at time of error:");
                                    debug(forest);
                                    throw new IllegalStateException(
                                        "for idle link X=" + node1.getData() + " to Y="
                                        + node2.getData() + " results in l(Y)=" + l2);
                                }
                                
                                node2.setKey(l2);
                            }

                        } else if (node2IsSource) {
                            // node2 is source so add to existing price changes?
                            if (!leftPriceIncreases.containsKey(idx1) && delta1 > 0) {
                                float p1 = (delta1 * eps);
                                leftPriceIncreases.put(idx1, p1);
                            }
                            if (rF.getBackwardLinksSourceRM().contains(idx1)) {
                                // saturated  node2(source) <--- node1(Left)
                                l2 -= (lt - l1);
                            
                                if (l2 < 0) {
                                    //X=N1  N2=SOURCE
                                    log.info("forest at time of error:");
                                    debug(forest);
                                    throw new IllegalStateException(
                                        "for saturated link source=" + node2.getData() 
                                            + " to X="
                                        + node1.getData() + " results in l(source)=" + l2);
                                }
                                
                                node2.setKey(l2);
                            } else {
                                assert(rF.getForwardLinksSourceRM().contains(idx1));
                                // idle  node2(source) --> node1(left)
                                l1 -= (lt - l2);
                                
                                if (l1 < 0) {
                                    //X=N1  N2=SOURCE
                                    log.info("forest at time of error:");
                                    debug(forest);
                                    throw new IllegalStateException(
                                        "for idle link source=" + node2.getData() + " to X="
                                        + node1.getData() + " results in l(X)=" + l1);
                                }
                                
                                node1.setKey(l1);
                            }
                        }
                        // end of node1 is left node
                    } else if (node1IsSource) {
                        assert (node2IsLeft);
                        if (!leftPriceIncreases.containsKey(idx2) && delta2 > 0) {
                            float p2 = ((lt - l2) * eps);
                            leftPriceIncreases.put(idx2, p2);
                        }
                        // node1 is source so add to existing price changes?
                        if (rF.getBackwardLinksSourceRM().contains(idx2)) {
                            // saturated  node1(source) <--- node2(Left)
                            l1 -= (lt - l2);
                            
                            if (l1 < 0) {
                                //N1=SOURCE  N2=X
                                log.info("forest at time of error:");
                                debug(forest);
                                throw new IllegalStateException(
                                    "for saturated link source=" + node1.getData() + " ro X="
                                    + node2.getData() + " results in l(source)=" + l1);
                            }
                            
                            node1.setKey(l1);
                        } else {
                            assert (rF.getForwardLinksSourceRM().contains(idx2));
                            // idle  node1(source) --> node2(left)
                            l2 -= (lt - l1);
                            
                            if (l2 < 0) {
                                //N1=SOURCE  N2=X
                                log.info("forest at time of error:");
                                debug(forest);
                                throw new IllegalStateException(
                                    "for idle link source=" + node1.getData() + " ro X="
                                    + node2.getData() + " results in l(X)=" + l2);
                            }
                            
                            node2.setKey(l2);
                        }
                    } else if (node1IsSink) {
                        assert(node2IsRight);
                        if (!rightPriceIncreases.containsKey(idx2) && delta2 > 0) {
                            float p2 = (delta2 * eps);
                            rightPriceIncreases.put(idx2, p2);
                        }
                        // node1 is sink so add existing price change?
                        if (rF.getBackwardLinksSinkRM().contains(idx2)) {
                            // "saturated" node2(right) <-- node1(sink)
                            l2 -= (lt - l1);
                            
                            if (l2 < 0) {
                                //N2=Y  N1=SINK
                                log.info("forest at time of error:");
                                debug(forest);
                                throw new IllegalStateException(
                                    "for saturated link Y=" + node2.getData() 
                                    + " to sink="
                                    + node1.getData() + " results in l(Y)=" + l2);
                            }
                            
                            node2.setKey(l2);
                        } else {
                            assert(rF.getForwardLinksSinkRM().contains(idx2));
                            // "idle" node2(right) --> node1(sink)
                            l1 -= (lt - l2);
                            
                            if (l1 < 0) {
                                //N2=Y  N1=SINK
                                log.info("forest at time of error:");
                                debug(forest);
                                throw new IllegalStateException(
                                    "for saturated link Y=" + node2.getData() 
                                    + " to sink="
                                    + node1.getData() + " results in l(sink)=" + l1);
                            }
                            
                            node1.setKey(l1);
                        }
                    } 
                }
                
                // if path lengths are zero,
                // store in output
                boolean doStore = true;
                for (PathNode node : path) {
                    int len = (int) node.getKey();
                    if (len != 0) {
                        doStore = false;
                        break;
                    }
                }
                if (doStore) {
                    
                    PathsAndPrices p = new PathsAndPrices();
                    p.path = path;
                    p.leftPriceIncreases = leftPriceIncreases;
                    p.rightPriceIncreases = rightPriceIncreases;
                    
                    output.add(p);                    
                }
            }
        }
        
        return output;
    }
    
    private List<LinkedList<PathNode>>
        findMaximalSetOfCompatiblePaths(FlowNetwork gFlow,
            List<PathsAndPrices> pathsList,
            TIntSet surplus, TIntSet deficit) {
       
        List<LinkedList<PathNode>> augPaths 
            = new ArrayList<LinkedList<PathNode>>();
        
        if (pathsList.isEmpty()) {
            return augPaths;
        }
        
        //pseudocode from pg 62 , Figure 8.2
        //   with input = length 0 adj list, the surplus
        //   list and the deficit list
        //   to create the maximal set of compatible paths.
        //
        //   state that needs to be tracked for a vertex:
        //    - visited (== marked)
        //    - identity (== Left or Right or Source or Sink)

        Map<PathNode, Set<PathNode>> pathLinksMap =
            createNewAdjacencyMap(pathsList, gFlow);
 
        //debug2(pathLinksMap);

        Set<LeftNode> surp = new HashSet<LeftNode>();
        Set<RightNode> def = new HashSet<RightNode>();
        makeSurplusAndDeficitSubSets(pathLinksMap, 
            surplus, deficit, surp, def);
        
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
    }
    
    private List<LinkedList<PathNode>> extractZeroLengthLinks(
        List<List<PathNode>> paths) {

        // each linked list is a pair of 2 path nodes whose link length = 0
        List<LinkedList<PathNode>> pathLists = new
            ArrayList<LinkedList<PathNode>>();
        
        for (int i = 0; i < paths.size(); ++i) {
            List<PathNode> path = paths.get(i);
            int n2 = path.size();
            for (int idx = 0; idx < (n2 - 1); ++idx) {
                PathNode node1 = path.get(idx);
                PathNode node2 = path.get(idx + 1);
                int l1 = (int) node1.getKey();
                int l2 = (int) node2.getKey();
                if (l2 != 0) {
                    continue;
                }
                LinkedList<PathNode> link = new LinkedList<PathNode>();
                link.add(node1);
                link.add(node2);
                pathLists.add(link);
            }
        }
        
        return pathLists;
    }
    
    private void augmentFlow(FlowNetwork gFlow, 
        List<LinkedList<PathNode>> cPaths) {

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
    }

    private Map<PathNode, Set<PathNode>> createNewAdjacencyMap(
        List<PathsAndPrices> pathsList,
        FlowNetwork gFlow) {
        
        PathNodes pathNodes = gFlow.getPathNodes();
        pathNodes.resetNodeExceptData(Long.MAX_VALUE);
   
        SinkNode sinkNode = pathNodes.getSinkNode();      
        
        SourceNode sourceNode = pathNodes.getSourceNode();
        
        Map<PathNode, Set<PathNode>> pathLinksMap 
            = new HashMap<PathNode, Set<PathNode>>();
        
        for (PathsAndPrices p : pathsList) {
            List<PathNode> path = p.path;
            for (int i = 0; i < path.size() - 1; ++i) {            
                PathNode node1 = path.get(i);
                PathNode node2 = path.get(i + 1);
            
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
                    node1 = pathNodes.getLeftNodes().get(index1);
                } else {
                    node1 = pathNodes.getRightNodes().get(index1);
                }
                if (node2 instanceof SinkNode) {
                    node2 = sinkNode;
                } else if (node2 instanceof SourceNode) {
                    node2 = sourceNode;
                } else if (node2 instanceof LeftNode) {
                    node2 = pathNodes.getLeftNodes().get(index2);
                } else {
                    node2 = pathNodes.getRightNodes().get(index2);
                }

                Set<PathNode> values = pathLinksMap.get(node1);
                if (values == null) {
                    values = new HashSet<PathNode>();
                    pathLinksMap.put(node1, values);
                }
                values.add(node2);
            }
        }

        return pathLinksMap;
    }
    
    /**
     * NOTE: the weights were condensed to make the min Heaps
     * faster, so corrections for that must be done before user
     * can use the flow network results
     * @return 
     */
    protected FlowNetwork getFinalFlowNetwork() {
        return finalFN;
    }
    
    private void makeSurplusAndDeficitSubSets(
        Map<PathNode, Set<PathNode>> pathLinksMap, 
        TIntSet surplus, TIntSet deficit, 
        Set<LeftNode> outputSurplus, 
        Set<RightNode> outputDeficit) {
  
        for (Entry<PathNode, Set<PathNode>> entry : pathLinksMap.entrySet()) {
            PathNode node1 = entry.getKey();
            Integer index1 = (Integer)node1.getData();
            int idx1 = index1.intValue();
            if (node1 instanceof LeftNode && 
                surplus.contains(idx1)) {
                outputSurplus.add((LeftNode)node1);
            } else if (node1 instanceof RightNode && 
                deficit.contains(idx1)) {
                outputDeficit.add((RightNode)node1);
            }
            for (PathNode node2 : entry.getValue()) {
                Integer index2 = (Integer)node2.getData();
                int idx2 = index2.intValue();
                if (node2 instanceof LeftNode && 
                    surplus.contains(idx2)) {
                    outputSurplus.add((LeftNode)node2);
                } else if (node2 instanceof RightNode && 
                    deficit.contains(idx2)) {
                    outputDeficit.add((RightNode)node2);
                }
            }
        }
    }
    
    private void debugPath(List<PathNode> path) {
        StringBuilder sb = new StringBuilder("path=");
        for (PathNode node : path) {
            String str = String.format("%s %s %s", 
                node.id, Long.toString(node.getKey()),
                node.getData().toString());
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
        if (g.getEdgeWeights().isEmpty()) {
            errors.append(" every vertex in g must have at least one edge.");
        } else {
            // every vertex must be part of at least one edge
            TIntSet xInEdges = new TIntHashSet();
            TIntSet yInEdges = new TIntHashSet();
            boolean lessThan1 = false;
            
            TObjectIntIterator<PairInt> iter2 
                = g.getEdgeWeights().iterator();
            for (int i = g.getEdgeWeights().size(); i-- > 0;) {
                iter2.advance();
                PairInt p = iter2.key();
                if (iter2.value() < 1) {
                    errors.append(
                    " edge wriths should be larger than 0.");
                }
                xInEdges.add(p.getX());
                yInEdges.add(p.getY());
            }
            if (!lessThan1) {
                if (g.getNLeft() > xInEdges.size() || g.getNRight() > yInEdges.size()) {
                    errors.append(" every vertex in g must have at least one edge.");
                }
            }
        }
        
        if (errors.length() > 0) {
            throw new IllegalArgumentException(errors.toString());
        }
    }

    private TIntIntMap singleNodesSolution(Graph g) {
    
        TIntIntMap m = new TIntIntHashMap();
        
        if (g.getEdgeWeights().containsKey(new PairInt(0, 0))) {
            m.put(0, 0);
        }
        
        throw new IllegalStateException("graph only had 1 left"
            + " and right node, but no edge between them.");
    }
  
    private void debug2(Map<PathNode, Set<PathNode>> pathLinksMap) {
        
        StringBuilder sb = new StringBuilder();
        
        for (Entry<PathNode, Set<PathNode>> entry : pathLinksMap.entrySet()) {
            
            PathNode node1 = entry.getKey();
            
            sb.append(String.format("%s %s", node1.id, node1.getData().toString()));
            sb.append("-->[");
            
            Set<PathNode> set = entry.getValue();
            for (PathNode node2 : set) {
                sb.append(String.format("%s %s", node2.id, 
                    node2.getData().toString()));
                sb.append(", ");
            }
            sb.append("]\n");
        }
        log.info(sb.toString());
    }

    private void raisePricesForMaximalSet(FlowNetwork gFlow, 
        List<LinkedList<PathNode>> cPaths, 
        List<PathsAndPrices> cachedPriceChanges) {
        
        /*
        Need a lookup structure to find the
        price maps in cachedPriceChanges for a path/
        
        The individual path lengths will be at most 4 (min cost edge
        plus possible zero length edge of sink and source arcs.
        
        identity for a lookup is class extending PathNode and the index.
        since sink and source indexes are unique within right and left,
        respectively, can simplify the key for a node being left or right and the index.
       
        The number of the vertices in the graph is assumed to be much less
        than (1<<31) - 1.
       
        to make a key Will use the index and make it negative sign for
        a RightNode.
        TODO: consider compressed keys.
        */
    
        TObjectIntMap<QuadInt> priceLookupMap = createPriceLookupMap(cachedPriceChanges);
        
        for (LinkedList path : cPaths) {
            
            Misc.reverse(path);
            
            QuadInt key = createKey(path);
                        
            if (!priceLookupMap.containsKey(key)) {
                throw new IllegalStateException("error in algorithm."
                    + "  partial paths need to be considered");
            }
            
            int idx = priceLookupMap.get(key);
            
            PathsAndPrices p = cachedPriceChanges.get(idx);
            
            TIntFloatMap leftP = p.leftPriceIncreases;
            
            TIntFloatIterator iter = leftP.iterator();
            for (int i = leftP.size(); i-- > 0;) {
                iter.advance();            
                float f = iter.value();
                gFlow.addToLeftPrice(iter.key(), f);
            }
            
            TIntFloatMap rightP = p.rightPriceIncreases;
            iter = rightP.iterator();
            for (int i = rightP.size(); i-- > 0;) {
                iter.advance();                
                float f = iter.value();
                gFlow.addToRightPrice(iter.key(), f);
            }
        }
        
    }
    
    private QuadInt createKey(List<PathNode> path) {
        
        QuadInt q = new QuadInt();
        
        for (int i = 0; i < path.size(); ++i) {
            PathNode node = path.get(i);
            int idx = ((Integer)node.getData()).intValue();
            if (node instanceof RightNode || node instanceof SinkNode) {
                idx *= -1;
            }
            switch(i) {
                case 0:
                    q.setA(idx);
                    break;
                case 1:
                    q.setB(idx);
                    break;
                case 2:
                    q.setC(idx);
                    break;
                case 3:
                    q.setD(idx);
                    break;
                default:
                    throw new IllegalStateException("Error in algorithm"
                        + ".  Expecting path length < 5");
            }
        }
        
        return q;
    }
    
    private TObjectIntMap<QuadInt> createPriceLookupMap(
        List<PathsAndPrices> cachedPriceChanges) {
        
        TObjectIntMap<QuadInt> map = new TObjectIntHashMap<QuadInt>();
        
        for (int i = 0; i < cachedPriceChanges.size(); ++i) {
            
            PathsAndPrices p = cachedPriceChanges.get(i);
            
            map.put(createKey(p.path), i);
        }
        
        return map;
    }

}
