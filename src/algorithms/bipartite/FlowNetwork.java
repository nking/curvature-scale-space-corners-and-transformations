package algorithms.bipartite;

import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntFloatHashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.logging.Logger;

/**
 * class representing a flow network for 
 * MinCostUnbalancedAssignment.java
   If the sink and source nodes and weights are set in the
   origin graph, sink and source logic is used in code.
 * 
 * @author nichole
 */
public class FlowNetwork {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    /*
     matchings M in Graph g are integral flows f in the 
     FlowNetwork.
        
     If allows for conditions which are not integral, one uses
     linear programming terms to replace the term slackness
     with "proper"ness on properties f and p where p is prices.
     - a node has a per unit price.
     - adopting the cost to dispose of a unit as the model p_d(v)
     - (the dispose cost is equal with opposite sign to acquire cost)
     - the cost of each arc is the difference between the
     disposal cost of the two encasing nodes.
     - the net cost is cp(v, w)
     - the benefits of the arc is bp(v,w) and the sum of it
     and cp(v,w) is 0.
     cp(v, w) := c(v, w) + pa(v) − pa(w) 
     = c(v, w) − pd(v) + pd(w) 
     bp(v, w) := b(v, w) − pa(v) + pa(w) 
     = b(v, w) + pd(v) − pd(w)
     
     cp(v, w) = c(v, w) − pd(v) + pd(w)

     - cost and net flow are related as:
     cp(f) = summation over arcs in flow network
     f(v,w) * cp(v,w)
     = summation over arcs in flow network
     f(v,w) * (c(v,w) − pd(v) + pd(w))
     = c(f) − |f| * pd(⊢) − pd(⊣)
     note that the value |f| is the total flow out of the
     source ⊢ and into the sink ⊣, 
     while flow is conserved at all other nodes.

     - each arc w/ f(v,w) = 0 and cp(x, y) ≥ 0
     and is idle
     and each arc w/ f(v,w) = 1 and cp(x, y) ≤ 0
     is saturated.
     (the pair of f and p are proper under those conditions)
     - arcs with fractional flow have zero net cost
     - (edges w/ 0 net cost are "tight") 
     */

    /**
     * <pre>
     * |-
     * </pre>
     */
    private final int sourceNode;
    
    /**
     * <pre>
     * -|
     * </pre>
     */
    private int sinkNode;

    /**
     * left (==X) vertices in graph
     */
    private final int nLeft;

    /**
     * right (==Y) vertices in graph
     */
    private final int nRight;

    /**
     * prices for left (author's adopted price to dispose of product).
     * The source node is at position nLeft.
     */
    private final float[] pLeft;

    /**
     * prices for right (author's adopted price to dispose of product).
     * The sink node is at position nRight.
     */
    private final float[] pRight;

    private int maxC = Integer.MIN_VALUE;
    
    private int minC = Integer.MAX_VALUE;

    /** forward arcs only in the flow graph
    may represent these differently soon.
    they are also known as "bipartite arcs".
    */
    private final TIntObjectHashMap<TIntSet> forwardArcs
        = new TIntObjectHashMap<TIntSet>();

    /** forward arcs from source to the left (a.k.a. X nodes).
     * If an arc is saturated, it will have flux=1 along this
     * arc.
    */
    private final TIntSet sourceForwardArcs = 
        new TIntHashSet();

    /** forward arcs from the right (a.k.a. Y nodes) to sink node
    If an arc is saturated, it will have flux=1 along this arc.
    */
    private final TIntSet sinkForwardArcs =
        new TIntHashSet();

    /**
     * per unit flow costs are the same as the edges given in G. 
     * key = index in left (==X==v) and index in right (==Y==w), 
     * value = cost of the arc. might
     * change to use sparse matrix format
     */
    private final TObjectIntHashMap<PairInt> c
        = new TObjectIntHashMap<PairInt>();

    /**
     * the flow on the arc. key = index in left (==X==v) and index in right
     * (==Y==w), value = number from 0 to 1 inclusive if it's a "pseudoflow":
     * the value 0 is "idle" and corresponds to an unmatched link in the
     * residual graph. the value 1 is "saturated" and corresponds to a matched
     * link in the residual graph (but in the residual graph, it would be
     * specified in format right to left). might change to use sparse matrix
     * format
     */
    private final TObjectFloatHashMap<PairInt> f
        = new TObjectFloatHashMap<PairInt>();
    
    /**
     * per unit dummy flow costs from source to X nodes
     * (initialized as 0).
     * key = index in left (==X==v) 
     * value = cost of the arc. might
     * change to use sparse matrix format
     */
    private final TIntIntHashMap sourceToLeftC =
        new TIntIntHashMap();
   
    /**
     * per unit dummy flow costs from sink to X nodes
     * (initialized as 0).
     * key = index of right node (== Y)
     * value = cost of the arc
     */
    private final TIntIntHashMap rightToSinkC =
        new TIntIntHashMap();

    /**
     * flux from source to x node
     */
    private final TIntFloatHashMap sourceToLeftF =
        new TIntFloatHashMap();
    
    /**
     * flux from y node to sink
     */
    private final TIntFloatHashMap rightToSinkF =
        new TIntFloatHashMap();
    
    // nodes placed here for convenience.  they're an adaptation of
    // the graph used in path algorithms.
    private PathNodes pathNodes = null;
    
    public FlowNetwork(Graph g, TIntIntMap m) {

        this.nLeft = g.getNLeft();
        this.nRight = g.getNRight();
        this.sourceNode = nLeft;
        this.sinkNode = nRight;
        pLeft = new float[nLeft + 1];
        pRight = new float[nRight + 1];
                
        TObjectIntIterator<PairInt> iter = g.getEdgeWeights().iterator();
        for (int i = g.getEdgeWeights().size(); i-- > 0;) {
            
            iter.advance();
            
            PairInt p = iter.key();

            int idx1 = p.getX();
            int idx2 = p.getY();

            TIntSet indexes2 = forwardArcs.get(idx1);
            if (indexes2 == null) {
                indexes2 = new TIntHashSet();
                forwardArcs.put(idx1, indexes2);
            }
            indexes2.add(idx2);

            int cost = g.getEdgeWeights().get(p);
            c.put(p, cost);
            f.put(p, 0);

            if (Math.abs(cost) > maxC) {
                maxC = Math.abs(cost);
            }
            if (Math.abs(cost) < minC) {
                minC = Math.abs(cost);
            }

            sourceForwardArcs.add(idx1);
            sinkForwardArcs.add(idx2);
            
            // set source and sink arcs to 0 and overwrite
            // those for matched vertexes in next block
            sourceToLeftF.put(idx1, 0);
            rightToSinkF.put(idx2, 0);
            
            // cost of dummy edge is 0
            sourceToLeftC.put(idx1, 0);
            rightToSinkC.put(idx2, 0);
        }
        
        TIntIntIterator iter2 = m.iterator();
        for (int ii = m.size(); ii-- > 0;) {
            iter2.advance();
            int leftIdx = iter2.key();
            int rightIdx = iter2.value();
            sourceToLeftF.put(leftIdx, 1);
            PairInt p = new PairInt(leftIdx, rightIdx);
            f.put(p, 1);
            rightToSinkF.put(rightIdx, 1);
        }

        assert (maxC > 1);

        /*
        For each vertex x in X, there is a left-dummy arc |- -> x, 
           directed from the source node |- to the node x. 
           The per-unit cost of a left-dummy arc is zero: 
           c(|-, x) := 0. 
        For each vertex y in Y , there is a right-dummy arc y -> -|, 
           directed from the node 
           y to the sink node -| and of cost zero: 
           c(y, -|) := 0.            
        */
        
    }
    
    public int getMaxC() {
        return maxC;
    }
    
    public int getMinC() {
        return minC;
    }

    /**
     * <pre>
     * cp(f) = c(f) - |f|*(pd(|-) - pd(-|))
     * </pre>
     * @return 
     */
    public double calcNetCost() {
        
        double cf = calcFluxCost();
        double flux = calcTotalFlow();
        
        double cp = cf - (flux * (pLeft[sourceNode] - pRight[sinkNode]));
        
        return cp;
    }
    
    /**
     * total flow out of the source (or equivalently,
     * into the sink or equivalently out of all bipartite arcs).
     * @return 
     */
    public double calcTotalFlow() {
        
        // excludes source and sink nodes too
        
        double sum = 0;
        
        TIntObjectIterator<TIntSet> iter3 = forwardArcs.iterator();
        for (int ii = forwardArcs.size(); ii-- > 0;) {
            iter3.advance();
            int idx1 = iter3.key();
            TIntSet idxs2 = iter3.value();
            
            TIntIterator iter4 = idxs2.iterator();
            while (iter4.hasNext()) {
                int idx2 = iter4.next();
                PairInt p = new PairInt(idx1, idx2);
                float unitFlow = f.get(p);
                sum += unitFlow;
            }
        }
        
        return sum;
    }
    
    /**
     * c(f) = summation over arcs of (f(v,w) * c(v,w))
     * @return 
     */
    public double calcFluxCost() {
        
        // includes source and sink nodes too, but their costs are 0
        
        double sum = 0;
        
        TIntObjectIterator<TIntSet> iter = forwardArcs.iterator();
        for (int i = forwardArcs.size(); i-- > 0;) {
            
            iter.advance();
            int idx1 = iter.key();
            TIntIterator iter2 = iter.value().iterator();
            while (iter2.hasNext()) {
                int idx2 = iter2.next();
                PairInt p = new PairInt(idx1, idx2);
                float unitFlow = f.get(p);
                float cost = c.get(p);
                sum += (unitFlow * cost);
            }
        }
        
        TIntIterator iter2 = sourceForwardArcs.iterator();
        while (iter2.hasNext()) {
            int idx1 = iter2.next();        
            float unitFlow = sourceToLeftF.get(idx1);
            float cost = sourceToLeftC.get(idx1);
            sum += (unitFlow * cost);
        }
        
        iter2 = sinkForwardArcs.iterator();
        while (iter2.hasNext()) {
            int idx1 = iter2.next();
            float unitFlow = rightToSinkF.get(idx1);
            float cost = rightToSinkC.get(idx1);
            sum += (unitFlow * cost);
        }
        
        return sum;
    }
    
    public float getSourceToLeftFlow(int idx) {
        return sourceToLeftF.get(idx);
    }
    public float getRightToSinkFlow(int idx) {
        return rightToSinkF.get(idx);
    }
    public int getSourceToLeftCost(int idx) {
        return sourceToLeftC.get(idx);
    }
    
    public float calcNetCost(int u, int v) {

        if (u == sourceNode) {
            return calcSourceNetCost(v);
        } else if (v == sinkNode) {
            return calcSinkNetCost(u);
        } 
       
        PairInt p = new PairInt(u, v);
        int cost = c.get(p);
        float pdX = pLeft[u];
        float pdY = pRight[v];
        
        float cp = cost - pdX + pdY;

        return cp;
    }
    
    public float getFlow(int u, int v) {

        if (u == sourceNode) {
            return sourceToLeftF.get(v);
        } else if (v == sinkNode) {
            return rightToSinkF.get(u);
        } 
        
        return f.get(new PairInt(u, v));
    }
    
    protected float calcSourceNetCost(int v) {

        int cost = sourceToLeftC.get(v);
        float pdX = pLeft[sourceNode];
        float pdY = pLeft[v];
       
        float cp = cost - pdX + pdY;

        return cp;
    }

    protected float calcSinkNetCost(int u) {

        int cost = rightToSinkC.get(u);
        float pdX = pRight[u];
        float pdY = pRight[sinkNode];
       
        float cp = cost - pdX + pdY;

        return cp;
    }

    public float calcNetCost(PairInt p) {
        return calcNetCost(p.getX(), p.getY());
    }
    
    /**
     * assertion, pg 44, I1.
     * The flux f on NG is a flow of value |f|=s
     * @param s
     * @return 
     */
    boolean assertFlowValue(int s) {
                
        //Should this be for the integral flow only?
        double flow = calcTotalFlow();
        
        log.info("assert I1: flow=" + flow + " s=" + s);
        
        return (Math.abs(flow - s) < 1);
    }
    
    boolean printFlowValueIncludingSrcSnk(int s) {
        
        log.info("s=" + s);
        
        double flow = 0;
        
        flow += calcTotalFlow();
        
        log.info("bipartite flow sum=" + flow);
        
        flow += calcTotalSourceFlow();
        
        log.info("bipartite + source flow sum=" + flow);
        
        flow += calcTotalSinkFlow();
        
        log.info("bipartite + source + sink flow sum=" + flow);
    
        return true;
    }
    
    /**
     * assert pg 52, I1'.
     * assert that flux f on NG is a flow of value |f|=s
     * @param nMatchingsHK the number of matched
     * nodes from hopcroft-karp before refine is invoked.
     * @return the number of surplus nodes which is the same
     * as the number of deficit nodes.  this may be smaller
     * than nMatchingsHK after an iteration of refine
     */
    boolean assertFlowValueIncludingSrcSnk(int nMatchingsHK) {
                
        double flow = 0;
        
        // NOTE: not completely sure about the intended assertion,
        // but I'm assuming that in refine, after zeroing
        // the bipartite matched flow then raising prices
        // and changing arc lengths, and then augemtning paths,
        // one has a matching
        // between the surplus and deficit nodes.
        // the matching might only be one node for example,
        // the the bipartite flow sum is 1, and the
        // that equals the number of hopcroftkarp matchings
        // minus the updated number of surplus matchings
        
        flow += calcTotalFlow();
        
        log.info("bipartite flow sum=" + flow);
        
        // nSurplus is the number of nodes, excluding the
        // sink where the flow into the node is larger
        // than the flow leaving the sink.
        
        TIntHashSet surplus = new TIntHashSet();
        getSurplusRightIndexes(surplus);
        assert(surplus.isEmpty());
        
        getSurplusLeftIndexes(surplus);
        
        int nSurplus = surplus.size();
        
        log.info("s=" + nMatchingsHK + " h=" + nSurplus);
                
        return (Math.abs(flow - (nMatchingsHK - nSurplus)) < 1);
    }
    
    public boolean printSurplusAndDeficit() {
    
        TIntHashSet surplus = new TIntHashSet();
        getSurplusLeftIndexes(surplus);
        
        TIntHashSet deficit = new TIntHashSet();
        getDeficitRightIndexes(deficit);
        
        StringBuilder sb = new StringBuilder("surplus=(");
        TIntIterator iter = surplus.iterator();
        while (iter.hasNext()) {
            int s = iter.next();
            sb.append(Integer.toString(s)).append(", ");
        }
        sb.append(")");
        log.info(sb.toString());
        
        sb = new StringBuilder("deficit=(");
        iter = deficit.iterator();
        while (iter.hasNext()) {
            int d = iter.next();
            sb.append(Integer.toString(d)).append(", ");
        }
        sb.append(")");
        log.info(sb.toString());
        
        return true;
    }
  
    double calcTotalSourceFlow() {
                
        double flow = 0;
        
        TIntIterator iter = sourceForwardArcs.iterator();
        while (iter.hasNext()) {
            int idx = iter.next();
            float unitFlow = sourceToLeftF.get(idx);
            flow += unitFlow;
        }
                
        return flow;
    }
    
    double calcTotalSinkFlow() {
                
        double flow = 0;
        
        TIntIterator iter = sinkForwardArcs.iterator();
        while (iter.hasNext()) {
            int idx = iter.next();
            float unitFlow = rightToSinkF.get(idx);
            flow += unitFlow;
        }
                
        return flow;
    }

    /**
     * assert pg 44, I2.
     * the prices at all nodes are multiples of eps
     * @param eps
     * @return 
     */
    boolean assertPricesAreQuantizedEps(float eps) {
        
        double tolerance = 0.01;
        
        for (int i = 0; i < pLeft.length; ++i) {
            float p = pLeft[i];
            float div = p/eps;
            double r = div - Math.floor(div);
            if (r > 0.05) {
                return false;
            }
        }
        
        return true;
    }
    
    /**
     * raise prices on all nodes so that every arc becomes eps-proper, 
     * for the resulting pseudoflow f
     * @param eps
     * @param q
     */
    public void raisePricesUntilEpsProper(float eps, int q) {
        
        float delta = (q - 1.f) * eps;
        
        float qEps = q * eps;
        
        log.fine("delta=" + delta + " qEps=" + qEps);
        
        //see Figure 7.4, pg 53
        // looks like the price raise is multiples of (q-1)*eps

        TIntIterator iter = sourceForwardArcs.iterator();
        // calc for source dummy arcs
        while (iter.hasNext()) {
            int xIdx = iter.next();
            float unitFlow = sourceToLeftF.get(xIdx);
            float cp = calcSourceNetCost(xIdx);
        
            if (unitFlow == 0) {
                // idle, cp > -qEps
                //cp = cost - pdX + pdY so inr pdY
                float count = 1.f;
                while (cp <= -qEps) {
                    pLeft[xIdx] += (count * delta);
                    cp = calcSourceNetCost(xIdx);
                    log.fine("raise price of left for source arc" 
                        + xIdx + " to cp=" + cp);
                }
            } else if (Math.abs(unitFlow - 1) < 0.01f) {
                // saturated
                //cp = cost - pdX + pdY so incr pdX
                float count = 1.f;
                while (cp > qEps) {
                    pLeft[sourceNode] += (count * delta);
                    cp = calcSourceNetCost(xIdx);
                    log.fine("lower price of left for source arc" 
                        + xIdx + " to cp=" + cp);
                }
            }
        }
        
        TIntObjectIterator<TIntSet> iter2 = forwardArcs.iterator();
        for (int ii = forwardArcs.size(); ii-- > 0;) {
            iter2.advance();
            int idx1 = iter2.key();
            TIntIterator iter3 = iter2.value().iterator();
            while (iter3.hasNext()) {
                int idx2 = iter3.next();
                PairInt p = new PairInt(idx1, idx2);
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
                if (unitFlow == 0) {
                    // idle, cp > -qEps
                    //cp = cost - pdX + pdY, so incr pdY
                    float count = 1.f;
                    while (cp <= -qEps) {
                        pRight[idx2] += (count * delta);
                        cp = calcNetCost(p);
                        log.fine("raise price of right " + idx2 +
                            " to cp=" + cp);
                    }
                } else if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated, cp <= qEps
                    //cp = cost - pdX + pdY, so incr pdX
                    while (cp > qEps) {
                        float count = 1.f;
                        pLeft[idx1] += (count * delta);
                        cp = calcNetCost(p);
                        log.fine("lower price of left " + idx1 +
                            " to cp=" + cp);
                    }
                }
            }
        }
        
        TIntIterator iter3 = sinkForwardArcs.iterator();
        // calc for sink dummy arcs
        while (iter3.hasNext()) {
            int yIdx = iter3.next();
            float unitFlow = rightToSinkF.get(yIdx);
            float cp = calcSinkNetCost(yIdx);
        
            if (unitFlow == 0) {
                // idle, cp > -qEps
                //cp = cost - pdX + pdY so inr pdY
                // where pdX is pRight and pdY is pRight[sink]
                float count = 1.f;
                while (cp <= -qEps) {
                    pRight[sinkNode] += (count * delta);
                    cp = calcSinkNetCost(yIdx);
                    log.fine("raise price of right for sink arc" 
                        + sinkNode + " to cp=" + cp);
                }
            } else if (Math.abs(unitFlow - 1) < 0.01f) {
                // saturated
                //cp = cost - pdX + pdY so incr pdX
                float count = 1.f;
                while (cp > qEps) {
                    pRight[yIdx] += (count * delta);
                    cp = calcSinkNetCost(yIdx);
                    log.fine("lower price of right for sink arc" 
                        + yIdx + " to cp=" + cp);
                }
            }
        }
    }
    
    /**
     * Every arc of NG, idle or saturated, is eps-proper.
     * all bipartite arcs are asserted, but not sink and source
     * @param eps
     * @return 
     */
    boolean integralBipartiteFlowIsEpsProper(float eps) {
        
        TIntObjectIterator<TIntSet> iter2 = forwardArcs.iterator();
        for (int ii = forwardArcs.size(); ii-- > 0;) {
            iter2.advance();
            int idx1 = iter2.key();
            TIntIterator iter3 = iter2.value().iterator();
            while (iter3.hasNext()) {
                int idx2 = iter3.next();
                                
                PairInt p = new PairInt(idx1, idx2);
                float unitFlow = f.get(p);
                int cp = (int)Math.ceil(calcNetCost(p));
                if (unitFlow == 0) {
                    // idle, cp > -epsT
                    if (cp <= -eps) {
                        return false;
                    }
                } else if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated
                    if (cp > eps) {
                        return false;
                    }
                }
            }
        }

        return true;
    }
    
    /**
     * assert pg 44 I3.
     * Every arc of NG, idle or saturated, is eps-proper.
     * all bipartite arcs and sink and source arcs should be eps-proper.
     * @param eps
     * @return 
     */
    boolean integralFlowIsEpsProper(float eps) {
        
        TIntObjectIterator<TIntSet> iter2 = forwardArcs.iterator();
        for (int ii = forwardArcs.size(); ii-- > 0;) {
            iter2.advance();
            int idx1 = iter2.key();
            TIntIterator iter3 = iter2.value().iterator();
            while (iter3.hasNext()) {
                int idx2 = iter3.next();
          
                PairInt p = new PairInt(idx1, idx2);
                float unitFlow = f.get(p);
                int cp = (int)Math.ceil(calcNetCost(p));
                if (unitFlow == 0) {
                    // idle, cp > -epsT
                    if (cp <= -eps) {
                        log.severe(String.format(
                            "NOT EPS-PROPER f=%.2f  cp=%.3f  eps=%.3f idx=%d to %d",
                            unitFlow, cp, eps, idx1, idx2));
                        return false;
                    }
                } else if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated
                    if (cp > eps) {
                        log.severe(String.format(
                            "NOT EPS-PROPER f=%.2f  cp=%.3f  eps=%.3f idx=%s td %d",
                            unitFlow, cp, eps, idx1, idx2));
                        return false;
                    }
                }
            }
        }

        TIntIterator iter = sourceForwardArcs.iterator();
        while (iter.hasNext()) {
            int idx1 = iter.next();
            float unitFlow = sourceToLeftF.get(idx1);
            int cp = (int)Math.ceil(calcSourceNetCost(idx1));
            if (unitFlow == 0) {
                // idle, cp > -epsT
                if (cp <= -eps) {
                    log.severe(String.format(
                        "NOT EPS-PROPER f=%.2f  cp=%.3f  eps=%.3f source to %d",
                        unitFlow, cp, eps, idx1));
                }
            } else if (Math.abs(unitFlow - 1) < 0.01f) {
                // saturated
                if (cp > eps) {
                    log.severe(String.format(
                        "NOT EPS-PROPER f=%.2f  cp=%.3f  eps=%.3f source to %d",
                        unitFlow, cp, eps, idx1));
                    return false;
                }
            }
        }
        
        TIntIterator iter1 = sinkForwardArcs.iterator();
        while (iter1.hasNext()) {
            int idx1 = iter1.next();
            float unitFlow = rightToSinkF.get(idx1);
            int cp = (int)Math.ceil(calcSinkNetCost(idx1));
            if (unitFlow == 0) {
                // idle, cp > -epsT
                if (cp <= -eps) {
                    log.severe(String.format(
                        "NOT EPS-PROPER f=%.2f  cp=%.3f  eps=%.3f %d to sink",
                        unitFlow, cp, eps, idx1));
                    return false;
                }
            } else if (Math.abs(unitFlow - 1) < 0.01f) {
                // saturated
                if (cp > eps) {
                    log.severe(String.format(
                        "NOT EPS-PROPER f=%.2f  cp=%.3f  eps=%.3f %d to sink",
                        unitFlow, cp, eps, idx1));
                    return false;
                }
            }
        }
        
        return true;
    }
    
    /**
     * Every arc of the network flow, idle or saturated, is proper.
     * @return 
     */
    public boolean integralFlowIsProper() {
        
        TIntObjectIterator<TIntSet> iter = forwardArcs.iterator();
        for (int ii = forwardArcs.size(); ii-- > 0;) {
            iter.advance();
            int idx1 = iter.key();
            TIntIterator iter2 = iter.value().iterator();
            while (iter2.hasNext()) {
                int idx2 = iter2.next();
                PairInt p = new PairInt(idx1, idx2);
                float unitFlow = f.get(p);
                int cp = (int)Math.ceil(calcNetCost(p));
                if (unitFlow == 0) {
                    // idle, cp >= 0
                    if (cp < -0.0f) {
                        return false;
                    }
                } else if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated
                    if (cp > 0.0f) {
                        return false;
                    }
                }
            }
        }

        return true;
    }
    
    /**
     * assert pg 44, I4.
     * Every saturated bipartite arc is eps-snug.
     * @param eps
     * @return 
     */
    boolean assertSaturatedBipartiteIsEpsSnug(float eps) {
        
        TIntObjectIterator<TIntSet> iter = forwardArcs.iterator();
        for (int ii = forwardArcs.size(); ii-- > 0;) {
            iter.advance();
            int idx1 = iter.key();
            
            if (idx1 == sourceNode) {
                continue;
            }
            
            TIntIterator iter3 = iter.value().iterator();
            while (iter3.hasNext()) {
                int idx2 = iter3.next();            
                if (idx2 == sinkNode) {
                    continue;
                }
                PairInt p = new PairInt(idx1, idx2);
                float unitFlow = f.get(p);
                int cp = (int)Math.ceil(calcNetCost(p));
                if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated
                    // snug is -eps < cp <= eps
                    if ((cp <= -eps) || (cp > eps)) {
                        log.severe(String.format(
                            "NOT EPS-SNUG f=%.2f  cp=%.3f  eps=%.3f idx=%d to %d",
                            unitFlow, cp, eps, idx1, idx2));
                        return false;
                    }
                }
            }
        }
        
        return true;
    }
    
    public void getSurplusLeftIndexes(TIntSet surplus) {

        for (int i = 0; i < nLeft; ++i) {
                        
            float flowInto = 0;
            float flowOutOf = 0;
            // flow into left from source
            if (sourceForwardArcs.contains(i)) {
                flowInto += sourceToLeftF.get(i);
            }
            TIntSet set = forwardArcs.get(i);
            if (set != null) {
                TIntIterator iter2 = set.iterator();
                while (iter2.hasNext()) {
                    int idx2 = iter2.next();
                    PairInt p = new PairInt(i, idx2);
                    flowOutOf += f.get(p);
                }
            }
            if (flowInto > flowOutOf) {
                surplus.add(i);
            }
        }
    }  
    
    public void getDeficitRightIndexes(TIntSet deficit) {

        // key = right, values = left
        TIntObjectHashMap<TIntSet> revMap = 
            createReverseMapOfForwardArcs();
        
        for (int i = 0; i < nRight; ++i) {
                        
            float flowInto = 0;
            float flowOutOf = 0;
            // flow out of right to sink
            if (sinkForwardArcs.contains(i)) {
                flowOutOf += rightToSinkF.get(i);
            }
            TIntSet set = revMap.get(i);
            if (set != null) {
                TIntIterator iter2 = set.iterator();
                while (iter2.hasNext()) {
                    int left = iter2.next();                
                    PairInt p = new PairInt(left, i);
                    flowInto += f.get(p);
                }
            }
            
            if (flowInto < flowOutOf) {
                deficit.add(i);
            }
        }
    }
    
    private void getSurplusRightIndexes(TIntSet surplus) {

        for (int i = 0; i < nRight; ++i) {
                        
            float flowInto = 0;
            float flowOutOf = 0;
            // flow into right from left node
            
            TIntObjectIterator<TIntSet> iter2 = forwardArcs.iterator();
            for (int ii = forwardArcs.size(); ii-- > 0;) {
                iter2.advance();
                int idx1 = iter2.key();
                TIntSet set = iter2.value();
                if (set.contains(i)) {
                    PairInt p = new PairInt(idx1, i);
                    flowInto += f.get(p);
                }
            }
            
            // flow from right to sink            
            if (sinkForwardArcs.contains(i)) {
                flowOutOf += rightToSinkF.get(i);
            }
            
            if (flowInto > flowOutOf) {
                surplus.add(i);
            }
        }
    }  
    
    /*
     the flow network uses forward arcs with "ceiling quantization".
     - the ceiling quantization is the choice of
     "eps-properness" for idle and saturated arcs whose
     net cost are equal to eps.
     those have been incorporated into integralFlowIsEpsProper
     - an idle arc v->w is eps-tight when -eps < cp(v,w) <= 0, 
     while a saturated arc v->w is eps-tight when 
     0 < cp(v,w) <= eps
     NOTE that eps-tight arcs are just barely eps-proper.
    */
    
    public TIntObjectHashMap<TIntSet> getForwardArcs() {
        return forwardArcs;
    }

    public TObjectIntHashMap<PairInt> getCosts() {
        return c;
    }

    public TObjectFloatHashMap<PairInt> getFlow() {
        return f;
    }

    public void setLeftPrice(int idx, float value) {
        pLeft[idx] = value;
    }

    public void setRightPrice(int idx, float value) {
        log.fine("add to right " + idx + " : " + pRight[idx] 
            + "  + " + value + " = " + (pRight[idx] + value));
        pRight[idx] = value;
    }
    
    public void addToLeftPrice(int idx, float value) {
        log.fine("add to left " + idx + " : " + pLeft[idx] 
            + " + " + value + " = " + (pLeft[idx] + value));
        pLeft[idx] += value;
    }

    public void addToSourcePrice(float value) {
        log.fine("add to source : " + pLeft[sourceNode] 
            + " + " + value + " = " + (pLeft[sourceNode] + value));
        pLeft[sourceNode] += value;
    }

    public void addToSinkPrice(float value) {
        log.fine("add to sink : " + pLeft[sinkNode] 
            + " + " + value + " = " + (pLeft[sinkNode] + value));
        pLeft[sinkNode] += value;
    }

    public void addToRightPrice(int idx, float value) {
        pRight[idx] += value;
    }
    
    public float getLeftPrice(int idx) {
        return pLeft[idx];
    }

    public float getRightPrice(int idx) {
        return pRight[idx];
    }

    /**
     * @return the number of left vertices
     */
    public int getNLeft() {
        return nLeft;
    }

    /**
     * @return the number of right vertices
     */
    public int getNRight() {
        return nRight;
    }
    
    public int getSinkNode() {
        return sinkNode;
    }
    
    public int getSourceNode() {
        return sourceNode;
    }

    public TIntSet getSourceForwardArcs() {
        return sourceForwardArcs;
    }

    public TIntSet getSinkForwardArcs() {
        return sinkForwardArcs;
    }

    public void augmentSourceToLeftFlowAndArc(int idx) {
        
        if (sourceToLeftC.containsKey(idx)) {
            float flow = sourceToLeftF.get(idx);
            if (Math.abs(flow) < 0.01) {
                // idle, so change to "saturated"
                sourceToLeftF.put(idx, 1);
                sourceForwardArcs.add(idx);
            } else if (Math.abs(flow - 1.) < 0.01) {
                // saturated, so change to "idle"
                sourceToLeftF.put(idx, 0);
                //TODO: consider not removing this
                //sourceForwardArcs.remove(index);
            } else {
                throw new IllegalStateException(
                "Error in algorithm.  not expecting"
                + " a fractional flow");
            }
        }
    }
    
    public void augmentRightToSinkFlowAndArc(int idx) {
            
        if (rightToSinkF.containsKey(idx)) {
            float flow = rightToSinkF.get(idx);
            if (Math.abs(flow) < 0.01) {
                // idle, so change to "saturated"
                rightToSinkF.put(idx, 1);
                sinkForwardArcs.add(idx);
            } else if (Math.abs(flow - 1.) < 0.01) {
                // saturated, so change to "idle"
                rightToSinkF.put(idx, 0);
                // TODO: consider not removing this:
                //sinkForwardArcs.remove(idx);
            } else {
                throw new IllegalStateException(
                    "Error in algorithm.  not expecting"
                    + " a fractional flow");
            }
        }
    }
    
    public void augmentFlowAndArc(int idx1, int idx2) {
        
        assert(forwardArcs.containsKey(idx1)
            && forwardArcs.get(idx1)
            .contains(idx2));
        
        PairInt p = new PairInt(idx1, idx2);
        
        float flow = f.get(p);
        if (Math.abs(flow - 1.) < 0.01) {
            // saturated, so change to "idle"
            f.put(p, 0);
        } else if (Math.abs(flow) < 0.01) {
            // idle, so change to "saturated"
            f.put(p, 1);
            
            // should the source and sink
            // arcs be reversed too?
            
        } else {
            throw new IllegalStateException(
                "Error in algorithm.  not expecting"
                + " a fractional flow");
        }
    }

    private TIntObjectHashMap<TIntSet> createReverseMapOfForwardArcs() {
        
        //key=right   values=left
        TIntObjectHashMap<TIntSet> revMap 
              = new TIntObjectHashMap<TIntSet>();

        TIntObjectIterator<TIntSet> iter2 = forwardArcs.iterator();
        for (int ii = forwardArcs.size(); ii-- > 0;) {
            iter2.advance();
            int leftIdx = iter2.key();
            TIntSet rightIdxs = iter2.value();
            TIntIterator iter3 = rightIdxs.iterator();
            while (iter3.hasNext()) {
                int rightIdx = iter3.next();
                TIntSet leftIdxs;
                if (!revMap.containsKey(rightIdx)) {
                    leftIdxs = new TIntHashSet();
                    revMap.put(rightIdx, leftIdxs);
                } else {
                    leftIdxs = revMap.get(rightIdx);
                } 
                leftIdxs.add(leftIdx);
            }
        }
        
        return revMap;
    }
        
    public int calculateNumberOfSaturatedArcs() {
    
        int n = 0;
        
        TIntObjectIterator<TIntSet> iter = forwardArcs.iterator();
        for (int ii = forwardArcs.size(); ii-- > 0;) {
            iter.advance();
            int idx1 = iter.key();
            TIntSet rightIdxs = iter.value();
            TIntIterator iter2 = rightIdxs.iterator();
            while (iter2.hasNext()) {
                int idx2 = iter2.next();
                PairInt p = new PairInt(idx1, idx2);
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
                if (Math.abs(unitFlow - 1) < 0.01f) {
                    n++;
                }
            }
        }
        
        return n;
    }

    public void printSaturatedLinks() {

        TIntObjectIterator<TIntSet> iter = forwardArcs.iterator();
        for (int ii = forwardArcs.size(); ii-- > 0;) {
            iter.advance();
            int idx1 = iter.key();
            TIntSet rightIdxs = iter.value();
            TIntIterator iter2 = rightIdxs.iterator();
            while (iter2.hasNext()) {
                int idx2 = iter2.next();
                PairInt p = new PairInt(idx1, idx2);
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
                if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated are matched arcs
                    log.info("saturated bipartite arc from " +
                        idx1 + " to " + idx2);
                }
            }
        }
        
    }
    
    public void printNetCosts() {

        StringBuilder sb = new StringBuilder();
        
        TIntObjectIterator<TIntSet> iter = forwardArcs.iterator();
        for (int ii = forwardArcs.size(); ii-- > 0;) {
            iter.advance();
            int idx1 = iter.key();
            TIntSet rightIdxs = iter.value();
            TIntIterator iter2 = rightIdxs.iterator();
            while (iter2.hasNext()) {
                int idx2 = iter2.next();
                PairInt p = new PairInt(idx1, idx2);
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);

                sb.append(String.format("%d to %d cp=%.2f f=%.2f\n",
                    idx1, idx2, cp, unitFlow));
            }
        }
        log.info(sb.toString());

        sb = new StringBuilder();

        TIntIterator iter2 = sourceForwardArcs.iterator();
        while (iter2.hasNext()) {
            int idx1 = iter2.next();

            float unitFlow = sourceToLeftF.get(idx1);
            float cp = calcSourceNetCost(idx1);
            sb.append(String.format("source to %d cp=%.2f f=%.2f\n",
                idx1, cp, unitFlow));
        }
        log.info(sb.toString());

        sb = new StringBuilder();
        iter2 = sinkForwardArcs.iterator();
        while (iter2.hasNext()) {
            int idx1 = iter2.next();
            float unitFlow = rightToSinkF.get(idx1);
            float cp = calcSinkNetCost(idx1);
            sb.append(String.format("%d to sink cp=%.2f f=%.2f\n",
                idx1, cp, unitFlow));
        }
        log.info(sb.toString());
    }

    public TIntIntMap extractMatches() {

        TIntIntMap m = new TIntIntHashMap();

        TIntObjectIterator<TIntSet> iter = forwardArcs.iterator();
        for (int ii = forwardArcs.size(); ii-- > 0;) {
            iter.advance();
            int idx1 = iter.key();
            TIntSet rightIdxs = iter.value();
            TIntIterator iter2 = rightIdxs.iterator();
            while (iter2.hasNext()) {
                int idx2 = iter2.next();
                PairInt p = new PairInt(idx1, idx2);
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
                if (Math.abs(unitFlow - 1) < 0.01f) {
                    m.put(idx1, idx2);
                }
            }
        }
        
        return m;
    }

    /**
     * given value, modify prices such that
     * //pd^~(v) = math.floor(pd(v) + value)
     * @param value
     */
    public void addToAllPrices(float value) {

        for (int i = 0; i < nLeft; ++i) {
            pLeft[i] = (float)Math.floor(pLeft[i] + value);
        }
        
        for (int i = 0; i < nRight; ++i) {
            pRight[i] = (float)Math.floor(pRight[i] + value);
        }
    }

    public void zeroTheMatchedBipartiteFlow(TIntSet outPutSurplus, 
        TIntSet outPutDeficit) {
        
        outPutSurplus.clear();
        outPutDeficit.clear();

        TIntObjectIterator<TIntSet> iter = forwardArcs.iterator();
        for (int ii = forwardArcs.size(); ii-- > 0;) {
            iter.advance();
            int idx1 = iter.key();
            TIntSet rightIdxs = iter.value();
            TIntIterator iter2 = rightIdxs.iterator();
            while (iter2.hasNext()) {
                int idx2 = iter2.next();
                PairInt p = new PairInt(idx1, idx2);
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
        
                if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated are matched arcs
                    outPutSurplus.add(idx1);
                    outPutDeficit.add(idx2);
                    f.put(p, 0);
                }
            }
        }
    }

    public void createPathNodes() {
        if (pathNodes == null) {
            pathNodes = new PathNodes(nLeft, nRight);
        }
    }
    
    public void resetPathNodes(long key) {
        if (pathNodes == null) {
            pathNodes = new PathNodes(nLeft, nRight);
        } else {
            pathNodes.resetNodeExceptData(key);
        }
    }
    
    public PathNodes getPathNodes() {
        return pathNodes;
    }
    
}
