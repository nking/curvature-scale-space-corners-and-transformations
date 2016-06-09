package algorithms.bipartite;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
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
    private final Map<Integer, Set<Integer>> forwardArcs
        = new HashMap<Integer, Set<Integer>>();

    /** forward arcs from source to the left (a.k.a. X nodes).
     * If an arc is saturated, it will have flux=1 along this
     * arc.
    */
    private final Set<Integer> sourceForwardArcs = 
        new HashSet<Integer>();

    /** forward arcs from the right (a.k.a. Y nodes) to sink node
    If an arc is saturated, it will have flux=1 along this arc.
    */
    private final Set<Integer> sinkForwardArcs =
        new HashSet<Integer>();

    /**
     * per unit flow costs are the same as the edges given in G. 
     * key = index in left (==X==v) and index in right (==Y==w), 
     * value = cost of the arc. might
     * change to use sparse matrix format
     */
    private final Map<PairInt, Integer> c
        = new HashMap<PairInt, Integer>();

    /**
     * the flow on the arc. key = index in left (==X==v) and index in right
     * (==Y==w), value = number from 0 to 1 inclusive if it's a "pseudoflow":
     * the value 0 is "idle" and corresponds to an unmatched link in the
     * residual graph. the value 1 is "saturated" and corresponds to a matched
     * link in the residual graph (but in the residual graph, it would be
     * specified in format right to left). might change to use sparse matrix
     * format
     */
    private final Map<PairInt, Float> f
        = new HashMap<PairInt, Float>();
    
    /**
     * per unit dummy flow costs from source to X nodes
     * (initialized as 0).
     * key = index in left (==X==v) 
     * value = cost of the arc. might
     * change to use sparse matrix format
     */
    private final Map<Integer, Integer> sourceToLeftC =
        new HashMap<Integer, Integer>();
   
    /**
     * per unit dummy flow costs from sink to X nodes
     * (initialized as 0).
     * key = index of right node (== Y)
     * value = cost of the arc
     */
    private final Map<Integer, Integer> rightToSinkC =
        new HashMap<Integer, Integer>();

    /**
     * flux from source to x node
     */
    private final Map<Integer, Float> sourceToLeftF =
        new HashMap<Integer, Float>();
    
    /**
     * flux from y node to sink
     */
    private final Map<Integer, Float> rightToSinkF =
        new HashMap<Integer, Float>();
    
    public FlowNetwork(Graph g, Map<Integer, Integer> m) {

        this.nLeft = g.getNLeft();
        this.nRight = g.getNRight();
        this.sourceNode = nLeft;
        this.sinkNode = nRight;
        pLeft = new float[nLeft + 1];
        pRight = new float[nRight + 1];
                
        for (Map.Entry<PairInt, Integer> entry : 
            g.getEdgeWeights().entrySet()) {

            PairInt p = entry.getKey();

            Integer index1 = Integer.valueOf(p.getX());
            Integer index2 = Integer.valueOf(p.getY());

            Set<Integer> indexes2 = forwardArcs.get(index1);
            if (indexes2 == null) {
                indexes2 = new HashSet<Integer>();
                forwardArcs.put(index1, indexes2);
            }
            indexes2.add(index2);

            Integer cost = g.getEdgeWeights().get(p);
            c.put(p, cost);
            f.put(p, Float.valueOf(0));

            if (Math.abs(cost.intValue()) > maxC) {
                maxC = Math.abs(cost.intValue());
            }
            if (Math.abs(cost.intValue()) < minC) {
                minC = Math.abs(cost.intValue());
            }

            sourceForwardArcs.add(index1);
            sinkForwardArcs.add(index2);
            
            // set source and sink arcs to 0 and overwrite
            // those for matched vertexes in next block
            sourceToLeftF.put(index1, Float.valueOf(0));
            rightToSinkF.put(index2, Float.valueOf(0));
            
            // cost of dummy edge is 0
            sourceToLeftC.put(index1, Integer.valueOf(0));
            rightToSinkC.put(index2, Integer.valueOf(0));
        }
        
        for (Entry<Integer, Integer> entry : m.entrySet()) {
            Integer leftIndex = entry.getKey();
            sourceToLeftF.put(leftIndex, Float.valueOf(1));
            PairInt p = new PairInt(leftIndex.intValue(),
                entry.getValue().intValue());
            f.put(p, Float.valueOf(1));
        }
        for (Integer rightIndex : m.values()) {
            rightToSinkF.put(rightIndex, Float.valueOf(1));
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
        
        for (Map.Entry<Integer, Set<Integer>> entry : forwardArcs.entrySet()) {
            Integer index1 = entry.getKey();
            for (Integer index2 : entry.getValue()) {
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
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
        
        for (Map.Entry<Integer, Set<Integer>> entry : forwardArcs.entrySet()) {
            Integer index1 = entry.getKey();
            for (Integer index2 : entry.getValue()) {
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
                float unitFlow = f.get(p);
                float cost = c.get(p);
                sum += (unitFlow * cost);
            }
        }
        
        for (Integer index1 : sourceForwardArcs) {
            float unitFlow = sourceToLeftF.get(index1);
            float cost = sourceToLeftC.get(index1);
            sum += (unitFlow * cost);
        }
        
        for (Integer index1 : sinkForwardArcs) {
            float unitFlow = rightToSinkF.get(index1);
            float cost = rightToSinkC.get(index1);
            sum += (unitFlow * cost);
        }
        
        return sum;
    }
    
    public float getSourceToLeftFlow(int idx) {
        return sourceToLeftF.get(Integer.valueOf(idx));
    }
    public float getRightToSinkFlow(int idx) {
        return rightToSinkF.get(Integer.valueOf(idx));
    }
    public int getSourceToLeftCost(int idx) {
        return sourceToLeftC.get(Integer.valueOf(idx));
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
            return sourceToLeftF.get(Integer.valueOf(v));
        } else if (v == sinkNode) {
            return rightToSinkF.get(Integer.valueOf(u));
        } 
        
        return f.get(new PairInt(u, v));
    }
    
    protected float calcSourceNetCost(int v) {

        int cost = sourceToLeftC.get(Integer.valueOf(v));
        float pdX = pLeft[sourceNode];
        float pdY = pLeft[v];
       
        float cp = cost - pdX + pdY;

        return cp;
    }

    protected float calcSinkNetCost(int u) {

        int cost = rightToSinkC.get(Integer.valueOf(u));
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
        
        List<Integer> surplus = new ArrayList<Integer>();
        getSurplusRightIndexes(surplus);
        assert(surplus.isEmpty());
        
        getSurplusLeftIndexes(surplus);
        
        int nSurplus = surplus.size();
        
        log.info("s=" + nMatchingsHK + " h=" + nSurplus);
                
        return (Math.abs(flow - (nMatchingsHK - nSurplus)) < 1);
    }
    
    public boolean printSurplusAndDeficit() {
    
        List<Integer> surplus = new ArrayList<Integer>();
        getSurplusLeftIndexes(surplus);
        
        List<Integer> deficit = new ArrayList<Integer>();
        getDeficitRightIndexes(deficit);
        
        StringBuilder sb = new StringBuilder("surplus=(");
        for (Integer s : surplus) {
            sb.append(s.toString()).append(", ");
        }
        sb.append(")");
        log.info(sb.toString());
        
        sb = new StringBuilder("deficit=(");
        for (Integer d : deficit) {
            sb.append(d.toString()).append(", ");
        }
        sb.append(")");
        log.info(sb.toString());
        
        return true;
    }
  
    double calcTotalSourceFlow() {
                
        double flow = 0;
        
        for (Integer index : sourceForwardArcs) {
            float unitFlow = sourceToLeftF.get(index);
            flow += unitFlow;
        }
                
        return flow;
    }
    
    double calcTotalSinkFlow() {
                
        double flow = 0;
        
        for (Integer index : sinkForwardArcs) {
            float unitFlow = rightToSinkF.get(index);
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

        // calc for source dummy arcs
        for (Integer xIndex : sourceForwardArcs) {
            float unitFlow = sourceToLeftF.get(xIndex);
            float cp = calcSourceNetCost(xIndex.intValue());
        
            if (unitFlow == 0) {
                // idle, cp > -qEps
                //cp = cost - pdX + pdY so inr pdY
                float count = 1.f;
                while (cp <= -qEps) {
                    pLeft[xIndex.intValue()] += (count * delta);
                    cp = calcSourceNetCost(xIndex.intValue());
                    log.fine("raise price of left for source arc" 
                        + xIndex + " to cp=" + cp);
                }
            } else if (Math.abs(unitFlow - 1) < 0.01f) {
                // saturated
                //cp = cost - pdX + pdY so incr pdX
                float count = 1.f;
                while (cp > qEps) {
                    pLeft[sourceNode] += (count * delta);
                    cp = calcSourceNetCost(xIndex.intValue());
                    log.fine("lower price of left for source arc" 
                        + xIndex + " to cp=" + cp);
                }
            }
        }
        
        for (Map.Entry<Integer, Set<Integer>> entry : 
            forwardArcs.entrySet()) {
            
            Integer index1 = entry.getKey();
            
            for (Integer index2 : entry.getValue()) {
                
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
                if (unitFlow == 0) {
                    // idle, cp > -qEps
                    //cp = cost - pdX + pdY, so incr pdY
                    float count = 1.f;
                    while (cp <= -qEps) {
                        pRight[index2.intValue()] += (count * delta);
                        cp = calcNetCost(p);
                        log.fine("raise price of right " + index2 +
                            " to cp=" + cp);
                    }
                } else if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated, cp <= qEps
                    //cp = cost - pdX + pdY, so incr pdX
                    while (cp > qEps) {
                        float count = 1.f;
                        pLeft[index1.intValue()] += (count * delta);
                        cp = calcNetCost(p);
                        log.fine("lower price of left " + index1 +
                            " to cp=" + cp);
                    }
                }
            }
        }
        
        // calc for sink dummy arcs
        for (Integer yIndex : sinkForwardArcs) {
            float unitFlow = rightToSinkF.get(yIndex);
            float cp = calcSinkNetCost(yIndex.intValue());
        
            if (unitFlow == 0) {
                // idle, cp > -qEps
                //cp = cost - pdX + pdY so inr pdY
                // where pdX is pRight and pdY is pRight[sink]
                float count = 1.f;
                while (cp <= -qEps) {
                    pRight[sinkNode] += (count * delta);
                    cp = calcSinkNetCost(yIndex.intValue());
                    log.fine("raise price of right for sink arc" 
                        + sinkNode + " to cp=" + cp);
                }
            } else if (Math.abs(unitFlow - 1) < 0.01f) {
                // saturated
                //cp = cost - pdX + pdY so incr pdX
                float count = 1.f;
                while (cp > qEps) {
                    pRight[yIndex.intValue()] += (count * delta);
                    cp = calcSinkNetCost(yIndex.intValue());
                    log.fine("lower price of right for sink arc" 
                        + yIndex + " to cp=" + cp);
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
        
        for (Map.Entry<Integer, Set<Integer>> entry : 
            forwardArcs.entrySet()) {
            
            Integer index1 = entry.getKey();
            
            for (Integer index2 : entry.getValue()) {
                
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
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
        
        for (Map.Entry<Integer, Set<Integer>> entry : 
            forwardArcs.entrySet()) {
            
            Integer index1 = entry.getKey();
            
            for (Integer index2 : entry.getValue()) {
                
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
                if (unitFlow == 0) {
                    // idle, cp > -epsT
                    if (cp <= -eps) {
    log.info(String.format
    ("NOT EPS-PROPER f=%.2f  cp=%.3f  eps=%.3f idx=%s to %s", 
    unitFlow, cp, eps, index1.toString(), index2.toString()));
                        return false;
                    }
                } else if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated
                    if (cp > eps) {
     log.info(String.format
    ("NOT EPS-PROPER f=%.2f  cp=%.3f  eps=%.3f idx=%s to %s", 
    unitFlow, cp, eps, index1.toString(), index2.toString()));
                        return false;
                    }
                }
            }
        }

        for (Integer index1 : sourceForwardArcs) {
            float unitFlow = sourceToLeftF.get(index1);
            float cp = calcSourceNetCost(index1.intValue());
            if (unitFlow == 0) {
                // idle, cp > -epsT
                if (cp <= -eps) {
    log.info(String.format
    ("NOT EPS-PROPER f=%.2f  cp=%.3f  eps=%.3f source to %s", 
    unitFlow, cp, eps, index1.toString()));
                }
            } else if (Math.abs(unitFlow - 1) < 0.01f) {
                // saturated
                if (cp > eps) {
    log.info(String.format
    ("NOT EPS-PROPER f=%.2f  cp=%.3f  eps=%.3f source to %s", 
    unitFlow, cp, eps, index1.toString()));
                    return false;
                }
            }
        }
        
        for (Integer index1 : sinkForwardArcs) {
            float unitFlow = rightToSinkF.get(index1);
            float cp = calcSinkNetCost(index1.intValue());
            if (unitFlow == 0) {
                // idle, cp > -epsT
                if (cp <= -eps) {
                    log.info(String.format
    ("NOT EPS-PROPER f=%.2f  cp=%.3f  eps=%.3f %s to sink", 
    unitFlow, cp, eps, index1.toString()));
                    return false;
                }
            } else if (Math.abs(unitFlow - 1) < 0.01f) {
                // saturated
                if (cp > eps) {
                    log.info(String.format
    ("NOT EPS-PROPER f=%.2f  cp=%.3f  eps=%.3f %s to sink", 
    unitFlow, cp, eps, index1.toString()));
                    return false;
                }
            }
        }
        return true;
    }
    /**
     * Every arc of the network flow, idle or saturated, is proper.
     * @param eps
     * @return 
     */
    public boolean integralFlowIsProper() {
        
        for (Map.Entry<Integer, Set<Integer>> entry : 
            forwardArcs.entrySet()) {
            
            Integer index1 = entry.getKey();
            
            for (Integer index2 : entry.getValue()) {
                
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
                if (unitFlow == 0) {
                    // idle, cp >= 0
                    if (cp < 0) {
                        return false;
                    }
                } else if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated
                    if (cp > 0) {
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
        
        for (Map.Entry<Integer, Set<Integer>> entry : 
            forwardArcs.entrySet()) {
            
            Integer index1 = entry.getKey();
            if (index1.intValue() == sourceNode) {
                continue;
            }
            for (Integer index2 : entry.getValue()) {
                if (index2.intValue() == sinkNode) {
                    continue;
                }
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
                if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated
                    // snug is -eps < cp <= eps
                    if ((cp <= -eps) || (cp > eps)) {
    log.info(String.format
    ("NOT EPS-SNUG f=%.2f  cp=%.3f  eps=%.3f idx=%s to %s", 
    unitFlow, cp, eps, index1.toString(), index2.toString()));
                        return false;
                    }
                }
            }
        }
        
        return true;
    }
    
    /**
     * populate given lists with the bipartite nodes of 
     * saturated arcs.
     * @param surplus
     * @param deficit 
     */
    public void getMatchedLeftRight(List<Integer> surplus, List<Integer> deficit) {

        for (Map.Entry<Integer, Set<Integer>> entry : 
            forwardArcs.entrySet()) {
            
            Integer index1 = entry.getKey();
            for (Integer index2 : entry.getValue()) {
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
                if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated are matched arcs
                    surplus.add(index1);
                    deficit.add(index2);
                }
            }
        }
    }

    public void getSurplusLeftIndexes(List<Integer> surplus) {

        for (int i = 0; i < nLeft; ++i) {
            
            Integer index = Integer.valueOf(i);        
            
            float flowInto = 0;
            float flowOutOf = 0;
            // flow into left from source
            if (sourceForwardArcs.contains(index)) {
                flowInto += sourceToLeftF.get(index);
            }
            Set<Integer> set = forwardArcs.get(index);
            if (set != null) {
                for (Integer index2 : set) {
                    PairInt p = new PairInt(index.intValue(), index2.intValue());
                    flowOutOf += f.get(p);
                }
            }
            if (flowInto > flowOutOf) {
                surplus.add(index);
            }
        }
    }  
    
    public void getDeficitRightIndexes(List<Integer> deficit) {

        // key = right, values = left
        Map<Integer, Set<Integer>> revMap =
            createReverseMapOfForwardArcs();
        
        for (int i = 0; i < nRight; ++i) {
            
            Integer index = Integer.valueOf(i);        
            
            float flowInto = 0;
            float flowOutOf = 0;
            // flow out of right to sink
            if (sinkForwardArcs.contains(index)) {
                flowOutOf += rightToSinkF.get(index);
            }
            Set<Integer> set = revMap.get(index);
            if (set != null) {
                for (Integer left : set) {
                    PairInt p = new PairInt(left.intValue(), index.intValue());
                    flowInto += f.get(p);
                }
            }
            
            if (flowInto < flowOutOf) {
                deficit.add(index);
            }
        }
    }  
    
    private void getSurplusRightIndexes(List<Integer> surplus) {

        for (int i = 0; i < nRight; ++i) {
            
            Integer index = Integer.valueOf(i);        
            
            float flowInto = 0;
            float flowOutOf = 0;
            // flow into right from left node
            for (Entry<Integer, Set<Integer>> entry :
                forwardArcs.entrySet()) {
                if (entry.getValue().contains(index)) {
                    PairInt p = new PairInt(entry.getKey().intValue(), 
                        index.intValue());
                    flowInto += f.get(p);
                }
            }
            
            // flow from right to sink            
            if (sinkForwardArcs.contains(index)) {
                flowOutOf += rightToSinkF.get(index);
            }
            
            if (flowInto > flowOutOf) {
                surplus.add(index);
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
    
    public Map<Integer, Set<Integer>> getForwardArcs() {
        return forwardArcs;
    }

    public Map<PairInt, Integer> getCosts() {
        return c;
    }

    public Map<PairInt, Float> getFlow() {
        return f;
    }

    public void setLeftPrice(int idx, float value) {
        pLeft[idx] = value;
    }

    public void setRightPrice(int idx, float value) {
        log.info("add to right " + idx + " : " + pRight[idx] + " + " + value + " = " + 
            (pRight[idx] + value));
        pRight[idx] = value;
    }
    
    public void addToLeftPrice(int idx, float value) {
        log.info("add to left " + idx + " : " + pLeft[idx] + " + " + value + " = " + 
            (pLeft[idx] + value));
        pLeft[idx] += value;
    }

    public void addToSourcePrice(float value) {
        log.info("add to source : " + pLeft[sourceNode] + " + " + value + " = " + 
            (pLeft[sourceNode] + value));
        pLeft[sourceNode] += value;
    }

    public void addToSinkPrice(float value) {
        log.info("add to sink : " + pLeft[sinkNode] + " + " + value + " = " + 
            (pLeft[sinkNode] + value));
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

    public Set<Integer> getSourceForwardArcs() {
        return sourceForwardArcs;
    }

    public Set<Integer> getSinkForwardArcs() {
        return sinkForwardArcs;
    }

    public void augmentSourceToLeftFlowAndArc(int idx) {
        
        Integer index = Integer.valueOf(idx);
        
        Float flow = sourceToLeftF.get(index);
        if (flow == null || (Math.abs(flow.floatValue()) < 0.01)) {
            // idle, so change to "saturated"
            sourceToLeftF.put(index, Float.valueOf(1));
            sourceForwardArcs.add(index);
        } else if (Math.abs(flow.floatValue() - 1.) < 0.01) {
            // saturated, so change to "idle"
            sourceToLeftF.put(index, Float.valueOf(0));
            //TODO: consider not removing this
            //sourceForwardArcs.remove(index);
        } else {
            throw new IllegalStateException(
                "Error in algorithm.  not expecting"
                + " a fractional flow");
        }
    }
    
    public void augmentRightToSinkFlowAndArc(int idx) {
        
        Integer index = Integer.valueOf(idx);
        
        Float flow = rightToSinkF.get(index);
        if (flow == null || (Math.abs(flow.floatValue()) < 0.01)) {
            // idle, so change to "saturated"
            rightToSinkF.put(index, Float.valueOf(1));
            sinkForwardArcs.add(index);
        } else if (Math.abs(flow.floatValue() - 1.) < 0.01) {
            // saturated, so change to "idle"
            rightToSinkF.put(index, Float.valueOf(0));
            // TODO: consider not removing this:
            //sinkForwardArcs.remove(index);
        } else {
            throw new IllegalStateException(
                "Error in algorithm.  not expecting"
                + " a fractional flow");
        }
    }
    
    public void augmentFlowAndArc(int idx1, int idx2) {
        
        assert(forwardArcs.containsKey(Integer.valueOf(idx1))
            && forwardArcs.get(Integer.valueOf(idx1))
            .contains(Integer.valueOf(idx2)));
        
        PairInt p = new PairInt(idx1, idx2);
        
        float flow = f.get(p);
        if (Math.abs(flow - 1.) < 0.01) {
            // saturated, so change to "idle"
            f.put(p, Float.valueOf(0));
        } else if (Math.abs(flow) < 0.01) {
            // idle, so change to "saturated"
            f.put(p, Float.valueOf(1));
            
            // should the source and sink
            // arcs be reversed too?
            
        } else {
            throw new IllegalStateException(
                "Error in algorithm.  not expecting"
                + " a fractional flow");
        }
    }

    private Map<Integer, Set<Integer>> 
        createReverseMapOfForwardArcs() {
        
        //key=right   values=left
        Map<Integer, Set<Integer>> revMap
            = new HashMap<Integer, Set<Integer>>();
    
        for (Entry<Integer, Set<Integer>> entry :
            forwardArcs.entrySet()) {
            
            Integer leftIndex = entry.getKey();
            Set<Integer> rightIndexes = entry.getValue();
            
            for (Integer rightIndex : rightIndexes) {
                Set<Integer> leftIndexes = revMap.get(rightIndex);
                if (leftIndexes == null) {
                    leftIndexes = new HashSet<Integer>();
                    revMap.put(rightIndex, leftIndexes);
                }
                leftIndexes.add(leftIndex);
            }
        }
        
        return revMap;
    }
        
    public int calculateNumberOfSaturatedArcs() {
    
        int n = 0;
        
        for (Map.Entry<Integer, Set<Integer>> entry : 
            forwardArcs.entrySet()) {
            
            Integer index1 = entry.getKey();
            for (Integer index2 : entry.getValue()) {
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
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

        for (Map.Entry<Integer, Set<Integer>> entry : 
            forwardArcs.entrySet()) {
            
            Integer index1 = entry.getKey();
            for (Integer index2 : entry.getValue()) {
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
                if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated are matched arcs
                    log.info("saturated bipartite arc from " +
                        index1 + " to " + index2);
                }
            }
        }
        
    }
    
    public void printNetCosts() {

        StringBuilder sb = new StringBuilder();
        
        for (Map.Entry<Integer, Set<Integer>> entry : 
            forwardArcs.entrySet()) {
            
            Integer index1 = entry.getKey();
            for (Integer index2 : entry.getValue()) {
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
                sb.append(String.format("%d to %d cp=%.2f f=%.2f\n",
                    index1.intValue(), index2.intValue(), cp, unitFlow));
            }
        }
        log.info(sb.toString());

        sb = new StringBuilder();
        for (Integer index1 : sourceForwardArcs) {
            float unitFlow = sourceToLeftF.get(index1);
            float cp = calcSourceNetCost(index1.intValue());
            sb.append(String.format("source to %d cp=%.2f f=%.2f\n",
                index1.intValue(), cp, unitFlow));
        }
        log.info(sb.toString());

        for (Integer index1 : sinkForwardArcs) {
            float unitFlow = rightToSinkF.get(index1);
            float cp = calcSinkNetCost(index1.intValue());
            sb.append(String.format("%d to sink cp=%.2f f=%.2f\n",
                index1.intValue(), cp, unitFlow));
        }
        log.info(sb.toString());
    }

    public Map<Integer, Integer> extractMatches() {
        
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();

        for (Map.Entry<Integer, Set<Integer>> entry : 
            forwardArcs.entrySet()) {
            
            Integer index1 = entry.getKey();
            for (Integer index2 : entry.getValue()) {
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
                if (Math.abs(unitFlow - 1) < 0.01f) {
                    m.put(index1, index2);
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

}
