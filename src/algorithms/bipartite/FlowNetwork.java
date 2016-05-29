package algorithms.bipartite;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * class representing a flow network for 
 * MinCostUnbalancedAssignment.java
   If the sink and source nodes and weights are set in the
   origin graph, sink and source logic is used in code.
 * 
 * @author nichole
 */
public class FlowNetwork {
    
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

            if (Math.abs(cost.intValue()) > maxC) {
                maxC = Math.abs(cost.intValue());
            }

            sourceForwardArcs.add(index1);
            sinkForwardArcs.add(index2);
            
            if (m.containsKey(index1) && m.get(index1).equals(index2)) {
                // create saturated arc, f(v,w)=1, net cost cp <= 0
                f.put(p, Float.valueOf(1));
                sourceToLeftF.put(index1, Float.valueOf(1));
                rightToSinkF.put(index2, Float.valueOf(1));
            } else {
                // create idle arc, f(v,w)=0, net cost cp >= 0
                f.put(p, Float.valueOf(0));                
                sourceToLeftF.put(index1, Float.valueOf(0));
                rightToSinkF.put(index2, Float.valueOf(0));
            }
            
            // cost of dummy edge is 0
            sourceToLeftC.put(index1, Integer.valueOf(0));
            rightToSinkC.put(index2, Integer.valueOf(0));
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
        
        return (Math.abs(flow - s) < 1);
    }
    
    /**
     * assert pg 52, I1'.
     * assert that flux f on NG is a flow of value |f|=s
     * @param s
     * @return 
     */
    boolean assertFlowValueIncludingSrcSnk(int s) {
                
        double flow = 0;
        
        for (Map.Entry<Integer, Set<Integer>> entry : forwardArcs.entrySet()) {
            Integer index1 = entry.getKey();
            for (Integer index2 : entry.getValue()) {
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
                float unitFlow = f.get(p);
                flow += unitFlow;
            }
        }
        
        for (Integer index1 : sourceForwardArcs) {
            float unitFlow = sourceToLeftF.get(index1);
            flow += unitFlow;
        }
        
        for (Integer index1 : sinkForwardArcs) {
            float unitFlow = rightToSinkF.get(index1);
            flow += unitFlow;
        }
        
        // calculate the number of source nodes == number of
        // deficit nodes
        List<Integer> surplus = new ArrayList<Integer>();        
        List<Integer> deficit = new ArrayList<Integer>();
        getMatchedLeftRight(surplus, deficit);
        assert(surplus.size() == deficit.size());
        int h = surplus.size();
                
        return (Math.abs(flow - (s - h)) < 1);
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
            if (r > 0.01) {
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
        
        //see Figure 7.4, pg 53
        // looks like the price raise is multiples of (q-1)*eps
        
        for (Map.Entry<Integer, Set<Integer>> entry : 
            forwardArcs.entrySet()) {
            
            Integer index1 = entry.getKey();
            
            for (Integer index2 : entry.getValue()) {
                
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
                if (unitFlow == 0) {
                    // idle, cp > -qEps
                    //cp = cost - pdX + pdY so inr pdY
                    float count = 1.f;
                    while (cp <= -qEps) {
                        pRight[index2.intValue()] += (count * delta);
                        cp = calcNetCost(p);
                        count += 1.f;
                    }
                } else if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated
                    //cp = cost - pdX + pdY so incr pdX
                    //TODO: a single price raise of 2*delta or 3*delta?
                    //      difficult to tell by pg 54 and 53
                    while (cp > qEps) {
                        float count = 1.f;
                        pLeft[index2.intValue()] += (count * delta);
                        cp = calcNetCost(p);
                        count += 1.f;
                    }
                }
            }
        }
        
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
                    count += 1.f;
                }
            } else if (Math.abs(unitFlow - 1) < 0.01f) {
                // saturated
                //cp = cost - pdX + pdY so incr pdX
                //TODO: a single price raise of 2*delta or 3*delta?
                //      difficult to tell by pg 54 and 53
                float count = 1.f;
                while (cp > qEps) {
                    pLeft[xIndex.intValue()] += (count * delta);
                    cp = calcSourceNetCost(xIndex.intValue());
                    count += 1.f;
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
                float count = 1.f;
                while (cp <= -qEps) {
                    pRight[yIndex.intValue()] += (count * delta);
                    cp = calcSinkNetCost(yIndex.intValue());
                    count += 1.f;
                }
            } else if (Math.abs(unitFlow - 1) < 0.01f) {
                // saturated
                //cp = cost - pdX + pdY so incr pdX
                //TODO: a single price raise of 2*delta or 3*delta?
                //      difficult to tell by pg 54 and 53
                float count = 1.f;
                while (cp > qEps) {
                    pRight[yIndex.intValue()] += (count * delta);
                    cp = calcSinkNetCost(yIndex.intValue());
                    count += 1.f;
                }
            }
        }
    }
    
    /**
     * assert pg 44 I3.
     * Every arc of NG, idle or saturated, is eps-proper.
     * all bipartite arcs should be eps-proper.
     * @param epsT
     * @return 
     */
    boolean integralFlowIsEpsProper(float epsT) {
        
        //NOTE that the flow is min-cost when eps < 1/6s
        //  and all integral arcs are "eps-proper"
        
        for (Map.Entry<Integer, Set<Integer>> entry : 
            forwardArcs.entrySet()) {
            
            Integer index1 = entry.getKey();
            
            for (Integer index2 : entry.getValue()) {
                
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
                float unitFlow = f.get(p);
                float cp = calcNetCost(p);
                if (unitFlow == 0) {
                    // idle, cp > -epsT
                    if (cp <= -epsT) {
                        return false;
                    }
                } else if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated
                    if (cp > epsT) {
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
        pRight[idx] = value;
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
}
