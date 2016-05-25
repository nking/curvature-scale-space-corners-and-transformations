package algorithms.bipartite;

import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
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
    int sourceNode = -1;
    
    /**
     * <pre>
     * -|
     * </pre>
     */
    int sinkNode = -1;

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
    * Note that the source node has x index nLeft and
    * that the sink node has y index nRight and both
    * are present in forwardArcs.
    */
    Map<Integer, Set<Integer>> forwardArcs
        = new HashMap<Integer, Set<Integer>>();

    public FlowNetwork(Graph g, Map<Integer, Integer> m) {

        this.nLeft = g.getNLeft();
        this.nRight = g.getNRight();
        pLeft = new float[nLeft + 1];
        pRight = new float[nRight + 1];

        for (Map.Entry<PairInt, Integer> entry : g.getEdgeWeights().entrySet()) {

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

            if (m.containsKey(index1) && m.get(index1).equals(index2)) {
                // create saturated arc, f(v,w)=1, net cost cp <= 0
                f.put(p, Float.valueOf(1));
            } else {
                // create idle arc, f(v,w)=0, net cost cp >= 0
                f.put(p, Float.valueOf(0));
            }
        }

        assert (maxC > 1);

        //Arrays.fill(pRight, maxC);
        // --- handle source and sink nodes ---
        this.sourceNode = nLeft;
        this.sinkNode = nRight;
        
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
        
        /*
        TODO: need to store the left and right stubs that
           are the result of an iteration of method refine(...)
        */
    }
    
    public int getMaxC() {
        return maxC;
    }

    boolean integralFlowIsEpsProper(float epsT) {
        
        //NOTE that the flow is min-cost when eps < 1/6s
        //  and all integral ars are "eps-proper"
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
                if (unitFlow == 0) {
                    // idle, cp > -epsT
                    if (cp <= -epsT) {
                        return false;
                    }
                } else if (unitFlow == 1.f) {
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
            if (index1.intValue() == sourceNode) {
                continue;
            }
            for (Integer index2 : entry.getValue()) {
                if (index2.intValue() == sinkNode) {
                    continue;
                }
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
        
        return sum;
    }
    
    public float calcNetCost(int u, int v) {

        PairInt p = new PairInt(u, v);

        int cost = c.get(p);
        float pdX = pLeft[u];
        float pdY = pRight[v];

        float cp = cost - pdX + pdY;

        return cp;
    }

    public float calcNetCost(PairInt p) {

        int cost = c.get(p);
        float pdX = pLeft[p.getX()];
        float pdY = pRight[p.getY()];

        float cp = cost - pdX + pdY;

        return cp;
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
    
    /**
     * per unit flow costs are the same as the edges given in G. 
     * key = index in left (==X==v) and index in right (==Y==w), 
     * value = cost of the arc. might
     * change to use sparse matrix format
     */
    Map<PairInt, Integer> c
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
    Map<PairInt, Float> f
        = new HashMap<PairInt, Float>();

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
}
