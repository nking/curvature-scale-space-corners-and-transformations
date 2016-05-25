package algorithms.bipartite;

import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * a residual digraph for use within refine method in
 * MinCostUnbalancedAssignment.java
 * 
 * @author nichole
 */
public class ResidualDigraph2 {
    
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
     * links X to Y (that is, left to right). These are "idle" arcs, f=0, in the
     * flow network N_G. They correspond to "unmarried" in the matched M graph.
     */
    private Map<Integer, Set<Integer>> forwardLinksRM
        = new HashMap<Integer, Set<Integer>>();

    /**
     * links Y to X (that is, right to left). These are "saturated" arcs, f=1,
     * in the flow network N_G. They correspond to "married" in the matched M
     * graph.
     */
    private Map<Integer, Integer> backwardLinksRM
        = new HashMap<Integer, Integer>();

    public ResidualDigraph2(FlowNetwork gFlow) {
    
        this.nLeft = gFlow.getNLeft();
        this.nRight = gFlow.getNRight();
        this.sourceNode = gFlow.getSourceNode();
        this.sinkNode = gFlow.getSinkNode();
        
        // this includes source and sink node arcs too
        
        for (Map.Entry<Integer, Set<Integer>> entry : 
            gFlow.getForwardArcs().entrySet()) {
            
            Integer index1 = entry.getKey();
            
            for (Integer index2 : entry.getValue()) {
               
                PairInt p = new PairInt(index1.intValue(), index2.intValue());
                float unitFlow = gFlow.getFlow().get(p);
                if (unitFlow == 0) {
                    // idle
                    Set<Integer> indexes2 = forwardLinksRM.get(index1);
                    if (indexes2 == null) {
                        indexes2 = new HashSet<Integer>();
                        forwardLinksRM.put(index1, indexes2);
                    }
                    indexes2.add(index2);
                } else if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated
                    backwardLinksRM.put(index2, index1);
                }
            }
        }
        
    }
    
    /*
    augmenting paths start at a surplus node and end at a 
    deficit node, and may visit source or sink nodes
    in between.  no node can be visited more than once 
    in an augmenting path.
       - augmenting path involves reversing link idle or
         saturated status.
    
    lengths of links:
       - quantization into units of eps is used for net cost:
         for a forward link: lp(v->w) = Math.ceil(cp(v, w)/eps)
         for a backward link: lp(w->v) = 1 - Math.ceil(cp(v, w)/eps)
         both results are >= 0
    
    raising prices:
       during main loop in refine, raising by eps the
         price gFlow.pRight at a right node in gFlow
         is complemented by 
           lowering by a value of 1 the length of any link
           in the residul digraph that leaves v and
           raising by a value of 1 the length of any link
           that enters v.
           (can see that is a + or - to cp in section above
            which would presumably be applied to the 
            flow network f??  looks like the authors only
            apply the change dynamically?)
          proof 7-6 on pg 52...
       so looks like link lenths might need to be a variable
       here
    
       need to add assert I5 on pg 52 to this class
    
       example of price changes on pg 53 and 54
    */
    
    /**
     * @return the forwardLinksRM
     */
    public Map<Integer, Set<Integer>> getForwardLinksRM() {
        return forwardLinksRM;
    }

    /**
     * @return the backwardLinksRM
     */
    public Map<Integer, Integer> getBackwardLinksRM() {
        return backwardLinksRM;
    }
}
