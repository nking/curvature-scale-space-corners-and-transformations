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
     * links X to Y (that is, left to right). These are "idle" 
     * arcs, f=0, in the flow network N_G. They correspond to 
     * "unmarried" in the matched M graph.
     */
    private Map<Integer, Set<Integer>> forwardLinksRM
        = new HashMap<Integer, Set<Integer>>();

    /**
     * links Y to X (that is, right to left). These are "saturated" 
     * arcs, f=1, in the flow network N_G. They correspond to "married" 
     * in the matched M graph.
     */
    private Map<Integer, Integer> backwardLinksRM
        = new HashMap<Integer, Integer>();

    /**
     * source to Left links are represented by this
     * set of Left indexes.  These are "saturated" 
     * Left vertexes at initialization from being
     * connected to a matched node.
     */
    private Set<Integer> forwardLinksSourceRM =
        new HashSet<Integer>();
    
    /**
     * Left to source links are represented by this set of
     * Left vertexes.  These are "idle" Left vertexes.
     */
    private Set<Integer> backwardLinksSourceRM =
        new HashSet<Integer>();
    
    /**
     * Right to sink links are represented by this set of
     * right vertexes.  These are "saturated" Right vertexes
     * at initialization, from being connected to a matched node.
     */
    private Set<Integer> forwardLinksSinkRM =
        new HashSet<Integer>();
    
    /**
     * sink to Right links are represented by this set of
     * Right vertexes.  These are "idle" Right vertexes.
     */
    private Set<Integer> backwardLinksSinkRM =
        new HashSet<Integer>();
    
    public ResidualDigraph2(FlowNetwork gFlow) {
    
        this.nLeft = gFlow.getNLeft();
        this.nRight = gFlow.getNRight();
        this.sourceNode = gFlow.getSourceNode();
        this.sinkNode = gFlow.getSinkNode();
       
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
        
        // see Figure 7.2 pg 49
        
        for (Integer xNode : gFlow.getSourceForwardArcs()) {
            float unitFlow = gFlow.getFlow(sourceNode, xNode.intValue());                
            if (unitFlow == 0) {
                // idle
                backwardLinksSourceRM.add(xNode);
            } else if (Math.abs(unitFlow - 1) < 0.01f) {
                // saturated
                forwardLinksSourceRM.add(xNode);
            }
        }

        for (Integer yNode : gFlow.getSinkForwardArcs()) {
            float unitFlow = gFlow.getFlow(yNode.intValue(), sinkNode);                
            if (unitFlow == 0) {
                // idle
                backwardLinksSinkRM.add(yNode);
            } else if (Math.abs(unitFlow - 1) < 0.01f) {
                // saturated
                forwardLinksSinkRM.add(yNode);
            }
        }
    }
  
    /**
     * @return the forwardLinksRM
     */
    public Map<Integer, Set<Integer>> getForwardLinksRM() {
        return forwardLinksRM;
    }

    /**
     * key is Right (==Y) index node, and value is
     * Left (==X) index node.
     * @return the backwardLinksRM
     */
    public Map<Integer, Integer> getBackwardLinksRM() {
        return backwardLinksRM;
    }

    /**
     * @return the sourceNode
     */
    public int getSourceNode() {
        return sourceNode;
    }

    /**
     * @return the sinkNode
     */
    public int getSinkNode() {
        return sinkNode;
    }
    
    /**
     * the arcs that are initialized from being connected
     * to matched nodes, that is Left nodes that
     * are in "saturated" arcs.
     * @return the forwardLinksSourceRM
     */
    public Set<Integer> getForwardLinksSourceRM() {
        return forwardLinksSourceRM;
    }

    /**
     * @return the backwardLinksSourceRM
     */
    public Set<Integer> getBackwardLinksSourceRM() {
        return backwardLinksSourceRM;
    }

    /**
     * the arcs that are initialized from being connected
     * to matched nodes, that is Right nodes that
     * are in "saturated" arcs.
     * @return the forwardLinksSinkRM
     */
    public Set<Integer> getForwardLinksSinkRM() {
        return forwardLinksSinkRM;
    }

    /**
     * @return the backwardLinksSinkRM
     */
    public Set<Integer> getBackwardLinksSinkRM() {
        return backwardLinksSinkRM;
    }

}
