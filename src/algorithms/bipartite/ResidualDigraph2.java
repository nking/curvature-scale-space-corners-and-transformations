package algorithms.bipartite;

import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

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
    private TIntObjectMap<TIntSet> forwardLinksRM
        = new TIntObjectHashMap<TIntSet>();

    /**
     * links Y to X (that is, right to left). These are "saturated" 
     * arcs, f=1, in the flow network N_G. They correspond to "married" 
     * in the matched M graph.
     */
    private TIntIntMap backwardLinksRM = new TIntIntHashMap();

    /**
     * source to Left links are represented by this
     * set of Left indexes.  These are "saturated" 
     * Left vertexes at initialization from being
     * connected to a matched node.
     */
    private TIntSet forwardLinksSourceRM =
        new TIntHashSet();
    
    /**
     * Left to source links are represented by this set of
     * Left vertexes.  These are "idle" Left vertexes.
     */
    private TIntSet backwardLinksSourceRM =
        new TIntHashSet();
    
    /**
     * Right to sink links are represented by this set of
     * right vertexes.  These are "saturated" Right vertexes
     * at initialization, from being connected to a matched node.
     */
    private TIntSet forwardLinksSinkRM =
        new TIntHashSet();
    
    /**
     * sink to Right links are represented by this set of
     * Right vertexes.  These are "idle" Right vertexes.
     */
    private TIntSet backwardLinksSinkRM =
        new TIntHashSet();
    
    public ResidualDigraph2(FlowNetwork gFlow) {
    
        this.nLeft = gFlow.getNLeft();
        this.nRight = gFlow.getNRight();
        this.sourceNode = gFlow.getSourceNode();
        this.sinkNode = gFlow.getSinkNode();
       
        TIntObjectIterator<TIntSet> iter = gFlow.getForwardArcs()
            .iterator();
        
        for (int i = gFlow.getForwardArcs().size(); i-- > 0;) {
            iter.advance();
            int idx1 = iter.key();
            TIntIterator iter2 = iter.value().iterator();
            while (iter2.hasNext()) {
                int idx2 = iter2.next();
                PairInt p = new PairInt(idx1, idx2);
                float unitFlow = gFlow.getFlow().get(p);
                if (unitFlow == 0) {
                    // idle
                    TIntSet indexes2 = forwardLinksRM.get(idx1);
                    if (indexes2 == null) {
                        indexes2 = new TIntHashSet();
                        forwardLinksRM.put(idx1, indexes2);
                    }
                    indexes2.add(idx2);
                } else if (Math.abs(unitFlow - 1) < 0.01f) {
                    // saturated
                    backwardLinksRM.put(idx2, idx1);
                }
            }
        }
        
        // see Figure 7.2 pg 49
        
        TIntIterator iter2 = gFlow.getSourceForwardArcs().iterator();
        while (iter2.hasNext()) {
            int xNode = iter2.next();
                    
            float unitFlow = gFlow.getFlow(sourceNode, xNode);                
            if (unitFlow == 0) {
                // idle
                backwardLinksSourceRM.add(xNode);
            } else if (Math.abs(unitFlow - 1) < 0.01f) {
                // saturated
                forwardLinksSourceRM.add(xNode);
            }
        }

        iter2 = gFlow.getSinkForwardArcs().iterator();
        while (iter2.hasNext()) {
            int yNode = iter2.next();
            float unitFlow = gFlow.getFlow(yNode, sinkNode);                
            if (unitFlow == 0) {
                // idle
                backwardLinksSinkRM.add(yNode);
            } else if (Math.abs(unitFlow - 1) < 0.01f) {
                // saturated
                forwardLinksSinkRM.add(yNode);
            }
        }
    }
    
    public int countOfForwardBipartiteLinks() {
        int n = 0;
        
        TIntObjectIterator<TIntSet> iter = forwardLinksRM
            .iterator();
        
        for (int i = forwardLinksRM.size(); i-- > 0;) {
            iter.advance();
            TIntSet set = iter.value();
            n += set.size();
        }
        return n;
    }
  
    /**
     * @return the forwardLinksRM
     */
    public TIntObjectMap<TIntSet> getForwardLinksRM() {
        return forwardLinksRM;
    }

    /**
     * key is Right (==Y) index node, and value is
     * Left (==X) index node.
     * @return the backwardLinksRM
     */
    public TIntIntMap getBackwardLinksRM() {
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
    public TIntSet getForwardLinksSourceRM() {
        return forwardLinksSourceRM;
    }

    /**
     * @return the backwardLinksSourceRM
     */
    public TIntSet getBackwardLinksSourceRM() {
        return backwardLinksSourceRM;
    }

    /**
     * the arcs that are initialized from being connected
     * to matched nodes, that is Right nodes that
     * are in "saturated" arcs.
     * @return the forwardLinksSinkRM
     */
    public TIntSet getForwardLinksSinkRM() {
        return forwardLinksSinkRM;
    }

    /**
     * @return the backwardLinksSinkRM
     */
    public TIntSet getBackwardLinksSinkRM() {
        return backwardLinksSinkRM;
    }

}
