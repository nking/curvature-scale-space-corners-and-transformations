package algorithms.graphs;

import algorithms.util.PairInt;
import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class CustomWatershedNode {

    protected PairInt key = null;
    
    /**
     * flag indicating whether node resolution has occurred (if true,
     * can just use first and only result in outgoing).
     */
    protected boolean resolved = false;

    protected PairInt[] outgoing = null;
    
    protected int nOutgoing = 0;

    /**
     * construct a node with a start capacity for edges of edgeCapacity.
     * @param nodeLocation
     */
    public CustomWatershedNode(PairInt nodeLocation) {
        this.key = nodeLocation;
        outgoing = new PairInt[1];
    }
    
    /**
     * construct a node with a start capacity for edges of edgeCapacity.
     * @param nodeLocation
     * @param edgeCapacity
     */
    public CustomWatershedNode(PairInt nodeLocation, int edgeCapacity) {
        this.key = nodeLocation;
        outgoing = new PairInt[edgeCapacity];
    }

    public void insertOutgoing(PairInt nodeLocation) {

        if ((nOutgoing + 1) > outgoing.length) {
            int expand = 2*nOutgoing;
            outgoing = Arrays.copyOf(outgoing, expand);
        }
        outgoing[nOutgoing] = nodeLocation;
        nOutgoing++;
    }

    public void setToResolved(PairInt nodeLocation) {
        
        resolved = true;
        
        outgoing[0] = nodeLocation;
        
        nOutgoing = 1;
        
        //TODO: compress if needed... there are at most 8 items in outgoing
    }
    
    public boolean isResolved() {
        return resolved;
    }
    
    public PairInt getResolved() {
        
        if (resolved) {
            return outgoing[0];
        }
        
        return null;
    }
}
