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
     * Note that if flag has been set to resolved, the single outgoing
     * connection may hold a value that it is not directly connected too,
     * but is the representative of its contiguous group.
     */
    protected boolean resolved = false;

    protected PairInt[] outgoing = null;
    
    protected int nOutgoing = 0;

    /**
     * construct a node with a start capacity for edges of edgeCapacity.
     * @param nodeLocation
     */
    public CustomWatershedNode(PairInt nodeLocation) {
        
        if (nodeLocation == null) {
            throw new IllegalStateException("nodeLocation cannot be null");
        }
        
        this.key = nodeLocation;
        outgoing = new PairInt[1];
    }
    
    /**
     * construct a node with a start capacity for edges of edgeCapacity.
     * @param nodeLocation
     * @param edgeCapacity
     */
    public CustomWatershedNode(PairInt nodeLocation, int edgeCapacity) {
        
        if (nodeLocation == null) {
            throw new IllegalStateException("nodeLocation cannot be null");
        }
        
        this.key = nodeLocation;
        outgoing = new PairInt[edgeCapacity];
    }

    public void insertOutgoing(PairInt nodeLocation) {
        
        if (resolved) {
            throw new IllegalStateException("this node has been set to resolved");
        }

        if ((nOutgoing + 1) > outgoing.length) {
            int expand = 2*nOutgoing;
            outgoing = Arrays.copyOf(outgoing, expand);
        }
        outgoing[nOutgoing] = nodeLocation;
        
        nOutgoing++;
    }

    public void setToResolved(PairInt nodeLocation) {
        
        if (resolved) {
            throw new IllegalStateException("this node has already been set to resolved");
        }
        
        if (nodeLocation == null) {
            throw new IllegalStateException("nodeLocation cannot be null");
        }
        
        resolved = true;
        
        if (nOutgoing == 0) {
            outgoing = new PairInt[]{nodeLocation};
        } else {
            outgoing[0] = nodeLocation;
            //TODO: compress if needed... there are at most 8 items in outgoing
        }
        
        nOutgoing = 1;        
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
    
    public int getConnectedNumber() {
        return nOutgoing;
    }
    
    public PairInt get(int nodeNumber) {
        
        if (nodeNumber < 0 || nodeNumber > (nOutgoing - 1)) {
            throw new IllegalArgumentException("nodeNumber is out of bounds");
        }
        
        return outgoing[nodeNumber];
    }

    public void reset(int nodeNumber, PairInt nodeValue) {
        
        if (nodeNumber < 0 || nodeNumber > (nOutgoing - 1)) {
            throw new IllegalArgumentException("nodeNumber is out of bounds");
        }
        
        if (nodeValue == null) {
            throw new IllegalArgumentException("nodeValue cannot be null");
        }
        
        outgoing[nodeNumber] = nodeValue;
    }
    
}
