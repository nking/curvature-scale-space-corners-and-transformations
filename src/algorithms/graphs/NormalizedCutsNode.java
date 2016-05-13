package algorithms.graphs;

/**
 *
 * @author nichole
 */
public class NormalizedCutsNode {
    
    private final int label;
    
    private int nCutsLabel = -1;
    
    public NormalizedCutsNode(int theLabel) {
        this.label = theLabel;
    }

    /**
     * @return the label
     */
    public int getLabel() {
        return label;
    }

    /**
     * @return the nCutsLabel
     */
    public int getNCutsLabel() {
        return nCutsLabel;
    }

    /**
     * @param nCutsLabel the nCutsLabel to set
     */
    public void setNCutsLabel(int nCutsLabel) {
        this.nCutsLabel = nCutsLabel;
    }
    
    
}
