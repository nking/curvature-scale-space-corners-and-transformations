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
    public int getnCutsLabel() {
        return nCutsLabel;
    }

    /**
     * @param nCutsLabel the nCutsLabel to set
     */
    public void setnCutsLabel(int nCutsLabel) {
        this.nCutsLabel = nCutsLabel;
    }
    
    
}
