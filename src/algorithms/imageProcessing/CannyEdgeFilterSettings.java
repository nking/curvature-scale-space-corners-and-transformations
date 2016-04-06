package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class CannyEdgeFilterSettings {
        
    private boolean doNormalizeByHistEqualization = false;
    
    private boolean useLineDrawingMode = false;

    public void setToNormalizeByHistogram() {
        this.doNormalizeByHistEqualization = true;
    }

    public void setUseLineDrawingMode() {
        this.useLineDrawingMode = true;
    }
    
    public boolean getNormalizeByHistogram() {
        return doNormalizeByHistEqualization;
    }
    
    public boolean getUseLineDrawingMode() {
        return useLineDrawingMode;
    }
    
}
