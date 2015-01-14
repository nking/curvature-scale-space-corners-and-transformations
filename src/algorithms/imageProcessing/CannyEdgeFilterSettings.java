package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class CannyEdgeFilterSettings {
    
    private boolean useOutdoorMode = false;
    
    private boolean doNotNormalizeByHistogram = false;
    
    private boolean useLineDrawingMode = false;
    
    private boolean overrideHighThreshold = false;
    
    private int[] shrinkToSize = null;
    
    /**
     * field only used if overrideHighThreshold is true.
     */
    private float highThreshold = 0.0f;
    
    public void setOverrideHighThreshold(float value) {
        
        overrideHighThreshold = true;
        
        highThreshold = value;
    }
    
    public boolean getOverrideHighThreshold() {
        return overrideHighThreshold;
    }

    public void setUseOutdoorMode() {
        this.useOutdoorMode = true;
    }

    public void setDoNotNormalizeByHistogram() {
        this.doNotNormalizeByHistogram = true;
    }

    public void setUseLineDrawingMode() {
        this.useLineDrawingMode = true;
    }

    /**
     * @return the highThreshold
     */
    public float getHighThreshold() {
        return highThreshold;
    }
    
    public boolean getUseOutdoorMode() {
        return useOutdoorMode;
    }
    
    public boolean getDoNotNormalizeByHistogram() {
        return doNotNormalizeByHistogram;
    }
    
    public boolean getUseLineDrawingMode() {
        return useLineDrawingMode;
    }
    
    public void setShrinkToSize(int xOffset, int yOffset, int width, int height) {
        shrinkToSize = new int[]{xOffset, yOffset, width, height};
    }
    
    public int[] getShrinkToSize() {
        return shrinkToSize;
    }
    
}
