package algorithms.imageProcessing.features;

public class FeatureMatcherSettings {

    private boolean useNormalizedFeatures = true;
    private String debugTag = "";
    private boolean debug = false;
    private boolean startWithBinnedImages = true;
    private boolean overrideWithCannySegmentation = false;
    
    //NOTE: not sure this belongs here.  temporarily placed until experimenting
    // w/ methods has found best practices
    private boolean use2ndDerivCorners = false;

    public void setToOverrideWithCannySegmentation() {
        overrideWithCannySegmentation = true;
    }
    public boolean doOverrideWithCannySegmentation() {
        return overrideWithCannySegmentation;
    }
    
    public void setToUse2ndDerivCorners() {
        use2ndDerivCorners = true;
    }
    
    public boolean doUse2ndDerivCorners() {
        return use2ndDerivCorners;
    }
    
    public boolean useNormalizedFeatures() {
        return useNormalizedFeatures;
    }

    public void setUseNormalizedFeatures(boolean useNormalizedFeatures) {
        this.useNormalizedFeatures = useNormalizedFeatures;
    }

    public String getDebugTag() {
        return debugTag;
    }

    public void setDebugTag(String debugTag) {
        this.debugTag = debugTag;
    }

    public boolean debug() {
        return debug;
    }

    public void setDebug(boolean debug) {
        this.debug = debug;
    }

    public boolean startWithBinnedImages() {
        return startWithBinnedImages;
    }

    public void setStartWithBinnedImages(boolean useBinnedImaged) {
        this.startWithBinnedImages = useBinnedImaged;
    }

    public FeatureMatcherSettings copy() {
        
        FeatureMatcherSettings s = new FeatureMatcherSettings();
        s.setDebug(debug);
        s.setDebugTag(debugTag);
        s.setStartWithBinnedImages(startWithBinnedImages);
        s.setUseNormalizedFeatures(useNormalizedFeatures);
        s.overrideWithCannySegmentation = overrideWithCannySegmentation;
        s.use2ndDerivCorners = use2ndDerivCorners;

        return s;
    }
}
