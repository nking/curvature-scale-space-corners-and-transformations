package algorithms.imageProcessing;

public class FeatureMatcherSettings {

    private boolean useNormalizedFeatures = true;
    private String debugTag = "";
    private boolean debug = false;
    private boolean startWithBinnedImages = true;

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

        return s;
    }
}
