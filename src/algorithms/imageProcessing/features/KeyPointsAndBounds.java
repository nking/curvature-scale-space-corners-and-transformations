package algorithms.imageProcessing.features;

import algorithms.imageProcessing.ImageSegmentation;
import algorithms.util.PairInt;
import java.util.List;
import java.util.Set;

/**
 * class containing groups of keypoints and details of the the bounding 
 * regions of each group.
 * 
 * @author nichole
 */
public class KeyPointsAndBounds {
    
    private final List<Set<PairInt>> keyPointGroups;
    
    private final ImageSegmentation.BoundingRegions boundingRegions;
    
    public KeyPointsAndBounds(List<Set<PairInt>> theKeyPointGroups,
        ImageSegmentation.BoundingRegions theBoundingRegions) {
        
        this.keyPointGroups = theKeyPointGroups;
        this.boundingRegions = theBoundingRegions;
    }

    /**
     * @return the keyPointGroups
     */
    public List<Set<PairInt>> getKeyPointGroups() {
        return keyPointGroups;
    }

    /**
     * @return the boundingRegions
     */
    public ImageSegmentation.BoundingRegions getBoundingRegions() {
        return boundingRegions;
    }
    
}
