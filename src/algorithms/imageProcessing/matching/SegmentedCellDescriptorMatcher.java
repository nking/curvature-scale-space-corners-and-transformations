package algorithms.imageProcessing.matching;

import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.features.IntensityClrFeatures2;
import algorithms.util.PairInt;
import java.util.List;
import java.util.Set;

/**
 * this class uses keypoints, their orientations,
 * and segmented cells to make a correspondence list.
 * 
 * @author nichole
 */
public class SegmentedCellDescriptorMatcher {
    
    private final ImageExt img1;
    private final ImageExt img2;
    private final GreyscaleImage gradient1;
    private final GreyscaleImage gradient2;
    private final double[] orientations1;
    private final double[] orientations2;
    private final List<Set<PairInt>> segmentedCells1;
    private final List<Set<PairInt>> segmentedCells2;
    private final List<Set<PairInt>> medialAxis1;
    private final List<Set<PairInt>> medialAxis2;
    private final RotatedOffsets rotatedOffsets;
    
    // depending upon other use, may externalize these, that is require them
    // to be given as arguments
    private final IntensityClrFeatures2 features1;
    private final IntensityClrFeatures2 features2;
    
    /**
     * note, RotatedOffsets caches the results of expensive transcendental
     * operations and is used in a shared manner, that is,
     * it is a singleton.  caution is needed in it's use because it 
     * contains unsynchronized unguarded cache variables that could be 
     * corrupted by multi-threaded use.  The design is meant to keep access
     * to it fast.
     * 
     * @param img1
     * @param img2
     * @param gradient1
     * @param gradient2
     * @param orientations1
     * @param orientations2
     * @param segmentedCells1
     * @param segmentedCells2
     * @param medialAxis1
     * @param medialAxis2
     * @param rotatedOffsets 
     */
    public SegmentedCellDescriptorMatcher(
        final ImageExt img1, final ImageExt img2,
        final GreyscaleImage gradient1, final GreyscaleImage gradient2,
        final double[] orientations1, final double[] orientations2,
        final List<Set<PairInt>> segmentedCells1, 
        final List<Set<PairInt>> segmentedCells2,
        final List<Set<PairInt>> medialAxis1, 
        final List<Set<PairInt>> medialAxis2,
        final RotatedOffsets rotatedOffsets) {
        
        this.img1 = img1;
        this.img2 = img2;
        this.gradient1 = gradient1;
        this.gradient2 = gradient2;
        this.orientations1 = orientations1;
        this.orientations2 = orientations2;
        this.segmentedCells1 = segmentedCells1;
        this.segmentedCells2 = segmentedCells2;
        this.rotatedOffsets = rotatedOffsets;
        this.medialAxis1 = medialAxis1;
        this.medialAxis2 = medialAxis2;
        
        this.features1 = new IntensityClrFeatures2(img1, 5, rotatedOffsets);
        this.features2 = new IntensityClrFeatures2(img2, 5, rotatedOffsets);
    }
    
    
}
