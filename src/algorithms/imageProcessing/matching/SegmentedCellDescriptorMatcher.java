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
    private final int[][] keypoints1;
    private final int[][] keypoints2;
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
     * @param keypoints1 
     * @param keypoints2 
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
        int[][] keypoints1, int[][] keypoints2,
        final double[] orientations1, final double[] orientations2,
        final List<Set<PairInt>> segmentedCells1, 
        final List<Set<PairInt>> segmentedCells2,
        final List<Set<PairInt>> medialAxis1, 
        final List<Set<PairInt>> medialAxis2,
        final RotatedOffsets rotatedOffsets) {
        
        assert(orientations1.length == keypoints1.length);
        assert(orientations2.length == keypoints2.length);
        
        this.img1 = img1;
        this.img2 = img2;
        this.gradient1 = gradient1;
        this.gradient2 = gradient2;
        this.keypoints1 = keypoints1;
        this.keypoints2 = keypoints2;
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
    
    /*
    there are many different ways to match these points, but the goal of the
    class is to use the segmentation information to match them using the
    group information for descriptor orientation and association with the 
    group location.
    
    The current goal is to find image 1 points as a single set of shape points
    so this class data for image 1 may change... noting requirements
    and what exists before making that decision...
    
    use case 1 (possibly the main and only use case):
       image 1 and one set of points in it (the template shape points).
       image 2 to be searched for the image 1 template.  image 2 has
       many sets of points (segmented cells, possibly filtered to
       be a small list of sets of points).
       the template shape might have projection transformation in the
       second image, so the solution should not assume stereo images.
       (note that have seen that an assumption of rough euclidean was necessary
       in other previous matching code before an epipolar allowed solution.
       this refined set of template points might make the rough first euclidean
       transformation unneccesary).
    
       The creation of a descriptor:
           - determine the rotation:
             - for a segmented cell boundary point:
               use medial axis points and the orientation to
               determine the integer rotation that points
               inward in the shape and is tangential to the 
               keypoint.
             - for a keypoint interior to the segmented cell:
               no further rotation correction necessary, but would
               only want to compare the keypoint to another interior
               keypoint in the other image.
               may want to try the rotation and 180 from it, and dither around 
               the keypoint in rotation angle and possibly small dither in x,y too
           - make the descriptor:
             use the intensity features object to extract that point and rotation
       To compare two descriptors:
          - use the static feature comparison method in the intensity features object.
       The stats of the comparisons then need to be evaluated to create a 
       a correspondence list.
    
       pausing here to consider global search pattern and how to store best
          results in evaluable steps.
    
       (1) using individual point matching:
          comparisons are each point in one image compared to each in other
             image.
          if iterate over sets, each point might be in 2 sets,
             so the number of comparisons is 
             --> n_template_keypoints * 2 * n_img2_keypoints <--
                 the interior point comparisons include dithers so are more than
                 a factor of 2 actually.
          various approaches to determining what is a match among those results exist,
          some using best results and nearest neighbors
             and excluding results in which there is a nearest neighbor
             within some distance of best answer.
          ==> this is the fastest solution.
              note that it only uses the group properties for
              orientation directions essentially.
          ransac helps to determine best consistent points within the implied
             geometry of a correspondence list.
          ** note that this approach should be less sensitive to scale differences
             as long as both objects (template and the true search match) are
             larger than the area of a descriptor at least...
             
       (2) OR ordered and aggregated (matches group points to another group's points:
          if use the keypoints on the boundaries only:
             - could make a local search just like the partial shape matcher
               but with color descriptors instead of chords and only using the
               keypoints.
               the runtime complexity is (review the code for this) 
                   probably m * n (where m and n are the number of points comparing).
             - that local search has to be wrapped in a pattern of aggregation
               combinations over segmented cells to find the best match.
               (did that with ShapeFinder, but using all boundary points and also
                using small (fast) color histograms).
               used floyd warshall pattern to be semi-complete.  that search runtime
                 is V^3 where V is the number of segmented cells in image 2.
               note that the same global search pattern but using a local 
               partial shape matcher style search but with hsv oriented 
               descriptors would make the local search much faster.
               - expecting that the inner keypoints if they exist, will be 
                 useful for either further evaluation of top results,
                 or further evaluation at the more detailed stage of the 
                 matching of individual aggregated cells with the template.
               - the adjacency properties and the finite dimensions of the template
                 shape are used to limit each search... can see need to partition
                 the image x,y space into 2D bins that are each searched with
                 a Floyd Warshall search (or other global search pattern).
                 
                 let V be the number of segmented cells in a bin plus extended area.
                 
                 for each bin, runtime is V^3 * (local matching stage)
                
                 the adjacency associations reduce the per bin runtime to be less
                     than V^3 for the outer runtime usually. 
                     can estimate min and max by assoc=1 for
                     each cell and assoc = all for each cell.
                (NOTE, the V^3 from Floyd Warshall pattern could be reduced using tabu, 
                scatter, or another pattern which is faster, but less definite)
    
                The local matching stage would be the descriptors matched using
                information present in their ordering (== a partial shape
                matcher tailored for keypoint descriptors rather than chords).  
                need to look at partial
                shape matching runtime, but think that was roughly m * n.
                note that the modified partial matcher has two different strict
                models - one in which the spacing between points is the same
                and the other being one in which the number of points in the
                shapes is the same.  so modifying the partial shape matcher
                to work with unevenly sampled keypoints does not look like it
                will be a faster solution after all...interpolation makes more
                points or exclusion reduces the points in a manner that is not
                clearly exlcuding the worst matchable points...then consider that
                the number of comparisons per point is not the single property of
                the chord difference, it's the descriptor size number of comparisons
                (which is 16 at smallest...).
    */
    
}
