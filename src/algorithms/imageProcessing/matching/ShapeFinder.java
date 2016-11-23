package algorithms.imageProcessing.matching;

import algorithms.util.PairIntArray;
import gnu.trove.map.TIntObjectMap;
import java.util.List;
import algorithms.imageProcessing.matching.PartialShapeMatcher.Result;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import gnu.trove.map.TObjectIntMap;
import java.util.Set;

/**
 * uses PartialShapeMatcher and search patterns to
 * find the best fitting match between a
 * group of adjacent segmented cells
 * to a template shape.  Any additional information
 * that helps limit the search should be done before
 * this stage and included or excluded in the
 * adjacency map.  
 *
 * Also note that the class uses internal cache which requires single
 * threaded use of this class.
 * 
 * @author nichole
 */
public class ShapeFinder {
    
    // partial shape matcher sampling distance
    private final int dp = 1;
    private final boolean useSameNSampl = false;

    private final PairIntArray bounds1;
    private final int[] ch1;
    private final float scale1;
    private final float sz1;
    private final List<Set<PairInt>> listOfSets2;
    private final List<OneDIntArray> listOfCH2s;
    private final float scale2;
    private final TObjectIntMap<int[]> keyIndexMap;
    private final TIntObjectMap<PairIntArray> indexBoundsMap;
    private final float intersectionLimit;
    
    /**
     * NOT READY FOR USE.
     * 
     * note that this is the data of one image in pyramid1 compared to the data
     * in one image of pyramid2 and it is assumed that one such pair will have
     * template and true match of nearly the same scale.  
     * On other words, this instance of ShapeFinder might not hold the true 
     * match at same scale, but another instance (holding other scles of the
     * pyramids) should if the true match is
     * present.
     * 
     * @param bounds1 clockwise ordered boundaries of template object.
     * @param ch1 1-D color histogram of template object
     * @param scale1 a bookkeeping number for the octave downsampling scale factor
     * with respect to the complete pyramid1.
     * @param listOfSets2 labeled point sets of dataset2 to be searched for template
     * object.listOfCH2s the color histograms for listOfSets2
     * @param scale2 a bookkeeping number for the octave downsampling scale factor
     * with respect to the complete pyramid2.
     * @param keyIndexMap an in/out variable for use in combination wtth
     * indexBoundsMap to cache the clockwise boundaries of aggregated adjacent
     * labeled cells.
     * The key = sorted indexes of listOfSets2, value = the lookup key for
     * indexBoundsMap,
     * @param indexBoundsMap an in/out variable for use in combination wtth
     * keyIndexMap to cache the clockwise boundaries of aggregated adjacent
     * labeled cells.
     * The key = the lookup key from keyIndexMap, value = clockwise ordered
     * bounding points of the aggregated sets of labels from the key in
     * indexBoundsMap,
     */
    public ShapeFinder(PairIntArray bounds1, int[] ch1, float scale1, float sz1,
        List<Set<PairInt>> listOfSets2, List<OneDIntArray> listOfCH2s,
        float scale2, TObjectIntMap<int[]> keyIndexMap,
        TIntObjectMap<PairIntArray>  indexBoundsMap,
        float intersectionLimit) {
        
        this.bounds1 = bounds1;
        this.ch1 = ch1;
        this.scale1 = scale1;
        this.sz1 = sz1;
        this.listOfSets2 = listOfSets2;
        this.listOfCH2s = listOfCH2s;
        this.scale2 = scale2;
        this.keyIndexMap = keyIndexMap;
        this.indexBoundsMap = indexBoundsMap;
        this.intersectionLimit = intersectionLimit;
    }
    
    public static class ShapeFinderResult extends Result {
        
        protected float intersection;
        protected final PairIntArray bounds1;
        protected final PairIntArray bounds2;
        
        public ShapeFinderResult(int n1, int n2, int offset, PairIntArray bounds1,
            PairIntArray bounds2) {
            
            super(n1, n2, offset);
            
            this.bounds1 = bounds1;
            this.bounds2 = bounds2;
        }
        
        public ShapeFinderResult(int n1, int n2, int nOriginal, 
            int offset, PairIntArray bounds1,
            PairIntArray bounds2) {
            
            super(n1, n2, nOriginal, offset);
            
            this.bounds1 = bounds1;
            this.bounds2 = bounds2;
        }
    }
    
    /**
     * NOT READY FOR USE.
     * uses PartialShapeMatcher and search patterns to
     * find the best fitting match between a
     * group of adjacent segmented cells
     * to a template shape. 
     *
     * The runtime is dependent upon which of the 3 search
     * patterns is kept in the end, so that will be
     * filled in here after implementation and testing.
     *
     * @return
     */
    public ShapeFinderResult findAggregated() {

        /*
        some notes that will be updated when return to the search
        task.
        
        basically, have a global grid search,
        wrapping local FloydWarshall searches where search
            is the best fitting bounds of aggregated labels sets and
            their color histograms
        
        step (1) O(N_cells) + O(N_cells_in_bin * log_2(N_cells_in_bin)):
           - 2D bins, that is x and y of size sz1.
           - sort the cells within each bin by x centroid and ties by y centroid
        step (2) search pattern, 3 are presented.
            The bigger picture is a sliding window of the 2D bin
            where the core of the search and aggregation is done
            within the 2D bin, but the aggregation can continue
            into the next bin to the right, next bin above,
            or above and to the right.
            
           this is the part that is changing.
           would like to use Dijkstra's within each bin to find the best 
             matching aggregation of cells, but unfortunately, the 
             aggregation of cell boundaries is not inductive, that is
             two wrong steps of aggregation might have better results than
             two correct steps.
           therefore, using the more time consuming but more complete
           Floyd-Warshal search within each 2D bin.
           Note that this might be improvable in the future.               
        */
        
        throw new UnsupportedOperationException("not yet implemented");
    }

}
