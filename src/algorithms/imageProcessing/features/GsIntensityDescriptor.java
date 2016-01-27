package algorithms.imageProcessing.features;

import algorithms.misc.MiscMath;
import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class GsIntensityDescriptor implements IntensityDescriptor {

    public static float sentinel = Float.MIN_VALUE;
    
    protected final float[] grey;
    
    protected float sumSquaredError = Float.NaN;
    
    protected boolean hasBeenNormalized = false;
    
    /**
     * the index within array grey that the central pixel
     * value is stored in.
     */
    protected final int centralIndex;
    
    /**
     * the number of columns across (it is the same as the number of rows)
     */
    protected final int nColumns;
    
    protected final int[] upperHalfIndexes;
    
    protected final int[] lowerHalfIndexes;
    
    public GsIntensityDescriptor(float[] intensities, int centralPixelIndex,
        int theNumberOfColumns) {
        this.grey = intensities;
        this.centralIndex = centralPixelIndex;
        this.nColumns = theNumberOfColumns;
        this.upperHalfIndexes = calculateUpperHalfIndexes();
        this.lowerHalfIndexes = calculateLowerHalfIndexes();
    }
    
    /**
     * apply a normalization to pixel values to try to reduce the differences 
     * due to images of the same region due to lighting or perspective 
     * for example.
     * The method invoked a second time does not change the internal values.
     */
    @Override
    public void applyNormalization() {
        
        if (hasBeenNormalized) {
            return;
        }
        
        /*        
        histogram equalization at the pre-processing stage of the entire image
        can stretch the range of values over the available range, but for
        images containing a small intersection of content that might not be
        a helpful operation.
        
        corrections at the block level for illumination probably need to 
        be derived at a larger level with knowledge of the illumination
        source...
        */
        
        float[] meanAndStDev = MiscMath.getAvgAndStDevIgnoreForSentinel(grey, 
            grey.length, sentinel);
        
        for (int i = 0; i < grey.length; ++i) {
            grey[i] -= meanAndStDev[0];
        }
        
        hasBeenNormalized = true;
    }

    @Override
    public boolean isNormalized() {
        return hasBeenNormalized;
    }
    
    @Override
    public float calculateSSD(IDescriptor otherDesc) {
        
        if (otherDesc == null) {
            throw new IllegalArgumentException("otherDesc cannot be null");
        }
        
        if (!(otherDesc instanceof GsIntensityDescriptor)) {
            throw new IllegalArgumentException(
            "otherDesc has to be type GsIntensityDescriptor");
        }
        
        GsIntensityDescriptor other = (GsIntensityDescriptor)otherDesc;
        
        if (this.grey.length != other.grey.length) {
            throw new IllegalArgumentException(
            "this and other arrays must have the same lengths");
        }
         
        float ssd = MiscMath.calculateSSD(grey, other.grey, sentinel);
                
        return ssd;
    }
    
    @Override
    public float calculateCosineSimilarity(IDescriptor otherDesc) {
        
        if (otherDesc == null) {
            throw new IllegalArgumentException("otherDesc cannot be null");
        }
        
        if (!(otherDesc instanceof GsIntensityDescriptor)) {
            throw new IllegalArgumentException(
            "otherDesc has to be type GsIntensityDescriptor");
        }
        
        GsIntensityDescriptor other = (GsIntensityDescriptor)otherDesc;
        
        if (this.grey.length != other.grey.length) {
            throw new IllegalArgumentException(
            "this and other arrays must have the same lengths");
        }
         
        float cSim = MiscMath.calculateCosineSimilarity(grey, other.grey, sentinel);
                
        return cSim;
    }
    
    /**
     * Determine the sum squared error within this descriptor using 
     * auto-correlation and the assumption that the value at the middle index 
     * is the value from the original central pixel.
     * Note that the value is persisted after one calculation.  Any use
     * of normalization should happen before this is first invoked.
     * (A function with a force calculation argument can be made if necessary though).
     * @return 
     */
    @Override
    public float sumSquaredError() {
        
        if (!Float.isNaN(sumSquaredError)) {
            return sumSquaredError;
        }
        
        float vc = grey[centralIndex];
        
        if (vc == sentinel) {
            throw new IllegalStateException(
            "ERROR: the central value for the array is somehow sentinel");
        }
        
        float sqErr = MiscMath.sumSquaredError(grey, sentinel, centralIndex);
            
        this.sumSquaredError = sqErr;
        
        return sumSquaredError;
    }
    
    /**
     * get the indexes of the array which represent the upper half of the
     * descriptor, that is, the half which are in the direction of the
     * orientation.
     * @return 
     */
    private int[] calculateUpperHalfIndexes() {
        
        /*
        for example:
             int cellDimension = 2;
             int range0 = 6;

         [5] (-6, 4)+2,+2  [11] (-4,-6)+2,+2  [17] (-2,-6)+2,+2  [23] (0,-6)+2,+2  [29] (2,-6)+2,+2  [35] (4,-6)+2,+2
         [4] (-6, 2)+2,+2  [10] (-4,-6)+2,+2  [16] (-2,-6)+2,+2  [22] (0,-6)+2,+2  [28] (2,-6)+2,+2  [34] (4,-6)+2,+2
         [3] (-6, 0)+2,+2  [9]  (-4,-6)+2,+2  [15] (-2,-6)+2,+2  [21] (0,-6)+2,+2  [27] (2,-6)+2,+2  [33] (4,-6)+2,+2
         [2] (-6,-2)+2,+2  [8]  (-4,-6)+2,+2  [14] (-2,-6)+2,+2  [20] (0,-6)+2,+2  [26] (2,-6)+2,+2  [32] (4,-6)+2,+2
         [1] (-6,-4)+2,+2  [7]  (-4,-6)+2,+2  [13] (-2,-6)+2,+2  [19] (0,-6)+2,+2  [25] (2,-6)+2,+2  [31] (4,-6)+2,+2
         [0] (-6,-6)+2,+2  [6]  (-4,-6)+2,+2  [12] (-2,-6)+2,+2  [18] (0,-6)+2,+2  [24] (2,-6)+2,+2  [30] (4,-6)+2,+2

         centralPixelIndex=21
         top half are indexes     3,4,5, 9,10,11, 15,16,17, 21,22,23, 27,28,29, 33,34,35
         bottom half are indexes  0,1,2, 6,7,8,   12,13,14, 18,19,20, 24,25,26, 30,31,32
        */
        
        int n = grey.length;
        int nDiv = n/nColumns;
        
        int[] indexes = new int[n];
        
        int startIdx = nColumns/2;
        
        int count = 0;
        for (int i = 0; i < nDiv; ++i) {
            int idx = startIdx + i*nDiv;
            for (int j = 0; j < nColumns; ++j) {
                indexes[count] = idx + j;
                count++;
            }
        }
        
        return indexes;
    }
    
    /**
     * get the indexes of the array which represent the upper half of the
     * descriptor, that is, the half which are in the direction of the
     * orientation.
     * @return 
     */
    public int[] getUpperHalfIndexes() {
        return upperHalfIndexes;
    }
    
    /**
     * get the indexes of the array which represent the upper half of the
     * descriptor, that is, the half which are in the direction of the
     * orientation.
     * @return 
     */
    private int[] calculateLowerHalfIndexes() {
        
        int n = grey.length;
        int nDiv = n/nColumns;
        
        int[] indexes = new int[n];
        
        int startIdx = 0;
        
        int count = 0;
        for (int i = 0; i < nDiv; ++i) {
            int idx = startIdx + i*nDiv;
            for (int j = 0; j < nColumns; ++j) {
                indexes[count] = idx + j;
                count++;
            }
        }
        
        return indexes;
    }

    /**
     * get the indexes of the array which represent the upper half of the
     * descriptor, that is, the half which are in the direction of the
     * orientation.
     * @return 
     */
    public int[] getLowerHalfIndexes() {
        return lowerHalfIndexes;
    }

    protected float[] getInternalArrayCopy() {
        return Arrays.copyOf(grey, grey.length);
    }

    @Override
    public int getCentralIndex() {
        return centralIndex;
    }

    @Override
    public String toString() {
        return Arrays.toString(grey);
    }
    
}
