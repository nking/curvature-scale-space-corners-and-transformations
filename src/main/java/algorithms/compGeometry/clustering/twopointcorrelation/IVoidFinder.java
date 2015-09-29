package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 * interface for classes that find the voids (that is, space between points) 
 * in the dataset
 * 
 * @author nichole
 */
public interface IVoidFinder {
    
    /**
     *
     * @param indexer
     * @throws TwoPointVoidStatsException
     */
    public void findVoids(AxisIndexer indexer) throws TwoPointVoidStatsException;
        
    /**
     *
     * @param xSortedIndex0
     * @param xSortedIndex1
     */
    public void processIndexedRegion(int xSortedIndex0, int xSortedIndex1);
    
    /**
     *
     * @param idx0
     * @param idx1
     */
    public void processIndexedPair(int idx0, int idx1);

    /**
     *
     * @return
     */
    public float[] getTwoPointDensities();
    
    /**
     *
     * @return
     */
    public float[] getTwoPointDensityErrors();
    
    /**
     *
     * @return
     */
    public int[] getPoint1();
    
    /**
     *
     * @return
     */
    public int[] getPoint2();
    
    /**
     *
     * @return
     */
    public int getNumberOfTwoPointDensities();
    
    /**
     *
     * @param sampling
     */
    public void setSampling(VoidSampling sampling);
    
    /**
     *
     * @param densities
     * @param point1Indexes
     * @param point2Indexes
     * @param xp
     * @param yp
     * @param xpe
     * @param ype
     * @return
     */
    public float[] calulateTwoPointDensityErrors(float[] densities, int[] point1Indexes, 
        int[] point2Indexes, float[] xp, float[] yp, float[] xpe, float[] ype);
        
    /**
     *
     * @return
     */
    public long approximateMemoryUsed();
    
    /**
     *
     * @param setDebugToTrue
     */
    public void setDebug(boolean setDebugToTrue);
}
