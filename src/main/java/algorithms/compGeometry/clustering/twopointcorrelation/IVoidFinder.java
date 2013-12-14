package algorithms.compGeometry.clustering.twopointcorrelation;

public interface IVoidFinder {
    
    public void findVoids(DoubleAxisIndexer indexer) throws TwoPointVoidStatsException;
        
    public void processIndexedRegion(int xSortedIndex0, int xSortedIndex1);
    
    public void processIndexedPair(int idx0, int idx1);
    
    
    public float[] getTwoPointDensities();
    
    public float[] getTwoPointDensityErrors();
    
    public int[] getPoint1();
    
    public int[] getPoint2();
    
    public int getNumberOfTwoPointDensities();
    
    public void setSampling(VoidSampling sampling);
    
    public float[] calulateTwoPointDensityErrors(float[] densities, int[] point1Indexes, 
        int[] point2Indexes, float[] xp, float[] yp, float[] xpe, float[] ype);
    
    public void releaseLargeVariables();
    
    public long approximateMemoryUsed();
    
    public void setDebug(boolean setDebugToTrue);
}
