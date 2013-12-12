package algorithms.compGeometry.clustering.twopointcorrelation;

public interface IVoidFinder {
    
    public void findVoids(DoubleAxisIndexer indexer) throws TwoPointVoidStatsException;
    
    public void processIndexedRegion(int xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi, boolean useCompleteSampling);
    
    public void processIndexedPair(int regionIndex0, int regionIndex1);
    
    
    public float[] getTwoPointDensities();
    
    public float[] getTwoPointDensityErrors();
    
    public int[] getPoint1();
    
    public int[] getPoint2();
    
    public int getNumberOfTwoPointDensities();
    
    public void setSampling(VoidSampling sampling);
    
    public float[] calulateTwoPointDensityErrors(float[] densities, int[] point1Indexes, 
        int[] point2Indexes, float[] xp, float[] yp, float[] xpe, float[] ype);
            
    
    public long approximateMemoryUsed();
    
    public void setDebug(boolean setDebugToTrue);
}
