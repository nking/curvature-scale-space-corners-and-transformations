package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

public class CompleteSamplingVoidFinder extends AbstractVoidFinder {
    
    protected int skipInterval = 1;
    
    protected int xSortedIdxLo = -1;
    
    protected int xSortedIdxHi = -1;
    
    protected int ySortedIdxLo = -1;
    
    protected int ySortedIdxHi = -1;
    
    public void constructLogger() {
        this.log = Logger.getLogger(this.getClass().getName());
    }   
    
    @Override
    protected void findVoidsImpl() {

        findVoidsByPairBounds();
    }
    
    @Override
    protected void initializeVariables() {
        
        super.initializeVariables();
        
        this.xSortedIdxLo = 0;
        
        this.xSortedIdxHi = indexer.getNXY();
        
        this.ySortedIdxLo = 0;
        
        this.ySortedIdxHi = indexer.getNXY();
    }
    
    public void setXSortedIdxLo(int idx) {
        this.xSortedIdxLo = idx;
    }
    
    public void setXSortedIdxHi(int idx) {
        this.xSortedIdxHi = idx;
    }
    
    public void setYSortedIdxLo(int idx) {
        this.ySortedIdxLo = idx;
    }
    
    public void setYSortedIdxHi(int idx) {
        this.ySortedIdxHi = idx;
    }
    
    public void setSkipInterval(int skip) {
        this.skipInterval = skip;
    }

    /**
     * find voids by looking for other points within bounds of pairs
     * 
     */
    protected void findVoidsByPairBounds() {
                
        // look at every pair of points and keep only if their bounds do not include other points
        //          * (ii,jj)
        //
        //    *(i,j)

        for (int uSortedXIndex = xSortedIdxLo; uSortedXIndex < xSortedIdxHi; uSortedXIndex+=skipInterval) {
            
            for (int vSortedXIndex = ySortedIdxLo; vSortedXIndex < ySortedIdxHi; vSortedXIndex+=skipInterval) {
            
                if (uSortedXIndex == vSortedXIndex) {
                    continue;
                }
                
                processIndexedRegion(uSortedXIndex, vSortedXIndex);
            }
        }
    }
}
