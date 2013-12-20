package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

public class SubsetSamplingVoidFinder extends AbstractVoidFinder {
        
    protected int sortedIdxLo = -1;
    
    protected int sortedIdxHi = -1;
    
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
        
        // if not set by user
        if (sortedIdxLo == -1) {
            
            this.sortedIdxLo = 0;

            this.sortedIdxHi = indexer.getNXY();
        }
    }
    
    /**
     * set the low index of the range of indexes.  the indexes are w.r.t. the array
     * indexer.sortedXIndexes.
     * 
     * @param idx
     */
    public void setSortedIdxLo(int idx) {
        this.sortedIdxLo = idx;
    }
    
    /**
     * set the high index of the range of indexes.  the indexes are w.r.t. the array
     * indexer.sortedXIndexes.
     * 
     * @param idx
     */
    public void setSortedIdxHi(int idx) {
        this.sortedIdxHi = idx;
    }

    /**
     * find voids by excluding pairs of points that contain other points within their bounds
     */
    protected void findVoidsByPairBounds() {
       
        for (int uSortedXIndex = sortedIdxLo; uSortedXIndex < sortedIdxHi; uSortedXIndex++) {
            
            for (int vSortedXIndex = (uSortedXIndex + 1); vSortedXIndex < sortedIdxHi; vSortedXIndex++) {
                                
                processIndexedRegion(uSortedXIndex, vSortedXIndex);
            }
        }
    }
}
