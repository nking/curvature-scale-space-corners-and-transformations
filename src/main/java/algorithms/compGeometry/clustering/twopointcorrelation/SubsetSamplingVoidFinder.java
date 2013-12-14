package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

public class SubsetSamplingVoidFinder extends AbstractVoidFinder {
        
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
        
        // if not set by user
        if (xSortedIdxLo == -1) {
            
            this.xSortedIdxLo = 0;

            this.xSortedIdxHi = indexer.getNXY();

            this.ySortedIdxLo = 0;

            this.ySortedIdxHi = indexer.getNXY();
        }
    }
    
    /**
     * set the low index of the range of x indexes.  the indexes are w.r.t. the array
     * indexer.sortedXIndexes.
     * 
     * @param idx
     */
    public void setXSortedIdxLo(int idx) {
        this.xSortedIdxLo = idx;
    }
    
    /**
     * set the high index of the range of x indexes.  the indexes are w.r.t. the array
     * indexer.sortedXIndexes.
     * 
     * @param idx
     */
    public void setXSortedIdxHi(int idx) {
        this.xSortedIdxHi = idx;
    }
    
    /**
     * set the low index of the range of y indexes.  the indexes are w.r.t. the array
     * indexer.sortedYIndexes.
     * 
     * @param idx
     */
    public void setYSortedIdxLo(int idx) {
        this.ySortedIdxLo = idx;
    }
    
    /**
     * set the high index of the range of y indexes.  the indexes are w.r.t. the array
     * indexer.sortedYIndexes.
     * 
     * @param idx
     */
    public void setYSortedIdxHi(int idx) {
        this.ySortedIdxHi = idx;
    }

    /**
     * find voids by looking for other points within bounds of pairs
     * 
     */
    protected void findVoidsByPairBounds() {
       
        for (int uSortedXIndex = xSortedIdxLo; uSortedXIndex < xSortedIdxHi; uSortedXIndex++) {
            
            for (int vSortedXIndex = ySortedIdxLo; vSortedXIndex < ySortedIdxHi; vSortedXIndex++) {

                // this is an index w.r.t. indexer.sortedYIndexes, so either need to provide
                //   another implementation of processIndexedRegion that is expecting that the
                //   2nd index is w.r.t. indexer.sortedYIndexes or convert the index to one
                //   relative to indexer.sortedXIndexes.  the later should be easier to maintain.
                
                float yValue = indexer.getY()[vSortedXIndex];
                
                int idx2 = indexer.findSortedXIndexesForY(yValue);
                
                if (uSortedXIndex == idx2) {
                    continue;
                }
                                
                processIndexedRegion(uSortedXIndex, idx2);
            }
        }
    }
}
