package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

/**
 * samples only the non-zero cells in a dataset to attempt 
 * 
 * @author nichole
 *
 */
public class SparseSamplingVoidFinder extends AbstractVoidFinder {
        
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

    /**
     * find voids by looking for other points within bounds of pairs
     * 
     */
    protected void findVoidsByPairBounds() {
       
        for (int uSortedXIndex = xSortedIdxLo; uSortedXIndex < xSortedIdxHi; uSortedXIndex++) {
            
            for (int vSortedXIndex = ySortedIdxLo; vSortedXIndex < ySortedIdxHi; vSortedXIndex++) {
            
                if (uSortedXIndex == vSortedXIndex) {
                    continue;
                }
                
                processIndexedRegion(uSortedXIndex, vSortedXIndex);
            }
        }
    }
}
