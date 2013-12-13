package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

public class CompleteSamplingVoidFinder extends AbstractVoidFinder {
    
    public void constructLogger() {
        this.log = Logger.getLogger(this.getClass().getName());
    }   
    
    @Override
    protected void findVoidsImpl() {

        findVoidsByPairBounds();
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

        for (int uSortedXIndex = 0; uSortedXIndex < indexer.getNumberOfPoints(); uSortedXIndex++) {
            
            for (int vSortedXIndex = (uSortedXIndex + 1); vSortedXIndex < indexer.getNumberOfPoints(); vSortedXIndex++) {
                
                processIndexedRegion(uSortedXIndex, vSortedXIndex);
            }
        }
    }
}
