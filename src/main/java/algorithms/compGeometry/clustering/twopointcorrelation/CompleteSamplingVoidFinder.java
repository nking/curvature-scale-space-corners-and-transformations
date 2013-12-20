package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

/**
 * RT complexity for findVoids()
 * 
 * findVoids() is < O(N^2), but it invokes processIndexedRegion which can be O(1) to O(N-1).
 * best time: < O(N^2)
 * worse time: < O(N^3)
 * 
 * @author nichole
 */
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

        // this is < O(N^2), but processIndexedRegion can be O(1) to O(N-1)
        for (int uSortedXIndex = 0; uSortedXIndex < indexer.getNumberOfPoints(); uSortedXIndex++) {
            
            for (int vSortedXIndex = (uSortedXIndex + 1); vSortedXIndex < indexer.getNumberOfPoints(); vSortedXIndex++) {
                
                processIndexedRegion(uSortedXIndex, vSortedXIndex);
            }
        }
    }
}
