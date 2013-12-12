package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

public class DivideAndConquerVoidFinder extends AbstractVoidFinder {

    @Override
    public void constructLogger() {
        this.log = Logger.getLogger(this.getClass().getName());
    }

    @Override
    protected void findVoidsImpl() {
        
        findVoids(0, indexer.nXY - 1, 0, indexer.nXY - 1);
    }
    

    /**
     * A divide and conquer approach to finding the rectangular areas
     * containing only two points.  it's a recursion with pattern
     * T(n) = 4T(n/2) + n  so the runtime is O(n^2).  
     * 
     * Note: It does not completely sample every pair of points.
     *
     * @param xIndexLo
     * @param xIndexHi
     * @param yIndexLo
     * @param yIndexHi
     */
    protected void findVoids(int xIndexLo, int xIndexHi,
        int yIndexLo, int yIndexHi) {
        
        boolean useCompleteSampling = (sampling.ordinal() == VoidSampling.COMPLETE.ordinal());
        
                                                                                 // cost     number of times
        if ((xIndexLo < xIndexHi) && (yIndexLo < yIndexHi)) {                    //

            int xIndexMid = (xIndexLo + xIndexHi)/2;                             //

            int yIndexMid = (yIndexLo + yIndexHi)/2;                             //

            findVoids(xIndexLo, xIndexMid, yIndexLo, yIndexMid);           // c4           N/2
            findVoids(xIndexLo, xIndexMid, yIndexMid + 1, yIndexHi);       // c5           N/2

            findVoids(xIndexMid + 1, xIndexHi, yIndexLo, yIndexMid);       // c6           N/2
            findVoids(xIndexMid + 1, xIndexHi, yIndexMid + 1, yIndexHi);   // c7           N/2

            processIndexedRegion(xIndexLo, xIndexHi, yIndexLo, yIndexHi, useCompleteSampling);  //              N/2
        }
    }
}
