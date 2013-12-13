package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

public class DivideAndConquerVoidFinder extends AbstractVoidFinder {

    @Override
    public void constructLogger() {
        this.log = Logger.getLogger(this.getClass().getName());
    }

    @Override
    protected void findVoidsImpl() {
        
        findVoids(0, indexer.nXY - 1);
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
    protected void findVoids(int xIndexLo, int xIndexHi) {
                
        if (xIndexLo < xIndexHi) {    

            int xIndexMid = (xIndexLo + xIndexHi)/2;         

            findVoids(xIndexLo, xIndexMid);      

            findVoids(xIndexMid + 1, xIndexHi);     

            processIndexedRegion(xIndexLo, xIndexHi); 
        }
    }
}
