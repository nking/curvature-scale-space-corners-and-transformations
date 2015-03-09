package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

/**
 * an implementation of IVoidFinder that uses a divide and conquer approach to sampling the
 * data.
 * 
 * @author nichole
 *
 */
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
     * <pre>
     * A divide and conquer approach to finding the rectangular areas
     * containing only two points.  
     *
     * runtime complexity:
     *     processIndexedRegion is a little worse than linear
     *     T(processIndexedRegion) = O(N)
     *
     *     then O(findVoids) = (c000 + c001) + (c005)*( O(N) ) + c003*T(findVoids) + c004*T(findVoids)
     *                       = 2T(N/2) + N
     *                       ==> this has a solution that is O(N lg2 N)
     * 
     * Note: It does not completely sample every pair of points.
     *</pre>
     * @param xIndexLo
     * @param xIndexHi
     */
    protected void findVoids(int xIndexLo, int xIndexHi) {
                
        if (xIndexLo < xIndexHi) {    

            int xIndexMid = (xIndexLo + xIndexHi) >> 1;         

            findVoids(xIndexLo, xIndexMid);      

            findVoids(xIndexMid + 1, xIndexHi);     

            processIndexedRegion(xIndexLo, xIndexHi); 
        }
    }
}
