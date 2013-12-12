package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

public class CompleteSamplingVoidFinder extends AbstractVoidFinder {

    protected int skipInterval = 1;
    
    public void constructLogger() {
        this.log = Logger.getLogger(this.getClass().getName());
    }   
    
    public void setSkipInterval(int skip) {
        this.skipInterval = skip;
    }
    
    @Override
    protected void findVoidsImpl() {

        findVoidsUsingDoubleIndexes(skipInterval);
    }

    protected void findVoidsUsingDoubleIndexes(int incr) {
        // N!/(2!(N-2)! * N!/(2!(N-2)!

        findVoidsUsingDoubleIndexes(0, indexer.getNXY() - 1, 0, indexer.getNXY() - 1, incr);
    }

    protected void findVoidsUsingDoubleIndexes(int xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi, int incr) {
        // N!/(2!(N-2)! * N!/(2!(N-2)!

        boolean useCompleteSampling = (sampling.ordinal() == VoidSampling.COMPLETE.ordinal());
        
        for (int i = xIndexLo; i < xIndexHi; i++) {
            if (debug) {
                log.info("findVoids i=" + i + "/" + indexer.getNXY());
            }
            for (int ii = (i + 1); ii < indexer.getNXY(); ii+=incr) {
                for (int j = yIndexLo; j < yIndexHi; j++) {
                    for (int jj = (j + 1); jj < indexer.getNXY(); jj+=incr) {
                        processIndexedRegion(i, ii, j, jj, useCompleteSampling);
                    }
                }
            }
        }
    }

}
