package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

public class RoughRangeSearchVoidFinder extends AbstractVoidFinder { 
    
    protected int xIdxLo = -1;
    
    protected int xIdxHi = -1;
    
    protected int yIdxLo = -1;
    
    protected int yIdxHi = -1;
    
    protected float binFactor = -1;
    
    protected int numberOfDivisions = -1;
    
    public void constructLogger() {
        this.log = Logger.getLogger(this.getClass().getName());
    }   
    
    public void setXIndexLow(int idx) {
        this.xIdxLo = idx;
    }
    
    public void setXIndexHigh(int idx) {
        this.xIdxHi = idx;
    }
    
    public void setYIndexLow(int idx) {
        this.yIdxLo = idx;
    }
    
    public void setYIndexHigh(int idx) {
        this.yIdxHi = idx;
    }
    
    public void setBinFactor(float bFactor) {
        this.binFactor = bFactor;
    }
    
    public void setNumberOfDivisions(int n) {
        this.numberOfDivisions = n;
    }
   
    @Override
    protected void findVoidsImpl() {
        
        if (xIdxLo == -1) {
            throw new IllegalStateException("xIdxLo is not set");
        }
        if (xIdxHi == -1) {
            throw new IllegalStateException("xIdxHi is not set");
        }
        if (yIdxLo == -1) {
            throw new IllegalStateException("yIdxLo is not set");
        }
        if (yIdxHi == -1) {
            throw new IllegalStateException("yIdxHi is not set");
        }
        if (numberOfDivisions == -1) {
            throw new IllegalStateException("numberOfDivisions is not set");
        }
        if (binFactor == -1) {
            throw new IllegalStateException("binFactor is not set");
        }
        
        findVoidsRoughRangeSearch(xIdxLo, xIdxHi, yIdxLo, yIdxHi, numberOfDivisions, binFactor);
    }

    /**
     * @param xIndexLo index given w.r.t. indexer.sortedXIndexes
     * @param xIndexHi index given w.r.t. indexer.sortedXIndexes
     * @param yIndexLo index given w.r.t. indexer.sortedYIndexes
     * @param yIndexHi index given w.r.t. indexer.sortedYIndexes
     * @param nDiv used to form the number of intervals = nPoints/nDiv
     *        which should usually be 2 or so
     * @param bFactor
     */
    protected void findVoidsRoughRangeSearch(int xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi, int nDiv, float bFactor) {

        boolean useCompleteSampling = (sampling.ordinal() == VoidSampling.COMPLETE.ordinal());
        
        //Divide y interval in half and execute the same size intervals in x over the full range
        int nYIntervals = (yIndexHi - yIndexLo) / nDiv;                     // cost     number of times

        for (int k = 0; k < nYIntervals; k++) {                             //               nDiv
            int binSz = (int)((k + 1) * bFactor);                           // c10

            int yLo = yIndexLo;
            while ((yLo + binSz) < yIndexHi) {                              //               nDiv

                int nXIntervals = (xIndexHi - xIndexLo)/ binSz;

                for (int j = 0; j < nXIntervals; j++) {                     //              nDiv/bfactor
                    int startX = xIndexLo + (j * binSz);
                    int endX = startX + binSz;
//System.out.println("processIndexedRegion: " + startX + ":" + endX + ":" + yLo + ":" + (yLo + binSz));
                    processIndexedRegion(startX, endX, yLo, yLo + binSz, useCompleteSampling);
                }
                yLo += binSz;
            }
        }
    }

}
