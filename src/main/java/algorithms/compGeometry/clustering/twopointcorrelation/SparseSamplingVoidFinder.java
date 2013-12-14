package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.Arrays;
import java.util.logging.Logger;

import algorithms.misc.Statistic;

/**
 * samples only the non-zero cells in a dataset and those below
 * avg + standard deviation to determine 2-point void densities.
 * 
 * @author nichole
 *
 */
public class SparseSamplingVoidFinder extends AbstractVoidFinder {
        
    protected Statistic statistic = null;
    
    public void constructLogger() {
        this.log = Logger.getLogger(this.getClass().getName());
    }   
    
    @Override
    protected void findVoidsImpl() {

        findVoidsInSparseSubsets();
    }
   
    
    public void setStatistic(Statistic stat) {
        this.statistic = stat;
    }

    /**
     * find voids in statistic cells that are not empty and do not have
     * counts above avg + stdev.
     * 
     */
    protected void findVoidsInSparseSubsets() {
       
        if (statistic == null) {
            throw new IllegalStateException("statistic cannot be null");
        } 
        
        float[] xcopy = Arrays.copyOf(indexer.getX(), indexer.getNumberOfPoints());

        Arrays.sort(xcopy);
                
        
        int[] xLoSortedXIndexes = findSortedXIndexesForX(statistic.getItemsX(), xcopy);
        
        int[] yLoSortedXIndexes = findSortedXIndexesForY(statistic.getItemsY());
        
        float[] xHiF = new float[xLoSortedXIndexes.length];
        
        float[] yHiF = new float[xLoSortedXIndexes.length];
        

        float xSz = statistic.getXSz();
        
        float ySz = statistic.getYSz();
        
        
        for (int i = 0; i < statistic.getItemsX().length; i++) {
            
            float xLo = statistic.getItemsX()[i];
            
            float xHi = xLo + xSz;
            
            xHiF[i] = xHi;
            
            float yLo = statistic.getItemsY()[i];
            
            float yHi = yLo + ySz;
            
            yHiF[i] = yHi;
        }
        
        int[] xHiSortedXIndexes = findSortedXIndexesForX(xHiF, xcopy);
        
        int[] yHiSortedXIndexes = findSortedXIndexesForY(yHiF);
       
        
        float upperLimit = statistic.getAverage() + statistic.getStandardDeviation();
        
        for (int i = 0; i < statistic.getNumberOfItems(); i++) {
            
            int item = statistic.getItems()[i];
            
            if (item != 0) {
                
                if (item <= upperLimit) {
                   
                    int xLo = xLoSortedXIndexes[i];
                    
                    int yLo = yLoSortedXIndexes[i];
                    
                    int xHi = xHiSortedXIndexes[i];
                    
                    int yHi = yHiSortedXIndexes[i];
                    
                    findVoidsInSubset(xLo, xHi, yLo, yHi);
                }
            }
        }
    }

    protected void findVoidsInSubset(int xLo, int xHi, int yLo, int yHi) {
        
        if (statistic == null) {
            throw new IllegalStateException("statistic cannot be null");
        }
        
        for (int uSortedXIndex = xLo; uSortedXIndex < xHi; uSortedXIndex++) {
            
            for (int vSortedXIndex = yLo; vSortedXIndex < yHi; vSortedXIndex++) {
                
                if (uSortedXIndex == vSortedXIndex) {
                    continue;
                }
                
                processIndexedRegion(uSortedXIndex, vSortedXIndex);
            }
        }
    }
    
    /**
     * return the indexes of xValues w.r.t. indexer.sortexXIndexes
     * 
     * @param xValues
     * @param sortedX
     * @return
     */
    public int[] findSortedXIndexesForX(float[] xValues, float[] sortedX) {
        
        int[] xIndexes = new int[xValues.length];
    
        for (int i = 0; i < xValues.length; i++) {
            
            int idx = findSortedXIndexesForX(xValues[i], sortedX);
                        
            xIndexes[i] = idx;
        }
        
        return xIndexes;
    }
    
    /**
     * return the index w.r.t. indexer.sortedXIndexes for the xValue.
     * the index frame of reference is needed for a subsequent use.
     * 
     * @param xValue
     * @param sortedX
     * @return
     */
    protected int findSortedXIndexesForX(float xValue, float[] sortedX) {
        
        int idx = Arrays.binarySearch(sortedX, xValue);
        
        if (idx < 0) {
            idx = -1*(idx + 1);
        }
        
        if (idx > (indexer.getNumberOfPoints() - 1)) {
            idx = indexer.getNumberOfPoints() - 1;
        }
        
        // now we have the index w.r.t. indexer.x
        
        // we need to find idx within the values of indexer.sortedXIndexes
        
        for (int i = 0; i < indexer.getSortedXIndexes().length; i++) {
            if (idx == indexer.getSortedXIndexes()[i]) {
                return i;
            }
        }
        
        throw new IllegalStateException("could not find index " + idx + " in indexer.sortedXIndexes");
    }
    
    /**
     * return the indexes w.r.t. indexer.sortedXIndexes for the yValues.
     * the index frame of reference is needed for a subsequent use.
     * 
     * @param yValues
     * @param sortedX
     * @return
     */
    protected int[] findSortedXIndexesForY(float[] yValues) {
                
        int[] indexes = new int[yValues.length];
    
        for (int i = 0; i < yValues.length; i++) {
            
            float yV = yValues[i];
            
            int idx = indexer.findSortedXIndexesForY(yV);
        
            indexes[i] = idx;
        }
        
        return indexes;
    }
}
