package algorithms.misc;

import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Random;
import java.util.logging.Logger;

import algorithms.compGeometry.clustering.twopointcorrelation.DoubleAxisIndexer;

public class DoubleAxisIndexerStats {

    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected float calculateBinSize(int numberOfBins, float minValue, float maxValue) {

        float divSz = (1.00001f*maxValue - minValue)/numberOfBins;
        
        return divSz;
    }
    
    /**
     * count the number of points within each cell and do basic statistics on those counts.
     * 
     * @param numberOfCellsInOneDimension number of cells across one dimension to divide the
     *  data into to count the number of points within each cell
     * @param indexer
     * @return statistic the basic statistics of the point counts within the cells
     */
    public Statistic calculateCellDensities(int numberOfCellsInOneDimension, DoubleAxisIndexer indexer) {

        //xmin, xmax, ymin, ymax 
        float[] xyMinMax = indexer.findXYMinMax();

        float xDivSz = calculateBinSize(numberOfCellsInOneDimension, xyMinMax[0], xyMinMax[1]);
                
        float yDivSz = calculateBinSize(numberOfCellsInOneDimension, xyMinMax[2], xyMinMax[3]);
        
        int[] items = new int[numberOfCellsInOneDimension*numberOfCellsInOneDimension];
        
        for (int i = 0; i < indexer.getX().length; i++) {
            
            float x = indexer.getX()[i];
            
            float y = indexer.getY()[i];
            
            int xCell = (int) ((x - xyMinMax[0])/xDivSz);
            
            int yCell = (int) ((y - xyMinMax[2])/yDivSz);
            
            int itemN = yCell*numberOfCellsInOneDimension + xCell;
            
            items[itemN] =  items[itemN] + 1;   
        }

        Statistic statistic = new Statistic(items);
        
        return statistic;
    }
    
    /**
     * looking at cells across the dataset to see if they have roughly equal numbers of points in them.
     * 
     * The division is by the x and y range of values.
     * 
     * @param numberOfCellsInOneDimension number of cells across one dimension to divide the
     *  data into to count the number of points within each cell
     * @param indexer
     * @return whether all of the counts within cells are similar within one standard deviation
     */
    public boolean allAreSame(int numberOfCellsInOneDimension, DoubleAxisIndexer indexer) {
        
        Statistic statistic = calculateCellDensities(numberOfCellsInOneDimension, indexer);
        
        if (statistic == null) {
            throw new IllegalArgumentException("could not calculate statistics from given indexer");
        }
                
        float stDev = statistic.getStandardDeviation();
        
        float avg = statistic.getAverage();
        
        for (int item : statistic.getItems()) {
            
            float diff = Math.abs(item - avg);
            
            if (diff > (1.0f * stDev)) {
            
                return false;
            }
        }
        
        return true;
    }
    
    /**
     * determine whether counts in each cell suggest that there are gaps in the data, that
     * is whether some cells have significantly fewer counts than others.
     * 
     * @param numberOfCellsInOneDimension
     * @param indexer
     * @return whether there are not cells with significantly fewer points in them than average
     */
    public boolean doesNotHaveLargeGaps(int numberOfCellsInOneDimension, DoubleAxisIndexer indexer) {
        
        Statistic statistic = calculateCellDensities(numberOfCellsInOneDimension, indexer);
        
        if (statistic == null) {
            throw new IllegalArgumentException("could not calculate statistics from given indexer");
        }
                
        float stDev = statistic.getStandardDeviation();
        
        float avg = statistic.getAverage();
        
        for (int item : statistic.getItems()) {
        
            if (item < (avg - stDev)) {
            
                return false;
            }
        }
        
        return true;
    }

    /**
     * choose a random cell within the data were a cell is a division of the data by
     * numberOfCellsInOneDimension for each dimension.  it returns the cell boundaries
     * in x and then y as indexes that are indexes relative to indexer.getX() and indexer.get().
     * 
     * @param numberOfCellsInOneDimension
     * @param indexer
     * @return an array holding the index boundaries of the cell 
     * as new int[]{xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi}
     */
    public int[] chooseARandomCell(int numberOfCellsInOneDimension, DoubleAxisIndexer indexer) {
        
        Random sr = null;
        
        try {
         
            sr = SecureRandom.getInstance("SHA1PRNG");
            
        } catch (NoSuchAlgorithmException e) {
            
            log.severe(e.getMessage());
            
            sr = new Random();
        }
        
        long seed = System.currentTimeMillis();
        seed = 1386620575944l;
        sr.setSeed(seed);
        
        int nXY = indexer.getNXY();
                                
        int xDivIndexesSz = nXY/numberOfCellsInOneDimension;

        int yDivIndexesSz = xDivIndexesSz;
        
        int anXIndex = (int)(sr.nextFloat()*(nXY - xDivIndexesSz));
        
        int anYIndex = (int)(sr.nextFloat()*(nXY - yDivIndexesSz));
                
        int xCell = (int) (anXIndex/xDivIndexesSz);

        int yCell = (int) (anYIndex/yDivIndexesSz);
        
        int xIndexLo = xCell*xDivIndexesSz;
  
        int xIndexHi = xIndexLo + xDivIndexesSz;
        
        int yIndexLo = yCell*yDivIndexesSz;
        
        int yIndexHi = yIndexLo + yDivIndexesSz;
        
        return new int[]{xIndexLo, xIndexHi, yIndexLo, yIndexHi};
    }
}
