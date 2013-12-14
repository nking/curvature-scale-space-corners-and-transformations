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
        
        if (numberOfCellsInOneDimension == 0) {
            throw new IllegalArgumentException("numberOfCellsInOneDimension must be larger than 0");
        }
        if (indexer == null) {
            throw new IllegalArgumentException("indexer cannot be null");
        }
        
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
            
            int itemN = yCell*numberOfCellsInOneDimension + xCell; // (0,0) (1,0)
            
            items[itemN] =  items[itemN] + 1;   
        }
        
        // calculate the start x and y values of the cells in items[].
        
        float[] xCells = new float[items.length];
        
        float[] yCells = new float[items.length];
        
        int count = 0;
        
        for (int i = 0; i < numberOfCellsInOneDimension; i++) {
            
            float yStart = xyMinMax[2] + i*yDivSz;
            
            if (yStart > xyMinMax[3]) {
                yStart = xyMinMax[3];
            }
            
            for (int j = 0; j < numberOfCellsInOneDimension; j++) {
                
                float xStart = xyMinMax[0] + j*xDivSz;
                
                if (xStart > xyMinMax[1]) {
                    xStart = xyMinMax[1];
                }
                
                xCells[count] = xStart;
                
                yCells[count] = yStart;
                
                count++;
            }
        }

        Statistic statistic = new Statistic(items, xCells, yCells, xDivSz, yDivSz);
        
        return statistic;
    }
    
    /**
     * return the fraction of cells whose counts are outside of average +- standardDeviation*factorStDev
     * 
     * @param numberOfCellsInOneDimension number of cells across one dimension to divide the
     *  data into to count the number of points within each cell
     * @param indexer
     * @param factorStDev the factor to test whether a cell's counts are within average += stDev*factorStDev.
     * factorStDev is usually a value of 2 or 3.  The value should be equal to TwoPointCorrelation.sigmaFactor.
     * @return whether all of the counts within cells are similar within one standard deviation
     */
    public float fractionOfCellsOutSideOfAvgTolerance(int numberOfCellsInOneDimension, DoubleAxisIndexer indexer, float factorStDev) {
        
        Statistic statistic = calculateCellDensities(numberOfCellsInOneDimension, indexer);
        
        if (statistic == null) {
            throw new IllegalArgumentException("could not calculate statistics from given indexer");
        }
                
        return fractionOfCellsOutSideOfAvgTolerance(statistic, factorStDev);
    }
    
    /**
     * return the fraction of cells whose counts are outside of average +- standardDeviation*factorStDev
     * 
     * @param statistic
     * @param factorStDev the factor to test whether a cell's counts are within average += stDev*factorStDev.
     * factorStDev is usually a value of 2 or 3.  The value should be equal to TwoPointCorrelation.sigmaFactor.
     * @return whether all of the counts within cells are similar within one standard deviation
     */
    public float fractionOfCellsOutSideOfAvgTolerance(Statistic statistic, float factorStDev) {
                
        if (statistic == null) {
            throw new IllegalArgumentException("statistic cannot be null");
        }
        
        int nAboveAvgPlusStDev = 0;
                
        float stDev = statistic.getStandardDeviation();
        
        float avg = statistic.getAverage();
        
        for (int item : statistic.getItems()) {
            
            float diff = Math.abs(item - avg);
            
            if (diff > (factorStDev * stDev)) {
            
                nAboveAvgPlusStDev++;
            }
        }
        
        return (float)nAboveAvgPlusStDev/(float)statistic.getNumberOfItems();
    }
    
    /**
     * determine whether counts in each cell suggest that there are gaps in the data, that
     * is whether some cells have significantly fewer counts than others.
     * 
     * @param numberOfCellsInOneDimension
     * @param indexer
     * @param factorStDev the factor to test whether a cell's counts are within average += stDev*factorStDev.
     * factorStDev is usually a value of 2 or 3.  The value should be equal to TwoPointCorrelation.sigmaFactor.
     * @return whether there are not cells with significantly fewer points in them than average
     */
    public boolean doesNotHaveLargeGaps(int numberOfCellsInOneDimension, DoubleAxisIndexer indexer, float factorStDev) {
        
        Statistic statistic = calculateCellDensities(numberOfCellsInOneDimension, indexer);
        
        if (statistic == null) {
            throw new IllegalArgumentException("could not calculate statistics from given indexer");
        }
                
        return doesNotHaveLargeGaps(statistic, factorStDev);
    }

    /**
     * determine whether counts in each cell suggest that there are gaps in the data, that
     * is whether some cells have significantly fewer counts than others.
     * 
     * @param statistic
     * @param factorStDev the factor to test whether a cell's counts are within average += stDev*factorStDev.
     * factorStDev is usually a value of 2 or 3.  The value should be equal to TwoPointCorrelation.sigmaFactor.
     * @return whether there are not cells with significantly fewer points in them than average
     */
    public boolean doesNotHaveLargeGaps(Statistic statistic, float factorStDev) {
                
        if (statistic == null) {
            throw new IllegalArgumentException("statistic cannot be null");
        }
                
        float stDev = statistic.getStandardDeviation();
        
        float avg = statistic.getAverage();
        
        for (int item : statistic.getItems()) {
        
            if ((item == 0) || (item < (avg - stDev*factorStDev)) ) {
            
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
    public float fractionOfCellsWithoutPoints(int numberOfCellsInOneDimension, DoubleAxisIndexer indexer) {
        
        Statistic statistic = calculateCellDensities(numberOfCellsInOneDimension, indexer);
        
        return fractionOfCellsWithoutPoints(statistic);
    }
    
    /**
     * determine whether counts in each cell suggest that there are gaps in the data, that
     * is whether some cells have significantly fewer counts than others.
     * 
     * @param statistic
     * @return whether there are not cells with significantly fewer points in them than average
     */
    public float fractionOfCellsWithoutPoints(Statistic statistic) {
        
        if (statistic == null) {
            throw new IllegalArgumentException("could not calculate statistics from given indexer");
        }
                
        int count = 0;
        
        for (int item : statistic.getItems()) {
            if (item == 0) {
                count++;
            }
        }
        
        return (float)count/statistic.getNumberOfItems();
    }

    /**
     * choose a random cell within the data were a cell is a division of the data by
     * numberOfCellsInOneDimension for each dimension.  it returns the cell boundaries
     * in x and then y as indexes that are indexes relative to 
     * {indexer.getSortedXIndexes(), indexer.getSortedXIndexes(), indexer.getSortedYIndexes(),
     * indexer.getSortedYIndexes()}
     * 
     * @param numberOfCellsInOneDimension
     * @param indexer
     * @return an array holding the index boundaries of the cell 
     * as new int[]{xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi} where the indexes should be treated
     * with respect to arrays {indexer.getSortedXIndexes(), indexer.getSortedXIndexes(), indexer.getSortedYIndexes(),
     * indexer.getSortedYIndexes()}
     */
    public int[] chooseARandomCell(int numberOfCellsInOneDimension, DoubleAxisIndexer indexer) {
        
        if (numberOfCellsInOneDimension < 1) {
            throw new IllegalArgumentException("numberOfCellsInOneDimension must be larger than 0");
        }
        
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
