package algorithms.misc;

import java.util.logging.Logger;

import algorithms.compGeometry.clustering.twopointcorrelation.DoubleAxisIndexer;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 * @author nichole
 */
public class DoubleAxisIndexerStatsTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void testCalculateCellDensities() throws Exception {
        
        DoubleAxisIndexerStats stats = new DoubleAxisIndexerStats();
        
        // uniform grid of data:
        
        float[] x = new float[] {
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8
        };
        float[] y = new float[] {
            0,  0,  0,  0,  0,  0,  0,  0,  0,
            1,  1,  1,  1,  1,  1,  1,  1,  1,
            2,  2,  2,  2,  2,  2,  2,  2,  2,
            3,  3,  3,  3,  3,  3,  3,  3,  3,
            4,  4,  4,  4,  4,  4,  4,  4,  4,
            5,  5,  5,  5,  5,  5,  5,  5,  5,
            6,  6,  6,  6,  6,  6,  6,  6,  6,
            7,  7,  7,  7,  7,  7,  7,  7,  7,
            8,  8,  8,  8,  8,  8,  8,  8,  8,
            9,  9,  9,  9,  9,  9,  9,  9,  9,
           10, 10, 10, 10, 10, 10, 10, 10, 10
        };
        
        int numberOfCellsInOneDimension = 2;
        
        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(x, y, x.length);
        
        Statistic stat = stats.calculateCellDensities(numberOfCellsInOneDimension, indexer);
        
        assertNotNull(stat);
        
        float expectedAvg = (float)(indexer.getNXY()/Math.pow(numberOfCellsInOneDimension, 2));
        
        assertTrue(Math.abs(expectedAvg - stat.getAverage()) < 0.5); 
                
    }

    public void testAllAreSame() throws Exception {
        
        DoubleAxisIndexerStats stats = new DoubleAxisIndexerStats();
        
        // uniform grid of data:
        float[] x = new float[] {
            0,  1,  2,  3, 
            0,  1,  2,  3,
            0,  1,  2,  3,  
            0,  1,  2,  3
        };
        float[] y = new float[] {
            0,  0,  0,  0,
            1,  1,  1,  1, 
            2,  2,  2,  2,
            3,  3,  3,  3
        };
        
        int numberOfCellsInOneDimension = 2;
        
        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(x, y, x.length);
     
        float fractionOutliers = stats.fractionOfCellsOutSideOfAvgTolerance
            (numberOfCellsInOneDimension, indexer, 2.5f);
        
        assertTrue(fractionOutliers < 0.1f);
    }

    public void testAllAreSame1() throws Exception {
        
        DoubleAxisIndexerStats stats = new DoubleAxisIndexerStats();
        
        // non-uniform grid of data:     
        float[] x = new float[] {
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9
        };
        float[] y = new float[] {
            0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
            1,  1,  1,  1,  1,  1,  1,  1, 1, 1,
            2,  2,  2,  2,  2,  2,  2,  2, 2, 2,
            3,  3,  3,  3,  3,  3,  3,  3, 3, 3,
            4,  4,  4,  2,  2,  2,  2,  2, 2, 2,
            5,  5,  5,  2,  2,  2,  2,  3, 3, 3,
            6,  6,  6,  6,  6,  6,  6,  6, 6, 6,
            7,  7,  7,  7,  7,  7,  7,  7, 7, 7,
            8,  8,  8,  8,  8,  8,  8,  8, 8, 8,
            9,  9,  9,  9,  9,  9,  9,  9, 9, 9
        };
        
        int numberOfCellsInOneDimension = 3;
        
        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(x, y, x.length);
     
        float fractionOutliers = stats.fractionOfCellsOutSideOfAvgTolerance
            (numberOfCellsInOneDimension, indexer, 2.5f);
        
        assertTrue(fractionOutliers > 0.1f);
    }
    
    public void testDoesNotHaveLargeGaps() throws Exception {
        
        DoubleAxisIndexerStats stats = new DoubleAxisIndexerStats();
                
        // non-uniform grid of data:     
        float[] x = new float[] {
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9,
            0,  1,  2,  3,  4,  5,  6, 7, 8, 9
        };
        float[] y = new float[] {
            0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
            1,  1,  1,  1,  1,  1,  1,  1, 1, 1,
            2,  2,  2,  2,  2,  2,  2,  2, 2, 2,
            3,  3,  3,  3,  3,  3,  3,  3, 3, 3,
            4,  4,  4,  2,  2,  2,  2,  2, 2, 2,
            5,  5,  5,  2,  2,  2,  2,  3, 3, 3,
            6,  6,  6,  6,  6,  6,  6,  6, 6, 6,
            7,  7,  7,  7,  7,  7,  7,  7, 7, 7,
            8,  8,  8,  8,  8,  8,  8,  8, 8, 8,
            9,  9,  9,  9,  9,  9,  9,  9, 9, 9
        };
        
        int numberOfCellsInOneDimension = 3;
        
        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(x, y, x.length);
     
        boolean same = stats.doesNotHaveLargeGaps(numberOfCellsInOneDimension, indexer, 2.5f);
        
        assertFalse(same);
    }
    
    public void testDoesNotHaveLargeGaps2() throws Exception {
        
        DoubleAxisIndexerStats stats = new DoubleAxisIndexerStats();
                
        float[] x = new float[] {
            0,  1,  2,  3, 
            0,  1,  2,  3,
            0,  1,  2,  3,  
            0,  1,  2,  3
        };
        float[] y = new float[] {
            0,  0,  0,  0,
            1,  1,  1,  1, 
            2,  2,  2,  2,
            3,  3,  3,  3
        };
        
        int numberOfCellsInOneDimension = 3;
        
        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(x, y, x.length);
     
        float fractionOutliers = stats.fractionOfCellsOutSideOfAvgTolerance
            (numberOfCellsInOneDimension, indexer, 2.5f);
        
        assertTrue(fractionOutliers < 0.1f);
    }
    
    public void testChooseARandomCell() throws Exception {
        
        DoubleAxisIndexerStats stats = new DoubleAxisIndexerStats();
        
        float[] x = new float[] {
            0,  1,  2,  3, 
            0,  1,  2,  3,
            0,  1,  2,  3,  
            0,  1,  2,  3
        };
        float[] y = new float[] {
            0,  0,  0,  0,
            1,  1,  1,  1, 
            2,  2,  2,  2,
            3,  3,  3,  3
        };
        
        int numberOfCellsInOneDimension = 2;
        
        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(x, y, x.length);
        
        int[] indexRanges = stats.chooseARandomCell(numberOfCellsInOneDimension, indexer);
        
        assertNotNull(indexRanges);
        
        assertTrue(indexRanges[0] < indexRanges[1]);
                        
        assertTrue(indexRanges[2] < indexRanges[3]);
                
    }
}
