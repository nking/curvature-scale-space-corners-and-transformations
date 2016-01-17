package algorithms.imageProcessing.features;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.List;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class BlobCornersEuclideanCalculator2Test extends TestCase {
    
    public BlobCornersEuclideanCalculator2Test() {
    }

    public void testCountFrequencyThenSortDesc() {
        
        List<List<FeatureComparisonStat>> matchedLists = new 
            ArrayList<List<FeatureComparisonStat>>();
        
        int c3 = 1;
        // making test data w/ different keys, and same number of values for subsequent odd i
        for (int i = 0; i < 10; ++i) {
            
            List<FeatureComparisonStat> list = new ArrayList<FeatureComparisonStat>();
            
            int nt = c3;
            c3++;
            if (c3 == 4) {
                c3 = 1;
            }
            
            int count = 0;
            while (count < nt) {
                list.add(new FeatureComparisonStat());
                count++;
            }
            
            matchedLists.add(list);
        }
        
        BlobCornersEuclideanCalculator2 instance = new BlobCornersEuclideanCalculator2();
        
        // x=size of lists of stats, y=number of blobs with that size
        List<PairInt> sortedFreq = new ArrayList<PairInt>();
        
        // x=size of lists of stats, y=number of blobs with that size
        List<PairInt> sortedSizes = new ArrayList<PairInt>();
        
        instance.countFrequencyThenSortDesc(matchedLists, sortedFreq, sortedSizes);
        
        assertTrue(sortedFreq.size() == 3);
        
        assertEquals(1, sortedFreq.get(0).getX());
        assertEquals(4, sortedFreq.get(0).getY());
        
        assertEquals(3, sortedFreq.get(1).getY());
        assertEquals(3, sortedFreq.get(2).getY());
        
        assertTrue((sortedFreq.get(1).getX() == 2) || (sortedFreq.get(1).getX() == 3));
        assertTrue((sortedFreq.get(2).getX() == 2) || (sortedFreq.get(2).getX() == 3));
        assertFalse(sortedFreq.get(1).getX() == sortedFreq.get(2).getX());
        /*
        2 - 3
        3 - 3
        1 - 4
        */
        assertTrue(sortedSizes.size() == 3);
        
        assertEquals(3, sortedSizes.get(0).getX());
        assertEquals(3, sortedSizes.get(0).getY());
        
        assertEquals(2, sortedSizes.get(1).getX());
        assertEquals(3, sortedSizes.get(1).getY());
        
        assertEquals(1, sortedSizes.get(2).getX());
        assertEquals(4, sortedSizes.get(2).getY());
    }
    
}
