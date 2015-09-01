package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.List;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class IntensityFeatureComparisonStatsTest extends TestCase {
    
    public IntensityFeatureComparisonStatsTest() {
    }
    
    public void testCompareTo() throws Exception {
        
        FixedSizeSortedVector<IntensityFeatureComparisonStats> vec = new 
            FixedSizeSortedVector<>(2, IntensityFeatureComparisonStats.class);
        
        IntensityFeatureComparisonStats cStats0 = new IntensityFeatureComparisonStats();
        List<FeatureComparisonStat> statList0 = new ArrayList<FeatureComparisonStat>();
        FeatureComparisonStat c0 = new FeatureComparisonStat();
        c0.setSumIntensitySqDiff(10);
        c0.setImg2PointIntensityErr(1);
        statList0.add(c0);
        cStats0.addAll(statList0);
        
        IntensityFeatureComparisonStats cStats1 = new IntensityFeatureComparisonStats();
        List<FeatureComparisonStat> statList1 = new ArrayList<FeatureComparisonStat>();
        FeatureComparisonStat c1 = new FeatureComparisonStat();
        c1.setSumIntensitySqDiff(100);
        c1.setImg2PointIntensityErr(10);
        statList1.add(c1);
        cStats1.addAll(statList1);
        
        IntensityFeatureComparisonStats cStats2 = new IntensityFeatureComparisonStats();
        List<FeatureComparisonStat> statList2 = new ArrayList<FeatureComparisonStat>();
        FeatureComparisonStat c2 = new FeatureComparisonStat();
        c2.setSumIntensitySqDiff(50);
        c2.setImg2PointIntensityErr(10);
        statList2.add(c2);
        cStats2.addAll(statList2);
        
        vec.add(cStats0);
        vec.add(cStats1);
        vec.add(cStats2);
        
        assertTrue(vec.getNumberOfItems() == 2);
        
        IntensityFeatureComparisonStats[] sorted = vec.getArray();
        
        IntensityFeatureComparisonStats sorted0 = sorted[0];
        assertTrue(sorted0.equals(cStats0));
        IntensityFeatureComparisonStats sorted1 = sorted[1];
        assertTrue(sorted1.equals(cStats2));
    }
}
