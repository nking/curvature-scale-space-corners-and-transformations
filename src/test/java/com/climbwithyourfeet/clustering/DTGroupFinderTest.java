package com.climbwithyourfeet.clustering;

import com.climbwithyourfeet.clustering.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class DTGroupFinderTest extends TestCase {
    
    public DTGroupFinderTest(String testName) {
        super(testName);
    }

    public void testCalculateGroups() throws Exception {
        
        Set<PairInt> points = getData0();
            
        DTGroupFinder<PairInt> finder = new DTGroupFinder<PairInt>();
        
        finder.setThreshholdFactor(1.0f);
        
        //float critSep = 2.f/(criticalDensity * threshholdFactor);
        
        float expectedCritSep = 2;
        
        float critDensity = 1.f/expectedCritSep;
        
        finder.calculateGroups(critDensity, points);
        
        assertTrue(finder.getNumberOfGroups() == 1);
        
        Set<PairInt> g0 = finder.getGroup(0);
        
        assertTrue(g0.size() == 6);
        
        assertTrue(points.contains(new PairInt(1, 2)));
        assertTrue(points.contains(new PairInt(2, 2)));
        assertTrue(points.contains(new PairInt(3, 2)));
        assertTrue(points.contains(new PairInt(2, 3)));
        assertTrue(points.contains(new PairInt(2, 5)));
        assertTrue(points.contains(new PairInt(3, 6)));
    }

    private Set<PairInt> getData0() {
        
        /*
         7                      @
         6          @
         5       @
         4          
         3       @
         2    @  @  @
         1
         0
           0  1  2  3  4  5  6  7
        
        */
        
        Set<PairInt> points = new HashSet<PairInt>();
        points.add(new PairInt(1, 2));
        points.add(new PairInt(2, 2));
        points.add(new PairInt(3, 2));
        points.add(new PairInt(2, 3));
        points.add(new PairInt(2, 5));
        points.add(new PairInt(3, 6));
        points.add(new PairInt(7, 7));
        
        return points;
    }
}
