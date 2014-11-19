package algorithms.imageProcessing;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PointMatcherTest {
    
    public PointMatcherTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testNumberOfMatches() {
        
        PointMatcher matcher = new PointMatcher();
        
        PairIntArray set1 = new PairIntArray();
        set1.add(10, 10);
        set1.add(20, 20);
        
        PairIntArray set2 = new PairIntArray();
        set2.add(47, 240);
        set2.add(100, 259);
        
        int transX = 125;
        int transY = 14;
        double transXTol = 10.3;
        double transYTol = 5.9;
        double rotation = 25*Math.PI/180.;
        double scale = 4;
        int centroidX1 = 100;
        int centroidY1 = 100;
        
        int nMatches = matcher.numberOfMatches(set1, set2, 
            transX, transY, transXTol, transYTol, rotation, scale, 
            centroidX1, centroidY1);
        
        assertTrue(nMatches == 2);
    }
    
    @Test
    public void testCalculateTranslation() {
        
        PointMatcher matcher = new PointMatcher();
        
        PairIntArray set1 = new PairIntArray();
        set1.add(10, 10);
        set1.add(20, 20);
        
        PairIntArray set2 = new PairIntArray();
        set2.add(47, 240);
        set2.add(100, 259);
        
        int transX = 125;
        int transY = 14;
        double transXTol = 10.3;
        double transYTol = 5.9;
        double rotation = 25*Math.PI/180.;
        double scale = 4;
        int centroidX1 = 100;
        int centroidY1 = 100;
        
        int[] transXY = matcher.calculateTranslation(set1, set2,
            transXTol, transYTol, rotation, scale, 
            centroidX1, centroidY1);
        
        assertTrue(Math.abs(transXY[0] - transX) < transXTol);
        
        assertTrue(Math.abs(transXY[1] - transY) < transYTol);
    }
    
}
