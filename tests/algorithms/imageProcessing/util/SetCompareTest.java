package algorithms.imageProcessing.util;

import algorithms.util.PairInt;
import java.security.SecureRandom;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SetCompareTest extends TestCase {
    
    public SetCompareTest() {
    }
    
    public void testCompare() throws Exception {
        /*
        public SetComparisonResults compare(Set<PairInt> expectedPoints,
        Set<PairInt> points) {
        */
        
        SecureRandom sr = new SecureRandom();
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nExpected = 1000;
        
        Set<PairInt> expectedPoints = new HashSet<PairInt>();
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        while (expectedPoints.size() < nExpected) {
            
            int x = sr.nextInt();
            int y = sr.nextInt();
            PairInt p = new PairInt(x, y);
            while (expectedPoints.contains(p)) {
                x = sr.nextInt();
                y = sr.nextInt();
                p = new PairInt(x, y);
            }
            
            expectedPoints.add(p);
            
            boolean include = true;
            if ((expectedPoints.size() % 10) == 0) {
                include = sr.nextBoolean();
            }
            
            if (include) {
                points.add(p);
            }
        }
        
        int nInside = points.size();
        
        int nOutside = (int)(sr.nextFloat()*(float)nExpected);
        
        for (int i = 0; i < nOutside; i++) {
            
            int x = sr.nextInt();
            int y = sr.nextInt();
            PairInt p = new PairInt(x, y);
            while (expectedPoints.contains(p)) {
                x = sr.nextInt();
                y = sr.nextInt();
                 p = new PairInt(x, y);
            }  
            
            points.add(p);
        }
        
        SetCompare s = new SetCompare();
        SetComparisonResults scr = s.compare(expectedPoints, points);
        
        assertTrue(scr.nExpectedPoints == nExpected);
        assertTrue(Math.abs(
            scr.numberMatchedDivExpected - (float)nInside/(float)nExpected)
            < 0.0000001);
        
        assertTrue(Math.abs(
            scr.numberOverrunDivExpected - (float)nOutside/(float)nExpected)
            < 0.0000001);
        
    }
}
