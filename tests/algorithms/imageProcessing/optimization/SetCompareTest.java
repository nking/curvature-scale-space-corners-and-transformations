package algorithms.imageProcessing.optimization;

import algorithms.imageProcessing.SkylineExtractor;
import algorithms.util.PairInt;
import algorithms.misc.MiscMath;
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
        
        SecureRandom sr = new SecureRandom();
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nExpected = 1000;
        int maxDimension = 2000;
        
        Set<PairInt> expectedPoints = new HashSet<PairInt>();
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        while (expectedPoints.size() < nExpected) {
            
            int x = sr.nextInt(maxDimension);
            int y = sr.nextInt(maxDimension);
            PairInt p = new PairInt(x, y);
            while (expectedPoints.contains(p)) {
                x = sr.nextInt(maxDimension);
                y = sr.nextInt(maxDimension);
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
        
        //--- find the border points (consider them all found) ---
        //xMin, xMax, yMin, yMax
        int[] imgMinMaxXY = MiscMath.findMinMaxXY(points);
        int imgWidth = imgMinMaxXY[1] + 1;
        int imgHeight = imgMinMaxXY[3] + 1;                
        Set<PairInt> embedded = new HashSet<PairInt>();
        Set<PairInt> borderPoints = new HashSet<PairInt>();
        SkylineExtractor.getEmbeddedAndBorderPoints(points,
            imgWidth, imgHeight, embedded, borderPoints);        
                
        int nInside = points.size();
        
        int nOutside = (int)(sr.nextFloat()*(float)nExpected);
        
        for (int i = 0; i < nOutside; i++) {
            
            int x = sr.nextInt(maxDimension);
            int y = sr.nextInt(maxDimension);
            PairInt p = new PairInt(x, y);
            while (expectedPoints.contains(p)) {
                x = sr.nextInt(maxDimension);
                y = sr.nextInt(maxDimension);
                p = new PairInt(x, y);
            }  
            
            points.add(p);
        }
       
        SetCompare s = new SetCompare();
        SetComparisonResults scr = s.compare(expectedPoints, points, 
            borderPoints, borderPoints);
        
        assertTrue(scr.nExpectedPoints == nExpected);
        assertTrue(Math.abs(
            scr.numberMatchedDivExpected - (float)nInside/(float)nExpected)
            < 0.0000001);
        
        assertTrue(Math.abs(
            scr.numberOverrunDivExpectedMatchedPoints - (float)nOutside/(float)nExpected)
            < 0.0000001);
        
        assertTrue(scr.nExpectedBorderPoints == borderPoints.size());
        assertTrue(Math.abs(
            scr.numberMatchedBorderDivExpected - 1)< 0.0000001);
        
    }
}
