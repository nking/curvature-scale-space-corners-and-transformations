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
        
        int maxDimension = 2000;
        
        Set<PairInt> expectedPoints = new HashSet<PairInt>();
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        int xMin = sr.nextInt(maxDimension/4) + 5;
        int xMax = xMin + sr.nextInt(3*maxDimension/4);
        
        int yMin = sr.nextInt(maxDimension/4) + 5;
        int yMax = yMin + sr.nextInt(3*maxDimension/4);
        
        long nExpected = (yMax - yMin + 1) * (xMax - xMin + 1);
        
        // make a solid block of points, and a set with some randomly missing.
        // do the same for border points
        long count = 0;
        for (int col = xMin; col <= xMax; ++col) {
            for (int row = yMin; row <= yMax; ++row) {
                
                PairInt p = new PairInt(col, row);
                
                expectedPoints.add(p);
                
                boolean add = true;
                // want roughly 10 percent not present:
                if ((count % (nExpected/10)) == 0) {
                    if (!sr.nextBoolean()) {
                        add = false;
                    }
                }
                if (add) {
                    points.add(p);
                }
                
                count++;
            }
        }
        
        Set<PairInt> expectedBorderPoints = new HashSet<PairInt>();
        for (int col = xMin; col <= xMax; ++col) {
            expectedBorderPoints.add(new PairInt(col, yMin));
            expectedBorderPoints.add(new PairInt(col, yMax));
        }
        for (int row = yMin; row <= yMax; ++row) {
            expectedBorderPoints.add(new PairInt(xMin, row));
            expectedBorderPoints.add(new PairInt(xMax, row));
        }
        
        int imgWidth = maxDimension;
        int imgHeight = maxDimension;
        
        Set<PairInt> embedded = new HashSet<PairInt>();
        Set<PairInt> borderPoints = new HashSet<PairInt>();
        SkylineExtractor.getEmbeddedAndBorderPoints(points,
            imgWidth, imgHeight, embedded, borderPoints);
        
        Set<PairInt> embedded2 = new HashSet<PairInt>();
        Set<PairInt> borderPoints2 = new HashSet<PairInt>();
        SkylineExtractor.getEmbeddedAndBorderPoints(expectedPoints,
            imgWidth, imgHeight, embedded2, borderPoints2);
        
        assertTrue(expectedBorderPoints.size() == borderPoints2.size());
        
        SetCompare s = new SetCompare();
        SetComparisonResults scr = s.compare(expectedPoints, points, 
            expectedBorderPoints, borderPoints);
        
        assertTrue(scr.nExpectedPoints == nExpected);
        assertTrue(Math.abs(
            scr.numberMatchedDivExpected - (float)points.size()/(float)nExpected)
            < 0.0000001);
        
        assertTrue(Math.abs(
            scr.numberOverrunDivExpectedMatchedPoints - 0)
            < 0.0000001);
        
        //System.out.println("scr.nExpectedBorderPoints=" + scr.nExpectedBorderPoints 
        //    + " expectedBorderPoints=" + expectedBorderPoints.size());
        
        //System.out.println("scr.numberMatchedBorderDivExpected=" + scr.numberMatchedBorderDivExpected 
        //    + " expected=" + ((float)borderPoints.size()/(float)expectedBorderPoints.size()));
            
        assertTrue(scr.nExpectedBorderPoints == expectedBorderPoints.size());
        assertTrue(Math.abs(
            scr.numberMatchedBorderDivExpected - 
                (float)borderPoints.size()/(float)expectedBorderPoints.size()) 
            < 0.03);
        
    }
}
