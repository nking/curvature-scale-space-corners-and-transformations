package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class EdgeExtractorForBlobBorderTest extends TestCase {
    
    public EdgeExtractorForBlobBorderTest() {
    }

    public void testExtractAndOrderTheBorder() {
        
        int width = 10;
        int height = 10;
        int xCenter = 50; 
        int yCenter = 50;
        
        Set<PairInt> rectangle10 = DataForTests.getRectangle(width, height, 
            xCenter, yCenter);
                
        int imageWidth = 200;
        int imageHeight = 200;
        boolean discardWhenCavityIsSmallerThanBorder = true;
        
        Set<PairInt> expected = new HashSet<PairInt>();
        for (int x = (xCenter - (width/2)); x < (xCenter + (width/2)); ++x) {
            expected.add(new PairInt(x, 45));
            expected.add(new PairInt(x, 54));
        }
        for (int y = (yCenter - (height/2)); y < (yCenter + (height/2)); ++y) {
            expected.add(new PairInt(45, y));
            expected.add(new PairInt(54, y));
        }
        assertTrue(expected.size() == 36);
        
        EdgeExtractorForBlobBorder instance = new EdgeExtractorForBlobBorder();
                
        PairIntArray result = instance.extractAndOrderTheBorder(
            rectangle10, imageWidth, imageHeight, 
            discardWhenCavityIsSmallerThanBorder);
        
        assertTrue(result.getN() == expected.size());
        
        for (int i = 0; i < result.getN(); ++i) {
            PairInt p = new PairInt(result.getX(i), result.getY(i));
            assertTrue(expected.remove(p));
        }
        
        assertTrue(expected.isEmpty());
    }

    public void testExtractAndOrderTheBorder0() {
        
        int width = 10;
        int height = 10;
        int xCenter = 50; 
        int yCenter = 50;
        
        Set<PairInt> rectangle10 = DataForTests.getRectangle(width, height, 
            xCenter, yCenter);
                
        int imageWidth = 200;
        int imageHeight = 200;
        boolean discardWhenCavityIsSmallerThanBorder = true;
        
        Set<PairInt> expected = new HashSet<PairInt>();
        for (int x = (xCenter - (width/2)); x < (xCenter + (width/2)); ++x) {
            expected.add(new PairInt(x, 45));
            expected.add(new PairInt(x, 54));
        }
        for (int y = (yCenter - (height/2)); y < (yCenter + (height/2)); ++y) {
            expected.add(new PairInt(45, y));
            expected.add(new PairInt(54, y));
        }
        assertTrue(expected.size() == 36);
        
        EdgeExtractorForBlobBorder instance = new EdgeExtractorForBlobBorder();
                
        PairIntArray result = instance.extractAndOrderTheBorder0(
            rectangle10, imageWidth, imageHeight, 
            discardWhenCavityIsSmallerThanBorder);
        
        assertTrue(result.getN() == expected.size());
        
        for (int i = 0; i < result.getN(); ++i) {
            PairInt p = new PairInt(result.getX(i), result.getY(i));
            assertTrue(expected.remove(p));
        }
        
        assertTrue(expected.isEmpty());
    }
    
}
