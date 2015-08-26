package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
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
    
    public void testExtractAndOrderTheBorder0_1() throws Exception {
        
        String fileName = "blobs1_236_123.dat";
        
        String filePath = ResourceFinder.findDirectory("testresources") + "/" 
            + fileName;
        
        Set<PairInt> blob = Misc.deserializeSetPairInt(filePath);
        
        assertNotNull(blob);
        assertFalse(blob.isEmpty());
        
        int[] minMaxXY = MiscMath.findMinMaxXY(blob);

        MiscDebug.plotPoints(blob, minMaxXY[1] + 1, minMaxXY[3] + 1,
            MiscDebug.getCurrentTimeFormatted());        

        boolean discardWhenCavityIsSmallerThanBorder = true;
        
        EdgeExtractorForBlobBorder instance = new EdgeExtractorForBlobBorder();
                
        PairIntArray result = instance.extractAndOrderTheBorder0(
            blob, minMaxXY[1] + 1, minMaxXY[3] + 1, 
            discardWhenCavityIsSmallerThanBorder);
        
        Image img3 = new Image(minMaxXY[1] + 1, minMaxXY[3] + 1);
        for (int j = 0; j < result.getN(); ++j) {
            int x = result.getX(j);
            int y = result.getY(j);
            if (j == 0 || (j == (result.getN() - 1))) {
                ImageIOHelper.addPointToImage(x, y, img3, 0, 200, 150, 0);
            } else {
                ImageIOHelper.addPointToImage(x, y, img3, 0, 255, 0, 0);
            }
        }
        MiscDebug.writeImageCopy(img3, "blob_ordered_perimeter_" + MiscDebug.getCurrentTimeFormatted() + ".png");

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        assertTrue(curveHelper.isAdjacent(result, 0, result.getN() - 1));

    }
}
