package algorithms.imageProcessing;

import algorithms.imageProcessing.features.HOGUtil;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class IntegralHistogramsTest extends TestCase {
    
    public IntegralHistogramsTest() {
    }
    
    public void testIntegralHistograms() {
        
        /*
        column major notation is used
        
        2  155  2   7
        1  3   16  42 
        0  2    5   9
           0    1   2
         */
        GreyscaleImage img = new GreyscaleImage(3, 3);
        img.setValue(0, 0, 2); img.setValue(0, 1, 3); img.setValue(0, 2, 155);
        img.setValue(1, 0, 5); img.setValue(1, 1, 16); img.setValue(1, 2, 2);
        img.setValue(2, 0, 9); img.setValue(2, 1, 42); img.setValue(2, 2, 7);
       
        IntegralHistograms sumTable = new IntegralHistograms();
        int[][] sHists = sumTable.create(img);
        assert(sHists.length == img.getNPixels());
        assert(sHists[0].length == 16);
        int minValue = 0;
        int maxValue = 255;
        int binWidth = 16;
        
        //TODO:  rewrite these
    }
}
