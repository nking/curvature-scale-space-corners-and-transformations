package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.Arrays;
import junit.framework.TestCase;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class HOGsTest extends TestCase {
    
    public HOGsTest() {
    }

    public void test0() throws IOException {
        
        int nBins = 9;
        
        // make a test with spiral.png
        String filePath = ResourceFinder.findFileInTestResources("spiral.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
        
        HOGs hogs = new HOGs(img);
        //hogs.setToDebug();
         
         // end of swirl
        int[] feature0 = new int[nBins];
        hogs.extractFeature(65, 122, feature0);
        
        // diagonal edge of swirl
        int[] feature1 = new int[nBins];
        hogs.extractFeature(114, 116, feature1);
    
        // middle of swirl
        int[] feature2 = new int[nBins];
        hogs.extractFeature(72, 68, feature2);
            
        // flat edge of solid swirl
        int[] feature3 = new int[nBins];
        hogs.extractFeature(75, 35, feature3);
        
        // in between arms in swirl
        int[] feature4 = new int[nBins];
        hogs.extractFeature(72, 54, feature4);
      
        System.out.println("0=" + Arrays.toString(feature0));
        System.out.println("1=" + Arrays.toString(feature1));
        System.out.println("2=" + Arrays.toString(feature2));
        System.out.println("3=" + Arrays.toString(feature3));
        System.out.println("4=" + Arrays.toString(feature4));
        
        /* with cell size=6X6:
            0=[73440, 0, 0, 0, 0, 0, 0, 0, 0]       // end of swirl
    [junit] 1=[0, 0, 0, 0, 73440, 0, 0, 0, 0]       // diagonal
    [junit] 2=[48404, 220, 0, 0, 24816, 0, 0, 0, 0] // middle
    [junit] 3=[73440, 0, 0, 0, 0, 0, 0, 0, 0]       // flat
    [junit] 4=[0, 0, 0, 0, 73440, 0, 0, 0, 0]       // in between arms
        */
        
        int orientation0 = 0;
        int orientation1 = 135;
        int orientation3 = 90;
        
        float intersection01 = hogs.intersection(feature0, orientation0, 
            feature1, orientation1);
        
        float intersection13 = hogs.intersection(feature1, orientation1, 
            feature3, orientation3);
        
        System.out.println("intersection 0:1=" + intersection01);
        System.out.println("intersection 1:3=" + intersection13);
        
        // assert that a keypoint on a flat edge of solid is similar
        assertTrue(intersection13 < 0.1);
        
        // assert that a keypoint on end of a point is different than one on
        //   the edge of a flat solid:
        assertTrue(intersection01 > 0.3);
    }
    
}
