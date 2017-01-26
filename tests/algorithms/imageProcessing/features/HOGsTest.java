package algorithms.imageProcessing.features;

import algorithms.imageProcessing.EdgeFilterProducts;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import com.climbwithyourfeet.clustering.util.MiscMath;
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
        
        HOGs hogs = new HOGs(img, 1);
        //hogs.setToDebug();
         
         // end of swirl
        int[] feature0 = new int[nBins];
        hogs.extractFeature(65, 122, feature0);
        
        // diagonal edge of swirl
        int[] feature1 = new int[nBins];
        hogs.extractFeature(114, 116, feature1);
    
        // middle of swirl
        int[] feature2 = new int[nBins];
        hogs.extractFeature(74, 68, feature2);
            
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
        
           with cell size=4X4
            0=[6076, 0, 0, 0, 7200, 5084, 0, 0, 0]
    [junit] 1=[9180, 9180, 0, 0, 0, 0, 0, 0, 0]
    [junit] 2=[3032, 10524, 4808, 0, 0, 0, 0, 0, 0]
    [junit] 3=[0, 0, 0, 0, 18360, 0, 0, 0, 0]
    [junit] 4=[0, 0, 0, 13504, 4856, 0, 0, 0, 0]
    [junit] intersection 0:1=0.39215687
    [junit] intersection 1:3=0.0
        
           with cell size=1X1
           0=[0, 0, 0, 0, 0, 2040, 0, 0, 0]
    [junit] 1=[0, 2040, 0, 0, 0, 0, 0, 0, 0]
    [junit] 2=[0, 2040, 0, 0, 0, 0, 0, 0, 0]
    [junit] 3=[0, 0, 0, 0, 2040, 0, 0, 0, 0]
    [junit] 4=[0, 0, 0, 0, 0, 0, 0, 0, 0]
        
        */
        
        int orientation0 = 90 + 20;  // shift seen is +1
        int orientation1 = 30;       // shift seen is -3
        int orientation3 = 90;       // shift is 0
        
        int expectedPeak0 = 4 + 1;
        int expectedPeak1 = 4 - 3;
        int expectedPeak3 = 4;
        
        assertEquals(expectedPeak0, MiscMath.findYMaxIndex(feature0));
        assertEquals(expectedPeak1, MiscMath.findYMaxIndex(feature1));
        assertEquals(expectedPeak3, MiscMath.findYMaxIndex(feature3));
        
        // --- redo for a more typical cell size ----
        hogs = new HOGs(img, 6);
        hogs.extractFeature(65, 122, feature0);
        hogs.extractFeature(114, 116, feature1);
        hogs.extractFeature(72, 68, feature2);
        hogs.extractFeature(75, 35, feature3);        
        hogs.extractFeature(72, 54, feature4);
        
        System.out.println("0=" + Arrays.toString(feature0));
        System.out.println("1=" + Arrays.toString(feature1));
        System.out.println("2=" + Arrays.toString(feature2));
        System.out.println("3=" + Arrays.toString(feature3));
        System.out.println("4=" + Arrays.toString(feature4));
        
        int orientation2 = 30;
      
        // diagonal=1, flat edge=3
        // middle=2, point at end=0
   
        float intersection13 = hogs.intersection(feature1, orientation1, 
            feature3, orientation3);
        
        float intersection23 = hogs.intersection(feature2, orientation2, 
            feature3, orientation3);
        
        System.out.println("intersection 1:3=" + intersection13);
        System.out.println("intersection 2:3=" + intersection23);
        
        // 1 and 3 should be similar
        assertTrue(intersection13 >= 0.5);
        
        // 2 and 3 should be different
        assertTrue(intersection23 < 0.4);
        
        // needs test w/ non-solid image and features
    }
    
    public void test1() throws IOException {
        
        //TODO: need to add a test for summing costs over a window
        //  of features
        
        int nBins = 9;
        
        // make a test with spiral.png
        String filePath = ResourceFinder.findFileInTestResources(
            "android_statues_01_sz1.jpg");
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
        
        HOGs hogs = new HOGs(img, 6);
        hogs.setToDebug();
         
        // middle of cupcake top
        int[] feature0 = new int[nBins];
        hogs.extractFeature(34, 100, feature0);
    
        int[] feature1 = new int[nBins];
        hogs.extractFeature(378, 94, feature1);
        
        System.out.println("0=" + Arrays.toString(feature0));
        System.out.println("1=" + Arrays.toString(feature1));
        
        int orientation0 = 110; 
        int orientation1 = 90;       // shift seen is -3
      
        float intersection01 = hogs.intersection(feature0, orientation0, 
            feature1, orientation1);
        
        float ssd01 = hogs.ssd(feature0, orientation0, 
            feature1, orientation1);
        
        // ssd should be large for higher difference
        System.out.println("intersection 0:1=" + intersection01 + " ssd=" +
            ssd01);
        
        // 0 and 1 should be different, that is, small intrsection
        //assertTrue(intersection01 < 0.4);
        
    }
}
