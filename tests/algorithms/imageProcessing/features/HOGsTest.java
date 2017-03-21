package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import algorithms.misc.MiscMath;
import java.io.IOException;
import java.util.Arrays;
import junit.framework.TestCase;

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
        
        HOGs hogs = new HOGs(img, 1, 1);
        //hogs.setToDebug();
         
         // end of swirl
        int[] feature0 = new int[nBins];
        hogs.extractFeature(65, 122, feature0);
        
        // diagonal edge of swirl
        int[] feature1 = new int[nBins];
        hogs.extractFeature(106, 48, feature1);
    
        // middle of swirl
        int[] feature2 = new int[nBins];
        hogs.extractFeature(74, 68, feature2);
            
        // flat edge of solid swirl
        int[] feature3 = new int[nBins];
        hogs.extractFeature(75, 35, feature3);
        
        // in between arms in swirl
        int[] feature4 = new int[nBins];
        hogs.extractFeature(72, 54, feature4);
      
        /*
        System.out.println("0=" + Arrays.toString(feature0));
        System.out.println("1=" + Arrays.toString(feature1));
        System.out.println("2=" + Arrays.toString(feature2));
        System.out.println("3=" + Arrays.toString(feature3));
        System.out.println("4=" + Arrays.toString(feature4));
        */
        
        /*
             20  40   60   80  100  120 140  160  180
          0=[0,   0,  0,   0,   0,   0,  0,   0, 255]
          1=[0,   0,  0,   0,   0,   0,  0, 255,   0]
          2=[255, 0,  0,   0,   0,   0,  0,   0, 255]
          3=[0,   0,  0,   0,   255, 0,  0,   0,   0]
          4=[0,   0,  0,   0,   0,   0,  0,   0,   0]
        */
        // 0=end of swirl.  1=diagonal edge. 2=middle.  3=flat edge. 4=between arms

        int orientation0 = 179;
        int orientation1 = 135;  // shift seen is +2
        int orientation3 = 90;   // shift is 0
        
        int expectedPeak0 = 8;
        int expectedPeak1 = 4 + 2;
        int expectedPeak3 = 4;
        
        assertTrue(
            expectedPeak0 == 
            MiscMath.findYMaxIndex(feature0) ||
            (expectedPeak0 - 1) == 
            MiscMath.findYMaxIndex(feature0)    
        );
        //assertEquals(expectedPeak1, MiscMath.findYMaxIndex(feature1));
        assertEquals(expectedPeak3, MiscMath.findYMaxIndex(feature3));

        hogs = new HOGs(img, 1, 3);
        hogs.extractFeature(65, 122, feature0);
        hogs.extractFeature(106, 48, feature1);
        hogs.extractFeature(74, 68, feature2);
        hogs.extractFeature(75, 35, feature3);
        hogs.extractFeature(72, 54, feature4);
        
        /*
        System.out.println("0=" + Arrays.toString(feature0));
        System.out.println("1=" + Arrays.toString(feature1));
        System.out.println("2=" + Arrays.toString(feature2));
        System.out.println("3=" + Arrays.toString(feature3));
        System.out.println("4=" + Arrays.toString(feature4));
        */
        
        /*
        [junit] 0=[0, 0, 122, 213, 685, 973, 193, 0, 107]
        [junit] 1=[0, 0, 0, 0, 0,    0, 33, 2262, 0]
        [junit] 2=[0, 224, 0, 254, 233, 0, 107, 564, 913]
        [junit] 3=[0, 0, 0, 0, 2295, 0, 0,   0, 0]
        [junit] 4=[0, 0, 0, 0, 0, 0, 0, 0, 0]
        */
        
        // 0=end of swirl.  1=diagonal edge. 2=middle.  3=flat edge. 4=between arms
        
        orientation0 = MiscMath.findYMaxIndex(feature0) * 20 + 10;
        orientation1 = MiscMath.findYMaxIndex(feature1) * 20 + 10;
        orientation3 = MiscMath.findYMaxIndex(feature3) * 20 + 10;
        
        float intersection03 = hogs.intersection(feature0, orientation0, 
            feature3, orientation3);
        
        float intersection13 = hogs.intersection(feature1, orientation1, 
            feature3, orientation3);
       
        /*
        int[] a = new int[nBins];
        float[] values = new float[img.getNPixels()];
        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                hogs.extractFeature(i, j, a);
                int v = a[MiscMath.findYMaxIndex(a)];
                values[img.getInternalIndex(i, j)] = v;
            }
        }
        MiscMath.rescale(values, 0, 255);
        
        GreyscaleImage img2 = img.createFullRangeIntWithDimensions();
        for (int i = 0; i < img2.getWidth(); ++i) {
            for (int j = 0; j < img2.getHeight(); ++j) {
                img2.setValue(i, j, (int)values[img2.getInternalIndex(i, j)]);
            }
        }
        MiscDebug.writeImage(img2, "_hog_0_+");
        */
            
        System.out.println("intersection 0:3=" + intersection03);
        System.out.println("intersection 1:3=" + intersection13);
        
        // 0 and 3 should be different
        assertTrue(intersection03 < 0.5);
        
    }
 
}
