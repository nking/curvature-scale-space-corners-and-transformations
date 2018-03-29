package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
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

    public void est0() throws IOException {
        
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
 
    public void test1() throws IOException {
        
        boolean useImg1 = true;
        boolean useImg2 = true;
        
        /* 
        226,160
        196:256, 190:130
        */
        int w = 60;
        int h = 60;
        int xc1 = 226;
        int yc1 = 160;
        int xc2 = 226;
        int yc2 = 102;
        
        int N_PIX_PER_CELL_DIM = 3;    //1;  3;
        int N_CELLS_PER_BLOCK_DIM = 3; //10; 3;
        int angle = 20;
        int nBins = 180/angle;
          
        Transformer transformer = new Transformer();
        TransformationParameters params = new TransformationParameters();
        params.setRotationInDegrees(90);
        params.setScale(1);
        params.setOriginX(0);
        params.setOriginY(0);
    
        // make a test with spiral.png
        String filePath = ResourceFinder.findFileInTestResources(
                "android_statues_objects.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
        
        GreyscaleImage img1 = null;
        GreyscaleImage img1_rot = null;
        GreyscaleImage img2 = null;
        GreyscaleImage img2_rot = null;
        HOGs hogs1 = null;
        HOGs hogs1_rot = null;
        HOGs hogs2 = null;
        HOGs hogs2_rot = null;
        int[] feature1 = null;
        int[] feature1_rot = null;
        int[] feature2 = null;
        int[] feature2_rot = null;
        
        int orientation1 = -1;
        int orientation1_rot = -1;
        int orientation2 = -1;
        int orientation2_rot = -1;
        float intersection_1 = -1;
        float intersection_2 = -1;
        float intersection_false_1 = -1;
        float intersection_false_2 = -1;
        float intersection_false_3 = -1;
        float intersection_false_4 = -1;
        float[] ssd_1 = null;
        float[] ssd_2 = null;
        float[] ssd_false_1 = null;
        float[] ssd_false_2 = null;
        float[] ssd_false_3 = null;
        float[] ssd_false_4 = null;
        
        if (useImg1) {
            img1 = img.subImage(xc1, yc1, w, h);    
            
            params.setTranslationX(0);
            params.setTranslationY(img1.getWidth() - 1);
            img1_rot = transformer.applyTransformation(img1, params, 
                img1.getHeight(), img1.getWidth());
        }
        
        if (useImg2) {
            img2 = img.subImage(xc2, yc2, w, h);
            params.setTranslationX(0);
            params.setTranslationY(img2.getWidth() - 1);
            img2_rot = transformer.applyTransformation(img2, params, 
                img2.getHeight(), img2.getWidth());
        }
        
        //MiscDebug.writeImage(img1, "img1");
        //MiscDebug.writeImage(img1_rot, "img1_rot");
        //MiscDebug.writeImage(img2, "img2");
        //MiscDebug.writeImage(img2_rot, "img2_rot");
        
        if (useImg1) {
            
            hogs1 = new HOGs(img1, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
            hogs1_rot = new HOGs(img1_rot, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
            
            //hogs1._printHistograms();
            //System.out.println("ROTATED: " + img1.getWidth() + "," + img1.getHeight());
            //hogs1_rot._printHistograms();
            
            feature1 = new int[nBins];
            hogs1.extractFeature((w/2), (h/2), feature1);
            feature1_rot = new int[nBins];
            hogs1_rot.extractFeature(w/2, h/2, feature1_rot);
            
            System.out.println(Arrays.toString(feature1));
            System.out.println(Arrays.toString(feature1_rot));
        
            orientation1 = MiscMath.findYMaxIndex(feature1) * angle + (angle/2);
            orientation1_rot = MiscMath.findYMaxIndex(feature1_rot) * angle + (angle/2);
        
            intersection_1 = hogs1.intersection(feature1, orientation1, 
            feature1_rot, orientation1_rot);
        
            ssd_1 = hogs1.diff(feature1, orientation1, feature1_rot, orientation1_rot);
            
            System.out.format("ang1=%d ang1rot90=%d\n", orientation1, orientation1_rot);
        
            System.out.format("intersection1=%.3f\n", intersection_1);
        
            System.out.format("ssd1=%.3f(%.3f)\n", ssd_1[0], ssd_1[1]);
        }
        
        if (useImg2) {
            hogs2 = new HOGs(img2, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
            hogs2_rot = new HOGs(img2_rot, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
            feature2 = new int[nBins];
            hogs2.extractFeature(w/2, h/2, feature2);
            feature2_rot = new int[nBins];
            hogs2_rot.extractFeature(w/2, h/2, feature2_rot);
            System.out.println(Arrays.toString(feature2));
            System.out.println(Arrays.toString(feature2_rot));
        
            orientation2 = MiscMath.findYMaxIndex(feature2) * angle + (angle/2);
            orientation2_rot = MiscMath.findYMaxIndex(feature2_rot) * angle + (angle/2);
        
            intersection_2 = hogs2.intersection(feature2, orientation2, 
                feature2_rot, orientation2_rot);
        
            ssd_2 = hogs2.diff(feature2, orientation2, feature2_rot, orientation2_rot);
        
            System.out.format("ang2=%d ang2rot90=%d\n", orientation2, orientation2_rot);
        
            System.out.format("intersection2=%.3f\n", intersection_2);
        
            System.out.format("ssd2=%.3f(%.3f)\n", ssd_2[0], ssd_2[1]);
        }
        
        if (useImg1 && useImg2) {
            
            intersection_false_1 = hogs1.intersection(feature1, orientation1,
                    feature2, orientation2);
            intersection_false_2 = hogs1.intersection(feature1, orientation1,
                    feature2_rot, orientation2_rot);
            intersection_false_3 = hogs1_rot.intersection(feature1_rot, orientation1_rot,
                    feature2, orientation2);
            intersection_false_4 = hogs1_rot.intersection(feature1_rot, orientation1_rot,
                    feature2_rot, orientation2_rot);

            ssd_false_1 = hogs1.diff(feature1, orientation1, feature2, orientation2);
            ssd_false_2 = hogs1.diff(feature1, orientation1, feature2_rot, orientation2_rot);
            ssd_false_3 = hogs1_rot.diff(feature1_rot, orientation1_rot, feature2, orientation2);
            ssd_false_4 = hogs1_rot.diff(feature1_rot, orientation1_rot, feature2_rot, orientation2_rot);

            System.out.format("false intersection=%.3f, %.3f, %.3f, %.3f\n",
                intersection_false_1, intersection_false_2,
                intersection_false_3, intersection_false_4);
            
            System.out.format("false ssds=%.3f(%.3f), %.3f(%.3f), %.3f(%.3f), %.3f(%.3f)\n",
                    ssd_false_1[0], ssd_false_1[1],
                    ssd_false_2[0], ssd_false_2[1],
                    ssd_false_3[0], ssd_false_3[1],
                    ssd_false_4[0], ssd_false_4[1]);
        }
    }
}
