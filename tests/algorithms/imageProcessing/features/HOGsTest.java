package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
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

    public void est0() throws IOException {
        
        int nBins = 9;
        int N_PIX_PER_CELL_DIM = 1;
        int N_CELLS_PER_BLOCK_DIM = 1;
        
        // make a test with spiral.png
        String filePath = ResourceFinder.findFileInTestResources("spiral.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
        
        HOGs hogs = new HOGs(img, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM);
        //hogs.setToDebug();
         
         // end of swirl
        int[] feature0 = new int[nBins];
        hogs.extractBlock(65, 122, feature0);
        
        // diagonal edge of swirl
        int[] feature1 = new int[nBins];
        hogs.extractBlock(106, 48, feature1);
    
        // middle of swirl
        int[] feature2 = new int[nBins];
        hogs.extractBlock(74, 68, feature2);
            
        // flat edge of solid swirl
        int[] feature3 = new int[nBins];
        hogs.extractBlock(75, 35, feature3);
        
        // in between arms in swirl
        int[] feature4 = new int[nBins];
        hogs.extractBlock(72, 54, feature4);
      
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
        hogs.extractBlock(65, 122, feature0);
        hogs.extractBlock(106, 48, feature1);
        hogs.extractBlock(74, 68, feature2);
        hogs.extractBlock(75, 35, feature3);
        hogs.extractBlock(72, 54, feature4);
        
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
                hogs.extractBlock(i, j, a);
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
        int dw = w-1;
        int dh = h-1;
        
        //int N_PIX_PER_CELL_DIM = 4;
        //int N_CELLS_PER_BLOCK_DIM = 2;
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
        //GreyscaleImage img3 = null;
        GreyscaleImage img4 = null;
        GreyscaleImage img5 = null;
        
        HOGs hogs1 = null;
        HOGs hogs1_rot = null;
        HOGs hogs2 = null;
        HOGs hogs2_rot = null;
        //HOGs hogs3 = null;
        HOGs hogs4 = null;
        HOGs hogs5 = null;
        
        int[] block1 = null;
        int[] block1_rot = null;
        int[] block2 = null;
        int[] block2_rot = null;
        int[] feature1 = null;
        int[] feature1_rot = null;
        int[] feature2 = null;
        int[] feature2_rot = null;
        //int[] feature3 = null;
        int[] feature4 = null;
        int[] feature5 = null;
                
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
        
        //int xc3 = 29;
        //int yc3 = 160;
        //img3 = img.subImage(xc3, yc3, w, h);
        int xc4 = 94;
        int yc4 = 160;
        img4 = img.subImage(xc4, yc4, w, h);
        int xc5 = 159;
        int yc5 = 160;
        img5 = img.subImage(xc5, yc5, w, h);
        
        
        /*
        MiscDebug.writeImage(img1, "img1");
        MiscDebug.writeImage(img1_rot, "img1_rot");
        MiscDebug.writeImage(img2, "img2");
        MiscDebug.writeImage(img2_rot, "img2_rot");
        */
        
        if (useImg1) {
            
            hogs1 = new HOGs(img1, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
            hogs1_rot = new HOGs(img1_rot, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
            
            //hogs1._printHistograms();
            //System.out.println("ROTATED: " + img1.getWidth() + "," + img1.getHeight());
            //hogs1_rot._printHistograms();
            
            block1 = new int[nBins];
            hogs1.extractBlock((w/2), (h/2), block1);
            block1_rot = new int[nBins];
            hogs1_rot.extractBlock(w/2, h/2, block1_rot);
            
            System.out.println(Arrays.toString(block1));
            System.out.println(Arrays.toString(block1_rot));
        
            orientation1 = MiscMath.findYMaxIndex(block1) * angle + (angle/2);
            orientation1_rot = MiscMath.findYMaxIndex(block1_rot) * angle + (angle/2);
        
            intersection_1 = hogs1.intersection(block1, orientation1, 
            block1_rot, orientation1_rot);
        
            ssd_1 = hogs1.diff(block1, orientation1, block1_rot, orientation1_rot);
            
            System.out.format("ang1=%d ang1rot90=%d\n", orientation1, orientation1_rot);
        
            System.out.format("intersection1=%.3f\n", intersection_1);
        
            System.out.format("ssd1=%.3f(%.3f)\n", ssd_1[0], ssd_1[1]);
        
            feature1 = hogs1.extractFeature(w/2, h/2, dw, dh);
            feature1_rot = hogs1_rot.extractFeature(w/2, h/2, dw, dh);
        }
        
        if (useImg2) {
            hogs2 = new HOGs(img2, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
            hogs2_rot = new HOGs(img2_rot, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
            block2 = new int[nBins];
            hogs2.extractBlock(w/2, h/2, block2);
            block2_rot = new int[nBins];
            hogs2_rot.extractBlock(w/2, h/2, block2_rot);
            System.out.println(Arrays.toString(block2));
            System.out.println(Arrays.toString(block2_rot));
        
            orientation2 = MiscMath.findYMaxIndex(block2) * angle + (angle/2);
            orientation2_rot = MiscMath.findYMaxIndex(block2_rot) * angle + (angle/2);
        
            intersection_2 = hogs2.intersection(block2, orientation2, 
                block2_rot, orientation2_rot);
        
            ssd_2 = hogs2.diff(block2, orientation2, block2_rot, orientation2_rot);
        
            System.out.format("ang2=%d ang2rot90=%d\n", orientation2, orientation2_rot);
        
            System.out.format("intersection2=%.3f\n", intersection_2);
        
            System.out.format("ssd2=%.3f(%.3f)\n", ssd_2[0], ssd_2[1]);
        
            feature2 = hogs2.extractFeature(w/2, h/2, dw, dh);
            feature2_rot = hogs2_rot.extractFeature(w/2, h/2, dw, dh);
        }
        
        //hogs3 = new HOGs(img3, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
        //feature3 = hogs3.extractFeature(w/2, h/2, dw, dh);
        hogs4 = new HOGs(img4, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
        feature4 = hogs4.extractFeature(w/2, h/2, dw, dh);
        hogs5 = new HOGs(img5, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
        feature5 = hogs5.extractFeature(w/2, h/2, dw, dh);
        
        if (useImg1 && useImg2) {
                        
            intersection_false_1 = hogs1.intersection(block1, orientation1,
                block2, orientation2);
            intersection_false_2 = hogs1.intersection(block1, orientation1,
                block2_rot, orientation2_rot);
            intersection_false_3 = hogs1_rot.intersection(block1_rot, orientation1_rot,
                block2, orientation2);
            intersection_false_4 = hogs1_rot.intersection(block1_rot, orientation1_rot,
                block2_rot, orientation2_rot);

            ssd_false_1 = hogs1.diff(block1, orientation1, block2, orientation2);
            ssd_false_2 = hogs1.diff(block1, orientation1, block2_rot, orientation2_rot);
            ssd_false_3 = hogs1_rot.diff(block1_rot, orientation1_rot, block2, orientation2);
            ssd_false_4 = hogs1_rot.diff(block1_rot, orientation1_rot, block2_rot, orientation2_rot);

            System.out.format("false intersection=%.3f, %.3f, %.3f, %.3f\n",
                intersection_false_1, intersection_false_2,
                intersection_false_3, intersection_false_4);
            
            System.out.format("false ssds=%.3f(%.3f), %.3f(%.3f), %.3f(%.3f), %.3f(%.3f)\n",
                    ssd_false_1[0], ssd_false_1[1],
                    ssd_false_2[0], ssd_false_2[1],
                    ssd_false_3[0], ssd_false_3[1],
                    ssd_false_4[0], ssd_false_4[1]);
            
            // ------ features ----------
            //need to compare the features that have similar rotations
            float inter2_pos = hogs1_rot.intersectionOfFeatures(feature1, feature4);
            float inter3_pos = hogs1_rot.intersectionOfFeatures(feature1, feature5);
            
            float[] diff2_pos = hogs1_rot.diffOfFeatures(feature1, feature4);
            float[] diff3_pos = hogs1_rot.diffOfFeatures(feature1, feature5);
            
            float inter1_neg = hogs1.intersectionOfFeatures(feature1, feature2_rot);
            float inter2_neg = hogs1_rot.intersectionOfFeatures(feature1_rot, 
                feature2);
            
            float[] diff1_neg = hogs1.diffOfFeatures(feature1, feature2_rot);
            float[] diff2_neg = hogs1_rot.diffOfFeatures(feature1_rot, 
                feature2);
            
            System.out.println("Comparing 'registered' block detector results for features:");
            System.out.format("    positive intersections=%.3f, %.3f\n", 
                inter2_pos, inter3_pos);
            System.out.format("    negative intersections=%.3f, %.3f\n", 
                inter1_neg, inter2_neg);
            
            System.out.format("    positive differences=(%.3f+-%.3f), (%.3f+-%.3f)\n", 
                diff2_pos[0], diff2_pos[1], diff3_pos[0], diff3_pos[1]);
            System.out.format("    negative differences=(%.3f+-%.3f), (%.3f+-%.3f)\n", 
                diff1_neg[0], diff1_neg[1], diff2_neg[0], diff2_neg[1]);
        }
    }
}
