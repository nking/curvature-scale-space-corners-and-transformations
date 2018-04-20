package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
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
public class HCPTTest extends TestCase {
    
    public HCPTTest() {
    }

    public void est0() throws IOException {
        
        int cellSize = 16;
        int nBins = 12;
        int maxDimension = 256;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        String filePath0 = ResourceFinder.findFileInTestResources(
            "android_statues_04.jpg");
        ImageExt img0 = ImageIOHelper.readImageExt(filePath0);
        int w0 = img0.getWidth();
        int h0 = img0.getHeight();
        int binFactor0 = (int) Math.ceil(Math.max(
            (float) w0 / maxDimension,
            (float) h0 / maxDimension));
        img0 = imageProcessor.binImage(img0, binFactor0);
        GreyscaleImage pt0 = imageProcessor
            .createCIELUVTheta_WideRangeLightness(img0, 255);
        
        //MiscDebug.writeImage(pt0, "_lab0_");
        
        String filePath1 = ResourceFinder.findFileInTestResources(
            "android_statues_02.jpg");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        int w1 = img1.getWidth();
        int h1 = img1.getHeight();
        int binFactor1 = (int) Math.ceil(Math.max(
            (float) w1 / maxDimension,
            (float) h1 / maxDimension));
        img1 = imageProcessor.binImage(img1, binFactor1);
        GreyscaleImage pt1 = imageProcessor
            .createCIELUVTheta_WideRangeLightness(img1, 255);
        
        //MiscDebug.writeImage(pt1, "_lab1_");
        
        HCPT hgs0 = new HCPT(pt0, 1, cellSize, nBins);
        
        HCPT hgs1 = new HCPT(pt1, 1, cellSize, nBins);
         
        int[] feature0 = new int[nBins];
        hgs0.extractBlock(93, 110, feature0);
        
        int[] feature1 = new int[nBins];
        hgs1.extractBlock(57, 65, feature1);
        
        System.out.println("0=" + Arrays.toString(feature0));
        System.out.println("1=" + Arrays.toString(feature1));
                
        float intersection01 = hgs0.intersection(feature0, feature1);
       
        System.out.println("intersection 0:1=" + intersection01);
                
        // 0 and 1 should be similar
        assertTrue(intersection01 > 0.5);
   
        
        feature1 = new int[nBins];
        hgs1.extractBlock(177, 59, feature1);
        intersection01 = hgs0.intersection(feature0, feature1);
        System.out.println("intersection 0:1=" + intersection01);
        assertTrue(intersection01 < 0.2);
        
        //-----
        
        feature0 = new int[nBins];
        hgs0.extractBlock(153, 91, feature0);
        
        feature1 = new int[nBins];
        hgs1.extractBlock(110, 55, feature1);
        intersection01 = hgs0.intersection(feature0, feature1);
        System.out.println("intersection 0:1=" + intersection01);
        
        feature0 = new int[nBins];
        hgs0.extractBlock(158, 107, feature0);
        
        feature1 = new int[nBins];
        hgs1.extractBlock(118, 71, feature1);
        intersection01 = hgs0.intersection(feature0, feature1);
        System.out.println("intersection 0:1=" + intersection01);
        
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
        //int angle = 20;
        int nBins = 8;//180/angle;
        float angle = 180.f/nBins;
          
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
        
        HCPT hogs1 = null;
        HCPT hogs1_rot = null;
        HCPT hogs2 = null;
        HCPT hogs2_rot = null;
        //HOGs hogs3 = null;
        HCPT hogs4 = null;
        HCPT hogs5 = null;
        
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
            
            hogs1 = new HCPT(img1, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
            hogs1_rot = new HCPT(img1_rot, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
            
            //hogs1._printHistograms();
            //System.out.println("ROTATED: " + img1.getWidth() + "," + img1.getHeight());
            //hogs1_rot._printHistograms();
            
            block1 = new int[nBins];
            hogs1.extractBlock((w/2), (h/2), block1);
            block1_rot = new int[nBins];
            hogs1_rot.extractBlock(w/2, (h/2)-1, block1_rot);
            
            System.out.println(Arrays.toString(block1));
            System.out.println(Arrays.toString(block1_rot));
        
            intersection_1 = hogs1.intersection(block1, 
                block1_rot);
        
            ssd_1 = hogs1.diff(block1, block1_rot);
                    
            System.out.format("intersection1=%.3f\n", intersection_1);
        
            System.out.format("ssd1=%.3f(%.3f)\n", ssd_1[0], ssd_1[1]);
        
            feature1 = hogs1.extractFeature(w/2, h/2, dw, dh);
            feature1_rot = hogs1_rot.extractFeature(w/2, (h/2)-1, dw, dh);
        }
        
        if (useImg2) {
            hogs2 = new HCPT(img2, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
            hogs2_rot = new HCPT(img2_rot, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
            block2 = new int[nBins];
            hogs2.extractBlock(w/2, h/2, block2);
            block2_rot = new int[nBins];
            hogs2_rot.extractBlock(w/2, (h/2)-1, block2_rot);
            System.out.println(Arrays.toString(block2));
            System.out.println(Arrays.toString(block2_rot));
        
            intersection_2 = hogs2.intersection(block2, 
                block2_rot);
        
            ssd_2 = hogs2.diff(block2, block2_rot);
                
            System.out.format("intersection2=%.3f\n", intersection_2);
        
            System.out.format("ssd2=%.3f(%.3f)\n", ssd_2[0], ssd_2[1]);
        
            feature2 = hogs2.extractFeature(w/2, h/2, dw, dh);
            feature2_rot = hogs2_rot.extractFeature(w/2, (h/2)-1, dw, dh);
        }
        
        //hogs3 = new HOGs(img3, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
        //feature3 = hogs3.extractFeature(w/2, h/2, dw, dh);
        hogs4 = new HCPT(img4, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
        feature4 = hogs4.extractFeature(w/2, h/2, dw, dh);
        hogs5 = new HCPT(img5, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
        feature5 = hogs5.extractFeature(w/2, h/2, dw, dh);
        
        if (useImg1 && useImg2) {
                        
            intersection_false_1 = hogs1.intersection(block1,
                block2);
            intersection_false_2 = hogs1.intersection(block1,
                block2_rot);
            intersection_false_3 = hogs1_rot.intersection(block1_rot,
                block2);
            intersection_false_4 = hogs1_rot.intersection(block1_rot,
                block2_rot);

            ssd_false_1 = hogs1.diff(block1, block2);
            ssd_false_2 = hogs1.diff(block1, block2_rot);
            ssd_false_3 = hogs1_rot.diff(block1_rot, block2);
            ssd_false_4 = hogs1_rot.diff(block1_rot, block2_rot);

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
