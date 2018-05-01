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
public class HOGRegionsTest extends TestCase {
    
    public HOGRegionsTest() {
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
        
        //TODO: finish unit test
    }
}
