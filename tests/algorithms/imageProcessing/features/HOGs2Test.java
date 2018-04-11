package algorithms.imageProcessing.features;

import algorithms.Rotate;
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
public class HOGs2Test extends TestCase {
    
    public HOGs2Test() {
    }

    public void test1() throws IOException {
        
        int w = 60;
        int h = 60;
        int xc1 = 226;
        int yc1 = 160;
        
        int dw = w-1;
        int dh = h-1;
        
        //int N_PIX_PER_CELL_DIM = 4;
        //int N_CELLS_PER_BLOCK_DIM = 2;
        int N_PIX_PER_CELL_DIM = 3;    //1;  3;
        int N_CELLS_PER_BLOCK_DIM = 3; //10; 3;
        //int angle = 20;
        int nBins = 8;//180/angle;
          
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
        HOGs hogs1 = null;
        HOGs hogs1_rot = null;
        
        img1 = img.subImage(xc1, yc1, w, h);    

        params.setTranslationX(0);
        params.setTranslationY(img1.getWidth() - 1);
        img1_rot = transformer.applyTransformation(img1, params, 
            img1.getHeight(), img1.getWidth());
        
        //MiscDebug.writeImage(img1, "img1");
        //MiscDebug.writeImage(img1_rot, "img1_rot");
        
        hogs1 = new HOGs(img1, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
        hogs1_rot = new HOGs(img1_rot, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, nBins);
        
        int[] histCenter = hogs1.extractHistogram(w/2, h/2);
        
        int[] histCenter_rot = hogs1_rot.extractHistogram((w/2), (h/2)-1);
        
        System.out.println("histCenter=    " + Arrays.toString(histCenter));
        int[] tmp = Arrays.copyOf(histCenter_rot, nBins);
        Rotate rotate = new Rotate();
        rotate.rotate(tmp, 4);
        System.out.println("histCenter_rot=" + Arrays.toString(tmp));
        
        assertTrue(Arrays.equals(histCenter, tmp));
    }
}
