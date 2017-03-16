package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ImageSegmentation2Test extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public ImageSegmentation2Test(String testName) {
        super(testName);
    }
    
    public void test0() throws Exception {

        String[] fileNames = new String[]{
           "seattle.jpg", 
         //  "tmp2.png", //snapshot from the paper
           "susan-in_plus.png", 
         /*  "lena.jpg",
           "campus_010.jpg", 
           "android_statues_01.jpg", 
           "android_statues_02.jpg", 
            "android_statues_03.jpg", 
            "android_statues_04.jpg"
          */  
        };

        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        for (String fileName : fileNames) {
            
            System.out.println("fileName=" + fileName);
        
            String filePath = ResourceFinder.findFileInTestResources(fileName);
        
            String fileNameRoot = fileName.substring(0, fileName.lastIndexOf("."));
            
            ImageExt img = ImageIOHelper.readImageExt(filePath);
            int w = img.getWidth();
            int h = img.getHeight();
            
            int binFactor = (int) Math.ceil(Math.max(
                (float) w / 256,
                (float) h / 256));

            ImageProcessor imageProcessor = new ImageProcessor();
            
            img = imageProcessor.binImage(img, binFactor);
            
            w = img.getWidth();
            h = img.getHeight();
            List<ImageExt> transformed = null;
            int selectIdx = -1;
            int minDimension = 512;
            if (w > minDimension || h > minDimension) {                
                MedianTransform mt = new MedianTransform();
                transformed = new ArrayList<ImageExt>();
                mt.<ImageExt>multiscalePyramidalMedianTransform2(img, transformed);
                for (int j = 0; j < transformed.size(); ++j) {
                    ImageExt tr = transformed.get(j);
                    if (selectIdx == -1) {
                        if (tr.getWidth() <= minDimension && tr.getHeight() <= minDimension) {
                            selectIdx = j;
                            img = transformed.get(selectIdx);
                            break;
                        }
                    }
                }                
            }
        
            ImageExt img2 = img.copyToImageExt();
            
            List<Set<PairInt>> segmentedCells = 
                imageSegmentation.createColorEdgeSegmentation(img,
                    fileNameRoot);
                        
            ImageIOHelper.addAlternatingColorPointSetsToImage(segmentedCells, 
                0, 0, 2, img2);
            
            MiscDebug.writeImage(img2, "_segmented_" + fileNameRoot);
            
        }
    }
    
}
