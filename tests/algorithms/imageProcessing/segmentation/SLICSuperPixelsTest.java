package algorithms.imageProcessing.segmentation;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SLICSuperPixelsTest extends TestCase {
    
    public SLICSuperPixelsTest() {
    }
    
    public void test0() throws Exception {
        
        String[] fileNames = new String[]{
            "android_statues_02.jpg", 
            "checkerboard_subimage.jpg",
            //"tmp3.png"
        };
        
        int[] kCells = new int[]{
            200,
            4
        };
        
        for (int i = 0; i < fileNames.length; ++i) {
                         
            String fileName = fileNames[i];
            
            System.out.println("fileName=" + fileName);

            String filePath = ResourceFinder.findFileInTestResources(fileName);

            String fileNameRoot = fileName.substring(0, fileName.lastIndexOf("."));

            ImageExt img = ImageIOHelper.readImageExt(filePath);

            SLICSuperPixels slic = new SLICSuperPixels(img, kCells[i]);

            slic.calculate();

            int[] labels = slic.getLabels();

            ImageIOHelper.addAlternatingColorLabelsToRegion(img, labels);

            MiscDebug.writeImage(img,  "_slic_" + fileNameRoot);
        }
    }
}
