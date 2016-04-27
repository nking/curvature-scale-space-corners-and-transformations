package algorithms.imageProcessing;

import algorithms.imageProcessing.features.BlobPerimeterCornerHelper;
import algorithms.imageProcessing.features.CornerRegion;
import algorithms.misc.Histogram;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static org.junit.Assert.*;

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
  //add the test images from their paper...should be findable on berkely site i think
           // "blox.gif", "lab.gif", "house.gif", "seattle.jpg", 
            "merton_college_I_001.jpg",
           // "susan-in_plus.png", "lena.jpg",
           // "campus_010.jpg", 
           //"android_statues_01.jpg", 
           //"android_statues_02.jpg", "android_statues_03.jpg", "android_statues_04.jpg"
        };

        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        for (String fileName : fileNames) {
            
            System.out.println("fileName=" + fileName);
        
            String filePath = ResourceFinder.findFileInTestResources(fileName);
        
            ImageExt img = ImageIOHelper.readImageExt(filePath);
        
            List<Set<PairInt>> segmentedCells = 
                imageSegmentation.createColorEdgeSegmentation(img);
            
            ImageExt img2 = ImageIOHelper.readImageExt(filePath);
            
            ImageIOHelper.addAlternatingColorPointSetsToImage(segmentedCells, 
                0, 0, 2, img2);
            
            MiscDebug.writeImage(img2, "_segmented_");
        }
    }
    
}
