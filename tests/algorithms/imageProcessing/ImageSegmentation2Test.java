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
import java.util.ArrayList;
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
           "seattle.jpg", 
           "tmp2.png",
           //"susan-in_plus.png", 
           "lena.jpg",
           "campus_010.jpg", 
           "android_statues_01.jpg", 
           "android_statues_02.jpg", 
            "android_statues_03.jpg", 
            "android_statues_04.jpg"
        };

        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        for (String fileName : fileNames) {
            
            System.out.println("fileName=" + fileName);
        
            String filePath = ResourceFinder.findFileInTestResources(fileName);
        
            String fileNameRoot = fileName.substring(0, fileName.lastIndexOf("."));
            
            ImageExt img = ImageIOHelper.readImageExt(filePath);
            
            int w = img.getWidth();
            int h = img.getHeight();
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
        
            List<Set<PairInt>> segmentedCells = 
                imageSegmentation.createColorEdgeSegmentation(img,
                    fileNameRoot);
            
            ImageExt img2 = ImageIOHelper.readImageExt(filePath);
            
            ImageIOHelper.addAlternatingColorPointSetsToImage(segmentedCells, 
                0, 0, 2, img2);
            
            MiscDebug.writeImage(img2, "_segmented_");
            
            /*
            TODO: here, need an evaluator to compare most f the content
             of descrptors with expected.
            
            the code that finds the best parameters to minimize the difference
            between expected and found will use centroids and number of points
            of segmented cell list.  (might add color terms later if that 
            helps move the evaluator towards better solution faster).
            
            */
        }
    }
    
}
