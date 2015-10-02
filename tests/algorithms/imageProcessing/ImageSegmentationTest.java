package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class ImageSegmentationTest extends TestCase {

    public ImageSegmentationTest(String testName) {
        super(testName);
    }
    
    public void testApplyColorSegmentation() throws Exception {
        
        //String fileName = "two_circles_color2.png";
        //String fileName = "two_circles_color.png";
        //String fileName = "cloudy_san_jose.jpg";
        //String fileName = "middlebury_cones_im2.png"; // a limitFrac of 0.1 works well
        //String fileName = "brown_lowe_2003_image1.jpg";
        //String fileName = "venturi_mountain_j6_0010.png";
        String fileName = "books_illum3_v6_695x555.png";
        //String fileName = "brown_lowe_2003_image1.jpg";
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        //HistogramEqualizationForColor hEq = new HistogramEqualizationForColor(img);
        //hEq.applyFilter();
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        List<Set<PairInt>> clusterSets = 
            imageSegmentation.calculateColorSegmentation3(img, /*0.1f,*/ true);
            //imageProcessor.calculateColorSegmentation2(img, true);
        
        int nPoints = count(clusterSets);
        
        int nExtraForDot = 0;
        
        Image img2 = new Image(img.getWidth(), img.getHeight());
        
        for (int i = 0; i < clusterSets.size(); ++i) {
            
            int[] rgb = ImageIOHelper.getNextRGB(i);
            
            Set<PairInt> set = clusterSets.get(i);
            
            ImageIOHelper.addToImage(set, 0, 0, img2, nExtraForDot, rgb[0], 
                rgb[1], rgb[2]);
        }
        String bin = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(bin + "/cluster.png", img2);
        
        
        //assertTrue(nPoints == img.getNPixels());
        
        boolean[] present = new boolean[img.getNPixels()];
        for (Set<PairInt> set : clusterSets) {
            for (PairInt p : set) {
                int pixIdx = img.getInternalIndex(p.getX(), p.getY());
                assertFalse(present[pixIdx]);
                present[pixIdx] = true;
            }
        }
        
        int z = 1;
    }
    
    private int count(List<Set<PairInt>> setList) {
        
        int c = 0;
        for (Set<PairInt> set : setList) {
            c += set.size();
        }
        
        return c;
    }
    
}
