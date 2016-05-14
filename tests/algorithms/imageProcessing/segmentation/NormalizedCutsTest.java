package algorithms.imageProcessing.segmentation;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class NormalizedCutsTest extends TestCase {
    
    public NormalizedCutsTest() {
    }

    public void testNormalizedCut() throws Exception {
            
        String[] fileNames = new String[]{
            //"color_squares.png",
            "android_statues_01.jpg"
        };
        int[] kCells = new int[] {
            //9,
            400
        };
        
        for (int i = 0; i < fileNames.length; ++i) {

            String fileName = fileNames[i];

            int kCell = kCells[i];

            System.out.println("fileName=" + fileName);

            String filePath = ResourceFinder.findFileInTestResources(fileName);

            String fileNameRoot = fileName.substring(0, fileName.lastIndexOf("."));

            ImageExt img = ImageIOHelper.readImageExt(filePath);

            SLICSuperPixels slic = new SLICSuperPixels(img, kCell);
            slic.calculate();
            int[] labels = slic.getLabels();

            System.out.println("have initial labels (super pixels)");
            System.out.flush();

            ImageExt img2 = img.copyToImageExt();
            ImageIOHelper.addAlternatingColorLabelsToRegion(img2, labels);
            MiscDebug.writeImage(img2,  "_slic_" + fileNameRoot);
            int w = img.getWidth();
            int h = img.getHeight();
            //int n = img.getNPixels();
    
            NormalizedCuts normCuts = new NormalizedCuts();
            int[] labels2 = normCuts.normalizedCut(img, labels);

            ImageIOHelper.addAlternatingColorLabelsToRegion(img2, labels2);
            MiscDebug.writeImage(img2,  "_ncuts_" + fileNameRoot);
            
            System.out.println("labels2=" + Arrays.toString(labels2));
        }
    }
    
}
