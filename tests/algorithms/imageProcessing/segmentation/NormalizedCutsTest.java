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
            
        // 5 X 5 image
        String fileName = "color_squares.png";
                 
        System.out.println("fileName=" + fileName);
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        String fileNameRoot = fileName.substring(0, fileName.lastIndexOf("."));
            
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        SLICSuperPixels slic = new SLICSuperPixels(img, 9);
        slic.calculate();
        int[] pixAssignments = slic.getLabels();
        
        System.out.println("have initial labels (super pixels)");
        System.out.flush();
        
        ImageExt img2 = img.copyToImageExt();
        ImageIOHelper.addAlternatingColorLabelsToRegion(img2, pixAssignments);
        MiscDebug.writeImage(img2,  "_slic_" + fileNameRoot);
        
        int w = img.getWidth();
        int h = img.getHeight();
        //int n = img.getNPixels();
        
        int nMaxLabel = Integer.MIN_VALUE;
        for (int i = 0; i < pixAssignments.length; ++i) {
            if (pixAssignments[i] > nMaxLabel) {
                nMaxLabel = pixAssignments[i];
            }
        }
        int[][] labels = new int[w][h];
        for (int i = 0; i < w; ++i) {
            labels[i] = new int[h];
        }
        for (int i = 0; i < pixAssignments.length; ++i) {
            int x = img.getCol(i);
            int y = img.getRow(i);
            labels[x][y] = pixAssignments[i];
        }
        
        NormalizedCuts normCuts = new NormalizedCuts();
        int[][] labels2 = normCuts.normalizedCut(img, labels);
        
        System.out.println("labels2=");
        for (int i = 0; i < labels2.length; ++i) {
            System.out.println(Arrays.toString(labels2[i]));
        }
        
        int z = 1;
    }
    
}
