package algorithms.imageProcessing.segmentation;

import algorithms.imageProcessing.*;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.MedianTransform;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.*;
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
            "color_squares.png"
        };
        int[] kCells = new int[] {
            9,
        };
        int[] numbersOfIterations = new int[] {
            1,
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
            
            int[] labels2 = null;
            for (int nIter = 0; nIter < numbersOfIterations[i]; ++nIter) {
    
                NormalizedCuts normCuts = new NormalizedCuts();
                labels2 = normCuts.normalizedCut(img, labels);
                labels = labels2;
                
                System.out.println("labels2=" + Arrays.toString(labels2));
                
                img2 = ImageIOHelper.readImageExt(filePath);
                ImageIOHelper.addAlternatingColorLabelsToRegion(img2, labels2);
                MiscDebug.writeImage(img2,  "_ncuts_" + fileNameRoot + "_" + nIter); 
            } 
            
            img2 = ImageIOHelper.readImageExt(filePath);
            LabelToColorHelper.applyLabels(img2, labels);
            MiscDebug.writeImage(img2,  "_ncuts_" + fileNameRoot + "_final"); 
        }
    }
   
}
