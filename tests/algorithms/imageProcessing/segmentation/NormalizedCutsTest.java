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

    public void testManyImages() throws Exception {

        ImageSegmentation imageSegmentation = new ImageSegmentation();

        ImageProcessor imageProcessor = new ImageProcessor();

        String[] fileNames = new String[] {
            "android_statues_01.jpg",
            "android_statues_02.jpg",
            "android_statues_03.jpg",
            "android_statues_04.jpg",
            "seattle.jpg",
            "stonehenge.jpg",
            "cloudy_san_jose.jpg",
            "patagonia_snowy_foreground.jpg",
            "mt_rainier_snowy_field.jpg",
            "brown_lowe_2003_image1.jpg",
            "brown_lowe_2003_image2.jpg",
            "venturi_mountain_j6_0001.png",
            "venturi_mountain_j6_0010.png",
            "campus_010.jpg",
            "campus_011.jpg",
            "merton_college_I_001.jpg",
            "merton_college_I_002.jpg",
            "arches.jpg",
            "stinson_beach.jpg",
            "norwegian_mtn_range.jpg",
            "halfdome.jpg",
            "halfdome2.jpg",
            "halfdome3.jpg",
            "costa_rica.jpg",
            "new-mexico-sunrise_w725_h490.jpg",
            "arizona-sunrise-1342919937GHz.jpg",
            "sky_with_rainbow.jpg",
            "sky_with_rainbow2.jpg",
            "books_illum3_v0_695x555.png",
            "books_illum3_v6_695x555.png",
            "klein_matterhorn_snowy_foreground.jpg",
            "30.jpg",
            "arches_sun_01.jpg",
            "stlouis_arch.jpg",
            "contrail.jpg",
            "checkerboard_01.jpg",
            "checkerboard_02.jpg" };

        for (int i = 0; i < fileNames.length; ++i) {

            String fileName = fileNames[i];

            System.out.println("fileName=" + fileName);

            String filePath = ResourceFinder.findFileInTestResources(fileName);

            String fileNameRoot = fileName.substring(0, fileName.lastIndexOf("."));

            ImageExt img = ImageIOHelper.readImageExt(filePath);

            int w = img.getWidth();
            int h = img.getHeight();
            float maxDimension = 128.f;
            int binFactor = (int) Math.ceil(Math.max((float) w / maxDimension,
                (float) h / maxDimension));
            
            MedianTransform mt = new MedianTransform();
            List<ImageExt> transformed = new ArrayList<ImageExt>();
            mt.<ImageExt>multiscalePyramidalMedianTransform2(img, transformed);
            ImageExt binnedImage = null;
            /*for (int j = 0; j < transformed.size(); ++j) {
                ImageExt tr = transformed.get(j);
                if (tr.getWidth() <= maxDimension && tr.getHeight() <= maxDimension) {
                    binnedImage = tr;
                    break;
                }
            }*/

            binnedImage = imageProcessor.binImage(img, binFactor);
 
            imageSegmentation.applySuperPixelsAndNormalizedCuts(binnedImage);

            MiscDebug.writeImage(binnedImage,  "_ncuts_" + fileNameRoot
                + "_final");
        }
    }
   
}
