package algorithms.imageProcessing.segmentation;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.optimization.segmentation.SData;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class NormalizedCutsTest extends TestCase {
    
    public NormalizedCutsTest() {
    }
    
    public void test0() throws Exception {
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        SData[] trainingData = getDetailedTrainingData();
        
        for (int i = 0; i < trainingData.length; ++i) {
        
            String filePath = trainingData[i].dirPath + "/" + 
                trainingData[i].imgFileName;
        
            ImageExt img = ImageIOHelper.readImageExt(filePath);
            
            System.out.println("filePath=" + filePath);
            
            imageSegmentation.applySuperPixelsAndNormalizedCuts(img);
            
            MiscDebug.writeImage(img,  "_ncuts_" + trainingData[i].imgFileName 
                + "_final"); 
        }
        
    }

    public void estNormalizedCut() throws Exception {
            
        String[] fileNames = new String[]{
            "color_squares.png",
            "android_statues_01.jpg"
        };
        int[] kCells = new int[] {
            9,
            400
        };
        int[] numbersOfIterations = new int[] {
            1,
            1
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
    
    public SData[] getDetailedTrainingData() throws IOException {
        
        String dir = ResourceFinder.findTmpDataDirectory() + 
            "/berkeleySegSubset";
        
        SData[] data = new SData[11];
        
        data[0] = new SData();
        data[0].dirPath = dir + "/101085";
        data[0].imgFileName = "101085.jpg";
        data[0].segFileName = "101085_1124.seg";
            
        data[1] = new SData();
        data[1].dirPath = dir + "/101087";
        data[1].imgFileName = "101087.jpg";
        data[1].segFileName = "101087_1126.seg";
        
        data[2] = new SData();
        data[2].dirPath = dir + "/126007";
        data[2].imgFileName = "126007.jpg";
        data[2].segFileName = "126007_1123.seg";
        
        data[3] = new SData();
        data[3].dirPath = dir + "/167062";
        data[3].imgFileName = "167062.jpg";
        data[3].segFileName = "167062_1112.seg";
        
        data[4] = new SData();
        data[4].dirPath = dir + "/216081";
        data[4].imgFileName = "216081.jpg";
        data[4].segFileName = "216081_1130.seg";
        
        data[5] = new SData();
        data[5].dirPath = dir + "/227092";
        data[5].imgFileName = "227092.jpg";
        data[5].segFileName = "227092_1121.seg";
        
        data[6] = new SData();
        data[6].dirPath = dir + "/229036";
        data[6].imgFileName = "229036.jpg";
        data[6].segFileName = "229036_1123.seg";
      
        data[7] = new SData();
        data[7].dirPath = dir + "/241004";
        data[7].imgFileName = "241004.jpg";
        data[7].segFileName = "241004_1108.seg";
        
        data[8] = new SData();
        data[8].dirPath = dir + "/37073";
        data[8].imgFileName = "37073.jpg";
        data[8].segFileName = "37073_1119.seg";
        
        data[9] = new SData();
        data[9].dirPath = dir + "/42049";
        data[9].imgFileName = "42049.jpg";
        data[9].segFileName = "42049_1123.seg";
        
        data[10] = new SData();
        data[10].dirPath = dir + "/62096";
        data[10].imgFileName = "62096.jpg";
        data[10].segFileName = "62096_1107.seg";
        
        return data;
    }
    
     public SData[] getLessDetailedTrainingData() throws IOException {
        
        String dir = ResourceFinder.findTmpDataDirectory() + 
            "/berkeleySegSubset";
        
        SData[] data = new SData[11];
        
        data[0] = new SData();
        data[0].dirPath = dir + "/101085";
        data[0].imgFileName = "101085.jpg";
        data[0].segFileName = "101085_1108.seg";
            
        data[1] = new SData();
        data[1].dirPath = dir + "/101087";
        data[1].imgFileName = "101087.jpg";
        data[1].segFileName = "101087_1109.seg";
        
        data[2] = new SData();
        data[2].dirPath = dir + "/126007";
        data[2].imgFileName = "126007.jpg";
        data[2].segFileName = "126007_1132.seg";
        
        data[3] = new SData();
        data[3].dirPath = dir + "/167062";
        data[3].imgFileName = "167062.jpg";
        data[3].segFileName = "167062_1123.seg";
        
        data[4] = new SData();
        data[4].dirPath = dir + "/216081";
        data[4].imgFileName = "216081.jpg";
        data[4].segFileName = "216081_1109.seg";
        
        data[5] = new SData();
        data[5].dirPath = dir + "/227092";
        data[5].imgFileName = "227092.jpg";
        data[5].segFileName = "227092_1130.seg";
        
        data[6] = new SData();
        data[6].dirPath = dir + "/229036";
        data[6].imgFileName = "229036.jpg";
        data[6].segFileName = "229036_1109.seg";
      
        data[7] = new SData();
        data[7].dirPath = dir + "/241004";
        data[7].imgFileName = "241004.jpg";
        data[7].segFileName = "241004_1108.seg";
        
        data[8] = new SData();
        data[8].dirPath = dir + "/37073";
        data[8].imgFileName = "37073.jpg";
        data[8].segFileName = "37073_1130.seg";
        
        data[9] = new SData();
        data[9].dirPath = dir + "/42049";
        data[9].imgFileName = "42049.jpg";
        data[9].segFileName = "42049_1115.seg";
        
        data[10] = new SData();
        data[10].dirPath = dir + "/62096";
        data[10].imgFileName = "62096.jpg";
        data[10].segFileName = "62096_1130.seg";
     
        return data;
    }
}
