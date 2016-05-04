package algorithms.imageProcessing.optimization.segmentation;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.optimization.segmentation.SegmentationNN.BData;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SegmentationNNTest extends TestCase {
    
    public SegmentationNNTest() {
    }
    
    public void est0() throws Exception {
        
        String dir = ResourceFinder.findTmpDataDirectory() + 
            "/berkeleySegSubset";
        
        String[][] data = getTrainingData();
        
        for (int i = 0; i < data.length; ++i) {
            String[] fileNames = data[i];
            String subDir = fileNames[0];
            String imgFileName = fileNames[1];
            for (int j = 2; j < fileNames.length; ++j) {
                plot(dir + "/" + subDir, imgFileName, fileNames[j]);
            }
        }
    }
    
    public void testNN() throws Exception {
    }
    
    public void testResults() throws Exception {
        
        boolean enabled = false;
        
        if (!enabled) {
            return;
        }
        
        BData[] data = getDetailedTrainingData();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        BerkeleySegmentationFileReader reader = new BerkeleySegmentationFileReader();
        
        for (int i = 0; i < data.length; ++i) {
            
            String rootName = data[i].imaFileName.split("\\.")[0];
            
            String imgFilePath = data[i].dirPath + "/" + data[i].imaFileName;
        
            String segFilePath = data[i].dirPath + "/" + data[i].segFileName;
                                        
            ImageExt img = ImageIOHelper.readImageExt(imgFilePath);
            
            List<Set<PairInt>> sets0 = 
                imageSegmentation.createColorEdgeSegmentation(img, rootName);
                        
            int[][] expectedCentroidsN = reader.readCentroids(segFilePath);
            
            ImageIOHelper.addAlternatingColorPointSetsToImage(sets0, 0, 0, 0, img);
        
            MiscDebug.writeImage(img, rootName + "_seg_to_train_");
            
            plot(data[i].dirPath, data[i].imaFileName, data[i].segFileName);
        }
    }
    
    public String[][] getTrainingData() {
        
        String[][] data = new String[11][];
        
        // detailed = 101085_1124.seg
        // least detailed = 101085_1108.seg
        data[0] = new String[]{
            "101085",
            "101085.jpg",
            "101085_1105.seg", "101085_1107.seg", "101085_1108.seg",
            "101085_1113.seg", "101085_1124.seg"
        };
        
        // detailed = 101087_1126.seg
        // least detailed = 101087_1109.seg
        data[1] = new String[]{
            "101087",
            "101087.jpg",
            "101087_1105.seg", "101087_1108.seg", "101087_1109.seg", 
            "101087_1126.seg"
        };
          
        // detailed = 126007_1123.seg
        // least detailed = 126007_1132.seg
        data[2] = new String[]{
            "126007",
            "126007.jpg",
            "126007_1103.seg", "126007_1104.seg", "126007_1107.seg", 
            "126007_1108.seg", "126007_1123.seg", "126007_1129.seg",
            "126007_1132.seg"
        };
        
        // detailed = 167062_1112.seg
        // least detailed = 167062_1123.seg
        data[3] = new String[]{
            "167062",
            "167062.jpg",
            "167062_1105.seg", "167062_1109.seg", "167062_1112.seg",
            "167062_1115.seg", "167062_1123.seg"
        };
        
        // detailed = 216081_1130.seg
        // least detailed = 216081_1109.seg
        data[4] = new String[]{
            "216081",
            "216081.jpg",
            "216081_1105.seg", "216081_1109.seg", "216081_1119.seg", 
            "216081_1124.seg", "216081_1130.seg"
        };
        
        // detailed = 227092_1121.seg
        // least detailed = 227092_1130.seg
        data[5] = new String[]{
            "227092",
            "227092.jpg",
            "227092_1108.seg", "227092_1115.seg", "227092_1121.seg",
            "227092_1124.seg", "227092_1130.seg"
        };
        
        // detailed = 229036_1123.seg
        // least detailed = 229036_1109.seg
        data[6] = new String[]{
            "229036",
            "229036.jpg",
            "229036_1104.seg", "229036_1105.seg", "229036_1109.seg", 
            "229036_1115.seg", "229036_1117.seg", "229036_1123.seg",
            "229036_1124.seg"
        };
        
        // detailed = 241004_1109.seg
        // least detailed = 241004_1108.seg
        data[7] = new String[]{
            "241004",
            "241004.jpg",
            "241004_1103.seg", "241004_1108.seg", "241004_1109.seg", 
            "241004_1122.seg", "241004_1124.seg"
        };
        
        // detailed = 37073_1119.seg
        // least detailed = 37073_1130.seg
        data[8] = new String[]{
            "37073",
            "37073.jpg",
            "37073_1109.seg", "37073_1119.seg", "37073_1123.seg", 
            "37073_1130.seg", "37073_1132.seg"
        };
        
        // detailed = 42049_1123.seg
        // least detailed = 42049_1115.seg
        data[9] = new String[]{
            "42049",
            "42049.jpg",
            "42049_1109.seg", "42049_1112.seg", "42049_1115.seg", 
            "42049_1119.seg", "42049_1123.seg"
        };
        
        // detailed = 62096_1107.seg
        // least detailed = 62096_1130.seg
        data[10] = new String[]{
            "62096",
            "62096.jpg",
            "62096_1105.seg", "62096_1109.seg", "62096_1127.seg", 
            "62096_1132.seg", "62096_1103.seg", "62096_1107.seg", 
            "62096_1123.seg", "62096_1130.seg"
        };
        
        return data;
    }
    
    public BData[] getDetailedTrainingData() throws IOException {
        
        String dir = ResourceFinder.findTmpDataDirectory() + 
            "/berkeleySegSubset";
        
        BData[] data = new BData[11];
        
        data[0] = new BData();
        data[0].dirPath = dir + "/101085";
        data[0].imaFileName = "101085.jpg";
        data[0].segFileName = "101085_1124.seg";
            
        data[1] = new BData();
        data[1].dirPath = dir + "/101087";
        data[1].imaFileName = "101087.jpg";
        data[1].segFileName = "101087_1126.seg";
        
        data[2] = new BData();
        data[2].dirPath = dir + "/126007";
        data[2].imaFileName = "126007.jpg";
        data[2].segFileName = "126007_1123.seg";
        
        data[3] = new BData();
        data[3].dirPath = dir + "/167062";
        data[3].imaFileName = "167062.jpg";
        data[3].segFileName = "167062_1112.seg";
        
        data[4] = new BData();
        data[4].dirPath = dir + "/216081";
        data[4].imaFileName = "216081.jpg";
        data[4].segFileName = "216081_1130.seg";
        
        data[5] = new BData();
        data[5].dirPath = dir + "/227092";
        data[5].imaFileName = "227092.jpg";
        data[5].segFileName = "227092_1121.seg";
        
        data[6] = new BData();
        data[6].dirPath = dir + "/229036";
        data[6].imaFileName = "229036.jpg";
        data[6].segFileName = "229036_1123.seg";
      
        data[7] = new BData();
        data[7].dirPath = dir + "/241004";
        data[7].imaFileName = "241004.jpg";
        data[7].segFileName = "241004_1108.seg";
        
        data[8] = new BData();
        data[8].dirPath = dir + "/37073";
        data[8].imaFileName = "37073.jpg";
        data[8].segFileName = "37073_1119.seg";
        
        data[9] = new BData();
        data[9].dirPath = dir + "/42049";
        data[9].imaFileName = "42049.jpg";
        data[9].segFileName = "42049_1123.seg";
        
        data[10] = new BData();
        data[10].dirPath = dir + "/62096";
        data[10].imaFileName = "62096.jpg";
        data[10].segFileName = "62096_1107.seg";
        
        return data;
    }
    
    public BData[] getLessDetailedTrainingData() throws IOException {
        
        String dir = ResourceFinder.findTmpDataDirectory() + 
            "/berkeleySegSubset";
        
        BData[] data = new BData[11];
        
        data[0] = new BData();
        data[0].dirPath = dir + "/101085";
        data[0].imaFileName = "101085.jpg";
        data[0].segFileName = "101085_1108.seg";
            
        data[1] = new BData();
        data[1].dirPath = dir + "/101087";
        data[1].imaFileName = "101087.jpg";
        data[1].segFileName = "101087_1109.seg";
        
        data[2] = new BData();
        data[2].dirPath = dir + "/126007";
        data[2].imaFileName = "126007.jpg";
        data[2].segFileName = "126007_1132.seg";
        
        data[3] = new BData();
        data[3].dirPath = dir + "/167062";
        data[3].imaFileName = "167062.jpg";
        data[3].segFileName = "167062_1123.seg";
        
        data[4] = new BData();
        data[4].dirPath = dir + "/216081";
        data[4].imaFileName = "216081.jpg";
        data[4].segFileName = "216081_1109.seg";
        
        data[5] = new BData();
        data[5].dirPath = dir + "/227092";
        data[5].imaFileName = "227092.jpg";
        data[5].segFileName = "227092_1130.seg";
        
        data[6] = new BData();
        data[6].dirPath = dir + "/229036";
        data[6].imaFileName = "229036.jpg";
        data[6].segFileName = "229036_1109.seg";
      
        data[7] = new BData();
        data[7].dirPath = dir + "/241004";
        data[7].imaFileName = "241004.jpg";
        data[7].segFileName = "241004_1108.seg";
        
        data[8] = new BData();
        data[8].dirPath = dir + "/37073";
        data[8].imaFileName = "37073.jpg";
        data[8].segFileName = "37073_1130.seg";
        
        data[9] = new BData();
        data[9].dirPath = dir + "/42049";
        data[9].imaFileName = "42049.jpg";
        data[9].segFileName = "42049_1115.seg";
        
        data[10] = new BData();
        data[10].dirPath = dir + "/62096";
        data[10].imaFileName = "62096.jpg";
        data[10].segFileName = "62096_1130.seg";
     
        return data;
    }
    
    private void plot(String dir, String imgFileName, String segFileName) throws Exception {
        
        String imgFilePath = dir + "/" + imgFileName;
        
        String segFilePath = dir + "/" + segFileName;
        
        String fileRootName = segFileName.split("\\.")[0];
                                
        ImageExt img = ImageIOHelper.readImageExt(imgFilePath);

        MiscDebug.writeImage(img, fileRootName);
        
        BerkeleySegmentationFileReader reader = new BerkeleySegmentationFileReader();
        
        List<Set<PairInt>> sets = reader.readFile(segFilePath);
        
        assertNotNull(sets);
                
        ImageIOHelper.addAlternatingColorPointSetsToImage(sets, 0, 0, 0, img);
        
        MiscDebug.writeImage(img, fileRootName + "_seg_");
    }
}
