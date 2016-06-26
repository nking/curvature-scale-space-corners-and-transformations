package algorithms.imageProcessing.optimization.segmentation;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.optimization.segmentation.DownhillSimplexSolver.SFit;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.Arrays;
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
    
    public void estBisectionSolver() throws Exception {
        
        boolean enabled = true;
        
        if (!enabled) {
            return;
        }
               
        SData[] data = getDetailedTrainingData();
        
        SequentialBisectorSolver solver = new SequentialBisectorSolver(true, 
            true, data);
        
        double diff = solver.solve();
        
        System.out.println(
            "hsv, reduceNoise=true, difference=" + diff +
            "  params=" + Arrays.toString(solver.getParameters()));
        
        // --------
        solver = new SequentialBisectorSolver(false, 
            true, data);
        
        diff = solver.solve();
        
        System.out.println(
            "cie, reduceNoise=true, difference=" + diff +
            "  parmaeters=" + Arrays.toString(solver.getParameters()));
        
        /*
        normalize = false:
            hsv, reduceNoise=true, difference=33097.4442013048  params=[6.0, 0.05, 1.5, 0.005]

            cie, reduceNoise=true, difference=2397540.564264881  parmaeters=[1.0, 2.5, 1.5, 0.005]
        */
    }
    
    public void estNelderMeadeSolver() throws Exception {
        
        boolean enabled = true;
        
        if (!enabled) {
            return;
        }
               
        SData[] data = getDetailedTrainingData();
        
        boolean useLowNoiseEdges = false;
        
        DownhillSimplexSolver solver = 
            new DownhillSimplexSolver(true, useLowNoiseEdges, data);
        
        SFit fit = solver.solve();
        double diff = fit.costF;
        
        System.out.println(
            "hsv, reduceNoise=" + useLowNoiseEdges + 
            " difference=" + diff +
            "  params=" + Arrays.toString(solver.getParameters()));
        
        // --------
        /*
        solver = new DownhillSimplexSolver(false, 
            useLowNoiseEdges, data);
        
        fit = solver.solve();
        diff = fit.costF;
        
        System.out.println(
            "cie, reduceNoise=" + useLowNoiseEdges
            + " difference=" + diff +
            "  parameters=" + Arrays.toString(solver.getParameters()));
        */
        /*
        normalize = false:
           hsv, nIter=3 bestFit=cost=2.7239726540785525E7 
              tLen=29 tColor=0.48000000019868216 
              tSmallMerge=0.03366666678339243 tR=1.5 
              reduceNoiseF=true colorSpace=1.0
            
           cie nIter=2 bestFit=cost=1.947937187996745E7 
             tLen=11 tColor=4.0 tSmallMerge=0.03500000014901161 
             tR=1.5 reduceNoiseF=true colorSpace=0.0
            
        */
    }
    
    public void testBenchmark() throws Exception {
        
        // simple test that fMeasure is larger
        //   for a better fit
        
        int dMax = 2;
        
        SData[] data = getDetailedTrainingData();
        
        String imgNumber = "101087";
        
        boolean useLowNoiseEdges = false;
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        BerkeleySegmentationFileReader reader = new BerkeleySegmentationFileReader();
        
        boolean reduceNoise = false;
                
        double[] fMeasures = new double[2];
        
        for (int tst = 0; tst < 2; ++tst) {
            int tLen;
            double tClr;
            double tR;  
            double tSMerge;
            if (tst == 0) {
                // expecting better results for these params.
                // which is a higher fMeasure
                tLen = 11;
                tClr = 0.135;
                tR = 2.25;
                tSMerge = 0.025;
            } else {
                tLen = 40;
                tClr = 0.35;
                tR = 0.5;
                tSMerge = 0.125;
            }
            
            for (int i = 0; i < data.length; ++i) {

                String rootName = data[i].imgFileName.split("\\.")[0];

                if (!rootName.contains(imgNumber)) {
                    continue;
                }

                String imgFilePath = data[i].dirPath + "/" + data[i].imgFileName;

                String segFilePath = data[i].dirPath + "/" + data[i].segFileName;

                ImageExt img = ImageIOHelper.readImageExt(imgFilePath);

                List<Set<PairInt>> modelSet = reader.readFile(segFilePath);

                List<PairIntArray> edges = imageSegmentation.extractEdges(img, 
                    reduceNoise, rootName);
                
                List<Set<PairInt>> dataSet = 
                    imageSegmentation.createColorEdgeSegmentation(img, 
                        edges,
                        1, (int)Math.round(tLen), tClr, tR, 
                        reduceNoise, tSMerge, rootName);
            
                // reduces the sets to perimeters for comparisons              
                SegmentationResults dataResults 
                    = new SegmentationResults(dataSet);
            
                SegmentationResults modelResults 
                    = new SegmentationResults(modelSet);
            
 //TODO: can see clearly that need to improve
 // or replace the perimeter algorithm
                
                MiscDebug.writeAlternatingColor(
                    img.createWithDimensions(), 
                    modelResults.getPerimeters(), "_model_" + tst + "_" + rootName);
                
                MiscDebug.writeAlternatingColor(
                    img.createWithDimensions(), 
                    dataResults.getPerimeters(), "_test_" + tst + "_" + rootName);
                
                double fMeasure = dataResults.evaluate(
                    modelResults, dMax);
                
                fMeasures[tst] = fMeasure;
                
                System.out.println("test=" + tst + " fMeasure=" + fMeasure);
                
            }
        }
        
    }
    
    public void estResults() throws Exception {
        
        boolean enabled = false;
        
        if (!enabled) {
            return;
        }
        
        SData[] data = getDetailedTrainingData();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        BerkeleySegmentationFileReader reader = new BerkeleySegmentationFileReader();
        
        for (int i = 0; i < data.length; ++i) {
            
            String rootName = data[i].imgFileName.split("\\.")[0];
            
            String imgFilePath = data[i].dirPath + "/" + data[i].imgFileName;
        
            String segFilePath = data[i].dirPath + "/" + data[i].segFileName;
                                        
            ImageExt img = ImageIOHelper.readImageExt(imgFilePath);
            
            List<Set<PairInt>> sets0 = 
                imageSegmentation.createColorEdgeSegmentation(img, rootName);
                        
            int[][] expectedCentroidsN = reader.readCentroids(segFilePath);
            
            ImageIOHelper.addAlternatingColorPointSetsToImage(sets0, 0, 0, 0, img);
        
            MiscDebug.writeImage(img, rootName + "_seg_to_train_");
            
            plot(data[i].dirPath, data[i].imgFileName, data[i].segFileName);
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
