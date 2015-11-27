package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class ImageSegmentationTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public ImageSegmentationTest(String testName) {
        super(testName);
    }
    
    public void test0() throws Exception {
        
        String[] fileNames = new String[2];
        
        fileNames[0] = "merton_college_I_002.jpg";
        fileNames[1] = "merton_college_I_001.jpg";
        fileNames[0] = "brown_lowe_2003_image1.jpg";
        fileNames[1] = "brown_lowe_2003_image2.jpg";
        //fileNames[0] = "venturi_mountain_j6_0001.png";
        //fileNames[1] = "venturi_mountain_j6_0010.png";
        //fileNames[0] = "books_illum3_v0_695x555.png";
        //fileNames[1] = "books_illum3_v6_695x555.png";
        //fileNames[0] = "campus_010.jpg";
        //fileNames[1] = "campus_011.jpg";
        //fileNames[0] = "checkerboard_01.jpg";
        //fileNames[1] = "checkerboard_02.jpg";
        
        for (String fileName : fileNames) {
            int idx = fileName.lastIndexOf(".");
            String fileNameRoot = fileName.substring(0, idx);        
            String filePath = ResourceFinder.findFileInTestResources(fileName);
            ImageExt img = ImageIOHelper.readImageExt(filePath);
            
            int binnedImageMaxDimension = 512;
            int binFactor = 
                (int) Math.ceil(Math.max((float) img.getWidth() / binnedImageMaxDimension, 
                (float) img.getHeight() / binnedImageMaxDimension));
            ImageProcessor imageProcessor = new ImageProcessor();
            img = imageProcessor.binImage(img, binFactor);
        
            ImageSegmentation imageSegmentation = new ImageSegmentation();
            
            //GreyscaleImage gsImg = imageSegmentation.createGreyscaleWithBWMask(img);
            GreyscaleImage gsImg = imageSegmentation.createGreyscale3(img);
            
            String bin = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(bin + "/" + fileNameRoot  + "_segmentation.png", gsImg);
            
            // canny edges are built into this detector:
            CurvatureScaleSpaceCornerDetector detector = new CurvatureScaleSpaceCornerDetector(gsImg);
            detector.findCorners();
            PairIntArray corners = detector.getCorners();
            List<PairIntArray> edges = detector.getEdgesInOriginalReferenceFrame();            
            Set<CornerRegion> cornerRegions = detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);
            
            ImageExt img2 = gsImg.copyToColorGreyscaleExt();
            ImageExt img3 = gsImg.copyToColorGreyscaleExt();
            ImageIOHelper.addAlternatingColorCurvesToImage(edges, img2);
            ImageIOHelper.addCornerRegionsToImage(cornerRegions, img2, 0, 0, 1, 255, 0, 0);        
            ImageIOHelper.writeOutputImage(bin + "/" + fileNameRoot + "_cornerRegions.png", img2);
            ImageIOHelper.addCurveToImage(corners, img3, 1, 255, 0, 0);            
            ImageIOHelper.writeOutputImage(bin + "/" + fileNameRoot + "_corners.png", img3);
            
            img = ImageIOHelper.readImageExt(filePath);
            BlobPerimeterHelper imgHelper = new BlobPerimeterHelper(img, fileNameRoot);
            imgHelper.createBinnedGreyscaleImage(binnedImageMaxDimension);
            imgHelper.applySegmentation(SegmentationType.GREYSCALE_HIST, true);
            Image img4 = imgHelper.getBinnedSegmentationImage(SegmentationType.GREYSCALE_HIST).copyToColorGreyscale();
            BlobCornerHelper blobCornerHelper = new BlobCornerHelper(imgHelper);
            List<List<CornerRegion>> cornerRegions2 = 
                blobCornerHelper.generatePerimeterCorners(SegmentationType.GREYSCALE_HIST, true);
           
            for (int i = 0; i < cornerRegions2.size(); ++i) {
                List<CornerRegion> crList = cornerRegions2.get(i);
                int[] rgb = ImageIOHelper.getNextRGB(i);
                ImageIOHelper.addCornerRegionsToImage(
                    crList, img4, 0, 0, 1, rgb[0], rgb[1], rgb[2]);
            }
            ImageIOHelper.writeOutputImage(bin + "/" + fileNameRoot + "_cornerregions2.png", img4);
            int z = 1;
        }
    }
    
    public void estApplyColorSegmentation() throws Exception {
        
        //String fileName = "two_circles_color2.png";
        //String fileName = "two_circles_color.png";
        //String fileName = "cloudy_san_jose.jpg";
        //String fileName = "middlebury_cones_im2.png"; // a limitFrac of 0.1 works well
        String fileName = "brown_lowe_2003_image2.jpg";
        //String fileName = "venturi_mountain_j6_0010.png";
        //String fileName = "books_illum3_v6_695x555.png";
        //String fileName = "brown_lowe_2003_image1.jpg";
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        //HistogramEqualizationForColor hEq = new HistogramEqualizationForColor(img);
        //hEq.applyFilter();
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        List<Set<PairInt>> clusterSets = 
            //imageSegmentation.calculateUsingPolarCIEXYAndFrequency(img, /*0.1f,*/ true);
            imageSegmentation.calculateUsingCIEXYAndClustering(img, true);
        
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
    
     public void estApplyColorSegmentation2() throws Exception {
        
         String[] fileNames1 = new String[]{
             "brown_lowe_2003_image1.jpg",
             "venturi_mountain_j6_0001.png",
             "books_illum3_v0_695x555.png"
         };
         String[] fileNames2 = new String[]{
             "brown_lowe_2003_image2.jpg",
             "venturi_mountain_j6_0010.png",
             "books_illum3_v6_695x555.png"
         };
         
         for (int i = 0; i < 1/*fileNames1.length*/; ++i) {
             String fileName1 = fileNames1[i];
             String fileName2 = fileNames2[i];
             
             int idx = fileName1.lastIndexOf(".");
             String fileNameRoot = fileName1.substring(0, idx);
             
             log.info("fileName=" + fileNameRoot);
             
             String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
             String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
             
             //HistogramEqualizationForColor hEq = new HistogramEqualizationForColor(img);
             //hEq.applyFilter();
             
             ImageSegmentation imageSegmentation = new ImageSegmentation();
             String bin = ResourceFinder.findDirectory("bin");
             
             for (int j = 0; j < 20; ++j) {
                 ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
                 ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
                 log.info("   j=" + j);
                 switch (j) {
                     case 0: {
                         int kBands = 2;
                         GreyscaleImage gsImg1 = img1.copyToGreyscale();
                         GreyscaleImage gsImg2 = img2.copyToGreyscale();
                         imageSegmentation.applyUsingKMPP(gsImg1, kBands);
                         imageSegmentation.applyUsingKMPP(gsImg2, kBands);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg2.copyToColorGreyscaleExt();
                         break;
                     }
                     case 1: {
                         int kBands = 3;
                         GreyscaleImage gsImg1 = img1.copyToGreyscale();
                         GreyscaleImage gsImg2 = img2.copyToGreyscale();
                         
                         /*float[] values1 = gsImg1.getFloatValues();
                         float[] errors1 = Errors.populateYErrorsBySqrt(values1);
                         HistogramHolder hist = Histogram.createSimpleHistogram(values1, errors1);
                         String str = hist.plotLogHistogram(fileNameRoot, "img_" + fileNameRoot + "_hist");
                         float maxValue = MiscMath.findMax(values1);*/
                         imageSegmentation.applyUsingKMPP(gsImg1, kBands);
                         imageSegmentation.applyUsingKMPP(gsImg2, kBands);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg2.copyToColorGreyscaleExt();
                         break;
                     }
                     case 2: {
                         int kBands = 8;
                         GreyscaleImage gsImg1 = img1.copyToGreyscale();
                         GreyscaleImage gsImg2 = img2.copyToGreyscale();
                         imageSegmentation.applyUsingKMPP(gsImg1, kBands);
                         imageSegmentation.applyUsingKMPP(gsImg2, kBands);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg2.copyToColorGreyscaleExt();
                         break;
                     }
                     case 3: {
                         int kBands = 2;
                         GreyscaleImage gsImg1 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenHistEq(img1, kBands, true);
                         GreyscaleImage gsImg2 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenHistEq(img2, kBands, true);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg2.copyToColorGreyscaleExt();
                         break;
                     }
                     case 4: {
                         int kBands = 3;
                         GreyscaleImage gsImg1 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenHistEq(img1, kBands, true);
                         GreyscaleImage gsImg2 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenHistEq(img2, kBands, true);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg2.copyToColorGreyscaleExt();
                         break;
                     }
                     case 5: {
                         int kBands = 8;
                         GreyscaleImage gsImg1 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenHistEq(img1, kBands, true);
                         GreyscaleImage gsImg2 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenHistEq(img2, kBands, true);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg2.copyToColorGreyscaleExt();
                         break;
                     }
                     case 6: {
                         int kBands = 2;
                         GreyscaleImage gsImg1 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenKMPPThenHistEq(img1, kBands, true);
                         GreyscaleImage gsImg2 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenKMPPThenHistEq(img2, kBands, true);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg2.copyToColorGreyscaleExt();
                         break;
                     }
                     case 7: {
                         int kBands = 3;
                         GreyscaleImage gsImg1 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenKMPPThenHistEq(img1, kBands, true);
                         GreyscaleImage gsImg2 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenKMPPThenHistEq(img2, kBands, true);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg2.copyToColorGreyscaleExt();
                         break;
                     }
                     case 8: {
                         int kBands = 8;
                         GreyscaleImage gsImg1 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenKMPPThenHistEq(img1, kBands, true);
                         GreyscaleImage gsImg2 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenKMPPThenHistEq(img2, kBands, true);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg2.copyToColorGreyscaleExt();
                         break;
                     }
                     case 9: {
                         int kBands = 2;
                         //TODO: may need a blur of 1
                         GreyscaleImage gsImg1 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenHistogram(img1, kBands);
                         GreyscaleImage gsImg2 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenHistogram(img2, kBands);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg2.copyToColorGreyscaleExt();
                         break;
                     }
                     case 10: {
                         int kBands = 3;
                         //TODO: may need a blur of 1
                         GreyscaleImage gsImg1 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenHistogram(img1, kBands);
                         GreyscaleImage gsImg2 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenHistogram(img2, kBands);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg2.copyToColorGreyscaleExt();
                         break;
                     }
                     case 11: {
                         int kBands = 8;
                         //TODO: may need a blur of 1
                         GreyscaleImage gsImg1 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenHistogram(img1, kBands);
                         GreyscaleImage gsImg2 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenHistogram(img2, kBands);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg2.copyToColorGreyscaleExt();
                         break;
                     }
                     case 12: {
                         int kBands = 2;
                         //TODO: may need a blur of 1
                         GreyscaleImage gsImg1 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenKMPP(img1, kBands);
                         GreyscaleImage gsImg2 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenKMPP(img2, kBands);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg2.copyToColorGreyscaleExt();
                         break;
                     }
                     case 13: {
                         int kBands = 3;
                         //TODO: may need a blur of 1
                         GreyscaleImage gsImg1 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenKMPP(img1, kBands);
                         GreyscaleImage gsImg2 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenKMPP(img2, kBands);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg2.copyToColorGreyscaleExt();
                         break;
                     }
                     case 14: {
                         int kBands = 8;
                         //TODO: may need a blur of 1
                         GreyscaleImage gsImg1 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenKMPP(img1, kBands);
                         GreyscaleImage gsImg2 = imageSegmentation
                             .applyUsingCIEXYPolarThetaThenKMPP(img2, kBands);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg2.copyToColorGreyscaleExt();
                         break;
                     }
                     case 15: {
                         List<Set<PairInt>> groups1 = 
                             imageSegmentation.calculateUsingCIEXYAndClustering(img1, true);
                         List<Set<PairInt>> groups2 = 
                             imageSegmentation.calculateUsingCIEXYAndClustering(img2, true);
                         int nExtraForDot = 0;
                         
                         ImageExt img11 = new ImageExt(img1.getWidth(), img1.getHeight());
                         for (int k = 0; k < groups1.size(); ++k) {
                             int[] rgb = ImageIOHelper.getNextRGB(k);
                             Set<PairInt> set = groups1.get(k);
                             ImageIOHelper.addToImage(set, 0, 0, img11, nExtraForDot, rgb[0],
                                 rgb[1], rgb[2]);
                         }
                         img1 = img11;
                         
                         ImageExt img22 = new ImageExt(img2.getWidth(), img2.getHeight());
                         for (int k = 0; k < groups2.size(); ++k) {
                             int[] rgb = ImageIOHelper.getNextRGB(k);
                             Set<PairInt> set = groups2.get(k);
                             ImageIOHelper.addToImage(set, 0, 0, img22, nExtraForDot, rgb[0],
                                 rgb[1], rgb[2]);
                         }
                         img2 = img22;
                         break;
                     }
                     case 16: {
                         GreyscaleImage gsImg1 = 
                             imageSegmentation.applyUsingPolarCIEXYAndFrequency(img1, true);
                         GreyscaleImage gsImg2 = 
                             imageSegmentation.applyUsingPolarCIEXYAndFrequency(img2, true);
                         img1 = gsImg1.copyToColorGreyscaleExt();
                         img2 = gsImg1.copyToColorGreyscaleExt();
                         break;
                     }
                     case 17: {
                         List<Set<PairInt>> groups1 = 
                             imageSegmentation.calculateUsingPolarCIEXYAndClustering(img1, true);
                         List<Set<PairInt>> groups2 = 
                             imageSegmentation.calculateUsingPolarCIEXYAndClustering(img2, true);
                         int nExtraForDot = 0;
                         
                         ImageExt img11 = new ImageExt(img1.getWidth(), img1.getHeight());
                         for (int k = 0; k < groups1.size(); ++k) {
                             int[] rgb = ImageIOHelper.getNextRGB(k);
                             Set<PairInt> set = groups1.get(k);
                             ImageIOHelper.addToImage(set, 0, 0, img11, nExtraForDot, rgb[0],
                                 rgb[1], rgb[2]);
                         }
                         img1 = img11;
                         
                         ImageExt img22 = new ImageExt(img2.getWidth(), img2.getHeight());
                         for (int k = 0; k < groups2.size(); ++k) {
                             int[] rgb = ImageIOHelper.getNextRGB(k);
                             Set<PairInt> set = groups2.get(k);
                             ImageIOHelper.addToImage(set, 0, 0, img22, nExtraForDot, rgb[0],
                                 rgb[1], rgb[2]);
                         }
                         img2 = img22;
                         break;
                     }
                     case 18: {
                         List<Set<PairInt>> groups1 = 
                             imageSegmentation.calculateUsingPolarCIEXYAndFrequency(img1, true);
                         List<Set<PairInt>> groups2 = 
                             imageSegmentation.calculateUsingPolarCIEXYAndFrequency(img2, true);
                         int nExtraForDot = 0;
                         
                         ImageExt img11 = new ImageExt(img1.getWidth(), img1.getHeight());
                         for (int k = 0; k < groups1.size(); ++k) {
                             int[] rgb = ImageIOHelper.getNextRGB(k);
                             Set<PairInt> set = groups1.get(k);
                             ImageIOHelper.addToImage(set, 0, 0, img11, nExtraForDot, rgb[0],
                                 rgb[1], rgb[2]);
                         }
                         img1 = img11;
                         
                         ImageExt img22 = new ImageExt(img2.getWidth(), img2.getHeight());
                         for (int k = 0; k < groups2.size(); ++k) {
                             int[] rgb = ImageIOHelper.getNextRGB(k);
                             Set<PairInt> set = groups2.get(k);
                             ImageIOHelper.addToImage(set, 0, 0, img22, nExtraForDot, rgb[0],
                                 rgb[1], rgb[2]);
                         }
                         img2 = img22;
                         break;
                     }
                     case 19: {
                         List<Set<PairInt>> groups1 = 
                             imageSegmentation.calculateUsingPolarCIEXYAndFrequency(img1, 0.2f, true);
                         List<Set<PairInt>> groups2 = 
                             imageSegmentation.calculateUsingPolarCIEXYAndFrequency(img2, 0.2f, true);
                         int nExtraForDot = 0;
                         
                         ImageExt img11 = new ImageExt(img1.getWidth(), img1.getHeight());
                         for (int k = 0; k < groups1.size(); ++k) {
                             int[] rgb = ImageIOHelper.getNextRGB(k);
                             Set<PairInt> set = groups1.get(k);
                             ImageIOHelper.addToImage(set, 0, 0, img11, nExtraForDot, rgb[0],
                                 rgb[1], rgb[2]);
                         }
                         img1 = img11;
                         
                         ImageExt img22 = new ImageExt(img2.getWidth(), img2.getHeight());
                         for (int k = 0; k < groups2.size(); ++k) {
                             int[] rgb = ImageIOHelper.getNextRGB(k);
                             Set<PairInt> set = groups2.get(k);
                             ImageIOHelper.addToImage(set, 0, 0, img22, nExtraForDot, rgb[0],
                                 rgb[1], rgb[2]);
                         }
                         img2 = img22;
                         break;
                     }
                     default: {
                         break;
                     }
                 }
                 
                 ImageIOHelper.writeOutputImage(bin + "/" + fileNameRoot + "_segmentation1_" + j + ".png", img1); 
                 ImageIOHelper.writeOutputImage(bin + "/" + fileNameRoot + "_segmentation2_" + j + ".png", img2); 
             }
             
         }                        
               
    }
    
    private int count(List<Set<PairInt>> setList) {
        
        int c = 0;
        for (Set<PairInt> set : setList) {
            c += set.size();
        }
        
        return c;
    }
    
}
