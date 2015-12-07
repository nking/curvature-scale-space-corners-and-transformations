package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.MultiArrayMergeSort;
import algorithms.misc.Histogram;
import algorithms.util.PairInt;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
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
    
    public void testNextSegmentation() throws Exception {
        
        String fileName1, fileName2;

        for (int i = 5; i < 6; ++i) {
            //fileName1 = "valve_gaussian.png";
            //fileName2 = "valve_gaussian.png";
            switch(i) {
                case 0: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    fileName2 = "brown_lowe_2003_image2.jpg";
                    break;
                }
                case 1: {
                    fileName1 = "venturi_mountain_j6_0001.png";
                    fileName2 = "venturi_mountain_j6_0010.png";
                    break;
                }
                case 2: {
                    fileName1 = "books_illum3_v0_695x555.png";
                    fileName2 = "books_illum3_v6_695x555.png";
                    break;
                }
                case 3: {
                    fileName1 = "campus_010.jpg";
                    fileName2 = "campus_011.jpg";
                    break;
                }
                case 4: {
                    fileName1 = "merton_college_I_001.jpg";
                    fileName2 = "merton_college_I_002.jpg";
                    break;
                }
                default: {
                    fileName1 = "checkerboard_01.jpg";
                    fileName2 = "checkerboard_02.jpg";
                    break;
                }
            }
            
            System.out.println("fileName1=" + fileName1);
          
            String bin = ResourceFinder.findDirectory("bin");
            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
            int idx = fileName1.lastIndexOf(".");
            String fileNameRoot = fileName1.substring(0, idx);

            GreyscaleImage img1 = ImageIOHelper.readImage(filePath1).copyToGreyscale();
            GreyscaleImage img2 = ImageIOHelper.readImage(filePath2).copyToGreyscale();
            
            ImageSegmentation imageSegmentation = new ImageSegmentation();
            GreyscaleImage segImg1 = imageSegmentation.createGreyscale5(img1);
            GreyscaleImage segImg2 = imageSegmentation.createGreyscale5(img2);
            
            String outPath1 = bin + "/seg_1_" + fileNameRoot +".png";
            String outPath2 = bin + "/seg_2_" + fileNameRoot +".png";
            ImageIOHelper.writeOutputImage(outPath1, segImg1);
            ImageIOHelper.writeOutputImage(outPath2, segImg2);
            
            // a look at finding lines with Hough transform
            ZhangSuenLineThinner lineThinner = new ZhangSuenLineThinner();
            lineThinner.applyFilter(segImg1);
            String outPath1_0 = bin + "/seg_1__lt_" + fileNameRoot + ".png";
            ImageIOHelper.writeOutputImage(outPath1_0, segImg1);
            
            algorithms.compGeometry.HoughTransform ht0 = 
                new algorithms.compGeometry.HoughTransform();
            Map<PairInt, Set<PairInt>> outputPolarCoordsPixMap = 
                ht0.calculateLineGivenEdges(segImg1);
            List<PairInt> outSortedKeys = ht0.sortByVotes(outputPolarCoordsPixMap);
            
            Image tmp1SegImg1 = segImg1.copyToColorGreyscale();
            
            for (int ii = 0; ii < outSortedKeys.size(); ++ii) {
                int[] c = ImageIOHelper.getNextRGB(ii);
                PairInt thetaRadius = outSortedKeys.get(ii);
                Set<PairInt> pixCoords = outputPolarCoordsPixMap.get(thetaRadius);
                
                for (PairInt p : pixCoords) {
                    tmp1SegImg1.setRGB(p.getX(), p.getY(), c[0], c[1], c[2]);
                }
            }
            ImageIOHelper.writeOutputImage(bin + "/seg_1_hough1_" + fileNameRoot + ".png", tmp1SegImg1);
            
            Map<PairInt, PairInt> pixToTRMap = ht0.createPixTRMapsFromSorted(
                outSortedKeys, outputPolarCoordsPixMap);
            Map<PairInt, Integer> trColorMap = new HashMap<PairInt, Integer>();
            Image tmp2SegImg1 = segImg1.copyToColorGreyscale();
            for (int col = 0; col < tmp2SegImg1.getWidth(); ++col) {
                for (int row = 0; row < tmp2SegImg1.getHeight(); ++row) {
                    int v = tmp2SegImg1.getR(col, row);
                    if (v > 0) {
                        PairInt p = new PairInt(col, row);
                        PairInt tr = pixToTRMap.get(p);
                        Integer cIndex = trColorMap.get(tr);
                        if (cIndex == null) {
                            cIndex = Integer.valueOf(trColorMap.size());
                            trColorMap.put(tr, cIndex);
                        }
                        int cIdx = cIndex.intValue();
                        int[] c = ImageIOHelper.getNextRGB(cIdx);
                        tmp2SegImg1.setRGB(col, row, c[0], c[1], c[2]);
                    }
                }
            }
            ImageIOHelper.writeOutputImage(bin + "/seg_1_hough2_" + fileNameRoot + ".png", tmp2SegImg1);
            
            int z0 = 1;
        }
    }
    
    public void est0() throws Exception {
        
        String[] fileNames = new String[2];
        
        /*
        tentatively, looks like ne segmentation algor + canny edge corners is 
        the right combination.
        */
                
        //fileNames[0] = "merton_college_I_002.jpg";
        //fileNames[1] = "merton_college_I_001.jpg";
        //fileNames[0] = "brown_lowe_2003_image1.jpg";
        //fileNames[1] = "brown_lowe_2003_image2.jpg";
        //fileNames[0] = "venturi_mountain_j6_0001.png";
        //fileNames[1] = "venturi_mountain_j6_0010.png";
        //fileNames[0] = "books_illum3_v0_695x555.png";
        //fileNames[1] = "books_illum3_v6_695x555.png";
        //fileNames[0] = "campus_010.jpg";
        //fileNames[1] = "campus_011.jpg";
        fileNames[0] = "checkerboard_02.jpg";
        fileNames[1] = "checkerboard_01.jpg";
        
        //fileNames[0] = "seattle.jpg";
        //fileNames[1] = "arches.jpg";
        //fileNames[0] = "stinson_beach.jpg";
        //fileNames[1] = "cloudy_san_jose.jpg";
        //fileNames[0] = "stonehenge.jpg";
        //fileNames[1] = "norwegian_mtn_range.jpg";
        //fileNames[0] = "halfdome.jpg";
        //fileNames[1] = "costa_rica.jpg";
        //fileNames[0] = "new-mexico-sunrise_w725_h490.jpg";
        //fileNames[1] = "arizona-sunrise-1342919937GHz.jpg";
        //fileNames[0] = "sky_with_rainbow.jpg";
        //fileNames[1] = "sky_with_rainbow2.jpg";
        //fileNames[0] = "patagonia_snowy_foreground.jpg";
        //fileNames[1] = "mt_rainier_snowy_field.jpg";  
        //fileNames[0] = "klein_matterhorn_snowy_foreground.jpg";
        //fileNames[1] = "30.jpg";
        //fileNames[0] = "arches_sun_01.jpg";
        //fileNames[1] = "stlouis_arch.jpg";
        //fileNames[0] = "contrail.jpg";
        
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
            
            GreyscaleImage gsImg = imageSegmentation.createGreyscale3(img);
            
            String bin = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(bin + "/" + fileNameRoot  + "_segmentation.png", gsImg);
            
            
            /*
            // canny edges are built into this detector:
            CurvatureScaleSpaceCornerDetector detector = new CurvatureScaleSpaceCornerDetector(
                imageProcessor.binImage(ImageIOHelper.readImageExt(filePath), binFactor));
            detector.findCorners();
            PairIntArray corners = detector.getCorners();
            List<PairIntArray> edges = detector.getEdgesInOriginalReferenceFrame();            
            Set<CornerRegion> cornerRegions = detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);
            
            //ImageExt img3 = gsImg.copyToColorGreyscaleExt();
            
            //ImageIOHelper.addAlternatingColorCurvesToImage(edges, img2);
            //ImageIOHelper.addCornerRegionsToImage(cornerRegions, img2, 0, 0, 1, 255, 0, 0);        
            //ImageIOHelper.writeOutputImage(bin + "/" + fileNameRoot + "_cornerRegions.png", img2);
            //ImageIOHelper.addCurveToImage(corners, img3, 1, 255, 0, 0);            
            //ImageIOHelper.writeOutputImage(bin + "/" + fileNameRoot + "_corners.png", img3);
            */
            
            img = ImageIOHelper.readImageExt(filePath);
            BlobPerimeterHelper imgHelper = new BlobPerimeterHelper(img, fileNameRoot);
            imgHelper.createBinnedGreyscaleImage(binnedImageMaxDimension);
            imgHelper.applySegmentation(SegmentationType.GREYSCALE_HIST, true);
            BlobCornerHelper blobCornerHelper = new BlobCornerHelper(imgHelper);
            List<List<CornerRegion>> cornerRegions2 = 
                blobCornerHelper.generatePerimeterCorners(SegmentationType.GREYSCALE_HIST, true);
           
            ImageExt img4 = gsImg.copyToColorGreyscaleExt();
            List<Set<PairInt>> blobs = imgHelper.getBlobs(SegmentationType.GREYSCALE_HIST, true);
            for (int i = 0; i < blobs.size(); ++i) {
                Set<PairInt> set = blobs.get(i);
                int[] rgb = ImageIOHelper.getNextRGB(i);
                ImageIOHelper.addToImage(set, 0, 0, img4, rgb[0], rgb[1], rgb[2]);
            }
            for (int i = 0; i < cornerRegions2.size(); ++i) {
                List<CornerRegion> crList = cornerRegions2.get(i);
                ImageIOHelper.addCornerRegionsToImage(crList, img4, 0, 0, 1, 255, 255, 255);
                ImageIOHelper.addCornerRegionsToImage(crList, img4, 0, 0, 0, 255, 0, 0);
            }            
            ImageIOHelper.writeOutputImage(bin + "/" + fileNameRoot + "_blob_cr.png", img4);
            
            if (fileNameRoot.contains("checkerboard")) {
                /*StringBuilder sb = new StringBuilder();
                for (List<CornerRegion> list : cornerRegions2) {
                    for (CornerRegion cr : list) {
                        sb.append(String.format("crImg.setValue(%d, %d);\n", 
                            Math.round(cr.getX()[cr.getKMaxIdx()]), 
                            Math.round(cr.getY()[cr.getKMaxIdx()]), 255));
                    }
                }
                log.info(sb.toString());*/
                GreyscaleImage crImg = new GreyscaleImage(img4.getWidth(), img4.getHeight(),
                    GreyscaleImage.Type.Bits32FullRangeInt);
                if (fileNameRoot.equals("checkerboard_01")) {
                    populateCleanedCR1(crImg);
                } else {
                    populateCleanedCR2(crImg);
                }
                imageProcessor.apply2DFFT(crImg, true);
                ImageIOHelper.writeOutputImage(bin + "/" + fileNameRoot + "_blob_cr_cleaned_fft.png", crImg);
                plotFFT(crImg, fileNameRoot + "_cleaned");
                
            }

            GreyscaleImage crImg = new GreyscaleImage(img4.getWidth(), img4.getHeight(),
                GreyscaleImage.Type.Bits32FullRangeInt);
            for (List<CornerRegion> list : cornerRegions2) {
                for (CornerRegion cr : list) {
                    crImg.setValue(cr.getX()[cr.getKMaxIdx()], 
                        cr.getY()[cr.getKMaxIdx()], 255);
                }
            }
            //ImageDisplayer.displayImage("img0", crImg);
            imageProcessor.apply2DFFT(crImg, true);
            //ImageDisplayer.displayImage("FFT of img0", crImg);
        
            ImageIOHelper.writeOutputImage(bin + "/" + fileNameRoot + "_blob_cr_fft.png", crImg);
            
            plotFFT(crImg, fileNameRoot);
            
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

    private void plotFFT(GreyscaleImage crImg, String fileNameRoot) throws IOException {

        int bn = 1;//8
        float[] xPoints = new float[crImg.getWidth() / bn];
        float[] yPoints = new float[crImg.getWidth() / bn];
        float xmn = 0;
        float xmx = crImg.getWidth() / bn;
        float ymn = Float.MAX_VALUE;
        float ymx = Float.MIN_VALUE;
        int row = 50;
        for (int i = 0; i < (crImg.getWidth() / bn) - 1; ++i) {
            int ii = bn * i;
            xPoints[i] = i;
            for (int k = 0; k < bn; ++k) {
                yPoints[i] += crImg.getValue(ii + k, row);
            }
            yPoints[i] /= bn;
            if (yPoints[i] < ymn) {
                ymn = yPoints[i];
            }
            if (yPoints[i] > ymx) {
                ymx = yPoints[i];
            }
        }

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(xmn, xmx, ymn, ymx,
            xPoints, yPoints, null, null, xPoints, yPoints,
            "X fft_" + fileNameRoot);

        xPoints = new float[crImg.getHeight() / bn];
        yPoints = new float[crImg.getHeight() / bn];
        xmn = 0;
        xmx = crImg.getHeight() / bn;
        ymn = Float.MAX_VALUE;
        ymx = Float.MIN_VALUE;
        int col = 50;
        for (int j = 0; j < (crImg.getHeight() / bn) - 1; ++j) {
            int jj = bn * j;
            xPoints[j] = j;
            for (int k = 0; k < bn; ++k) {
                yPoints[j] += crImg.getValue(jj + k, row);
            }
            yPoints[j] /= bn;
            if (yPoints[j] < ymn) {
                ymn = yPoints[j];
            }
            if (yPoints[j] > ymx) {
                ymx = yPoints[j];
            }
        }

        plotter.addPlot(xmn, xmx, ymn, ymx,
            xPoints, yPoints, null, null, xPoints, yPoints,
            "y fft_" + fileNameRoot);

        plotter.writeFile(fileNameRoot + "_blob_cr_fft");

    }

    private void populateCleanedCR2(GreyscaleImage crImg) {
        crImg.setValue(84, 146);
        crImg.setValue(36, 146);
        crImg.setValue(35, 99);
        crImg.setValue(85, 99);
        crImg.setValue(83, 97);
        crImg.setValue(37, 97);
        //crImg.setValue(36, 77); artifact
        //crImg.setValue(35, 65); artifact
        crImg.setValue(36, 50);
        //crImg.setValue(45, 49); artifact
        //crImg.setValue(67, 50); artifact
        crImg.setValue(83, 51);
        crImg.setValue(83, 193);
        crImg.setValue(37, 194);
        crImg.setValue(36, 148);
        crImg.setValue(84, 148);
        crImg.setValue(131, 99);
        crImg.setValue(86, 98);
        crImg.setValue(86, 51);
        //crImg.setValue(105, 52); artifact
        crImg.setValue(131, 53);
        crImg.setValue(131, 192);
        crImg.setValue(86, 193);
        crImg.setValue(86, 146);
        crImg.setValue(131, 146);
        crImg.setValue(130, 52);
        //crImg.setValue(105, 51); artifact
        crImg.setValue(86, 50);
        crImg.setValue(86, 4);
        //crImg.setValue(104, 5); artifact
        //crImg.setValue(119, 6); artifact
        crImg.setValue(131, 7);
        crImg.setValue(176, 54);
        crImg.setValue(133, 53);
        crImg.setValue(132, 13);
        crImg.setValue(134, 6);
        //crImg.setValue(154, 7); artifact
        //crImg.setValue(168, 8); artifact
        crImg.setValue(176, 9);
        crImg.setValue(130, 145);
        crImg.setValue(86, 144);
        crImg.setValue(87, 99);
        //crImg.setValue(92, 100); artifact
        crImg.setValue(131, 101);
        crImg.setValue(176, 144);
        //crImg.setValue(138, 145); artifact
        crImg.setValue(133, 146);
        crImg.setValue(132, 101);
        crImg.setValue(177, 101);
        crImg.setValue(175, 99);
        crImg.setValue(133, 98);
        crImg.setValue(134, 54);
        //crImg.setValue(155, 55); artifact
        crImg.setValue(176, 56);
        //crImg.setValue(219, 89); artifact
        crImg.setValue(218, 101);
        crImg.setValue(178, 100);
        crImg.setValue(177, 56);
        //crImg.setValue(187, 55); artifact
        //crImg.setValue(190, 56); artifact
        crImg.setValue(220, 58);
        crImg.setValue(175, 189);
        //crImg.setValue(168, 190); artifact
        crImg.setValue(134, 191);
        crImg.setValue(133, 148);
        crImg.setValue(176, 147);
        crImg.setValue(218, 189);
        crImg.setValue(178, 190);
        crImg.setValue(177, 146);
        crImg.setValue(219, 146);
        crImg.setValue(218, 56);
        //crImg.setValue(190, 55); artifact
        crImg.setValue(179, 54);
        crImg.setValue(179, 10);
        //crImg.setValue(190, 11); artifact
        //crImg.setValue(203, 12); artifact
        crImg.setValue(219, 14);
        crImg.setValue(218, 144);
        crImg.setValue(178, 144);
        //crImg.setValue(177, 116); artifact
        crImg.setValue(178, 102);
        //crImg.setValue(201, 101); artifact
        //crImg.setValue(204, 102); artifact
        crImg.setValue(219, 103);
        crImg.setValue(259, 144);
        crImg.setValue(221, 145);
        crImg.setValue(221, 101);
        crImg.setValue(260, 103);
        crImg.setValue(260, 58);
        //crImg.setValue(248, 57); artifact
        crImg.setValue(222, 56);
        crImg.setValue(222, 12);
        //crImg.setValue(237, 13); artifact
        //crImg.setValue(250, 14); artifact
        //crImg.setValue(253, 15); artifact
        crImg.setValue(261, 16);
    }

    private void populateCleanedCR1(GreyscaleImage crImg) {
        crImg.setValue(123, 148);
        crImg.setValue(78, 148);
        crImg.setValue(77, 102);
        //crImg.setValue(110, 101); artifact
        //crImg.setValue(113, 102); artifact
        crImg.setValue(124, 103);
        crImg.setValue(123, 55);
        crImg.setValue(78, 54);
        crImg.setValue(78, 8);
        crImg.setValue(124, 9);
        crImg.setValue(169, 101);
        crImg.setValue(125, 102);
        crImg.setValue(125, 55);
        crImg.setValue(138, 56);
        crImg.setValue(170, 57);
        //crImg.setValue(170, 96); artifact
        crImg.setValue(75, 147);
        crImg.setValue(31, 147);
        crImg.setValue(30, 102);
        crImg.setValue(50, 101);
        crImg.setValue(53, 102);
        crImg.setValue(76, 103);
        crImg.setValue(75, 53);
        crImg.setValue(31, 53);
        crImg.setValue(30, 17);
        crImg.setValue(31, 14);
        crImg.setValue(32, 7);
        crImg.setValue(53, 8);
        crImg.setValue(76, 9);
        crImg.setValue(214, 56);
        crImg.setValue(171, 56);
        crImg.setValue(171, 10);
        //crImg.setValue(206, 11); artifact
        crImg.setValue(215, 12);
        crImg.setValue(213, 148);
        crImg.setValue(171, 148);
        crImg.setValue(171, 102);
        crImg.setValue(215, 103);
        //crImg.setValue(214, 112); artifact
        //crImg.setValue(215, 136); artifact
        crImg.setValue(168, 193);
        crImg.setValue(125, 194);
        crImg.setValue(125, 148);
        crImg.setValue(169, 149);
        crImg.setValue(122, 101);
        crImg.setValue(79, 100);
        crImg.setValue(78, 56);
        crImg.setValue(108, 55);
        crImg.setValue(111, 56);
        crImg.setValue(123, 57);
        crImg.setValue(259, 102);
        crImg.setValue(215, 101);
        //crImg.setValue(215, 66); artifact
        crImg.setValue(217, 57);
        crImg.setValue(260, 58);
        crImg.setValue(168, 55);
        crImg.setValue(135, 54);
        crImg.setValue(125, 53);
        crImg.setValue(125, 10);
        //crImg.setValue(141, 9); artifact
        //crImg.setValue(144, 10); artifact
        crImg.setValue(169, 11);
        crImg.setValue(258, 193);
        crImg.setValue(216, 193);
        crImg.setValue(216, 148);
        crImg.setValue(258, 148);
        crImg.setValue(122, 193);
        crImg.setValue(79, 193);
        crImg.setValue(78, 150);
        crImg.setValue(123, 150);
        crImg.setValue(168, 147);
        crImg.setValue(126, 147);
        crImg.setValue(125, 104);
        crImg.setValue(169, 104);
        //crImg.setValue(214, 66); artifact
        crImg.setValue(213, 101);
        crImg.setValue(172, 101);
        crImg.setValue(172, 57);
        crImg.setValue(215, 58);
        crImg.setValue(213, 192);
        crImg.setValue(170, 191);
        crImg.setValue(170, 150);
        crImg.setValue(214, 150);
        crImg.setValue(258, 147);
        crImg.setValue(217, 147);
        crImg.setValue(216, 104);
        crImg.setValue(258, 103);
        crImg.setValue(258, 56);
        crImg.setValue(217, 56);
        crImg.setValue(217, 12);
        crImg.setValue(258, 13);
    }

}
