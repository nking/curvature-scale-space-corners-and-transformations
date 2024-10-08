package algorithms.imageProcessing;

import algorithms.misc.Histogram;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertEquals;

public class OtsuThresholdingTest extends TestCase {
    
    public void testCalculateBinaryThreshold256() throws Exception {
        
        String[] fileNames = new String[]{
            "android_statues_01.jpg",
            /*"android_statues_02.jpg",
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
            "checkerboard_02.jpg"*/
        };
        
        for (String fileName : fileNames) {
       
            int idx = fileName.lastIndexOf(".");
            String fileNameRoot = fileName.substring(0, idx);

            String filePath = ResourceFinder.findFileInTestResources(fileName);
            ImageExt img = ImageIOHelper.readImageExt(filePath);
            
            GreyscaleImage gsImg0 = img.copyToGreyscale();
            GreyscaleImage gsImg = gsImg0.copyImage();
        
            OtsuThresholding thrshFinder = new OtsuThresholding();
            
            int nLevels = 3;
            
            thrshFinder.applyMultiLevelThreshold256(gsImg, nLevels);
            
            MiscDebug.writeImage(gsImg, "_adaptive_n3_" + fileNameRoot);
        }
    }
    
    public void testCalculateBinaryThreshold2D() throws Exception {
        
        String[] fileNames = new String[]{
            "susan-in_plus.png",
            "house.gif",
            "android_statues_01.jpg",
            
            "seattle.jpg",
            "stonehenge.jpg",
            "cloudy_san_jose.jpg",
            "patagonia_snowy_foreground.jpg",
            "mt_rainier_snowy_field.jpg",
            "brown_lowe_2003_image1.jpg",
            "venturi_mountain_j6_0010.png",
            "campus_010.jpg",
            "merton_college_I_001.jpg",
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
            "klein_matterhorn_snowy_foreground.jpg",
            "30.jpg",
            "arches_sun_01.jpg",
            "stlouis_arch.jpg",
            "contrail.jpg",
            "checkerboard_01.jpg"
        };
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        for (String fileName : fileNames) {
       
            int idx = fileName.lastIndexOf(".");
            String fileNameRoot = fileName.substring(0, idx);

            String filePath = ResourceFinder.findFileInTestResources(fileName);
            GreyscaleImage img = ImageIOHelper.readImageAsGreyscaleFullRange(
                filePath);
            
            int w = img.getWidth();
            int h = img.getHeight();
            double[][] imgD0 = new double[w][];
            double[][] imgD1 = new double[w][];
            for (int i = 0; i < w; ++i) {
                imgD0[i] = new double[h];
                imgD1[i] = new double[h];
                for (int j = 0; j < h; ++j) {
                    imgD0[i][j] = img.getValue(i, j);
                    imgD1[i][j] = img.getValue(i, j)/255;
                }
            }
        
            OtsuThresholding thrshFinder = new OtsuThresholding();
            
            double thrsh2D256 = thrshFinder.calculateBinaryThreshold2D(imgD0,
                256);
           
            double thrsh1D256 = thrshFinder.calculateBinaryThreshold256(img);
            
            assertTrue(thrsh2D256 > 0);
            assertTrue(thrsh2D256 < 255);
            
            assertTrue(thrsh1D256 > 0);
            assertTrue(thrsh1D256 < 255);
            
            double[][] imgDCp = imageProcessor.copy(imgD0);
            GreyscaleImage img1 = img.copyImage();
            
            for (int i = 0; i < img.getNPixels(); ++i) {
                if (img.getValue(i) > thrsh2D256) {
                    img.setValue(i, 255);
                } else {
                    img.setValue(i, 0);
                }
                if (img.getValue(i) > thrsh1D256) {
                    img1.setValue(i, 255);
                } else {
                    img1.setValue(i, 0);
                }
            }
            
            MiscDebug.writeImage(img, "_twoD_" + fileNameRoot);
            
            MiscDebug.writeImage(img1, "_oneD_" + fileNameRoot);
        
            AdaptiveThresholding at = new AdaptiveThresholding();
            at.applyAdaptiveThresholdImage(imgDCp, 15, 0.2, 255.);
            
            MiscDebug.writeImage(imgDCp, "img_adap_thrsh_" + fileNameRoot
                + ".png");
            
            System.out.format("thresh 2D 256=%.3f, thresh 1D 256 = %.3f\n", 
                    (float)thrsh2D256, (float)thrsh1D256);
        }
    }
    
    public void testHistogramInt() throws Exception {
        
        OtsuThresholding thrshFinder = new OtsuThresholding();
        
        GreyscaleImage img = new GreyscaleImage(16, 16);
        for (int i = 0; i < img.getNPixels(); ++i) {
            img.setValue(i, i);
        }
        int[] h = thrshFinder.createHistogram(img, 0, 255, 256);        
        for (int i = 0; i < 256; ++i) {
            assertEquals(1, h[i]);
        }
        
        img = new GreyscaleImage(16, 16);
        img.fill(255);
        h = thrshFinder.createHistogram(img, 0, 255, 256);
        for (int i = 0; i < 256; ++i) {
            if (i == 255) {
                assertEquals(256, h[i]);
            } else {
                assertEquals(0, h[i]);
            }
        }
        
        img = new GreyscaleImage(16, 16);
        img.fill(0);
        h = thrshFinder.createHistogram(img, 0, 255, 256);
        for (int i = 0; i < 256; ++i) {
            if (i == 0) {
                assertEquals(256, h[i]);
            } else {
                assertEquals(0, h[i]);
            }
        }
        
        // set 2 pixel values up to 256, e.g. {0,0, 2,2, 4,4, ...} 
        img = new GreyscaleImage(16, 16);
        for (int i = 0; i < img.getNPixels(); i += 2) {
            img.setValue(i, i);
            img.setValue(i + 1, i);
        }
        h = thrshFinder.createHistogram(img, 0, 255, 256);
        for (int i = 0; i < img.getNPixels(); i += 2) {
            assertEquals(2, h[i]);
        }
        for (int i = 1; i < img.getNPixels(); i += 2) {
            assertEquals(0, h[i]);
        }
    }
}
