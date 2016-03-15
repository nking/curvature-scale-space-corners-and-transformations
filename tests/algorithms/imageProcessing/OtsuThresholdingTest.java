package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;

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
            
            GreyscaleImage gsImg0 = img.copyBlueToGreyscale();
            GreyscaleImage gsImg = gsImg0.copyImage();
        
            OtsuThresholding thrshFinder = new OtsuThresholding();
            
            int thrsh0 = thrshFinder.calculateBinaryThreshold256(gsImg0);
            thrshFinder.applyMultiLevelThreshold256(gsImg, 2);
           
            assertTrue(thrsh0 > 0);
            assertTrue(thrsh0 < 255);
            
            for (int i = 0; i < gsImg0.getNPixels(); ++i) {
                if (gsImg0.getValue(i) > thrsh0) {
                    gsImg0.setValue(i, 255);
                } else {
                    gsImg0.setValue(i, 0);
                }
            }
            
            MiscDebug.writeImage(gsImg0, "oneD_" + fileNameRoot);
            MiscDebug.writeImage(gsImg, "twoD_" + fileNameRoot);
        }
    }
    
}
