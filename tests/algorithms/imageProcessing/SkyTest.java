package algorithms.imageProcessing;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.Sky.SkyObject;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.awt.Color;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SkyTest extends TestCase {
    
    public SkyTest() {
    }
    
    private String[] fileNames = new String[] {
        "seattle.jpg",
        "stonehenge.jpg",
        "costa_rica.jpg",
        "cloudy_san_jose.jpg",
        "patagonia_snowy_foreground.jpg",
        "mt_rainier_snowy_field.jpg",
        "brown_lowe_2003_image1.jpg",
        "brown_lowe_2003_image2.jpg",
        "venturi_mountain_j6_0001.png",
        "venturi_mountain_j6_0010.png",
        "arches.jpg",
        "stinson_beach.jpg",
        "norwegian_mtn_range.jpg",
        "halfdome.jpg",
        "halfdome2.jpg",
        "halfdome3.jpg",
        "new-mexico-sunrise_w725_h490.jpg",
        "arizona-sunrise-1342919937GHz.jpg",
        "sky_with_rainbow.jpg",
        "sky_with_rainbow2.jpg",
        "klein_matterhorn_snowy_foreground.jpg",
        "arches_sun_01.jpg",
        "stlouis_arch.jpg",
        "contrail.jpg"
    };
    
    private String[] sunFileNames = new String[] {
        "costa_rica.jpg",
        "arizona-sunrise-1342919937GHz.jpg",
        "arches_sun_01.jpg",
        "stlouis_arch.jpg",
    };
    
    private String[] rainbowFileNames = new String[] {
        "sky_with_rainbow.jpg",  // bright
        "sky_with_rainbow2.jpg" // dark
    };
    
    public void estFindSun() throws Exception {
        
        /*
        String filePath = ResourceFinder.findFileInTestResources("tmp2.png");
        ImageExt img0 = ImageIOHelper.readImageExt(filePath);

        float[] hsb = new float[3];
        for (int y = 0; y < 120; y+=10) {
            for (int x = 0; x < 95; x+=10) {
                Color.RGBtoHSB(img0.getR(x, y), img0.getG(x, y), img0.getB(x, y), 
                    hsb);
                System.out.format("(%d,%d) h,s,v=%.2f %.2f %.2f\n", 
                    x, y, hsb[0], hsb[1], hsb[2]);
            }
            for (int x = 520; x < img0.getWidth(); x+=10) {
                Color.RGBtoHSB(img0.getR(x, y), img0.getG(x, y), img0.getB(x, y), 
                    hsb);
                System.out.format("(%d,%d) h,s,v=%.2f %.2f %.2f\n", 
                    x, y, hsb[0], hsb[1], hsb[2]);
            }
        }
        if (true) {
            return;
        }
        */
        
        int maxDimension = 256;//512;

        String fileName1 = "";

        for (int i = 0; i < sunFileNames.length; ++i) {

            fileName1 = sunFileNames[i];
            
            int idx = fileName1.lastIndexOf(".");
            String fileName1Root = fileName1.substring(0, idx);

            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img = ImageIOHelper.readImageExt(filePath1);

            ImageProcessor imageProcessor = new ImageProcessor();

            int w1 = img.getWidth();
            int h1 = img.getHeight();

            int binFactor1 = (int) Math.ceil(Math.max(
                (float) w1 / maxDimension,
                (float) h1 / maxDimension));

            img = imageProcessor.binImage(img, binFactor1);
            //MiscDebug.writeImage(img, "_"  + fileName1Root);
        
            GreyscaleImage[] lma = imageProcessor.createLCHForLUV(img);
            for (int k = 0; k < lma.length; ++k) {
                MiscDebug.writeImage(lma[k], "_"
                    + fileName1Root + "_lma_" + k + "_");
            }
            
            Sky sky = new Sky(img);
            SkyObject obj = sky.findSun();
            
            assertNotNull(obj);
            int[] xyCenter = obj.xyCenter;
            if (fileName1Root.contains("costa")) {
                assertTrue(Math.abs(xyCenter[0] - 22) < 3);
                assertTrue(Math.abs(xyCenter[1] - 85) < 3);
            } else if (fileName1Root.contains("arizona")) {
                assertTrue(Math.abs(xyCenter[0] - 117) < 3);
                assertTrue(Math.abs(xyCenter[1] - 118) < 3);
            } else if (fileName1Root.contains("stlouis")) {
                assertTrue(Math.abs(xyCenter[0] - 136) < 3);
                assertTrue(Math.abs(xyCenter[1] - 34) < 3);
            } else if (fileName1Root.contains("arches_sun")) {
                assertTrue(Math.abs(xyCenter[0] - 155) < 3);
                assertTrue(Math.abs(xyCenter[1] - 80) < 3);
            }
            
            ImageIOHelper.addCurveToImage(obj.points, img, 
                0, 0, 255, 0);
            
            ImageIOHelper.addPointToImage(obj.xyCenter[0], obj.xyCenter[1], 
                img, 0, 255, 0, 0);
            
            MiscDebug.writeImage(img, "_" + fileName1Root + "_SUN_");
        
        }
    }
    
    public void testFindRainbows() throws Exception {
       
        int maxDimension = 256;//512;

        String fileName1 = "";

        for (int i = 0; i < rainbowFileNames.length; ++i) {

            fileName1 = rainbowFileNames[i];
            
            int idx = fileName1.lastIndexOf(".");
            String fileName1Root = fileName1.substring(0, idx);

            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img = ImageIOHelper.readImageExt(filePath1);

            ImageProcessor imageProcessor = new ImageProcessor();

            int w1 = img.getWidth();
            int h1 = img.getHeight();

            int binFactor1 = (int) Math.ceil(Math.max(
                (float) w1 / maxDimension,
                (float) h1 / maxDimension));

            img = imageProcessor.binImage(img, binFactor1);
            MiscDebug.writeImage(img, "_"  + fileName1Root);
        
            GreyscaleImage[] lma = imageProcessor.createLCHForLUV(img);
            for (int k = 0; k < lma.length; ++k) {
                MiscDebug.writeImage(lma[k], "_"
                    + fileName1Root + "_lma_" + k + "_");
            }
            
            Sky sky = new Sky(img);
            List<SkyObject> objs = sky.findRainbows();
                  
            assertNotNull(objs);
            
            for (int k = 0; k < objs.size(); ++k) {
                int[] clr = ImageIOHelper.getNextRGB(k);
                ImageIOHelper.addCurveToImage(objs.get(i).points, img, 
                    0, clr[0], clr[1], clr[2]);
            }
            
            System.out.println("size=" + objs.size());
            
            if (fileName1Root.contains("rainbow2")) {
                //assertTrue(Math.abs(xyCenter[0] - 22) < 3);
                //assertTrue(Math.abs(xyCenter[1] - 85) < 3);
            } else {
            }
            
            MiscDebug.writeImage(img, "_" + fileName1Root 
                + "_RAINBOW_");
        }
    }

    public void est0() throws Exception {

        int maxDimension = 256;//512;

        String fileName1 = "";

        //for (int i = 2; i < 3; ++i) {
        for (int i = 0; i < fileNames.length; ++i) {

            fileName1 = fileNames[i];
            
            int idx = fileName1.lastIndexOf(".");
            String fileName1Root = fileName1.substring(0, idx);

            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img = ImageIOHelper.readImageExt(filePath1);

            ImageProcessor imageProcessor = new ImageProcessor();

            int w1 = img.getWidth();
            int h1 = img.getHeight();

            int binFactor1 = (int) Math.ceil(Math.max(
                (float) w1 / maxDimension,
                (float) h1 / maxDimension));

            img = imageProcessor.binImage(img, binFactor1);
            //MiscDebug.writeImage(img, "_"  + fileName1Root);
                
            Sky sky = new Sky(img);
            //sky._printGsRegions0();
            sky._printGsRegions1();
            //sky._printPtRegions1();
            //sky._printPtRegions1();
            
            //GreyscaleImage[] lma = imageProcessor.createLCHForLUV(img);
            //for (int k = 0; k < lma.length; ++k) {
            //    MiscDebug.writeImage(lma[k], "_" + 
            //        fileName1Root + "_lma_" + k + "_");
            //}
        }
    }
    
}
