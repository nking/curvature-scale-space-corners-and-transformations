package algorithms.imageProcessing.features.mser;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.misc.MiscDebug;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class MSERTest extends TestCase {
    
    public void test0() throws IOException {
        
        String[] files = new String[] {
            "android_statues_01_sz1.jpg",
            "android_statues_02.jpg",
            "android_statues_03_sz3.jpg",
            "android_statues_04_sz1.jpg"
        };
            
        for (String file : files) {
        
            File fl = ResourceFinder.findFileInTmpData(file);
        
            Image img = 
                ImageIOHelper.readImage(fl.getPath());
            
            MSER mser = new MSER();
           
            List<List<Region>> regions = mser.findRegions(img.copyToGreyscale2());
            
            Image img1 = ImageIOHelper.readImage(fl.getPath());
            
            Image img2 = img1.copyImage();
            
            int nExtraDot = 0;
                        
            for (int i = 0; i < regions.get(0).size(); ++i) {
                regions.get(0).get(i).drawEllipse(img1, nExtraDot, 
                    255, 255, 255);
            }

            for (int i = 0; i < regions.get(1).size(); ++i) {
                regions.get(1).get(i).drawEllipse(img2, nExtraDot, 
                    255, 255, 255);
            }
            
            MiscDebug.writeImage(img1, file + "_out_0_");
            
            MiscDebug.writeImage(img2, file + "_out_1_");
        
            System.out.println("num extracted regions=" +
                (regions.get(0).size() + regions.get(1).size()));
        
            //---------- looking at centroids ----
            nExtraDot = 1;
            Canonicalizer cr = new Canonicalizer();

            GreyscaleImage gs1 = img1.copyToGreyscale2();
            int width = gs1.getWidth();
            int height = gs1.getHeight();
            int[] data = new int[width * height];
            for (int i = 0; i < gs1.getNPixels(); ++i) {
                data[i] = gs1.getValue(i);
            }
            
            Image img_1 = img2.copyImage();
            
            PairIntArray xyCens = cr.extractRegionXYCenters(regions.get(1), 
                data, img2.getWidth(), img2.getHeight());

            ImageIOHelper.addCurveToImage(xyCens, img_1, 
                nExtraDot, 255, 0, 0);

            MiscDebug.writeImage(img_1, file + "_out_xy1_");


            Image img_2 = img2.copyImage();

            xyCens = cr.calculateIntensityCentroids(regions.get(1), 
                data, img2.getWidth(), img2.getHeight());

            ImageIOHelper.addCurveToImage(xyCens, img_2, 
                nExtraDot, 255, 0, 0);

            MiscDebug.writeImage(img_2, file + "_out_xy2_");

        }
    }
    
    public void est1() throws IOException {
        
        int width = 256;
        int height = 171;
        
        String file = "android_statues_02.jpg";
                    
        File fl = ResourceFinder.findFileInTmpData(file);
        
        int[] data = new int[width * height];
        File flTxt = ResourceFinder.findFileInTmpData("andr02.txt");
        BufferedReader in = null;
        int count = 0;
        try {
            in = new BufferedReader(new FileReader(flTxt));
            String line = in.readLine();
            while (line != null) {
                line = line.trim();
                data[count] = Integer.valueOf(line).intValue();
                count++;
                line = in.readLine();
            }
        } finally {
            if (in != null) {
                in.close();
            }
        }
        
        MSER mser = new MSER();
        
        List<List<Region>> regions = mser.findRegions(data, width, height);
            
        Image img2 = ImageIOHelper.readImage(fl.getPath());

        Region.drawEllipses(img2, regions, 0);

        MiscDebug.writeImage(img2, file + "_out_txt_");
    
        int nR = regions.get(0).size() + regions.get(1).size();
        
        System.out.println("num extracted regions=" + nR);
    
        assertEquals(317, nR);
        
    }
    
}
