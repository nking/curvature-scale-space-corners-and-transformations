package algorithms.imageProcessing.features.mser;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.misc.MiscDebug;
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
    
    public void est0() throws IOException {
        
        String[] files = new String[] {
            //"android_statues_01_sz1.jpg",
            "android_statues_02.jpg",
           // "android_statues_03_sz3.jpg",
           // "android_statues_04_sz1.jpg"
        };
            
        for (String file : files) {
        
            File fl = ResourceFinder.findFileInTmpData(file);
        
            GreyscaleImage img = 
                ImageIOHelper.readImageAsGrayScaleAvgRGB(fl.getPath());
            
            MSER mser = new MSER();
           
            List<List<Region>> regions = mser.findRegions(img);
            
            Image img2 = ImageIOHelper.readImage(fl.getPath());
            
            Region.drawEllipses(img2, regions, 1);
            
            MiscDebug.writeImage(img2, file + "_out_");
        }
    }
    
    public void test1() throws IOException {
        
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

        Region.drawEllipses(img2, regions, 1);

        MiscDebug.writeImage(img2, file + "_out_");
    }
    
}
