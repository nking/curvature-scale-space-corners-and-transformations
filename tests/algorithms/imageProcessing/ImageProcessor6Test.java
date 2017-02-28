package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import algorithms.util.VeryLongBitString;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ImageProcessor6Test extends TestCase {

    public ImageProcessor6Test(String testName) {
        super(testName);
    }
    
    public void test0() throws Exception {
      
        /* pt values:
        red = 0 - 18
        orange = 18 - 40
        yellow = 41 - 60ish
        green = 61 - 106
        blue = 107 - 192
        purple = 193 - 255
        */
        
        String filePath = ResourceFinder.findFileInTestResources("colors.png");
        ImageExt img0 = ImageIOHelper.readImageExt(filePath);

        ImageProcessor imageProcessor = new ImageProcessor();
        GreyscaleImage ptImg = imageProcessor.createCIELUVTheta(img0, 255);
        
        for (int y = 10; y < 11; y+=2) {
            for (int x = 0; x < img0.getWidth(); x+=2) {
                int v = ptImg.getValue(x, y);
                System.out.format("(%d,%d) v=%d\n", x, y, v);
            }
        }
    }
    
}
