package algorithms.imageProcessing.features.mser;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import gnu.trove.set.TIntSet;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MSER3Test extends TestCase {
    
    public MSER3Test() {
    }

    public void test0() throws Exception {

        /*
        An example run through of MSER with default settings and this as
        input image array:
             pixel values        pixel indexes
              0    1    2        
        0    250  200  129       00 01 02
        1    200  200  129       03 04 05
        2    150  150  129       06 07 08
        */
        int[] img = new int[]{250, 200, 129, 200, 200, 129, 150, 150, 129};
        MSER mser = new MSER();
        List<List<Region>> regionsList = mser.findRegions(img, 3, 3);
        
    }
    
}
