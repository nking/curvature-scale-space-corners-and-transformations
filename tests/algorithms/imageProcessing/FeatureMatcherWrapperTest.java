package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.Collection;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertTrue;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class FeatureMatcherWrapperTest extends TestCase {
    
    public FeatureMatcherWrapperTest() {
    }
    
    public void test0() throws Exception {
        
        String fileName1, fileName2;
        
        for (int i = 0; i < 4;/*3;*/ ++i) {
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
                default: {
                    fileName1 = "campus_010.jpg";
                    fileName2 = "campus_011.jpg";
                    break;
                }
            }
            runCorrespondenceList(fileName1, fileName2);
        }
    }

    private void runCorrespondenceList(String fileName1, String fileName2) 
        throws Exception {
        
        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);
        idx = fileName2.lastIndexOf(".");
        String fileName2Root = fileName2.substring(0, idx);
        
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
        
        FeatureMatcherWrapper wrapper = new FeatureMatcherWrapper(img1, img2,
            fileName1Root);
        
        CorrespondenceList cl = wrapper.matchFeatures();
        
        assertNotNull(cl);
        
        Collection<PairInt> m1 = cl.getPoints1();
        Collection<PairInt> m2 = cl.getPoints2();
        GreyscaleImage gsImg1 = img1.copyToGreyscale();
        GreyscaleImage gsImg2 = img2.copyToGreyscale();
        MiscDebug.plotCorners(gsImg1, m1, "1_" + fileName1Root  + "_matched", 2);
        MiscDebug.plotCorners(gsImg2, m2, "2_" + fileName2Root + "_matched", 2);
                
    }
}
