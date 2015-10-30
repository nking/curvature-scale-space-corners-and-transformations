package algorithms.imageProcessing;

import algorithms.misc.MiscMath;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SegmentedImageHelperTest extends TestCase {
    
    public SegmentedImageHelperTest() {
    }

    public void testApplyEqualizationIfNeededByComparison() throws Exception {
        
        /*
        making 2 images, one with values 0 to 126 and the other with
        values 127 to 255
        */
        
        int w = 10;
        int h = 10;
        
        ImageExt img1 = new ImageExt(w, h);
        ImageExt img2 = new ImageExt(w, h);
        populateWithValues(img1, 0, 126);
        populateWithValues(img2, 127, 255);
        
        SegmentedImageHelper imgHelper1 = new SegmentedImageHelper(img1);
        SegmentedImageHelper imgHelper2 = new SegmentedImageHelper(img2);
        
        boolean applied = imgHelper1.applyEqualizationIfNeededByComparison(
            imgHelper2);
        
        assertTrue(applied);
        
        GreyscaleImage gsImg1 = imgHelper1.getGreyscaleImage();
        GreyscaleImage gsImg2 = imgHelper2.getGreyscaleImage();

        int min1, max1;
        if (gsImg1.is64Bit) {
            min1 = (int)MiscMath.findMinForByteCompressed(gsImg1.aL, gsImg1.len, 
                gsImg1.itemByteLength);
            max1 = (int)MiscMath.findMaxForByteCompressed(gsImg1.aL, gsImg1.len, 
                gsImg1.itemByteLength);
        } else {
            min1 = MiscMath.findMinForByteCompressed(gsImg1.a, gsImg1.len, 
                gsImg1.itemByteLength);
            max1 = MiscMath.findMaxForByteCompressed(gsImg1.a, gsImg1.len, 
                gsImg1.itemByteLength);
        }
        
        int min2, max2;
        if (gsImg2.is64Bit) {
            min2 = (int)MiscMath.findMinForByteCompressed(gsImg2.aL, gsImg2.len, 
                gsImg2.itemByteLength);
            max2 = (int)MiscMath.findMaxForByteCompressed(gsImg2.aL, gsImg2.len, 
                gsImg2.itemByteLength);
        } else {
            min2 = MiscMath.findMinForByteCompressed(gsImg2.a, gsImg2.len, 
                gsImg2.itemByteLength);
            max2 = MiscMath.findMaxForByteCompressed(gsImg2.a, gsImg2.len, 
                gsImg2.itemByteLength);
        }
               
        assertEquals(255, max1);
        assertEquals(255, max2);
        assertEquals(0, min1);
        assertEquals(1, min2);

    }

    private void populateWithValues(ImageExt img, int minV, int maxV) {
        
        int lastV = minV - 1;
        
        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                int v = lastV + 1;
                if (v > maxV) {
                    v = maxV;
                }
                img.setRGB(i, j, v, v, v);
                lastV = v;
            }
        }
    }

}
