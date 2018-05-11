package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ImageProcessor6Test extends TestCase {
    
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
        GreyscaleImage ptImg = imageProcessor
            .createCIELUVTheta_WideRangeLightness(img0, 255);
        /*
        for (int y = 10; y < 11; y+=2) {
            for (int x = 0; x < img0.getWidth(); x+=2) {
                int v = ptImg.getValue(x, y);
                System.out.format("(%d,%d) v=%d\n", x, y, v);
            }
        }*/
        
        // ---- the above was a table for building Sky.java filters
        
        /*
        for the "H" of LCH made from LUV colorspace, wanting to
        te test that same color at larger brightness is same theta
        and that different colors are different theta.
        The above table print out already shows the later.
        */
        
        img0 = new ImageExt(10, 10);
        img0.setRGB(1, 1, 10, 150, 10);
        img0.setRGB(1, 2, 10, 180, 10);
        img0.setRGB(1, 3, 10, 210, 10);
        img0.setRGB(1, 4, 10, 250, 10);
        
        img0.setRGB(5, 1, 150, 10, 10);
        img0.setRGB(5, 2, 180, 10, 10);
        img0.setRGB(5, 3, 210, 10, 10);
        img0.setRGB(5, 4, 250, 10, 10);
   
        ptImg = imageProcessor
            .createCIELUVTheta_WideRangeLightness(img0, 255);
        
        assertTrue(Math.abs(ptImg.getValue(1, 1) - ptImg.getValue(1, 2)) <= 1);
        assertTrue(Math.abs(ptImg.getValue(1, 1) - ptImg.getValue(1, 3)) <= 1);
        assertTrue(Math.abs(ptImg.getValue(1, 1) - ptImg.getValue(1, 4)) <= 1);
        assertTrue(Math.abs(ptImg.getValue(5, 1) - ptImg.getValue(5, 2)) <= 1);
        assertTrue(Math.abs(ptImg.getValue(5, 1) - ptImg.getValue(5, 3)) <= 1);
        assertTrue(Math.abs(ptImg.getValue(5, 1) - ptImg.getValue(5, 4)) <= 1);
        
        /*
        System.out.println("greens:" + ptImg.getValue(1, 1) + " " +
            ptImg.getValue(1, 2) + " " + ptImg.getValue(1, 3) + " " +
            ptImg.getValue(1, 4));
        
        System.out.println("reds: " + ptImg.getValue(5, 1) + " " +
            ptImg.getValue(5, 2) + " " + ptImg.getValue(5, 3) + " " +
            ptImg.getValue(5, 4));
        */
        
        GreyscaleImage[] lch = imageProcessor.createLCHForLUV(img0);
        
        /*
        for "C" of LCH expect magnitude changes in green and red to be incresing
        magnitudes
        */
        
        int prevG = lch[1].getValue(1, 1);
        int prevR = lch[1].getValue(5, 1);
        for (int i = 2; i <= 4; ++i) {
            int vG = lch[1].getValue(1, i);
            int vR = lch[1].getValue(5, i);
            assertTrue(vG > prevG);
            assertTrue(vR > prevR);
            prevG = vG;
            prevR = vR;
        }
        
        /*
        System.out.println("greens:" + lch[1].getValue(1, 1) + " " +
            lch[1].getValue(1, 2) + " " + lch[1].getValue(1, 3) + " " +
            lch[1].getValue(1, 4));
        
        System.out.println("reds: " + lch[1].getValue(5, 1) + " " +
            lch[1].getValue(5, 2) + " " + lch[1].getValue(5, 3) + " " +
            lch[1].getValue(5, 4));
        */
        
        int a = lch[1].getValue(6, 4);
        imageProcessor.blur(lch[1], SIGMA.ONE);
        
        assertTrue(lch[1].getValue(6, 4) > a);
       
        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;
        for (int k = 0; k < dxs.length; ++k) {
            if (dxs[k] == 0) { continue;}
            int x2 = 5 + dxs[k];
            int y2 = 5 + dys[k];
            assertTrue(img0.getRGB(x2, y2) == 0);
        }
        
        GreyscaleImage rImg = img0.copyRedToGreyscale();
        
        imageProcessor.blur(img0, SIGMA.ONE, 0, 255);
        
        for (int k = 0; k < dxs.length; ++k) {
            if (dxs[k] == 0) { continue;}
            int x2 = 5 + dxs[k];
            int y2 = 5 + dys[k];
            assertTrue(img0.getRGB(x2, y2) > 0);
        }
        
        for (int k = 0; k < dxs.length; ++k) {
            if (dxs[k] == 0) { continue;}
            int x2 = 5 + dxs[k];
            int y2 = 5 + dys[k];
            assertTrue(rImg.getValue(x2, y2) == 0);
        }
        /*
        for (int y = 0; y < rImg.getWidth(); y++) {
            for (int x = 0; x < rImg.getWidth(); x++) {
                int v = rImg.getValue(x, y);
                System.out.format("before (%d,%d) v=%d\n", x, y, v);
            }
        }*/
        GreyscaleImage rImg2 = null;
        rImg2= imageProcessor.downSample(rImg, 5, 5, 0, 255);
        assertEquals(5, rImg2.getWidth());
        assertEquals(5, rImg2.getHeight());
        
        for (int y = 0; y < rImg2.getWidth(); y++) {
            for (int x = 0; x < rImg2.getWidth(); x++) {
                int v = rImg2.getValue(x, y);
                System.out.format("resampled (%d,%d) v=%d\n", x, y, v);
            }
        }
        //assertTrue(rImg2.getValue(0, 0) > 0 && rImg2.getValue(0,0) <
        //    10./4.);
        //assertTrue(rImg2.getValue(4, 0) == 0);
        //assertTrue(rImg2.getValue(2, 2) > 0 && rImg2.getValue(2,2) <
        //    100./4.);
        
        /*
        rImg2 = imageProcessor
            .downSample(rImg, 20, 20, 0, 255);
        assertEquals(20, rImg2.getWidth());
        assertEquals(20, rImg2.getHeight());
        for (int y = 0; y < rImg2.getWidth(); y++) {
            for (int x = 0; x < rImg2.getWidth(); x++) {
                int v = rImg2.getValue(x, y);
                System.out.format("resampled (%d,%d) v=%d\n", x, y, v);
            }
        }
        //assertTrue(rImg2.getValue(0, 0) > 0 && rImg2.getValue(0,0) <
        //    10./4.);
        assertTrue(rImg2.getValue(3, 0) == 0);
        //assertTrue(rImg2.getValue(2, 2) > 0 && rImg2.getValue(2,2) <
        //    100./4.);
        */
        
        /*
        img0.setRGB(5, 1, 150, 10, 10);
        img0.setRGB(5, 2, 180, 10, 10);
        img0.setRGB(5, 3, 210, 10, 10);
        img0.setRGB(5, 4, 250, 10, 10);
        */
        
        float[] f0 = imageProcessor.copyToFloat(rImg);
        assertEquals(150.f, f0[(1 * rImg.getWidth()) + 5]);
        assertEquals(0.f, f0[(1 * rImg.getWidth()) + 6]);
        
        int[] i0 = imageProcessor.copyToInt(rImg);
        assertEquals(150, i0[(1 * rImg.getWidth()) + 5]);
        assertEquals(0, i0[(1 * rImg.getWidth()) + 6]);
        
        double[][] d1 = imageProcessor.copy(rImg);
        assertEquals(150., d1[5][1]);
        assertEquals(0., d1[6][1]);
        
        float[][] fRowMaj = imageProcessor.copyToRowMajor(rImg);
        assertEquals(150.f, fRowMaj[1][5]);
        assertEquals(0.f, fRowMaj[1][6]);
    
    }
    
    
}
