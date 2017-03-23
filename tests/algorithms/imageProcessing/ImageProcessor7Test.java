package algorithms.imageProcessing;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ImageProcessor7Test extends TestCase {

    public ImageProcessor7Test(String testName) {
        super(testName);
    }
    
    public void testKernels() {
        
        ImageExt img0 = new ImageExt(10, 10);
        for (int i = 4; i < 7; ++i) {
            for (int j = 4; j < 7; ++j) {
                img0.setRGB(i, j, 10, 150, 10);
            }
        }
        
        GreyscaleImage img = img0.copyToGreyscale2();
        ImageProcessor imageProcessor = new ImageProcessor();
        img = imageProcessor.applyLaplacianKernel(img, 0, 1000);
        assertEquals(0, img.getValue(0, 6));
        assertEquals(0, img.getValue(1, 6));
        assertEquals(0, img.getValue(2, 6));
        assertEquals(0, img.getValue(9, 6));
        assertTrue(img.getValue(4, 4) != 0);
       
        float[][] img1 = imageProcessor.copyToFloat2D(
            img0.copyToGreyscale2());
        
        IKernel kernel = new Laplacian();
        Kernel kernelXY = kernel.getKernel();
        float norm = kernel.getNormalizationFactor();

        imageProcessor.applyKernel(img1, kernelXY, norm);
        
        GreyscaleImage img2 = imageProcessor.applyKernel(
            img0.copyToGreyscale2(), kernelXY, norm, 0, 1000);
        
        assertEquals(0, img2.getValue(9, 7));
        assertEquals(0, img2.getValue(8, 7));
        assertTrue(img2.getValue(4, 4) != 0);
       
        assertEquals(0.f, img1[9][7]);
        assertEquals(0.f, img1[8][7]);
        assertTrue(img1[7][7] != 0.f);
       
        // ------
        img1 = imageProcessor.copyToFloat2D(
            img0.copyToGreyscale2());
         
        imageProcessor.applyKernel1D(img1, new float[]{0.5f, 0.f, -0,5f},
            true);
        imageProcessor.applyKernel1D(img1, new float[]{0.5f, 0.f, -0,5f},
            false);
        assertEquals(0.f, img1[0][6]);
        assertEquals(0.f, img1[1][6]);
        assertEquals(0.f, img1[2][6]);
        assertEquals(0.f, img1[9][6]);
        assertTrue(img1[7][6] != 0.f);
        
        // -----
        img2 = img0.copyToGreyscale2();
        imageProcessor.applyKernel1D(img2, new float[]{0.5f, 0.f, -0,5f},
            true, 0, 255);
        imageProcessor.applyKernel1D(img2, new float[]{0.5f, 0.f, -0,5f},
            false, 0, 255);
        assertEquals(0, img2.getValue(0, 6));
        assertEquals(0, img2.getValue(1, 6));
        assertEquals(0, img2.getValue(2, 6));
        assertEquals(0, img2.getValue(9, 6));
        assertTrue(img2.getValue(7, 6) != 0);
       
        int[][] img4 = imageProcessor.copyToInt2D(img0.copyToGreyscale2());
        imageProcessor.applyKernelTwo1Ds(img4, new float[]{0.5f, 0.f, -0,5f});
        assertEquals(0, img4[0][6]);
        assertEquals(0, img4[1][6]);
        assertEquals(0, img4[2][6]);
        assertEquals(0, img4[9][6]);
        assertTrue(img4[7][6] != 0);
    
        img4 = imageProcessor.copyToInt2D(img0.copyToGreyscale2());
        imageProcessor.applyKernel1D(img4, new float[]{0.5f, 0.f, -0,5f}, true);
        imageProcessor.applyKernel1D(img4, new float[]{0.5f, 0.f, -0,5f}, false);
        assertEquals(0, img4[0][6]);
        assertEquals(0, img4[1][6]);
        assertEquals(0, img4[2][6]);
        assertEquals(0, img4[9][6]);
        assertTrue(img4[7][6] != 0);
        
        img2 = img0.copyToGreyscale2().copyToFullRangeIntImage();
        imageProcessor.applyKernel1D(img2, new 
            float[]{0.5f, 0.f, -0,5f},
            true);
        imageProcessor.applyKernel1D(img2, new float[]{0.5f, 0.f, -0,5f},
            false);
        assertEquals(0, img2.getValue(0, 6));
        assertEquals(0, img2.getValue(1, 6));
        assertEquals(0, img2.getValue(2, 6));
        assertEquals(0, img2.getValue(9, 6));
        assertTrue(img2.getValue(7, 6) != 0);
        
        
        img2 = img0.copyToGreyscale2().copyToFullRangeIntImage();
        imageProcessor.applySecondDerivGaussian(img2, SIGMA.ZEROPOINTSEVENONE,
            0, 1000);
        assertEquals(0, img2.getValue(7, 6));
        
        
        img2 = img0.copyToGreyscale2().copyToFullRangeIntImage();
        int[] img5 = imageProcessor.performSecondDerivGaussian(img2,
            SIGMA.ZEROPOINTSEVENONE);
        assertEquals(0, img5[(6 * 10) + 7]);
        
    }
    
}
