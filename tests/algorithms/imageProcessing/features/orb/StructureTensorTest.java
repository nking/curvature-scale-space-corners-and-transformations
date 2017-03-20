package algorithms.imageProcessing.features.orb;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.StructureTensor;
import algorithms.misc.MiscDebug;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class StructureTensorTest extends TestCase {
    
    public StructureTensorTest() {
    }
    
    public void test00() {
        
        /*
        this test is from the scipy source code.
        (add licensing here)
        https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/corner.py
        
            --------
            >>> from skimage.feature import structure_tensor
            >>> square = np.zeros((5, 5))
            >>> square[2, 2] = 1
            >>> Axx, Axy, Ayy = structure_tensor(square, sigma=0.1)
            >>> Axx
            array([[ 0.,  0.,  0.,  0.,  0.],
                   [ 0.,  1.,  0.,  1.,  0.],
                   [ 0.,  4.,  0.,  4.,  0.],
                   [ 0.,  1.,  0.,  1.,  0.],
                   [ 0.,  0.,  0.,  0.,  0.]])
            """
        
        */
       
        int sz = 5;
                
        float[][] img = new float[sz][];
        for (int i = 0; i < img.length; ++i) {
            img[i] = new float[sz];
        }
        img[2][2] = 1.f;

        StructureTensor tensor = new StructureTensor(img, 
            0.f, true);
        
        float[][] Axx = tensor.getDXSquared();
        float[][] Ayy = tensor.getDYSquared();
        float[][] Axy = tensor.getDXDY();
        
        
        //String str = MiscDebug.getPrintRowMajor(Axx, "Axx=");
        //System.out.println(str);
        
        assertTrue(Math.abs((Axx[2][1]/Axx[1][1]) - 4.) < 0.001);
        assertTrue(Math.abs((Axx[2][1]/Axx[3][1]) - 4.) < 0.001);
        assertTrue(Math.abs((Axx[2][3]/Axx[1][3]) - 4.) < 0.001);
        assertTrue(Math.abs((Axx[2][3]/Axx[3][3]) - 4.) < 0.001);        
        
        for (int i = 0; i < Axx.length; ++i) {
            assertEquals(0.f, Axx[i][0]);
            assertEquals(0.f, Axx[i][2]);
            assertEquals(0.f, Axx[i][4]);
        }
        for (int j = 0; j < Axx[0].length; ++j) {
            assertEquals(0.f, Axx[0][j]);
            assertEquals(0.f, Axx[4][j]);
        }
        
        //str = MiscDebug.getPrintRowMajor(Axy, "Axy=");
        //System.out.println(str);
        
        assertEquals(0.0625f , Axy[1][1]);
        assertEquals(-0.0625f , Axy[3][1]);
        assertEquals(-0.0625f , Axy[1][3]);
        assertEquals(0.0625f , Axy[3][3]);
        for (int i = 0; i < Axy.length; ++i) {
            assertEquals(0.f, Math.abs(Axy[i][0]));
            assertEquals(0.f, Math.abs(Axy[i][2]));
            assertEquals(0.f, Math.abs(Axy[i][4]));
        }
        for (int j = 0; j < Axy[0].length; ++j) {
            assertEquals(0.f, Math.abs(Axy[0][j]));
            assertEquals(0.f, Math.abs(Axy[4][j]));
        }
        
        //str = MiscDebug.getPrintRowMajor(Ayy, "Ayy=");
        //System.out.println(str);
        
        assertTrue(Math.abs((Ayy[1][2]/Ayy[1][1]) - 4.) < 0.001);
        assertTrue(Math.abs((Ayy[1][2]/Ayy[1][3]) - 4.) < 0.001);
        assertTrue(Math.abs((Ayy[3][2]/Ayy[3][1]) - 4.) < 0.001);
        assertTrue(Math.abs((Ayy[3][2]/Ayy[3][3]) - 4.) < 0.001);        
        for (int i = 0; i < Axy.length; ++i) {
            assertEquals(0.f, Math.abs(Ayy[0][i]));
            assertEquals(0.f, Math.abs(Ayy[2][i]));
            assertEquals(0.f, Math.abs(Ayy[4][i]));
        }
        for (int j = 0; j < Axy[0].length; ++j) {
            assertEquals(0.f, Math.abs(Ayy[j][0]));
            assertEquals(0.f, Math.abs(Ayy[j][4]));
        }
    }
    
    public void est0() {
        
        int sz = 10;
        
        GreyscaleImage image = new GreyscaleImage(sz, sz);
        
        float[][] img = new float[sz][];
        for (int i = 0; i < img.length; ++i) {
            img[i] = new float[sz];
        }
        for (int i = 2; i < 8; ++i) {
            for (int j = 2; j < 8; ++j) {
                img[i][j] = 1.f;
            }
        }
        
        StructureTensor tensorComponents = new 
            StructureTensor(img, 1, false);
        
        float[][] detA = tensorComponents.getDeterminant();

        float[][] traceA = tensorComponents.getTrace();
        
        //float[][] hc = orb.cornerHarris(img, detA, traceA);
        
        //orb.debugPrint("hc=", hc);
        
        /*
         >>> from skimage.feature import corner_harris, corner_peaks
        >>> import numpy as np
        >>> square3 = np.zeros([10, 10])
        >>> square3[2:8, 2:8] = 1
        >>> square3.astype(int)
        array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
               [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
               [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
               [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
               [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
               [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
        
        
        */
        
       
    }
    
    
}
