package algorithms.imageProcessing.features.orb;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.StructureTensor;
import algorithms.imageProcessing.StructureTensorR;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscDebug;
import algorithms.util.FormatArray;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class StructureTensorTest extends TestCase {
    
    public StructureTensorTest() {
    }

    static double[][] convertIntToDouble(float[][] a) {
        double[][] c = new double[a.length][];
        int i, j;
        for (i = 0; i < a.length; ++i) {
            c[i] = new double[a[0].length];
            for (j = 0; j < a[0].length; ++j) {
                c[i][j] = a[i][j];
            }
        }
        return c;
    }

    public void testDel() {
        System.out.println("testDel");
        float[][] gX = new float[3][3];
        gX[0] = new float[]{1f, 0, 1f};
        gX[1] = new float[]{0f, 1, 0f};
        gX[2] = new float[]{1f, 0, 1f};

        float[][] gY = new float[3][3];
        gY[0] = new float[]{1f, 0, 1f};
        gY[1] = new float[]{0f, 1, 0f};
        gY[2] = new float[]{1f, 0, 1f};

        double fNormX0 = MatrixUtil.frobeniusNorm(convertIntToDouble(gX));
        double fNormY0 = MatrixUtil.frobeniusNorm(convertIntToDouble(gY));

        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.applySobelX(gX);
        imageProcessor.applySobelY(gY);

        double fNormX = MatrixUtil.frobeniusNorm(convertIntToDouble(gX));
        double fNormY = MatrixUtil.frobeniusNorm(convertIntToDouble(gY));

        System.out.printf("gX=\n%s\n", FormatArray.toString(gX, "%.3e"));
        System.out.printf("gY=\n%s\n", FormatArray.toString(gY, "%.3e"));
        System.out.printf("fNormX=%.3e\n", fNormX);
        System.out.printf("fNormY=%.3e\n", fNormY);
        System.out.printf("fNormX0=%.3e\n", fNormX0);
        System.out.printf("fNormY0=%.3e\n", fNormY0);
        System.out.printf("fNormX0/fNormX=%.3e\n", fNormX0/fNormX);
        System.out.printf("fNormY0/fNormY=%.3e\n", fNormY0/fNormY);

        //MatrixUtil.multiply(gX, norm);
        //MatrixUtil.multiply(gY, norm);
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
        
        float norm = 4.f;
        
        float[][] Axx = tensor.getDXSquared();
        float[][] Ayy = tensor.getDYSquared();
        float[][] Axy = tensor.getDXDY();
        
        float[][] detA = tensor.getDeterminant();
        float[][] traceA = tensor.getTrace();
        
        for (int i = 0; i < detA.length; ++i) {
            for (int j = 0; j < detA[0].length; ++j) {
                assertEquals(0.0f, detA[i][j]);
            }
        }
        
        //String str = MiscDebug.getPrintRowMajor(detA, "detA=");
        //System.out.println(str);
        //str = MiscDebug.getPrintRowMajor(traceA, "traceA=");
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
        
        assertEquals(0.0625f * norm * norm, Axy[1][1]);
        assertEquals(-0.0625f * norm * norm, Axy[3][1]);
        assertEquals(-0.0625f * norm * norm, Axy[1][3]);
        assertEquals(0.0625f * norm * norm, Axy[3][3]);
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
    
    public void test0() {

        /*
         >>> from skimage.feature import corner_harris, corner_peaks
        >>> import numpy as np
        >>> square3 = np.zeros([10, 10])
        >>> square3[2:8, 2:8] = 1
        >>> square3.astype(int)
        array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
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
        double[][] img = MatrixUtil.zeros(11, 10);
        int i;
        int j;
        for (i = 3; i < 9; ++i) {
            for (j = 2; j < 8; ++j) {
                img[i][j] = 1;
            }
        }

        ImageProcessor imageProcessor = new ImageProcessor();
        int[][] keypoints = imageProcessor.calcHarrisCorners(img);

        int[][] expected = new int[4][2];
        expected[0] = new int[]{3, 2};
        expected[1] = new int[]{3, 7};
        expected[2] = new int[]{8, 2};
        expected[3] = new int[]{8, 7};

        assertEquals(expected.length, keypoints.length);

        for (i = 0; i < keypoints.length; ++i) {
            for (j = 0; j < 2; ++j) {
                assertEquals(expected[i][j], keypoints[i][j]);
            }
        }
       
    }
    
    
}
