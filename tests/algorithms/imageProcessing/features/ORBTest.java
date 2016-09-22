package algorithms.imageProcessing.features;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageDisplayer;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.features.ORB.TwoDFloatArray;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ORBTest extends TestCase {
    
    public ORBTest() {
    }
    
    public void estPeakLocalMax() {
        
        /*
        using a test embedded in scipy code.
        
        (see their copyright in the ORB.java class documentation).
        
        image = np.zeros((5, 5))
        image[2][2] = 1.
        image[3][2] = 1.
        image[2][3] = 1.
        image[3][3] = 1.
        min_distance=1
        size = 2 * min_distance + 1
        image_max = ndi.maximum_filter(image, size=size, mode='constant')
        print image_max
        [[ 0.  0.  0.  0.  0.]
         [ 0.  1.  1.  1.  1.]
         [ 0.  1.  1.  1.  1.]
         [ 0.  1.  1.  1.  1.]
         [ 0.  1.  1.  1.  1.]]
        mask = image == image_max
        >>> print image
        [[ 0.  0.  0.  0.  0.]
         [ 0.  0.  0.  0.  0.]
         [ 0.  0.  1.  1.  0.]
         [ 0.  0.  1.  1.  0.]
         [ 0.  0.  0.  0.  0.]]
        >>> print mask
        [[ True  True  True  True  True]
         [ True False False False False]
         [ True False  True  True False]
         [ True False  True  True False]
         [ True False False False False]]
        
        mask &= image > 0.1
        print mask
        [[False False False False False]
         [False False False False False]
         [False False  True  True False]
         [False False  True  True False]
         [False False False False False]]
        */
        
        /*
        array([
           [ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  1.,  1.,  0.],
           [ 0.,  0.,  1.,  1.,  0.],
           [ 0.,  0.,  0.,  0.,  0.]])
    >>> peak_local_max(response)
    array([[2, 2],
           [2, 3],
           [3, 2],
           [3, 3]])
        */
        
        ORB orb = new ORB(4);
        
        float[][] img = new float[5][];
        for (int i = 0; i < img.length; ++i) {
            img[i] = new float[5];
        }
        for (int i = 2; i < 4; ++i) {
            for (int j = 2; j < 4; ++j) {
                img[i][j] = 1.f;
            }
        }
        
        TIntList outputKeypoints0 = new TIntArrayList();
        TIntList outputKeypoints1 = new TIntArrayList();
        orb.peakLocalMax(img, 1, 0.1f, outputKeypoints0, outputKeypoints1);
        
        //System.out.println("peakRowCols=" + peakRowCols.toString());
        
        /*
        3, 3, 3, 2, 2, 3, 2, 2
        array([[2, 2],
           [2, 3],
           [3, 2],
           [3, 3]])
        */
        
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(2, 2));
        expected.add(new PairInt(2, 3));
        expected.add(new PairInt(3, 2));
        expected.add(new PairInt(3, 3));
        
        for (int i = 0; i < outputKeypoints0.size(); ++i) {
            int x = outputKeypoints0.get(i);
            int y = outputKeypoints1.get(i);
            PairInt p = new PairInt(x, y);
            assertTrue(expected.remove(p));
        }
    }
    
    public void estCornerPeaks() {
        
        ORB orb = new ORB(10);
        
        float[][] img = new float[12][];
        for (int i = 0; i < img.length; ++i) {
            img[i] = new float[12];
        }
        for (int i = 3; i < 9; ++i) {
            for (int j = 3; j < 9; ++j) {
                img[i][j] = 1.f;
            }
        }
        
        /*
        using a test embedded in scipy code.
        
        print square
        [[ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
         [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
         [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
         [ 0.  0.  0.  1.  1.  1.  1.  1.  1.  0.  0.  0.]
         [ 0.  0.  0.  1.  1.  1.  1.  1.  1.  0.  0.  0.]
         [ 0.  0.  0.  1.  1.  1.  1.  1.  1.  0.  0.  0.]
         [ 0.  0.  0.  1.  1.  1.  1.  1.  1.  0.  0.  0.]
         [ 0.  0.  0.  1.  1.  1.  1.  1.  1.  0.  0.  0.]
         [ 0.  0.  0.  1.  1.  1.  1.  1.  1.  0.  0.  0.]
         [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
         [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
         [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]]

        from skimage.morphology import octagon
        from skimage.feature import (corner_fast, corner_peaks, corner_orientations)
        square = np.zeros((12, 12))
        square[3:9, 3:9] = 1
        square.astype(int)
        */
        
        float[][] cf = orb.cornerFast(img, 9, 0.15f);
        
        /*
        using a test embedded in scipy code.
        
        cf = corner_fast(square, 9); cf

        [[  0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.]
         [  0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.]
         [  0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.]
         [  0.   0.   0.  11.  10.   9.   9.  10.  11.   0.   0.   0.]
         [  0.   0.   0.  10.   9.   0.   0.   9.  10.   0.   0.   0.]
         [  0.   0.   0.   9.   0.   0.   0.   0.   9.   0.   0.   0.]
         [  0.   0.   0.   9.   0.   0.   0.   0.   9.   0.   0.   0.]
         [  0.   0.   0.  10.   9.   0.   0.   9.  10.   0.   0.   0.]
         [  0.   0.   0.  11.  10.   9.   9.  10.  11.   0.   0.   0.]
         [  0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.]
         [  0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.]
         [  0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.]]
        */
        assertEquals(11.f, cf[3][3]);
        assertEquals(10.f, cf[3][4]);
        assertEquals(9.f, cf[3][5]);
        assertEquals(9.f, cf[3][6]);
        assertEquals(10.f, cf[3][7]);
        assertEquals(11.f, cf[3][8]);
       
        assertEquals(10.f, cf[4][3]);
        assertEquals(9.f, cf[4][4]);
        assertEquals(9.f, cf[4][7]);
        assertEquals(10.f, cf[4][8]);
        
        assertEquals(9.f, cf[5][3]);
        assertEquals(9.f, cf[5][8]);
        
        assertEquals(9.f, cf[6][3]);
        assertEquals(9.f, cf[6][8]);
        
        assertEquals(10.f, cf[7][3]);
        assertEquals(9.f, cf[7][4]);
        assertEquals(9.f, cf[7][7]);
        assertEquals(10.f, cf[7][8]);
        
        assertEquals(11.f, cf[8][3]);
        assertEquals(10.f, cf[8][4]);
        assertEquals(9.f, cf[8][5]);
        assertEquals(9.f, cf[8][6]);
        assertEquals(10.f, cf[8][7]);
        assertEquals(11.f, cf[8][8]);        
        
        /*        
        min_distance = 1
        thresholdRel = 0.1
        size = 2 * min_distance + 1
        image_max = ndi.maximum_filter(cf, size=size, mode='constant')
        print image_max
        [[  0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.]
         [  0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.]
         [  0.   0.  11.  11.  11.  10.  10.  11.  11.  11.   0.   0.]
         [  0.   0.  11.  11.  11.  10.  10.  11.  11.  11.   0.   0.]
         [  0.   0.  11.  11.  11.  10.  10.  11.  11.  11.   0.   0.]
         [  0.   0.  10.  10.  10.   9.   9.  10.  10.  10.   0.   0.]
         [  0.   0.  10.  10.  10.   9.   9.  10.  10.  10.   0.   0.]
         [  0.   0.  11.  11.  11.  10.  10.  11.  11.  11.   0.   0.]
         [  0.   0.  11.  11.  11.  10.  10.  11.  11.  11.   0.   0.]
         [  0.   0.  11.  11.  11.  10.  10.  11.  11.  11.   0.   0.]
         [  0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.]
         [  0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.]]

        mask = cf == image_max
        print mask
        print mask
        [[ True  True  True  True  True  True  True  True  True  True  True  True]
         [ True  True  True  True  True  True  True  True  True  True  True  True]
         [ True  True False False False False False False False False  True  True]
         [ True  True False  True False False False False  True False  True  True]
         [ True  True False False False False False False False False  True  True]
         [ True  True False False False False False False False False  True  True]
         [ True  True False False False False False False False False  True  True]
         [ True  True False False False False False False False False  True  True]
         [ True  True False  True False False False False  True False  True  True]
         [ True  True False False False False False False False False  True  True]
         [ True  True  True  True  True  True  True  True  True  True  True  True]
         [ True  True  True  True  True  True  True  True  True  True  True  True]]

         mask &= cf > max(thresholds)
        >>> print mask
        [[False False False False False False False False False False False False]
         [False False False False False False False False False False False False]
         [False False False False False False False False False False False False]
         [False False False  True False False False False  True False False False]
         [False False False False False False False False False False False False]
         [False False False False False False False False False False False False]
         [False False False False False False False False False False False False]
         [False False False False False False False False False False False False]
         [False False False  True False False False False  True False False False]
         [False False False False False False False False False False False False]
         [False False False False False False False False False False False False]
         [False False False False False False False False False False False False]]
        */    
        
        TIntList keypoints0 = new TIntArrayList();
        TIntList keypoints1 = new TIntArrayList();
        
        orb.cornerPeaks(cf, 1, keypoints0, keypoints1);
        
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(3, 3));
        expected.add(new PairInt(3, 8));
        expected.add(new PairInt(8, 3));
        expected.add(new PairInt(8, 8));
        
        /*
        expecting
        array([[3, 3],
               [3, 8],
               [8, 3],
               [8, 8]])
        
        */
        
        for (int i = 0; i < keypoints0.size(); ++i) {
            int x = keypoints0.get(i);
            int y = keypoints1.get(i);
            PairInt p = new PairInt(x, y);
            assertTrue(expected.remove(p));
        } 
        assertTrue(expected.isEmpty());
        
        double[] orientation = orb.cornerOrientations(img, keypoints0, keypoints1);
        
        assertEquals(keypoints0.size(), orientation.length);
        Map<PairInt, Integer> expectedOrientations = new HashMap<PairInt, Integer>();
        expectedOrientations.put(new PairInt(3, 3), Integer.valueOf(45));
        expectedOrientations.put(new PairInt(3, 8), Integer.valueOf(135));
        expectedOrientations.put(new PairInt(8, 3), Integer.valueOf(-45));
        expectedOrientations.put(new PairInt(8, 8), Integer.valueOf(-135));
        
        //System.out.println(Arrays.toString(orientation));
        
        for (int i = 0; i < keypoints0.size(); ++i) {
            int x = keypoints0.get(i);
            int y = keypoints1.get(i);
            PairInt p = new PairInt(x, y);
            Integer expAng = expectedOrientations.get(p);
            double ang = orientation[i] * 180./Math.PI;
            //System.out.println("found=" + ang + " expected=" + expAng);
            assertTrue(Math.abs(expAng.intValue() - ang) < 0.001);
        }
        
        /*
        repeating a few steps:
        cp = corner_peaks(cf, min_distance=1)
        >>> print cp
        [[3 3]
         [3 8]
         [8 3]
         [8 8]]
        mask = skimage.feature.orb._mask_border_keypoints(cf.shape, cp, distance=3);
        
        >>> print mask
        [ True  True  True  True]
        >>> print cp[mask]
        [[3 3]
         [3 8]
         [8 3]
         [8 8]]
        cp = cp[mask]
        orientations = skimage.feature.orb.corner_orientations(square, cp, OFAST_MASK)
        
        print orientations
        [ 0.78539816  2.35619449 -0.78539816 -2.35619449]
        
        harris_response = corner_harris(square, method='k', k=0.04);
        >>> print harris_response.shape
        (12, 12)
        
        responses = harris_response[cp[:, 0], cp[:, 1]]
        >>> print responses
        [ 21.4776577  21.4776577  21.4776577  21.4776577]
        
        */
    }
    
    public void estTensor() {
        
        /*
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
        */
        
        ORB orb = new ORB(100);
        
        int sz = 5;
        
        float[][] img = new float[sz][];
        for (int i = 0; i < img.length; ++i) {
            img[i] = new float[sz];
        }
        img[2][2] = 1.f;
        
        TwoDFloatArray[] tensorComponents = orb.structureTensor(img, 0.1f);
        
        float[][] axx = tensorComponents[0].a;
        //for (int i = 0; i < axx.length; ++i) {
        //    System.out.println("axx[" + i + "]=" + Arrays.toString(axx[i]));
        //}
        
        assertTrue(Math.abs(axx[0][0]) < 0.01);
        assertTrue(Math.abs(axx[1][0]) < 0.01);
        assertTrue(Math.abs(axx[2][0]) < 0.01);
        assertTrue(Math.abs(axx[3][0]) < 0.01);
        assertTrue(Math.abs(axx[4][0]) < 0.01);
        
        assertTrue(Math.abs(axx[0][1]) < 0.01);
        assertTrue(Math.abs(axx[4][1]) < 0.01);
        
        assertTrue(Math.abs(axx[0][2]) < 0.01);
        assertTrue(Math.abs(axx[1][2]) < 0.01);
        assertTrue(Math.abs(axx[2][2]) < 0.01);
        assertTrue(Math.abs(axx[3][2]) < 0.01);
        assertTrue(Math.abs(axx[4][2]) < 0.01);
        
        assertTrue(Math.abs(axx[0][3]) < 0.01);
        assertTrue(Math.abs(axx[4][3]) < 0.01);
        
        assertTrue(Math.abs(axx[0][4]) < 0.01);
        assertTrue(Math.abs(axx[1][4]) < 0.01);
        assertTrue(Math.abs(axx[2][4]) < 0.01);
        assertTrue(Math.abs(axx[3][4]) < 0.01);
        assertTrue(Math.abs(axx[4][4]) < 0.01);
        
        /*
        array([[ 0.,  0.,  0.,  0.,  0.],
               [ 0.,  1.,  0.,  1.,  0.],
               [ 0.,  4.,  0.,  4.,  0.],
               [ 0.,  1.,  0.,  1.,  0.],
               [ 0.,  0.,  0.,  0.,  0.]])
        */
        
        int factor = 16;
        //System.out.println(Math.round(axx[1][1])/factor);
        //System.out.println(Math.round(axx[2][1])/factor);
        //System.out.println(Math.round(axx[3][1])/factor);
        assertTrue((Math.round(axx[1][1])/factor - 1.) < 0.01);
        assertTrue((Math.round(axx[2][1])/factor - 4.) < 0.01);
        assertTrue((Math.round(axx[3][1])/factor - 1.) < 0.01);
    }
    
    public void estCornerHarris() {
        
        ORB orb = new ORB(100);
        
        int sz = 10;
        
        float[][] img = new float[sz][];
        for (int i = 0; i < img.length; ++i) {
            img[i] = new float[sz];
        }
        for (int i = 2; i < 8; ++i) {
            for (int j = 2; j < 8; ++j) {
                img[i][j] = 1.f;
            }
        }
        
        float[][] hc = orb.cornerHarris(img);
        
        orb.debugPrint("hc=", hc);
        
        /*
         >>> from skimage.feature import corner_harris, corner_peaks
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
        >>> corner_peaks(corner_harris(square3), min_distance=1)
        array([[2, 2],
               [2, 7],
               [7, 2],
               [7, 7]])
        
        */
        
        TIntList keypoints0 = new TIntArrayList();
        TIntList keypoints1 = new TIntArrayList();
        
        orb.cornerPeaks(hc, 1, keypoints0, keypoints1);
        
        //System.out.println("coords=" + coords);
        
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(2, 2));
        expected.add(new PairInt(2, 7));
        expected.add(new PairInt(7, 2));
        expected.add(new PairInt(7, 7));
        
        for (int i = 0; i < keypoints0.size(); ++i) {
            int x = keypoints0.get(i);
            int y = keypoints1.get(i);
            PairInt p = new PairInt(x, y);
            assertTrue(expected.remove(p));
        } 
        assertTrue(expected.isEmpty());
    }
    
    public void testKeypoints_1() throws Exception {
        
        String fileName = "susan-in_plus.png";   
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        //ImageExt img = ImageIOHelper.readImageExt(filePath);
        ImageExt img = img0.copyToImageExt();
        
        //NOTE: the ridges are picked up well with reduced threshold
        
        ORB orb = new ORB(500);
        //orb.overrideFastN(12);
        orb.overrideFastThreshold(0.01f);
        
        orb.detectAndExtract(img);
        
        TIntList keypoints0 = orb.getAllKeyPoints0();
        TIntList keypoints1 = orb.getAllKeyPoints1();
        TFloatList scales = orb.getAllScales();
        TFloatList responses = orb.getAllHarrisResponses();
        TDoubleList orientations = orb.getAllOrientations();
        //System.out.println("keypoints0=" + keypoints0.toString());
        //System.out.println("keypoints1=" + keypoints1.toString());
        //System.out.println("scales=" + scales.toString());
        //System.out.println("responses=" + responses.toString());
        //System.out.println("orientations=" + orientations.toString());
        
        for (int i = 0; i < keypoints0.size(); ++i) {
            int y = keypoints0.get(i);
            int x = keypoints1.get(i);
            ImageIOHelper.addPointToImage(x, y, img0, 2, 255, 0, 0);
        }
        MiscDebug.writeImage(img0, "orb_keypoints");
        
    }
}
