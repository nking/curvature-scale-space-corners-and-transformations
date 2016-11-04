package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageDisplayer;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.StructureTensor;
import algorithms.imageProcessing.features.ORB.Descriptors;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import algorithms.util.TwoDFloatArray;
import algorithms.util.VeryLongBitString;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
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
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        TIntList outputKeypoints0 = new TIntArrayList();
        TIntList outputKeypoints1 = new TIntArrayList();
        imageProcessor.peakLocalMax(img, 1, 0.1f, outputKeypoints0, outputKeypoints1);
        
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
        
        TDoubleList orientation = 
            orb.cornerOrientations(img, keypoints0, keypoints1);
        
        assertEquals(keypoints0.size(), orientation.size());
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
            double ang = orientation.get(i) * 180./Math.PI;
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
        orb.overrideToUseSmallestPyramid();
        
        int sz = 5;
        
        float[][] img = new float[sz][];
        for (int i = 0; i < img.length; ++i) {
            img[i] = new float[sz];
        }
        img[2][2] = 1.f;
        
        StructureTensor tensorComponents = new 
            StructureTensor(img, 0.1f, false);
                
        float[][] axx = tensorComponents.getDXSquared();
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
        orb.overrideToUseSmallestPyramid();
        
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
        
        StructureTensor tensorComponents = new 
            StructureTensor(img, 1, false);
        
        float[][] detA = tensorComponents.getDeterminant();

        float[][] traceA = tensorComponents.getTrace();
        
        float[][] hc = orb.cornerHarris(img, detA, traceA);
        
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
    
    public void estKeypoints_1() throws Exception {
        
        String fileName = "susan-in_plus.png";  
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        //ImageExt img = ImageIOHelper.readImageExt(filePath);
        ImageExt img = img0.copyToImageExt();
        
        //NOTE: the ridges are picked up well with reduced threshold
        
        ORB orb = new ORB(500);
        orb.overrideToUseSmallestPyramid();
        //orb.overrideFastN(12);
        orb.overrideFastThreshold(0.01f);
        orb.overrideToNotCreateDescriptors();
        
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
        MiscDebug.writeImage(img0, "orb_keypoints_01");
        
    }
    
    public void estKeypoints_2() throws Exception {
        
        String fileName = "susan-in_plus.png";  
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        //ImageExt img = ImageIOHelper.readImageExt(filePath);
        ImageExt img = img0.copyToImageExt();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(img, SIGMA.getValue(SIGMA.ONE), 0, 255);
        
        //NOTE: the ridges are picked up well with reduced threshold
        
        ORB orb = new ORB(500);
        orb.overrideToUseSmallestPyramid();
        //orb.overrideFastN(12);
        //orb.overrideFastThreshold(0.01f);
        orb.overrideToNotCreateDescriptors();
        orb.overrideToAlsoCreate1stDerivKeypoints();
        
        orb.detectAndExtract(img);
        
        TIntList keypoints0 = orb.getAllKeyPoints0();
        TIntList keypoints1 = orb.getAllKeyPoints1();
        
        for (int i = 0; i < keypoints0.size(); ++i) {
            int y = keypoints0.get(i);
            int x = keypoints1.get(i);
            ImageIOHelper.addPointToImage(x, y, img0, 2, 255, 0, 0);
        }
        MiscDebug.writeImage(img0, "orb_keypoints_02");
        
    }
    
    public void estKeypoints_3() throws Exception {
        
        //NOTE: can see there may still be errors in the code
        // or need to allow relaxation of border distance
        
        String fileName = "checkerboard_01.jpg";   
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        //ImageExt img = ImageIOHelper.readImageExt(filePath);
        ImageExt img = img0.copyToImageExt();
                
        ORB orb = new ORB(500);
        orb.overrideToUseSmallestPyramid();
        //orb.overrideFastN(12);
        //orb.overrideFastThreshold(0.01f);
        orb.overrideToAlsoCreate1stDerivKeypoints();
        
        orb.detectAndExtract(img);
        
        TIntList keypoints0 = orb.getAllKeyPoints0();
        TIntList keypoints1 = orb.getAllKeyPoints1();
        
        for (int i = 0; i < keypoints0.size(); ++i) {
            int y = keypoints0.get(i);
            int x = keypoints1.get(i);
            ImageIOHelper.addPointToImage(x, y, img0, 2, 255, 0, 0);
        }
        MiscDebug.writeImage(img0, "orb_keypoints_03");
        
    }
    
    public void estKeypoints_4() throws Exception {
        
        String fileName = "blox.gif";
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        //ImageExt img = ImageIOHelper.readImageExt(filePath);
        ImageExt img = img0.copyToImageExt();
        
        //NOTE: the ridges are picked up well with reduced threshold
        
        ORB orb = new ORB(500);
        orb.overrideToUseSmallestPyramid();
        //orb.overrideFastThreshold(0.01f);
        orb.detectAndExtract(img);
        orb.overrideToNotCreateDescriptors();
        orb.overrideToAlsoCreate1stDerivKeypoints();
        
        TIntList keypoints0 = orb.getAllKeyPoints0();
        TIntList keypoints1 = orb.getAllKeyPoints1();
        
        for (int i = 0; i < keypoints0.size(); ++i) {
            int y = keypoints0.get(i);
            int x = keypoints1.get(i);
            ImageIOHelper.addPointToImage(x, y, img0, 2, 255, 0, 0);
        }
        MiscDebug.writeImage(img0, "orb_keypoints_04");        
    }
    
    public void estKeypoints_5() throws Exception {
        
        //NOTE: can see there may still be errors in the code
        // or need to allow relaxation of border distance
        
        String fileName = "android_statues_01.jpg";   
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        //ImageExt img = ImageIOHelper.readImageExt(filePath);
        ImageExt img = img0.copyToImageExt();
                
        ORB orb = new ORB(500);
        orb.overrideToUseSmallestPyramid();
        //orb.overrideFastN(12);
        orb.overrideFastThreshold(0.01f);
        orb.overrideToNotCreateDescriptors();
        
        orb.detectAndExtract(img);
        
        TIntList keypoints0 = orb.getAllKeyPoints0();
        TIntList keypoints1 = orb.getAllKeyPoints1();
        
        for (int i = 0; i < keypoints0.size(); ++i) {
            int y = keypoints0.get(i);
            int x = keypoints1.get(i);
            ImageIOHelper.addPointToImage(x, y, img0, 2, 255, 0, 0);
        }
        MiscDebug.writeImage(img0, "orb_keypoints_05");
        
    }
    
    public void estKeypoints_6() throws Exception {
        
        //NOTE: can see there may still be errors in the code
        // or need to allow relaxation of border distance
        
        String fileName = "android_statues_01.jpg";   
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        //ImageExt img = ImageIOHelper.readImageExt(filePath);
        ImageExt img = img0.copyToImageExt();
                
        ORB orb = new ORB(500);
        orb.overrideToUseSmallestPyramid();
        //orb.overrideFastN(12);
        orb.overrideToNotCreateDescriptors();
        orb.overrideToAlsoCreate1stDerivKeypoints();
        
        orb.detectAndExtract(img);
        
        TIntList keypoints0 = orb.getAllKeyPoints0();
        TIntList keypoints1 = orb.getAllKeyPoints1();
        
        for (int i = 0; i < keypoints0.size(); ++i) {
            int y = keypoints0.get(i);
            int x = keypoints1.get(i);
            ImageIOHelper.addPointToImage(x, y, img0, 2, 255, 0, 0);
        }
        MiscDebug.writeImage(img0, "orb_keypoints_06");
        
    }
    
     public void estKeypoints_7() throws Exception {
        
        //NOTE: can see there may still be errors in the code
        // or need to allow relaxation of border distance
        
        String fileName = "house.gif";   
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        //ImageExt img = ImageIOHelper.readImageExt(filePath);
        ImageExt img = img0.copyToImageExt();
                
        ORB orb = new ORB(500);
        orb.overrideToUseSmallestPyramid();
        //orb.overrideFastN(12);
        orb.overrideFastThreshold(0.01f);
        orb.overrideToNotCreateDescriptors();
        
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
        MiscDebug.writeImage(img0, "orb_keypoints_07");
        
    }
     
    public void estKeypoints_8() throws Exception {
        
        //NOTE: can see there may still be errors in the code
        // or need to allow relaxation of border distance
        
        String fileName = "house.gif";   
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        //ImageExt img = ImageIOHelper.readImageExt(filePath);
        ImageExt img = img0.copyToImageExt();
                
        ORB orb = new ORB(500);
        orb.overrideToUseSmallestPyramid();
        //orb.overrideFastN(12);
        //orb.overrideFastThreshold(0.01f);
        orb.overrideToNotCreateDescriptors();
        orb.overrideToAlsoCreate1stDerivKeypoints();
        
        orb.detectAndExtract(img);
        
        TIntList keypoints0 = orb.getAllKeyPoints0();
        TIntList keypoints1 = orb.getAllKeyPoints1();
        
        for (int i = 0; i < keypoints0.size(); ++i) {
            int y = keypoints0.get(i);
            int x = keypoints1.get(i);
            ImageIOHelper.addPointToImage(x, y, img0, 2, 255, 0, 0);
        }
        MiscDebug.writeImage(img0, "orb_keypoints_08");
        
    }
     
    public void estKeypoints_9() throws Exception {
        
        //NOTE: can see there may still be errors in the code
        // or need to allow relaxation of border distance
        
        String fileName = "android_statues_02.jpg";   
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        //ImageExt img = ImageIOHelper.readImageExt(filePath);
        ImageExt img = img0.copyToImageExt();
                
        ORB orb = new ORB(2000);
        orb.overrideToUseSmallestPyramid();
        //orb.overrideFastN(12);
        orb.overrideFastThreshold(0.01f);
        orb.overrideToNotCreateDescriptors();
        
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
        MiscDebug.writeImage(img0, "orb_keypoints_09");
        
    }
    
    public void estKeypoints_10() throws Exception {
        
        //NOTE: can see there may still be errors in the code
        // or need to allow relaxation of border distance
        
        String fileName = "android_statues_02.jpg";   
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        //ImageExt img = ImageIOHelper.readImageExt(filePath);
        ImageExt img = img0.copyToImageExt();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(img, SIGMA.getValue(SIGMA.TWO), 0, 255);
                
        ORB orb = new ORB(2000);
        orb.overrideToUseSmallestPyramid();
        //orb.overrideFastN(12);
        orb.overrideToNotCreateDescriptors();
        orb.overrideToAlsoCreate1stDerivKeypoints();
        
        orb.detectAndExtract(img);
        
        TIntList keypoints0 = orb.getAllKeyPoints0();
        TIntList keypoints1 = orb.getAllKeyPoints1();
        
        for (int i = 0; i < keypoints0.size(); ++i) {
            int y = keypoints0.get(i);
            int x = keypoints1.get(i);
            ImageIOHelper.addPointToImage(x, y, img0, 2, 255, 0, 0);
        }
        MiscDebug.writeImage(img0, "orb_keypoints_10");
        
    }
   
    public void estKeypoints_11() throws Exception {
        
        //NOTE: can see there may still be errors in the code
        // or need to allow relaxation of border distance
        
        String fileName = "lab.gif";   
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        //ImageExt img = ImageIOHelper.readImageExt(filePath);
        ImageExt img = img0.copyToImageExt();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(img, SIGMA.getValue(SIGMA.TWO), 0, 255);
                
        ORB orb = new ORB(2000);
        orb.overrideToUseSmallestPyramid();
        //orb.overrideFastN(12);
        orb.overrideToNotCreateDescriptors();
        
        orb.detectAndExtract(img);
        
        TIntList keypoints0 = orb.getAllKeyPoints0();
        TIntList keypoints1 = orb.getAllKeyPoints1();
        
        for (int i = 0; i < keypoints0.size(); ++i) {
            int y = keypoints0.get(i);
            int x = keypoints1.get(i);
            ImageIOHelper.addPointToImage(x, y, img0, 2, 255, 0, 0);
        }
        MiscDebug.writeImage(img0, "orb_keypoints_11");
        
    }
   
    public void estKeypoints_12() throws Exception {
        
        //NOTE: can see there may still be errors in the code
        // or need to allow relaxation of border distance
        
        String fileName = "lab.gif";   
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        //ImageExt img = ImageIOHelper.readImageExt(filePath);
        ImageExt img = img0.copyToImageExt();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(img, SIGMA.getValue(SIGMA.TWO), 0, 255);
                
        ORB orb = new ORB(2000);
        orb.overrideToUseSmallestPyramid();
        //orb.overrideFastN(12);
        orb.overrideToNotCreateDescriptors();
        orb.overrideToAlsoCreate1stDerivKeypoints();
        
        orb.detectAndExtract(img);
        
        TIntList keypoints0 = orb.getAllKeyPoints0();
        TIntList keypoints1 = orb.getAllKeyPoints1();
        
        for (int i = 0; i < keypoints0.size(); ++i) {
            int y = keypoints0.get(i);
            int x = keypoints1.get(i);
            ImageIOHelper.addPointToImage(x, y, img0, 2, 255, 0, 0);
        }
        MiscDebug.writeImage(img0, "orb_keypoints_12");
        
    }
   
    public void estDescriptors() throws IOException {
        
        /*
        2 test images:
           one a figure 8 of length 32 and other same image rotated by 90.
        
        positive test:
           test that descriptors of the center match each other
        negative test:
           test that descriptor on edge does not match well a descriptor
              from the center of the figure.
        */
        
        Image img = getFigureEight();
        int w = img.getWidth();
        int h = img.getHeight();
        
        Transformer tr = new Transformer();
        TransformationParameters params = new TransformationParameters();
        params.setOriginX(w/2);
        params.setOriginY(w/2);
        params.setRotationInDegrees(90);        
        Image img90 = tr.applyTransformation(img, params, w, h);
        
        //MiscDebug.writeImage(img90, "_rotated_");
        
        ORB orb = new ORB(2);
        orb.overrideToUseSmallestPyramid();
        orb.detectAndExtract(img);
        List<PairInt> kp0 = orb.getAllKeyPoints();
        Descriptors desc0 = orb.getAllDescriptors();

        ORB orb1 = new ORB(2);
        orb1.detectAndExtract(img90);
        List<PairInt> kp90 = orb1.getAllKeyPoints();
        Descriptors desc90 = orb1.getAllDescriptors();
            
        int[][] matches = ORB.matchDescriptors(
            desc0.descriptors, desc90.descriptors, 
            kp0, kp90);
        
        // center is 31, 32    side=17,34
        //           32, 33         34,47      
     
        int m11 = matches[0][0];
        int m12 = matches[0][1];
        
        int m21 = matches[1][0];
        int m22 = matches[1][1];
        
        assertTrue(
            (m11 == 1 && m12 == 1) ||
            (m11 == 0 && m12 == 0));
               
        assertTrue(
            (m21 == 1 && m22 == 1) ||
            (m21 == 0 && m22 == 0));
        
        int[][] costMatrix = ORB.calcDescriptorCostMatrix(
            desc0.descriptors, desc90.descriptors);

        int c00 = costMatrix[0][0];
        int c11 = costMatrix[1][1];
        
        int c01 = costMatrix[0][1];
        int c10 = costMatrix[1][0];
        
        assertEquals(0, c00);
        assertEquals(0, c11);
        
        assertTrue(c01 > 0);
        assertTrue(c10 > 0);
    
        int np = 8;
        orb = new ORB(np);
        orb.detectAndExtract(img);
        kp0 = orb.getAllKeyPoints();
        desc0 = orb.getAllDescriptors();
    
        // --- looking at same image as img, but size 85% and 70%
        for (int i = 0; i < 2; ++i) {
            Image img2;
            if (i == 0) {
                img2 = getFigureEight1();
            } else {
                img2 = getFigureEight2();
            }
            ORB orb2 = new ORB(np);
            orb2.detectAndExtract(img2);
            List<PairInt> kp2 = orb2.getAllKeyPoints();
            Descriptors desc2 = orb2.getAllDescriptors();
            List<TDoubleList> or2 = orb2.getOrientationsList();

            int[][] costMatrix2 = ORB.calcDescriptorCostMatrix(
                desc0.descriptors, desc2.descriptors);

            boolean foundCenterCostIsZero = false;
            
            for (int j = 0; j < costMatrix2.length; ++j) {
                for (int k = 0; k < costMatrix2[j].length; ++k) {
                    int c = costMatrix2[j][k];
                    //System.out.println(
                    //    "i=" + i + 
                    //    String.format(" %s:%s ", kp0.get(j),
                    //    kp2.get(k)) +
                    //    " c[" + j +"]["+k+"] c=" + c);
                    
                    if ((Math.abs(kp0.get(j).getX() - 32) < 3)
                        && (Math.abs(kp0.get(j).getY() - 32) < 3)
                        && (Math.abs(kp2.get(k).getX() - 32) < 3)
                        && (Math.abs(kp2.get(k).getY() - 32) < 3)) {
                        if (c == 0) {
                            foundCenterCostIsZero = true;
                        }
                    }
                }
            }
            
            assertTrue(foundCenterCostIsZero);
        }       
    }
    
    public void testHSVDescriptors() throws IOException {
        
        /*
        test an image with 2 color rectangles having very different
        HSV colors.
        
        (1) matching the original image descriptors with copy of image
            rotated by 90 degrees.
        (2) matching original image descriptors with copy of image 
            rotated by 90 and then scaled down to 40 percent.
        (3) matching original image descriptors with copy of image
            rotated gy 90 and then scaled up to 140 percent.
        */
        
        // ----- test (1) ----------
        Image img = getColorRectangles();
        int w = img.getWidth();
        int h = img.getHeight();
        
        Transformer tr = new Transformer();
        TransformationParameters params = new TransformationParameters();
        params.setOriginX(w/2);
        params.setOriginY(w/2);
        params.setRotationInDegrees(90);        
        Image img90 = tr.applyTransformation(img, params, w, h);
        
        MiscDebug.writeImage(img, "_rect_");
        MiscDebug.writeImage(img90, "_rect_rotated_");
        
        int nPyramidImages = 3;//ORB.estimateNumberOfDefaultScales(w, h);
        int np = 8 * nPyramidImages;
        
        
        PairIntArray expected0 = new PairIntArray(8);
        expected0.add(20, 30);
        expected0.add(50, 30);
        expected0.add(20, 50);
        expected0.add(50, 50);
        expected0.add(90, 70);
        expected0.add(110, 70);
        expected0.add(90, 110);
        expected0.add(110, 110);
        PairIntArray expected90 = tr.applyTransformation(params, 
            expected0);
        
        ORB orb = new ORB(np);
        orb.overrideToUseSmallestPyramid();
        orb.overrideToCreateHSVDescriptors();
        orb.detectAndExtract(img);
        
        List<PairInt> kp0 = orb.getAllKeyPoints();
        List<TIntList> yList = orb.getKeyPoint0List();
        List<TIntList> xList = orb.getKeyPoint1List();
        List<TDoubleList> orList = orb.getOrientationsList();
        List<TFloatList> scalesList = orb.getScalesList();
        Descriptors dH0 = orb.getDescriptorsH().get(0);
        Descriptors dS0 = orb.getDescriptorsS().get(0);
        Descriptors dV0 = orb.getDescriptorsV().get(0);
        
        assertTrue(nPyramidImages >= yList.size());
        assertTrue(yList.size() == xList.size());
        assertTrue(yList.get(0).size() == dH0.descriptors.length);
        assertTrue(yList.get(0).size() == dS0.descriptors.length);
        assertTrue(yList.get(0).size() == dV0.descriptors.length);
        assertTrue(orList.get(0).size() == dV0.descriptors.length);
        
        ORB orb1 = new ORB(np);
        orb1.overrideToUseSmallestPyramid();
        orb1.overrideToCreateHSVDescriptors();
        orb1.detectAndExtract(img90);
        
        List<PairInt> kp90 = orb1.getAllKeyPoints();
        List<TIntList> yList90 = orb1.getKeyPoint0List();
        List<TIntList> xList90 = orb1.getKeyPoint1List();        
        List<TDoubleList> orList90 = orb1.getOrientationsList();
        Descriptors dH90 = orb1.getDescriptorsH().get(0);
        Descriptors dS90 = orb1.getDescriptorsS().get(0);
        Descriptors dV90 = orb1.getDescriptorsV().get(0);
        
        assertTrue(nPyramidImages >= yList90.size());
        assertTrue(yList90.size() == xList90.size());
        assertTrue(yList90.get(0).size() == dH90.descriptors.length);
        assertTrue(yList90.get(0).size() == dS90.descriptors.length);
        assertTrue(yList90.get(0).size() == dV90.descriptors.length);
        assertTrue(orList90.get(0).size() == dV90.descriptors.length);
        
        // key=expected index, value=orb index
        TIntIntMap p1IndexMap = new TIntIntHashMap();
        TIntIntMap p2IndexMap = new TIntIntHashMap();
        for (int i = 0; i < xList.get(0).size(); ++i) {
            int y = yList.get(0).get(i);
            int x = xList.get(0).get(i);
            for (int j = 0; j < expected0.getN(); ++j) {
                int x1 = expected0.getX(j);
                int y1 = expected0.getY(j);
                if (Math.abs(x - x1) < 2 && Math.abs(y - y1) < 2) {
                    PairInt p = new PairInt(x, y);
                    p1IndexMap.put(j, i);
                    break;
                }
            }
        }
        for (int i = 0; i < xList90.get(0).size(); ++i) {
            int y = yList90.get(0).get(i);
            int x = xList90.get(0).get(i);
            for (int j = 0; j < expected90.getN(); ++j) {
                int x1 = expected90.getX(j);
                int y1 = expected90.getY(j);
                if (Math.abs(x - x1) < 2 && Math.abs(y - y1) < 2) {
                    PairInt p = new PairInt(x, y);
                    p2IndexMap.put(j, i);
                    break;
                }
            }
        }
        assertEquals(expected0.getN(), p1IndexMap.size());
        assertEquals(expected90.getN(), p2IndexMap.size());
        assertEquals(p1IndexMap.size(), p2IndexMap.size());
        
        //to make debugging easier, will rearrange data so that index 0
        //is matched to index 0, etc
        int np0_3 = p1IndexMap.size();
        int np90_3 = p2IndexMap.size();
        Descriptors dH0_3 = new Descriptors();
        dH0_3.descriptors = new VeryLongBitString[np0_3];
        Descriptors dS0_3 = new Descriptors();
        dS0_3.descriptors = new VeryLongBitString[np0_3];
        Descriptors dV0_3 = new Descriptors();
        dV0_3.descriptors = new VeryLongBitString[np0_3];
        List<PairInt> kp0_3 = new ArrayList<PairInt>(np0_3);
        
        Descriptors dH90_3 = new Descriptors();
        dH90_3.descriptors = new VeryLongBitString[np90_3];
        Descriptors dS90_3 = new Descriptors();
        dS90_3.descriptors = new VeryLongBitString[np90_3];
        Descriptors dV90_3 = new Descriptors();
        dV90_3.descriptors = new VeryLongBitString[np90_3];
        List<PairInt> kp90_3 = new ArrayList<PairInt>(np90_3);
        
        for (int i = 0; i < p1IndexMap.size(); ++i) {
            int idx1 = p1IndexMap.get(i);
            int idx2 = p2IndexMap.get(i);
            
            dH0_3.descriptors[i] = dH0.descriptors[idx1];
            dS0_3.descriptors[i] = dS0.descriptors[idx1];
            dV0_3.descriptors[i] = dV0.descriptors[idx1];
            kp0_3.add(kp0.get(idx1));
            
            dH90_3.descriptors[i] = dH90.descriptors[idx2];
            dS90_3.descriptors[i] = dS90.descriptors[idx2];
            dV90_3.descriptors[i] = dV90.descriptors[idx2];
            kp90_3.add(kp90.get(idx2));
        }
        
        {
            MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();
            TransformationParameters params0 = tc.calulateEuclidean(
                kp0_3.get(0).getX(), kp0_3.get(0).getY(),
                kp0_3.get(1).getX(), kp0_3.get(1).getY(),
                kp90_3.get(0).getX(), kp90_3.get(0).getY(),
                kp90_3.get(1).getX(), kp90_3.get(1).getY(),
                0, 0);
            assertTrue(Math.abs(params0.getScale() - 1) < 0.001);
        }
        
        TFloatList scales1 = new TFloatArrayList();
        scales1.add(orb.getScalesList().get(0).get(0));
        TFloatList scales2 = new TFloatArrayList();
        scales2.add(orb1.getScalesList().get(0).get(0));
        List<Descriptors> descH1 = new ArrayList<Descriptors>();
        descH1.add(dH0_3);
        List<Descriptors> descS1 = new ArrayList<Descriptors>();
        descS1.add(dS0_3);
        List<Descriptors> descV1 = new ArrayList<Descriptors>();
        descV1.add(dV0_3);
        
        List<Descriptors> descH2 = new ArrayList<Descriptors>();
        descH2.add(dH90_3);
        List<Descriptors> descS2 = new ArrayList<Descriptors>();
        descS2.add(dS90_3);
        List<Descriptors> descV2 = new ArrayList<Descriptors>();
        descV2.add(dV90_3);
        
        List<TIntList> keypointsX1 = new ArrayList<TIntList>();
        List<TIntList> keypointsY1 = new ArrayList<TIntList>();
        keypointsX1.add(new TIntArrayList());
        keypointsY1.add(new TIntArrayList());
        for (int i = 0; i < kp0_3.size(); ++i) {
            int x = kp0_3.get(i).getX();
            int y = kp0_3.get(i).getY();
            keypointsX1.get(0).add(x);
            keypointsY1.get(0).add(y);
        }
        assertEquals(expected0.getN(), kp0_3.size());
        assertEquals(expected90.getN(), kp90_3.size());
        
        List<TIntList> keypointsX2 = new ArrayList<TIntList>();
        List<TIntList> keypointsY2 = new ArrayList<TIntList>();
        keypointsX2.add(new TIntArrayList());
        keypointsY2.add(new TIntArrayList());
        for (int i = 0; i < kp90_3.size(); ++i) {
            int x = kp90_3.get(i).getX();
            int y = kp90_3.get(i).getY();
            keypointsX2.get(0).add(x);
            keypointsY2.get(0).add(y);
        }
        
        CorrespondenceList cor = ORB.matchDescriptors2(
            scales1, scales2,
            descH1, descS1, descV1,
            descH2, descS2, descV2,
            keypointsX1, keypointsY1,
            keypointsX2, keypointsY2, 2);
        
        for (int i = 0; i < cor.getPoints1().size(); ++i) {
            PairInt p1 = cor.getPoints1().get(i);
            PairInt p2 = cor.getPoints2().get(i);
            int idx1 = kp0_3.indexOf(p1);
            int idx2 = kp90_3.indexOf(p2);
            assert(idx1 == idx2);
        }
        assertEquals(expected0.getN(), cor.getPoints1().size());
        
        // ------ use all scale descriptors for multiscale test -----
        for (int i = 1; i < scalesList.size(); ++i) {
            scales1.add(scalesList.get(i).get(0));
            descH1.add(orb.getDescriptorsH().get(i));
            descS1.add(orb.getDescriptorsS().get(i));
            descV1.add(orb.getDescriptorsV().get(i));
            keypointsX1.add(yList.get(i));
            keypointsY1.add(xList.get(i));           
        }
        for (int i = 1; i < orList90.size(); ++i) {
            scales2.add(orb1.getScalesList().get(i).get(0));
            descH2.add(orb1.getDescriptorsH().get(i));
            descS2.add(orb1.getDescriptorsS().get(i));
            descV2.add(orb1.getDescriptorsV().get(i));
            keypointsX2.add(yList90.get(i));
            keypointsY2.add(xList90.get(i));
        }
        
        cor = ORB.matchDescriptors2(
            scales1, scales2,
            descH1, descS1, descV1,
            descH2, descS2, descV2,
            keypointsX1, keypointsY1,
            keypointsX2, keypointsY2, 2);
        
        for (int i = 0; i < cor.getPoints1().size(); ++i) {
            PairInt p1 = cor.getPoints1().get(i);
            PairInt p2 = cor.getPoints2().get(i);
            int idx1 = kp0_3.indexOf(p1);
            int idx2 = kp90_3.indexOf(p2);
            assert(idx1 == idx2);
        }
        assertEquals(expected0.getN(), cor.getPoints1().size());
                   
        // ----- use all scale descriptors for the multiscale test which
        //       returns the top result and those within tolerance of it
        List<CorrespondenceList> corList = ORB.matchDescriptors2(
            scales1, scales2,
            descH1, descS1, descV1,
            descH2, descS2, descV2,
            keypointsX1, keypointsY1,
            keypointsX2, keypointsY2, 2, 0.1f);
        cor = corList.get(0);
        for (int i = 0; i < cor.getPoints1().size(); ++i) {
            PairInt p1 = cor.getPoints1().get(i);
            PairInt p2 = cor.getPoints2().get(i);
            int idx1 = kp0_3.indexOf(p1);
            int idx2 = kp90_3.indexOf(p2);
            assert(idx1 == idx2);
        }
        assertEquals(expected0.getN(), cor.getPoints1().size());         
        
        // ---- scramble the order of points in the dataset2 frame for just
        //      the fist size data and repeat the test.
        //      (this makes sure that the best answer is found even when
        //      order of points does not start w/ true matches)
        {
            Descriptors h2 = descH2.get(0);
            Descriptors s2 = descS2.get(0);
            Descriptors v2 = descV2.get(0);
            TIntList x2 = keypointsX2.get(0);
            TIntList y2 = keypointsY2.get(0);
            // change the order of points.  reverse
            int nSep = x2.size() >> 1;
            for (int idx = 0; idx < nSep; ++idx) {
                int idx2 = x2.size() - 1 - idx;
                
                VeryLongBitString swap1 = h2.descriptors[idx];
                h2.descriptors[idx] = h2.descriptors[idx2];
                h2.descriptors[idx2] = swap1;
                
                swap1 = s2.descriptors[idx];
                s2.descriptors[idx] = s2.descriptors[idx2];
                s2.descriptors[idx2] = swap1;
                
                swap1 = v2.descriptors[idx];
                v2.descriptors[idx] = v2.descriptors[idx2];
                v2.descriptors[idx2] = swap1;
                
                int swap2 = x2.get(idx);
                x2.set(idx, x2.get(idx2));
                x2.set(idx2, swap2);
                
                swap2 = y2.get(idx);
                y2.set(idx, y2.get(idx2));
                y2.set(idx2, swap2);
            }
        }
        corList = ORB.matchDescriptors2(
            scales1, scales2,
            descH1, descS1, descV1,
            descH2, descS2, descV2,
            keypointsX1, keypointsY1,
            keypointsX2, keypointsY2, 2, 0.1f);
        cor = corList.get(0);
        for (int i = 0; i < cor.getPoints1().size(); ++i) {
            PairInt p1 = cor.getPoints1().get(i);
            PairInt p2 = cor.getPoints2().get(i);
            int idx1 = kp0_3.indexOf(p1);
            int idx2 = kp90_3.indexOf(p2);
            assert(idx1 == idx2);
        }
        assertEquals(expected0.getN(), cor.getPoints1().size()); 
        
        // (2) ------- copy the rotated image and scale it by 0.7 -------
        ImageProcessor imageProcessor = new ImageProcessor();
        int w2 = Math.round(0.7f * img90.getWidth());
        int h2 = Math.round(0.7f * img90.getHeight());
        Image r90s70Img = imageProcessor.bilinearDownSampling(
            img90, w2, h2, 0, 255);
        
        TransformationParameters paramsr90s70 =
            new TransformationParameters();
        paramsr90s70.setOriginX(0);
        paramsr90s70.setOriginY(0);
        paramsr90s70.setScale(0.7f);
        
        PairIntArray expectedr90s70 = tr.applyTransformation(
            paramsr90s70, expected90);
        
        ORB orb7 = new ORB(np);
        orb7.overrideToUseSmallestPyramid();
        orb7.overrideToCreateHSVDescriptors();
        orb7.detectAndExtract(r90s70Img);
        
        List<TIntList> yListr90s70 = orb7.getKeyPoint0List();
        List<TIntList> xListr90s70 = orb7.getKeyPoint1List();        
        List<TDoubleList> orListr90s70 = orb7.getOrientationsList();
        List<Descriptors> dHr90s70 = orb7.getDescriptorsH();
        List<Descriptors> dSr90s70 = orb7.getDescriptorsS();
        List<Descriptors> dVr90s70 = orb7.getDescriptorsV();
        TFloatList scalesListr90s70 = new TFloatArrayList();
        for (int i = 0; i < yListr90s70.size(); ++i) {
            scalesListr90s70.add(orb7.getScalesList().get(i).get(0));
        }
            
        p2IndexMap.clear();
        for (int i = 0; i < xListr90s70.get(0).size(); ++i) {
            int y = yListr90s70.get(0).get(i);
            int x = xListr90s70.get(0).get(i);
            for (int j = 0; j < expectedr90s70.getN(); ++j) {
                int x1 = expectedr90s70.getX(j);
                int y1 = expectedr90s70.getY(j);
                if (Math.abs(x - x1) < 2 && Math.abs(y - y1) < 2) {
                    PairInt p = new PairInt(x, y);
                    p2IndexMap.put(j, i);
                    break;
                }
            }
        }
        assertEquals(expectedr90s70.getN(), p2IndexMap.size());
        assertEquals(p1IndexMap.size(), p2IndexMap.size());
        List<PairInt> ordered2 = new ArrayList<PairInt>();
        for (int i = 0; i < expectedr90s70.getN(); ++i) {
            int idx = p2IndexMap.get(i);
            ordered2.add(
                new PairInt(xListr90s70.get(0).get(idx),
                yListr90s70.get(0).get(idx)));
        }
        
        // order the first scale data in the second set of data for easier debugging
        {
            np90_3 = p2IndexMap.size();
            Descriptors dH90s70_3 = new Descriptors();
            dH90s70_3.descriptors = new VeryLongBitString[np90_3];
            Descriptors dS90s70_3 = new Descriptors();
            dS90s70_3.descriptors = new VeryLongBitString[np90_3];
            Descriptors dV90s70_3 = new Descriptors();
            dV90s70_3.descriptors = new VeryLongBitString[np90_3];
            TIntList yr90s70 = new TIntArrayList();
            TIntList xr90s70 = new TIntArrayList();     
            for (int i = 0; i < p2IndexMap.size(); ++i) {
                int idx2 = p2IndexMap.get(i);
                dH90s70_3.descriptors[i] = dHr90s70.get(0).descriptors[idx2];
                dS90s70_3.descriptors[i] = dSr90s70.get(0).descriptors[idx2];
                dV90s70_3.descriptors[i] = dVr90s70.get(0).descriptors[idx2];
                xr90s70.add(xListr90s70.get(0).get(idx2));
                yr90s70.add(yListr90s70.get(0).get(idx2));
            }
            dHr90s70.set(0, dH90s70_3);
            dSr90s70.set(0, dS90s70_3);
            dVr90s70.set(0, dV90s70_3);
            xListr90s70.set(0, xr90s70);
            yListr90s70.set(0, yr90s70);
            
            MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();
            TransformationParameters params0 = tc.calulateEuclidean(
                kp0_3.get(0).getX(), kp0_3.get(0).getY(),
                kp0_3.get(1).getX(), kp0_3.get(1).getY(),
                xListr90s70.get(0).get(0), yListr90s70.get(0).get(0),
                xListr90s70.get(0).get(1), yListr90s70.get(0).get(1),
                0, 0);
            assertTrue(Math.abs(params0.getScale() - 0.7) < 0.02); 
            assertTrue(Math.abs(params0.getRotationInDegrees()- 90) < 0.2); 
        }
        
        corList = ORB.matchDescriptors2(
            scales1, scalesListr90s70,
            descH1, descS1, descV1,
            dHr90s70, dSr90s70, dVr90s70,
            keypointsX1, keypointsY1,
            xListr90s70, yListr90s70, 2, 
            0.1f);
        cor = corList.get(0);
        for (int i = 0; i < cor.getPoints1().size(); ++i) {
            PairInt p1 = cor.getPoints1().get(i);
            PairInt p2 = cor.getPoints2().get(i);
            int idx1 = kp0_3.indexOf(p1);
            int idx2 = ordered2.indexOf(p2);
            assert(idx1 == idx2);
        }
        assertEquals(expected0.getN(), cor.getPoints1().size()); 
    
        // ----- reverse the order to compare a smaller to larger set of images ----
        corList = ORB.matchDescriptors2(
            scalesListr90s70, scales1,
            dHr90s70, dSr90s70, dVr90s70,
            descH1, descS1, descV1,
            xListr90s70, yListr90s70, 
            keypointsX1, keypointsY1,
            2, 0.1f);
        cor = corList.get(0);
        for (int i = 0; i < cor.getPoints1().size(); ++i) {
            PairInt p1 = cor.getPoints1().get(i);
            PairInt p2 = cor.getPoints2().get(i);
            int idx1 = kp0_3.indexOf(p2);
            int idx2 = ordered2.indexOf(p1);
            assert(idx1 == idx2);
        }
        assertEquals(expected0.getN(), cor.getPoints1().size()); 
    
        
    }
    
    private Image getColorRectangles() {
        
        /*
        image dimension 126 X 126
        
        rectangle of magentish HSB: 0.83,<1,1  RGB:(255,48,255)
           x:20 to 50, y:30 to 50
        
        rectangle of silver    HSB: 0.,0,0.75  RGB:(192,192,192) 
           x:90 to 110, y:70 to 110
        
        */
        int w = 126;
        int h = 126;
        
        Image img = new Image(w, h);
        
        for (int i = 20; i < 50; ++i) {
            for (int j = 30; j < 50; ++j) {
                img.setRGB(i, j, 255,48,255);
            }
        }
        
        for (int i = 90; i < 110; ++i) {
            for (int j = 70; j < 110; ++j) {
                img.setRGB(i, j, 192,192,192);
            }
        }
        
        return img;
    }
    
    private Image getFigureEight() throws IOException {
        
        String fileName = "test_img.png";  
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        
        return img0;
    }
    private Image getFigureEight1() throws IOException {
        
        String fileName = "test_img_1.png";  
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        
        return img0;
    }
    private Image getFigureEight2() throws IOException {
        
        String fileName = "test_img_2.png";  
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img0 = ImageIOHelper.readImageAsGrayScale(filePath);
        
        return img0;
    }
}
