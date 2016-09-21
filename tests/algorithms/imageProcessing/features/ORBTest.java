package algorithms.imageProcessing.features;

import algorithms.util.PairInt;
import gnu.trove.list.TIntList;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ORBTest extends TestCase {
    
    public ORBTest() {
    }
    
    public void testPeakLocalMax() {
        
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
        
        TIntList peakRowCols = orb.peakLocalMax(img, 1, 0.1f);
        
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
        
        for (int i = 0; i < peakRowCols.size(); i += 2) {
            int x = peakRowCols.get(i);
            int y = peakRowCols.get(i + 1);
            PairInt p = new PairInt(x, y);
            assertTrue(expected.remove(p));
        }
    }
    
    public void testCornerPeaks() {
        
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
        
        
        TIntList coords = orb.cornerPeaks(cf, 1);
        
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
        
        for (int i = 0; i < coords.size(); i += 2) {
            int x = coords.get(i);
            int y = coords.get(i + 1);
            PairInt p = new PairInt(x, y);
            assertTrue(expected.remove(p));
        }         
    }
}
