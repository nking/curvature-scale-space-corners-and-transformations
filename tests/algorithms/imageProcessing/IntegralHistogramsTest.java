package algorithms.imageProcessing;

import algorithms.imageProcessing.features.HOGUtil;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class IntegralHistogramsTest extends TestCase {
    
    public IntegralHistogramsTest() {
    }
    
    public void testIntegralHistograms() {
        
        /*
        column major notation is used
        
        2  155  2   7
        1  3   16  42 
        0  2    5   9
           0    1   2
         */
        GreyscaleImage img = new GreyscaleImage(3, 3);
        img.setValue(0, 0, 2); img.setValue(0, 1, 3); img.setValue(0, 2, 155);
        img.setValue(1, 0, 5); img.setValue(1, 1, 16); img.setValue(1, 2, 2);
        img.setValue(2, 0, 9); img.setValue(2, 1, 42); img.setValue(2, 2, 7);
       
        IntegralHistograms sumTable = new IntegralHistograms();
        int[][] sHists = sumTable.create(img);
        assert(sHists.length == img.getNPixels());
        assert(sHists[0].length == 16);
        int minValue = 0;
        int maxValue = 255;
        int binWidth = 16;
        
        /*
        column major notation is used
        2  155  2   7      6  7  8
        1  3   16  42      3  4  5
        0  2    5   9      0  1  2
           0    1   2
        pixIdx = (row * width) + col
        0,0 = 0  1,0 = 1  2,0 = 2
        0,1 = 3  1,1 = 4  2,1 = 5
        0,2 = 6  1,2 = 7  2,2 = 8
        */
        int[][] expected = new int[img.getNPixels()][16];
        expected[0] = new int[]{1, 0, 0, 0,
                                0, 0, 0, 0,
                                0, 0, 0, 0,
                                0, 0, 0, 0};
        expected[1] = new int[]{2, 0, 0, 0,
                                0, 0, 0, 0,
                                0, 0, 0, 0,
                                0, 0, 0, 0};
        expected[2] = new int[]{3, 0, 0, 0,
                                0, 0, 0, 0,
                                0, 0, 0, 0,
                                0, 0, 0, 0};
        expected[3] = new int[]{2, 0, 0, 0,
                                0, 0, 0, 0,
                                0, 0, 0, 0,
                                0, 0, 0, 0};
        expected[4] = new int[]{3, 1, 0, 0,    
                                0, 0, 0, 0,
                                0, 0, 0, 0,
                                0, 0, 0, 0};
        expected[5] = new int[]{4, 1, 1, 0,
                                0, 0, 0, 0,
                                0, 0, 0, 0,
                                0, 0, 0, 0};
        expected[6] = new int[]{2, 0, 0, 0,
                                0, 0, 0, 0,
                                0, 1, 0, 0,
                                0, 0, 0, 0};
        expected[7] = new int[]{4, 1, 0, 0,
                                0, 0, 0, 0,
                                0, 1, 0, 0,
                                0, 0, 0, 0};
        expected[8] = new int[]{6, 1, 1, 0,
                                0, 0, 0, 0,
                                0, 1, 0, 0,
                                0, 0, 0, 0};
       
        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                int pixIdx = img.getInternalIndex(i, j);
                assertTrue(Arrays.equals(expected[pixIdx], 
                    sHists[pixIdx]));
            }
        }
        
        /*
        column major notation is used
        2  155  2   7      6  7  8
        1  3   16  42      3  4  5
        0  2    5   9      0  1  2
           0    1   2
        pixIdx = (row * width) + col
        0,0 = 0  1,0 = 1  2,0 = 2
        0,1 = 3  1,1 = 4  2,1 = 5
        0,2 = 6  1,2 = 7  2,2 = 8
        */
       
        int expectedSum;
        int[] outputN = new int[1];
        int[] output = new int[16];
        int w = img.getWidth();
        int h = img.getHeight();
        int[] expectedWindow = new int[16];
        
        // window=5
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 0, 0, 5, output, outputN);
        HOGUtil.add(expectedWindow, sHists[8]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(9, outputN[0]);
        
        // window=3
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 0, 0, 3, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[4]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(4, outputN[0]);
        
        // window=1
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 0, 0, 1, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[0]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(1, outputN[0]);
        
        /*
        column major notation is used
        2  155  2   7      6  7  8
        1  3   16  42      3  4  5
        0  2    5   9      0  1  2
           0    1   2
        pixIdx = (row * width) + col
        0,0 = 0  1,0 = 1  2,0 = 2
        0,1 = 3  1,1 = 4  2,1 = 5
        0,2 = 6  1,2 = 7  2,2 = 8
        */
        
        // window=5
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 0, 1, 5, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[8]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(9, outputN[0]);
        
        // window=3
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 0, 1, 3, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[7]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(6, outputN[0]);
        
        // window=1
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 0, 1, 1, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[3]);
        HOGUtil.subtract(expectedWindow, sHists[0]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(1, outputN[0]);
        
        /*
        column major notation is used
        2  155  2   7      6  7  8
        1  3   16  42      3  4  5
        0  2    5   9      0  1  2
           0    1   2
        pixIdx = (row * width) + col
        0,0 = 0  1,0 = 1  2,0 = 2
        0,1 = 3  1,1 = 4  2,1 = 5
        0,2 = 6  1,2 = 7  2,2 = 8
        */
        
        // window=5
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 0, 2, 5, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[8]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(9, outputN[0]);
       
        // window=3
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 0, 2, 3, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[7]);
        HOGUtil.subtract(expectedWindow, sHists[1]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(4, outputN[0]);
        
        // window=1
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 0, 2, 1, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[6]);
        HOGUtil.subtract(expectedWindow, sHists[3]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(1, outputN[0]);
        
        /*
        column major notation is used
        2  155  2   7      6  7  8
        1  3   16  42      3  4  5
        0  2    5   9      0  1  2
           0    1   2
        pixIdx = (row * width) + col
        0,0 = 0  1,0 = 1  2,0 = 2
        0,1 = 3  1,1 = 4  2,1 = 5
        0,2 = 6  1,2 = 7  2,2 = 8
        */
        
        // window=5
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 1, 0, 5, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[8]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(9, outputN[0]);
        
        // window=3
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 1, 0, 3, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[5]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(6, outputN[0]);
        
        // window=1
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 1, 0, 1, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[1]);
        HOGUtil.subtract(expectedWindow, sHists[0]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(1, outputN[0]);
        
        /*
        column major notation is used
        2  155  2   7      6  7  8
        1  3   16  42      3  4  5
        0  2    5   9      0  1  2
           0    1   2
        pixIdx = (row * width) + col
        0,0 = 0  1,0 = 1  2,0 = 2
        0,1 = 3  1,1 = 4  2,1 = 5
        0,2 = 6  1,2 = 7  2,2 = 8
        */
        
        // window=5
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 1, 1, 5, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[8]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(9, outputN[0]);
        
        // window=3
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 1, 1, 3, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[8]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(9, outputN[0]);
        
        // window=1
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 1, 1, 1, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[4]);
        HOGUtil.subtract(expectedWindow, sHists[3]);
        HOGUtil.subtract(expectedWindow, sHists[1]);
        HOGUtil.add(expectedWindow, sHists[0]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(1, outputN[0]);
        
        /*
        column major notation is used
        2  155  2   7      6  7  8
        1  3   16  42      3  4  5
        0  2    5   9      0  1  2
           0    1   2
        pixIdx = (row * width) + col
        0,0 = 0  1,0 = 1  2,0 = 2
        0,1 = 3  1,1 = 4  2,1 = 5
        0,2 = 6  1,2 = 7  2,2 = 8
        */
        
        // window=5
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 1, 2, 5, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[8]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(9, outputN[0]);
        
        // window=3
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 1, 2, 3, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[8]);
        HOGUtil.subtract(expectedWindow, sHists[2]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(6, outputN[0]);
        
        // window=1
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 1, 2, 1, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[7]);
        HOGUtil.subtract(expectedWindow, sHists[6]);
        HOGUtil.subtract(expectedWindow, sHists[4]);
        HOGUtil.add(expectedWindow, sHists[3]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(1, outputN[0]);
        
        /*
        column major notation is used
        2  155  2   7      6  7  8
        1  3   16  42      3  4  5
        0  2    5   9      0  1  2
           0    1   2
        pixIdx = (row * width) + col
        0,0 = 0  1,0 = 1  2,0 = 2
        0,1 = 3  1,1 = 4  2,1 = 5
        0,2 = 6  1,2 = 7  2,2 = 8
        */
        // window=5
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 2, 0, 5, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[8]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(9, outputN[0]);
        
        // window=3
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 2, 0, 3, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[5]);
        HOGUtil.subtract(expectedWindow, sHists[3]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(4, outputN[0]);
        
        // window=1
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 2, 0, 1, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[2]);
        HOGUtil.subtract(expectedWindow, sHists[1]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(1, outputN[0]);
        
        /*
        column major notation is used
        2  155  2   7      6  7  8
        1  3   16  42      3  4  5
        0  2    5   9      0  1  2
           0    1   2
        pixIdx = (row * width) + col
        0,0 = 0  1,0 = 1  2,0 = 2
        0,1 = 3  1,1 = 4  2,1 = 5
        0,2 = 6  1,2 = 7  2,2 = 8
        */
        // window=5
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 2, 1, 5, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[8]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(9, outputN[0]);
        
        // window=3
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 2, 1, 3, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[8]);
        HOGUtil.subtract(expectedWindow, sHists[6]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(6, outputN[0]);
        
        // window=1
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 2, 1, 1, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[5]);
        HOGUtil.subtract(expectedWindow, sHists[4]);
        HOGUtil.subtract(expectedWindow, sHists[2]);
        HOGUtil.add(expectedWindow, sHists[1]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(1, outputN[0]);
        
        /*
        column major notation is used
        2  155  2   7      6  7  8
        1  3   16  42      3  4  5
        0  2    5   9      0  1  2
           0    1   2
        pixIdx = (row * width) + col
        0,0 = 0  1,0 = 1  2,0 = 2
        0,1 = 3  1,1 = 4  2,1 = 5
        0,2 = 6  1,2 = 7  2,2 = 8
        */
        // window=5
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 2, 2, 5, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[8]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(9, outputN[0]);
        
        // window=3
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 2, 2, 3, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[8]);
        HOGUtil.subtract(expectedWindow, sHists[6]);
        HOGUtil.subtract(expectedWindow, sHists[2]);
        HOGUtil.add(expectedWindow, sHists[0]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(4, outputN[0]);
        
        // window=1
        sumTable.extractWindowFromIntegralHistograms(sHists, 
            w, h, 2, 2, 1, output, outputN);
        Arrays.fill(expectedWindow, 0);
        HOGUtil.add(expectedWindow, sHists[8]);
        HOGUtil.subtract(expectedWindow, sHists[7]);
        HOGUtil.subtract(expectedWindow, sHists[5]);
        HOGUtil.add(expectedWindow, sHists[4]);
        assertTrue(Arrays.equals(expectedWindow, output));
        assertEquals(1, outputN[0]);
        
    }
}
