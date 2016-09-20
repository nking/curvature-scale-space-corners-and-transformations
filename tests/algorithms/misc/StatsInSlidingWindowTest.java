package algorithms.misc;

import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class StatsInSlidingWindowTest extends TestCase {
    
    public StatsInSlidingWindowTest() {
    }
    
    public void testCalculateMaximum() {
        
        float[][] input = new float[5][5];
        float[][] output = new float[5][5];
        float[][] expected = new float[5][5];
        for (int i = 0; i < input.length; ++i) {
            input[i] = new float[5];
            output[i] = new float[5];
            expected[i] = new float[5];
            for (int j = 0; j < input[i].length; ++j) {
                input[i][j] = 4 - j + i;
            }
        }
        
        /*
        4 3 2 1 0 
        5 4 3 2 1
        6 5 4 3 2
        7 6 5 4 3
        8 7 6 5 4
        */
        expected[0] = new float[]{6, 5, 4, 3, 2};
        expected[1] = new float[]{7, 6, 5, 4, 3};
        expected[2] = new float[]{8, 7, 6, 5, 4};
        expected[3] = new float[]{8, 7, 6, 5, 4};
        expected[4] = new float[]{8, 7, 6, 5, 4};
        
        StatsInSlidingWindow sw = new StatsInSlidingWindow();
        
        sw.calculateMaximum(input, output, 3, 3);
        
        for (int i = 0; i < input.length; ++i) {
            assertTrue(Arrays.equals(expected[i], output[i]));
        }
    }
    
    public void testCalculateMaximum2() {
        
        float[][] input = new float[6][6];
        float[][] output = new float[6][6];
        float[][] expected = new float[6][6];
        for (int i = 0; i < input.length; ++i) {
            input[i] = new float[6];
            output[i] = new float[6];
            expected[i] = new float[6];
            for (int j = 0; j < input[i].length; ++j) {
                input[i][j] = 5 - j + i;
            }
        }
        
        int xWindow = 4;
        int yWindow = 4;
        
        /*
        5  4 3 2 1 0 
        6  5 4 3 2 1
        7  6 5 4 3 2
        8  7 6 5 4 3
        9  8 7 *6 5 4
        10 9 8 7 6 5
        */
        expected[0] = new float[]{8, 7, 6, 5, 4, 3};
        expected[1] = new float[]{9, 8, 7, 6, 5, 4};
        expected[2] = new float[]{10, 9, 8, 7, 6, 5};
        expected[3] = new float[]{10, 9, 8, 7, 6, 5};
        expected[4] = new float[]{10, 9, 8, 7, 6, 5};
        expected[5] = new float[]{10, 9, 8, 7, 6, 5};
        
        StatsInSlidingWindow sw = new StatsInSlidingWindow();
        
        sw.calculateMaximum(input, output, xWindow, yWindow);
        
        for (int i = 0; i < input.length; ++i) {
            assertTrue(Arrays.equals(expected[i], output[i]));
        }
        
        //------
    }
}
