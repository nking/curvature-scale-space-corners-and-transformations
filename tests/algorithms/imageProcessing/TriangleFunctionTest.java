package algorithms.imageProcessing;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TriangleFunctionTest extends TestCase {
    
    public TriangleFunctionTest() {
    }
    
    public void testCalculate() {
        
        int d = 10;
        
        int n = d*d;
        
        GreyscaleImage input = new GreyscaleImage(d, d, 
            GreyscaleImage.Type.Bits32FullRangeInt);
        
        for (int x = 0; x < d; ++x) {
            for (int y = 0; y < d; ++y) {
                int pixIdx = (y * d) + x;
                input.setValue(pixIdx, y * 256);
            }
        }
        
        /*
        input:
        0  0  0  0  0  0  0  0  0  0 
        1  1  1  1  1  1  1  1  1  1
        2  2  2  2  2  2  2  2  2  2
        3  3  3  3 *3* 3  3  3  3  3
        4  4  4  4  4  4  4  4  4  4
        5  5  5  5  5  5  5  5  5  5
        ...
         // add:
        //    c_(j + 1,k) = (1/4)*c_(j,k-(2^j)) + (1/2)*c_(j,k) + (1/4)*c_(j,k+(2^j))
        // subtract:
        //    w_(j+1,k) = (-1/4)*c_(j,k-(2^j)) + (1/2)*c_(j,k) - (1/4)*c_(j,k+(2^j))
        */ 
        
        int s = 1;
        
        TriangleFunction tf = new TriangleFunction();
        GreyscaleImage inputSubtr = tf.subtractLevels(input, s);
                
        GreyscaleImage inputAdd = tf.calculateNextLevel(input, s);
        
        int expectedAdd = 3 * 256;
        int expectedSubtr = 0;
         
        assertEquals(expectedSubtr, inputSubtr.getValue(3, 3));
        
        assertEquals(expectedAdd, inputAdd.getValue(3, 3));
        
    }
   
}
