package thirdparty.libquantum;

import java.util.Arrays;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class ShorTest extends TestCase {
    
    public ShorTest() {
    }
    
    public void testRun() {
                
        int number = 15;
        
        final int[] expected = new int[]{3, 5};
        
        int nTests = 100;
        int nCorrect = 0;
        
        for (int i = 0; i < nTests; ++i) {
            
            Shor shor = new Shor(15, 8);
        
            int[] factors = shor.run();
            
            if (factors != null && factors.length == 2) {
                Arrays.sort(factors);
                if (Arrays.equals(expected, factors)) {
                    nCorrect++;
                }
            } 
        }
        
        System.out.println("nTests=" + nTests + " nCorect=" + nCorrect);
        
    }
}
