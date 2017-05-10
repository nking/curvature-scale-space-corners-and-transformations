package algorithms.quantum;

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
         
        /*
 from "Experimental demonstration of Shor’s algorithm with quantum entanglement:
Consider N=15: if we choose a=2, the quantum routine
finds r=4, and the prime factors are given by the nontrivial
greatest common divisor of C
r/2±1 and N, i.e. 3
and 5; similarly if we choose the next possible co-prime,
C=4, we find the order r=2, yielding the same factors
        */
        
        int number = 15;
        
        Shor shor = new Shor(number, 1494367752772L);
        
        int[] expResult = new int[]{3, 5};
        int[] result = shor.run();
        
        //Arrays.sort(result);
        
        //assertArrayEquals(expResult, result);
    }

    public void estRandomGCD() {
        
        int number = 15;
        
        Shor shor = new Shor(number, 1494367040443L);
        
        int[] expResult = new int[]{3, 5};
        int[] result = shor.randomGCD(number);
        Arrays.sort(result);
        
        assertArrayEquals(expResult, result);        
    }
    
}
