package thirdparty.libquantum;

import algorithms.misc.NumberTheory;
import java.util.Arrays;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ShorTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public ShorTest() {
    }
    
    public void testRun() {
        
        log.info("testRun");
        
        /*
         from "Experimental demonstration of Shor’s algorithm with 
        quantum entanglement":
        Consider N=15: if we choose a=2, the quantum routine
        finds r=4, and the prime factors are given by the nontrivial
        greatest common divisor of C
        r/2±1 and N, i.e. 3
        and 5; similarly if we choose the next possible co-prime,
        C=4, we find the order r=2, yielding the same factors
        */
        
                
        int number = 15;
        
        final int[] expected = new int[]{3, 5};
        
        int nTests = 100;
        int nCorrect = 0;
        
        int maxNRetries = 0;
        int nRetries = 0;
        
        for (int i = 0; i < nTests; ++i) {
            
            Shor shor = new Shor(15, 8);
        
            int[] factors = shor.run();
            
            if (factors != null && factors.length == 2) {
                Arrays.sort(factors);
                if (Arrays.equals(expected, factors)) {
                    nCorrect++;
                    if (nRetries > maxNRetries) {
                        maxNRetries = nRetries;
                    }
                    nRetries = 0;
                } else {
                    nRetries++;
                }
            } 
        }
        
        log.info("nTests=" + nTests + " nCorrect=" + nCorrect +
            " max number of retries between fails=" + maxNRetries);
    }
    
    public void testEuclidModularMethods() {
        
        log.info("testEuclidModularMethods");
        
        long[] factors;
                                      //14 * x = 30 % 100
                                      //       = 30
        //GreatestCommonDenominator.gcdModularLinearEqnSolver(14, 30, 100);
        //  given (a, b, n) solves for x in the equation a * x ≡ b (mod n) (which is actually (a*x) mod n = b)
        // where d is the gcd of number n and d|b (that is, d divides b).
        // finds the smallest gcd for which a*x + b*y = d.
        // The equation may have zero, one, or more than one such solution.
        long[] mcd = NumberTheory.gcdModularLinearEqnSolver(78, 99, 99); // (b,y): (-77,3), (-33,7)
        //System.out.println("mcd = " + Arrays.toString(mcd)); //[66, 0, 33]
        //assertEquals(mcd[0], 3L);
        
        factors = NumberTheory.extendedEuclid(78, 99);
        //System.out.println("ext euc factors = " + Arrays.toString(factors));
        assertEquals(factors[0], 3);
        assertEquals(factors[1], 14);
        
    }
    
    public void testRun2() {
        
        log.info("testRun2");
        
        int[] factors;
        
        int nTests = 10;
        
        for (int i = 0; i < nTests; ++i) {
            
            Shor shor = new Shor(99, 14);//32768
            //Shor shor = new Shor(12345); //cofactors 5, 2469
            //Shor shor = new Shor(2469);//cofactors 3, 823
            
            factors = shor.run();

            log.info("*factors = " + Arrays.toString(factors));
        }
    }
    
}
