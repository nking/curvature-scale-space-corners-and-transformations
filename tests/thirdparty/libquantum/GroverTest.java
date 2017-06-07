package thirdparty.libquantum;

import algorithms.misc.MiscMath;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class GroverTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public GroverTest() {
    }
    
    public void testRun() {
        
        log.info("testRun");
               
        int number = 3;
        int nbits = 3;
                
        int nTests = 1;
        
        for (int i = 0; i < nTests; ++i) {
            
            Grover grover = new Grover();
            
            grover.run(number, nbits);
        }        
    }
    
    public void testRun_list() {
        
        log.info("testRun_list");
               
        int number = 3;
        int nbits = 3;
        //int[] list = new int[]{0,1,2,3,4,5,6,7};
        //int[] list = new int[]{1,2,3,0,1,4,2,7};
         
        Grover grover = new Grover();
            
        // (2*nBits) + 4
        nbits = 14;
        int[] list = new int[]{1,2,3,4,11,14,12,17};
        
        int r = grover.run(number, nbits, list);
        System.out.println("search=" + number + " found=" + r);
        
        number = 5;
        r = grover.run(number, nbits, list);
        System.out.println("search=" + number + " found=" + r);
    }
    
    public void testSetBits() {
        
        System.out.println("testSetBits");
        
        // 1 1 0 1
        //int setBits = (1 << 0) + (1 << 2) + (1 << 3);
        
        // 0 1 1 0
        //int setBits = (1 << 1) + (1 << 2);
        
        // 1 1 1 1
        int setBits = (1 << 0) + (1 << 1) + (1 << 2) + (1 << 3);
        
        /*
        qubits
        0000
    [junit] 0001
    [junit] 0100
    [junit] 0101
    [junit] 1000
    [junit] 1001
    [junit] 1100
    [junit] 1101
        */
        
        System.out.println("setBits=" + setBits);
        TIntSet set = new TIntHashSet();
        int nBits = MiscMath.numberOfBits(setBits);
        for (int i = 0; i < nBits; ++i) {
            if ((setBits & (1 << i)) != 0) {
                set.add(1 << i);
            }
        }
        
        // (2*nBits) + 1
        int width = 9;
        for (int i = 0; i < 10;++i) {
            int number = i;
        
            System.out.println("searching for " + number);
                 
            Grover grover = new Grover();
                
            int r = grover.run(number, width, setBits);
            
            System.out.println("number=" + number + " set=" +
                (set.contains(number)) + " r=" + r);
            if (set.contains(number)) {
                assertEquals(number, r);
            } else {
                assertFalse(r == number);    
            }
        }
        
    }
}
