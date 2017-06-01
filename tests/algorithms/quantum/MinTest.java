package thirdparty.libquantum;

import java.util.logging.Logger;
import junit.framework.TestCase;
import algorithms.quantum.Grover2;

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
    
}
