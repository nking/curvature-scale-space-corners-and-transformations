package thirdparty.libquantum;

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
         
        Grover grover = new Grover();
            
        int[] list = new int[]{0,1,2,3,4,5,6,7};
        
        //int[] list = new int[]{1,2,3,0,1,4,2,7};
        
        grover.run(number, nbits, list);
          
    }
    
    
}
