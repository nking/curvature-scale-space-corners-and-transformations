package algorithms.quantum;

import algorithms.misc.Misc;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Random;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MinTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public MinTest() {
    }
    
    public void testHashCode() {
        
        int length = 100000;//256;
        
        int n = 1000000;
        
        TIntObjectMap<TIntSet> idxVMap = new TIntObjectHashMap<TIntSet>(length);
        
        Min min = new Min();
        
        Random rng = Misc.getSecureRandom();
        
        int nCollisions = 0;
        
        for (int i = 0; i < n; ++i) {
            int v = rng.nextInt(length);
            int idx = min.hashToIndex(v, length);
            assertTrue(idx < length);
            assertTrue(idx >= 0);
            if (idxVMap.containsKey(idx) && idxVMap.get(idx).contains(v)) {
                nCollisions++;
            }
            TIntSet set = idxVMap.get(idx);
            if (set == null) {
                set = new TIntHashSet();
                idxVMap.put(idx, set);
            }
            set.add(v);
        }
                
        System.out.println("nCollisons=" + nCollisions 
            + " idxs.size=" + idxVMap.size() + " nC/N=" 
            + ((double)nCollisions/n));
        
        assertTrue(idxVMap.size() > 0.9*length);
    }
    
    public void testRun() {
        
        log.info("testRun");
                       
        int nbits = 3;
         
        Min min = new Min();
            
        //int[] list = new int[]{0,1,2,3,4,5,6,7};
        
        int[] list = new int[]{1,2,3,0,1,4,2,7};
        
        int m = min.run(nbits, list);        
    
        System.out.println("min=" + m);
        
        assertEquals(3, m);
    }
    
}
