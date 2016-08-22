package algorithms.imageProcessing.matching;

import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SequenceTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass()
        .getName());
    
    public SequenceTest() {
    }

    public void estIntersects() {
        
        Sequence s1, s2;
        int n1 = 187;
        int n2 = 231;
        int offset;
        
        // 7:19 to 21
        // 11:23 to 29
        offset = 12;
        s1 = new Sequence(n1, n2, offset);
        s1.startIdx1 = 7;
        s1.startIdx2 = 19;
        s1.stopIdx2 = 21;
        assertEquals(12, s1.getOffset());
        
        s2 = new Sequence(n1, n2, offset);
        s2.startIdx1 = 11;
        s2.startIdx2 = 23;
        s2.stopIdx2 = 29;
        assertEquals(12, s2.getOffset());
        
        assertFalse(s1.intersects(s2));
        
        
        // 17:54 to 62  n1=187  n2=231
        // 15:54 to 64
        
        s1 = new Sequence(n1, n2, 54-17);
        s1.startIdx1 = 17;
        s1.startIdx2 = 54;
        s1.stopIdx2 = 62;
        assertEquals(54 - 17, s1.getOffset());
        
        s2 = new Sequence(n1, n2, 54-15);
        s2.startIdx1 = 15;
        s2.startIdx2 = 54;
        s2.stopIdx2 = 64;
        assertEquals(54 - 15, s2.getOffset());
        
        assertTrue(s1.intersects(s2));
        
        //43:0 to 8
        //44:0 to 7
        s1 = new Sequence(n1, n2, -43);
        s1.startIdx1 = 43;
        s1.startIdx2 = 0;
        s1.stopIdx2 = 8;
        assertEquals(0 - 43, s1.getOffset());
        
        s2 = new Sequence(n1, n2, -44);
        s2.startIdx1 = 44;
        s2.startIdx2 = 0;
        s2.stopIdx2 = 9;
        assertEquals(0 - 44, s2.getOffset());
        
        assertTrue(s1.intersects(s2));
        
    }
    
    public void estMerge() {
        
        Sequence s1, s2;
        Sequence[] merged;
        int n1 = 187;
        int n2 = 231;
        
        // 2:11  34
        // 13:22 67
        s1 = new Sequence(n1, n2, 9);
        s1.startIdx1 = 2;
        s1.startIdx2 = 11;
        s1.stopIdx2 = 34;
        
        s2 = new Sequence(n1, n2, 9);
        s2.startIdx1 = 13;
        s2.startIdx2 = 22;
        s2.stopIdx2 = 67;
        
        //43:0 to 8
        //44:0 to 7
        s1 = new Sequence(n1, n2, -43);
        s1.startIdx1 = 43;
        s1.startIdx2 = 0;
        s1.stopIdx2 = 8;
        
        s2 = new Sequence(n1, n2, -44);
        s2.startIdx1 = 44;
        s2.startIdx2 = 0;
        s2.stopIdx2 = 9;

    }
 
}
