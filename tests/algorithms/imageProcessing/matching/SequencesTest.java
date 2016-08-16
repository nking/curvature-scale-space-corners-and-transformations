package algorithms.imageProcessing.matching;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SequencesTest extends TestCase {
    
    public SequencesTest() {
    }
    
    public void testIsConsistentClockwise() {
        
        Sequences sequences = new Sequences();
        
        int n1 = 187;
        int n2 = 283;
        Sequence s;
        
        s = new Sequence(n1, n2, 0);
        s.startIdx1 = 0;
        s.startIdx2 = 0;
        s.stopIdx2 = 7;
        sequences.add(s);
        
        s = new Sequence(n1, n2, 1);
        s.startIdx1 = 0;
        s.startIdx2 = 1;
        s.stopIdx2 = 10;
        sequences.add(s);
        
        assertFalse(sequences.isConsistentClockwise());
        assertFalse(Sequences.isConsistentClockwise(
            sequences.getSequences()));
        
        sequences.getSequences().clear();
        s = new Sequence(n1, n2, 0);
        s.startIdx1 = 26;
        s.startIdx2 = 26;
        s.stopIdx2 = 33;
        sequences.add(s);
        
        s = new Sequence(n1, n2, 1);
        s.startIdx1 = 27;
        s.startIdx2 = 28;
        s.stopIdx2 = 31;
        sequences.add(s);
        
        assertFalse(sequences.isConsistentClockwise());
        assertFalse(Sequences.isConsistentClockwise(
            sequences.getSequences()));
        
        // ---------
        sequences.getSequences().clear();
        s = new Sequence(n1, n2, 1);
        s.startIdx1 = 33;
        s.startIdx2 = 34;
        s.stopIdx2 = 38;
        sequences.add(s);
        
        s = new Sequence(n1, n2, 1);
        s.startIdx1 = 46;
        s.startIdx2 = 47;
        s.stopIdx2 = 51;
        sequences.add(s);
        
        assertTrue(sequences.isConsistentClockwise());
        assertTrue(Sequences.isConsistentClockwise(
            sequences.getSequences()));
        
       // --------------------
       
    }
}
