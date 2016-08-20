package algorithms.imageProcessing.matching;

import java.io.Console;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SequenceTest extends TestCase {
    
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
        
        merged = s1.merge(s2);
        assertNotNull(merged);
        assertEquals(1, merged.length);
        assertEquals(2, merged[0].startIdx1);
        assertEquals(11, merged[0].startIdx2);
        assertEquals(67, merged[0].stopIdx2);
        merged = s2.merge(s1);
        assertNotNull(merged);
        assertEquals(1, merged.length);
        assertEquals(2, merged[0].startIdx1);
        assertEquals(11, merged[0].startIdx2);
        assertEquals(67, merged[0].stopIdx2);

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

        merged = s1.merge(s2);
        assertNull(merged);
        merged = s2.merge(s1);
        assertNull(merged);
    }
    
    public void estMerge2() throws Exception {
        
        /* trying to test all of the merge branch logic
           by making many valid random sequences within
           same n1, n2 space and merging them.
           eventually, should have 3 or 4 sequences, 
           two of which are the start and stop sentinel 
           sequences.
                          idx1   idx1   idx2   idx2
        idx1: 0 1 2 3      0      0      3      3
        idx2: 3 0 1 2      1      1      4      4
                           1      1      0      0
                           2      3      1      2
        */
        
        SecureRandom sr = SecureRandom.getInstance(
            "SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1471293393476L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        System.out.flush();
        
        int n1 = 50;
        int n2 = n1;
        int nSeq = 10 * n2;
        
        int nTests = 1;
        for (int nTest = 0; nTest < nTests; ++nTest) {
            
            int offset12 = sr.nextInt((n1 - 2)/2);
            
            List<Sequence> seqs = 
                new ArrayList<Sequence>();
            
            for (int i = 0; i < nSeq; ++i) {
                Sequence s0 = createSequence(n1, n2, sr, offset12);
                seqs.add(s0);
            }
            
            Sequence.mergeSequences(seqs);
            
            System.out.println("final nMerged=" + seqs.size());
            
            for (Sequence s : seqs) {
                System.out.println("SEQ " + s);
            }
              
            assertTrue(seqs.size() <= 4);
            boolean found0 = false;
            boolean found1 = false;
            for (Sequence seq : seqs) {
                if (seq.isStartSentinel()) {
                    found0 = true;
                }
                if (seq.isStopSentinel()) {
                    found1 = true;
                }
            }
            assertTrue(found0);
            assertTrue(found1);          
        }
        
        // TODO: test for larger n1, n2 with mixed offsets
        
    }

    public void testMerge3() throws Exception {
        
        /* trying to test all of the merge branch logic
           by making many valid random sequences within
           same n1, n2 space and merging them.
           eventually, should have 3 or 4 sequences, 
           two of which are the start and stop sentinel 
           sequences.
                          idx1   idx1   idx2   idx2
        idx1: 0 1 2 3      0      0      3      3
        idx2: 3 0 1 2      1      1      4      4
                           1      1      0      0
                           2      3      1      2
        */
        
        SecureRandom sr = SecureRandom.getInstance(
            "SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1471293393476L;
        sr.setSeed(seed);
        System.out.println("SEED3=" + seed);
        System.out.flush();
        
        int n1 = 50;
        int n2 = n1;
        int nSeq = 10 * n2;
        
        int nTests = 1;
        for (int nTest = 0; nTest < nTests; ++nTest) {
            
            int[] offsets = new int[4];
            for (int j = 0; j < offsets.length; ++j) {
                offsets[j] = sr.nextInt((n1 - 2)/2);
            }
            
            List<Sequence> seqs = 
                new ArrayList<Sequence>();
            
            int j = 0;
            for (int i = 0; i < nSeq; ++i) {
                Sequence s0 = createSequence(n1, n2, sr, offsets[j]);
                seqs.add(s0);
                j++;
                if (j > (offsets.length - 1)) {
                    j = 0;
                }
            }
            
            LinkedHashSet<Sequence> seqs2 = 
                new LinkedHashSet<Sequence>();
            
            int nIter = 0;
            
            boolean didMerge = false;
            do {
                
                if (seqs.size() == 1) {
                    // break because loop below won't store it
                    break;
                }
                
                // sort seqs
                Collections.sort(seqs, new SequenceComparator4());
                
                if (seqs.size() <= 3) {
                    boolean found0 = false;
                    boolean found1 = false;
                    for (Sequence seq : seqs) {
                        if (seq.isStartSentinel()) {
                            found0 = true;
                        }
                        if (seq.isStopSentinel()) {
                            found1 = true;
                        }
                    }
                    if (found0 && found1) {
                        break;
                    }
                }
                
                System.out.println("seqs.size=" + seqs.size());
                for (Sequence s : seqs) {
                    System.out.println(nIter + ": SEQ " + s);
                }
                
                didMerge = false;
                Iterator<Sequence> iter = seqs.iterator();
                Sequence prev = iter.next();
                boolean prevMerged = false;
                while (iter.hasNext()) {
                    Sequence s = iter.next();
                    
                    Sequence[] merged = prev.merge(s);
                    if (merged == null) {
                        seqs2.add(prev);
                        seqs2.add(s);
                        prevMerged = false;
                        prev = s;
                    } else {
                        didMerge = true;
                        seqs2.remove(prev);
                        seqs2.remove(s);
                        for (Sequence st : merged) {
                            seqs2.add(st);
                            prev = st;
                        }
                    }
                }

                seqs.clear();
                seqs.addAll(seqs2);
                seqs2.clear();
                
                nIter++;
            
            } while (didMerge);
            
            System.out.println("final nMerged=" + seqs.size());
            
            for (Sequence s : seqs) {
                System.out.println("SEQ " + s);
            }
              
            /*
            assertTrue(seqs.size() <= 4);
            boolean found0 = false;
            boolean found1 = false;
            for (Sequence seq : seqs) {
                if (seq.isStartSentinel()) {
                    found0 = true;
                }
                if (seq.isStopSentinel()) {
                    found1 = true;
                }
            }
            assertTrue(found0);
            assertTrue(found1); 
            */
        }
       
    }

    private Sequence createSequence(int n1, int n2, 
        SecureRandom sr, int offset) {
        
        /*
        n2 >= 0
        offset>=0
        offset=stopIdx2-stopIdx1=startIdx2-startIdx1
        len=stopIdx2-startIdx2=stopIdx1-startIdx1
           stopIdx2  range: len-1+offset  : n1-1 + len + offset
           startIdx2 range: offset        : n1-1 + offset
        
           stopIdx1  range: len-1         : n1-1 + len
           startIdx1 range:   0           : n1-1
       
        len = randomInt(n1-1)
        
        the above will be missing the last numbers in 
        idx1
        */
        
        //len = stopIdx2 - startIdx2 + 1
        int len = sr.nextInt(n1) + 1;
        
        Sequence s = new Sequence(n1, n2, offset);
        s.sumDiffs = 0;
        s.startIdx1 = sr.nextInt(n1);
        int stopIdx1 = s.startIdx1 + len - 1;
        if (stopIdx1 >= n1) {
            // truncate an rewrite length len
            len -= (stopIdx1 - (n1 - 1));
            if (len <= 0) {
                // unlikely to be caught in endless loop
                return createSequence(n1, n2, sr, offset);
            }
            stopIdx1 = s.startIdx1 + (len - 1);
        }
        s.startIdx2 = s.startIdx1 + offset;
        if (s.startIdx2 >= n2) {
            s.startIdx2 -= n2;
        }
        s.stopIdx2 = s.startIdx2 + len - 1;

        if (s.stopIdx2 < n2) {
            assert(s.length() == len);
            return s;
        } else {
            /* needs to return 2 sequences to
               keep idx2 ranges increasing in value
            example:
            39:42 to 2  len=11, n1=n2=50
               should be    
            39:42 to 49  len=8
            47:0  to 2
               since the code is handling the mergine,
               and there are many random sequences,
               including the 2nd portion,
               will just truncate the range here
               instead of returning 2 sequences. 
            */
            len -= (s.stopIdx2 - (n2-1));
            if (len <= 0) {
                // unlikely to be caught in endless loop
                return createSequence(n1, n2, sr, offset);
            }
            assert(s.length() == len);
            return s;
        }     
    }
    
}
