package algorithms.imageProcessing.matching;

import java.io.Console;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
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
        
        // 7:19 to 21
        // 11:23 to 29
        s1 = new Sequence(n1, n2);
        s1.startIdx1 = 7;
        s1.startIdx2 = 19;
        s1.stopIdx2 = 21;
        assertEquals(12, s1.calcOffset12());
        
        s2 = new Sequence(n1, n2);
        s2.startIdx1 = 11;
        s2.startIdx2 = 23;
        s2.stopIdx2 = 29;
        assertEquals(12, s2.calcOffset12());
        
        assertFalse(s1.intersects(s2));
        
        
        // 17:54 to 62  n1=187  n2=231
        // 15:54 to 64
        s1 = new Sequence(n1, n2);
        s1.startIdx1 = 17;
        s1.startIdx2 = 54;
        s1.stopIdx2 = 62;
        assertEquals(54 - 17, s1.calcOffset12());
        
        s2 = new Sequence(n1, n2);
        s2.startIdx1 = 15;
        s2.startIdx2 = 54;
        s2.stopIdx2 = 64;
        assertEquals(54 - 15, s2.calcOffset12());
        
        assertTrue(s1.intersects(s2));
        
        //43:0 to 8
        //44:0 to 7
        s1 = new Sequence(n1, n2);
        s1.startIdx1 = 43;
        s1.startIdx2 = 0;
        s1.stopIdx2 = 8;
        assertEquals(0 - 43, s1.calcOffset12());
        
        s2 = new Sequence(n1, n2);
        s2.startIdx1 = 44;
        s2.startIdx2 = 0;
        s2.stopIdx2 = 9;
        assertEquals(0 - 44, s2.calcOffset12());
        
        assertTrue(s1.intersects(s2));
        
        n1 = 50;
        n2 = 60;
        s2 = new Sequence(n1, n2);
        s2.startIdx1 = 44;
        s2.startIdx2 = 0;
        s2.stopIdx2 = 9;
        assertEquals(n1 - 44, s2.calcOffset12());
    }
    
    public void estMerge() {
        
        Sequence s1, s2;
        Sequence[] merged;
        int n1 = 187;
        int n2 = 231;
        
        // 2:11  34
        // 13:22 67
        s1 = new Sequence(n1, n2);
        s1.startIdx1 = 2;
        s1.startIdx2 = 11;
        s1.stopIdx2 = 34;
        
        s2 = new Sequence(n1, n2);
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
        s1 = new Sequence(n1, n2);
        s1.startIdx1 = 43;
        s1.startIdx2 = 0;
        s1.stopIdx2 = 8;
        
        s2 = new Sequence(n1, n2);
        s2.startIdx1 = 44;
        s2.startIdx2 = 0;
        s2.stopIdx2 = 9;

        merged = s1.merge(s2);
        assertNull(merged);
        merged = s2.merge(s1);
        assertNull(merged);
    }
    
    public void testMerge2() throws Exception {
        
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
        //seed = 1471239341724L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        System.out.flush();
        
        int n1 = 50;
        int n2 = n1;
        int nSeq = 10 * n2;
        
        int nTests = 1;
        for (int nTest = 0; nTest < nTests; ++nTest) {
            
            int offset12 = sr.nextInt((n1 - 2)/2);
            //if (sr.nextBoolean()) {
            //    offset12 *= -1;
            //}
            
            // hard wire while debugging
            offset12 = 0;
            
            List<Sequence> seqs = 
                new ArrayList<Sequence>();
            
            for (int i = 0; i < nSeq; ++i) {
                Sequence s0 = createSequence(n1, n2, sr, offset12);
                seqs.add(s0);
            }
            
            LinkedHashSet<Sequence> seqs2 = 
                new LinkedHashSet<Sequence>();
            
            int nIter = 0;
            
            // TODO: the sequences need to be sorted
            // by startIdx1, then startIdx2
            
            boolean didMerge = false;
            do {
                
                if (seqs.size() == 1) {
                    // break because loop below won't store it
                    break;
                }
                
                // sort seqs
                Collections.sort(seqs, new SequenceComparator4());
                
                if (seqs.size() == 2) {
                    if (seqs.get(0).isStartSentinel() &&
                        seqs.get(1).isStopSentinel()) {
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
            
            if (seqs.size() == 24) {
               int z = 1;
            }
                    Sequence[] merged = prev.merge(s);
                    if (merged == null) {
                        if (!prevMerged) {
                            seqs2.add(prev);
                        }
                        seqs2.add(s);
                        prevMerged = false;
                        prev = s;
                    } else {
                        didMerge = true;
                        prevMerged = true;
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
            
            //TODO: assert seqs2 is 3 or 4 in size
            // and has the start sentinel and a 
            // stop sentinel
        }
        
        // TODO: test for larger n where the calc offset
        // may need revision
        
    }

    private Sequence createSequence(int n1, int n2, 
        SecureRandom sr, int offset) {
        
        if (offset < 0) {
            return createSequence1(n1, n2, sr, offset);
        }
        
        return createSequence0(n1, n2, sr, offset);
    }
    
    // tailored for positive offsets
    private Sequence createSequence0(int n1, int n2, 
        SecureRandom sr, int offset) {
       
        /*
        n2 >= 0
        offset>=0
        offset=stopIdx2-stopIdx1=startIdx2-startIdx1
        len=stopIdx2-startIdx2=stopIdx1-startIdx1
           stopIdx2  range: len+offset  : n1-1
           startIdx2 range: offset      : n1-1-len
        
           stopIdx1  range: len         : n1-1-offset
           startIdx1 range:   0         : n1-1-len-offset 
       
        len = randomInt(n1-1-offset)
        */
        int len = sr.nextInt(n1 - 2 - offset) + 1;
        
        Sequence s = new Sequence(n1, n2);
        s.startIdx1 = sr.nextInt(n1 - len - offset);
        int stopIdx1 = s.startIdx1 + len;
        s.startIdx2 = s.startIdx1 + offset;
        s.stopIdx2 = s.startIdx2 + len;
        s.fractionOfWhole = (float)len/(float)n1;
        s.absAvgSumDiffs = 0;
        
        return s;
    }
    
    // tailored for offset < 0
    private Sequence createSequence1(int n1, int n2, 
        SecureRandom sr, int offset) {
        
        // calc ranges for offset<0
        throw new UnsupportedOperationException("not yet impl");
    }
    
}
