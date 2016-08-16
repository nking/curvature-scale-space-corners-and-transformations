package algorithms.imageProcessing.matching;

import java.lang.reflect.Field;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.logging.Logger;
import sun.misc.Unsafe;

/**
 <pe>
The rules for sequence fields are
presented here.  note that for
best use as an api, they should be
encapsulated by Sequences.java or other
means.
  The idx1 axis range is 0 to n1-1.
  The idx2 axis range is 0 to n2-1.
  The field offset is always g.e. 0.
  The field startIdx1 is in range of n1.
  The fields startIdx2 and stopIdx2 are
  in range of n2. 
  stopIdx2 is g.e. startIdx2.
  The length of the sequence is stopIdx2 - startIdx2 + 1.
  The implied stopIdx1 is g.e. startIdx1.
Note that because the indexes are on
  closed curves, the indexes at n1-1 are
  adjacent to those at index 0.
  The wrap-around is handled by cutting the
  sequence into pieces that obey the above rules
  and by adding a sentinel sequence when the
  wrap around has occurred.
  The sentinel sequence is setting the startIdx2
  and stopIdx2 to value n2 (and the corresponding
  startIdx1 calculated from the offset).

    For example, a set of sequences that have
    offset=1, n1=n2=4:
                      start  stop  start   stop
                       idx1   idx1   idx2   idx2
      idx1: 0 1 2 3      0      2      1      3  --- seq 0
      idx2: 1 2 3 0      3      3      4      4  --- seq 2
                         3      3      0      0  --- seq 1

 Another example, written in format written by toString():
       offset=9, n1=n2=50

       SEQ (0:9 to 49, f=0.8200 d=0.0000  n1=50, n2=50)
       SEQ (41:0 to 8, f=0.0800 d=NaN  n1=50, n2=50)
       SEQ (41:50 to 50, f=0.0000 d=0.0000  n1=50, n2=50)
 </pre>
 * @author nichole
 */
public class Sequence {
    
    private Logger log = Logger.getLogger(
        this.getClass().getName());

    int startIdx1;
    int startIdx2 = -1;
    int stopIdx2 = -1;
    float absAvgSumDiffs;
    float fractionOfWhole;

    private final int n1;
    private final int n2;
    private final int offset;
    
    public Sequence(int nIndexes1, int nIndexes2,
        int offset12) {
        n1 = nIndexes1;
        n2 = nIndexes2;
        offset = offset12;
    }
    
    public int getN1() {
        return n1;
    }
    
    public int getN2() {
        return n2;
    }
    
    public int getStartIdx1() {
        return startIdx1;
    }

    public int getStartIdx2() {
        return startIdx2;
    }

    public int getStopIdx2() {
        return stopIdx2;
    }

    /**
     * TODO: this method needs many tests
     * 
     * @param sTest
     * @return 
     */
    public boolean intersects(Sequence sTest) {
        
        if (sTest.getN1() != n1 || sTest.getN2() != n2) {
            throw new IllegalArgumentException("the n1 or n2 "
            + " of sTest have to be same as this n1 and n2");
        }
        
        if ((sTest.isStartSentinel() || sTest.isStopSentinel())
            && (isStartSentinel() || isStopSentinel())) {
            return true;
        }
           
        int sTestStopIdx1 = sTest.startIdx1 +
            (sTest.stopIdx2 - sTest.startIdx2);
         
        int stopIdx1 = startIdx1
            + (stopIdx2 - startIdx2);
        
        //TODO: might need corrections for 
        // wrap around here
        
        //       s s
        //  t t      t t
        // test that first indexes don't intersect
        if ((sTestStopIdx1 < startIdx1)
            || (sTest.startIdx1 > stopIdx1)) {
            // test that second indexes don't intersect
            if ((sTest.stopIdx2 < startIdx2)
                || (sTest.startIdx2 > stopIdx2)) {
                return false;
            }
        }
       
        return true;
    }
    
    /**
     * this method works best if sequences have the same offset
     * due to internal re-ordering of the sequences.
     * @param sequences 
     */
    public static void mergeSequences(List<Sequence> sequences) {

        LinkedHashSet<Sequence> seqs2 = 
            new LinkedHashSet<Sequence>();
        
        int nIter = 0;
        
        int nIterMax = 2 * sequences.size();
        
        boolean didMerge = false;
        do {

            if (sequences.size() == 1) {
                // break because loop below won't store it
                break;
            }

            // sort sequences
            Collections.sort(sequences, new SequenceComparator4());

            if (sequences.size() <= 3) {
                boolean found0 = false;
                boolean found1 = false;
                for (Sequence seq : sequences) {
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

            System.out.println("seqs.size=" + 
                sequences.size());
            for (Sequence s : sequences) {
                System.out.println(nIter + ": SEQ " + s);
            }

            didMerge = false;
            Iterator<Sequence> iter = sequences.iterator();
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

            sequences.clear();
            sequences.addAll(seqs2);
            seqs2.clear();

            nIter++;

        } while (didMerge);
    }
    
    public int getOffset() {
        return offset;
    }

    /**
     NOTE: not ready for use.  needs alot of tests.... 
     
     merge this sequence with mergeFrom if they have the
     same offsets from idx1 to idx2.
       
       <pre>
        Note, in PartialShapeMatcher, aggregation 
        first proceeds by startIdx1 and only
        up until n1-1.  
        
        For idx2, to handle wrap around, need to
        make sentinel sequences of value n2 and 0.
          
       For example, written in format written by toString():
       offset=9, n1=n2=50

       SEQ (0:9 to 49, f=0.8200 d=0.0000  n1=50, n2=50)
       SEQ (41:0 to 8, f=0.0800 d=NaN  n1=50, n2=50)
       SEQ (41:50 to 50, f=0.0000 d=0.0000  n1=50, n2=50)
       where the later is the stop sentinel sequence.
       </pre>
     * @param mergeFrom
     * @return 
     */
    public Sequence[] merge(Sequence mergeFrom) {
        
        if (mergeFrom.getN1() != n1 || mergeFrom.getN2() != n2) {
            throw new IllegalArgumentException("the n1 or n2 "
            + " of mergeFrom have to be same as this n1 and n2");
        }
        
        if (this.isStopSentinel() && 
            mergeFrom.isStartSentinel()) {
            // return both, they should be kept, but
            // cannot be merged
            log.info("*sentinels:" + this + "\n " + mergeFrom);
            return new Sequence[]{this.copy(), mergeFrom};
        } else if (mergeFrom.isStopSentinel() &&
            isStartSentinel()) {
            log.info("**sentinels:" + mergeFrom + "\n " + this);
            return new Sequence[]{mergeFrom, this.copy()};
        } else if (mergeFrom.stopIdx2 == (n2 - 1) 
            && isStartSentinel()) {
           
            // this sequence should be followed by a 
            // sequence indicating it's the end sentinel, then
            // can return all 3 sequences
            /*
            example w/ n1=n2=50 and offset=3
            48 :  1   49
            47 : 50   50
            47 :  0    0
            */
            
            // merge this and mergeFrom, then add a stop sentinel
            Sequence mergeInto = copy();
            
            int len0 = mergeInto.length();
            float f0 = mergeInto.fractionOfWhole;
            float d0 = mergeInto.absAvgSumDiffs;
            float d0Tot = d0 * len0;
            int nAdded = mergeFrom.stopIdx2 - 
                mergeInto.stopIdx2;
            len0 += nAdded;
            float d1Tot = mergeFrom.absAvgSumDiffs * mergeFrom.length();
            d0Tot += (d1Tot/(float)nAdded);
            mergeInto.stopIdx2 = mergeFrom.stopIdx2;
            mergeInto.fractionOfWhole = (float)len0/(float)n1;
            mergeInto.absAvgSumDiffs = d0Tot/(float)len0;
            assert(mergeInto.length() == len0);     
            
            Sequence endS = new Sequence(n1, n2, offset);
            endS.startIdx1 = n1 - offset;
            endS.startIdx2 = n2;
            endS.stopIdx2 = n2;
            // in order to not add to the class Sequences stats, 
            // the end sentinel sequence will have values
            // that do not affect the result, but this
            // needs to be enforced elsewhere too.
            endS.fractionOfWhole = 0;
            endS.absAvgSumDiffs = 0;
            log.info("***sentinels: \n" + 
                mergeFrom + "\n " + this
                + "\n  ==> " + mergeInto 
                + "\n  ==> " + endS);
            return new Sequence[]{mergeInto, endS};
        } else if (stopIdx2 == (n2 - 1) && 
            mergeFrom.isStartSentinel()) {
            // this sequence should be followed by a 
            // sequence indicating it's the end sentinel, then
            // can return all 3 sequences
            /*
            example w/ n1=n2=50 and offset=3
            48 :  1   49
            47 : 50   50 
            47 :  0    0 
            */
            
            // merge this and mergeFrom, then add a stop sentinel
            Sequence mergeInto = mergeFrom.copy();
            
            int len0 = mergeInto.length();
            float f0 = mergeInto.fractionOfWhole;
            float d0 = mergeInto.absAvgSumDiffs;
            float d0Tot = d0 * len0;
            int nAdded = this.stopIdx2 - mergeInto.stopIdx2;
   
            len0 += nAdded;
            float d1Tot = this.absAvgSumDiffs * this.length();
            d0Tot += (d1Tot/(float)nAdded);
            mergeInto.stopIdx2 = this.stopIdx2;
            mergeInto.fractionOfWhole = (float)len0/(float)n1;
            mergeInto.absAvgSumDiffs = d0Tot/(float)len0;
            assert(mergeInto.length() == len0);
            
            Sequence endS = new Sequence(n1, n2, offset);
            endS.startIdx1 = n1 - offset;
            endS.startIdx2 = n2;
            endS.stopIdx2 = n2;
            // in order to not add to the class Sequences stats, 
            // the end sentinel sequence will have values
            // that do not affect the result, but this
            // needs to be enforced elsewhere too.
            endS.fractionOfWhole = 0;
            endS.absAvgSumDiffs = 0;
            log.info("****sentinels:" + 
                this + " \n " + endS 
                + "\n " + mergeFrom);
            return new Sequence[]{this.copy(), endS, mergeFrom};
        } else if (mergeFrom.stopIdx2 == (n2 - 1) &&
            isStopSentinel()) {
            // this is a proper sentinel boundary, so return
            // both in order
            log.info("existing sentinels:" + 
                mergeFrom + " \n " + "\n " + this);
            return new Sequence[]{mergeFrom, this.copy()};
        } else if ((stopIdx2 == (n2 - 1)) && mergeFrom.isStopSentinel()) {
            // this is a proper sentinel boundary, so return 
            // both in order
            log.info("*existing sentinels:" + 
                this + " \n " +  "\n " + mergeFrom);
            return new Sequence[]{this.copy(), mergeFrom};
        }
       
        Sequence mergeInto;
        
        // make mergeInto.startIdx1 < mergeFrom.startIdx1
        if (mergeFrom.startIdx1 < startIdx1) {
            mergeInto = mergeFrom;
            mergeFrom = this.copy();
        } else {
            mergeInto = this.copy();
        }
                
        if (mergeInto.getOffset() != 
            mergeFrom.getOffset()) {
            log.info("NO MERGE. offsets=" + 
                mergeInto.getOffset() + " " +
                mergeFrom.getOffset());
            return null;
        }
        
        if (mergeInto.isStopSentinel() || mergeFrom.isStopSentinel()) {
            return new Sequence[]{mergeInto, mergeFrom};
        }
        
        if (mergeInto.equals(mergeFrom)) {
            // these are same ranges, so let invoker remove mergeFrom
            log.fine("same sequence: " + " \n " + 
                mergeInto + " \n " + mergeFrom);
            return new Sequence[]{mergeInto};
        }

        // a rule is that for a single sequence,
        // the implied stops should be in range n1 too        
        int mIStopIdx1 = mergeInto.startIdx1 +
            (mergeInto.stopIdx2 - mergeInto.startIdx2);
        
        int mFStopIdx1 = mergeFrom.startIdx1 +
            (mergeFrom.stopIdx2 - mergeFrom.startIdx2);

        assert(mIStopIdx1 <= n1);
        assert(mFStopIdx1 <= n1);
        
        // can merge if they are adjacent or 
        // intersecting.

//if (!(mergeFrom.stopIdx2 >= mergeFrom.startIdx2)) {
    log.info("****CHECK: " + 
    "\n  mergeFrom=" + mergeFrom +
    "\n  mergeIntp=" + mergeInto
    + "\n   mIStopIdx1=" + mIStopIdx1);
//}

        assert(mergeInto.stopIdx2 >= mergeInto.startIdx2);
        assert(mergeFrom.stopIdx2 >= mergeFrom.startIdx2);
        
        // -- check for adjacent (before intersects filter) --
        if ((mIStopIdx1 + 1) == mergeFrom.startIdx1) {

            StringBuilder sb = new StringBuilder("*merge ")
                .append(mergeInto).append("\n into ")
                .append(mergeFrom);
            
            int len0 = mergeInto.length();
            float f0 = mergeInto.fractionOfWhole;
            float d0 = mergeInto.absAvgSumDiffs;
            float d0Tot = d0 * len0;
            int nAdded = mergeFrom.stopIdx2 - 
                mergeInto.stopIdx2;
            len0 += nAdded;
            
            float d1Tot = mergeFrom.absAvgSumDiffs * mergeFrom.length();
            d0Tot += (d1Tot/(float)nAdded);
            
            mergeInto.stopIdx2 = mergeFrom.stopIdx2;
            mergeInto.fractionOfWhole = (float)len0/(float)n1;
            mergeInto.absAvgSumDiffs = d0Tot/(float)len0;
            assert(mergeInto.length() == len0);
            
            sb.append("\n => ").append(mergeInto.toString());
            
            // if stopIdx2 == n2-1, create a sentinel too
            if (mergeInto.stopIdx2 == (n2 - 1)) {
                /*
                example w/ n1=n2=50 and offset=3
                48 :  1   49
                47 : 50   50
                47 :  0    0
                */
                Sequence endS = new Sequence(n1, n2, offset);
                endS.startIdx1 = n1 - offset;
                endS.startIdx2 = n2;
                endS.stopIdx2 = n2;
                // in order to not add to the class Sequences stats, 
                // the end sentinel sequence will have values
                // that do not affect the result, but this
                // needs to be enforced elsewhere too.
                endS.fractionOfWhole = 0;
                endS.absAvgSumDiffs = 0;
                
                sb.append("\n  + sentinel=" + endS);
                
                return new Sequence[]{mergeInto, endS};
            }
            
            log.info(sb.toString());
            
            return new Sequence[]{mergeInto};
        }  
        if (!mergeInto.intersects(mergeFrom)) {
            log.info("NO MERGE: not intersecting: " + 
                " \n" + mergeInto + "\n " + mergeFrom);
            return null;
        }
        
        /*
        Note, in PartialShapeMatcher, aggregation 
        first proceeds by startIdx1 and only
        up until n1-1.  
        
        For idx2, to handle wrap around, need to
        make sentinel sequences of value n2 and 0.
          
        examples w/ n1=4; n2=4;
        
                          start  stop  start   stop
                          idx1   idx1   idx2   idx2
        idx1: 0 1 2 3      0      3      0      3
        idx2: 0 1 2 3         
        
        a wrap around 
        idx2 offset=1:     
                          idx1   idx1   idx2   idx2
        idx1: 0 1 2 3      0      2      1      3
        idx2: 1 2 3 0      3      3      4      4
                           3      3      0      0
        
        */

        // --handle the remaining cases near the sentinels--
        
        // mergeInto.startIdx1 < mergeFrom.startIdx1
        // and the case where both are sentinels already handled
        if (mergeInto.isStopSentinel()) {
            // the next would need to be start sentinel
            // and that case has already been handled
            log.info("NO MERGE: has a stop sentinel: " + 
                " \n" + mergeInto + "\n " + mergeFrom);
            return null;
        }
        
        // all cases involving sentinel should be handled
        // by now and so is the adjacent rnage case,
        // so what remains in merging overlapping regions
        // that do not cross the sentinels
        
        // they're ordered by startIdx1.  startIdx2 should be
        // increasing also
        StringBuilder sb = new StringBuilder("**merge ")
            .append(mergeInto).append("\n into ")
            .append(mergeFrom);

        if (mergeInto.stopIdx2 >= mergeFrom.stopIdx2) {
            // mergeInto includes all of mergeFrom already
            log.fine("NO MERGE: already includes other: " + 
                " \n" + mergeInto + "\n " + mergeFrom);
            return new Sequence[]{mergeInto};
        }
        int len0 = mergeInto.length();
        float f0 = mergeInto.fractionOfWhole;
        float d0 = mergeInto.absAvgSumDiffs;
        float d0Tot = d0 * len0;
        int nAdded = mergeFrom.stopIdx2 - mergeInto.stopIdx2;
        len0 += nAdded;
        
        float d1Tot = mergeFrom.absAvgSumDiffs * mergeFrom.length();
        d0Tot += (d1Tot/(float)nAdded);
            
        mergeInto.stopIdx2 = mergeFrom.stopIdx2;
        mergeInto.fractionOfWhole = (float)len0/(float)n1;
        mergeInto.absAvgSumDiffs = d0Tot/(float)len0;
        assert(mergeInto.length() == len0);
        
        sb.append("\n => ").append(mergeInto.toString());
        log.info(sb.toString());

        return new Sequence[]{mergeInto};
    }

    public boolean isStartSentinel() {
        return (startIdx2 == 0);
    }
    public boolean isStopSentinel() {
        return (stopIdx2 == n2);
    }

    public int length() {
        return stopIdx2 - startIdx2 + 1;
    }
    
    @Override
    public boolean equals(Object obj) {

        if (obj == null || !(obj instanceof Sequence)) {
            return false;
        }

        Sequence other = (Sequence) obj;

        if (startIdx1 == other.startIdx1
            && startIdx2 == other.startIdx2
            && stopIdx2 == other.stopIdx2 
            && n1 == other.n1 && n2 == other.n2) {
            return true;
        }

        return false;
    }

    protected static int fnv321aInit = 0x811c9dc5;
    protected static int fnv32Prime = 0x01000193;

    protected int fnvHashCode(int i0, int i1,
        int i2, int i3, int i4) {
        
        /*
         * hash = offset_basis
         * for each octet_of_data to be hashed
         *     hash = hash xor octet_of_data
         *     hash = hash * FNV_prime
         * return hash
         *
         * Public domain:  http://www.isthe.com/chongo/src/fnv/hash_32a.c
         */
        int hash = 0;
        int sum = fnv321aInit;
        sum ^= i0;
        sum *= fnv32Prime;        
        sum ^= i1;
        sum *= fnv32Prime;
        sum ^= i2;
        sum *= fnv32Prime;
        sum ^= i3;
        sum *= fnv32Prime;
        sum ^= i4;
        sum *= fnv32Prime;
        
        hash = sum;
        
        return hash;
    }
    
    @Override
    public int hashCode() {
        
        int hash = fnvHashCode(this.startIdx1,
            this.startIdx2, this.stopIdx2, this.n1,
            this.n2);

        return hash;
    }
    
    public void transpose() {
        
        //TODO: this may need revision.
        // the idx1 axis can wrap around and the
        // idx2 axis does not.
        // for consistency, may need to break
        // a sequence into 2 to keep that pattern
        // and may need to merge sequences after
        // these changes.
                
        int n = stopIdx2 - startIdx2;
        int tmpStartIdx2 = startIdx1;
        startIdx1 = startIdx2;
        startIdx2 = tmpStartIdx2;
        stopIdx2 = startIdx2 + n;
            
        String errMsg = "jvm not allowing this method to "
            + " change 2 variables.  The code needs to be edited "
            + " and recompiled if this is the case.";
        
        try {
            // reset n1 and n2 which are final fields
            
            Field f = Unsafe.class.getDeclaredField("theUnsafe");
            f.setAccessible(true);
            Unsafe unsafe = (Unsafe) f.get(null);

            int tmpN1 = n1;
            int tmpN2 = n2;

            Field f2 = this.getClass().getDeclaredField("n1");
            long offset = unsafe.objectFieldOffset(f2);
            unsafe.putObject(this, offset, tmpN2);

            f2 = this.getClass().getDeclaredField("n2");
            offset = unsafe.objectFieldOffset(f2);
            unsafe.putObject(this, offset, tmpN1);

            assert (n1 == tmpN2);
            assert (n2 == tmpN1);

        } catch (NoSuchFieldException ex) {
            log.severe(errMsg + " : " + ex.getMessage());
        } catch (SecurityException ex) {
            log.severe(errMsg + " : " + ex.getMessage());
        } catch (IllegalArgumentException ex) {
            log.severe(errMsg + " : " + ex.getMessage());
        } catch (IllegalAccessException ex) {
            log.severe(errMsg + " : " + ex.getMessage());
        }
        
    }
    
    boolean assertNoWrapAround() {
        
        int stopIdx1 = startIdx1 + (stopIdx2 - startIdx2);
    
        return (stopIdx1 < n1);
    }
    
    public Sequence copy() {
        Sequence cp = new Sequence(n1, n2, offset);
        cp.startIdx1 = startIdx1;
        cp.startIdx2 = startIdx2;
        cp.stopIdx2 = stopIdx2;
        cp.absAvgSumDiffs = absAvgSumDiffs;
        cp.fractionOfWhole = fractionOfWhole;
        return cp;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format(
            "(%d:%d to %d, f=%.4f d=%.4f  n1=%d, n2=%d offset=%d)",
            startIdx1, startIdx2, stopIdx2,
            fractionOfWhole, absAvgSumDiffs, n1, n2,
            offset));
        return sb.toString();
    }
}
