package algorithms.imageProcessing.matching;

import java.lang.reflect.Field;
import java.util.logging.Logger;
import sun.misc.Unsafe;

/**
 *
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
    
    public Sequence(int nIndexes1, int nIndexes2) {
         n1 = nIndexes1;
         n2 = nIndexes2;
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
        
        assert(sTestStopIdx1 < n1);
        assert(stopIdx1 < n1);
                
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
        
    public int calcOffset12() {
        
        /* TODO: this might need to be revised.
           trying to correct the offset for
           correspondence near the sentinels
           idx1: 0 1 2 3 
           idx2: 3 0 1 2
           offset is -1, but would be naively calculated
           as 3. 
           idx1: 1 2 3 4 0  
           idx2: 3 4 0 1 2    5/2ceil=3-1=2 
        */
        
        int offset0 = startIdx2 - startIdx1;
        int offset1 = (n1 - startIdx1) + startIdx2;
        if (Math.abs(offset0) <= Math.abs(offset1)) {
            return offset0;
        }
        
        return offset1;
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
        
        a wrap around 
        idx2 offset=-1
                          idx1   idx1   idx2   idx2
        idx1: 0 1 2 3      0      0      3      3
        idx2: 3 0 1 2      1      1      4      4
                           1      1      0      0
                           2      3      1      2
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
            /*               idx1   idx1   idx2   idx2
            idx1: 0 1 2 3      0      2      1      3
            idx2: 1 2 3 0      3      3      4      4
                               3      3      0      0  
            */
            Sequence endS = new Sequence(n1, n2);
            endS.startIdx1 = 1 + mergeFrom.startIdx1 
                + (mergeFrom.stopIdx2 - mergeFrom.startIdx2);
            endS.startIdx2 = n2;
            endS.stopIdx2 = n2;
            // in order to not add to the class Sequences stats, 
            // the end sentinel sequence will have values
            // that do not affect the result, but this
            // needs to be enforced elsewhere too.
            endS.fractionOfWhole = 0;
            endS.absAvgSumDiffs = 0;
            log.info("***sentinels:" + 
                mergeFrom + " \n " + endS + "\n " + this);
            return new Sequence[]{mergeFrom, endS, this.copy()};
        } else if (stopIdx2 == (n2 - 1) && 
            mergeFrom.isStartSentinel()) {
            // this sequence should be followed by a 
            // sequence indicating it's the end sentinel, then
            // can return all 3 sequences
            /*               idx1   idx1   idx2   idx2
            idx1: 0 1 2 3      0      2      1      3
            idx2: 1 2 3 0      3      3      4      4
                               3      3      0      0  
            */
            Sequence endS = new Sequence(n1, n2);
            endS.startIdx1 = 1 + startIdx1 + (stopIdx2 - startIdx2);
            endS.startIdx2 = n2;
            endS.stopIdx2 = n2;
            // in order to not add to the class Sequences stats, 
            // the end sentinel sequence will have values
            // that do not affect the result, but this
            // needs to be enforced elsewhere too.
            endS.fractionOfWhole = 0;
            endS.absAvgSumDiffs = 0;
            log.info("****sentinels:" + 
                this + " \n " + endS + "\n " + mergeFrom);
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
                
        if (mergeInto.calcOffset12() != 
            mergeFrom.calcOffset12()) {
            log.info("NO MERGE. offsets=" + 
                mergeInto.calcOffset12() + " " +
                mergeFrom.calcOffset12());
            return null;
        }
        
        if (mergeInto.equals(mergeFrom)) {
            // these are same ranges, so let invoker remove mergeFrom
            log.info("same sequence: " + " \n " + 
                mergeInto + " \n " + mergeFrom);
            return new Sequence[]{mergeInto};
        }
        
        int mIStopIdx1 = mergeInto.startIdx1 +
            (mergeInto.stopIdx2 - mergeInto.startIdx2);
        
        int mFStopIdx1 = mergeFrom.startIdx1 +
            (mergeFrom.stopIdx2 - mergeFrom.startIdx2);
        
        assert(mFStopIdx1 < n1 && mIStopIdx1 < n1);

        // can merge if they are adjacent or 
        // intersecting.
        
        /*
                          idx1   idx1   idx2   idx2
        idx1: 0 1 2 3      0      2      1      3
        idx2: 1 2 3 0      3      3      4      4
                           3      3      0      0        
                          idx1   idx1   idx2   idx2
        idx1: 0 1 2 3      0      0      3      3
        idx2: 3 0 1 2      1      1      4      4
                           1      1      0      0
                           2      3      1      2
        */

        assert(mergeInto.stopIdx2 >= mergeInto.startIdx2);
        assert(mergeFrom.stopIdx2 >= mergeFrom.startIdx2);
        
        // -- check for adjacent (before intersects filter) --
        if ((mIStopIdx1 + 1) == mergeFrom.startIdx1) {
           
            StringBuilder sb = new StringBuilder("*merge ")
                .append(mergeInto).append("\n into ")
                .append(mergeFrom);
            
            int len0 = mergeInto.stopIdx2 - mergeInto.startIdx2 + 1;
            float f0 = mergeInto.fractionOfWhole;
            float d0 = mergeInto.absAvgSumDiffs;
            float d0Tot = d0 * len0;
            int nAdded = mergeFrom.stopIdx2 - 
                mergeInto.stopIdx2;
            len0 += nAdded;
            mergeInto.stopIdx2 = mergeFrom.stopIdx2;
            mergeInto.fractionOfWhole = (float)len0/(float)n1;
            mergeInto.absAvgSumDiffs = d0Tot/(float)len0;
            
            sb.append("\n => ").append(mergeInto.toString());
            log.info(sb.toString());
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
        
        a wrap around 
        idx2 offset=-1
                          idx1   idx1   idx2   idx2
        idx1: 0 1 2 3      0      0      3      3
        idx2: 3 0 1 2      1      1      4      4
                           1      1      0      0
                           2      3      1      2
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
        mergeInto.stopIdx2 = mergeFrom.stopIdx2;
        mergeInto.fractionOfWhole = (float)len0/(float)n1;
        mergeInto.absAvgSumDiffs = d0Tot/(float)len0;

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
        Sequence cp = new Sequence(n1, n2);
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
            "(%d:%d to %d, f=%.4f d=%.4f  n1=%d, n2=%d)",
            startIdx1, startIdx2, stopIdx2,
            fractionOfWhole, absAvgSumDiffs, n1, n2));
        return sb.toString();
    }
}
