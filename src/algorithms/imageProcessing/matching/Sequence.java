package algorithms.imageProcessing.matching;

import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.logging.Logger;

/**
 <pre>
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

        example, n1=n2=50, offset=3
            0  : 3 ...
            46 : 49  49
            47 : 0   1
            49 : 2   2

        example, n1=50;n2=55, offset=3
            0  : 3 ...    note, startIdx2 is never 0, or n2-1
            46 : 49  49
            47 : 50  50
            48 : 51  51
            49 : 52  52
        example, n1=50;n2=55, offset=47
             0 : 47  49
             3 : 50
            49 : 46  46

   NOTE that much of the logic below relies
    on n1 leq n2, so if use transpose() to
    reverse correspondence, should not use
    methods such as merge afterwards...
    TODO: change code to create a corresp
    list at transpose stage...
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
        //assert(n1 <= n2);
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
        
        // order by startIdx1
        Sequence s0, s1;
        // make mergeInto.startIdx1 < mergeFrom.startIdx1
        if (sTest.startIdx1 <= startIdx1) {
            s0 = sTest;
            s1 = this;
        } else {
            s0 = this;
            s1 = sTest;
        }
        
        int s0Sentinel1 = (n2 - 1 - offset);
        if (s0Sentinel1 < s1.startIdx1) {
            // can compare idx2 alone to see does not
            // intersect
            if ((s0.getStopIdx2() < s1.startIdx2)
                || (s0.startIdx2 > s1.getStopIdx2())) {
                return false;
            } else {
                return true;
            }
        }

        //TODO: might need corrections for
        // wrap around here

        //       s s
        //  t t      t t
        // test that first indexes don't intersect
        if ((s0.getStopIdx1() < s1.startIdx1)
            || (s0.startIdx1 > s1.getStopIdx1())) {
            // test that second indexes don't intersect
            if ((s0.getStopIdx2() < s1.startIdx2)
                || (s0.startIdx2 > s1.getStopIdx2())) {
                return false;
            }
        }

        return true;
    }

    public int getStopIdx1() {
        int stopIdx1 = startIdx1 + length() - 1;
        if (stopIdx1 >= n1) {
            stopIdx1 -= n1;
        }
        return stopIdx1;
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
        *
        example, n1=n2=50, offset=3
            0  : 3 ...
            46 : 49  49
            47 : 0   1
            49 : 2   2
        example, n1=50;n2=55, offset=3
            0  : 3 ...    note, startIdx2 is never 0, or n2-1
            46 : 49  49
            47 : 50  50
            48 : 51  51
            49 : 52  52
        example, n1=50;n2=55, offset=47
             0 : 47  49
             3 : 50
            49 : 46  46
       </pre>
     * @param mergeFrom
     * @return
     */
    public Sequence[] merge(Sequence mergeFrom) {

        if (mergeFrom.getN1() != n1 || mergeFrom.getN2() != n2) {
            throw new IllegalArgumentException("the n1 or n2 "
            + " of mergeFrom have to be same as this n1 and n2");
        }

        if (getOffset() != mergeFrom.getOffset()) {
            log.info("NO MERGE. offsets=" +
                getOffset() + " " +
                mergeFrom.getOffset());
            return null;
        }

        if (this.isStopSentinel() || mergeFrom.isStopSentinel()) {
            log.info("*sentinels:" + this + "\n " + mergeFrom);
            return null;
        }

        Sequence mergeInto;

        // make mergeInto.startIdx1 < mergeFrom.startIdx1
        if (mergeFrom.startIdx1 < startIdx1) {
            mergeInto = mergeFrom;
            mergeFrom = this.copy();
        } else {
            mergeInto = this.copy();
        }

        if (mergeInto.startIdx2 > mergeFrom.startIdx2) {
            return null;
        }

        if (mergeInto.equals(mergeFrom)) {
            // these are same ranges, so let invoker remove mergeFrom
            log.fine("same sequence: " + " \n " +
                mergeInto + " \n " + mergeFrom);
            return new Sequence[]{mergeInto};
        }

        // a rule is that for a single sequence,
        // the implied stops should be in range n1 too
        int mIStopIdx1 = mergeInto.getStopIdx1();

        int mFStopIdx1 = mergeFrom.getStopIdx1();

//if (!(mergeFrom.stopIdx2 >= mergeFrom.startIdx2)) {
    log.info("****CHECK: " +
    "\n  mergeFrom=" + mergeFrom +
    "\n  mergeInto=" + mergeInto
    + "\n   mIStopIdx1=" + mIStopIdx1
    + "\n   mFStopIdx1=" + mFStopIdx1);
//}

        assert(mIStopIdx1 <= n1);
        assert(mFStopIdx1 <= n1);

        // can merge if they are adjacent or
        // intersecting.

        assert(mergeInto.stopIdx2 >= mergeInto.startIdx2);
        assert(mergeFrom.stopIdx2 >= mergeFrom.startIdx2);

        // -- check for adjacent (before intersects filter) --
        if (
            // idx2 has to have room before n2 for a merge
            (mergeInto.getStopIdx2() <
            (n2 - 1 + (mergeFrom.stopIdx2 - mergeInto.stopIdx2)))
            &&
            (mIStopIdx1 + 1) == mergeFrom.startIdx1) {

            StringBuilder sb = new StringBuilder("*merge ")
                .append(mergeInto).append("\n into ")
                .append(mergeFrom);

            int len0 = mergeInto.length();

            mergeInto.stopIdx2 = mergeFrom.stopIdx2;

            float f0 = mergeInto.fractionOfWhole;
            float d0 = mergeInto.absAvgSumDiffs;
            float d0Tot = d0 * len0;
            int nAdded = mergeInto.length() - len0;
            len0 = mergeInto.length();

            float d1Tot = mergeFrom.absAvgSumDiffs * mergeFrom.length();
            d0Tot += (d1Tot/(float)nAdded);

            mergeInto.fractionOfWhole = (float)len0/(float)n1;
            mergeInto.absAvgSumDiffs = d0Tot/(float)len0;
            assert(mergeInto.length() == len0);

            sb.append("\n => ").append(mergeInto.toString());

            log.info(sb.toString());

            assert(mergeInto.length() <= n1);
            
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
        */

        // they're ordered by startIdx1.  startIdx2 should be
        // increasing also
        StringBuilder sb = new StringBuilder("**merge ")
            .append(mergeFrom).append("\n into ")
            .append(mergeInto);

        int len0 = mergeInto.length();

        mergeInto.stopIdx2 = mergeFrom.stopIdx2;

        float f0 = mergeInto.fractionOfWhole;
        float d0 = mergeInto.absAvgSumDiffs;
        float d0Tot = d0 * len0;
        int nAdded = (mergeInto.length() - len0);
        if (nAdded == 0) {
            sb.append("\n => ").append(mergeInto.toString());
            log.info(sb.toString());

            assert(mergeInto.length() <= n1);
            
            return new Sequence[]{mergeInto};
        }

        len0 = mergeInto.length();

        float d1Tot = mergeFrom.absAvgSumDiffs * mergeFrom.length();
        d0Tot += (d1Tot/(float)nAdded);

        mergeInto.fractionOfWhole = (float)len0/(float)n1;
        mergeInto.absAvgSumDiffs = d0Tot/(float)len0;
        assert(mergeInto.length() == len0);

        sb.append("\n ==> ").append(mergeInto.toString());
        log.info(sb.toString());

        assert(mergeInto.length() <= n1);
        
        return new Sequence[]{mergeInto};
    }

    /**
     used to correct the bounds of a sequence
     which has been aggregated from min diffs,
     for example.
     the idx2 values have not been corrected for
     wrap around so their differences should
     give the correct length and offsets to
     start with.
     note that startIdx1 should always be in n1 range
     at this point already too and n1 is leq n2.
     @param s
     @return
    */
    public static Sequence[] parse(Sequence s) {

        int startIdx1 = s.startIdx1;
        int len = s.length();
        int offset = s.getOffset();
        int stopIdx1 = s.getStopIdx1();

        int startIdx2 = s.startIdx2;
        int stopIdx2 = s.stopIdx2;

        // if location idx1 + offset == n2 is within
        // s range of idx1, then need at least one split

        int t1 = s.n2 - offset;
        if ((t1 >= startIdx1) && (t1 <= stopIdx1)) {
            // s0 startIdx1 to t1-1
            // s1 t1 to stopIdx1

            Sequence[] seqs = new Sequence[2];
            seqs[0] = new Sequence(s.n1, s.n2, offset);
            seqs[0].startIdx1 = startIdx1;
            int len0 = t1 - startIdx1;
            seqs[0].startIdx2 = startIdx1 + offset;
            seqs[0].stopIdx2 = seqs[0].startIdx2 + len0 - 1;

            float f0 = (float)len0/(float)len;
            seqs[0].absAvgSumDiffs = s.absAvgSumDiffs * f0;
            seqs[0].fractionOfWhole = (float)len0/(float)s.n1;

            seqs[1] = new Sequence(s.n1, s.n2, offset);
            seqs[1].startIdx1 = t1;
            int len1 = stopIdx1 - seqs[1].startIdx1 + 1;
            seqs[1].startIdx2 = 0;
            seqs[1].stopIdx2 = seqs[1].startIdx2 + len1 - 1;

            seqs[1].absAvgSumDiffs = s.absAvgSumDiffs
                * (float)len1/(float)len;
            seqs[1].fractionOfWhole = (float)len1/(float)s.n1;

            System.out.println("parsed s=" + s +
                "\n into " + seqs[0] +
                "\n and  " + seqs[1]);

            assert(len0 + len1 == len);

            return seqs;
        }
        
        // look for idx2 being out of range of n2
        if (s.startIdx2 < s.n2 && s.stopIdx2 >= s.n2) {
            
            Sequence[] seqs = new Sequence[2];
            seqs[0] = new Sequence(s.n1, s.n2, offset);
            seqs[0].startIdx1 = startIdx1;
            seqs[0].startIdx2 = startIdx2;
            seqs[0].stopIdx2 = s.n2 - 1;
            int len0 = seqs[0].length();
            
            float f0 = (float)len0/(float)len;
            seqs[0].absAvgSumDiffs = s.absAvgSumDiffs * f0;
            seqs[0].fractionOfWhole = (float)len0/(float)s.n1;

            seqs[1] = new Sequence(s.n1, s.n2, offset);
            seqs[1].startIdx1 = seqs[0].getStopIdx1() + 1;
            seqs[1].startIdx2 = 0;
            int len1 = stopIdx1 - seqs[1].startIdx1 + 1;
            seqs[1].stopIdx2 = seqs[1].startIdx2 + len1 - 1;
        
            seqs[1].absAvgSumDiffs = s.absAvgSumDiffs
                * (float)len1/(float)len;
            seqs[1].fractionOfWhole = (float)len1/(float)s.n1;

            System.out.println("*parsed s=" + s +
                "\n into " + seqs[0] +
                "\n and  " + seqs[1]);

            assert(len0 + len1 == len);

            return seqs;
            
        } else if (s.startIdx2 >= s.n2) {
            
            assert(s.stopIdx2 >= s.n2);
            startIdx2 -= s.n2;
            stopIdx2 -= s.n2;
            assert(startIdx2 < s.n2);
            if (stopIdx2 < s.n2) {
                Sequence[] seqs = new Sequence[1];
                seqs[0] = new Sequence(s.n1, s.n2, offset);
                seqs[0].startIdx1 = startIdx1;
                seqs[0].startIdx2 = startIdx2;
                seqs[0].stopIdx2 = stopIdx2;
                int len0 = seqs[0].length();
                System.out.println("*parsed s=" + s +
                "\n into " + seqs[0]);
                assert(len0 == len);
            
                return seqs;
            }
            
            throw new IllegalStateException(
                "ERROR: Does this case "
               + " happen?  think not because reading the matrix"
               + " dimensions will wrap at most once"
            );            
        }

        return new Sequence[]{s};
    }

    public boolean precedes(Sequence sTest) {
        //TODO: this may need corrections for wrap around
        if (sTest.getStopIdx1() == (startIdx1 - 1)) {
            return true;
        }
        return false;
    }

    public boolean proceeds(Sequence sTest) {
        //TODO: this may need corrections for wrap around
        if (getStopIdx1() == (sTest.getStartIdx1() + 1)) {
            return true;
        }
        return false;
    }

    public boolean isStartSentinel() {
        return (startIdx1 == 0);
    }
    public boolean isStopSentinel() {
        return (getStopIdx1() == (n1 - 1));
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

    /**
     * transpose the sequence to swap reference frames idx1
     * with idx2.  more than one sequence is possibly
     * returned to keep a format following class rules
     * (see class documentation).
     * @return
     */
    public Sequence transpose() {

        Sequence s0 = new Sequence(n2, n1, n1 - offset);

        int len = length();

        s0.startIdx1 = startIdx2;
        s0.startIdx2 = startIdx1;
        s0.stopIdx2 = startIdx1 + len - 1;
        //int stopIdx1 = stopIdx2;

        s0.absAvgSumDiffs = absAvgSumDiffs;
        s0.fractionOfWhole = fractionOfWhole;
        return s0;
    }

    boolean assertNoWrapAround() {

        int stopIdx1 = startIdx1 + length() - 1;

        if (stopIdx1 > n1) {
            log.info("WRAPPED? stopIdx1=" + stopIdx1
                + " : " + toString());
        }

        return (stopIdx1 <= n1);
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
            "(%d:%d to %d, f=%.4f d=%.4f  n1=%d, n2=%d offset=%d len=%d)",
            startIdx1, startIdx2, stopIdx2,
            fractionOfWhole, absAvgSumDiffs, n1, n2,
            offset, length()));
        return sb.toString();
    }

}
