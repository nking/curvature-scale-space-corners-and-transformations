package algorithms.imageProcessing.matching;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;

import java.util.Arrays;

public class Match implements Comparable<Match> {

    /**
     * start index of diagonal block in A1 (i.e. p).  the matched range is starts1[i] through stops1[],
     * inclusive.
     * the indices are encoded into longs with assumption that the number of points in a curve
     * fits within 8 bits.
     */
    protected final long[] starts1;
    /**
     * the sizes of the matched blocks.
     * stops1 is calculated as starts1[i] + blocks[i] - 1.
     * the matched range is starts1[i] through stops1[], inclusive.
     */
    protected final long[] blocks;

    /**
     * the offset from starts1 index.  the A2 index (i.e. q) is calculated as starts1[i] + offset[i]
     * and individual points may need mod n2 for wrap around.
     */
    protected final long[] offsets2;
    /**
     * sum of the absolute value of chord differences between the range (stops1[i] through starts1[i]) - range (stops2[i] through starts2[i])
     */
    protected final double[] diffChordSums;

    /**
     * the index of the last item added to the arrays.
     */
    protected int lastIdx = -1;

    protected int BITS_PER_NUM = 8;

    // machine precision to use for equals of maxChordSum.  note that float and double representations are used throughout code.
    final static double tol = 1E-10;

    protected double diffChordSum = 0.;
    protected double maxChordSum = Double.NEGATIVE_INFINITY;
    protected int mLen = 0;
    protected final int n1;
    protected final int n2;
    protected final int nMaxMatchable;

    /**
     *
     * @param n1
     * @param n2
     */
    public Match(int n1, int n2) {
        if (n1 < 1 || n2 < 1) {
            throw new IllegalArgumentException(
                    String.format("n1 and n2 must be greater than 0. received n1=%, n2=%d", n1, n2));
        }
        if (n1 > (1 << 16) || n2 > (1 << 16)) {
            throw new IllegalArgumentException("n1 and n2 must be less than " + (1<<16));
        }
        this.n1 = n1;
        this.n2 = n2;
        this.nMaxMatchable = Math.min(n1, n2);
        // num bits need to represent a 0-based index
        if (nMaxMatchable > 256) {
            BITS_PER_NUM = 16;
        } else {
            BITS_PER_NUM = 8;
        }
        int len = (int)Math.ceil((float)(Math.max(n1, n2))/BITS_PER_NUM);
        this.starts1 = new long[len];
        this.blocks = new long[len];
        this.offsets2 = new long[len];
        this.diffChordSums = new double[nMaxMatchable];
    }

    /**
     * class to hold the result matches in a format easy to use and reverse p,with q, and scal if needed.
     */
    public static class Points {
        int[] pIdxs;
        int[] qIdxs;
        /**
         * sum of the chord differences w.r.t. the compared reference frames in A1, and A2
         */
        double chordDiffSum;
        /**
         * number of matching points in the compared reference frames of A1 and A2.
         */
        int mLen;

        public Points(Match m) {
            TIntList rStarts1 = new TIntArrayList();
            TIntList rStops1 = new TIntArrayList();
            TIntList rStarts2 = new TIntArrayList();
            TIntList rStops2 = new TIntArrayList();
            for (int i = 0; i < m.starts1.length; ++i) {
                long start1Row = m.starts1[i];
                long offsetRow = m.offsets2[i];
                long blockRow = m.blocks[i];
                long tmps, tmpo, tmpb;
                for (int colNum = 0; colNum < m.BITS_PER_NUM; ++colNum) {
                    // mask from MSB to bit colNum * BITS_PER_NUM
                    long x = (colNum + 1) * m.BITS_PER_NUM;
                    tmps = start1Row & ((1L << x) - 1);
                    tmpo = offsetRow & ((1L << x) - 1);
                    tmpb = blockRow & ((1L << x) - 1);
                    long y = colNum  * m.BITS_PER_NUM;
                    tmps >>= y;
                    tmpo >>= y;
                    tmpb >>= y;
                    rStarts1.add((int)tmps);
                    rStops1.add((int)(tmps + tmpb - 1));
                    rStarts2.add((int)(tmps + tmpo));
                    rStops2.add((int)(tmps + tmpb + tmpo - 1));
                    if ((i * m.BITS_PER_NUM + colNum) == m.lastIdx) {
                        break;
                    }
                }
            }
            TIntList t1 = new TIntArrayList();
            TIntList t2 = new TIntArrayList();
            for (int i = 0; i < rStarts1.size(); ++i) {
                for (int idx1 = rStarts1.get(i), idx2 = rStarts2.get(i); idx1 <= rStops1.get(i); ++idx1, ++idx2) {
                    t1.add(idx1);
                    t2.add(idx2 % m.n2);
                }
            }
            this.pIdxs = t1.toArray();
            this.qIdxs = t2.toArray();
            this.chordDiffSum = m.diffChordSum;
            this.mLen = m.mLen;
            assert(this.pIdxs.length == this.mLen);
        }

        public void scale(int factor1, int factor2) {
            for (int i = 0; i < pIdxs.length; ++i) {
                pIdxs[i] *= factor1;
                qIdxs[i] *= factor2;
            }
        }
        public void interchange() {
            for (int i = 0; i < pIdxs.length; ++i) {
                if (pIdxs[i] != qIdxs[i]) {
                    pIdxs[i] ^= qIdxs[i];
                    qIdxs[i] ^= pIdxs[i];
                    pIdxs[i] ^= qIdxs[i];
                }
            }
        }
    }

    public Match copy() {
        Match t = new Match(this.n1, this.n2);
        System.arraycopy(starts1, 0, t.starts1, 0, starts1.length);
        System.arraycopy(blocks, 0, t.blocks, 0, blocks.length);
        System.arraycopy(offsets2, 0, t.offsets2, 0, offsets2.length);
        System.arraycopy(diffChordSums, 0, t.diffChordSums, 0, diffChordSums.length);
        t.mLen = this.mLen;
        t.maxChordSum = this.maxChordSum;
        t.diffChordSum = this.diffChordSum;
        t.lastIdx = this.lastIdx;
        return t;
    }

    /**
     * add a match range to the interval
     * @param idx the index of diagonal position with respect to A1.
     * @param offset2
     * @param blockSize
     * @param diffChordSum
     * @param matchNumber the number of this match with respect to the matched sequence it is in.  This
     *                    parameter is needed for use with backtracking.
     */
    public void add(final int idx, final int offset2, final int blockSize, final double diffChordSum,
                    final int matchNumber) {

        if (matchNumber == nMaxMatchable) {
            throw new IllegalStateException("nMaxMatchable is " + nMaxMatchable
                    + " but you are attempting to insert more than that.  ");
        }
        if (matchNumber > lastIdx + 1) {
            throw new IllegalStateException(String.format("Error in algorithm: lastIdx=%d and match_number to insert is %d",
                    lastIdx, matchNumber));
        }
        if (idx > (1<<BITS_PER_NUM) || idx + offset2 > (1<<BITS_PER_NUM)) {
            throw new IllegalArgumentException("expecting curve indices to be < " + (1 << BITS_PER_NUM));
        }
        if (blockSize > (1<<BITS_PER_NUM)) {
            throw new IllegalArgumentException("expecting blockSize to be < " + (1 << BITS_PER_NUM));
        }
        lastIdx = matchNumber;
        int r = lastIdx / BITS_PER_NUM;
        int c = lastIdx % BITS_PER_NUM;
        this.starts1[r] |= ((long)idx << c * BITS_PER_NUM);
        this.blocks[r] |= ((long)blockSize << c * BITS_PER_NUM);
        this.offsets2[r] |= ((long)offset2 << c * BITS_PER_NUM);
        this.diffChordSums[lastIdx] = diffChordSum;
        this.mLen += blockSize;
        this.diffChordSum += diffChordSum;
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof Match)) {
            return false;
        }
        Match other = (Match)obj;
        if (n1 != other.n1) return false;
        if (n2 != other.n2) return false;
        if (BITS_PER_NUM != other.BITS_PER_NUM) return false;
        if (mLen != other.mLen) return false;
        if (lastIdx != other.lastIdx) return false;
        // no need to compare maxChordSum as that will be reset upong use, though could reconsider this.
        if (nMaxMatchable != other.nMaxMatchable) return false;
        if (Math.abs(diffChordSum - other.diffChordSum) > tol) return false;
        if (!Arrays.equals(starts1, other.starts1)) return false;
        if (!Arrays.equals(blocks, other.blocks)) return false;
        if (!Arrays.equals(offsets2, other.offsets2)) return false;
        if (!listEquals(diffChordSums,other.diffChordSums)) return false;
        return true;
    }

    private boolean listEquals(double[] a, double[] b) {
        if (a.length != b.length) return false;
        for (int i = 0; i < a.length; ++i) {
            if (Math.abs(a[i] - b[i]) > tol) return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = fnvHashCode();
        return hash;
    }
    //Public domain:  http://www.isthe.com/chongo/src/fnv/hash_32a.c
    protected static int fnv321aInit = 0x811c9dc5;
    protected static int fnv32Prime = 0x01000193;
    protected int fnvHashCode() {
        int sum = fnv321aInit;

        for (int i = 0; i < starts1.length; ++i) {
            sum ^= starts1[i];
            sum *= fnv32Prime;
            sum ^= blocks[i];
            sum *= fnv32Prime;
            sum ^= offsets2[i];
            sum *= fnv32Prime;
        }
        sum ^= diffChordSums.hashCode();
        sum *= fnv32Prime;

        sum ^= lastIdx;
        sum *= fnv32Prime;

        sum ^= mLen;
        sum *= fnv32Prime;

        sum ^= nMaxMatchable;
        sum *= fnv32Prime;

        sum ^= n1;
        sum *= fnv32Prime;

        sum ^= n2;
        sum *= fnv32Prime;

        sum ^= BITS_PER_NUM;
        sum *= fnv32Prime;

        sum ^= Double.hashCode(diffChordSum);
        sum *= fnv32Prime;
        sum ^= Double.hashCode(maxChordSum);
        sum *= fnv32Prime;

        return sum;
    }

    /**
     * compare to objects using thire Salukwdze distances.
     * for best results, make sure that all objects being compared have the same maxChordSum and
     * nMaxMatchable.
     * @param other the object to be compared.
     * @return
     */
    @Override
    public int compareTo(Match other) {
        double maxDiffChord = Math.max(maxChordSum, other.maxChordSum);
        int maxLen = Math.max(nMaxMatchable, other.nMaxMatchable);
        double d1 = calcSalukDist(diffChordSum, maxDiffChord, mLen, maxLen);
        double d2 = calcSalukDist(other.diffChordSum, maxDiffChord, other.mLen, maxLen);
        return Double.compare(d1, d2);
    }

    double calcSalukDist(double compChord, double maxChord, int length, int maxMatchable) {
        double d;
        if (maxChord == 0 || Double.isInfinite(maxChord)) {
            d = 0;
        } else {
            d = compChord / maxChord;
        }
        double f = 1. - ((double)length/(double)maxMatchable);
        return f*f + d*d;
    }
}
