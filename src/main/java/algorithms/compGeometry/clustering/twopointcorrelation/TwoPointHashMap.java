package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.Arrays;
import java.util.logging.Logger;

/**
 * A holder for two-point identities, where the identities are the indexes
 * of the indexer internal arrays.  N is the size of the dataset, that is
 * the indexer.nXY.
 *
 * TwoPointHashMap can only be used if
 *     (N^2) < Integer.MAX_VALUE
 *     ==> N < 46340
 *
 * For datasets with N larger than 46340, the binary search tree impl should be
 * used instead.
 * TODO:  if memory does not need to be conserved, could consider making many
 *    instances of TwoPointHashMap (Math.ceiling(N/46340.) instances)
 *    to divide N into ranges kept within each TwoPointHashMap instance.
 *
 * This data structure is optimized for inserts and comparisons.
 *
 * Runtime complexity:
 *    inserts are O(1), just the cost of the hash function, so asymptotically constant for large N
 *    comparisons are also O(1)
 *
 * Space complexity:
 *    O(N)
 *
 * @author nichole
 */
class TwoPointHashMap implements ITwoPointIdentity {

    protected int[] a0;
    protected int[] a1;
    protected int n = 0;

    protected final int nDimen;

    public static final int nMax = 46340;

    protected static int emptyVal = -1;

    protected transient Logger log = null;

    /**
     * constructor with data size as argument.  Note that indexerNXY has to be
     * smaller than 46340.
     *
     * @param indexerNXY number of points in the original dataset
     */
    TwoPointHashMap(int indexerNXY) {

        if (indexerNXY > 46340) {
            throw new IllegalArgumentException("please choose another impl of ITwoPointIdentity."
                + " this one can store values if the dataset size is <= 46340");
        }

        log = Logger.getLogger(this.getClass().getName());

        int nt = indexerNXY*indexerNXY;

        log.fine("creating arrays of size " + nt + " for data size " + indexerNXY);

        a0 = new int[nt];
        a1 = new int[nt];

        this.nDimen = indexerNXY;

        Arrays.fill(a0, emptyVal);
        Arrays.fill(a1, emptyVal);
    }

    @Override
    public long approximateMemoryUsed() {

        String arch = System.getProperty("sun.arch.data.model");

        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;

        int nbits = (is32Bit) ? 32 : 64;

        int arrayRefBits = 32;

        int overheadBytes = 16;

        // 4 ints and 2 int arrays
        long sumBits = 2*arrayRefBits + 4*nbits + 2*a0.length*nbits;

        long sumBytes = (sumBits/8) + overheadBytes;

        long padding = (sumBytes % 8);

        sumBytes += padding;

        return sumBytes;
    }

    /**
     * check if combination is already stored, if not add it and return true, else
     * return false
     *
     * @param index0
     * @param index1
     * @return
     */
    @Override
    public boolean storeIfDoesNotContain(int index0, int index1) {

        // order the indexes to avoid double counting.
        int i0, i1;
        if (index0 < index1) {
            i0 = index0;
            i1 = index1;
        } else {
            i0 = index1;
            i1 = index0;
        }

        int hash = hash(i0, i1);

        if (a0[hash] != emptyVal) {
            return false;
        }

        a0[hash] = index0;
        a1[hash] = index1;
        n++;

        return true;
    }

    int hash(int i0, int i1) {

        int i = (nDimen*i0) + i1;

        return i;
    }

}
