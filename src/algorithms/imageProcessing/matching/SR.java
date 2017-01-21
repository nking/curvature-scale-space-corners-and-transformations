package algorithms.imageProcessing.matching;

/**
 *
 * @author nichole
 */
class SR implements Comparable<SR> {
    
    // range of match is idx1 thru idx2, inclusive and if idx2 < idx1,
    //   the interval has wrapped around the closed curve
    int startIdx1;
    int stopIdx1;
    int offsetIdx2;
    int row;
    int mLen;
    int nMax;
    double diffChordSum;
    double maxChordSum;
    boolean chordSumNeedsUpdate = true;

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof SR)) {
            return false;
        }
        SR other = (SR)obj;
        return (other.startIdx1 == startIdx1 && other.stopIdx1 == stopIdx1 &&
            other.mLen == mLen && other.diffChordSum == diffChordSum &&
            other.row == row && other.offsetIdx2 == offsetIdx2);
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
        sum ^= startIdx1;
        sum *= fnv32Prime;
        sum ^= stopIdx1;
        sum *= fnv32Prime;
        sum ^= offsetIdx2;
        sum *= fnv32Prime;
        sum ^= row;
        sum *= fnv32Prime;
        sum ^= mLen;
        sum *= fnv32Prime;
        sum ^= Double.hashCode(maxChordSum);
        sum *= fnv32Prime;
        return sum;
    }

    @Override
    public int compareTo(SR other) {
        //to handle a changing maxDiffChordSum, will attempt to
        // use the largest on both SR instances
        double maxDiffChord = Math.max(maxChordSum, other.maxChordSum);
        double d1 = calcSalukDist(diffChordSum, maxDiffChord, mLen, nMax);
        double d2 = calcSalukDist(other.diffChordSum, maxDiffChord, 
            other.mLen, other.nMax);
        if (d1 < d2) {
            return -1;
        } else if (d1 > d2) {
            return 1;
        }
        return 0;
    }

    public void setChordSumNeedsUpdate(boolean needsUpdate) {
        chordSumNeedsUpdate = needsUpdate;
    }

    double calcSalukDist() {
        return calcSalukDist(diffChordSum, maxChordSum, 
            (stopIdx1 - startIdx1 + 1), nMax);
    }

    double calcSalukDist(double compChord, double maxChord,
        int length, int maxMatchable) {
        double d = compChord/maxChord;
        double f = 1. - ((double)length/(double)maxMatchable);
        return f*f + d*d;
    }
}
