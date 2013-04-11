package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.Arrays;

/**
 * a holder for StringLite objects
 *
 * @author nichole
 */
class StringArrayLite {

    protected StringLite[] tstr;
    protected int[] sumTStrChars;
    protected int nTStr = 0;

    StringArrayLite() {
        tstr = new StringLite[100];
        sumTStrChars = new int[100];
        nTStr = 0;
    }

    StringArrayLite(int initialCapacity) {
        tstr = new StringLite[initialCapacity];
        sumTStrChars = new int[initialCapacity];
        nTStr = 0;
    }

    public long approximateMemoryUsed() {

        int n = (tstr == null) ? 0 : tstr.length;

        long sumBytes = 8 + 16 + (n*(16 + 8 + 8));
        long sumBits = (n*32) + 32;

        sumBytes += (sumBits/8);

        // amount of padding needed to make it a round 8 bytes
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
    protected boolean storeIfDoesNotContain(int index0, int index1) {

        // order the indexes to avoid double counting.
        int i0, i1;
        if (index0 < index1) {
            i0 = index0;
            i1 = index1;
        } else {
            i0 = index1;
            i1 = index0;
        }

        byte[] identity = createIdentity(index0, index1);

        int sumChars = sum(identity);

        for (int i = 0; i < nTStr; i++) {
            int tsum = sumTStrChars[i];
            if (sumChars == tsum) {
                if (tstr[i].equals(identity)) {
                    return false;
                }
            }
        }

        expandIfNeeded(identity.length);

        sumTStrChars[nTStr] = sumChars;
        tstr[nTStr] = new StringLite(identity);
        nTStr++;

        return true;
    }

    protected void expandIfNeeded(int nCharacters) {
        if ( (nTStr + 1) > tstr.length) {
            tstr = Arrays.copyOf(tstr, nTStr + 100);
            sumTStrChars = Arrays.copyOf(sumTStrChars, nTStr + 100);
        }
    }

    protected int sum(byte[] bytes) {
        int sum = 0;
        for (int i = 0; i < bytes.length; i++) {
            sum += bytes[i];
        }
        return sum;
    }

     static byte[] createIdentity(int index0, int index1) {

        int nB = 8; // 2 integers

        byte[] bytes = new byte[nB];

        System.arraycopy( writeIntegerToBytes(index0, 4), 0, bytes, 0, 4);
        System.arraycopy( writeIntegerToBytes(index1, 4), 0, bytes, 4, 4);

        return bytes;
    }

    /**
     * write the integer num to 4 bytes
     *
     * @param num
     * @param numBytes
     * @return
     */
    protected static byte[] writeIntegerToBytes(int num, int numBytes) {

        /*
         *  byte          int
         *  -----        -----
         *  0            0
         *  127          127    0111 1111
         * -128          128    1000 0000
         * -1            255    1111 1111
         */

        byte[] bytes = new byte[numBytes];

        for (int i = 0; i < numBytes; i++) {
            int shift = i * 8;
            int a = (num >> shift) & 255;
            byte b = (byte) a;

            // write in reverse order as most numbers will be small, and this results
            //   in comparing low order digits first instead of zeros
            int index = i;
            bytes[index] = b;
        }

        return bytes;
    }

}
