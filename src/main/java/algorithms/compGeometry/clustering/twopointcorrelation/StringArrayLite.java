package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.Arrays;

/**
 * a holder for StringLite objects
 *
 * @author nichole
 */
class StringArrayLite implements ITwoPointIdentity {

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

    @Override
    public long approximateMemoryUsed() {

        /*
        stack word size is 32 or 64 bits.
        Stack holds:  local variables

        Heap holds:  object references (32 bits) and arrays (32 bits x nItems X itemSize)

         StringLite[] tstr;
         int[] sumTStrChars;
         int nTStr;
                                                     32-bit platform          64-bit platform
                                                Field   Size on          Field   Size on
           Java types                           size    stack            size    stack
           boolean                              32       32              32      64
           byte                                 32       32              32      64
           char                                 32       32              32      64
           short                                32       32              32      64
           int                                  32       32              32      64
           ﬂoat                                 32       32              32      64
           reference                            32       32              64      64
           array reference                      32       32              32      32
           returnAddress                        32       32              64      64
           long                                 64       64              64      128
           double                               64       64              64      128
           from Table II of http://users.elis.ugent.be/~leeckhou/papers/SPE06.pdf
        */

        String arch = System.getProperty("sun.arch.data.model");

        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;


        int n = (tstr == null) ? 0 : tstr.length;

        int nbits = (is32Bit) ? 32 : 64;

        int arrayRefBits = 32;

        int overheadBytes = 16;

        // StringLite size:
        // byte[] bytes  size is 32 *  8 * nbits
        // int nBytes    size is nbits
        long stringLiteSizeInBits = (arrayRefBits * 8 * 32) + nbits;
        long stringLiteSizeInBytes = (stringLiteSizeInBits/8) + overheadBytes;
        long padding = (stringLiteSizeInBytes % 8);
        stringLiteSizeInBytes += padding;


        long sumBits = (arrayRefBits * n * (stringLiteSizeInBytes*8)) + (arrayRefBits * n * 32) + nbits;

        long sumBytes = (sumBits/8) + overheadBytes;

        padding = (sumBytes % 8);

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

        //TODO:  use 2 integer arrays in StringArrayLite instead
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
