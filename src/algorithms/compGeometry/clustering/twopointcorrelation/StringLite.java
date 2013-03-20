package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 * a fast simple character holder.
 *
 * The hashcode is tailored
 * to hold the following 11 ascii characters:   0-9 ' '
 * and byte[] bytes is expected to hold up to 20 characters.
 *
 * This is for use with TwoPointVoidStats to speed up identity checks.
 *
 * @author nichole
 */
class StringLite {

    protected final byte[] bytes;
    protected int nBytes = 0;

    /**
     * constructor.  value should be less than 33 characters.
     *
     * @param value
     */
    StringLite(byte[] value) {
        this.bytes = value;
        nBytes = value.length;
    }

    public int length() {
        return nBytes;
    }

    /**
     * fast check that contents are the same.  assumes that they are never null
     *
     * @param other
     * @return
     */
    public boolean equals(StringLite other) {
        if (other == null) {
            return false;
        }
        return equals(other.bytes);
    }
    /**
     * fast check that contents are the same. assumes that they are never null
     *
     * @param other
     * @return
     */
    public boolean equals(byte[] other) {
        if (other == null) {
            return false;
        }
        if (bytes.length != other.length) {
            return false;
        }
        for (int i = 0; i < nBytes; i++) {
            if (bytes[i] != other[i]) {
                return false;
            }
        }
        return true;
    }

    @Override
    public boolean equals(Object o) {
        if (o instanceof StringLite) {
            return equals((StringLite)o);
        } else if (o instanceof byte[]) {
            return equals(o);
        } else {
            return false;
        }
    }

    protected int hash = 0;

    @Override
    public int hashCode() {

        // like strings, would like to have same code for same content

        // bytes holds up to 20 characters of 11 possible symbols

        if (bytes == null) {
            return hash;
        }

        return fnvHashCode();
    }

    protected int simpleHashCode() {

        if (hash == 0) {
            int sum = 0;
            for (int i = 0; i < bytes.length; i++) {
                sum += bytes[i];
            }
            hash = sum;
        }

        return hash;
    }

    protected static int fnv321aInit = 0x811c9dc5;
    protected static int fnv32Prime = 0x01000193;

    protected int fnvHashCode() {

        /*
         * hash = offset_basis
         * for each octet_of_data to be hashed
         *     hash = hash xor octet_of_data
         *     hash = hash * FNV_prime
         * return hash
         *
         * Public domain:  http://www.isthe.com/chongo/src/fnv/hash_32a.c
         */

        if (hash == 0) {

            int sum = fnv321aInit;

            for (int i = 0; i < bytes.length; i++) {
                // xor the bottom with the current octet.
	            sum ^= bytes[i];

                // multiply by the 32 bit FNV magic prime mod 2^32
                sum *= fnv32Prime;
            }
            hash = sum;
        }

        return hash;
    }

    static byte[] createIdentity(int index0, int index1) {

        int nB = 8; // 2 integers

        byte[] bytes = new byte[nB];

        System.arraycopy( writeIntegerToBytes(index0, 4), 0, bytes, 0, 4);
        System.arraycopy( writeIntegerToBytes(index1, 4), 0, bytes, 4, 4);

        return bytes;
    }

    /**
     * write the unsigned integer num to 4 bytes
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
            //   in comparing low order digits first
            int index = i;
            bytes[index] = b;
        }

        return bytes;
    }

}
