package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 * a fast simple character holder for use with TwoPointVoidStats to speed up identity checks.
 *
 * @author nichole
 */
class StringLite {

    protected final byte[] bytes;
    protected int nBytes = 0;

    /**
     * constructor
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

        if (bytes == null) {
            return hash;
        }

        return fnvHashCode();
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

}
