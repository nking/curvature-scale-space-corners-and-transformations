package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 * a fast simple character holder.
 *
 * The hashcode is tailored
 * to hold the following 11 ascii characters:   0-9 ' '
 * and char[] chars is expected to hold up to 20 characters.
 *
 * This is for use with TwoPointVoidStats to speed up identity checks.
 *
 * @author nichole
 */
public class StringLite {

    protected final char[] chars;
    protected int nChars = 0;

    /**
     * constructor.  value should be less than 33 characters.
     *
     * @param value
     */
    public StringLite(char value[]) {
        this.chars = value;
        nChars = value.length;
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
        return equals(other.chars);
    }
    /**
     * fast check that contents are the same. assumes that they are never null
     *
     * @param other
     * @return
     */
    public boolean equals(char[] other) {
        if (other == null) {
            return false;
        }
        if (chars.length != other.length) {
            return false;
        }
        for (int i = 0; i < nChars; i++) {
            if (chars[i] != other[i]) {
                return false;
            }
        }
        return true;
    }

    @Override
    public boolean equals(Object o) {
        if (o instanceof StringLite) {
            return equals((StringLite)o);
        } else {
            return false;
        }
    }

    protected int hash = 0;

    @Override
    public int hashCode() {

        // like strings, would like to have same code for same content

        // chars holds up to 20 characters of 11 possible symbols

        if (chars == null) {
            return hash;
        }

        return fnvHashCode();
    }

    protected int simpleHashCode() {

        if (hash == 0) {
            int sum = 0;
            for (int i = 0; i < chars.length; i++) {
                sum += chars[i];
            }
            hash = sum;
        }

        return hash;
    }

    protected int fnv321aInit = 0x811c9dc5;
    protected int fnv32Prime = 0x01000193;

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

            for (int i = 0; i < chars.length; i++) {
                // xor the bottom with the current octet.
                // chars[i] is 16 bits, but StringLite should be holding only
                //      low number ascii characters,
                //      so resulting value of chars[i] is same as extracted for 8 bit
	            sum ^= chars[i];

                // multiply by the 32 bit FNV magic prime mod 2^32
                sum *= fnv32Prime;
            }
            hash = sum;
        }

        return hash;
    }
}
