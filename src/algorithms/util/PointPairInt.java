package algorithms.util;

/**
 *
 * @author nichole
 */
public class PointPairInt {

    private final PairInt key;
    
    private PairInt value;

    public PointPairInt(PairInt theKey, PairInt theValue) {

        if (theKey == null) {
            throw new IllegalArgumentException("theKey cannot be null");
        }
        if (theValue == null) {
            throw new IllegalArgumentException("theValue cannot be null");
        }

        this.key = theKey;
        this.value = theValue;
    }

    public PairInt getKey() {
        return key;
    }

    public PairInt getValue() {
        return value;
    }

    @Override
    public boolean equals(Object obj) {

        if (!(obj instanceof PointPairInt)) {
            return false;
        }

        PointPairInt other = (PointPairInt) obj;

        return (other.getKey().equals(key) && other.getValue().equals(value));
    }

    @Override
    public int hashCode() {

        int hash = fnvHashCode(this.key, this.value);

        return hash;
    }

    protected static int fnv321aInit = 0x811c9dc5;
    protected static int fnv32Prime = 0x01000193;

    protected int fnvHashCode(PairInt i0, PairInt i1) {

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

        // xor the bottom with the current octet.
        sum ^= i0.getX();

        // multiply by the 32 bit FNV magic prime mod 2^32
        sum *= fnv32Prime;

        sum ^= i0.getY();
        sum *= fnv32Prime;

        sum ^= i1.getX();
        sum *= fnv32Prime;

        sum ^= i1.getY();
        sum *= fnv32Prime;

        hash = sum;

        return hash;
    }
}
