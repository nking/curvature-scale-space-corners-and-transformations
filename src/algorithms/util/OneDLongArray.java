package algorithms.util;

import java.util.Arrays;

/**
 * a class to enclose a primitive int array and provide an
 * equals and hashcode which will have the same identity for two
 * instances that have the same ordered content.
 * @author nichole
 */
public class OneDLongArray {
    public long[] a;
    public OneDLongArray(long[] t) {
        a = t;
    }
    
    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof OneDLongArray)) {
            return false;    
        }
        OneDLongArray other = (OneDLongArray)obj;
        if (a.length != other.a.length) {
            return false;
        }
        
        for (int i = 0; i < a.length; ++i) {
            if (a[i] != other.a[i]) {
                return false;
            }
        }
        
        return true;
    }

    @Override
    public int hashCode() {
        
        int hash = fnvHashCode();

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

        int mask = (1 << 31) - 1;
        int shift = 32;
        int sum = fnv321aInit;
        long b0, b1;
        
        for (int i = 0; i < a.length; ++i) {

            b0 = a[i] & mask;
            b1 = (a[i] >> shift) & mask;
            
            // xor the bottom with the current octet.
            sum ^= b0;
            // multiply by the 32 bit FNV magic prime mod 2^32
            sum *= fnv32Prime;
            
            sum ^= b1;
            sum *= fnv32Prime;
        }
        
        return sum;
    }

    @Override
    public String toString() {
        return Arrays.toString(a);
    }
}
