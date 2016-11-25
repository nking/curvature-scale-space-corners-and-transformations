package algorithms.util;

import java.util.Arrays;

/**
 * a class to enclose a primitive int array and provide an
 * equals and hashcode which will have the same identity for two
 * instances that have the same ordered content.
 * @author nichole
 */
public class OneDIntArray {
    public int[] a;
    public OneDIntArray(int[] t) {
        a = t;
    }
    
    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof OneDIntArray)) {
            return false;    
        }
        OneDIntArray other = (OneDIntArray)obj;
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

        int sum = fnv321aInit;
        
        for (int i = 0; i < a.length; ++i) {

            // xor the bottom with the current octet.
            sum ^= a[i];

            // multiply by the 32 bit FNV magic prime mod 2^32
            sum *= fnv32Prime;
        }
        
        return sum;
    }

    @Override
    public String toString() {
        return Arrays.toString(a);
    }
}
