package com.climbwithyourfeet.clustering.util;

/**
 adapted from 
 * http://github.com/nking/curvature-scale-space-corners-and-transformations/blob/master/src/algorithms/util/PairFloat.java
 * under MIT License (MIT), Nichole King 2014

 * @author nichole
 */
public class PairFloat {
    
    /**
     *
     */
    protected final float x;

    /**
     *
     */
    protected final float y;
    
    /**
     *
     * @param xPoint
     * @param yPoint
     */
    public PairFloat(float xPoint, float yPoint) {
        x = xPoint;
        y = yPoint;
    }

    /**
     *
     * @return
     */
    public float getX() {
        return x;
    }

    /**
     *
     * @return
     */
    public float getY() {
        return y;
    }

    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof PairFloat)) {
            return false;    
        }
        
        PairFloat other = (PairFloat)obj;
        
        return (x == other.getX()) && (y == other.getY());
    }

    @Override
    public int hashCode() {
        
        int hash = fnvHashCode(Float.floatToIntBits(this.x), 
            Float.floatToIntBits(this.y));

        return hash;
    }

    /**
     *
     */
    protected static int fnv321aInit = 0x811c9dc5;

    /**
     *
     */
    protected static int fnv32Prime = 0x01000193;

    /**
     *
     * @param i0
     * @param i1
     * @return
     */
    protected int fnvHashCode(int i0, int i1) {

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
        sum ^= i0;

        // multiply by the 32 bit FNV magic prime mod 2^32
        sum *= fnv32Prime;
        
        sum ^= i1;
        
        sum *= fnv32Prime;
        
        hash = sum;

        return hash;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder("(");
        sb.append(Float.toString(x)).append(",").append(Float.toString(y)).append(")");
        
        return sb.toString();
    }
    
}