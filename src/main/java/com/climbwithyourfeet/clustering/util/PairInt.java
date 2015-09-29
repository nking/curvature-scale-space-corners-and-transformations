package com.climbwithyourfeet.clustering.util;

/**
 adapted from 
 * http://github.com/nking/curvature-scale-space-corners-and-transformations/blob/master/src/algorithms/util/PairInt.java
 * under MIT License (MIT), Nichole King 2014

 * @author nichole
 */
public class PairInt {
    
    /**
     *
     */
    protected final int x;

    /**
     *
     */
    protected final int y;
    
    /**
     *
     * @param xPoint
     * @param yPoint
     */
    public PairInt(int xPoint, int yPoint) {
        x = xPoint;
        y = yPoint;
    }

    /**
     *
     * @return
     */
    public int getX() {
        return x;
    }

    /**
     *
     * @return
     */
    public int getY() {
        return y;
    }

    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof PairInt)) {
            return false;    
        }
        
        PairInt other = (PairInt)obj;
        
        return (x == other.getX()) && (y == other.getY());
    }

    @Override
    public int hashCode() {
        
        int hash = fnvHashCode(this.x, this.y);

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
    
    /**
     *
     * @return
     */
    public PairInt copy() {
        return new PairInt(x, y);
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder("(");
        sb.append(Integer.toString(x)).append(",").append(Integer.toString(y)).append(")");
        
        return sb.toString();
    }
    
}
