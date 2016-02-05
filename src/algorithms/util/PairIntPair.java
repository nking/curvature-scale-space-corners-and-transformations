package algorithms.util;

/**
 *
 * @author nichole
 */
public class PairIntPair {
    
    private final int p1X;
    
    private final int p1Y;
    
    private final int p2X;
    
    private final int p2Y;
    
    public PairIntPair (int x1, int y1, int x2, int y2) {
        this.p1X = x1;
        this.p1Y = y1;
        this.p2X = x2;
        this.p2Y = y2;
    }
    
    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof PairIntPair)) {
            return false;    
        }
        
        PairIntPair other = (PairIntPair)obj;
        
        if ((p1X == other.p1X) && (p1Y == other.p1Y) && (p2X == other.p2X) && 
            (p2Y == other.p2Y)) {
            return true;
        }
        
        if ((p1X == other.p2X) && (p1Y == other.p2Y) && (p2X == other.p1X) && 
            (p2Y == other.p1Y)) {
            return true;
        }
                
        return false;
    }
    
    @Override
    public int hashCode() {
        
        int x1, y1, x2, y2;
        
        if (p1X < p2X) {
            x1 = p1X;
            y1 = p1Y;
            x2 = p2X;
            y2 = p2Y;
        } else {
            x1 = p2X;
            y1 = p2Y;
            x2 = p1X;
            y2 = p1Y;
        }
        
        int hash = fnvHashCode(x1, y1, x2, y2);

        return hash;
    }

    private static int fnv321aInit = 0x811c9dc5;
    private static int fnv32Prime = 0x01000193;

    protected int fnvHashCode(int i0, int i1, int i2, int i3) {

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
        
        sum ^= i2;
        
        sum *= fnv32Prime;
        
        sum ^= i3;
        
        sum *= fnv32Prime;
        
        hash = sum;

        return hash;
    }

    /**
     * @return the x1
     */
    public int getX1() {
        return p1X;
    }

    /**
     * @return the y1
     */
    public int getY1() {
        return p1Y;
    }

    /**
     * @return the x2
     */
    public int getX2() {
        return p2X;
    }

    /**
     * @return the y2
     */
    public int getY2() {
        return p2Y;
    }
}
