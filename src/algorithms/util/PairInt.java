package algorithms.util;

/**
 *
 * @author nichole
 */
public class PairInt {
    
    private int x = Integer.MIN_VALUE;
    private int y = Integer.MIN_VALUE;
    
    public PairInt() {
    }
    public PairInt(int xPoint, int yPoint) {
        x = xPoint;
        y = yPoint;
    }
    public void setX(int xPoint) {
        x = xPoint;
    }
    public void setY(int yPoint) {
        y = yPoint;
    }
    public int getX() {
        return x;
    }
    public int getY() {
        return y;
    }

    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof PairInt)) {
            return false;    
        }
        
        PairInt other = (PairInt)obj;
        
        if ((x == other.getX()) && (y == other.getY())) {
            return true;
        }
        
        return false;
    }

    @Override
    public int hashCode() {
        
        int hash = fnvHashCode(this.x, this.y);

        return hash;
    }

    protected static int fnv321aInit = 0x811c9dc5;
    protected static int fnv32Prime = 0x01000193;

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
        
        StringBuilder sb = new StringBuilder("x=");
        sb.append(Integer.toString(x)).append(" y=").append(Integer.toString(y));
        
        return sb.toString();
    }
    
}
