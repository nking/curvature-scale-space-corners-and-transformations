package algorithms.util;

/**
 *
 * @author nichole
 */
public class TrioInt {
    
    private int x = Integer.MIN_VALUE;
    private int y = Integer.MIN_VALUE;
    private int z = Integer.MIN_VALUE;
    
    public TrioInt() {
    }
    public TrioInt(int xPoint, int yPoint, int zPoint) {
        x = xPoint;
        y = yPoint;
        z = zPoint;
    }
    public void setX(int xPoint) {
        x = xPoint;
    }
    public void setY(int yPoint) {
        y = yPoint;
    }
    public void setZ(int zPoint) {
        z = zPoint;
    }
    public int getX() {
        return x;
    }
    public int getY() {
        return y;
    }
    public int getZ() {
        return z;
    }
    
    public TrioInt copy() {
        TrioInt c = new TrioInt(x, y, z);
        return c;
    }

    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof TrioInt)) {
            return false;    
        }
        
        TrioInt other = (TrioInt)obj;
        
        return (x == other.getX()) && (y == other.getY()) && (z == other.getZ());
    }

    @Override
    public int hashCode() {
        
        int hash = fnvHashCode(this.x, this.y, this.z);

        return hash;
    }

    protected static int fnv321aInit = 0x811c9dc5;
    protected static int fnv32Prime = 0x01000193;

    protected int fnvHashCode(int i0, int i1, int i2) {

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
        
        hash = sum;

        return hash;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder("(");
        sb.append(Integer.toString(x)).append(",")
            .append(Integer.toString(y)).append(")")
            .append(Integer.toString(z)).append(")");
        
        return sb.toString();
    }
    
}
