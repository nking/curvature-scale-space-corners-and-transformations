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
        
        //TODO: revisit this...
        
        int hash = 7;
        hash = 11 * hash + this.x;
        hash = 11 * hash + this.y;
        return hash;
    }
    
}
