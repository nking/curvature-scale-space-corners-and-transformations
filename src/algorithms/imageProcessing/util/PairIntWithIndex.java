package algorithms.imageProcessing.util;

/**
 *
 * @author nichole
 */
public class PairIntWithIndex extends com.climbwithyourfeet.clustering.util.PairInt {

    int pixIdx;

    public PairIntWithIndex(int xPoint, int yPoint, int thePixIndex) {
        super(xPoint, yPoint);
        pixIdx = thePixIndex;
    }
    
    public int getPixIndex() {
        return pixIdx;
    }

    @Override
    public boolean equals(Object obj) {

        if (!(obj instanceof com.climbwithyourfeet.clustering.util.PairInt)) {
            return false;
        }

        com.climbwithyourfeet.clustering.util.PairInt other
            = (com.climbwithyourfeet.clustering.util.PairInt) obj;

        return (x == other.getX()) && (y == other.getY());
    }

    @Override
    public int hashCode() {

        int hash = fnvHashCode(this.x, this.y);

        return hash;
    }

    @Override
    public com.climbwithyourfeet.clustering.util.PairInt copy() {
         return new PairIntWithIndex(x, y, pixIdx);
    }

    @Override
    public String toString() {

        StringBuilder sb = new StringBuilder(super.toString());
        sb.append(" pixIdx=").append(Integer.toString(pixIdx));

        return sb.toString();
    }

}
