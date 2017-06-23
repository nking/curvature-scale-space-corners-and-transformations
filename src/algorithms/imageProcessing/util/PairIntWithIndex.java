package algorithms.imageProcessing.util;

import algorithms.util.PairInt;

/**
 *
 * @author nichole
 */
public class PairIntWithIndex extends PairInt {

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

        if (!(obj instanceof PairInt)) {
            return false;
        }

        PairInt other = (PairInt) obj;

        return (getX() == other.getX()) && (getY() == other.getY());
    }

    @Override
    public int hashCode() {

        int hash = fnvHashCode(getX(), getY());

        return hash;
    }

    @Override
    public PairInt copy() {
         return new PairIntWithIndex(getX(), getY(), pixIdx);
    }

    @Override
    public String toString() {

        StringBuilder sb = new StringBuilder(super.toString());
        sb.append(" pixIdx=").append(Integer.toString(pixIdx));

        return sb.toString();
    }

}
