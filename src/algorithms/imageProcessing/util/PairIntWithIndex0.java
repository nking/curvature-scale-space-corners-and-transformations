package algorithms.imageProcessing.util;

import algorithms.util.PairInt;

/**
 *
 * @author nichole
 */
public class PairIntWithIndex0 extends algorithms.util.PairInt {

    int pixIdx;

    public PairIntWithIndex0(int xPoint, int yPoint, int thePixIndex) {
        super(xPoint, yPoint);
        pixIdx = thePixIndex;
    }
    
    public int getPixIndex() {
        return pixIdx;
    }

    @Override
    public boolean equals(Object obj) {

        if (!(obj instanceof algorithms.util.PairInt)) {
            return false;
        }

        PairInt other = (PairInt) obj;

        return (this.getX() == other.getX()) && (this.getY() == other.getY());
    }

    @Override
    public int hashCode() {

        int hash = fnvHashCode(this.getX(), this.getY());

        return hash;
    }

    @Override
    public PairInt copy() {
         return new PairIntWithIndex0(getX(), getY(), pixIdx);
    }

    @Override
    public String toString() {

        StringBuilder sb = new StringBuilder(super.toString());
        sb.append(" pixIdx=").append(Integer.toString(pixIdx));

        return sb.toString();
    }

}
