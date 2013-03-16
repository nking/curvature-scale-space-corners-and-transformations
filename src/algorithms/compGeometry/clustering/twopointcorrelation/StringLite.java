package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 * a fast simple character holder
 *
 * @author nichole
 */
public class StringLite {

    protected final char[] chars;
    protected int nChars = 0;

    public StringLite(char value[]) {
        this.chars = value;
        nChars = value.length;
    }

    /**
     * fast check that contents are the same.  assumes that they are never null
     *
     * @param other
     * @return
     */
    public boolean equals(StringLite other) {
        if (other == null) {
            return false;
        }
        return equals(other.chars);
    }
    /**
     * fast check that contents are the same. assumes that they are never null
     *
     * @param other
     * @return
     */
    public boolean equals(char[] other) {
        if (other == null) {
            return false;
        }
        if (chars.length != other.length) {
            return false;
        }
        for (int i = 0; i < nChars; i++) {
            if (chars[i] != other[i]) {
                return false;
            }
        }
        return true;
    }

    @Override
    public boolean equals(Object o) {
        if (o instanceof StringLite) {
            return equals((StringLite)o);
        } else {
            return false;
        }
    }

    protected int hash = 0;

    @Override
    public int hashCode() {
        // like strings, would like to have same code for same content

        if (chars == null) {
            return hash;
        }
        if (hash == 0) {
            int sum = 0;
            for (int i = 0; i < chars.length; i++) {
                sum += chars[i]*11;
            }
            hash = sum;
        }
        return hash;
    }

}
