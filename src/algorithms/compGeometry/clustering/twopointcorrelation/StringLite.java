package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 * a fast simple character holder.
 *
 * The hashcode is tailored
 * to hold the following 15 ascii characters:   0-9 ' ' '.' 'e' '-' 'f'
 * and char[] chars is expected to hold up to 32 characters.
 *
 * This is for use with TwoPointVoidStats to speed up identity checks.
 *
 * It's use for more than 32 characters is not recommended and use with
 * character sets that extend beyond values 67108863 is not recommended.
 *
 * @author nichole
 */
public class StringLite {

    protected final char[] chars;
    protected int nChars = 0;

    /**
     * constructor.  value should be less than 33 characters.
     *
     * @param value
     */
    public StringLite(char value[]) {
        // value should be less than 33 in length.  not enforced.
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

        // StringLite was created for use in TwoPointVoidStats specifically to hold
        // only 15 symbols: 0-9 ' ' '.' 'e' '-' 'f'
        // and only up to 32 of those symbols.
        //
        // so to hold each combination uniquely within an integer,
        // can divide the integer into 32 slots
        //    0 to 67108863            is slot '0'
        //    67108863 to 134217726    is slot '1'
        //    134217726 to 201326589   is slot '2'
        //    ...
        //    2080374753 to 2147483616 is slot '32'
        // where delta = 67108863 which is 2147483647/32
        //

        if (chars == null) {
            return hash;
        }
        if (hash == 0) {
            int sum = 0;
            for (int i = 0; i < chars.length; i++) {
                sum += chars[i]*67108863;
            }
        }

        return hash;
    }
}
