package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.Arrays;
import java.util.Locale;

/**
 * a holder for StringLite objects
 *
 * @author nichole
 */
public class StringArrayLite {

    protected StringLite[] tstr;
    protected int[] sumTStrChars;
    protected int nTStr = 0;

    public StringArrayLite() {
        tstr = new StringLite[100];
        sumTStrChars = new int[100];
        nTStr = 0;
    }

    public StringArrayLite(int initialCapacity) {
        tstr = new StringLite[initialCapacity];
        sumTStrChars = new int[initialCapacity];
        nTStr = 0;
    }

    /**
     * check if combination is already stored, if not add it and return true, else
     * return false
     *
     * @param index0
     * @param index1
     * @return
     */
    protected boolean storeIfDoesNotContain(int index0, int index1) {

        // order the indexes to avoid double counting.
        int i0, i1;
        if (index0 < index1) {
            i0 = index0;
            i1 = index1;
        } else {
            i0 = index1;
            i1 = index0;
        }

        char[] strFmted = String.format(Locale.US, "%d %d", i0, i1).toCharArray();

        int sumChars = sumAsciiChars(strFmted);
        for (int i = 0; i < nTStr; i++) {
            int tsum = sumTStrChars[i];
            if (sumChars == tsum) {
                if (tstr[i].equals(strFmted)) {
                    return false;
                }
            }
        }

        expandIfNeeded(strFmted.length);

        sumTStrChars[nTStr] = sumChars;
        tstr[nTStr] = new StringLite(strFmted);
        nTStr++;

        return true;
    }

    protected void expandIfNeeded(int nCharacters) {
        if ( (nTStr + 1) > tstr.length) {
            tstr = Arrays.copyOf(tstr, nTStr + 100);
            sumTStrChars = Arrays.copyOf(sumTStrChars, nTStr + 100);
        }
    }
    
    protected int sumAsciiChars(char[] str) {
        int sum = 0;
        for (int i = 0; i < str.length; i++) {
            sum += str[i];
        }
        return sum;
    }
}
