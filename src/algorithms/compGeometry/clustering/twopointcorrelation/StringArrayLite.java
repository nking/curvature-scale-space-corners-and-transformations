package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.Arrays;

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

    protected boolean containsInTStr(char[] str) {
        int sumChars = sumAsciiChars(str);
        for (int i = 0; i < nTStr; i++) {
            int tsum = sumTStrChars[i];
            if (sumChars == tsum) {
                if (tstr[i].equals(str)) {
                    return true;
                }
            }
        }
        return false;
    }
    protected int sumAsciiChars(char[] str) {
        int sum = 0;
        for (int i = 0; i < str.length; i++) {
            sum += str[i];
        }
        return sum;
    }
    protected void storeInTStr(StringLite str) {
        if ( (nTStr + 1) > tstr.length) {
            tstr = Arrays.copyOf(tstr, nTStr + 100);
            sumTStrChars = Arrays.copyOf(sumTStrChars, nTStr + 100);
        }
        sumTStrChars[nTStr] = sumAsciiChars(str.chars);
        tstr[nTStr] = str;
        nTStr++;
    }
}
