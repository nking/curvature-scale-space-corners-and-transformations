package algorithms.strings;

import algorithms.util.PairIntArray;
import java.util.Arrays;
import junit.framework.TestCase;


/**
 *
 * @author nichole
 */
public class StringEditDistanceTest extends TestCase {
    
    public static String reverse(String str) {
        int n = str.length();
        char[] r = new char[n];
        for (int i = 0;  i < n; ++i) {
            r[i] = str.charAt(n-i-1);
        }
        return new String(r);
    }
    
    String a = "rosettacode"; // 6 matches, m=11
    String b = "raisethysword"; // n=13
    int nExpected = 8;
        
    /*
    //          01234567890123456
    //             ... ....     <=== better match
    //                  ...  ....
    String a = "acadefcqdefabqdef";
    //String b =       "zdefgqdef";
    //                   ... ....
    //             ... ....
    String b =   "zdefgqdef"; // <=== better match for having same single gap
    //            012345678
    // 10,8  9,7  8,6  7,5    5,3  4,2  3,1
    int nExpected = 10;
    */
    /*String a = "abcd";
    String b = "zxplmn";
    int nExpected = 6;*/
    
    public StringEditDistanceTest() {
    }
    
    public void testCalculateWithWagnerFischer() {
        
        PairIntArray outIndexes = new PairIntArray();
        StringEditDistance sed = new StringEditDistance();
        int nEdits = sed.calculateWithWagnerFischer(a, b, outIndexes);
        
        System.out.println("nEdits=" + nEdits);
        assertEquals(nExpected, nEdits);
        
        System.out.println("outIndexes=" + outIndexes.toString());
        
        char[] result = new char[outIndexes.getN()];
        for (int i = 0; i < result.length; ++i) {
            result[i] = a.charAt(outIndexes.getX(i));
        }   
        System.out.println("result=" + Arrays.toString(result));
    }
    
    public void testCalculateWithModifiedWagnerFischer() {
        
        System.out.println("testCalculateWithModifiedWagnerFischer");
        
        StringEditDistance sed = new StringEditDistance();
        int nEdits = sed.calculateWithModifiedWagnerFischer(a, b);
        System.out.println("nEdits=" + nEdits);
        assertEquals(nExpected, nEdits);
        
    }
}
