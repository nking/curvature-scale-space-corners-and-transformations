package algorithms.strings;

import algorithms.util.PairIntArray;
import gnu.trove.list.TCharList;
import java.util.Arrays;
import junit.framework.TestCase;


/**
 *
 * @author nichole
 */
public class StringEditDistanceTest extends TestCase {
    
    String a = "rosettacode";   // 6 matches, m=11
    String b = "raisethysword"; // n=13
    int nExpected = 8;
        
 
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
        //System.out.println("nEdits=" + nEdits);
        assertEquals(nExpected, nEdits);
        
    }
    
    public void testHirschberg() {
        System.out.println("testHirschberg");
        
        //a = "AGTACGCA";
        //b = "TATGC";
        
        StringEditDistance sed = new StringEditDistance();
        TCharList[] xyAligned = sed.hirschbergOptimal(
        //    x.toCharArray(), y.toCharArray());
            a.toCharArray(), b.toCharArray());
        System.out.println("*"+ Arrays.toString(xyAligned[0].toArray()));
        System.out.println("*"+ Arrays.toString(xyAligned[1].toArray()));
        /*
        W = AGTACGCA
        Z = --TATGC-
        */
    }
    
    /*
    String x = "AGTA";
    String y = "ATA";
    
    "AGTACGCA"
    "TATGC"
    
    
    */
}
