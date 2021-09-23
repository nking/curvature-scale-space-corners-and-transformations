/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package algorithms.strings;

import algorithms.util.PairIntArray;
import junit.framework.TestCase;


/**
 *
 * @author nichole
 */
public class StringEditDistanceTest extends TestCase {
    
    public StringEditDistanceTest() {
    }
    
    public void testCalculateWithWagnerFischer() {
        String a = "rosettacode";
        String b = "raisethysword";
        
        PairIntArray outIndexes = new PairIntArray();
        StringEditDistance sed = new StringEditDistance();
        int nEdits = sed.calculateWithWagnerFischer(a, b, outIndexes);
        
        assertEquals(8, nEdits);
        
        char[] result = new char[outIndexes.getN()];
        for (int i = 0; i < result.length; ++i) {
            result[i] = a.charAt(outIndexes.getX(i));
        } 
        
        
    }
}
