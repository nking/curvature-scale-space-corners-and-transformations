/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package algorithms.strings;

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
        
        StringEditDistance sed = new StringEditDistance();
        int nEdits = sed.calculateWithWagnerFischer(a, b, true);
        
        assertEquals(8, nEdits);
        
    }
}
