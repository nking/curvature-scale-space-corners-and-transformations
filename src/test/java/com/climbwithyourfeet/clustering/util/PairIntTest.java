package com.climbwithyourfeet.clustering.util;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PairIntTest extends TestCase {
    
    public PairIntTest() {
    }
    
    public void testGetSetX() {
        int xPoint = 3;
        int yPoint = 100;
        PairInt instance = new PairInt(xPoint, yPoint);
        
        assertTrue(instance.getX() == xPoint);
        assertTrue(instance.getY() == yPoint);
    }
    
}
