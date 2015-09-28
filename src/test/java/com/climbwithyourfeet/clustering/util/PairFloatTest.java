package com.climbwithyourfeet.clustering.util;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PairFloatTest extends TestCase {
    
    public PairFloatTest(String testName) {
        super(testName);
    }

    public void testGetSetX() {
        float xPoint = 3;
        float yPoint = 100;
        
        PairFloat instance = new PairFloat();
        instance.setX(xPoint);
        instance.setY(yPoint);
        
        assertTrue(instance.getX() == xPoint);
        assertTrue(instance.getY() == yPoint);
    }
    
}
