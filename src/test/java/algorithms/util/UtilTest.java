package algorithms.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;

import algorithms.compGeometry.clustering.twopointcorrelation.AxisIndexer;
import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator;
import algorithms.compGeometry.clustering.twopointcorrelation.TwoPointCorrelation;

import junit.framework.TestCase;

public class UtilTest extends TestCase {
    
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    public void test() throws Exception {
        
        long availMem = Util.getAvailableHeapMemory();
        assertTrue(availMem > 0);
        
        // for coverage:
        Util u = new Util();
    }
}
