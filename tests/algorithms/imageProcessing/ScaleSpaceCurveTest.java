package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class ScaleSpaceCurveTest {
    
    public ScaleSpaceCurveTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
/*
    public ScaleSpaceCurve(float theSigma, PairIntArray curve, 
        boolean curveIsClosed) {
    */
    @Test
    public void test() {
        
        /*
             *  *     3
          *        *  2
             *        1
                *     0
          0  1  2  3 
        */
        PairIntArray xy = new PairIntArray();
        xy.add(0, 2);
        xy.add(1, 1);
        xy.add(2, 0);
        xy.add(3, 2);
        xy.add(2, 3);
        xy.add(0, 3);
        xy.set(5, 1, 3);
        
        ScaleSpaceCurve curve = new ScaleSpaceCurve(0.4f, xy, true);
        
        assertTrue(curve.getSigma() == 0.4f);
        assertTrue(curve.curveIsClosed());
        
        assertTrue(curve.getX(0) == 0 && curve.getY(0) == 2);
        assertTrue(curve.getX(5) == 1 && curve.getY(5) == 3);
        
        assertTrue(curve.getSize() == 6);
        
        assertTrue(curve.getT().length == 6);
        
        curve.setK(0, -0.5f);
        assertTrue(curve.getK(0) == -0.5f);
        
        PairIntArray xyR = curve.getXYCurve();
        assertTrue(xyR.getN() == xy.getN());
        for (int i = 0; i < xyR.getN(); i++) {
            assertTrue(xyR.getX(i) == xy.getX(i));
            assertTrue(xyR.getY(i) == xy.getY(i));
        }
        
        assertTrue(curve.getKIsZeroIdxSize() == 0);
        
        
        curve.setK(2, 0.5f);
        assertTrue(curve.getK(2) == 0.5f);
       
        curve.addKIsZeroIdx(1, xy.getX(1), xy.getY(1));
        curve.compressKIsZeroIdx();
        
        int[] idxes = curve.getKIsZeroIdx();
        assertTrue(idxes.length == 1);
        assertTrue(idxes[0] == 1);
        
        assertTrue(curve.getKIsZeroX().length == 1);
        assertTrue(curve.getKIsZeroY().length == 1);
        assertTrue(curve.getKIsZeroX()[0] == xy.getX(1));
        assertTrue(curve.getKIsZeroY()[0] == xy.getY(1));
    }

}
