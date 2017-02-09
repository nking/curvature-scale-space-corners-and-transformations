package algorithms.compGeometry;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ParabolaLeastSquaresTest extends TestCase {
    
    public ParabolaLeastSquaresTest() {
    }
    
    public void testSolve() {
        ParabolaLeastSquares pLS = new ParabolaLeastSquares();
        
        pLS.accumulate(0, -4);
        pLS.accumulate(1, -1);
        pLS.accumulate(2, 4);
        pLS.accumulate(3, 11);
        pLS.accumulate(4, 20);
        
        //x    y    x^2    x^3    x^4    xy    yx^2
        assertEquals(10., pLS.moments_[0]);
        assertEquals(30., pLS.moments_[1]);
        assertEquals(120., pLS.moments_[5]);
        assertEquals(30., pLS.moments_[2]);
        assertEquals(434., pLS.moments_[6]);
        assertEquals(100., pLS.moments_[3]);
        assertEquals(354., pLS.moments_[4]);
        
        float[] coeffs = pLS.solve();
        
        assertTrue(Math.abs(coeffs[2] - -4) < 0.01);
        assertTrue(Math.abs(coeffs[1] - 2) < 0.01);
        assertTrue(Math.abs(coeffs[0] - 1) < 0.01);
    }
}
