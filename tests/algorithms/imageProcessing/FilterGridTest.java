package algorithms.imageProcessing;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class FilterGridTest extends TestCase {
    
    public FilterGridTest() {
    }

    public void testFiltergrid() {
        
        int nCols = 8;
        int nRows = 4;
        
        FilterGrid fg = new FilterGrid();
        FilterGrid.FilterGridProducts products = fg.filtergrid(nRows, nCols);
        
        assertNotNull(products);
        
        assertEquals(nRows, products.getRadius().length);
        assertEquals(nCols, products.getRadius()[0].length);
        
        assertEquals(nRows, products.getU1().length);
        assertEquals(nCols, products.getU1()[0].length);
        
        assertEquals(nRows, products.getU2().length);
        assertEquals(nCols, products.getU2()[0].length);
        
        double[][] expected = new double[nRows][];
        expected[0] = new double[]{0., 0.125, 0.25, 0.375, 0.5, 0.375, 0.25, 0.125};
        expected[1] = new double[]{0.25, 0.2795085, 0.35355339, 0.45069391,
            0.55901699, 0.45069391, 0.35355339, 0.2795085};
        expected[2] = new double[]{0.5, 0.5153882, 0.55901699, 0.625, 
            0.70710678, 0.625, 0.55901699, 0.5153882};
        expected[3] = new double[]{0.25, 0.2795085, 0.35355339, 0.45069391,
            0.55901699, 0.45069391, 0.35355339, 0.2795085};
        
        for (int i = 0; i < expected.length; ++i) {
            for (int j = 0; j < expected[i].length; ++j) {
                assertTrue(
                    Math.abs(expected[i][j] - products.getRadius()[i][j])
                    < 0.1);
            }
        }
        
        expected = new double[nRows][];
        expected[0] = new double[]{0., 0.125, 0.25, 0.375, -0.5, -0.375, -0.25, -0.125};
        expected[1] = new double[]{0., 0.125, 0.25, 0.375, -0.5, -0.375, -0.25, -0.125};
        expected[2] = new double[]{0., 0.125, 0.25, 0.375, -0.5, -0.375, -0.25, -0.125};
        expected[3] = new double[]{0., 0.125, 0.25, 0.375, -0.5, -0.375, -0.25, -0.125};
        
        for (int i = 0; i < expected.length; ++i) {
            for (int j = 0; j < expected[i].length; ++j) {
                assertTrue(
                    Math.abs(expected[i][j] - products.getU1()[i][j])
                    < 0.1);
            }
        }
        
        expected = new double[nRows][];
        expected[0] = new double[]{0., 0., 0., 0., 0., 0., 0., 0};
        expected[1] = new double[]{0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
        expected[2] = new double[]{-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5};
        expected[3] = new double[]{-0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25};
        
        for (int i = 0; i < expected.length; ++i) {
            for (int j = 0; j < expected[i].length; ++j) {
                assertTrue(
                    Math.abs(expected[i][j] - products.getU2()[i][j])
                    < 0.1);
            }
        }
    }
    
    public void testPolarFiltergrid() {
        
        int nCols = 8;
        int nRows = 4;
        
        PolarFilterGrid g = new PolarFilterGrid();
        PolarFilterGrid.FilterGridProducts products = g.filtergrid(nRows, nCols);
        
        assertNotNull(products);
        
        assertEquals(nRows, products.getRadius().length);
        assertEquals(nCols, products.getRadius()[0].length);
        
        assertEquals(nRows, products.getCosTheta().length);
        assertEquals(nCols, products.getCosTheta()[0].length);
        
        assertEquals(nRows, products.getSinTheta().length);
        assertEquals(nCols, products.getSinTheta()[0].length);
        
        // assert radius
        double[][] expected = new double[nRows][];
        expected[0] = new double[]{1., 0.125, 0.25, 0.375, 0.5, 0.375, 0.25, 0.125};
        expected[1] = new double[]{0.25, 0.2795085, 0.35355339, 0.45069391,
            0.55901699, 0.45069391, 0.35355339, 0.2795085};
        expected[2] = new double[]{0.5, 0.5153882, 0.55901699, 0.625, 
            0.70710678, 0.625, 0.55901699, 0.5153882};
        expected[3] = new double[]{0.25, 0.2795085, 0.35355339, 0.45069391,
            0.55901699, 0.45069391, 0.35355339, 0.2795085};
        
        for (int i = 0; i < expected.length; ++i) {
            for (int j = 0; j < expected[i].length; ++j) {
                assertTrue(
                    Math.abs(expected[i][j] - products.getRadius()[i][j])
                    < 0.1);
            }
        }
       
        // assert theta
        expected = new double[nRows][];
        expected[0] = new double[]{-0., -0., -0., -0., -3.14159265, -3.14159265,
            -3.14159265, -3.14159265};
        expected[1] = new double[]{-1.57079633, -1.10714872, -0.78539816,
            -0.5880026, -2.67794504, -2.55359005, -2.35619449, -2.03444394};
        expected[2] = new double[]{1.57079633, 1.32581766, 1.10714872, 
            0.92729522, 2.35619449, 2.21429744, 2.03444394,  1.81577499};
        expected[3] = new double[]{1.57079633, 1.10714872, 0.78539816,
            0.5880026, 2.67794504, 2.55359005, 2.35619449, 2.03444394};
        
        for (int i = 0; i < expected.length; ++i) {
            for (int j = 0; j < expected[i].length; ++j) {
                double expTh = expected[i][j];
                double se = Math.sin(expTh);
                double ce = Math.cos(expTh);
                assertTrue(se - products.getSinTheta()[i][j] < 0.1);
                assertTrue(ce - products.getCosTheta()[i][j] < 0.1);
            }
        }
        
    }
    
}
