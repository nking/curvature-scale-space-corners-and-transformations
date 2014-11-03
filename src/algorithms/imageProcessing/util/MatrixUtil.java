package algorithms.imageProcessing.util;

import Jama.Matrix;

/**
 *
 * @author nichole
 */
public class MatrixUtil {
    
    public static double[] multiply(double[][] m, double[] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mcols = m.length;

        int mrows = m[0].length;

        int ncols = n.length;
        
        if (mcols != ncols) {
            throw new IllegalArgumentException(
                "the number of columns in m must equal the number of rows in n");
        }
        
        double[] c = new double[mrows];

        int cCol = 0;
        
        for (int row = 0; row < mrows; row++) {
                        
            for (int col = 0; col < mcols; col++) {
                
                c[cCol] += (m[col][row] * n[col]);
            }
            
            cCol++;        
        }

        return c;
    }
    
    public static double[][] dot(Matrix m1, Matrix m2) {
        
        if (m1 == null) {
            throw new IllegalArgumentException("m1 cannot be null");
        }
        if (m2 == null) {
            throw new IllegalArgumentException("m2 cannot be null");
        }
        if (m1.getArray().length != m2.getArray()[0].length) {
            throw new IllegalArgumentException(
                "the number of columns in m1 != number of rows in m2");
        }
         
        int cCols = m2.getArray().length;
        int cRows = m1.getArray()[0].length;        
        
        // m1 dot m2
        double[][] m = new double[cCols][cRows];
        for (int i = 0; i < m.length; i++) {
            m[i] = new double[cRows];
        }
               
        for (int cRow = 0; cRow < cRows; cRow++) {        
            for (int cCol = 0; cCol < cCols; cCol++) {
                for (int col = 0; col < m1.getArray().length; col++) {
                    m[cCol][cRow] += (m1.get(col, cRow) * m2.get(cCol, col));
                }
            }
        }

        return m;
    }
    
}
