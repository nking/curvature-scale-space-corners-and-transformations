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
        if (m1.getColumnDimension() != m2.getRowDimension()) {
            throw new IllegalArgumentException(
                "the number of columns in m1 != number of rows in m2");
        }
        
        int cCols = m2.getColumnDimension();
        int cRows = m1.getRowDimension();
        
        // m1 dot m2
        double[][] m = new double[cRows][cCols];
        for (int i = 0; i < cRows; i++) {
            m[i] = new double[cCols];
        }
        /*
        0     1     0     2  1
        1000  100  10     3  0
                          4  0
        
        0*2    + 1*3   + 0*4     0*1    +  1*0   +  0*0
        1000*2 + 100*3 + 10*4    1000*1 +  100*0 + 10*0
        */
                
        for (int colAdd = 0; colAdd < m2.getColumnDimension(); colAdd++) {
            for (int cRow = 0; cRow < cRows; cRow++) {
                for (int col = 0; col < m1.getColumnDimension(); col++) {
                
                    // a[0][0]  b[0][0]
                    // a[0][1]  b[1][0]
                    // a[0][2]  b[2][0]
                    //
                    // a[1][0]  b[0][0]
                    //
                    // a[1][0]  b[0][0]
                    //
                    double a = m1.get(cRow, col);
                    double b = m2.get(col, colAdd);
                    m[cRow][colAdd] += (a * b);
                }
            }
        }

        return m;
    }
    
}
