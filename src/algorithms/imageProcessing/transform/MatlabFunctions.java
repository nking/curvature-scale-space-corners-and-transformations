package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import no.uib.cipr.matrix.DenseMatrix;

/**
 *
 * @author nichole
 */
public class MatlabFunctions {
    
    
    /**
     * implements matlab function and parameters for 
     *    repmat(m(row,:), nRowsFactor, nColsFactor)
     * e.g. extract a single row from the m.numRows present in m,
     * then replicate into the output matrix that single row 
     * nRowsFactor times along rows then replicate all of the output
     * columns, nColsFactor times along columns.
     * the result is a matrix of size nRowsFactor X nColsFactor * m.numCols().
     * @param m
     * @param row
     * @param nRowsFactor the factor to use in replicating row 'row' from m.
     * @param nColsFactor the factor to use in replicating 'row' from m along columns.
     * @return matrix of size nRowsFactor X (nColsFactor * m.numCols())
     */
    public static DenseMatrix exRowRepl(DenseMatrix m, int row, int nRowsFactor, int nColsFactor) {

        //repmat(m(3,:),1,3) <--- replicate row 2, 1 time by row, 3 times by columns
        
        int nRows = nRowsFactor;
        int nCols = m.numColumns() * nColsFactor;
        
        DenseMatrix out = new DenseMatrix(nRows, nCols);
        // replicate the row along rows
        for (int c = 0; c < m.numColumns(); ++c) {
            double v = m.get(row, c);
            for (int r = 0; r < nRows; ++r) {
                out.set(r, c, v);
            }
        }
        
        /*  nRowsFactor=2  nColsFactor=3
         0 1 2   0 1 2   0 1 2
         0 1 2
        */
        
        // replicate along columns
        for (int r = 0; r < nRows; ++r) {
            for (int cf = 0; cf < nColsFactor; ++cf) {
                int c2 = m.numColumns() * cf;
                for (int c = 0; c < m.numColumns(); ++c) {
                    double v = out.get(r, c);
                    out.set(r, c2 + c, v);
                }
            }
        }

        return out;
    }
   
    /**
     * implements matlab function and parameters for 
     *    repmat(m(row,:)', nRowsFactor, nColsFactor)
     * e.g. extract a single row from the m.numRows present in m,
     * (this is a 1 X m.numCols vector),  
     * then transpose that vector (result is m.numCols X 1),
     * then replicate into the output matrix that single column 
     * nRowsFactor times along rows then replicate all of the output
     * columns, nColsFactor times along columns.
     * the result is a matrix of size (nRowsFactor * m.numCols) X nColsFactor.
     * @param m
     * @param row
     * @param nRowsFactor the factor to use in replicating row 'row' from m.
     * @param nColsFactor the factor to use in replicating 'row' from m along columns.
     * @return matrix of size (nRowsFactor * m.numCols) X nColsFactor
     */
    public static DenseMatrix exRowConjRepl(DenseMatrix m, int row, int nRowsFactor, int nColsFactor) {
    
    // 0) repmat(X2(3,:)',1,3).   X2 size is 3 X nData.
        //    extract row 2 for all columns.   this is 1 X nData of all "1's"
        //    transpose X2 extract.  result is nData X 1 matrix.
        //    replicate it to new matrix once along rows, and 3 times along
        //    columns.
            
        DenseMatrix at = exRowRepl(m, row, nColsFactor, nRowsFactor);
        
        DenseMatrix out = MatrixUtil.transpose(at);
        
        return out;
    }
    
    public static double[][] piecewiseMult(DenseMatrix m1, DenseMatrix m2) {
        
        if (m1.numRows() != m2.numRows() || m1.numColumns() != m2.numColumns()) {
            throw new IllegalArgumentException("m1 and m2 must be same size");
        }
        
        double[][] out = new double[m1.numRows()][];
        for (int r = 0; r < m1.numRows(); ++r) {
            out[r] = new double[m1.numColumns()];
            for (int c = 0; c < m1.numColumns(); ++c) {
                out[r][c] = m1.get(r, c) * m2.get(r, c);
            }
        }
        
        return out;
    }
    
    /**
     * implement in one part of matlab's function 'sum' to return the sum of 
     * the elements of A along the first array dimension whose size does not 
     * equal 1.  A is a matrix, so sum(A) returns a row vector containing the 
     * sum of each column.
     * @param a
     * @return a 1 X a[0].length array which is the sum of each column of a across
     * all rows.
     */
    public static double[] sumEachColumn(double[][] a) {
        
        double[] out = new double[a[0].length];
        for (int r = 0; r < a.length; ++r) {
            for (int c = 0; c < a[r].length; ++c) {
                out[c] += a[r][c];
            }
        }
        
        return out;
    }
    

}
