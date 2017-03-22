package algorithms.imageProcessing.util;

import algorithms.imageProcessing.ImageProcessor;
import algorithms.misc.Complex;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.TDoubleSet;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TDoubleHashSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.VectorEntry;
import no.uib.cipr.matrix.sparse.FlexCompColMatrix;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;
import no.uib.cipr.matrix.sparse.SparseVector;

/**
 *
 * @author nichole
 */
public class MatrixUtil {
    
    public static double[] add(double[] m, double[] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        if (m.length != n.length) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        
        int len = m.length;
     
        double[] c = new double[len];
        
        for (int i = 0; i < len; i++) {
            c[i] = m[i] + n[i];
        }

        return c;
    }
    
    public static float[] add(float[] m, float[] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        if (m.length != n.length) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        
        int len = m.length;
     
        float[] c = new float[len];
        
        for (int i = 0; i < len; i++) {
            c[i] = m[i] + n[i];
        }

        return c;
    }
    
    public static float[][] subtract(float[][] m, float[][] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        if (m.length != n.length) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        if (m[0].length != n[0].length) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        
        float[][] c = new float[m.length][];

        for (int i = 0; i < m.length; ++i) {
            c[i] = new float[m[0].length];
            for (int j = 0; j < m[0].length; ++j) {
                c[i][j] -= m[i][j] - n[i][j];
            }
        }

        return c;
    }
    
    public static float[][] add(float[][] m, float[][] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        if (m.length != n.length) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        if (m[0].length != n[0].length) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        
        float[][] c = new float[m.length][];

        for (int i = 0; i < m.length; ++i) {
            c[i] = new float[m[0].length];
            for (int j = 0; j < m[0].length; ++j) {
                c[i][j] = m[i][j] + n[i][j];
            }
        }

        return c;
    }

    public static float[] subtract(float[] m, float[] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        if (m.length != n.length) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        
        int len = m.length;
     
        float[] c = new float[len];
        
        for (int i = 0; i < len; i++) {
            c[i] = m[i] - n[i];
        }

        return c;
    }
    
    public static DenseMatrix subtract(DenseMatrix m, DenseMatrix n) {

        if (m == null) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        if (m.numRows() != n.numRows() || m.numColumns() != n.numColumns()) {
            throw new IllegalArgumentException("m and n must be same length");
        }
        
        DenseMatrix output = new DenseMatrix(m.numRows(), m.numColumns());
        
        for (int i = 0; i < m.numRows(); ++i) {
            for (int j = 0; j < m.numColumns(); ++j) {
                double v0 = m.get(i, j);
                double v1 = n.get(i, j);
                output.set(i, j, v0 - v1);
            }
        }
        
        return output;
    }
    
    public static void add(int[] m, int n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int len = m.length;
             
        for (int i = 0; i < len; i++) {
            m[i] += n;
        }
    }
    
    /**
     * performs eigenvalue decomposition, checks the matrix products and returns
     * those real vectors that pass the real product tests.
     * Also returns their eigenvalues.
     * new Object[]{rightEigenVector, eigenValues] as
     * DenseMatrix and TDoubleList.
     * Note that the results are also filtered for uniqueness.
     * @param m
     * @return 
     */
    public static Object[] eigenWithErrorFilter(DenseMatrix m) {
        
        EVD evd;
        try {
            evd = EVD.factorize(m);
        } catch (NotConvergedException ex) {
            Logger.getLogger(MatrixUtil.class.getName()).log(
                Level.SEVERE, null, ex);
            return null;
        }
        
        DenseMatrix rightEigenVectors = evd.getRightEigenvectors();
        
        //System.out.println("a=\n" + m.toString());
        //System.out.println("left e=\n" + evd.getLeftEigenvectors().toString());
        //System.out.println("right e=\n" + rightEigenVectors.toString());
        
        //System.out.println("leftT * M * right=\n" +
        //    MatrixUtil.multiply(leftEigenVectors.transpose(),
        //        MatrixUtil.multiply(m, rightEigenVectors)));
        
        double[] eigenValues = evd.getRealEigenvalues();
        
        DenseMatrix d = new DenseMatrix(eigenValues.length, eigenValues.length);
        for (int i = 0; i < eigenValues.length; ++i) {
            d.set(i, i, eigenValues[i]);
        }
        
        /*
        The right eigenvector v(j) of A satisfies
                          A * v(j) = lambda(j) * v(j)
        where lambda(j) is its eigenvalue.
        
        NOTE: to create a method which returns the imaginary portion of a
        complex solution, the left check is:
        The left eigenvector u(j) of A satisfies
                       u(j)**H * A = lambda(j) * u(j)**H
        where u(j)**H denotes the conjugate-transpose of u(j).
        u is the left eigen vector and lambda is the eigen values.
        */
        // A * V
        DenseMatrix check0_right = MatrixUtil.multiply(m, rightEigenVectors);
        // D * V
        DenseMatrix check1_right = MatrixUtil.multiply(rightEigenVectors, d);
                
        //System.out.println("check0_right=\n" + check0_right);
        //System.out.println("check1_right=\n" + check1_right);

        // get columns of vectors passing product test
        TIntSet columns = check(check0_right, check1_right);
        
        // filter for unique eigen vectors and values
        TIntSet rm = new TIntHashSet();
        TDoubleSet exists = new TDoubleHashSet();
        TIntIterator iter = columns.iterator();
        while (iter.hasNext()) {
            int col = iter.next();
            double ev = eigenValues[col];
            if (exists.contains(ev)) {
                rm.add(col);
            } else {
                exists.add(ev);
            }
        }
        columns.removeAll(rm);
        
        DenseMatrix rightVectors2 = new DenseMatrix(rightEigenVectors.numRows(),
            columns.size());
        TDoubleList eigenValues2 = new TDoubleArrayList(columns.size());
        int col2 = 0;
        
        iter = columns.iterator();
        while (iter.hasNext()) {
            int col = iter.next();            
            eigenValues2.add(eigenValues[col]);
            
            for (int row = 0; row < rightEigenVectors.numRows(); ++row) {
                rightVectors2.set(row, col2, 
                    rightEigenVectors.get(row, col));
            }
            
            col2++;
        }
        
        return new Object[]{rightVectors2, eigenValues2};
    }
    
    // this is adapted from JAMA matrix test http://math.nist.gov/javanumerics/jama/
    // Jama-1.0.3.zip
    private static TIntSet check(DenseMatrix m, DenseMatrix n) {
        
        TIntSet columns = new TIntHashSet();
            
        double eps = Math.pow(2.0, -52.0);
        
        if (m.norm(Matrix.Norm.Frobenius) == 0. && 
            n.norm(Matrix.Norm.Frobenius) < 10 * eps) {
            return columns;
        }
        if (n.norm(Matrix.Norm.Frobenius) == 0. && 
            m.norm(Matrix.Norm.Frobenius) < 10 * eps) {
            return columns;
        }
        
        double mNorm = Math.max(m.norm(Matrix.Norm.One), 
            n.norm(Matrix.Norm.Frobenius));
        double c2 = 1000 * eps * mNorm;
        
        DenseVector d = new DenseVector(m.numRows());
        
        for (int j = 0; j < m.numColumns(); ++j) {
            for (int i = 0; i < m.numRows(); ++i) {
                double v0 = m.get(i, j);
                double v1 = n.get(i, j);
                double diff = v0 - v1;
                d.set(i, diff);
            }
            double c = d.norm(Vector.Norm.One);            
            if (c < c2) {
                columns.add(j);
            }
        }
        
        return columns;
    }
    
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
    
    public static float[] multiply(float[][] m, float[] n) {

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
        
        float[] c = new float[mrows];

        int cCol = 0;
        
        for (int row = 0; row < mrows; row++) {
                        
            for (int col = 0; col < mcols; col++) {
                
                c[cCol] += (m[col][row] * n[col]);
            }
            
            cCol++;        
        }

        return c;
    }
    
    public static void multiply(int[] m, int factor) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int len = m.length;
                
        for (int i = 0; i < len; i++) {                        
            m[i] *= factor;
        }
    }
    
    
    public static void multiply(float[][] a, float m) {

        if (a == null || a.length == 0) {
            throw new IllegalArgumentException("a cannot be null or empty");
        }
        
        int mcols = a.length;

        int mrows = a[0].length;
        
        for (int col = 0; col < mcols; col++) {
            for (int row = 0; row < mrows; row++) {
                a[col][row] *= m;
            }            
        }
    }
    
    public static void multiply(int[] m, int[] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        if (m.length != n.length) {
            throw new IllegalArgumentException("m must be the same size as n");
        }
        
        int len = m.length;
                
        for (int i = 0; i < len; i++) {                        
            m[i] *= n[i];
        }
    }
    
    public static double[][] multiply(double[][] m, double[][] n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mrows = m.length;

        int mcols = m[0].length;

        int nrows = n.length;
        
        int ncols = n[0].length;
        
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of columns in m must equal the number of rows in n");
        }
        
        /*
        a b c      p0 p1 p2
        d e f      p3 p4 p5
                   p6 p7 p8        
        a*p0+... a*p a*p
        d*p0+... d*p d*p
        */
        
        double[][] c = new double[mrows][];
        
        for (int row = 0; row < mrows; row++) {
            c[row] = new double[ncols];
            for (int ncol = 0; ncol < ncols; ncol++) {
                double sum = 0;                
                for (int mcol = 0; mcol < mcols; mcol++) {
                    sum += (m[row][mcol] * n[mcol][ncol]);                    
                }
                c[row][ncol] = sum;
            }            
        }

        return c;
    }
    
    public static DenseMatrix multiply(
        Matrix m, Matrix n) {

        if (m == null || m.numRows() == 0 || m.numColumns() == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.numRows() == 0 || n.numColumns() == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mrows = m.numRows();

        int mcols = m.numColumns();

        int nrows = n.numRows();
        
        int ncols = n.numColumns();
        
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of columns in m must equal the number of rows in n");
        }
        
        /*
        a b c      p0 p1 p2
        d e f      p3 p4 p5
                   p6 p7 p8        
        a*p0+... a*p a*p
        d*p0+... d*p d*p
        */
        
        no.uib.cipr.matrix.DenseMatrix c = new DenseMatrix(mrows, ncols);
        
        for (int row = 0; row < mrows; row++) {
            for (int ncol = 0; ncol < ncols; ncol++) {
                double sum = 0;                
                for (int mcol = 0; mcol < mcols; mcol++) {
                    sum += (m.get(row, mcol) * n.get(mcol, ncol));                    
                }
                c.set(row, ncol, sum);
            }            
        }

        return c;
    }
    
    public static long countNodes(FlexCompRowMatrix a) {
        
        long count = 0;
        Iterator<MatrixEntry> iter = a.iterator();
        while (iter.hasNext()) {
            iter.next();
            ++count;
        }
        return count;
    }

    public static long countNodes(FlexCompColMatrix a) {
        
        long count = 0;
        Iterator<MatrixEntry> iter = a.iterator();
        while (iter.hasNext()) {
            iter.next();
            ++count;
        }
        return count;
    }
    
    /**
     * for sparse matrices m and n (which are same size), subtract them sparsely
     * and return the sparse result.
     * @param m
     * @param n
     * @return 
     */
    public static FlexCompRowMatrix sparseMatrixSubtract(FlexCompRowMatrix m, 
        FlexCompRowMatrix n) {

        if (m == null || m.numRows() == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.numRows() == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mrows = m.numRows();

        int mcols = m.numColumns();

        if (mrows != n.numRows() || mcols != n.numColumns()) {
            throw new IllegalArgumentException(
                "m and n must be the same size");
        }
        
        FlexCompRowMatrix c = new FlexCompRowMatrix(mrows, mcols);
        
        Set<PairInt> processed = new HashSet<PairInt>();
        
        long count = 0;
        
        Iterator<MatrixEntry> iter = m.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            int col = entry.column();
            int row = entry.row();
            double v = entry.get();
            
            double v2 = n.get(row, col);
            
            double result = v - v2;
            
            if (result != 0) {
            
                c.set(row, col, v - v2);
                
                count++;
            }
            
            processed.add(new PairInt(row, col));
           
        }
        
        iter = n.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            int col = entry.column();
            int row = entry.row();
            
            PairInt p = new PairInt(row, col);
            if (processed.contains(p)) {
                continue;
            }
            
            double v2 = entry.get();
            
            double v = m.get(row, col);
            
            double result = v - v2;
            
            if (result != 0) {
            
                c.set(row, col, v - v2);
                
                count++;
            }
            
            //processed.add(p);
        }
        
        //System.out.println(count + " nodes in output sparse matrix");
        
        return c;
    }
    
    /**
     * for sparse matrices m and n (which are same size), subtract them sparsely
     * and return the sparse result.
     * @param m
     * @param n
     * @return 
     */
    public static FlexCompColMatrix sparseMatrixSubtract(FlexCompColMatrix m, 
        FlexCompColMatrix n) {
        
        //TODO: consider use parameterization here with FlexCompRowMatrix method

        if (m == null || m.numRows() == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.numRows() == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mrows = m.numRows();

        int mcols = m.numColumns();

        if (mrows != n.numRows() || mcols != n.numColumns()) {
            throw new IllegalArgumentException(
                "m and n must be the same size");
        }
        
        FlexCompColMatrix c = new FlexCompColMatrix(mrows, mcols);
        
        Set<PairInt> processed = new HashSet<PairInt>();
        
        long count = 0;
        
        Iterator<MatrixEntry> iter = m.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            int col = entry.column();
            int row = entry.row();
            double v = entry.get();
            
            double v2 = n.get(row, col);
            
            double result = v - v2;
            
            if (result != 0) {
            
                c.set(row, col, v - v2);
                
                count++;
            }
            
            processed.add(new PairInt(row, col));
           
        }
        
        iter = n.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            int col = entry.column();
            int row = entry.row();
            
            PairInt p = new PairInt(row, col);
            if (processed.contains(p)) {
                continue;
            }
            
            double v2 = entry.get();
            
            double v = m.get(row, col);
            
            double result = v - v2;
            
            if (result != 0) {
            
                c.set(row, col, v - v2);
                
                count++;
            }
            
            //processed.add(p);
        }
        
        //System.out.println(count + " nodes in output sparse matrix");
        
        return c;
    }
    
    /**
     * multiply m by n and returns result in a sparse matrix.
     * iterates over nearly full size of nRows X nCols, but only writes non-zero values
     * to a sparse matrix result.
     * @param m
     * @param n
     * @return 
     */
    public static FlexCompRowMatrix sparseMatrixMultiply(
        FlexCompRowMatrix m, FlexCompRowMatrix n) {

        if (m == null || m.numRows() == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.numRows() == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mrows = m.numRows();

        int mcols = m.numColumns();

        int nrows = n.numRows();
        
        int ncols = n.numColumns();
        
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of columns in m must equal the number of rows in n");
        }
        
        /*
        a b c      p0 p1 p2
        d e f      p3 p4 p5
                   p6 p7 p8        
        a*p0+... a*p a*p
        d*p0+... d*p d*p
        */
        
        long count = 0;
        
        Iterator<VectorEntry> iter = null;
        
        FlexCompRowMatrix c = new FlexCompRowMatrix(mrows, ncols);
                        
        for (int row = 0; row < mrows; row++) {
            
            for (int ncol = 0; ncol < ncols; ncol++) {
                
                double sum = 0;
                
                iter = m.getRow(row).iterator();
                while (iter.hasNext()) {
                    VectorEntry entry = iter.next();
                    int mcol = entry.index();
                    double vm = entry.get();
                    double vn = n.get(mcol, ncol);
                    sum += (vm * vn);
                }
                //for (int mcol = 0; mcol < mcols; mcol++) {
                //    sum += (m.get(row, mcol) * n.get(mcol, ncol));                    
                //}
                if (sum != 0) {
                    c.set(row, ncol, sum);
                    count++;
                }
            } 
        }

        //System.out.println(count + " nodes in output sparse matrix");
        
        return c;
    }
    
     /**
     * multiply m by n and returns result in a sparse matrix.
     * iterates over nearly full size of nRows X nCols, but only writes non-zero values
     * to a sparse matrix result.
     * @param m
     * @param n
     * @return 
     */
    public static FlexCompColMatrix sparseMatrixMultiply(
        FlexCompColMatrix m, FlexCompColMatrix n) {

        if (m == null || m.numRows() == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (n == null || n.numRows() == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mrows = m.numRows();

        int mcols = m.numColumns();

        int nrows = n.numRows();
        
        int ncols = n.numColumns();
        
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of columns in m must equal the number of rows in n");
        }
        
        /*
        a b c      p0 p1 p2
        d e f      p3 p4 p5
                   p6 p7 p8        
        a*p0+b*p3  ...  a*p1    a*p2
        d*p0+e*p3  ...  d*p1    d*p2
        */
        
        long count = 0;
        
        Iterator<VectorEntry> iter = null;
        
        FlexCompColMatrix c = new FlexCompColMatrix(mrows, ncols);
                            
        for (int ncol = 0; ncol < ncols; ncol++) {
        
            SparseVector nColVec = n.getColumn(ncol);
            
            for (int row = 0; row < mrows; row++) {
                
                double sum = 0;
                
                int colCount = 0;
                iter = nColVec.iterator();
                while (iter.hasNext()) {
                    VectorEntry entry = iter.next();
                    int mcol = entry.index();
                    double vn = entry.get();
                    double vm = m.get(row, colCount);
                    sum += (vm * vn);
                    colCount++;
                }
                if (sum != 0) {
                    c.set(row, ncol, sum);
                    count++;
                }
            } 
        }

        //System.out.println(count + " nodes in output sparse matrix");
        
        return c;
    }
    
    public static double[][] convertToRowMajor(DenseMatrix a) {
        int nc = a.numColumns();
        int nr = a.numRows();
        double[][] out = new double[nr][];
        for (int i = 0; i < nr; ++i) {
            out[i] = new double[nc];
            for (int j = 0; j < nc; ++j) {
                out[i][j] = a.get(i, j);
            }
        }
        return out;
    }
    
    public static double[][] dot(DenseMatrix m1, DenseMatrix m2) {
        
        if (m1 == null) {
            throw new IllegalArgumentException("m1 cannot be null");
        }
        if (m2 == null) {
            throw new IllegalArgumentException("m2 cannot be null");
        }
        int cCols = m2.numColumns();
        int cRows = m1.numRows();
        
        if (m1.numColumns() != m2.numRows()) {
            throw new IllegalArgumentException(
                "the number of columns in m1 != number of rows in m2");
        }
        
        // m1 dot m2
        double[][] m = new double[cRows][cCols];
        for (int i = 0; i < cRows; i++) {
            m[i] = new double[cCols];
        }
        
        /*
        t00  t01  t02       x1  x2  x3  x4
        t10  t11  t12       y1  y2  y3  y4
        t20  t21  t22       1    1   1   1
        
        row=0, col=0:nCols0  times and plus col=0, row=0:nRows1 --> stored in row, row + (cAdd=0)
        row=1, col=0:nCols0  times and plus col=0, row=0:nRows1 --> stored in row, row + (cAdd=0)
                
        row=0, col=0:nCols0  times and plus col=(cAdd=1), row=0:nRows1 --> stored in row, row + (cAdd=0)
        */
        
        for (int colAdd = 0; colAdd < m2.numColumns(); colAdd++) {
            for (int row = 0; row < m1.numRows(); ++row) {
                for (int col = 0; col < m1.numColumns(); col++) {
                    double a = m1.get(row, col);
                    double b = m2.get(col, colAdd);
                    m[row][colAdd] += (a * b);
                }
            }
        }

        return m;
    }
        
    /**
     * apply dot operator to m1 and m2 which are formatted using same as 
     * DenseMatrix, that is row major [row][col].
     * @param m1
     * @param m2
     * @return 
     */
    public static double[][] dot(double[][] m1, double[][] m2) {
        
        if (m1 == null) {
            throw new IllegalArgumentException("m1 cannot be null");
        }
        if (m2 == null) {
            throw new IllegalArgumentException("m2 cannot be null");
        }
        int cCols = m2[0].length;
        int cRows = m1.length;
        
        if (m1[0].length != m2.length) {
            throw new IllegalArgumentException(
                "the number of columns in m1 != number of rows in m2");
        }
        
        // m1 dot m2
        double[][] m = new double[cRows][cCols];
        for (int i = 0; i < cRows; i++) {
            m[i] = new double[cCols];
        }
        
        for (int colAdd = 0; colAdd < m2[0].length; colAdd++) {
            for (int row = 0; row < m1.length; ++row) {
                for (int col = 0; col < m1[0].length; col++) {
                    double a = m1[row][col];
                    double b = m2[col][colAdd];
                    m[row][colAdd] += (a * b);
                }
            }
        }

        return m;
    }
    
    public static float[][] transpose(float[][] m) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int mRows = m.length;
        int mCols = m[0].length;
        
        float[][] t = new float[mCols][];
        for (int i = 0; i < mCols; i++) {
            t[i] = new float[mRows];
        }
        
        for (int i = 0; i < mRows; i++) {
            for (int j = 0; j < mCols; j++) {
                t[j][i] = m[i][j];
            }
        }
        
        return t;
    }
    
    public static DenseMatrix transpose(DenseMatrix m) {
        
        int mRows = m.numRows();
        int mCols = m.numColumns();
        
        double[][] t = new double[mRows][];
        for (int i = 0; i < mRows; i++) {
            t[i] = new double[mCols];
            for (int j = 0; j < mCols; j++) {
                t[i][j] = m.get(i, j);
            }
        } 
        
        double[][] transposed = transpose(t);
        
        DenseMatrix mT = new DenseMatrix(transposed);
        
        return mT;
    }
    
    public static double[][] transpose(double[][] m) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int mRows = m.length;
        int mCols = m[0].length;
        
        double[][] t = new double[mCols][];
        for (int i = 0; i < mCols; i++) {
            t[i] = new double[mRows];
        }
        
        for (int i = 0; i < mRows; i++) {
            for (int j = 0; j < mCols; j++) {
                t[j][i] = m[i][j];
            }
        }
        
        return t;
    }
    
    public static Complex[][] swapXandY(Complex[][] m) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int mRows = m.length;
        int mCols = m[0].length;
        
        Complex[][] t = new Complex[mCols][];
        for (int i = 0; i < mCols; i++) {
            t[i] = new Complex[mRows];
        }
        
        for (int i = 0; i < mRows; i++) {
            for (int j = 0; j < mCols; j++) {
                t[j][i] = m[i][j];
            }
        }
        
        return t;
    }
    
    /**
     * multiply matrix m by the transpose of n
     * @param p
     * @param n
     * @return 
     */
    public static double[][] multiplyByTranspose(double[][] p, double[][] n) {

        if (p == null || p.length == 0) {
            throw new IllegalArgumentException("p cannot be null or empty");
        }
        if (n == null || n.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int prows = p.length;

        int pcols = p[0].length;

        int nrows = n.length;
        
        int ncols = n[0].length;
        
        if (pcols != ncols) {
            throw new IllegalArgumentException(
                "the number of columns in p must equal the number of cols in n "
                + "for multiplication of m by transpose of n");
        }
        
        /*
        example:  m is p0 p1 p2
                       p3 p4 p5
                       p6 p7 p8
        
                  n is x1  y1  1
                       x2  y2  1
        
        multiply m by transpose of n:
        
        p0 p1 p2     x1  y2
        p3 p4 p5     y1  y2
        p6 p7 p8      1   1
        
        (p0*x1 + p1*y1 + p2*1)  (p0*x2 + p1*y2 + p2*1)
        (p3*x1 + p4*y1 + p5*1)  (p3*x2 + p4*y2 + p5*1)
        (p6*x1 + p7*y1 + p8*1)  (p6*x2 + p7*y2 + p8*1)
        */
        
        double[][] c = new double[prows][];
       
        for (int row = 0; row < prows; row++) {
            c[row] = new double[nrows];
            for (int nrow = 0; nrow < nrows; nrow++) {
                double sum = 0;                
                for (int mcol = 0; mcol < pcols; mcol++) {
                    sum += (p[row][mcol] * n[nrow][mcol]);                    
                }
                c[row][nrow] = sum;
            }            
        }

        return c;
    }
    

    public static double[] multiply(Matrix a, double[] b) {
        
        if (a == null || a.numRows() == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (b == null || b.length == 0) {
            throw new IllegalArgumentException("n cannot be null or empty");
        }
        
        int mrows = a.numRows();

        int mcols = a.numColumns();

        int nrows = b.length;
        
        if (mcols != nrows) {
            throw new IllegalArgumentException(
                "the number of cols in a must equal the length of b");
        }
        
        double[] c = new double[mrows];

        int cCol = 0;
        
        /*
        a0 1 2     0
                   1
                   2
        */
        
        for (int row = 0; row < mrows; row++) {
            for (int col = 0; col < mcols; col++) {
                c[cCol] += (a.get(row, col) * b[col]);
            }
            cCol++;        
        }

        return c;
    }
    
    public static void multiply(Matrix a, double b) {
        
        if (a == null || a.numRows() == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        Iterator<MatrixEntry> iter = a.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            entry.set(entry.get() * b);
        }
        
    }

    public static double trace(double[][] d) {
        
        double sum = 0;
        
        int n = Math.min(d.length, d[0].length);
        for (int i = 0; i < n; ++i) {
            sum += d[i][i];
        }
        
        return sum;
    }
    
    public static class EigenValuesAndVectors {
        private final double[] eigenValues;
        private final double[][] eigenVectors;
        public EigenValuesAndVectors(int nComponents) {
            eigenValues = new double[nComponents];
            eigenVectors = new double[2][];
        }
        public void setEigenValueAndVector(int index, double value, double[] vector) {
            eigenValues[index] = value;
            eigenVectors[index] = Arrays.copyOf(vector, vector.length);
        }
        public double[] getEigenVector(int index) {
            return eigenVectors[index];
        }
        public double getEigenValue(int index) {
            return eigenValues[index];
        }
        public int getNumberOfComponents() {
            return eigenValues.length;
        }
        public int getLengthOfVector() {
            return eigenVectors[0].length;
        }
    }
    
    public static double[][] copy(double[][] a) {
        
        double[][] m = new double[a.length][];
        
        for (int i = 0; i < a.length; ++i) {
            int n0 = a[i].length;
            m[i] = new double[n0];
            System.arraycopy(a[i], 0, m[i], 0, n0);
        }
        
        return m;
    }
    
    /**
     * calculate the dot product of a.transpose with a.  the result is in the
     * format as DenseMatrix, row major double[][] is [row][col].
     * @param a
     * @return 
     */
    public static double[][] arrayTransposeDotArray(double[] a) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length == 0) {
            throw new IllegalArgumentException(
                "a length must be > 0");
        }
        
        int n = a.length;
        
        /*
        a0    a0 a1 a2
        a1
        a2
        */
        
        double[][] m = new double[n][n];
        for (int i = 0; i < n; i++) {
            m[i] = new double[n];
        }
        
        for (int row = 0; row < n; ++row) {
            double v1 = a[row];
            for (int col = 0; col < n; col++) {
                double v2 = a[col];
                m[row][col] = (v1 * v2);
            }
        }

        return m;
    }
    
    /**
     * given a matrix with rows being features and columns being data, scale
     * the data to a mean of zero and stdev of 1 for each feature.
     * @param a
     * @return 
     */
    public static DenseMatrix scaleToUnitStandardDeviation(DenseMatrix a) {
        
        int n = a.numColumns();
        int nRows = a.numRows();
                
        double[] mean = new double[nRows];
        for (int i = 0; i < n; ++i) {            
            for (int j = 0; j < nRows; ++j) {
                mean[j] += a.get(j, i);
            }
        }
        for (int j = 0; j < nRows; ++j) {
            mean[j] /= (double)n;
        }
        
        double[] diff = new double[nRows];
        for (int i = 0; i < n; ++i) {            
            for (int j = 0; j < nRows; ++j) {
                double d = a.get(j, i) - mean[j];
                diff[j] += (d * d);
            }
        }
        for (int j = 0; j < nRows; ++j) {
            diff[j] /= (double)n;//((double)n - 1.);
        }
        
        double[] scaleFactors = new double[nRows];
        for (int j = 0; j < nRows; ++j) {
            scaleFactors[j] = Math.sqrt(1./diff[j]);
        }
        
        /*
        transformation for no rotation:  x*s - xc*s
           (data - mean) * scaleFactor
        */        
        DenseMatrix a2 = new DenseMatrix(nRows, a.numColumns());
        for (int i = 0; i < n; ++i) {            
            for (int j = 0; j < nRows; ++j) {
                double centered = a.get(j, i) - mean[j];
                double v = centered * scaleFactors[j];
                a2.set(j, i, v);
            }
        }
        
        return a2;
    }
    
    /**
     * given a matrix with rows being features and columns being data, scale
     * the data to a mean of zero and stdev of 1 for each feature.
     * @param a
     * @return 
     */
    public static DenseMatrix scaleToUnitStandardDeviation2(DenseMatrix a) {
        
        int n = a.numColumns();
        int nRows = a.numRows();
                
        double[] mean = new double[nRows];
        for (int i = 0; i < n; ++i) {            
            for (int j = 0; j < nRows; ++j) {
                mean[j] += a.get(j, i);
            }
        }
        for (int j = 0; j < nRows; ++j) {
            mean[j] /= (double)n;
        }
        
        double[] diff = new double[nRows];
        for (int i = 0; i < n; ++i) {            
            for (int j = 0; j < nRows; ++j) {
                double d = a.get(j, i) - mean[j];
                diff[j] += Math.abs(d);
            }
        }
        for (int j = 0; j < nRows; ++j) {
            diff[j] /= (n - 1.);
        }
        
        double sqrtTwo = Math.sqrt(2.);
        double[] scaleFactors = new double[nRows];
        for (int j = 0; j < nRows; ++j) {
            scaleFactors[j] = sqrtTwo/diff[j];
        }
        
        /*
        transformation for no rotation:  x*s - xc*s
           (data - mean) * scaleFactor
        */        
        DenseMatrix a2 = new DenseMatrix(nRows, a.numColumns());
        for (int i = 0; i < n; ++i) {            
            for (int j = 0; j < nRows; ++j) {
                double centered = a.get(j, i) - mean[j];
                double v = centered * scaleFactors[j];
                a2.set(j, i, v);
            }
        }
        
        return a2;
    }
    
    public static double[] extractRawPitchRollFromRotation(DenseMatrix rotMatrix) {
        
        double yaw = Math.atan2(rotMatrix.get(1, 0), rotMatrix.get(0, 0));
        
        double pitch = Math.atan2(-rotMatrix.get(2, 0),
            Math.sqrt(rotMatrix.get(2, 1)*rotMatrix.get(2, 1) +
                rotMatrix.get(2, 2)*rotMatrix.get(2, 2)));
         
        double roll = Math.atan2(rotMatrix.get(2, 1), rotMatrix.get(2, 2));
        
        return new double[]{yaw, pitch, roll};
    }
    
    public static DenseMatrix calculateRotationMatrix(double yaw,
        double pitch, double roll) {
                
        DenseMatrix rot = new DenseMatrix(3, 3);
        
        return calculateRotationMatrix(rot, yaw, pitch, roll);
    }
    
    protected static DenseMatrix calculateRotationMatrix(
        DenseMatrix m, double yaw, double pitch, double roll) {
                        
        m.set(0, 0, Math.cos(yaw) * Math.cos(pitch));
        m.set(0, 1, 
            (Math.cos(yaw)*Math.sin(pitch)*Math.sin(roll) - 
            Math.sin(yaw)*Math.cos(roll)));
        m.set(0, 2, 
            (Math.cos(yaw)*Math.sin(pitch)*Math.cos(roll) + 
            Math.sin(yaw)*Math.sin(roll)));
        
        m.set(1, 0, Math.sin(yaw) * Math.cos(pitch));
        m.set(1, 1, 
            (Math.sin(yaw)*Math.sin(pitch)*Math.sin(roll) + 
            Math.cos(yaw)*Math.cos(roll)));
        m.set(1, 2, 
            (Math.sin(yaw)*Math.sin(pitch)*Math.cos(roll) - 
            Math.cos(yaw)*Math.sin(roll)));
        
        m.set(2, 0, -Math.sin(pitch));
        m.set(0, 0, Math.cos(pitch) * Math.sin(roll));
        m.set(0, 0, Math.cos(pitch) * Math.cos(roll));
        
        return m;
    }
    
    public static String printMatrix(double[][] params) {
        
        if (params == null) {
            return "";
        }
        
        StringBuffer sb = new StringBuffer();
        for (int row = 0; row < params.length; row++) {
            sb.append("row ").append(Integer.toString(row)).append(":");
            for (int col = 0; col < params[0].length; col++) {
                sb.append(" ").append(Double.toString(params[row][col]));
            }
            sb.append("\n");
        }
        
        return sb.toString();
    }

   /**
     * given a one dimensional matrix of integer class numbers, modify the
     * matrix to start class number at zero and use intervals of value 1 and
     * return the number of classes.
     * For example, a matrix containing [1, 2, 5] would be transformed to
     * [0, 1, 2] and return the value 3.
     * @param classes
     * @return 
     */
    public static int transformToZeroBasedClasses(DenseMatrix classes) {
        
        Set<Integer> set = new HashSet<Integer>();        
        
        for (int i = 0; i < classes.numColumns(); ++i) {
            double v = classes.get(0, i);
            set.add(Integer.valueOf((int)Math.round(v)));
        }
                
        List<Integer> sorted = new ArrayList<Integer>(set);
        Collections.sort(sorted);
        
        // make a map for O(1) look ups
        Map<Integer, Integer> transformIndexMap = new HashMap<Integer, Integer>();
        for (int i = 0; i < sorted.size(); ++i) {
            transformIndexMap.put(sorted.get(i), Integer.valueOf(i));
        }
        
        for (int i = 0; i < classes.numColumns(); ++i) {
            double v = classes.get(0, i);
            int key = (int)Math.round(v);
            int v2 = transformIndexMap.get(Integer.valueOf(key));
            classes.set(0, i, v2);
        }
        
        return set.size();
    }
    
    /**
     * extract the null space from the SVD matrix.
     * The null space the non-zero vector x in which 
     * all solutions to the original matrix A * x = 0.
     * 
     * @param svd
     * @return 
     */
    public static DenseMatrix nullSpace(SVD svd) {
        
        /*
        SVD is the result of A = U * SIGMA * VT
        
            A is a real m-by-n matrix.
            SIGMA is an m-by-n matrix which is zero except for its min(m,n) 
                 diagonal elements.
            U is an m-by-m orthogonal matrix.
            VT (V transposed) is an n-by-n orthogonal matrix.
            
            The diagonal elements of SIGMA are the singular values of A; 
                they are real and non-negative, and are returned in 
                descending order. 
            The first min(m,n) columns of U and V are the left and right 
                singular vectors of A.

            The routine returns VT, not V.
        
        The (right) null space of A is the columns of V corresponding to
        singular values equal to zero. 
        The left null space of A is the rows of U corresponding to
        singular values equal to zero 
        (or the columns of U corresponding to singular values equal to
        zero, transposed).
        */
        
        double tol = singularThreshold(svd);
        
        //NXN
        DenseMatrix vT = svd.getVt();
        
        int N = vT.numColumns();
        
        // null space nRows = M (which is size of U)
        // null space nCols N - the number of items in s above tol.
        
        // s sometimes has fewer than N items if the last are all zeroes
        double[] s = svd.getS();
                
        // indexes in s that are zero or < tol
        TIntList zeroIdxs = new TIntArrayList();
        
        for (int i = 0; i < N; ++i) {
            if (i >= (s.length)) {
                zeroIdxs.add(i);
            } else {
                if (s[i] < tol) {
                    zeroIdxs.add(i);
                }
            }
        }
        
        DenseMatrix nullSpace = new DenseMatrix(N, zeroIdxs.size());
        
        /*
        v0=Type = dense real , numRows = 5 , numCols = 5
         -0.194   0.298  -0.473   0.087  -0.801
          0.388  -0.597  -0.323   0.621  -0.058
         -0.452   0.084   0.585   0.655  -0.134
         -0.710  -0.131  -0.525   0.086   0.443
          0.322   0.728  -0.233   0.413   0.376
                           /|\
                            |
                           z0      z1      z2
        */
        
        int r = 0;
        for (int i0 = 0; i0 < zeroIdxs.size(); ++i0) {
            int vCol = zeroIdxs.get(i0);
            // extract column vCol from vT and store it as column r in nullSpace
            for (int vRow = 0; vRow < vT.numRows(); ++vRow) {
                double v = vT.get(vRow, vCol);
                nullSpace.set(vRow, r, v);
            }
            ++r;
        }
        
        return nullSpace;
    }
    
    /**
     * Returns a reasonable threshold for singular values.<br><br>
     *
     * tol = max (size (A)) * largest sigma * eps;
     * 
     * The implementation is adapted from the EJML project
     * https://github.com/lessthanoptimal/ejml/blob/69baa142637e2adf45d90722cf785fabe3d74fe0/main/dense64/src/org/ejml/ops/SingularOps.java
     * 
     * which has copyright:
     * Copyright (c) 2009-2014, Peter Abeles. All Rights Reserved.
     *
     * This file is part of Efficient Java Matrix Library (EJML).
     *
     * Licensed under the Apache License, Version 2.0 (the "License");
     * you may not use this file except in compliance with the License.
     * You may obtain a copy of the License at
     *
     *   http://www.apache.org/licenses/LICENSE-2.0
     *
     * Unless required by applicable law or agreed to in writing, software
     * distributed under the License is distributed on an "AS IS" BASIS,
     * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
     * See the License for the specific language governing permissions and
     * limitations under the License.
     * 
     * @param svd
     * @return 
     */
    public static double singularThreshold(SVD svd) {
        
        double largest = 0;
        
        double w[] = svd.getS();

        int N = w.length;

        for( int j = 0; j < N; j++ ) {
            if( w[j] > largest)
                largest = w[j];
        }

        int M = Math.max(svd.getU().numColumns(), svd.getU().numRows());
        
        return M * largest * Math.pow(2,-52);
    }
    
    public static DenseMatrix extractAColumn(DenseMatrix m, int column) {
        
        DenseMatrix a = new DenseMatrix(m.numRows(), 1);
        for (int i = 0; i < m.numRows(); ++i) {
            a.set(i, 0, m.get(i, column));
        }
        
        return a;
    }
    
    public static DenseMatrix extractARow(DenseMatrix m, int row) {
        
        DenseMatrix a = new DenseMatrix(1, m.numColumns());
        for (int i = 0; i < m.numColumns(); ++i) {
            a.set(0, i, m.get(row, i));
        }
        
        return a;
    }
    
    /**
     * using cofactors and minors of the matrix, return the determinant.
     * in practice one can use any row as the primary set of cofactors or
     * any column.  this method may be optimized in the future, but for now,
     * uses the first column as the cofactors.
     *
     * e.g.    | 1  -5  2 |         | 3 4 |         | 7 4 |         | 7 3 |
     *         | 7   3  4 |  =  1 * | 1 5 |  +  5 * | 2 5 |  +  2 * | 2 1 |  = 11 
        + 135 + 2 = 148
     *         | 2   1  5 |
     */
    public static double determinant(Matrix m) {

        double[][] a = no.uib.cipr.matrix.Matrices.getArray(m);
        
        return determinant(a);
    }
    
    /**
     * using cofactors and minors of the matrix, return the determinant.
     * in practice one can use any row as the primary set of cofactors or
     * any column.  this method may be optimized in the future, but for now,
     * uses the first column as the cofactors.
     *
     * e.g.    | 1  -5  2 |         | 3 4 |         | 7 4 |         | 7 3 |
     *         | 7   3  4 |  =  1 * | 1 5 |  +  5 * | 2 5 |  +  2 * | 2 1 |  = 11 
        + 135 + 2 = 148
     *         | 2   1  5 |
     */
    public static double determinant(double[][] m) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        if (m.length != m[0].length) {
            throw new IllegalArgumentException("m must be a square");
        }
        if (m.length == 1) {
            return m[0][0];
        } else if (m.length == 2) {
            double s = ( m[0][0]*m[1][1] ) - ( m[0][1]*m[1][0] );
            return s;
        } else {
            double s = 0.0;
            // use 1st row as cofactors and minors
            for (int i = 0; i < m.length; i++) {

                double[][] n = copyExcept(m, i, 0);
                
                double tmp = m[i][0] * determinant(n);
                                
                if ((i & 1) == 0) {
                    s +=  tmp;
                } else {
                    s -=  tmp;
                }
            }
            return s;
        }
    }
    
    /**
     * create copy of matrix m except row and col
     * @param m
     * @param i
     * @param i0
     * @return
     */
    private static double[][] copyExcept(double[][] m, int col, int row) {

        double[][] n = new double[m.length - 1][m.length - 1];

        int nr = 0;
        int nc = 0;

        for (int mCol = 0; mCol < m.length; mCol++) {
            if (mCol == col) {
                continue;
            }

            n[nc] = new double[m.length - 1];
            
            nr = 0;
            for (int mRow = 0; mRow < m[0].length; mRow++) {
                if (mRow == row) {
                    continue;
                }

                n[nc][nr] = m[mCol][mRow];
                nr++;
            }
            nc++;
        }

        return n;
    }

    /**
     * find the equation for which A * A^(-1) = the identity matrix
     *
     *             1
     * A^(-1) =  ------ C^(T)  where C_ij = cofactor of a_ij
     *            det A
     *
     * @param m
     * @return
     */
    public static DenseMatrix inverse(DenseMatrix m) {

        ImageProcessor imageProcessor = new ImageProcessor();
            
        double[][] m2 = imageProcessor.copyToDouble2D(m);
           
        double[][] invM2 = inverse(m2);
        
        DenseMatrix invM = new DenseMatrix(invM2);
        
        return invM;
    }
    
    /**
     * find the equation for which A * A^(-1) = the identity matrix
     *
     *             1
     * A^(-1) =  ------ C^(T)  where C_ij = cofactor of a_ij
     *            det A
     *
     * @param m
     * @return
     */
    public static double[][] inverse(double[][] m) {

        // create cofactor of matrix:
        double[][] cofactor = createCofactor(m);

        double[][] cofactorTransposed = transpose(cofactor);

        double det = determinant(m);

        multiply(cofactorTransposed, 1./det);

        return cofactorTransposed;
    }
    
    public static void multiply(double[][] m, double factor) {

        int nrows = m.length;
        int ncols = m[0].length;

        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                m[i][j] = factor*m[i][j];
            }
        }
    }
    
    /*
         * e.g.    | 1  -5  2 |         | 3 4 |         | 7 4 |         | 7 3 |
         *         | 7   3  4 |  =  1 * | 1 5 |  +  5 * | 2 5 |  +  2 * | 2 1 |  = 11 + 135 + 2 =
 148
         *         | 2   1  5 |
         *
         *          3 4     7  4    7  3
         *          1 5     2  5    2  1
         * 
         *         -5 2     1  2    1 -5
         *          1 5     2  5    2  1
         *
         *         -5 2     1  2    1 -5
         *          3 4     7  4    7  3
    */
    public static double[][] createCofactor(double[][] m) {

        int ncols = m.length;
        int nrows = m[0].length;

        double[][] cofactor = new double[ncols][nrows];

        for (int i = 0; i < ncols; i++) {
            
            cofactor[i] = new double[nrows];

            boolean si = ((i & 1) == 1); // sign is -

            for (int j = 0; j < nrows; j++) {

                boolean sj = ((j & 1) == 1); // sign is -

                double[][] n = copyExcept(m, i, j);

                double cfctr = determinant(n);

                if (si ^ sj) { // XOR if either is 1 but not both
                    cfctr = -1*cfctr;
                }

                cofactor[i][j] = cfctr;
            }
         }
        return cofactor;
    }
    
    
}
