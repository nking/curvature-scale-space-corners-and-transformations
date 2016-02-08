package algorithms.imageProcessing.util;

import algorithms.MultiArrayMergeSort;
import algorithms.misc.Complex;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import org.ejml.simple.*;

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
    
    public static void add(int[] m, int n) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int len = m.length;
             
        for (int i = 0; i < len; i++) {
            m[i] += n;
        }
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
    
    public static double[][] dot(SimpleMatrix m1, SimpleMatrix m2) {
        
        if (m1 == null) {
            throw new IllegalArgumentException("m1 cannot be null");
        }
        if (m2 == null) {
            throw new IllegalArgumentException("m2 cannot be null");
        }
        if (m1.numCols() != m2.numRows()) {
            throw new IllegalArgumentException(
                "the number of columns in m1 != number of rows in m2");
        }
        
        int cCols = m2.numCols();
        int cRows = m1.numRows();
        
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
        
        for (int colAdd = 0; colAdd < m2.numCols(); colAdd++) {
            for (int row = 0; row < m1.numRows(); ++row) {
                for (int col = 0; col < m1.numCols(); col++) {
                    double a = m1.get(row, col);
                    double b = m2.get(col, colAdd);
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
    
    public static Complex[][] transpose(Complex[][] m) {

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
    
    /**
     * find best separation of classes given the data and return the
     * transformation.
     * adapted from http://sebastianraschka.com/Articles/2014_python_lda.html
     * 
     * To apply the results to a data matrix, use W.T.dot(X.T).T where W is
     * the returned matrix.  
     * Note, X internally has been scaled to unit standard deviation.
     * @param data
     * @param classes
     * @return 
     */
    public static SimpleMatrix createLDATrasformation(SimpleMatrix data,
        SimpleMatrix classes) {
        
        int n = data.numCols();
        
        if (classes.numCols() != n) {
            throw new IllegalArgumentException(
                "data and classes must have some number of data columns");
        }
                
        SimpleMatrix normData = scaleToUnitStandardDeviation(data);
       
        int nRows = normData.numRows();
        
        assert(nRows == data.numRows());
        
        // transforms from integer classes to zero based counting with delta of 1
        // for example:  [1, 2, 5, ...] becomes [0, 1, 2, ...]
        int nClasses = trasformToZeroBasedClasses(classes);
        
        SimpleMatrix w = createLDATrasformation2(normData, classes);
        
        return w;
    }
    
    /**
     * find best separation of classes given the already scaled data 
     * and the zero based classes and return the
     * transformation.
     * adapted from http://sebastianraschka.com/Articles/2014_python_lda.html
     * 
     * To apply the results to a data matrix, use W.T.dot(X.T).T where W is
     * the returned matrix.  
     * @param normData matrix already scaled to unit standard deviation (mean = 0
     * and standard deviation of the mean = 1).
     * @param classes the classes of the matrix data written such that the
     * classes are integers and the smallest value is 0 and the difference
     * between values is one (for example, the zero based classes of
     * [1, 2, 5] would be [0, 1, 2].
     * @return 
     */
    public static SimpleMatrix createLDATrasformation2(SimpleMatrix normData,
        SimpleMatrix classes) {
        
        int n = classes.numCols();
        
        if (classes.numCols() != n) {
            throw new IllegalArgumentException(
                "data and classes must have some number of data columns");
        }
        
        //Logger log = Logger.getLogger(MatrixUtil.class.getName());
               
        int nRows = normData.numRows();
        
        int nClasses = classes.numCols();
        
        double[][] mean = new double[nClasses][nRows];
        int[][] count = new int[nClasses][nRows];
        for (int k = 0; k < nClasses; ++k) {
            mean[k] = new double[nRows];
            count[k] = new int[nRows];
        }
        for (int i = 0; i < n; ++i) {
            int k = (int)Math.round(classes.get(0, i));
            for (int j = 0; j < nRows; ++j) {
                mean[k][j] += normData.get(j, i);
                count[k][j]++;
            }
        }
        
        for (int k = 0; k < nClasses; ++k) {
            for (int j = 0; j < nRows; ++j) {
                mean[k][j] /= (double)count[k][j];
            }
            //log.info(String.format("mean vector class %d = %s", k, 
            //    Arrays.toString(mean[k])));
        }
        
        // ----- calculate the scatter within each class --------
        double[][] scatWithinClass = new double[nRows][nRows];
        for (int j = 0; j < nRows; ++j) {
            scatWithinClass[j] = new double[nRows];
        }
            
        for (int k = 0; k < nClasses; ++k) {
            
            double[][] scatWithinClassPerClass = new double[nRows][nRows];
            for (int j = 0; j < nRows; ++j) {
                scatWithinClassPerClass[j] = new double[nRows];
            }
             
            for (int i = 0; i < n; ++i) {
                if (((int)Math.round(classes.get(0, i))) != k) {
                    continue;
                }
                double[] diff = new double[nRows];
                for (int j = 0; j < nRows; ++j) {
                    diff[j] = normData.get(j, i) - mean[k][j];
                }
                // calc (row-mv).dot((row-mv).T)
                double[][] dTdotD = arrayTransposeDotArray(diff);
                
                // add to scatWithinClassPerClass
                for (int j = 0; j < nRows; ++j) {
                    for (int jj = 0; jj < nRows; ++jj) {
                        scatWithinClassPerClass[j][jj] += dTdotD[j][jj];
                    }
                    //log.info(Arrays.toString(scatWithinClassPerClass[j]));
                }
            }
            
            // add scatWithinClassPerClass to scatWithinClass
            for (int j = 0; j < nRows; ++j) {
                for (int jj = 0; jj < nRows; ++jj) {
                    scatWithinClass[j][jj] += scatWithinClassPerClass[j][jj];
                }
                //log.info(Arrays.toString(scatWithinClass[j]));
            }
        }
        
        double[] overallMean = new double[nRows];
        for (int k = 0; k < nClasses; ++k) {
            for (int j = 0; j < nRows; ++j) {
                overallMean[j] += mean[k][j];
            }
        }
        for (int j = 0; j < nRows; ++j) {
            overallMean[j] /= (double)nClasses;
        }
        
        // ------- calculate the between class scatter matrix -------
        double[][] scatBetweenClass = new double[nRows][nRows];
        for (int j = 0; j < nRows; ++j) {
            scatBetweenClass[j] = new double[nRows];
        }
        
        // mean is length nClasses by nRows
        // overall_mean is length nRows
        
        double[] diff = new double[nRows];
        
        for (int k = 0; k < nClasses; ++k) {
            
            int countN = 0;
            Arrays.fill(diff, 0);
            for (int j = 0; j < nRows; ++j) {                                
                diff[j] = mean[k][j] - overallMean[j];
                countN += count[k][j];
            }
                       
            // calc (mean_vec - overall_mean).dot((mean_vec - overall_mean).T)
            double[][] dTdotD = arrayTransposeDotArray(diff);
            
            //S_B += n * (mean_vec - overall_mean).dot((mean_vec - overall_mean).T)
            // add to scatBetweenClass
            for (int j = 0; j < nRows; ++j) {
                for (int jj = 0; jj < nRows; ++jj) {
                    scatBetweenClass[j][jj] += (countN * dTdotD[j][jj]);
                }
                //log.info(Arrays.toString(scatBetweenClass[j]));
            }
        }
        
        //for (int j = 0; j < nRows; ++j) {
        //    log.info(Arrays.toString(scatBetweenClass[j]));
        //}
        
        // --- solve for the generalized eigenvalue, S_W^-1 * S_B ----
        SimpleMatrix sw = new SimpleMatrix(scatWithinClass);
        SimpleMatrix sb = new SimpleMatrix(scatBetweenClass);
        
        double[][] invSWDotSB = dot(sw.invert(), sb);
        SimpleMatrix m = new SimpleMatrix(invSWDotSB);
        SimpleEVD evd = m.eig();
        
        int nEigen = evd.getNumberOfEigenvalues();
        
        int[] indexes = new int[nEigen];
        double[] eigenValues = new double[nEigen];
        double[][] eigenVectors = new double[nEigen][nRows];
        for (int i = 0; i < nEigen; ++i) {
            eigenValues[i] = evd.getEigenvalue(i).getMagnitude();
            SimpleMatrix ev = evd.getEigenVector(i);
            eigenVectors[i] = new double[ev.numRows()];
            for (int j = 0; j < ev.numRows(); ++j) {
                eigenVectors[i][j] = ev.get(j, 0);
            }
            indexes[i] = i;
            //log.info("eigenvalue " + i + " = " + evd.getEigenvalue(i));
            //log.info("eigenvector=" + evd.getEigenVector(i));
        }
        
        // ----- sort by decreasing eigen value -----
        MultiArrayMergeSort.sortByDecr(eigenValues, indexes);
        
        SimpleMatrix w = new SimpleMatrix(2, nRows);
        for (int i = 0; i < 2; ++i) {
            int index = indexes[i];
            double[] eV = eigenVectors[index];
            for (int col = 0; col < eV.length; ++col) {
                w.set(i, col, eV[col]);
            }
        }
        
        //log.info("w=" + w);
        
        return w;
    }
    
    /**
     * calculate the dot product of a.transpose with a.  the result is in the
     * format that SimpleMatrix uses, double[][] is [row][col].
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
    public static SimpleMatrix scaleToUnitStandardDeviation(SimpleMatrix a) {
        
        int n = a.numCols();
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
        SimpleMatrix a2 = new SimpleMatrix(nRows, a.numCols());
        for (int i = 0; i < n; ++i) {            
            for (int j = 0; j < nRows; ++j) {
                double centered = a.get(j, i) - mean[j];
                double v = centered * scaleFactors[j];
                a2.set(j, i, v);
            }
        }
        
        return a2;
    }
    
    
    public static double[] extractRawPitchRollFromRotation(SimpleMatrix rotMatrix) {
        
        double yaw = Math.atan(rotMatrix.get(1, 0)/ rotMatrix.get(0, 0));
        
        double pitch = Math.atan(-rotMatrix.get(2, 0)/
            Math.sqrt(rotMatrix.get(2, 1)*rotMatrix.get(2, 1) +
                rotMatrix.get(2, 2)*rotMatrix.get(2, 2)));
         
        double roll = Math.atan(rotMatrix.get(2, 1)/ rotMatrix.get(2, 2));
        
        return new double[]{yaw, pitch, roll};
    }
    
    public static SimpleMatrix calculateRotationMatrix(double yaw,
        double pitch, double roll) {
                
        SimpleMatrix rot = new SimpleMatrix(3, 3);
        
        return calculateRotationMatrix(rot, yaw, pitch, roll);
    }
    
    protected static SimpleMatrix calculateRotationMatrix(
        SimpleMatrix m, double yaw, double pitch, double roll) {
                        
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

    private static int trasformToZeroBasedClasses(SimpleMatrix classes) {
        
        Set<Integer> set = new HashSet<Integer>();        
        
        for (int i = 0; i < classes.numCols(); ++i) {
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
        
        for (int i = 0; i < classes.numCols(); ++i) {
            double v = classes.get(0, i);
            int key = (int)Math.round(v);
            int v2 = transformIndexMap.get(Integer.valueOf(key));
            classes.set(0, i, v2);
        }
        
        return set.size();
    }
}
