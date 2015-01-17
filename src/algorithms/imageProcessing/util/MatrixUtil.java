package algorithms.imageProcessing.util;

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
    
    public static void multiply(int[] m, int factor) {

        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int len = m.length;
                
        for (int i = 0; i < len; i++) {                        
            m[i] *= factor;
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
            c[row] = new double[mcols];
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
        0     1     0     2  1
        1000  100  10     3  0
                          4  0
        
        0*2    + 1*3   + 0*4     0*1    +  1*0   +  0*0
        1000*2 + 100*3 + 10*4    1000*1 +  100*0 + 10*0
        */
                
        for (int colAdd = 0; colAdd < m2.numCols(); colAdd++) {
            for (int cRow = 0; cRow < cRows; cRow++) {
                for (int col = 0; col < m1.numCols(); col++) {
                
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
}
