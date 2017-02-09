package algorithms.compGeometry;

import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.util.PairInt;
import algorithms.util.PolynomialFitter;
import java.util.Arrays;
import java.util.Set;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.QR;
import no.uib.cipr.matrix.UpperTriangDenseMatrix;

/**
 
 y = ax^2 + bx + c

 least squares fit is minimizing the sum of square of residuals

 summation from i=1 to n points of
     (y_i - (ax_i^2 + bx_i + c))^2

 to minimize, need each partial deriv = 0,
 that is d(eqn)/da = 0, etc

 the 3 partial derivatives are summarized as

     a * sum_i=1toN (x_i)^4 + b * sum_i=1toN (x_i)^3
         + c * sum_i=1toN (x_i)^2
         = sum_i=1toN (x_i)^2 * y_i

     a * sum_i=1toN (x_i)^3 + b * sum_i=1toN (x_i)^2
         + c * sum_i=1toN (x_i)
         = sum_i=1toN (x_i) * y_i

     a * sum_i=1toN (x_i)^2 + b * sum_i=1toN (x_i)^1
         + c * n
         = sum_i=1toN y_i

     given a set of points, can calc moments
       x    y    x^2    x^3    x^4    xy    yx^2
     (when placed in a matrix, is known as Vandermonde)

       sum each of those for the points

       then can solve the 3 equations for a, b, c using linear algebra.
       
       

 * @author nichole
 */
public class ParabolaLeastSquares {
 
    //x    y    x^2    x^3    x^4    xy    yx^2
    protected final double[] moments_;
    protected int n;
    
    private boolean debug = false;
    
    public ParabolaLeastSquares() {
        moments_ = new double[7];
        n = 0;
    }
    
    public void setToDebug() {
        debug = true;
    }
    
    public void accumulate(double x, double y) {
        ++n;
        //x    y    x^2    x^3    x^4    xy    yx^2
        moments_[0] += x;
        moments_[1] += y;
        double tmp = x * x;
        moments_[2] += tmp;
        tmp *= x;
        moments_[3] += tmp;
        tmp *= x;
        moments_[4] += tmp;
        moments_[5] += x * y;
        moments_[6] += y * x * x;
    }
    
    double[][] getMatrixA() {
        
        double[][] A = new double[3][3];
        A[0] = new double[]{moments_[4], moments_[3], moments_[2] };
        A[1] = new double[]{moments_[3], moments_[2], moments_[0] };
        A[2] = new double[]{moments_[2], moments_[0], n};   
        
        return A;
    }
    
    double[][] getRHS() {
        double[][] Y = new double[3][1];
        Y[0] = new double[]{moments_[6]};
        Y[1] = new double[]{moments_[5]};
        Y[2] = new double[]{moments_[1]};
        return Y;
    }
    
    double[] getRHSVector() {
        double[] Y = new double[]{moments_[6], moments_[5], moments_[1]};
        return Y;
    }
    
    public float[] solve() {
    
        /*
        a * sum_i=1toN (x_i)^4 + b * sum_i=1toN (x_i)^3
             + c * sum_i=1toN (x_i)^2
             = sum_i=1toN (x_i)^2 * y_i

         a * sum_i=1toN (x_i)^3 + b * sum_i=1toN (x_i)^2
             + c * sum_i=1toN (x_i)
             = sum_i=1toN (x_i) * y_i

         a * sum_i=1toN (x_i)^2 + b * sum_i=1toN (x_i)^1
             + c * n
             = sum_i=1toN y_i
        
        //x    y    x^2    x^3    x^4    xy    yx^2
        
        a * moments_[4] + b * moments_[3] + c * moments_[2] 
            = moments_[6]

        a * moments_[3] + b * moments_[2] + c * moments_[0]
             = moments_[5]
        
        a * moments_[2] + b * moments_[0] + c*area 
             = moments_[1]
        
        |a|
        |b| = inv of coeff above * right sides 
        |c|
        
        */
        
        // helpful also was the Go solution of http://rosettacode.org/wiki/Polynomial_Fitting
        
        double[][] A = getMatrixA();
        
        double[][] Y = getRHS();
        
        Matrix mA = new DenseMatrix(A);
        
        if (debug) {
            System.out.println("row0=" + Arrays.toString(A[0]));
            System.out.println("row1=" + Arrays.toString(A[1]));
            System.out.println("row2=" + Arrays.toString(A[2]));
        }
        
        QR qr = QR.factorize(mA);
        DenseMatrix q = qr.getQ();
        UpperTriangDenseMatrix r = qr.getR();
        
        Matrix qT = q.transpose();
        double[] qTY = MatrixUtil.multiply(qT, getRHSVector());
        float[] c = new float[3];
        for (int i = (3 - 1); i >= 0; i--) {
            c[i] = (float)qTY[i];
            for (int j = (i + 1); j < 3; j++) {
                c[i] -= c[j] * r.get(i, j);
            }
            c[i] /= r.get(i, i);
        }
        
        if (debug) {
            
            System.out.println("abc=" + Arrays.toString(c));
        
            for (int i = 0; i < 3; ++i) {
                double s = 0;
                for (int j = 0; j < 3; ++j) {
                    s += c[j] * A[i][j];
                }
                System.out.println("sum=" + s + " rhs=" + Y[i][0]);
            }
        }
        //qr_res  = test of QR factorization, norm1(Q*R-A)/(n*eps)
        
        return c;
    }
    
    public static String plotFit(float[] coefficients, Set<PairInt> points, 
        int plotXMax, int plotYMax, int plotNumber, String plotLabel) {
    
        return PolynomialFitter.plotFit(coefficients, points, plotXMax, plotYMax, 
            plotNumber, plotLabel);
    }
    
    public static double calcResiduals(float[] coefficients, Set<PairInt> points) {
        
        return PolynomialFitter.calcResiduals(coefficients, points);
    }
}
