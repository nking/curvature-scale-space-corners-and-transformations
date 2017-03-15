package algorithms.util;

import algorithms.MultiArrayMergeSort;
import algorithms.misc.MiscMath;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.QRDecomposition;
import org.ejml.simple.SimpleMatrix;

/**
 * methods associated with fitting a 2nd order polynomial curve.
 * 
 * @author nichole
 */
public class PolynomialFitter {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * solve for 2nd order curve for a random sample of 1000 points from
     * (dataX, dataY).
     * 
     * @param points
     * @return 2nd order polynomial coefficients if solved, else null
     */
    public float[] solveAfterRandomSampling(Set<PairInt> points) {
        
        SecureRandom sr = new SecureRandom();
        
        sr.setSeed(System.currentTimeMillis());
        
        return solveAfterRandomSampling(points, sr);
    }
    
    /**
     * solve for 2nd order curve for a random sample of 1000 points from
     * (dataX, dataY).
     * 
     * @param points
     * @param sr instance of secure random to use for generating random numbers
     * @return 2nd order polynomial coefficients if solved, else null.
     * the coefficients are used in y = c0*1 + c1*x[i] + c2*x[i]*x[i]
     */
    protected float[] solveAfterRandomSampling(Set<PairInt> points,
        SecureRandom sr) {
        
        int n = (points.size() > 2500) ? 2500 : points.size();
        
        List<PairInt> tmp = new ArrayList<PairInt>(points);
        
        int[] indexes = new int[n];
        if (n != points.size()) {
            MiscMath.chooseRandomly(sr, indexes, points.size());
        } else {
            for (int i = 0; i < n; i++) {
                indexes[i] = i;
            }
        }
        
        float[] xP = new float[n];
        float[] yP = new float[xP.length];
        int i = 0;
        for (int idx : indexes) {
            PairInt p = tmp.get(idx);
            xP[i] = p.getX();
            yP[i] = p.getY();
            i++;
        }
        tmp = null;
        
        // y = c0*1 + c1*x[i] + c2*x[i]*x[i]
        float[] coeff = solve(xP, yP);
        
        return coeff;
    }
    
    /**
     * solve for 2nd order curves.   This uses a Vandermonde matrix and QR 
     * decomposition.  Note that the solver is not ideal for a mostly
     * vertical polynomial with an extremely large gap separating two
     * groups of points.
     * 
     * @param dataX
     * @param dataY
     * @return 2nd order polynomial coefficients if solved, else null.
     * The coefficients are used in y = c0*1 + c1*x[i] + c2*x[i]*x[i]
     */
    public float[] solve(float[] dataX, float[] dataY) {
        
        // adapted from the Go solution of http://rosettacode.org/wiki/Polynomial_Fitting
        
        /*
        polyDegree unknowns and m data
        
        solving for the coefficients of y = c0*1 + c1*x[i] + c2*x[i]*x[i]
        
        The Vandermonde matrix fills the matrix with the representation of the
        x powers of data, that is the "monomials" 1, x[i], x[i]^2, ... x[i]^(n-1)
        
        | 1  x[i]    (x[i])^2   |  | c0 |    | y[i]   |
        | 1  x[i+1]  (x[i+1])^2 |  | c1 | =  | y[i+1] |
        | ..  ...      ...      |  | c2 |    | ...    |
        
        X * c = Y
        
        then c = Y / X is ==> c = X^-1 * Y
        
        for over determined systems, QR decomposition (or SVD or normal equations) 
        can be used.
            find Q, R such that X = Q^T * R (QR decomposition)
            then solve R*c = Q*Y
        
        */
        int m = dataX.length;
        
        int polyDegree = 2;
        
        int n = polyDegree + 1;
        
        SimpleMatrix y = new SimpleMatrix(m, 1);
        
        SimpleMatrix x = new SimpleMatrix(m, n);
        
        for (int i = 0; i < m; i++) {
            y.set(i, 0, dataY[i]);
            float ip = 1.f;
            for (int j = 0; j < n; j++) {
                x.set(i, j, ip);
                ip *= dataX[i];
            }
        }
        
        QRDecomposition<DenseMatrix64F> qrDecomp = DecompositionFactory.qr(m, n);
        if (!qrDecomp.decompose(x.getMatrix())) {
            return null; 
        }
       
        //nCols=X.length, nRows=X.length
        DenseMatrix64F q = qrDecomp.getQ(null, false);
        
        //nCols=3, nRows=X.length
        DenseMatrix64F r = qrDecomp.getR(null, false);
        
        //nCols=X.length, nRows=X.length
        SimpleMatrix qq = new SimpleMatrix(q);
        
        //nCols=1, nRows=X.length
        SimpleMatrix qty = qq.transpose().mult(y);
        
        float[] c = new float[n];

        for (int i = (n - 1); i >= 0; i--) {
            c[i] = (float)qty.get(i, 0);
            for (int j = (i + 1); j < n; j++) {
                c[i] -= c[j] * r.get(i, j);
            }
            c[i] /= r.get(i, i);
        }
        
        return c;
    }
   
    /**
     * plot the points and 2nd order curve given the coefficients.
     * note that if there are more than 1000 points, only 1000 of
     * them are 
     * @param coefficients
     * @param points
     * @param plotXMax
     * @param plotYMax
     * @param plotNumber
     * @param plotLabel
     * @return 
     */
    public static String plotFit(float[] coefficients, Set<PairInt> points, 
        int plotXMax, int plotYMax, int plotNumber, String plotLabel) {
                
        // shape the rainbow points into a more even ribbon
        
        float[] xP = new float[points.size()];
        float[] yP = new float[xP.length];
        int i = 0;
        for (PairInt p : points) {
            float x = p.getX();
            xP[i] = x;
            yP[i] = p.getY();
            i++;
        }
        
        float[] minXYMaxXY = determineGoodEndPoints(coefficients, points);
        float xMin = minXYMaxXY[0];
        float xMax = minXYMaxXY[2];
        
        int nCurve = 100;
        float dx = (xMax - xMin)/(float)nCurve;
        float[] xPoly = new float[nCurve];
        float[] yPoly = new float[xPoly.length];
        for (i = 0; i < nCurve; i++) {
            xPoly[i] = xMin + i*dx;
            yPoly[i] = coefficients[0] + coefficients[1]*xPoly[i] 
                + coefficients[2]*xPoly[i]*xPoly[i];
        }
        
        try {
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0, 
                plotXMax, 0, plotYMax);
            plotter.addPlot(xP, yP, xPoly, yPoly, plotLabel);
            
            String fileName = plotter.writeFile(Integer.valueOf(plotNumber));
            
            return fileName;
            
        } catch (IOException e) {
            Logger.getLogger(PolynomialFitter.class.getName()).severe(e.getMessage());
        }
        
        return null;
    }
    
    /**
     * determine good end points for the polynomial fit to the points.
     * runtime complexity is O(N*lg_2(N)) + O(N)
     * @param coefficients
     * @param points
     * @return float[]{minX, YForMinX, maxX, YFoMaxX}
     */
    public static float[] determineGoodEndPoints(float[] coefficients,
        Set<PairInt> points) {
        
        /*
        determine good endpoints for the polynomial solution.
        it depends upon the orientation of the rainbow polynomial.
        
        should be able to determine the average to the 100 or so median
        values as the middle of the rainbow,
        then would find the smallest residual points from the model
        polynomial that are located furthest from the median location.
        */
        
        // sort points by x then y
        int[] x = new int[points.size()];
        int[] y = new int[x.length];
        int i = 0;
        for (PairInt p : points) {
            x[i] = p.getX();
            y[i] = p.getY();
            i++;
        }
        //O(N*lg_2(N))
        MultiArrayMergeSort.sortBy1stArgThen2nd(x, y);
        
        // average of central 10 or so median
        int mid = points.size() >> 1;
        float medianX;
        float medianY;
        if (points.size() < 10) {
            medianX = x[points.size()/2];
            medianY = y[points.size()/2];
        } else {
            double sumX = 0;
            double sumY = 0;
            int nMed = 5;
            for (i = (mid - nMed); i <= (mid + nMed); i++) {
                sumX += x[i];
                sumY += y[i];
            }
            medianX = (float)sumX/(float)(2*nMed);
            medianY = (float)sumY/(float)(2*nMed);
        }
        
        // find the furthest points that have the smallest residuals on each side
        int minX = -1;
        int yForMinX = -1;
        int maxX = -1;
        int yForMaxX = -1;
        for (int half = 0; half < 2; half++) {
            double minResid = Double.MAX_VALUE;
            int minResidIdx = -1;
            double distFromMedianSq = Double.MIN_VALUE;
            if (half == 0) {
                for (i = 0; i < mid; ++i) {
                    float yV = coefficients[0] * (coefficients[1] * x[i]) +
                        (coefficients[2] * x[i] * x[i]);
                    double resid = Math.abs(y[i] - yV);

                    double distX = x[i] - medianX;
                    double distY = y[i] - medianY;

                    double distSq = (distX * distX) + (distY * distY);

                    if ((resid < minResid) && (distSq >= distFromMedianSq)) {
                        minResid = resid;
                        minResidIdx = i;
                        distFromMedianSq = distSq;
                    }
                }
                minX = x[minResidIdx];
                yForMinX = y[minResidIdx];
            } else {
                for (i = (points.size() - 1); i > mid; --i) {
                    float yV = coefficients[0] * (coefficients[1] * x[i]) +
                        (coefficients[2] * x[i] * x[i]);
                    double resid = Math.abs(y[i] - yV);

                    double distX = x[i] - medianX;
                    double distY = y[i] - medianY;

                    double distSq = (distX * distX) + (distY * distY);

                    if ((resid < minResid) && (distSq >= distFromMedianSq)) {
                        minResid = resid;
                        minResidIdx = i;
                        distFromMedianSq = distSq;
                    }
                }
                maxX = x[minResidIdx];
                yForMaxX = y[minResidIdx];
            }
        } 
        
        return new float[]{minX, yForMinX, maxX, yForMaxX};
    }
    
    /**
     * calculate the square root of the sum of the squared differences between 
     * a 2nd order polygon defined by the given coefficients and the given 
     * points.
     * Note that if coefficients or points are null or empty, it returns
     * a result of infinity.
     * 
     * @param coefficients
     * @param points
     * @return 
     */
    public static double calcResiduals(float[] coefficients, Set<PairInt> points) {
        
        if (points == null || points.isEmpty()) {
            return Double.POSITIVE_INFINITY;
        }
        
        if (coefficients == null || (coefficients.length != 3)) {
            return Double.POSITIVE_INFINITY;
        }
                
        //float[] minXYMaxXY = determineGoodEndPoints(coefficients,points);
                
        //double xMin = minXYMaxXY[0];
        
        // these can be large, so use abs value instead of sum of squares
        
        double sum = 0;
        
        for (PairInt p : points) {
            
            double x = p.getX();
            
            double y = p.getY();         
            
            double yPoly = coefficients[0] + (coefficients[1]*x) 
                + (coefficients[2]*x*x);
            
            double diff = y - yPoly;
            
            sum += Math.abs(diff);
        }
        
        sum /= (double)points.size();
        
        return sum;
    }

    /**
     * calculate the square root of the sum of the squared differences between 
     * a 2nd order polygon defined by the given coefficients and the given 
     * points.
     * Note that if coefficients or points are null or empty, it returns
     * a result of infinity.
     * 
     * @param coefficients
     * @param points
     * @return 
     */
    public double[] calcResidualsForAvg(float[] coefficients, Set<PairInt> points) {
        
        if (points == null || points.isEmpty()) {
            return new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        }
        
        if (coefficients == null || (coefficients.length != 3)) {
            return new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        }
                
        //float[] minXYMaxXY = determineGoodEndPoints(coefficients,points);
                
        //double xMin = minXYMaxXY[0];
        
        // these can be large, so use abs value instead of sum of squares
        
        double sum = 0;
        
        for (PairInt p : points) {
            
            double x = p.getX();
            
            double y = p.getY();            
            
            double yPoly = coefficients[0] + (coefficients[1]*x) 
                + (coefficients[2]*x*x);
            
            double diff = y - yPoly;
            
            sum += diff;
        }
        
        double avg = sum/(double)points.size();
        
        sum = 0;
        
        for (PairInt p : points) {
            
            double x = p.getX();
            
            double y = p.getY();
            
            double yPoly = coefficients[0] + (coefficients[1]*x) 
                + (coefficients[2]*x*x);
            
            double diff = (y - yPoly) - avg;
            
            sum += (diff * diff);
        }
        
        double stDev = Math.sqrt(sum/((double)points.size() - 1.0));
                
        return new double[]{avg, stDev};
    }
    
    private static class IntArrayWrapper {
        int[] a;
        public IntArrayWrapper(int[] values) {
            a = Arrays.copyOf(values, values.length);
        }
    }
}
