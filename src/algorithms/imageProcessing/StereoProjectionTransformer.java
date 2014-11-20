package algorithms.imageProcessing;

import algorithms.util.PairFloatArray;
import Jama.*;
import algorithms.imageProcessing.util.MatrixUtil;
import java.util.Arrays;
import java.util.logging.Logger;

/**
 * class to solve for the epipoles for two images with stereo projection
 * and apply the solution.
 * 
 * Following:
 * 
 * IEEE TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE, VOL. 19, 
 * NO. 6, JUNE 1997
 * "In Defense of the Eight-Point Algorithm" by Richard Hartley

 * The 8-point algorithm matrix represents epipolar geometry, and can be
 * used with data from cameras whose characteristics are not known to solve
 * up to the projective transformation.
 * 
 * Some definitions:
    u^T*v represents the inner product
    u*v^T is a matrix
    the norm of a vector f is the square root of the sum or squares of its 
        entries.
    u_1 is the (x,y) points from image 1 and u_2 are the matched (x,y) points 
        from image 2.
    
    the fundamental matrix is defined:
        u_2^T * F * u_1 = 0  where u are the x,y points in images _1 and _2

    u_1 = (x_1, y_1, 1)^T
    u_2 = (x_2, y_2, 1)^T

    x_1*x_2*F_1_1 + x_1*y_2*F_2_1 + x_1*F_3_1 + y_1*x_2*F_1_2 + y_1*y_2*F_2_2
        + y_1*F_3_2 + x_2*F_1_3 + y_2*F_2_3 + F_3_3 = 0

    A * f = 0

    where A = x_1*x_2, x_1*y_2, x_1, y_1*x_2, y_1*y_2, y_1, x_2, y_2, 1

    To avoid the trivial scale, ||f|| = 1 where f is the norm of f

    And we need least squares fits because the set may be over determined 
    and not have a zero solution.

    we want the vector f that minimizes ||A*f|| subject to the constraint
    that ||f|| = f^T*f = 1

    the solution is the unit eigenvector corresponding to the smallest 
    eigenvalue of A^T*A.

    Since A^T*A is semi-definite and symmetric, all of its eigenvectors
    are real and positive or zero.
    This eigenvector is what he calls the least eigenvector of A^T*A and
    it is found via the Jacobi algorithm or Singular Value Decomposition.

    The solved for matrix will in general not have rank 2 and needs to, so
    further corrections are necessary:
        matrix F is replaced by F' that minimizes the Frobenius norm
        ||F - F'|| subject to the condition det F' = 0.
        A convenient method of doing this is to use the Singular Value
        Decomposition (SVD).
           let F = U*D*V^T be the SVD of F, where D is diagonal matrix 
           D = diag(r, s, t) satisfying r >= s >= t.
           let F' = U*diag(r, s, 0)*V^T.

 (1) Transforming the coordinates:

     Normalization for isotropic scaling.

     utrans = T * u ==> u = utrans * inv(T)

     u_2^T * F * u_1 = 0 
     
        becomes   utrans_2^T * inv(T_2) * F * inv(T_1) * utrans_1 = 0

        and inv(T_2) * F * inv(T_1) is the fundamental matrix for
        utrans_2 < -- > utrans_1 which when found, will be subsequently
        denormalized.

    a) points are translated so that their centroid is at the origin.
    b) points are then scaled so that the average distance from the
       origin is sqrt(2)
    c) the transformation is applied to each of the 2 images separately.
     
 (2) build matrix A with the normalized x,y points
 
 (3) compute linear least square solution to the least eigenvector of f.
     solve A = U * D * V^T   for A*f = [..x...]*f = 0
     A has rank 8.  f has rank 2.
     calculate [U,D,V] from svd(A)
 
 (4) make the fundamental matrix have a rank of 2
     by performing a svd and then reconstructing with the two largest 
     singular values.
         [U,D,V] = svd(F,0);
         F = U * diag([D(1,1) D(2,2) 0]) * V'; 
 
 (5) denormalize the fundamental matrix
     The related part of the normalization equation: inv(T_2) * F * inv(T_1)
     so denormalizing is:
     
         F = (T_1)^T * F * T_2
   
 (6) estimate the error in the fundamental matrix by calculting epipolar
     lines for points in image 1 and find their nearest points in image 2 
     and measure the perpendicular distance from the epipolar line for
     those nearest points.
     
 * @author nichole
 */
public class StereoProjectionTransformer {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    private Matrix leftXY = null;
        
    private Matrix rightXY = null;
        
    private Matrix fundamentalMatrix = null;
    
    private double[] leftEpipole = null;
    
    private double[] rightEpipole = null;
    
    /**
     * each row is an epipolar line in the right image.
     * Each column corresponds to a point in leftXY and rightXY which are
     * in the same column.
     */
    private double[][] epipolarLinesInRight = null;
    
    /**
     * each row is an epipolar line in the left image.
     * Each column corresponds to a point in leftXY and rightXY which are
     * in the same column.
     */
    private double[][] epipolarLinesInLeft = null;
        
    public void calculateEpipolarProjection(PairFloatArray refactorLeftXY, 
        PairFloatArray refactorRightXY) {
        
        if (refactorLeftXY == null) {
            throw new IllegalArgumentException("refactorLeftXY cannot be null");
        }
        if (refactorRightXY == null) {
            throw new IllegalArgumentException("refactorRightXY cannot be null");
        }
        if (refactorLeftXY.getN() != refactorRightXY.getN()) {
            throw new IllegalArgumentException(
                "refactorLeftXY and refactorRightXY must be same size");
        }
        
        if (refactorLeftXY.getN() < 8) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the 8-point problem requires 8 or more points." 
                + " refactorLeftXY.n=" + refactorLeftXY.getN());
        }
        
        leftXY = rewriteInto3ColumnMatrix(refactorLeftXY);
        
        rightXY = rewriteInto3ColumnMatrix(refactorRightXY);
        
        fundamentalMatrix = calculateFundamentalMatrix(leftXY, rightXY);
     
        // 2 x 3 matrix of leftEpipole in column 0 and rightEpipole in column 1.
        double[][] leftRightEpipoles = calculateEpipoles(fundamentalMatrix);
        
        leftEpipole = leftRightEpipoles[0];
        
        rightEpipole = leftRightEpipoles[1];
        
        epipolarLinesInRight = calculateRightEpipolarLines();
        
        epipolarLinesInLeft = calculateLeftEpipolarLines();
       
        /*
        compute the perpendicular errors:
        
        transform points from the first image to get the equipolar lines
           in the second image.
        
        find the closest points in the 2nd image to the epipolar lines 
           store the difference
        
        do the same for the othe image
        
        equipolar lines:
           aVector = F^T * u_1
           then aVector^T*u_2 = aVector_1*x_2 + aVector_2*y_2 + aVector_3 = 0
        
        plot the differences, calc stats, determine inliers.
        use iterative method of choosing 8 or 8-best and error inspection to
        create a better solution (terminate when set of inliers does not 
        change).
        
        */
    }
    
    protected Matrix calculateFundamentalMatrix(Matrix leftXY, 
        Matrix rightXY) {
        
        //x is xy[0], y is xy[1], xy[2] is all 1's
        NormalizedXY normalizedXY1 = normalize(leftXY);
        
        NormalizedXY normalizedXY2 = normalize(rightXY);        
        
        //build the fundamental matrix
        Matrix aMatrix = new Matrix(createFundamentalMatrix(
            normalizedXY1.getXy(), normalizedXY2.getXy()));
               
        /*
        compute linear least square solution:
            solve A = U * D * V^T   for A*f = [..x...]*f = 0
            A has rank 8.  f has rank 2.
            calculate [U,D,V] from svd(A)
        */       
        SingularValueDecomposition svd = aMatrix.svd();
        
        // creates U as 9 x nXY1 matrix
        //         D as length 9 array
        //         V as 9 x 9 matrix
        
        // f is the last column of V.  it's got 9 items in it.
        double[] f = svd.getV().getArray()[svd.getV().getArray().length - 1];
        
        // reshape it to 3x3
        double[][] ff = new double[3][3];
        for (int i = 0; i < 3; i++) {
            ff[i] = new double[3];
        }
        for (int i = 0; i < 3; i++) {
            ff[i][0] = f[(i * 3) + 0];
            ff[i][1] = f[(i * 3) + 1];
            if (i == 2 && (f.length >= 9)) { // npoints == 8
                ff[i][2] = f[(i * 3) + 2];
            }
        }
        Matrix fMatrix = new Matrix(ff);
        
        /* make the fundamental matrix have a rank of 2
           by performing a svd and then reconstructing with the two largest 
           singular values.
              [U,D,V] = svd(F,0);
              F = U * diag([D(1,1) D(2,2) 0]) * V';
        */
        svd = fMatrix.svd();
        
        // creates U as 3 x 3 matrix
        //         D as length 3 array
        //         V as 3 x 3 matrix
        
        Matrix d = svd.getS();
        
        // remove the smallest singular value from D, making it rank 2
        double[] keep = new double[]{d.get(0, 0), d.get(1, 1), d.get(2, 2)};
        Arrays.sort(keep);
        d.set(0, 0, keep[2]);
        d.set(1, 1, keep[1]);
        d.set(2, 2, 0);
        
        /*
        multiply the terms:
             F = dot(U, dot(diag(D),V))
        */
        double[][] dDotV = MatrixUtil.dot(d, svd.getV());
                                
        // 3x3
        Matrix theFundamentalMatrix = new Matrix(MatrixUtil.dot(svd.getU(), 
            new Matrix(dDotV)));
                        
        /*
        denormalize
            F = (T_1)^T * F * T_2  
            where T_1 is normalizedXY1.getNormalizationMatrix();
            and T2 is normalizedXY2.getNormalizationMatrix();
        */
        
        // 3x3
        Matrix t1Transpose = normalizedXY1.getNormalizationMatrix().transpose();
        Matrix t2 = normalizedXY2.getNormalizationMatrix();
        
        Matrix denormFundamentalMatrix = t1Transpose.times(
            theFundamentalMatrix.times(t2));
        
        log.info("F=");
        for (int j = 0; j < denormFundamentalMatrix.getArray()[0].length; j++) {
            StringBuilder sb = new StringBuilder();            
            for (int i = 0; i < denormFundamentalMatrix.getArray().length; i++) {
                sb.append(denormFundamentalMatrix.get(i, j)).append(" ");                
            }            
            log.info(sb.toString());
        }
        
        return denormFundamentalMatrix;
    }   

    /**
     normalize the x,y coordinates as recommended by Hartley 1997 and return
     the matrix and coordinates.
       
     * @param xyPair
     * @return 
     */
    private NormalizedXY normalize(Matrix xy) {
        
        /*
        utrans = T * u ==> u = utrans * inv(T)
        
            utrans_2^T * inv(T_2) * F * inv(T_1) * utrans_1
        
        format the tensors T_1 and T_2 such that the applied translation
        and scaling have the effect of:
       
        a) points are translated so that their centroid is at the origin.
        b) points are then scaled so that the average distance from the
           origin is sqrt(2)
        c) the transformation is applied to each of the 2 images separately.
        */
       
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        double[] centroidXY = curveHelper.calculateXYCentroids(xy);
        
        double mean = 0;
        int n = xy.getArray()[0].length;
        for (int i = 0; i < n; i++) {
            double diffX = xy.get(0, i) - centroidXY[0];
            double diffY = xy.get(1, i) - centroidXY[1];
            double dist = Math.sqrt((diffX * diffX) + (diffY * diffY));
            mean += dist;
        }
        
        mean /= (double)n;
        
        /*
        mean * factor = sqrt(2)        
        */
        double scaleFactor = Math.sqrt(2)/mean;
        
        /*
        x_1_0   y_1_0  1          t00   0   0
        x_1_1   y_1_1  1            0  t11  0
        x_1_2   y_1_2  1          t02  t12  1
        x_1_3   y_1_3  1
        
        x_1_0*t00 + y_1_0*t01 + 1*t02 = (x_1_0 - cX) * s = x_1_0 * s - cX * s
                         0
             => t01 = 0
             => t00 = scaleFactor 
             => t02 = -scaleFactor * centroidXY[0]
        
        x_1_0*t10 + y_1_0*t11 + 1*t12 = (y_1_0 - cY) * s = y_1_0 * s - cY * s
            0
             => t10 = 0
             => t11 = scaleFactor 
             => t12 = -scaleFactor * centroidXY[1]
        */
        
        double[][] t = new double[3][];
        t[0] = new double[]{scaleFactor, 0,           -scaleFactor * centroidXY[0]};
        t[1] = new double[]{0,           scaleFactor, -scaleFactor * centroidXY[1]};
        t[2] = new double[]{0,           0,           1};
        Matrix tMatrix = new Matrix(t);
                
        Matrix normXY = new Matrix(MatrixUtil.dot(xy, tMatrix));
              
        NormalizedXY normalizedXY = new NormalizedXY();
        normalizedXY.setCentroidXY(centroidXY);
        normalizedXY.setNormMatrix(tMatrix);
        normalizedXY.setXy(normXY);
        
        return normalizedXY;
    }
    
    protected Matrix rewriteInto3ColumnMatrix(PairFloatArray xyPairs) {
        
        // rewrite xyPairs into a matrix of size 3 X xy.getN();
        // column 0 is x
        // column 1 is y
        // column 2 is all 1's
        double[][] xyPoints = new double[3][xyPairs.getN()];
        for (int i = 0; i < 3; i++) {
            xyPoints[i] = new double[xyPairs.getN()];
        }
        for (int i = 0; i < xyPairs.getN(); i++) {
            xyPoints[0][i] = xyPairs.getX(i);
            xyPoints[1][i] = xyPairs.getY(i);
            xyPoints[2][i] = 1;
        }
        
        Matrix xy = new Matrix(xyPoints);
        
        return xy;
    }
  
    /**
     * @param normXY1 a matrix of size 3 x nPoints, where 1st column is x,
     * second is y.
     * @param normXY2 a matrix of size 3 x nPoints, where 1st column is x,
     * second is y.
     * @return 
     */
    private double[][] createFundamentalMatrix(Matrix normXY1, 
        Matrix normXY2) {
        
        if (normXY1 == null) {
            throw new IllegalArgumentException("normXY1 cannot be null");
        }
        if (normXY2 == null) {
            throw new IllegalArgumentException("normXY2 cannot be null");
        }
        if (normXY1.getArray()[0].length != normXY2.getArray()[0].length) {
            throw new IllegalArgumentException(
            "the number of columns in normXY1 != number of rows in normXY2");
        }
        
        int nXY1 = normXY1.getArray()[0].length;
        
        /*
        (2) each row in matrix A:
            x_1*x_2, x_1*y_2, x_1, y_1*x_2, y_1*y_2, y_1, x_2, y_2, 1
        */
        double[][] a = new double[9][nXY1];
        for (int i = 0; i < 9; i++) {
            a[i] = new double[nXY1];
        }
        for (int i = 0; i < nXY1; i++) {
            double x1 = normXY1.get(0, i);
            double x2 = normXY2.get(0, i);
            double y1 = normXY1.get(1, i);
            double y2 = normXY2.get(1, i);
            a[0][i] = x1 * x2;
            a[1][i] = x1 * y2;
            a[2][i] = x1;
            a[3][i] = y1 * x2;
            a[4][i] = y1 * y2;
            a[5][i] = y1;
            a[6][i] = x2;
            a[7][i] = y2;
            a[8][i] = 1;
        }
        
        return a;
    }

    /**
     * calculate the epipoles of the fundamental matrix and return them as
     * an array with left epipole in column 0 and right epipole in column 1.
     * @param fundamentalMatrix
     * @return 
     */
    private double[][] calculateEpipoles(Matrix fundamentalMatrix) {
        /*
        epipoles: 
             [U,D,V] = svd(denormalized FundamentalMatrix);
             e1 = last column of V divided by it's last item
             e2 = last column of U divided by it's last item
        */
        SingularValueDecomposition svdE = fundamentalMatrix.svd();
        
        double[] e1 = svdE.getV().getArray()[2];
        double e1Div = e1[2];
        for (int i = 0; i < e1.length; i++) {
            e1[i] /= e1Div;
        }
        double[] e2 = svdE.getU().getArray()[2];
        double e2Div = e2[2];
        for (int i = 0; i < e2.length; i++) {
            e2[i] /= e2Div;
        }
        
        StringBuilder sb = new StringBuilder("e1=");
        for (int i = 0; i < e1.length; i++) {
            sb.append(e1[i]).append(" ");
        }
        log.info(sb.toString());
        
        sb = new StringBuilder("e2=");
        for (int i = 0; i < e2.length; i++) {
            sb.append(e2[i]).append(" ");
        }
        log.info(sb.toString());
        
        double[][] e = new double[2][];
        e[0] = e1;
        e[1] = e2;
        
        return e;
    }

    private double[][] calculateRightEpipolarLines() {
        
        /* calculate right epipolar lines
        F * leftPoint
        */
        int nPoints = leftXY.getArray()[0].length;
        double[][] rightEpipolarLines = new double[3][nPoints];
        for (int i = 0; i < nPoints; i++) {
            rightEpipolarLines[i] = new double[3];
        }
        for (int i = 0; i < nPoints; i++) {
            
            double[] leftPoint = new double[]{leftXY.get(0, i), 
                leftXY.get(1, i), 1};
            
            rightEpipolarLines[i] = MatrixUtil.multiply(
                fundamentalMatrix.getArray(), leftPoint);
        }
 
        return rightEpipolarLines;
    }
    
    private double[][] calculateLeftEpipolarLines() {
        
        /* calculate left epipolar lines
        F^T * rightPoint
        */
        
        Matrix fundamentalMatrixTranspose = fundamentalMatrix.transpose();
        
        int nPoints = rightXY.getArray()[0].length;
        double[][] leftEpipolarLines = new double[3][nPoints];
        for (int i = 0; i < nPoints; i++) {
            leftEpipolarLines[i] = new double[3];
        }
        for (int i = 0; i < nPoints; i++) {
            
            double[] rightPoint = new double[]{rightXY.get(0, i), 
                rightXY.get(1, i), 1};
            
            leftEpipolarLines[i] = MatrixUtil.multiply(
                fundamentalMatrixTranspose.getArray(), 
                rightPoint);
        }
 
        return leftEpipolarLines;
    }
    
    //TODO: put this in aspect
    void drawEpipolarLinesOnImage(Image image, double[][] epipolarLines) {
        
            /*
            [ a b c ]
            y = - (a/b) * x - (c/b)
            */
        
    }
    
    //TODO: put this in aspect
    void drawPointsOnImage(Image image, Matrix xyPoints) {
        
    }
    
    private static class NormalizedXY {

        /**
         * 3 dimensional matrix, with column 0 being x, column 1 being y,
         * and the last column is place holder 1's
         */
        private Matrix xy = null;
        
        private double[] centroidXY = null;
        
        private Matrix normalizationMatrix = null;

        /**
         * @return the centroidXY
         */
        public double[] getCentroidXY() {
            return centroidXY;
        }

        /**
         * @param centroidXY the centroidXY to set
         */
        public void setCentroidXY(double[] centroidXY) {
            this.centroidXY = centroidXY;
        }

        /**
         * @return the factor
         */
        public Matrix getNormalizationMatrix() {
            return normalizationMatrix;
        }

        /**
         * @param matrix holding the scale and offsets to apply to x, y
         */
        public void setNormMatrix(Matrix normMatrix) {
            this.normalizationMatrix = normMatrix;
        }

        /**
         * @return the xy
         */
        public Matrix getXy() {
            return xy;
        }

        /**
         * @param xy the xy to set
         */
        public void setXy(Matrix xy) {
            this.xy = xy;
        }
    }
    
}
