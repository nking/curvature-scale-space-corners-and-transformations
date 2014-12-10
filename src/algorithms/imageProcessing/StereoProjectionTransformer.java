package algorithms.imageProcessing;

import algorithms.util.PairFloatArray;
import Jama.*;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import thirdparty.HungarianAlgorithm;

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
         F = U * diag([D(1,1) D(2,2) 0]) * V^T; 
 
 (5) denormalize the fundamental matrix
     The related part of the normalization equation: inv(T_2) * F * inv(T_1)
     so denormalizing is:
     
         F = (T_1)^T * F * T_2
   
 (6) estimate the error in the fundamental matrix by calculating epipolar
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
    private Matrix epipolarLinesInRight = null;
    
    /**
     * each row is an epipolar line in the left image.
     * Each column corresponds to a point in leftXY and rightXY which are
     * in the same column.
     */
    private Matrix epipolarLinesInLeft = null;
   
    public void calculateEpipolarProjection(PairFloatArray pointsLeftXY, 
        PairFloatArray pointsRightXY) {
        
        if (pointsLeftXY == null) {
            throw new IllegalArgumentException("refactorLeftXY cannot be null");
        }
        if (pointsRightXY == null) {
            throw new IllegalArgumentException("refactorRightXY cannot be null");
        }
        if (pointsLeftXY.getN() != pointsRightXY.getN()) {
            throw new IllegalArgumentException(
                "refactorLeftXY and refactorRightXY must be same size");
        }
        
        if (pointsLeftXY.getN() < 8) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the 8-point problem requires 8 or more points." 
                + " refactorLeftXY.n=" + pointsLeftXY.getN());
        }
        
        /*
        note that the Jama matrix convention is new Matrix(mRows, nCols)
        */
        
        leftXY = rewriteInto3ColumnMatrix(pointsLeftXY);
        
        rightXY = rewriteInto3ColumnMatrix(pointsRightXY);
        
        fundamentalMatrix = calculateFundamentalMatrix(leftXY, rightXY)
            .transpose();
     
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
    
    /**
     * NOTE: this method should only be used for comparison.  Prefer 
     * calculateEpipolarProjection().
     * 
     * @param pointsLeftXY
     * @param pointsRightXY 
     */
    public void calculateEpipolarProjectionWithoutNormalization(
        PairFloatArray pointsLeftXY, PairFloatArray pointsRightXY) {
        
        if (pointsLeftXY == null) {
            throw new IllegalArgumentException("refactorLeftXY cannot be null");
        }
        if (pointsRightXY == null) {
            throw new IllegalArgumentException("refactorRightXY cannot be null");
        }
        if (pointsLeftXY.getN() != pointsRightXY.getN()) {
            throw new IllegalArgumentException(
                "refactorLeftXY and refactorRightXY must be same size");
        }
        
        log.warning("NOTE:  consider using calculateEpipolarProjection instead");
        
        if (pointsLeftXY.getN() < 8) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the 8-point problem requires 8 or more points." 
                + " refactorLeftXY.n=" + pointsLeftXY.getN());
        }
        
        /*
        note that the Jama matrix convention is new Matrix(mRows, nCols)
        */
        
        leftXY = rewriteInto3ColumnMatrix(pointsLeftXY);
        
        rightXY = rewriteInto3ColumnMatrix(pointsRightXY);
        
        fundamentalMatrix = calculateFundamentalMatrixWithoutNormalization(
            leftXY, rightXY).transpose();
     
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
        
        return calculateFundamentalMatrix(normalizedXY1, normalizedXY2);
    }
    
    Matrix calculateFundamentalMatrix(NormalizedXY normalizedXY1, 
        NormalizedXY normalizedXY2) {
        
        //build the fundamental matrix
        Matrix aMatrix = new Matrix(createFundamentalMatrix(
            normalizedXY1.getXy(), normalizedXY2.getXy()));

        /*
        compute linear least square solution:
            solve A = U * D * V^T   for A*f = [..x...]*f = 0
            A has rank 8.  f has rank 2.
        
        calculate [U,D,V] from svd(A):
           result has mRows = number of data points
                      nCols = 9
        */
        SingularValueDecomposition svd = aMatrix.svd();

        // creates U as 9 x nXY1 matrix
        //         D as length 9 array
        //         V as 9 x 9 matrix
        
        // mRows = 9; nCols = 9
        Matrix V = svd.getV();
        
        assert(V.getColumnDimension() == 9);
        assert(V.getRowDimension() == 9);
        
        // reshape V to 3x3
        
        int vNCols = V.getColumnDimension();
        
        double[][] ff = new double[3][3];
        for (int i = 0; i < 3; i++) {
            ff[i] = new double[3];
            ff[i][0] = V.get((i * 3) + 0, vNCols - 1);
            ff[i][1] = V.get((i * 3) + 1, vNCols - 1);
            ff[i][2] = V.get((i * 3) + 2, vNCols - 1);
        }
        Matrix fMatrix = new Matrix(ff);
        
        /* make the fundamental matrix have a rank of 2
        by performing a svd and then reconstructing with the two largest 
        singular values.
            [U,D,V] = svd(F,0);
        
        From [U,D,V] we create:
            F = U * diag([D(1,1) D(2,2) 0]) * V^T, where V^T is V transposed.
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
             F = dot(U, dot(diag(D),V^T))
        */               
        Matrix dDotV = d.times(svd.getV().transpose());
        
        // 3x3        
        Matrix theFundamentalMatrix = svd.getU().times(dDotV);
        
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
        
        denormFundamentalMatrix = denormFundamentalMatrix.times(
            1./denormFundamentalMatrix.get(2, 2));
        
        return denormFundamentalMatrix;
    }   
    
    /**
     * calculate the fundamental matrix without normalization.  Note, this
     * method should not be used, but is present for comparing solutions
     * during tests.
     * @param matchedXY1
     * @param matchedXY2
     * @return 
     */
    Matrix calculateFundamentalMatrixWithoutNormalization(Matrix matchedXY1, 
        Matrix matchedXY2) {
        
        //build the fundamental matrix
        Matrix aMatrix = new Matrix(createFundamentalMatrix(
            matchedXY1, matchedXY2));

        /*
        compute linear least square solution:
            solve A = U * D * V^T   for A*f = [..x...]*f = 0
            A has rank 8.  f has rank 2.
        
        calculate [U,D,V] from svd(A):
           result has mRows = number of data points
                      nCols = 9
        */
        SingularValueDecomposition svd = aMatrix.svd();

        // creates U as 9 x nXY1 matrix
        //         D as length 9 array
        //         V as 9 x 9 matrix
        
        // mRows = 9; nCols = 9
        Matrix V = svd.getV();
        int vNCols = V.getColumnDimension();
        assert(V.getColumnDimension() == 9);
        assert(V.getRowDimension() == 9);
        // reshape it to 3x3
        double[][] ff = new double[3][3];
        for (int i = 0; i < 3; i++) {
            ff[i] = new double[3];
            ff[i][0] = V.get((i * 3) + 0, vNCols - 1);
            ff[i][1] = V.get((i * 3) + 1, vNCols - 1);
            ff[i][2] = V.get((i * 3) + 2, vNCols - 1);
        }
        Matrix fMatrix = new Matrix(ff);
        
        /* make the fundamental matrix have a rank of 2
           by performing a svd and then reconstructing with the two largest 
           singular values.
              [U,D,V] = svd(F,0);
        
           then from [U,D,V], create F:
              F = U * diag([D(1,1) D(2,2) 0]) * V^T;
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
             F = dot(U, dot(diag(D),V^T))
        */
        Matrix dDotV = d.times(svd.getV().transpose());
        
        // 3x3        
        Matrix theFundamentalMatrix = svd.getU().times(dDotV);  
        
        theFundamentalMatrix = theFundamentalMatrix.times(
            1./theFundamentalMatrix.get(2, 2));
        
        return theFundamentalMatrix;
    } 

    /**
     normalize the x,y coordinates as recommended by Hartley 1997 and return
     the matrix and coordinates.
       
     * @param xyPair
     * @return 
     */
    NormalizedXY normalize(Matrix xy) {
        
        /*
        uTransposed = T * u 
        uTransposed * inv(T) = u
        
                uTransposed_2^T * inv(T_2) * F * inv(T_1) * uTransposed_1
        
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
        x_1_0  x_1_1  x_1_2  x_1_3
        y_1_0  y_1_1  y_1_2  y_1_3
        1      1      1      1
        
        t00     t01(=0)  t02
        t10(=0) t11      t12
        0        0        1
        
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
                
        Matrix normXY = new Matrix(MatrixUtil.dot(tMatrix, xy));
              
        NormalizedXY normalizedXY = new NormalizedXY();
        normalizedXY.setCentroidXY(centroidXY);
        normalizedXY.setNormMatrix(tMatrix);
        normalizedXY.setXy(normXY);
        
        return normalizedXY;
    }
    
    /**
     * write a matrix of size mRows = 3, nCols = xyPairs.getN()
     * @param xyPairs
     * @return 
     */
    public Matrix rewriteInto3ColumnMatrix(PairFloatArray xyPairs) {
        
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
        
        // matrix of size mRows x nCols
        
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
    double[][] createFundamentalMatrix(Matrix normXY1, Matrix normXY2) {
        
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
        
        // mRows = 484;  nCols = 9
        
        /*
        (2) each row in matrix A:
            x_1*x_2,  x_1*y_2,  x_1,  y_1*x_2,  y_1*y_2,  y_1,  x_2,  y_2,  1
        */
        double[][] a = new double[nXY1][9];
        for (int i = 0; i < nXY1; i++) {
            a[i] = new double[9];
            double x1 = normXY1.get(0, i);
            double x2 = normXY2.get(0, i);
            double y1 = normXY1.get(1, i);
            double y2 = normXY2.get(1, i);
            a[i][0] = x1 * x2;
            a[i][1] = x1 * y2;
            a[i][2] = x1;
            a[i][3] = y1 * x2;
            a[i][4] = y1 * y2;
            a[i][5] = y1;
            a[i][6] = x2;
            a[i][7] = y2;
            a[i][8] = 1;
        }
        
        return a;
    }

    /**
     * calculate the epipoles of the fundamental matrix and return them as
     * an array with left epipole in column 0 and right epipole in column 1.
     * @param fundamentalMatrix
     * @return 
     */
    double[][] calculateEpipoles(Matrix fundamentalMatrix) {
        
        /*
        The representation of lines in homogeneous projective coordinates
        is:   line a*x + b*y + c = 0
            | a |
            | b |
            | c |
        The line can be rewritten in slope, intercept form:
            y = intercept + slope * x
              = -(c/b) - slope*(a/b)*x
        
        written as homogenization form of lines:
            | -a/b |
            | -c/b |
        
        
        From u_2^T * F * u_1 = 0 
        
        The right epipole, 
            e_right = 
        
        
        epipoles: 
             [U,D,V] = svd(denormalized FundamentalMatrix);
             e1 = last column of V divided by it's last item
             e2 = last column of U divided by it's last item
                
        */
        SingularValueDecomposition svdE = fundamentalMatrix.svd();
        Matrix V = svdE.getV().transpose();
        double[] e1 = V.getArray()[2];
        double e1Div = e1[2];
        for (int i = 0; i < e1.length; i++) {
            e1[i] /= e1Div;
        }
        Matrix U = svdE.getU();
        double[] e2 = U.getArray()[2];
        double e2Div = e2[2];
        for (int i = 0; i < e2.length; i++) {
            e2[i] /= e2Div;
        }
        
        StringBuilder sb = new StringBuilder("e1=");
        for (int i = 0; i < e1.length; i++) {
            sb.append(e1[i]).append(" ");
        }
        log.fine(sb.toString());
        
        sb = new StringBuilder("e2=");
        for (int i = 0; i < e2.length; i++) {
            sb.append(e2[i]).append(" ");
        }
        log.fine(sb.toString());
        
        double[][] e = new double[2][];
        e[0] = e1;
        e[1] = e2;
        
        return e;
    }

    private Matrix calculateRightEpipolarLines() {
        
        /* calculate right epipolar lines
        F * leftPoint
        */
        
        Matrix m = fundamentalMatrix.times(leftXY);
        
        return m;
    }
    
    private Matrix calculateLeftEpipolarLines() {
        
        /* calculate left epipolar lines
        F^T * rightPoint
        */
        
        Matrix fundamentalMatrixTranspose = fundamentalMatrix.transpose();
        
        Matrix m = fundamentalMatrixTranspose.times(rightXY);
        
        return m;
    }

    /**
    calculate the distance of points in set 2 from the epipolar lines and 
    estimate the average and standard deviation and count the number within 
    tolerance.
     
     * @param tolerance
     * @return 
     */
    public StereoProjectionTransformerFit evaluateFitForRight(double tolerance) {
        
        return evaluateFitInRightImageForMatchedPoints(leftXY, rightXY, 
            tolerance);
    }
    
    public static class NormalizedXY {

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
         * @param normMatrix holding the scale and offsets to apply to x, y
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
    
    public double[] getLeftEpipole() {
        return leftEpipole;
    }
    public double[] getRightEpipole() {
        return rightEpipole;
    }
    public Matrix getEpipolarLinesInRight() {
        return epipolarLinesInRight;
    }
    public Matrix getEpipolarLinesInLeft() {
        return epipolarLinesInLeft;
    }
    
    public PairIntArray getEpipolarLineInLeft(int imgWidth, int imgHeight, 
        int pointNumber) {
        
        return getEpipolarLine(epipolarLinesInLeft, imgWidth, imgHeight, 
             pointNumber);
    }
    
    public PairIntArray getEpipolarLineInRight(int imgWidth, int imgHeight, 
        int pointNumber) {
        
        return getEpipolarLine(epipolarLinesInRight, imgWidth, imgHeight, 
            pointNumber);
    }
     
    PairIntArray getEpipolarLine(Matrix epipolarLines, int imgWidth, 
        int imgHeight, int pointNumber) {
        
        int n = imgWidth/10;
        
        PairIntArray line = new PairIntArray(n);
        
        double a = epipolarLines.get(0, pointNumber);
        double b = epipolarLines.get(1, pointNumber);
        double c = epipolarLines.get(2, pointNumber);
        boolean isHoriz = (Math.abs(a/b) <= 1.0);
        if (isHoriz) {
            for (int x = 0; x < imgWidth; x++) {
                //y = - (a/b) * x - (c/b)
                double y = (c + (a * (double)x)) / (-b);
                line.add(x, (int) Math.round(y));
            }
        } else {
            for (int y = 0; y < imgHeight; y++) {
                //y = - (a/b) * x - (c/b)
                //y+(c/b) = - (a/b) * x 
                // ==> x = (-b/a) * (y+(c/b)) = y*(-b/a) - (c/a)
                double x = -(c + (b * (double)y))/a;
                line.add((int) Math.round(x), y);
            }
        }
  
        return line;
    }
    
    private double[] getEpipolarLineYEndpoints(Matrix epipolarLines, int xBegin,
        int xEnd, int pointNumber) {
        
        double a = epipolarLines.get(0, pointNumber);
        double b = epipolarLines.get(1, pointNumber);
        double c = epipolarLines.get(2, pointNumber);
        
        //y = - (a/b) * x - (c/b)
        double yBegin = (c + (a * (double)xBegin)) / (-b);
        
        double yEnd = (c + (a * (double)xEnd)) / (-b);
  
        return new double[]{yBegin, yEnd};
    }
    
    private double[] getEpipolarLineYEndpoints(Matrix epipolarLines, float xBegin,
        float xEnd, int pointNumber) {
        
        double a = epipolarLines.get(0, pointNumber);
        double b = epipolarLines.get(1, pointNumber);
        double c = epipolarLines.get(2, pointNumber);
        
        //y = - (a/b) * x - (c/b)
        double yBegin = (c + (a * xBegin)) / (-b);
        
        double yEnd = (c + (a * xEnd)) / (-b);
  
        return new double[]{yBegin, yEnd};
    }
    
    Matrix calculateEpipolarRightLines(Matrix points) {
        return fundamentalMatrix.times(points);
    }
    
     Matrix calculateEpipolarLeftLines(Matrix points) {
        return fundamentalMatrix.transpose().times(points);
    }
    
     /**
      * evaluate the fit as the distance of points in the left image from the
      * projected right epipolar lines if the distance is within tolerance.
      * Note that to include all points regardless of tolerance, set tolerance
      * to a very high number such as Double.MAX_VALUE.
      * 
      * @param leftImagePoints
      * @param rightImagePoints
      * @param tolerance
      * @return 
      */
    public StereoProjectionTransformerFit evaluateFit(
        PairFloatArray leftImagePoints, PairFloatArray rightImagePoints, 
        double tolerance) {
        
        /*
        comparing the points in left image with the epipolar lines projected
        from the right image points.
        */
        
        Matrix theRightPoints = rewriteInto3ColumnMatrix(rightImagePoints);
        
        Matrix theLeftEpipolarLines = calculateEpipolarLeftLines(theRightPoints);
        
        int[] minMaxLineXEndpoints = getMinMaxX(leftImagePoints);
        float lineX0 = minMaxLineXEndpoints[0];
        float lineX1 = minMaxLineXEndpoints[1];
        
        int nPoints1 = leftImagePoints.getN();
        int nPoints2 = rightImagePoints.getN();
               
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        double[] diffs = new double[nPoints1];
        int nMatched = 0;
        long diffSum = 0;
            
        boolean useBipartiteMinCost = true;
        
        if (useBipartiteMinCost) {
            
            // === using an optimal bipartite min cost matching ~O(N^4) ====
            // make the cost matrix
            float[][] diffsAsCost = new float[nPoints2][nPoints1];
            
            // the algorithm modifies diffsAsCost, so make a copy
            float[][] diffsAsCostCopy = new float[nPoints2][nPoints1];

            for (int i = 0; i < nPoints2; i++) {

                diffsAsCost[i] = new float[nPoints1];
                diffsAsCostCopy[i] = new float[nPoints1];

                double[] epipolarLineYEndPoints = getEpipolarLineYEndpoints(
                    theLeftEpipolarLines, lineX0, lineX1, i);

                float lineY0 = (float)epipolarLineYEndPoints[0];
                float lineY1 = (float)epipolarLineYEndPoints[1];

                for (int j = 0; j < leftImagePoints.getN(); j++) {

                    double dist = curveHelper.distanceFromPointToALine(
                        lineX0, lineY0, lineX1, lineY1, 
                        leftImagePoints.getX(j), leftImagePoints.getY(j));

                    diffsAsCost[i][j] = (float)dist;
                    diffsAsCostCopy[i][j] = (float)dist;
                }
            }

            HungarianAlgorithm b = new HungarianAlgorithm();
            int[][] match = b.computeAssignments(diffsAsCostCopy);

            assert(match.length == nPoints2);
            
log.info("set1 n=" + nPoints1 + " set2 n=" + nPoints2);

            
            for (int i = 0; i < match.length; i++) {
                int idx2 = match[i][0];
                int idx1 = match[i][1];
                if (idx1 == -1 || idx2 == -1) {
                    continue;
                }
                double diff = diffsAsCost[idx2][idx1];
                if (diff < tolerance) {
                    diffs[nMatched] = diff;
                    diffSum += diff;
                    nMatched++;
                }
            }
            
        } else {
        
            // ==== using a greedy non-optimal min cost matching ~O(N^2) ======
            Set<Integer> chosen = new HashSet<Integer>();
            for (int i = 0; i < nPoints2; i++) {
                double[] epipolarLineYEndPoints = getEpipolarLineYEndpoints(
                    theLeftEpipolarLines, lineX0, lineX1, i);
                float lineY0 = (float)epipolarLineYEndPoints[0];
                float lineY1 = (float)epipolarLineYEndPoints[1];
                double minDiff = Double.MAX_VALUE;
                int minIdx = -1;
                for (int j = 0; j < leftImagePoints.getN(); j++) {
                    if (chosen.contains(Integer.valueOf(j))) {
                        continue;
                    }
                    double dist = curveHelper.distanceFromPointToALine(
                        lineX0, lineY0, lineX1, lineY1, 
                        leftImagePoints.getX(j), leftImagePoints.getY(j));
                    if (dist < minDiff) {
                        minDiff = dist;
                        minIdx = j;
                    }
                }
                if ((minDiff < Double.MAX_VALUE) && (minDiff < tolerance)) {
                    diffs[nMatched] = minDiff;
                    diffSum += minDiff;
                    nMatched++;
                    chosen.add(Integer.valueOf(minIdx));
                }
            }
        }
        
        double avgDist = diffSum/(double)nMatched;
        
        diffSum = 0;
        for (int i = 0; i < nMatched; i++) {
            double d = diffs[i] - avgDist;
            diffSum += (d * d);
        }
        double stdDevDist = Math.sqrt(diffSum/((double)nMatched - 1));
        
        StereoProjectionTransformerFit fit = new StereoProjectionTransformerFit(
            nMatched, tolerance, avgDist, stdDevDist);
        
        return fit;
    }

    /**
      * evaluate the fit in the right image as the distance of points there
      * from the projected epipolar lines (created from the left points).
      * The distances are only added if they are within tolerance.
      * Note that to include all points regardless of tolerance, set tolerance
      * to a very high number such as Double.MAX_VALUE.
      * 
      * @param leftImagePoints
      * @param rightImagePoints
      * @param tolerance
      * @return 
      */
    public StereoProjectionTransformerFit evaluateFitInRightImage(
        PairFloatArray leftImagePoints, PairFloatArray rightImagePoints, 
        double tolerance) {
        
        Matrix theLeftPoints = rewriteInto3ColumnMatrix(leftImagePoints);
        
        Matrix theRightPoints = rewriteInto3ColumnMatrix(rightImagePoints);
       
        return evaluateFitInRightImage(theLeftPoints, theRightPoints, tolerance);
    }
    
    /**
      * evaluate the fit in the right image as the distance of points there
      * from the projected epipolar lines (created from the left points).
      * The distances are only added if they are within tolerance.
      * Note that to include all points regardless of tolerance, set tolerance
      * to a very high number such as Double.MAX_VALUE.
      * 
      * @param leftImagePoints
      * @param rightImagePoints
      * @param tolerance
      * @return 
      */
    public StereoProjectionTransformerFit 
        evaluateFitInRightImageForMatchedPoints(
        PairFloatArray leftImagePoints, PairFloatArray rightImagePoints,
        double tolerance) {
        
        Matrix theLeftPoints = rewriteInto3ColumnMatrix(leftImagePoints);
        
        Matrix theRightPoints = rewriteInto3ColumnMatrix(rightImagePoints);
       
        return evaluateFitInRightImageForMatchedPoints(theLeftPoints, 
            theRightPoints, tolerance);
    }
        
    public StereoProjectionTransformerFit evaluateRightForMatchedAndStoreOutliers(
        PairFloatArray leftImagePoints, PairFloatArray rightImagePoints,
        float factor, double minRemoval, LinkedHashSet<Integer> skipIndexes) {
        
        Matrix theLeftPoints = rewriteInto3ColumnMatrix(leftImagePoints);
        
        Matrix theRightPoints = rewriteInto3ColumnMatrix(rightImagePoints);
       
        return evaluateRightForMatchedAndStoreOutliers(theLeftPoints, 
            theRightPoints, factor, minRemoval, skipIndexes);
    }

    /**
      * evaluate the fit in the right image as the distance of points there
      * from the projected epipolar lines (created from the left points).
      * The distances are only added if they are within tolerance.
      * Note that to include all points regardless of tolerance, set tolerance
      * to a very high number such as Double.MAX_VALUE.
      * 
      * @param leftPoints
      * @param rightPoints
      * @param tolerance
      * @return 
      */
    public StereoProjectionTransformerFit evaluateFitInRightImage(
        Matrix leftPoints, Matrix rightPoints, 
        double tolerance) {
       
        Matrix theRightEpipolarLines = calculateEpipolarRightLines(leftPoints);
      
        int[] minMaxLineXEndpoints = getMinMaxX(rightPoints);
        float lineX0 = minMaxLineXEndpoints[0];
        float lineX1 = minMaxLineXEndpoints[1];
        
        int nPointsLeft = leftPoints.getColumnDimension();
        int nPointsRight = rightPoints.getColumnDimension();
      
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        double[] diffs = new double[nPointsLeft];
        int nMatched = 0;
        long diffSum = 0;
            
        boolean useBipartiteMinCost = true;
        
        if (useBipartiteMinCost) {
            
            // === using an optimal bipartite min cost matching ~O(N^4) ====
            // make the cost matrix
            float[][] diffsAsCost = new float[nPointsLeft][nPointsRight];

            // the algorithm modifies diffsAsCost, so make a copy
            float[][] diffsAsCostCopy = new float[nPointsLeft][nPointsRight];

            for (int i = 0; i < nPointsLeft; i++) {

                diffsAsCost[i] = new float[nPointsRight];
                diffsAsCostCopy[i] = new float[nPointsRight];

                double[] epipolarLineYEndPoints = getEpipolarLineYEndpoints(
                    theRightEpipolarLines, lineX0, lineX1, i);

                float lineY0 = (float)epipolarLineYEndPoints[0];
                float lineY1 = (float)epipolarLineYEndPoints[1];

                for (int j = 0; j < nPointsRight; j++) {

                    double dist = curveHelper.distanceFromPointToALine(
                        lineX0, lineY0, lineX1, lineY1, 
                        (float)rightPoints.get(0, j), 
                        (float)rightPoints.get(1, j));

                    diffsAsCost[i][j] = (float)dist;
                    diffsAsCostCopy[i][j] = (float)dist;
                }
            }

            HungarianAlgorithm b = new HungarianAlgorithm();
            int[][] match = b.computeAssignments(diffsAsCostCopy);

            assert(match.length == nPointsLeft);
            
            for (int i = 0; i < match.length; i++) {
                int idxLeft = match[i][0];
                int idxRight = match[i][1];
                if (idxLeft == -1 || idxRight == -1) {
                    continue;
                }
                double diff = diffsAsCost[idxLeft][idxRight];
                if (diff < tolerance) {
                    diffs[nMatched] = diff;
                    diffSum += diff;
                    nMatched++;
                }
            }
            
        } else {
        
            // ==== using a greedy non-optimal min cost matching ~O(N^2) ======
            Set<Integer> chosen = new HashSet<Integer>();
            for (int i = 0; i < nPointsLeft; i++) {
                double[] epipolarLineYEndPoints = getEpipolarLineYEndpoints(
                    theRightEpipolarLines, lineX0, lineX1, i);
                float lineY0 = (float)epipolarLineYEndPoints[0];
                float lineY1 = (float)epipolarLineYEndPoints[1];
                double minDiff = Double.MAX_VALUE;
                int minIdx = -1;
                for (int j = 0; j < nPointsRight; j++) {
                    if (chosen.contains(Integer.valueOf(j))) {
                        continue;
                    }
                    double dist = curveHelper.distanceFromPointToALine(
                        lineX0, lineY0, lineX1, lineY1, 
                        (float)rightPoints.get(0, j), 
                        (float)rightPoints.get(1, j));
                    if (dist < minDiff) {
                        minDiff = dist;
                        minIdx = j;
                    }
                }
                if ((minDiff < Double.MAX_VALUE) && (minDiff < tolerance)) {
                    diffs[nMatched] = minDiff;
                    diffSum += minDiff;
                    nMatched++;
                    chosen.add(Integer.valueOf(minIdx));
                }
            }
        }
        
        log.fine("n epipolar lines=" + nPointsLeft + " right n=" + nPointsRight 
            + " nMatched=" + nMatched);
        
        double avgDist = diffSum/(double)nMatched;
        
        diffSum = 0;
        for (int i = 0; i < nMatched; i++) {
            double d = diffs[i] - avgDist;
            diffSum += (d * d);
        }
        double stdDevDist = Math.sqrt(diffSum/((double)nMatched - 1));
        
        StereoProjectionTransformerFit fit = new StereoProjectionTransformerFit(
            nMatched, tolerance, avgDist, stdDevDist);
        
        return fit;
    }
    
    /**
      * evaluate the fit in the right image as the distance of points there
      * from the projected epipolar lines (created from the left points).
      * The distances are only added if they are within tolerance.
      * Note that to include all points regardless of tolerance, set tolerance
      * to a very high number such as Double.MAX_VALUE.
      * 
      * @param leftPoints
      * @param rightPoints
      * @param tolerance
      * @return 
      */
    public StereoProjectionTransformerFit evaluateFitInRightImageForMatchedPoints(
        Matrix leftPoints, Matrix rightPoints, 
        double tolerance) {
               
        if (leftPoints == null) {
            throw new IllegalArgumentException("leftPoints cannot be null");
        }
        
        if (rightPoints == null) {
            throw new IllegalArgumentException("rightPoints cannot be null");
        }
         
        int nPointsLeft = leftPoints.getColumnDimension();
        int nPointsRight = rightPoints.getColumnDimension();
        
        if (nPointsLeft != nPointsRight) {
            throw new IllegalArgumentException("point lists must have same "
                + " length since these are matched. " +
                " For unmatched lists, use evaluateFitInRightImage()");
        }
      
        Matrix theRightEpipolarLines = calculateEpipolarRightLines(leftPoints);
      
        int[] minMaxLineXEndpoints = getMinMaxX(rightPoints);
        float lineX0 = minMaxLineXEndpoints[0];
        float lineX1 = minMaxLineXEndpoints[1];
       
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        double[] diffs = new double[nPointsLeft];
        int nMatched = 0;
        long diffSum = 0;
                    
        for (int i = 0; i < nPointsLeft; i++) {

            double[] epipolarLineYEndPoints = getEpipolarLineYEndpoints(
                theRightEpipolarLines, lineX0, lineX1, i);

            float lineY0 = (float)epipolarLineYEndPoints[0];
            float lineY1 = (float)epipolarLineYEndPoints[1];

            double dist = curveHelper.distanceFromPointToALine(
                lineX0, lineY0, lineX1, lineY1, 
                (float)rightPoints.get(0, i), 
                (float)rightPoints.get(1, i));

            if (dist < tolerance) {
                diffs[nMatched] = dist;
                diffSum += dist;
                nMatched++;
            }
        }
        
        double avgDist = diffSum/(double)nMatched;
        
        diffSum = 0;
        for (int i = 0; i < nMatched; i++) {
            double d = diffs[i] - avgDist;
            diffSum += (d * d);
        }
        double stdDevDist = Math.sqrt(diffSum/((double)nMatched - 1));
        
        StereoProjectionTransformerFit fit = new StereoProjectionTransformerFit(
            nMatched, tolerance, avgDist, stdDevDist);
        
        return fit;
    }
    
    /**
     * evaluate the fit of projection by calculating the difference between
     * the right points and the epipolar lines that are calculated from the
     * left points.  After the average distance and standard deviation are
     * learned for the right points, those that are greater than 
     * average distance + factor * standard deviation are added to the list
     * skipIndexes.  Note that any points already in skipIndexes are
     * skipped during the evaluation.
     * @param leftPoints input matched points from left image, excepting those
     * with indexes in skipIndexes.
     * @param rightPoints input matched points from right image, excepting those
     * with indexes in skipIndexes.
     * @param factor
     * @param minRemoval minimum distance above which an outlier distance
     * can be applied
     * @param skipIndexes existing indexes to skip in evaluation and to append 
     * to with outliers from this evaluation.
     * @return 
     */
    public StereoProjectionTransformerFit evaluateRightForMatchedAndStoreOutliers(
        Matrix leftPoints, Matrix rightPoints, 
        float factor, double minRemoval, LinkedHashSet<Integer> skipIndexes) {
               
        if (leftPoints == null) {
            throw new IllegalArgumentException("leftPoints cannot be null");
        }
        
        if (rightPoints == null) {
            throw new IllegalArgumentException("rightPoints cannot be null");
        }
        
        if (skipIndexes == null) {
            throw new IllegalArgumentException("skipIndexes cannot be null");
        }
         
        int nPointsLeft = leftPoints.getColumnDimension();
        int nPointsRight = rightPoints.getColumnDimension();
        
        if (nPointsLeft != nPointsRight) {
            throw new IllegalArgumentException("point lists must have same "
                + " length since these are matched. " +
                " For unmatched lists, use evaluateFitInRightImage()");
        }
      
        Matrix theRightEpipolarLines = calculateEpipolarRightLines(leftPoints);
      
        int[] minMaxLineXEndpoints = getMinMaxX(rightPoints);
        float lineX0 = minMaxLineXEndpoints[0];
        float lineX1 = minMaxLineXEndpoints[1];
       
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        double[] diffs = new double[nPointsLeft];
        int nMatched = 0;
        long diffSum = 0;
                    
        for (int i = 0; i < nPointsLeft; i++) {

            double[] epipolarLineYEndPoints = getEpipolarLineYEndpoints(
                theRightEpipolarLines, lineX0, lineX1, i);

            float lineY0 = (float)epipolarLineYEndPoints[0];
            float lineY1 = (float)epipolarLineYEndPoints[1];

            double dist = curveHelper.distanceFromPointToALine(
                lineX0, lineY0, lineX1, lineY1, 
                (float)rightPoints.get(0, i), 
                (float)rightPoints.get(1, i));

            diffs[i] = dist;
            
            if (!skipIndexes.contains(Integer.valueOf(i))) {
                diffSum += dist;
                nMatched++;
            }
        }
        
        double avgDist = diffSum/(double)nMatched;
        
        diffSum = 0;
        for (int i = 0; i < nPointsLeft; i++) {
            if (!skipIndexes.contains(Integer.valueOf(i))) {
                double d = diffs[i] - avgDist;
                diffSum += (d * d);
            }
        }
        
        double stdDevDist = Math.sqrt(diffSum/((double)nMatched - 1));
        
        double comp = avgDist + factor * stdDevDist;
        
        if (comp > minRemoval) {
                    
            for (int i = 0; i < nPointsLeft; i++) {
                Integer idx = Integer.valueOf(i);
                if (!skipIndexes.contains(idx)) {
                    if (diffs[i] > comp) {
                        skipIndexes.add(idx);
                    }
                }
            }
            
        }
        
        StereoProjectionTransformerFit fit = new StereoProjectionTransformerFit(
            nMatched, comp, avgDist, stdDevDist);
        
        return fit;
    }

    public PairIntArray getRightXYInt() {
        return getXYInt(rightXY);
    }
    
    public PairIntArray getLeftXYInt() {                         
        return getXYInt(leftXY);
    }
    
    public int getNumberOfMatches() {
        return leftXY.getColumnDimension();
    }
    
    public PairFloatArray getRightXYFloat() {                         
        return getXYFloat(rightXY);
    }
    
    public PairFloatArray getLeftXYFloat() {
        return getXYFloat(leftXY);
    }
    
    private PairIntArray getXYInt(Matrix leftOrRightXY) {
                
        int nPoints = leftOrRightXY.getColumnDimension();
         
        PairIntArray out = new PairIntArray();
        
        for (int i = 0; i < nPoints; i++) {
            
            float xP = (float) leftOrRightXY.get(0, i);
            float yP = (float) leftOrRightXY.get(1, i);
           
            out.add(Math.round(xP), Math.round(yP));
        }
        
        return out;
    }
    
    private PairFloatArray getXYFloat(Matrix leftOrRightXY) {
                
        int nPoints = leftOrRightXY.getColumnDimension();
         
        PairFloatArray out = new PairFloatArray();
        
        for (int i = 0; i < nPoints; i++) {
            
            float xP = (float) leftOrRightXY.get(0, i);
            float yP = (float) leftOrRightXY.get(1, i);
           
            out.add(xP, yP);
        }
        
        return out;
    }
    
    private int[] getMinMaxX(PairFloatArray points) {
        
        int nPoints = points.getN();
         
        // estimate endpoints of the epipolar line as the min and max x of right
        int xBegin = Integer.MAX_VALUE;
        int xEnd = Integer.MIN_VALUE;
        for (int i = 0; i < nPoints; i++) {
            double xP = points.getX(i);
            int x = (int)Math.round(xP);
            if (x < xBegin) {
                xBegin = x;
            }
            if (x > xEnd) {
                xEnd = x;
            }
        }
        
        return new int[]{xBegin, xEnd};
    }
    
    private int[] getMinMaxX(Matrix points) {
        
        int nPoints = points.getColumnDimension();
         
        // estimate endpoints of the epipolar line as the min and max x of right
        int xBegin = Integer.MAX_VALUE;
        int xEnd = Integer.MIN_VALUE;
        for (int i = 0; i < nPoints; i++) {
            double xP = points.get(0, i);
            int x = (int)Math.round(xP);
            if (x < xBegin) {
                xBegin = x;
            }
            if (x > xEnd) {
                xEnd = x;
            }
        }
        
        return new int[]{xBegin, xEnd};
    }
}
