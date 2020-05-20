package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.features.FeatureComparisonStat;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.matrix.MatrixUtil;
import algorithms.util.PairFloatArray;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
 * class to solve for the epipoles for two images with stereo projection
 * and apply the solution.
 *
 * <pre>
 * The fundamental matrix is the projective solution for transformation
 * between 2 images of the same objects in pixel coordinates.
 * Present below is the solution for having 7 matched points between images
 * and the solution for having 8 or more matched points between the images.
 * Both use numerical conditioning and recipes suggested by Hartley
 * (see reference below).
 *
 * Following:
 *
 * IEEE TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE, VOL. 19,
 * NO. 6, JUNE 1997
 * "In Defense of the Eight-Point Algorithm" by Richard Hartley

 * The 8-point algorithm matrix represents epipolar geometry, and can be
 * used with data from cameras whose characteristics are not known to solve
 * up to the projective transformation.
 * ...a simple transformation (translation and scaling) of the points in 
 * the image before formulating the linear equations leads to an enormous 
 * improvement in the condition of the problem and hence of the stability 
 * of the result...
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
     
    u_1   = x_1_i    x_1_i+1  ...    
            y_1_i    y_1_i+1  ...
            1        1        ...
 
    u_2^T = x_2_i    y_2_i    1
            x_2_i+1  y_2_i+1  1
            ...      ...      ...

    x_1*x_2*F_1_1 + x_1*y_2*F_2_1 + x_1*F_3_1 + y_1*x_2*F_1_2 + y_1*y_2*F_2_2
        + y_1*F_3_2 + x_2*F_1_3 + y_2*F_2_3 + F_3_3 = 0

    A * f = 0

    where A = x_1*x_2, x_1*y_2, x_1, y_1*x_2, y_1*y_2, y_1, x_2, y_2, 1
             (which is each element of column 0 of u_1 dotted separately with
             row 0 of u2_T)
 
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

 The 7-point algorithm is also implemented below and is similar to the
 8-point solution except that is solves for the null space of the fundamental
 matrix and results in one or 3 solutions which can for some geometries
 be reduced to a single solution.
 The nullspace is where Ax=0 in reduced echelon, that is, the free variable rows.
 The normalization and denormalization steps before and following the solution,
 are the same as in the 8-point solution.

 this from comments in the Hartley & Zisserman matlab code vgg_F_from_7pts_2img:
   Solutions for the 7-points are pruned by requirement that
   scalars s in all equations s * cross(e1,u_1) == F*u_2 are positive.
   In case of multiple solutions, F has one dimension
   more such that F(:,:,n) is the n-th solution.

 </pre>
 NOTE:
For "7-point" correspondences, consider implementing MLESAC.
     "MLESAC: A new robust estimator with application to estimating image geometry"
     by P. H. S. Torr and A. Zisserman
     1996 http://www.robots.ox.ac.uk/~vgg/publications/papers/torr00.pdf
 
 * @author nichole
 */
public class EpipolarTransformer {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    private static double eps = 1e-12;

    /**
     * calculate the fundamental matrix for the given matched left and
     * right correspondence of 8 or more matched points.
     *
     * @param pointsLeftXY
     * @param pointsRightXY
     * @return
     */
    public DenseMatrix calculateEpipolarProjection(
        PairIntArray pointsLeftXY,  PairIntArray pointsRightXY) {

        if (pointsLeftXY == null) {
            throw new IllegalArgumentException("pointsLeftXY cannot be null");
        }
        if (pointsRightXY == null) {
            throw new IllegalArgumentException("pointsRightXY cannot be null");
        }
        if (pointsLeftXY.getN() != pointsRightXY.getN()) {
            throw new IllegalArgumentException(
                "pointsLeftXY and pointsRightXY must be same size");
        }

        if (pointsLeftXY.getN() == 7) {
            throw new IllegalArgumentException(
                "for 7 points, use calculateEpipolarProjectionFor7Points");
        }

        if (pointsLeftXY.getN() < 7) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the algorithms require 7 or more points."
                + " pointsLeftXY.n=" + pointsLeftXY.getN());
        }

        return calculateEpipolarProjection(
            rewriteInto3ColumnMatrix(pointsLeftXY),
            rewriteInto3ColumnMatrix(pointsRightXY));
    }

    /**
     * calculate the epipolar projection for a set of 8 or more matched points.
     *
     * @param theLeftXY
     * @param theRightXY
     * @return
     */
    public DenseMatrix calculateEpipolarProjection(
        DenseMatrix theLeftXY, DenseMatrix theRightXY) {

        if (theLeftXY == null) {
            throw new IllegalArgumentException("theLeftXY cannot be null");
        }
        if (theRightXY == null) {
            throw new IllegalArgumentException("refactorRightXY cannot be null");
        }
        if (theLeftXY.numColumns()!= theRightXY.numColumns()) {
            throw new IllegalArgumentException(
                "theLeftXY and theRightXY must be same size");
        }

        if (theLeftXY.numColumns() == 7) {
            throw new IllegalArgumentException(
                "for 7 points, use calculateEpipolarProjectionFor7Points");
        }

        if (theLeftXY.numColumns() < 7) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the algorithms require 7 or more points."
                + " refactorLeftXY.n=" +theLeftXY.numColumns());
        }

        //the matrix convention is [mRows][nCols]

        DenseMatrix fundamentalMatrix = (DenseMatrix)
            calculateFundamentalMatrix(theLeftXY, theRightXY).transpose();

        return fundamentalMatrix;
    }

    protected DenseMatrix calculateFundamentalMatrix(DenseMatrix leftXY,
        DenseMatrix rightXY) {

        //x is xy[0], y is xy[1], xy[2] is all 1's
        NormalizedXY normalizedXY1 = normalize(leftXY);

        NormalizedXY normalizedXY2 = normalize(rightXY);

        return calculateFundamentalMatrix(normalizedXY1, normalizedXY2);
    }

    /*
    Following the 7-point algorithm by R. Hartley and A. Zisserman, 
    in their book "Multiple View Geometry in Computer Vision"

    (1) SVD of matrix A (as is done in 8-point algorithm)
        giving a matrix of rank 7
    (2) The homogeneous system AX = 0 : X is the null space of matrix A.
        The system is nullable because rank 7 < number of columns, 9.

        The nullable system must have a solution other than trivial where
        |A| = 0.

        There should be 9-7=2 linearly independent vectors u1, u2, ... , un-r
        that span the null space of A.

        The right null space of A reduced by SVD is then 2D and the last
        2 columns of V can be extracted and reshaped to [3x3] as F1 and F2.

        A linear convex combination of F1 and F2 form the estimate of F.

        F = α*F1 + (1 − α)*F2  where α is between 0 and 1

        The eigenvalues of F are possible only if the determinant of F is 0.
        The determinant of F is a polynomial function, the characteristic
        polynomial whose degree is the order of the matrix, which is 3 in this
        case. Therefore, the answer(s) to determinant(F) = 0 requires the cubic
        roots of the equation.

        det A = 0 ==> det(α*F1 + (1 − α)*F2) = 0

        because det(F1 + F2) != det(F1) + det(F2), have to step through the
        determinant of the sums, and group the terms by a^3, a^2, a^1, and a^0
        and then solve for the cubic roots as the values of 'a'.

   The matrices multiplied and summed:

    a*ff1[0][0] + (1-a)*ff2[0][0]   a*ff1[0][1] + (1-a)*ff2[0][1]   a*ff1[0][2] + (1-a)*ff2[0][2]
    a*ff1[1][0] + (1-a)*ff2[1][0]   a*ff1[1][1] + (1-a)*ff2[1][1]   a*ff1[1][2] + (1-a)*ff2[1][2]
    a*ff1[2][0] + (1-a)*ff2[2][0]   a*ff1[2][1] + (1-a)*ff2[2][1]   a*ff1[2][2] + (1-a)*ff2[2][2]

    The terms are further grouped below in methods
       calculateCubicRoot...OrderCoefficientFor7Point(ff1, ff2)

    After the cubic root(s) are solved, they are back substituted into :
        Fi = a(i) * FF{1} + (1-a(i)) * FF{2};
    to get the solutions Fi which may be one or 3 solutions.
    */

    /**
     * calculate the epipolar projection for a set of matched points that are
     * at 7 points in length.
     *
     * @param pointsLeftXY
     * @param pointsRightXY
     * @return
     */
    public List<DenseMatrix> calculateEpipolarProjectionFor7Points(
        PairIntArray pointsLeftXY, PairIntArray pointsRightXY) {

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

        if (pointsLeftXY.getN() != 7) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the 7-point problem requires 7 points."
                + " pointsLeftXY.n=" + pointsLeftXY.getN());
        }

        return calculateEpipolarProjectionFor7Points(
            rewriteInto3ColumnMatrix(pointsLeftXY),
            rewriteInto3ColumnMatrix(pointsRightXY));
    }

    /**
     * calculate the epipolar projection for a set of matched points that are
     * at 7 points in length.
     *
     * @param theLeftXY
     * @param theRightXY
     * @return
     */
    @SuppressWarnings({"unchecked"})
    public List<DenseMatrix> calculateEpipolarProjectionFor7Points(
        DenseMatrix theLeftXY, DenseMatrix theRightXY) {

        if (theLeftXY == null) {
            throw new IllegalArgumentException("refactorLeftXY cannot be null");
        }
        if (theRightXY == null) {
            throw new IllegalArgumentException("refactorRightXY cannot be null");
        }
        if (theLeftXY.numRows() != theRightXY.numRows()) {
            throw new IllegalArgumentException(
                "theLeftXY and theRightXY must be same size");
        }
        if (theLeftXY.numColumns() != theRightXY.numColumns()) {
            throw new IllegalArgumentException(
                "theLeftXY and theRightXY must be same size");
        }

        if (theLeftXY.numColumns() != 7) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the 7-point problem requires 7 points."
                + " theLeftXY.n=" + theLeftXY.numColumns());
        }

        //x is xy[0], y is xy[1], xy[2] is all 1's
        NormalizedXY normalizedXY1 = normalize(theLeftXY);

        NormalizedXY normalizedXY2 = normalize(theRightXY);

        double[][] m = createFundamentalMatrix(
            normalizedXY1.getXy(), normalizedXY2.getXy());

        DenseMatrix aMatrix = new DenseMatrix(m);
        SVD svd = null;

        try {
            svd = SVD.factorize(aMatrix);
        } catch (Throwable t) {
            System.err.println(t.getMessage());
            return null;
        }

        // nCols = 9
        DenseMatrix nullSpace = algorithms.imageProcessing.util.MatrixUtil.nullSpace(svd);


        double[][] ff1 = new double[3][3];
        double[][] ff2 = new double[3][3];
        for (int i = 0; i < 3; i++) {

            ff1[i] = new double[3];
            ff1[i][0] = nullSpace.get((i * 3) + 0, 0);
            ff1[i][1] = nullSpace.get((i * 3) + 1, 0);
            ff1[i][2] = nullSpace.get((i * 3) + 2, 0);

            ff2[i] = new double[3];
            ff2[i][0] = nullSpace.get((i * 3) + 0, 1);
            ff2[i][1] = nullSpace.get((i * 3) + 1, 1);
            ff2[i][2] = nullSpace.get((i * 3) + 2, 1);
        }

        DenseMatrix[] solutions = solveFor7Point(ff1, ff2);

        //denormalize:  F = (T_1)^T * F * T_2
        //    T_1 is normalizedXY1.getNormalizationMatrix();
        //    T2 is normalizedXY2.getNormalizationMatrix();

        List<DenseMatrix> denormalizedSolutions = new ArrayList<DenseMatrix>();

        DenseMatrix t1Transpose = (DenseMatrix) normalizedXY1
            .getNormalizationMatrix().transpose();
        DenseMatrix t2 = normalizedXY2.getNormalizationMatrix();

        for (DenseMatrix solution : solutions) {

            if (solution == null) {
                continue;
            }
            
            DenseMatrix validated = validateSolution(solution, 
                normalizedXY1.getXy(), normalizedXY2.getXy());

            DenseMatrix denormFundamentalMatrix =
                MatrixUtil.multiply(t1Transpose,
                    MatrixUtil.multiply(solution, t2));

            double s = 1./(denormFundamentalMatrix.get(2, 2) + eps);
            MatrixUtil.multiply(denormFundamentalMatrix, s);

            denormFundamentalMatrix = (DenseMatrix) denormFundamentalMatrix.transpose();

            
            

            if (validated != null) {
                denormalizedSolutions.add(validated);
            }
            
            //denormalizedSolutions.add(denormFundamentalMatrix);
        }

        return denormalizedSolutions;
    }

    /**
     * The validation of the 7-point algorithm follows source code adapted
     * from this site and license:
     *
     * based upon code within  www.robots.ox.ac.uk/~vgg/hzbook/code/
        MIT License
        License for
        "MATLAB Functions for Multiple View Geometry"

        Copyright (c) 1995-2005 Visual Geometry Group
        Department of Engineering Science
        University of Oxford
        http://www.robots.ox.ac.uk/~vgg/
        Permission is hereby granted, free of charge, to any person obtaining a
        * copy of this software and associated documentation files
        * (the "Software"), to deal in the Software without restriction,
        * including without limitation the rights to use, copy, modify, merge,
        * publish, distribute, sublicense, and/or sell copies of the Software,
        * and to permit persons to whom the Software is furnished to do so,
        * subject to the following conditions:

        The above copyright notice and this permission notice shall be included
        * in all copies or substantial portions of the Software.

        The software is provided "as is", without warranty of any kind, express
        * or implied, including but not limited to the warranties of
        * merchantability, fitness for a particular purpose and noninfringement.
        * In no event shall the authors or copyright holders be liable for any
        * claim, damages or other liability, whether in an action of contract,
        * tort or otherwise, arising from, out of or in connection with the
        * software or the use or other dealings in the software.

       vgg_multiview/vgg_F_from_7pts_2img.m

       The method "signs_OK" validates the solution matrices:

        for i = 1:length(a)
          Fi = a(i)*FF{1} + (1-a(i))*FF{2};
          %for n = 1:7, disp(norm(x(:,n,1)'*Fi*x(:,n,2))), end  % test code
          if signs_OK(Fi,x1,x2)
            F = cat(3, F, Fi);
          end
        end

        return

        %%%%%%%%%%%%%%%%%%%%%%%%%

        % Checks sign consistence of F and x
        function OK = signs_OK(F,x1,x2)
        [u,s,v] = svd(F');
        e1 = v(:,3);
        l1 = vgg_contreps(e1)*x1;
        s = sum( (F*x2) .* l1 );
        OK = all(s>0) | all(s less than 0);
        return

    More on the subject is present in "Cheirality in Epipolar Geometry" by
    Werner & Pajdla, 2000 regarding realizability of two images.
    http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.32.9013&rep=rep1&type=pdf

    Very clear paper on cheirality:
    * http://users.cecs.anu.edu.au/~hartley/Papers/cheiral/revision/cheiral.pdf
    * The cheirality of a point is whether it lies in front of or behind a given
    * camera.  It's used to  distinguish between four different possible scene
    * reconstructions from two views.
    * A transform is cheirality-reversing for a given point if it swaps the
    * point from the front to the back of the camera, or vice-versa.
    * Otherwise it is called cheirality-preserving.
    */
    @SuppressWarnings({"unchecked"})
    private DenseMatrix validateSolution(DenseMatrix solution, DenseMatrix leftXY,
        DenseMatrix rightXY) {
        
        /*        
        NOTE: vgg_contreps of a 3X1 vector e1 is
            Y = [0      e1(3)  -e1(2)
                -e1(3)  0      e1(1)
                 e1(2) -e1(1)  0];
        NOTE: '.*' is matlab notation to operate on each field
        NOTE: the ' is mathematica syntax for conjugate transpose, a.k.a.
               the Hermitian. it's a matrix with signs reversed for imaginary
                components of complex numbers and then the matrix transposed.
                
        function OK = signs_OK(F,x1,x2)
        [u,s,v] = svd(F');
        e1 = v(:,3);
        l1 = vgg_contreps(e1)*x1;
        s = sum( (F*x2) .* l1 );
        OK = all(s>0) | all(s<0);

        (F*x2) .* l1 ==>  (solution * rightXY) .* (testE1 * leftXY)
        
        'sum' is a matlab function to sum for each column

        'all' is a function that returns '1' is all items are non-zero, else
            returns 0
        */
        DenseMatrix solutionHermitian = MatrixUtil.transpose(solution);
        double[][] leftRightEpipoles = calculateEpipoles(solutionHermitian);

        // 3 columns (x,y,1):
        double[] testE1 = leftRightEpipoles[0];
        //vgg_contreps of a 3X1 vector e1 is
        //    Y = [0      e1(3)  -e1(2)
        //        -e1(3)  0      e1(1)
        //         e1(2) -e1(1)  0];
        //  [x y ] [x0 x1 x2 x3 x4 x5 x6]
        //         [y0 y1 y2 y3 y4 y5 y6]
        double[][] contrepsE1 = new double[3][3];
        contrepsE1[0] = new double[]{0, testE1[2], -testE1[1]};
        contrepsE1[1] = new double[]{-testE1[2], 0, testE1[0]};
        contrepsE1[2] = new double[]{testE1[1], -testE1[0], 0};
        
        //l1 = e1contreps*x1;
        
        // 3 X 7
        double[][] l1 = MatrixUtil.multiply(contrepsE1, MatrixUtil.convertToRowMajor(leftXY));

        // s = sum( (F*x2) .* l1 );
        // .* is matlab notation 
        // https://www.mathworks.com/help/matlab/ref/times.html
        
        // 3 X 7
        DenseMatrix fx2 = MatrixUtil.multiply(solution, rightXY);
        
        // 3 x 7
        DenseMatrix fx2l1 = fx2.copy();
        for (int row = 0; row < testE1.length; ++row){
            for (int col = 0; col < fx2.numColumns(); ++col) {
                double value = fx2l1.get(row, col) * l1[row][col];
                fx2l1.set(row, col, value);
            }
        }
        
        // matlab sum function:
        //    If A is a matrix, then sum(A) returns a row vector containing 
        //    the sum of each column.
        double[] sum = new double[fx2l1.numColumns()];
        for (int row = 0; row < fx2l1.numRows(); ++row){
            for (int col = 0; col < fx2l1.numColumns(); ++col) {
                double value = fx2l1.get(row, col);
                sum[col] += value;
            }
        }

        int nGTZ = 0;
        int nLTZ = 0;
        for (int i = 0; i < sum.length; ++i) {
            if (sum[i] == 0) {
                return null;
            } else if (sum[i] > 0) {
                nGTZ++;
            } else {
                nLTZ++;
            }
        }

        if ((nGTZ > 0 && nLTZ == 0) || (nLTZ > 0 && nGTZ == 0)) {
            return solution;
        }
        return null;
    }

    DenseMatrix[] solveFor7Point(double[][] ff1, double[][] ff2) {

        //solve for the roots of equation a0 * x^3 + a1 * x^2 + a2 * x + a3 = 0;
        
        double a0 = calculateCubicRoot3rdOrderCoefficientFor7Point(ff1, ff2);
        double a1 = calculateCubicRoot2ndOrderCoefficientFor7Point(ff1, ff2);
        double a2 = calculateCubicRoot1stOrderCoefficientFor7Point(ff1, ff2);
        double a3 = calculateCubicRoot0thOrderCoefficientFor7Point(ff1, ff2);

        double[] roots = MiscMath.solveCubicRoots(a0, a1, a2, a3);

        double[][] m = new double[3][];
        for (int i = 0; i < 3; i++) {
            m[i] = new double[3];
        }

        DenseMatrix[] solutions = new DenseMatrix[roots.length];

        for (int i = 0; i < roots.length; i++) {

            //Fi = a(i)*FF{1} + (1-a(i))*FF{2};

            double a = roots[i];

            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    m[row][col] = a*ff1[row][col] + (1. - a)*ff2[row][col];
                }
            }

            solutions[i] = new DenseMatrix(m);
        }

        return solutions;
    }

    private double calculateCubicRoot3rdOrderCoefficientFor7Point(
        double[][] ff1, double[][] ff2) {

        double b = ff1[0][0];
        double e = ff1[1][0];
        double h = ff1[2][0];
        double c = ff1[0][1];
        double f = ff1[1][1];
        double i = ff1[2][1];
        double d = ff1[0][2];
        double g = ff1[1][2];
        double j = ff1[2][2];

        double k = ff2[0][0];
        double n = ff2[1][0];
        double q = ff2[2][0];
        double l = ff2[0][1];
        double o = ff2[1][1];
        double r = ff2[2][1];
        double m = ff2[0][2];
        double p = ff2[1][2];
        double s = ff2[2][2];

        double sum = h*g*c + h*o*d + h*o*m + h*p*l + i*e*d + i*g*k
            + i*n*m + i*p*b + j*e*l + j*k*o + j*n*c + q*f*d
            + q*f*m + q*g*l + q*p*c + r*e*m + r*g*b + r*n*d
            + r*p*k + s*b*o + s*e*c + s*k*f + s*n*l + b*j*f
            - h*f*d - h*f*m - h*g*l - h*p*c - i*e*m - i*g*b
            - i*n*d - i*p*k - j*b*o - j*e*c - j*k*f - j*n*l
            - q*g*c - q*o*d - q*o*m - q*p*l - r*e*d - r*g*k
            - r*n*m - r*p*b - s*b*f - s*e*l - s*k*o - s*n*c;

        return sum;
    }

    private double calculateCubicRoot2ndOrderCoefficientFor7Point(double[][] ff1,
        double[][] ff2) {

        double b = ff1[0][0];
        double e = ff1[1][0];
        double h = ff1[2][0];
        double c = ff1[0][1];
        double f = ff1[1][1];
        double i = ff1[2][1];
        double d = ff1[0][2];
        double g = ff1[1][2];
        double j = ff1[2][2];

        double k = ff2[0][0];
        double n = ff2[1][0];
        double q = ff2[2][0];
        double l = ff2[0][1];
        double o = ff2[1][1];
        double r = ff2[2][1];
        double m = ff2[0][2];
        double p = ff2[1][2];
        double s = ff2[2][2];

        double sum = h*f*m + h*g*l + h*p*c + i*e*m + i*n*d + i*p*k
            + i*p*k + j*b*o + j*k*f + j*n*l + j*n*l + q*g*c
            + q*o*d + q*o*d + q*o*m + q*o*m + q*o*m + q*p*l
            + q*p*l + q*p*l + r*e*d + r*g*k + r*g*k + r*n*m
            + r*n*m + r*n*m + r*p*b + r*p*b + s*o*k + s*b*f
            + s*e*l + s*e*l + s*k*o + s*k*o + s*n*c + s*n*c
            - h*o*d - h*o*m - h*o*m - h*p*l - h*p*l - i*g*k
            - i*n*m - i*n*m - i*p*b - j*e*l - j*k*o - j*n*c
            - j*o*k - q*f*d - q*f*m - q*f*m - q*g*l - q*g*l
            - q*p*c - q*p*c - r*e*m - r*e*m - r*g*b - r*n*d
            - r*n*d - r*p*k - r*p*k - r*p*k - s*b*o - s*b*o
            - s*e*c - s*k*f - s*k*f - s*n*l - s*n*l - s*n*l;

        return sum;
    }

    private double calculateCubicRoot1stOrderCoefficientFor7Point(double[][] ff1,
        double[][] ff2) {

        /*
        f1 =
         b c d
         e f g
         h i j
        f2 =
         k l m
         n o p
         q r s
        */

        double b = ff1[0][0];
        double e = ff1[1][0];
        double h = ff1[2][0];
        double c = ff1[0][1];
        double f = ff1[1][1];
        double i = ff1[2][1];
        double d = ff1[0][2];
        double g = ff1[1][2];
        double j = ff1[2][2];

        double k = ff2[0][0];
        double n = ff2[1][0];
        double q = ff2[2][0];
        double l = ff2[0][1];
        double o = ff2[1][1];
        double r = ff2[2][1];
        double m = ff2[0][2];
        double p = ff2[1][2];
        double s = ff2[2][2];

        double sum = h*o*m + h*p*l + i*n*m + j*o*k + q*f*m + q*g*l + q*p*c
            + r*e*m + r*n*d + r*p*k + r*p*k + r*p*k + s*b*o + s*k*f
            + s*n*l + s*n*l + s*n*l
            - i*p*k - j*n*l - q*o*d - q*o*m - q*o*m - q*o*m
            - q*p*l - q*p*l - q*p*l - r*g*k - r*n*m - r*n*m
            - r*n*m - r*p*b - s*e*l - s*k*o - s*n*c - s*o*k
            - s*o*k;

        return sum;
    }

    private double calculateCubicRoot0thOrderCoefficientFor7Point(
        double[][] ff1, double[][] ff2) {

        /*
        f1 =
         b c d
         e f g
         h i j
        f2 =
         k l m
         n o p
         q r s
        */

        double k = ff2[0][0];
        double n = ff2[1][0];
        double q = ff2[2][0];
        double l = ff2[0][1];
        double o = ff2[1][1];
        double r = ff2[2][1];
        double m = ff2[0][2];
        double p = ff2[1][2];
        double s = ff2[2][2];

        double sum = q * o * m + q * p * l + r * n * m + s * o * k - r*p*k
            - s*n*l;

        return sum;
    }

    @SuppressWarnings({"unchecked"})
    DenseMatrix calculateFundamentalMatrix(NormalizedXY normalizedXY1,
        NormalizedXY normalizedXY2) {

        //build the fundamental matrix
        double[][] m = createFundamentalMatrix(normalizedXY1.getXy(),
            normalizedXY2.getXy());

        /*
        compute linear least square solution:
            solve A = U * D * V^T   for A*f = [..x...]*f = 0
            A has rank 8.  f has rank 2.

        calculate [U,D,V] from svd(A):
        */
        DenseMatrix aMatrix = new DenseMatrix(m);

        SVD svd = null;
        try {
            svd = SVD.factorize(aMatrix);
        } catch (NotConvergedException e) {
            log.severe(e.getMessage());
            return null;
        }

        //A is nData rows X 3 columns

        DenseMatrix V = (DenseMatrix) svd.getVt().transpose();

        // creates U as nXY1 x nXY1 matrix  (M X M)
        //         D as length 3 array      (vector of len N)
        //         V as 9 x 9 matrix        (N*N X N*N)

        // mRows = 9; nCols = 9

        // reshape V to 3x3  (just the last column)

        int vNCols = V.numColumns();

        double[][] ff = new double[3][3];
        for (int i = 0; i < 3; i++) {
            ff[i] = new double[3];
            ff[i][0] = V.get((i * 3) + 0, vNCols - 1);
            ff[i][1] = V.get((i * 3) + 1, vNCols - 1);
            ff[i][2] = V.get((i * 3) + 2, vNCols - 1);
        }
        DenseMatrix fMatrix = new DenseMatrix(ff);

        /* make the fundamental matrix have a rank of 2
        by performing a svd and then reconstructing with the two largest
        singular values.
            [U,D,V] = svd(F,0);

        From [U,D,V] we create:
            F = U * diag([D(1,1) D(2,2) 0]) * V^T, where V^T is V transposed.
        */
        try {
            svd = SVD.factorize(fMatrix);
        } catch (NotConvergedException e) {
            log.severe(e.getMessage());
            return null;
        }

        // creates U as 3 x 3 matrix
        //         D as length 3 array
        //         V as 3 x 3 matrix

        double[] sDiag = svd.getS();

        //F = U * diag([D(1,1) D(2,2) 0]) * V^T, where V^T is V transposed.

        // keep the largest 2 values in sDiag to make the diagonal rank 2
        DenseMatrix d = new DenseMatrix(3, 3);
        if (sDiag.length > 0) {
            d.set(0, 0, sDiag[0]);
        }
        if (sDiag.length > 1) {
            d.set(1, 1, sDiag[1]);
        }

        V = svd.getVt();

        /*
        multiply the terms:
             F = dot(U, dot(diag(D),V^T))
        */
        DenseMatrix dDotV = MatrixUtil.multiply(d, V);

        // 3x3
        DenseMatrix theFundamentalMatrix = MatrixUtil.multiply(svd.getU(), dDotV);

        //System.out.println("fm=" + theFundamentalMatrix.toString());

        DenseMatrix denormFundamentalMatrix =
            denormalizeTheFundamentalMatrix(theFundamentalMatrix,
                normalizedXY1, normalizedXY2);

        return denormFundamentalMatrix;
    }

    @SuppressWarnings({"unchecked"})
    DenseMatrix denormalizeTheFundamentalMatrix(
        DenseMatrix normalizedFundamentalMatrix,
        NormalizedXY normalizedLeftXY, NormalizedXY normalizedRightXY) {

        /*
        denormalize
            F = (T_1)^T * F * T_2
            where T_1 is normalizedXY1.getNormalizationMatrix();
            and T2 is normalizedXY2.getNormalizationMatrix();
        */

        DenseMatrix t1Transpose = (DenseMatrix) normalizedLeftXY
            .getNormalizationMatrix().transpose();
        DenseMatrix t2 = normalizedRightXY.getNormalizationMatrix();

        DenseMatrix denormFundamentalMatrix =
            MatrixUtil.multiply(t1Transpose,
                MatrixUtil.multiply(normalizedFundamentalMatrix, t2));

        double factor = 1./(denormFundamentalMatrix.get(2, 2) + eps);
        MatrixUtil.multiply(denormFundamentalMatrix, factor);

        return denormFundamentalMatrix;
    }

    /**
     normalize the x,y coordinates as recommended by Hartley 1997 and return
     the matrix and coordinates.
     does not modify the state of this transformer instance.

     the normalized coordinates have an origin = centroid of
     the points.
     the mean square distance of the normalized points from the
     origin is sqrt(2) pixels.


     * @param xyPair
     * @return
     */
    @SuppressWarnings({"unchecked"})
    NormalizedXY normalize(DenseMatrix xy) {

        /*
        format points such that the applied translation
        and scaling have the effect of:

        a) points are translated so that their centroid is at the origin.
        b) points are then scaled so that the average distance from the
           origin is sqrt(2)
        c) the transformation is applied to each of the 2 images separately.
        */

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        int n = xy.numColumns();

        //x is xy[0], y is xy[1], xy[2] is all 1's
        double cen0 = 0;
        double cen1 = 0;
        for (int i = 0; i < n; ++i) {
            cen0 += xy.get(0, i);
            cen1 += xy.get(1, i);
        }
        cen0 /= (double)n;
        cen1 /= (double)n;

        double mean = 0;

        for (int i = 0; i < n; i++) {
            double diffX = xy.get(0, i) - cen0;
            double diffY = xy.get(1, i) - cen1;
            double dist = Math.sqrt((diffX * diffX) + (diffY * diffY));
            mean += dist;
        }

        mean /= (double)n;

        /*
        mean * factor = sqrt(2)
        */
        double scaleFactor = Math.sqrt(2)/mean;

        DenseMatrix tMatrix = createScaleTranslationMatrix(scaleFactor, cen0, cen1);

        DenseMatrix normXY = new DenseMatrix(MatrixUtil.dot(tMatrix, xy));

        NormalizedXY normalizedXY = new NormalizedXY();
        normalizedXY.setCentroidXY(new double[]{cen0, cen1});
        normalizedXY.setNormMatrix(tMatrix);
        normalizedXY.setXy(normXY);

        return normalizedXY;
    }

    /**
     * create a matrix to be applied on the left side of the dot operator
     * with a matrix of points to transform the points by scale and translation.
     * @param scale
     * @param centroidX
     * @param centroidY
     * @return
     */
    protected DenseMatrix createScaleTranslationMatrix(double scale,
        double centroidX, double centroidY) {

        /*
        scale, rotate, then translate.
        let xc = x centroid
        let yc = y centroid
        let sX = scaleFactor in x
        let sY = scaleFactor in y
            here, sX = sY
        let r = rotation
            here, rotation = 0, so math.cos(0)=1, math.sin(0)=0
        let tX = x translation
            here, tX = -xc*sX
        let tY = y translation
            here, tY = -yc*sY

        transformation equations for x are for translating the points to the
        centroid, then scaling them
           x_transformed = xc*s + ((x - xc)*s) + tX = x*s - xc*s
           y_transformed = yc*s + ((y - yc)*s) + tY = y*s - yc*s

        matrix for xy matrix has format xy[row][col]
            where row=0 is x and row=1 is y, and row=2 is placeholder of value 1
        x[0]  x[1]  x[2]  x[3]
        y[0]  y[1]  y[2]  y[3]
        1      1      1     1

        Formatting the scale and translation into a matrix that can be used
        with dot operator to transform the points.
        Because the xy points have the x and y along rows, this new transformation
        matrix must be used on the left side of the dot operation.

           t dot xy = tranformed xy

        t00     t01      t02
        t10     t11      t12
        0        0        1

        t00*x[0] + t01*y[0] + t02*1 = x[0]*s - xc*s
                         0
             => t00 = s
             => t01 = 0
             => t02 = -xc*s

        t10*x[0] + t11*y[0] + t12*1 = y[0]*s - yc*s
            0
             => t10 = 0
             => t11 = s
             => t12 = -yc*s
        */

        double[][] t = new double[3][];
        t[0] = new double[]{scale,       0,     -centroidX*scale};
        t[1] = new double[]{0,           scale, -centroidY*scale};
        t[2] = new double[]{0,           0,           1};
        DenseMatrix tMatrix = new DenseMatrix(t);

        return tMatrix;
    }

    /**
     * write a matrix of size mRows = 3, nCols = xyPairs.getN()
     * @param xyPairs
     * @return
     */
    public DenseMatrix rewriteInto3ColumnMatrix(PairFloatArray xyPairs) {

        DenseMatrix xy = new DenseMatrix(3, xyPairs.getN());

        // rewrite xyPairs into a matrix of size 3 X xy.getN();
        // row 0 is x
        // row 1 is y
        // row 2 is all 1's
        for (int i = 0; i < xyPairs.getN(); i++) {
            xy.set(0, i, xyPairs.getX(i));
            xy.set(1, i, xyPairs.getY(i));
            xy.set(2, i, 1);
        }

        return xy;
    }

    /**
     * write a matrix of size mRows = 3, nCols = xyPairs.getN()
     * @param xyPairs
     * @return
     */
    public DenseMatrix rewriteInto3ColumnMatrix(List<PairInt> xyPairs) {

        DenseMatrix xy = new DenseMatrix(3, xyPairs.size());

        // rewrite xyPairs into a matrix of size 3 X xy.getN();
        // row 0 is x
        // row 1 is y
        // row 2 is all 1's
        for (int i = 0; i < xyPairs.size(); i++) {
            PairInt p = xyPairs.get(i);
            xy.set(0, i, p.getX());
            xy.set(1, i, p.getY());
            xy.set(2, i, 1);
        }

        return xy;
    }

    /**
     * write a matrix of size mRows = 3, nCols = xyPairs.getN()
     * @param xyPairs
     * @return
     */
    public DenseMatrix rewriteFirstItemInto3ColumnMatrix(List<List<PairInt>> xyPairs) {

        DenseMatrix xy = new DenseMatrix(3, xyPairs.size());

        // rewrite xyPairs into a matrix of size 3 X xy.getN();
        // row 0 is x
        // row 1 is y
        // row 2 is all 1's
        for (int i = 0; i < xyPairs.size(); i++) {
            List<PairInt> points = xyPairs.get(i);
            PairInt p = points.get(0);
            xy.set(0, i, p.getX());
            xy.set(1, i, p.getY());
            xy.set(2, i, 1);
        }

        return xy;
    }

    /**
     * write a matrix of size mRows = 3, nCols = xyPairs.getN()
     * @param xyPairs
     * @return
     */
    public DenseMatrix rewriteInto3ColumnMatrix(PairIntArray xyPairs) {

        DenseMatrix xy = new DenseMatrix(3, xyPairs.getN());

        // rewrite xyPairs into a matrix of size 3 X xy.getN();
        // row 0 is x
        // row 1 is y
        // row 2 is all 1's
        for (int i = 0; i < xyPairs.getN(); i++) {
            xy.set(0, i, xyPairs.getX(i));
            xy.set(1, i, xyPairs.getY(i));
            xy.set(2, i, 1);
        }

        return xy;
    }

    /**
     * @param normXY1 a matrix of size 3 x nPoints, where 1st column is x,
     * second is y.
     * @param normXY2 a matrix of size 3 x nPoints, where 1st column is x,
     * second is y.
     * @return
     */
    double[][] createFundamentalMatrix(DenseMatrix normXY1,
        DenseMatrix normXY2) {

        if (normXY1 == null) {
            throw new IllegalArgumentException("normXY1 cannot be null");
        }
        if (normXY2 == null) {
            throw new IllegalArgumentException("normXY2 cannot be null");
        }
        if (normXY1.numColumns() != normXY2.numColumns()) {
            throw new IllegalArgumentException(
            "the number of columns in normXY1 != number of cols in normXY2");
        }

        int nXY1 = normXY1.numColumns();

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
    @SuppressWarnings({"unchecked"})
    double[][] calculateEpipoles(DenseMatrix fundamentalMatrix) {

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

        epipoles:
             [U,D,V] = svd(denormalized FundamentalMatrix);
             e1 = last column of V divided by it's last item
             e2 = last column of U divided by it's last item

        */
        SVD svdE;
        try {
            svdE = SVD.factorize(fundamentalMatrix);
        } catch (NotConvergedException ex) {
            Logger.getLogger(EpipolarTransformer.class.getName())
                .log(Level.SEVERE, null, ex);
            return null;
        }
        DenseMatrix V = (DenseMatrix) svdE.getVt().transpose();
        double[] e1 = new double[V.numColumns()];
        double e1Div = V.get(2, 2);
        for (int i = 0; i < e1.length; i++) {
            e1[i] = V.get(i, 2)/e1Div;
        }
        DenseMatrix U = svdE.getU();
        double[] e2 = new double[U.numColumns()];
        double e2Div = U.get(2, 2);
        for (int i = 0; i < e2.length; i++) {
            e2[i] = U.get(i, 2)/e2Div;
        }

        double[][] e = new double[2][];
        e[0] = e1;
        e[1] = e2;

        return e;
    }

    private EpipolarFeatureTransformationFit combineErrors(EpipolarTransformationFit
        distanceErrors, EpipolarFeatureTransformationFit featureErrors) {

        /*
        in order to have the distance errors count as much as the SSD errors,
        need to scale them up or SSD down by a factor.

        will use the descriptor size and the average of the maximum SSD that
        could be calculated and the minimum and make a factor for the
        distances of SSDFactor/tolerance.
        */

        // sum square diffs / size = (d0*d0) + (d3*d3).../n
        double maxSSD = 255. * 255.;

        double distScaleFactor = (maxSSD/2.)/distanceErrors.getTolerance();

        Map<Integer, Double> indexSSDErrorsMap = new HashMap<Integer, Double>();
        Map<Integer, FeatureComparisonStat> indexFeatureMap =
            new HashMap<Integer, FeatureComparisonStat>();
        for (int i = 0; i < featureErrors.getInlierIndexes().size(); ++i) {
            Integer index = featureErrors.getInlierIndexes().get(i);
            indexSSDErrorsMap.put(index, featureErrors.getErrors().get(i));
            indexFeatureMap.put(index,
                featureErrors.getFeatureComparisonStats().get(i));
        }

        List<Integer> outputInliers = new ArrayList<Integer>();
        List<Double> outputDistances = new ArrayList<Double>();
        List<FeatureComparisonStat> fcs = new ArrayList<FeatureComparisonStat>();

        for (int i = 0; i < distanceErrors.getInlierIndexes().size(); ++i) {
            Integer index = distanceErrors.getInlierIndexes().get(i);
            Double ssd = indexSSDErrorsMap.get(index);
            if (ssd != null) {

                outputInliers.add(index);

                Double dist = distanceErrors.getErrors().get(i);
                double cost = dist.doubleValue() * ssd.doubleValue();
                outputDistances.add(Double.valueOf(cost));

                fcs.add(indexFeatureMap.get(index));
            }
        }
        double costTerm2 = 1./(double)outputDistances.size();
        for (int i = 0; i < outputDistances.size(); ++i) {
            double err = outputDistances.get(i).doubleValue() * costTerm2 * costTerm2
                * distScaleFactor;
            outputDistances.set(i, Double.valueOf(err));
        }

        EpipolarFeatureTransformationFit fit =
            new EpipolarFeatureTransformationFit(
            distanceErrors.getFundamentalMatrix(),
            outputInliers, fcs,
            distanceErrors.getErrorType(), outputDistances,
            distanceErrors.getTolerance());

        return fit;
    }

    /**
     * extract a single row from the 3 rows present,
     * transpose, then replicate that single column.
     * the result should be m.columns X 3
     * @param m
     * @param row
     * @return matrix of size nData X 3 (== m.numColumns X 3)
     */
    private DenseMatrix exRowTRepl(DenseMatrix m, int row) {

         //repmat(X2(3,:)',1,3)
        //   extract row 2 for all columns (== 1 X nData)
        //   and transpose that (== nData X 1)
        //   then replicate that column twice more

        DenseMatrix repl = exRowRepl(m, row);

        return MatrixUtil.transpose(repl);
    }
    
    /**
     * extract a single row from the 3 rows present,
     * then replicate that single column.
     * the result should be 3 X m.columns
     * @param m
     * @param row
     * @return matrix of size 3 X nData (== 3 X m.numColumns)
     */
    private DenseMatrix exRowRepl(DenseMatrix m, int row) {

         //repmat(X2(3,:),1,3)
        //   extract row 2 for all columns (== 1 X nData)
        //   then replicate that row twice more

        int nCols = m.numColumns();

        DenseMatrix out = new DenseMatrix(3, nCols);

        for (int col = 0; col < m.numColumns(); ++col) {
            double v = m.get(row, col);
            for (int row2 = 0; row2 < 3; ++row2) {
                out.set(row2, col, v);
            }
        }

        return out;
    }

    private double[] sumMult(DenseMatrix m1, DenseMatrix m2) {

        int n = m1.numColumns();
        
        double[] out = new double[n];
        for (int i = 0; i < n; ++i) {
            for (int row = 0; row < m1.numRows(); ++row) {
                out[i] += (m1.get(row, i) * m2.get(row, i));
            }
        }
        
        return out;
    }

    public static class NormalizedXY {

        /**
         * 3 dimensional matrix, with column 0 being x, column 1 being y,
         * and the last column is place holder 1's
         */
        private DenseMatrix xy = null;

        private double[] centroidXY = null;

        private DenseMatrix normalizationMatrix = null;

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
        public DenseMatrix getNormalizationMatrix() {
            return normalizationMatrix;
        }

        /**
         * @param normMatrix holding the scale and offsets to apply to x, y
         */
        public void setNormMatrix(DenseMatrix normMatrix) {
            this.normalizationMatrix = normMatrix;
        }

        /**
         * @return the xy
         */
        public DenseMatrix getXy() {
            return xy;
        }

        /**
         * @param xy the xy to set
         */
        public void setXy(DenseMatrix xy) {
            this.xy = xy;
        }
    }

    public PairIntArray getEpipolarLine(DenseMatrix epipolarLines, int imgWidth,
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

    /**
     * evaluate fit for already matched point lists
     * @param fm
     * @param leftPoints
     * @param rightPoints
     * @param tolerance
     * @return
     */
    public EpipolarTransformationFit calculateEpipolarDistanceError(
        DenseMatrix fm, DenseMatrix leftPoints, DenseMatrix rightPoints,
        double tolerance) {

        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        if (leftPoints == null) {
            throw new IllegalArgumentException("leftPoints cannot be null");
        }
        if (rightPoints == null) {
            throw new IllegalArgumentException("rightPoints cannot be null");
        }
        int nRows = leftPoints.numRows();
        if (nRows != rightPoints.numRows()) {
            throw new IllegalArgumentException("matrices must have same number of rows");
        }

        //2D point (x,y) and line (a, b, c): dist=(a*x + b*y + c)/sqrt(a^2 + b^2)

        PairFloatArray distances = calculateDistancesFromEpipolar(fm,
            leftPoints, rightPoints);

        List<Double> errors = new ArrayList<Double>();

        List<Integer> inlierIndexes = new ArrayList<Integer>();

        for (int i = 0; i < distances.getN(); ++i) {

            float leftPtD = distances.getX(i);

            float rightPtD = distances.getY(i);

            float dist = (float)Math.sqrt(leftPtD*leftPtD + rightPtD*rightPtD);

            if (dist > tolerance) {
                continue;
            }

            inlierIndexes.add(Integer.valueOf(i));

            errors.add(Double.valueOf(dist));
        }

        EpipolarTransformationFit fit = null;

        if (errors.size() > 0) {
            fit = new EpipolarTransformationFit(fm, inlierIndexes,
                ErrorType.DIST_TO_EPIPOLAR_LINE, errors, tolerance);
        } else {
            fit = new EpipolarTransformationFit(fm, new ArrayList<Integer>(),
                ErrorType.DIST_TO_EPIPOLAR_LINE, new ArrayList<Double>(), tolerance);
        }

        fit.setNMaxMatchable(leftPoints.numColumns());

        return fit;
    }

    /**
     * evaluate fit for already matched point lists
     * @param fm
     * @param leftPoints
     * @param rightPoints
     * @param tolerance
     * @return
     */
    public EpipolarTransformationFit calculateEpipolarDistanceErrorThenFilter(
        DenseMatrix fm, DenseMatrix leftPoints, DenseMatrix rightPoints,
        double tolerance) {

        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        if (leftPoints == null) {
            throw new IllegalArgumentException("leftPoints cannot be null");
        }
        if (rightPoints == null) {
            throw new IllegalArgumentException("rightPoints cannot be null");
        }
        int nRows = leftPoints.numRows();
        if (nRows != rightPoints.numRows()) {
            throw new IllegalArgumentException("matrices must have same number of rows");
        }

        //2D point (x,y) and line (a, b, c): dist=(a*x + b*y + c)/sqrt(a^2 + b^2)

        PairFloatArray distances = calculateDistancesFromEpipolar(fm,
            leftPoints, rightPoints);

        List<Double> errors = new ArrayList<Double>();

        List<Integer> inlierIndexes = new ArrayList<Integer>();

        for (int i = 0; i < distances.getN(); ++i) {

            float leftPtD = distances.getX(i);

            float rightPtD = distances.getY(i);

            float dist = (float)Math.sqrt(leftPtD*leftPtD + rightPtD*rightPtD);

            if (dist > tolerance) {
                continue;
            }

            inlierIndexes.add(Integer.valueOf(i));

            errors.add(Double.valueOf(dist));
        }

        filterForDegenerate(leftPoints, inlierIndexes, errors);
        filterForDegenerate(rightPoints, inlierIndexes, errors);

        EpipolarTransformationFit fit = null;

        if (errors.size() > 0) {
            fit = new EpipolarTransformationFit(fm, inlierIndexes,
                ErrorType.DIST_TO_EPIPOLAR_LINE, errors, tolerance);
        } else {
            fit = new EpipolarTransformationFit(fm, new ArrayList<Integer>(),
                ErrorType.DIST_TO_EPIPOLAR_LINE, new ArrayList<Double>(), tolerance);
        }

        fit.setNMaxMatchable(leftPoints.numColumns());

        return fit;
    }

    /**
     * find the distance of the given points from their respective projective
     * epipolar lines.
     * @param fm
     * @param rightEpipolarLines
     * @param leftEpipolarLines
     * @param leftPoints
     * @param rightPoints
     * @param tolerance
     * @return
     */
    public PairFloatArray calculateDistancesFromEpipolar(
        DenseMatrix fm, DenseMatrix matchedLeftPoints,
        DenseMatrix matchedRightPoints) {

        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        if (matchedLeftPoints == null) {
            throw new IllegalArgumentException("matchedLeftPoints cannot be null");
        }
        if (matchedRightPoints == null) {
            throw new IllegalArgumentException("rightPoints cannot be null");
        }
        int nRows = matchedLeftPoints.numRows();
        if (nRows != matchedRightPoints.numRows()) {
            throw new IllegalArgumentException("matrices must have same number of rows");
        }

        /*
        u_2^T * F * u_1 = 0  where u are the x,y points in images _1 and _2
        u_1 = (x_1, y_1, 1)^T
        u_2 = (x_2, y_2, 1)^T
        */

        int n = matchedLeftPoints.numColumns();

        PairFloatArray distances = new PairFloatArray(n);

        DenseMatrix rightEpipolarLines = MatrixUtil.multiply(fm, matchedLeftPoints);

        DenseMatrix leftEpipolarLines = MatrixUtil.multiply(fm.transpose(),
            matchedRightPoints);

        float[] output = new float[2];

        for (int i = 0; i < matchedLeftPoints.numColumns(); i++) {

            calculatePerpDistFromLines(matchedLeftPoints,
                matchedRightPoints, 
                rightEpipolarLines, leftEpipolarLines,
                i, i, output);

            distances.add(output[0], output[1]);
        }

        return distances;
    }

    public void calculatePerpDistFromLines(DenseMatrix leftPoints,
        DenseMatrix rightPoints, 
        DenseMatrix epipolarLinesFromLeft,
        DenseMatrix epipolarLinesFromRight, 
        int leftIdx, int rightIdx,
        float[] output) {
        
        double a = epipolarLinesFromLeft.get(0, leftIdx);
        double b = epipolarLinesFromLeft.get(1, leftIdx);
        double c = epipolarLinesFromLeft.get(2, leftIdx);

        double aplusb = Math.sqrt((a*a) + (b*b));

        //dist = (a*x + b*y + c)/sqrt(a^2 + b^2)

        double x = rightPoints.get(0, rightIdx);
        double y = rightPoints.get(1, rightIdx);

        double d = (a*x + b*y + c)/aplusb;

        // find the reverse distance by projection:
        double aRev = epipolarLinesFromRight.get(0, rightIdx);
        double bRev = epipolarLinesFromRight.get(1, rightIdx);
        double cRev = epipolarLinesFromRight.get(2, rightIdx);

        double xL = leftPoints.get(0, leftIdx);
        double yL = leftPoints.get(1, leftIdx);

        double dRev = (aRev*xL + bRev*yL + cRev)/
            Math.sqrt((aRev*aRev + bRev*bRev));

        output[0] = (float)dRev;
        output[1] = (float)d;
    }

    public EpipolarTransformationFit calculateErrorThenFilter(DenseMatrix fm,
        DenseMatrix x1, DenseMatrix x2, ErrorType errorType, double tolerance) {

        if (errorType.equals(ErrorType.SAMPSONS)) {
            return calculateSampsonsErrorThenFilter(fm, x1, x2, tolerance);
        } else {
            return calculateEpipolarDistanceErrorThenFilter(fm, x1, x2, tolerance);
        }
    }

    public EpipolarTransformationFit calculateError(DenseMatrix fm,
        DenseMatrix x1, DenseMatrix x2, ErrorType errorType, double tolerance) {

        if (errorType.equals(ErrorType.SAMPSONS)) {
            return calculateSampsonsError(fm, x1, x2, tolerance);
        } else {
            return calculateEpipolarDistanceError(fm, x1, x2, tolerance);
        }
    }

    /**
     Return the "algebraic distance" needed for use in calculating
     Sampson's distance.
     The algebraic distance has no geometrical significance,
     (it isn't the perpendicular distance).
     * The topic is discussed in
     "The Development and Comparison of Robust Methods
      for Estimating the Fundamental Matrix" by Torr and Murray, 1997
      as Equation (2) on page 274.
      The return here contains 2 matrices which contain parts of equation
      2, formatted for use in calculating Sampson's distance.
     <pre>
     The method is adapted from code from the book
     "Multiple View Geometry in Computer Vision" by 
     Hartley and Zisserman, 2004.
       
     obtained from 
     http://www.robots.ox.ac.uk/~vgg/hzbook/

        license is MIT license

        Acknowledgements: These functions are written by: David Capel,
        Andrew Fitzgibbon, Peter Kovesi, Tomas Werner, Yoni Wexler,
        and Andrew Zisserman
     </pre>
     * @param fm 3 X 3 matrix
     * @param x1 3 X nData matrix
     * @param x2 3 X nData matrix
     * @return DX matrix of size nData X 3 is element 0
     *         DY matrix of size nData X 3 is element 1
     */
    private DenseMatrix[] algebraicDistance(DenseMatrix fm, DenseMatrix x1,
        DenseMatrix x2) {

        // N = size(X1,2);
        int n = x1.numColumns();

        //repmat(X2(3,:)',1,3)
        //   extract row 2 for all columns (== 1 X nData)
        //   and transpose that (== nData X 1)
        //   then replicate that column twice more
        //   to make nData X 3 matrix
        //X1' is nData X 3

        // x coord of x2
        DenseMatrix x2RowT0 = exRowTRepl(x2, 0);

        // y coord of x2
        DenseMatrix x2RowT1 = exRowTRepl(x2, 1);

        DenseMatrix x2RowT2 = exRowTRepl(x2, 2);

        DenseMatrix x1T = MatrixUtil.transpose(x1);

        //nData X 3
        // x1 times 1 = x1_x, x1_y, 1
        DenseMatrix dX0 = algorithms.imageProcessing.util.MatrixUtil.multiplyPointwise(
            x1T, x2RowT2);
        assert(dX0.numRows() == n);
        assert(dX0.numColumns() == 3);

        //nData X 3
        DenseMatrix dX1 = new DenseMatrix(x2.numColumns(), 3);
        assert(dX1.numColumns() == 3);
        assert(dX1.numRows() == n);
        
        //nData X 3
        // x1 times -x cooord of x2 = x1_x * x2_x, x1_y * x2_x, x2_x
        DenseMatrix dX2 = algorithms.imageProcessing.util.MatrixUtil.multiplyPointwise(
            x1T, x2RowT0);
        MatrixUtil.multiply(dX2, -1);
        assert(dX2.numColumns() == 3);
        assert(dX2.numRows() == n);
        
        //nData X 3
        DenseMatrix dY0 = new DenseMatrix(x2.numColumns(), 3);
        assert(dY0.numColumns() == 3);
        assert(dY0.numRows() == n);
        
        //nData X 3
        // x1 times 1 = x1_x, x1_y, 1
        DenseMatrix dY1 = algorithms.imageProcessing.util.MatrixUtil.multiplyPointwise(
            x1T, x2RowT2);
        assert(dY1.numColumns() == 3);
        assert(dY1.numRows() == n);
        
        //nData X 3
        // x1 times -y cooord x2 = of x1_x * y2_x, x1_y * y2_x, y2_x
        DenseMatrix dY2 = algorithms.imageProcessing.util.MatrixUtil.multiplyPointwise(
            x1T, x2RowT1);
        MatrixUtil.multiply(dY2, -1);
        assert(dY2.numColumns() == 3);
        assert(dY2.numRows() == n);
        
        // h = reshape(H',9,1);
        //reshape draws all from column 0 to copy to column 0 of output, etc
        DenseMatrix fmT = MatrixUtil.transpose(fm);
        double[] h2 = new double[9];
        int count = 0;
        for (int col = 0; col < fmT.numColumns(); ++col) {
            for (int row = 0; row < fmT.numRows(); ++row) {
                h2[count] = fmT.get(row, col);
                count++;
            }
        }
        assert(count == 9);

        // d = [Dx * h , Dy * h]';
        //     [nData X 3, nData X 3, nData X 3] * 9 X 1, ...
        //     results in [nData X 3 matrix, nData X 3 matrix]^T

        //nData X 3
        DenseMatrix dXAll = new DenseMatrix(n, 3);
        //nData X 3
        DenseMatrix dYAll = new DenseMatrix(n, 3);
        assert(dX0.numRows() == n);
        // dot each row of dX0 by first 3 elements of h2, etc
        for (int row = 0; row < dX0.numRows(); ++row) {
            double sumdX0 = 0;
            double sumdX1 = 0;
            double sumdX2 = 0;
            double sumdY0 = 0;
            double sumdY1 = 0;
            double sumdY2 = 0;
            for (int col = 0; col < dX0.numColumns(); ++col) {
                sumdX0 = (dX0.get(row, col) * h2[col]);
                sumdX1 = (dX1.get(row, col) * h2[col + 3]);
                sumdX2 = (dX2.get(row, col) * h2[col + 6]);
                sumdY0 = (dY0.get(row, col) * h2[col]);
                sumdY1 = (dY1.get(row, col) * h2[col + 3]);
                sumdY2 = (dY2.get(row, col) * h2[col + 6]);
            }
            // row, col=0 for dX0
            // row, col=1 for dX1
            // row, col=2 for dX2
            dXAll.set(row, 0, sumdX0);
            dXAll.set(row, 1, sumdX1);
            dXAll.set(row, 2, sumdX2);
            dYAll.set(row, 0, sumdY0);
            dYAll.set(row, 1, sumdY1);
            dYAll.set(row, 2, sumdY2);
        }
        
        return new DenseMatrix[]{dXAll, dYAll};
    }

    /**
     * calculate the Sampson's error for the correspondence and given
     * fundamental matrix.
     * 
     <pre>
     The method is adapted from code from the book
     "Multiple View Geometry in Computer Vision" by 
     Hartley and Zisserman, 2004.
       
     obtained from 
     http://www.robots.ox.ac.uk/~vgg/hzbook/

        license is MIT license

        Acknowledgements: These functions are written by: David Capel,
        Andrew Fitzgibbon, Peter Kovesi, Tomas Werner, Yoni Wexler,
        and Andrew Zisserman
     </pre>
     
     * @param fm
     * @param x1
     * @param x2
     * @param tolerance .001
     * @return 
     */
    public EpipolarTransformationFit calculateSampsonsError(DenseMatrix fm,
        DenseMatrix x1, DenseMatrix x2, double tolerance) {

        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        if (x1 == null) {
            throw new IllegalArgumentException("x1 cannot be null");
        }
        if (x2 == null) {
            throw new IllegalArgumentException("x2 cannot be null");
        }
        if (fm.numRows() != 3 || fm.numColumns() != 3) {
            throw new IllegalArgumentException("fm should have 3 rows and 3 columns");
        }
        if (x1.numRows() != 3 || x2.numRows() != 3) {
            throw new IllegalArgumentException("x1 and x2 must "
                + "have 3 rows");
        }
        if (x1.numColumns() != x2.numColumns()) {
            throw new IllegalArgumentException("x1 and x2 must be same sizes");
        }

        int n = x1.numColumns();

        /*        
        geometric error of the final solution or the 7-point sample trial,
        can be approximated by Sampson's error:
             (x2_i * F * x1_i^T)^2                 (x2_i * F * x1_i^T)^2
           ---------------------------------  +  ---------------------------
             (F*x1_i^T)_x^2 + (F*x1_i^T)_y^2     (x2_i*F)_x^2 + (x2_i*F)_y^2
        */
        
        // 3 X nData
        // inverse of x and y coordinates
        DenseMatrix p1 = exRowRepl(x1, 2); 
        DenseMatrix p2 = exRowRepl(x2, 2);
        assert(p1.numRows() == 3);
        assert(p2.numRows() == 3);
        assert(p1.numColumns() == n);
        assert(p2.numColumns() == n);
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < n; ++col) {
                double v = p1.get(row, col)/(x1.get(row, col) + eps);
                p1.set(row, col, v);
                v = p2.get(row, col)/(x2.get(row, col) + eps);
                p2.set(row, col, v);
            }
        }
        
        //DX matrix of size nData X 3, DY matrix of size nData X 3
        DenseMatrix[] alg = algebraicDistance(fm, p1, p2);
        assert(alg[0].numRows() == n);
        assert(alg[0].numColumns() == 3);
        assert(alg[1].numRows() == n);
        assert(alg[1].numColumns() == 3);
        
        DenseMatrix g1 = new DenseMatrix(4, n);
        DenseMatrix g2 = new DenseMatrix(4, n);
        for (int i = 0; i < n; ++i) {
            double v1 = fm.get(0, 0) - p2.get(0, i) * fm.get(2, 0);
            double v2 = fm.get(1, 0) - p2.get(1, i) * fm.get(2, 0);
            g1.set(0, i, v1);
            g2.set(0, i, v2);
            v1 = fm.get(0, 1) - p2.get(0, i) * fm.get(2, 1);
            v2 = fm.get(1, 1) - p2.get(1, i) * fm.get(2, 1);
            g1.set(1, i, v1);
            g2.set(1, i, v2);
            double v = -1 * p1.get(0, i) * fm.get(2, 0) 
                - p1.get(1, i) * fm.get(2, 1) - fm.get(2, 2);
            g1.set(2, i, v);
            g2.set(3, i, v);
        }
        
        // nData X 1        
        double[] magG1 = sumMult(g1, g1);
        double[] magG2 = sumMult(g2, g2);
        assert(magG1.length == n);
        assert(magG2.length == n);
        for (int i = 0; i < n; ++i) {
            magG1[i] = Math.sqrt(magG1[i]);
            magG2[i] = Math.sqrt(magG2[i]);
        }
        // nData X 1
        double[] magG1G2 = sumMult(g1, g2);
        
        double[] alpha = new double[n];
        for (int i = 0; i < n; ++i) {
            double v1 = magG1G2[i];
            double v2 = magG1[i] * magG2[i];
            double v = v1/(v2 + eps);
            if (!(v < 1.)) {
                System.out.format("v1=%.3f, v2=%.3f, v=%.3f\n", 
                        (float)v1, (float)v2, (float)v);
            }
            assert(v < 1.);
            alpha[i] = Math.acos(v);
        }
        
        /*DX, which is alg[0]
        // X1 times 1 * h[0:3]
        // zeroes
        // X1 times -x cooord of X2 * h[6:9]
        
        DY, which is alg[1]
        // zeroes
        // X1 times 1    * h[3:6]
        // X1 times -y cooord of X2 * h[6:9]
        
        p1,p2 are inverse of x and y coords
        
        g1, g2 are composed of x and y inverse times fm
        */
        
        List<Integer> outputInliers = new ArrayList<Integer>();
        List<Double> outputDistances = new ArrayList<Double>();
        
        for (int i = 0; i < n; ++i) {
            double d1 = alg[0].get(i, 0)/(magG1[i] + eps);
            double d2 = alg[1].get(i, 1)/(magG2[i] + eps);
            double d = (
                (d1 * d1) + (d2 * d2) 
                - (2. * d1 * d2 * Math.cos(alpha[i])/Math.sin(alpha[i])));
            
            //System.out.println("d=" + d);
            
            if (d < tolerance) {
                outputInliers.add(Integer.valueOf(i));
                outputDistances.add(Double.valueOf(d));
            }
        }
      
        EpipolarTransformationFit fit = new EpipolarTransformationFit(fm,
            outputInliers, ErrorType.SAMPSONS, outputDistances, tolerance);

        fit.setNMaxMatchable(x1.numColumns());

        return fit;
    }

    //follow errors w/ filter for degeneracy
    public EpipolarTransformationFit calculateSampsonsErrorThenFilter(DenseMatrix fm,
        DenseMatrix x1, DenseMatrix x2, double tolerance) {

        EpipolarTransformationFit fit = calculateSampsonsError(
            fm, x1, x2, tolerance);
        
        List<Integer> outputInliers = fit.getInlierIndexes();
        List<Double> outputDistances = fit.getErrors();

        filterForDegenerate(x1, outputInliers, outputDistances);
        filterForDegenerate(x2, outputInliers, outputDistances);

        EpipolarTransformationFit fit2 = new EpipolarTransformationFit(fm,
            outputInliers, ErrorType.SAMPSONS, outputDistances, tolerance);

        fit2.setNMaxMatchable(x1.numColumns());

        return fit2;
    }

    private void filterForDegenerate(DenseMatrix xy1,
        List<Integer> outputInliers, List<Double> outputDistances) {

        Map<PairInt, List<Integer>> pointIndexes = new HashMap<PairInt, List<Integer>>();

        for (int i = 0; i < outputInliers.size(); ++i) {

            int idx = outputInliers.get(i);

            int x1 = (int)Math.round(xy1.get(0, idx));
            int y1 = (int)Math.round(xy1.get(1, idx));

            PairInt p1 = new PairInt(x1, y1);

            List<Integer> oIndexes = pointIndexes.get(p1);
            if (oIndexes == null) {
                oIndexes = new ArrayList<Integer>();
                pointIndexes.put(p1, oIndexes);
            }
            oIndexes.add(Integer.valueOf(i));
        }

        List<Integer> remove = new ArrayList<Integer>();

        for (Entry<PairInt, List<Integer>> entry : pointIndexes.entrySet()) {

            List<Integer> oIndexes = entry.getValue();
            if (oIndexes.size() < 2) {
                continue;
            }

            double minError = Double.MAX_VALUE;
            int minIdx = -1;
            for (Integer index : oIndexes) {
                int idx = index.intValue();
                double error = outputDistances.get(idx);
                if (error < minError) {
                    minError = error;
                    minIdx = idx;
                }
            }

            assert(minIdx > -1);

            for (Integer index : oIndexes) {
                int idx = index.intValue();
                if (idx == minIdx) {
                    continue;
                }
                remove.add(Integer.valueOf(idx));
            }
        }
        Collections.sort(remove);

        for (int i = (remove.size() - 1); i > -1; --i) {
            int idx = remove.get(i);
            outputDistances.remove(idx);
            outputInliers.remove(idx);
        }
    }

    /**
    calculate the 4 possible projection matrices from the essential matrix.
    * Note that the essential matrix is the transformation matrix between points
    */
    /*
    public DenseMatrix[] calculatePFromEssential(DenseMatrix essentialMatrix) {
    }
    */

}
