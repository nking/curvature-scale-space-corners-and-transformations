package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.features.FeatureComparisonStat;
import algorithms.imageProcessing.features.FeatureMatcher;
import algorithms.imageProcessing.features.IntensityClrFeatures;
import algorithms.imageProcessing.features.KeyPointsAndBounds;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.util.PairFloatArray;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Logger;
import org.ejml.simple.*;

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

 The 7-point algorithm is also implemented below and is similar to the
 8-point solution except that is solves for the null space of the fundamental
 matrix and results in one or 3 solutions which can for some geometries
 be reduced to a single solution.
 The normalization and denormalization steps before and following the solution,
 are the same as in the 8-point solution.
 * 

    right epipolar lines:
        fm * leftXY
    left epipolar lines:
        fm^T * rightXY
   
 * </pre>
 *
 * @author nichole
 */
public class EpipolarTransformer {

    private Logger log = Logger.getLogger(this.getClass().getName());

    /**
     * calculate the fundamental matrix for the given matched left and
     * right correspondence of 8 or more matched points.
     * 
     * @param pointsLeftXY
     * @param pointsRightXY
     * @return 
     */
    public SimpleMatrix calculateEpipolarProjection(
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
    public SimpleMatrix calculateEpipolarProjection(
        SimpleMatrix theLeftXY, SimpleMatrix theRightXY) {

        if (theLeftXY == null) {
            throw new IllegalArgumentException("theLeftXY cannot be null");
        }
        if (theRightXY == null) {
            throw new IllegalArgumentException("refactorRightXY cannot be null");
        }
        if (theLeftXY.numCols()!= theRightXY.numCols()) {
            throw new IllegalArgumentException(
                "theLeftXY and theRightXY must be same size");
        }

        if (theLeftXY.numCols() == 7) {
            throw new IllegalArgumentException(
                "for 7 points, use calculateEpipolarProjectionFor7Points");
        }

        if (theLeftXY.numCols() < 7) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the algorithms require 7 or more points."
                + " refactorLeftXY.n=" +theLeftXY.numCols());
        }

        //the matrix convention is [mRows][nCols]

        SimpleMatrix fundamentalMatrix = calculateFundamentalMatrix(theLeftXY, 
            theRightXY).transpose();

        return fundamentalMatrix;
    }

    protected SimpleMatrix calculateFundamentalMatrix(SimpleMatrix leftXY,
        SimpleMatrix rightXY) {

        //x is xy[0], y is xy[1], xy[2] is all 1's
        NormalizedXY normalizedXY1 = normalize(leftXY);

        NormalizedXY normalizedXY2 = normalize(rightXY);

        return calculateFundamentalMatrix(normalizedXY1, normalizedXY2);
    }

    /*
    for 7-point algorithm:

    (1) SVD of matrix A (as is done in 8-point algorithm)
        giving a matrix of rank 7
    (2) The homogeneous system AX = 0 is called the null space of matrix A.
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
    public List<SimpleMatrix> calculateEpipolarProjectionFor7Points(
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
    public List<SimpleMatrix> calculateEpipolarProjectionFor7Points(
        SimpleMatrix theLeftXY, SimpleMatrix theRightXY) {

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
        if (theLeftXY.numCols() != theRightXY.numCols()) {
            throw new IllegalArgumentException(
                "theLeftXY and theRightXY must be same size");
        }

        if (theLeftXY.numCols() != 7) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the 7-point problem requires 7 points."
                + " theLeftXY.n=" + theLeftXY.numCols());
        }

        //x is xy[0], y is xy[1], xy[2] is all 1's
        NormalizedXY normalizedXY1 = normalize(theLeftXY);

        NormalizedXY normalizedXY2 = normalize(theRightXY);

        double[][] m = createFundamentalMatrix(
            normalizedXY1.getXy(), normalizedXY2.getXy());

        SimpleMatrix aMatrix = new SimpleMatrix(m);
        SimpleSVD<SimpleMatrix> svd = aMatrix.svd();
        SimpleMatrix nullSpace = svd.nullSpace();

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

        SimpleMatrix[] solutions = solveFor7Point(ff1, ff2);

        //denormalize:  F = (T_1)^T * F * T_2
        //    T_1 is normalizedXY1.getNormalizationMatrix();
        //    T2 is normalizedXY2.getNormalizationMatrix();

        List<SimpleMatrix> denormalizedSolutions = new ArrayList<SimpleMatrix>();

        SimpleMatrix t1Transpose = normalizedXY1.getNormalizationMatrix().transpose();
        SimpleMatrix t2 = normalizedXY2.getNormalizationMatrix();

        for (SimpleMatrix solution : solutions) {

            if (solution == null) {
                continue;
            }

            SimpleMatrix denormFundamentalMatrix = t1Transpose.mult(
                solution.mult(t2));

            denormFundamentalMatrix = denormFundamentalMatrix.scale(
                1./denormFundamentalMatrix.get(2, 2));

            denormFundamentalMatrix = denormFundamentalMatrix.transpose();

            /*
            SimpleMatrix validated = validateSolution(denormFundamentalMatrix);

            if (validated != null) {
                denormalizedSolutions.add(validated);
            }*/
                denormalizedSolutions.add(denormFundamentalMatrix);
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
    private SimpleMatrix validateSolution(SimpleMatrix solution, SimpleMatrix leftXY,
        SimpleMatrix rightXY) {

        /*
        function OK = signs_OK(F,x1,x2)
        [u,s,v] = svd(F');
        e1 = v(:,3);
        l1 = vgg_contreps(e1)*x1;
        s = sum( (F*x2) .* l1 );
        OK = all(s>0) | all(s<0);

        (F*x2) .* l1 ==>  (solution * rightXY) .* (testE1 * leftXY)
        '.*' is mattlab notation to operate on each field

        'sum' is a matlab function to sum for each column

        'all' is a function that returns '1' is all items are non-zero, else
            returns 0
        */

        double[][] leftRightEpipoles = calculateEpipoles(solution);

        double[] testE1 = leftRightEpipoles[0];

        SimpleMatrix l1 = leftXY.copy();
        for (int row = 0; row < testE1.length; ++row){
            for (int col = 0; col < l1.numCols(); ++col) {
                double value = testE1[row] * l1.get(row, col);
                l1.set(row, col, value);
            }
        }

        double[] sum = new double[l1.numCols()];
        SimpleMatrix t1 = solution.mult(rightXY);
        for (int row = 0; row < testE1.length; ++row){
            for (int col = 0; col < t1.numCols(); ++col) {
                double value = l1.get(row, col) * t1.get(row, col);
                t1.set(row, col, value);
                sum[col] += value;
            }
        }

        for (int i = 0; i < sum.length; ++i) {
            if (sum[i] == 0) {
                return null;
            }
        }

        return solution;
    }

    SimpleMatrix[] solveFor7Point(double[][] ff1, double[][] ff2) {

        double a0 = calculateCubicRoot3rdOrderCoefficientFor7Point(ff1, ff2);
        double a1 = calculateCubicRoot2ndOrderCoefficientFor7Point(ff1, ff2);
        double a2 = calculateCubicRoot1stOrderCoefficientFor7Point(ff1, ff2);
        double a3 = calculateCubicRoot0thOrderCoefficientFor7Point(ff1, ff2);

        double[] roots = MiscMath.solveCubicRoots(a0, a1, a2, a3);

        double[][] m = new double[3][];
        for (int i = 0; i < 3; i++) {
            m[i] = new double[3];
        }

        SimpleMatrix[] solutions = new SimpleMatrix[roots.length];

        for (int i = 0; i < roots.length; i++) {

            //Fi = a(i)*FF{1} + (1-a(i))*FF{2};

            double a = roots[i];

            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    m[row][col] = a*ff1[row][col] + (1. - a)*ff2[row][col];
                }
            }

            solutions[i] = new SimpleMatrix(m);
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
    SimpleMatrix calculateFundamentalMatrix(NormalizedXY normalizedXY1,
        NormalizedXY normalizedXY2) {

        //build the fundamental matrix
        double[][] m = createFundamentalMatrix(normalizedXY1.getXy(),
            normalizedXY2.getXy());

        /*
        compute linear least square solution:
            solve A = U * D * V^T   for A*f = [..x...]*f = 0
            A has rank 8.  f has rank 2.

        calculate [U,D,V] from svd(A):
           result has mRows = number of data points
                      nCols = 9
        */
        SimpleMatrix aMatrix = new SimpleMatrix(m);
        SimpleSVD<SimpleMatrix> svd = aMatrix.svd();
        SimpleMatrix V = svd.getV();

        // creates U as 9 x nXY1 matrix
        //         D as length 9 array
        //         V as 9 x 9 matrix

        // mRows = 9; nCols = 9

        // reshape V to 3x3

        int vNCols = V.numCols();

        double[][] ff = new double[3][3];
        for (int i = 0; i < 3; i++) {
            ff[i] = new double[3];
            ff[i][0] = V.get((i * 3) + 0, vNCols - 1);
            ff[i][1] = V.get((i * 3) + 1, vNCols - 1);
            ff[i][2] = V.get((i * 3) + 2, vNCols - 1);
        }
        SimpleMatrix fMatrix = new SimpleMatrix(ff);

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

        SimpleMatrix d = svd.getW();

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
        SimpleMatrix dDotV = d.mult(svd.getV().transpose());

        // 3x3
        SimpleMatrix theFundamentalMatrix = svd.getU().mult(dDotV);

        SimpleMatrix denormFundamentalMatrix =
            denormalizeTheFundamentalMatrix(theFundamentalMatrix,
                normalizedXY1, normalizedXY2);

        return denormFundamentalMatrix;
    }

    @SuppressWarnings({"unchecked"})
    SimpleMatrix denormalizeTheFundamentalMatrix(
        SimpleMatrix normalizedFundamentalMatrix,
        NormalizedXY normalizedLeftXY, NormalizedXY normalizedRightXY) {

        /*
        denormalize
            F = (T_1)^T * F * T_2
            where T_1 is normalizedXY1.getNormalizationMatrix();
            and T2 is normalizedXY2.getNormalizationMatrix();
        */

        SimpleMatrix t1Transpose = normalizedLeftXY.getNormalizationMatrix().transpose();
        SimpleMatrix t2 = normalizedRightXY.getNormalizationMatrix();

        SimpleMatrix denormFundamentalMatrix = t1Transpose.mult(
            normalizedFundamentalMatrix.mult(t2));

        denormFundamentalMatrix = denormFundamentalMatrix.scale(
            1./denormFundamentalMatrix.get(2, 2));

        return denormFundamentalMatrix;
    }

    /**
     normalize the x,y coordinates as recommended by Hartley 1997 and return
     the matrix and coordinates.
     does not modify the state of this transformer instance.
     * @param xyPair
     * @return
     */
    @SuppressWarnings({"unchecked"})
    NormalizedXY normalize(SimpleMatrix xy) {

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
        int n = xy.numCols();
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

        SimpleMatrix tMatrix = createScaleTranslationMatrix(scaleFactor, 
            centroidXY[0], centroidXY[1]);
        
        SimpleMatrix normXY = new SimpleMatrix(MatrixUtil.dot(tMatrix, xy));
        
        NormalizedXY normalizedXY = new NormalizedXY();
        normalizedXY.setCentroidXY(centroidXY);
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
    protected SimpleMatrix createScaleTranslationMatrix(double scale, 
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
        SimpleMatrix tMatrix = new SimpleMatrix(t);

        return tMatrix;
    }

    /**
     * write a matrix of size mRows = 3, nCols = xyPairs.getN()
     * @param xyPairs
     * @return
     */
    public SimpleMatrix rewriteInto3ColumnMatrix(PairFloatArray xyPairs) {

        SimpleMatrix xy = new SimpleMatrix(3, xyPairs.getN());

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
    public SimpleMatrix rewriteInto3ColumnMatrix(List<PairInt> xyPairs) {

        SimpleMatrix xy = new SimpleMatrix(3, xyPairs.size());

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
    public SimpleMatrix rewriteFirstItemInto3ColumnMatrix(List<List<PairInt>> xyPairs) {

        SimpleMatrix xy = new SimpleMatrix(3, xyPairs.size());

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
    public SimpleMatrix rewriteInto3ColumnMatrix(PairIntArray xyPairs) {

        SimpleMatrix xy = new SimpleMatrix(3, xyPairs.getN());

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
    double[][] createFundamentalMatrix(SimpleMatrix normXY1,
        SimpleMatrix normXY2) {

        if (normXY1 == null) {
            throw new IllegalArgumentException("normXY1 cannot be null");
        }
        if (normXY2 == null) {
            throw new IllegalArgumentException("normXY2 cannot be null");
        }
        if (normXY1.numCols() != normXY2.numCols()) {
            throw new IllegalArgumentException(
            "the number of columns in normXY1 != number of cols in normXY2");
        }

        int nXY1 = normXY1.numCols();

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
    double[][] calculateEpipoles(SimpleMatrix fundamentalMatrix) {

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
        SimpleSVD<SimpleMatrix> svdE = fundamentalMatrix.svd();
        SimpleMatrix V = svdE.getV().transpose();
        double[] e1 = new double[V.numCols()];
        double e1Div = V.get(2, 2);
        for (int i = 0; i < e1.length; i++) {
            e1[i] = V.get(i, 2)/e1Div;
        }
        SimpleMatrix U = svdE.getU();
        double[] e2 = new double[U.numCols()];
        double e2Div = U.get(2, 2);
        for (int i = 0; i < e2.length; i++) {
            e2[i] = U.get(i, 2)/e2Div;
        }

        double[][] e = new double[2][];
        e[0] = e1;
        e[1] = e2;

        return e;
    }

    private EpipolarFeatureTransformationFit calculateFeatureErrors(SimpleMatrix fm, 
        SimpleMatrix x1, SimpleMatrix x2, 
        IntensityClrFeatures features1, IntensityClrFeatures features2, 
        KeyPointsAndBounds keyPointsAndBounds1, int bmaIndex1, 
        KeyPointsAndBounds keyPointsAndBounds2, int bmaIndex2,
        GreyscaleImage rImg1, GreyscaleImage gImg1, GreyscaleImage bImg1, 
        GreyscaleImage rImg2, GreyscaleImage gImg2, GreyscaleImage bImg2, 
        double tolerance, boolean useHalfDescriptors) {
        
        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        if (x1 == null) {
            throw new IllegalArgumentException("x1 cannot be null");
        }
        if (x2 == null) {
            throw new IllegalArgumentException("x2 cannot be null");
        }
        if (fm.numRows() != 3 || fm.numCols() != 3) {
            throw new IllegalArgumentException("fm should have 3 rows and 3 columns");
        }
        if (x1.numRows() != x2.numRows() || x1.numCols() != x2.numCols()) {
            throw new IllegalArgumentException("x1 and x2 must be same sizes");
        }
        
        FeatureMatcher matcher = new FeatureMatcher();
        
        List<Integer> outputInliers = new ArrayList<Integer>();
        List<Double> outputDistances = new ArrayList<Double>();
        List<FeatureComparisonStat> fcs = new ArrayList<FeatureComparisonStat>();
        
        int n = x1.numCols();
        
        for (int col = 0; col < n; ++col) {
            int xPt1 = (int)Math.round(x1.get(0, col));
            int yPt1 = (int)Math.round(x1.get(1, col));
            int xPt2 = (int)Math.round(x2.get(0, col));
            int yPt2 = (int)Math.round(x2.get(1, col));
            
            FeatureComparisonStat stat;
            if (useHalfDescriptors) {
                stat = matcher.matchHalfDescriptors(
                    features1, features2,
                    keyPointsAndBounds1, bmaIndex1, keyPointsAndBounds2, bmaIndex2,
                    xPt1, yPt1, xPt2, yPt2,
                    rImg1, gImg1, bImg1, rImg2, gImg2, bImg2);
            } else {
                stat = matcher.matchDescriptors(
                    features1, features2, xPt1, yPt1, xPt2, yPt2,
                    rImg1, gImg1, bImg1, rImg2, gImg2, bImg2);
            }
            
            if (stat == null || 
                (stat.getSumIntensitySqDiff() > stat.getImg2PointIntensityErr())) {
                continue;
            }
            
            outputDistances.add(Double.valueOf(stat.getSumIntensitySqDiff()));
            outputInliers.add(Integer.valueOf(col));
            fcs.add(stat);
        }
        
        EpipolarFeatureTransformationFit fit = new EpipolarFeatureTransformationFit(
            fm, outputInliers, fcs, ErrorType.SAMPSONS, outputDistances, tolerance);
    
        fit.setNMaxMatchable(x1.numCols());
        
        return fit;
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
        
        EpipolarFeatureTransformationFit fit = new EpipolarFeatureTransformationFit(
            distanceErrors.getFundamentalMatrix(),
            outputInliers, fcs,
            distanceErrors.getErrorType(), outputDistances, 
            distanceErrors.getTolerance());
            
        return fit;
    }

    public static class NormalizedXY {

        /**
         * 3 dimensional matrix, with column 0 being x, column 1 being y,
         * and the last column is place holder 1's
         */
        private SimpleMatrix xy = null;

        private double[] centroidXY = null;

        private SimpleMatrix normalizationMatrix = null;

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
        public SimpleMatrix getNormalizationMatrix() {
            return normalizationMatrix;
        }

        /**
         * @param normMatrix holding the scale and offsets to apply to x, y
         */
        public void setNormMatrix(SimpleMatrix normMatrix) {
            this.normalizationMatrix = normMatrix;
        }

        /**
         * @return the xy
         */
        public SimpleMatrix getXy() {
            return xy;
        }

        /**
         * @param xy the xy to set
         */
        public void setXy(SimpleMatrix xy) {
            this.xy = xy;
        }
    }

    public PairIntArray getEpipolarLine(SimpleMatrix epipolarLines, int imgWidth,
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
        SimpleMatrix fm, SimpleMatrix leftPoints, SimpleMatrix rightPoints, 
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
        int nRows = leftPoints.getMatrix().getNumRows();
        if (nRows != rightPoints.getMatrix().getNumRows()) {
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

        fit.setNMaxMatchable(leftPoints.numCols());

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
        SimpleMatrix fm, SimpleMatrix leftPoints, SimpleMatrix rightPoints, 
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
        int nRows = leftPoints.getMatrix().getNumRows();
        if (nRows != rightPoints.getMatrix().getNumRows()) {
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

        fit.setNMaxMatchable(leftPoints.numCols());

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
    PairFloatArray calculateDistancesFromEpipolar(
        SimpleMatrix fm, SimpleMatrix matchedLeftPoints,
        SimpleMatrix matchedRightPoints) {

        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        if (matchedLeftPoints == null) {
            throw new IllegalArgumentException("matchedLeftPoints cannot be null");
        }
        if (matchedRightPoints == null) {
            throw new IllegalArgumentException("rightPoints cannot be null");
        }
        int nRows = matchedLeftPoints.getMatrix().getNumRows();
        if (nRows != matchedRightPoints.getMatrix().getNumRows()) {
            throw new IllegalArgumentException("matrices must have same number of rows");
        }

        int n = matchedLeftPoints.numCols();

        PairFloatArray distances = new PairFloatArray(n);

        SimpleMatrix rightEpipolarLines = fm.mult(matchedLeftPoints);

        SimpleMatrix leftEpipolarLines = fm.transpose().mult(matchedRightPoints);

        for (int i = 0; i < matchedLeftPoints.numCols(); i++) {

            double a = rightEpipolarLines.get(0, i);
            double b = rightEpipolarLines.get(1, i);
            double c = rightEpipolarLines.get(2, i);

            double aplusb = Math.sqrt((a*a) + (b*b));

            double xL = matchedLeftPoints.get(0, i);
            double yL = matchedLeftPoints.get(1, i);

            //dist = (a*x + b*y + c)/sqrt(a^2 + b^2)

            double x = matchedRightPoints.get(0, i);
            double y = matchedRightPoints.get(1, i);

            double d = (a*x + b*y + c)/aplusb;

            // find the reverse distance by projection:
            double aRev = leftEpipolarLines.get(0, i);
            double bRev = leftEpipolarLines.get(1, i);
            double cRev = leftEpipolarLines.get(2, i);

            double dRev = (aRev*xL + bRev*yL + cRev)/
                Math.sqrt((aRev*aRev + bRev*bRev));

            distances.add((float)dRev, (float)d);
        }

        return distances;
    }
    
    public EpipolarTransformationFit calculateErrorThenFilter(SimpleMatrix fm,
        SimpleMatrix x1, SimpleMatrix x2, ErrorType errorType, double tolerance) {
         
        if (errorType.equals(ErrorType.SAMPSONS)) {
            return calculateSampsonsErrorThenFilter(fm, x1, x2, tolerance);
        } else {
            return calculateEpipolarDistanceErrorThenFilter(fm, x1, x2, tolerance);
        }
    }
    
    public EpipolarTransformationFit calculateError(SimpleMatrix fm,
        SimpleMatrix x1, SimpleMatrix x2, ErrorType errorType, double tolerance) {
         
        if (errorType.equals(ErrorType.SAMPSONS)) {
            return calculateSampsonsError(fm, x1, x2, tolerance);
        } else {
            return calculateEpipolarDistanceError(fm, x1, x2, tolerance);
        }
    }
    
    public EpipolarFeatureTransformationFit calculateError(SimpleMatrix fm,
        SimpleMatrix x1, SimpleMatrix x2, 
        IntensityClrFeatures features1, IntensityClrFeatures features2,
        KeyPointsAndBounds kpab1, int bmaIndex1,
        KeyPointsAndBounds kpab2, int bmaIndex2,
        GreyscaleImage rImg1, GreyscaleImage gImg1, GreyscaleImage bImg1,
        GreyscaleImage rImg2, GreyscaleImage gImg2, GreyscaleImage bImg2,
        ErrorType errorType, double tolerance, boolean useHalfDescriptors) {
         
        EpipolarTransformationFit distanceErrors;
        
        if (errorType.equals(ErrorType.SAMPSONS)) {            
            distanceErrors = calculateSampsonsError(fm, x1, x2, tolerance);
        } else {
            distanceErrors = calculateEpipolarDistanceError(fm, x1, x2, tolerance);
        }
        
        EpipolarFeatureTransformationFit featureErrors = calculateFeatureErrors(
            fm, x1, x2, features1, features2,
            kpab1, bmaIndex1, kpab2, bmaIndex2,
            rImg1, gImg1, bImg1, rImg2, gImg2, bImg2,
            tolerance, useHalfDescriptors);
        
        EpipolarFeatureTransformationFit combinedFit = combineErrors(distanceErrors,
            featureErrors);
        combinedFit.setNMaxMatchable(x1.numCols());
        
        return combinedFit;
    }
    
    public EpipolarTransformationFit calculateSampsonsError(SimpleMatrix fm,
        SimpleMatrix x1, SimpleMatrix x2, double tolerance) {
        
        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        if (x1 == null) {
            throw new IllegalArgumentException("x1 cannot be null");
        }
        if (x2 == null) {
            throw new IllegalArgumentException("x2 cannot be null");
        }
        if (fm.numRows() != 3 || fm.numCols() != 3) {
            throw new IllegalArgumentException("fm should have 3 rows and 3 columns");
        }
        if (x1.numRows() != x2.numRows() || x1.numCols() != x2.numCols()) {
            throw new IllegalArgumentException("x1 and x2 must be same sizes");
        }
        
        /*
        geometric error of the final solution or the 7-point sample trial,
        can be approximated by Sampon's error:
             (x2_i * F * x1_i^T)^2                 (x2_i * F * x1_i^T)^2
           ---------------------------------  +  ---------------------------
             (F*x1_i^T)_x^2 + (F*x1_i^T)_y^2     (x2_i*F)_x^2 + (x2_i*F)_y^2
        
        reference?
        */
        
        List<Integer> outputInliers = new ArrayList<Integer>();
        List<Double> outputDistances = new ArrayList<Double>();
                
        SimpleMatrix fmT = fm.transpose();
                
        int n = x1.numCols();
        
        // 1 x n
        SimpleMatrix x2tFx1 = new SimpleMatrix(1, n);
        
        boolean extractRow = false;
        for (int col = 0; col < n; ++col) {
            // 3 x 1 ==> T ==> 1 x 3
            SimpleMatrix x2T_i = x2.extractVector(extractRow, col).transpose();
            
            // x2T_i * F is 1X3 * 3X3 = 1X3
            SimpleMatrix x2T_iF = x2T_i.mult(fm);
            
            // 3 x 1
            SimpleMatrix x1_i = x1.extractVector(extractRow, col);
            
            // x2T_iF * x1_i is 1X3 * 3X1 = 1X1
            SimpleMatrix result = x2T_iF.mult(x1_i);
            
            x2tFx1.set(0, col, result.get(0, 0));            
        }
        
        //Fx1 = F * x1  is 3X3 * 3Xn = 3Xn
        SimpleMatrix fx1 = fm.mult(x1);
        
        //Ftx2 = F' * x2  is 3x3 * 3xn = 3xn
        SimpleMatrix ftx2 = fmT.mult(x2);
        
        for (int col = 0; col < n; ++col) {
            double x2tFx1_j = x2tFx1.get(0, col);
            double a = x2tFx1_j * x2tFx1_j;
            
            double t1 = fx1.get(0, col);
            t1 *= t1;
            
            double t2 = fx1.get(1, col);
            t2 *= t2;
            
            double t3 = ftx2.get(0, col);
            t3 *= t3;
            
            double t4 = ftx2.get(1, col);
            t4 *= t4;
            
            double b = t1 + t2 + t3 + t4;
            
            double error = a/b;
                        
            if (error < tolerance) {
                outputInliers.add(Integer.valueOf(col));
                outputDistances.add(Double.valueOf(error));
            }
        }
        
        EpipolarTransformationFit fit = new EpipolarTransformationFit(fm,
            outputInliers, ErrorType.SAMPSONS, outputDistances, tolerance);
    
        fit.setNMaxMatchable(x1.numCols());
        
        return fit;
    }
    
    //follow errors w/ filter for degeneracy
    public EpipolarTransformationFit calculateSampsonsErrorThenFilter(SimpleMatrix fm,
        SimpleMatrix x1, SimpleMatrix x2, double tolerance) {
        
        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        if (x1 == null) {
            throw new IllegalArgumentException("x1 cannot be null");
        }
        if (x2 == null) {
            throw new IllegalArgumentException("x2 cannot be null");
        }
        if (fm.numRows() != 3 || fm.numCols() != 3) {
            throw new IllegalArgumentException("fm should have 3 rows and 3 columns");
        }
        if (x1.numRows() != x2.numRows() || x1.numCols() != x2.numCols()) {
            throw new IllegalArgumentException("x1 and x2 must be same sizes");
        }
        
        /*
        geometric error of the final solution or the 7-point sample trial,
        can be approximated by Sampon's error:
             (x2_i * F * x1_i^T)^2                 (x2_i * F * x1_i^T)^2
           ---------------------------------  +  ---------------------------
             (F*x1_i^T)_x^2 + (F*x1_i^T)_y^2     (x2_i*F)_x^2 + (x2_i*F)_y^2
        
        reference?
        */
        
        List<Integer> outputInliers = new ArrayList<Integer>();
        List<Double> outputDistances = new ArrayList<Double>();
                
        SimpleMatrix fmT = fm.transpose();
                
        int n = x1.numCols();
        
        // 1 x n
        SimpleMatrix x2tFx1 = new SimpleMatrix(1, n);
        
        boolean extractRow = false;
        for (int col = 0; col < n; ++col) {
            // 3 x 1 ==> T ==> 1 x 3
            SimpleMatrix x2T_i = x2.extractVector(extractRow, col).transpose();
            
            // x2T_i * F is 1X3 * 3X3 = 1X3
            SimpleMatrix x2T_iF = x2T_i.mult(fm);
            
            // 3 x 1
            SimpleMatrix x1_i = x1.extractVector(extractRow, col);
            
            // x2T_iF * x1_i is 1X3 * 3X1 = 1X1
            SimpleMatrix result = x2T_iF.mult(x1_i);
            
            x2tFx1.set(0, col, result.get(0, 0));            
        }
        
        //Fx1 = F * x1  is 3X3 * 3Xn = 3Xn
        SimpleMatrix fx1 = fm.mult(x1);
        
        //Ftx2 = F' * x2  is 3x3 * 3xn = 3xn
        SimpleMatrix ftx2 = fmT.mult(x2);
        
        for (int col = 0; col < n; ++col) {
            double x2tFx1_j = x2tFx1.get(0, col);
            double a = x2tFx1_j * x2tFx1_j;
            
            double t1 = fx1.get(0, col);
            t1 *= t1;
            
            double t2 = fx1.get(1, col);
            t2 *= t2;
            
            double t3 = ftx2.get(0, col);
            t3 *= t3;
            
            double t4 = ftx2.get(1, col);
            t4 *= t4;
            
            double b = t1 + t2 + t3 + t4;
            
            double error = a/b;
                        
            if (error < tolerance) {
                outputInliers.add(Integer.valueOf(col));
                outputDistances.add(Double.valueOf(error));
            }
        }
        
        filterForDegenerate(x1, outputInliers, outputDistances);
        filterForDegenerate(x2, outputInliers, outputDistances);
        
        EpipolarTransformationFit fit = new EpipolarTransformationFit(fm,
            outputInliers, ErrorType.SAMPSONS, outputDistances, tolerance);
    
        fit.setNMaxMatchable(x1.numCols());
        
        return fit;
    }
    
    private void filterForDegenerate(SimpleMatrix xy1,
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
    * that are in camera coordinates rather than in pixel coordinates.
    the method is adapted from: https://github.com/jesolem/PCV
    Their code is licensed under  BSD license (2-clause "Simplified BSD License").
    */
    /*
    public SimpleMatrix[] calculatePFromEssential(SimpleMatrix essentialMatrix) {
    }
    */

}
