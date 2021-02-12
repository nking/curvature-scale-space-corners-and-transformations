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
import java.util.Arrays;
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
    vectors are treated as columns unless noted otherwise.
    vector as a row uses notation u^T.
    u^T*v represents the inner product.
    u*v^T is a matrix.
    u_1 is the (x,y) points from image 1 and u_2 are the matched (x,y) points
        from image 2.

    the fundamental matrix is defined:
        u_2^T * F * u_1 = 0  
        where u are the x,y points in images _1 and _2
        and F is a 3 × 3 matrix of rank 2
        
    An explanation of the derivation of the fundamental matrix can be found in
    Zhengyou Zhang. "Determining the Epipolar Geometry and its Uncertainty: A
    Review". RR-2927, INRIA. 1996. ffinria-00073771
    equations (1) and (2) 
    
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
    and f is a nine-vector containing the entries of the matrix F
 
    To avoid the trivial scale, ||f|| = 1 where f is the norm of f

    the rank of f is then 8.
    if A is longer than 8, the system is over-determined (over specified) and
    so must be solved using least squares.  the set may be over determined
    and not have a zero solution.

    we want the vector f that minimizes ||A*f|| subject to the constraint
    that ||f|| = f^T*f = 1

    the solution is the unit eigenvector corresponding to the smallest
    eigenvalue of A^T*A.

    Since A^T*A is semi-definite and symmetric, all of its eigenvectors
    are real and positive or zero.
    This eigenvector is what he calls the least eigenvector of A^T*A and
    it is found via the Jacobi algorithm or Singular Value Decomposition.
    NOTE: A^T*A is V * D * D^T * V^T
      and A*A^T is U * D * D^T * U^T
      adn that for the SVD of A^T*A and that of A*A^T the U and V vectors equal one another.

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
     1) The points are translated so that their centroid is at
        the origin.
     2) The points are then scaled so that the average distance
        from the origin is equal to 2 .
     3) This transformation is applied to each of the two images independently.

     NOTE: if needed to use non-isotropic scaling (e.g. rectangular pixels, etc):
     transform the points so that
     1) Their centroid is at the origin.
     2) The principal moments are both equal to unity
      
     
     scaling the coordinate so that the homogeneous coordinates are on the 
     average equal to unity will improve the condition of the matrix A^T*A.

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
     singular values (similar to dimensionality reduction)
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
            Util.rewriteInto3ColumnMatrix(pointsLeftXY),
            Util.rewriteInto3ColumnMatrix(pointsRightXY));
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

        DenseMatrix fundamentalMatrix = algorithms.matrix.MatrixUtil.transpose(
            calculateFundamentalMatrix(theLeftXY, theRightXY));

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
            Util.rewriteInto3ColumnMatrix(pointsLeftXY),
            Util.rewriteInto3ColumnMatrix(pointsRightXY));
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
        // These are the last 2 columns of V
        DenseMatrix nullSpace = algorithms.imageProcessing.util.MatrixUtil.nullSpace(svd);


        double[][] ff1 = new double[3][3];
        double[][] ff2 = new double[3][3];
        for (int i = 0; i < 3; i++) {

            ff1[i] = new double[3];
            ff2[i] = new double[3];
            
            ff1[i][0] = nullSpace.get((i * 3) + 0, 0);
            ff1[i][1] = nullSpace.get((i * 3) + 1, 0);
            ff1[i][2] = nullSpace.get((i * 3) + 2, 0);
            ff2[i][0] = nullSpace.get((i * 3) + 0, 1);
            ff2[i][1] = nullSpace.get((i * 3) + 1, 1);
            ff2[i][2] = nullSpace.get((i * 3) + 2, 1);
        }

        DenseMatrix[] solutions = solveFor7Point(ff1, ff2);
        //DenseMatrix[] solutions = solveFor7Point2(ff1, ff2);
        
        //denormalize:  F = (T_1)^T * F * T_2
        //    T_1 is normalizedXY1.getNormalizationMatrix();
        //    T2 is normalizedXY2.getNormalizationMatrix();

        List<DenseMatrix> denormalizedSolutions = new ArrayList<DenseMatrix>();

        DenseMatrix t1Transpose = MatrixUtil.transpose(normalizedXY1.getNormalizationMatrix());
        DenseMatrix t2 = normalizedXY2.getNormalizationMatrix();

        for (DenseMatrix solution : solutions) {

            DenseMatrix validated = validateSolution(solution, 
                normalizedXY1.getXy(), normalizedXY2.getXy());
            
            if (validated == null) {
                continue;
            }
                            
            DenseMatrix denormFundamentalMatrix
                = MatrixUtil.multiply(t1Transpose,
                MatrixUtil.multiply(validated, t2));

            double s = 1. / (denormFundamentalMatrix.get(2, 2) + eps);
            MatrixUtil.multiply(denormFundamentalMatrix, s);

            denormFundamentalMatrix = (DenseMatrix) denormFundamentalMatrix.transpose();

            denormalizedSolutions.add(denormFundamentalMatrix);
            
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
    
    /**
     * This method is ported from the matlab code from the book
     * @Book{Hartley2004,
            author = "Hartley, R.~I. and Zisserman, A.",
            title = "Multiple View Geometry in Computer Vision",
            edition = "Second",
            year = "2004",
            publisher = "Cambridge University Press, ISBN: 0521540518"
        }
     * https://www.robots.ox.ac.uk/~vgg/hzbook/code/
     * The code is available using the MIT license.
     * 
     * code VGG_SINGF_FROM_FF.m
     * 
     * @param ff1  3X3
     * @param ff2  3X3
     * @return 
     */
    DenseMatrix[] solveFor7Point2(double[][] ff1, double[][] ff2) {

        double[][][] d = new double[2][2][2];
        for (int i1 = 0; i1 < 2; ++i1) {
            d[i1] = new double[2][2];
            for (int i2 = 0; i2 < 2; ++i2) {
                d[i1][i2] = new double[2];
            }
        }
        
        double[][] f1, f2, f3;
        double[][] tempF = new double[3][3];
        for (int i = 0; i < 3; ++i) {
            tempF[i] = new double[3];
        }
        
        for (int i1 = 0; i1 < 2; ++i1) {
            if (i1 == 0) {
                f1 = ff1;
            } else {
                f1 = ff2;
            }
            tempF[0][0] = f1[0][0];
            tempF[0][1] = f1[1][0];
            tempF[0][2] = f1[2][0];
            for (int i2 = 0; i2 < 2; ++i2) {
                if (i2 == 0) {
                    f2 = ff1;
                } else {
                    f2 = ff2;
                }
                tempF[1][0] = f2[0][1];
                tempF[1][1] = f2[1][1];
                tempF[1][2] = f2[2][1];
                for (int i3 = 0; i3 < 2; ++i3) {
                    if (i3 == 0) {
                        f3 = ff1;
                    } else {
                        f3 = ff2;
                    }
                    tempF[2][0] = f3[0][2];
                    tempF[2][1] = f3[1][2];
                    tempF[2][2] = f3[2][2];
                    //d[i1][i2][i3] = det([F{i1}(:,1) F{i2}(:,2) F{i3}(:,3)]);
                    d[i1][i2][i3] = MatrixUtil.determinant(tempF);
                }
            }
        }
        
        //-D(2,1,1)+D(1,2,2)+D(1,1,1)+D(2,2,1)+D(2,1,2)-D(1,2,1)-D(1,1,2)-D(2,2,2)
        double a0 = -d[1][0][0] + d[0][1][1] + d[0][0][0] + d[1][1][0]
           + d[1][0][1] - d[0][1][0] - d[0][0][1] - d[1][1][1];
        
        //D(1,1,2)-2*D(1,2,2)-2*D(2,1,2)+D(2,1,1)-2*D(2,2,1)+D(1,2,1)+3*D(2,2,2)
        double a1 = d[0][0][1] - 2.*d[0][1][1] - 2.*d[1][0][1] + d[1][0][0]
           - 2.*d[1][1][0] + d[0][1][0] + 3.*d[1][1][1];
                  
        //D(2,2,1)+D(1,2,2)+D(2,1,2)-3*D(2,2,2)
        double a2 = d[1][1][0] + d[0][1][1] + d[1][0][1] - 3.*d[1][1][1];
                
        //D(2,2,2)]
        double a3 = d[1][1][1];
                
        //solve for the roots of equation a0 * x^3 + a1 * x^2 + a2 * x + a3 = 0;
        
        double _a0 = calculateCubicRoot3rdOrderCoefficientFor7Point(ff1, ff2);
        double _a1 = calculateCubicRoot2ndOrderCoefficientFor7Point(ff1, ff2);
        double _a2 = calculateCubicRoot1stOrderCoefficientFor7Point(ff1, ff2);
        double _a3 = calculateCubicRoot0thOrderCoefficientFor7Point(ff1, ff2);
        
        double[] roots = MiscMath.solveCubicRoots(a0, a1, a2, a3);
        
        double[] _roots = MiscMath.solveCubicRoots(_a0, _a1, _a2, _a3);

        
        System.out.println("a0=" + a0);
        System.out.println("a1=" + a1);
        System.out.println("a2=" + a2);
        System.out.println("a3=" + a3);
        System.out.println("_a0=" + _a0);
        System.out.println("_a1=" + _a1);
        System.out.println("_a2=" + _a2);
        System.out.println("_a3=" + _a3);
        
        System.out.println("roots=" + Arrays.toString(roots));
        System.out.println("_roots=" + Arrays.toString(_roots));
        System.out.flush();
        
        double[][] m = new double[3][];
        for (int i = 0; i < 3; i++) {
            m[i] = new double[3];
        }

        DenseMatrix[] solutions = new DenseMatrix[roots.length];

        for (int i = 0; i < roots.length; i++) {

            //Fi = a(i)*FF{1} + (1-a(i))*FF{2};

            double a = roots[i];
            double aa = 1. - a;
            
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    m[row][col] = (a*ff1[row][col]) + (aa*ff2[row][col]);
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

        //aMatrix is m x n  (== nData X 9)
        // U   is  m X m     the left singular vectors, column-wise. Not available for partial decompositions
        // S   is  min(m, n) the singular values (stored in descending order)
        // V^T is  n X n     the right singular vectors, row-wise. Not available for partial decompositions
        SVD svd = null;
        try {
            svd = SVD.factorize(aMatrix);
        } catch (NotConvergedException e) {
            log.severe(e.getMessage());
            return null;
        }

        //A is nData rows X 3 columns

        DenseMatrix V = algorithms.matrix.MatrixUtil.transpose(svd.getVt());

        /*
        DenseMatrix V = MatrixUtil.transpose(svd.getVt());

        System.out.println("A=" + aMatrix.toString());
        System.out.println("U=" + svd.getU().toString());
        System.out.println("S=" + Arrays.toString(svd.getS()));
        System.out.println("V^T=" + svd.getVt().toString());
        System.out.println("V=" + V.toString());
        System.out.flush();
        */
        
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
        d = (DenseMatrix)d.zero();
        if (sDiag.length > 0) {
            d.set(0, 0, sDiag[0]);
        }
        if (sDiag.length > 1) {
            d.set(1, 1, sDiag[1]);
        }

        /*
        multiply the terms:
             F = dot(U, dot(diag(D),V^T))
        */
        DenseMatrix dDotV = MatrixUtil.multiply(d, svd.getVt());

        // 3x3
        DenseMatrix theFundamentalMatrix = MatrixUtil.multiply(svd.getU(), dDotV);

        //System.out.println("fm before de-normalization=" + theFundamentalMatrix.toString());

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

        DenseMatrix t1Transpose = 
            algorithms.matrix.MatrixUtil.transpose(normalizedLeftXY
            .getNormalizationMatrix());
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

        /*
        double[][] t = new double[3][];
        t[0] = new double[]{scale,       0,     -centroidX*scale};
        t[1] = new double[]{0,           scale, -centroidY*scale};
        t[2] = new double[]{0,           0,           1};
        DenseMatrix tMatrix = new DenseMatrix(t);
        
        xy is size 3 X nData
        
        (x_0*scale-centroidX*scale) ...for i=1 to n
        (y_0*scale-centroidY*scale)
        (1)
        */
        
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
     * @param normXY1 a matrix of size 3 x nPoints, where 1st row is x,
     * second is y.
     * @param normXY2 a matrix of size 3 x nPoints, where 1st row is x,
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
        DenseMatrix V = algorithms.matrix.MatrixUtil.transpose(svdE.getVt());
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
    calculate the 4 possible projection matrices from the essential matrix.
    * Note that the essential matrix is the transformation matrix between points
    */
    /*
    public DenseMatrix[] calculatePFromEssential(DenseMatrix essentialMatrix) {
    }
    */

}
