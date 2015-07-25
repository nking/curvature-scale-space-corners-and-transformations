package algorithms.imageProcessing;

import algorithms.util.PairFloatArray;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;
import thirdparty.HungarianAlgorithm;
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
 * Matching the points before using this can be done with PointMatcher.
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
 * </pre>
 * 
 * @author nichole
 */
public class StereoProjectionTransformer {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    private SimpleMatrix leftXY = null;
        
    private SimpleMatrix rightXY = null;
        
    private SimpleMatrix fundamentalMatrix = null;
    
    private double[] leftEpipole = null;
    
    private double[] rightEpipole = null;
    
    /**
     * each row is an epipolar line in the right image.
     * Each column corresponds to a point in leftXY and rightXY which are
     * in the same column.
     */
    private SimpleMatrix epipolarLinesInRight = null;
    
    /**
     * calculate the epipolar projection for a set of 8 or more unmatched
     * points.  It uses PointMatcher, RANSACSolver and PointPartitioner
     * to obtain a solution for part of the images that best fits all of
     * the points.
     */
    private SimpleMatrix epipolarLinesInLeft = null;
    
    private boolean fitIsBetter(StereoProjectionTransformerFit bestFit, 
        StereoProjectionTransformerFit compareFit) {
        
        if (compareFit == null) {
            return false;
        }
        if (bestFit == null) {
            return true;
        }
       
        long bestN = bestFit.getNMatches();
        long compareN = compareFit.getNMatches();
        
        if (compareN == 0) {
            return false;
        }
        if (bestN == 0) {
            return true;
        }
        
        double bestStat = (double)bestN/(double)bestFit.getNMaxMatchable();
        
        double compareStat = (double)compareN/(double)compareFit.getNMaxMatchable();
        
        if (compareStat > bestStat) {
            
            return true;
            
        } else if (compareStat == bestStat) {
            
            if (!Double.isNaN(compareFit.getMeanDistance()) && (
                compareFit.getMeanDistance()
                < bestFit.getMeanDistance())) {
                
                return true;
                
            } else if (compareFit.getMeanDistance()
                == bestFit.getMeanDistance()) {
                
                if (compareFit.getStDevFromMean()< bestFit.getStDevFromMean()) {
                    
                    return true;
                }
            }
        }
        
        return false;
    }
    
    public SimpleMatrix calculateEpipolarProjectionForPerfectlyMatched(
        PairFloatArray pointsLeftXY,  PairFloatArray pointsRightXY) {
        
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
        
        if (pointsLeftXY.getN() == 7) {
            throw new IllegalArgumentException(
                "for 7 points, use calculateEpipolarProjectionFor7Points");
        }
        
        if (pointsLeftXY.getN() < 7) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the algorithms require 7 or more points." 
                + " refactorLeftXY.n=" + pointsLeftXY.getN());
        }
        
        return calculateEpipolarProjectionForPerfectlyMatched(
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
    public SimpleMatrix calculateEpipolarProjectionForPerfectlyMatched(
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
        
        leftXY = theLeftXY;
        
        rightXY = theRightXY;
        
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
        
        return fundamentalMatrix;
    }
    
    /**
     * NOTE: this method should only be used for comparison.  Prefer 
     * calculateEpipolarProjection().
     * 
     * @param pointsLeftXY
     * @param pointsRightXY 
     */
    public SimpleMatrix calculateEpipolarProjectionWithoutNormalization(
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
        
        //the matrix convention is [mRows][nCols]
        
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
        
        leftXY = theLeftXY;
        
        rightXY = theRightXY;
        
        //x is xy[0], y is xy[1], xy[2] is all 1's
        NormalizedXY normalizedXY1 = normalize(leftXY);
        
        NormalizedXY normalizedXY2 = normalize(rightXY);      
        
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
            
            SimpleMatrix validated = validateSolution(denormFundamentalMatrix);
            
            if (validated != null) {
                denormalizedSolutions.add(validated);
            }
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
    private SimpleMatrix validateSolution(SimpleMatrix solution) {
                            
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
        
       denormFundamentalMatrix = denormFundamentalMatrix.scale(
            1./denormFundamentalMatrix.get(2, 2));
        
        return denormFundamentalMatrix;
    }
    
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
    SimpleMatrix calculateFundamentalMatrixWithoutNormalization(
        SimpleMatrix matchedXY1, SimpleMatrix matchedXY2) {
        
        //build the fundamental matrix
        SimpleMatrix aMatrix = new SimpleMatrix(createFundamentalMatrix(
            matchedXY1, matchedXY2));

        /*
        compute linear least square solution:
            solve A = U * D * V^T   for A*f = [..x...]*f = 0
            A has rank 8.  f has rank 2.
        
        calculate [U,D,V] from svd(A):
           result has mRows = number of data points
                      nCols = 9
        */
        SimpleSVD<SimpleMatrix> svd = aMatrix.svd();

        // creates U as 9 x nXY1 matrix
        //         D as length 9 array
        //         V as 9 x 9 matrix
        
        // mRows = 9; nCols = 9
        SimpleMatrix V = svd.getV();
        int vNCols = V.numCols();
        
        // reshape it to 3x3
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
        
           then from [U,D,V], create F:
              F = U * diag([D(1,1) D(2,2) 0]) * V^T;
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
        
        theFundamentalMatrix = theFundamentalMatrix.scale(
            1./theFundamentalMatrix.get(2, 2));
        
        return theFundamentalMatrix;
    } 

    /**
     normalize the x,y coordinates as recommended by Hartley 1997 and return
     the matrix and coordinates.
     does not modify the state of this transformer instance.
     * @param xyPair
     * @return 
     */
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
        SimpleMatrix tMatrix = new SimpleMatrix(t);
                
        SimpleMatrix normXY = new SimpleMatrix(MatrixUtil.dot(tMatrix, xy));
              
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
    public static SimpleMatrix rewriteInto3ColumnMatrix(PairFloatArray xyPairs) {
        
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
        
        SimpleMatrix xy = new SimpleMatrix(xyPoints);
        
        return xy;
    }
    
    /**
     * write a matrix of size mRows = 3, nCols = xyPairs.getN()
     * @param xyPairs
     * @return 
     */
    public static SimpleMatrix rewriteInto3ColumnMatrix(PairIntArray xyPairs) {
        
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
        
        SimpleMatrix xy = new SimpleMatrix(xyPoints);
        
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

    private SimpleMatrix calculateRightEpipolarLines() {
        
        /* calculate right epipolar lines
        F * leftPoint
        */
        
        SimpleMatrix m = fundamentalMatrix.mult(leftXY);
        
        return m;
    }
    
    private SimpleMatrix calculateLeftEpipolarLines() {
        
        //calculate left epipolar lines:  F^T * rightPoint
        
        SimpleMatrix fundamentalMatrixTranspose = fundamentalMatrix.transpose();
        
        SimpleMatrix m = fundamentalMatrixTranspose.mult(rightXY);
        
        return m;
    }

    /**
     * NOT YET READY FOR USE.
     * using the homography matrix bestFitProj, project the points in one set 
     * to epipolar lines in the other set and find the closest match that
     * works best in both frames.  The cost function for the match is the
     * perpendicular distance of a point to the epipolar line.
     * Currently, no information about the points position parallel to the
     * line is used.
     * 
     * @param fm projective transformation coefficients between left and right
     * @param unMatchedLeft
     * @param unMatchedRight
     * @param tolerance
     * @param outputMatchedLeft the output matched points for the left image.
     * the indexes for this are the same as the indexes for outputMatchedRight.
     * @param outputMatchedRight the output matched points for the right
     * image.  the indexes are the same as those in outputMatchedLeft.
     * @return 
     */
    public StereoProjectionTransformerFit evaluateFitAndMatchPoints(
        SimpleMatrix fm, PairIntArray unMatchedLeft, PairIntArray 
        unMatchedRight, double tolerance,
        PairIntArray outputMatchedLeft, PairIntArray outputMatchedRight) {
        
        SimpleMatrix in1 = rewriteInto3ColumnMatrix(unMatchedLeft);
        SimpleMatrix in2 = rewriteInto3ColumnMatrix(unMatchedRight);
        
        return evaluateFitAndMatchPoints(fm, in1, in2, tolerance, unMatchedLeft, 
            unMatchedRight);
    }
    
    private StereoProjectionTransformerFit evaluateFitAndMatchPoints(
        SimpleMatrix fm, 
        SimpleMatrix unMatchedLeft, SimpleMatrix unMatchedRight, 
        double tolerance,
        PairIntArray outputMatchedLeft, PairIntArray outputMatchedRight) {
        
        PairIntArray matchedIndexes = new PairIntArray();
        
        SimpleMatrix theRightEpipolarLines = fm.mult(unMatchedLeft);
        
        SimpleMatrix theLeftEpipolarLines = fm.transpose().mult(unMatchedRight);
        
        StereoProjectionTransformerFit fit =
            evaluateFitForUnmatchedOptimal(fm, theRightEpipolarLines,
            theLeftEpipolarLines, unMatchedLeft, unMatchedRight, tolerance,
            matchedIndexes);
        
        for (int i = 0; i < matchedIndexes.getN(); i++) {
            int idx1 = matchedIndexes.getX(i);
            int idx2 = matchedIndexes.getY(i);
            outputMatchedLeft.add(
                (int)Math.round(unMatchedLeft.get(0, idx1)),
                (int)Math.round(unMatchedLeft.get(1, idx1)));
            outputMatchedRight.add(
                (int)Math.round(unMatchedRight.get(0, idx2)),
                (int)Math.round(unMatchedRight.get(1, idx2)));
        }
        
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
    
    public double[] getLeftEpipole() {
        return leftEpipole;
    }
    public double[] getRightEpipole() {
        return rightEpipole;
    }
    public SimpleMatrix getEpipolarLinesInRight() {
        return epipolarLinesInRight;
    }
    public SimpleMatrix getEpipolarLinesInLeft() {
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
     
    PairIntArray getEpipolarLine(SimpleMatrix epipolarLines, int imgWidth, 
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
    
    SimpleMatrix calculateEpipolarRightLines(SimpleMatrix points) {
        return fundamentalMatrix.mult(points);
    }
    
    SimpleMatrix calculateEpipolarLeftLines(SimpleMatrix points) {
        return fundamentalMatrix.transpose().mult(points);
    }

    public StereoProjectionTransformerFit evaluateFitForUnmatched(
        SimpleMatrix fm, SimpleMatrix leftPoints, SimpleMatrix rightPoints, 
        double tolerance, PairIntArray matchedIndexes) {
        
        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        
        if (leftPoints == null) {
            throw new IllegalArgumentException("leftPoints cannot be null");
        }
        
        if (rightPoints == null) {
            throw new IllegalArgumentException("rightPoints cannot be null");
        }
        
        SimpleMatrix theRightEpipolarLines = fm.mult(leftPoints);
       
        SimpleMatrix theLeftEpipolarLines = fm.transpose().mult(rightPoints);
        
        //2D point (x,y) and line (a, b, c): dist=(a*x + b*y + c)/sqrt(a^2 + b^2)
        
        /*
        for each left point, it's match in the right must be one of the closest
        to the right epipolar line and vice versa.
        
        one needs a bipartite matching (optimal matching) for each right
        epipolar line to the best reciprocating point within a certain distance,
        else no match.
        
        leftxy[0] has re[0] which finds rightxy[10],rightxy[12], rightxy[51]
           within tolerance.
           rightxy[10] has le[10] which > tolerance from leftxy[0].
           rightxy[12] has le[12] which < tolerance from leftxy[0].
           rightxy[51] has le[51] which > tolerance from leftxy[0].
           ==> so leftxy[0] would be matched with rightxy[12]
           and put into a left skip list to avoid repeating and
           put into a right skip list.
        
        this approach is however, a greedy match.  a point visited after 
        leftxy[0] might be the correct match to rightxy[12] and leftxy[0]
        might be noise.
        
        so an optimal match is possibly needed for some datasets.  the runtime
        complexity for an optimal match is large, so it should be avoided
        when possible (~O(N^4)).
       
        so, the distance stored for the match leftxy[i] to rightxy[j] would
        be the 2 added and this
        */
        
        return evaluateFitForUnmatchedOptimal(fm, theRightEpipolarLines,
            theLeftEpipolarLines, leftPoints, rightPoints, tolerance,
            matchedIndexes);
    }
  
    /**
     * 
     * @param fm
     * @param rightEpipolarLines
     * @param leftEpipolarLines
     * @param leftPoints
     * @param rightPoints
     * @param tolerance
     * @param matchedIndexes x is index from leftPoints.  y is index from 
     * rightPoints.
     * @return 
     */
    StereoProjectionTransformerFit evaluateFitForUnmatchedOptimal(
        SimpleMatrix fm, 
        SimpleMatrix rightEpipolarLines, SimpleMatrix leftEpipolarLines,
        SimpleMatrix leftPoints, SimpleMatrix rightPoints, double tolerance,
        PairIntArray matchedIndexes) {
        
        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        if (rightEpipolarLines == null) {
            throw new IllegalArgumentException(
            "rightEpipolarLines cannot be null");
        }
        if (leftEpipolarLines == null) {
            throw new IllegalArgumentException(
            "leftEpipolarLines cannot be null");
        }
        if (leftPoints == null) {
            throw new IllegalArgumentException("leftPoints cannot be null");
        }
        if (rightPoints == null) {
            throw new IllegalArgumentException("rightPoints cannot be null");
        }
        
        float[] diffs = matchOptimallyAndCalcResiduals(fm, 
            rightEpipolarLines, leftEpipolarLines,
            leftPoints, rightPoints, tolerance, matchedIndexes);
                    
        float[] avgAndStdDev = MiscMath.getAvgAndStDev(diffs);
        
        StereoProjectionTransformerFit fit = null;
        if (diffs.length > 0) {
            fit = new StereoProjectionTransformerFit(fm,
                diffs.length, tolerance, 
                avgAndStdDev[0], avgAndStdDev[1]);
        } else {
            fit = new StereoProjectionTransformerFit(fm,
                diffs.length, tolerance, 
                Double.MAX_VALUE, Double.MAX_VALUE);
        }
        
        int nMaxMatchable = leftPoints.numCols();
        if (rightPoints.numCols() < nMaxMatchable) {
            nMaxMatchable = rightPoints.numCols();
        }
        
        fit.setNMaxMatchable(nMaxMatchable);
        
        return fit;
    }
    
    float[] matchOptimallyAndCalcResiduals(
        SimpleMatrix fm, 
        SimpleMatrix rightEpipolarLines, SimpleMatrix leftEpipolarLines,
        SimpleMatrix leftPoints, SimpleMatrix rightPoints, double tolerance,
        PairIntArray outputMatchedIndexes) {
        
        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        if (rightEpipolarLines == null) {
            throw new IllegalArgumentException(
            "rightEpipolarLines cannot be null");
        }
        if (leftEpipolarLines == null) {
            throw new IllegalArgumentException(
            "leftEpipolarLines cannot be null");
        }
        if (leftPoints == null) {
            throw new IllegalArgumentException("leftPoints cannot be null");
        }
        if (rightPoints == null) {
            throw new IllegalArgumentException("rightPoints cannot be null");
        }
        
        //2D point (x,y) and line (a, b, c): dist=(a*x + b*y + c)/sqrt(a^2 + b^2)
        
        /*
        for each left point, it's match in the right must be one of the closest
        to the right epipolar line and vice versa.
        
        one needs a bipartite matching (optimal matching) for each right
        epipolar line to the best reciprocating point within a certain distance,
        else no match.
        
        leftxy[0] has re[0] which finds rightxy[10],rightxy[12], rightxy[51]
           within tolerance.
           rightxy[10] has le[10] which > tolerance from leftxy[0].
           rightxy[12] has le[12] which < tolerance from leftxy[0].
           rightxy[51] has le[51] which > tolerance from leftxy[0].
           ==> so leftxy[0] would be matched with rightxy[12]
           and put into a left skip list to avoid repeating and
           put into a right skip list.
        
        this approach is however, a greedy match.  a point visited after 
        leftxy[0] might be the correct match to rightxy[12] and leftxy[0]
        might be noise.
        
        so an optimal match is possibly needed for some datasets.  the runtime
        complexity for an optimal match is large, so it should be avoided
        when possible (~O(N^4)).
       
        so, the distance stored for the match leftxy[i] to rightxy[j] would
        be the 2 added and this
        */
                
        int nPoints1 = rightEpipolarLines.numCols();
        int nPoints2 = leftEpipolarLines.numCols();
        
        if ((nPoints1 == 0) || (nPoints2 == 0)) {
            return new float[0];
        }
        
        float[][] diffsAsCost = new float[nPoints1][nPoints2];
        // the algorithm modifies diffsAsCost, so make a copy
        float[][] diffsAsCostCopy = new float[nPoints1][nPoints2];
        
        for (int i = 0; i < leftPoints.numCols(); i++) {

            diffsAsCost[i] = new float[nPoints2];
            diffsAsCostCopy[i] = new float[nPoints2];
        
            double a = rightEpipolarLines.get(0, i);
            double b = rightEpipolarLines.get(1, i);
            double c = rightEpipolarLines.get(2, i);
            
            double aplusb = Math.sqrt((a*a) + (b*b));
            
            double xL = leftPoints.get(0, i);
            double yL = leftPoints.get(1, i);
        
            //dist = (a*x + b*y + c)/sqrt(a^2 + b^2)
            
            for (int j = 0; j < rightPoints.numCols(); j++) {
                
                double x = rightPoints.get(0, j);
                double y = rightPoints.get(1, j);
                
                double d = (a*x + b*y + c)/aplusb;
                
                // find the reverse distance by projection:
                double aRev = leftEpipolarLines.get(0, j);
                double bRev = leftEpipolarLines.get(1, j);
                double cRev = leftEpipolarLines.get(2, j);
                
                double dRev = (aRev*xL + bRev*yL + cRev)/
                    Math.sqrt((aRev*aRev + bRev*bRev));
                
                float dist = (float)Math.sqrt(d*d + dRev*dRev);

                if (dist > 1.414*tolerance) {
                    dist = Float.MAX_VALUE;
                }
                diffsAsCost[i][j] = dist;
                diffsAsCostCopy[i][j] = dist;
            }
        }
        
        boolean transposed = false;
        if (nPoints1 > nPoints2) {
            diffsAsCostCopy = MatrixUtil.transpose(diffsAsCostCopy);
            transposed = true;
        }
        
        HungarianAlgorithm b = new HungarianAlgorithm();
        int[][] match = b.computeAssignments(diffsAsCostCopy);
        
        // count the number of matches before tolerance filter
        int count = 0;
        for (int i = 0; i < match.length; i++) {
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            count++;
        }
        
        // at this point, have matches between 1 --> 2 and a tolerance.
        // could return a tolerance filtered
        // reduced match list and the diffs for use with other methods
        // (for example, an invoker with mathing goals could use the
        // match list
        
        // x is index from leftPoints.  y is index from rightPoints.                
        float[] diffs = new float[count];
        
        int nMatched = 0;
        
        for (int i = 0; i < match.length; i++) {
           
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
                        
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }
            
            float dist = diffsAsCost[idx1][idx2];
            
            if (dist > tolerance) {
                continue;
            }
            
            diffs[nMatched] = dist;
            
            outputMatchedIndexes.add(idx1, idx2);
            
            nMatched++;
        }
        
        diffs = Arrays.copyOf(diffs, nMatched);
        
        return diffs;
    }
    
    /**
     * evaluate fit for matched point sets and the fundamental matrix
     * @param fm
     * @param matchedLeftPoints
     * @param matchedRightPoints
     * @param tolerance
     * @return 
     */
    public StereoProjectionTransformerFit evaluateFitForAlreadyMatched(
        SimpleMatrix fm, SimpleMatrix matchedLeftPoints, 
        SimpleMatrix matchedRightPoints, double tolerance) {
        
        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        
        if (matchedLeftPoints == null) {
            throw new IllegalArgumentException(
            "matchedLeftPoints cannot be null");
        }
        
        if (matchedRightPoints == null) {
            throw new IllegalArgumentException(
            "matchedRightPoints cannot be null");
        }
        
        SimpleMatrix theRightEpipolarLines = fm.mult(matchedLeftPoints);
       
        SimpleMatrix theLeftEpipolarLines = fm.transpose().mult(matchedRightPoints);
        
        //2D point (x,y) and line (a, b, c): dist=(a*x + b*y + c)/sqrt(a^2 + b^2)
      
        return evaluateFitForAlreadyMatched(fm, theRightEpipolarLines,
            theLeftEpipolarLines, matchedLeftPoints, matchedRightPoints, 
            tolerance);
    }
    
    /**
     * evaluate fit for already matched point lists
     * @param fm
     * @param rightEpipolarLines
     * @param leftEpipolarLines
     * @param leftPoints
     * @param rightPoints
     * @param tolerance
     * @return 
     */
    StereoProjectionTransformerFit evaluateFitForAlreadyMatched(
        SimpleMatrix fm, 
        SimpleMatrix rightEpipolarLines, SimpleMatrix leftEpipolarLines,
        SimpleMatrix leftPoints, SimpleMatrix rightPoints, double tolerance) {
        
        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        if (rightEpipolarLines == null) {
            throw new IllegalArgumentException(
            "rightEpipolarLines cannot be null");
        }
        if (leftEpipolarLines == null) {
            throw new IllegalArgumentException(
            "leftEpipolarLines cannot be null");
        }
        if (leftPoints == null) {
            throw new IllegalArgumentException("leftPoints cannot be null");
        }
        if (rightPoints == null) {
            throw new IllegalArgumentException("rightPoints cannot be null");
        }
        
        float[] diffs = new float[leftPoints.numCols()];
        
        int nMatched = 0;
        
        List<Integer> inlierIndexes = new ArrayList<Integer>();
        
        for (int i = 0; i < leftPoints.numCols(); i++) {

            double a = rightEpipolarLines.get(0, i);
            double b = rightEpipolarLines.get(1, i);
            double c = rightEpipolarLines.get(2, i);
            
            double aplusb = Math.sqrt((a*a) + (b*b));
            
            double xL = leftPoints.get(0, i);
            double yL = leftPoints.get(1, i);
        
            //dist = (a*x + b*y + c)/sqrt(a^2 + b^2)
            
            double x = rightPoints.get(0, i);
            double y = rightPoints.get(1, i);
                
            double d = (a*x + b*y + c)/aplusb;
                
            // find the reverse distance by projection:
            double aRev = leftEpipolarLines.get(0, i);
            double bRev = leftEpipolarLines.get(1, i);
            double cRev = leftEpipolarLines.get(2, i);

            double dRev = (aRev*xL + bRev*yL + cRev)/
                Math.sqrt((aRev*aRev + bRev*bRev));

            float dist = (float)Math.sqrt(d*d + dRev*dRev);

            if (dist > tolerance) {
                continue;
            }
            
            inlierIndexes.add(Integer.valueOf(i));
            
            diffs[nMatched] = dist;
            
            nMatched++;
        }
        
        diffs = Arrays.copyOf(diffs, nMatched);
        
        float[] avgAndStdDev = MiscMath.getAvgAndStDev(diffs);
        
        StereoProjectionTransformerFit fit = null;
        if (diffs.length > 0) {
            fit = new StereoProjectionTransformerFit(fm, 
                diffs.length, tolerance, 
                avgAndStdDev[0], avgAndStdDev[1]);
        } else {
            fit = new StereoProjectionTransformerFit(fm, 
                diffs.length, tolerance, 
                Double.MAX_VALUE, Double.MAX_VALUE);
        }
        
        fit.setNMaxMatchable(leftPoints.numCols());
        
        fit.setInlierIndexes(inlierIndexes);
        
        return fit;
    }

    public PairIntArray getRightXYInt() {
        return getXYInt(rightXY);
    }
    
    public PairIntArray getLeftXYInt() {                         
        return getXYInt(leftXY);
    }
    
    public int getNumberOfMatches() {
        return leftXY.numCols();
    }
    
    public PairFloatArray getRightXYFloat() {                         
        return getXYFloat(rightXY);
    }
    
    public PairFloatArray getLeftXYFloat() {
        return getXYFloat(leftXY);
    }
    
    private PairIntArray getXYInt(SimpleMatrix leftOrRightXY) {
                
        int nPoints = leftOrRightXY.numCols();
         
        PairIntArray out = new PairIntArray();
        
        for (int i = 0; i < nPoints; i++) {
            
            float xP = (float) leftOrRightXY.get(0, i);
            float yP = (float) leftOrRightXY.get(1, i);
           
            out.add(Math.round(xP), Math.round(yP));
        }
        
        return out;
    }
    
    private PairFloatArray getXYFloat(SimpleMatrix leftOrRightXY) {
                
        int nPoints = leftOrRightXY.numCols();
         
        PairFloatArray out = new PairFloatArray();
        
        for (int i = 0; i < nPoints; i++) {
            
            float xP = (float) leftOrRightXY.get(0, i);
            float yP = (float) leftOrRightXY.get(1, i);
           
            out.add(xP, yP);
        }
        
        return out;
    }
    
    /**
     * calculate the epipolar projection among the given points with the
     * assumption that some of the matches are not true matches.  It uses
     * a RANSAC approach and subsets to find inliers (true, low error matches)
     * and the calculates the fundamental matrix from those.
     * 
     * @param matchedLeftXY
     * @param matchedRightXY
     * @param outputLeftXY
     * @param outputRightXY
     * @return
     * @throws NoSuchAlgorithmException 
     */
    public StereoProjectionTransformerFit calculateEpipolarProjection(
        PairFloatArray matchedLeftXY, PairFloatArray matchedRightXY,
        PairFloatArray outputLeftXY, PairFloatArray outputRightXY) 
        throws NoSuchAlgorithmException {
        
        RANSACSolver ransacSolver = new RANSACSolver();
        
        return ransacSolver.calculateEpipolarProjection(matchedLeftXY, 
            matchedRightXY, outputLeftXY, outputRightXY);
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
        
        SimpleSVD<SimpleMatrix> svd = essentialMatrix.svd();
        
        SimpleMatrix v = svd.getV();
        double val = svd.getU().dot(v);
        if (val < 0) {
            v = v.scale(-1);
        }
       
        SimpleMatrix diag = new SimpleMatrix();
        diag.set(0, 0, 1);
        diag.set(1, 1, 1);
        double[][] dDotV = MatrixUtil.dot(diag, v);
        SimpleMatrix d = new SimpleMatrix(dDotV);
        
        SimpleMatrix essentialMatrix2 = svd.getU().mult(d);
        
        SimpleMatrix w = new SimpleMatrix(3, 3);
        w.set(0, 1, -1);
        w.set(1, 0, 1);
        w.set(2, 2, 1);
        
        SimpleMatrix u = svd.getU();
        
        SimpleMatrix uwv = w.mult(v);
        uwv = u.mult(uwv);
        
        SimpleMatrix uwTransposev = u.mult(w.transpose().mult(v));
        
        SimpleMatrix uLastCol = u.extractVector(false, 2);
        
        SimpleMatrix pt1 = u.mult(uwv).transpose();
        SimpleMatrix pt2 = u.mult(uwTransposev).transpose();
        
        //dot(U,dot(W,V)).T,U[:,2])
        SimpleMatrix p1;//pt1,  uLastCol
        
        //(dot(U,dot(W,V)).T,-U[:,2])
        SimpleMatrix p2;//pt1, -uLastCol
        
        //(dot(U,dot(W.T,V)).T,U[:,2])
        SimpleMatrix p3;//pt2, uLastCol
        
        //(dot(U,dot(W.T,V)).T,-U[:,2])
        SimpleMatrix p4;//pt2, -uLastCol
        
    }*/
}
