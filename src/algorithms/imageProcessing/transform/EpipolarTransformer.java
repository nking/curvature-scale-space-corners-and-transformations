package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.List;
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
 * For the 8-point algorithm:
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
    NOTE:
        SVD(A).U == SVD(AA^T).U == SVD(AA^T).V
        SVD(A).V == SVD(A^TA).V == SVD(A^TA).U

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
     I talso allows use of error calculations such as Sampson's.

 (2) build matrix A with the normalized x,y points

 (3) compute linear least square solution to the least eigenvector of f:
     solve A = U * D * V^T   for A*f = [..x...]*f = 0
     A has rank 8.  F has rank 2.
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
     NOTE: This is best done using normalized coordinates and fundamental matrix,
     after step (4) and before step (5).

For the 7-point algorithm:
references are:
  the the Hartley & Zisserman matlab code vgg_F_from_7pts_2img
from a version of http://www.robots.ox.ac.uk/~vgg/hzbook/code/ which is part
of the supplementary material for their book "Multiple View Geometry in Computer Vision
Second Edition"
   and Hartley, R. I. (1994a). Projective reconstruction and invariants from 
multiple images. PAMI, 16(10):1036–1041
   and Torr, P. H. S. and Murray, D. (1997). 
"The development and comparison of robust meth- ods for estimating the 
fundamental matrix. International Journal of Computer Vision", 24(3):271–300.

 The 7-point algorithm solves for the null space of the fundamental
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

 NOTE: the epipolar lines in right image are the projections of the left image
 points  = F * x1
 The epipolar lines in left image are the projections of the  right image 
 points  = F^T * x2
 </pre>
 
 <pre>
    The plunder-dl scoring can be used for comparison between different models.
    for example, comparing results of the 7-point and 8-point 
    solutions or comparing 7-point projection to 6-point affine, etc.

    plunder-dl is from equation 33 of
    Torr, Zisserman, & Maybank 1996, 
    “Robust Detection of Degenerate Configurations whilst Estimating 
    the Fundamental Matrix"
    https://www.robots.ox.ac.uk/~phst/Papers/CVIU97/m.ps.gz
     EQN 33: PL = DOF + (4*n_o + n_i dimension of model)
                   where n_i = number of inliers
                   n_o = number of outliers
                   DOF = 7 for this solver
    n=7               PL = DOF + 4*n_o + n_i* (model_dimension)
         ni=7, no=0   PL = 7   + 0     + 0 * md
         ni=5, no=2   PL = 7   + 8     + 8 * md
         ni=4, no=3   PL = 7   + 12    + 28 * md
    PLUNDER stands for Pick Least UNDEgenerate Randomly, Description Length

    For nPoints=8, model_dimension = 1.
    for nPoints=7 amd only 1 solution in the cubic constraints, model_dimension=2,
    else for nPoints=7, model_dimension = 3.
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
     * calculate the epipolar projection for a set of 8 or more matched points.
     * NOTE that for best results, the method should be given unit standard
     * normalized coordinates.
     * @param leftXY left image correspondence in format of
     * double array with x points on row 0, y points on row 1, 
     *     and 1's on row 2.  the number of columns is the number of data points.
     * @param rightXY right image correspondence in format of 
     * double array with x points on row 0, y points on row 1, 
     *     and 1's on row 2.  the number of columns is the number of data points.
     * @return
     */
    public DenseMatrix calculateEpipolarProjection(
        DenseMatrix leftXY, DenseMatrix rightXY) {

        if (leftXY.numColumns() != rightXY.numColumns() || leftXY.numRows() != rightXY.numRows()) {
            throw new IllegalArgumentException(
                "leftXY and rightXY must be same size");
        }

        if (leftXY.numColumns() == 7) {
            throw new IllegalArgumentException(
                "for 7 points, use calculateEpipolarProjectionFor7Points");
        }

        if (leftXY.numColumns() < 8) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the algorithms requires 8 or more points.");
        }
        
        double[][] m = createFundamentalMatrix(leftXY, rightXY);
        
        /*
        compute linear least square solution:
            solve A = U * D * V^T   for A*f = [..x...]*f = 0
            A has rank 8.  f has rank 2.

        calculate [U,D,V] from svd(A):
        */
        DenseMatrix aMatrix = new DenseMatrix(m);
        
        System.out.printf("matrix A dimensions = %d x %d\n", m.length, m[0].length);
        
        //aMatrix is m x n  (== nData X 9)
        // U   is  m X m = nData X nData the left singular vectors, **column-wise**
        // S   is  min(m, n) the singular values (stored in descending order)
        // V^T is  n X n = 9x9    the right singular vectors, **row-wise**
        SVD svd = null;
        DenseMatrix vT = null;
        DenseMatrix u = null;
        double[] sDiag = null;
        
        try {
            svd = SVD.factorize(aMatrix);
            vT = svd.getVt();
        } catch (NotConvergedException e) {
        
            double[][] aTa = MatrixUtil.multiply(MatrixUtil.transpose(m), m);
            //SVD(A).U == SVD(AA^T).U == SVD(AA^T).V
            //SVD(A).V == SVD(A^TA).V == SVD(A^TA).U 
            //SVD(A) eigenvalues are the same as sqrt( SVD(AA^T) eigenvalues )
            //    and sqrt( SVD(A^TA) eigenvalues )
            try { 
                svd = SVD.factorize(new DenseMatrix(aTa));
                vT = svd.getVt();
            } catch (NotConvergedException ex) {
                Logger.getLogger(EpipolarTransformer.class.getName()).log(Level.SEVERE, null, ex);
                return null;
            }
        }
        
        /*        
        set f to be the eigenvector associated with the smallest eigenvalue
        (which is the last row of V^T or the last column of V).
         the smallest eigenvalue determines the plane of closest fit.
         */
      
        int n = vT.numRows();
        assert(n == 9);
        assert(vT.numColumns() == 9);

        // dimensions of V are nxn and n=9.  smallet eigenvector is last row of v^T and A
        double[][] ff = new double[3][3];
        for (int i = 0; i < 3; i++) {
            ff[i] = new double[3];
          
            ff[i][0] = vT.get(n - 1, (i * 3) + 0);
            ff[i][1] = vT.get(n - 1, (i * 3) + 1);
            ff[i][2] = vT.get(n - 1, (i * 3) + 2);
        }
        DenseMatrix fMatrix = new DenseMatrix(ff);

        /* make the fundamental matrix have a rank of 2
        by performing svd and reconstruction with the two largest
        singular values.
            [U,D,V] = svd(F,0);
        (a.k.a. dimension reduction.  
        see Chap 11 of book "Mining of Massive Datasets" 
        by Jure Leskovec, Anand Rajaraman, Jeff Ullman
        http://www.mmds.org/)

        From [U,D,V] we create:
            F = U * diag([D(1,1) D(2,2) 0]) * V^T, where V^T is V transposed.
        */
        
        vT = null;
        u = null;
        sDiag = null;
        svd = null;
        try {
            svd = SVD.factorize(fMatrix);
            vT = svd.getVt();
            u = svd.getU();
            sDiag = svd.getS();
        } catch (NotConvergedException e) {
            double[][] aTa = MatrixUtil.multiply(MatrixUtil.transpose(ff), ff);
            double[][] aaT = MatrixUtil.multiply(ff, MatrixUtil.transpose(ff));
            //SVD(A).U == SVD(AA^T).U == SVD(AA^T).V
            //SVD(A).V == SVD(A^TA).V == SVD(A^TA).U 
            //SVD(A) eigenvalues are the same as sqrt( SVD(AA^T) eigenvalues )
            //    and sqrt( SVD(A^TA) eigenvalues )
            try { 
                svd = SVD.factorize(new DenseMatrix(aTa));
                vT = svd.getVt();
                sDiag = svd.getS();
                sDiag[0] = Math.sqrt(sDiag[0]);
                sDiag[1] = Math.sqrt(sDiag[1]);
                
                svd = SVD.factorize(new DenseMatrix(aaT));
                u = svd.getU();
            } catch (NotConvergedException ex) {
                Logger.getLogger(EpipolarTransformer.class.getName()).log(Level.SEVERE, null, ex);
                return null;
            }
        }

        // creates U as 3 x 3 matrix
        //         D as length 3 array
        //         V as 3 x 3 matrix

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

        // dimension reduction in this case zeroes out instead of reducing the
        // sizes of the matrices.  if wanted to reduce the size:
        //   remove the last column of U and the last row of V and
        //   the last item in the diagonal of S
                        
        /*
        multiply the terms:
             F = dot(U, dot(diag(D),V^T))
        */
        DenseMatrix dDotV = MatrixUtil.multiply(d, vT);

        // 3x3 with rank 2
        DenseMatrix fm = MatrixUtil.multiply(u, dDotV);
 
        return fm;
    }

    /*
    the 7-point algorithm.
    references are:
      the the Hartley & Zisserman matlab code vgg_F_from_7pts_2img
    from a version of http://www.robots.ox.ac.uk/~vgg/hzbook/code/ which is part
    of the supplementary material for their book "Multiple View Geometry in Computer Vision
    Second Edition"
       and Hartley, R. I. (1994a). Projective reconstruction and invariants from 
    multiple images. PAMI, 16(10):1036–1041
       and Torr, P. H. S. and Murray, D. (1997). 
    "The development and comparison of robust meth- ods for estimating the 
    fundamental matrix. International Journal of Computer Vision", 24(3):271–300.

    (1) Transform the coordinates using unit normal standaridization.

    (2) build matrix A with the normalized x,y points.
    (3) The homogeneous system AX = 0 : X is the null space of matrix A.
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
     * calculate the epipolar projection for 7 correspondences and filter
     * with a chirality check.  returns a list of the filtered solutions.
     * NOTE that for best results, the method should be given unit standard
     * normalized coordinates.
     * references are:
        (1) the the Hartley & Zisserman matlab code vgg_F_from_7pts_2img
        from a version of http://www.robots.ox.ac.uk/~vgg/hzbook/code/ which is part
        of the supplementary material for their book "Multiple View Geometry in Computer Vision
        Second Edition"
        (2) Section IVa of Hartley, R. I. (1994a). Projective reconstruction and invariants from 
        multiple images. PAMI, 16(10):1036–1041
        (3) Torr, P. H. S. and Murray, D. (1997). 
        "The development and comparison of robust methods for estimating the 
        fundamental matrix. International Journal of Computer Vision", 24(3):271–300.
     * @param leftXY
     * @param rightXY
     * 
     * @return
     */
    public List<DenseMatrix> calculateEpipolarProjectionFor7Points(
        DenseMatrix leftXY, DenseMatrix rightXY) {

        final int nSet = 7;
        
        if (leftXY.numRows() != rightXY.numRows() || leftXY.numColumns() != rightXY.numColumns()) {
            throw new IllegalArgumentException("leftXY and rightXY must have same dimensions");
        }
        if (leftXY.numColumns() != nSet) {
            throw new IllegalArgumentException("leftXY.numColumns must be at least 7");
        }
        if (leftXY.numRows() != 3) {
            throw new IllegalArgumentException("leftXY.numRows must be 3");
        }
        
        double[][] m = createFundamentalMatrix(leftXY, rightXY);

        DenseMatrix aMatrix = new DenseMatrix(m);
        SVD svd = null;
        DenseMatrix vT = null;
        try {
            svd = SVD.factorize(aMatrix);
            vT = svd.getVt();
        } catch (NotConvergedException e) {
            double[][] aTa = MatrixUtil.multiply(MatrixUtil.transpose(m), m);
            //SVD(A).U == SVD(AA^T).U == SVD(AA^T).V
            //SVD(A).V == SVD(A^TA).V == SVD(A^TA).U 
            //SVD(A) eigenvalues are the same as sqrt( SVD(AA^T) eigenvalues )
            //    and sqrt( SVD(A^TA) eigenvalues )
            try { 
                svd = SVD.factorize(new DenseMatrix(aTa));
                vT = svd.getVt();
            } catch (NotConvergedException ex) {
                Logger.getLogger(EpipolarTransformer.class.getName()).log(Level.SEVERE, null, ex);
                return null;
            }
        }

        // U   is 7x7
        // V^T is 9X9
        // A*V = U*s
        
        //for i <=r:
        //    A*v_i = σ*u_i
        //for i >r:
        //    A*v_i = 0 and A^T*u_i = 0

        int n = vT.numRows();
        assert(n == 9);
        assert(svd.getU().numRows() == nSet);
                
        double[][] ff1 = new double[3][3];
        double[][] ff2 = new double[3][3];
        for (int i = 0; i < 3; i++) {

            ff1[i] = new double[3];
            ff2[i] = new double[3];
            
            ff1[i][0] = vT.get(n - 1, (i * 3) + 0);
            ff1[i][1] = vT.get(n - 1, (i * 3) + 1);
            ff1[i][2] = vT.get(n - 1, (i * 3) + 2);
            
            ff2[i][0] = vT.get(n - 2, (i * 3) + 0);
            ff2[i][1] = vT.get(n - 2, (i * 3) + 1);
            ff2[i][2] = vT.get(n - 2, (i * 3) + 2);
        }
        
        DenseMatrix[] solutions = solveFor7Point(ff1, ff2);
        
        EpipolarTransformationFit fit = null;
        
        List<DenseMatrix> validatedFM = new ArrayList<DenseMatrix>();
        
        for (DenseMatrix solution : solutions) {

            // chirality (cheirality) check
            DenseMatrix validated = validateSolution(solution, leftXY, rightXY);
            
            if (validated == null) {
                continue;
            }
            
            validatedFM.add(solution);
        }
        
        return validatedFM;
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
    * 
    * note: alternate spelling is chirality.
    */
    @SuppressWarnings({"unchecked"})
    DenseMatrix validateSolution(DenseMatrix solution, DenseMatrix leftXY,
        DenseMatrix rightXY) {
        
        /*        
        NOTE: vgg_contreps of a 3X1 vector e1 is
            Y = [0      e1(3)  -e1(2)
                -e1(3)  0      e1(1)
                 e1(2) -e1(1)  0];
        (looks like the skew-symmetric of e1)
        NOTE: '.*' is matlab notation to operate on each field
        NOTE: the ' is mathematica syntax for conjugate transpose, a.k.a.
               the Hermitian. it's a matrix with signs reversed for imaginary
                components of complex numbers and then the matrix transposed.
                
        function OK = signs_OK(F,x1,x2)
        [u,s,v] = svd(F'); where F' is the conjugate transpose of F
        e1 = v(:,3); (NLK this is == the svd(F).u
        l1 = vgg_contreps(e1) * x1;
        s = sum( (F*x2) .* l1 );
        OK = all(s>0) | all(s<0);

        (F*x2) .* l1 ==>  (solution * rightXY) .* (testE1 * leftXY)
               
        'sum' is a matlab function to sum for each column

        'all' is a function that returns '1' is all items are non-zero, else
            returns 0
        
        //NLK: without having transposed F:
         // e2^T*F = 0  and F*e1 = 0
         // F^T * x2 = leftEpipolarLines, that is lines in left image 
                       due to epipolar projection of right points
         // F * x1 = rightEpipolarLines, that is lines in right image 
                       due to epipolar projection of left points
         // e1 = last column of U / last item of that column
         //    = the left image position of the epipolar projection of the right camera center
         // e2 = last row of V / last item of that row
         //    = the right image position of the epipolar projection of the left camera center
        */
                
        // row 0 = e1 = right camera center projected into left image reference frame
        // row 1 = e2 = left camera center projected into right image reference frame
        double[][] leftRightEpipoles = calculateEpipoles(solution);

        // 3 columns (x,y,1):
        double[] testE1 = leftRightEpipoles[0];  // NOTE: this is normalized
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
        
        // 3X3 * 3X7 
        //  3 X 7
        DenseMatrix fX2 = MatrixUtil.multiply(solution, rightXY);
        
        // 3 x 7
        double[][] fX2l1 = MatrixUtil.elementwiseMultiplication(
            MatrixUtil.convertToRowMajor(fX2), l1);
        
        // matlab sum function:
        //    If A is a matrix, then sum(A) returns a row vector containing 
        //    the sum of each column.
        double[] sum = new double[fX2l1[0].length];
        for (int row = 0; row < fX2l1.length; ++row){
            for (int col = 0; col < fX2l1[0].length; ++col) {
                sum[col] += fX2l1[row][col];
            }
        }

        //'all' is a function that returns '1' is all items are non-zero, else
        //    returns 0
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
    
    /**
     * Determine the camera matrices of an image pair up to a scene 
     * homography, given their fundamental matrix using algorithm
     * of Hartley and Zisserman 2004.
     * <pre>
     * If x2'*F*x1 = 0 for any pair of image points x1 and x2,
         then the camera matrices of the image pair are 
         P1 = [I 0] (as 3x4 matrix) and P2 = vgg_P_from_F(F), up to a scene homography.
     * </pre>
     * algorithm follows source code for vgg_P_from_F.m 
     * adapted from this site and license:
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

       
        %%%%%%%%%%%%%%%%%%%%%%%%%
     * @param f fundamental matrix in format 3x3
     * @return camera matrices in format 3x4
    */
    @SuppressWarnings({"unchecked"})
    public DenseMatrix pFromF(DenseMatrix f) {
        
        /*
         // e2^T*F = 0  and F*e1 = 0
         // F^T * x2 = leftEpipolarLines, that is lines in left image 
                       due to epipolar projection of right points
         // F * x1 = rightEpipolarLines, that is lines in right image 
                       due to epipolar projection of left points
         // e1 = last column of U / last item of that column
         //    = the left image position of the epipolar projection of the right camera center
         // e2 = last row of V / last item of that row
         //    = the right image position of the epipolar projection of the left camera center
        */
        
        /*        
        NOTE: vgg_contreps of a 3X1 vector e1 is
            Y = [0      e1(3)  -e1(2)
                -e1(3)  0      e1(1)
                 e1(2) -e1(1)  0];
        (looks like the skew-symmetric of e1)
              
        function P = vgg_P_from_F(F)
        [U,S,V] = svd(F);
        e = U(:,3);
        P = [-vgg_contreps(e)*F e];
        return
        */
        SVD svdE;
        DenseMatrix u = null;
        DenseMatrix vT = null;
        double[] sDiag = null;
        try {
            svdE = SVD.factorize(f);
            vT = svdE.getVt();
            u = svdE.getU();
            sDiag = svdE.getS();
        } catch (NotConvergedException e) {
            double[][] fm = MatrixUtil.convertToRowMajor(f);
            double[][] aTa = MatrixUtil.multiply(MatrixUtil.transpose(fm), fm);
            double[][] aaT = MatrixUtil.multiply(fm, MatrixUtil.transpose(fm));
            //SVD(A).U == SVD(AA^T).U == SVD(AA^T).V
            //SVD(A).V == SVD(A^TA).V == SVD(A^TA).U 
            //SVD(A) eigenvalues are the same as sqrt( SVD(AA^T) eigenvalues )
            //    and sqrt( SVD(A^TA) eigenvalues )
            try { 
                svdE = SVD.factorize(new DenseMatrix(aTa));
                vT = svdE.getVt();
                sDiag = svdE.getS();
                sDiag[0] = Math.sqrt(sDiag[0]);
                sDiag[1] = Math.sqrt(sDiag[1]);
                
                svdE = SVD.factorize(new DenseMatrix(aaT));
                u = svdE.getU();
            } catch (NotConvergedException ex) {
                Logger.getLogger(EpipolarTransformer.class.getName()).log(Level.SEVERE, null, ex);
                return null;
            }
        }
        
        assert(u.numColumns() == 3);
        assert(u.numRows() == 3);
        assert(vT.numColumns() == 3);
        assert(vT.numRows() == 3);
        // e1 = last column of U but not normalized by the last item of that column
        double[] e1 = new double[u.numRows()];
        for (int i = 0; i < e1.length; i++) {
            e1[i] = u.get(i, 2);
        }
        
        //vgg_contreps of a 3X1 vector e1 is
        //    Y = [0      e1(3)  -e1(2)
        //        -e1(3)  0      e1(1)
        //         e1(2) -e1(1)  0];
        double[][] contrepsE1 = new double[3][3];
        contrepsE1[0] = new double[]{0, e1[2], -e1[1]};
        contrepsE1[1] = new double[]{-e1[2], 0, e1[0]};
        contrepsE1[2] = new double[]{e1[1], -e1[0], 0};
        
        /*NOTE:  contrepsE1 * e1 = 0
            and  e1^T * contrepsE1 = 0
        */
        
        double[][] P = MatrixUtil.copy(contrepsE1);
        MatrixUtil.multiply(P, -1);
        
        // 3X3
        P = MatrixUtil.multiply(P, MatrixUtil.convertToRowMajor(f));
        
        // append e onto last column of P
        double[][] out = new double[P.length][P[0].length + 1];
        for (int row = 0; row < P.length; ++row){
            out[row] = new double[P[0].length + 1];
            System.arraycopy(P[row], 0, out[row], 0, P[0].length);
            out[row][P[0].length] = e1[row];
        }

        return new DenseMatrix(out);
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
     packages the x, y coordinates into a data structure that is used for 
     normalized coordinates, but does not perform normalization.
     The transformation included as a centroid of 0 and scale of 1.
     
     * @param xy
     * @return
     */
    @SuppressWarnings({"unchecked"})
    NormalizedXY formatForAlreadyNormalized(DenseMatrix xy) {

        double scale = 1;
        double cenX = 0;
        double cenY = 0;   
        DenseMatrix tMatrix = createScaleTranslationMatrix(scale, cenX, cenY);

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
        
        NormalizedXY normalizedXY = new NormalizedXY();
        normalizedXY.setCentroidXY(new double[]{cenX, cenY});
        normalizedXY.setNormMatrix(tMatrix);
        normalizedXY.setXy(xy);

        return normalizedXY;
    }

    /**
     normalize the x,y coordinates as recommended by Hartley 1997 and return
     the matrix and coordinates.
     does not modify the state of this transformer instance.

     the normalized coordinates have an origin = centroid of
     the points.
     the mean square distance of the normalized points from the
     origin is sqrt(2) pixels.


     * @param xy
     * @return
     */
    @SuppressWarnings({"unchecked"})
    public static NormalizedXY normalize(DenseMatrix xy) {

        /*
        format points such that the applied translation
        and scaling have the effect of:

        a) points are translated so that their centroid is at the origin.
        b) points are then scaled so that the average distance from the
           origin is sqrt(2)
        c) the transformation is applied to each of the 2 images separately.
        */

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

        double stDev = 0;

        // using a euclidean distance, chosen for expected use on images with square pixels
        for (int i = 0; i < n; i++) {
            double diffX = xy.get(0, i) - cen0;
            double diffY = xy.get(1, i) - cen1;
            double dist = (diffX * diffX) + (diffY * diffY);
            stDev += dist;
        }

        stDev = Math.sqrt(stDev/(n - 1.));

        //stDev * factor = sqrt(2)
        double scaleFactor = Math.sqrt(2)/stDev;
        
        DenseMatrix tMatrix = createScaleTranslationMatrix(scaleFactor, cen0, cen1);
        //DenseMatrix tMatrix = createScaleTranslationMatrix(1, 0, 0);

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
                        
        DenseMatrix normXY = new DenseMatrix(MatrixUtil.multiply(tMatrix, xy));

        NormalizedXY normalizedXY = new NormalizedXY();
        normalizedXY.setCentroidXY(new double[]{cen0, cen1});
        normalizedXY.setNormMatrix(tMatrix);
        normalizedXY.setXy(normXY);

        return normalizedXY;
    }

    /**
     * create a matrix to be applied on the left side of the dot operator
     * with a matrix of points to transform the points by scale and translation.
     * @param scale (this is actually sqrt(2)/stDev)
     * @param centroidX
     * @param centroidY
     * @return
     */
    protected static DenseMatrix createScaleTranslationMatrix(double scale,
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
     * 
     * @param fm the normalized fundamental matrix, that is the solution for
     * the fundamental metrix using the normalized correspondence.
     * @param tMatrix1 the transformation matrix used to perform normalization
     * on the left correspondence.  The matrix includes terms for the x and y 
     * centroids and scalings.
     * @param tMatrix2 the transformation matrix used to perform normalization
     * on the right correspondence.  The matrix includes terms for the x and y 
     * centroids and scalings.
     * @return 
     */
    public static DenseMatrix denormalizeTheFundamentalMatrix(DenseMatrix fm,
        DenseMatrix tMatrix1, DenseMatrix tMatrix2) {
        
        //denormalize:  F = (T_1)^T * F * T_2
        //    T_1 is normalizedXY1.getNormalizationMatrix();
        //    T2 is normalizedXY2.getNormalizationMatrix();
        DenseMatrix t1Transpose = MatrixUtil.transpose(tMatrix1);
        DenseMatrix t2 = tMatrix2;

        DenseMatrix denormFundamentalMatrix = MatrixUtil.multiply(t1Transpose,
            MatrixUtil.multiply(fm, t2));

        double s = 1. / (denormFundamentalMatrix.get(2, 2) + eps);
        MatrixUtil.multiply(denormFundamentalMatrix, s);

        denormFundamentalMatrix = (DenseMatrix) denormFundamentalMatrix.transpose();
        
        return denormFundamentalMatrix;
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
     The left epipole is e1 = the left image position of the epipolar projection 
     of the right camera center.  It is the last column of U / last item of that 
     column.
     The right epipole e2 = the right image position of the epipolar projection 
     of the left camera center.  It is the last row of V / last item of that row.
     * @param fundamentalMatrix
     * @return a matrix holding whose first row holds e1 and 2nd row holds e2.
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

        NOTE:
            SVD(A).U == SVD(AA^T).U == SVD(AA^T).V
            SVD(A).V == SVD(A^TA).V == SVD(A^TA).U
        
        coords of epipole e’ (==e2) w.r.t. left image coords (==img1) is where 
        the right camera is w.r.t. left image coords=
            svd(F^T*F).u[2]/svd(F^T*F).u[2][2]
        */
                
        SVD svdE;
        DenseMatrix u = null;
        DenseMatrix vT = null;
        double[] sDiag = null;
        try {
            svdE = SVD.factorize(fundamentalMatrix);
            vT = svdE.getVt();
            u = svdE.getU();
            sDiag = svdE.getS();
        } catch (NotConvergedException e) {
            double[][] fm = MatrixUtil.convertToRowMajor(fundamentalMatrix);
            double[][] aTa = MatrixUtil.multiply(MatrixUtil.transpose(fm), fm);
            double[][] aaT = MatrixUtil.multiply(fm, MatrixUtil.transpose(fm));
            //SVD(A).U == SVD(AA^T).U == SVD(AA^T).V
            //SVD(A).V == SVD(A^TA).V == SVD(A^TA).U 
            //SVD(A) eigenvalues are the same as sqrt( SVD(AA^T) eigenvalues )
            //    and sqrt( SVD(A^TA) eigenvalues )
            try { 
                svdE = SVD.factorize(new DenseMatrix(aTa));
                vT = svdE.getVt();
                sDiag = svdE.getS();
                sDiag[0] = Math.sqrt(sDiag[0]);
                sDiag[1] = Math.sqrt(sDiag[1]);
                
                svdE = SVD.factorize(new DenseMatrix(aaT));
                u = svdE.getU();
            } catch (NotConvergedException ex) {
                Logger.getLogger(EpipolarTransformer.class.getName()).log(Level.SEVERE, null, ex);
                return null;
            }
        }
        
        /*
            The left epipole is e1 = last column of U / last item of that column
         It is  the left image position of the epipolar projection of the right camera center
         The right epipole e2 = last row of V / last item of that row
         It is the right image position of the epipolar projection of the left camera center
        */
        
        assert(u.numColumns() == 3);
        assert(u.numRows() == 3);
        assert(vT.numColumns() == 3);
        assert(vT.numRows() == 3);
        // e1 = last column of U / last item of that column
        // e2 = last row of V / last item of that row
        double[] e1 = new double[u.numRows()];
        double e1Div = u.get(2, 2);
        for (int i = 0; i < e1.length; i++) {
            e1[i] = u.get(i, 2)/e1Div;
        }
        
        double[] e2 = new double[vT.numColumns()];
        double e2Div = vT.get(2, 2);
        for (int i = 0; i < e2.length; i++) {
            e2[i] = vT.get(2, i)/e2Div;
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
    
    public static DenseMatrix extractIndices(DenseMatrix m, List<Integer> inlierIndexes) {
        DenseMatrix out = new DenseMatrix(m.numRows(), inlierIndexes.size());
        int r = 0;
        for (int i = 0; i < inlierIndexes.size(); ++i) {
            int idx = inlierIndexes.get(i);
            for (int j = 0; j < m.numRows(); ++j) {
                out.add(j, r, m.get(j, idx));
            }
            r++;
        }
        return out;
    }

    /**
     * if assume gaussian errors and chi-squared statistics, for 1 degree of 
     * freedom (i.e. fitting a line, fundamental matrix, d^2 = 3.84*(st.dev^2)
     * @param standardDeviation
     * @return 
     */
    public static double estimateToleranceForDOF1(double standardDeviation) {
        double d = Math.sqrt(3.84*standardDeviation*standardDeviation);
        return d;
    }
    
    /**
     * if assume gaussian errors and chi-squared statistics, for 2 degree of freedom (i.e. fitting
     * a line, fundamental matrix, d^2 = 5.99*(st.dev^2)
     * @param standardDeviation
     * @return 
     */
    public static double estimateToleranceForDOF2(double standardDeviation) {
        double d = Math.sqrt(5.99*standardDeviation*standardDeviation);
        return d;
    }
    /**
     * if assume gaussian errors and chi-squared statistics, for 2 degree of freedom (i.e. fitting
     * a line, fundamental matrix, d^2 = 7.82*(st.dev^2)
     * @param standardDeviation
     * @return 
     */
    public static double estimateToleranceForDOF3(double standardDeviation) {
        double d = Math.sqrt(7.82*standardDeviation*standardDeviation);
        return d;
    }

    /**
    calculate the 4 possible projection matrices from the essential matrix.
    * Note that the essential matrix is the transformation matrix between points
    */
    /*
    public DenseMatrix[] calculatePFromEssential(DenseMatrix essentialMatrix) {
    }
    */

    private String _toString(DenseMatrix a, String decimalFormat) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < a.numRows(); ++i) {
            for (int j = 0; j < a.numColumns(); ++j) {
                sb.append(String.format(decimalFormat, a.get(i, j)));
                if (j < (a.numColumns() - 1)) {
                    sb.append(", ");
                }
            }
            sb.append("\n");
        }
        return sb.toString();
    }
}
