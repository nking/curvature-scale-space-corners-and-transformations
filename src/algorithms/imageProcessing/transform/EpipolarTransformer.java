package algorithms.imageProcessing.transform;

import algorithms.SubsetChooser;
import algorithms.misc.PolynomialRootSolver;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath;
import algorithms.misc.MiscMath0;
import algorithms.statistics.Standardization;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;
import algorithms.matrix.MatrixUtil.SVDProducts;
import algorithms.util.FormatArray;

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
    u^T*v represents the inner product (= dot product).  result is a scalar.
    u*v^T represents the outer product. result is a matrix.
    u_1 is the (x,y) points from image 1 and u_2 are the matched (x,y) points
        from image 2.

     Epipolar Plane:
        An epipolar plane is defined by the principal points o1, and o2 of the cameras and the point p.
     Epipoles:
         The projection of one camera center onto the other camera's image plane is an epipole.
         (the epipole can lie outside the image boundaries)
         in the Essential Matrix: e2 ~ T and e1 ~ R*T up to a scalar factor.
         in general: e1 is the left null space of E and e2 is the right null space of E.
         e2^T*E = 0
         E*e1 = 0
         and e1 = SVD(E).u <-- last column, then divided by last item
         e2 = SVD(E).v <-- last row, then divided by last item
     Epipolar line of p:
         intersection of the epipolar place of p with an image place is the epipolar line of p.
         l2 ~ E*x1
         l1 ~ E^T*x2

         l*e=0 for (l1,e1) and also (l2,e2)
         l^T*x=0 for (l1,x1) and also (l2,x2)

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
             A is nData X 9
    and f is a nine-vector (size 9X1) containing the entries of the matrix F
 
    To avoid the trivial scale, ||f|| = 1 where f is the norm of f

    the rank of A is 8, but noise and other errors make it 9.
    if A is longer than 8, the system is over-determined (over specified) and
    so must be solved using least squares.  the set may be over determined
    and not have a zero solution.

    we want the vector f that minimizes ||A*f|| subject to the constraint
    that ||f|| = f^T*f = 1

    the solution is the unit eigenvector corresponding to the smallest
    eigenvalue of A^T*A.

    Since A^T*A is positive semi-definite and symmetric, all of its eigenvectors
    are real and positive or zero.
    This eigenvector is what he calls the least eigenvector of A^T*A and
    it is found via the Jacobi algorithm or Singular Value Decomposition.
    NOTE:
        SVD(A).U == SVD(A^T).V == SVD(AA^T).U == SVD(AA^T).V
        SVD(A).V == SVD(A^T).U == SVD(A^TA).V == SVD(A^TA).U

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

         F = (T_2)^T * F * T_1

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

  note:
    because translation is not a linear transformation (see Strang Chap 7)
       one has to keep it as a separate transformation matrix when
       performing operations on a sequence of matrices such as inverse
       and transpose operations.

    translation matrix: inverse changes the signs of the translation elements,
        but not the diagonal

    rotation matrix: inverse is the transpose of rotation matrix.
     
    scaling matrix: inverse is performed on each element, that is, the reciprocal.

    (A*B*C)^-1 = (C^-1) * (B^-1) * (A^-1)

    also, when A * A^(-1) = I, one can use:
                    1
        A^(-1) =  ------ C^(T)  where C_ij = cofactor of a_ij
                   det A

  NOTE: the normalization suggested by Hartley is explored further in
  Chojnacki et al. 2003,  "Revisiting Hartley’s Normalized Eight-Point Algorithm"
  
  Some excerpts:
   
     xc is the centroid of x coordinates.
     yc is the centroid of y coordinates.
     s is the root mean square distance of (x-xc,y-yc) to origin divided by sqrt(2)
    
     u1_normalized = T1 * u1 
     u2_normalized = T2 * u2 

     denormalized F = T2^T * F_normalized * T1

     denormalized u1 = T1^-1 * u1_normalized and similar foru2

     FM_normalized = inverse(transpose(T2)) * FM * inverse(T1)

     denormalized FM = transpose(T2) * FM_normalized * T1 

     u2^T * FM * u1 = u2_normalized^T * FM_normalized * u1_normalized = residual
  
              | 1/s   0  0 |   | 1  0  -xc |   | 1/s    0   -xc/s |
          T = |  0  1/s  0 | * | 0  1  -yc | = |   0  1/s   -yc/s |
              |  0    0  1 |   | 0  0   1  |   |   0    0      1  |

                 | 1  0  xc |   | s  0   0 |   | s   0  xc |
          T^-1 = | 0  1  yc | * | 0  s   0 | = | 0   s  yc |
                 | 0  0   1 |   | 0  0   1 |   | 0   0   1 |
   
                        | 1  0  xc |   | s^2  0    0 |   | 1  0  0 |
          T^-1 * T^-T = | 0  1  yc | * | 0   s^2   0 | * | 0  1  0 |
                        | 0  0   1 |   | 0    0    1 |   | xc yc 1 |

                        | s^2 + xc^2   xc*yc       xc |
                      = | yc*xc        s^2 + yc^2  yc |
                        | xc           yc          1  |
        
                                  | s  0   0 |   | 1  0  0 |   |  s   0   0 |
        from that can see  T^-T = | 0  s   0 | * | 0  1  0 | = |  0   s   0 |
                                  | 0  0   1 |   | xc yc 1 |   | xc  yc   1 |
        
              |  1    0    0 |   | 1/s   0  0 |   |   1/s    0    0 |
        T^T = |  0    1    0 | * |  0  1/s  0 | = |     0   1/s   0 |
              |-xc  -yc    1 |   |  0    0  1 |   | -xc/s  -yc/s  1 |
        
                   ( | 1/s   0  0 | )^-1     |  1    0    0 |   |  s  0   0 |
        (T^T)^-1 = ( |  0  1/s  0 | )     *  |  0    1    0 | = |  0  s   0 |
                   ( |  0    0  1 | )        | xc   yc    1 |   | xc  yc  1 |
        
        can see that (T^-1)^T = (T^T)^-1        
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
 <pre>
 A summary of epipolar geometry is in chapter 5 of 
 Ma, Soatto, Kosecká, and Sastry 2012, "An Invitation to 3-D Vision".
 
    e2 when normalized by 3rd coord is in coord space of left image and
       it is the location of the right camera center.
    e2 is the last column of svd.u
    e2 is the left null space of F
    (e2^T*F = 0  e2^T*E = 0)
      ==> e2~T  where T is translation and ~ is equality up to a scale factor

    e1 is the last row of svd.vt
    e1 is the right null space of F
    (F*e1 = 0  E*e1 = 0)
      ==> e1~R^T*T  where R is rotation and T is translation
    l2 = E*x1
    l1 = E^T*x2
    (x1^T*l1=0 and l1^T*e1=0)
    (x2^T*l2=0 and l2^T*e2=0)
 </pre>
 
 NOTE:
For "7-point" correspondences, consider implementing MLESAC.
     "MLESAC: A new robust estimator with application to estimating image geometry"
     by P. H. S. Torr and A. Zisserman
     1996 http://www.robots.ox.ac.uk/~vgg/publications/papers/torr00.pdf
 
 TODO: implement Nister's 5-point algorithm to determine the essential matrix 
 * for the case where arguments are the camera intrinsic parameters and the correspondence.
 * @author nichole
 */
public class EpipolarTransformer {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    private static double eps = 1e-12;

    /**
     * calculate the epipolar projection.  all normalizations and de-normalizations are handled internally.
     * if the camera calibration is known, the essential matrix is returned.
     * @param leftXY the left image coordinates of the feature correspondences
     * @param rightXY the right image coordinates of the feature correspondences
     * @param calibrated whether or not the camera calibration is known.  if true, the essential matrix is returned.
     * @return the fundamental matrix (or essential matrix) of the epipolar projection, else null if SVD failed.
     */
    public double[][] calculateEpipolarProjection2(double[][] leftXY, double[][] rightXY,
                                                   boolean calibrated) throws IOException {

        if (leftXY.length != 3) {
            throw new IllegalArgumentException("leftXY.length must be 3");
        }
        if (rightXY.length != 3) {
            throw new IllegalArgumentException("rightXY.length must be 3");
        }
        int n = leftXY[0].length;
        if (rightXY[0].length != n) {
            throw new IllegalArgumentException("leftXY and rightXY must have same dimensions");
        }
        if (n < 7) {
            throw new IllegalArgumentException("not enough points.  currently, this method only calculates the epipolar projection" +
                    "matrix for 7 or more points");
        }

        if (n == 7) {
            return calculateEpipolarProjectionUsing7Points(leftXY, rightXY);
        }

        //leftXY = MatrixUtil.copy(leftXY);
        //rightXY = MatrixUtil.copy(rightXY);

        double[][] a = createFundamentalMatrix(leftXY, rightXY);

        SVDProducts svd0 = null;
        try {
            // we only need V in first SVD operation, so can let the method use a or a^T*a
            svd0 = MatrixUtil.performSVD(a);
        } catch (NotConvergedException ex) {
            //Logger.getLogger(EpipolarTransformer.class.getName()).log(Level.SEVERE, ex.getMessage(), ex);
            return null;
        }
        int ns = 9;
        double[][] vT = svd0.vT;
        assert(vT[0].length == ns);
        double[][] ff = new double[3][3];
        int i;
        for (i = 0; i < 3; i++) {
            ff[i] = new double[3];
            ff[i][0] = vT[ns - 1][(i * 3) + 0];
            ff[i][1] = vT[ns - 1][(i * 3) + 1];
            ff[i][2] = vT[ns - 1][(i * 3) + 2];
        }

        // enforce rank=2
        SVD svd = null;
        try {
            // we need to use "U" and "V" below
            svd = SVD.factorize(new DenseMatrix(ff));
        } catch (NotConvergedException e) {
            Logger.getLogger(EpipolarTransformer.class.getName()).log(Level.SEVERE, null, e);
            return null;
        }
        // keep the largest 2 values in sDiag to make the diagonal rank 2
        double[][] d = MatrixUtil.zeros(3,3);
        if (calibrated) {
            //from Serge Belongie lectures from Computer Vision II, CSE 252B, USSD
            // refers to Ma, Soatto, Kosecka, and Sastry 2003
            // "An Invitation to 3D Vision From Images to Geometric Models", p. 114
            // if this is solving the Essential matrix instead of the Fundamental
            //    Matrix, the diagonal is
            // d = [lambda, lambda, 0] where lambda = (svd.s[0] + svd.s[1])/2.
            // MASKS Theorem 5.9
            double sig = (svd.getS()[0] + svd.getS()[1])/2.;
            d[0][0] = sig;
            d[1][1] = sig;
            //NOTE: to normalize E, one can use sig = 1 here. see p. 122 of MASKS chapt 5.
        } else {
            d[0][0] = svd.getS()[0];
            d[1][1] = svd.getS()[1];
        }

        /*
        multiply the terms:
             F = dot(U, dot(diag(D),V^T))
        */
        vT = MatrixUtil.convertToRowMajor(svd.getVt());
        double[][] dDotV = MatrixUtil.multiply(d, vT);

        // 3x3 with rank 2
        double[][] fm = MatrixUtil.multiply(MatrixUtil.convertToRowMajor(svd.getU()), dDotV);

        return fm;
    }

    /**
     * calculate the epipolar projection given 7 pairs of points.  all normalizations and de-normalizations are handled internally.
     * if the camera calibration is known, the essential matrix is returned.
     * @param leftXY the left image coordinates of the feature correspondences
     * @param rightXY the right image coordinates of the feature correspondences
     * @return the fundamental matrix (or essential matrix) of the epipolar projection.
     */
    public double[][] calculateEpipolarProjectionUsing7Points(double[][] leftXY, double[][] rightXY)
    throws IOException {

        if (leftXY.length != 3) {
            throw new IllegalArgumentException("leftXY.length must be 3");
        }
        if (rightXY.length != 3) {
            throw new IllegalArgumentException("rightXY.length must be 3");
        }
        int n = 7;
        if (rightXY[0].length != n || leftXY[0].length != n) {
            throw new IllegalArgumentException("leftXY and rightXY must have 7 points");
        }

        double[][] a = createFundamentalMatrix(leftXY, rightXY);

        SVDProducts svd = null;
        try {
            svd = MatrixUtil.performSVD(a);
        } catch (NotConvergedException ex) {
            Logger.getLogger(EpipolarTransformer.class.getName()).log(Level.SEVERE, null, ex);
            return null;
        }
        int ns = 9;
        assert(svd.vT[0].length == ns);
        // dimensions of V are nxn and n=9.  smallest eigenvector is last row of v^T and A
        double[][] ff1 = new double[3][3];
        double[][] ff2 = new double[3][3];
        int i;
        for (i = 0; i < 3; i++) {
            ff1[i] = new double[3];
            ff2[i] = new double[3];

            ff1[i][0] = svd.vT[ns - 2][(i * 3) + 0];
            ff1[i][1] = svd.vT[ns - 2][(i * 3) + 1];
            ff1[i][2] = svd.vT[ns - 2][(i * 3) + 2];

            ff2[i][0] = svd.vT[ns - 1][(i * 3) + 0];
            ff2[i][1] = svd.vT[ns - 1][(i * 3) + 1];
            ff2[i][2] = svd.vT[ns - 1][(i * 3) + 2];
        }

        //det A = 0 ==> det(α*F1 + (1 − α)*F2) = 0

        /*
        to solve for 'a' as the roots of a polynomial, write out the multiplication
        including that in the determinant in detail, then group by the powers of a to
        calculate the polynomial coefficients.
        see docs/miscNotes/fundamental_cube_roots.txt
         */

        double[] coeffs = calcOrderCoeffs(ff1, ff2);

        //double[] coeffs3 = calculateCubicRootCoefficientsHartley(ff1, ff2);
        //System.out.printf("Hartley's coeffs:\n%s\n", FormatArray.toString(coeffs3, "%.3e"));
        //System.out.printf("coeffs current:\n%s\n", FormatArray.toString(coeffs, "%.3e"));

        double eps = 1E-4;
        double[] roots;
        if (Math.abs(coeffs[3]) < 1E-7) {
            roots = MiscMath.solveCubicRoots(coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
        } else {
            try {
                roots = PolynomialRootSolver.solveForRealUsingMPSolve(coeffs, eps);
            } catch (IOException e) {
                log.severe("PolynomialRootSolver error: " + e.getMessage());
                roots = MiscMath.solveCubicRoots(coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
            }
        }

        // MASKS Appendix 6.A: usually only one of the solved 'a's is consistent with det(F1 + a*F2) = 0

        return filterForZeroDeterminant(ff1, ff2, roots, eps);
    }

    private double[] calculateCubicRootCoefficientsHartley(double[][] ff1, double[][] ff2) {

        // follow Hartley's vgg_singF_from_FF.m to calc the coeffs:
        // which is rewritten by imkaywu.github.io for opencv
        //TODO: fill in those 2 references here
        double[][][] d = new double[2][2][2];
        double[][] dTmp = MatrixUtil.zeros(3, 3);
        for (int i1 = 0; i1 < 2; ++i1) {
            d[i1] = new double[2][2];
            for (int i2 = 0; i2 < 2; ++i2) {
                d[i1][i2] = new double[2];
                for (int i3 = 0; i3 < 2; ++i3) {
                    if (i1 == 0) {
                        dTmp[0] = MatrixUtil.extractColumn(ff1, 0) ;
                    } else {
                        dTmp[0] = MatrixUtil.extractColumn(ff2, 0) ;
                    }
                    if (i2 == 0) {
                        dTmp[1] = MatrixUtil.extractColumn(ff1, 1) ;
                    } else {
                        dTmp[1] = MatrixUtil.extractColumn(ff2, 1) ;
                    }
                    if (i3 == 0) {
                        dTmp[2] = MatrixUtil.extractColumn(ff1, 2) ;
                    } else {
                        dTmp[2] = MatrixUtil.extractColumn(ff2, 2) ;
                    }
                    dTmp = MatrixUtil.transpose(dTmp);
                    d[i1][i2][i3] = MatrixUtil.determinant(dTmp);
                }
            }
        }

        double c3 = -d[1][0][0]+d[0][1][1]+d[0][0][0]+d[1][1][0]+d[1][0][1]-d[0][1][0]-d[0][0][1]-d[1][1][1];
        double c2 = d[0][0][1]-2*d[0][1][1]-2*d[1][0][1]+d[1][0][0]-2*d[1][1][0]+d[0][1][0]+3*d[1][1][1];
        double c1 = d[1][1][0]+d[0][1][1]+d[1][0][1]-3*d[1][1][1];
        double c0 = d[1][1][1];

        return new double[]{c3, c2, c1, c0};
    }

    // det(α*F1 + (1 − α)*F2) = 0
    private double[][] filterForZeroDeterminant(double[][] f1, double[][] f2, double[] roots, double eps) {
        int i;
        int j;
        int c = 0;
        double[][] f;
        double[][] t1;
        double[][] t2;
        double detF;
        double[][] solns = MatrixUtil.zeros(3*roots.length, 3);
        for (i = 0; i < roots.length; ++i) {
            t1 = MatrixUtil.copy(f1);
            t2 = MatrixUtil.copy(f2);
            MatrixUtil.multiply(t1, roots[i]);
            MatrixUtil.multiply(t2, 1. - roots[i]);
            f = MatrixUtil.pointwiseAdd(t1, t2);
            detF = MatrixUtil.determinant(f);
            if (Math.abs(detF) < eps) {
                for (j = 0; j < 3; ++j) {
                    solns[c * 3 + j] = f[j];
                    System.arraycopy(f[j], 0, solns[c * 3 + j], 0, f[j].length);
                }
                ++c;
            }
        }
        // number of solutions is c
        if (c < roots.length) {
            solns = MatrixUtil.copySubMatrix(solns, 0, 3*c - 1, 0, 2);
        }
        return solns;
    }

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
         * @param calibrated if true, solves for the Essential Matrix, else solves for the
         * Fundamental matrix.  The difference is only in the diagonal used for
         * dimension reduction.
         * @return
         */
    public DenseMatrix calculateEpipolarProjection(
        DenseMatrix leftXY, DenseMatrix rightXY, boolean calibrated) {

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
        double[][] aTa = MatrixUtil.createATransposedTimesA(m);
        //DenseMatrix aMatrix = new DenseMatrix(m);
        DenseMatrix aTaMatrix = new DenseMatrix(aTa);

        //System.out.printf("matrix A dimensions = %d x %d\n", m.length, m[0].length);
        
        //aMatrix is m x n  (== nData X 9)
        // U   is  m X m = nData X nData the left singular vectors, **column-wise**
        // S   is  min(m, n) the singular values (stored in descending order)
        // V^T is  n X n = 9x9    the right singular vectors, **row-wise**

        //A^T * A is 9x9

        SVDProducts svd = null;
        try {
            svd = MatrixUtil.performSVD(aTaMatrix);
        } catch (NotConvergedException e) {            
            Logger.getLogger(EpipolarTransformer.class.getName()).log(Level.SEVERE, null, e);
            return null;
        }

        //when condition number is large, the 2 smallest eigenvalues are close to on another
        //and that makes their eigenvectors sensitive to small perturbations.
        double conditionNumber = svd.s[0]/svd.s[svd.s.length-2];
        log.fine("conditionNumber=" + conditionNumber);
        
        /*        
        set f to be the eigenvector associated with the smallest eigenvalue
        (which is the last row of V^T or the last column of V).
         the smallest eigenvalue determines the plane of closest fit.
         */
      
        int n = svd.vT.length;
        assert(n == 9);
        assert(svd.vT[0].length == 9);

        // dimensions of V are nxn and n=9.  smallest eigenvector is last row of v^T and A
        double[][] ff = new double[3][3];
        for (int i = 0; i < 3; i++) {
            ff[i] = new double[3];
          
            ff[i][0] = svd.vT[n - 1][(i * 3) + 0];
            ff[i][1] = svd.vT[n - 1][(i * 3) + 1];
            ff[i][2] = svd.vT[n - 1][(i * 3) + 2];
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
        
        SVD svd2 = null;
        try {
            svd2 = SVD.factorize(fMatrix);
        } catch (NotConvergedException e) {            
            Logger.getLogger(EpipolarTransformer.class.getName()).log(Level.SEVERE, null, e);
            return null;
        }
       
        // creates U as 3 x 3 matrix
        //         D as length 3 array
        //         V as 3 x 3 matrix

        //F = U * diag([D(1,1) D(2,2) 0]) * V^T, where V^T is V transposed.

        // keep the largest 2 values in sDiag to make the diagonal rank 2
        DenseMatrix d = new DenseMatrix(3, 3);
        d = (DenseMatrix)d.zero();
        if (calibrated) {
            //from Serge Belongie lectures from Computer Vision II, CSE 252B, USSD
            // refers to Ma, Soatto, Kosecka, and Sastry 2003
            // "An Invitation to 3D Vision From Images to Geometric Models", p. 114
            // if this is solving the Essential matrix instead of the Fundamental 
            //    Matrix, the diagonal is 
            // d = [lambda, lambda, 0] where lambda = (svd.s[0] + svd.s[1])/2.
            // MASKS Theorem 5.9
            double sig = (svd2.getS()[0] + svd2.getS()[1])/2.;
            d.set(0, 0, sig);
            d.set(1, 1, sig);
            //NOTE: to normalize E, one can use sig = 1 here. see p. 122 of MASKS chapt 5.
        } else {
            if (svd2.getS().length > 0) {
                d.set(0, 0, svd2.getS()[0]);
            }
            if (svd2.getS().length > 1) {
                d.set(1, 1, svd2.getS()[1]);
            }
        }        
        
        // dimension reduction in this case zeroes out instead of reducing the
        // sizes of the matrices.  if wanted to reduce the size:
        //   remove the last column of U and the last row of V and
        //   the last item in the diagonal of S
                        
        /*
        multiply the terms:
             F = dot(U, dot(diag(D),V^T))
        */
        DenseMatrix dDotV = MatrixUtil.multiply(d, svd2.getVt());

        // 3x3 with rank 2
        DenseMatrix fm = MatrixUtil.multiply(svd2.getU(), dDotV);
 
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
     * NOT READY FOR USE YET
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

        // m is [nData X 9]
        double[][] m = createFundamentalMatrix(leftXY, rightXY);

        // [9X9]
        double[][] aTa = MatrixUtil.multiply(MatrixUtil.transpose(m), m);

        // S   is  min(m, n) the singular values (stored in descending order)
        // V^T is  n X n = 9x9    the right singular vectors, **row-wise**

        SVDProducts svd = null;
        try {
            svd = MatrixUtil.performSVD(aTa);
        } catch (NotConvergedException e) {
            Logger.getLogger(EpipolarTransformer.class.getName()).log(Level.SEVERE, null, e);
            return null;
        }

        //when condition number is large, the 2 smallest eigenvalues are close to on another
        //and that makes their eigenvectors sensitive to small perturbations.
        double conditionNumber = svd.s[0]/svd.s[svd.s.length-2];
        log.fine("conditionNumber=" + conditionNumber);

        //for i <=r:
        //    A*v_i = σ*u_i
        //for i >r:
        //    A*v_i = 0 and A^T*u_i = 0

        int n = svd.vT.length;
        assert(n == 9);

        // last singular value being 0, the last column of V is a solution to the nullspace, that is A*x=0

        double[][] ff1 = new double[3][3];
        double[][] ff2 = new double[3][3];
        for (int i = 0; i < 3; i++) {

            ff1[i] = new double[3];
            ff2[i] = new double[3];
            
            ff1[i][0] = svd.vT[n - 2][(i * 3) + 0];
            ff1[i][1] = svd.vT[n - 2][(i * 3) + 1];
            ff1[i][2] = svd.vT[n - 2][(i * 3) + 2];
            
            ff2[i][0] = svd.vT[n - 1][(i * 3) + 0];
            ff2[i][1] = svd.vT[n - 1][(i * 3) + 1];
            ff2[i][2] = svd.vT[n - 1][(i * 3) + 2];
        }
        
        DenseMatrix[] solutions = solveFor7Point(ff1, ff2);

        List<DenseMatrix> validatedFM = new ArrayList<DenseMatrix>();
        for (DenseMatrix solution : solutions) {
            //System.out.printf("validating FM=\n%s\n", FormatArray.toString(solution, "%.3e"));
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
        e1 = v(:,3); <== this is e2 of svd(F)
        l1 = vgg_contreps(e1) * x1;
        s = sum( (F*x2) .* l1 );
        OK = all(s>0) | all(s<0);

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
         // e2 = last row of V / last item of that row
        */

        // 3X3 * 3X7 = 3 X 7
        DenseMatrix fX2 = MatrixUtil.multiply(solution, rightXY);

        // row 0 = e1 = right camera center projected into left image reference frame
        // row 1 = e2 = left camera center projected into right image reference frame
        double[][] leftRightEpipoles = calculateEpipoles(solution);

        // 3 columns (x,y,1):
        //vgg_contreps of a 3X1 vector e1 is
        //    Y = [0      e1(3)  -e1(2)
        //        -e1(3)  0      e1(1)
        //         e1(2) -e1(1)  0];
        //  [x y ] [x0 x1 x2 x3 x4 x5 x6]
        //         [y0 y1 y2 y3 y4 y5 y6]
        double[][] contrepsE2 = MatrixUtil.skewSymmetric(leftRightEpipoles[1]);
        MatrixUtil.multiply(contrepsE2, -1);

        // sum of ( (F * x2) .*  [e2]_x*x1 )
        // =      ( l2 .*  [e2]_x*x1 )

        // from MASKS, have
        //    l2 ~ F*x1
        //    l1 ~ F^T*x2
        //     l^T*x=0 for (l1,x1) and also (l2,x2)
        // and l*e=0 for (l1,e1) and also (l2,e2)

        // 3 X 7
        double[][] l1 = MatrixUtil.multiply(contrepsE2, MatrixUtil.convertToRowMajor(leftXY));

        // 3 x 7
        double[][] fX2l1 = MatrixUtil.pointwiseMultiplication(MatrixUtil.convertToRowMajor(fX2), l1);

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
        EpipolarTransformer tr = new EpipolarTransformer();
        
        SVDProducts svdE = null;
        try {
            svdE = MatrixUtil.performSVD(f);
        } catch (NotConvergedException ex) {
            Logger.getLogger(EpipolarTransformer.class.getName()).log(Level.SEVERE, null, ex);
            return null;
        }
        
        assert(svdE.u[0].length == 3);
        assert(svdE.u.length == 3);
        assert(svdE.vT[0].length == 3);
        assert(svdE.vT.length == 3);
        // e2 = last column of U but not normalized by the last item of that column
        double[] e2 = new double[svdE.u.length];
        for (int i = 0; i < e2.length; i++) {
            e2[i] = svdE.u[i][2];
        }
        
        //-1* vgg_contreps of a 3X1 vector e is the skewSymmetric
        //    Y = [0      e(3)  -e(2)
        //        -e(3)  0      e(1)
        //         e(2) -e(1)  0];
        double[][] contrepsE2 = new double[3][3];
        contrepsE2[0] = new double[]{0, -e2[2], e2[1]};
        contrepsE2[1] = new double[]{e2[2], 0, -e2[0]};
        contrepsE2[2] = new double[]{-e2[1], e2[0], 0};
        
        // 3X3
        double[][] P = MatrixUtil.copy(contrepsE2);
        P = MatrixUtil.multiply(P, MatrixUtil.convertToRowMajor(f));
        
        // append e onto last column of P
        double[][] out = new double[P.length][P[0].length + 1];
        for (int row = 0; row < P.length; ++row){
            out[row] = new double[P[0].length + 1];
            System.arraycopy(P[row], 0, out[row], 0, P[0].length);
            out[row][P[0].length] = e2[row];
        }

        return new DenseMatrix(out);
    }

    DenseMatrix[] solveFor7Point(double[][] ff1, double[][] ff2) {

        //solve for the roots of equation a0 * x^3 + a1 * x^2 + a2 * x + a3 = 0;

        double[] coeffs = calcOrderCoeffs(ff1, ff2);
        double[] roots = null;//MiscMath.solveCubicRoots(coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
        if (Math.abs(coeffs[3]) < 1E-7) {
            roots = MiscMath.solveCubicRoots(coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
        } else {
            try {
                roots = PolynomialRootSolver.solveForRealUsingMPSolve(coeffs, eps);
            } catch (IOException e) {
                log.severe("PolynomialRootSolver error: " + e.getMessage());
                roots = MiscMath.solveCubicRoots(coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
            }
        }

        double[][] m = new double[3][];
        for (int i = 0; i < 3; i++) {
            m[i] = new double[3];
        }

        DenseMatrix[] solutions = new DenseMatrix[roots.length];

        double detF;
        int c = 0;
        for (int i = 0; i < roots.length; i++) {
            //Fi = a(i)*FF{1} + (1-a(i))*FF{2};
            double a = roots[i];
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    m[row][col] = a*ff1[row][col] + (1. - a)*ff2[row][col];
                }
            }
            // MASKS Appendix 6.A: usually only one of the solved 'a's is consistent with det(F1 + a*F2) = 0
            detF = MatrixUtil.determinant(m);
            if (Math.abs(detF) < eps) {
                solutions[c] = new DenseMatrix(m);
                ++c;
            }
        }
        double eps = 1E-4;
        if (c < roots.length) {
            solutions = Arrays.copyOfRange(solutions,0, c);
        }

        return solutions;
    }

    // see docs/miscNotes/fundamental_cube_roots.txt
    private double[] calcOrderCoeffs( double[][] ff1, double[][] ff2) {

        double f100 = ff1[0][0];
        double f110 = ff1[1][0];
        double f120 = ff1[2][0];
        double f101 = ff1[0][1];
        double f111 = ff1[1][1];
        double f121 = ff1[2][1];
        double f102 = ff1[0][2];
        double f112 = ff1[1][2];
        double f122 = ff1[2][2];

        double f200 = ff2[0][0];
        double f210 = ff2[1][0];
        double f220 = ff2[2][0];
        double f201 = ff2[0][1];
        double f211 = ff2[1][1];
        double f221 = ff2[2][1];
        double f202 = ff2[0][2];
        double f212 = ff2[1][2];
        double f222 = ff2[2][2];

        double c3 = f100*f111*f122
                - f100*f111*f222
                - f100*f112*f121
                + f100*f112*f221
                + f100*f121*f212
                - f100*f122*f211
                + f100*f211*f222
                - f100*f212*f221
                - f101*f110*f122
                + f101*f110*f222
                + f101*f112*f120
                - f101*f112*f220
                - f101*f120*f212
                + f101*f122*f210
                - f101*f210*f222
                + f101*f212*f220

                + f102*f110*f121
                - f102*f111*f120
                + f102*f111*f220
                + f102*f120*f211
                - f102*f121*f210
                + f102*f210*f221
                - f102*f211*f220

                - f110*f102*f221
                + f110*f122*f201
                - f110*f202*f121
                + f110*f202*f221
                - f110*f222*f201

                + f111*f120*f202
                + f111*f200*f222
                - f111*f200*f122
                - f111*f202*f220

                - f112*f120*f201
                + f112*f121*f200
                - f112*f221*f200
                + f112*f201*f220

                + f120*f201*f212
                - f120*f202*f211

                + f121*f202*f210
                - f121*f212*f200

                - f122*f201*f210
                + f122*f211*f200

                - f200*f211*f222
                + f200*f212*f221

                + f201*f210*f222
                - f201*f212*f220

                - f202*f210*f221
                + f202*f211*f220;

        double c2 = f100*f111*f222
                - f100*f112*f221
                - f100*f121*f212
                + f100*f122*f211
                - 2 * f100*f211*f222
                + 2 * f100*f212*f221

                - f101*f110*f222
                + f101*f112*f220
                + f101*f120*f212
                - f101*f122*f210
                + 2 * f101*f210*f222
                - 2 * f101*f212*f220

                + f102*f110*f221
                - f102*f111*f220
                + f102*f121*f210
                - f102*f120*f211
                - 2 * f102*f210*f221
                + 2 * f102*f211*f220

                + f110*f121*f202
                - f110*f122*f201
                + 2 * f110*f201*f222
                - 2 * f110*f202*f221

                - f111*f120*f202
                + f111*f122*f200
                - 2 * f111*f200*f222
                + 2 * f111*f202*f220

                + f112*f120*f201
                - f112*f121*f200
                + 2 * f112*f200*f221
                - 2 * f112*f201*f220

                - 2 * f120*f201*f212
                + 2 * f120*f202*f211

                - 2 * f122*f200*f211
                + 2 * f121*f200*f212
                - 2 * f121*f202*f210

                + 2 * f122*f201*f210

                + 3 * f200*f211*f222
                - 3 * f200*f212*f221

                - 3 * f201*f210*f222
                + 3 * f201*f212*f220

                + 3 * f202*f210*f221
                - 3 * f202*f211*f220;

        double c1 = f100*f211*f222
                - f100*f221*f212

                - f101*f210*f222
                + f101*f212*f220

                + f102*f210*f221
                - f102*f211*f220

                + f200*f111*f222
                - f200*f112*f221
                - f200*f121*f212
                + f200*f211*f122

                - 3 * f200*f211*f222

                + 3 * f200*f221*f212

                - f201*f110*f222
                + f201*f112*f220
                + f201*f120*f212
                - f201*f210*f122

                + 3 * f201*f210*f222

                - 3 * f201*f212*f220
                + f202*f110*f221
                - f202*f120*f211

                - f202*f111*f220
                + f202*f121*f210

                - 3 * f202*f210*f221

                + 3 * f202*f211*f220;

        double c0 = f202*f210*f221
                - f201*f210*f222
                - f202*f211*f220
                + f200*f211*f222
                - f200*f212*f221
                + f201*f212*f220;

        return new double[]{c3, c2, c1, c0};
    }

    @SuppressWarnings({"unchecked"})
    DenseMatrix denormalizeTheFundamentalMatrix(
        DenseMatrix normalizedFundamentalMatrix,
        NormalizedXY normalizedLeftXY, NormalizedXY normalizedRightXY) {

        return denormalizeTheFundamentalMatrix(normalizedFundamentalMatrix, 
            normalizedLeftXY.getNormalizationMatrices(), 
            normalizedRightXY.getNormalizationMatrices());
       
    }

    /**
     normalize the x,y coordinates as recommended by Hartley 1997 and return
     the matrix and coordinates.
     does not modify the state of this transformer instance.

     the normalized coordinates have an origin = centroid of
     the points.
     the mean square distance of the normalized points from the
     origin is sqrt(2) pixels.

     normalized xy = T * xy

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

        double scale = 0;
        double diffX, diffY;
        for (int i = 0; i < n; i++) {
            diffX = xy.get(0, i) - cen0;
            diffY = xy.get(1, i) - cen1;
            scale += (diffX*diffX + diffY*diffY);
        }
        //scale: root mean square distance of (x,y) to origin is sqrt(2)
        scale = Math.sqrt(scale/(2.*(n-1)));
        
        NormalizationTransformations tMatrices = createScaleTranslationMatrices(
            scale, cen0, cen1);
       
        DenseMatrix normXY = MatrixUtil.multiply(new DenseMatrix(tMatrices.t), xy);
        
        NormalizedXY normalizedXY = new NormalizedXY();
        normalizedXY.setNormMatrices(tMatrices);
        normalizedXY.setXy(normXY);

        return normalizedXY;
    }
    
   /**
     normalize the x,y coordinates as recommended by Hartley 1997 and return
     the matrix and coordinates.
     does not modify the state of this transformer instance.

     the normalized coordinates have an origin = centroid of
     the points and a scale equal to the standard deviation of the
     mean subtracted coordinates.

     * @param xy
     * @return
     */
    @SuppressWarnings({"unchecked"})
    public static NormalizedXY normalizeUsingUnitStandard(DenseMatrix xy) {

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

        double scale = 0;
        double diffX, diffY;
        for (int i = 0; i < n; i++) {
            diffX = xy.get(0, i) - cen0;
            diffY = xy.get(1, i) - cen1;
            scale += (diffX*diffX + diffY*diffY);
        }
        scale = Math.sqrt(scale/(n-1.));
        
        NormalizationTransformations tMatrices = createScaleTranslationMatrices(
            scale, cen0, cen1);
       
        DenseMatrix normXY = MatrixUtil.multiply(new DenseMatrix(tMatrices.t), xy);
        
        NormalizedXY normalizedXY = new NormalizedXY();
        normalizedXY.setNormMatrices(tMatrices);
        normalizedXY.setXy(normXY);

        return normalizedXY;
    }

    /**
     perform unit standard normalization of the points xy to result in a mean of 0 for xy and a standard deviation of 1,
     and return the mean of x, mean of y, st. dev of x and st dev of y.  note that no normalization is performed on the z coordinates.
     note that the normalization alters the input.
     * @param xy data points in format [3 X n].  e.g. xy[0][0] is x, xy[1][0] is y, xy[2][0] is z for point 0.
     * @return mean of x, mean of y, st. dev of x and st dev of y
     */
    @SuppressWarnings({"unchecked"})
    public static double[] normalizeUsingUnitStandard(double[][] xy) {

        if (xy.length != 3) {
            throw new IllegalArgumentException("xy.length must be 3");
        }
        int n = xy[0].length;

        double[] mS = new double[4];

        int i;
        int j;
        for (i = 0; i < n; ++i) {
            for (j = 0; j < 2; ++j) {
                mS[j] += xy[j][i];
            }
        }

        for (j = 0; j < 2; ++j) {
            mS[j] /= n;
        }

        double d;
        for (i = 0; i < n; ++i) {
            for (j = 0; j < 2; ++j) {
                d = (xy[j][i] - mS[j]);
                mS[j + 2] += (d * d);
            }
        }
        for (j = 2; j < 4; ++j) {
            mS[j] = Math.sqrt(mS[j]/(n - 1.0));
        }

        double[][] t = new double[3][];
        t[0] = new double[]{1./mS[2],       0,     -mS[0]/mS[2]};
        t[1] = new double[]{0,           1./mS[3], -mS[1]/mS[3]};
        t[2] = new double[]{0,           0,           1};

        double[] chk = MatrixUtil.multiplyMatrixByColumnVector(t, MatrixUtil.extractColumn(xy, 0));

        double[][] xyN = MatrixUtil.multiply(t, xy);
        for (i = 0; i < 3; ++i) {
            System.arraycopy(xyN[i], 0, xy[i], 0, xyN[i].length);
        }
        return mS;
    }

    /**
     * create 2 scale and translation matrices to normalize homogeneous 2D coordinates
     * by multiplying on the left side, and to de-normalize the normalized
     * coordinates by multiplying on the left side.
     * The normalization matrix transforms x by (x - centroidX)/scale and similarly transforms y.
     * The de-normalization matrix can be used to multiply by scale and add the centroid.
     * @param scale
     * @param centroidX
     * @param centroidY
     * @return
     */
    protected static NormalizationTransformations createScaleTranslationMatrices(
        double scale,
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

        t * xy = transformed xy

        to reverse the operation, that is de-normalize, need division by scale, 
        then addition of centroid    
        
        Revisiting Hartley’s Normalized Eight-Point Algorithm
        Chojnacki et al. 2003

        inverse changes the order of operations,
        inverse translation matrix: inverse changes the signs of the translation elements,
                but not the diagonal.
        inverse rotation matrix: inverse is the transpose of rotation matrix.
        inverse scaling matrix: inverse is performed on each element, that is, the reciprocal.

        denormalized x2 = transpose(normalized x2) * transpose(T2^1)
        denormalized x1 = (T1^-1) * (normalized x1)
        denormalized FM = transpose(T2) * FM_normalized * T1

                | 1/s   0  0 |   | 1  0  -xc |
            T = |  0  1/s  0 | * | 0  1  -yc |
                |  0    0  1 |   | 0  0   1  |

                   | 1  0  xc |   | s  0   0 | = | s  0  xc |
            T^-1 = | 0  1  yc | * | 0  s   0 |   | 0  s  yc |
                   | 0  0   1 |   | 0  0   1 |   | 0  0  1  |

            transposed(T^-1) = | s   0    0 |
                               | 0   s    0 |
                               | xc  yc   1 |

                          | 1  0  xc |   | s^2  0    0 |   | 1  0  0 |
            T^-1 * T^-T = | 0  1  yc | * | 0   s^2   0 | * | 0  1  0 |
                          | 0  0   1 |   | 0    0    1 |   | xc yc 1 |

                          | s^2 + xc^2   xc*yc       xc |
                        = | yc*xc        s^2 + yc^2  yc |
                          | xc           yc          1  |
        */
        
        double[][] t = new double[3][];
        t[0] = new double[]{1./scale,       0,     -centroidX/scale};
        t[1] = new double[]{0,           1./scale, -centroidY/scale};
        t[2] = new double[]{0,           0,           1};

        double[][] d = new double[3][];
        d[0] = new double[]{scale,       0,     centroidX};
        d[1] = new double[]{0,           scale, centroidY};
        d[2] = new double[]{0,           0,           1};
        
        NormalizationTransformations nt = new NormalizationTransformations();
        nt.t = t;
        nt.tDenorm = d;
        
        return nt;
    }
    
    /**
     * 
     * @param normalizedFundamentalMatrix the normalized fundamental matrix, that is the solution for
     * the fundamental matrix using the normalized correspondence.
     * @param leftNT the transformations created for normalization and de-normalization
     * of the left image coordinates.
     * @param rightNT the transformations created for normalization and de-normalization
     * of the right image coordinates.
     * 
     * @return 
     */
    public static DenseMatrix denormalizeTheFundamentalMatrix(
        DenseMatrix normalizedFundamentalMatrix,
        NormalizationTransformations leftNT, NormalizationTransformations rightNT) {

        //NOTE: denormalized x2 = transpose(normalized x2) * transpose(T2^1)
        //      denormalized x1 = (T1^-1) * (normalized x1)
        //      denormalized FM = transpose(T2) * FM_normalized * T1
        
        DenseMatrix t2T = new DenseMatrix(MatrixUtil.transpose(rightNT.t));

        DenseMatrix t1 = new DenseMatrix(leftNT.t);

        DenseMatrix denormFundamentalMatrix =
            MatrixUtil.multiply(t2T,
                    MatrixUtil.multiply(normalizedFundamentalMatrix, t1));
        
        double factor = 1./(denormFundamentalMatrix.get(2, 2) + eps);
        MatrixUtil.multiply(denormFundamentalMatrix, factor);

        return denormFundamentalMatrix;
    }

    /**
     * a = xLeft ⊗ xRight
     * @param normXY1 a matrix of size 3 x nPoints, where 1st row is x,
     * second is y.
     * @param normXY2 a matrix of size 3 x nPoints, where 1st row is x,
     * second is y.
     * @return fundamental matrix of size [nPoints x 9]
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
        double x1, x2, y1, y2;
        double[][] a = new double[nXY1][9];
        for (int i = 0; i < nXY1; i++) {
            a[i] = new double[9];
            x1 = normXY1.get(0, i);
            x2 = normXY2.get(0, i);
            y1 = normXY1.get(1, i);
            y2 = normXY2.get(1, i);
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
     * a = xLeft ⊗ xRight
     * @param normXY1 a matrix of size 3 x nPoints, where 1st row is x,
     * second is y.
     * @param normXY2 a matrix of size 3 x nPoints, where 1st row is x,
     * second is y.
     * @return fundamental matrix of size [nPoints x 9]
     */
    double[][] createFundamentalMatrix(double[][] normXY1,
                                       double[][] normXY2) {

        if (normXY1 == null) {
            throw new IllegalArgumentException("normXY1 cannot be null");
        }
        if (normXY2 == null) {
            throw new IllegalArgumentException("normXY2 cannot be null");
        }
        int nXY1 = normXY1[0].length;
        if (nXY1 != normXY2[0].length) {
            throw new IllegalArgumentException(
                    "the number of columns in normXY1 != number of cols in normXY2");
        }

        /*
        (2) each row in matrix A:
            x_1*x_2,  x_1*y_2,  x_1,  y_1*x_2,  y_1*y_2,  y_1,  x_2,  y_2,  1
        */
        double x1, x2, y1, y2;
        double[][] a = new double[nXY1][9];
        for (int i = 0; i < nXY1; i++) {
            a[i] = new double[9];
            x1 = normXY1[0][i];
            x2 = normXY2[0][i];
            y1 = normXY1[1][i];
            y2 = normXY2[1][i];
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
     *Epipolar Plane:
     *       An epipolar plane is defined by the principal points o1, and o2 of the cameras and the point p.
     *    Epipoles:
     *       In the epipolar plane, the point where the camera principal point intersects it image is the
     *       epipole.  The location of the epipole in the other image can be thought of as the projection
     *       of the other camera center onto its image.
     *       (the epipole can lie outside the image boundaries)
     *       e1 is the left null space of E and e2 is the right null space of E.
     *         e2^T*E = 0
     *         E*e1 = 0
     *         and e1 = SVD(E).u <-- last column, then divided by last item
     *             e2 = SVD(E).v <-- last row, then divided by last item
     *         for the essential matrix: e2 ~ T and e1 ~ R*T up to a scalar factor.

     * @param fundamentalMatrix
     * @return a matrix holding whose first row holds e1 and 2nd row holds e2
     */
    @SuppressWarnings({"unchecked"})
    double[][] calculateEpipoles(DenseMatrix fundamentalMatrix) {

        /*
        NOTE:
            SVD(A).U == SVD(A^T).V == SVD(AA^T).U == SVD(AA^T).V
            SVD(A).V == SVD(A^T).U == SVD(A^TA).V == SVD(A^TA).U
        
        epipoles:
             [U,D,V] = svd(denormalized FundamentalMatrix);
             e1 = last column of V divided by it's last item
             e2 = last column of U divided by it's last item

        e2 when normalized by 3rd coord is in coord space of left image and
               it is the location of the right camera center.
        */
             
        SVD svdE = null;
        try {
            svdE = SVD.factorize(fundamentalMatrix);
        } catch (NotConvergedException e) {
            Logger.getLogger(EpipolarTransformer.class.getName()).log(Level.SEVERE, null, e);
            return null;
        }

        // e1 = last column of U / last item of that column
        // e2 = last row of V / last item of that row
        double[] e1 = new double[3];
        double e1Div = svdE.getU().get(2, 2);
        for (int i = 0; i < e1.length; i++) {
            e1[i] = svdE.getU().get(i, 2)/e1Div;
        }
        
        double[] e2 = new double[3];
        double e2Div = svdE.getVt().get(2, 2);
        for (int i = 0; i < e2.length; i++) {
            e2[i] = svdE.getVt().get(2, i)/e2Div;
        }

        double[][] e = new double[2][];
        e[0] = e1;
        e[1] = e2;

        return e;
    }
    
    /**
     * convert the fundamental matrix into the essential matrix.
     * @param k1 camera 1 intrinsics matrix.
     * @param k2 camera 2 intrinsics matrix.
     * @param fm the fundamental matrix in the reference frame of the original image space
     * of pixel coordinates.
     * @return the essential matrix in the reference frame of normalized image
     * coordinates (origin is optical center of the image).
     */
    public static double[][] createEssentialFromFundamentalMatrix(
        double[][] k1, double[][] k2, double[][] fm) throws NotConvergedException {
        
        // E = K2^T * F * K1
        //   then SVD(E) to reform E:
        //   E = U * [1 0 0] * V^T
        //           [0 1 0]
        //           [0 0 0]
        double[][] k2T = MatrixUtil.transpose(k2);
        
        double[][] em = MatrixUtil.multiply(k2T, fm);
        em = MatrixUtil.multiply(em, k1);
        
        System.out.printf("E before:\n%s\n", FormatArray.toString(em, "%.3e"));
        
        SVDProducts svd = MatrixUtil.performSVD(new DenseMatrix(em));
        double[][] d = MatrixUtil.createIdentityMatrix(3);
        d[2][2] = 0;
        em = MatrixUtil.multiply(svd.u, d);
        em = MatrixUtil.multiply(em, svd.vT);
        
        System.out.printf("E:\n%s\n", FormatArray.toString(em, "%.3e"));
        
        return em;
    }
    
    /**
     * convert the essential matrix into the fundamental matrix using 
     * F= K2^-T * E * K1^-1
     * @param k1 camera 1 intrinsics matrix.
     * @param k2 camera 2 intrinsics matrix.
     * @param em the essential matrix in the reference frame of normalized image
     * coordinates (origin is optical center of the image)
     * @return the fundamental matrix .
     */
    public static double[][] createFundamentalFromEssentialMatrix(
        double[][] k1, double[][] k2, double[][] em) {
        
        // F= K2^-T * E * K1^-1
        double[][] k1Inv = Camera.createIntrinsicCameraMatrixInverse(k1);
        double[][] k2Inv = Camera.createIntrinsicCameraMatrixInverse(k2);
        double[][] k2InvT = MatrixUtil.transpose(k2Inv);
        
        double[][] fm = MatrixUtil.multiply(k2InvT, em);
        fm = MatrixUtil.multiply(fm, k1Inv);
        
        return fm;
    }
    
    public static class NormalizationTransformations {
        /**
         * a normalization matrix which subtracts a centroid and divides by a scale factor
         * to each of x and y.
         * The matrix is 3X3.  can use as t*xy where xy is size 3XN.
         */
        public double[][] t = null;
        
        /**
         * a de-normalization matrix which multiplies by a scale factor and adds 
         * a centroid to each of x and y
         * to each of x and y.
         * the matrix is 3X3. can use as t*xy where xy is size 3XN.
         */
        public double[][] tDenorm = null;
    } 
    
    public static class NormalizedXY {

        /**
         * 3 dimensional matrix, with column 0 being x, column 1 being y,
         * and the last column is place holder 1's
         */
        private DenseMatrix xy = null;
        
        private NormalizationTransformations normalizationMatrices = null;

        /**
         * @return the factor
         */
        public NormalizationTransformations getNormalizationMatrices() {
            return normalizationMatrices;
        }

        /**
         * @param normMatrices holding the scale and centroid offsets to apply to x, y
         */
        public void setNormMatrices(NormalizationTransformations normMatrices) {
            this.normalizationMatrices = normMatrices;
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
                out.set(j, r, m.get(j, idx));
            }
            r++;
        }
        return out;
    }

    public static double[][] extractIndices(double[][] m, List<Integer> inlierIndexes) {
        double[][] out = MatrixUtil.zeros(m.length, inlierIndexes.size());
        int r = 0;
        int i, j;
        for (i = 0; i < inlierIndexes.size(); ++i) {
            int idx = inlierIndexes.get(i);
            for (j = 0; j < m.length; ++j) {
                out[j][r] = m[j][idx];
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
