package algorithms.imageProcessing.transform;

import algorithms.dimensionReduction.CURDecomposition;
import algorithms.dimensionReduction.CURDecomposition.CUR;
import algorithms.imageProcessing.features.RANSACSolver;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.transform.Camera.CameraParameters;
import algorithms.imageProcessing.transform.Camera.CameraProjection;
import static algorithms.imageProcessing.transform.CameraPose.eps;
import algorithms.matrix.LinearEquations;
import algorithms.matrix.MatrixUtil;
import algorithms.matrix.MatrixUtil.SVDProducts;
import algorithms.util.FormatArray;
import java.util.Arrays;
import java.util.List;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;

/**
 * given correspondence between two images calculate the camera
 * parameters as intrinsic and extrinsic parameters,
 * and the real world position.
 * 
 * Euler rotations:
        
        about z-axis (yaw):           about x-axis (roll):       about the y-axis (pitch):
            | cos φ   -sin φ    0 |    |    1       0       0 |  |  cos ψ    0  sin ψ |
            | sin φ    cos φ    0 |    |    0   cos θ   sin θ |  |      0    1      0 |
            |     0        0    1 |    |    0  -sin θ   cos θ |  | -sin ψ    0  cos ψ |        
        
 * useful reading:
 * <pre>
 * http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
 * add other references here
 * </pre>
 * @author nichole
 */
public class Reconstruction {
    
    /*
    // notes from Serge Belongie lectures from Computer Vision II, CSE 252B, USSD
    
    // homogeneous repr of a point is x_vec = (x, y, 1)^T
    // equation f a line is ell = a*x + b*y + c = 0;
    
    // line rewritten in homogeneous coordinatrs is x_vec^T * ell.
    
    // general conic in 3 dimensions is a*x^2 + b*x*y + c*y^2 + d*x*z + e*y*z + f*z^2 = 0.
    //     rewritten using 2D homogenouse coords, quadratice form: x_vec^T * C * x_vec = 0
    //                  [a   b/2   d/2 ]
    //        where C = [b/2   c   c/2 ]
    //                  [d/2 c/2     f ]
    //        C has 6 variable, 5 DOF, so need 5 points
    //     can then reformat x_vec^T * C * x_vec = 0 into
    //        the "design matrix" * "the carrier vector" = 0
    //               A * c = 0
    //        c = SVD(A).V^T[n-1], the eigenvector assoc w/ smallest eigenvalue.
    
    Parallel lines lines intersect in points at infinity (also known as ideal points) 
    and these points have the form (x, y, 0)^⊤
    the set of ideal points (i.e., points at infinity) is the set of points 
    where parallel lines intersect.
    
    The intersection of two lines is given by their vector cross product.
    
    The line passing through any two points is given by their cross product.
    
    The line at infinity is (0,0,1)^⊤ which can be seen by taking the cross
        product of 2 points at infinity.  e.g. (x1,y1,0)^⊤ cross (x2,y2,0)^⊤ = (0,0,1)^⊤
     
        The two circular points are defined as
           I = (1, i, 0)^⊤
           J = (1, −i, 0)^T
        The circular points lie on l∞, along with all other ideal points. All 
        circles intersect l∞ at points I and J.
    
        Recall the duality between points and lines:
           x2 = H*x1,  l2 = (H^−⊤)*l1
    
        The dual of a conic C is the set of lines satisfying:
           l^⊤ * C∗ * l = 0
           where C∗ is the adjoint in this case, so C∗ ∼ C^−1  
        Dual conics tranform under homography H as:
           C∗′ = H * C∗ * H^⊤
    
        The “conic dual to the circular points” is defined as
           C ∞∗ = I*J^⊤ + J*I^⊤
                ~ [ 1  0  0 ]
                  [ 0  1  0 ]
                  [ 0  0  0 ]
    
    Stratified Reconstruction (cahpters 8 and 6):
    
    Notes on stratified reconstruction in the 2D case (lec 6):
        Each stratum represents a different level of reconstruction we may wish 
        to obtain, namely projective, affine and Euclidean. 
        A general 2D homography can be decomposed into three components:
           H = H_p * H_a * H_e 
        which are the projective, affine and euclidean components.
              
          H = [  I   0 ] * [ K    0 ] * [ R   T ]
              [ v^T  1 ]   [ 0^T  1 ]   [ 0^T 1 ]
                 H_p           H_a         H_e
    
          H_e is a 2D rigid transformation
          H_a is an affine trnsformation
          H_p is a projective transformation known as an “elation.”
               v^T affets the line at infinity, l∞ = (0, 0, 1)^⊤.
               ** Only H_p can map l∞ to a finite line in the image plane, or vice versa.
               suppose the image of l∞ is some line l = (a, b, c)^⊤, then the
               following matrix H will send l back to infinity:
                                [ 1  0  0 ]
                      H = H_a * [ 0  1  0 ]
                                [ a  b  c ]
                          where where H_a is any affine transformation
          
       The key to the affine upgrade is the behavior of the line at infinity. 
       The counterpart to this for the Euclidean upgrade is the behavior of the 
       “circular points.”
    
       Chap 6, Section 3 has more about euclidean upgrades in transformation
          for circular points and conics.
    
    Notes from Lec 8:
        Levels of reconstruction of a scene: projective, affine and euclidean components.
    
        If the camera intrinisc parameters, K, are not known, then only the
        projective reconstruction is possible, but this can be upgraded to 
        affine (parallelism preserved) and Euclidean (parallelism and orthogonality preserved) reconstructions.
    
        Given a set of point correspondences between two views {(x′1,x′2)} 
        we can get the Projective Structure X_p.
        where x' = x * K, that is, x = x' * K^-1.  x' are the image coordinates as pixels.
           first calculate F from the point correspondences, then use F to
           get the projection matrix for camera 2, P2, then triangulate to get
           X_p (which is the 3D point in projective space).
           caveat is that K, R, T from in P2 will not be unique.
           The canonical choice for these two projection matrices is:
               P1 = [ I | 0 ]
               P2 = [ ([T']_x)^T  F  | T' ]
           where T' = K*T and ||T'|| = 1
           and recall that T′ ∼ e2
    
           The triangulation uses DLT and SVD.
    
           for the case where there is no noise:
                X_P is SVD(M).V^T[last row]
           for the case where there is noise, the SVD solution is the initial
                values for a non-linear optimization method.
    
    NOTE: to solve affine reconstruction for the case of pure translation, 
    see Example 6.6 of Ma, Soatto, Kosecká, & Shankar Sastry book
    "Invitation to 3D".
    For the case of pure rotation, see Example 6.10.
    
    */
    
     /**
     * given 2 sets of correspondence from 2 different images taken from
     * 2 cameras whose intrinsic and extrinsic parameters are known,
     * determine the world scene coordinates of the correspondence points.
     * This method simply uses triangulation on each correspondence pair.
     * <pre>
     * following CMU lectures of Kris Kitani at 
     * http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
     * 
     * </pre>
     * @param camera1 image 1 camera matrices of intrinsic and extrinsic parameters.
     * the size is 3 x 4.
     * @param camera2 image 2 camera matrices of intrinsic and extrinsic parameters.
     * the size is 3 x 4.
     * @param x1 the image 1 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @param x2 the image 2 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @return the world scene coordinates and the intrinsic and extrinsic
     * camera matrices (the later were given to the code, but are convenient to return in results).
     * @throws no.uib.cipr.matrix.NotConvergedException 
     */
    public static ReconstructionResults calculateReconstruction(
        CameraParameters camera1, CameraParameters camera2,
        double[][] x1, double[][] x2) throws NotConvergedException {
        
        if (x1.length != 3 || x2.length != 3) {
            throw new IllegalArgumentException("x1.length must be 3 and so must x2.length");
        }
        int n = x1[0].length;
        if (x2[0].length != n) {
            throw new IllegalArgumentException("x1 and x2 must be same dimensions");
        }
        
        /*
        http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
        
        (3) For each point correspondence, compute the point X in 3D space (triangularization)
        */
                        
        double[][] XW = new double[4][n];
        for (int i = 0; i < 4; ++i) {
            XW[i] = new double[n];
        }
        double[] XWPt;
        
        double[][] x1Pt = new double[3][1];
        double[][] x2Pt = new double[3][1];
        int i, j, ii;
        for (i = 0; i < 3; ++i) {
            x1Pt[i] = new double[1];
            x2Pt[i] = new double[1];
        }
                    
        for (i = 0; i < n; ++i) {
            for (ii = 0; ii < 3; ++ii) {
                x1Pt[ii][0] = x1[ii][i];
                x2Pt[ii][0] = x2[ii][i];
            }
            //length is 4
            XWPt = Triangulation.calculateWCSPoint(camera1, camera2, x1Pt, x2Pt);
            for (ii = 0; ii < 4; ++ii) {
                XW[ii][i] = XWPt[ii];
            } 
        }
                
        ReconstructionResults rr = new ReconstructionResults();
        rr.XW = XW;

        return rr;
    }
    
    /**
     * This is also called Projective Structure From Motion for the
     * Two-camera case.   it's a distorted version of euclidean 3d.
     * 
     * NOTE that because the camera calibration, that is, intrinsic parameters,
     * are not known, only the projective reconstruction is possible,
     * but this can be upgraded to 
     affine (parallelism preserved) and Euclidean (parallelism and orthogonality preserved) 
     reconstructions.
     To upgrade to an affine projection, need 3 vanishing points
     (see Section 9.2.2 of Belongie lec 9).
     To directly upgrade from projective to euclidean projection, need
     5 ground truth points in general position, that is, no 4 points
     are coplanar (see Section 9.3 of Belongie lec 9).
     * NOTE: this solution is fine for cases with no noise, otherwise, the
     * results should be the initial values for a non-linear optimization method.
     
     * The method uses the fundamental matrix to compute the camera matrices P1, P2
     * and then uses triangulation on each correspondence pair.
     * 
     * <pre>
     * following CMU lectures of Kris Kitani at 
     * http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
     * 
     * other references:
     * Sect 7.2.1 of Szeliski 2010
     * </pre>
     * @param x1 the image 1 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * NOTE: since intrinsic parameters are not known, users of this method should 
     * presumably center the coordinates in some manner 
     * (e.g. subtract the image center or centroid of points) since internally
     * an identity matrix is used for K.  
     * @param x2 the image 2 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * NOTE: since intrinsic parameters are not known, users of this method should 
     * presumably center the coordinates in some manner 
     * (e.g. subtract the image center or centroid of points).
     * @return the estimated projections P1 and P2 and the objects locations as 3-D points;
     */
    public static ProjectionResults calculateProjectiveReconstruction(
        double[][] x1, double[][] x2) throws NotConvergedException {
                        
        if (x1.length != 3 || x2.length != 3) {
            throw new IllegalArgumentException("x1.length must be 3 and so must x2.length");
        }
        int n0 = x1[0].length;
        if (x2[0].length != n0) {
            throw new IllegalArgumentException("x1 and x2 must be same dimensions");
        }
        if (n0 < 7) {
            throw new IllegalArgumentException("need at least 7 points for the uncalirated camera (s) 2 view solution");
        }
        
        /*
        Similar to CameraPose.calculateUsingEssentialMatrix()
        except for use of uncalibrated points and use of a different diagonal
        in constrcting the fundamental matrix.
                        
        http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
        
        (1) compute fundamental matrix FM from the correspondence x1, x2
        (2) compute the camera matrices P1, P2 from FM.
        (3) For each point correspondence, compute the point X in 3D space (triangularization)
        
        see also notes above from notes from Serge Belongie lectures from Computer Vision II, CSE 252B, USSD
        */
        
        DenseMatrix x1M = new DenseMatrix(x1);
        DenseMatrix x2M = new DenseMatrix(x2);
        
        EpipolarTransformer.NormalizedXY normXY1 = EpipolarTransformer.normalize(x1M);
        EpipolarTransformer.NormalizedXY normXY2 = EpipolarTransformer.normalize(x2M);
        DenseMatrix leftM = normXY1.getXy();
        DenseMatrix rightM = normXY2.getXy();
        
        double tolerance = 3.84; //3.84 5.99 7.82        
        boolean useToleranceAsStatFactor = true;
        ErrorType errorType = ErrorType.SAMPSONS;
        EpipolarTransformationFit fitR = null;
        
        boolean reCalcIterations = true;
        
        //EpipolarTransformer tr = new EpipolarTransformer();
        
        /*
        DenseMatrix normalizedFM = tr.calculateEpipolarProjection(leftM, rightM);
        DenseMatrix vNFM = tr.validateSolution(normalizedFM, leftM, rightM);
        
        Distances distances = new Distances();
        if (useToleranceAsStatFactor) {
            fitR = distances.calculateError2(normalizedFM, leftM, rightM,
                    errorType, tolerance);
        } else {
            fitR = distances.calculateError(normalizedFM, leftM, rightM,
                    errorType, tolerance);
        }
        */
        
        boolean calibrated = false;
                
        RANSACSolver solver = new RANSACSolver();
        fitR = solver.calculateEpipolarProjection(
            leftM, rightM, errorType, useToleranceAsStatFactor, tolerance,
                reCalcIterations, calibrated);
        
        DenseMatrix fm = EpipolarTransformer.denormalizeTheFundamentalMatrix(
            fitR.getFundamentalMatrix(), 
            normXY1.getNormalizationMatrices(),
            normXY2.getNormalizationMatrices());
        
        double[][] _fm = MatrixUtil.convertToRowMajor(fm);
        
        x1M = extractIndices(x1M, fitR.inlierIndexes);
        x2M = extractIndices(x2M, fitR.inlierIndexes);
        x1 = MatrixUtil.convertToRowMajor(x1M);
        x2 = MatrixUtil.convertToRowMajor(x2M);
        
        int n = x1[0].length;
        
        System.out.println("RANSAC fit=" + fitR.toString());
        
        //(2) compute the camera matrices P1, P2 from FM.
        
        MatrixUtil.SVDProducts svdF = MatrixUtil.performSVD(_fm);
        
        assert(svdF.u[0].length == 3 && svdF.u.length == 3);
        
        double detU = MatrixUtil.determinant(svdF.u);
        double detV = MatrixUtil.determinant(svdF.vT);
        
        System.out.printf("SVD.u=\n%s\n", FormatArray.toString(svdF.u, "%.3e"));
        System.out.printf("SVD.s=\n%s\n", FormatArray.toString(svdF.s, "%.3e"));
        System.out.printf("SVD.vT=\n%s\n", FormatArray.toString(svdF.vT, "%.3e"));
        System.out.printf("det(SVD.u)=%.2f\n", detU);
        System.out.printf("det(SVD.vT)=%.2f\n", detV);
        
        //H = +U * R_Z(+-90)^T * sigma_hat * V^T
        double[][] r90T = new double[3][3];
        r90T[0] = new double[]{0, 1, 0};
        r90T[1] = new double[]{-1, 0, 0};
        r90T[2] = new double[]{0, 0, 1};
        
        double[][] h = MatrixUtil.multiply(
            MatrixUtil.multiplyByDiagonal(
                MatrixUtil.multiply(svdF.u, r90T),
                new double[]{svdF.s[0], svdF.s[1], svdF.s[2]}),
            svdF.vT);
        
        // last column in u is the second epipole and is the direction of vector t
        double[] t1 = MatrixUtil.extractColumn(svdF.u, 2);
        double[] e2 = t1;
        double[] e1 = svdF.vT[2];
        
        System.out.printf("e1=%s\n", FormatArray.toString(e1, "%.3e"));
        System.out.printf("e2=%s\n", FormatArray.toString(e2, "%.3e"));
         
        // camera matrix P for left image = [I | 0 ]
        double[][] camera1 = MatrixUtil.zeros(3, 4);
        for (int i = 0; i < 3; ++i) {
            camera1[i][i] = 1;
        }
        
        // Szeliksi 2010 eqn 7.34:  P0 =[I|0] and P0 =[H|e],
        //    and then he finishes w/ triangulation
        // kitani lecture has P2 = [ [e1]_x * F | e2 ]
        // Hartley & Zisserman: has P2 = [ -[e2]_x^T * F | e2 ]
        //    note that slide 59 of http://16720.courses.cs.cmu.edu/lec/sfm.pdf
        //    also uses the Hartley & Zisserman: version of P2 and refers to proof in 
        //    Forsyth & Ponce Sec 8.3
        // Belongie uses P2 = [ [e2]_x^T * F | e2 ] and so does Ma et al. textbook (msks, "Invitation to 3D"
        // 
        
        // NOTE: not necessary to normalize the epipoles by the last value
        
        double[][] camera2 = MatrixUtil.zeros(3, 4);
        for (int i = 0; i < 3; ++i) {
            camera2[i] = new double[4];
            System.arraycopy(h[i], 0, camera2[i], 0, 3);
            camera2[i][3] = e2[i];
        }
        
        System.out.printf("Szeliski P2\n%s\n", FormatArray.toString(camera2, "%.3e"));
        
        CameraProjection P1 = new CameraProjection(camera1);
        CameraProjection P2 = new CameraProjection(camera2);
        
        double[][] XW;
        double[] XWPt;
        XW = new double[4][n];
        for (int i = 0; i < 4; ++i) {
            XW[i] = new double[n];
        }
        
        double[][] x1Pt = new double[3][1];
        double[][] x2Pt = new double[3][1];
        int i, j;
        for (i = 0; i < 3; ++i) {
            x1Pt[i] = new double[1];
            x2Pt[i] = new double[1];
        }
        
        for (i = 0; i < n; ++i) {
            for (j = 0; j < 3; ++j) {
                x1Pt[j][0] = x1[j][i];
                x2Pt[j][0] = x2[j][i];
            }
            XWPt = Triangulation.calculateWCSPoint(P1, P2, x1Pt, x2Pt);
            for (j = 0; j < 4; ++j) {
                XW[j][i] = XWPt[j];
            }
        }
        
        System.out.printf("Szeliski reconstruction\n%s\n", FormatArray.toString(XW, "%.3e"));
        
        ProjectionResults rr = new ProjectionResults();
        rr.XW = XW;
        rr.projectionMatrices = MatrixUtil.zeros(3*2, 4);
        for (i = 0; i < 3; ++i) {
            System.arraycopy(P1.getP()[i], 0, rr.projectionMatrices[i], 0, 4);
        }
        for (i = 0; i < 3; ++i) {
            System.arraycopy(P2.getP()[i], 0, rr.projectionMatrices[3 + i], 0, 4);
        }
        return rr;
    }
   
    /**
     * TODO: proof read the algorithm and write test for this.
     * for the case of un-calibrated cameras viewing the same scene features,
     * recover the 3-D coordinates in WCS and the projection matrices 
     * from pairs of corresponding
     * un-calibrated image points, that is, points in the image reference frame in pixels.
     * 
     * The method implements the Sturm & Triggs 1996 algorithm: 
     "a method for the recovery of projective shape and motion from multiple 
     images of a scene by the factorization of a matrix containing the images 
     of all points in all views. This factorization is only possible when the
     image points are correctly scaled. The major technical contribution of 
     this paper is a practical method for the recovery of these scalings, 
     using only fundamental matrices and epipoles estimated from the image data."
     "[it is a] closed form solutions, not iterative bundle-adjustment..."
     * <pre>
     * references:
     * 
     * Sturm and Triggs 1996, 
    "A Factorization Based Algorithm for Multi-Image Projective Structure and Motion"
     https://link.springer.com/content/pdf/10.1007/3-540-61123-1_183.pdf
    
    see also proj_recons_fsvd.m from http://lear.inrialpes.fr/people/triggs/src/
    which has a very liberal copyright in the file COPYRIGHT
    Copyright Bill Triggs (http://www.inrialpes.fr/movi/people/Triggs),
    INRIA (http://www.inria.fr) and CNRS (http://www.cnrs.fr),
    1995-2002. All rights reserved.

    You may use and distribute [*] this work with or without modification,
    for any purpose and without fee or royalty, subject to the following
    conditions:
       (see file COPYRIGHT)
    
     * </pre>
     * 
     * NOTE: Sturm & Triggs 1996 state in their code, "% The projective output 
     frame is numerically well-conditioned, but otherwise *completely* 
     arbitrary. It has *no* relation to any Euclidean frame.
     * 
     * @param x the image coordinates of feature correspondences in 2 or more
     * images.  format is 2 X (nImages * nFeatures) where row 0 holds the x-coordinates
     * and row 1 holds the y-coordinates and each image's features are given
     * before the next and the features are ordered in the same manner within
     * all images.
     * for example: row 0 = [img_0_feature_0, ... img_0_feature_n-1, ... img_m-1_feature_0,...
     *     img_m-1_feature_n-1].
     * @param mImages the number of images in x.
     * @return the estimated projections P1 and P2 and the objects locations as 3-D points;
     */
    public static ProjectionResults calculateProjectiveReconstruction(
        double[][] x, int mImages) throws NotConvergedException {
                        
        if (x.length != 2) {
            throw new IllegalArgumentException("x.length must be 2");
        }
        if ((x[0].length % mImages) != 0) {
            throw new IllegalArgumentException("x must have a multiple of mImages as the number of columns");
        }
        int nFeatures = x[0].length / mImages;
        
        // need at least 7 points in each image for the point version of fundamental
        // matrix.
        // not implementing the line version as Triggs in another paper states
        //    that they may be more affected by outliers
        if (nFeatures < 7) {
            throw new IllegalArgumentException("need at least 7 points per image");
        }
        
        /*
        3.3 Outline of the Algorithm
        The complete algorithm is composed of the following steps.
        i. Normalize the image coordinates, by applying transformations Ti.
        2. Estimate the fundamental matrices and epipoles with the method of [Har95].
        3. Determine the scale factors Aip using equation (3).
        4. Build the rescaled measurement matrix W.
        5. Balance W by column-wise and "triplet-of-rows"-wise scalar mutliplications.
        6. Compute the SVD of the balanced matrix W.
        7. From the SVD, recover projective motion and shape.
        8. Adapt projective motion, to account for the normalization transformations Ti of
        step 1.
        */
                
        // following proj_recons_fsvd.m
        //    pairs of image sets can be formed either by using the first
        //    image as x1 for all images, or chaining them all together.
        // i.e. [(0,1), (0, 2), (0,3)] or [(0,1), (1,2), (2,3)].
        // choosing the later here.
        
        RANSACSolver ransac = new RANSACSolver();
        double[] e12 = new double[3];
        EpipolarTransformationFit fit;
        double[][] fm;
        ErrorType errorType = ErrorType.DIST_TO_EPIPOLAR_LINE;
        boolean useToleranceAsStatFactor = true;
        double tolerance = 3.84;
        // TODO: estimate this:
        boolean recalcIterations = (nFeatures > 100 || (mImages > 10 && nFeatures > 20));
        boolean calibrated = false;
        
        int i, j;
        // image pairs extracted from x:
        double[][] x1, x2;
        // the normalization centroidX, centroidY, and sigma values applied to x1 and x2:
        double[] tt = new double[3*mImages];
        
        // use a reference depth of 1 for first image's features, as have no measurements for any depths to bootstrap from.
        double[][] lambdas = MatrixUtil.zeros(mImages, nFeatures);
        Arrays.fill(lambdas[0], 1.);
        
        x1 = extractAndNormalize(x, 0, nFeatures, tt);
        DenseMatrix x1M, x2M;
        x1M = new DenseMatrix(x1);
        double[] x1p = new double[3];
        double[] x2p = new double[3];
        double[] x2e, tmp;
        double tmp2, x2esq;
        
        // format x into shape W (3*mImages X nFeatures):
        //  row 0:2 = image 1 points where row 0 is the x coordinates, row 1 is the y coordinates
        //  row 3:5 = image 2 points
        double[][] w = MatrixUtil.zeros(3*mImages, nFeatures);
        System.arraycopy(x1[0], 0, w[0], 0, nFeatures);
        System.arraycopy(x1[1], 0, w[1], 0, nFeatures);
        System.arraycopy(x1[2], 0, w[2], 0, nFeatures);
        
        for (i = 1; i < mImages; ++i) {
            
            x2 = extractAndNormalize(x, i, nFeatures, tt);
            x2M = new DenseMatrix(x2);
            
            System.arraycopy(x2[0], 0, w[i*3], 0, nFeatures);
            System.arraycopy(x2[1], 0, w[i*3+1], 0, nFeatures);
            System.arraycopy(x2[2], 0, w[i*3+2], 0, nFeatures);
            
            fit = ransac.calculateEpipolarProjection(x1M, x2M, errorType, 
                useToleranceAsStatFactor, tolerance, recalcIterations, calibrated);
            
            fm = MatrixUtil.convertToRowMajor(fit.getFundamentalMatrix());
            
            /*
            TODO: consider keeping only the inliers in a future version that handles
            occlusion.  by imputation or applied factorization or other means...
            x1M = extractIndices(x1M, fitR.inlierIndexes);
            x2M = extractIndices(x2M, fitR.inlierIndexes);
            x1 = MatrixUtil.convertToRowMajor(x1M);
            x2 = MatrixUtil.convertToRowMajor(x2M);
            int nFeaturesI = x1[0].length;
            */
            
            // note: e12 is not normalized by last component
            calculateLeftEpipole(fit.getFundamentalMatrix(), e12);
            
            for (j = 0; j < nFeatures; ++j) {
                extractColumn(x1, j, x1p);
                extractColumn(x2, j, x2p);
                
                //xe = cross(x2_j, e12);  // same as x2 cross -e21
                x2e = MatrixUtil.crossProduct(x2p, e12);
                
                //lambda(i,j) = lambda(x1,j) * abs((x1_j' * FM * xe) / (xe' *xe));
                //    note: epipolar line2 = x1_j' * FM
                x2esq = MatrixUtil.innerProduct(x2e, x2e);
                
                tmp = MatrixUtil.multiplyRowVectorByMatrix(x1p, fm);
                
                tmp2 = MatrixUtil.innerProduct(tmp, x2e);
                lambdas[i][j] = lambdas[i-1][j] * Math.abs(tmp2/x2esq);
            }
            
            x1 = x2;
            x1M = x2M;
        }
        
        /*
        4. Build the rescaled measurement matrix W.
        5. Balance W by column-wise and "triplet-of-rows"-wise scalar mutliplications.
        
        as stated in Sturm & Triggs 1996 Section 3.2, the balancing of the
        rescaled measurement matrix by Q's then P's can be replaced by
        balancing the m x n matrix lambdas instead of W because of the simplification
        of working with normalized image coordinates Q.
        The balance operations are demonstrated in proj_recons_fsvd.m
        */
        double eps = 1.e-11;
        double[] lambdaj;
        int k;
        // authors find 2 iterations is heuristically enough:
        for (i = 0; i < 2; ++i) {
            for (j = 0; j < nFeatures; ++j) {
                lambdaj = MatrixUtil.extractColumn(lambdas, j);
                tmp2 = MatrixUtil.lPSum(lambdaj, 2);
                if (Math.abs(tmp2) < eps) { 
                    tmp2 = eps;
                }
                for (k = 0; k < mImages; ++k) {
                    lambdas[k][j] = lambdaj[k] / tmp2;
                }
            }
            for (k = 0; k < mImages; ++k) {
                tmp2 = MatrixUtil.lPSum(lambdas[k], 2);
                if (Math.abs(tmp2) < eps) { 
                    tmp2 = eps;
                }
                MatrixUtil.multiply(lambdas[k], 1./tmp2);
            }
        }
        
        // rescale the image points
        for (i = 0; i < mImages; ++i) {
            for (j = 0; j < nFeatures; ++j) {
                w[3*i + 0][j] *= lambdas[i][j];
                w[3*i + 1][j] *= lambdas[i][j];
                w[3*i + 2][j] *= lambdas[i][j];
            }
        }
        
        /*
        6. Compute the SVD of the balanced matrix W.
        7. From the SVD, recover projective motion and shape.
        8. Adapt projective motion, to account for the normalization transformations Ti of
        step 1.
        */
        
        // if the number of images is larger than 10 or the number of features
        //   is greater than 30, will use cur decomposition
        double[][] u, vT;
        double[][] s;
        //double[][] wRescaled;
        if (mImages > 10 || nFeatures > 30) {
            CUR cur = CURDecomposition.calculateDecomposition(w, 4);
            SVDProducts curSVD = cur.getApproximateSVD();
            u = curSVD.u;
            vT = curSVD.vT;
            s = curSVD.sigma;
            
            // Sturm & Triggs divide by the largest eigenvalue:
            MatrixUtil.multiply(s, 1./s[0][0]);
            
            //wRescaled = cur.getResult();
        } else {
            SVDProducts svd = MatrixUtil.performSVD(w);
            
            //reduce rank to 4
            u = MatrixUtil.copySubMatrix(svd.u, 0, svd.u.length-1, 0, 3);
            vT = MatrixUtil.copySubMatrix(svd.vT, 0, 3, 0, svd.vT[0].length-1);
            
            s = MatrixUtil.zeros(4, 4);
            s[0][0] = svd.s[0];
            s[1][1] = svd.s[1];
            s[2][2] = svd.s[2];
            s[3][3] = svd.s[3];
            
            /*double[][] sqrts4 = MatrixUtil.zeros(4, 4);
            sqrts4[0][0] = Math.sqrt(svd.s[0]);
            sqrts4[1][1] = Math.sqrt(svd.s[1]);
            sqrts4[2][2] = Math.sqrt(svd.s[2]);
            sqrts4[3][3] = Math.sqrt(svd.s[3]);

            u = MatrixUtil.multiply(u4, sqrts4);
            vT = MatrixUtil.multiply(sqrts4, vT4);
            wRescaled = MatrixUtil.multiply(u, vT);
            */
            
            // Sturm & Triggs divide by the largest eigenvalue:
            MatrixUtil.multiply(s, 1./s[0][0]);
        }
        
        //Ps = fliplr(U(:,1:4));             // 3*MImages X 4
        //Xs = flipud(S(1:4,1:4)*V(:,1:4)'); // 4 X 4*mImages
        double[][] ps = MatrixUtil.copy(u);
        MatrixUtil.flipLR(ps);
        double[][] XW = MatrixUtil.multiply(s, vT);
        MatrixUtil.flipUD(XW);
   
        //denormalize ps.  = ps * T^-1
        // tt is the normalization centroidX, centroidY, and sigma values applied to x1 and x2:
        /*
                 | 1  0  xc |   | s  0   0 |   | s   0  xc |
          T^-1 = | 0  1  yc | * | 0  s   0 | = | 0   s  yc |
                 | 0  0   1 |   | 0  0   1 |   | 0   0   1 |
        */
        double[][] tInv = MatrixUtil.zeros(3, 3);
        tInv[2][2] = 1;
        double ts, txc, tyc;
        double[] p0, p1, p2;
        double[][] p = MatrixUtil.zeros(3, 4);
        for (i = 0; i < mImages; ++i) {
            txc = tt[3*i];
            tyc = tt[3*i + 1];
            ts = tt[3*i + 2];
            tInv[0][0] = ts;
            tInv[1][1] = ts;
            tInv[0][2] = txc;
            tInv[1][2] = tyc;
            
            p0 = ps[3*mImages];
            p1 = ps[3*mImages + 1];
            p2 = ps[3*mImages + 2];
            System.arraycopy(p0, 0, p[0], 0, 4);
            System.arraycopy(p1, 0, p[1], 0, 4);
            System.arraycopy(p2, 0, p[2], 0, 4);
            
            p = MatrixUtil.multiply(p, tInv);
            
            System.arraycopy(p[0], 0, ps[3*mImages], 0, 4);
            System.arraycopy(p[1], 0, ps[3*mImages + 1], 0, 4);
            System.arraycopy(p[2], 0, ps[3*mImages + 2], 0, 4);
        }
        
        ProjectionResults rr = new ProjectionResults();
        rr.XW = XW;
        rr.projectionMatrices = ps;
        
        return rr;
    }
    
     /**
     * given correspondence between two images in image coordinates calculate 
     * the extrinsic camera parameters and the 3-D points.
     * 
     * This method calculates the essential matrix and uses the SVD of it to
     * extract the translation and possible rotation matrices which are
     * filtered to find the best while calculating triangulation for each point.
     * 
     * Note that the absolute translation between the two cameras can never be 
     * recovered from pure image measurements alone, regardless of how many 
     * cameras or points are used as ground control points are
     * needed.
     * <pre>
     * following CMU lectures of Kris Kitani at 
     http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
     Szeliski 2010, Chapter 7, and eqn (7.25).
     Ma, Soatto, Kosecká, and Sastry 2012, "An Invitation to 3-D Vision", pg 121 
     * </pre>
     * @param k1 intrinsic camera matrix for image 1 in units of pixels.
     * @param k2 intrinsic camera matrix for image 2 in units of pixels.
     * @param x1 the image 1 set of correspondence points in image reference frame.  
     * format is 3 x N where N is the number of points.
     * @param x2 the image 2 set of correspondence points in image reference frame.  
     * format is 3 x N where N is the number of points.
     * @return 
     * @throws no.uib.cipr.matrix.NotConvergedException 
     */
    public static ReconstructionResults calculateUsingEssentialMatrix(
        double[][] k1, double[][] k2,
        double[][] x1, double[][] x2) throws NotConvergedException {
                
        if (x1.length != 3 || x2.length != 3) {
            throw new IllegalArgumentException("x1.length must be 3 and so must x2.length");
        }
        int n = x1[0].length;
        if (x2[0].length != n) {
            throw new IllegalArgumentException("x1 and x2 must be same dimensions");
        }
        
        /*
        http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
        
        (1) compute fundamental mat5rix FM from the correspondence x1, x2
        (2) compute the camera matrices P1, P2 from FM.
        (3) For each point correspondence, compute the point X in 3D space (triangularization)
        */
        
        double[][] k1IntrInv = Camera.createIntrinsicCameraMatrixInverse(k1);
        double[][] k2IntrInv = Camera.createIntrinsicCameraMatrixInverse(k2);
        
        // the direction of the points is calculated by K^-1 * x
        double[][] x1Direction = MatrixUtil.multiply(k1IntrInv, x1);
        double[][] x2Direction = MatrixUtil.multiply(k2IntrInv, x2);
                
        DenseMatrix x1M = new DenseMatrix(x1Direction);
        DenseMatrix x2M = new DenseMatrix(x2Direction);
        
        // normalizing by unit standard to improve results of epipolar solution:
        EpipolarTransformer.NormalizedXY normXY1 
            = EpipolarTransformer.normalizeUsingUnitStandard(x1M);
        EpipolarTransformer.NormalizedXY normXY2 
            = EpipolarTransformer.normalizeUsingUnitStandard(x2M);
        
        DenseMatrix leftM = normXY1.getXy();
        DenseMatrix rightM = normXY2.getXy();
        
        double tolerance = 3.84; //3.84 5.99 7.82        
        boolean useToleranceAsStatFactor = true;
        ErrorType errorType = ErrorType.SAMPSONS;
        EpipolarTransformationFit fitR = null;
        boolean reCalcIterations = true;
        
        DenseMatrix normalizedE;
        
        /*EpipolarTransformer tr = new EpipolarTransformer();        
        normalizedE = tr.calculateEpipolarProjection(leftM, rightM);
        DenseMatrix vNFM = tr.validateSolution(normalizedFM, leftM, rightM);
        
        Distances distances = new Distances();
        if (useToleranceAsStatFactor) {
            fitR = distances.calculateError2(vNFM, leftM, rightM,
                    errorType, tolerance);
        } else {
            fitR = distances.calculateError(vNFM, leftM, rightM,
                    errorType, tolerance);
        }*/
        
        //if true, solves for the Essential Matrix, else solves
        //for the Fundamental Matrix.  The difference is in the diagonal used for
        //dimension reduction.
        boolean coordsAreInCameraRefFrame = true;
        
        RANSACSolver solver = new RANSACSolver();
        fitR = solver.calculateEpipolarProjection(
            leftM, rightM, errorType, useToleranceAsStatFactor, tolerance,
                reCalcIterations, coordsAreInCameraRefFrame);
        
        System.out.println("RANSAC fit=" + fitR.toString());
        
        normalizedE = fitR.getFundamentalMatrix();
        
        // this is now back in the reference frame of the x1Direction and x2Direction
        DenseMatrix essentialM = EpipolarTransformer.denormalizeTheFundamentalMatrix(
            normalizedE, normXY1.getNormalizationMatrices(),
            normXY2.getNormalizationMatrices());
        
        //=========
        double[][] _essentialMatrix = MatrixUtil.convertToRowMajor(essentialM);
        
        MatrixUtil.SVDProducts svdE = MatrixUtil.performSVD(_essentialMatrix);
        
        assert(svdE.u[0].length == 3 && svdE.u.length == 3);
        
        double detU = MatrixUtil.determinant(svdE.u);
        double detV = MatrixUtil.determinant(svdE.vT);
        
        System.out.printf("SVD.u=\n%s\n", FormatArray.toString(svdE.u, "%.3e"));
        System.out.printf("SVD.s=\n%s\n", FormatArray.toString(svdE.s, "%.3e"));
        System.out.printf("SVD.vT=\n%s\n", FormatArray.toString(svdE.vT, "%.3e"));
        System.out.printf("det(SVD.u)=%.2f\n", detU);
        System.out.printf("det(SVD.vT)=%.2f\n", detV);
        
        /*
        Szeliski 2010 chap 7:
        Once an estimate for the essential matrix E has been recovered, 
        the direction of the translation vector t can be estimated. 
        
        Note that the absolute distance between the two cameras can never 
        be recovered from pure image measurements alone without knowledge 
        about absolute camera and point positions or distances, often called ground 
        control points in photogrammetry.
        */
        
        // det(R)=1 is a proper rotation matrix.  rotation angles are counterclockwise.
        //           it's a special orthogonal matrix and provides the
        //           defining matrix representation of the group of proper n-dimensional rotations, denoted
        //           by SO(n). http://scipp.ucsc.edu/~haber/ph251/rotreflect_17.pdf
        // det(R)=-1 is an improper rotation matrix representing rotations that
        //           require mirrors.
        //           The most general improper rotation matrix is a product of a proper rotation by an
        //           angle θ about some axis nˆ and a mirror reflection through a plane that passes through
        //           the origin and is perpendicular to nˆ.  NOTE: nˆ is determined by
        //           the right hand rule.
        
        /*
        Sect 7.2 of Szeliski 2010 eqn (7.25) introduces
        R3 and R4 constructed from -U as 2 more rotation possibilities to be tested
        and that is necessary in some cases where det(R) would otherwise be -1
        (reflection).
        ...we only know both E and tˆup to a sign. Furthermore, the matrices U and V
        are not guaranteed to be rotations (you can flip both their signs and 
        still get a valid SVD).   
        For this reason, we have to generate all four possible rotation matrices
        
        R = +-U * R_Z(+-90)^T * V^T
           and keep the 2 whose determinant = 1
        */
        double[][] r1 = MatrixUtil.zeros(3, 3);
        double[][] r2 = MatrixUtil.zeros(3, 3);
        double[][] u = MatrixUtil.zeros(3, 3);
        populateWithDet1Rs(svdE, r1, r2, u);
        
        // last column in u is the second epipole and is the direction of vector t
        double[] t1 = MatrixUtil.extractColumn(u, 2);
        double[] t2 = Arrays.copyOf(t1, t1.length);
        MatrixUtil.multiply(t2, -1); 
        
        System.out.printf("R1=\n%s\n", FormatArray.toString(r1, "%.3e"));
        System.out.printf("R2=\n%s\n", FormatArray.toString(r2, "%.3e"));
        System.out.printf("t1=\n%s\n", FormatArray.toString(t1, "%.3e"));
        System.out.printf("t2=\n%s\n", FormatArray.toString(t2, "%.3e"));
        
        //then of the 4 possible choices find the one with largest number of positive Z.
          
        //NOTE: the last column vector in u is the smallest
        //    eigenvector.  it is epipole2, that is, the right image position 
        //    of the epipolar projection of the left camera center.
        //    it's int the left null space of E.
                        
        x1M = extractIndices(new DenseMatrix(x1Direction), fitR.inlierIndexes);
        x2M = extractIndices(new DenseMatrix(x2Direction), fitR.inlierIndexes);
        x1 = MatrixUtil.convertToRowMajor(x1M);
        x2 = MatrixUtil.convertToRowMajor(x2M);
        
        double[][] rSelected = MatrixUtil.zeros(3, 3);
        double[] tSelected = new double[3];
        double[][] XW = MatrixUtil.zeros(4, x1[0].length);
        bestInCheiralityTest(x1, x2, k1, k2, r1, r2, t1, t2, rSelected, tSelected, XW);  
        
        ReconstructionResults rr = new ReconstructionResults();
        rr.XW = XW;
        rr.k1ExtrRot = MatrixUtil.createIdentityMatrix(3);
        rr.k1ExtrTrans = new double[]{0, 0, 0};
        rr.k1Intr = k1;
        rr.k2ExtrRot = rSelected;
        rr.k2ExtrTrans = tSelected;
        rr.k2Intr = k2;
        rr.essentialMatrix = _essentialMatrix;
        rr.svd = svdE;
        rr.fundamentalMatrix = null;

        return rr;        
    }

    /**
     * NOTE: not ready for use yet.
     * 
     * TODO: proof read the algorithm and write test for this.
     * for the case where the cameras are viewing small, distant scenes,
     * recover the 3-D coordinates in WCS and the rotation matrices 
     * from pairs of corresponding
     * un-calibrated image points, that is, points in the image reference frame in pixels.
     * assumes an orthographic camera model.
     * can use the orthographic camera model when
     *    (the average distance of an object from the camera) 
     *     .geq. 10*(the average width of the object (measured along the optical axis of the camera).
     * <pre>
     * references:
     * 
     * lecture 16 notes from Serge Belongie lectures from Computer Vision II, CSE 252B, USSD
     * http://www-cse.ucsd.edu/classes/sp04/cse252b/notes/lec16/lec16.pdf
     * 
     * lectures of Deva Ramanan at http://16720.courses.cs.cmu.edu/lec/sfm.pdf
     * .:w
     * 
     * Tomasi & Kanade 1991, "Shape and motion from image streams under 
     * orthography: a factorization method", International journal of computer vision 
     * 
     *  Morita and Kanade 1997 for solving Q.
         T. Morita and T. Kanade, A Sequential Factorization Method for Recovering Shape and Motion
         from Image Streams, Pattern Analysis and Machine Intelligence, IEEE Transactions on, vol. 19,
         no.8, pp.858-867, Aug 1997  (1994?)
         
     * Higham, 1988, “Computing a Nearest Symmetric Positive Semidefinite Matrix,” 
     *    Linear Algebra and Appl., 103:103-118, 1988
     * 
     * a great summary of the above:
     * http://note.sonots.com/SciSoftware/Factorization.html#cse252b
     * http://note.sonots.com/?plugin=attach&refer=SciSoftware%2FFactorization&openfile=Factorization.pdf
     * 
     * and a derivation of the geometry of the tracking equation:
     * Birchfield 1997, "Derivation of Kanade-Lucas-Tomasi Tracking Equation"
     * http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.185.413&rep=rep1&type=pdf
     * </pre>
     * NOTE: could overload this method to enable handling of occlusion 
     * following Section 5 of Tomasi & Kanade 1991, but might want to alter the
     * algorithm to use geometric median in place of centroid so that the
     * "centers" are not as affected by removing or adding a point.
     * NOTE: comments from Poelman & Kanade 1992:
     * Orthographic projection does not account for the apparent change in size 
     * of an object as it moves toward or away from the camera, nor the different 
     * angle from which an object is viewed as it moves parallel to the image plane.
     * NOTE: consider implementing Section 3.3 Sequential Factorization Algorithm
     * from the Morita & Kanade 1997 paper (1994?)
     * @param x the image coordinates of feature correspondences in 2 or more
     * images.  format is 2 X (nImages * nFeatures) where row 0 holds the x-coordinates
     * and row 1 holds the y-coordinates and each image's features are given
     * before the next and the features are ordered in the same manner within
     * all images.
     * for example: row 0 = [img_0_feature_0, ... img_0_feature_n-1, ... img_m-1_feature_0,...
     *     img_m-1_feature_n-1].
     * @param mImages the number of images in x.
     * @return the estimated projections P1 and P2 and the objects locations as 3-D points;
     */
    public static OrthographicProjectionResults calculateAffineReconstruction(
        double[][] x, int mImages) throws NotConvergedException {
                                
        if (x.length != 2) {
            throw new IllegalArgumentException("x.length must be 2");
        }
        if ((x[0].length % mImages) != 0) {
            throw new IllegalArgumentException("x must have a multiple of mImages as the number of columns");
        }
        int nFeatures = x[0].length / mImages;
        
        //2mn >= 8m + 3n – 12
        if ((2*mImages * nFeatures) < (8*mImages +3*nFeatures - 12)) {
            throw new IllegalArgumentException("more points are necessary:"
                + "2 * mImages * nFeatures >= 8 * mImages + 3 * nFeatures – 12."
                +  "\narguments: mImages=" + mImages + " nFeatures=" + nFeatures
                + "\n" + (2*mImages * nFeatures) + " < " +
                    (8*mImages +3*nFeatures - 12));
        }
        // for mImages=2, need 4 features
        
        // Tomasi & Kanade 1992: An image stream can be represented by the
        //  2*F X P measurment matrix W of image coordinates of P points
        // tracked through F frames.
        
        // create matrix W composed of U and V 
        //     where U is rows of each image's x coordinates (size is mImages X nFeatures).
        //     where V is rows of each image's y coordinates (size is mImages X nFeatures).
        // create matrix t which holds the centroids of each row of W
        // create matrix WC = W - t
        
        // row[0] = image_0_x[0], image_0_x[1], image_0_x[2], ...
        // row[1] = image_1_x[0], image_1_x[1], image_1_x[2], ...
        // ...
        // row[mImages+0] = image_0_y[0], image_0_y[1], image_0_y[2], ...
        double[][] w = MatrixUtil.zeros(2*mImages, nFeatures);
        // t points to the camera's focal point
        double[] t = new double[2*mImages];
        int i, j;
        int m, n, rowU, rowV, col;
        for (m = 0; m < mImages; ++m) {
            rowU = m;
            rowV = mImages + m;
            for (n = 0; n < nFeatures; ++n) {
                col = nFeatures*m + n;
                w[rowU][n] = x[0][col];
                t[rowU] += x[0][col];
                w[rowV][n] = x[1][col];
                t[rowV] += x[1][col];
            }
            t[rowU] /= (double)nFeatures;
            t[rowV] /= (double)nFeatures;
        }
        
        //registered measurement matrix:
        double[][] wC = MatrixUtil.copy(w);
        for (i = 0; i < t.length; ++i) {
            for (n = 0; n < wC[i].length; ++n) {
                wC[i][n] -= t[i];
            }
        }
        
        // see Fig 3.1 of Tomasi & Kanade 1991 or Fig 2. of Belongie lecture notes
        
        // under orthography, this coordinate centered matrix has rank 3
        // and can be factored into the product of 2 matrices, R and S:
        // the camera rotation matrix R (size 2FX3) and 
        // the shape matrix S (size 3XP) which is shape in a coordinate system 
        // attached to the object centroid.
        //The two components of the camera translation along the image plane
        // are computed as averages of the rows of W
                
        // the registered measurement matrix is highly rank deficient
        // eqn (3.11)
        SVDProducts svd = MatrixUtil.performSVD(wC);
        
        // Tomasi & Kanade 1992,eqn (3.12): 1st 3 columns of U, upper 3X3 of S, and first 3 rows of V^T
        //U3 is 2FX3 where F is the number of image frames
        //S3 is 3X3
        //VT3 is 3XP where P is the number of points, that is features, per image
        double[][] u3 = MatrixUtil.copySubMatrix(svd.u, 0, svd.u.length-1, 0, 2);
        double[][] s3 = MatrixUtil.zeros(3, 3);
        s3[0][0] = svd.s[0];
        s3[1][1] = svd.s[1];
        s3[2][2] = svd.s[2];
        double[][] sqrts3 = MatrixUtil.copy(s3);
        sqrts3[0][0] = Math.sqrt(sqrts3[0][0]);
        sqrts3[1][1] = Math.sqrt(sqrts3[1][1]);
        sqrts3[2][2] = Math.sqrt(sqrts3[2][2]);
        double[][] vT3 = MatrixUtil.copySubMatrix(svd.vT, 0, 2, 0, svd.vT[0].length-1);
        
        // if the ratio between the 3rd and 4th largest singular value of the registered measurement matrix
        //   is large, then the noise portion of the full decomposition
        //   can be ignored (the noise portion is the block partitions not 
        //   copied to U3, sqrtS3 and vT3).
        // there is more about this in Morita & Kanade 1994/1997
        double sRatio = svd.s[2]/svd.s[3];
        System.out.printf("svd.s[2]/svd.s[3]=%.3e\n", sRatio);
        
        //Tamasi & Kanade
        // (2*mImages)X3
        double[][] rC = MatrixUtil.multiply(u3, sqrts3);
        // 3XnFeatures
        double[][] sC = MatrixUtil.multiply(sqrts3, vT3);
        
        System.out.printf("rC=\n%s\n", FormatArray.toString(rC, "%.4e"));
        
        System.out.printf("wC=\n%s\n", FormatArray.toString(wC, "%.4e"));
        System.out.printf("rC*sC=\n%s\n", FormatArray.toString(
            MatrixUtil.multiply(rC, sC), 
            "%.4e"));
        
        // wC = rC * sC  
        //rC and sC are linear translations of the true rotation matrix R and
        //  the true shape matrix S, respectively.
        // Morita and Kanade: the decomposition is not completely unique.  
        //     it's unique only up to an affine transformation
        
        //Tomasi & Kanade eqn (3.15) and Belongie Section 16.4.4 (c)
        // metric constraints:  
        // note that R is composed of rows of unit vectors.
        // note that the first F rows in R are orthogonal to the last F rows in R.
        
        /*
        NOTE: can make a rotation matrix orthonormal:
        svd = MatrixUtil.performSVD(rotationMAtrix);
        ortho = MatrixUtil.multiply(svd.u, MatrixUtil.transpose(svd.vT));
        detR = MatrixUtil.determinant(ortho);
        assert(Math.abs(detR - 1.)<1.e-7);

        if no translation, can use procrustes to get difference in rotation:
           double[][] ar = Rotation.procrustesAlgorithmForRotation(rot1, _r2);
        */
        
        double[][] q = solveForTransformationToOrthoNormal(rC);
        
       
        //Q is an affine transformation which transforms rC into R in motion space
        //   and the inverse of Q transforms sC into S in the shape space
        
        // finding Q is called "Metric Transformation"
        
        // rC size is  (2*mImages)X3
        // sC size is 3XnFeatures
        double[][] r2 = MatrixUtil.multiply(rC, q);
        
        assertDotProductMetrics(r2, mImages);
        
        double[][] s2 = MatrixUtil.multiply(MatrixUtil.pseudoinverseRankDeficient(q), sC);
        
        // assert that wC is the same as wC2
        System.out.printf("r2*s2=\n%s\n", FormatArray.toString(MatrixUtil.multiply(r2, s2), 
            "%.4e"));
        
        /*
        from Tomasi & Kanade 1992:
        If desired, align the first camera reference system with the world 
        reference system by forming the products R*R_0 and R_0^T*S, 
        where the orthonormal matrix R_0 = [i1 j1 k1] rotates the first camera 
        reference system into the identity matrix
        
        i0x i0y i0z    *  r0ix  r0jx  r0kx   = (i0 dot r0i)  (i0 dot r0j) (i0 dot r0k)
        i1x i1y i1z       r0iy  r0jy  r0ky
        i2x i2y i2z       r0iz  r0jz  r0kz
        i3x i3y i3z
        ...
        j0x j0y j0z
        j1x j1y j1z
        */
        
        System.out.printf("r2=\n%s\n", FormatArray.toString(r2, "%.4e"));
 
        /* r2 * R0 = I
                R0 = inv(r2)
        */
        double[][] rFirst = new double[3][3];
        rFirst[0] = Arrays.copyOf(r2[0], r2[0].length);
        rFirst[1] = Arrays.copyOf(r2[mImages], r2[mImages].length);
        rFirst[2] = MatrixUtil.crossProduct(rFirst[0], rFirst[1]);
        double[][] r0 = MatrixUtil.pseudoinverseFullColumnRank(rFirst);
        
        System.out.printf("r0= \n%s\n", FormatArray.toString(r0,"%.4e"));
        
        System.out.printf("chk==1: \n%s\n", FormatArray.toString(
            MatrixUtil.multiply(rFirst, r0),"%.4e"));
        
        // with orthographic, can only recover rotation, not translation
        //(2*mImages)X3
        //apply to indiv rotation matrices
         
        // r has the i and j direction and k=i cross j.
        // create a stack of rotation matrices, one per image.
        double[][] rotStack = MatrixUtil.zeros(3*mImages, 3);
        double[][] r3 = new double[2*mImages][];//(2*mImages)X3
        double[][] rTmp = new double[3][];
        //r2 size is (2*mImages)X3
        for (i = 0; i < mImages; ++i) {
            rTmp[0] = Arrays.copyOf(r2[i], r2[i].length);
            rTmp[1] = Arrays.copyOf(r2[mImages + i], r2[mImages + i].length);
            rTmp[2] = MatrixUtil.crossProduct(rTmp[0], rTmp[1]);
            rTmp = MatrixUtil.multiply(rTmp, r0);
            for (j = 0; j < 3; ++j) {
                rotStack[i*3 + j] = rTmp[j]; 
            }
            r3[i] = rTmp[0];
            r3[i + mImages] = rTmp[1];
        }
        
        double[][] shape = MatrixUtil.multiply(rFirst, s2);
        
        System.out.printf("**r3=\n%s\n", FormatArray.toString(r3, 
            "%.4e"));
        System.out.printf("**rot stack=\n%s\n", FormatArray.toString(rotStack, 
            "%.4e"));
        
        System.out.printf("shape=\n%s\n", FormatArray.toString(shape, 
            "%.4e"));
        
        System.out.printf("t=\n%s\n", FormatArray.toString(t, 
            "%.4e"));
        
        //assert eqn (3.7) of Tomasi & Kanade:
        // original measurement matrix 
        //    W = R*X + t*(e_p)^T 
        //        where t = the vector of centroids a_0, a_1,...a_(F-1),b_0,b_1,...b_(F-1)
        //         and e_P^T is a vector of P 1's.
        double[] ep = new double[nFeatures];
        Arrays.fill(ep, 1);
        double[][] tep = MatrixUtil.outerProduct(t, ep);
        System.out.printf("W = R*X + t*(e_p)^T=\n%s\n", 
            FormatArray.toString(MatrixUtil.elementwiseAdd(
                MatrixUtil.multiply(r3, shape), tep), 
            "%.4e"));
        System.out.printf("orig W = \n%s\n", FormatArray.toString(w, 
            "%.4e"));
        
        OrthographicProjectionResults results = new OrthographicProjectionResults();
        results.XW = shape;
        results.rotationMatrices = rotStack;
                
        return results;
    }
    
    /**
     * a look at enforcing orthonormal rotation
     * 
     * @param x
     * @param mImages
     * @return
     * @throws NotConvergedException 
     */
    static OrthographicProjectionResults _DoNotUseThisCalculateAffineReconstruction(
        double[][] x, int mImages) throws NotConvergedException {
                                
        if (x.length != 2) {
            throw new IllegalArgumentException("x.length must be 2");
        }
        if ((x[0].length % mImages) != 0) {
            throw new IllegalArgumentException("x must have a multiple of mImages as the number of columns");
        }
        int nFeatures = x[0].length / mImages;
        
        //2mn >= 8m + 3n – 12
        if ((2*mImages * nFeatures) < (8*mImages +3*nFeatures - 12)) {
            throw new IllegalArgumentException("more points are necessary:"
                + "2 * mImages * nFeatures >= 8 * mImages + 3 * nFeatures – 12."
                +  "\narguments: mImages=" + mImages + " nFeatures=" + nFeatures
                + "\n" + (2*mImages * nFeatures) + " < " +
                    (8*mImages +3*nFeatures - 12));
        }
        // for mImages=2, need 4 features
        
        // Tomasi & Kanade 1992: An image stream can be represented by the
        //  2*F X P measurment matrix W of image coordinates of P points
        // tracked through F frames.
        
        // create matrix W composed of U and V 
        //     where U is rows of each image's x coordinates (size is mImages X nFeatures).
        //     where V is rows of each image's y coordinates (size is mImages X nFeatures).
        // create matrix t which holds the centroids of each row of W
        // create matrix WC = W - t
        
        // row[0] = image_0_x[0], image_0_x[1], image_0_x[2], ...
        // row[1] = image_1_x[0], image_1_x[1], image_1_x[2], ...
        // ...
        // row[mImages+0] = image_0_y[0], image_0_y[1], image_0_y[2], ...
        //NOTE: adding the Z axis=1
        double[][] w = MatrixUtil.zeros(3*mImages, nFeatures);
        double[] t = new double[3*mImages];
        int i, j;
        int m, n, rowX, rowY, rowZ, col;
        for (m = 0; m < mImages; ++m) {
            rowX = m;
            rowY = mImages + m;
            rowZ = 2*mImages + m;
            for (n = 0; n < nFeatures; ++n) {
                col = nFeatures*m + n;
                w[rowX][n] = x[0][col];
                t[rowX] += x[0][col];
                w[rowY][n] = x[1][col];
                t[rowY] += x[1][col];
                w[rowZ][n] = 1;
                t[rowZ] += 1;
            }
            t[rowX] /= (double)nFeatures;
            t[rowY] /= (double)nFeatures;
            t[rowZ] /= (double)nFeatures;
        }
        
        //registered measurement matrix:
        double[][] wC = MatrixUtil.copy(w);
        for (i = 0; i < t.length; ++i) {
            for (n = 0; n < wC[i].length; ++n) {
                wC[i][n] -= t[i];
            }
        }
        
        // see Fig 3.1 of Tomasi & Kanade 1991 or Fig 2. of Belongie lecture notes
        
        // under orthography, this coordinate centered matrix has rank 3
        // and can be factored into the product of 2 matrices, R and S:
        // the camera rotation matrix R (size 2FX3) and 
        // the shape matrix S (size 3XP) which is shape in a coordinate system 
        // attached to the object centroid.
        //The two components of the camera translation along the image plane
        // are computed as averages of the rows of W
                
        // the registered measurement matrix is highly rank deficient
        // eqn (3.11)
        SVDProducts svd = MatrixUtil.performSVD(wC);
        
        // Tomasi & Kanade 1992,eqn (3.12): 1st 3 columns of U, upper 3X3 of S, and first 3 rows of V^T
        //U3 is 2FX3 where F is the number of image frames
        //S3 is 3X3
        //VT3 is 3XP where P is the number of points, that is features, per image
        double[][] u3 = MatrixUtil.copySubMatrix(svd.u, 0, svd.u.length-1, 0, 2);
        double[][] s3 = MatrixUtil.zeros(3, 3);
        s3[0][0] = svd.s[0];
        s3[1][1] = svd.s[1];
        s3[2][2] = svd.s[2];
        double[][] sqrts3 = MatrixUtil.copy(s3);
        sqrts3[0][0] = Math.sqrt(sqrts3[0][0]);
        sqrts3[1][1] = Math.sqrt(sqrts3[1][1]);
        sqrts3[2][2] = Math.sqrt(sqrts3[2][2]);
        double[][] vT3 = MatrixUtil.copySubMatrix(svd.vT, 0, 2, 0, svd.vT[0].length-1);
        
        // if the ratio between the 3rd and 4th largest singular value of the registered measurement matrix
        //   is large, then the noise portion of the full decomposition
        //   can be ignored (the noise portion is the block partitions not 
        //   copied to U3, sqrtS3 and vT3).
        // there is more about this in Morita & Kanade 1994/1997
        double sRatio = svd.s[2]/svd.s[3];
        System.out.printf("svd.s[2]/svd.s[3]=%.3e\n", sRatio);
        
        //Tamasi & Kanade
        // (2*mImages)X3
        double[][] rC = MatrixUtil.multiply(u3, sqrts3);
        // 3XnFeatures
        double[][] sC = MatrixUtil.multiply(sqrts3, vT3);
        
        System.out.printf("rC=\n%s\n", FormatArray.toString(rC, "%.4e"));
        
        System.out.printf("w=\n%s\n", FormatArray.toString(w, "%.4e"));
        System.out.printf("wC=\n%s\n", FormatArray.toString(wC, "%.4e"));
        System.out.printf("rC*sC=\n%s\n", FormatArray.toString(
            MatrixUtil.multiply(rC, sC), 
            "%.4e"));
        
        // wC = rC * sC  
        //rC and sC are linear translations of the true rotation matrix R and
        //  the true shape matrix S, respectively.
        // Morita and Kanade: the decomposition is not completely unique.  
        //     it's unique only up to an affine transformation
        
        /*
        NLK: the rotation matrices aren't orthormal, so one needs to apply a linear
        transformation to rC that makes it orthonormal while applying the
        inverse transformation to sC to maintain the value of wC.
        
        wC = rC * sC

        let rC' be a single rotation matrix formed from the x-row of a frame in rC,
            the y-row of the same frame in rC, and the cross product of the x and y rows.
        rCO is obtained from SVD(rC').U*(SVD(rC').VT)^T)

        let R be the orthogonal rCO matrices reformatted into the x, y row format of rC
        
        since wC = R * Z * z^-1 * sC
              wC = R * Z^-1 * sC
              pseudoInv(R)*wC = pseudoInv(R)*R * Z^-1 * sC
              pseudoInv(R)*wC * pseudoInv(sC) = Z^-1 if inv(R)*R=I
               
              then can use Z^-1 to transform sC into s
        
        (3) Can continue with the rest of the Tomasi & Kanade and Morita & Kanade
            algorithm, eqn (3.15) below...
        */
        
        //sC is 3XnFeatures
        
        double[][] invR;
        
        double[][] rot = new double[3*mImages][3];
        double[][] r3 = new double[2*mImages][3];
        double[][] shape = new double[3*mImages][nFeatures];
        
        // 3X3
        double[][] rCP = new double[3][];
        double[] tmp;
        //3X3
        double[][] rCO;
        SVDProducts svdRC;
        double[][] si;
        double[][] wCij = new double[3][];//3XnFeatures
        for (int ii = 0; ii < mImages; ++ii) {
            rCP[0] = rC[ii];
            rCP[1] = rC[ii + mImages];
            rCP[2] = rC[ii + 2*mImages];
            tmp = MatrixUtil.crossProduct(rCP[0], rCP[1]);
            System.out.printf("compare %s to %s\n", 
                FormatArray.toString(rCP[2], "%.3e"), FormatArray.toString(tmp, "%.3e"));
    
            System.out.printf("rC_%d=\n%s\n", ii, FormatArray.toString(rCP, "%.4e"));
            double detR = MatrixUtil.determinant(rCP);
            System.out.printf("det(rC_%d)=%.4e\n", ii, detR);
            
            svdRC = MatrixUtil.performSVD(rCP);
            rCO = MatrixUtil.multiply(svdRC.u, MatrixUtil.transpose(svdRC.vT));
            System.out.printf("r_uvtt=\n%s\n", FormatArray.toString(rCO, "%.4e"));
            detR = MatrixUtil.determinant(rCO);
            System.out.printf("det(r_uvtt_%d)=%.4e\n", ii, detR);
           
            System.out.printf("rC_%d * (rC_%d)^T=\n%s\n", ii, ii,
                FormatArray.toString(MatrixUtil.multiply(rCP, MatrixUtil.transpose(rCP)), "%.4e"));
            
            System.out.printf("r_uvtt_%d * (r_uvtt_%d)^T=\n%s\n", ii, ii,
                FormatArray.toString(MatrixUtil.multiply(rCO, MatrixUtil.transpose(rCO)), "%.4e"));
            
            System.out.flush();
           
            for (j = 0; j < 3; ++j) {
                rot[ii*3 + j] = rCO[j];
            }
            r3[ii] = Arrays.copyOf(rCO[0], rCO[0].length);
            r3[ii + mImages] = Arrays.copyOf(rCO[1], rCO[1].length);
            
            /*
            block[0] = rCP[0] times sC = 1XnFeatures
            block[1] = rCP[1] times sC
            block[2] = rCp[2] times sC 
            
            rC*z * (z^-1)*sC = 3X3*3X3*3XnFeatures = 3XnFeatures = block of wC
            rOrth * (z^-1)*sC= wC_i_j
            (z^-1)*sC= inv(rOth)*wC_i_j
            store in blocks of trans
            */     
            invR = MatrixUtil.pseudoinverseFullColumnRank(rCO);
            
            // 3XnFeatures
            wCij[0] = wC[ii];
            wCij[1] = wC[ii + mImages];
            wCij[2] = wC[ii + 2*mImages];
            
            //3XnFeatures
            si = MatrixUtil.multiply(invR, wCij);
            
            //double[][] trans = new double[3*mImages][nFeatures];
            for (j = 0; j < 3; ++j) {
                shape[ii*3 + j] = si[j];
            }
            
            // TODO: need to consider the rotation and translation origins now
            System.out.printf("for image%d have shape=\n%s\n", ii, FormatArray.toString(si, "%.4e"));
        }
        
        System.out.printf("rot stack=\n%s\n", FormatArray.toString(rot, 
            "%.4e"));
        
        System.out.printf("shape=\n%s\n", FormatArray.toString(shape, 
            "%.4e"));
        
        System.out.printf("t=\n%s\n", FormatArray.toString(t, 
            "%.4e"));
        // shape size is
        //r3 = new double[2*mImages][3];
        double[][] _si = new double[3][nFeatures];
        double[][] _ri = new double[2][3];
        double[][] _rs;
        for (i = 0; i < mImages; ++i) {
            _si[0] = shape[i*mImages + i];//3XnFeatures
            _si[1] = shape[i*mImages + i + 1];
            _si[2] = shape[i*mImages + i + 2]; 
            _ri[0] = r3[i*mImages + 0];//2X3
            _ri[1] = r3[i*mImages + 1];
            _rs = MatrixUtil.multiply(_ri, _si);//2XnFeatures
            for (int k = 0; k < nFeatures; ++k) {
                _rs[0][k] += t[i*mImages + i];
                _rs[1][k] += t[i*mImages + i + 1];
            }
            
            System.out.printf("W = R*X + t*(e_p)^T=\n%s\n", 
                FormatArray.toString(_rs, "%.4e"));
        }
        System.out.printf("original W = \n%s\n", 
            FormatArray.toString(w, "%.4e"));
        
        OrthographicProjectionResults results = new OrthographicProjectionResults();
        results.XW = shape;
        results.rotationMatrices = rot;
                
        return results;
           
    }
    
    /*
    TODO: proof read and write test for this.
    from Szeliski 2010 and Poelman & Kanade 1992 (year?  a few published translations with different years):
    Para-perspective provides a more accurate projection model than scaled 
    orthography, without incurring the added complexity of per-pixel perspective 
    division, which invalidates traditional factoriza- tion methods 
    (Poelman and Kanade 1997).
    
    Scaled orthographic projection, sometimes referred to as "weak perspective", 
    accounts for the scaling effect of an object as it moves towards and away 
    from the camera. Paraperspective projection, first introduced by Ohta in 
    [4] and named by Aloimonos in [1], accounts for the scaling effect as well 
    as the different angle from which an object is viewed as it moves in a 
    direction parallel to the image plane.

    */
          
    /**
     * for the case where the cameras are viewing small, distant scenes,
     * recover the 3-D coordinates in WCS and the projection matrices 
     * from pairs of corresponding
     * un-calibrated image points, that is, points in the image reference frame in pixels.
     * assumes a para-perspective camera model.
     *
     * Input set of P feature point coordinates (x_f_p,y_f_p) , for each of 
     * the F frames of the image sequence. From this information, our goal is 
     * to recover the estimated shape of the object, given by the position 
     * s_P, of every point, and the estimated motion of the
     camera, given by iHat_f, jHat_f, kHat_f for each frame in the sequence. 
     Rather than recover iHat_f in world coordinates, we generally recover 
     the three.eparate components tHat_f dot iHat_f, tHat_f dot jHat_f,
     tHat_f dot kHat_f.
 
      
     <pre>
      references:
      
     Poelman & Kanade 1997 (1994), "A Paraperspective Factorization Method for Shape 
     and Motion Recovery" 
     
     Description from Poelman & Kanade:
      
     Each feature point p that we track corresponds to a single world point, 
      located at position s. in some fixed world coordinate system.

      Each image f was taken at some specific camera orientation, which we 
      describe by the orthonormal unit vectors i_f, j_f and k_f 
      where kf_ points along the camera's line of sight, 
      i_f corresponds to the camera image plane's x-axis, 
      and j_f corresponds to the camera image's y-axis.
      
      t_f is a vector pointing from the origin of the fixed world coordinate system
      to the camera's focal plane.  it's the position of the camera in each fram f.
      
     </pre>
     * 
     * @param x the image coordinates of feature correspondences in 2 or more
     * images.  format is 2 X (nImages * nFeatures) where row 0 holds the x-coordinates
     * and row 1 holds the y-coordinates and each image's features are given
     * before the next and the features are ordered in the same manner within
     * all images.
     * for example: row 0 = [img_0_feature_0, ... img_0_feature_n-1, ... img_m-1_feature_0,...
     *     img_m-1_feature_n-1].
     * @param mImages the number of images in x.
     * @return the estimated projections P1 and P2 and the objects locations as 3-D points;
     */
    public static ParaperspectiveProjectionResults calculateParaperspectiveReconstruction(
        double[][] x, int mImages) throws NotConvergedException {
                        
        if (x.length != 2) {
            throw new IllegalArgumentException("x.length must be 2");
        }
        if ((x[0].length % mImages) != 0) {
            throw new IllegalArgumentException("x must have a multiple of mImages as the number of columns");
        }
        int nFeatures = x[0].length / mImages;
        
        ///2mn >= 8m + 3n – 12
        if ((2*mImages * nFeatures) < (8*mImages +3*nFeatures - 12)) {
            throw new IllegalArgumentException("more points are necessary:"
                + "2 * mImages * nFeatures >= 8 * mImages + 3 * nFeatures – 12."
                +  "\narguments: mImages=" + mImages + " nFeatures=" + nFeatures
                + "\n" + (2*mImages * nFeatures) + " < " +
                    (8*mImages +3*nFeatures - 12));
        }
        // for mImages=2, need 4 features
        
        //NOTE: assumes the camera's focal length = 1
        
        // create measurement matrix W composed of U and V 
        //     where U is rows of each image's x coordinates (size os mImages X nFeatures).
        //     where V is rows of each image's y coordinates (size os mImages X nFeatures).
        // create matrix t which holds the centroids of each row of W
        // create registered measurement matrix WC = W - t
        double[][] w = MatrixUtil.zeros(2*mImages, nFeatures);
        // t points to the camera's focal point
        double[] t = new double[2*mImages];
        int i, j;
        int m, n, rowU, rowV, col;
        for (m = 0; m < mImages; ++m) {
            rowU = m;
            rowV = mImages + m;
            for (n = 0; n < nFeatures; ++n) {
                col = nFeatures*m + n;
                w[rowU][n] = x[0][col];
                t[rowU] += x[0][col];
                w[rowV][n] = x[1][col];
                t[rowV] += x[1][col];
            }
            t[rowU] /= (double)nFeatures;
            t[rowV] /= (double)nFeatures;
        }
        
        //registered measurement matrix:
        double[][] wC = MatrixUtil.copy(w);
        for (i = 0; i < t.length; ++i) {
            for (n = 0; n < wC[i].length; ++n) {
                wC[i][n] -= t[i];
            }
        }
        
        SVDProducts svd = MatrixUtil.performSVD(wC);
        //U3 is 2FX3 where F is the number of image frames
        //S3 is 3X3
        //VT3 is 3XP where P is the number of points, that is features, per image
        double[][] u3 = MatrixUtil.copySubMatrix(svd.u, 0, svd.u.length-1, 0, 2);
        double[][] s3 = MatrixUtil.zeros(3, 3);
        s3[0][0] = svd.s[0];
        s3[1][1] = svd.s[1];
        s3[2][2] = svd.s[2];
        double[][] sqrts3 = MatrixUtil.copy(s3);
        sqrts3[0][0] = Math.sqrt(sqrts3[0][0]);
        sqrts3[1][1] = Math.sqrt(sqrts3[1][1]);
        sqrts3[2][2] = Math.sqrt(sqrts3[2][2]);
        double[][] vT3 = MatrixUtil.copySubMatrix(svd.vT, 0, 2, 0, svd.vT[0].length-1);
        
        // if the ratio between the 3rd and 4th largest singular value of the registered measurement matrix
        //   is large, then the noise portion of the full decomposition
        //   can be ignored (the noise portion is the block partitions not 
        //   copied to U3, sqrtS3 and vT3).
        // there is more about this in Morita & Kanade 1994/1997
        double sRatio = svd.s[2]/svd.s[3];
        System.out.printf("svd.s[2]/svd.s[3]=%.3e\n", sRatio);
                
        // (2*mImages)X3
        double[][] mC = MatrixUtil.multiply(u3, sqrts3);
        // 3XnFeatures
        double[][] sC = MatrixUtil.multiply(sqrts3, vT3);
        
        System.out.printf("wC=\n%s\n", FormatArray.toString(wC, "%.4e"));
        
        // assert that wC is the same as wC2
        System.out.printf("mC*sC=\n%s\n", FormatArray.toString(MatrixUtil.multiply(mC, sC), 
            "%.4e"));
        
        /*
        ------------------------------------------------------------------
        Paraperspective Normalization
        ------------------------------------------------------------------
        
         3 constraints:

         eqn(15) of paper:
                 |m_f|^2/(1+x_f^2) - |n_f|^2/(1+y_f^2) = 0
         eqn(17) of paper:
                 m_f dot n_f = x_f * y*f * 0.5 * ( |m_f|^2/(1+x_f^2) + |n_f|^2/(1+y_f^2) )
         eqn(18) of paper:
                 |m_0|=1

         those are 2*F + 1 equations as metric constraints

         from the SVD of the registered measurement matrix, there is M and S
         `M is size 2*F X 3
         `M = vectorized( m_0, m_1, ... m_{F-1}, n_0, n_1, ... n_{F-1},

         let M = `M*A where A is a 3X3 matrix, and as before, Q = symmetric matrix, but Q=A^T*A.

         Equations (15), (17), and (18) give us 2F+ 1 equations,
         We compute the 3 X 3 matrix A such that M = `M*A best satisfies these metric constraints
         in the least sum-of-squares error sense.

         This is a simple problem because the constraints are linear in the 6 unique elements
         of the symmetric 3 x 3 matrix Q = A^TA.
        
         Thus we compute Q by solving the overconstrained linear system of 2F + 1 equations
         in 6 variables defined by the metric constraints, ...

             [ q1  q2  q3 ]
         Q = [ q2  q4  q5 ]
             [ q3  q5  q6 ]

         m_f = `m_f * Q = [`mf0  `mf1  `mf2] * [ q1  q2  q3 ]
                                               [ q2  q4  q5 ]
                                               [ q3  q5  q6 ]
                        = [ (q1*`mf0 + q2*`mf1 + q3*`mf2 )  (q2*`mf0 + q4*`mf1 + q5*`mf2 )  (q3*`mf0 + q5*`mf1 + q6*`mf2 ) ]

         |vector| is the magnitude of a vector = square root of the sum of squares of its components.
        
         as an aside, in case can simplify any future steps with this:
         and Q*Q = q1*q1 + q2*q2 + q3*q3   q1*q2 + q2*q4 + q3*q5   q1*q3 + q2*q5 + q3*q6
                   q1*q2 + q2*q4 + q3*q5   q2*q2 + q4*q4 + q5*q5   q2*q3 + q4*q5 + q5*q6
                   q1*q3 + q2*q5 + q3*q6   q2*q3 + q4*q5 + q5*q6   q3*q3 + q5*q5 + q6*q6
                 = Q_col0 dot Q_col0   Q_col0 dot Q_col1  Q_col0 dot Q_col2
                   Q_col0 dot Q_col1   Q_col1 dot Q_col1  Q_col1 dot Q_col2
                   Q_col0 dot Q_col2   Q_col1 dot Q_col2  Q_col2 dot Q_col2

        NOTE: below are the expanded details of the multiplication.
        The result can be rewritten using other notation:
        since need the dot product of vector `m_f*Q with itself, can use the inner product
        of the vector multiplied by its transpose:
            |m_f|^2 = `m_f*Q*Q^T*`m_f^T where m_f is a 1X3 vector and Q is a 3X3 symmetric matrix
                    = `m_f*Q^2*`m_f^T
                    = [`m_f*Q^2[:][0] `m_f*Q^2[:][1] `m_f*Q^2[:][2]] * `m_f^T
                    = ['m_f[0]*`m_f*Q^2[:][0] + 'm_f[1]*`m_f*Q^2[:][1] + 'm_f[2]*`m_f*Q^2[:][2]]


         expand |m_f|^2/(1+x_f^2) :
             (1/(1+x_f^2)) * [ (q1*`mf0 + q2*`mf1 + q3*`mf2 )^2 + (q2*`mf0 + q4*`mf1 + q5*`mf2 )^2 + (q3*`mf0 + q5*`mf1 + q6*`mf2 )^2 ]
             (1/(1+x_f^2)) * [ (q1*`mf0*q1*`mf0 + q1*`mf0*q2*`mf1 + q1*`mf0*q3*`mf2)
                              + (q2*`mf1*q1*`mf0 + q2*`mf1*q2*`mf1 + q2*`mf1*q3*`mf2)
                              + (q3*`mf2*q1*`mf0 + q3*`mf2*q2*`mf1 + q3*`mf2*q3*`mf2)
                              + (q2*`mf0*q2*`mf0 + q2*`mf0*q4*`mf1 + q2*`mf0*q5*`mf2 )
                              + (q4*`mf1*q2*`mf0 + q4*`mf1*q4*`mf1 + q4*`mf1*q5*`mf2 )
                              + (q5*`mf2*q2*`mf0 + q5*`mf2*q4*`mf1 + q5*`mf2*q5*`mf2 )
                              + (q3*`mf0*q3*`mf0 + q3*`mf0*q5*`mf1 + q3*`mf0*q6*`mf2 )
                              + (q5*`mf1*q3*`mf0 + q5*`mf1*q5*`mf1 + q5*`mf1*q6*`mf2 )
                              + (q6*`mf2*q3*`mf0 + q6*`mf2*q5*`mf1 + q6*`mf2*q6*`mf2 ) ]

             (1/(1+x_f^2))  * [ (q1*`mf0*q1*`mf0 + q1*`mf0*q2*`mf1 + q1*`mf0*q3*`mf2)
                              + (q1*`mf0*q2*`mf1 + q2*`mf1*q2*`mf1 + q2*`mf1*q3*`mf2)
                              + (q1*`mf0*q3*`mf2 + q2*`mf1*q3*`mf2 + q3*`mf2*q3*`mf2)
                              + (q2*`mf0*q2*`mf0 + q2*`mf0*q4*`mf1 + q2*`mf0*q5*`mf2 )
                              + (q2*`mf0*q4*`mf1 + q4*`mf1*q4*`mf1 + q4*`mf1*q5*`mf2 )
                              + (q2*`mf0*q5*`mf2 + q4*`mf1*q5*`mf2 + q5*`mf2*q5*`mf2 )
                              + (q3*`mf0*q3*`mf0 + q3*`mf0*q5*`mf1 + q3*`mf0*q6*`mf2 )
                              + (q3*`mf0*q5*`mf1 + q5*`mf1*q5*`mf1 + q5*`mf1*q6*`mf2 )
                              + (q3*`mf0*q6*`mf2 + q6*`mf2*q5*`mf1 + q6*`mf2*q6*`mf2 ) ]
        
                          (1/(1+x_f^2))  * [ q1*q1*`mf0*`mf0 + q1*q2*`mf0*`mf1 + q1*q3*`mf0*`mf2
                              + q1*q2*`mf0*`mf1 + q2*q2*`mf1*`mf1 + q2*q3*`mf1*`mf2
                              + q1*q3*`mf0*`mf2 + q2*q3*`mf1*`mf2 + q3*q3*`mf2*`mf2
                              + q2*q2*`mf0*`mf0 + q2*q4*`mf0*`mf1 + q2*q5*`mf0*`mf2
                              + q2*q4*`mf0*`mf1 + q4*q4*`mf1*`mf1 + q4*q5*`mf1*`mf2
                              + q2*q5*`mf0*`mf2 + q4*q5*`mf1*`mf2 + q5*q5*`mf2*`mf2
                              + q3*q3*`mf0*`mf0 + q3*q5*`mf0*`mf1 + q3*q6*`mf0*`mf2
                              + q3*q5*`mf0*`mf1 + q5*q5*`mf1*`mf1 + q5*q6*`mf1*`mf2
                              + q3*q6*`mf0*`mf2 + q6*q5*`mf2*`mf1 + q6*q6*`mf2*`mf2 ]

             (1/(1+x_f^2))  * [ q1*q1*`mf0*`mf0 + q2*q2*`mf0*`mf0 + q3*q3*`mf0*`mf0
                              + q1*q2*`mf0*`mf1 + q2*q4*`mf0*`mf1 + q3*q5*`mf0*`mf1   <== 1st 3 lines of addition are:
                              + q1*q3*`mf0*`mf2 + q2*q5*`mf0*`mf2 + q3*q6*`mf0*`mf2       [ `mf0*`mf0  `mf0*`mf1  `mf0*`mf2 ] * [Q^2_col0]

                              + q1*q2*`mf0*`mf1 + q2*q4*`mf0*`mf1 + q3*q5*`mf0*`mf1   <== 2nd 3 lines of addition are:
                              + q2*q2*`mf1*`mf1 + q4*q4*`mf1*`mf1 + q5*q5*`mf1*`mf1       [ `mf0*`mf1  `mf1*`mf1  `mf1*`mf2 ] * [Q^2_col1]
                              + q2*q3*`mf1*`mf2 + q4*q5*`mf1*`mf2 + q5*q6*`mf1*`mf2

                              + q1*q3*`mf0*`mf2 + q2*q5*`mf0*`mf2 + q3*q6*`mf0*`mf2   <== 3rd 3 lines of addition are:
                              + q2*q3*`mf1*`mf2 + q4*q5*`mf1*`mf2 + q5*q6*`mf1*`mf2       [ `mf0*`mf2  `mf1*`mf2  `mf2*`mf2 ] * [Q^2_col2]
                              + q3*q3*`mf2*`mf2 + q5*q5*`mf2*`mf2 + q6*q6*`mf2*`mf2 ]

             let z0 = [ `mf0*`mf0  `mf0*`mf1  `mf0*`mf2 ]
                 z1 = [ `mf0*`mf1  `mf1*`mf1  `mf1*`mf2 ]
                 z2 = [ `mf0*`mf2  `mf1*`mf2  `mf2*`mf2 ]
             (1/(1+x_f^2))  * [ z0 z1 z2] * [Q^2_col0]
                                            [Q^2_col1]
                                            [Q^2_col2]

         reminder of Q*Q:
                 = q1*q1 + q2*q2 + q3*q3   q1*q2 + q2*q4 + q3*q5   q1*q3 + q2*q5 + q3*q6
                   q1*q2 + q2*q4 + q3*q5   q2*q2 + q4*q4 + q5*q5   q2*q3 + q4*q5 + q5*q6
                   q1*q3 + q2*q5 + q3*q6   q2*q3 + q4*q5 + q5*q6   q3*q3 + q5*q5 + q6*q6
                 = Q_col0 dot Q_col0   Q_col0 dot Q_col1  Q_col0 dot Q_col2
                   Q_col0 dot Q_col1   Q_col1 dot Q_col1  Q_col1 dot Q_col2
                   Q_col0 dot Q_col2   Q_col1 dot Q_col2  Q_col2 dot Q_col2
            
             let z0 = [ `mf0*`mf0  `mf0*`mf1  `mf0*`mf2 ]
                 z1 = [ `mf1*`mf0  `mf1*`mf1  `mf1*`mf2 ]
                 z2 = [ `mf2*`mf0  `mf2*`mf1  `mf2*`mf2 ]
             let y0 = [ `nf0*`nf0  `nf0*`nf1  `nf0*`nf2 ]
                 y1 = [ `nf1*`nf0  `nf1*`nf1  `nf1*`nf2 ]
                 y2 = [ `nf2*`nf0  `nf2*`nf1  `nf2*`nf2 ]
             let w0 = [ `mf0*`nf0  `mf0*`nf1  `mf0*`nf2 ]
                 w1 = [ `mf1*`nf0  `mf1*`nf1  `mf1*`nf2 ]
                 w2 = [ `mf2*`nf0  `mf2*`nf1  `mf2*`nf2 ]
             let c2 = (x_f*y_f*0.5)

        factoring the constraints to separate Q unknowns from m and n knowns:
         eqn(15) of paper:
                 |m_f|^2/(1+x_f^2) - |n_f|^2/(1+y_f^2) = 0

             (1/(1+x_f^2)) * [ z0 z1 z2] * [Q^2_col0]  -  (1/(1+y_f^2)) * [ y0 y1 y2] * [Q^2_col0] = 0
                                           [Q^2_col1]                                   [Q^2_col1]
                                           [Q^2_col2]                                   [Q^2_col2]

             [ (z0/(1+x_f^2) -y0/(1+y_f^2))  (z1/(1+x_f^2) -y1/(1+y_f^2))  (z2/(1+x_f^2) -y2/(1+y_f^2)) ] * [Q^2_col0] = 0
                                                                                                            [Q^2_col1]
                                                                                                            [Q^2_col2]

         eqn(17) of paper:
                 m_f dot n_f = x_f*y_f*0.5 * ( |m_f|^2/(1+x_f^2) + |n_f|^2/(1+y_f^2) )

                 [ w0 w1 w2] * [Q^2_col0] - c2 * [ z0 z1 z2] * [Q^2_col0]  -  c2 * [ y0 y1 y2] * [Q^2_col0] = 0
                               [Q^2_col1]                      [Q^2_col1]                        [Q^2_col1]
                               [Q^2_col2]                      [Q^2_col2]                        [Q^2_col2]

                 [ (w0 -z0*c2 -y0*c2)  (w1 -z1*c2 -y1*c2)  (w2 -z2*c2 -y2*c2)] * [Q^2_col0] = 0
                                                                                 [Q^2_col1]
                                                                                 [Q^2_col2]
         eqn(18) of paper:
                 |m_0|=1

                 square to use the same factorization by Q^2?

             let v0 = [ `m_0[0]*`mf0  `m_0[0]*`m_0[1]  `m_0[0]*`m_0[2] ]
                 v1 = [ `m_0[0]*`mf1  `m_0[1]*`m_0[1]  `m_0[1]*`m_0[2] ]
                 v2 = [ `m_0[0]*`mf2  `m_0[1]*`m_0[2]  `m_0[2]*`m_0[2] ]
             |m_0|^2 = [ v0 v1 v2] * [Q^2_col0] = 1
                                     [Q^2_col1]
                                     [Q^2_col2]
        */
        
        double[][] g = new double[2*mImages + 1][9];
        double[] z0 = new double[3];
        double[] z1 = new double[3];
        double[] z2 = new double[3];
        double[] y0 = new double[3];
        double[] y1 = new double[3];
        double[] y2 = new double[3];
        double[] w0 = new double[3];
        double[] w1 = new double[3];
        double[] w2 = new double[3];
        double c2, xf, yf, divXf, divYf;
        
        for (i = 0; i < mImages; ++i) {
            xf = t[i];
            yf = t[mImages + i];
            
            z0[0] = mC[i][0] * mC[i][0]; 
            z0[1] = mC[i][0] * mC[i][1];
            z0[2] = mC[i][0] * mC[i][2];
            
            z1[0] = mC[i][1] * mC[i][0]; 
            z1[1] = mC[i][1] * mC[i][1];
            z1[2] = mC[i][1] * mC[i][2];
            
            z2[0] = mC[i][2] * mC[i][0]; 
            z2[1] = mC[i][2] * mC[i][1];
            z2[2] = mC[i][2] * mC[i][2];
            
            y0[0] = mC[mImages + i][0] * mC[mImages + i][0];
            y0[1] = mC[mImages + i][0] * mC[mImages + i][1];
            y0[2] = mC[mImages + i][0] * mC[mImages + i][2];
            
            y1[0] = mC[mImages + i][1] * mC[mImages + i][0];
            y1[1] = mC[mImages + i][1] * mC[mImages + i][1];
            y1[2] = mC[mImages + i][1] * mC[mImages + i][2];
            
            y2[0] = mC[mImages + i][2] * mC[mImages + i][0];
            y2[1] = mC[mImages + i][2] * mC[mImages + i][1];
            y2[2] = mC[mImages + i][2] * mC[mImages + i][2];
            
            w0[0] = mC[i][0] * mC[mImages + i][0];
            w0[1] = mC[i][0] * mC[mImages + i][1];
            w0[2] = mC[i][0] * mC[mImages + i][2];
            
            w1[0] = mC[i][1] * mC[mImages + i][0];
            w1[1] = mC[i][1] * mC[mImages + i][1];
            w1[2] = mC[i][1] * mC[mImages + i][2];
            
            w2[0] = mC[i][2] * mC[mImages + i][0];
            w2[1] = mC[i][2] * mC[mImages + i][1];
            w2[2] = mC[i][2] * mC[mImages + i][2];
            
            c2 = 0.5*xf*yf;
            
            divXf = 1./(1.+xf*xf);
            divYf = 1./(1.+yf*yf);
            
            // length 9
            // eqn(15) of paper is the 1st mImages rows of g
            //   [ (z0/(1+x_f^2) -y0/(1+y_f^2))  (z1/(1+x_f^2) -y1/(1+y_f^2))  (z2/(1+x_f^2) -y2/(1+y_f^2)) ] ... = 0
            g[i] = new double[9]; // g is (2*mImages + 1) X 9;
            for (j = 0; j < 3; ++j) { 
                g[i][j] = (z0[j]*divXf - y0[j]*divYf);                
                g[i][3+j] = (z1[j]*divXf - y1[j]*divYf);
                g[i][6+j] = (z2[j]*divXf - y2[j]*divYf);                
            }
            
            //eqn(17) of paper is the 2nd mImages rows of g
            //  [ (w0 -z0*c2 -y0*c2)  (w1 -z1*c2 -y1*c2)  (w2 -z2*c2 -y2*c2)] ... = 0
            g[mImages + i] = new double[9]; // g is (2*mImages + 1) X 9;
            for (j = 0; j < 3; ++j) { 
                g[mImages + i][j] = (w0[j] - z0[j]*c2 - y0[j]*c2);                
                g[mImages + i][3+j] = (w1[j] - z1[j]*c2 - y1[j]*c2);
                g[mImages + i][6+j] = (w2[j] - z2[j]*c2 - y2[j]*c2);               
            }
        }
        
        //let v0 = [ `m_0[0]*`mf0  `m_0[0]*`m_0[1]  `m_0[0]*`m_0[2] ]
        //    v1 = [ `m_0[0]*`mf1  `m_0[1]*`m_0[1]  `m_0[1]*`m_0[2] ]
        //    v2 = [ `m_0[0]*`mf2  `m_0[1]*`m_0[2]  `m_0[2]*`m_0[2] ]
        // length 9
        // [ v0 v1 v2 ] = 1
        i = 0;
        double[] v0 = new double[] {
           mC[i][0] * mC[i][0], mC[i][0] * mC[i][1], mC[i][0] * mC[i][2]};
        double[] v1 = new double[] {
            mC[i][1] * mC[i][0], mC[i][1] * mC[i][1], mC[i][1] * mC[i][2]};
        double[] v2 = new double[] {
            mC[i][2] * mC[i][0], mC[i][2] * mC[i][1], mC[i][2] * mC[i][2]};
        g[2*mImages] = new double[9]; // g is (2*mImages + 1) X 9;
        for (j = 0; j < 3; ++j) { 
            g[2*mImages][j] = v0[j];                
            g[2*mImages][3+j] = v1[j];
            g[2*mImages][6+j] = v2[j];               
        }
        
        double[] c = new double[2*mImages + 1];
        c[2*mImages] = 1;
        
        double[][] gInv = MatrixUtil.pseudoinverseFullColumnRank(g);
        double[] iVector = MatrixUtil.multiplyMatrixByColumnVector(gInv, c);
        assert(iVector.length == 9);
        
        // 3X3
        double[][] ell = new double[3][3];
        ell[0] = new double[]{iVector[0], iVector[1], iVector[2]};
        ell[1] = new double[]{iVector[1], iVector[3], iVector[4]};
        ell[2] = new double[]{iVector[2], iVector[4], iVector[5]};

        // Q can be determined :
        //   as the square root of ell,
        //   or with the Cholesky decomposition
        //   or with eigendecomposition

        double eps = 1e-16;//1.e-11; eps close to zero within machine precision to perturb the matrix to smallest eigenvalue of eps
        double[][] lPSD = MatrixUtil.nearestPositiveSemidefiniteToASymmetric(ell, eps);
        EVD evd2 = EVD.factorize(new DenseMatrix(lPSD));
        double[] eig = evd2.getRealEigenvalues();
        System.out.printf("eig(L_PSD)=\n%s\n", FormatArray.toString(eig, "%.4e"));
        double[][] aMinusPSD = MatrixUtil.elementwiseSubtract(ell, lPSD);
        double dist1 = MatrixUtil.frobeniusNorm(aMinusPSD);
        
        boolean ipd = MatrixUtil.isPositiveDefinite(lPSD);
        assert(ipd);
        
        //decompose Q = L * (sigma+) * L^T;  Q is size 3X3
        double[][] q = LinearEquations.choleskyDecompositionViaMTJ(lPSD);
        
        // (2*mImages)X3
        double[][] _M = MatrixUtil.multiply(mC, q);
        // 3XnFeatures
        double[][] _S = MatrixUtil.multiply(MatrixUtil.pseudoinverseRankDeficient(q), sC);
        
        System.out.printf("_M=\n%s\n", FormatArray.toString(_M, "%.4e"));
        
        // assert that wC is the same as _M*_S
        System.out.printf("_M*_S=\n%s\n", FormatArray.toString(MatrixUtil.multiply(_M, _S), 
            "%.4e"));
                
        /*
        --------------------------------
        Paraperspective Motion Recovery
        --------------------------------
        eqn(19) :
         `i_f = z_f*m_f + x_f*`k_f

         `j_f = z_f*n_f + y_f*`k_f


         Since the `i, `j, `k  produced must be orthonormal,
            they can be written as functions of only 3 rotational variables.
            We can then view the problem as, for each frame f,
              solving an overconstrained system of 6 equations
              (the expansion of (19) to each of its vector components)
              in 4 variables (the 3 rotational variables and zs).

              `i_f[0] = z_f * m_f[0] + x_f * `k_f[0]
              `i_f[1] = z_f * m_f[1] + x_f * `k_f[1]
              `i_f[2] = z_f * m_f[2] + x_f * `k_f[2]

              `j_f[0] = z_f * n_f[0] + y_f * `k_f[0]
              `j_f[1] = z_f * n_f[1] + y_f * `k_f[1]
              `j_f[2] = z_f * n_f[2] + y_f * `k_f[2]

          using the equalities of `k_f terms from eqn (19):
            `k_f[0]: (same for vector components [1] and [2]):

                (1/x_f)*(`i_f[0] - z_f * m_f[0]) = `k_f[0]

                (1/y_f)*(`j_f[0] - z_f * n_f[0]) = `k_f[0]

                (1/x_f)*(`i_f[0] - z_f * m_f[0]) = (1/y_f)*(`j_f[0] - z_f * n_f[0])
                `i_f[0] = z_f * m_f[0] + (x_f/y_f)*(`j_f[0] - z_f * n_f[0])
                        = `j_f[0]*(x_f/y_f) + z_f * m_f[0] - (x_f/y_f)*(z_f * n_f[0])
               generalized for each component::
                `i_f = `j_f*(x_f/y_f) + z_f*m_f - (x_f/y_f)*(z_f*n_f)
        
         rewriting `i_f, `j_f, and `k_f in terms of `j_f, m_f, n_f:
            `i_f = `j_f*(x_f/y_f) + z_f*m_f - (x_f/y_f)*(z_f*n_f)
            `j_f
            `k_f = `j_f*(1/y_f) - z_f*n_f*(1/y_f)

           dot product metrics:
            idoti: `i_f dot `i_f = 1
            jdotj: `j_f dot `j_f = 1
            idotj: `i_f dot `j_f = 0

            idoti: `i_f dot `i_f = 1:
                (`j_f dot `j_f)*(x_f/y_f)^2 + (m_f dot m_f)*(z_f)^2 - (n_f dot n_f)*(z_f)^2*(x_f/y_f)^2 = 1
                (1)*(x_f/y_f)^2 + (m_f dot m_f)*(z_f)^2 - (n_f dot n_f)*(z_f)^2*(x_f/y_f)^2 = 1
                (z_f)^2*((m_f dot m_f)-(n_f dot n_f)*(x_f/y_f)^2) + (x_f/y_f)^2 = 1
                (z_f)^2*((m_f[0]^2 + m_f[1]^2 + m_f[2]^2)-(n_f[0]^2 + n_f[1]^2 + n_f[2]^2)*(x_f/y_f)^2) + (x_f/y_f)^2 = 1
                (z_f)^2 = (1-(x_f/y_f)^2) / ((m_f[0]^2 + m_f[1]^2 + m_f[2]^2)-(n_f[0]^2 + n_f[1]^2 + n_f[2]^2)*(x_f/y_f)^2)
*          ====> can solve for (z_f)^2 from this
           ==>   Instead, using eqn (15) averages: 
                  zf^2 = (1/2)*( ((1+xf^2)/(mf*mf)) + ((1+yf^2)/(nf*nf)))

            jdotj:`j_f dot `j_f = 1:

            idotj:`i_f dot `j_f = 0:
                `j_f*(x_f/y_f) + z_f*m_f - (x_f/y_f)*(z_f*n_f) dot j_f = 0
                (`j_f dot `j_f)*(x_f/y_f) + (j_f dot m_f)*z_f - (j_f dot nf)*(z_f*x_f/y_f) = 0
                 (j_f dot m_f)*z_f - (j_f dot n_f)*(z_f*x_f/y_f) + (x_f/y_f) = 0
               factor this out for each component:
                (j_f[0]*m_f[0])*z_f - (j_f[0]*n_f[0])*(z_f*x_f/y_f) = -(x_f/y_f)
                j_f[0]*(m_f[0]*z_f - (n_f[0]*z_f*x_f/y_f)) = -(x_f/y_f)
                j_f[0] = -(x_f/y_f)/ (m_f[0]*z_f - (n_f[0]*z_f*x_f/y_f))
*            ===> can solve for j_f from this

            rewrite `i_f from eqn (19) using the equalities of `k_f terms from above:
               `i_f[0] = `j_f[0]*(x_f/y_f) + m_f[0]*z_f - n_f[0]*(z_f*x_f/y_f)
               `i_f[1] = `j_f[1]*(x_f/y_f) + m_f[1]*z_f - n_f[1]*(z_f*x_f/y_f)
               `i_f[2] = `j_f[2]*(x_f/y_f) + m_f[2]*z_f - n_f[2]*(z_f*x_f/y_f)
*            ===> can solve for `i_f from this

          using orthogonality of `i_f,`j_f,`k_f to define `k_f:
               `k_f = `i_f cross `j_f
               `k_f[0] = `i_f[1]*`j_f[2] - `i_f[2]*`j_f[1]
               `k_f[1] = `i_f[2]*`j_f[0] - `i_f[0]*`j_f[2]
               `k_f[2] = `i_f[0]*`j_f[1] - `i_f[1]*`j_f[0]
*            ===> can solve for `k_f from this
              uneeded details:
               `k_f[0] = `i_f[1]*`j_f[2] - `i_f[2]*`j_f[1]
                       = `j_f[2]*(`j_f[1]*(x_f/y_f) + m_f[1]*z_f - n_f[1]*z_f*(x_f/y_f))
                          - `j_f[1]*(`j_f[2]*(x_f/y_f) + m_f[2]*z_f - n_f[2]*z_f*(x_f/y_f))
               `k_f[1] = `i_f[2]*`j_f[0] - `i_f[0]*`j_f[2]
                       = `j_f[0]*(`j_f[2]*(x_f/y_f) + m_f[2]*z_f - n_f[2]*z_f*(x_f/y_f))
                          - `j_f[2]*(`j_f[0]*(x_f/y_f) + m_f[0]*z_f - n_f[0]*(z_f*x_f/y_f))
               `k_f[2] = `i_f[0]*`j_f[1] - `i_f[1]*`j_f[0]
                       = `j_f[1]*(`j_f[0]*(x_f/y_f) + m_f[0]*z_f - n_f[0]*(z_f*x_f/y_f))
                          - `j_f[0]*(`j_f[1]*(x_f/y_f) + m_f[1]*z_f - n_f[1]*z_f*(x_f/y_f))
    */    
        
        /*
        solve for z_f
           zf^2 = (1/2)*( ((1+xf^2)/(mf*mf)) + ((1+yf^2)/(nf*nf)))
        
        solve for `j_f:
            j_f[0] = -(x_f/y_f)/ (m_f[0]*z_f - (n_f[0]*z_f*x_f/y_f))
                  ...
        
        solve for `i_f:
            `i_f[0] = `j_f[0]*(x_f/y_f) + m_f[0]*z_f - n_f[0]*(z_f*x_f/y_f)
             i_f[1] = `j_f[1]*(x_f/y_f) + m_f[1]*z_f - n_f[1]*(z_f*x_f/y_f)
             i_f[2] = `j_f[2]*(x_f/y_f) + m_f[2]*z_f - n_f[2]*(z_f*x_f/y_f)
        
        solve for k_f:
            `k_f[0] = `i_f[1]*`j_f[2] - `i_f[2]*`j_f[1]
            `k_f[1] = `i_f[2]*`j_f[0] - `i_f[0]*`j_f[2]
            `k_f[2] = `i_f[0]*`j_f[1] - `i_f[1]*`j_f[0]
        */
                
        double[] zf = new double[mImages];
        double[] mf, nf, tmp1, tmp2, tmp3;
        double xDivY, mDotM, nDotN, tmp4, tmp5;
        double[][] _MCameraOrientation2D = new double[2*mImages][3];
        for (i = 0; i < mImages; ++i) {
            xf = t[i];
            yf = t[mImages + i];
            xDivY = xf/yf;
            mf = _M[i];
            nf = _M[mImages + i];
            mDotM = 0;
            nDotN = 0;
            for (j=0; j<mf.length; ++j) {
                mDotM += mf[i] * mf[i];
                nDotN += nf[i] * nf[i];
            }
            //for z_f^2 take an average of eqn (15).  it's always positive
            tmp4 = 1 + xf * xf;
            tmp4 /= mDotM;
            tmp5 = 1 + yf * yf;
            tmp5 /= nDotN;
            zf[i] = (tmp4 + tmp5) / 2.;
            zf[i] = Math.sqrt(zf[i]);

            //j_f[0] = -(x_f/y_f)/ (m_f[0]*z_f - (n_f[0]*z_f*x_f/y_f))
            tmp1 = Arrays.copyOf(mf, mf.length);
            MatrixUtil.multiply(tmp1, zf[i]);
            tmp2 = Arrays.copyOf(nf, nf.length);
            MatrixUtil.multiply(tmp2, zf[i]);
            MatrixUtil.multiply(tmp2, xDivY);
            // j_f
            _MCameraOrientation2D[mImages + i] = new double[3];
            for (j = 0; j < 3; ++j) {
                _MCameraOrientation2D[mImages + i][j] = -xDivY/(tmp1[j] - (tmp2[j]/yf));
            }
            
            //`i_f = `j_f*(x_f/y_f) + m_f*z_f - n_f*(z_f*x_f/y_f)
            tmp1 = Arrays.copyOf(_MCameraOrientation2D[mImages + i], _MCameraOrientation2D[mImages + i].length);
            MatrixUtil.multiply(tmp1, xDivY);
            tmp2 = Arrays.copyOf(mf, mf.length);
            MatrixUtil.multiply(tmp2, zf[i]);
            tmp3 = Arrays.copyOf(nf, nf.length);
            MatrixUtil.multiply(tmp3, zf[i]);
            MatrixUtil.multiply(tmp3, xDivY);
            // i_f
            _MCameraOrientation2D[i] = new double[3];
            for (j = 0; j < 3; ++j) {
                _MCameraOrientation2D[i][j] = tmp1[j] + tmp2[j] - tmp3[j];
            }
        }
        
        System.out.printf("_MCameraOrientation2D=\n%s\n", FormatArray.toString(_MCameraOrientation2D, "%.4e"));
        
        // _M: rank is 3, mRows is >= 4, nCols=3
        // _MCameraOrientation2D: rank is 3, mRows is >= 4, nCols=3

        // _M * A = _MCameraOrientation2D
        // A = pseudoInv(_M)*_MCameraOrientation2D;
        // A^-1 = pseudoInv(_MCameraOrientation2D)*_M 
        
        double[][] q2 = solveForTransformationToOrthoNormal(_MCameraOrientation2D);
                
        _MCameraOrientation2D = MatrixUtil.multiply(_MCameraOrientation2D, q2);
        
        assertDotProductMetrics(_MCameraOrientation2D, mImages);
                
        /* eqn (3)
           u_f_p = m_f dot s_p + x_f
           and 
           v_f_p = n_f dot s_p + y_f
              where u_f_p and v_f_p are the original feature coordinates per image
                 in the measurement matrix
              where x_f and y_f are each image's centroid of all features.
           the registered measurement matrix holds
             `u_f_p = u_f_p - x_f
             `v_f_p = v_f_p - y_f
        
        u_f_p = m_f dot s_p + x_f
        `u_f_p = m_f dot s_p  <== decomposition gives motion and shape
         u_f_p/m_f[0] = s_p
        */
        
        //Looking at overall transformation to get to camera orientation
        //    and applying inverse to _S before inv(q2)
        
        double[][] _S2 = MatrixUtil.multiply(MatrixUtil.pseudoinverseRankDeficient(q2), _S);
        
        System.out.printf("after orthonormalization: _MCameraOrientation2D=\n%s\n", 
            FormatArray.toString(_MCameraOrientation2D, "%.4e"));
        
        System.out.printf("_S2=inv(q2)*_S=\n%s\n", FormatArray.toString(_S2, "%.4e"));
        
        System.out.printf("_MCameraOrientation2D*_S2=\n%s\n", 
            FormatArray.toString(MatrixUtil.multiply(_MCameraOrientation2D, _S2), "%.4e"));
                  
        /*
        from Tomasi & Kanade 1992:
        If desired, align the first camera reference system with the world
        reference system by forming the products R*R_0 and R_0^T*S,
        where the orthonormal matrix R_0 = [i1 j1 k1] rotates the ^Lfirst camera
        reference system into the identity matrix

        i0x i0y i0z    *  r0ix  r0jx  r0kx   = (i0 dot r0i)  (i0 dot r0j) (i0 dot r0k)
        i1x i1y i1z       r0iy  r0jy  r0ky
        i2x i2y i2z       r0iz  r0jz  r0kz
        i3x i3y i3z
        ...
        j0x j0y j0z
        j1x j1y j1z
        */
        

        double[][] rFirst = new double[3][3];
        rFirst[0] = Arrays.copyOf(_MCameraOrientation2D[0], _MCameraOrientation2D[0].length);
        rFirst[1] = Arrays.copyOf(_MCameraOrientation2D[mImages], _MCameraOrientation2D[mImages].length);
        rFirst[2] = MatrixUtil.crossProduct(rFirst[0], rFirst[1]);
        
        double[][] r0 = MatrixUtil.pseudoinverseFullColumnRank(rFirst);
        
        System.out.printf("r0= \n%s\n", FormatArray.toString(r0,"%.4e"));
        
        System.out.printf("chk==1: \n%s\n", FormatArray.toString(
            MatrixUtil.multiply(rFirst, r0),"%.4e"));
                  
        // multiply _M2 by r0
        double[][] rotStack = MatrixUtil.zeros(3*mImages, 3);
        double[][] _M3 = new double[3*mImages][];//(2*mImages)X3
        double[][] rTmp = new double[3][];
        for (i = 0; i < mImages; ++i) {
            rTmp[0] = Arrays.copyOf(_MCameraOrientation2D[i], _MCameraOrientation2D[i].length);
            rTmp[1] = Arrays.copyOf(_MCameraOrientation2D[mImages + i], _MCameraOrientation2D[mImages + i].length);
            rTmp[2] = MatrixUtil.crossProduct(rTmp[0], rTmp[1]);
            rTmp = MatrixUtil.multiply(rTmp, r0);
            for (j = 0; j < 3; ++j) {
                rotStack[i*3 + j] = rTmp[j]; 
            }
            _M3[i] = rTmp[0];
            _M3[i + mImages] = rTmp[1];
            _M3[i + 2*mImages] = rTmp[2];
        }
        
        double[][] _S3 = MatrixUtil.multiply(rFirst, _S2);
        
        System.out.printf("_M3=\n%s\n", FormatArray.toString(_M3,  "%.4e"));
        System.out.printf("rot stack=\n%s\n", FormatArray.toString(rotStack,  "%.4e"));
        System.out.printf("shape=\n%s\n", FormatArray.toString(_S3, "%.4e"));
        System.out.printf("t=\n%s\n", FormatArray.toString(t, 
            "%.4e"));
        /*
        Poelman & Kanade, last paragraph, Sect 3.4:
        All that remain to be computed are the translations for each frame. 
        We calculate the depth z_f from (15).  
        Once we know we x_f, y_f, z_f, `i_f, `j_f, `k_f 
        we can calculate `t_f using (4) and (5).
        
        eqn(4):
            z_f = -t_f dot k_f
        eqn (5):
            x_f = (-1/z_f)*(t_f dot i_f)
            and
            y_f = (-1/z_f)*(t_f dot j_f)

        use z_f equalities:
            z_f = -t_f dot k_f
                = -tf[0]*kf[0] + -tf[1]*kf[1] + -tf[2]*kf[2]

            z_f = (-1/x_f) * (t_f dot i_f)
                = -(1/xf)*tf[0]*if[0] + -(1/xf)*tf[1]*if[1] + -(1/xf)*tf[2]*if[2]

            z_f = (-1/y_f) * (t_f dot j_f)
                = -(1/yf)*tf[0]*jf[0] + -(1/yf)*tf[1]*jf[1] + -(1/yf)*tf[2]*jf[2]
            
        factor:
           tf[0]           tf[1]           tf[2]         const
           -------------------------------------------
            kf[0]           kf[1]          kf[2]         -zf
            (1/xf)*if[0]    (1/xf)*if[1]   (1/xf)*if[2]  -zf
            (1/yf)*jf[0]    (1/yf)*jf[1]   (1/yf)*jf[2]  -zf
        */
        double[][] trans = new double[mImages][3];
        double[][] tf = MatrixUtil.zeros(3, 3);
        double[][] g3Inv;
        double[] i3Vector;
        double[] cs = new double[3];
        for (i = 0; i < mImages; ++i) {
            xf = t[i];
            yf = t[mImages + i];
            
            // i_f[i] is _M2[i]
            // j_f[i] is _M2[mImages + i]
            // k_f[i] is _M2[2*mImages + i]
            for (j = 0; j < 3; ++j) {
                tf[0][j] = _M3[2*mImages + i][j];
                tf[1][j] = _M3[i][j]/xf;
                tf[2][j] = _M3[mImages + i][j]/yf;
            }
            
            Arrays.fill(cs, -zf[i]);
            SVDProducts svdi = MatrixUtil.performSVD(tf);
            g3Inv = MatrixUtil.pseudoinverseRankDeficient(tf);
            i3Vector = MatrixUtil.multiplyMatrixByColumnVector(g3Inv, cs);
            assert(i3Vector.length == 3);
            
            trans[i] = new double[]{i3Vector[0], i3Vector[1], i3Vector[2]};
        }
        
        System.out.printf("trans=\n%s\n", FormatArray.toString(trans, "%.4e"));
        
        ParaperspectiveProjectionResults results = new ParaperspectiveProjectionResults();
        results.XW = _S;
        results.rotationMatrices = rotStack;
        results.translationVectors = trans;
          
        return results;        
    }
    
    private static DenseMatrix extractIndices(DenseMatrix m, List<Integer> inlierIndexes) {
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

    private static double[] gT(double[] a, double[] b) {
        double[] gT = new double[]{
            a[0]*b[0], 
            a[0]*b[1] + a[1]*b[0],
            a[0]*b[2] + a[2]*b[0],
            a[1]*b[1],
            a[1]*b[2] + a[2]*b[1],
            a[2]*b[2]
        };
        return gT;
    }

    private static double[][] extractAndNormalize(double[][] x, int imageNumber, 
        int nFeatures, double[] outputNorm) {
        
        double[][] xN = new double[3][nFeatures];
        xN[0] = new double[nFeatures];
        xN[1] = new double[nFeatures];
        xN[2] = new double[nFeatures];
        Arrays.fill(xN[2], 1);
        
        int imageIdx = (x.length/nFeatures)*imageNumber;
        
        double cen0 = 0;
        double cen1 = 0;
        for (int i = 0; i < nFeatures; ++i) {
            cen0 += x[0][imageIdx + i];
            cen1 += x[0][imageIdx + i];
        }
        cen0 /= (double)nFeatures;
        cen1 /= (double)nFeatures;

        double scale = 0;
        double diffX, diffY;
        for (int i = 0; i < nFeatures; ++i) {
            diffX = x[0][imageIdx + i] - cen0;
            diffY = x[0][imageIdx + i] - cen1;
            scale += (diffX*diffX + diffY*diffY);
        }
        scale = Math.sqrt(scale/(2.*(nFeatures - 1.)));
        // to use std dev instead: scale = Math.sqrt(scale/(n-1.));
        
        outputNorm[imageNumber*3 + 0] = cen0;
        outputNorm[imageNumber*3 + 1] = cen1;
        outputNorm[imageNumber*3 + 2] = scale;
        
        for (int i = 0; i < nFeatures; ++i) {
            xN[0][i] = (x[0][imageIdx + i] - cen0)/scale;
            xN[1][i] = (x[1][imageIdx + i] - cen1)/scale;
        }
        
        return xN;
    }

    private static void calculateLeftEpipole(DenseMatrix fundamentalMatrix,
        double[] outputE01) throws NotConvergedException {
        
        if (outputE01.length != 3) {
            throw new IllegalArgumentException("outputE01 length must be 3");
        }
        
        SVDProducts svdE = MatrixUtil.performSVD(fundamentalMatrix);
        
        /*
         The left epipole is e1 = last column of U / last item of that column
         It is  the left image position of the epipolar projection of the right camera center
         The right epipole e2 = last row of V / last item of that row
         It is the right image position of the epipolar projection of the left camera center
        */
        
        assert(svdE.u[0].length == 3);
        assert(svdE.u.length == 3);
        assert(svdE.vT[0].length == 3);
        assert(svdE.vT.length == 3);
        
        //double e1Div = svdE.vT[2][2];
        for (int i = 0; i < 3; i++) {
            outputE01[i] = svdE.vT[2][i];///e1Div;
        }
    }

    private static void extractColumn(double[][] x, int idx, double[] outputPoint) {
        outputPoint[0] = x[0][idx];
        outputPoint[1] = x[1][idx];
        outputPoint[2] = x[2][idx];
    }

    /**
     * 
     * @param motion 2*F X 3 matrix of rotation for images
     * where the first F rows are the first rows of each image's rotation matrix
     * and the second F rows are the second rows of each image's rotation matrix.
     * the third row of rotation is recreated when needed using the cross
     * product of the first 2 rows of each image's rotation matrix;
     * @return
     * @throws NotConvergedException 
     */
    private static double[][] solveForTransformationToOrthoNormal(double[][] motion) throws NotConvergedException {
         /*
        The rows of R represent the orientations of the horizontal and vertical camera
        reference axes throughout the stream, 
        while the columns of S are the coordinates of the P feature
        points with respect to their centroid.
        
    motion = [ i_hat_0[0] i_hat_0[2] i_hat_0[2] ] where i_hat_f are unit vectors
             [ i_hat_1[0] i_hat_1[2] i_hat_1[2] ]
             [   ..._m-1 ...                    ]
             [ j_hat_0[0] j_hat_0[2] j_hat_0[2] ]
             [ j_hat_1[0] j_hat_1[2] j_hat_1[2] ] where j_hat_f are unit vectors
             [  ..._m-1 ... ]
        
    shape = [s_C_1  ...  s_C_m] is [3XP] where P is the number of points per image frame
        NOTE that the summation over columns of sC = 0 (they are centered w.r.t. image points)
        
        R = motion*Q
        S = (Q^-1)*shape
        
          eqn (1)  (`i_f)^T * Q * Q^T * (`i_f) = 1   [dimensions 1X3 * 3X3 * 3X1 = 1]
          eqn (2)  (`j_f)^T * Q * Q^T * (`j_f) = 1
          eqn (3)  (`i_f)^T * Q * Q^T * (`j_f) = 0
        
        find Q as a 3 × 3 matrix, non-singular matrix
        
        //Morita and Kanade
        L = Q^T * Q
        solve the linear system of equations for L 
            and use Cholesky decomposition to get Q.
            Correct the decomposition to enforce L to be positive definite
            symmetric.
        
        See notes reference Morita and Kanade for solving Q.
         T. Morita and T. Kanade, A Sequential Factorization Method for Recovering Shape and Motion
         from Image Streams, Pattern Analysis and Machine Intelligence, IEEE Transactions on, vol. 19,
         no.8, pp.858-867, Aug 1997  (1994?)
        
        http://note.sonots.com/SciSoftware/Factorization.html#cse252b
        
        L = [ l1 l2 l3 ]
            [ l2 l4 l5 ]
            [ l3 l5 l6 ]     

        the knowns are `i_f and `j_f, so we are solving for the 6 unknowns in L.

        expand the terms:

        eqn(1):
        if_0  if_1  if_2 ] * [ l1 l2 l3 ] * [ if_0 ] = 1
                             [ l2 l4 l5 ]   [ if_1 ]
                             [ l3 l5 l6 ]   [ if_2 ]
        eqn(2):
        jf_0  jf_1  jf_2 ] * [ l1 l2 l3 ] * [ jf_0 ] = 1
                             [ l2 l4 l5 ]   [ jf_1 ]
                             [ l3 l5 l6 ]   [ jf_2 ]
        eqn(3):
        if_0  if_1  if_2 ] * [ l1 l2 l3 ] * [ jf_0 ] = 0
                             [ l2 l4 l5 ]   [ jf_1 ]
                             [ l3 l5 l6 ]   [ jf_2 ]
        
        eqn(1):
        l1*if_0*if_0 + l2*if_1*if_0 + l3*if_2*if_0 + l2*if_0*if_1 + l4*if_1*if_1 + l5*if_2*if_1 + l3*if_0*if_2 + l5*if_1*if_2 + l6*if_2*if_2 = 1
        eqn(2):
        l1*jf_0*jf_0 + l2*jf_1*jf_0 + l3*jf_2*jf_0 + l2*jf_0*jf_1 + l4*jf_1*jf_1 + l5*jf_2*jf_1 + l3*jf_0*jf_2 + l5*jf_1*jf_2 + l6*jf_2*jf_2 = 1
        eqn(3):
        l1*if_0*jf_0 + l2*if_1*jf_0 + l3*if_2*jf_0 + l2*if_0*jf_1 + l4*if_1*jf_1 + l5*if_2*jf_1 + l3*if_0*jf_2 + l5*if_1*jf_2 + l6*if_2*jf_2 = 1
        
        factor out the L terms, linearly
        l1           l2                    l3                    l4           l5                    l6         const
        (if0*if0)    (if1*if0 + if0*if1)   (if2*if0 + if0*if2)   (if1*if1)    (if2*if1 + if1*if2)   (if2*if2)   1
        (jf0*jf0)    (jf1*jf0 + jf0*jf1)   (jf2*jf0 + jf0*jf2)   (jf1*jf1)    (jf2*jf1 + jf1*jf2)   (jf2*jf2)   1
        (if0*jf0)    (if1*jf0 + if0*jf1)   (if2*jf0 + if0*jf2)   (if1*jf1)    (if2*jf1 + if1*jf2)   (if2*jf2)   0

        since the terms in the rows have a similar pattern, can write the equation more concisely using
        a function to generate them:
           g(a,b) = [a0*b0         ]
                    [a0*b1 + a1*b0 ]
                    [a0*b2 + a2*b0 ]
                    [a1*b1         ]
                    [a1*b2 + a2*b1 ]
                    [a2*b2         ]

        G size is 3*F X 6
        L is a vector of length 6
        c size is 3*F
        
        the G = [ g(i_0, i_0)^T       ]   L_vectorized = [l1]    c = [2*F rows of 1]
                [ ...each row thru F  ]                  [l2]        [F rows of 0  ]
                [ g(j_0, j_0)^T       ]                  [l3]
                [ ...each row thru F  ]                  [l4]
                [ g(i_0, j_0)^T       ]                  [l5]
                [ ...each row thru F  ]                  [l6]

        G*L_vectorized = c ==>  L_vectorized = G^-1 * c
        */
        
         int mImages = motion.length/2;
        
        int i, j;
        double[] c = new double[3*mImages];
        for (i = 0; i < 2*mImages; ++i) {
            c[i] = 1;
        }
        
        //g is 3F X 6
        double[][] g = new double[3*mImages][6];
        for (i = 0; i < 2*mImages; ++i) {
            g[i] = gT(motion[i], motion[i]);
            assert(g[i].length == 6);
        }
        for (i = 2*mImages, j=0; i < 3*mImages; ++i, j++) {
            g[i] = gT(motion[j], motion[mImages + j]);
        }
                
        // rank of G is 4
        //G * L = c
        // G^T * G * L = G^T * c
        // L = (G^T*G)^-1*G^T * c
        //     pseudoinverse of G is (G^T*G)^-1*G^T
        double[][] gInv = MatrixUtil.pseudoinverseRankDeficient(g);
        
        // 6X1
        double[] lVector = MatrixUtil.multiplyMatrixByColumnVector(gInv, c);
        assert(lVector.length == 6);
        
        // 3X3
        double[][] ell = new double[3][3];
        ell[0] = new double[]{lVector[0], lVector[1], lVector[2]};
        ell[1] = new double[]{lVector[1], lVector[3], lVector[4]};
        ell[2] = new double[]{lVector[2], lVector[4], lVector[5]};
        
        // Q can be determined :
        //   as the square root of ell,
        //   or with the Cholesky decomposition
        //   or with eigendecomposition
        
        // enforcing positive definiteness of L (which is ell here).
        double eps = 1e-16;//1.e-11; eps close to zero within machine precision to perturb the matrix to smallest eigenvalue of eps
        double[][] lPSD = MatrixUtil.nearestPositiveSemidefiniteToASymmetric(ell, eps);
        /*EVD evd2 = EVD.factorize(new DenseMatrix(lPSD));
        double[] eig = evd2.getRealEigenvalues();
        double[][] aMinusPSD = MatrixUtil.elementwiseSubtract(ell, lPSD);
        double dist1 = MatrixUtil.frobeniusNorm(aMinusPSD);
        */
        boolean ipd = MatrixUtil.isPositiveDefinite(lPSD);
        assert(ipd);
       
        //decompose Q = L * (sigma+) * L^T;  Q is size 3X3
        double[][] q = LinearEquations.choleskyDecompositionViaMTJ(lPSD);
        
        return q;
    }

    private static void assertDotProductMetrics(double[][] motion,
        int mImages) {
        double[] tmp1, tmp2, tmp3;
        int i, j;
        double tmp4, tmp5, tmp6, tmp7;
        for (i = 0; i < mImages; ++i) {
            tmp1 = Arrays.copyOf(motion[i], motion[i].length);
            tmp2 = Arrays.copyOf(motion[mImages + i], motion[mImages + i].length);
            tmp3 = MatrixUtil.crossProduct(tmp1, tmp2);
            tmp4 = 0; // i dot i = 1
            tmp5 = 0; // j dot j = 1
            tmp6 = 0; // i dot j = 0
            tmp7 = 0; // k dot k = 1
            for (j = 0; j < 3; ++j) {
                tmp4 += (tmp1[j]*tmp1[j]);
                tmp5 += (tmp2[j]*tmp2[j]);
                tmp6 += (tmp1[j]*tmp2[j]);
                tmp7 += (tmp3[j]*tmp3[j]);
            }
            System.out.printf("tmps:%.4e, %.4e, %.4e, %.4e\n", tmp4, tmp5, tmp7, tmp6);
            /* if had perfect data:
            assert(Math.abs(tmp4 - 1) < 1e-2);
            assert(Math.abs(tmp5 - 1) < 1e-2);
            assert(Math.abs(tmp7 - 1) < 1e-2);
            assert(Math.abs(tmp6) < 1e-2);
            */
        }
    }
    
    public static class ProjectionResults {
        /**
         * world coordinate system points in matrix of size 4 X nFeatures.
         * The points are stacked along columns sequentially.
         */
        public double[][] XW;
        
        /**
         * the projection matrices stacked along rows for each image.
         * so projection for image 0 will be in rows [0, 3);
         * projection for image 1 will be in rows [3, 6), etc.
         * This matrix's size is 3*nImages X 4
         */
        public double[][] projectionMatrices;
    }
    
    public static class OrthographicProjectionResults {
        /**
         * world coordinate system points
         */
        public double[][] XW;
        
        /**
         * the rotation matrices stacked along rows for each image.
         * so rotation for image 0 will be in rows [0, 3);
         * rotation for image 1 will be in rows [3, 6), etc.
         */
        public double[][] rotationMatrices;
    }
    
    public static class ParaperspectiveProjectionResults {
        /**
         * world coordinate system points
         */
        private double[][] XW;
        
        /**
         * the rotation matrices (as extrinsic parameters) stacked along rows for each image.
         * so rotation for image 0 will be in rows [0, 3);
         * rotation for image 1 will be in rows [3, 6), etc.
         */
        private double[][] rotationMatrices;
        
        /**
         * the translation vectors (as extrinsic parameters) 
         * stacked along rows for each image.
         * 
         */
        private double[][] translationVectors;

        /**
         * @return the XW
         */
        public double[][] getXW() {
            return XW;
        }

        /**
         * @param XW the XW to set
         */
        public void setXW(double[][] XW) {
            this.XW = XW;
        }

        /**
         * @return the rotationMatrices (as extrinsic parameters) 
         * stacked along rows for each image.
         * so rotation for image 0 will be in rows [0, 3);
         * rotation for image 1 will be in rows [3, 6), etc.
         */
        public double[][] getRotationStack() {
            return rotationMatrices;
        }

        /**
         * @param rotationMatrices the rotationMatrices to set.
         * the rotation matrices (as extrinsic parameters) stacked along rows for each image.
         * so rotation for image 0 will be in rows [0, 3);
         * rotation for image 1 will be in rows [3, 6), etc.
         */
        public void setRotationStack(double[][] rotationMatrices) {
            this.rotationMatrices = rotationMatrices;
        }

        /**
         * @return the translationVectors
         */
        public double[][] getTranslationVectorStack() {
            return translationVectors;
        }

        /**
         * @param translationVectors the translationVectors to set
         */
        public void setTranslationVectorStack(double[][] translationVectors) {
            this.translationVectors = translationVectors;
        }
        
        public double[][] getExtrinsicProjection(int imageNumber) {
            double[][] p = new double[3][4];
            for (int i = 0; i < 3; ++i) {
                p[i] = new double[4];
                System.arraycopy(rotationMatrices[imageNumber*3 + i], 0, p[i], 0, 3);
                p[i][4] = translationVectors[imageNumber][i];
            }
            return p;
        }
    }

    public static class ReconstructionResults {
        double[][] XW;
        double[][] k1Intr;
        double[][] k2Intr;
        double[][] k1ExtrRot;
        double[] k1ExtrTrans;
        double[][] k2ExtrRot;
        double[] k2ExtrTrans;
        double[][] essentialMatrix;
        SVDProducts svd;
        double[][] fundamentalMatrix;
        
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("XW=\n");
            if (XW != null) {
                sb.append(FormatArray.toString(XW, "%.4e"));
            }
            sb.append("k1 intrinsic=\n");
            if (k1Intr != null) {
                sb.append(FormatArray.toString(k1Intr, "%.4e"));
            }
            sb.append("k1 extrinsic rotation=\n");
            if (k1ExtrRot != null) {
                sb.append(FormatArray.toString(k1ExtrRot, "%.4e"));
            }
            sb.append("k1 extrinsic translation=\n");
            if (k1ExtrTrans != null) {
                sb.append(FormatArray.toString(k1ExtrTrans, "%.4e"));
                sb.append("\n");
            }
            sb.append("k2 intrinsic=\n");
            if (k2Intr != null) {
                sb.append(FormatArray.toString(k2Intr, "%.4e"));
            }
            sb.append("k2 extrinsic rotation=\n");
            if (k2ExtrRot != null) {
                sb.append(FormatArray.toString(k2ExtrRot, "%.4e"));
            }
            sb.append("k2 extrinsic translation=\n");
            if (k2ExtrTrans != null) {
                sb.append(FormatArray.toString(k2ExtrTrans, "%.4e"));
                sb.append("\n");
            }
            return sb.toString();
        }
    }
    
      /**
     * among the 4 rotation and translation combinations from R1, R1, T1, and T2, 
     * select the one with the largest number of projected Z coordinates which are
     * positive, that is, in front of both cameras.
     * NOTE that inaccuracies in this chirality are larger for points further 
     * away from the cameras and closer to the plane at infinity.
     * NOTE that the determinants of R1 and R2 should have already been checked to be +1.
     * @param x1 image 1 portion of the correspondence pairs.
     * @param x2 image 2 portion of the correspondence pairs.
     * @param k1 intrinsic camera matrix for camera 1
     * @param k2 intrinsic camera matrix for camera 2
     * @param R1 rotation matrix whose determinant is +1
     * @param R2 rotation matrix whose determinant is +1
     * @param t1 translation vector (the direction between camera centers)
     * @param t2 translation vector (the direction between camera centers)
     * @param rSelected output variable holding the R1 or R2, whichever was the 
     * first found as a valid solution.
     * @param tSelected output variable holding the t1 or t2, whichever was the 
     * first found as a valid solution.
     * @param outputX the real world coordinates of the projection of x1 and x2 using
     * triangulation. else null if no valid solution was found
     */
    private static void bestInCheiralityTest(double[][] x1, double[][] x2, 
        double[][] k1, double[][] k2,
        double[][] R1, double[][] R2, double[] t1, double[] t2, 
        double[][] rSelected, double[] tSelected, double[][] outputX) {
    
        int n = x1[0].length;
        
        if (outputX.length != 4) {
            throw new IllegalArgumentException("outputX.length must be 4");
        }
        if (outputX[0].length != n) {
            throw new IllegalArgumentException("outputX[0].length must be the same as x1[0].length");
        }
        
        // for this model, for the first image, the camera extrinsics are
        //    R = I and t = [0], which leaves all rotation and translation in
        //    the 2nd camera extrinsics w.r.t. the first.
        double[][] k1ExtrRot = MatrixUtil.createIdentityMatrix(3);
        double[] k1ExtrTrans = new double[3];
        
        // save the first that pass the tests for Z>=0.
        double[][] bestR = null;
        double[] bestT = null;
        double[][] bestXW = null;
        String bestLabel = null;
        int bestNPosZ = Integer.MIN_VALUE;
        
        double[][] XW;
        double[] XWPt;
        String label = null;
        
        XWPt = new double[4];
        XW = new double[4][n];
        for (int i = 0; i < 4; ++i) {
            XW[i] = new double[n];
        }
            
        double[][] rTst = null;
        double[] tTst = null;
        double[][] x1Pt = new double[3][1];
        double[][] x2Pt = new double[3][1];
        int i, j, ii;
        for (i = 0; i < 3; ++i) {
            x1Pt[i] = new double[1];
            x2Pt[i] = new double[1];
        }
        
        int nPosZ; 
        
        for (j = 0; j < 4; ++j) {
            switch(j) {
                case 0: {
                    label = "R1, T1";
                    rTst = R1;
                    tTst = t1;
                    break;
                }
                case 1: {
                    label = "R1, T2";
                    rTst = R1;
                    tTst = t2;
                    break;
                }
                case 2: {
                    label = "R2, T1";
                    rTst = R2;
                    tTst = t1;
                    break;
                }
                default: {                    
                    label = "R2, T2";
                    rTst = R2;
                    tTst = t2;
                    break;
                }
            }
            nPosZ = 0;
            for (i = 0; i < n; ++i) {
                for (ii = 0; ii < 3; ++ii) {
                    x1Pt[ii][0] = x1[ii][i];
                    x2Pt[ii][0] = x2[ii][i];
                }
                //
                XWPt = Triangulation.calculateWCSPoint(
                    k1, k1ExtrRot, k1ExtrTrans, 
                    k2, rTst, tTst, 
                    x1Pt, x2Pt);
                if (XWPt[2] >= 0) {
                    nPosZ++;
                }
                for (ii = 0; ii < 4; ++ii) {
                    XW[ii][i] = XWPt[ii];
                } 
            }
            if (nPosZ > bestNPosZ) {
                bestNPosZ = nPosZ;
                bestR = rTst;
                bestT = tTst;
                bestLabel = label;
                bestXW = MatrixUtil.copy(XW);
            }
        }
        
        if (bestR == null) {
            return;
        }
        
        // copy into output variables:
        for (i = 0; i < bestR.length; ++i) {
            System.arraycopy(bestR[i], 0, rSelected[i], 0, bestR[i].length);
        }
        System.arraycopy(bestT, 0, tSelected, 0, bestT.length);
        
        System.out.println("choosing solution: " + bestLabel);
        //double estimatedRotY = Math.atan(R[0][2]/R[0][0]) * (180./Math.PI);
        double estimatedRotZ = Math.atan(-bestR[1][0]/bestR[1][1]) * (180./Math.PI);
        System.out.printf("estimated rotation in degrees about z axis from R=%.2f\n", estimatedRotZ);
        //System.out.printf("X_WCS=\n%s\n", FormatArray.toString(bestXW, "%.3e"));
        System.out.flush();
        
        for (i = 0; i < XW.length; ++i) {
            System.arraycopy(bestXW[i], 0, outputX[i], 0, bestXW[i].length);
        }
        
    }

    static void populateWithDet1Rs(MatrixUtil.SVDProducts svdE, 
        double[][] r1Out, double[][] r2Out, double[][] uOut) {
        
        //Szeliski 2010, eqn (7.25)
        
        // R_Z+90 and R_Z_-90 from 
        // Ma, Soatto, Kosecká, and Sastry 2012, "An Invitation to 3-D Vision", pg 121
        
        //R_z_90^T  = [ [0, 1, 0], [0, -1, 0], [0, 0, 1] ]
        //R_z_-90^T = [ [0, -1, 0], [0, 1, 0], [0, 0, 1] ]
         
        //R_z_90
        double[][] r90T = new double[3][3];
        r90T[0] = new double[]{0, 1, 0};
        r90T[1] = new double[]{-1, 0, 0};
        r90T[2] = new double[]{0, 0, 1};
        double[][] r90NegT = MatrixUtil.transpose(r90T);
        
        double[][] u = svdE.u;
        double[][] uNeg = MatrixUtil.copy(u);
        MatrixUtil.multiply(uNeg, -1);
        
        double[][] rr1 = MatrixUtil.multiply(MatrixUtil.multiply(u, r90T), svdE.vT);
        double[][] rr2 = MatrixUtil.multiply(MatrixUtil.multiply(uNeg, r90T), svdE.vT);
        double[][] rr3 = MatrixUtil.multiply(MatrixUtil.multiply(u, r90NegT), svdE.vT);
        double[][] rr4 = MatrixUtil.multiply(MatrixUtil.multiply(uNeg, r90NegT), svdE.vT);
        
        double det1 = MatrixUtil.determinant(rr1);
        double det2 = MatrixUtil.determinant(rr2);
        double det3 = MatrixUtil.determinant(rr3);
        double det4 = MatrixUtil.determinant(rr4);
        
        System.out.printf("det(r1,r2,r3,r4)=%.3e,%.3e,%.3e,%.3e\n", det1, det2, det3, det4);
        System.out.printf("r1:\n%s\n", FormatArray.toString(rr1, "%.4e"));
        System.out.printf("r2:\n%s\n", FormatArray.toString(rr2, "%.4e"));
        System.out.printf("r3:\n%s\n", FormatArray.toString(rr3, "%.4e"));
        System.out.printf("r4:\n%s\n", FormatArray.toString(rr4, "%.4e"));
        
        boolean useUPos = true;
        
        int i;
        if (Math.abs(det1 - 1.) < eps) {
            System.out.printf("using +U\n");
            for (i = 0; i < 3; ++i) {
                System.arraycopy(rr1[i], 0, r1Out[i], 0, rr1[i].length);
            }
        } else if (Math.abs(det2 - 1.) < eps) {
            System.out.printf("using -U\n");
            useUPos = false;
            for (i = 0; i < 3; ++i) {
                System.arraycopy(rr2[i], 0, r1Out[i], 0, rr2[i].length);
            }
        } else {
            throw new IllegalStateException("neither rotation matrix is SO(3)");
        }
        if (Math.abs(det3 - 1.) < eps) {
            if (!useUPos) {
                throw new IllegalStateException("expecting to need +U");
            }
            System.out.printf("using +U\n");
            for (i = 0; i < 3; ++i) {
                System.arraycopy(rr3[i], 0, r2Out[i], 0, rr3[i].length);
            }
        } else if (Math.abs(det4 - 1.) < eps) {
            if (useUPos) {
                throw new IllegalStateException("expecting to need -U");
            }
            System.out.printf("using -U\n");
            for (i = 0; i < 3; ++i) {
                System.arraycopy(rr4[i], 0, r2Out[i], 0, rr4[i].length);
            }
        } else {
            throw new IllegalStateException("neither rotation matrix is SO(3)");
        }
        
        if (useUPos) {
            for (i = 0; i < 3; ++i) {
                System.arraycopy(u[i], 0, uOut[i], 0, u[i].length);
            }
        } else {
            for (i = 0; i < 3; ++i) {
                System.arraycopy(uNeg[i], 0, uOut[i], 0, uNeg[i].length);
            }
        }
    }
}
