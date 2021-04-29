package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.features.RANSACSolver;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.transform.Camera.CameraExtrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraIntrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraParameters;
import algorithms.imageProcessing.transform.Camera.CameraProjection;
import algorithms.matrix.MatrixUtil;
import algorithms.matrix.MatrixUtil.SVDProducts;
import algorithms.util.FormatArray;
import java.util.Arrays;
import java.util.List;
import no.uib.cipr.matrix.DenseCholesky;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SymmDenseEVD;

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
    */
            
     /**
     * given correspondence between two images calculate the camera
     * parameters, and the real world position.
     * 
     * <pre>
     * following CMU lectures of Kris Kitani at 
     * http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
     * 
     * add other references here
     * </pre>
     * @param camera1 image 1 camera matrices of intrinsic and extrinsic parameters.
     * the size is 3 x 4.
     * @param camera2 image 2 camera matrices of intrinsic and extrinsic parameters.
     * the size is 3 x 4.
     * @param x1 the image 1 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @param x2 the image 2 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @return 
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
        double[] XWPt = new double[4];
        
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
            //
            XWPt = Triangulation.calculateWCSPoint(
                camera1, camera2, x1Pt, x2Pt);
            for (ii = 0; ii < 4; ++ii) {
                XW[ii][i] = XWPt[ii];
            } 
        }
                
        ReconstructionResults rr = new ReconstructionResults();
        rr.XW = XW;

        return rr;
    }
    
    /**
     * recover the 3-D coordinates in WCS from pairs of corresponding
     * un-calibrated image points, that is, points in the image reference frame in pixels.
     * 
     * This is also called Projective Structure From Motion for the
     * Two-camera case.
     * 
     * NOTE that because the camera calibration, that is, intrinsic parameters,
     * are not known, only the projective reconstruction is possible,,
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
    public static ReconstructionResults2 calculateProjectiveReconstruction(
        double[][] x1, double[][] x2) throws NotConvergedException {
                        
        if (x1.length != 3 || x2.length != 3) {
            throw new IllegalArgumentException("x1.length must be 3 and so must x2.length");
        }
        int n0 = x1[0].length;
        if (x2[0].length != n0) {
            throw new IllegalArgumentException("x1 and x2 must be same dimensions");
        }
        
        /*
        http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
        
        (1) compute fundamental mat5rix FM from the correspondence x1, x2
        (2) compute the camera matrices P1, P2 from FM.
        (3) For each point correspondence, compute the point X in 3D space (triangularization)
        
        see also notes above frpm notes from Serge Belongie lectures from Computer Vision II, CSE 252B, USSD
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
        boolean reCalcIterations = false;
        
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
        RANSACSolver solver = new RANSACSolver();
        fitR = solver.calculateEpipolarProjection(
            leftM, rightM, errorType, useToleranceAsStatFactor, tolerance,
                reCalcIterations, false);
        
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
        
        MatrixUtil.SVDProducts svd = MatrixUtil.performSVD(_fm);
        double[] s = Arrays.copyOf(svd.s, svd.s.length);
        System.out.printf("SVD.s=\n%s\n", FormatArray.toString(s, "%.3e"));
        
        assert(svd.u[0].length == 3 && svd.u.length == 3);

        // Szeliski 2010 eqn 7.30.reform dropping the last singular value.
        s[2] = 0;
        _fm = MatrixUtil.multiplyByDiagonal(svd.u, s);
        _fm = MatrixUtil.multiply(_fm, svd.vT);
        svd = MatrixUtil.performSVD(_fm);
        double[] e2 = MatrixUtil.extractColumn(svd.u, 2);
        double[] e1 = svd.vT[2];
        
        double detV = MatrixUtil.determinant(svd.vT);
        
        System.out.printf("SVD.u=\n%s\n", FormatArray.toString(svd.u, "%.3e"));
        System.out.printf("SVD.s=\n%s\n", FormatArray.toString(s, "%.3e"));
        System.out.printf("SVD.vT=\n%s\n", FormatArray.toString(svd.vT, "%.3e"));
        System.out.printf("det(SVD.u)=%.2f\n", MatrixUtil.determinant(svd.u));
        System.out.printf("det(SVD.vT)=%.2f\n", detV);
                   
        // form the ambiguous homography:
        //   ~H = U * (R_90)^T * ~Sigma * V^T
        //    where ~sigma is the singular value matrix with the smallest value 
        //  replaced by a reasonable alternative (say, the middle value)
       
        // camera matrix P for left image = [I | 0 ]
        double[][] camera1 = MatrixUtil.zeros(3, 4);
        for (int i = 0; i < 3; ++i) {
            camera1[i][i] = 1;
        }
        
        double[][] r90 = new double[3][3];
        r90[0] = new double[]{0, -1, 0};
        r90[1] = new double[]{1, 0, 0};
        r90[2] = new double[]{0, 0, 1};        
        double[][] r90T = MatrixUtil.transpose(r90);
        s = Arrays.copyOf(svd.s, svd.s.length);
        s[2] = s[1]; // should this be 0?
        
        double[][] homog = MatrixUtil.multiply(svd.u, r90T);
        homog = MatrixUtil.multiplyByDiagonal(homog, s);
        homog = MatrixUtil.multiply(homog, svd.vT);
        
        e2 = MatrixUtil.extractColumn(svd.u, 2);
        e1 = svd.vT[2];
        
        System.out.printf("e1=%s\n", FormatArray.toString(e1, "%.3e"));
        System.out.printf("e2=%s\n", FormatArray.toString(e2, "%.3e"));
         
        // Szeliksi 2010 eqn 7.34:  P0 =[I|0] and P0 =[H|e],
        //    and then he finishes w/ triangulation
        // kitani lecture has P2 = [ [e1]_x * F | e2 ]
        // Hartley & Zisserman: has P2 = [ -[e2]_x^T * F | e2 ]
        //    note that slide 59 of http://16720.courses.cs.cmu.edu/lec/sfm.pdf
        //    also uses the Hartley & Zisserman: version of P2 and refers to proof in 
        //    Forsyth & Ponce Sec 8.3
        // Belongie uses P2 = [ [e2]_x^T * F | e2 ]
        
        // NOTE: not necessary to normalize the epipoles by the last value
        
        double[][] camera2S = MatrixUtil.zeros(3, 4);
        for (int i = 0; i < 3; ++i) {
            camera2S[i] = new double[4];
            System.arraycopy(homog[i], 0, camera2S[i], 0, 3);
            camera2S[i][3] = e2[i];
        }
        
        double[][] e1SkewSym = MatrixUtil.skewSymmetric(e1);
        double[][] e1F = MatrixUtil.multiply(e1SkewSym, _fm);
        double[][] camera2K = MatrixUtil.zeros(3, 4);
        for (int i = 0; i < 3; ++i) {
            camera2K[i] = new double[4];
            System.arraycopy(e1F[i], 0, camera2K[i], 0, 3);
            camera2K[i][3] = e2[i];
        }
        double[][] e2SkewSymT = MatrixUtil.transpose(MatrixUtil.skewSymmetric(e2));
        MatrixUtil.multiply(e2SkewSymT, -1);
        double[][] e2TF = MatrixUtil.multiply(e2SkewSymT, _fm);
        double[][] camera2H = MatrixUtil.zeros(3, 4);
        for (int i = 0; i < 3; ++i) {
            camera2H[i] = new double[4];
            System.arraycopy(e2TF[i], 0, camera2H[i], 0, 3);
            camera2H[i][3] = e2[i];
        }
        
        e2SkewSymT = MatrixUtil.transpose(MatrixUtil.skewSymmetric(e2));
        e2TF = MatrixUtil.multiply(e2SkewSymT, _fm);
        double[][] camera2B = MatrixUtil.zeros(3, 4);
        for (int i = 0; i < 3; ++i) {
            camera2B[i] = new double[4];
            System.arraycopy(e2TF[i], 0, camera2B[i], 0, 3);
            camera2B[i][3] = e2[i];
        }
        
        System.out.printf("Szeliski P2\n%s\n", FormatArray.toString(camera2S, "%.3e"));
        
        System.out.printf("Kitani P2\n%s\n", FormatArray.toString(camera2K, "%.3e"));
        
        System.out.printf("Hartley & Zisserman P2\n%s\n", FormatArray.toString(camera2H, "%.3e"));
        
        System.out.printf("Belongie P2\n%s\n", FormatArray.toString(camera2B, "%.3e"));
        
        CameraProjection P1 = new CameraProjection(camera1);
        CameraProjection P2S = new CameraProjection(camera2S);
        CameraProjection P2K = new CameraProjection(camera2K);
        CameraProjection P2H = new CameraProjection(camera2H);
        CameraProjection P2B = new CameraProjection(camera2B);
        
        double[][] XWS, XWK, XWH, XWB;
        double[] XWPt;
        
        XWPt = new double[4];
        XWS = new double[4][n];
        XWK = new double[4][n];
        XWH = new double[4][n];
        XWB = new double[4][n];
        for (int i = 0; i < 4; ++i) {
            XWS[i] = new double[n];
            XWK[i] = new double[n];
            XWH[i] = new double[n];
            XWB[i] = new double[n];
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
            XWPt = Triangulation.calculateWCSPoint(P1, P2S, x1Pt, x2Pt);
            for (j = 0; j < 4; ++j) {
                XWS[j][i] = XWPt[j];
            }
            XWPt = Triangulation.calculateWCSPoint(P1, P2K, x1Pt, x2Pt);
            for (j = 0; j < 4; ++j) {
                XWK[j][i] = XWPt[j];
            }
            XWPt = Triangulation.calculateWCSPoint(P1, P2H, x1Pt, x2Pt);
            for (j = 0; j < 4; ++j) {
                XWH[j][i] = XWPt[j];
            }
            XWPt = Triangulation.calculateWCSPoint(P1, P2B, x1Pt, x2Pt);
            for (j = 0; j < 4; ++j) {
                XWB[j][i] = XWPt[j];
            }
        }
        
        System.out.printf("Szeliski reconstruction\n%s\n", FormatArray.toString(XWS, "%.3e"));
        
        System.out.printf("Kitani reconstruction\n%s\n", FormatArray.toString(XWK, "%.3e"));
        
        System.out.printf("Hartley & Zisserman reconstruction\n%s\n", FormatArray.toString(XWH, "%.3e"));
        
        System.out.printf("Belongie reconstruction\n%s\n", FormatArray.toString(XWB, "%.3e"));
        
        ReconstructionResults2 rr = new ReconstructionResults2();
        rr.P1 = P1;
        rr.P2 = P2K;
        rr.XW = XWK;
        return rr;
    }

     /**
     * given correspondence between two images calculate the extrinsic camera
     * parameters, and the real world position.  This is called reconstruction
     * with calibrated cameras where the calibration refers to having the intrinsic
     * camera parameters.
     * 
     * <pre>
     * following CMU lectures of Kris Kitani at 
     * http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
     * 
     * add other references here
     * </pre>
     * @param k1 intrinsic camera matrix for image 1 in units of pixels.
     * @param k2 intrinsic camera matrix for image 2 in units of pixels.
     * @param x1 the image 1 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @param x2 the image 2 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @return 
     */
    public static ReconstructionResults calculateReconstruction(
        CameraIntrinsicParameters k1, CameraIntrinsicParameters k2,
        double[][] x1, double[][] x2) throws NotConvergedException {
        
        double[][] outputXW = MatrixUtil.zeros(4, x1[0].length);
        CameraExtrinsicParameters[] cps = 
            CameraPose.calculateUsingEssentialMatrix(
            k1.getIntrinsic(), k2.getIntrinsic(), x1, x2, outputXW);
                
        ReconstructionResults rr = new ReconstructionResults();
        rr.XW = outputXW;
        rr.k1ExtrRot = cps[0].getRotation();
        rr.k1ExtrTrans = cps[0].getTranslation();
        rr.k1Intr = k1.getIntrinsic();
        rr.k2ExtrRot = cps[1].getRotation();
        rr.k2ExtrTrans = cps[11].getTranslation();;
        rr.k2Intr = k2.getIntrinsic();

        return rr;
    }
    
    /**
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
     * .
     * Tomasi & Kanade 1991, "Shape and motion from image streams under 
     * orthography: a factorization method", International journal of computer vision 
     * 
     * Higham, 1988, “Computing a Nearest Symmetric Positive Semidefinite Matrix,” 
     *    Linear Algebra and Appl., 103:103-118, 1988
     * 
     * a great summary of the above:
     * http://note.sonots.com/SciSoftware/Factorization.html#cse252b
     * http://note.sonots.com/?plugin=attach&refer=SciSoftware%2FFactorization&openfile=Factorization.pdf
     * </pre>
     * NOTE: could overload this method to enable handling of occlusion 
     * following Section 5 of Tomasi & Kanade 1991, but might want to alter the
     * algorithm to use geometric median in place of centroid so that the
     * "centers" are not as affected by removing or adding a point.
     * NOTE: comments from Poelman & Kanade 1992:
     * Orthographic projection does not account for the apparent change in size 
     * of an object as it moves toward or away from the camera, nor the different 
     * angle from which an object is viewed as it moves parallel to the image plane.
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
                + "2 * mImages * nFeatures >= 8 * mImages + 3 * nFeatures – 12");
        }
        // for mImages=2, need 4 features
        
        // create matrix W composed of U and V 
        //     where U is rows of each image's x coordinates (size os mImages X nFeatures).
        //     where V is rows of each image's y coordinates (size os mImages X nFeatures).
        // create matrix t which holds the centroids of each row of W
        // create matrix WC = W - t
        
        double[][] w = MatrixUtil.zeros(2*mImages, nFeatures);
        double[] t = new double[2*mImages];
        int i, j;
        int m, n, xCol, vRow;
        for (m = 0; m < mImages; ++m) {
            vRow = mImages + m;
            for (n = 0; n < nFeatures; ++n) {
                xCol = m * nFeatures + n;
                w[m][n] = x[0][xCol];
                t[m] += x[0][xCol];
                w[vRow][n] = x[1][xCol];
                t[vRow] += x[1][xCol];
            }
            t[m] /= (double)nFeatures;
            t[vRow] /= (double)nFeatures;
        }
        
        //registered measurement matrix:
        double[][] wC = MatrixUtil.copy(w);
        for (i = 0; i < t.length; ++i) {
            for (n = 0; n < nFeatures; ++n) {
                wC[i][n] -= t[i];
            }
        }
        
        SVDProducts svd = MatrixUtil.performSVD(wC);
        double[][] u3 = MatrixUtil.copySubMatrix(svd.u, 0, svd.u.length, 0, 3);
        double[][] s3 = MatrixUtil.zeros(3, 3);
        s3[0][0] = svd.s[0];
        s3[1][1] = svd.s[1];
        s3[2][2] = svd.s[2];
        double[][] sqrts3 = MatrixUtil.zeros(3, 3);
        s3[0][0] = Math.sqrt(svd.s[0]);
        s3[1][1] = Math.sqrt(svd.s[1]);
        s3[2][2] = Math.sqrt(svd.s[2]);
        double[][] vT3 = MatrixUtil.copySubMatrix(svd.vT, 0, 3, 0, svd.vT[0].length);
        
        // if this is large, then the noise contribution can be ignored (cholesky not necessary)
        double sRatio = svd.s[2]/svd.s[3];
        System.out.printf("svd.s[2]/svd.s[3]=%.3e\n", sRatio);
                
        //Tamasi & Kanade
        // (2*mImages)X3
        double[][] rC = MatrixUtil.multiply(u3, sqrts3);
        // 3XnFeatures
        double[][] sC = MatrixUtil.multiply(sqrts3, vT3);
        
        // see Fig 3.1 of Tomasi & Kanade 1991 or Fig 2. of Belongie lecture notes
        
        // Belongie Section 16.4.4 (c)
        // Seo Step 3 - Metric Constraints
        /*
        The rows of R represent the orientations of the horizontal and vertical camera
        reference axes throughout the stream, 
        while the colums of S are the coordinates of the P feature
        points with respect to their centroid.
        
        rC = [ i_C_1^T ]
             [   ...   ]
             [ i_C_m^T ]
             [ j_C_1^T ]
             [   ...   ]
             [ j_C_m^T ]
        
        sC = [s_C_1  ...  s_C_m]
        */
        
        // constraints: enforce image axes to be orthonogal and length 1
        //    that is, the the rows of rC must have unit norm
        //    and the i_f's of rC must be perpendicular to the j_f’s where f is the
        //    image number(== frame number).
        /*
          eqn (1)  (`i_f)^T * Q * Q^T * (`i_f) = 1
          eqn (2)  (`j_f)^T * Q * Q^T * (`j_f) = 1
          eqn (3)  (`i_f)^T * Q * Q^T * (`j_f) = 0
        where Q is a 3 × 3 matrix
        
        L = Q*Q^T and solve the linear system of equations for L 
            and use Cholesky decomposition to get Q.
            Correct the decomposition to enforce L to be positive definite
            symmetric.
        
        Seo notes reference Morita and Kanade for solving Q.
         T. Morita and T. Kanade, A Sequential Factorization Method for Recovering Shape and Motion
         from Image Streams, Pattern Analysis and Machine Intelligence, IEEE Transactions on, vol. 19,
         no.8, pp.858-867, Aug 1997    
        
        http://note.sonots.com/SciSoftware/Factorization.html#cse252b
        
          eqn (1)  (`i_f)^T * Q * Q^T * (`i_f) = 1
          eqn (2)  (`j_f)^T * Q * Q^T * (`j_f) = 1
          eqn (3)  (`i_f)^T * Q * Q^T * (`j_f) = 0

        let L = Q*Q^T.  it's symmetric and square:
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
        
        (if0*l1 + if1*l2 + if2*l3)*if0 + (if0*l2 + if1*l4 + if2*l5)*if1 + (if0*l3 + if1*l5 + if2*l6)*if2   = 1
        (jf0*l1 + jf1*l2 + jf2*l3)*jf0 + (jf0*l2 + jf1*l4 + jf2*l5)*jf1 + (jf0*l3 + jf1*l5 + jf2*l6)*if2   = 1
        (if0*l1 + if1*l2 + if2*l3)*jf0 + (if0*l2 + if1*l4 + if2*l5)*jf1 + (if0*l3 + if1*l5 + if2*l6)*jf2   = 0

        rewriting:
        l1*if0*if0 + l2*if1*if0 + l3*if2*if0 + l2*if0*if1 + l4*if1*if1 + l5*if2*if1 + l3*if0*if2 + l5*if1*if2 + l6*if2*if2   = 1
        l1*jf0*jf0 + l2*jf1*jf0 + l3*jf2*jf0 + l2*jf0*jf1 + l4*jf1*jf1 + l5*jf2*jf1 + l3*jf0*jf2 + l5*jf1*jf2 + l6*jf2*jf2   = 1
        l1*if0*jf0 + l2*if1*jf0 + l3*if2*jf0 + l2*if0*jf1 + l4*if1*jf1 + l5*if2*jf1 + l3*if0*jf2 + l5*if1*jf2 + l6*if2*jf2   = 0

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

        the G = [ g(i_0, i_0)^T       ]   L_vectorized = [l1]    c = [2*F rows of 1]
                [ ...each row thru F  ]                  [l2]        [F rows of 0  ]
                [ g(j_0, j_0)^T       ]                  [l3]
                [ ...each row thru F  ]                  [l4]
                [ g(i_0, j_0)^T       ]                  [l5]
                [ ...each row thru F  ]                  [l6]

        G*L_vectorized = c ==>  L_vectorized = G^-1 * c
        */
        
        double[] c = new double[3*mImages];
        for (i = 0; i < 2*mImages; ++i) {
            c[i] = 1;
        }
        
        //g is 3F × 6
        double[][] g = new double[3*mImages][6];
        for (i = 0; i < 2*mImages; ++i) {
            g[i] = gT(rC[i], rC[i]);
            assert(g[i].length == 6);
        }
        j = 0;
        for (i = 2*mImages; i < 3*mImages; ++i) {
            g[i] = gT(rC[j], rC[mImages + j]);
            j++;
        }
        double[][] gInv = MatrixUtil.pseudoinverseRankDeficient(g);
        
        // 6X1
        double[] iVector = MatrixUtil.multiplyMatrixByColumnVector(gInv, c);
        assert(iVector.length == 6);
        
        // 3X3
        double[][] ell = new double[3][3];
        ell[0] = new double[]{iVector[0], iVector[1], iVector[2]};
        ell[1] = new double[]{iVector[1], iVector[3], iVector[4]};
        ell[2] = new double[]{iVector[2], iVector[4], iVector[5]};
        
        // Q can be determined :
        //   as the squaare root of ell,
        //   or with the Cholesky decomposition
        //   or with eigendecomposition
        
        /*
        // enforcing positive definiteness of L (which is ell here).
        ell = MatrixUtil.elementwiseAdd(ell, MatrixUtil.transpose(ell));
        MatrixUtil.multiply(ell, 0.5);
        
        Not using this as ell is written above as symmetric.
        The correction might be for use with the eigenvalue decomposition.
        If one uses the general EVD instead of the symmetric EVD method, 
        the left and right eigenvector matrices are a little different
        so forming one out of the average of both would be useful as input
        for forming Q below.  The correction's not necessary for the symmetric EVD method.
        */
        
        double eps = 1e-5;
        
        SymmDenseEVD evd = SymmDenseEVD.factorize(new DenseMatrix(ell));
        double[][] ellSigmaSqrt = MatrixUtil.zeros(evd.getEigenvalues().length, evd.getEigenvalues().length);
        for (i = 0; i < ellSigmaSqrt.length; ++i) {
            if (ellSigmaSqrt[i][i] < 0) {
                // replace with very small value
                ellSigmaSqrt[i][i] = eps;
            } else {
                ellSigmaSqrt[i][i] = Math.sqrt(ellSigmaSqrt[i][i]);
            }
        }
        double[][] lEig = MatrixUtil.convertToRowMajor(evd.getEigenvectors());
        // 3X3
        double[][] q = MatrixUtil.multiply(lEig, ellSigmaSqrt);
        
        // rC size is  (2*mImages)X3
        // sC size is 3XnFeatures
        double[][] r = MatrixUtil.multiply(rC, q);
        double[][] s = MatrixUtil.multiply(MatrixUtil.pseudoinverseRankDeficient(q), sC);
        
        // Align the first camera reference system with the world reference system
        //
        double[] i1 = Arrays.copyOf(r[0], r[0].length);
        double i1Norm = MatrixUtil.lPSum(i1, 2);
        MatrixUtil.multiply(i1, 1./i1Norm);
        double[] j1 = Arrays.copyOf(r[mImages], r[mImages].length);
        double j1Norm = MatrixUtil.lPSum(j1, 2);
        MatrixUtil.multiply(j1, 1./j1Norm);
        double[] k1 = MatrixUtil.crossProduct(i1, j1);
        double k1Norm = MatrixUtil.lPSum(k1, 2);
        MatrixUtil.multiply(k1, 1./k1Norm);
        
        double[][] r0 = new double[3][3];
        for (i = 0; i < 3; ++i) {
            r0[i] = new double[]{i1[i], j1[i], k1[i]};
        }
        
        // with orthographic, can only reocver rotation, not translation
        //(2*mImages)X3
        r = MatrixUtil.multiply(r, r0);
        // s now holds the world reference frame coordinates
        // 3XnFeatures
        s = MatrixUtil.multiply(MatrixUtil.transpose(r0), s);
        
        // r has the i and j direction and k=i cross j.
        // create a stack of rotation matrices, one per image.
        double[][] rotStack = MatrixUtil.zeros(3*mImages, 3);
        double[] ic, jc, kc;
        for (i = 0; i < mImages; ++i) {
            ic = r[i];
            jc = r[mImages + i];
            kc = MatrixUtil.crossProduct(ic, jc);
            for (j = 0; j < 3; ++j) {
                rotStack[i*3 + j][0] = ic[j]; 
                rotStack[i*3 + j][1] = jc[j];
                rotStack[i*3 + j][2] = kc[j];
            }
        }
        
        OrthographicProjectionResults results = new OrthographicProjectionResults();
        results.XW = s;
        results.rotationMatrices = rotStack;
                
        return results;
    }
    
    /*
    from Szeliski 2010 and Poelman & Kanade 1992:
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
     * recover the 3-D coordinates in WCS and the projection matrices 
     * from pairs of corresponding
     * un-calibrated image points, that is, points in the image reference frame in pixels.
     * assumes a para-perspective camera model.
     *
     * <pre>
     * references:
     * 
     * Poelman & Kanade 1992, "A Paraperspective Factorization Method for Shape 
     * and Motion Recovery" 
     * 
     * </pre>
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
    public static ProjectionResults calculateParaperspectiveReconstruction(
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
                + "2 * mImages * nFeatures >= 8 * mImages + 3 * nFeatures – 12");
        }
        // for mImages=2, need 4 features
        
        /*
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
        
        paused here
        */
        
        throw new UnsupportedOperationException("not yet finished");
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
    
    public static class OrthographicProjectionResults {
        /**
         * world coordinate system points
         */
        double[][] XW;
        
        /**
         * the rotation matrices stacked along rows for each image.
         * so rotation for image 0 will be in rows [0, 3);
         * rotation for image 1 will be in rows [3, 6), etc.
         */
        double[][] rotationMatrices;
    }

    public static class ReconstructionResults {
        double[][] XW;
        double[][] k1Intr;
        double[][] k2Intr;
        double[][] k1ExtrRot;
        double[] k1ExtrTrans;
        double[][] k2ExtrRot;
        double[] k2ExtrTrans;
        
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
    
    public static class ReconstructionResults2 {
        double[][] XW;
        CameraProjection P1;
        CameraProjection P2;
        
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("XW=\n");
            if (XW != null) {
                sb.append(FormatArray.toString(XW, "%.4e"));
            }
            sb.append("projection matrix for camera 1=\n");
            if (P1 != null) {
                sb.append(FormatArray.toString(P1.getP(), "%.4e"));
            }
            sb.append("projection matrix for camera 2=\n");
            if (P2 != null) {
                sb.append(FormatArray.toString(P2.getP(), "%.4e"));
            }
            return sb.toString();
        }
    }
}
