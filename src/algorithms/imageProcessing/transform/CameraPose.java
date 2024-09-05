package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Camera.CameraExtrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraIntrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraPoseParameters;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import no.uib.cipr.matrix.*;

import java.io.IOException;
import java.util.Arrays;

/**
 given a set of features in 1 image in image coordinates and world coordinate space,
 estimate the intrinsic and extrinsic camera parameters.
 This is also called "geometric camera calibration".

 TODO: write overloaded methods to use quaternion rotation.
 see
 T. Barfoot, et al., Pose estimation using linearized rotations and quaternion algebra,
 Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049
 
 * TODO: consider solving with M-estimators.
 * see http://research.microsoft.com/en- us/um/people/zhang/INRIA/Publis/Tutorial-Estim/node24.html
   or robust MM-estimator, or Least trimmed squares (LTS)
 * 
 * From Tumurbaatar, and Kim 2019, Sensors (Basel). 2019 Apr; 19(8): 1905.
 * "Comparative Study of Relative-Pose Estimations from a Monocular Image 
 * Sequence in Computer Vision and Photogrammetry"
 * We show that homography-based approaches are more accurate than essential-matrix 
 * or relative orientation–based approaches under noisy conditions.
 *
 * @author nichole
 */
public class CameraPose {
    
    public static double eps = 1e-7;
    
    /**
     * given an image and matched image features coordinates with real world coordinates, estimate the
     * camera matrix intrinsic and extrinsic parameters using a pin-hole camera model.
     *
     * This is also known as estimating the Motion.
     * This method uses DLT and should be followed by non-linear optimization
     * to improve the parameter estimates.
     <pre>
      references:
        https://www.ipb.uni-bonn.de/html/teaching/msr2-2020/sse2-13-DLT.pdf
        slides from a lecture titled "Photogrammetry & Robotics Lab,
         Camera Calibration: Direct Linear Transform" by Cyrill Stachniss
     see also http://www.ipb.uni-bonn.de/photogrammetry-i-ii/
     and
       Kitani lecture notes http://www.cs.cmu.edu/~16385/s17/Slides/11.3_Pose_Estimation.pdf
     </pre>
     <pre>

     Regarding translation, p3 is included in the results.  p3 is the last column in the projection
     matrix calculated internally.  (2) and (4) outlined below are what you should consider using
     to estimate the translation, depending upon your system's use of translate then rotate or vice versa.
         If the user is assuming translate then rotate: X_c = R * (X_wcs - t).
             (1) The projection matrix constructed would be [R | -R*t]
             where the last column is -R*t, R is rotation, t is translation,
             XW is object in real world coordinate frame, X_c is the object location seen in
             the camera reference frame.
             In this case, one would extract the translation
             using t = -1*(R^-1)*p3.  Note that when properly formed, R^-1 = R^T because rotation is orthogonal and unitary.
             (2) For the context of X_im = K * X_c, we have P = [K*R | -K*R*t]
             where K is the intrinsic parameter matrix for the camera.
             In this case, one would extract the translation
             using t = -1 * R^-1 * K^-1 * p3.
         If the user is assuming rotate then translate, X_c = R * X_wcs + t.
             (3) The projection matrix constructed would be [R | t].
             In this case, one would extract the translation
             using t = p3.
             (4) For the context of X_im = K * X_c, we have P = [K*R | K*t].
             In this case, one would extract the translation
             using t = K^-1 * p3.

     This method returns case (2) translation in the CameraExtrinsics field and assumes that x are given
     in image coordinates.

     see also CameraCalibration.solveForHomography(...)
     </pre>
     * @param x the image coordinates of the features in pixels in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 6 features are needed.
     * NOTE x and X should both be distortion-free or both should be distorted.
     * @param X the world coordinates of the features in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 6 features are needed.
     * NOTE x and X should both be distortion-free or both should be distorted.
     * @return 
     */
    public static CameraPoseParameters calculatePoseFromXXW(double[][] x, double[][] X)
        throws NotConvergedException {

        if (x.length != 3) {
            throw new IllegalArgumentException("x.length must be 3");
        }
        if (X.length != 3) {
            throw new IllegalArgumentException("X.length must be 3");
        }
        int n = x[0].length;
        if (n < 6) {
            throw new IllegalArgumentException("x must have at least 6 correspondences");
        }
        if (X[0].length != n) {
            throw new IllegalArgumentException("the number of columns in X must be the same as in x");
        }

        // from x = P * X.
        double[][] p = calculatePFromXXW(x, X);

        //TODO: refine p to minimize the reproduction error ||x - P*X||^2

        return calculatePoseFromP(p);
    }

    /**
     * given the camera projection matrix (i.e. K*[R|t]), estimate the
     * camera matrix intrinsic and extrinsic parameters using a pin-hole camera model.
     *
     * This is also known as estimating the Motion.
     * This method uses DLT and should be followed by non-linear optimization
     * to improve the parameter estimates.
     <pre>
     references:
     https://www.ipb.uni-bonn.de/html/teaching/msr2-2020/sse2-13-DLT.pdf
     slides from a lecture titled "Photogrammetry & Robotics Lab,
     Camera Calibration: Direct Linear Transform" by Cyrill Stachniss
     see also http://www.ipb.uni-bonn.de/photogrammetry-i-ii/
     and
     Kitani lecture notes http://www.cs.cmu.edu/~16385/s17/Slides/11.3_Pose_Estimation.pdf

     Zhang, chap 2 of camera calibration book:
     https://people.cs.rutgers.edu/~elgammal/classes/cs534/lectures/CameraCalibration-book-chapter.pdf
     </pre>
     <pre>

     Regarding translation, p3 is included in the results.  p3 is the last column in the projection
     matrix calculated internally.  (2) and (4) outlined below are what you should consider using
     to estimate the translation, depending upon your system's use of translate then rotate or vice versa.
     If the user is assuming translate then rotate: X_c = R * (X_wcs - t).
     (1) The projection matrix constructed would be [R | -R*t]
     where the last column is -R*t, R is rotation, t is translation,
     XW is object in real world coordinate frame, X_c is the object location seen in
     the camera reference frame.
     In this case, one would extract the translation
     using t = -1*(R^-1)*p3.  Note that when properly formed, R^-1 = R^T because rotation is orthogonal and unitary.
     (2) For the context of X_im = K * X_c, we have P = [K*R | -K*R*t]
     where K is the intrinsic parameter matrix for the camera.
     In this case, one would extract the translation
     using t = -1 * R^-1 * K^-1 * p3.
     If the user is assuming rotate then translate, X_c = R * X_wcs + t.
     (3) The projection matrix constructed would be [R | t].
     In this case, one would extract the translation
     using t = p3.
     (4) For the context of X_im = K * X_c, we have P = [K*R | K*t].
     In this case, one would extract the translation
     using t = K^-1 * p3.

     This method returns case (2) translation in the CameraExtrinsics field and assumes that x are given
     in image coordinates.

     </pre>
     * @param p the camera projection matrix, i.e. K*[R|t] where K is the intrinsic camera matrix, R is rotation,
     *          matrix, t is translation vector.
     * @return
     */
    public static CameraPoseParameters calculatePoseFromP(double[][] p)
            throws NotConvergedException {

        if (p.length != 3) {
            throw new IllegalArgumentException("p.length must be 3");
        }
        if (p[0].length != 4) {
            throw new IllegalArgumentException("p[0].length must be 4");
        }

        // P = K * [ R | t ] = [K*R | K*t]  = [K*R | -K*R*C] where t=-R*C
        // let M = K*R
        // P = [M | M*C]

        double[] p3 = MatrixUtil.extractColumn(p, 3);

        double[][] M = MatrixUtil.copySubMatrix(p, 0, 2, 0, 2);

        double detM = MatrixUtil.determinant(M);
        detM = Math.round(detM * 1E11)/1E11;
        /*if (detM == 0) {
            System.out.println("camera center is at infinity");
        } else {
            System.out.println("finite projective camera");
        }*/

        // method from zhang chap2 of camera calibration book
        /* K2 = M*M^T = K*K^T
              = [ ku=(alpha^2 + gamma^2 + u0^2)  kc=(u0*v0 + c+alpha)     u0 ]
                [ kc=(u0*v0 + c+alpha)           kv=(alphav^2 + v0^2)     v0 ]
                [ u0                              v0                      1  ]
           normalize K2 so that K2[2][2]=1

          u0 = k2[0][2]
          v0 = k2[1][2]
          beta = sqrt(kv - v0^2)
          gamma = (kc - u0*v0)/beta
          alpha = sqrt(ku - u0^2 - gamma^2)

          alpha, beta > 0

          then form K from those
          R = K^-1 * M
          t = K^-1 * (last column of P)
         */
        double[][] k2 = MatrixUtil.multiply(M, MatrixUtil.transpose(M));
        MatrixUtil.multiply(k2, 1./k2[2][2]);

        double u0 = k2[0][2];
        double v0 = k2[1][2];
        double beta = Math.sqrt(k2[1][1] - v0*v0);
        double gamma = (k2[0][1] - u0*v0)/beta;
        double alpha = Math.sqrt(k2[0][0] - u0*u0 - gamma*gamma);
        double[][] kEst1 = new double[][]{
                {alpha, gamma, u0}, {0, beta, v0}, {0, 0, 1}
        };
        double[][] kEst1Inv = Camera.createIntrinsicCameraMatrixInverse(kEst1);
        double[][] rEst1 = MatrixUtil.multiply(kEst1Inv, M);
        double[] tEst1 = MatrixUtil.multiplyMatrixByColumnVector(kEst1Inv, p3);

        // scale _t and _r by:
        double lambda1_1 = 1./MatrixUtil.lPSum(MatrixUtil.extractColumn(rEst1, 0), 2);
        double lambda1_2 = 1./MatrixUtil.lPSum(MatrixUtil.extractColumn(rEst1, 1), 2);
        MatrixUtil.multiply(tEst1, lambda1_1);
        MatrixUtil.multiply(rEst1, lambda1_1);

        // because the intrinsic camera matrix K has a positive diagonal, we can rewrite K*R as:
        //  (K * D) * (D^-1 * R)
        // to enforce the positive diagonal, though this should be checked above when constructing kEst1
        double[] d = new double[3];
        for (int i = 0; i < 3; ++i) {
            d[i] = (kEst1[i][i] > 0) ? 1 : -1;
        }
        kEst1 = MatrixUtil.multiplyByDiagonal(kEst1, d);
        rEst1 = MatrixUtil.multiplyDiagonalByMatrix(d, rEst1);

        /*
        // ==========================
        //(1) method from kitani lecture
        RQ rqDecomp = RQ.factorize(new DenseMatrix(M));
        double[][] kEst = MatrixUtil.convertToRowMajor(rqDecomp.getR());
        double[][] rEst = MatrixUtil.convertToRowMajor(rqDecomp.getQ());

        // because the intrinsic camera matrix K has a positive diagonal, we can rewrite K*R as:
        //  (K * D) * (D^-1 * R)
        // to enforce the positive diagonal:
        d = new double[3];
        for (int i = 0; i < 3; ++i) {
            d[i] = (kEst[i][i] > 0) ? 1 : -1;
        }
        kEst = MatrixUtil.multiplyByDiagonal(kEst, d);
        rEst = MatrixUtil.multiplyDiagonalByMatrix(d, rEst);
        MatrixUtil.multiply(kEst, 1./kEst[2][2]);

        //(2) estimate the camera position C
        // use assumption P*C=0 to solve for the null space of P, but we lose scale
        MatrixUtil.SVDProducts svdP = MatrixUtil.performSVD(p);
        double[] C = Arrays.copyOf(svdP.vT[svdP.vT.length - 1], 3);
        MatrixUtil.multiply(C, 1./C[2]);

        // t = -R*C
        double[] tEst = MatrixUtil.multiplyMatrixByColumnVector(rEst, C);
        MatrixUtil.multiply(tEst, -1);
        */// ===========================

        CameraExtrinsicParameters extrinsics = new CameraExtrinsicParameters();
        extrinsics.setRotation(rEst1);
        extrinsics.setTranslation(tEst1);

        CameraIntrinsicParameters intrinsics = new CameraIntrinsicParameters();
        intrinsics.setIntrinsic(kEst1);

        CameraPoseParameters camera = new CameraPoseParameters(intrinsics, extrinsics, p3);

        return camera;
    }

    private static double meanOfAbs(double[] a) {
        double sum = 0;
        for (int i = 0; i < a.length; ++i) {
            sum += Math.abs(a[i]);
        }
        return sum/a.length;
    }

    /**
     * given a single image which has matched image positions x with real world positions X,
     * estimate the camera matrix P using a camera model x = P * X and DLT.
     * @param x
     * @param X
     * @return
     * @throws NotConvergedException
     */
    public static double[][] calculatePFromXXW(double[][] x, double[][] X) throws NotConvergedException {
        if (x.length != 3) {
            throw new IllegalArgumentException("x.length must be 3");
        }
        if (X.length != 3 && X.length != 4) {
            throw new IllegalArgumentException("X.length must be 3 or 4");
        }
        int n = x[0].length;
        if (n < 6) {
            throw new IllegalArgumentException("x must have at least 6 correspondences");
        }
        if (X[0].length != n) {
            throw new IllegalArgumentException("the number of columns in X must be the same as in x");
        }

        boolean useNormConditioning = true;

        double[][] tx = null;
        double[][] tX = null;
        double[][] txInv = null;
        if (useNormConditioning) {
            if (X.length != 4) {
                double[][] tmp = new double[4][X[0].length];
                Arrays.fill(tmp[3], 1);
                for (int i = 0; i < X.length; ++i) {
                    System.arraycopy(X[i], 0, tmp[i], 0, X[i].length);
                }
                X = tmp;
            } else {
                X = MatrixUtil.copy(X);
            }
            x = MatrixUtil.copy(x);
            tx = EpipolarNormalizationHelper.unitStandardNormalize(x);
            txInv = EpipolarNormalizationHelper.inverseT(tx);
            tX = EpipolarNormalizationHelper.unitStandardNormalize(X);
        }

        int i, j;

        // 2*n X 12
        double xi, yi, Xi, Yi, Zi;
        double[][] ell = new double[2*n][12];
        for (i = 0; i < n; ++i) {
            xi = x[0][i];
            yi = x[1][i];
            Xi = X[0][i];
            Yi = X[1][i];
            Zi = X[2][i];
            ell[2*i]     = new double[]{Xi, Yi, Zi, 1, 0, 0, 0, 0, -xi*Xi, -xi*Yi, -xi*Zi, -xi};
            ell[2*i + 1] = new double[]{0, 0, 0, 0, Xi, Yi, Zi, 1, -yi*Xi, -yi*Yi, -yi*Zi, -yi};
            //ell[2*i]     = new double[]{Xi, Yi, Zi, 1, 0, 0, 0, 0, xi*Xi, xi*Yi, xi*Zi, xi};
            //ell[2*i + 1] = new double[]{0, 0, 0, 0, Xi, Yi, Zi, 1, yi*Xi, yi*Yi, yi*Zi, yi};
        }

        // SVD(ell).V is 12 X 12.  SVD(ell).U is 2n X 2n
        MatrixUtil.SVDProducts svd = MatrixUtil.performSVD(MatrixUtil.createATransposedTimesA(ell));

        // vT is 12X12.  last row in vT is the eigenvector for the smallest eigenvalue
        // it is also the epipole e1, defined as the right nullspace
        double[] xOrth = svd.vT[svd.vT.length - 1];

        // subject to ||x|| = 1
        xOrth = MatrixUtil.normalizeLP(xOrth, 2);

        // assert that ell * xOrth ~ 0
        //double[] chk = MatrixUtil.multiplyMatrixByColumnVector(ell, xOrth);
        //System.out.printf("check that A*x=0:%s\n", FormatArray.toString(chk, "%.3e"));

        // reshape into 3 X 4
        double[][] P2 = MatrixUtil.zeros(3,4);
        System.arraycopy(xOrth, 0, P2[0], 0, 4);
        System.arraycopy(xOrth, 4, P2[1], 0, 4);
        System.arraycopy(xOrth, 8, P2[2], 0, 4);

        // enforce det(R) > 0
        double detR = MatrixUtil.determinant(MatrixUtil.copySubMatrix(P2, 0,2, 0, 2));
        detR = (detR>0) ? 1 : -1;
        MatrixUtil.multiply(P2, detR);

        if (useNormConditioning) {
            P2 = MatrixUtil.multiply(txInv, P2);
            P2 = MatrixUtil.multiply(P2, tX);
        }

        return P2;
    }

    /**
     * calculate the camera extrinsic parameters R and T, given a list of coordinates of objects in the world reference frame
     * corresponding to a list of their coordinates in an image frame and given the camera intrinsic parameters.
     * The method is ported from github repositories holding the Bouguet Matlab Toolbox code.
     <pre>
     The Bouguet toolbox webpage is currently at http://robots.stanford.edu/cs223b04/JeanYvesCalib/
     and states that the source code is freely available.
     The github repositories with forked Bouguet Matlab code do not have license
     information.  Those references are
     https://github.com/fragofer/TOOLBOX_calib
     and
     https://github.com/hunt0r/Bouguet_cam_cal_toolbox
     and the methods adapted from are
     compute_extrinsic_init.m, normalize_pixel.m, compute_homography.m,
     compute_extrinsic_refine.m,  project_points2.m, rigid_motion.m
     </pre>
     * @param intrinsics
     * @param x objects in image coordinate reference frame.  size [3Xn].  if given [2Xn], will stack a row of 1's onto it
     *          internally.
     * @param X objects in world coordinate reference frame.  size [3Xn]
     * @param useBouguetForRodrigues if true, uses only the Bouguet algoirthms for Rodrigues rotation matrices and vectors
     * @return
     * @throws NotConvergedException
     * @throws IOException
     */
    public static CameraExtrinsicParameters calculatePoseUsingBouguet(
            Camera.CameraIntrinsicParameters intrinsics, double[][] x,
            double[][] X, boolean refine, boolean useBouguetForRodrigues) throws NotConvergedException, IOException {

        if (x.length != 2 && x.length != 3) {
            throw new IllegalArgumentException("x.length must be 3 or 2");
        }
        if (X.length != 3) {
            throw new IllegalArgumentException("X.length must be 3");
        }
        int n = x[0].length;
        if (n < 4) {
            throw new IllegalArgumentException("x must have at least 4 correspondences");
        }
        if (X[0].length != n) {
            throw new IllegalArgumentException("the number of columns in X must be the same as in x");
        }

        int i, j;

        if (x.length == 3) {
            x = MatrixUtil.copy(x);
        } else {
            double[][] x2 = MatrixUtil.zeros(3, n);
            System.arraycopy(x[0], 0, x2[0], 0, n);
            System.arraycopy(x[1], 0, x2[1], 0, n);
            Arrays.fill(x2[2], 1);
            x = x2;
        }
        X = MatrixUtil.copy(X);

        //xn = normalize_pixel(x_kk,fc,cc,kc,alpha_c);
        double[][] xc = Camera.pixelToCameraCoordinates(x, intrinsics);

        // x = x_kk
        // X = X_kk
        //X_mean = mean(X_kk')';  [3 X 1]
        double[] XMean = new double[X.length];
        for (i = 0; i < X.length; ++i) {
            XMean[i] = MiscMath0.getAvgAndStDev(X[i])[0];
        }

        //Y = X_kk - (X_mean*ones(1,Np));  [3Xn]
        double[][] Y = new double[X.length][];
        for (i = 0; i < X.length; ++i) {
            Y[i] = MatrixUtil.subtract(X[i], XMean[i]);
        }

        //YY = Y*Y';  [3 X n][n X 3] = [3X3]
        double[][] YY = MatrixUtil.multiply(Y, MatrixUtil.transpose(Y));

        //[U,S,V] = svd(YY); [3X3] for U and V
        SVD svd = SVD.factorize(new DenseMatrix(YY));

        //r = S(3,3)/S(2,2);
        double r = svd.getS()[2]/svd.getS()[1];

        CameraExtrinsicParameters soln;

        if ((r < 1e-3)|| (n < 5)) { // test of planarity
            // planar structure
            //Transform the plane to bring it in the Z=0 plane:
            soln = bouguetPoseInitPlanar(intrinsics, xc, X, svd.getVt(), XMean, useBouguetForRodrigues);
        } else {
            //%fprintf(1,'Non planar structure detected: r=%f\n',r);
            soln = bouguetPoseInitNonPlanar(xc, X, useBouguetForRodrigues);
        }

        if (refine) {
            // could return more intermediate arrays such as JJ
            //this needs image coordinates because internally it is projecting X to camera then image and comparing that to x
            soln = bouguetPoseRefine(soln, intrinsics, x, X, useBouguetForRodrigues);
        }

        //computation of the homography (not useful in the end)
        if (true) {
            //H = [Rckk(:,1:2)Tckk];  // [3X2][3X1]
            double[][] H = new double[3][];
            for (i = 0; i < 3; ++i) {
                H[i] = new double[]{soln.getRotation()[i][0], soln.getRotation()[i][1], soln.getTranslation()[i]};
            }
            //% Computes the reprojection error in pixels:
            //x = project_points2(X_kk, omckk, Tckk, fc, cc, kc, alpha_c);
            ProjectedPoints pp = bouguetProjectPoints2(X, soln.getRodriguesVector(), soln.getTranslation(),
                    intrinsics, useBouguetForRodrigues);

            //ex = x_kk - x;
            double[][] ex = MatrixUtil.pointwiseSubtract(x, pp.getXEstAs3XN());

            double[][] err = MatrixUtil.copySubMatrix(ex, 0, 1, 0, ex[0].length - 1);
            double[] xMeanStdv = MiscMath0.getAvgAndStDev(err[0]);
            double[] yMeanStdv = MiscMath0.getAvgAndStDev(err[1]);
            System.out.printf("x err=%s\n", FormatArray.toString(xMeanStdv, "%.4e"));
            System.out.printf("y err=%s\n", FormatArray.toString(yMeanStdv, "%.4e"));

            //% Converts the homography in pixel units:
            //KK = [fc(1) alpha_c * fc(1) cc(1);
            //0 fc(2) cc(2);
            //0 0 1];
            //H = KK * H;
            H = MatrixUtil.multiply(intrinsics.getIntrinsic(), H);
            System.out.printf("H from [R_3X2 | t]=\n%s\n", FormatArray.toString(H, "%.3e"));
        }

        //can return [omckk,Tckk,Rckk,H,x,ex,JJ] if return JJ from refine
        return soln;
    }

    public static class ProjectedPoints {

        /**
         * [2 X n] projected points xEst = R*X+T, where R = rodrigues(om), X is world coordinates of object, and T is translation
         */
        public double[][] xEst;

        /**
         * [2*n X 3] derivatives of XP w.r.t. rotation vector om
         */
        public double[][] dxdom;

        /**
         * [2*n X 3] derivatives of XP w.r.t. translation vector
         */
        public double[][] dxdT;

        public double[][] getXEstAs3XN() {
            double[][] x2 = new double[3][];
            x2[0] = Arrays.copyOf(xEst[0], xEst[0].length);
            x2[1] = Arrays.copyOf(xEst[1], xEst[1].length);
            x2[2] = new double[xEst[0].length];
            Arrays.fill(x2[2], 1);
            return x2;
        }

        /**
         * [2*n X 2] derivatives of XP w.r.t. camera focal length
         */
        public double[][] dxdF;

        /**
         * [2*n X 2] derivatives of XP w.r.t. camera principal point.
         * Not all methods produce output for this.
         */
        public double[][] dxdC = null;

        /**
         * [2*n X 4] derivatives of XP w.r.t. camera distortion coefficients
         */
        public double[][] dxdK;

        /**
         * [2*n X 1] derivatives of XP w.r.t. camera skew coefficient between x and y pixel
         * (alpha = 0 <=> square pixels results in derivative of 0 also).
         * Not all methods produce output for this.
         */
        public double[] dxdAlpha = null;

        /**
         * [2*n X 3] derivatives of XP w.r.t. the real world point.
         * Not all methods produce output for this.
         */
        public double[][] dxdX = null;
    }

    /**
     * Projects a 3D structure onto the image plane.
     * Bouguet toolbox code project_points2.m
     *
     * @param X 3D structure in the world coordinate frame (3xN matrix for N points)
     * @param om rotation vector (3x1 vector) between world coordinate frame and camera reference frame.
     * @param t translation vector (3x1 vector) between world coordinate frame and camera reference frame.
     * @param intrinsics camera intrinsic parameters
*    @param useBouguetForRodrigues if true, uses only the Bouguet algorithms for Rodrigues rotation matrices and vectors
     * @return [xp, dxpdom, dxpdT] where xp are the Projected pixel coordinates (2xN matrix for N points)
     * dxpdom are the Derivatives of xp with respect to om ((2N)x3 matrix), and
     * dxpdT are the derivatives of xp with respect to T ((2N)x3 matrix).
     */
    @SuppressWarnings({"fallthrough"})
    public static ProjectedPoints bouguetProjectPoints2(double[][] X, double[] om, double[] t,
         CameraIntrinsicParameters intrinsics, boolean useBouguetForRodrigues) {

        //[m,n] = size(X);
        int m = X.length;
        int n = X[0].length;

        if (m != 3) {
            throw new IllegalArgumentException("X.length should be 3");
        }
        if (om.length != 3) {
            throw new IllegalArgumentException("om.length should be 3");
        }
        if (t.length != 3) {
            throw new IllegalArgumentException("t.length should be 3");
        }

        /*
        %Definitions:
            %Let P be a point in 3D of coordinates X in the world reference frame (stored in the matrix X)
            %The coordinate vector of P in the camera reference frame is:
                Xc = R*X + T
                    %where R is the rotation matrix corresponding to the rotation vector om:
                    R = rodrigues(om);
            %call x, y and z the 3 coordinates of Xc: x = Xc(1); y = Xc(2); z = Xc(3);
            %The pinehole projection coordinates of
                P is [a;b] where a=x/z and b=y/z.
            %call r^2 = a^2 + b^2.
            %The distorted point coordinates are:
                xd = [xx;yy] where:
            %
            %xx = a * (1 + kc(1)*r^2 + kc(2)*r^4 + kc(5)*r^6)      +      2*kc(3)*a*b + kc(4)*(r^2 + 2*a^2);
            %yy = b * (1 + kc(1)*r^2 + kc(2)*r^4 + kc(5)*r^6)      +      kc(3)*(r^2 + 2*b^2) + 2*kc(4)*a*b;
            %
            %The left terms correspond to radial distortion (6th degree), the right terms correspond to tangential distortion
            %
            %Finally, conversion into pixel coordinates: The final pixel coordinates vector xp=[xxp;yyp] where:
            %
            %xxp = f(1)*(xx + alpha*yy) + c(1)
            %yyp = f(2)*yy + c(2)
            %
            %
            %NOTE: About 90 percent of the code takes care of computing the Jacobian matrices
            %
            %
            %Important function called within that program:
            %
            %rodrigues.m: Computes the rotation matrix corresponding to a rotation vector
            %
            %rigid_motion.m: Computes the rigid motion transformation of a given structure
         */

        //[Y,dYdom,dYdT] = rigid_motion(X,om,T);
        ProjectedPoints pRM = bouguetRigidMotion(X, om, t, useBouguetForRodrigues); // in camera reference frame
        double[][] Y = pRM.xEst; // [3 X n]
        double[][] dYdom = pRM.dxdom; // [3*n X 3]
        double[][] dYdT = pRM.dxdT;   // [3*n X 3]

        //inv_Z = 1./Y(3,:);  [1Xn]
        double[] invZ = new double[n];
        int i, j;
        for (i = 0; i < n; ++i) {
            invZ[i] = 1./Y[2][i];
        }

        //x = (Y(1:2,:) .* (ones(2,1) * inv_Z)) ;
        double[][] x = new double[2][];
        x[0] = MatrixUtil.pointwiseMultiplication(Y[0], invZ);
        x[1] = MatrixUtil.pointwiseMultiplication(Y[1], invZ);

        //     ([1Xn] dot [1Xn])^T. [nX1] * [1X3] = [nX3]
        //bb = (-x(1,:) .* inv_Z)'*ones(1,3);
        //cc = (-x(2,:) .* inv_Z)'*ones(1,3);
        double[] tmp1 = new double[n];
        double[] tmp2 = new double[n];
        for (i = 0; i < n; ++i) {
            tmp1[i] = -x[0][i] * invZ[i];
            tmp2[i] = -x[1][i] * invZ[i];
        }
        // [n X 3]
        double[][] bb = MatrixUtil.outerProduct(tmp1, new double[]{1, 1, 1});
        double[][] cc = MatrixUtil.outerProduct(tmp2, new double[]{1, 1, 1});

        //dxdom = zeros(2*n,3);
        //                       [nX1][1X3]       . [3*n->n  X 3]       [nX3] . [[3*n->n  X 3]
        //dxdom(1:2:end,:) = ((inv_Z')*ones(1,3)) .* dYdom(1:3:end,:) + bb .* dYdom(3:3:end,:);
        //dxdom(2:2:end,:) = ((inv_Z')*ones(1,3)) .* dYdom(2:3:end,:) + cc .* dYdom(3:3:end,:);

        // [nX3]
        double[][] d0 = MatrixUtil.outerProduct(invZ, new double[]{1, 1, 1});  //((inv_Z')*ones(1,3))
        double[][] d1 = new double[n][]; // dYdom(1:3:end,:)
        double[][] d2 = new double[n][]; // dYdom(2:3:end,:)
        double[][] d3 = new double[n][]; // dYdom(3:3:end,:)
        for (i = 0; i < n; ++i) {
            d1[i] = Arrays.copyOf(dYdom[i*3], 3);
            d2[i] = Arrays.copyOf(dYdom[i*3 + 1], 3);
            d3[i] = Arrays.copyOf(dYdom[i*3 + 2], 3);
        }
        // [nX3]
        double[][] dd1a = MatrixUtil.pointwiseMultiplication(d0, d1);//((inv_Z')*ones(1,3)) .* dYdom(1:3:end,:)
        double[][] dd1b = MatrixUtil.pointwiseMultiplication(bb, d3);//bb .* dYdom(3:3:end,:);
        double[][] dd1 = MatrixUtil.pointwiseAdd(dd1a, dd1b);

        double[][] dd2a = MatrixUtil.pointwiseMultiplication(d0, d2);//((inv_Z')*ones(1,3)) .* dYdom(2:3:end,:)
        double[][] dd2b = MatrixUtil.pointwiseMultiplication(cc, d3);//cc .* dYdom(3:3:end,:);
        double[][] dd2 = MatrixUtil.pointwiseAdd(dd2a, dd2b);

        // [2*n X 3]
        double[][] dxdom = new double[2*n][];
        for (i = 0; i < n; ++i) {
            dxdom[2*i] = Arrays.copyOf(dd1[i], dd1[i].length);
            dxdom[2*i + 1] = Arrays.copyOf(dd2[i], dd2[i].length);
        }

        // dYdT is [3*n X 3]
        //dxdT = zeros(2*n,3);
        //dxdT(1:2:end,:) = ((inv_Z')*ones(1,3)) .* dYdT(1:3:end,:) + bb .* dYdT(3:3:end,:);
        //dxdT(2:2:end,:) = ((inv_Z')*ones(1,3)) .* dYdT(2:3:end,:) + cc .* dYdT(3:3:end,:);
        d1 = new double[n][]; // dYdT(1:3:end,:)
        d2 = new double[n][]; // dYdT(2:3:end,:)
        d3 = new double[n][]; // dYdT(3:3:end,:)
        for (i = 0; i < n; ++i) {
            d1[i] = Arrays.copyOf(dYdT[i*3], 3);
            d2[i] = Arrays.copyOf(dYdT[i*3 + 1], 3);
            d3[i] = Arrays.copyOf(dYdT[i*3 + 2], 3);
        }
        // [nX3]
        dd1a = MatrixUtil.pointwiseMultiplication(d0, d1);//((inv_Z')*ones(1,3)) .* dYdT(1:3:end,:)
        dd1b = MatrixUtil.pointwiseMultiplication(bb, d3);//bb .* dYdT(3:3:end,:)
        dd1 = MatrixUtil.pointwiseAdd(dd1a, dd1b);

        dd2a = MatrixUtil.pointwiseMultiplication(d0, d2);//((inv_Z')*ones(1,3)) .* dYdT(2:3:end,:)
        dd2b = MatrixUtil.pointwiseMultiplication(cc, d3);//cc .* dYdT(3:3:end,:)
        dd2 = MatrixUtil.pointwiseAdd(dd2a, dd2b);

        //[2*n X 3]
        double[][] dxdT = new double[2*n][];
        for (i = 0; i < n; ++i) {
            dxdT[2*i] = Arrays.copyOf(dd1[i], dd1[i].length);
            dxdT[2*i + 1] = Arrays.copyOf(dd2[i], dd2[i].length);
        }

        double[] k = intrinsics.getRadialDistortionCoeffs();
        if (k == null) {
            k = new double[]{0, 0};
        }

        //% Add distortion:
        //r2 = x(1,:).^2 + x(2,:).^2;
        double[] r2 = new double[n];
        for (i = 0;i < n; ++i) {
            r2[i] = x[0][i]*x[0][i] + x[1][i]*x[1][i];
        }

        //dxdom is [2*n X 3]
        //        [nX3]
        //dr2dom = 2*((x(1,:)')*ones(1,3)) .* dxdom(1:2:end,:) + 2*((x(2,:)')*ones(1,3)) .* dxdom(2:2:end,:);
        double[][] d01 = MatrixUtil.outerProduct(x[0], new double[]{1, 1, 1});//2*((x(1,:)')*ones(1,3))
        double[][] d02 = MatrixUtil.outerProduct(x[1], new double[]{1, 1, 1});//2*((x(2,:)')*ones(1,3))
        MatrixUtil.multiply(d01, 2);
        MatrixUtil.multiply(d02, 2);
        d1 = new double[n][]; // dxdom(1:2:end,:)   [nX3]
        d2 = new double[n][]; // dxdom(2:2:end,:)   [nX3]
        for (i = 0; i < n; ++i) {
            d1[i] = Arrays.copyOf(dxdom[i*2], 3);
            d2[i] = Arrays.copyOf(dxdom[i*2 + 1], 3);
        }
        // [nX3]
        dd1a = MatrixUtil.pointwiseMultiplication(d01, d1);//2*((x(1,:)')*ones(1,3)) .* dxdom(1:2:end,:)
        dd1b = MatrixUtil.pointwiseMultiplication(d02, d2);//2*((x(2,:)')*ones(1,3)) .* dxdom(2:2:end,:)
        double[][] dr2dom = MatrixUtil.pointwiseAdd(dd1a, dd1b);

        // dxdT is [2*n X 3]
        //dr2dT = 2*((x(1,:)')*ones(1,3)) .* dxdT(1:2:end,:) + 2*((x(2,:)')*ones(1,3)) .* dxdT(2:2:end,:);
        for (i = 0; i < n; ++i) {
            d1[i] = Arrays.copyOf(dxdT[i*2], 3);
            d2[i] = Arrays.copyOf(dxdT[i*2 + 1], 3);
        }
        dd1a = MatrixUtil.pointwiseMultiplication(d01, d1);//2*((x(1,:)')*ones(1,3)) .* dxdom(1:2:end,:)
        dd1b = MatrixUtil.pointwiseMultiplication(d02, d2);//2*((x(2,:)')*ones(1,3)) .* dxdom(2:2:end,:)
        double[][] dr2dT = MatrixUtil.pointwiseAdd(dd1a, dd1b);

        //r4 = r2.^2;
        double[] r4 = new double[r2.length];
        for (i = 0;i < r2.length; ++i) {
            r4[i] = r2[i]*r2[i];
        }

        // [nX1][1X3]  . [nX3]
        //dr4dom = 2*((r2')*ones(1,3)) .* dr2dom;
        //dr4dT = 2*((r2')*ones(1,3)) .* dr2dT;
        d0 = MatrixUtil.outerProduct(r2, new double[]{1, 1, 1});
        MatrixUtil.multiply(d0, 2);
        double[][] dr4dom = MatrixUtil.pointwiseAdd(d0, dr2dom);
        double[][] dr4dT = MatrixUtil.pointwiseAdd(d0, dr2dT);

        //r6 = r2.^3;
        double[] r6 = new double[r2.length];
        for (i = 0;i < r2.length; ++i) {
            r4[i] = r2[i]*r4[i];
        }

        //dr6dom = 3*((r2'.^2)*ones(1,3)) .* dr2dom;
        //dr6dT = 3*((r2'.^2)*ones(1,3)) .* dr2dT;
        d0 = MatrixUtil.outerProduct(r4, new double[]{1, 1, 1});
        MatrixUtil.multiply(d0, 3);
        double[][] dr6dom = MatrixUtil.pointwiseAdd(d0, dr2dom);
        double[][] dr6dT = MatrixUtil.pointwiseAdd(d0, dr2dT);

         //% Radial distortion:
        //cdist = 1   + k(1) * r2   + k(2) * r4   + k(5) * r6;
        //dcdistdom = k(1) * dr2dom + k(2) * dr4dom + k(5) * dr6dom; // [nX3]
        //dcdistdT  = k(1) * dr2dT  + k(2) * dr4dT  + k(5) * dr6dT;  // [nX3]
        double[] cdist = new double[n];
        Arrays.fill(cdist, 1);
        double[][] dcdistdom = null;
        double[][] dcdistdT = null;

        if (k != null) {
            dcdistdom = MatrixUtil.zeros(n, 3);
            dcdistdT = MatrixUtil.zeros(n, 3);
            double[] tmp;
            double[][] tmp3;
            switch(k.length) {
                case 5 : {
                    tmp = Arrays.copyOf(r6, r6.length);
                    MatrixUtil.multiply(tmp, k[4]);
                    cdist = MatrixUtil.add(cdist, tmp);

                    tmp3 = MatrixUtil.copy(dr6dom);
                    MatrixUtil.multiply(tmp3, k[4]);
                    dcdistdom = MatrixUtil.pointwiseAdd(dcdistdom, tmp3);

                    tmp3 = MatrixUtil.copy(dr6dT);
                    MatrixUtil.multiply(tmp3, k[4]);
                    dcdistdT = MatrixUtil.pointwiseAdd(dcdistdT, tmp3);
                    // fall through
                }
                case 2 : {
                    tmp = Arrays.copyOf(r4, r4.length);
                    MatrixUtil.multiply(tmp, k[1]);
                    cdist = MatrixUtil.add(cdist, tmp);

                    tmp3 = MatrixUtil.copy(dr4dom);
                    MatrixUtil.multiply(tmp3, k[1]);
                    dcdistdom = MatrixUtil.pointwiseAdd(dcdistdom, tmp3);

                    tmp3 = MatrixUtil.copy(dr4dT);
                    MatrixUtil.multiply(tmp3, k[1]);
                    dcdistdT = MatrixUtil.pointwiseAdd(dcdistdT, tmp3);
                    // fall through
                }
                case 1 : {
                    tmp = Arrays.copyOf(r2, r2.length);
                    MatrixUtil.multiply(tmp, k[0]);
                    cdist = MatrixUtil.add(cdist, tmp);

                    tmp3 = MatrixUtil.copy(dr2dom);
                    MatrixUtil.multiply(tmp3, k[0]);
                    dcdistdom = MatrixUtil.pointwiseAdd(dcdistdom, tmp3);

                    tmp3 = MatrixUtil.copy(dr2dT);
                    MatrixUtil.multiply(tmp3, k[0]);
                    dcdistdT = MatrixUtil.pointwiseAdd(dcdistdT, tmp3);
                    break;
                }
                default :
                    throw new IllegalStateException("k.length not handled for " + k.length);
            }
        }

        //[n X 5]
        //dcdistdk = [ r2' r4' zeros(n,2) r6'];
        double[][] dcdistdk = new double[5][];
        dcdistdk[0] = Arrays.copyOf(r2, r2.length);
        dcdistdk[1] = Arrays.copyOf(r4, r4.length);
        dcdistdk[2] = new double[n];
        dcdistdk[3] = new double[n];
        dcdistdk[4] = Arrays.copyOf(r6, r6.length);
        dcdistdk = MatrixUtil.transpose(dcdistdk);

        // [2Xn] . [ [2x1]*[1Xn] ]
        //xd1 = x .* (ones(2,1)*cdist);
        double[][] tmp3 = MatrixUtil.outerProduct(new double[]{1, 1}, cdist);
        double[][] xd1 = MatrixUtil.pointwiseMultiplication(x, tmp3);

        //dxd1dom = zeros(2*n,3);
        //dxd1dom(1:2:end,:) = (x(1,:)'*ones(1,3)) .* dcdistdom; // [nX1][1X3] . [nX3]
        //dxd1dom(2:2:end,:) = (x(2,:)'*ones(1,3)) .* dcdistdom;
        d01 = MatrixUtil.outerProduct(x[0], new double[]{1, 1, 1});//(x(1,:)'*ones(1,3)) // [nX1][1X3]=nX3
        d02 = MatrixUtil.outerProduct(x[1], new double[]{1, 1, 1});//(x(2,:)'*ones(1,3))
        d1 = MatrixUtil.pointwiseMultiplication(d01, dcdistdom); // [nX3]
        d2 = MatrixUtil.pointwiseMultiplication(d01, dcdistdom);
        double[][] dxd1dom = new double[2*n][]; // [2nX3]
        for (i = 0; i < n; ++i) {
            dxd1dom[2*i] = Arrays.copyOf(d1[i], d1[i].length);
            dxd1dom[2*i + 1] = Arrays.copyOf(d2[i], d2[i].length);
        }
        //coeff = (reshape([cdist;cdist],2*n,1)*ones(1,3)); // cdist is [1Xn]
        //dxd1dom = dxd1dom + coeff.* dxdom;
        tmp1 = new double[2*n];
        for (i = 0; i < n; ++i) {
            tmp1[i] = cdist[i];
            tmp1[i + n] = cdist[i];
        }
        //[2*n X 3]
        double[][] coeff = MatrixUtil.outerProduct(tmp1, new double[]{1, 1, 1});
        tmp3 = MatrixUtil.pointwiseMultiplication(coeff, dxdom);
        dxd1dom = MatrixUtil.pointwiseAdd(dxd1dom, tmp3);

        //dxd1dT = zeros(2*n,3);
        //dxd1dT(1:2:end,:) = (x(1,:)'*ones(1,3)) .* dcdistdT; // [nX3] .* [nX3]
        //dxd1dT(2:2:end,:) = (x(2,:)'*ones(1,3)) .* dcdistdT;
        //dxd1dT = dxd1dT + coeff.* dxdT;
        d1 = MatrixUtil.pointwiseMultiplication(d01, dcdistdT); // [nX3]
        d2 = MatrixUtil.pointwiseMultiplication(d02, dcdistdT);
        double[][] dxd1dT = MatrixUtil.zeros(2*n, 3); // [2nX3]
        for (i = 0; i < n; ++i) {
            dxd1dT[2*i] = d1[i];
            dxd1dT[2*i + 1] = d2[i];
        }
        tmp3 = MatrixUtil.pointwiseMultiplication(coeff, dxdT);
        dxd1dT = MatrixUtil.pointwiseAdd(dxd1dT, tmp3);

        //dxd1dk = zeros(2*n,5);
        //dxd1dk(1:2:end,:) = (x(1,:)'*ones(1,5)) .* dcdistdk;
        //dxd1dk(2:2:end,:) = (x(2,:)'*ones(1,5)) .* dcdistdk;
        d01 = MatrixUtil.outerProduct(x[0], new double[]{1, 1, 1, 1, 1});
        d02 = MatrixUtil.outerProduct(x[1], new double[]{1, 1, 1, 1, 1});
        d1 = MatrixUtil.pointwiseMultiplication(d01, dcdistdk); // [nX5] [nX5]
        d2 = MatrixUtil.pointwiseMultiplication(d02, dcdistdk);
        double[][] dxd1dk = new double[2*n][]; // [2nX3]
        for (i = 0; i < n; ++i) {
            dxd1dk[2*i] = Arrays.copyOf(d1[i], d1[i].length);
            dxd1dk[2*i + 1] = Arrays.copyOf(d2[i], d2[i].length);
        }

        // excluding tangential distortion.  lines 161 - 191

        double[] deltaX = new double[2];

        //xd2 = xd1 + delta_x;
        //[2Xn]
        double[][] xd2 = MatrixUtil.copy(xd1);
        double[][] dxd2dom = MatrixUtil.copy(dxd1dom);// + ddelta_xdom ;
        double[][] dxd2dT = MatrixUtil.copy(dxd1dT);//+ ddelta_xdT;
        double[][] dxd2dk = MatrixUtil.copy(dxd1dk);// + ddelta_xdk ;

        //% Add Skew:
        double alpha = intrinsics.getIntrinsic()[0][1];

        //xd3 = [xd2(1,:) + alpha*xd2(2,:);xd2(2,:)];
        tmp1 = Arrays.copyOf(xd2[1], xd2[1].length);
        MatrixUtil.multiply(tmp1, alpha);
        double[][] xd3 = new double[2][]; // [2Xn]
        xd3[0] = MatrixUtil.add(xd2[0], tmp1);
        xd3[1] = Arrays.copyOf(xd2[1], xd2[1].length);

        //% Compute: dxd3dom, dxd3dT, dxd3dk, dxd3dalpha
        //dxd3dom = zeros(2*n,3);
        //dxd3dom(1:2:2*n,:) = dxd2dom(1:2:2*n,:) + alpha*dxd2dom(2:2:2*n,:);
        //dxd3dom(2:2:2*n,:) = dxd2dom(2:2:2*n,:);
        // dxd2dom is [2*n X 3]
        double[][] dxd3dom = new double[2*n][];
        for (i = 0; i < n; ++i) {
            tmp1 = Arrays.copyOf(dxd2dom[2*i + 1], dxd2dom[2*i + 1].length);
            MatrixUtil.multiply(tmp1, alpha);
            dxd3dom[2*i] = MatrixUtil.add(dxd2dom[2*i + 0], tmp1);
            dxd3dom[2*i + 1] = Arrays.copyOf(dxd2dom[2*i + 1], dxd2dom[2*i + 1].length);
        }

        //dxd3dT = zeros(2*n,3);
        //dxd3dT(1:2:2*n,:) = dxd2dT(1:2:2*n,:) + alpha*dxd2dT(2:2:2*n,:);
        //dxd3dT(2:2:2*n,:) = dxd2dT(2:2:2*n,:);
        double[][] dxd3dT = new double[2*n][];
        for (i = 0; i < n; ++i) {
            tmp1 = Arrays.copyOf(dxd2dT[2*i + 1], dxd2dT[2*i + 1].length);
            MatrixUtil.multiply(tmp1, alpha);
            dxd3dT[2*i] = MatrixUtil.add(dxd2dT[2*i + 0], tmp1);
            dxd3dT[2*i + 1] = Arrays.copyOf(dxd2dT[2*i + 1], dxd2dT[2*i + 1].length);
        }

        //dxd3dk = zeros(2*n,5);
        //dxd3dk(1:2:2*n,:) = dxd2dk(1:2:2*n,:) + alpha*dxd2dk(2:2:2*n,:);
        //dxd3dk(2:2:2*n,:) = dxd2dk(2:2:2*n,:);
        double[][] dxd3dk = new double[2*n][];
        for (i = 0; i < n; ++i) {
            tmp1 = Arrays.copyOf(dxd2dk[2*i + 1], dxd2dk[2*i + 1].length);
            MatrixUtil.multiply(tmp1, alpha);
            dxd3dk[2*i] = MatrixUtil.add(dxd2dk[2*i + 0], tmp1);
            dxd3dk[2*i + 1] = Arrays.copyOf(dxd2dk[2*i + 1], dxd2dk[2*i + 1].length);
        }

        //dxd3dalpha = zeros(2*n,1);
        //dxd3dalpha(1:2:2*n,:) = xd2(2,:)';
        // xd2 is [2*n X 1]
        double[] dxd3dalpha = new double[2*n];
        for (i = 0; i < n; ++i) {
            dxd3dalpha[i*2] = xd2[1][i];
        }

        //% Pixel coordinates:
        //if length(f)>1,
        boolean focalXYSame = Math.abs(intrinsics.getIntrinsic()[0][0] - intrinsics.getIntrinsic()[1][1]) < 1e-3;

        double[] c = new double[]{intrinsics.getIntrinsic()[0][2], intrinsics.getIntrinsic()[1][2]};

        double[][] xp;
        double[][] dxpdom, dxpdT, dxpdk, dxpdf;
        double[] dxpdalpha;

        if (!focalXYSame) {
            double[] f = new double[]{intrinsics.getIntrinsic()[0][0], intrinsics.getIntrinsic()[1][1]};
            //xp = xd3. * (f(:) *ones(1, n))+c(:)*ones(1, n);
            tmp1 = new double[n];
            Arrays.fill(tmp1, 1);
            tmp3 = MatrixUtil.outerProduct(f, tmp1);
            double[][] tmp4 = MatrixUtil.outerProduct(c, tmp1);
            xp = MatrixUtil.pointwiseAdd(MatrixUtil.pointwiseMultiplication(xd3, tmp3), tmp4);

            //coeff = reshape(f(:)*ones(1, n), 2 * n, 1);
            double[] coeff2 = MatrixUtil.stack(tmp3);  //[2*n X 1]
            //dxpdom = (coeff * ones(1, 3)). * dxd3dom;
            dxpdom = MatrixUtil.outerProduct(coeff2, new double[]{1, 1, 1}); // [2*n X 3]
            dxpdom = MatrixUtil.pointwiseMultiplication(dxpdom, dxd3dom);

            //dxpdT = (coeff * ones(1, 3)). * dxd3dT;
            dxpdT = MatrixUtil.outerProduct(coeff2, new double[]{1, 1, 1});
            dxpdT = MatrixUtil.pointwiseMultiplication(dxpdT, dxd3dT); // [2*n X 3]

            //dxpdk = (coeff * ones(1, 5)). * dxd3dk;  [2*n X 5]
            dxpdk = MatrixUtil.outerProduct(coeff2, new double[]{1, 1, 1, 1, 1});
            dxpdk = MatrixUtil.pointwiseMultiplication(dxpdk, dxd3dk);

            //dxpdalpha = (coeff). * dxd3dalpha;  //[2*n X 1] [2*n X 1]
            dxpdalpha = MatrixUtil.pointwiseMultiplication(coeff2, dxd3dalpha);
            //dxpdf = zeros(2 * n, 2);
            //dxpdf(1:2:end, 1) =xd3(1,:)';
            //dxpdf(2:2:end, 2) =xd3(2,:)';
            dxpdf = MatrixUtil.zeros(2*n, 2); //[2*nX2] for f as an array
            for (i = 0; i < n; ++i) {
                dxpdf[2*i][0] = xd3[0][i];
                dxpdf[2*i + 1][1] = xd3[1][i];
            }

        } else {

            //xp = f * xd3 + c * ones(1, n); // [2Xn]
            double f = intrinsics.getIntrinsic()[0][0];
            tmp1 = new double[n];
            Arrays.fill(tmp1, 1);
            tmp3 = MatrixUtil.outerProduct(c, tmp1);
            xp = MatrixUtil.copy(xd3);
            MatrixUtil.multiply(xp, f);
            xp = MatrixUtil.pointwiseAdd(xp, tmp3);

            //dxpdom = f * dxd3dom;  [2*n X 3]
            dxpdom = MatrixUtil.copy(dxd3dom);
            MatrixUtil.multiply(dxpdom, f);

            //dxpdT = f * dxd3dT;  [2*n X 3]
            dxpdT = MatrixUtil.copy(dxd3dT);
            MatrixUtil.multiply(dxpdT, f);

            //dxpdk = f * dxd3dk;  [2*n X 5]
            dxpdk = MatrixUtil.copy(dxd3dk);
            MatrixUtil.multiply(dxpdk, f);

            //dxpdalpha = f. * dxd3dalpha;
            dxpdalpha = Arrays.copyOf(dxd3dalpha, dxd3dalpha.length);
            MatrixUtil.multiply(dxpdalpha, f);

            // xd3 is [2Xn].  dxpdf for scalar f is [2*n X 1]
            //dxpdf = xd3(:);
            dxpdf = new double[1][];
            dxpdf[0] = MatrixUtil.stack(xd3);
            dxpdf = MatrixUtil.transpose(dxpdf);
        }

        //dxpdc = zeros(2*n,2);
        //dxpdc(1:2:end,1) = ones(n,1);
        //dxpdc(2:2:end,2) = ones(n,1);
        double[][] dxpdc = MatrixUtil.zeros(2*n, 2);
        for (i = 0; i < n; ++i) {
            dxpdc[2*i][0] = 1;
            dxpdc[2*i + 1][1] = 1;
        }

        // arrays that could be returned:
        // arrays xp,dxpdom,dxpdT,dxpdf,dxpdc,dxpdk,dxpdalpha

        ProjectedPoints pp = new ProjectedPoints();
        pp.xEst = xp;
        pp.dxdom = dxpdom;
        pp.dxdT = dxpdT;
        pp.dxdF = dxpdf;
        pp.dxdC = dxpdc;
        pp.dxdK = dxpdk;
        pp.dxdAlpha = dxpdalpha;

        return pp;
    }

    /**
     * Computes the rigid motion transformation Y = R*X+T, where R = rodrigues(om).
     * <pre>
     *     rigid_motion.m
     *     TODO: put other Bouguet toolbox references here
     * </pre>
     * @param X
     * @param om
     * @param T
     * @param useBouguetForRodrigues if true, uses only the Bouguet algoirthms for Rodrigues rotation matrices and vectors
     * @return Y = R*X+T, where R = rodrigues(om).  returns
     * Y: 3D coordinates of the structure points in the camera reference frame (3xN matrix for N points)
     * %        dYdom: Derivative of Y with respect to om ((3N)x3 matrix)
     * %        dYdT: Derivative of Y with respect to T ((3N)x3 matrix)
     */
    static ProjectedPoints bouguetRigidMotion(double[][] X, double[] om, double[] T, boolean useBouguetForRodrigues) {

        if (X.length != 3) {
            throw new IllegalArgumentException("X.length should be 3");
        }
        if (om.length != 3) {
            throw new IllegalArgumentException("om.length should be 3");
        }
        if (T.length != 3) {
            throw new IllegalArgumentException("T.length should be 3");
        }

        boolean passive = true;

        //[R, dRdom] = rodrigues(om);
        Rotation.RodriguesRotation rRot = Rotation.createRotationRodriguesBouguet(om, passive);
        double[][] r;
        if (useBouguetForRodrigues) {
            r = rRot.r;
        } else {
            r = Rotation.createRotationRodriguesFormula(om, passive);
        }

        //[m,n] = size(X);
        int m = X.length;
        int n = X[0].length;

        //Y = R*X + repmat(T,[1 n]);
        double[][] Y = MatrixUtil.multiply(r, X);
        int i, j;
        double tmp;
        for (j = 0; j < m; ++j) {
            tmp = T[j];
            for (i = 0; i < n; ++i) {
                Y[j][i] += tmp;
            }
        }

        //if nargout > 1
        //dYdR = zeros(3 * n, 9);
        //dYdT = zeros(3 * n, 3);
        double[][] dYdR = MatrixUtil.zeros(3*n, 9);
        double[][] dYdT = MatrixUtil.zeros(3*n, 3);

        //dYdR(1:3:end, 1:3:end) =X';  //dYdR[i*3]     cols 0, 3, 6
        //dYdR(2:3:end, 2:3:end) =X';  //dYdR[i*3 + 1] cols 1, 4, 7
        //dYdR(3:3:end, 3:3:end) =X';  //dYdR[i*3 + 2] cols 2, 5, 8
        int ii;
        int[] c0 = new int[]{0, 3, 6};
        int[] c1 = new int[]{1, 4, 7};
        int[] c2 = new int[]{2, 5, 8};
        for (i = 0; i < n; ++i) {
            for (ii = 0; ii < 3; ++ii) {
                dYdR[i * 3][c0[ii]] = X[ii][i];
                dYdR[i*3 + 1][c1[ii]] = X[ii][i];
                dYdR[i*3 + 2][c2[ii]] = X[ii][i];
            }
        }

        //dYdT(1:3:end, 1) =ones(n, 1);
        //dYdT(2:3:end, 2) =ones(n, 1);
        //dYdT(3:3:end, 3) =ones(n, 1);
        for (i = 0; i < n; ++i) {
            dYdT[i * 3][0] = 1;
            dYdT[i * 3 + 1][1] = 1;
            dYdT[i * 3 + 2][2] = 1;
        }

        //dYdom = dYdR * dRdom;
        //        [3*n X 9] [9X3] = [3*n X 3]
        double[][] dYdom = MatrixUtil.multiply(dYdR, rRot.dRdR);

        ProjectedPoints pp = new ProjectedPoints();
        pp.xEst = Y; //[3 X n]
        pp.dxdom = dYdom; //[3*n X 3]
        pp.dxdT = dYdT;  // [3*n X 3]

        return pp;
    }

    /**
     *
     * @param init
     * @param intrinsics
     * @param xi
     * @param X
     * @param useBouguetsRodrigues if true,uses the Bouguet algorithms for Rodrigues Rotation matrix and vector,
     *                             else, uses the other Rotation.java Rodrigues methods.
     *                             Recommend using false at this time.
     * @return
     * @throws NotConvergedException
     */
    public static CameraExtrinsicParameters bouguetPoseRefine(CameraExtrinsicParameters init,
          CameraIntrinsicParameters intrinsics, double[][] xi, double[][] X,
          boolean useBouguetsRodrigues) throws NotConvergedException {

        if (xi.length != 2 && xi.length != 3) {
            throw new IllegalArgumentException("xi length must be 3 or 2");
        }
        if (X.length != 3) {
            throw new IllegalArgumentException("X length must be3");
        }

        if (init == null || init.getRodriguesVector() == null || init.getTranslation() == null) {
            throw new IllegalArgumentException("inital solution must have the Rodrigues rotation vector" +
                    " and translation");
        }

        int i, j;

        //if size(m,1)<3,
        //    m = [m;ones(1,Np)];
        //end;
        //if size(M,1)<3,
        //    M = [M;ones(1,Np)];
        //end;
        //m = m ./ (ones(3,1)*m(3,:));
        //M = M ./ (ones(3,1)*M(3,:));

        //TODO: make this an argument or consider if it should be less than infinity
        double threshCond = Double.POSITIVE_INFINITY;
        int MaxIter = 20;

        //% Initialization:
        //omckk = omc_init;
        //Tckk = Tc_init;
        double[] omckk = Arrays.copyOf(init.getRodriguesVector(), init.getRodriguesVector().length);
        double[] Tckk = Arrays.copyOf(init.getTranslation(), init.getTranslation().length);

        int n = xi[0].length;

        //[2 X n]
        double[][] xkk = MatrixUtil.copySubMatrix(xi, 0, 1, 0, n - 1);

        //% Final optimization (minimize the reprojection error in pixel):
        //% through Gradient Descent:

        //param = [omckk;Tckk];  [6X1]
        double[] param = new double[6];
        System.arraycopy(omckk, 0, param, 0, 3);
        System.arraycopy(Tckk, 0, param, 3, 3);

        //change = 1;
        //iter = 0;
        double change = 1;
        int iter = 0;

        //%keyboard;
        //%fprintf(1,'Gradient descent iterations: ');

        ProjectedPoints pp;

        double[][] JJ = MatrixUtil.zeros(2*n, 6);

        while ((change > 1e-10) && (iter < MaxIter)) {

            //%fprintf(1,'%d...',iter+1);
            //[x,dxdom,dxdT] = project_points2(X_kk,omckk,Tckk,fc,cc,kc,alpha_c);
            pp = bouguetProjectPoints2(X, omckk, Tckk, intrinsics, useBouguetsRodrigues); // these are in image reference frame
            double[][] x = pp.xEst;  //[2 X n]
            double[][] dxdom = pp.dxdom; // [2*n X 3]
            double[][] dxdT = pp.dxdT; // [2*n X 3]

            //ex = x_kk - x; //both are [2Xn]
            double[][] ex = MatrixUtil.pointwiseSubtract(xkk, x);

            //%keyboard;

            //JJ = [dxdom dxdT];  // [2*n X 3]  [2*n X 3] = [2*n X 6]
            for (i = 0; i < 2*n; ++i) {
                System.arraycopy(dxdom[i], 0, JJ[i], 0, dxdom[i].length);
                System.arraycopy(dxdT[i], 0, JJ[i], dxdom[i].length, dxdT[i].length);
            }

            //Condition number of a matrix is the ratio of the largest singular value of that matrix to the smallest singular value.
            //if cond(JJ) > thresh_cond,
            if (MatrixUtil.conditionNumber(JJ) > threshCond) {
                change = 0;
            } else {
                //JJ2 = JJ'*JJ;  //[6 X 2*n][2*n X 6] = [6X6]
                double[][] JJ2 = MatrixUtil.createATransposedTimesA(JJ);

                //param_innov = inv(JJ2)*(JJ')*ex(:); //  [6X6][6 X 2*n] * [2*nX1] = [6 X 1]
                double[] paramInnov = MatrixUtil.multiplyMatrixByColumnVector(
                        MatrixUtil.multiply(MatrixUtil.pseudoinverseFullColumnRank(JJ2), MatrixUtil.transpose(JJ)),
                        MatrixUtil.stack(ex)
                );

                //           [6X1] + [6X1]
                //param_up = param + param_innov;
                double[] paramUp = MatrixUtil.add(param, paramInnov);

                //change = norm(param_innov)/norm(param_up);
                change = MatrixUtil.lPSum(paramInnov, 2) / MatrixUtil.lPSum(paramUp, 2);

                //param = param_up;
                System.arraycopy(paramUp, 0, param, 0, paramUp.length);
                //iter = iter + 1;
                ++iter;

                //omckk = param(1:3);
                System.arraycopy(param, 0, omckk, 0, omckk.length);

                //Tckk = param(4:6);
                System.arraycopy(param, 3, Tckk, 0, Tckk.length);

            }// end if
        } // end while

        System.out.printf("bouguet refine iter=%d\n", iter);
        //%fprintf(1,'\n');

        boolean passive = true;

        //Rckk = rodrigues(omckk);
        double[][] r;
        if (useBouguetsRodrigues) {
            r = Rotation.createRotationRodriguesBouguet(omckk, passive).r;
        } else {
            r = Rotation.createRotationRodriguesFormula(omckk, passive);
        }

        //can return [omckk,Tckk,Rckk,JJ]
        CameraExtrinsicParameters extr = new CameraExtrinsicParameters(r, omckk, Tckk);
        return extr;
    }

    /**
     * calc a rotation (ambiguous) and translation between the measurements of a point in the real world.
     * the lists xC and X are correspondences of image and object.
     *
     * https://github.com/fragofer/TOOLBOX_calib/
     * compute_extrinsic_init.m
     * @param intrinsics
     * @param xc objects in camera coordinates
     * @param X objects in real world coordinates
     * @param vT formed from the SVD of X (hence, the last column is orthogonal to X)
     * @param XMean
     * @param useBouguetForRodrigues if true, uses only the Bouguet algoirthms for Rodrigues rotation matrices and vectors
     * @return
     * @throws NotConvergedException
     * @throws IOException
     */
    static CameraExtrinsicParameters bouguetPoseInitPlanar(
            Camera.CameraIntrinsicParameters intrinsics, double[][] xc,
            double[][] X, DenseMatrix vT, double[] XMean, boolean useBouguetForRodrigues) throws NotConvergedException, IOException {

        if (xc.length != 2 && xc.length != 3) {
            throw new IllegalArgumentException("xc length must be 3 or 2");
        }
        if (X.length != 2 && X.length != 3 && X.length != 4) {
            throw new IllegalArgumentException("X length must be 2, 3, or 4");
        }

        boolean passive = true;

        int i, j;

        //if size(m,1)<3,
        //    m = [m;ones(1,Np)];
        //end;
        //if size(M,1)<3,
        //    M = [M;ones(1,Np)];
        //end;
        //m = m ./ (ones(3,1)*m(3,:));
        //M = M ./ (ones(3,1)*M(3,:));

        if (xc.length == 2) {
            double[][] x2 = new double[3][];
            x2[0] = Arrays.copyOf(xc[0], xc[0].length);
            x2[1] = Arrays.copyOf(xc[1], xc[1].length);
            x2[2] = new double[xc[0].length];
            Arrays.fill(x2[2], 1);
            xc = x2;
        } else {
            xc = MatrixUtil.copy(xc);
            // normalize by last coordinate just in case not performed already:
            for (i = 0; i < xc[0].length; ++i) {
                for (j = 0; j < xc.length; ++j) {
                    xc[j][i] /= xc[xc.length - 1][i];
                }
            }
        }

        if (X.length == 2) {
            double[][] x2 = new double[3][];
            x2[0] = Arrays.copyOf(X[0], X[0].length);
            x2[1] = Arrays.copyOf(X[1], X[1].length);
            x2[2] = new double[X[0].length];
            Arrays.fill(x2[2], 1);
            X = x2;
        }

        if (X.length == 4) {
            for (i = 0; i < X[0].length; ++i) {
                for (j = 0; j < X.length; ++j) {
                    X[j][i] /= X[X.length - 1][i];
                }
            }
        }

        // planar structure
        //Transform the plane to bring it in the Z=0 plane:

        int n = xc[0].length;

        //R_transform = V';  [3X3]
        double[][] Rtransform = MatrixUtil.convertToRowMajor(vT); // orthogonal to X...
        //%norm(R_transform(1:2,3))

        //if norm(R_transform(1:2,3)) < 1e-6,
        //    R_transform = eye(3);
        //end;
        double norm0 = MatrixUtil.lPSum(new double[]{Rtransform[0][2], Rtransform[1][2]},2);
        if (norm0 < 1e-6) {
            Rtransform = MatrixUtil.createIdentityMatrix(3);
        }

        //if det(R_transform) < 0, R_transform = -R_transform; end;
        if (MatrixUtil.determinant(Rtransform) < 0) {
            MatrixUtil.multiply(Rtransform, -1);
        }

        // [3X3][3 X 1] = [3 X 1]
        //T_transform = -(R_transform)*X_mean;
        double[] Ttransform = MatrixUtil.multiplyMatrixByColumnVector(Rtransform, XMean);
        MatrixUtil.multiply(Ttransform, -1);

        //   [3X3] [3Xn] = [3Xn];   [3X1][1Xn] = 3Xn
        //X_new = R_transform*X_kk + T_transform*ones(1,Np);
        double[] ones = new double[n];
        Arrays.fill(ones, 1);
        double[][] t2 = MatrixUtil.outerProduct(Ttransform, ones);
        double[][] Xnew = MatrixUtil.multiply(Rtransform, X);
        Xnew = MatrixUtil.pointwiseAdd(Xnew, t2);

        //% Compute the planar homography:

        //H = compute_homography(xn,X_new(1:2,:));
        //NLK: replace Xnew[2] with 1's because we are giving the method only the first
        // 2 rows of Xnew, then compute_homography.m when receiving Xnew of length 2,
        // appends a row of 1's in the Matlab code.
        Arrays.fill(Xnew[2], 1);
        double[][] H = CameraCalibration.solveForHomographyBouget(xc, Xnew);

        //% De-embed the motion parameters from the homography:
        // Matlab norm of a vector is a euclidean norm
        //sc = mean([norm(H(:,1));norm(H(:,2))]);
        norm0 = MatrixUtil.lPSum(MatrixUtil.extractColumn(H, 0), 2);
        double norm1 = MatrixUtil.lPSum(MatrixUtil.extractColumn(H, 1), 2);
        double sc = (norm0 + norm1)/2.;

        //H = H/sc;
        MatrixUtil.multiply(H, 1./sc);

        //u1 = H(:,1);
        //u1 = u1 / norm(u1);
        double[] u1 = MatrixUtil.extractColumn(H, 0);
        MatrixUtil.multiply(u1, 1./MatrixUtil.lPSum(u1, 2));

        //u2 = H(:,2) - dot(u1,H(:,2)) * u1;
        //u2 = u2 / norm(u2);
        double[] tu2 = MatrixUtil.extractColumn(H, 1);
        double d = MatrixUtil.dot(u1, tu2);
        double[] tu1 = Arrays.copyOf(u1, u1.length);
        MatrixUtil.multiply(tu1, d);
        double[] u2 = new double[u1.length];
        MatrixUtil.pointwiseSubtract(tu2, tu1, u2);

        //u3 = cross(u1,u2);
        double[] u3 = MatrixUtil.crossProduct(u1, u2);
        //RRR = [u1 u2 u3];
        double[][] RRR = new double[3][];
        RRR[0] = u1;
        RRR[1] = u2;
        RRR[2] = u3;
        RRR = MatrixUtil.transpose(RRR);
        RRR = Rotation.orthonormalizeUsingSVD(RRR);

        //omckk = rodrigues(RRR);
        double[] omckk;
        if (useBouguetForRodrigues) {
            Rotation.RodriguesRotation rRot = Rotation.extractRotationVectorRodriguesBouguet(RRR, passive);
            omckk = rRot.om;
        } else {
            omckk = Rotation.extractRotationVectorRodrigues(RRR);
        }

        //%omckk = rodrigues([H(:,1:2) cross(H(:,1),H(:,2))]);
        //Rckk = rodrigues(omckk);
        double[][] Rckk;
        if (useBouguetForRodrigues) {
            Rotation.RodriguesRotation rRot2 = Rotation.createRotationRodriguesBouguet(omckk, passive);
            Rckk = rRot2.r;
        } else {
            Rckk = Rotation.createRotationRodriguesFormula(omckk, passive);
        }

        //Tckk = H(:,3);
        double[] Tckk = MatrixUtil.extractColumn(H, 2);
        System.out.printf("T_transform of X from its origin=\n%s\n", FormatArray.toString(Ttransform, "%.4e"));
        System.out.printf("Tckk derived from homography between x and X_origin =\n%s\n", FormatArray.toString(Tckk, "%.4e"));

        //%If Xc = Rckk * X_new + Tckk, then Xc = Rckk * R_transform * X_kk + Tckk + T_transform
        //NLK: Xc = Rckk * (R_transform * X_kk + T_transform) + Tckk
        //Tckk = Tckk + Rckk* T_transform;
        Tckk = MatrixUtil.add(Tckk, MatrixUtil.multiplyMatrixByColumnVector(Rckk, Ttransform));
        System.out.printf("Tckk += Rckk* T_transform = \n%s\n", FormatArray.toString(Tckk, "%.4e"));
        //Rckk = Rckk * R_transform;
        Rckk = MatrixUtil.multiply(Rckk, Rtransform);
        //omckk = rodrigues(Rckk);
        if (useBouguetForRodrigues) {
            omckk = Rotation.extractRotationVectorRodriguesBouguet(Rckk, passive).om;
            //Rckk = rodrigues(omckk);
            Rckk = Rotation.createRotationRodriguesBouguet(omckk, passive).r;
        } else {
            omckk = Rotation.extractRotationVectorRodrigues(Rckk);
            // this should be the same.  TODO: follow up on simplifying this method w.o. losing accuracy though
            Rckk = Rotation.createRotationRodriguesFormula(omckk, passive);
        }

        return new CameraExtrinsicParameters(Rckk, omckk, Tckk);
    }

    /**
     *
     * @param xc
     * @param X
     * @param useBouguetForRodrigues if true, uses only the Bouguet algoirthms for Rodrigues rotation matrices and vectors
     * @return
     * @throws NotConvergedException
     */
    static CameraExtrinsicParameters bouguetPoseInitNonPlanar(double[][] xc,
            double[][] X, boolean useBouguetForRodrigues) throws NotConvergedException {

        if (xc.length != 2 && xc.length != 3) {
            throw new IllegalArgumentException("xc length must be 3 or 2");
        }
        if (X.length != 2 && X.length != 3 && X.length != 4) {
            throw new IllegalArgumentException("X length must be 2, 3, or 4");
        }

        int i, j;

        //if size(m,1)<3,
        //    m = [m;ones(1,Np)];
        //end;
        //if size(M,1)<3,
        //    M = [M;ones(1,Np)];
        //end;
        //m = m ./ (ones(3,1)*m(3,:));
        //M = M ./ (ones(3,1)*M(3,:));

        if (xc.length == 2) {
            double[][] x2 = new double[3][];
            x2[0] = Arrays.copyOf(xc[0], xc[0].length);
            x2[1] = Arrays.copyOf(xc[1], xc[1].length);
            x2[2] = new double[xc[0].length];
            Arrays.fill(x2[2], 1);
            xc = x2;
        } else {
            xc = MatrixUtil.copy(xc);
            // normalize by last coordinate just in case not performed already:
            for (i = 0; i < xc[0].length; ++i) {
                for (j = 0; j < xc.length; ++j) {
                    xc[j][i] /= xc[xc.length - 1][i];
                }
            }
        }

        if (X.length == 2) {
            double[][] x2 = new double[3][];
            x2[0] = Arrays.copyOf(X[0], X[0].length);
            x2[1] = Arrays.copyOf(X[1], X[1].length);
            x2[2] = new double[X[0].length];
            Arrays.fill(x2[2], 1);
            X = x2;
        }

        if (X.length == 4) {
            for (i = 0; i < X[0].length; ++i) {
                for (j = 0; j < X.length; ++j) {
                    X[j][i] /= X[X.length - 1][i];
                }
            }
        }

        //% Computes an initial guess for extrinsic parameters (works for general 3d structure, not planar!!!):
        //% The DLT method is applied here!!

        int n = xc[0].length;

        //J = zeros(2*Np,12);
        double[][] J = new double[2*n][];

        //xX = (ones(3,1)*xn(1,:)).*X_kk;
        //yX = (ones(3,1)*xn(2,:)).*X_kk;
        //J(1:2:end,[1 4 7]) = -X_kk';
        //J(2:2:end,[2 5 8]) = X_kk';
        //J(1:2:end,[3 6 9]) = xX';
        //J(2:2:end,[3 6 9]) = -yX';
        //J(1:2:end,12) = xn(1,:)';
        //J(2:2:end,12) = -xn(2,:)';
        //J(1:2:end,10) = -ones(Np,1);
        //J(2:2:end,11) = ones(Np,1);

        for (i = 0; i < n; ++i) {
            J[i * 2] = new double[]{-X[0][i], 0, (xc[0][i] * X[0][i]), -X[1][i],
                    0, xc[0][i] * X[1][i], -X[2][i], 0, xc[0][i] * X[2][i], -1, 0, xc[0][i]};
            J[i * 2 + 1] = new double[]{0, X[0][i], -xc[1][i] * X[0][i], 0, X[1][i],
                    -xc[1][i] * X[1][i], 0, X[2][i],
                    -xc[1][i] * X[2][i], 0, 1, -xc[1][i]};
        }

        //JJ = J'*J; [12 X 12]
        double[][] JTJ = MatrixUtil.createATransposedTimesA(J);
        //[U,S,V] = svd(JJ);
        SVD svd = SVD.factorize(new DenseMatrix(JTJ));

        //RR = reshape(V(1:9,12),3,3);
        double[][] Vt = MatrixUtil.convertToRowMajor(svd.getVt()); //[12 X 12]
        double[] orth = Vt[11];
        // reshape fills each column first, then next column, etc, so will fill by rows then transpose
        double[][] RR = new double[3][];
        RR[0] = Arrays.copyOfRange(orth, 0, 3);
        RR[1] = Arrays.copyOfRange(orth, 3, 6);
        RR[2] = Arrays.copyOfRange(orth, 6, 9);
        RR = MatrixUtil.transpose(RR);

        //if det(RR) < 0,
        //        V(:,12) = -V(:,12);
        //RR = -RR;
        //end;
        if (MatrixUtil.determinant(RR) < 0) {
            MatrixUtil.multiply(orth, -1);
            MatrixUtil.multiply(RR, -1);
        }

        //[Ur,Sr,Vr] = svd(RR);
        svd = SVD.factorize(new DenseMatrix(RR));
        double[][] Ur = MatrixUtil.convertToRowMajor(svd.getU());
        double[][] Vrt = MatrixUtil.convertToRowMajor(svd.getVt());

        //Rckk = Ur*Vr';
        double[][] Rckk = MatrixUtil.multiply(Ur, Vrt);

        //sc = norm(V(1:9,12)) / norm(Rckk(:));
        double norm0 = MatrixUtil.lPSum(Arrays.copyOfRange(orth, 0, 9), 2);
        double norm1 = MatrixUtil.lPSum(MatrixUtil.stack(Rckk), 2);
        double sc = norm0/norm1;

        boolean passive = true;

        //Tckk = V(10:12,12)/sc;
        double[] Tckk = Arrays.copyOfRange(orth, 9, 12);
        MatrixUtil.multiply(Tckk, 1./sc);
        //omckk = rodrigues(Rckk);
        double[] omckk;
        if (useBouguetForRodrigues) {
            omckk = Rotation.extractRotationVectorRodriguesBouguet(Rckk, passive).om;
            //Rckk = rodrigues(omckk);
            Rckk = Rotation.createRotationRodriguesBouguet(omckk, passive).r;
        } else {
            omckk = Rotation.extractRotationVectorRodrigues(Rckk);
            Rckk = Rotation.createRotationRodriguesFormula(omckk, passive);
        }

        return new CameraExtrinsicParameters(Rckk, omckk, Tckk);
    }
    
    /**
     * NOT YET IMPLEMENTED.
     * given n 3D-to-2D point correspondences, estimates the pose 
     * of a calibrated camera (a.k.a. P-n-P) with computational complexity O(n)
     * using the Moreno-Noguer et al. 2007 non-iterative algorithm.
     * This could be followed by non-linear optimization
     * to improve the parameter estimates.
     * <pre>
     * references:
     * Moreno-Noguer, Lepetite, & Fua 2007, "Accurate Non-Iterative O(n) Solution to the PnP Problem"
     * Szeliski 2010 draft of "Computer Vision: Algorithms and Applications"
     * </pre>
     * @param intrinsics
     * @param x the image coordinates of the features in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 3 features are needed to 
     * calculate the extrinsic parameters.
     * NOTE x and X should both be distortion-free or both should be distorted.
     * @param X the world coordinates of the features in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 3features are needed to 
     * calculate the extrinsic parameters.
     * NOTE x and X should both be distortion-free or both should be distorted.
     @return 
     */
    /*public static CameraExtrinsicParameters calculatePoseUsingPNP(
        Camera.CameraIntrinsicParameters intrinsics, double[][] x,
        double[][] X) throws NotConvergedException {
                
        if (x.length != 3) {
            throw new IllegalArgumentException("x.length must be 3");
        }
        if (X.length != 3) {
            throw new IllegalArgumentException("X.length must be 3");
        }
        int n = x[0].length;
        
        if (n < 4) {
            throw new IllegalArgumentException("x must have at least 4 correspondences");
        }
        if (X[0].length != n) {
            throw new IllegalArgumentException("the number of columns in X must be the same as in x");
        }
        
        // Szeliski 2010 refers to perspective-n-point-problem (PnP) references  
        //   (Haralick, Lee, Ottenberg et al. 1994; Quan and Lan 1999; Moreno-Noguer, Lepetit, and Fua 2007)
        
        //port the c++ impl of  Moreno-Noguer, Lepetit, and Fua (2007)  here?
        //https://github.com/cvlab-epfl/EPnP/tree/master/cpp   
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    */


}
