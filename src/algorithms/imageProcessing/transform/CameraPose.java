package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Camera.CameraExtrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraIntrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraPoseParameters;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import no.uib.cipr.matrix.*;

import java.util.Arrays;

/**
 given a set of features in image coordinates and world coordinate space with
  known camera intrinsic parameters, estimate the camera pose, that is
  extract the camera extrinsic parameters.
 <em>See also PNP.java</em>
 
 * TODO: consider solving with M-estimators.
 * see http://research.microsoft.com/en- us/um/people/zhang/INRIA/Publis/Tutorial-Estim/node24.html
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
     * given a set of features in image space and world coordinate space,
     * estimate the camera pose, that is
     * extract the camera extrinsic parameters of rotation and translation and also the camera
     * intrinsic parameters.
     * This is also known as estimating the Motion.
     * This method uses DLT and could be followed by non-linear optimization
     *      * to improve the parameter estimates.
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
     While reading this, keep in mind that this method needs x in image reference frame (units os pixels).
     The case of camera coordinates is explained also, butshould be applied to calculatePoseUsingDLT().

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
     * @param x the image coordinates of the features in pixels in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 6 features are needed.
     * NOTE x and X should both be distortion-free or both should be distorted.
     * @param X the world coordinates of the features in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 6 features are needed.
     * NOTE x and X should both be distortion-free or both should be distorted.
     * @return 
     */
    public static CameraPoseParameters calculatePoseAndKUsingDLT(double[][] x, double[][] X)
        throws NotConvergedException {
                
        double[][] p = calculatePFromXXW(x, X);

        /*
        MatrixUtil.SVDProducts svdP = MatrixUtil.performSVD(p);
        double[] c = Arrays.copyOf(svdP.vT[svdP.vT.length - 1], svdP.vT[0].length);
        // assert P*c = 0
        double[] check0 = MatrixUtil.multiplyMatrixByColumnVector(p, c);
        System.out.printf("check that P*c=0:%s\n", FormatArray.toString(check0, "%.3e"));
        */

        /*
        Szeliski Sect 6.2.1
        Since K is by convention upper-triangular 
        (see the discussion in Section 2.1.5), both K and R can be obtained 
        from the front 3 ⇥ 3 sub-matrix of P using RQ factorization 
        (Golub and Van Loan 1996)
        */

        // [3X3]
        double[][] M = MatrixUtil.copySubMatrix(p, 0, 2, 0, 2);
        double[][] invM = MatrixUtil.pseudoinverseFullColumnRank(M);
        QR qr = QR.factorize(new DenseMatrix(invM)); // Q=rot^T R=K^-1
        double[][] rot = MatrixUtil.transpose(MatrixUtil.convertToRowMajor(qr.getQ()));
        double[][] k = MatrixUtil.pseudoinverseFullColumnRank(MatrixUtil.convertToRowMajor(qr.getR()));
        MatrixUtil.multiply(k, 1./k[2][2]);

        double[][] rzpi = MatrixUtil.createIdentityMatrix(3);
        rzpi[0][0] = -1;
        rzpi[1][1] = -1;

        rot = MatrixUtil.multiply(rzpi, rot);
        k = MatrixUtil.multiply(k, rzpi);

        // TODO: presumably, should then orthonormalize rot, but then need to apply
        // similar multiplication to right multiply of k
        // see Rotation.orthonormalizeUsingSVD or Rotation.orthonormalizeUsingSkewCayley.

        // this assumes xc=R*(xw-t) instead of xc=R*xw + t
        // calculates: last column of P = P[*][3] = -K*R*X_0 where X_0 is projection center of camera.
        // X_0 = -1 * R^-1 * K^-1 * P[*][3]
        double[] projectionCenter = MatrixUtil.multiplyMatrixByColumnVector(invM,
                MatrixUtil.extractColumn(p, 3));
        MatrixUtil.multiply(projectionCenter, -1);

        double[] p3 = MatrixUtil.extractColumn(p, 3);

        CameraExtrinsicParameters extrinsics = new CameraExtrinsicParameters();
        extrinsics.setRotation(rot);
        extrinsics.setTranslation(projectionCenter);

        CameraIntrinsicParameters intrinsics = new CameraIntrinsicParameters();
        intrinsics.setIntrinsic(k);

        CameraPoseParameters camera = new CameraPoseParameters(intrinsics, extrinsics, p3);
        
        return camera;
    }

    /**
     * calculate projection matrix P from x = P * X.
     <pre>
     x = P * X where x and X are homogeneous coordinates (normalized so that last item = 1).
             p11*X[0] + p12*X[1] + p13*X[2] + p14
     x[0] = -------------------------------------
             p31*X[0] + p32*X[1] + p32*X[2] + p34

             p21*X[0] + p22*X[1] + p23*X[2] + p24
     x[1] = -------------------------------------
             p31*X[0] + p32*X[1] + p32*X[2] + p34

     rewrite in terms of factoring p members for DLT:
     -p11*X[0] - p12*X[1] - p13*X[2] - p14                                       + x[0]*(p31*X[0] + p32*X[1] + p32*X[2] + p34) = 0
                                           -p21*X[0] - p22*X[1] - p23*X[2] - p24 + x[1]*(p31*X[0] + p32*X[1] + p32*X[2] + p34) = 0

    A =  [ -X[0], -X[1], -X[2], -1,   0,        0,    0,   0,  x[0]*X[0], x[0]*X[1], x[0]*X[2], x[0] ]   *  [p11, p12, p13, p14, p21, p22, ...]^T
         [ 0,      0,     0,  0,  -X[0], -X[1], -X[2], -1,  x[1]*X[0], x[1]*X[1], x[1]*X[2], x[1]    ]

      A * p = 0

     and svd for least squares fit.
     </pre>
     * @param x coordinates in camera or image reference frame of the objects in X.  need at least 6 points
     * @param X coordinates of objects in world reference frame.
     * @return
     */
    public static double[][] calculatePFromXXW(double[][] x, double[][] X) throws NotConvergedException {
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

        x = MatrixUtil.copy(x);
        X = MatrixUtil.copy(X);

        int i, j;

        // normalize by last coordinate just in case nor performed already:
        for (i = 0; i < x[0].length; ++i) {
            for (j = 0; j < x.length; ++j) {
                x[j][i] /= x[x.length - 1][i];
            }
        }
        for (i = 0; i < X[0].length; ++i) {
            for (j = 0; j < X.length; ++j) {
                X[j][i] /= X[X.length - 1][i];
            }
        }

        // 2*n X 12
        double xi, yi, Xi, Yi, Zi;
        double[][] ell = new double[2*n][12];
        for (i = 0; i < n; ++i) {
            xi = x[0][i];
            yi = x[1][i];
            Xi = X[0][i];
            Yi = X[1][i];
            Zi = X[2][i];
            ell[2*i]     = new double[]{-Xi, -Yi, -Zi, -1, 0, 0, 0, 0, xi*Xi, xi*Yi, xi*Zi, xi};
            ell[2*i + 1] = new double[]{0, 0, 0, 0, -Xi, -Yi, -Zi, -1, yi*Xi, yi*Yi, yi*Zi, yi};
        }

        // SVD.V is 12 X 12.  SVD.U is 2n X 2n
        MatrixUtil.SVDProducts svd = MatrixUtil.performSVD(ell);

        // vT is 12X12.  last row in vT is the eigenvector for the smallest eigenvalue
        // it is also the epipole e1, defined as the right nullspace
        double[] xOrth = svd.vT[svd.vT.length - 1];

        // reshape into 3 X 4
        double[][] P2 = MatrixUtil.zeros(3,4);
        System.arraycopy(xOrth, 0, P2[0], 0, 4);
        System.arraycopy(xOrth, 4, P2[1], 0, 4);
        System.arraycopy(xOrth, 8, P2[2], 0, 4);

        return P2;
    }

    /**
     * given a set of features in image space and world coordinate space with
     * known camera intrinsic parameters, estimate the camera pose, that is
     * extract the camera extrinsic parameters.
     * calibrating the camera extrinsic parameters is a.k.a. 
     * perspective-n-point-problem where n is the number of features (a.k.a. points).
     * It's also called planar homography decomposition.
     * This method uses planar homography DLT and could be followed by non-linear optimization
     * to improve the parameter estimates.
     * Note that the projective matrix assumed is P = [K*R|K*t] from  x_im = K * X_c = K * R * X_wcs + t.
     If you instead are using the convention x_im = K * X_c = K * R * (X_wcs - t),
     the projection matrix would contain P = [K*R | -K*R*t] and so to correct the translation to
     your convention, you can calculate t2 = -1 * pseudoInv(R) * t where t is the returned translation from this method.
     * Note that for a proper rotation matrix that is orthogonal and unitary, one can use R^T for inv(R).
     * <pre>
     * references:
     * Ma, Chen, & Moore 2003 "Camera Calibration: a USU Implementation"
     * http://www.cs.cmu.edu/~16385/s17/Slides/11.3_Pose_Estimation.pdf
     * Zhang 1999, "Flexible Camera Calibration By Viewing a Plane From Unknown Orientations"
     * Szeliski 2010 draft of "Computer Vision: Algorithms and Applications"
     * </pre>
     * @param intrinsics
     * @param x the image coordinates of the features in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 3 features are needed to 
     * calculate the extrinsic parameters.
     * NOTE x and X should both be distortion-free or both should be distorted.
     * @param X the world coordinates of the features in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 3 features are needed to
     * calculate the extrinsic parameters.
     * NOTE x and X should both be distortion-free or both should be distorted.
     * @return 
     */
    public static CameraExtrinsicParameters calculatePoseUsingCameraCalibration(
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
        
        int nImages = 1;

        x = MatrixUtil.copy(x);
        X = MatrixUtil.copy(X);

        int i, j;

        // normalize by last coordinate just in case nor performed already:
        for (i = 0; i < x[0].length; ++i) {
            for (j = 0; j < x.length; ++j) {
                x[j][i] /= x[x.length - 1][i];
            }
        }
        for (i = 0; i < X[0].length; ++i) {
            for (j = 0; j < X.length; ++j) {
                X[j][i] /= X[X.length - 1][i];
            }
        }
        
        // following Ma et al. 2003
        double[][] h = CameraCalibration.solveForHomography(x, X);

        CameraExtrinsicParameters kExtr = CameraCalibration.solveForExtrinsic(
            intrinsics, h);
        
        return kExtr;
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
