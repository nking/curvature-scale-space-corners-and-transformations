package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Camera.CameraExtrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraIntrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraMatrices;
import algorithms.matrix.MatrixUtil;
import algorithms.matrix.MatrixUtil.SVDProducts;
import algorithms.misc.MiscMath;
import algorithms.misc.MiscMath0;
import algorithms.misc.PolynomialRootSolver;
import algorithms.util.FormatArray;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import no.uib.cipr.matrix.*;

/**
 * estimate the camera intrinsic and extrinsic parameters using 3 images
 * or more of the same objects with different camera poses.
 * Following the algorithm of Ma, Chen, & Moore 2003 in 
 * "Camera Calibration: a USU Implementation" available as a preprint
 * at arXiv  https://arxiv.org/pdf/cs/0307072
 * Note that Ma et al. 2003 algorithm is based upon Zhang 1999 
 * "Flexible Camera Calibration By Viewing a Plane From Unknown Orientations"
 * available at https://www.microsoft.com/en-us/research/wp-content/uploads/2016/11/zhan99.pdf
 * 
 * <pre>
 * "one image observed by a camera can provide 2 constraints about this camera’s 
 * intrinsic parameters that are regarded to be unchanged here. 
 * With 3 images observed by the same camera, 6 constraints are established and 
 * we are able to recover the 5 intrinsic parameters. Once the intrinsic 
 * parameters are known, we can estimate the extrinsic parameters, 
 * the distortion coefficients (k1,k2), and put every initial guess of these 
 * parameters into some nonlinear optimization routine to get the final estimations. 
 * </pre>
 <pre>
   Note, if need to estimate the instrinsic camera for initial conditions, one can generally
   start with:
   MASKS Algorithm 11.6, step 1:
   Guess a calibration matrix K by choosing the optical center at the center of the image,
   assuming the pixels to be square, and guessing the focal length f. For example, for
   an image plane of size (Dx X Dy) pixels, a typical guess is
       | f O Dx/2
   K = | 0 f Dy/2
       | 0 0 1
   with f = k X Dx, where k is typically chosen in the interval [0.5, 2].
 </pre>
 * @author nichole
 */
public class CameraCalibration {
    
    public static final double eps = 1e-7;
    public static final double eps2 = 1e-5;
    
    private static final Level LEVEL = Level.INFO;//Level.FINEST;
    private static final Logger log;
    static {
        log = Logger.getLogger(CameraCalibration.class.getSimpleName());
        log.setLevel(LEVEL);
    }

    /**
     * for a camera at rest and images taken of a geometric plane in motion,
     * estimate the intrinsic and extrinsic camera parameters.
     <pre>
     The method follows Zhang 1999 (
     Z. Zhang. A flexible new technique for camera calibration.
     IEEE Trans. Pattern Analysis and Machine Intelligence, 22(11):1330–1334, 2000.)
     and Zhang 2004, Chap 2 Camera Calibration
      in "Emerging Topics in Computer Vision" by Medioni & Kang

     </pre>
     * @param n n is the number of points in each image which is the
              same for all images.  n is number of features.
     * @param coordsI  holds the homogeneous image coordinates in pixels of
               features present in all images ordered in the same
               manner and paired with features in coordsW.
               It is a 2-dimensional double array of format
               3 X (N*n) where N is the number of images.
               the first row is the x coordinates, the second row
               is the y coordinates, and the third row is "1"'s.
               The columns hold each image in order and within each image's
               columns are the features presented in the same order in each image.
               In Table 1 of Ma, Chen, & Moore 2003 "Camera Calibration"
               these are the (u_d, v_d) pairs.
    <pre>
    e.g. coordsI = [ [xim1[0], xim1[1], ...xim1[n-1], xim2[0], xim2[1],...xim2[n-1],..., ]
                    [yim1[0], yim1[1], ...yim1[n-1], yim2[0], yim2[1],...yim2[n-1],...]
                    [1, 1, ...1, 1, 1,...1,...]]
    </pre>
     * @param coordsW holds the homogeneous world coordinates of features, ordered
               by the same features in the images.
               the first row is the X coordinates, the second row
               is the Y coordinates, and the third row is 1's 
               (Z_w = 0, the scale factor is lost in the homography).
               It is a 2-dimensional double array of format
               3 X n
       @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
    f(r) = 1 +k1*r^2 + k2*r^4 if true,
    else use model #3 f(r) = 1 +k1*r + k2*r^2.
    * @throws Exception if there is an error in use of MPSolver during the
     * removal of radial distortion, a generic exception is thrown with the
     * error message from the MPSolve documentation.
     * @return camera intrinsic parameters, extrinsic parameters, and radial
     * distortion coefficients
     */
    public static CameraMatrices estimateCameraPlanar(final int n, double[][] coordsI,
                                                      double[][] coordsW, boolean useR2R4) throws NotConvergedException,
        Exception {
        
        if (coordsI.length != 3) {
            throw new IllegalArgumentException("coordsI must have 3 rows.");
        }
        int nImages = coordsI[0].length/n;
        if (coordsI[0].length != nImages*n) {
            throw new IllegalArgumentException("coordsI must have nImages * n features as the number of columns");
        }
        if (coordsW.length != 3) {
            throw new IllegalArgumentException("coordsW must have 3 rows.");
        }

        //TODO: add optimization steps consider using PNP solve for pose

        boolean passive = true;
        
        int i;

        //TODO: consider normalizing coordsI by coordsI[2][*] if the last row
        //      is not already 1's

        boolean useNormConditioning = true;
        
        //(1) for each image: invoke homography solver using SVD on a (2N)X9 DLT
        //    where the homographies are the projections of the 3D points onto the images
        //    and the model plane has Z=0 (hence the (2N)X9 DLT instead of (2N)X12 DLT)
        double[][] h = solveForHomographies(coordsI, coordsW, n, nImages, useNormConditioning);

        //(2) using all homographies, solve for the camera intrinsic parameters
        // this is where at least 3 points are needed per image to equal the number of unknown intrinsic parameters.
        CameraIntrinsicParameters kIntr = solveForIntrinsicPlanar(h);

        double[] u = new double[n*nImages];
        double[] v = new double[n*nImages];
        List<CameraExtrinsicParameters> extrinsics = new ArrayList<>();
        int j, k;
        for (i = 0, j=i, k=i ; i < nImages; ++i, j+=3, k+=n) {

            double[][] hI = MatrixUtil.copySubMatrix(h, j, j+2, 0, 2);

            CameraExtrinsicParameters extr = solveForExtrinsic(kIntr, hI);

            double[] om = Rotation.extractRotationVectorRodriguesBouguet(extr.getRotation()).rotVec;

            extrinsics.add(new CameraExtrinsicParameters(extr.getRotation(), om,
                    extr.getTranslation()));

            // project the world coordinates into image frame, needed to estimate radial distortion
            double[][] x = MatrixUtil.multiply(hI, coordsW);
            System.arraycopy(x[0], 0, u, i*n, n);
            System.arraycopy(x[1], 0, v, i*n, n);
        }

        CameraMatrices cameraMatrices = new CameraMatrices();
        cameraMatrices.setIntrinsics(kIntr);
        cameraMatrices.getExtrinsics().addAll(extrinsics);

        double[] kRadial = solveForRadialDistortion(coordsI, u, v, cameraMatrices, useR2R4);

        kIntr.setRadialDistortionCoeffs(kRadial);
        kIntr.setUseR2R4(useR2R4);

        // (5) optimization to improve the parameter estimates

        return cameraMatrices;
    }

    /**
     * estimate the intrinsic and extrinsic camera parameters.
     To refine the returned values, follow with PNP.solveForPose(...).
     * The method uses double[] rVec = Rotation.extractRotationVectorRodrigues(r);
     * to extract the included rotation vectors.
     * TODO: this method could be overloaded to use an indicator variable for a feature presence in
     * an image, etc.
     * @param n n is the number of points in each image which is the
    same for all images.  n is number of features.
     * @param coordsI  holds the image coordinates in pixels of
    features present in all images ordered in the same
    manner and paired with features in coordsW.
    It is a 2-dimensional double array of format
    3 X (N*n) where N is the number of images.
    the first row is the x coordinates, the second row
    is the y coordinates, and the third row is "1"'s.
    The columns hold each image in order and within each image's
    columns are the features presented in the same order in each image.
    In Table 1 of Ma, Chen, & Moore 2003 "Camera Calibration"
    these are the (u_d, v_d) pairs.
     * @param coordsW holds the world coordinates of features, ordered
    by the same features in the images.
    the first row is the X coordinates, the second row
    is the Y coordinates, and the third row is 1's
    (Z_w = 0, the scale factor is lost in the homography).
    It is a 2-dimensional double array of format
    3 X n
     @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
     f(r) = 1 +k1*r^2 + k2*r^4 if true,
     else use model #3 f(r) = 1 +k1*r + k2*r^2.
      * @throws Exception if there is an error in use of MPSolver during the
     * removal of radial distortion, a generic exception is thrown with the
     * error message from the MPSolve documentation.
     * @return camera intrinsic parameters, extrinsic parameters, and radial
     * distortion coefficients
     */
    public static CameraMatrices estimateCamera3D(final int n, double[][] coordsI,
            double[][] coordsW, boolean useR2R4) throws NotConvergedException,
            Exception {

        //TODO: add refinement stages

        if (coordsI.length != 3) {
            throw new IllegalArgumentException("coordsI must have 3 rows.");
        }
        int nImages = coordsI[0].length/n;
        if (coordsI[0].length != nImages*n) {
            throw new IllegalArgumentException("coordsI must have nImages * n features as the number of columns");
        }
        if (coordsW.length != 3) {
            throw new IllegalArgumentException("coordsW must have 3 rows.");
        }

        boolean passive = true;

        int i;

        List<CameraExtrinsicParameters> extrinsics = new ArrayList<>();
        // average the intrinsic matrices.
        double[][] k = new double[3][3];
        for (i = 0; i < nImages; ++i) {
            double[][] x = MatrixUtil.copySubMatrix(coordsI, 0, 2, i*n, (i+1)*n-1);
            double[][] p = CameraPose.calculatePFromXXW(x, coordsW);
            Camera.CameraPoseParameters pose = CameraPose.calculatePoseFromP(p);

            k = MatrixUtil.pointwiseAdd(k, pose.getIntrinsicParameters().getIntrinsic());

            double[] om = Rotation.extractRotationVectorRodriguesBouguet(pose.getExtrinsicParameters().getRotation()).rotVec;

            extrinsics.add(new CameraExtrinsicParameters(pose.getExtrinsicParameters().getRotation(), om,
                    pose.getExtrinsicParameters().getTranslation()));
        }
        MatrixUtil.multiply(k, 1./nImages);
        CameraIntrinsicParameters kIntr = new CameraIntrinsicParameters(k);

        CameraMatrices cameraMatrices = new CameraMatrices();
        cameraMatrices.setIntrinsics(kIntr);
        cameraMatrices.getExtrinsics().addAll(extrinsics);

        // (4) estimate the radial distortion coefficients
        //     NOTE: There are a couple of radial functions which are commonly
        //           used.
        //           (a) f_r = 1 + k_1*r^2 + k_2*r^4
        //               which is Eqn #4 of Table 2 of Ma et al. 2004.
        //               Ma et al. 2004 statistics for it were among the best
        //               so they use it in their algorithms.
        //           (b) f_r = 1 + k_1*r + k_2*r^2
        //               which is Eqn #3 of Table 2 of Ma et al. 2004.
        //               it's a lower order function so may be a better choice
        //               for some data.
        //               Ma, Chen, Moore 2003 prefer (b) to (a) because:
        //                  (I) Low order fitting, better for fixed-point implementation
        //                  (II) Explicit inverse function with no numerical iterations
        //                  (IV) Better accuracy than radial distortion model (a)
        //
        //    NOTE: The Ma et al. 2003 paper found that the non-linear optimization
        //          works just as well with initial estimates of 0 for the radial
        //          distortion coefficients, so one can exclude this step.
        //
        //    NOTE: Ma, Soatta, Kosecka, & Sastry (year 2012? 2004?) have
        //          specified f(r) in terms of a center of radial distortion
        //          which is not necessarily the image centering Section (3.3.3)

        // Using the estimated intrinsic and extrinsic parameters,
        //  we can get the ideal projected image points Along with the real
        //  observed image points, we can estimate the
        // two distortion coefficients (k1,k2)

        double[][] xW = addDimension1(coordsW);

        // calculating the projected world coordinates using eqn (17).
        //    the homography transforms the world reference coordinates to the
        //    image reference frame (which are w.r.t. the corner of the image,
        //    not the center).
        double[] u = new double[n*nImages];
        double[] v = new double[n*nImages];
        for (i = 0; i < nImages; ++i) {
            // recalculating P here, but could instead save it from above and reuse here
            // P = K * [R | t]
            CameraExtrinsicParameters extr = cameraMatrices.getExtrinsics().get(i);
            double[][] r = extr.getRotation();
            double[] t = extr.getTranslation();
            double[][] p = new double[][] {
                    {r[0][0], r[0][1], r[0][2], t[0]},
                    {r[1][0], r[1][1], r[1][2], t[1]},
                    {r[2][0], r[2][1], r[2][2], t[2]}
            };
            p = MatrixUtil.multiply(kIntr.getIntrinsic(), p);
            double[][] x = MatrixUtil.multiply(p, xW);
            System.arraycopy(x[0], 0, u, i*n, n);
            System.arraycopy(x[1], 0, v, i*n, n);
        }

        double[] kRadial = solveForRadialDistortion(coordsI, u, v, cameraMatrices, useR2R4);

        kIntr.setRadialDistortionCoeffs(kRadial);
        kIntr.setUseR2R4(useR2R4);

        // (5) optimization to improve the parameter estimates

        // ============ iterate over the above steps after non-linear optimization
        //  for extrinsic parameters.
        double[][] cI;
        CameraExtrinsicParameters kExtr;
        CameraExtrinsicParameters extrinsic;
        int nMaxIter = 100;
        /*
        for (i = 0; i < nImages; ++i) {

            cI = MatrixUtil.copySubMatrix(coordsI, 0, 2, n*i, n*(i + 1)-1);

            kExtr = cameraMatrices.getExtrinsics().get(i);

            // improve the extrinsic parameter estimates:
            extrinsic = PNP.solveForPose(cI, coordsW, kIntr,
                kExtr, nMaxIter);

            cameraMatrices.getExtrinsics().set(i, extrinsic);
        }*/

        return cameraMatrices;
    }

    /**
     * given a matrix of length n, return same matrix with a row of 1s appended to it
     * @param a
     * @return
     */
    public static double[][] addDimension1(double[][] a) {
        int nR = a.length;
        double[][] out = new double[nR+1][a[0].length];
        for (int j = 0; j < nR; ++j) {
            System.arraycopy(a[j], 0, out[j], 0, a[j].length);
        }
        Arrays.fill(out[nR], 1);
        return out;
    }

    /**
     * for a given set of feature coordinates in image reference frame and in
     * world coordinate system, calculates the homography following the 
     * algorithm in Ma et al. 2003.
     * scaleFactor * [u v 1] (col) = K * [r1  r2  t] (col) * [X_W  Y_W 1] (col).
     * this is also in sect 2.4 of camera calibration book by zhang
     * as homography between a model plane and its image (in camera calibration w/ 2d objects: plane based techniques).
     *
     * see also CameraPose.calculatePFromXXW(...)
     *
     * @param coordsI holds the image coordinates in pixels of features present in image i as format [3 X nPoints].
     *                Only the first 2 dimensions are used, so if the 3rd dimension (z axis) is present, it is
     *                the responsibility of the invoker to have normalized the 1st 2 dimensions by the 3rd.
     * @param coordsW holds the world coordinates of features present in image 1 corresponding
               to the same features and order of coordsI_i as format [3 X nPoints].
     *                Only the first 2 dimensions are used, so if the 3rd dimension (z axis) is present, it is
     *                the responsibility of the invoker to have normalized the 1st 2 dimensions by the 3rd.
     * @return the homography, projection matrix
     */
    public static double[][] solveForHomography(double[][] coordsI, double[][] coordsW, boolean useNormConditioning)
            throws NotConvergedException {
        
        if (coordsI.length < 2 || coordsI.length > 3) {
            throw new IllegalArgumentException("coordsI must have 2 or 3 rows.");
        }
        if (coordsW.length < 2 || coordsW.length > 3) {
            throw new IllegalArgumentException("coordsW must have 2 or 3 rows.");
        }
        int n = coordsI[0].length;
        if (coordsW[0].length != n) {
            throw new IllegalArgumentException("coordsW must have same number of columns as coordsI.");
        }
        
        /*
        creates matrix L, and finds the solution to x as orthogonal to L by using the SVD(L)
           to find the eigenvector belonging to the smallest eigenvalue.
          -reformats x into 3x3 H to return
        */
        
        // Section 6.1 of Ma et al. 2003
        
        /*
          H =   [ h11 h12 h13 ]
                [ h21 h22 h23 ]
                [ h31 h32 h33 ]
          H^T = [ h11 h21 h31 ]
                [ h12 h22 h32 ]
                [ h13 h23 h33 ]

          Let h_i be the ith row of H:
              h_i = [h_i_1]^T = [h_i_1  h_i_2  h_i_3]
                    [h_i_2]
                    [h_i_3]
    */

        double[][] tI = null;
        double[][] tW = null;
        if (useNormConditioning) {
            coordsI = MatrixUtil.copy(coordsI);
            coordsW = MatrixUtil.copy(coordsW);
            tI = EpipolarNormalizationHelper.unitStandardNormalize(coordsI);
            tW = EpipolarNormalizationHelper.unitStandardNormalize(coordsW);
        }

        // 2*n X 9
        double u, v, X, Y;
        double[][] ell = new double[2*n][9];
        for (int i = 0; i < n; ++i) {
            u = coordsI[0][i];
            v = coordsI[1][i];
            X = coordsW[0][i];
            Y = coordsW[1][i];
            // eqn(15) of Ma et al. 2003
            ell[2*i]     = new double[]{X, Y, 1, 0, 0, 0, -u*X, -u*Y, -u};
            ell[2*i + 1] = new double[]{0, 0, 0, X, Y, 1, -v*X, -v*Y, -v};

            //NOTE: Zhang chap in Camera Calibration book uses:
            //ell[2*i]     = new double[]{X, Y, 1, 0, 0, 0, u*X, u*Y, u};
            //ell[2*i + 1] = new double[]{0, 0, 0, X, Y, 1, v*X, v*Y, v};
        }
        coordsI = null;
        coordsW = null;

        MatrixUtil.SVDProducts svd = MatrixUtil.performSVD(ell);

        // vT is 9X9.  last row in vT is the eigenvector for the smallest eigenvalue
        double[] xOrth = svd.vT[svd.vT.length - 1];

        // subject to ||x|| = 1
        xOrth = MatrixUtil.normalizeLP(xOrth, 2);
        if (xOrth[xOrth.length - 1] < 0) {
            MatrixUtil.multiply(xOrth, -1);
        }

        //h00, h01, h02, h10, h11, h12, h20,h21,h22

        double[][] h = new double[3][3];
        for (int i = 0; i < 3; i++) {
            System.arraycopy(xOrth, (i * 3), h[i], 0, 3);
        }

        svd = null;

        if (useNormConditioning) {
            h = MatrixUtil.multiply(EpipolarNormalizationHelper.inverseT(tI), h);
            h = MatrixUtil.multiply(h, tW);
        }

        return h;
    }

    /**
     * for a given set of feature coordinates in image reference frame and in
     * world coordinate system, calculates the homography following the
     * algorithm in Wetzstein "EE 267 Virtual Reality
     * Course Notes: 6-DOF Pose Tracking with the VRduino".
     *
     * This algorithm uses a passive right-hand transformation system.
     * It flips the z-coordinate system so that observer is in the origin looking down the negative z axis towards
     * the moving object.  The algorithm also sets h[2][2] to 1 in solving for pose.
     *
     * @param coordsC holds the feature coordinates in camera reference frame image i as format [3 X nPoints].
     *                Only the first 2 dimensions are used, so if the 3rd dimension (z axis) is present, it is
     *                the responsibility of the invoker to have normalized the 1st 2 dimensions by the 3rd.
     *                Note that the method should work similarly if input is feature coordinates in image frame instead,
     *                but the -z should be considered afterward when using the homography.
     * @param coordsW holds the world coordinates of features present in image 1 corresponding
    to the same features and order of coordsC_i as format [3 X nPoints].
     *                Only the first 2 dimensions are used, so if the 3rd dimension (z axis) is present, it is
     *                the responsibility of the invoker to have normalized the 1st 2 dimensions by the 3rd.
     * @return the homography, projection matrix
     */
    public static double[][] solveFor8PointHomography(double[][] coordsC, double[][] coordsW, boolean useNormConditioning)
            throws NotConvergedException {

        if (coordsC.length < 2 || coordsC.length > 3) {
            throw new IllegalArgumentException("coordsC must have 2 or 3 rows.");
        }
        if (coordsW.length < 2 || coordsW.length > 3) {
            throw new IllegalArgumentException("coordsW must have 2 or 3 rows.");
        }
        int n = coordsC[0].length;
        if (coordsW[0].length != n) {
            throw new IllegalArgumentException("coordsW must have same number of columns as coordsC.");
        }

        /*
          H =   [ h11 h12 h13 ]
                [ h21 h22 h23 ]
                [ h31 h32 h33 ]
          H^T = [ h11 h21 h31 ]
                [ h12 h22 h32 ]
                [ h13 h23 h33 ]
        */

        double[][] tC = null;
        double[][] tW = null;
        if (useNormConditioning) {
            coordsC = MatrixUtil.copy(coordsC);
            coordsW = MatrixUtil.copy(coordsW);
            tC = EpipolarNormalizationHelper.unitStandardNormalize(coordsC);
            tW = EpipolarNormalizationHelper.unitStandardNormalize(coordsW);
        }

        // 2*n X 9
        double u, v, X, Y;
        double[][] ell = new double[2*n][8];
        for (int i = 0; i < n; ++i) {
            u = coordsC[0][i];
            v = coordsC[1][i];
            X = coordsW[0][i];
            Y = coordsW[1][i];
            ell[2*i]     = new double[]{X, Y, 1, 0, 0, 0, -u*X, -u*Y};
            ell[2*i + 1] = new double[]{0, 0, 0, X, Y, 1, -v*X, -v*Y};
        }

        MatrixUtil.SVDProducts svd = MatrixUtil.performSVD(ell);

        // vT is 9X9.  last row in vT is the eigenvector for the smallest eigenvalue
        double[] xOrth = svd.vT[svd.vT.length - 1];

        // subject to ||x|| = 1
        xOrth = MatrixUtil.normalizeLP(xOrth, 2);
        if (xOrth[xOrth.length - 1] < 0) {
            MatrixUtil.multiply(xOrth, -1);
        }

        //h00, h01, h02, h10, h11, h12, h20,h21,h22

        double[][] h = new double[3][3];
        for (int i = 0; i < 3; i++) {
            System.arraycopy(xOrth, (i * 3), h[i], 0, 3);
        }
        h[2][2] = 1;

        if (useNormConditioning) {
            h = MatrixUtil.multiply(EpipolarNormalizationHelper.inverseT(tC), h);
            h = MatrixUtil.multiply(h, tW);
        }

        return h;
    }

    /**
     * solve for the planar homography between world coordinate objects xW and the imaged objects' coordinates.
     * x ~ H*xW where "~" means equal up to a non zero scalar factor.
     *
     * this method is similar to solveForHomography(), but has additional normalization and a refinement step,
     * all ported from github repositories holding the Bouguet Matlab Toolbox code.
     * The Bouguet's toolbox web page implies that the source is freely available.
     * The github repositories with the Bouguet Matlab code that the individual authors have modified do not have license
     * information.  Those references are here and the method this is adapted from.
     * <pre>
     *     http://robots.stanford.edu/cs223b04/JeanYvesCalib/
     *     compute_homography.m in
     *     https://github.com/fragofer/TOOLBOX_calib
     *     and
     *     https://github.com/hunt0r/Bouguet_cam_cal_toolbox
     * </pre>
     *
     * @param x0 image coordinates
     * @param xW0 world coordinates
     * @return
     * @throws NotConvergedException
     */
    public static double[][] solveForHomographyBouget(double[][] x0, double[][] xW0) throws NotConvergedException {

        if (x0.length != 3) {
            throw new IllegalArgumentException("x0 must have 3 rows.");
        }
        if (xW0.length != 3) {
            throw new IllegalArgumentException("xW0 must have 3 rows.");
        }
        int n = x0[0].length;
        if (xW0[0].length != n) {
            throw new IllegalArgumentException("xW0 must have same number of columns as x0.");
        }

        /*
         First computes an initial guess for the homography through quasi-linear method.
         Then, if the total number of points is larger than 4, optimize the solution by minimizing
         the reprojection error (in the least squares sense)
         */

        int i, j;

        double[][] x = MatrixUtil.copy(x0);
        double[][] xW = MatrixUtil.copy(xW0);

        //m = m ./ (ones(3,1)*m(3,:));
        //M = M ./ (ones(3,1)*M(3,:));
        for (i = 0; i < x[0].length; ++i) {
            for (j = 0; j < x.length; ++j) {
                x[j][i] /= x[x.length - 1][i];
            }
        }
        for (i = 0; i < xW[0].length; ++i) {
            for (j = 0; j < xW.length; ++j) {
                xW[j][i] /= xW[xW.length - 1][i];
            }
        }

        //Prenormalization of point coordinates (very important):
        // (Affine normalization)

        double[] ax = x[0];
        double[] ay = x[1];

        double mxx = MiscMath.getAvgAndStDev(ax)[0];
        double myy = MiscMath.getAvgAndStDev(ay)[0];
        ax = MatrixUtil.subtract(ax, mxx);
        ay = MatrixUtil.subtract(ay, myy);

        double scxx = meanOfAbs(ax);
        double scyy = meanOfAbs(ay);

        double[][] hNorm = new double[3][];
        hNorm[0] = new double[]{1./scxx, 0., -mxx/scxx};
        hNorm[1] = new double[]{0, 1./scyy, -myy/scyy};
        hNorm[2] = new double[]{0, 0, 1};

        double[][] invHNorm = new double[3][];
        invHNorm[0] = new double[]{scxx, 0, mxx};
        invHNorm[1] = new double[]{0, scyy, myy};
        invHNorm[2] = new double[]{0, 0, 1};

        //mn = Hnorm*m;
        //[3 X 3] * [3 X n]
        double[][] xn = MatrixUtil.multiply(hNorm, x);

        // 2*n X 9
        double u, v, X, Y, Z;
        double[][] ell = new double[2*n][9];
        for (i = 0; i < n; ++i) {
            u = xn[0][i];
            v = xn[1][i];
            X = xW[0][i];
            Y = xW[1][i];
            Z = xW[2][i];
            // last rows are "1" for both xn and xW
            // eqn(15) of Ma et al. 2003
            ell[2*i]     = new double[]{X, Y, Z, 0, 0, 0, -u*X, -u*Y, -u*Z};
            ell[2*i + 1] = new double[]{0, 0, 0, X, Y, Z, -v*X, -v*Y, -v*Z};
        }

        if (n > 4) {
            // SVD(A).V == SVD(A^T*A).V
            ell = MatrixUtil.createATransposedTimesA(ell);
        }

        SVD svd = SVD.factorize(new DenseMatrix(ell));
        double[][] vT = MatrixUtil.convertToRowMajor(svd.getVt());
        // vT is 9X9.  last row in vT is the eigenvector for the smallest eigenvalue
        double[] xOrth = vT[vT.length - 1];
        //h00, h01, h02, h10, h11, h12,h20,h21,h22
        MatrixUtil.multiply(xOrth, 1./xOrth[xOrth.length - 1]);
        //TODO: consider normalization instead:
        //xOrth = MatrixUtil.normalizeLP(xOrth, 2);
        //if (xOrth[xOrth.length - 1] < 0) {
        //    MatrixUtil.multiply(xOrth, -1);
        //}

        // Hrem = reshape(hh,3,3)';
        // Matlab reshape fills along columns, but this is transposed, so fill rows
        double[][] h = new double[3][3];
        h[0] = Arrays.copyOfRange(xOrth, 0, 3);
        h[1] = Arrays.copyOfRange(xOrth, 3, 6);
        h[2] = Arrays.copyOfRange(xOrth, 6, 9);

        h = MatrixUtil.multiply(invHNorm, h);

        if (true) {
            // a quick look at errors
            double[][] xEst = MatrixUtil.multiply(h, xW);
            for (i = 0; i < xEst[0].length; ++i) {
                for (j = 0; j < xEst.length; ++j) {
                    xEst[j][i] /= xEst[xEst.length - 1][i];
                }
            }
            double[][] err = MatrixUtil.pointwiseSubtract(x, xEst);
            err = MatrixUtil.copySubMatrix(err, 0, 1, 0, err[0].length - 1);
            double[] xMeanStdv = MiscMath0.getAvgAndStDev(err[0]);
            double[] yMeanStdv = MiscMath0.getAvgAndStDev(err[1]);
            System.out.printf("x err=%s\n", FormatArray.toString(xMeanStdv, "%.4e"));
            System.out.printf("y err=%s\n", FormatArray.toString(yMeanStdv, "%.4e"));

            double[][] xEst2 = MatrixUtil.multiply(h, xW0);
            for (i = 0; i < xEst2[0].length; ++i) {
                for (j = 0; j < xEst2.length; ++j) {
                    xEst2[j][i] /= xEst2[xEst.length - 1][i];
                }
            }
            double[][] err2 = MatrixUtil.pointwiseSubtract(x0, xEst2);
            err2 = MatrixUtil.copySubMatrix(err2, 0, 1, 0, err2[0].length - 1);
            double[] xMeanStdv2 = MiscMath0.getAvgAndStDev(err2[0]);
            double[] yMeanStdv2 = MiscMath0.getAvgAndStDev(err2[1]);
            System.out.printf("x0 err=%s\n", FormatArray.toString(xMeanStdv2, "%.4e"));
            System.out.printf("y0 err=%s\n", FormatArray.toString(yMeanStdv2, "%.4e"));
        }

        if (n <= 4) {
            return h;
        }

        // refinement to improve solution
        //hhv = reshape(H',9,1);
        double[] hhv = MatrixUtil.stack(MatrixUtil.transpose(h));

        //hhv = hhv(1:8);
        hhv = Arrays.copyOf(hhv, 8);

        for (int iter = 1; iter <= 10; ++iter) {

            //mrep = H * M;
            double[][] mrep = MatrixUtil.multiply(h, xW);

            double[][] J = MatrixUtil.zeros(2*n,8);

            //MMM = (M ./ (ones(3,1)*mrep(3,:)));   // ./ = elementwise division
            double[][] MMM = new double[3][n];
            for (i = 0; i < 3; ++i) {
                MMM[i] = MatrixUtil.pointwiseDivision(xW[i], mrep[2]);
            }

            //J(1:2:2*Np,1:3) = -MMM';
            //J(2:2:2*Np,4:6) = -MMM';
            for (j = 0; j < n; ++j) {
                //J[2*j] = new double[]{-MMM[0][j], -MMM[1][j], -MMM[2][j],   0,0,0,    0,0};
                System.arraycopy(new double[]{-MMM[0][j], -MMM[1][j], -MMM[2][j]},
                        0, J[2*j], 0, 3);
                //J[2*j + 1] = new double[]{0, 0, 0, -MMM[0][j], -MMM[1][j], -MMM[2][j], 0,0};
                System.arraycopy(new double[]{-MMM[0][j], -MMM[1][j], -MMM[2][j]},
                        0, J[2*j], 3, 3);
            }

            //mrep = mrep ./ (ones(3,1)*mrep(3,:));
            for (i = 0; i < 3; ++i) {
                mrep[i] = MatrixUtil.pointwiseDivision(mrep[i], mrep[2]);
            }

            // [2 * n]
            //m_err = m(1:2,:) - mrep(1:2,:);
            //m_err = m_err(:); // <=== this stacks each column after the previous
            double[] merr = new double[2*n];
            int c = 0;
            for (i = 0; i < n; ++i) {
                for (j = 0; j < 2; ++j) {
                    merr[c] = xn[j][i] - mrep[j][i];
                    c++;
                }
            }

            //MMM2 = (ones(3,1)*mrep(1,:)) .* MMM;
            //MMM3 = (ones(3,1)*mrep(2,:)) .* MMM;
            //     .* is elementwise multiplication
            double[][] MMM2 = new double[3][]; // [3 X n]
            for (i = 0; i < 3; ++i) {
                MMM2[i] = MatrixUtil.pointwiseMultiplication(mrep[0], MMM[i]);
            }
            double[][] MMM3 = new double[3][]; // [3 X n]
            for (i = 0; i < 3; ++i) {
                MMM3[i] = MatrixUtil.pointwiseMultiplication(mrep[1], MMM[i]);
            }

            //J(1:2:2*Np,7:8) = MMM2(1:2,:)';
            //J(2:2:2*Np,7:8) = MMM3(1:2,:)';
            for (j = 0; j < n; ++j) {
                J[2*j][6] = MMM2[0][j];
                J[2*j][7] = MMM2[1][j];
                J[2*j + 1][6] = MMM3[0][j];
                J[2*j + 1][7] = MMM3[1][j];
            }

            //MMM = (M ./ (ones(3,1)*mrep(3,:)))';
            for (i = 0; i < 3; ++i) {
                MMM[i] = MatrixUtil.pointwiseDivision(xW[i], mrep[2]);
            }
            // [n X 3]
            MMM = MatrixUtil.transpose(MMM);

            //                ([8 X 2*n]*[2*n X 8])= [8 X 8];
            //                                                J'*m_err = [8 X 2*n] * [2 * n X 1] = [8 X 1]
            //hh_innov  = inv(J'*J)*J'*m_err; // [8X1]
            double[][] t1 = MatrixUtil.pseudoinverseFullColumnRank(MatrixUtil.createATransposedTimesA(J));
            double[] t2 = MatrixUtil.multiplyMatrixByColumnVector(MatrixUtil.transpose(J), merr);
            double[] hhInov = MatrixUtil.multiplyMatrixByColumnVector(t1, t2);

            // length 8:
            //hhv_up = hhv - hh_innov;
            double[] hhvUp = new double[8];
            MatrixUtil.pointwiseSubtract(hhv, hhInov, hhvUp);

            //H_up = reshape([hhv_up;1],3,3)';
            // reshape writes into columns, but the transpose means write into rows
            double[][] hUp = new double[3][];
            hUp[0] = Arrays.copyOfRange(hhvUp, 0, 3);
            hUp[1] = Arrays.copyOfRange(hhvUp, 3, 6);
            hUp[2] = new double[]{hhvUp[6], hhvUp[7], 1};

            // %norm(m_err)
            // %norm(hh_innov)

            //hhv = hhv_up;
            hhv = hhvUp;

            //H = H_up;
            h = hUp;
        }

        if (true) {
            // a quick look at errors
            double[][] xEst = MatrixUtil.multiply(h, xW);
            for (i = 0; i < xEst[0].length; ++i) {
                for (j = 0; j < xEst.length; ++j) {
                    xEst[j][i] /= xEst[xEst.length - 1][i];
                }
            }
            double[][] err = MatrixUtil.pointwiseSubtract(x, xEst);
            err = MatrixUtil.copySubMatrix(err, 0, 1, 0, err[0].length - 1);
            double[] xMeanStdv = MiscMath0.getAvgAndStDev(err[0]);
            double[] yMeanStdv = MiscMath0.getAvgAndStDev(err[1]);
            System.out.printf("x err=%s\n", FormatArray.toString(xMeanStdv, "%.4e"));
            System.out.printf("y err=%s\n", FormatArray.toString(yMeanStdv, "%.4e"));
        }

        return h;
    }

    private static double meanOfAbs(double[] a) {
        double sum = 0;
        for (int i = 0; i < a.length; ++i) {
            sum += Math.abs(a[i]);
        }
        return sum/a.length;
    }

    /**
     * estimate the camera intrinsic parameters from the image homographies.
     * @param h H as (3*NImages)x3 homography, projection matrices
              where each image homography is stacked row-wise
     * @return the camera intrinsic parameters.
     */
    static CameraIntrinsicParameters solveForIntrinsicPlanar(double[][] h) throws NotConvergedException {
        
        if (h[0].length != 3) {
            throw new IllegalArgumentException("h must have 3 columns");
        }
        int nImages = h.length/3;
        if (nImages < 3) {
            throw new IllegalArgumentException("we need at least 3 images for the planar solution");
        }

        log.log(LEVEL, "h=\n%s\n", FormatArray.toString(h, "%.3e"));

        // Section 6.3 of Ma et al. 2003
        
        /*
          H =   [ h11 h12 h13 ]
                [ h21 h22 h23 ]
                [ h31 h32 h33 ]
          H^T = [ h11 h21 h31 ]
                [ h12 h22 h32 ]
                [ h13 h23 h33 ]

        Let h_i be the ith column vector of H:
              h_i = [h_i_1]^T = [h_i_1  h_i_2  h_i_3]
                    [h_i_2]
                    [h_i_3]
        
        - for each H:
                    form a matrix V_i_j out of the first 2 columns of each H matrix
                    and stack them by rows, into a matrix called V
              - perform SVD(V) to get right singular vector of V associated with the smallest singular value
                as the solution to b.
              - b holds the contents of the upper right triangle of B
                where B = A^-T * A^-1 known as the absolute conic.
              - the intrinsic parameters are extracted from combinations of the solved
                for B and other coefficients.
        
        b = [B11, B12, B22, B13, B23, B33]^T
        */

        int n = h.length/3;

        // 2*nImages X 6
        double[][] v = new double[2*n][6];
        double[] v22 = new double[6];
        double h11, h12, h13, h21, h22, h23;
        //Vij = [hi1*hj1, hi1*hj2 + hi2*hj1, hi2*hj2, hi3*hj1 + hi1*hj3, hi3*hj2 + hi2*hj3, hi3*hj3]T
        for (int i = 0; i < n; ++i) {
            // h_i is the ith column vector of H
            // h11 = column 0 of h, first element: h[0][0]
            // h12 = column 0 of h, 2nd element:   h[1][0]
            // h13 = column 0 of h, 3rd element:   h[2][0]
            // h21 = column 1 of h, first element: h[0][1]
            // h22 = column 1 of h, 2nd element:   h[1][1]
            // h23 = column 1 of h, 3rd element:   h[2][1]
            h11 = h[0+3*i][0]; h12 = h[1+3*i][0]; h13 = h[2+3*i][0];
            h21 = h[0+3*i][1]; h22 = h[1+3*i][1]; h23 = h[2+3*i][1];

            //h11 = h[0+3*i][0]; h12 = h[0+3*i][1]; h13 = h[0+3*i][2];
            //h21 = h[1+3*i][0]; h22 = h[1+3*i][1]; h23 = h[1+3*i][2];

            //V12 = [h11*h21, h11*h22 + h12*h21, h12*h22, h13*h21 + h11*h23,
            //       h13*h22 + h12*h23, h13*h23]T
            v[2*i] = new double[]{
                    h11*h21, h11*h22 + h12*h21, h12*h22, h13*h21 + h11*h23,
                    h13*h22 + h12*h23, h13*h23
            };

            //V11 = [
            // h11*h11 ,
            // h11*h12 + h12*h11,
            // h12*h12,
            // h13*h11 + h11*h13,
            // h13*h12 + h12*h13,
            // h13*h13]T
            // 2nd row is V11 - V12
            /*v[2*i + 1] = new double[]{
                h11*h11           - v[2*i][0],
                h11*h12 + h12*h11 - v[2*i][1],
                h12*h12           - v[2*i][2],
                h13*h11 + h11*h13 - v[2*i][3],
                h13*h12 + h12*h13 - v[2*i][4],
                h13*h13           - v[2*i][5]
            };*/

            //Vij = [hi1*hj1, hi1*hj2 + hi2*hj1, hi2*hj2, hi3*hj1 + hi1*hj3, hi3*hj2 + hi2*hj3, hi3*hj3]T
            v22[0] = h21*h21;
            v22[1] = h21*h22 + h22*h21;
            v22[2] = h22*h22;
            v22[3] = h23*h21 + h21*h23;
            v22[4] = h23*h22 + h22*h23;
            v22[5] = h23*h23;

            // v11 - v22; Zhang 99 eqn(8)
            v[2*i + 1] = new double[]{
                    h11*h11           - v22[0],
                    h11*h12 + h12*h11 - v22[1],
                    h12*h12           - v22[2],
                    h13*h11 + h11*h13 - v22[3],
                    h13*h12 + h12*h13 - v22[4],
                    h13*h13           - v22[5]
            };

        }

        //Vb = 0 and b = [B11, B12, B22, B13, B23, B33]^T
        SVDProducts svd = MatrixUtil.performSVD(MatrixUtil.createATransposedTimesA(v));

        if (svd.rank < 6) {
            System.out.printf("warning, rank < 6 for the design matrix of planar scene and image points\n");
        }
        // vT is 6X6.  last row in vT is the eigenvector for the smallest eigenvalue for full rank
        // it's the solution for the right null space.   there are problems when rank < 6
        double[] b = svd.vT[5];

        double[][] kIntr = null;
        double[][] B = new double[][]{{b[0], b[1], b[3]}, {b[1], b[2], b[4]}, {b[3], b[4], b[5]}};
        if (true) {
            boolean isPD = MatrixUtil.isPositiveDefinite(B);
            B = MatrixUtil.nearestPositiveSemidefiniteToA(B, 1E-9);
            DenseCholesky chol = new DenseCholesky(B.length, false);
            chol = chol.factor(new LowerSPDDenseMatrix(new DenseMatrix(B)));
            LowerTriangDenseMatrix _cholL = chol.getL();
            double[][] cholL = Matrices.getArray(_cholL);
            double[][] cholLT = MatrixUtil.transpose(cholL);
            kIntr = MatrixUtil.inverse(cholLT);
            MatrixUtil.multiply(kIntr, cholL[2][2]);
        } else {

            //       0    1    2    3    4    5
            //b = [B11, B12, B22, B13, B23, B33]^T
            log.log(LEVEL, String.format("b=%s\n", FormatArray.toString(b, "%.3e")));
            double B11 = b[0];
            double B12 = b[1];
            double B22 = b[2];
            double B13 = b[3];
            double B23 = b[4];
            double B33 = b[5];
            //Zhang 99 Appendix B; Ma et al. 2003 eqn (26)
            double v0 = (B12 * B13 - B11 * B23) / (B11 * B22 - B12 * B12);
            double lambda = B33 - ((B13 * B13 + v0 * (B12 * B13 - B11 * B23)) / B11);
            double alpha = Math.sqrt(lambda / B11);
            double beta = Math.sqrt(lambda * B11 / (B11 * B22 - B12 * B12));
            double gamma = -B12 * alpha * alpha * beta / lambda;
            double u0 = (gamma * v0 / beta) - (B13 * alpha * alpha / lambda);
            //u0 = (gamma*v0/alpha) - (B13*alpha*alpha/lambda);

            log.log(LEVEL, String.format("v0=%.4e (exp=220.866)\n", v0));
            log.log(LEVEL, String.format("lambda=%.4e\n", lambda));
            log.log(LEVEL, String.format("alpha=%.4e\n", alpha));
            log.log(LEVEL, String.format("beta=%.4e\n", beta));
            log.log(LEVEL, String.format("gamma=%.4e\n", gamma));
            log.log(LEVEL, String.format("u0=%.4e\n", u0));

            kIntr = Camera.createIntrinsicCameraMatrix(
                    alpha, beta, u0, v0, gamma);
        }

        // enforce element [2][2]=1
        MatrixUtil.multiply(kIntr, 1./kIntr[2][2]);

        //B = lambda * K^-T * K  [3X3] = [3X3]*[3X3]
        //K^T*K^-1*B = lambda
        double[][] kTKInvB = MatrixUtil.multiply(MatrixUtil.multiply(MatrixUtil.transpose(kIntr),
                MatrixUtil.inverse(kIntr)), B);

        double lambda1_1 = 1./MatrixUtil.lPSum(MatrixUtil.extractColumn(kTKInvB, 0), 2);
        double lambda1_2 = 1./MatrixUtil.lPSum(MatrixUtil.extractColumn(kTKInvB, 1), 2);
        // or mean or median or column 0 and column 1

        CameraIntrinsicParameters intrinsics = new CameraIntrinsicParameters();
        intrinsics.setIntrinsic(kIntr);
        intrinsics.setLambda1(lambda1_1);
        intrinsics.setLambda2(lambda1_2);

        return intrinsics;
    }

    /**
     * estimate the camera intrinsic parameters from the image homographies.
     * @param h H as (3*NImages)x3 homography, projection matrices
    where each image homography is stacked row-wise
     * @return the camera intrinsic parameters.
     */
    static CameraIntrinsicParameters _solveForIntrinsicPlanar(double[][] h) throws NotConvergedException {

        if (h[0].length != 3) {
            throw new IllegalArgumentException("h must have 3 columns");
        }
        int nImages = h.length/3;
        if (nImages < 3) {
            throw new IllegalArgumentException("we need at least 3 images for the planar solution");
        }

        log.log(LEVEL, "h=\n%s\n", FormatArray.toString(h, "%.3e"));

        // Section 6.3 of Ma et al. 2003

        /*
          H =   [ h11 h12 h13 ]
                [ h21 h22 h23 ]
                [ h31 h32 h33 ]
          H^T = [ h11 h21 h31 ]
                [ h12 h22 h32 ]
                [ h13 h23 h33 ]

        Let h_i be the ith column vector of H:
              h_i = [h_i_1]^T = [h_i_1  h_i_2  h_i_3]
                    [h_i_2]
                    [h_i_3]

        - for each H:
                    form a matrix V_i_j out of the first 2 columns of each H matrix
                    and stack them by rows, into a matrix called V
              - perform SVD(V) to get right singular vector of V associated with the smallest singular value
                as the solution to b.
              - b holds the contents of the upper right triangle of B
                where B = A^-T * A^-1 known as the absolute conic.
              - the intrinsic parameters are extracted from combinations of the solved
                for B and other coefficients.

        b = [B11, B12, B22, B13, B23, B33]^T
        */

        int n = h.length/3;

        // 2*nImages X 6
        double[][] v = new double[2*n][6];
        double[] v22 = new double[6];
        double h11, h12, h13, h21, h22, h23;
        //Vij = [hi1*hj1, hi1*hj2 + hi2*hj1, hi2*hj2, hi3*hj1 + hi1*hj3, hi3*hj2 + hi2*hj3, hi3*hj3]T
        for (int i = 0; i < n; ++i) {
            // h_i is the ith column vector of H
            // h11 = column 0 of h, first element: h[0][0]
            // h12 = column 0 of h, 2nd element:   h[1][0]
            // h13 = column 0 of h, 3rd element:   h[2][0]
            // h21 = column 1 of h, first element: h[0][1]
            // h22 = column 1 of h, 2nd element:   h[1][1]
            // h23 = column 1 of h, 3rd element:   h[2][1]
            h11 = h[0+3*i][0]; h12 = h[1+3*i][0]; h13 = h[2+3*i][0];
            h21 = h[0+3*i][1]; h22 = h[1+3*i][1]; h23 = h[2+3*i][1];

            //h11 = h[0+3*i][0]; h12 = h[0+3*i][1]; h13 = h[0+3*i][2];
            //h21 = h[1+3*i][0]; h22 = h[1+3*i][1]; h23 = h[1+3*i][2];

            //V12 = [h11*h21, h11*h22 + h12*h21, h12*h22, h13*h21 + h11*h23,
            //       h13*h22 + h12*h23, h13*h23]T
            v[2*i] = new double[]{
                    h11*h21, h11*h22 + h12*h21, h12*h22, h13*h21 + h11*h23,
                    h13*h22 + h12*h23, h13*h23
            };

            //V11 = [
            // h11*h11 ,
            // h11*h12 + h12*h11,
            // h12*h12,
            // h13*h11 + h11*h13,
            // h13*h12 + h12*h13,
            // h13*h13]T
            // 2nd row is V11 - V12
            /*v[2*i + 1] = new double[]{
                h11*h11           - v[2*i][0],
                h11*h12 + h12*h11 - v[2*i][1],
                h12*h12           - v[2*i][2],
                h13*h11 + h11*h13 - v[2*i][3],
                h13*h12 + h12*h13 - v[2*i][4],
                h13*h13           - v[2*i][5]
            };*/

            //Vij = [hi1*hj1, hi1*hj2 + hi2*hj1, hi2*hj2, hi3*hj1 + hi1*hj3, hi3*hj2 + hi2*hj3, hi3*hj3]T
            v22[0] = h21*h21;
            v22[1] = h21*h22 + h22*h21;
            v22[2] = h22*h22;
            v22[3] = h23*h21 + h21*h23;
            v22[4] = h23*h22 + h22*h23;
            v22[5] = h23*h23;

            // v11 - v22; Zhang 99 eqn(8)
            v[2*i + 1] = new double[]{
                    h11*h11           - v22[0],
                    h11*h12 + h12*h11 - v22[1],
                    h12*h12           - v22[2],
                    h13*h11 + h11*h13 - v22[3],
                    h13*h12 + h12*h13 - v22[4],
                    h13*h13           - v22[5]
            };

        }

        //Vb = 0 and b = [B11, B12, B22, B13, B23, B33]^T
        SVDProducts svd = MatrixUtil.performSVD(v);

        if (svd.rank < 6) {
            System.out.printf("warning, rank < 6 for the design matrix of planar scene and image points\n");
        }
        // vT is 6X6.  last row in vT is the eigenvector for the smallest eigenvalue for full rank
        // it's the solution for the right null space.   there are problems when rank < 6
        double[] b = svd.vT[svd.rank-1];

        /*
        this isnt necessary, but could be used:
        // B is symmetric
            double[][] B = new double[][]{{b[0], b[1], b[3]}, {b[1], b[2], b[4]}, {b[3], b[4], b[5]}};

            boolean isPD = MatrixUtil.isPositiveDefinite(B);
            if (!isPD && !m1) {
                MatrixUtil.multiply(b, -1);
                m1 = true;
                isPD = MatrixUtil.isPositiveDefinite(B);
            }
            if (!isPD) {
                B = MatrixUtil.nearestPositiveSemidefiniteToA(B, 1E-9);
                int _t1 = 1;
            }

            // use cholesky decomposition
            // L = Solve(L · L^T = B)
            // A = (L^-1)^T * L[2][2]
            // where A is kIntr
            DenseCholesky chol = new DenseCholesky(B.length, false);
            chol = chol.factor(new LowerSPDDenseMatrix(new DenseMatrix(B)));
            LowerTriangDenseMatrix _cholL = chol.getL();
            double[][] cholL = Matrices.getArray(_cholL);
            double[][] cholLT = MatrixUtil.transpose(cholL);
            double[][] kIntr0 = MatrixUtil.inverse(cholLT);
            MatrixUtil.multiply(kIntr0, cholL[2][2]);

            //       0    1    2    3    4    5
            //b = [B11, B12, B22, B13, B23, B33]^T
            log.log(LEVEL, String.format("b=%s\n", FormatArray.toString(b, "%.3e")));
            double B11 = b[0];
            double B12 = b[1];
            double B22 = b[2];
            double B13 = b[3];
            double B23 = b[4];
            double B33 = b[5];
         */

        //       0    1    2    3    4    5
        //b = [B11, B12, B22, B13, B23, B33]^T
        log.log(LEVEL, String.format("b=%s\n", FormatArray.toString(b, "%.3e")));
        double B11 = b[0];
        double B12 = b[1];
        double B22 = b[2];
        double B13 = b[3];
        double B23 = b[4];
        double B33 = b[5];
        //Zhang 99 Appendix B; Ma et al. 2003 eqn (26)
        double v0 = (B12*B13 - B11*B23)/(B11*B22 - B12*B12);
        double lambda = B33 - ((B13*B13 + v0*(B12*B13 - B11*B23))/B11);
        double alpha = Math.sqrt(lambda/B11);
        double beta = Math.sqrt( lambda*B11 / (B11*B22 - B12*B12) );
        double gamma = -B12*alpha*alpha*beta / lambda;
        double u0 = (gamma*v0/beta) - (B13*alpha*alpha/lambda);
        //u0 = (gamma*v0/alpha) - (B13*alpha*alpha/lambda);

        log.log(LEVEL, String.format("v0=%.4e (exp=220.866)\n", v0));
        log.log(LEVEL, String.format("lambda=%.4e\n", lambda));
        log.log(LEVEL, String.format("alpha=%.4e\n", alpha));
        log.log(LEVEL, String.format("beta=%.4e\n", beta));
        log.log(LEVEL, String.format("gamma=%.4e\n", gamma));
        log.log(LEVEL, String.format("u0=%.4e\n", u0));

        double[][] kIntr = Camera.createIntrinsicCameraMatrix(
                alpha, beta, u0, v0, gamma);

        //B = lambda * K^-T * K  [3X3] = [3X3]*[3X3]
        //K^T*K^-1*B = lambda
        double[][] B = new double[][]{
                {B11, B12, B13},
                {B12, B22, B23},
                {B13, B23, B33}
        };
        double[][] kTKInvB = MatrixUtil.multiply(MatrixUtil.multiply(MatrixUtil.transpose(kIntr),
                MatrixUtil.inverse(kIntr)), B);

        double lambda1_1 = 1./MatrixUtil.lPSum(MatrixUtil.extractColumn(kTKInvB, 0), 2);
        double lambda1_2 = 1./MatrixUtil.lPSum(MatrixUtil.extractColumn(kTKInvB, 1), 2);
        // or mean or median or column 0 and column 1

        CameraIntrinsicParameters intrinsics = new CameraIntrinsicParameters();
        intrinsics.setIntrinsic(kIntr);
        intrinsics.setLambda1(lambda1_1);
        intrinsics.setLambda2(lambda1_2);

        return intrinsics;
    }

    /**
     * following Ma et al. 2003
     * estimate the extrinsic parameters from the image of the absolute conic.
     * @param kIntr camera intrinsic parameters
     * @param h homography for the projection for an image. at least 5 points should have been used
     *          to generate h.
     * @return
     */
    static Camera.CameraExtrinsicParameters solveForExtrinsic(
            CameraIntrinsicParameters kIntr, double[][] h) throws NotConvergedException {

        // notes from Serge Belongie lectures from Computer Vision II, CSE 252B, USSD
        // homogeneous repr of a point is x_vec = (x, y, 1)^T
        // equation f a line is ell = a*x + b*y + c = 0;
        // line rewritten in homogeneous coordinates is x_vec^T * ell.
        // general conic in 3 dimensions is a*x^2 + b*x*y + c*y^2 + d*x*z + e*y*z + f*z^2 = 0.
        //     rewritten using 2D homogeneous coordinates, quadratic form: x_vec^T * C * x_vec = 0
        //                  [a   b/2   d/2 ]
        //        where C = [b/2   c   c/2 ]
        //                  [d/2 c/2     f ]
        //        C has 6 variable, 5 DOF, so need 5 points
        //     can then reformat x_vec^T * C * x_vec = 0 into
        //        the "design matrix" * "the carrier vector" = 0
        //               A * c = 0
        //        c = SVD(A).V^T[n-1], the eigenvector assoc w/ smallest eigenvalue.
        //
        //        there are 3 cases for the smallest eigenvalue of SVD(A):
        //          (1) SVD(A).s[5] == 0, and n=5, then a conic exists that
        //              fits the data exactly
        //          (2) SVD(A).s[5] >=0, and n > 5, then the value is the goodness of fit
        //          (3) n < 5, the conic is undetermined and requires other means to solve.
        //

        // points at infinity, a.k.a. ideal points, have the form (x, y, 0)^T
        // the line at infinity is (0, 0, 1)^T.

        boolean passive = true;

        double[][] aInv = Camera.createIntrinsicCameraMatrixInverse(kIntr.getIntrinsic());// K^-1

        //h_i is the ith column vector of H
        //r1 = λ * A^−1 * h1
        //r2 = λ * A^−1 * h2
        //r3 = r1×r2
        // t = λ * A^−1 * h3

        double[] h1 = MatrixUtil.extractColumn(h, 0);
        double[] h2 = MatrixUtil.extractColumn(h, 1);
        double[] h3 = MatrixUtil.extractColumn(h, 2);

        double[] r1 = MatrixUtil.multiplyMatrixByColumnVector(aInv, h1);
        double[] r2 = MatrixUtil.multiplyMatrixByColumnVector(aInv, h2);
        double[] t = MatrixUtil.multiplyMatrixByColumnVector(aInv, h3);

        //λ = 1/||(A^−1)*(h1)||_2 = 1/||(A^−1)*(h2)||_2
        double lambda1_1 = 1./MatrixUtil.lPSum(r1, 2);
        double lambda1_2 = 1./MatrixUtil.lPSum(r2, 2);
        //double scaleFactor = 2./(Math.sqrt(sumOfSquares(h1)) + Math.sqrt(sumOfSquares(h1)));
        log.log(LEVEL, String.format("lambda1=%.3e, lambda2=%.3e\n", lambda1_1, lambda1_2));

        MatrixUtil.multiply(r1, lambda1_1);
        MatrixUtil.multiply(r2, lambda1_1);
        MatrixUtil.multiply(t, lambda1_1);

        double[] r3 = MatrixUtil.crossProduct(r1, r2);
        // r1, r2, and r3 are columns of R
        double[][] r = MatrixUtil.zeros(3, 3);
        for (int row = 0; row < 3; ++row) {
            r[row][0] = r1[row];
            r[row][1] = r2[row];
            r[row][2] = r3[row];
        }

        // orthonormalization of r:
        SVDProducts svd = MatrixUtil.performSVD(r);
        r = MatrixUtil.multiply(svd.u, svd.vT);

        double detR = MatrixUtil.determinant(r);
        detR = Math.round(detR*1E11)/1E11;
        assert(Math.abs(detR - 1) < 1E-7);

        CameraExtrinsicParameters kExtr = new Camera.CameraExtrinsicParameters();
        kExtr.setRotation(r);
        kExtr.setTranslation(t);
        kExtr.setRodriguesVector(Rotation.extractRotationVectorRodriguesBouguet(r).rotVec);

        return kExtr;
    }

    /**
     calculate the camera rotation and translation given camera coordinates and world coordinates
     of features.
     Note that coordsC = (intrinsicCamera)^-1 * coordsI where coordsI are the feature coordinates in
     the image frame in units of pixels.
     * following Wetzstein "EE 267 Virtual Reality
     *      * Course Notes: 6-DOF Pose Tracking with the VRduino"
     * estimate the extrinsic parameters from the features given positions in camera coordinates
     * and WCS real world coordinates.
     *
     * The pose results can be improved by following it with use of non-linear Levenberg-Marquardt or
     * other optimization method. see Appendix A of Wetzstein reference.
     *
     * Also, one can follow with radial distortion corrections.
     *
     * @param coordsC holds the feature coordinates in camera reference frame image i as format [3 X nPoints].
     *                Only the first 2 dimensions are used, so if the 3rd dimension (z axis) is present, it is
     *                the responsibility of the invoker to have normalized the 1st 2 dimensions by the 3rd.
     *                Note that the method should work similarly if input is feature coordinates in image frame instead,
     *                but the -z should be considered afterward when using the homography.
     * @param coordsW holds the world coordinates of features present in image 1 corresponding
    to the same features and order of coordsC_i as format [3 X nPoints].
     *                Only the first 2 dimensions are used, so if the 3rd dimension (z axis) is present, it is
     *                the responsibility of the invoker to have normalized the 1st 2 dimensions by the 3rd.
     * @return the pose, that is camera rotation and translation
     * @return 
     */
    static Camera.CameraExtrinsicParameters solveForExtrinsicPlanarWetzstein(
        double[][] coordsC, double[][] coordsW) throws NotConvergedException {

        boolean passive = true;

        boolean useNormConditioning = false;
        double[][] h = solveFor8PointHomography(coordsC, coordsW, useNormConditioning);
        
        double[] h1 = MatrixUtil.extractColumn(h, 0);
        double[] h2 = MatrixUtil.extractColumn(h, 1);
        double[] h3 = MatrixUtil.extractColumn(h, 2);

        double lambda1 = MatrixUtil.lPSum(h1, 2);
        double lambda2 = MatrixUtil.lPSum(h2, 2);
        double lambda = 2./(lambda1 + lambda2);
        //double scaleFactor = 2./(Math.sqrt(sumOfSquares(h1)) + Math.sqrt(sumOfSquares(h1)));
        log.log(LEVEL, String.format("lambda1=%.3e, lambda2=%.3e, lambda=%.3e\n", lambda1, lambda2, lambda));

        double[] t = new double[] {lambda * h[0][2], lambda * h[1][2], -lambda};

        double[] r1 = Arrays.copyOf(h1, h1.length);
        MatrixUtil.multiply(r1, 1./lambda1);

        // r1 orthogonal to r2:
        double r1DotH2 = MatrixUtil.dot(r1, h2);
        double[] r2 = new double[] {
                h2[0] - (r1[0] * r1DotH2),
                h2[1] - (r1[1] * r1DotH2),
                -h2[2] - (r1[2] * r1DotH2)
        };

        // r3 orthogonal to r1 and r2:
        double[] r3 = MatrixUtil.crossProduct(r1, r2);
        // r1, r2, and r3 are columns of R
        double[][] r = MatrixUtil.zeros(3, 3);
        for (int row = 0; row < 3; ++row) {
            r[row][0] = r1[row];
            r[row][1] = r2[row];
            r[row][2] = r3[row];
        }

        //TODO: consider converting the rotation matrix to a quaternion or Euler angles (eqns 39, 37)

        // further orthonormalization of r:
        SVDProducts svd = MatrixUtil.performSVD(r);
        r = MatrixUtil.multiply(svd.u, svd.vT);

        double detR = MatrixUtil.determinant(r);
        detR = Math.round(detR*1E11)/1E11;
        assert(Math.abs(detR - 1) < 1E-7);

        CameraExtrinsicParameters kExtr = new Camera.CameraExtrinsicParameters();
        kExtr.setRotation(r);
        kExtr.setTranslation(t);
        kExtr.setRodriguesVector(Rotation.extractRotationVectorRodriguesBouguet(r).rotVec);

        return kExtr;
    }

    /**
    apply radial distortion to distortion-free camera centered coordinates using 
    the algorithm of Ma, Chen & Moore (which is Ma et al. 2003) for the
    distortion function expressed as f(r) = 1 + k1*r + k2*r^2 (which is 
    equation #3 in Table 2 of Ma et al. 2004).
    In terms of the variables outlined below, the algorithm input is
    (x, y), k1, k2, and cameraIntrinsics and the output is (x_d, y_d).
    TODO: consider overloading this method to implement equation #4. 
    <pre>
    Ma, Chen & Moore 2004, "Rational Radial Distortion Models of Camera Lenses 
    with Analytical Solution for Distortion Correction."
    International Journal of Information Acquisition · June 2004    
    
    defining variables:
        K            :  camera intrinsics matrix
        (u_d, v_d)   :  Distorted image point in pixel
        (u, v)       :  Distortion-free image point in pixel
        (x_d, y_d)   :  [x_d, y_d, 1]^T = K^−1[u_d, v_d, 1]^T
        (x, y)       :  [x, y, 1]^T = K^−1[u, v, 1]^T
        r_d          :  r_d^2 = x_d^2 + y_d^2
        r            :  r^2 = x^2 + y^2
        k1, k2, ...  :  Radial distortion coefficients

    and the projection equation variables:
    [X_c,Y_c,Z_c]^T denotes a point in the camera frame which is related to the 
        corresponding point 
    [X_w, Y_w, Z_w]^T in the world reference frame 
        by
    P_c = R P_w + t
    R is thr rotation matrix
    t is the translation vector
    
    lambda * [u] = K [R|t] [X_w] = K [X_c]
             [v]           [Y_w]     [Y_c]
             [1]           [Z_w]     [Z_c]
                           [1  ]
    
    ==> To Apply Radial Distortion, given coefficients k1, k2 and coordinates
      (eqn 7) r_d = r * f(r) = r*(1 + k1*r + k2*r^2 + k3*r^3 + ...)
      (eqn 8) using only 2 coeffs: f(r) = (1 + k1*r + k2*r^2)
    
      Ma et al. 2003 distortion model:
         x_d = x * f(r)
         y_d = y * f(r)
    
    ==> Radial Undistortion:
        solve for cubic roots
        https://en.wikipedia.org/wiki/Cubic_equation#Reduction_to_a_depressed_cubic
        though authors use Pearson's 1983 version of 
           Handbook of Applied Mathematics: Selected Results and Methods
    
        from k2*r^3 +k1*r^2 +r - r_d = 0
        solve 
           r_bar^3 + r_bar*p + q = 0
        where 
          r_bar = r + (a/3)
          a = k1/k2
          b = 1/k2
          c = −r_d/k2
          p = b − (a^2/3)
          q = (2a^3)/27 − ab/3 + c
        then use depressed cubic root to solve for r_bar.
        r = r_bar = (a/3)
    
        After r is determined, (u,v) can be calculated from (Eqn 5) which is
           u_d - u_0 = (u−u_0) * f(r)
           v_d − v_0 = (v−v_0) * f(r)
          
    Useful reading is also:
    Drap et al, "An Exact Formula for Calculating Inverse Radial Lens Distortions"
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4934233/
        Radial distortion is caused by the spherical shape of the lens, 
        whereas tangential distortion is caused by the decentering and 
        non-orthogonality of the lens components with respect to the 
        optical axis
        ...Barrel distortion can be physically present in small focal length 
        systems, while larger focal lengths can result in pincushion distortion 
        ...
        barrel distortion corresponds to a negative value of k1.
        pincushion distortion to a positive value of k1.
         
    NOTE that Drap et al. also apply the camera intrinsics, but they use a 
    focal length of "1".
    * 
    </pre>
    @param xC distortion-free camera centered coordinates.  format is 3XN where N is the
    number of points.
    In terms of Table 1 of Chen et al. 2004, this is a double array of (x, y).
    @param k1 first radial distortion coefficient
    @param k2 second radial distortion coefficient
    @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
    f(r) = 1 +k1*r^2 + k2*r^4 if true,
    else use model #3 f(r) = 1 +k1*r + k2*r^2.
    @return  distorted camera centered coordinates in format 3XN where N is
    the number of points.
    In terms of Table 1 of Chen et al. 2004, this is a double array of (x_d, y_d).
    */
    static double[][] applyRadialDistortion(double[][] xC, double k1, double k2,
        boolean useR2R4) {
        
        if (xC.length != 3) {
            throw new IllegalArgumentException("xC.length must be 3");
        }
                
        double[][] distorted = MatrixUtil.copy(xC);
        
        //double r, r2;
        double c, signx, signy, c2p1, divc2p1, xm, ym, x, x2, y, y2;
        double[] c2s = new double[2];
        int i;
        
        for (i = 0; i < distorted[0].length; ++i) {
            //r2 = distorted[0][i]*distorted[0][i] + distorted[1][i]*distorted[1][i];
            //r = Math.sqrt(r2);
            
            // following Ma et al. 2004 Table 2,column 3 for model #3:
            /*
            for more information on rewriting r^2 in terms of only x or y, see
            Eqn 3.1 of Boas "Mathematical Methods in the Physical Sciences".
            where r^2 is (x-x0)^2 + (y-y0)^2
            let _x=x-x0  this is the notation here for camera coordinate frame
            let c^2 = (_y)^2/(_x)^2
            let c2p1 = (1 + c^2)
            then r^2 = _x^2 + _x^2*c^2 = _x^2*(1 + c^2)
                      = _x^2*c2p1
            also, for r^2 in terms of _y^2:
            let divc2p1 = ((1/c^2) + 1)
            then r^2 = _y^2*((1/c^2) + 1)
                     = _y^2*divc2p1
            */
            x = distorted[0][i];
            x2 = x*x;
            y = distorted[1][i];
            y2 = y*y;
            
            calculateC2s(x, y, c2s);
            c2p1 = c2s[0];
            divc2p1 = c2s[1];
            
            if (useR2R4) {
                /* model #4
                for x:
                 magnitude of distortion = xd-x
                                         = _x*(k1*_x^2*c2p1 + k2*_x^4*(c2p1^2))
               for y:
                 magnitude of distortion = yd-y
                                       = _y*(k1*_y^2*divc2p1 + k2*_y^4*(divc2p1^2))
                */
                xm = x*(k1*c2p1*x2 + k2*c2p1*c2p1*x2*x2);
                ym = y*(k1*divc2p1*y2 + k2*divc2p1*divc2p1*y2*y2);
            } else {
                /* model #3:
                magnitude of distortion = xd-x
                                       = _x*(k1*r + k2*r^2)
                                       = _x*(signx*k1*_x*sqrt(c2p1) + k2*_x^2*c2p1)
               magnitude of distortion = yd-y
                                       = _y*(signy*k1*_y*sqrt(divc2p1) + k2*_y^2*divc2p1)
                */
                signx = (x < 0) ? -1 : 1;
                signy = (y < 0) ? -1 : 1;
                // model #3, Table 2, column 2:
                xm = x*(signx*k1*x*Math.sqrt(c2p1) + k2*x*x*c2p1);
                ym = y*(signy*k1*y*Math.sqrt(divc2p1) + k2*y*y*divc2p1);
            }
  
            log.log(LEVEL, String.format(
                "%d) distort fx=%.4f, fy=%.4f for (x,y)=(%.3f,%.3f) => (xd,yd)=(%.3f,%.3f)\n", 
                i, (xm/distorted[0][i])+1, (ym/distorted[1][i])+1,
                distorted[0][i], distorted[1][i], distorted[0][i]+xm, distorted[1][i]+ym));
            
            distorted[0][i] += xm;
            distorted[1][i] += ym;
        }
                
        return distorted;
    }
    
    /**
    apply radial distortion to distortion-free camera centered coordinates using 
    the algorithm of Ma, Chen & Moore (which is Ma et al. 2003) for the
    distortion function expressed as f(r) = 1 + k1*r + k2*r^2 (which is 
    equation #3 in Table 2 of Ma et al. 2004).
    In terms of the variables outlined below, the algorithm input is
    (x, y), k1, k2, and cameraIntrinsics and the output is (x_d, y_d).
    TODO: consider overloading this method to implement equation #4. 
    <pre>
    Ma, Chen & Moore 2004, "Rational Radial Distortion Models of Camera Lenses 
    with Analytical Solution for Distortion Correction."
    International Journal of Information Acquisition · June 2004    
    
    defining variables:
        K            :  camera intrinsics matrix
        (u_d, v_d)   :  Distorted image point in pixel
        (u, v)       :  Distortion-free image point in pixel
        (x_d, y_d)   :  [x_d, y_d, 1]^T = K^−1[u_d, v_d, 1]^T
        (x, y)       :  [x, y, 1]^T = K^−1[u, v, 1]^T
        r_d          :  r_d^2 = x_d^2 + y_d^2
        r            :  r^2 = x^2 + y^2
        k1, k2, ...  :  Radial distortion coefficients

    and the projection equation variables:
    [X_c,Y_c,Z_c]^T denotes a point in the camera frame which is related to the 
        corresponding point 
    [X_w, Y_w, Z_w]^T in the world reference frame 
        by
    P_c = R P_w + t
    R is thr rotation matrix
    t is the translation vector
    
    lambda * [u] = K [R|t] [X_w] = K [X_c]
             [v]           [Y_w]     [Y_c]
             [1]           [Z_w]     [Z_c]
                           [1  ]
    
    ==> To Apply Radial Distortion, given coefficients k1, k2 and coordinates
      (eqn 7) r_d = r * f(r) = r*(1 + k1*r + k2*r^2 + k3*r^3 + ...)
      (eqn 8) using only 2 coeffs: f(r) = (1 + k1*r + k2*r^2)
    
      Ma et al. 2003 distortion model:
         x_d = x * f(r)
         y_d = y * f(r)
    
    ==> Radial Undistortion:
        solve for cubic roots
        https://en.wikipedia.org/wiki/Cubic_equation#Reduction_to_a_depressed_cubic
        though authors use Pearson's 1983 version of 
           Handbook of Applied Mathematics: Selected Results and Methods
    
        from k2*r^3 +k1*r^2 +r - r_d = 0
        solve 
           r_bar^3 + r_bar*p + q = 0
        where 
          r_bar = r + (a/3)
          a = k1/k2
          b = 1/k2
          c = −r_d/k2
          p = b − (a^2/3)
          q = (2a^3)/27 − ab/3 + c
        then use depressed cubic root to solve for r_bar.
        r = r_bar = (a/3)
    
        After r is determined, (u,v) can be calculated from (Eqn 5) which is
           u_d - u_0 = (u−u_0) * f(r)
           v_d − v_0 = (v−v_0) * f(r)
          
    Useful reading is also:
    Drap et al, "An Exact Formula for Calculating Inverse Radial Lens Distortions"
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4934233/
        Radial distortion is caused by the spherical shape of the lens, 
        whereas tangential distortion is caused by the decentering and 
        non-orthogonality of the lens components with respect to the 
        optical axis
        ...Barrel distortion can be physically present in small focal length 
        systems, while larger focal lengths can result in pincushion distortion 
        ...
        barrel distortion corresponds to a negative value of k1.
        pincushion distortion to a positive value of k1.
         
    NOTE that Drap et al. also apply the camera intrinsics, but they use a 
    focal length of "1".
    * 
    </pre>
    @param xCPt distortion-free camera centered coordinates for a point.  The
    length is 3.
    In terms of Table 1 of Chen et al. 2004, this is a double array of (x, y).
    @param k1 first radial distortion coefficient
    @param k2 second radial distortion coefficient
    @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
    f(r) = 1 +k1*r^2 + k2*r^4 if true,
    else use model #3 f(r) = 1 +k1*r + k2*r^2.
    @return  distorted camera centered coordinates in format 3XN where N is
    the number of points.
    In terms of Table 1 of Chen et al. 2004, this is a double array of (x_d, y_d).
    */
    static double[] applyRadialDistortion(double[] xCPt, double k1, double k2,
        boolean useR2R4) {
        
        double[] distorted = new double[3];
        
        applyRadialDistortion(xCPt, k1, k2, useR2R4, distorted);
        
        return distorted;
    }
    
    /**
    apply radial distortion to distortion-free camera centered coordinates using 
    the algorithm of Ma, Chen & Moore (which is Ma et al. 2003) for the
    distortion function expressed as f(r) = 1 + k1*r + k2*r^2 (which is 
    equation #3 in Table 2 of Ma et al. 2004).
    In terms of the variables outlined below, the algorithm input is
    (x, y), k1, k2, and cameraIntrinsics and the output is (x_d, y_d).
    TODO: consider overloading this method to implement equation #4. 
    <pre>
    Ma, Chen & Moore 2004, "Rational Radial Distortion Models of Camera Lenses 
    with Analytical Solution for Distortion Correction."
    International Journal of Information Acquisition · June 2004    
    
    defining variables:
        K            :  camera intrinsics matrix
        (u_d, v_d)   :  Distorted image point in pixel
        (u, v)       :  Distortion-free image point in pixel
        (x_d, y_d)   :  [x_d, y_d, 1]^T = K^−1[u_d, v_d, 1]^T
        (x, y)       :  [x, y, 1]^T = K^−1[u, v, 1]^T
        r_d          :  r_d^2 = x_d^2 + y_d^2
        r            :  r^2 = x^2 + y^2
        k1, k2, ...  :  Radial distortion coefficients

    and the projection equation variables:
    [X_c,Y_c,Z_c]^T denotes a point in the camera frame which is related to the 
        corresponding point 
    [X_w, Y_w, Z_w]^T in the world reference frame 
        by
    P_c = R P_w + t
    R is thr rotation matrix
    t is the translation vector
    
    lambda * [u] = K [R|t] [X_w] = K [X_c]
             [v]           [Y_w]     [Y_c]
             [1]           [Z_w]     [Z_c]
                           [1  ]
    
    ==> To Apply Radial Distortion, given coefficients k1, k2 and coordinates
      (eqn 7) r_d = r * f(r) = r*(1 + k1*r + k2*r^2 + k3*r^3 + ...)
      (eqn 8) using only 2 coeffs: f(r) = (1 + k1*r + k2*r^2)
    
      Ma et al. 2003 distortion model:
         x_d = x * f(r)
         y_d = y * f(r)
    
    ==> Radial Undistortion:
        solve for cubic roots
        https://en.wikipedia.org/wiki/Cubic_equation#Reduction_to_a_depressed_cubic
        though authors use Pearson's 1983 version of 
           Handbook of Applied Mathematics: Selected Results and Methods
    
        from k2*r^3 +k1*r^2 +r - r_d = 0
        solve 
           r_bar^3 + r_bar*p + q = 0
        where 
          r_bar = r + (a/3)
          a = k1/k2
          b = 1/k2
          c = −r_d/k2
          p = b − (a^2/3)
          q = (2a^3)/27 − ab/3 + c
        then use depressed cubic root to solve for r_bar.
        r = r_bar = (a/3)
    
        After r is determined, (u,v) can be calculated from (Eqn 5) which is
           u_d - u_0 = (u−u_0) * f(r)
           v_d − v_0 = (v−v_0) * f(r)
          
    Useful reading is also:
    Drap et al, "An Exact Formula for Calculating Inverse Radial Lens Distortions"
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4934233/
        Radial distortion is caused by the spherical shape of the lens, 
        whereas tangential distortion is caused by the decentering and 
        non-orthogonality of the lens components with respect to the 
        optical axis
        ...Barrel distortion can be physically present in small focal length 
        systems, while larger focal lengths can result in pincushion distortion 
        ...
        barrel distortion corresponds to a negative value of k1.
        pincushion distortion to a positive value of k1.
         
    NOTE that Drap et al. also apply the camera intrinsics, but they use a 
    focal length of "1".
    * 
    </pre>
    @param xCPt distortion-free camera centered coordinates for a point.  The
    length is 3.
    In terms of Table 1 of Chen et al. 2004, this is a double array of (x, y).
    @param k1 first radial distortion coefficient
    @param k2 second radial distortion coefficient
    @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
    f(r) = 1 +k1*r^2 + k2*r^4 if true,
    else use model #3 f(r) = 1 +k1*r + k2*r^2.
    @param outputDistorted  distorted camera centered coordinates in format 3XN where N is
    the number of points.
    In terms of Table 1 of Chen et al. 2004, this is a double array of (x_d, y_d).
    */
    static void applyRadialDistortion(double[] xCPt, double k1, double k2,
        boolean useR2R4, double[] outputDistorted) {
        
        if (xCPt.length != 3) {
            throw new IllegalArgumentException("xCPt.length must be 3");
        }
        if (outputDistorted.length != 3) {
            throw new IllegalArgumentException("outputDistorted.length must be 3");
        }
                        
        //double r, r2;
        double signx, signy, c2p1, divc2p1, xm, ym, x, x2, y, y2;
        double[] c2s = new double[2];
        
        //r2 = distorted[0][i]*distorted[0][i] + distorted[1][i]*distorted[1][i];
        //r = Math.sqrt(r2);

        // following Ma et al. 2004 Table 2,column 3 for model #3:
        /*
        for more information on rewriting r^2 in terms of only x or y, see
        Eqn 3.1 of Boas "Mathematical Methods in the Physical Sciences".
        where r^2 is (x-x0)^2 + (y-y0)^2
        let _x=x-x0  this is the notation here for camera coordinate frame
        let c^2 = (_y)^2/(_x)^2
        let c2p1 = (1 + c^2)
        then r^2 = _x^2 + _x^2*c^2 = _x^2*(1 + c^2)
                  = _x^2*c2p1
        also, for r^2 in terms of _y^2:
        let divc2p1 = ((1/c^2) + 1)
        then r^2 = _y^2*((1/c^2) + 1)
                 = _y^2*divc2p1
        */
        x = xCPt[0];
        x2 = x*x;
        y = xCPt[1];
        y2 = y*y;

        calculateC2s(x, y, c2s);
        c2p1 = c2s[0];
        divc2p1 = c2s[1];

        if (useR2R4) {
            /* model #4
            for x:
             magnitude of distortion = xd-x
                                     = _x*(k1*_x^2*c2p1 + k2*_x^4*(c2p1^2))
           for y:
             magnitude of distortion = yd-y
                                   = _y*(k1*_y^2*divc2p1 + k2*_y^4*(divc2p1^2))
            */
            xm = x*(k1*c2p1*x2 + k2*c2p1*c2p1*x2*x2);
            ym = y*(k1*divc2p1*y2 + k2*divc2p1*divc2p1*y2*y2);
        } else {
            /* model #3:
            magnitude of distortion = xd-x
                                   = _x*(k1*r + k2*r^2)
                                   = _x*(signx*k1*_x*sqrt(c2p1) + k2*_x^2*c2p1)
           magnitude of distortion = yd-y
                                   = _y*(signy*k1*_y*sqrt(divc2p1) + k2*_y^2*divc2p1)
            */
            signx = (x < 0) ? -1 : 1;
            signy = (y < 0) ? -1 : 1;
            // model #3, Table 2, column 2:
            xm = x*(signx*k1*x*Math.sqrt(c2p1) + k2*x*x*c2p1);
            ym = y*(signy*k1*y*Math.sqrt(divc2p1) + k2*y*y*divc2p1);
        }

        log.log(LEVEL, String.format(
            "distort fx=%.4f, fy=%.4f for (x,y)=(%.3f,%.3f) => (xd,yd)=(%.3f,%.3f)\n", 
            (xm/x)+1, (ym/y)+1,
            x, y, x+xm, y+ym));

        outputDistorted[0] = xCPt[0] + xm;
        outputDistorted[1] = xCPt[1] + ym;
        outputDistorted[2] = xCPt[2];// <=== presumably this coordinate is '1'
    }
    
    /**
    remove radial distortion from image.
    The algorithm follows Ma, Chen, & Moore 2004.
    In terms of the variables outlined in comments below, the algorithm input is
    distorted points as a double array of (x_d, x_d), and the radial distortion
    coefficients k1, k2.  The output is a double array of (x, y).
    Choices for the distortion function are models #3 and #4.
    <pre>
    Ma, Chen & Moore 2004, "Rational Radial Distortion Models of Camera Lenses 
    with Analytical Solution for Distortion Correction."
    International Journal of Information Acquisition · June 2004    
    
    defining variables:
        K            :  camera intrinsics matrix
        (u_d, v_d)   :  Distorted image point in pixel
        (u, v)       :  Distortion-free image point in pixel
        (x_d, y_d)   :  [x_d, y_d, 1]^T = K^−1[u_d, v_d, 1]^T
        (x, y)       :  [x, y, 1]^T = K^−1[u, v, 1]^T
        r_d          :  r_d^2 = x_d^2 + y_d^2
        r            :  r^2 = x^2 + y^2
        k1, k2, ...  :  Radial distortion coefficients
      
    Useful reading is also:
    Drap et al, "An Exact Formula for Calculating Inverse Radial Lens Distortions"
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4934233/
        Radial distortion is caused by the spherical shape of the lens, 
        whereas tangential distortion is caused by the decentering and 
        non-orthogonality of the lens components with respect to the 
        optical axis
        ...Barrel distortion can be physically present in small focal length 
        systems, while larger focal lengths can result in pincushion distortion 
        ...
        barrel distortion corresponds to a negative value of k1.
        pincushion distortion to a positive value of k1.
         
    </pre>
    @param xC distorted points in the camera reference frame, presumably 
    already center subtracted.  format is 3XN where N is the
    number of points.  These are (x_d, x_d) pairs in terms of Table 1 in Ma et al. 2004.
    @param k1 first radial distortion coefficient
    @param k2 second radial distortion coefficient
    @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
        f(r) = 1 +k1*r^2 + k2*r^4 if true,
        else use model #3 f(r) = 1 +k1*r + k2*r^2.
        note that if rCoeffs is null or empty, no radial distortion is removed.
    @return undistorted points in the camera reference frame.  Format is 3XN where N is the
    number of points.  These are (x, y) pairs in terms of Table 1 in Ma et al. 2004.
    @throws no.uib.cipr.matrix.NotConvergedException
    */
    static double[][] removeRadialDistortion(double[][] xC, double k1, double k2,
        boolean useR2R4) throws NotConvergedException, IOException {
        
        if (useR2R4) {
            return removeRadialDistortion4(xC, k1, k2);
        }
        return removeRadialDistortion(xC, k1, k2);
    }
    
    /**
    remove radial distortion from image.
    The algorithm follows Ma, Chen, & Moore 2004 and is for model #3,
    * f(r) = 1 + k1*r + k2*r^2.
    In terms of the variables outlined in comments below, the algorithm input is
    distorted points as a double array of (x_d, x_d), and the radial distortion
    coefficients k1, k2.  The output is a double array of (x, y).
    <pre>
    Ma, Chen & Moore 2004, "Rational Radial Distortion Models of Camera Lenses 
    with Analytical Solution for Distortion Correction."
    International Journal of Information Acquisition · June 2004    
    
    defining variables:
        K            :  camera intrinsics matrix
        (u_d, v_d)   :  Distorted image point in pixel
        (u, v)       :  Distortion-free image point in pixel
        (x_d, y_d)   :  [x_d, y_d, 1]^T = K^−1[u_d, v_d, 1]^T
        (x, y)       :  [x, y, 1]^T = K^−1[u, v, 1]^T
        r_d          :  r_d^2 = x_d^2 + y_d^2
        r            :  r^2 = x^2 + y^2
        k1, k2, ...  :  Radial distortion coefficients
      
    Useful reading is also:
    Drap et al, "An Exact Formula for Calculating Inverse Radial Lens Distortions"
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4934233/
        Radial distortion is caused by the spherical shape of the lens, 
        whereas tangential distortion is caused by the decentering and 
        non-orthogonality of the lens components with respect to the 
        optical axis
        ...Barrel distortion can be physically present in small focal length 
        systems, while larger focal lengths can result in pincushion distortion 
        ...
        barrel distortion corresponds to a negative value of k1.
        pincushion distortion to a positive value of k1.
    </pre>
    @param xC distorted points in the camera reference frame, presumably 
    already center subtracted.  format is 3XN where N is the
    number of points.  These are (x_d, x_d) pairs in terms of Table 1 in Ma et al. 2004.
    @param k1 first radial distortion coefficient
    @param k2 second radial distortion coefficient
    @return undistorted points in the camera reference frame.  Format is 3XN where N is the
    number of points.  These are (x, y) pairs in terms of Table 1 in Ma et al. 2004.
    @throws no.uib.cipr.matrix.NotConvergedException
    */
    static double[][] removeRadialDistortion(double[][] xC, double k1, double k2) 
        throws NotConvergedException, IOException {
        
        if (xC.length != 3) {
            throw new IllegalArgumentException("xC.length must be 3");
        }
          
        // output array that will be corrected
        double[][] corrected = MatrixUtil.copy(xC);
        
        if (Math.abs(k1) < eps) {
            if (Math.abs(k2) < eps) {
                return corrected;
            }
            throw new IllegalArgumentException("k1 must be non-zero for this algorithm.  "
                + "To solve for the 2nd order term only, use useR2R4=true");
        }
        
        double[] coeffs = new double[]{k2, k1, 1, 0};
        
        if (Math.abs(k2) < eps) {
            // k2 is 0 so solve for reduced order of coeffs:
            coeffs = new double[]{k1, 1, 0};
        }
        double[] rBar;
        double rd, r, fr, theta;
        double p, q, a, b, c, signx, signy;
        //double c2p1, signx, signy, fx, fy;
        int i;
        for (i = 0; i < xC[0].length; ++i) {
            rd = Math.sqrt(corrected[0][i]*corrected[0][i] + corrected[1][i]*corrected[1][i]);
            if (Math.abs(rd) < eps) {
                continue;
            }
            coeffs[coeffs.length - 1] = -rd;
            
            /*
            Ma et al. 2004, reduction to cubic root
            k2*r^3 +k1*r^2 + r - r_d = 0
            
            NOTE: checked the coefficient factoring
            a = k1/k2, b = 1/k2, and c = -rd/k2
            let rbar = r + (a/3)
            then r = rbar - (a/3)
            let p = b - (a*a/3)
            let q = 2*(a*a*a/27) - (a*b/3) + c
   
            rbar^3 + p*rbar + q = 0
            */
            //a = k1/k2;
            //b = 1./k2;
            //c = -rd/k2;
            //p = b - (a*a/3);
            //q = 2*(a*a*a/27) - (a*b/3) + c;
      
            log.log(LEVEL, String.format("\ni=%d\n",i));
            
            //rBar = CubicRootSolver.solve(coeffs);
            rBar = PolynomialRootSolver.solveForRealUsingMPSolve(coeffs, eps2);
            if (rBar != null)
                log.log(LEVEL, String.format("rBar=%s\n", FormatArray.toString(rBar, "%.4e")));
            if (rBar == null || rBar.length == 0) {
                r = rd;
                // check solution: 
                //  k2*r^3 +k1*r^2 +r - r_d = 0
                double chk = k2 * r * r * r + k1 * r * r + r - rd;
                log.log(LEVEL, String.format("chk 0=%.4e\n", chk));
                //assert(Math.abs(chk) < tol);
            } else {
                log.log(LEVEL, String.format("rBar=%s\n", FormatArray.toString(rBar, "%.4f")));
                r = rBar[0];
                if (r < 0 && rBar.length > 1) {
                    for (int ii = 1; ii < rBar.length; ++ii) {
                        r = rBar[ii];
                        if (r > 0) {
                            break;
                        }
                    }
                }
                // check solution: 
                //  k2*r^3 +k1*r^2 +r - r_d = 0
                //double chk = k2 * r * r * r + k1 * r * r + r - rd;
                //log.log(LEVEL, String.format("chk 0==%.4e\n", chk));
                //assert(Math.abs(chk) < tol);
            }
                       
            // remove radial distortion
            // we have r now and have ud, vd
            // model #3: fr = 1 + k1*r + k2*r*r;
            
            // eqn (5) from Ma et al. 2004
            // (ud - u0) = (u - u0) * fr
            // where (ud, vd) are the real observed image points.
            // since the "apply radial distortion" is to the points in the camera
            // reference frame, the removal must be also in order for the radial
            // coefficients to be of the right scale.
            //
            // so one can look at eqn (4) instead:
            //    xd = x*fr
            //    where x is the undistorted in the camera reference frame.
            // solving for x:
            //    x = xd/fr
            //fr = 1 + k1*r + k2*r*r;
            theta = Math.atan2(corrected[1][i], corrected[0][i]);
            
            log.log(LEVEL, String.format("(xd,yd)=(%.4f,%.4f)  (x,y)=(%.4f,%.4f)\n",
                corrected[0][i], corrected[1][i], 
                r*Math.cos(theta), r*Math.sin(theta)));
            
            corrected[0][i] = (r*Math.cos(theta));
            corrected[1][i] = (r*Math.sin(theta));
        }
                
        return corrected;
    }
    
    /**
    remove radial distortion from image.
    The algorithm follows Ma, Chen, & Moore 2004 and is for model #4,
    * f(r) = 1 + k1*r^2 + k2*r^4.
    In terms of the variables outlined in comments below, the algorithm input is
    distorted points as a double array of (x_d, x_d), and the radial distortion
    coefficients k1, k2.  The output is a double array of (x, y).
    <pre>
    Ma, Chen & Moore 2004, "Rational Radial Distortion Models of Camera Lenses 
    with Analytical Solution for Distortion Correction."
    International Journal of Information Acquisition · June 2004    
    
    defining variables:
        K            :  camera intrinsics matrix
        (u_d, v_d)   :  Distorted image point in pixel
        (u, v)       :  Distortion-free image point in pixel
        (x_d, y_d)   :  [x_d, y_d, 1]^T = K^−1[u_d, v_d, 1]^T
        (x, y)       :  [x, y, 1]^T = K^−1[u, v, 1]^T
        r_d          :  r_d^2 = x_d^2 + y_d^2
        r            :  r^2 = x^2 + y^2
        k1, k2, ...  :  Radial distortion coefficients
      
    Useful reading is also:
    Drap et al, "An Exact Formula for Calculating Inverse Radial Lens Distortions"
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4934233/
        Radial distortion is caused by the spherical shape of the lens, 
        whereas tangential distortion is caused by the decentering and 
        non-orthogonality of the lens components with respect to the 
        optical axis
        ...Barrel distortion can be physically present in small focal length 
        systems, while larger focal lengths can result in pincushion distortion 
        ...
        barrel distortion corresponds to a 
            negative value of k1.
            present in small focal length systems.
        pincushion distortion to a 
            positive value of k1.
            present in larger focal length systems
            *
    </pre>
    @param xC distorted points in the camera reference frame, presumably 
    already center subtracted.  format is 3XN where N is the
    number of points.  These are (x_d, x_d) pairs in terms of Table 1 in Ma et al. 2004.
    @param k1 first radial distortion coefficient
    @param k2 second radial distortion coefficient
    @return undistorted points in the camera reference frame.  Format is 3XN where N is the
    number of points.  These are (x, y) pairs in terms of Table 1 in Ma et al. 2004.
    @throws no.uib.cipr.matrix.NotConvergedException
    */
    static double[][] removeRadialDistortion4(double[][] xC, double k1, double k2) 
        throws NotConvergedException, IOException {
        
        if (xC.length != 3) {
            throw new IllegalArgumentException("xC.length must be 3");
        }
        
        // 5th order root finding
        // a polynomial of degree n has at most n real or complex roots
                
        double[][] corrected = MatrixUtil.copy(xC);
        
        if (Math.abs(k1) < eps) {
            if (Math.abs(k2) < eps) {
                return corrected;
            }
            throw new IllegalArgumentException("k1 must be non-zero for this algorithm.");
        }
        
        log.log(LEVEL, String.format("distorted coords=\n%s\n", FormatArray.toString(xC, "%.3f")));
        
        //model#4: k2*r^5 + k1*r^3 + r - rd = 0.
        //         (k2*(c2p1^2))*_x^5 + (k1*c2p1)*_x^3 + _x - xd = 0
        //         (k2*(divc2p1^2))*_y^5 + (k1*divc2p1)*_y^3 + _y - yd = 0
        
        double rd, r, theta, x, y, c, c2p1, divc2p1, fx, fy;
        int i;
        double[] rootsX, rootsY;
        double[] coeffsX = new double[]{0, 0, 0, 0, 1, 0};
        double[] coeffsY = new double[]{0, 0, 0, 0, 1, 0};
        double[] c2s = new double[2];
        if (Math.abs(k2) < eps) {
            // k2 is 0 so solve for reduced order of coeffs:
            coeffsX = new double[]{0, 0, 1, 0};
            coeffsY = new double[]{0, 0, 1, 0};
        }
        
        for (i = 0; i < xC[0].length; ++i) {
            rd = Math.sqrt(corrected[0][i]*corrected[0][i] + corrected[1][i]*corrected[1][i]);
            if (Math.abs(rd) < eps) {
                continue;
            }
            
            //for more information on rewriting r^2 in terms of only x or y, see
            //Eqn 3.1 of Boas "Mathematical Methods in the Physical Sciences".
            //c2p1 = (y/x)^2 + 1;
            //divc2p1 = (x/y)^2 + 1;
            
            //handle cases where x=0 or y=0 
            calculateC2s(corrected[0][i], corrected[1][i], c2s);
            c2p1 = c2s[0];
            divc2p1 = c2s[1];
            
            log.log(LEVEL, String.format("\ni=%d\n",i));
            // solve for x in (k2*(c2p1^2))*_x^5 + (k1*c2p1)*_x^3 + _x - xd = 0
            if (coeffsX.length == 6) {
                coeffsX[0] = k2 * c2p1 * c2p1;
                coeffsX[2] = k1 * c2p1;
            } else {
                coeffsX[0] = k1 * c2p1;
            }
            coeffsX[coeffsX.length - 1] = -corrected[0][i];
            
            if (Math.abs(corrected[0][i]) < eps) {
                x = corrected[0][i];
            } else {
                rootsX = PolynomialRootSolver.solveForRealUsingMPSolve(coeffsX, eps2);
                if (rootsX == null || rootsX.length == 0) {
                    x = corrected[0][i];
                    log.log(LEVEL, String.format("rootsX=null\n"));
                } else {
                    log.log(LEVEL, String.format("poly rootsX=%s\n", FormatArray.toString(rootsX, "%.4f")));
                    if (corrected[0][i] < 0) {
                        x = rootsX[0];
                    } else {
                        x = rootsX[rootsX.length - 1];
                    }
                }
            }
            //xd = x*(1 + k1*c2p1*x2 + k2*c2p1*c2p1*x2*x2)
            fx = (k2*(c2p1*c2p1))*Math.pow(x, 4) + (k1*c2p1)*Math.pow(x, 2) + 1;
            // xd = x*fx
            double chkX = fx * x;
            
            // solve for y in (k2*(divc2p1^2))*_y^5 + (k1*divc2p1)*_y^3 + _y - yd = 0
            if (coeffsY.length == 6) {
                coeffsY[0] = k2 * divc2p1 * divc2p1;
                coeffsY[2] = k1 * divc2p1;
            } else {
                coeffsY[0] = k1 * divc2p1;
            }
            coeffsY[coeffsY.length - 1] = -corrected[1][i];
            
            if (Math.abs(corrected[1][i]) < eps) {
                y = corrected[1][i];
            } else {
                rootsY = PolynomialRootSolver.solveForRealUsingMPSolve(coeffsY, eps2);
                if (rootsY == null || rootsY.length == 0) {
                    y = corrected[1][i];
                    log.log(LEVEL, String.format("rootsY=null\n"));
                } else {
                    log.log(LEVEL, String.format("poly rootsY=%s\n", FormatArray.toString(rootsY, "%.4f")));
                    if (corrected[1][i] < 0) {
                        y = rootsY[0];
                    } else {
                        y = rootsY[rootsY.length - 1];
                    }
                }
            }
            fy = (k2*(divc2p1*divc2p1))*Math.pow(y, 4) + (k1*divc2p1)*Math.pow(y, 2) + 1;
            // yd = y*fy
            double chkY = fy * y;
            log.log(LEVEL, String.format("(xd,yd)=(%.4f,%.4f)  (x,y)=(%.4f,%.4f)\n",
                corrected[0][i], corrected[1][i], x, y));
            log.log(LEVEL, String.format("fx=%.4f, fy=%.4f, checkX:%.4f==%.4f? checkY:%.4f==%.4f?\n",
                fx, fy, chkX, xC[0][i], chkY, xC[1][i]));
            
            corrected[0][i] = x;
            corrected[1][i] = y;
            
        }
                
        return corrected;
    }

    /**
     * calculate the projection of world features in coordsW by the
     * homography h into the image plane, storing the results in ud and vd.
     * The method follows eqn (17) of Ma, Chen, & Moore 2003 "Camera Calibration".
     * @param coordsW the coordinates of the features in world reference frame.
     *    size is 3 X n.
     * @param h the homography.  size is (nImages*3) X 3
     * @param u output projected image x coordinates for all images.
     *     length is (n*nImages).
     * @param v output projected image y coordinated for all images.
     *     length is  (n*nImages).
     */
    static void calculateProjected(double[][] coordsW, double[][] h, 
        double[] u, double[] v) {
        
        // n is the number of features
        int n = coordsW[0].length;
        int nImages = h.length/3;
        
        //u, v are 1 X (n*nImages)
        //h is nImages*3 X 3
        //coordsW is 3 X n
        
        // eqn (17) denom = h[2][0]*X_w + h[2][1]*Y_w + h[2][2]
        //          ud = (h[0][0]*X_w + h[0][1]*Y_w + h[0][2])/denom
        //          vd = (h[1][0]*X_w + h[1][1]*Y_w + h[1][2])/denom
        
        double[] xw1 = new double[3];
        xw1[2] = 1;
        double denom;
        double[] h0, h1, h2;
        int i, j;
        for (i = 0; i < nImages; ++i) {
            h0 = h[i*3 + 0];
            h1 = h[i*3 + 1];
            h2 = h[i*3 + 2];
            for (j = 0; j < n; ++j) { // n features
                xw1[0] = coordsW[0][j];
                xw1[1] = coordsW[1][j];
                denom = MatrixUtil.innerProduct(h2, xw1);
                u[i*n + j] = MatrixUtil.innerProduct(h0, xw1);
                u[i*n + j] /= denom;
                v[i*n + j] = MatrixUtil.innerProduct(h1, xw1);
                v[i*n + j] /= denom;
            }
        }        
    }

    /**
     * 
     * @param uvD (ud, vd) are the Real observed distorted image points.
     * uvD holds the features in each image in pixel coordinates ordered 
     * such that all features of one image are followed by all features
     * of the next image.  the x-axis coordinates are in row 0.
     * the y-axis coordinates are in row 1.  the third row is all 1's.
               It is a 2 dimensional double array of size
               3 X (N*n) where N is the number of images.
               In Table 1 of Ma, Chen, & Moore 2003 "Camera Calibration"
               these are the (u_d, v_d) pairs.
     * @param u projections of the WCS feature x coordinates into the image
     * reference frame.  array length is n*nImages
     * @param v projections of the WCS feature y coordinates into the image
     * reference frame.  array length is n*nImages
     * @param cameraMatrices data structure holding the camera intrinsic parameters
     * and the extrinsic parameter matrices for each image.
     * @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
    f(r) = 1 +k1*r^2 + k2*r^4 if true,
    else use model #3 f(r) = 1 +k1*r + k2*r^2.
    * @throws Exception if there is an error in use of MPSolver during the
     * removal of radial distortion, a generic exception is thrown with the
     * error message from the MPSolve documentation.
     * @return 
     */
    static double[] solveForRadialDistortion(double[][] uvD, 
        double[] u, double[] v, 
        CameraMatrices cameraMatrices, boolean useR2R4) 
        throws NotConvergedException, Exception {
        
        int nImages = cameraMatrices.getExtrinsics().size();
        int nFeatures = u.length/nImages;
                
        /* 
        (ud, vd) are Real observed distorted image points in image reference frame.
        (u, v) Ideal projected undistorted image points in image reference frame (the projection of
               the points from the world reference frame to the image reference frame).
        [x,y,1] = A^-1 * [u, v, 1] are transformed to camera reference frame.

        eqn (5) of Ma, Chen, & Moore 2004, "Rational Radial Distortion..."
           ud-u0 = (u-u0)*f(r)
           vd-v0 = (v-v0)*f(r)

        eqn (8) of Ma, Chen & Moore 2003, "Camera Calibration..."
           ud = u + (u−u0)*f_r
           vd = v + (v−v0)*f_r

        eqn (11) of Zhang 1998, "Flexible Camera Calibration ..."
           ud = u + (u−u0)*[k1*r + k2*r^2]
           vd = v + (v−v0)*[k1*r + k2*r^2]

        (5) and (8) use equation #3 or #4 of Ma et al. 2004 Table 2
            #3: f_r = 1 + k_1*r + k_2*r^2
                    = 1 + k_1*(x^2 + y^2)^-1/2 + k_2*(x^2 + y^2)
            #4: f_r = 1 + k_1*r^2 + k_2*r^4
                    = 1 + k_1*(x^2 + y^2) + k_2*(x^2 + y^2)^2

        factored out eqn (5) of Ma, Chen, & Moore 2004:
           ud-u0 = (u-u0)*(1 + k_1*r + k_2*r^2)
                 = (u-u0) + (u-u0)*(k_1*r + k_2*r^2)
           ud = u + (u-u0)*(k_1*r + k_2*r^2)
           ud-u = (u-u0)*(k_1*r + k_2*r^2)
         is the same as eqn (11) of Zhang 1998.

        Given n points in nImages, we can stack all equations together
        to obtain totally 2Nn equations in matrix form as
        Dk = d, where k = [k1, k2]^T .

        The linear least-square solutions for k is k = (D^T*D)^−1*D^T*d = pseudoInv(D)*d.

          if choose #3:

               k1                        k2                    const
              ----------------------------------------------------------
        D = [ (u-u0)*sqrt(x^2 + y^2)    (u-u0)*(x^2 + y^2) ]   d = [ ud - u ]
            [ (v-v0)*sqrt(x^2 + y^2)    (v-v0)*(x^2 + y^2) ]       [ vd - v ]

          if choose #4:

               k1                   k2                       const
              ----------------------------------------------------------
        D = [ (u-u0)*(x^2 + y^2)    (u-u0)*(x^2 + y^2)^2 ]   d = [ ud - u ]
            [ (v-v0)*(x^2 + y^2)    (v-v0)*(x^2 + y^2)^2 ]       [ vd - v ]
        
       The linear least-square solutions for k is k = (D^T*D)^−1*D^T*d.
                                                      (2X2nN * 2nNX2)^-1 * (2X2nN) * (2nNX1)
                                                      (2X2)              * (2X2nN) * (2nNX1)
                                                      (2X2nN) * (2nNX1) = 2X1
        */
        
        int i, j;
        double ui, vi, udi, vdi, xi, yi, xi2, yi2;
        double signx, signy, c2p1, divc2p1;
        double[] c2s = new double[2];
        double ud0 = cameraMatrices.getIntrinsics().getIntrinsic()[0][2];
        double vd0 = cameraMatrices.getIntrinsics().getIntrinsic()[1][2];
        double u0 = ud0; 
        double v0 = vd0;
        double[][] xy;
        double[][] dM = new double[2*nFeatures*nImages][2];
        double[] dV = new double[2*nFeatures*nImages];
        for (i = 0; i < nImages; ++i) {
            xy = MatrixUtil.copySubMatrix(uvD, 0, 2, nFeatures*i, nFeatures*(i + 1)-1);
            xy = Camera.pixelToCameraCoordinates(xy, cameraMatrices.getIntrinsics());
            for (j = 0; j < nFeatures; ++j) {
                ui = u[nFeatures*i + j];
                vi = v[nFeatures*i + j];
                udi = uvD[0][nFeatures*i + j];
                vdi = uvD[1][nFeatures*i + j];
                xi = xy[0][j];
                yi = xy[1][j];
                xi2 = xi*xi;
                yi2 = yi*yi;
        
                
                calculateC2s(xi, yi, c2s);
                c2p1 = c2s[0];
                divc2p1 = c2s[1];
                                
                // e.g. nFeatures=3
                //i:0 j:0          idx=0, idy=1
                //i:0 j:1          idx=2, idy=3
                //i:0 j:2          idx=4, idy=5
                //i:1 j:0  idx=6, idy=7
                //i:1 j:1  idx=8, idy=9
      
                if (useR2R4) {
                    dM[2*nFeatures*i + 2*j] = new double[]{
                        (ui-u0)*c2p1*xi2, (ui-u0)*c2p1*c2p1*xi2*xi2};
                    dM[2*nFeatures*i + 2*j + 1] = new double[]{
                        (vi-v0)*divc2p1*yi2, (vi-v0)*divc2p1*divc2p1*yi2*yi2};
                } else {
                    signx = (xi < 0) ? -1 : 1;
                    signy = (yi < 0) ? -1 : 1;
                    dM[2*nFeatures*i + 2*j] = new double[]{
                        (ui-u0)*signx*xi*Math.sqrt(c2p1), (ui-u0)*xi2*c2p1};
                    dM[2*nFeatures*i + 2*j + 1] = new double[]{
                        (vi-v0)*signy*yi*Math.sqrt(divc2p1), (vi-v0)*yi2*divc2p1};
                }
                dV[2*nFeatures*i + 2*j] = udi - ui;
                dV[2*nFeatures*i + 2*j + 1] = vdi - vi;
            }
        }
           
        //k = (D^T*D)^−1*D^T*d = pseudoInv(D) * d
        double[][] dInv = MatrixUtil.pseudoinverseFullColumnRank(dM);
        double[] k = MatrixUtil.multiplyMatrixByColumnVector(dInv, dV);
                
        return k;
    }

    static double[][] solveForHomographies(double[][] coordsI, 
        double[][] coordsW, int n, int nImages, boolean useNormConditioning) throws NotConvergedException {
        
        // n is the number of features
        
        double[][] h = MatrixUtil.zeros(nImages*3, 3);
        
        double[][] g, cI;
        int i, i2;
        for (i = 0; i < nImages; ++i) {
            cI = MatrixUtil.copySubMatrix(coordsI, 0, 2, n*i, n*(i + 1)-1);
            
            //3X3  and contains intrinsic camera information from the point-to-point mappings
            g = solveForHomography(cI, coordsW, useNormConditioning);

            /*
            h image 0// 3X3 in rows 0:3
            h image 1// 3X3 in rows 3:6
            h image 2// 3X3 in rows 6:9
            */
            for (i2 = 0; i2 < 3; ++i2) {
                System.arraycopy(g[i2], 0, h[i*3 + i2], 0, g[i2].length);
            }
        }
        
        return h;
    }

    static List<CameraExtrinsicParameters> solveForExtrinsics(
        CameraIntrinsicParameters kIntr, double[][] h, int nImages) 
        throws NotConvergedException, Exception {
        
        List<CameraExtrinsicParameters> list = new ArrayList<CameraExtrinsicParameters>();
        //(3) for each image homography and inverse intrinsic parameter matrix,
        //    estimate the extrinsic parameters for the pose of the camera for that image.
        CameraExtrinsicParameters rtExtr;
        
        int i;
        double[][] g;
        
        for (i = 0; i < nImages; ++i) {
            //                               ri,      rf,ci,cf
            g = MatrixUtil.copySubMatrix(h, 3*i, 3*i + 2, 0, 2);
            
            //NOTE: here, internal to solveForExtrinsic() would be a different place 
            // where one could remove radial distortion.
            // the method forms the image of the "absolute conic"
            rtExtr = solveForExtrinsic(kIntr, g);
            
            list.add(rtExtr);
        }
        
        return list;
    }

    /**
     * calculate c2p1 and divc2p1
     * where c = y/x
     * and c2p1 = c*c + 1
     * and divc2p1 = (1./(c*c)) + 1.
     * the method handles cases where x or y are 0.
     * @param x
     * @param y
     * @param output array of length 2 to return []{c2p1, divc2p1};
     */
    private static void calculateC2s(double x, double y, double[] output) {
        double c, c2p1, divc2p1;
        if (Math.abs(x) < eps) {
            c2p1 = 0;
        } else {
            if (Math.abs(y) < eps) {
                // r^2 = _x^2
                c2p1 = 1;
            } else {
                c = y / x;
                c2p1 = c * c + 1;
            }
        }
        if (Math.abs(y) < eps) {
            divc2p1 = 0;
        } else {
            if (Math.abs(x) < eps) {
                // r^2 = _y^2
                divc2p1 = 1;
            } else {
                c = y / x;
                divc2p1 = (1. / (c * c)) + 1;
            }
        }
        output[0] = c2p1;
        output[1] = divc2p1;
    }
}
