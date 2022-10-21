package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Camera.CameraExtrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraIntrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraParameters;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import no.uib.cipr.matrix.RQ;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;

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
     * given a set of features in image space and world coordinate space with
     * known camera intrinsic parameters, estimate the camera pose, that is
     * extract the camera extrinsic parameters of rotation and translation.
     * calibrating the camera extrinsic parameters is a.k.a. 
     * perspective-n-point-problem where n is the number of features (a.k.a. points).
     * This method uses DLT and could be followed by non-linear optimization
     * to improve the parameter estimates.
     * This is also known as estimating the Motion.
     * <pre>
     * references:
     *  Kitani lecture notes http://www.cs.cmu.edu/~16385/s17/Slides/11.3_Pose_Estimation.pdf
     *
     *  a few lines of code from source provided by book authors:
     *  @Book{Hartley2004,
     *     author = "Hartley, R.~I. and Zisserman, A.",
     *     title = "Multiple View Geometry in Computer Vision",
     *     edition = "Second",
     *     year = "2004",
     *     publisher = "Cambridge University Press, ISBN: 0521540518"
     * }
     * https://www.robots.ox.ac.uk/~vgg/hzbook/
     * vgg_multiview/vgg_KR_from_P.m
     * </pre>
     * @param x the images coordinates of the features in pixels in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 3 features are needed to 
     * calculate the extrinsic parameters.
     * NOTE x and X should both be distortion-free or both should be distorted.
     * @param X the world coordinates of the features in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 3 features are needed to
     * calculate the extrinsic parameters.
     * NOTE x and X should both be distortion-free or both should be distorted.
     * @return 
     */
    public static CameraParameters calculatePoseUsingDLT(double[][] x, double[][] X) 
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
        
        // normalize by last coordinate:
        /*for (int i = 0; i < x[0].length; ++i) {
            x[0][i] /= x[2][i];
            x[1][i] /= x[2][i];
        }*/
                
        // 2*n X 12       
        double xi, yi, Xi, Yi, Zi;
        double[][] ell = new double[2*n][12];
        for (int i = 0; i < n; ++i) {
            xi = x[0][i];
            yi = x[1][i];
            Xi = X[0][i];
            Yi = X[1][i];
            Zi = X[2][i];
            // http://www.cs.cmu.edu/~16385/s17/Slides/11.3_Pose_Estimation.pdf
            ell[2*i]     = new double[]{Xi, Yi, Zi, 1, 0, 0, 0, 0, -xi*Xi, -xi*Yi, -xi*Zi, -xi};
            ell[2*i + 1] = new double[]{0, 0, 0, 0, Xi, Yi, Zi, 1, -yi*Xi, -yi*Yi, -yi*Zi, -yi};
        }

        // SVD.V is 12 X 12.  SVD.U is 2n X 2n
        MatrixUtil.SVDProducts svd = MatrixUtil.performSVD(ell);
        
        // vT is 12X12.  last row in vT is the eigenvector for the smallest eigenvalue
        // it is also the epipole e1, defined as the right nullspace
        double[] xOrth = svd.vT[svd.vT.length - 1];
                
        // reshape into 3 X 4
        double[][] P2 = MatrixUtil.zeros(3, 4);
        System.arraycopy(xOrth, 0, P2[0], 0, 4);
        System.arraycopy(xOrth, 4, P2[1], 0, 4);
        System.arraycopy(xOrth, 8, P2[2], 0, 4);

        MatrixUtil.SVDProducts svdP = MatrixUtil.performSVD(P2);
        double[] c = Arrays.copyOf(svdP.vT[svdP.vT.length - 1], svdP.vT[0].length);

        // assert P2*c = 0
        double[] check0 = MatrixUtil.multiplyMatrixByColumnVector(P2, c);
        System.out.printf("check that P2*c=0:%s\n", FormatArray.toString(check0, "%.3e"));
        
        /*
        Szeliski Sect 6.2.1
        Since K is by convention upper-triangular 
        (see the discussion in Section 2.1.5), both K and R can be obtained 
        from the front 3 ⇥ 3 sub-matrix of P using RQ factorization 
        (Golub and Van Loan 1996)
        */

        // [3X3]
        double[][] M = MatrixUtil.copySubMatrix(P2, 0, 2, 0, 2);
        RQ rq = RQ.factorize(new DenseMatrix(M));
        
        System.out.printf("RQ.R=intr\n%s\n", FormatArray.toString(
            MatrixUtil.convertToRowMajor(rq.getR()), "%.3e"));
        System.out.printf("RQ.Q=rot=\n%s\n", FormatArray.toString(
            MatrixUtil.convertToRowMajor(rq.getQ()), "%.3e"));
        
        double[][] kIntr = MatrixUtil.convertToRowMajor(rq.getR());
        MatrixUtil.multiply(kIntr, 1./kIntr[2][2]);

        double[][] kExtrRot = MatrixUtil.convertToRowMajor(rq.getQ());
        System.out.printf("  decomposed into intrinsic=\n   %s\n", FormatArray.toString(kIntr, "%.3e"));
        System.out.printf("  decomposed into extrinsic rotation=\n   %s\n", FormatArray.toString(kExtrRot, "%.3e"));
        System.out.printf("  decomposed into extrinsic translation=\n   %s\n", FormatArray.toString(c, "%.3e"));

        if (kIntr[0][0] < 0) {
            // from vgg_multiview/vgg_KR_from_P.m
            double[][] diag = MatrixUtil.createIdentityMatrix(3);
            diag[0][0] = -1;
            diag[1][1] = -1;
            kIntr = MatrixUtil.multiply(kIntr, diag);
            kExtrRot = MatrixUtil.multiply(diag, kExtrRot);
            System.out.printf("  decomposed into intrinsic=\n   %s\n", FormatArray.toString(kIntr, "%.3e"));
            System.out.printf("  decomposed into extrinsic rotation=\n   %s\n", FormatArray.toString(kExtrRot, "%.3e"));
        }

        // from vgg_multiview/vgg_KR_from_P.m
        // t = -P(:,1:N)\P(:,end) = -M\P(:,end) = inv(-M)*P
        double[][] invM = MatrixUtil.copy(M);
        MatrixUtil.multiply(invM, -1);
        invM = MatrixUtil.pseudoinverseRankDeficient(invM);
        double[] t = MatrixUtil.extractColumn(P2, 3);
        t = MatrixUtil.multiplyMatrixByColumnVector(invM, t);

        System.out.printf("t=\n%s\n", FormatArray.toString(t, "%.3e"));
        System.out.printf("c=\n%s\n", FormatArray.toString(c, "%.3e"));

        CameraExtrinsicParameters extrinsics = new CameraExtrinsicParameters();
        extrinsics.setRotation(kExtrRot);
        extrinsics.setTranslation(c);
        
        CameraIntrinsicParameters intrinsics = new CameraIntrinsicParameters();
        intrinsics.setIntrinsic(kIntr);
        CameraParameters camera = new CameraParameters(intrinsics, extrinsics);
        
        return camera;
    }
    
    /**
     * given a set of features in image space and world coordinate space with
     * known camera intrinsic parameters, estimate the camera pose, that is
     * extract the camera extrinsic parameters.
     * calibrating the camera extrinsic parameters is a.k.a. 
     * perspective-n-point-problem where n is the number of features (a.k.a. points).
     * It's also called planar homography decomposition.
     * This method uses DLT and could be followed by non-linear optimization
     * to improve the parameter estimates.
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
        if (n < 3) {
            throw new IllegalArgumentException("x must have at least 3 correspondences");
        }
        if (X[0].length != n) {
            throw new IllegalArgumentException("the number of columns in X must be the same as in x");
        }
        
        int nImages = 1;
        
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
