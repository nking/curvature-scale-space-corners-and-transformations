package algorithms.imageProcessing.transform;

/**
 given a set of features in image space and world coordinate space with
  known camera intrinsic parameters, estimate the camera pose, that is
  extract the camera extrinsic parameters from the pose.
 * 
 * TODO: consider solving with M-estimators.
 * see http://research.microsoft.com/en- us/um/people/zhang/INRIA/Publis/Tutorial-Estim/node24.html
 * 
 * @author nichole
 */
public class CameraPose {
    
    /**
     * given a set of features in image space and world coordinate space with
     * known camera intrinsic parameters, estimate the camera pose, that is
     * extract the camera extrinsic parameters from the pose.
     * calibrating the camera extrinsic parameters is a.k.a. 
     * perspective-n-point-problem where n is the number of features (a.k.a. points).
     * 
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
     * @param X the world coordinates of the features in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 3features are needed to 
     * calculate the extrinsic parameters.
     * @return 
     */
    public static Camera.CameraMatrices calculatePose(
        Camera.CameraIntrinsicParameters intrinsics, double[][] x,
        double[][] X) {
        
        editing... considering the versions which need 3 points (given K, solves R,T) and one that needs 6 points (solves for K, R, T)
        
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
        
        /*
        Szeliski:
        (1) use DLT to solve the projection
            The resulting algorithm is called the direct linear transform (DLT) and is commonly attributed to Sutherland (1974).
        (2) use non-linear least squares on equations 6.33 and 6.44 to imrpove 
            the parameter estimates.
        (3) extract K and R from the projection using RQ-decomposition
        */
    }
}
