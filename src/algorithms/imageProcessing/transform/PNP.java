package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Camera.CameraExtrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraIntrinsicParameters;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import no.uib.cipr.matrix.NotConvergedException;

/**
 * NOT READY FOR USE.
 * iterative non-linear optimization using Levenberg-Marquardt algorithm
 * to minimize the re-projection error of perspective projection
 * in "perspective-n-point".
 * 
 * L-M is guaranteed to converge to eventually find an improvement, 
 * because an update with a sufficiently small magnitude and a negative scalar 
 * product with the gradient is guaranteed to do so.
 * 
 * TODO: consider implementing the Szeliski 2010 chapter 6 equations (6.44)-(6.47)
 * 
 * NOTE: this implementation uses euler rotation angles and singularity safe
 * updates, but future versions could consider using quaternions.
 * This from Szeliski 2010:
   ..."Quaternions, on the other hand, are better if you want to keep track of 
   a smoothly moving camera, since there are no discontinuities in the 
   representation. It is also easier to interpolate between rotations and to 
   chain rigid transformations (Murray, Li, and Sastry 1994; Bregler and Malik 1998).
   My usual preference is to use quaternions, but to update their estimates 
   using an incremental rotation, as described in Section 6.2.2."
 * @author nichole
 */
public class PNP {
    
    /**
     * NOT YET TESTED.
     * given the fixed intrinsic camera calibration for the single camera
     * used in all images,
     * and given the initial estimates
     * of extrinsic camera parameters for each image, use the 
     * Levenberg-Marquardt algorithm to improve the extrinsic camera parameters rotation and translation 
     * by minimizing the re-projection error using feature 
     * measurements and their world coordinates.
     * NOTE: this is a dense Hessian matrix solver, not the sparse Hessian
     * matrix bundle adjustment which is in BundleAdjustment.java.
     * <pre>
     References:
     
     Gordon Wetzstein lecture notes, Stanford University, EE 267 Virtual Reality, 
     "Course Notes: 6-DOF Pose Tracking with the VRduino",
       https://stanford.edu/class/ee267/notes/ee267_notes_tracking.pdf
     
     Danping Zou lecture notes, Shanghai Jiao Tong University,
     EE382-Visual localization & Perception, “Lecture 08- Nonlinear least square & RANSAC”
     http://drone.sjtu.edu.cn/dpzou/teaching/course/lecture07-08-nonlinear_least_square_ransac.pdf

     Szeliski 2010, "Computer Vision: Algorithms and Applications", Chapter 6
     
     * </pre>
     * @param imageC observed features in one image followed by the same observed
     * features in the next image, etc. format of double array with
     * first dimension length 3 (for x, y, z) and the second dimension length
     * of nFeatures * mImages. 
     * @param worldC world scene features in format of double array with
     * first dimension length 3 (for x, y, z) and the second dimension length
     * of nFeatures.
     * @param kIntr
     * @param kExtrs
     * @param nMaxIter
     *
     * @throws Exception if there is an error in use of MPSolver during the
     * removal of radial distortion, a generic exception is thrown with the
     * error message from the MPSolve documentation.
     * @return 
     * @throws no.uib.cipr.matrix.NotConvergedException 
     */
    public static List<CameraExtrinsicParameters> solveForPose(double[][] imageC, double[][] worldC, 
        CameraIntrinsicParameters kIntr, 
        List<CameraExtrinsicParameters> kExtrs, final int nMaxIter)
        throws NotConvergedException, Exception {
        
        if (imageC.length != 3) {
            throw new IllegalArgumentException("imageC.length must be 3 for the x,y,z dimensions");
        }
        if (worldC.length != 3) {
            throw new IllegalArgumentException("worldC.length must be 3 for the x,y,z dimensions");
        }
        
        int nFeatures = worldC[0].length;
        
        if (nFeatures < 4) {
            throw new IllegalArgumentException("worldC[0].length must be at least 4");
        }
        if ((imageC[0].length % nFeatures) != 0) {
            throw new IllegalArgumentException("imageC[0].length must be an integer multiple of worldC[0].length");
        }
        
        int nImages = imageC[0].length/nFeatures;
        
        if (kExtrs.size() != nImages) {
            throw new IllegalArgumentException("from imageC[0].length, "
               + "nImages is expected to be " + nImages
               + " therefore kExtrs.size() should be " + nImages + " also.");
        }

        List<CameraExtrinsicParameters> refined = new ArrayList<CameraExtrinsicParameters>();
        
        double[][] imageCI = MatrixUtil.zeros(3, nFeatures);
        
        CameraExtrinsicParameters r;
        
        for (int i = 0; i < nImages; ++i) {
            
            MatrixUtil.copySubMatrix(imageC, 0, 2, i*nFeatures, i*nFeatures+nFeatures-1, imageCI);
            
            r = solveForPose(imageCI, worldC, kIntr, kExtrs.get(i), nMaxIter);
            
            refined.add(r);
        }
        return refined;
    }
    
    /**
     * NOT YET TESTED.
     * given the fixed intrinsic camera calibration and the initial estimates
     * of extrinsic camera parameters as lists of the same size, each item being
     * for a single image, use the 
     * Levenberg-Marquardt algorithm to improve the extrinsic camera parameters rotation and translation 
     * by minimizing the re-projection error using feature 
     * measurements and their world coordinates.
     * NOTE: this is a dense Hessian matrix solver, not the sparse Hessian
     * matrix bundle adjustment which is in BundleAdjustment.java.
     * <pre>
     References:
     
     Gordon Wetzstein lecture notes, Stanford University, EE 267 Virtual Reality, 
     "Course Notes: 6-DOF Pose Tracking with the VRduino",
       https://stanford.edu/class/ee267/notes/ee267_notes_tracking.pdf
     
     Danping Zou lecture notes, Shanghai Jiao Tong University,
     EE382-Visual localization & Perception, “Lecture 08- Nonlinear least square & RANSAC”
     http://drone.sjtu.edu.cn/dpzou/teaching/course/lecture07-08-nonlinear_least_square_ransac.pdf

     Szeliski 2010, "Computer Vision: Algorithms and Applications", Chapter 6
     
     * </pre>
     * @param imageC
     * @param worldC
     * @param kIntrs
     * @param kExtrs
     * @param nMaxIter
     * @throws Exception if there is an error in use of MPSolver during the
     * removal of radial distortion, a generic exception is thrown with the
     * error message from the MPSolve documentation.
     * @return 
     * @throws no.uib.cipr.matrix.NotConvergedException 
     */
    public static List<CameraExtrinsicParameters> solveForPose(double[][] imageC, double[][] worldC, 
        List<CameraIntrinsicParameters> kIntrs, 
        List<CameraExtrinsicParameters> kExtrs, final int nMaxIter)
        throws NotConvergedException, Exception {
        
        if (kIntrs.size() != kExtrs.size()) {
            throw new IllegalArgumentException("the sizes of Kintrs and kExtrs must be the same");
        }
        
        if (imageC.length != 3) {
            throw new IllegalArgumentException("imageC.length must be 3 for the x,y,z dimensions");
        }
        if (worldC.length != 3) {
            throw new IllegalArgumentException("worldC.length must be 3 for the x,y,z dimensions");
        }
        
        int nFeatures = worldC[0].length;
        
        if (nFeatures < 4) {
            throw new IllegalArgumentException("worldC[0].length must be at least 4");
        }
        if ((imageC[0].length % nFeatures) != 0) {
            throw new IllegalArgumentException("imageC[0].length must be an integer multiple of worldC[0].length");
        }
        
        int nImages = imageC[0].length/nFeatures;
        
        if (kExtrs.size() != nImages) {
            throw new IllegalArgumentException("from imageC[0].length, "
               + "nImages is expected to be " + nImages
               + " therefore kExtrs.size() should be " + nImages + " also.");
        }
        if (kIntrs.size() != nImages) {
            throw new IllegalArgumentException("from imageC[0].length, "
               + "nImages is expected to be " + nImages
               + " therefore kIntrs.size() should be " + nImages + " also.");
        }

        List<CameraExtrinsicParameters> refined = new ArrayList<CameraExtrinsicParameters>();
        
        double[][] imageCI = MatrixUtil.zeros(3, nFeatures);
        
        CameraExtrinsicParameters r;
        
        for (int i = 0; i < nImages; ++i) {
            
            MatrixUtil.copySubMatrix(imageC, 0, 2, i*nFeatures, i*nFeatures+nFeatures-1, imageCI);
            
            r = solveForPose(imageCI, worldC, kIntrs.get(i), kExtrs.get(i), nMaxIter);
            
            refined.add(r);
        }
        return refined;
    }
    
    /**
     * NOT YET TESTED.
     * given the camera intrinsic calibration  
     * and initial estimates for extrinsic camera parameters for an image, use the 
     * Levenberg-Marquardt algorithm to improve the extrinsic camera parameters rotation and translation 
     * estimates by minimizing the re-projection error using feature 
     * measurements and their world coordinates.
     * NOTE: this is a dense Hessian matrix solver, not the sparse Hessian
     * matrix bundle adjustment which is in BundleAdjustment.java.
     * <pre>
     References:
     
     Gordon Wetzstein lecture notes, Stanford University, EE 267 Virtual Reality, 
     "Course Notes: 6-DOF Pose Tracking with the VRduino",
       https://stanford.edu/class/ee267/notes/ee267_notes_tracking.pdf
     
     Danping Zou lecture notes, Shanghai Jiao Tong University,
     EE382-Visual localization & Perception, “Lecture 08- Nonlinear least square & RANSAC”
     http://drone.sjtu.edu.cn/dpzou/teaching/course/lecture07-08-nonlinear_least_square_ransac.pdf

     Szeliski 2010, "Computer Vision: Algorithms and Applications", Chapter 6
     
     * </pre>
     * @param imageC
     * @param worldC
     * @param kIntr
     * @param kExtr
     * @param nMaxIter
     * @throws Exception if there is an error in use of MPSolver during the
     * removal of radial distortion, a generic exception is thrown with the
     * error message from the MPSolve documentation.
     * @return 
     * @throws no.uib.cipr.matrix.NotConvergedException 
     */
    public static CameraExtrinsicParameters solveForPose(double[][] imageC, double[][] worldC, 
        CameraIntrinsicParameters kIntr, 
        CameraExtrinsicParameters kExtr, final int nMaxIter)
        throws NotConvergedException, Exception {

        boolean useBouguetsRodrigues = false;
        if (true) {
            return CameraPose.bouguetPoseRefine(kExtr, kIntr, imageC, worldC, useBouguetsRodrigues);
        }

        if (imageC.length != 3) {
            throw new IllegalArgumentException("imageC.length must be 3 for the x,y,z dimensions");
        }
        if (worldC.length != 3) {
            throw new IllegalArgumentException("worldC.length must be 3 for the x,y,z dimensions");
        }
        
        // number of features
        int n = worldC[0].length;
        
        if (n < 4) {
            throw new IllegalArgumentException("worldC[0].length must be at least 4");
        }
        if (imageC[0].length != n) {
            throw new IllegalArgumentException("imageC[0].length must equal worldC[0].length");
        }
        
        double[] b = new double[2*n];
        final double[][] xn = Camera.pixelToCameraCoordinates(imageC, kIntr);
        for (int i = 0; i < n; ++i) {
            xn[0][i] /= xn[2][i];
            xn[1][i] /= xn[2][i];
            //TODO: follow up on signs from pixelToCameraCoordinates.
            b[2*i] = xn[0][i];
            b[2*i + 1] = xn[1][i];
        }
                
        //TODO: consider adding constraints suggested in Seliski 2010:
        // u_0 and v_0 are close to half the image lengths and widths, respectively.
        // the angle between 2 image axes is close to 90.
        // the focal lengths along both axes are greater than 0.
      
        // ==========  initialize ===================
        
        // extract pose as (theta_x, theta_y, theta_z, t_x, t_y, t_z)
        double[][] rTest = MatrixUtil.copy(kExtr.getRotation());
        double[] thetasTest = Rotation.extractThetaFromZYX(rTest);
        double[] tTest = Arrays.copyOf(kExtr.getTranslation(), 4);
        
        CameraExtrinsicParameters outExtr = new CameraExtrinsicParameters();
        outExtr.setRotation(MatrixUtil.copy(kExtr.getRotation()));
        outExtr.setTranslation(Arrays.copyOf(kExtr.getTranslation(), 4));
        
//TODO: check the signs of the last row.  wetzstein uses a convention of looking down
//the negative z-axis.
        
        // equation (19).  size is 1 X 9
        double[] h = new double[9];
        populateHWithRT(h, rTest, tTest);
        
        // eqn (20) of Wetzstein.  length is 2*N
        // project the world coordinates to the camera coord frame, using R and T in h:
        double[] fgp = transformToCamera(worldC, h);
        
        // in camera reference frame, subtract the projected world points from the observations:
        // length is 2*N
        double[] bMinusFGP = MatrixUtil.subtract(b, fgp);   
                
// the projected x and y have opposite signs than they should      
        System.out.printf("observed in camera ref frame=\n%s\n",
            FormatArray.toString(b, "%.3e"));
        System.out.printf("world feature projected to camera ref frame=\n%s\n",
            FormatArray.toString(fgp, "%.3e"));
        System.out.printf("(observed-projected) in camera ref frame=\n%s\n",
            FormatArray.toString(bMinusFGP, "%.3e"));
        
        // sum the squares to evaluate the re-projection error:
        double fPrev = evaluateObjective(bMinusFGP);
        double fTest = Double.POSITIVE_INFINITY;
        
        // size is (2N) X 6
        double[][] j = calculateJ(worldC, h, thetasTest);
        
        // size is 6 X (2N)
        double[][] jT = MatrixUtil.transpose(j);
        // size is 6 X 6
        double[][] jTJ = MatrixUtil.multiply(jT, j);
        
        double lambda = maxDiag(jTJ);
        System.out.printf("max diag of Hessian lambda=%.7e\n", lambda);

        //factor to raise or lower lambda.  
        //   consider using the eigenvalue spacing of J^T*J (Transtrum & Sethna, "Improvements to the Levenberg-Marquardt algorithm for nonlinear least-squares minimization")
        final double lambdaF = 2;
                
        // deltaP is the array of length 6 holding the steps of change for theta and translation.
        double[] deltaP = new double[6];
        
        // deltaTheta is used to extract the 1st 3 elements of deltaPM 
        double[] deltaTheta = new double[3];
        // deltaT is used to extract the last 3 elements of deltaPLM
        double[] deltaT = new double[3];
                
        double eps = 1E-12;
        
        double gainRatio;
        
        double[] gradient = MatrixUtil.multiplyMatrixByColumnVector(jT, bMinusFGP);
        
        final double tolP = 1.e-3;
        final double tolG = 1.e-3;
        
        int nIter = 0;
        //begin loop of tentative changes
        while ((nIter < nMaxIter) && (Math.abs(fPrev - fTest) >= eps)) {
            
            if (nIter == 0) {
                initDeltaPWithQu(deltaP);
            } else {
                calculateDeltaPLMSzeliski(jTJ, lambda, gradient, deltaP);
            }

            nIter++;
            
            // p is (theta_x, theta_y, theta_z, t_x, t_y, t_z)
            System.arraycopy(deltaP, 0, deltaTheta, 0, 3);
            System.arraycopy(deltaP, 3, deltaT, 0, 3);
            updateT(tTest, deltaT);
            updateRTheta(rTest, thetasTest, deltaTheta);

            populateHWithRT(h, rTest, tTest);
            
            // eqn (20) of Wetzstein.  length is 2*N
            // project the world coordinates to the camera coord frame, using R and T in h:
            fgp = transformToCamera(worldC, h);
            // in camera reference frame, subtract the projected world points from the observations:
            bMinusFGP = MatrixUtil.subtract(b, fgp);   
            // sum the squares:
            fTest = evaluateObjective(bMinusFGP);
            
            System.out.printf("nIter=%d) fPrev=%.11e, fTest=%.11e  diff=%.11e\n", nIter, fPrev, fTest,
                (fPrev-fTest));
            System.out.flush();

            // ===== calculate step ========
            // eqns (24-34) of Wetzstein
            j = calculateJ(worldC, h, thetasTest); // 2NX6
            jT = MatrixUtil.transpose(j); // 6X2N
            jTJ = MatrixUtil.multiply(jT, j); // 6X6
           
            //gradient is the local direction of steepest ascent
            // 6X1
            gradient = MatrixUtil.multiplyMatrixByColumnVector(jT, bMinusFGP);
            
            System.out.printf("\nj^T*(b-fgp)=%s\n", FormatArray.toString(gradient, "%.3e"));
            System.out.printf("deltaP=%s\n", FormatArray.toString(deltaP, "%.3e"));
        
            gainRatio = calculateGainRatio(fTest, fPrev, deltaP, lambda, 
                gradient, eps);
            
            if (gainRatio <= Double.NEGATIVE_INFINITY) {
                break;
            }
            
            System.out.printf("lambda=%.6e\ngainRatio=%.6e\nfPrev=%.11e, fTest=%.11e\n", 
                lambda, gainRatio, fPrev, fTest);
            System.out.flush();

            /*
            for large values of lambda, the update is a very steep descent and
            deltaP is very small.
            If the damping term is small the approach is a nearly linear problem.
            
            NOTE: the damping term is used like a factor in the perturbation
            of a symmetric matrix.  see:
                https://nhigham.com/2021/02/16/diagonally-perturbing-a-symmetric-matrix-to-make-it-positive-definite/
            */
            if (gainRatio > 0) {
                // near the minimimum, which is good.
                // decrease lambda
                lambda /= lambdaF;
                assert(fTest < fPrev);
                fPrev = fTest;
                
                System.out.printf("fPrev updated to fTest: fPrev=%.11e fTest=%.11e\n", fPrev, fTest);
                System.out.flush();
                
                System.out.printf("pointwise diff in rotation:\n%s\n", 
                    FormatArray.toString(MatrixUtil.pointwiseSubtract(outExtr.getRotation(), rTest),
                    "%.11e"));
            
                outExtr.setRotation(MatrixUtil.copy(rTest));
                outExtr.setTranslation(Arrays.copyOf(tTest, 3));
                
                // ======= stopping conditions ============
                //   step length vanishes:  deltaPLM --> 0
                //   gradient of f(x) vanishes: -J^T * (b - fgp) --> 0
                //MatrixUtil.multiply(gradientCheck, -1.);            
                if (isNegligible(deltaP, tolP) || isNegligible(gradient, tolG)) {
                    break;
                }
                
            } else {
                // increase lambda to reduce step length and get closer to 
                // steepest descent direction
                lambda *= lambdaF;
            }
            System.out.printf("new lambda=%.6e\n", lambda);    
        }
        
        return outExtr;
         
        /*
        Initialization: A = J^T*J, lambda=max{a_i_i}
        • Repeat until the step length vanishes, || delta x_LM || --> 0, 
          or the gradient of f(x) vanishes, del f = -J^T * r --> 0
          a) Solve (A + lambda*I)*(delta x) = J^T*r to get delta x_LM
          b) x = x + (delta x_LM)
          c) Adjust the damping parameter by checking the gain ratio
             1. gain > 0 Good approximation, decrease the damping parameter
             2. gain <= 0 Bad approximation, increase the damping parameter
        
        where r_i(x) = y_i - h_i(x) where y_i are the measurements
           and h_i(x) are the projections.
        where delta x_LM = = (J^T*J + lambda*I)^-1 * J^T * r
        where f(x) = (1/2) * || r(x) ||^2
        where rho = (f(x + (delta x_LM)) - f(x)) / (ell(delta x_LM))
        
        where the incremental of the objective function predicted by the 
           linear model is given by
              ell(delta x) = - (delta x)^T * J^T * J * (delta x) - 2*(delta x)^T * J^T * r
        where the incremental predicted by the LM step is computed as
              ell(delta x_LM) = - (delta x)^T * (lambda * (delta x_LM) + J^T*r)
        */
    }
    
    /**
     * use the homography matrix to the project the world reference points
     * to 2D point coordinates in the camera frame.  eqn (20).
     * @param worldC
     * @param h
     * @return 
     */
    static double[] transformToCamera(double[][] worldC, double[] h) {
        
        int n = worldC[0].length;
        
        double[] f = new double[2*n];
        
        int i, j;
        double X, Y, s;
        for (i = 0; i < n; ++i) {
            X = worldC[0][i];
            Y = worldC[1][i];
            s = h[6]*X + h[7]*Y + h[8];
            
            f[2*i] = (h[0]*X + h[1]*Y + h[2])/s;
            f[2*i+1] = ((h[3]*X + h[4]*Y + h[5])/s); 
        }
        
        return f;
    }
    
    /**
     * evaluate the objective (|| b − f (g (p)) ||_2)^2
     * eqn (21)
     * @param bMinusFGP array b - f(g(p))
     * @return 
     */
    static double evaluateObjective(double[] bMinusFGP) {
        
        int n = bMinusFGP.length;
        
        double sum = 0;
        int i;
        for (i = 0; i < n; ++i) {
            sum += (bMinusFGP[i] * bMinusFGP[i]);
        }
        
        return sum;
    }
    
    /**
     * calculate J for extrinsic parameters in solving pose
     <pre> 
     Sect 6.1 of George Wetzstein Stanford Course Notes 
     for EE 267 Virtual Reality, 6-DOF Pose Tracking with the VRduino  
     </pre>
     <pre>
     J = J_f * J_g where p is the parameter vector [thetas, translations]
                   where J_g = dh/dp 
                         where h is the 2-D projection matrix of size 3x3 as
                         2 columns of rotation and last column is translation
                   and 
                   where J_f = df/dh 
                         where f is the world point transformed by the homography h.
         J = df/dp
     </pre>
     * @param worldC
     * @param h
     * @param thetas
     * @return matrix size (2N) X 6
     */
    static double[][] calculateJ(double[][] worldC, double[] h, double[] thetas) {
        //(2N) X 9
        double[][] jF = calculateJF(worldC, h);
        //9 X 6
        double[][] jG = calculateJG(thetas);
        // (2N)X6)
        double[][] j = MatrixUtil.multiply(jF, jG);
        return j;
    }
    
    /**
     <pre>
     J_f = df/dh 
         where f is the world point transformed by the homography h.
     </pre>
     <pre> 
     Sect 6.1 of George Wetzstein Stanford Course Notes 
     for EE 267 Virtual Reality, 6-DOF Pose Tracking with the VRduino  
     </pre>
     * @param worldC
     * @param h
     * @return a (2*n)X9 matrix
     */
    static double[][] calculateJF(double[][] worldC, double[] h) {
        
        int n = worldC[0].length;
        
        if (n < 4) {
            throw new IllegalArgumentException("need at least 4 features in worldC");
        }
        
        //TODO: assert worldC[2][*] = 1 for local device frame
        assert(Math.abs(worldC[2][0] - 1.) < 1e-5);
        
        // (2*n) X 9
        double[][] jF = MatrixUtil.zeros(2*n, 9);
        int i, j;
        double x, y, s1, s2, d, dsq;
        for (i = 0; i < n; ++i) {
            x = worldC[0][i];
            y = worldC[1][i];
            d = h[6] * x + h[7] * y + h[8];
            s1 = h[0] * x + h[1] * y + h[2];
            s2 = h[3] * x + h[4] * y + h[5];
            dsq = s1*s1;
            
            jF[2*i][0] = x/d;
            jF[2*i][1] = y/d;
            jF[2*i][2] = 1./d;
            jF[2*i][6] = -(s1/dsq)*x;
            jF[2*i][7] = -(s1/dsq)*y;
            jF[2*i][8] = -(s1/dsq);
            
            jF[2*i + 1][3] = x/d;
            jF[2*i + 1][4] = y/d;
            jF[2*i + 1][5] = 1./d;
            jF[2*i + 1][6] = -(s2/dsq)*x;
            jF[2*i + 1][7] = -(s2/dsq)*y;
            jF[2*i + 1][8] = -(s2/dsq);
        }
        return jF;
    }
    
    static void calculateJF(double[] worldC, double[] h, double[][] outJF) {
        
        //TODO: assert worldC[2] = 1 for local device frame
        //assert(Math.abs(worldC[2] - 1.) < 1e-5);
        
        assert(outJF.length == 2);
        assert(outJF[0].length == 9);
        
        int i, j;
        double x, y, s1, s2, d, dsq;
        x = worldC[0];
        y = worldC[1];
        d = h[6] * x + h[7] * y + h[8];
        s1 = h[0] * x + h[1] * y + h[2];
        s2 = h[3] * x + h[4] * y + h[5];
        dsq = s1*s1;

        outJF[0][0] = x/d;
        outJF[0][1] = y/d;
        outJF[0][2] = 1./d;
        outJF[0][6] = -(s1/dsq)*x;
        outJF[0][7] = -(s1/dsq)*y;
        outJF[0][8] = -(s1/dsq);

        outJF[1][3] = x/d;
        outJF[1][4] = y/d;
        outJF[1][5] = 1./d;
        outJF[1][6] = -(s2/dsq)*x;
        outJF[1][7] = -(s2/dsq)*y;
        outJF[1][8] = -(s2/dsq);
    }
    
    /**
    <pre>
     J_g = dh/dp 
         where h is the 2-D projection matrix of size 3x3 as
         2 columns of rotation and last column is translation
    </pre>
    <pre> 
     Sect 6.1 of George Wetzstein Stanford Course Notes 
     for EE 267 Virtual Reality, 6-DOF Pose Tracking with the VRduino  
     </pre>
     * @param thetas
     * @return 9X6 matrix
     */
    static double[][] calculateJG(double[] thetas) {
                        
        // 9 X 6
        double[][] jG = new double[9][6];//MatrixUtil.zeros(9, 6);
        double cx, cy, cz, sx, sy, sz;
        cx = Math.cos(thetas[0]);
        cy = Math.cos(thetas[1]);
        cz = Math.cos(thetas[2]);
        sx = Math.sin(thetas[0]);
        sy = Math.sin(thetas[1]);
        sz = Math.sin(thetas[2]);
        
        jG[0] = new double[]{-cx*sy*sz, -sy*cz-sx*cy*sz, -cy*sz-sx*sy*cz,
           0, 0, 0};
        
        jG[1] = new double[]{sx*sz, 0, -cx*cz, 0, 0, 0};
        
        jG[2] = new double[]{0, 0, 0, 1, 0, 0};
        
        jG[3] = new double[]{cx*sy*cz, -sy*sz+sx*cy*cz, cy*cz-sx*sy*sz, 0, 0, 0};
        
        jG[4] = new double[]{-sx*cz, 0, -cx*sz, 0, 0, 0};
        
        jG[5] = new double[]{0, 0, 0, 0, 1, 0};
        
        jG[6] = new double[]{-sx*sy, cx*cy, 0, 0, 0, 0};
        
        jG[7] = new double[]{-cx, 0, 0, 0, 0, 0};
        
        jG[8] = new double[]{0, 0, 0, 0, 0, -1};
        
        return jG;
    }
    
    /**
     * following Szeliski 2010, Chap 6, eqn (6.18)
     * @param jTJ J^T * J.  size is 6X6.
     * @param lambda
     * @param jTBFG J^T * (B-F(g(p))). size is 6X1
     * @param outDeltaP calculated step length
     * @throws no.uib.cipr.matrix.NotConvergedException 
     */
    private static void calculateDeltaPLMSzeliski(double[][] jTJ, 
        double lambda, double[] jTBFG, double[] outDeltaP) throws NotConvergedException {
                
        //                        [6X6]                   * [6X1] = [6X1]
        //delta p = pseudoInv(J^T*J + lambda*diag(J^T*J)) * J^T*BFG
        
        int i, j;
        // J^T J + λ diag(J^TJ)     
        // [6X6]
        double[][] a = MatrixUtil.copy(jTJ);
        for (i = 0; i < 6; ++i) {
            a[i][i] += (lambda*(jTJ[i][i]));
        }
        //[6X6]
        double[][] aInv = MatrixUtil.pseudoinverseRankDeficient(a);
                
        //[6X6] * [6X1] = [6X1]
        MatrixUtil.multiplyMatrixByColumnVector(aInv, jTBFG, outDeltaP);        
    }
        
    /**
     * from lecture notes of Danping Zou
     * http://drone.sjtu.edu.cn/dpzou/teaching/course/lecture07-08-nonlinear_least_square_ransac.pdf
     *                          6X6  * (6 * (2N)) * (2NX1) = 6 X (2N) * (2NX1) = 6X1
     * calculate the step as (J^T*J + lambda*I)^-1 * J^T * (b-f(g(p))
     * @param jTJ 
     * @param lambda
     * @param jTBF
     * @return an array of length 6 
     * @throws NotConvergedException 
     */
    private static double[] calculateDeltaPLM(double[][] jTJ, 
        double lambda, double[] jTBF) throws NotConvergedException {
        
        // (J^T*J + lambda*I)^-1 * J^T * (b-f(g(p))
        double[][] identity = MatrixUtil.createIdentityMatrix(6);
        MatrixUtil.multiply(identity, lambda);
        // 6 X 6
        double[][] a = MatrixUtil.pointwiseAdd(jTJ, identity);
        double[][] aInv = MatrixUtil.pseudoinverseFullColumnRank(a);
        
        double[] step = MatrixUtil.multiplyMatrixByColumnVector(aInv, jTBF);
        
        return step;
    }
    
    /**
     * gain = (f(p) - f(p + delta p)) / ell(delta p)
             where ell(delta p) is (delta p)^T * (lambda * (delta p)) + J^T * ( b - fgp))
       gain = (f - fPrev) / ( (delta p)^T * (lambda * (delta p) + J^T * ( b - fgp)) )
     * @param fNew
     * @param fPrev
     * @param deltaP
     * @param lambda
     * @param jTBFG j^T*(b-f(g(p))). size is 6X1
     * @return 
     */
    private static double calculateGainRatio(double fNew, double fPrev, 
        double[] deltaP, double lambda, double[] jTBFG,
        double eps) {
        
        // (M. Lourakis, A. Argyros: SBA: A Software Package For Generic
        // Sparse Bundle Adjustment. ACM Trans. Math. Softw. 36(1): (2009))
        //  gain ratio = ( fPrev - fNew) / ( deltaP^T * (lambda * deltaP + J^T*fPrev) )
        
        
        //   1X6         *    ([1X1]*[6X1]         + [6X1])           = [1X1]
        //(delta p LM)^T *  (lambda * (delta p LM) + J^T * (b - fgp))
        
             
        //      1X6          *            ( 6X1   +   6 X (2N) * (2NX1) )
        //      1X6                        6X1 
        //  1X1
        //(delta p LM)^T * (lambda * (delta p LM) + J^T * (b - fgp))
        double[] denom = Arrays.copyOf(deltaP, deltaP.length);
        for (int i = 0; i < denom.length; ++i) {
            denom[i] *= lambda;
        }
        
        denom = MatrixUtil.add(denom, jTBFG);
        
        double d = MatrixUtil.innerProduct(deltaP, denom);
        
        if (Math.abs(d) < eps) {
            return Double.NEGATIVE_INFINITY;
        }
        
        double gain = (fPrev - fNew)/d;
        
        return gain;
    }

    private static boolean isNegligible(double[] c, double eps) {
        for (int i = 0; i < c.length; ++i) {
            if (Math.abs(c[i]) > eps) {
                return false;
            }
        }
        return true;
    }

    private static void populateHWithRT(double[] h, double[][] r, double[] t) {
        h[0] = r[0][0];
        h[1] = r[0][1];
        h[2] = t[0];
        h[3] = r[1][0];
        h[4] = r[1][1];
        h[5] = t[1];
        // TODO: review the code for right-handedness vs left-handedness.
        //    Wetzstein convention looks down the negative z-axis.
        //    (see Rotation.java notes). may need transposition...
        //h[6] = -r[2][0];
        //h[7] = -r[2][1];
        //h[8] = -t[2];
        h[6] = r[2][0];
        h[7] = r[2][1];
        h[8] = t[2];
    }
    
    /**
     * 
     * @param h input and output variable homography to be updated with given
     * rotation angles and translations.
     * @param thetas
     * @param t 
     */
    private static void updateHWithThetaT(double[] h, double[] thetas, double[] t) {
        // equation (19).  size is 1 X 9
        
        double cx = Math.cos(thetas[0]);
        double cy = Math.cos(thetas[1]);
        double cz = Math.cos(thetas[2]);
        double sx = Math.sin(thetas[0]);
        double sy = Math.sin(thetas[1]);
        double sz = Math.sin(thetas[2]);
        
        h[0] = cy*cz - sx*sy*sz;
        h[1] = -cx*sz;
        h[2] = t[0];
        h[3] = cy*sz + sx*sy*cz;
        h[4] = cx*cz;
        h[5] = t[1];
        h[6] = cx*sy;
        h[7] = -sx;
        h[9] = -t[2];
    }

    /**
     * update translation t by deltaT
     */
    private static void updateT(double[] t, double[] deltaT) {
        // from Danping Zou lecture notes, Shanghai Jiao Tong University,
        // EE382-Visual localization & Perception, “Lecture 08- Nonlinear least square & RANSAC”
        // http://drone.sjtu.edu.cn/dpzou/teaching/course/lecture07-08-nonlinear_least_square_ransac.pdf

        // parameter perturbations for a vector are 
        //     x + delta x
        int i;
        
        // vector perturbation for translation:
        for (i = 0; i < 3; ++i) {
            t[i] += deltaT[i];
        }
    }
    
    private static void updateT(double[] t, double[] deltaT, double[][] r) {
        // t_local(t + deltaT) = t + r*deltaT
        int i;
        
        double[] rdt = MatrixUtil.multiplyMatrixByColumnVector(r, deltaT);
        
        // vector perturbation for translation:
        for (i = 0; i < 3; ++i) {
            t[i] += rdt[i];
        }
    }
    
    /**
     * update rotation matrix r with deltaTheta.
     * the approach used avoids singularities and the need to restore 
     * the constraint afterwards (i.e., constraint restoration is built in).
     * 
     * @param r
     * @param thetas input and output array holding euler rotation angles 
     *    theta_x, theta_y, theta_
     * @param deltaTheta
     */
    private static void updateRTheta(double[][] r, double[] thetas, double[] deltaTheta) {
        // from Danping Zou lecture notes, Shanghai Jiao Tong University,
        // EE382-Visual localization & Perception, “Lecture 08- Nonlinear least square & RANSAC”
        // http://drone.sjtu.edu.cn/dpzou/teaching/course/lecture07-08-nonlinear_least_square_ransac.pdf
        // parameter perturbations for a Lie group such as rotation are:
        //     R * (I - [delta x]_x) where [delta x]_x is the skew-symetric matrix of delta_x 
        
        // T. Barfoot, et al. 2010, 
        // Pose estimation using linearized rotations and quaternion algebra, 
        // Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049
        // eqn (31) for updating rotation:
        // C(theta) = C(deltaPhi) * previous C
        //     where C is rotation matrix r
        //           calculated as C(theta) = 
        //     and deltaPhi = sTheta * deltaTheta
        
        //double[][] out;// = MatrixUtil.zeros(3, 3);
        double[] qUpdated = Rotation.applySingularitySafeRotationPerturbation(thetas, deltaTheta);
        
        // [4X4]
        double[][] qR = Rotation.createRotationMatrixFromQuaternion4(qUpdated);
        qR = MatrixUtil.transpose(qR);
        
        // rotation is [0:2, 0:2] of qR  
                
        // update in-out variable r
        for (int i = 0; i < 3; ++i) {
            System.arraycopy(qR[i], 0, r[i], 0, 3);
        }

        // ---- update theta ----        
        //extracting theta from the updated rotation would keep the theta
        //    vector consistent with the rotation matrix,
        //    but the value is a little different that updating theta with delta theta
        //    by addition.
        double[] thetaExtracted = Rotation.extractThetaFromZYX(r);
        System.arraycopy(thetaExtracted, 0, thetas, 0, thetas.length);
        
        /*
        for (int i = 0; i < 3; ++i) {
            thetas[i] += deltaTheta[i];
        }*/
        
    }

    /**
     * initialize the first parameter steps with small values suggested
     * by Qu.
     * @param deltaP 
     */
    private static void initDeltaPWithQu(double[] deltaP) {
        /*Qu thesis eqn (3.38)
        
        delta thetas ~ 1e-8
        delta translation ~1e-5
        delta focus ~ 1
        delta kRadial ~ 1e-3
        delta x ~ 1e-8
        */
        int i;
        for (i = 0; i < 3; ++i) {
            // delta theta
            deltaP[i] = 1e-8;
        }
        for (i = 0; i < 3; ++i) {
            // delta translation
            deltaP[3+i] = 1e-5;
        }
        
    }

    private static double maxDiag(double[][] jTJ) {
        double max = Double.NEGATIVE_INFINITY;
        //Arrays.fill(max, Double.NEGATIVE_INFINITY);
        for (int i = 0; i < jTJ.length; ++i) {
            if (jTJ[i][i] > max) {
                max = jTJ[i][i];
            }
        }
        return max;
    }
}
