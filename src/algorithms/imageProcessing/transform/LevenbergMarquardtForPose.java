package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Camera.CameraExtrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraIntrinsicParameters;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import java.util.Arrays;
import no.uib.cipr.matrix.NotConvergedException;

/**
 * iterative non-linear optimization using Levenberg-Marquardt algorithm
 * to minimize the re-projection error of perspective projection.
 * 
 * L-M is guaranteed to converge to eventually find an improvement, 
 * because an update with a sufficiently small magnitude and a negative scalar 
 * product with the gradient is guaranteed to do so.
 * 
 * TODO: consider implementing the Szeliski 2010 chapter 6 equations (6.44)-(6.47)
 * 
 * NOTE: consider implementing the version of L-M by Barfoot et al. 2010
 * which utilizes the sparseness of the block-diagonal structure in the
 * hessian approximation by Engles, Stewenius, & Nister "Bundle Adjustment Rules".
 * 
 * @author nichole
 */
public class LevenbergMarquardtForPose {
    
    /**
     * given initial camera calibration for a single camera 
     * and initial estimates for extrinsic parameters for each image, use the 
     * Levenberg-Marquardt algorithm to improve the rotation and translation 
     * estimates by minimizing the re-projection error using feature 
     * measurements in world coordinates and in each image.
     * <pre>
     References:
     
     Gordon Wetzstein lecture notes, Stanford University, EE 267 Virtual Reality, 
     "Course Notes: 6-DOF Pose Tracking with the VRduino",
       https://stanford.edu/class/ee267/notes/ee267_notes_tracking.pdf
     
     Danping Zou lecture notes, Shanghai Jiao Tong University,
     EE382-Visual localization & Perception, “Lecture 08- Nonlinear least square & RANSAC”
     http://drone.sjtu.edu.cn/dpzou/teaching/course/lecture07-08-nonlinear_least_square_ransac.pdf

     Szeliski 2010, "Computer Vision: Algorithms and Applications", Chapter 6
     
     NOTE: this algorithm is vulnerable to singularities induced by the
     changes of the rotation angles.
     TODO: implement the sparse block version algorithm of Barfoot et al. 2010
     in this same class.
     * </pre>
     * @param imageC
     * @param worldC
     * @param kIntr
     * @param kExtr
     * @param kRadial
     * @param nMaxIter
     * @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
        f(r) = 1 +k1*r^2 + k2*r^4 if true,
        else use model #3 f(r) = 1 +k1*r + k2*r^2.
        note that if rCoeffs is null or empty, no radial distortion is removed.
     * @throws Exception if there is an error in use of MPSolver during the
     * removal of radial distortion, a generic exception is thrown with the
     * error message from the MPSolve documentation.
     * @return 
     * @throws no.uib.cipr.matrix.NotConvergedException 
     */
    public static CameraExtrinsicParameters solveForPose(double[][] imageC, double[][] worldC, 
        CameraIntrinsicParameters kIntr, CameraExtrinsicParameters kExtr, 
        double[] kRadial, final int nMaxIter, boolean useR2R4) 
        throws NotConvergedException, Exception {
        
        // number of features
        int n = imageC[0].length;
        
        if (n < 4) {
            throw new IllegalArgumentException("imageC[0].length must be at least 4");
        }
        if (worldC[0].length != n) {
            throw new IllegalArgumentException("imageC[0].length must equal worldC[0].length");
        }
        
        double[] b = new double[2*n];
        final double[][] xn = Camera.pixelToCameraCoordinates(imageC, kIntr, 
            kRadial, useR2R4);
        for (int i = 0; i < n; ++i) {
            xn[0][i] /= xn[2][i];
            xn[1][i] /= xn[2][i];
            b[2*i] = xn[0][i];
            b[2*i + 1] = xn[1][i];
        }
        
        //TODO: consider adding constraints suggestded in Seliski 2010:
        // u_0 and v_0 are close to half the image lengths and widths, respectively.
        // the angle between 2 image axes is close to 90.
        // the focal lengths along both axes are greater than 0.
        
        // extract pose as (theta_x, theta_y, theta_z, t_x, t_y, t_z)
        double[][] r0 = kExtr.getRotation();
        double[][] r = MatrixUtil.copy(r0);
        double[] thetas = Rotation.extractRotation(r);
        double[] t = kExtr.getTranslation();
        
        // equation (19).  size is 1 X 9
        double[] h = new double[]{r[0][0], r[0][1], t[0], 
            r[1][0], r[1][1], t[1], -r[2][0], -r[2][1], -t[2]};
      
        // equation 20.  length is 2*N
        double[] fgp;
        
        // length is 2*N
        double[] bMinusFGP;
        
        // size is (2N) X 6
        double[][] j = calculateJ(worldC, h, thetas);
        
        // size is 6 X (2N)
        double[][] jT = MatrixUtil.transpose(j);
        // size is 6 X 6
        double[][] jTJ = MatrixUtil.multiply(jT, j);
        
        double lambda = maxDiag(jTJ);
        double lambdaF = 2;
        
        // the residual, that is, the evaulation of the objective which is the 
        //   re-projection error.
        double f = Double.POSITIVE_INFINITY;
        double fPrev;
        
        // deltaPLM is the array of length 6 holding the steps of change for theta and translation.
        double[] deltaPLM;
        // deltaTheta is used to extract the 1st 3 elements of deltaPM 
        double[] deltaTheta = new double[3];
        // deltaT is used to extract the last 3 elements of deltaPLM
        double[] deltaT = new double[3];
                
        double eps = 1.e-5;
        
        // gain ratio:
        //gain = (f(p + delta p LM) - f(p)) / ell(delta p LM)
        //     where ell(delta p LM) is (delta p LM)^T * (lambda * (delta p LM)) + J^T * ( b - fgp))
        //gain = (f - fPrev) / (delta p LM)^T * (lambda * (delta p LM)) + J^T * ( b - fgp))
        double gainRatio;
        
        double[] stepLengthCheck;
        double[] gradientCheck;
        final double tolP = 1.e-3;
        final double tolG = 1.e-3;
                    
        /*
        init: f=infinity
              lambda = max(trace(J^T*J) which means calc J for initial conditions
        loop:
            set fPrev = f
            calc fgp = map(worldC, h);
            set f = b-fgp
            calc j from (worldC, h, thetas);
            calc deltaPM from (j, lambda, b-fgp)
            check stopping conditions: jT-bMinusFGP small or deltaPM small
            update thetas, t from (deltaPM)
            update h from (thetas, t)
            adjust lambda from (f, fPrev, deltaPM, b-fgp)
        */
        
        int doUpdate = 0;
        int nIter = 0;
        while (nIter < nMaxIter) {
            
            nIter++;
                   
            fPrev = f;
            // eqn (20) of Wetzstein.  length is 2*N
            // project the world coordinates to the camera coord frame, using R and T in h:
            fgp = map(worldC, h);
            // in camera reference frame, subtract the projected world points from the observations:
            bMinusFGP = MatrixUtil.subtract(b, fgp);            
            f = evaluateObjective(bMinusFGP);
            
            // ===== calculate step ========
            // eqns (24-34) of Wetzstein
            j = calculateJ(worldC, h, thetas); //(2N) X 6
            jT = MatrixUtil.transpose(j);
            jTJ = MatrixUtil.multiply(jT, j);
           
            //gradient is the local direction of steepest ascent
            
            gradientCheck = MatrixUtil.multiplyMatrixByColumnVector(jT, bMinusFGP);
            deltaPLM = calculateDeltaPLMSzeliski(jTJ, lambda, gradientCheck);
            
            System.out.printf("j^T*(b-fgp)=%s\n", FormatArray.toString(gradientCheck, "%.3e"));
            System.out.printf("deltaP=%s\n", FormatArray.toString(deltaPLM, "%.3e"));
        
            // ======= stopping conditions ============
            //   step length vanishes:  deltaPLM --> 0
            //   gradient of f(x) vanishes: -J^T * (b - fgp) --> 0
            stepLengthCheck = deltaPLM;
            //MatrixUtil.multiply(gradientCheck, -1.);            
            if (isNegligible(stepLengthCheck, tolP) || !isNegligible(gradientCheck, tolG)) {
                break;
            }
            
             //TODO: consider whether to accept or reject step?
             // comments from: https://arxiv.org/pdf/1201.5885
             //           Transtruma & Sethna 2012
            // the qualitative effect of the damping term is to modify the 
            // eigenvalues of the matrix JT J + λDT D to be at least λ
            
            // ======= revise parameters =======
             
            // ====== accept or reject changes and change lambda ======
            if (nIter == 0) {
                doUpdate = 1;   
            } else {
                // gain ratio:
                //    gain = (f(p + delta p LM) - f(p)) / ell(delta p LM)
                //         where ell(delta p LM) is (delta p LM)^T * (lambda * (delta p LM)) + J^T * ( b - fgp))
                //    gain = (f - fPrev) / (delta p LM)^T * (lambda * (delta p LM)) + J^T * ( b - fgp))
                gainRatio = calculateGainRatio(f/2., fPrev/2., deltaPLM, lambda, 
                    gradientCheck, eps);
                System.out.printf("lambda=%.4e, gainRatio=%.4e\n", lambda, gainRatio);
                if (gainRatio > 0) {
                    doUpdate = 1;
                    // near the minimimum, which is good.
                    // decrease lambda
                    lambda /= lambdaF;
                } else {
                    doUpdate = 0;
                    // increase lambda to reduce step length and get closer to 
                    // steepest descent direction
                    lambda *= lambdaF;
                }
                System.out.printf("new lambda=%.4e\n", lambda);
            }
            if (doUpdate == 1) {
                 // add deltaPLM to p, where p is (theta_x, theta_y, theta_z, t_x, t_y, t_z)
                 
                System.arraycopy(deltaPLM, 0, deltaTheta, 0, 3);
                System.arraycopy(deltaPLM, 3, deltaT, 0, 3);
                
                updateT(t, deltaT);
                
                //NOTE: potential singularity in this method:
                updateRTheta(r, thetas, deltaTheta);
                             
                //updateHWithThetaT(h, thetas, t);
                
                updateHWithRT(h, r, t);
            }            
        }
        
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
        
        // ===== create rotation matrix from thetas
        CameraExtrinsicParameters extrinsic = new CameraExtrinsicParameters();
        extrinsic.setRotation(Rotation.createEulerRotationMatrix(thetas[0], 
            thetas[1], thetas[2]));
        extrinsic.setTranslation(t);
        
        return extrinsic;
    }
    
    /**
     * given initial camera calibration and extrinsic parameters, use the 
     * Levenberg-Marquardt
     * algorithm to improve the rotation and translation estimates
     * by minimizing the re-projection error.
     * This algorithm assumes one camera calibration for all images, and
     * individual poses for each image.
     * This method exploits some of the properties of sparse matrices to
     * reduce computation following Barfoot et al. 2010.
     <pre>
     References:
     
     T. Barfoot, et al., "Pose estimation using linearized rotations and 
     quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049
     </pre>
     @param coordsI  holds the image coordinates in pixels of
               features present in all images ordered in the same
               manner and paired with features in coordsW.
               It is a 2 dimensional double array of format
               3 X (N*n) where N is the number of images.
               the first row is the x coordinates, the second row
               is the y coordinates, and the third row is "1"'s.
               The columns hold each image in order and within each image's
               columns are the features presented in the same order in each image.
               In Table 1 of Ma, Chen, & Moore 2003 "Camera Calibration"
               these are the (u_d, v_d) pairs.
     * @param kIntr
     * @param kExtr
     * @param kRadial
     * @param nMaxIter
     * @param coordsW holds the world coordinates of features, ordered
               by the same features in the images.
               the first row is the X coordinates, the second row
               is the Y coordinates, and the third row is 1's 
               (Z_w = 0, the scale factor is lost in the homography).
               It is a 2 dimensional double array of format
               3 X n
     * @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
        f(r) = 1 +k1*r^2 + k2*r^4 if true,
        else use model #3 f(r) = 1 +k1*r + k2*r^2.
        note that if rCoeffs is null or empty, no radial distortion is removed.
     * @throws Exception if there is an error in use of MPSolver during the
     * removal of radial distortion, a generic exception is thrown with the
     * error message from the MPSolve documentation.
     * @return 
     * @throws no.uib.cipr.matrix.NotConvergedException 
     */
    public static CameraExtrinsicParameters solveForPoseUsingSparse(
        double[][] imageC, double[][] worldC, 
        CameraIntrinsicParameters kIntr, CameraExtrinsicParameters kExtr, 
        double[] kRadial, final int nMaxIter, boolean useR2R4) 
        throws NotConvergedException, Exception {
        
        int n = imageC[0].length;
        
        if (n < 4) {
            throw new IllegalArgumentException("imageC[0].length must be at least 4");
        }
        if (worldC[0].length != n) {
            throw new IllegalArgumentException("imageC[0].length must equal worldC[0].length");
        }
        
        /*
        The algortihm operates on pairs of reference frames u and v 
        (it was desgined for stereoimaging from similar cameras), so will
        either use the pattern of reference frame u always being that of image 0
        or reference frame always being sequential pairs.
        
        Regarding the use of applying it to pairs of images rather than 
        looking at the sizes of structures from forming them using all images
        at once:
        
        A quote from Engels, Stewénius, & Nistér, 2006, "Bundle Adjustment Rules"
        <pre>
            However, several researchers have noted (Fitzgibbon and 
            Zisserman, 1998, Nistér, 2001, Pollefeys, 1999) that in the 
            application of camera tracking, performing bundle adjustment each 
            time a new frame has been added to the estimation can prevent the 
            tracking process from failing over time. 
            Thus, bundle adjustment can over time in a sequential estimation 
            process have a much more dramatic impact than mere accuracy 
            improvement, since it improves the initialization for future 
            estimates, which can ultimately enable success in cases that would 
            otherwise miserably fail.
        </pre>
        */
        see https://cseweb.ucsd.edu/classes/fa04/cse252c/manmohan1.pdf
        and add notes here (slide 10 and on)
        
        see http://users.ics.forth.gr/~lourakis/sba/PRCV_colloq.pdf
        NOTE that Engels et al. use J_f=[J_P J_C] which is reverse order
        to that in  Barfoot et al. 2010, hence the different ordering in 
        the final Schur complement reduced equation.
        
        Engels [ H_PP    H_PC ] [ dP ] = [ b_P ]
               [ H_PC^T  H_CC ] [ dC ]   [ b_C ]
        
        From Barfoot et al. 2010:        
        The objective function J on eqn (75a) is the sum of the
        squares of the errors in reprojection in u and in v,
        so one should be able to minimize J over more than 2 images similarly
        as errors add in quadrature.
                
        --> The pose state variables, dx_1, can be determined directly from (82) 
              [A_11 - A_12*A_22^-1*A_12^T] * dx_1 = [b1 - A_12*A_22^-1*b2]
        
              (where dx_1 is [dtransl drot])
        
              A * x = b
                  NOTE: A is called the RCM (Reduced Camea Matrix)
        
              Can be solved without inverting A since it is a sparse matrix!
              https://www.cvg.ethz.ch/teaching/3dvision/2014/slides/class06eth14.pdf
              double[] x = LinearEquations.solveXFromLUDecomposition(a, b);
              Sparse Matrix factorization:
                  (1) LU factorization: A = L*U where original equation is A*x = b
                      (a) use LU decomposition to solve L,U = luDecomp(A)
                      (b) let c = U*x and rewrite A*x=b as L*c = b.
                            then solve c from forward elimination
                      (c) solve x from back substitution in c = U*x
                  (2) QR factorization: A = Q*R
                  (3) Cholesky factorization:A = L*L^T
                  Some details for the above 3 are in Triggs et al. 1999/2000 Appendix B and near page 23
                  ("Bundle Ajustment – A Modern Synthesis",
                  Springer-Verlag, pp.298–372, 2000, 
                  Lecture Notes in Computer Science (LNCS), 
                  10.1007/3-540-44480-7_21. inria-00590128)
                  https://hal.inria.fr/file/index/docid/590128/filename/Triggs-va99.pdf
              Or Iterative methods:
                  (1) Conjugate Gradient
                  (2) Gauss-Seidel
        
              consider Google Ceres  https://code.google.com/p/ceres-solver/
              or Lourakis SBA
        
             more suggestions in http://users.ics.forth.gr/~lourakis/sba/PRCV_colloq.pdf
               (1) Store as dense, decompose with ordinary linear algebra ◦
                    M. Lourakis, A. Argyros: SBA: A Software Package For Generic 
                       Sparse Bundle Adjustment. ACM Trans. Math. Softw. 36(1): (2009) ◦ 
                    C. Engels, H. Stewenius, D. Nister: Bundle Adjustment Rules. 
                       Photogrammetric Computer Vision (PCV), 2006. 
               (2) Store as sparse, factorize with sparse direct solvers ◦ 
                    K. Konolige: Sparse Sparse Bundle Adjustment. BMVC 2010: 1-11
               (3) Store as sparse, use conjugate gradient methods memory efficient, 
                    iterative, precoditioners necessary! ◦ 
                    S. Agarwal, N. Snavely, S.M. Seitz, R. Szeliski: 
                       Bundle Adjustment in the Large. ECCV (2) 2010: 29-42 ◦ 
                    M. Byrod, K. Astrom: Conjugate Gradient Bundle Adjustment. 
                       ECCV (2) 2010: 114-127 
               (4) Avoid storing altogether ◦ 
                    C. Wu, S. Agarwal, B. Curless, S.M. Seitz: Multicore Bundle 
                       Adjustment. CVPR 2011: 30 57-3064 ◦ 
                    M. Lourakis: Sparse Non-linear Least Squares Optimization 
                       for Geometric Vision. ECCV (2) 2010: 43-56
        
              and computer vision library in java http://boofcv.org/index.php?title=Example_Sparse_Bundle_Adjustment
              there is also a java binding for SBA http://seinturier.fr/jorigin/jsba.html
        
        --> the point position variables, dx_2, found inexpensively through 
            backsubstitution. 
              [A_12^T]*[dx_1] + [A_22]*[dx_2] = b2
              [A_22]*[dx_2] = b2 - [A_12^T]*[dx_1]
              A_22 is invertible
              dx_2 = (A_22)^-1 * ( b2 - [A_12^T]*[dx_1] )
        --> After solving for dx we update r using (76a), 
            p_j using (76c), 
            and C using (31). 
        --> We iterate through the entire procedure until dx is sufficiently 
            small. Once converged, the covariance of the state estimate is 
            provided by A^-1
        */
        
        //TODO: consider adding constraints suggestded in Seliski 2010:
        // u_0 and v_0 are close to half the image lengths and widths, respectively.
        // the angle between 2 image axes is close to 90.
        // the focal lengths along both axes are greater than 0.
        
        // extract pose as (theta_x, theta_y, theta_z, t_x, t_y, t_z)
        double[][] r0 = kExtr.getRotation();
        double[][] r = MatrixUtil.copy(r0);
        double[] thetas = Rotation.extractRotation(r);
        double[] t = kExtr.getTranslation();
        
        /*
        NOTE: the Barfoot et al. 2010 notation here is different from the
        Ma, Chen, & Moore 2003, 2004 notation.
        
        u, v: transform coordinates to camera reference frame
        
        p: project the world coordinates to the u camera reference frame
           using the current rotation and translation estimates
        */
        
        
                
        throw new UnsupportedOperationException("not yet finished");
    }
    
    /**
     * map the homography matrix to the projected 2D point coordinates of the 
     * world reference points eqn (20).
     * @param worldC
     * @param h
     * @return 
     */
    static double[] map(double[][] worldC, double[] h) {
        
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
        int i, j;
        double r;
        for (i = 0; i < n; ++i) {
            sum += (bMinusFGP[i] * bMinusFGP[i]);
        }
        
        return sum;
    }
    
    /**
     * 
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
     * 
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
    
    /**
     * 
     * @param thetas
     * @return 9X6 matrix
     */
    static double[][] calculateJG(double[] thetas) {
                        
        // 9 X 9
        double[][] jG = new double[9][6];//MatrixUtil.zeros(9, 6);
        int i, j;
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

    private static double maxDiag(double[][] a) {
        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < a.length; ++i) {
            if (a[i][i] > max) {
                max = a[i][i];
            }
        }
        return max;
    }

    /**
     * following Szeliski 2010, Chap 6, eqn (6.18)
     * @param jTJ J^T * J.  size is 6X6.
     * @param jT.  size is 6X(2*nFeatures)
     * @param lambda
     * @param jTBFG J^T * (B-F(g(p))). size is 6X1
     * @return 
     */
    private static double[] calculateDeltaPLMSzeliski(double[][] jTJ, 
        double lambda, double[] jTBFG) throws NotConvergedException {
        
        //                 inv(6X6)                       6X2N * 2N
        //delta p = pseudoInv(J^T*J + lambda*diag(J^T*J)) * J^T*BFG
        
        int i, j;
        // J^T J + λ diag(J^TJ)     
        // 6 X 6
        double[][] a = MatrixUtil.copy(jTJ);
        for (i = 0; i < 6; ++i) {
            a[i][i] += (lambda*(jTJ[i][i]));
        }
        double[][] aInv = MatrixUtil.pseudoinverseFullRank(a);
                
        double[] step = MatrixUtil.multiplyMatrixByColumnVector(aInv, jTBFG);
        
        return step;
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
        double[][] a = MatrixUtil.elementwiseAdd(jTJ, identity);
        double[][] aInv = MatrixUtil.pseudoinverseFullRank(a);
        
        double[] step = MatrixUtil.multiplyMatrixByColumnVector(aInv, jTBF);
        
        return step;
    }
    
    /**
     * gain = (f(p + delta p LM) - f(p)) / ell(delta p LM)
             where ell(delta p LM) is (delta p LM)^T * (lambda * (delta p LM)) + J^T * ( b - fgp))
       gain = (f - fPrev) / (delta p LM)^T * (lambda * (delta p LM)) + J^T * ( b - fgp))
     * @param f
     * @param fPrev
     * @param deltaPLM
     * @param lambda
     * @param jTBFG j^T*(b-f(g(p))). size is 6X1
     * @return 
     */
    private static double calculateGainRatio(double f, double fPrev, 
        double[] deltaPLM, double lambda, double[] jTBFG,
        double eps) {
             
        //      1X6          *            ( 6X1   +   6 X (2N) * (2NX1) )
        //      1X6                        6X1 
        //  1X1
        //(delta p LM)^T * (lambda * (delta p LM) + J^T * (b - fgp))
        double[] pt1 = Arrays.copyOf(deltaPLM, deltaPLM.length);
        MatrixUtil.multiply(pt1, lambda);
        
        double ell = MatrixUtil.innerProduct(deltaPLM, jTBFG);
        
        if (Math.abs(ell) < eps) {
            return Double.POSITIVE_INFINITY;
        }
        
        double gain = (f - fPrev)/ell;
        
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

    private static void updateHWithRT(double[] h, double[][] r, double[] t) {
        //h = new double[]{r[0][0], r[0][1], t[0], 
       //     r[1][0], r[1][1], t[1], -r[2][0], -r[2][1], -t[2]};
        h[0] = r[0][0];
        h[1] = r[0][1];
        h[2] = t[0];
        h[3] = r[1][0];
        h[4] = r[1][1];
        h[5] = t[1];
        h[6] = -r[2][0];
        h[7] = -r[2][1];
        h[8] = -t[2];
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
     * update t by deltaT
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
    
    /**
     * update rotation matrix r and the theta vector with deltaTheta.
     * NOTE: this method has a potential error from singularity from 90 degree 
     * angles (near zenith and nadir, for example).
     * @param thetas input and output array holding euler rotation angles 
     *    theta_x, theta_y, theta_
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
        
        // potential problem changing an axis when a theta=90
        double[][] sTheta = Rotation.sTheta(thetas);
        double[] deltaPhi = MatrixUtil.multiplyMatrixByColumnVector(sTheta, deltaTheta);
        
        // The Barfoot et al. paper uses
        //    rotation matrix formed from rZ * rY * rX (yaw, pitch, and roll)
        //    which is the same convention used by Wetzstein
        
        double[][] rDeltaPhi = Rotation.calculateRotationZYX(deltaPhi);
        //double[][] rDeltaPhi = Rotation.calculateRotationXYZ(deltaPhi);
        
        double[][] r2 = MatrixUtil.multiply(rDeltaPhi, r);
        
        // update r
        for (int i = 0; i < 3; ++i) {
            System.arraycopy(r2, 0, r, 0, 3);
        }
        
        // T. Barfoot, et al. 2012, 
        // Pose estimation using linearized rotations and quaternion algebra, 
        // Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049
        // eqn (28) for updating thetas:
        // theta = previous theta + sTheta^-1 * deltaPhi
        
        // rotation matrix: inverse is the transpose of rotation matrix.
        double[] sTInvDP = MatrixUtil.multiplyMatrixByColumnVector(
            MatrixUtil.transpose(sTheta), deltaPhi);
        
        double[] thetas2 = MatrixUtil.add(thetas, sTInvDP);
        
        for (int i = 0; i < 3; ++i) {
            System.arraycopy(thetas2, 0, thetas, 0, 3);
        }
        
    }

}
