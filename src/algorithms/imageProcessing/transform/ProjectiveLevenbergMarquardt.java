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
 * TODO: consider implementing the Szeliski 2010 chapter 6 equations (6.44)-(6.47)
 * @author nichole
 */
public class ProjectiveLevenbergMarquardt {
    
    /**
     * given initial camera calibration and extrinsic parameters, use the 
     * Levenberg-Marquardt
     * algorithm to improve those values by minimizing the re-projection error.
     * This version of the algorithm uses details such as elements of the
     * Jacobian provided in the lecture notes 
     * of Gordon Wetzstein at Stanford University,
       EE 267 Virtual Reality, "Course Notes: 6-DOF Pose Tracking with the VRduino",
       https://stanford.edu/class/ee267/notes/ee267_notes_tracking.pdf.
       
       NOTE:the invoking code could wrap this in a larger iteration as feedback to 
       optimize the camera intrinsic parameters and radial distortion parameters
       also.
     * @param imageC
     * @param worldC
     * @param kIntr
     * @param kExtr
     * @param kRadial
     * @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
        f(r) = 1 +k1*r^2 + k2*r^4 if true,
        else use model #3 f(r) = 1 +k1*r + k2*r^2.
        note that if rCoeffs is null or empty, no radial distortion is removed.
     * @return 
     */
    public static CameraExtrinsicParameters solve(double[][] imageC, double[][] worldC, 
        CameraIntrinsicParameters kIntr, CameraExtrinsicParameters kExtr, 
        double[] kRadial, final int nMaxIter, boolean useR2R4) throws NotConvergedException {
        
        int n = imageC[0].length;
        
        if (n < 4) {
            throw new IllegalArgumentException("imageC[0].length must be at least 4");
        }
        if (worldC[0].length != n) {
            throw new IllegalArgumentException("imageC[0].length must equal worldC[0].length");
        }
        
        // TODO: remove this temporary change
        kRadial = null;
        
        double[] b = new double[2*n];
        final double[][] xn = Camera.pixelToCameraCoordinates(imageC, kRadial, 
            kIntr.getIntrinsic(), useR2R4);
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
        double[][] r = kExtr.getRotation();
        double[] thetas = Rotation.extractRotation(r);
        double[] t = kExtr.getTranslation();
        
        // equation (19).  size is 1 X 9
        double[] h = new double[]{r[0][0], r[0][1], t[0], 
            r[1][0], r[1][1], t[1], -r[2][0], -r[2][1], -t[2]};
       
        // equation 20.  length is 2*N
        double[] fgp = map(worldC, h);
        
        // length is 2*N
        double[] bMinusFGP = MatrixUtil.subtract(b, fgp);
        
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
        
        // from a different lecture:: delta p LM: (J^T*J + lambda*I)^-1 * J^T * (b-f(g(p)))
        // length is 6
        double[] deltaPLM;
        
        //from eqn (22) :           inv( J^T J + λ diag(J^TJ)) · (b−f)
        double[] deltaPLM22;
        
        double eps = 1.e-5;
        
        // gain ratio:
        //gain = (f(p + delta p LM) - f(p)) / ell(delta p LM)
        //     where ell(delta p LM) is (delta p LM)^T * (lambda * (delta p LM)) + J^T * ( b - fgp))
        //gain = (f - fPrev) / (delta p LM)^T * (lambda * (delta p LM)) + J^T * ( b - fgp))
        double gainRatio;
        
        // loop stopping conditions:
        //   step length vanishes:  deltaPLM --> 0
        //   gradient of f(x) vanishes: -J^T * (b - fgp) --> 0
        double[] stepLengthCheck;
        double[] gradientCheck;
        final double tolP = 1.e-3;
        final double tolG = 1.e-3;
                    
        /*
        init: f=infinity
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
        
        int nIter = 0;
        while (nIter < nMaxIter) {
            
            nIter++;
                   
            fPrev = f;
            fgp = map(worldC, h);
            bMinusFGP = MatrixUtil.subtract(b, fgp);            
            f = evaluateObjective(bMinusFGP);
            
            // ===== calculate step ========
            j = calculateJ(worldC, h, thetas); //(2N) X 6
            jT = MatrixUtil.transpose(j);
            jTJ = MatrixUtil.multiply(jT, j);
    editing here: deltaPM not correct        
            deltaPLM22 = calculateDeltaPLM(jTJ, lambda, bMinusFGP);        
            // compare with other step function:
            deltaPLM = calculateDeltaPLM(jTJ, jT, lambda, bMinusFGP);
            
            System.out.printf("delta0=%s\ndelta1=%s\n", FormatArray.toString(deltaPLM, "%.3e"),
                FormatArray.toString(deltaPLM22, "%.3e"));
        
            // ======= stopping conditions ============
            stepLengthCheck = deltaPLM22;
            gradientCheck = MatrixUtil.multiplyMatrixByColumnVector(jT, bMinusFGP);
            MatrixUtil.multiply(gradientCheck, -1.);            
            if (isNegligible(stepLengthCheck, tolP) || !isNegligible(gradientCheck, tolG)) {
                break;
            }
            
             //TODO: consider whether to accept or reject step?
             // comments from: https://arxiv.org/pdf/1201.5885
             //           Transtruma & Sethna 2012
            // the qualitative effect of the damping term is to modify the 
            // eigenvalues of the matrix JT J + λDT D to be at least λ
            
            
            // ======= revise parameters =======
            
            // add deltaPM22 to p, which is (theta_x, theta_y, theta_z, t_x, t_y, t_z)
            updateBySteps(thetas, t, deltaPLM22);
            
            updateH(h, thetas, t);
            
            // ====== change lambda ======
            if (nIter > 1) {
                gainRatio = calculateGainRatio(f/2., fPrev/2.,
                    deltaPLM22, lambda, jT, bMinusFGP, eps);
                if (gainRatio > 0) {
                    // near the minimimum, which is good.
                    // decrease lambda
                    lambda /= lambdaF;
                } else {
                    // increase lambda
                    lambda *= lambdaF;
                }
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
     * map the homography matrix to the projected 2D point coordinates of the 
     * world reference points eqn (20).
     * @param worldC
     * @param h
     * @return 
     */
    static double[] map(double[][] worldC, double[] h) {
        
        int n = worldC[0].length;
        
        double[] f = new double[2*n];
        
        double sum = 0;
        int i, j;
        double X, Y, s;
        for (i = 0; i < n; ++i) {
            X = worldC[0][i];
            Y = worldC[1][i];
            s = h[6] * X + h[7] * Y + h[8];
            
            f[2*i] = (h[0] * X + h[1] * Y + h[2])/s;
            f[2*i+1] = ((h[3] * X + h[4] * Y + h[5])/s);            
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
        // (@N)X6)
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
     *                          6X6  * (6 * (2N)) * (2NX1) = 6 X (2N) * (2NX1) = 6X1
     * calculate the step as (J^T*J + lambda*I)^-1 * J^T * (b-f(g(p))
     * @param jTJ 
     * @param lambda
     * @param jT
     * @return an array of length 6 
     * @throws NotConvergedException 
     */
    private static double[] calculateDeltaPLM(double[][] jTJ, double[][] jT, 
        double lambda, double[] bMinusF) throws NotConvergedException {
        
        // (J^T*J + lambda*I)^-1 * J^T * (b-f(g(p))
        double[][] identity = MatrixUtil.createIdentityMatrix(6);
        MatrixUtil.multiply(identity, lambda);
        // 6 X 6
        double[][] a = MatrixUtil.elementwiseAdd(jTJ, identity);
        double[][] aInv = MatrixUtil.pseudoinverseFullRank(a);
        
        double[][] pt1 = MatrixUtil.multiply(aInv, jT);
        double[] step = MatrixUtil.multiplyMatrixByColumnVector(pt1, bMinusF);
        
        return step;
    }
    
    /**
     * calculate the step as inv( J^T J + λ diag(J^TJ)) · (b−f)
     * @param jTJ
     * @param lambda
     * @param bMinusF
     * @return an array of length 6
     * @throws NotConvergedException 
     */
    private static double[] calculateDeltaPLM(double[][] jTJ, 
        double lambda, double[] bMinusF) throws NotConvergedException {
       
        //inv( J^T J + λ diag(J^TJ)) · (b−f)
        
        int i, j;
        // J^T J + λ diag(J^TJ)     
        // 6 X 6
        double[][] a = MatrixUtil.copy(jTJ);
        for (i = 0; i < 6; ++i) {
            a[i][i] += (lambda*(jTJ[i][i]));
        }
        
        //6X6
        double[][] aInv = MatrixUtil.pseudoinverseFullRank(a);
              
        //512
        double[] step = MatrixUtil.multiplyMatrixByColumnVector(aInv, bMinusF);
        
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
     * @param jT
     * @param bMinusFGP
     * @return 
     */
    private static double calculateGainRatio(double f, double fPrev, 
        double[] deltaPLM, double lambda, double[][] jT, double[] bMinusFGP,
        double eps) {
             
        //      1X6          *            ( 6X1   +   6 X (2N) * (2NX1) )
        //      1X6                        6X1 
        //  1X1
        //(delta p LM)^T * (lambda * (delta p LM) + J^T * (b - fgp))
        double[] pt1 = Arrays.copyOf(deltaPLM, deltaPLM.length);
        MatrixUtil.multiply(pt1, lambda);
        double[] pt2 = MatrixUtil.multiplyMatrixByColumnVector(jT, bMinusFGP);
        pt2 = MatrixUtil.add(pt1, pt2);
        
        double ell = MatrixUtil.innerProduct(deltaPLM, pt2);
        
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

    /**
     * 
     * @param h input and output variable homography to be updated with given
     * rotation angles and translations.
     * @param thetas
     * @param t 
     */
    private static void updateH(double[] h, double[] thetas, double[] t) {
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
     * update thetas and t by deltaPLM22
     * @param thetas input and output array holding euler rotation angles 
     *    theta_x, theta_y, theta_
     * @param t input and output array holding translation vector in x,y,z. 
     * @param deltaPLM22 
     */
    private static void updateBySteps(double[] thetas, double[] t, double[] deltaPLM22) {
        int i;
        for (i = 0; i < 3; ++i) {
            thetas[i] += deltaPLM22[i];
        }
        for (i = 0; i < 3; ++i) {
            t[i] += deltaPLM22[i+3];
        }
    }
}
