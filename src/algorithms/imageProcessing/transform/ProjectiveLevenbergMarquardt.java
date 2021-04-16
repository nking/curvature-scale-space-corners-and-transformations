package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Camera.CameraExtrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraIntrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraMatrices;
import algorithms.matrix.MatrixUtil;

/**
 * iterative non-linear optimization using Levenberg-Marquardt algorithm
 * to minimize the re-projection error of perspective projection.
 * 
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
       https://stanford.edu/class/ee267/notes/ee267_notes_tracking.pdf
     * @param imageC
     * @param worldC
     * @param kIntr
     * @param kExtr
     * @param kRadial
     * @return 
     */
    public static CameraMatrices solve(double[][] imageC, double[][] worldC, 
        CameraIntrinsicParameters kIntr, CameraExtrinsicParameters kExtr, double[] kRadial) {
        
        int n = imageC[0].length;
        
        if (n < 4) {
            throw new IllegalArgumentException("imageC[0].length must be at least 4");
        }
        if (worldC[0].length != n) {
            throw new IllegalArgumentException("imageC[0].length must equal worldC[0].length");
        }
        
        //TODO: consider adding constraints suggestded in Seliski 2010:
        // u_0 and v_0 are close to half the image lengths and widths, respectively.
        // the angle between 2 image axes is close to 90.
        // the focal lengths along both axes are greater than 0.
        
        // extract pose as (theta_x, theta_y, theta_z, t_x, t_y, t_z)
        double[][] r = kExtr.getRotation();
        double[] thetas = Rotation.extractRotation(r);
        double[] t = kExtr.getTranslation();
        
        // equation (19)
        double[] h = new double[]{r[0][0], r[0][1], t[0], r[1][0], r[1][1], t[1], -r[2][0], -r[2][1], t[2]};
        
        /*
        Initialization: A = J^T*J, lambda=max{a_i_i}
        • Repeat until the step length vanishes, || delta x_LM || --> 0, 
          or the gradient of f(x) vanishes, del f = -J^T * r --> 0
          a) Solve (A + lambda*I)*(delta x) = J^T*r to get delta x_LM
          b) x = x + (delta x_LM)
          c) Adjust the damping parameter by checking the gain ratio
             1. rho > 0 Good approximation, decrease the damping parameter
             2. rho <= 0 Bad approximation, increase the damping parameter
        
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
        
        throw new UnsupportedOperationException("not yet complete");
    }
    
    /**
     * 
     * @param worldC
     * @param h
     * @return a (2*n)X9 matrix
     */
    static double[][] calculateJF(double[][] worldC, double[] h) {
        
        int n = worldC[0].length;
        
 //TODO: assert worldC[2][*] = be 0 for local device frame
        
        // (2*n) X 9
        double[][] jF = MatrixUtil.zeros(2*n, 9);
        int i, j;
        double x, y, s1, s2,s1sq;
        for (i = 0; i < n; ++i) {
            x = worldC[0][i];
            y = worldC[1][i];
            s1 = h[6] * x + h[7] * y + h[8];
            s2 = h[0] * x + h[1] * y + h[2];
            s1sq = s1*s1;
            
            jF[0][i] = x/s1;
            jF[1][i] = y/s1;
            jF[2][i] = 1./s1;
            //jF[3][i] = 0; jF[4][i] = 0; jF[5][i] = 0;
            jF[6][i] = -(s2/s1sq)*x;
            jF[7][i] = -(s2/s1sq)*y;
            jF[8][i] = -(s2/s1sq);
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
}
