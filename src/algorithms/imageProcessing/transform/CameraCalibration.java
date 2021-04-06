package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.CubicRootSolver;
import algorithms.misc.PolynomialRootSolver;
import java.util.Arrays;
import no.uib.cipr.matrix.NotConvergedException;

/**
 * estimate the camera intrinsic and extrinsic parameters using 3 images
 * or more of the same objects with different camera poses.
 * Following the algorithm os Ma, Chen, & Moore 2003 in 
 * "Camera Calibration: a USU Implementation" available as a preprint
 * at arXiv.
 * <pre>
 * "one image observed by a camera can provide 2 constraints about this camera’s 
 * intrinsic parameters that are regarded to be unchanged here. 
 * With 3 images observed by the same camera, 6 constraints are established and 
 * we are able to recover the 5 intrinsic parameters. Once the intrinsic 
 * parameters are known, we can estimate the extrinsic parameters, 
 * the distortion coefficients (k1,k2), and put every initial guess of these 
 * parameters into some nonlinear optimization routine to get the final estimations. 
 * </pre>
 * @author nichole
 */
public class CameraCalibration {
    
    /*
    need data for N images where data is n features in each
    image as image coordinates and world coordinates.
    At least N=3 images are required.    
    */
    
    
    /**
    apply radial distortion to distortion-free camera centered coordinates using 
    the algorithm of Ma, Chen & Moore (which is Ma et al. 2003) for the
    distortion function expressed as f(r) = 1 + k1*r + k2*r^2;
    In terms of the variables outlined below, the algorithm input is
    (x, y), k1, k2, and cameraIntrinsics and the output is (x_d, y_d).
    
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
    
    ==> To Apply Radial Distorion, given coefficients k1, k2 and coordinates
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
    TODO: implement Ma et al. 2004 Tabe 2, #6.
    </pre>
    @param xC distortion-free camera centered coordinates.  format is 3XN where N is the
    number of points.
    In terms of Table 1 of Chen et al. 2004, this is a double array of (x, y).
    @param k1 first radial distortion coefficient
    @param k2 second radial distortion coefficient
    @return  distorted camera centered coordinates in format 3XN where N is
    the number of points.
    In terms of Table 1 of Chen et al. 2004, this is a double array of (x_d, y_d).
    */
    public static double[][] applyRadialDistortion(double[][] xC, double k1, double k2) {
        
        if (xC.length != 3) {
            throw new IllegalArgumentException("xC.length must be 3");
        }
        
        double[][] distorted = MatrixUtil.copy(xC);
        
        double r, r2, fr, signx, signy, c2p1, fx, fy;
        int i;
                
        for (i = 0; i < distorted[0].length; ++i) {
            r2 = distorted[0][i]*distorted[0][i] + distorted[1][i]*distorted[1][i];
            r = Math.sqrt(r2);
            //f(r) = (1 + k1*r + k2*r^2)
            fr = 1 + k1*r + k2*r2;
            //distorted[0][i] *= fr;
            //distorted[1][i] *= fr;
            
            // following Ma et al. 2004 Table 2,column 3 for model #3:
            // where c = y_d/x_d = y/x
            // f(x) = (1 + k1*math.sqrt(1+c^2)*x*sign(x) + k2*(1+c^2)*x^2)
            c2p1 = Math.pow(distorted[1][i]/distorted[0][i], 2.) + 1;
            signx = (distorted[0][i] < 0) ? -1 : 0;
            signy = (distorted[1][i] < 0) ? -1 : 0;
            fx = 1 + (k1*Math.sqrt(c2p1)*distorted[0][i]*signx) + (k2*c2p1*distorted[0][i]*distorted[0][i]);
            fy = 1 + (k1*Math.sqrt(c2p1)*distorted[1][i]*signy) + (k2*c2p1*distorted[1][i]*distorted[1][i]);
            distorted[0][i] *= fx;
            distorted[1][i] *= fy;
        }
                
        return distorted;
    }
    
    /**
    apply radial distortion correction to points in the camera reference frame.
    This method is used by pixelToCameraCoordinates().
    The algorithm follows Ma, Chen, & Moore 2004 to correct the distortion
    estimated as f(r) = 1 + k1*r + k2*r^2;
    In terms of the variables outlined in comments below, the algorithm input is
    distorted points as a double array of (x_d, x_d), and the radial distortion
    coefficients k1, k2.  The output is a a double array of (x, y).
    
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
    
    ==> To Apply Radial Distorion, given coefficients k1, k2 and coordinates
      (eqn 7) r_d = r * f(r) = r*(1 + k1*r + k2*r^2 + k3*r^3 + ...)
      (eqn 8) using only 2 coeffs: f(r) = (1 + k1*r + k2*r^2)
    
      Ma et al. 2003 distortion model:
         x_d = x * f(r)
         y_d = y * f(r)
    
    ==> Radial Undistortion, Section 2.2 of Ma et al. 2004:
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
         
    </pre>
    TODO: implement Ma et al. 2004 Tabe 2, #6.
    @param xC distorted points in the camera reference frame, preseumably 
    already center subtracted.  format is 3XN where N is the
    number of points.  These are (x_d, x_d) pairs in terms of Table 1 in Ma et al. 2004.
    @param k1 first radial distortion coefficient
    @param k2 second radial distortion coefficient
    @return undistorted points in the camera reference frame.  Format is 3XN where N is the
    number of points.  These are (x, y) pairs in terms of Table 1 in Ma et al. 2004.
    @throws no.uib.cipr.matrix.NotConvergedException
    */
    static double[][] removeRadialDistortion(double[][] xC, double k1, double k2) throws NotConvergedException {
        
        if (xC.length != 3) {
            throw new IllegalArgumentException("xC.length must be 3");
        }
                
        double[][] corrected = MatrixUtil.copy(xC);
        
        /*
        from k2*r^3 +k1*r^2 + r - r_d = 0
        solve 
           r_bar^3 + r_bar*p + q = 0
        where 
          r_bar = r + (a/3)  <=== r = r_bar - (a/3)
          a = k1/k2
          b = 1/k2
          c = −r_d/k2
          p = b − (a^2/3)
          q = (2a^3)/27 − ab/3 + c
        */
        
        double a = k1/k2;
        double a2 = a*a;
        double a3 = a2*a;
        double b = 1./k2;
        double rd, c, p, q, r, fr;
        double tol = 1e-5;
        int i;
        for (i = 0; i < xC[0].length; ++i) {
            rd = Math.sqrt(corrected[0][i]*corrected[0][i] + corrected[1][i]*corrected[1][i]);
            if (Math.abs(rd) < tol) {
                continue;
            }
            c = -rd/k2;
            p = b - (a2/3.);
            q = (2.*a3/27.) - (a*b/3.) + c;
            
            // Ma et al. 2004 expect rd = 0, r=0 when cubic root discriminant = 0.
            //   so consider whether to continue use cubic root for that case
            
            double[] rBar = CubicRootSolver.solveUsingDepressedCubic(p, q);
            if (rBar == null || rBar.length == 0) {
                //k2*r^3 +k1*r^2 + r - r_d = 0
                rBar = PolynomialRootSolver.realRoots(new double[]{k2, k1, 1, -rd});
                if (rBar == null || rBar.length == 0) {
                    //TODO: consider how to handle this case
                    r = 0;
                }
                r = rBar[0];
            } else if (rBar.length == 1) {
                r = rBar[0] - (a/3.);
                // check solution: 
                //  k2*r^3 +k1*r^2 +r - r_d = 0
                double chk = k2 * r * r * r + k1 * r * r + r - rd;
                assert(Math.abs(chk) < tol);
            } else {
                assert(rBar.length == 3);
                if ((Math.pow(q/2., 2) + Math.pow(p/3., 3)) < tol) {
                    //from Ma et al 2004:
                    //if ∆<0,then there are three solutions [Pearson, 1983]. 
                    //In general, the middle one is what we need, since the 
                    //first root is at a negative radius and the third lies 
                    //beyond the positive turning point 
                    Arrays.sort(rBar);
                }
                r = rBar[1] - (a/3.);
                // check solution: 
                //  k2*r^3 +k1*r^2 +r - r_d = 0
                double chk = k2 * r * r * r + k1 * r * r + r - rd;
                assert(Math.abs(chk) < tol);
            }
                        
            fr = 1 + k1*r + k2*r*r;
                  
            // eqn (5) from Ma et al. 2004
            //u = u0 + ((ud - u0)/fr);
            //v = v0 + ((vd - v0)/fr);
            // but points are currently in image reference frame as xd, yd
            // and the radial distortion w.r.t image reference frame
            // so presumably should use eqn (4) instead
            // eqn 4): xd = x *f(r) ==> x = x_d/f(r)
            corrected[0][i] /= fr;
            corrected[1][i] /= fr;
        }
                
        return corrected;
    }
}
