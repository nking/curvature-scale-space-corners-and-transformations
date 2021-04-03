package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.CubicRootSolver;
import algorithms.misc.PolynomialRootSolver;
import algorithms.util.PolynomialFitter;
import java.util.Arrays;
import no.uib.cipr.matrix.NotConvergedException;

/**
 * utility methods for camera intrinsic and extrinsic matrices.
 * 
 * @author nichole
 */
public class Camera {
    
    /**
     *  create camera intrinsic matrix k with assumptions of square pixels
     * and no skew.  the focal length and optical centers should be in units of pixels.
     * NOTE that given the field of view (FOV) and the image dimensions,
     * one can roughly estimate the focal length as (image width/2) / tan(FOV/2).
     * @param focalLength focal length of camera in units of pixels.
     * @param xC x coordinate of camera optical center in image pixel coordinates.
     * @param yC y coordinate of camera optical center in image pixel coordinates.
     * @return intrinsic camera matrix in units of pixels.
     */
    public static double[][] createIntrinsicCameraMatrix(double focalLength,
        double xC, double yC) {
        
        double[][] k = new double[3][3];
        k[0] = new double[]{-focalLength, 0, xC};
        k[1] = new double[]{0, -focalLength, yC};
        k[2]= new double[]{0, 0, 1};
        
        return k;
    }
    
    /**
     *  create the inverse of camera intrinsic matrix k with assumptions of square pixels
     * and no skew.  the focal length and optical centers should be in units of pixels.
     * NOTE that given the field of view (FOV) and the image dimensions,
     * one can roughly estimate the focal length as (image width/2) / tan(FOV/2).
     * @param focalLength focal length of camera in units of pixels.
     * @param xC x coordinate of camera optical center in image pixel coordinates.
     * @param yC y coordinate of camera optical center in image pixel coordinates.
     * @return intrinsic camera matrix in units of pixels.
     */
    public static double[][] createIntrinsicCameraMatrixInverse(double focalLength,
        double xC, double yC) {
        
        double[][] k = new double[3][3];
        k[0] = new double[]{-1./focalLength, 0, -xC};
        k[1] = new double[]{0, -1./focalLength, -yC};
        k[2]= new double[]{0, 0, 1};
        
        return k;
    }
    
    /**
     *  create the inverse of camera intrinsic matrix k.
     * The multiplicative elements such as focal length,
     * are inverted, and the translation elements are multiplied by -1.
     * then the matrix is transposed.
     * @param kIntr intrinsic camera matrix.
     * @return intrinsic camera matrix inverse.
     */
    public static double[][] createIntrinsicCameraMatrixInverse(double[][] kIntr) {
        
        /*
        double[][] kInv = new double[3][3];
        kInv[0] = new double[]{1./kIntr[0][0], 0, -1*kIntr[0][2]};
        kInv[1] = new double[]{0, 1./kIntr[1][1], -1*kIntr[1][2]};
        kInv[2]= new double[]{0, 0, 1};
        */
        double[][] kInv = MatrixUtil.copy(kIntr);
        int i, j;
        double tol = 1e-7;
        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 2; ++j) {
                if (Math.abs(kIntr[i][j]) > tol) {
                    kInv[i][j] = 1./kInv[i][j];
                }
            }
        }
        kInv[0][2] *= -1;                    
        kInv[1][2] *= -1;
        
        return kInv;
    }
   
    /**
     * 
     * @param k camera intrinsics matrix of size 3 x 3.
     * @param r camera extrinsics rotation matrix of size 3 x 3.
     * @param t camera extrinsics translation vector of length 2.
     * @return the camera matrix resulting from intrinsic and extrinsic parameters.
     * the size is 3 x 4.
     */
    public static double[][] createCamera(double[][] k, double[][] r, double[] t) {
        if (k.length != 3 || k[0].length != 3) {
            throw new IllegalArgumentException("k must be 3 x 3");
        }
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r must be 3 x 3");
        }
        if (t.length != 3) {
            throw new IllegalArgumentException("t must be length 3");
        }
        
        /*
            4x4     
        [ R  -R*t ]
        [ 0   1   ]
        
        P = K * R * [I | -t]
        
        alternately, can write as P = K * [ R | -R*t]
        */
        double[] rt = MatrixUtil.multiplyMatrixByColumnVector(r, t);
        
        double[][] kExtr = new double[3][4];
        for (int i = 0; i < 3; ++i) {
            kExtr[i] = new double[4];
            System.arraycopy(r[i], 0, kExtr[i], 0, 3);
            kExtr[i][3] = rt[i];
        }
        
        double[][] p = MatrixUtil.multiply(k, kExtr);
        
        return p;
    }

    /**
     * not ready for use.  a quick rough method to estimate the 3D homogeneous point
     * from the 2D-homogenous point and this inverse camera matrix, with caveat 
     * about missing information on the last dimension.  
     * One should use reconstruction methods instead of this method.
     * to use:
     * <pre>
     * double[][] X = MatrixUtil.multiply(cameraInv, x);
     * then divide each column by the 3rd row.
     * </pre>
     * 
     * @param k camera intrinsics matrix of size 3 x 3.
     * @param r camera extrinsics rotation matrix of size 3 x 3.
     * @param t camera extrinsics translation vector of length 2.
     * @return the inverse camera matrix resulting from intrinsic and extrinsic parameters.
     * the size is 4x3
     */
    public static double[][] createCameraInverse(double[][] k, double[][] r, double[] t) {
        if (k.length != 3 || k[0].length != 3) {
            throw new IllegalArgumentException("k must be 3 x 3");
        }
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r must be 3 x 3");
        }
        if (t.length != 3) {
            throw new IllegalArgumentException("t must be length 3");
        }
        
        /*
        translation matrix: inverse changes the signs of the translation elements, but not the diagonal.
        rotation matrix: inverse is the transpose of rotation matrix.
        scaling matrix: inverse is performed on each element, that is, the reciprocal.
        */
                
        double[] tInv = Arrays.copyOf(t, t.length);
        tInv[0] *= -1;
        tInv[1] *= -1;
                
        double[] rTInv = MatrixUtil.multiplyMatrixByColumnVector(r, tInv);
        
        double[][] kInv = Camera.createIntrinsicCameraMatrixInverse(k);
        
        /*           
        inverse of   K * R * [I | -t]             
            
        is  | r  | r*tInv ]^T  * kInv
        */
        
        double[][] cInv = new double[3][4];
        for (int i = 0; i < 3; ++i) {
            cInv[i] = new double[4];
            System.arraycopy(r[i], 0, cInv[i], 0, 3);
            cInv[i][3] = rTInv[i];
        }
        cInv = MatrixUtil.transpose(cInv);
        
        cInv = MatrixUtil.multiply(cInv, kInv);
        
        return cInv;
    }
    
    /**
     * applies radial distortion to distortion-free camera centered coordinates
     * then multiplies by the camera intrinsics to result in distorted coordinates 
     * in the image reference frame in units of pixels.
     * In terms of Table 1 of Ma et al. 2004, the input is a double array of (x, y)
     * and the output is a double array of (u_d, v_d).
     * Also useful reading is NVM Tools by Alex Locher
    https://github.com/alexlocher/nvmtools.git
     * @param xC distortion-free camera centered coordinates. 
     * format is 3XN for N points.  
     * In terms of Table 1 of Ma et al. 2004, this is a double array of (x, y).
     * @param rCoeffs radial distortion vector of length 2 or radial and tangential
     * distortion vector of length 5.  can be null to skip lens distortion correction.
     * @param focalLength focal length of camera in units of pixels.
     * @param centerX x coordinate of camera optical center in image pixel coordinates.
     * @param centerY y coordinate of camera optical center in image pixel coordinates.
     * @return pixels in the reference frame of image with distortion applied.
     * In terms of Table 1 of Ma et al. 2004, this is a double array of (u_d, v_d)
     */
    public static double[][] cameraToPixelCoordinates(double[][] xC, double[] rCoeffs,
        double focalLength, double centerX, double centerY) {
        
        // http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/parameters.html
        
        double[][] cc = MatrixUtil.copy(xC);
        for (int i = 0; i < xC[0].length; ++i) {
            // normalized pinhole projection X_c/Z_c and 
            cc[0][i] /= xC[2][i];
            cc[1][i] /= xC[2][i];
        }
        
        if (rCoeffs != null) {
            // input and output cc are in camera reference frame
            cc = applyRadialDistortion(cc, rCoeffs[0], rCoeffs[1]);
            //TODO: consider implementing the higher order terms to include tangential
        }
         
        focalLength = Math.abs(focalLength);
        
        double[][] cameraIntr = Camera.createIntrinsicCameraMatrix(focalLength, centerX, centerY);
                       
        cc = MatrixUtil.multiply(cc, cameraIntr);
        
        return cc;
    }
    
    /* converts pixel coordinates to camera coordinates by transforming them to camera 
    reference frame then removing radial distortion.
    The radial distortion removal follows Ma et al. 2004.
    The input in terms of Table 1 of Ma et al. 2004 is a double array of (u_d, v_d)
    and the output is a double array of (x, y).
    Also useful reading is NVM Tools by Alex Locher
    https://github.com/alexlocher/nvmtools.git
    
     * @param xC points in the camera centered reference frame. 
     * format is 3XN for N points.  
     * @param rCoeffs radial distortion vector of length 2 or radial and tangential
     * distortion vector of length 5.  can be null to skip lens distortion correction.
     * @param focalLength focal length of camera in units of pixels.
     * @param centerX x coordinate of camera optical center in image pixel coordinates.
     * @param centerY y coordinate of camera optical center in image pixel coordinates.
     * @return pixels in the reference frame of 
     */
    public static double[][] pixelToCameraCoordinates(double[][] x, double[] rCoeffs,
        double focalLength, double centerX, double centerY) throws NotConvergedException {
        
        // http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/parameters.html
        
        focalLength = Math.abs(focalLength);
        
        double[][] cameraIntrInv = Camera.createIntrinsicCameraMatrixInverse(
            focalLength, centerX, centerY);
        
        // put x into camera coordinates reference frame:
        double[][] pix = MatrixUtil.multiply(cameraIntrInv, x);
        
        if (rCoeffs != null) {
            pix = removeRadialDistortion(pix, rCoeffs[0], rCoeffs[1]);
            //TODO: consider implementing the higher order terms to include tangential
        }
                
        return pix;
    }
    
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
        discriminator delta = (q/2)^2 + (p/3)^3
           delta .gt. 0 there is 1 real root
           delta .eq. 0 has a multiple root
           delta .lt. 0 there are 3 real roots and the middle one is what is 
                     needed since the first root is at a negative radius 
                     and the third lies beyond the positive turning point
    
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
        
        double r, r2, fr;
        int i;
                
        for (i = 0; i < distorted[0].length; ++i) {
            r2 = distorted[0][i]*distorted[0][i] + distorted[1][i]*distorted[1][i];
            r = Math.sqrt(r2);
            //f(r) = (1 + k1*r + k2*r^2)
            fr = 1 + k1*r + k2*r2;
            distorted[0][i] *= fr;
            distorted[1][i] *= fr;
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
        discriminator delta = (q/2)^2 + (p/3)^3
           delta .gt. 0 there is 1 real root
           delta .eq. 0 has a multiple root
           delta .lt. 0 there are 3 real roots and the middle one is what is 
                     needed since the first root is at a negative radius 
                     and the third lies beyond the positive turning point
    
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
    @param xC distorted points in the camera reference frame.  format is 3XN where N is the
    number of points.  These are (x_d, x_d) pairs in terms of Table 1 in Ma et al. 2004.
    @param k1 first radial distortion coefficient
    @param k2 second radial distortion coefficient
    @return undistorted points in the camera reference frame.  Format is 3XN where N is the
    number of points.  These are (x, y) pairs in terms of Table 1 in Ma et al. 2004.
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
        discriminator delta = (q/2)^2 + (p/3)^3
           delta .gt. 0 there is 1 real root
           delta .eq. 0 has a multiple root
           delta .lt. 0 there are 3 real roots and the middle one is what is 
                     needed since the first root is at a negative radius 
                     and the third lies beyond the positive turning point
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
            c = -rd/k2;
            p = b - (a2/3.);
            q = (2.*a3/27.) - (a*b/3.) + c;
            
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
