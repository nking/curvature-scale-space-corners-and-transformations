package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Camera.CameraExtrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraIntrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraMatrices;
import algorithms.matrix.MatrixUtil;
import algorithms.matrix.MatrixUtil.SVDProducts;
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
    
    /**
     * 
     * @param n n is the number of points in each image which is the
              same for all images.
     * @param coordsI  holds the image coordinates in pixels of
               features present in all images ordered in the same
               manner and paired with features in coordsW.
               It is a 2 dimensional double array of format
               3 X (N*n) where N is the number of images.
               the first row is the x coordinates, the second row
               is the y coordinates, and the third row is "1"'s.
     * @param coordsW holds the world coordinates of features, ordered
               by the same features in the images.
               the first row is the X coordinates, the second row
               is the Y coordinates, and the third row is assumed to
               be zero as the scale factor is lost in the homography.
               It is a 2 dimensional double array of format
               3 X n
     * @return camera intrinsic parameters, extrinsic parameters, and radial
     * distortion coefficients
     */
    public static CameraMatrices estimateCamera(int n, double[][] coordsI, 
        double[][] coordsW) throws NotConvergedException {
        
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
        
        //TODO: consider normalizing coordsI by coordsI[2][*] if the last row
        //      is not already 1's
        
        //(1) for each image: invoke homography solver using SVD
        double[][] h = MatrixUtil.zeros(nImages*3, 3);
        
        double[][] g, cI;
        int i, i2;
        for (i = 0; i < nImages; ++i) {
            cI = MatrixUtil.copySubMatrix(coordsI, 0, 3, i, i+n);
            g = solveForHomography(cI, coordsW);
            
            /*
            h image 0
            h image 1
            h image 2
            */
            for (i2 = 0; i2 < 3; ++i2) {
                System.arraycopy(g[i2], 0, h[i+i2], 0, g[i2].length);
            }
        }
        
        //(2) for all homographies, solve for the camera intrinsic parameters
        CameraIntrinsicParameters kIntr = solveForIntrinsic(h);
        
        CameraMatrices cameraMatrices = new CameraMatrices();
        cameraMatrices.setIntrinsics(kIntr);
        
        //(3) for each image homography and inverse intrinsic parameter matrix,
        //    estimate the extrinsic parameters for the pose of the camera for that image.
        Camera.CameraExtrinsicParameters kExtr;
        
        for (i = 0; i < nImages; ++i) {
            cI = MatrixUtil.copySubMatrix(h, 3*i, 3+3*i, 0, 3);
            kExtr = solveForExtrinsic(kIntr, cI);
            cameraMatrices.addExtrinsics(kExtr);
        }
        
        // (4) estimate the radial distortion coefficients
        
        // (5) optimization to improve the parameter estimates
        
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    /**
     * for a given set of feature coordinates in image reference frame and in
     * world coordinate system, calculates the homography following the 
     * algorithm in Ma et al. 2003.
     * @param coordsI holds the image coordinates in pixels of features present in image i
     * @param coordsW holds the world coordinates of features present in image 1 corresponding
               to the same features and order of coordsI_i
     * @return the homography, projection matrix
     */
    static double[][] solveForHomography(double[][] coordsI, double[][] coordsW) throws NotConvergedException {
        
        if (coordsI.length != 3) {
            throw new IllegalArgumentException("coordsI must have 3 rows.");
        }
        if (coordsW.length != 3) {
            throw new IllegalArgumentException("coordsW must have 3 rows.");
        }
        int n = coordsI.length;
        if (coordsW[0].length != n) {
            throw new IllegalArgumentException("coordsW must have same number of rows as coordsI.");
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

          Let h_i be the ith column vector of H:
              h_i = ]h_i_1]^T = [h_i_1  h_i_2  h_i_3]
                    [h_i_2]
                    [h_i_3]
        */
        
        // 2*n X 9       
        double u, v, X, Y;
        double[][] ell = new double[2*n][9];
        for (int i = 0; i < n; ++i) {
            u = coordsI[0][i];
            v = coordsI[1][i];
            X = coordsW[0][i];
            Y = coordsW[1][i];
            ell[2*i] = new double[]{X, Y, 1, 0, 0, 0, -u*X, -u*Y, -u};
            ell[2*i + 1] = new double[]{0, 0, 0, X, Y, 1, -v*X, -v*Y, -v};
        }
        
        MatrixUtil.SVDProducts svd = MatrixUtil.performSVD(ell);
        
        // vT is 9X9.  last row in vT is the eigenvector for the smallest eigenvalue
        double[] xOrth = svd.vT[svd.vT.length - 1];
        
        double[][] h = new double[3][3];
        for (int i = 0; i < 3; i++) {
            h[i] = new double[3];
            h[i][0] = xOrth[(i * 3) + 0];
            h[i][1] = xOrth[(i * 3) + 1];
            h[i][2] = xOrth[(i * 3) + 2];
        }
        
        return h;
    }
    
    /**
     * estimate the camera intrinsic parameters from the image homographies.
     * @param h H as (3*NImages)x3 homography, projection matrices
              where each image homography is stacked row-wise
     * @return the camera intrinsic parameters.
     */
    static CameraIntrinsicParameters solveForIntrinsic(double[][] h) throws NotConvergedException {
        
        if (h[0].length != 3) {
            throw new IllegalArgumentException("h must have 3 columns");
        }
        
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
        
        //Vij = [hi1*hj1, hi1*hj2 + hi2*hj1, hi2*hj2, hi3*hj1 + hi1*hj3, hi3*hj2 + hi2*hj3, hi3*hj3]T 
        for (int i = 0; i < n; ++i) {
            // h_i is the ith column vector of H
            // h11 = column 0 of h, first element: h[0][0]
            // h12 = column 0 of h, 2nd element:   h[1][0]
            // h13 = column 0 of h, 3rd element:   h[2][0]
            // h21 = column 1 of h, first element: h[0][1]
            // h22 = column 1 of h, 2nd element:   h[1][1]
            // h23 = column 1 of h, 3rd element:   h[2][1]
            
            //V12 = [h11*h21, h11*h22 + h12*h21, h12*h22, h13*h21 + h11*h23, 
            //       h13*h22 + h12*h23, h13*h23]T 
            v[2*i] = new double[]{
                h[0+3*i][0]*h[0+3*1][1],
                h[0+3*i][0]*h[1+3*1][1] + h[1+3*i][0]*h[0+3*1][1],
                h[1+3*i][0]*h[1+3*1][1],
                h[2+3*i][0]*h[0+3*1][1] + h[0+3*i][0]*h[2+3*1][1],
                h[2+3*i][0]*h[1+3*1][1] + h[1+3*i][0]*h[2+3*1][1],
                h[2+3*i][0]*h[2+3*1][1]
            };
            
            //V11 = [
            // h11*h11 , 
            // h11*h12 + h12*h11, 
            // h12*h12, 
            // h13*h11 + h11*h13, 
            // h13*h12 + h12*h13, 
            // h13*h13]T
            v[2*i + 1] = new double[]{
                h[0+3*i][0]*h[0+3*1][0] - v[2*i][0],
                h[0+3*i][0]*h[1+3*1][0] + h[1+3*i][0]*h[0+3*1][0] - v[2*i][1],
                h[1+3*1][0]*h[1+3*1][0] - v[2*i][2],
                h[2+3*i][0]*h[0+3*1][0] + h[0+3*i][0]*h[2+3*1][0] - v[2*i][3],
                h[2+3*i][0]*h[1+3*1][0] + h[1+3*i][0]*h[2+3*1][0] - v[2*i][4],
                h[2+3*i][0]*h[2+3*i][0] - v[2*i][5]
            };
        }
        
        //Vb = 0 and b = [B11, B12, B22, B13, B23, B33]^T
        SVDProducts svd = MatrixUtil.performSVD(v);
        
        // vT is 9X9.  last row in vT is the eigenvector for the smallest eigenvalue
        double[] b = svd.vT[svd.vT.length - 1];
        
        /*
        double[][] B = MatrixUtil.zeros(3, 3);
        B[0][0] = b[0];
        B[0][1] = b[1];
        B[1][1] = b[2];
        B[0][2] = b[3];
        B[1][2] = b[4];
        B[2][2] = b[5];
        */
        
        //       0    1    2    3    4    5
        //b = [B11, B12, B22, B13, B23, B33]^T
        
        //Ma et al. 2003 eqn (26)
        // v0 = (B12*B13 - B11*B23)/(B11*B22 - B12*B12)
        // lambda = B33 - (B13*B13 - v0*(B12*B13 - B11*B23))/B11
        // alpha = sqrt( lambda/B11 )
        // beta = sqrt( lambda*B11 / (B11*B22 - B12*B12) )
        // gamma = -B12*alpha*alpha*beta / lambda
        // u0 = (gamma*v0/beta) - B13*alpha*alpha/lambda
        
        double v0 = (b[1]*b[3] - b[0]*b[4])/(b[0]*b[2] - b[1]*b[1]);
        double lambda = b[5] - (b[3]*b[3] - v0*(b[1]*b[3] - b[0]*b[4]))/b[0];
        double alpha = Math.sqrt(lambda / b[0]);
        double beta = Math.sqrt( lambda*b[0] / (b[0]*b[2] - b[1]*b[1]) );
        double gamma = -b[1]*alpha*alpha*beta / lambda;
        double u0 = (gamma*v0/beta) - b[3]*alpha*alpha/lambda;
        
        double[][] kIntr = Camera.createIntrinsicCameraMatrix(
            alpha, beta, u0, v0, gamma);
            
        CameraIntrinsicParameters intrinsics = new CameraIntrinsicParameters();
        intrinsics.setIntrinsic(kIntr);
        intrinsics.setLambda(lambda);
        
        return intrinsics;
    }
    
    
    /**
     * estimate the extrinsic parameters
     * @param kIntr camera intrinsic parameters
     * @param h homography for the projection for an image
     * @return 
     */
    private static Camera.CameraExtrinsicParameters solveForExtrinsic(
        CameraIntrinsicParameters kIntr, double[][] h) throws NotConvergedException {
        
        double[][] kIntrInv = Camera.createIntrinsicCameraMatrixInverse(kIntr.getIntrinsic());
        
        //h_i is the ith column vector of H
        //r1 = λ * A^−1 * h1
        //r2 = λ * A^−1 * h2
        //r3 = r1×r2
        // t = λ * A^−1 * h3
        
        double[] h1 = MatrixUtil.extractColumn(h, 0);
        double[] h2 = MatrixUtil.extractColumn(h, 1);
        double[] h3 = MatrixUtil.extractColumn(h, 2);
        
        //λ = 1/||A−1h1||2 = 1/||A−1h2||2
        double lambda = 1./sumOfSquares(MatrixUtil.multiplyMatrixByColumnVector(kIntrInv, h1));
        double lambda2 = 1./sumOfSquares(MatrixUtil.multiplyMatrixByColumnVector(kIntrInv, h2));
        
        double[][] r1M = MatrixUtil.copy(kIntrInv);
        MatrixUtil.multiply(r1M, lambda);
        double[][] r2M = MatrixUtil.copy(r1M);
        double[][] tM = MatrixUtil.copy(r1M);
        
        double[] r1 = MatrixUtil.multiplyMatrixByColumnVector(r1M, h1);
        double[] r2 = MatrixUtil.multiplyMatrixByColumnVector(r2M, h2);
        double[] t = MatrixUtil.multiplyMatrixByColumnVector(tM, h3);
        
        double[] r3 = MatrixUtil.crossProduct(r1, r2);
        
        double[][] r = MatrixUtil.zeros(3, 3);
        for (int row = 0; row < 3; ++row) {
            r[row][0] = r1[row];
            r[row][1] = r2[row];
            r[row][2] = r3[row];
        }
        
        SVDProducts svd = MatrixUtil.performSVD(r);
        
        r = MatrixUtil.multiply(svd.u, svd.vT);
        
        CameraExtrinsicParameters kExtr = new Camera.CameraExtrinsicParameters();
        kExtr.setRotation(r);
        kExtr.setTranslation(t);
        
        return kExtr;
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
    TODO: implement Ma et al. 2004 Table 2, #6.
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

    private static double sumOfSquares(double[] m) {
        double sum = 0;
        for (int i = 0; i < m.length; ++i) {
            sum += (m[i]*m[i]);
        }
        return sum;
    }

}
