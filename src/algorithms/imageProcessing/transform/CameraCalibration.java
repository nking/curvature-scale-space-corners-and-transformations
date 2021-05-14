package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Camera.CameraExtrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraIntrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraMatrices;
import algorithms.matrix.MatrixUtil;
import algorithms.matrix.MatrixUtil.SVDProducts;
import algorithms.misc.CubicRootSolver;
import algorithms.misc.PolynomialRootSolver;
import algorithms.util.FormatArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import no.uib.cipr.matrix.NotConvergedException;

/**
 * estimate the camera intrinsic and extrinsic parameters using 3 images
 * or more of the same objects with different camera poses.
 * Following the algorithm os Ma, Chen, & Moore 2003 in 
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
               The columns hold each image in order and within each image's
               columns are the features presented in the same order in each image.
               In Table 1 of Ma, Chen, & Moore 2003 "Camera Calibration"
               these are the (u_d, v_d) pairs.
     * @param coordsW holds the world coordinates of features, ordered
               by the same features in the images.
               the first row is the X coordinates, the second row
               is the Y coordinates, and the third row is 1's 
               (Z_w = 0, the scale factor is lost in the homography).
               It is a 2 dimensional double array of format
               3 X n
       @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
    f(r) = 1 +k1*r^2 + k2*r^4 if true,
    else use model #3 f(r) = 1 +k1*r + k2*r^2 if true.
     * @return camera intrinsic parameters, extrinsic parameters, and radial
     * distortion coefficients
     */
    public static CameraMatrices estimateCamera(final int n, double[][] coordsI, 
        double[][] coordsW, boolean useR2R4) throws NotConvergedException {
        
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
        
        int i;
        
        // the normalization centroidX, centroidY, and sigma values applied to x1 and x2:
        //double[] tt = new double[3*nImages];        
        //coordsI = extractAndNormalize(coordsI, 0, n, tt);
        
        //TODO: consider normalizing coordsI by coordsI[2][*] if the last row
        //      is not already 1's
        
        //(1) for each image: invoke homography solver using SVD on a (2N)X9 DLT
        //    where the homographies are the projections of the 3D points onto the images
        //    and the model plane has Z=0 (hence the (2N)X9 DLT instead of (2N)X12 DLT)
        double[][] h = solveForHomographies(coordsI, coordsW, n, nImages);
        
        
        //(2) using all homographies, solve for the camera intrinsic parameters
        // this is where at least 3 points are needed per image to euqal the number of unknown intrinsic parameters.
        CameraIntrinsicParameters kIntr = solveForIntrinsic(h);
        
        CameraMatrices cameraMatrices = new CameraMatrices();
        cameraMatrices.setIntrinsics(kIntr);
        
        List<CameraExtrinsicParameters> extrinsics = solveForExtrinsics(kIntr, h, nImages);
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
        //          which is not necessarily the image centerin Section (3.3.3)
        
        // Using the estimated intrinsic and extrinsic parameters, 
        //  we can get the ideal projected image points Along with the real 
        //  observed image points, we can estimate the 
        // two distortion coefficients (k1,k2)
        
        // calculating the projected world coordinates using eqn (17)
        double[] u = new double[n*nImages];
        double[] v = new double[n*nImages];
        calculateProjected(coordsW, h, u, v);
                
        double[] kRadial = solveForRadialDistortion(coordsI, u, v, cameraMatrices,
            useR2R4);
        
        // (5) optimization to improve the parameter estimates
        
        // ============ iterate over the above steps after non-linear optimization
        //  for extrinsic parameters.
        double[][] cI;
        CameraExtrinsicParameters kExtr;
        CameraExtrinsicParameters extrinsic;
        int nMaxIter = 100;
        for (i = 0; i < nImages; ++i) {
            
            cI = MatrixUtil.copySubMatrix(coordsI, 0, 2, n*i, n*(i + 1)-1);
            
            kExtr = cameraMatrices.getExtrinsics().get(i);
            
            // improve the extrinsic parameter estimates:
            extrinsic = ProjectiveLevenbergMarquardt.solve(cI, coordsW, kIntr, 
                kExtr, kRadial, nMaxIter, useR2R4);
            
            cameraMatrices.getExtrinsics().set(i, extrinsic);
        }        
            
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
        int n = coordsI[0].length;
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

          Let h_i be the ith row of H:
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
            // eqn(15) of Ma et al. 2003
            ell[2*i]     = new double[]{X, Y, 1, 0, 0, 0, -u*X, -u*Y, -u};
            ell[2*i + 1] = new double[]{0, 0, 0, X, Y, 1, -v*X, -v*Y, -v};
        }
        
        MatrixUtil.SVDProducts svd = MatrixUtil.performSVD(ell);
        
        // vT is 9X9.  last row in vT is the eigenvector for the smallest eigenvalue
        double[] xOrth = svd.vT[svd.vT.length - 1];
        //h00, h01, h02, h10, h11, h12,h20,h21,h22
        
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
        
        System.out.printf("h=\n%s\n", FormatArray.toString(h, "%.3e"));
        
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
        System.out.printf("b=%s\n", FormatArray.toString(b, "%.3e"));
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
        
        System.out.printf("v0=%.4e (exp=220.866)\n", v0);
        System.out.printf("lambda=%.4e\n", lambda);
        System.out.printf("alpha=%.4e\n", alpha);
        System.out.printf("beta=%.4e\n", beta);
        System.out.printf("gamma=%.4e\n", gamma);
        System.out.printf("u0=%.4e\n", u0);
        
        double[][] kIntr = Camera.createIntrinsicCameraMatrix(
            alpha, beta, u0, v0, gamma);
            
        CameraIntrinsicParameters intrinsics = new CameraIntrinsicParameters();
        intrinsics.setIntrinsic(kIntr);
        intrinsics.setLambda(lambda);
        
        return intrinsics;
    }
    
    
    /**
     * following Ma et al. 2003
     * estimate the extrinsic parameters from the image of the absolute conic.
     * @param kIntr camera intrinsic parameters
     * @param h homography for the projection for an image
     * @return 
     */
    static Camera.CameraExtrinsicParameters solveForExtrinsic(
        CameraIntrinsicParameters kIntr, double[][] h) throws NotConvergedException {
        
        // notes from Serge Belongie lectures from Computer Vision II, CSE 252B, USSD
        // homogeneous repr of a point is x_vec = (x, y, 1)^T
        // equation f a line is ell = a*x + b*y + c = 0;
        // line rewritten in homogeneous coordinatrs is x_vec^T * ell.
        // general conic in 3 dimensions is a*x^2 + b*x*y + c*y^2 + d*x*z + e*y*z + f*z^2 = 0.
        //     rewritten using 2D homogenouse coords, quadratice form: x_vec^T * C * x_vec = 0
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
        //          (1) SVD(A).s[5] == 0, and n=5, then a conic xists that
        //              fits the data exactly
        //          (2) SVD(A).s[5] >=0, and n > 5, then the value is the goodness of fit
        //          (3) n < 5, the conic is undetermined and requires other means to solve.
        //        
        
        // points at infinity, a.k.a. ideal points, have the form (x, y, 0)^T
        // the line at infinity is (0, 0, 1)^T.
        
        double[][] aInv = Camera.createIntrinsicCameraMatrixInverse(kIntr.getIntrinsic());
        
        //h_i is the ith column vector of H
        //r1 = λ * A^−1 * h1
        //r2 = λ * A^−1 * h2
        //r3 = r1×r2
        // t = λ * A^−1 * h3
        
        double[] h1 = MatrixUtil.extractColumn(h, 0);
        double[] h2 = MatrixUtil.extractColumn(h, 1);
        double[] h3 = MatrixUtil.extractColumn(h, 2);
        
        //λ = 1/||A−1h1||2 = 1/||A−1h2||2
        double lambda = 1./sumOfSquares(MatrixUtil.multiplyMatrixByColumnVector(aInv, h1));
        double lambda2 = 1./sumOfSquares(MatrixUtil.multiplyMatrixByColumnVector(aInv, h2));
        //double scaleFactor = 2./(Math.sqrt(sumOfSquares(h1)) + Math.sqrt(sumOfSquares(h1)));
        System.out.printf("lambda1=%.3e, lambda2=%.3e\n", lambda, lambda2);
        
        double[][] lambdaAInv = MatrixUtil.copy(aInv);
        MatrixUtil.multiply(lambdaAInv, lambda);
        
        
        double[] r1 = MatrixUtil.multiplyMatrixByColumnVector(lambdaAInv, h1);
        double[] r2 = MatrixUtil.multiplyMatrixByColumnVector(lambdaAInv, h2);
        double[] t = MatrixUtil.multiplyMatrixByColumnVector(lambdaAInv, h3);
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
     * following Wetzstein "EE 267 Virtual Reality Course Notes: 6-DOF Pose 
     * Tracking with the VRduino"
     * estimate the extrinsic parameters
     * @param kIntr camera intrinsic parameters
     * @param kRadial an array holding the 2 radial distortion coefficients, or a null array
     * @param coordsI holds the image coordinates in pixels of features present in image i
     * @param coordsW holds the world coordinates of features present in image 1 corresponding
               to the same features and order of coordsI_i
     * @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
    f(r) = 1 +k1*r^2 + k2*r^4 if true,
    else use model #3 f(r) = 1 +k1*r + k2*r^2 if true.
    * note that if rCoeffs is null or empty, no radial distortion is applied.
     * @return 
     * @throws no.uib.cipr.matrix.NotConvergedException 
     */
    static CameraExtrinsicParameters solveForExtrinsic2(
        CameraIntrinsicParameters kIntr, double[] kRadial, 
        double[][] coordsI, double[][] coordsW, boolean useR2R4) throws NotConvergedException {
        
        if (coordsI.length != 3) {
            throw new IllegalArgumentException("coordsI length should be 3");
        }
        int n = coordsI[0].length;
        if (coordsW[0].length != n) {
            throw new IllegalArgumentException("coordsW must have same number of rows as coordsI.");
        }
        if (n < 4) {
            throw new IllegalArgumentException("coordsI[0].length cannot be less than 4");
        }
        
        // convert to camera centered coordinates
        double[][] xc = Camera.pixelToCameraCoordinates(coordsI, kRadial, 
            kIntr.getIntrinsic(), useR2R4);
        
        // normalize by last coordinate:
        for (int i = 0; i < xc[0].length; ++i) {
            xc[0][i] /= xc[2][i];
            xc[1][i] /= xc[2][i];
        }
        
        // 2*n X 8 
        // making it 2*n X 9 by adding the -u and -v columns onto end of rowsto use SVD instead of psuedoinverse
        double u, v, X, Y;
        double[][] a = new double[2*n][8];
        for (int i = 0; i < n; ++i) {
            u = xc[0][i];
            v = xc[1][i];
            X = coordsW[0][i];
            Y = coordsW[1][i];
            // eqn(10) of lecture notes   same homography as above, but missing the last column for -u and -v
            a[2*i]     = new double[]{X, Y, 1, 0, 0, 0, -u*X, -u*Y, -u};
            a[2*i + 1] = new double[]{0, 0, 0, X, Y, 1, -v*X, -v*Y, -v};
        }
        
        MatrixUtil.SVDProducts svd = MatrixUtil.performSVD(a);
        
        // vT is 8X8  last row in vT is the eigenvector for the smallest eigenvalue
        // making if 9x9 now
        double[] h = svd.vT[svd.vT.length - 1];
        
        // check the order here
        double[] hcol1 = new double[]{h[0], h[3], h[6]};
        double[] hcol2 = new double[]{h[1], h[4], h[7]};
        double[] hcol3 = new double[]{h[2], h[5], h[8]};
        
        double s = 2./(Math.sqrt(sumOfSquares(hcol1)) + Math.sqrt(sumOfSquares(hcol2)));

        // translational component of the camera pose:
        double[] t = Arrays.copyOf(hcol3, hcol3.length);
        MatrixUtil.multiply(t, s);
        t[2] *= -1;
        
        // estimate rotation from the homography
        
        double denom = sumOfSquares(hcol1);
        
        // column 1:
        double[] r1 = Arrays.copyOf(hcol1, hcol1.length);
        MatrixUtil.multiply(r1, 1./denom);
        r1[2] *= -1;
        
        //extract the second column of the rotation matrix r2 from the homography, 
        //but we have to make sure that it is orthogonal to the first column. 
        //We can enforce that as follows
        double[] r2 = Arrays.copyOf(hcol2, hcol2.length);
        r2[2] *= -1;
        
        double r1doth2 = MatrixUtil.innerProduct(r1, hcol2);
        for (int i = 0; i < r2.length; ++i) {
            r2[i] -= (r1[i] * r1doth2);
        }
        r2 = MatrixUtil.normalizeL2(r2);
        
        // check that r1 dot r2 = 0;
        double chk = MatrixUtil.innerProduct(r1, r2);
        assert(Math.abs(chk) < 1.e-3);
        
        // r3 is r1 cross r2
        double[] r3 = MatrixUtil.crossProduct(r1, r2);
        
        
        double[][] r = MatrixUtil.zeros(3, 3);
        for (int row = 0; row < 3; ++row) {
            r[row][0] = r1[row];
            r[row][1] = r2[row];
            r[row][2] = r3[row];
        }
        
        // unit norm:
        SVDProducts svd2 = MatrixUtil.performSVD(r);
        r = MatrixUtil.multiply(svd2.u, svd2.vT);
        
        CameraExtrinsicParameters kExtr = new Camera.CameraExtrinsicParameters();
        kExtr.setRotation(r);
        kExtr.setTranslation(t);
        
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
    </pre>
    @param xC distortion-free camera centered coordinates.  format is 3XN where N is the
    number of points.
    In terms of Table 1 of Chen et al. 2004, this is a double array of (x, y).
    @param k1 first radial distortion coefficient
    @param k2 second radial distortion coefficient
    @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
    f(r) = 1 +k1*r^2 + k2*r^4 if true,
    else use model #3 f(r) = 1 +k1*r + k2*r^2 if true.
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
        
        double r, r2, fr;
        //double signx, signy, c2p1, fx, fy;
        int i;
        
        for (i = 0; i < distorted[0].length; ++i) {
            r2 = distorted[0][i]*distorted[0][i] + distorted[1][i]*distorted[1][i];
            r = Math.sqrt(r2);
            //f(r) = (1 + k1*r + k2*r^2)
            if (useR2R4) {
                fr = 1 + k1*r2 + k2*r2*r2;
            } else {
                fr = 1 + k1*r + k2*r2;
            }
            // x_d = x*f_r:
            distorted[0][i] *= fr;
            distorted[1][i] *= fr;
            
            // following Ma et al. 2004 Table 2,column 3 for model #3:
            // where c = y_d/x_d = y/x
            // f(x) = (1 + k1*math.sqrt(1+c^2)*x*sign(x) + k2*(1+c^2)*x^2)
            //c2p1 = Math.pow(distorted[1][i]/distorted[0][i], 2.) + 1;
            //signx = (distorted[0][i] < 0) ? -1 : 0;
            //signy = (distorted[1][i] < 0) ? -1 : 0;
            //fx = 1 + (k1*Math.sqrt(c2p1)*distorted[0][i]*signx) + (k2*c2p1*distorted[0][i]*distorted[0][i]);
            //fy = 1 + (k1*Math.sqrt(c2p1)*distorted[1][i]*signy) + (k2*c2p1*distorted[1][i]*distorted[1][i]);
            //distorted[0][i] *= fx;
            //distorted[1][i] *= fy;
        }
                
        return distorted;
    }
    
    /**
    apply radial distortion correction to points in the camera reference frame.
    This method is used by pixelToCameraCoordinates().
    The algorithm follows Ma, Chen, & Moore 2004 to correct the distortion
    estimated as f(r) = 1 + k1*r + k2*r^2 (eqn #3 in Table 2 of Ma et al. 2004)
    In terms of the variables outlined in comments below, the algorithm input is
    distorted points as a double array of (x_d, x_d), and the radial distortion
    coefficients k1, k2.  The output is a a double array of (x, y).
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
    
    ==> To Apply Radial Distorion, given coefficients k1, k2 and coordinates
      (eqn 7) r_d = r * f(r) 
      (eqn 8) using only 2 coeffs: f(r) = (1 + k1*r + k2*r^2)
              This is equation #3 in Table 2 of Ma et al. 2004
    
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
    @param xC distorted points in the camera reference frame, presumably 
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
        from k2*r^3 + k1*r^2 + r - r_d = 0
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
        //double c2p1, signx, signy, fx, fy;
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
                                    
            // eqn (5) from Ma et al. 2004
            //u = u0 + ((ud - u0)/fr);
            //v = v0 + ((vd - v0)/fr);
            // but points are currently in image reference frame as xd, yd
            // and the radial distortion w.r.t image reference frame
            // so presumably should use eqn (4) instead
            // eqn 4): xd = x *f(r) ==> x = x_d/f(r)
            fr = 1 + k1*r + k2*r*r;
            corrected[0][i] /= fr;
            corrected[1][i] /= fr;
            
            
            // following Ma et al. 2004 Table 2,column 3 for model #3:
            // where c = y_d/x_d = y/x
            // f(x) = (1 + k1*math.sqrt(1+c^2)*x*sign(x) + k2*(1+c^2)*x^2)
            //c2p1 = Math.pow(corrected[1][i]/corrected[0][i], 2.) + 1;
            //signx = (corrected[0][i] < 0) ? -1 : 0;
            //signy = (corrected[1][i] < 0) ? -1 : 0;
            //fx = 1 + (k1*Math.sqrt(c2p1)*corrected[0][i]*signx) + (k2*c2p1*corrected[0][i]*corrected[0][i]);
            //fy = 1 + (k1*Math.sqrt(c2p1)*corrected[1][i]*signy) + (k2*c2p1*corrected[1][i]*corrected[1][i]);
            //corrected[0][i] /= fx;
            //corrected[1][i] /= fy;
        }
                
        return corrected;
    }
    
    /**
    apply radial distortion correction to points in the camera reference frame.
    This method is used by pixelToCameraCoordinates4().
    The algorithm follows Ma, Chen, & Moore 2004 to correct the distortion
    estimated as f(r) = 1 + k1*r^2 + k2*r^4 (eqn #4 in Table 2 of Ma et al. 2004)
    In terms of the variables outlined in comments below, the algorithm input is
    distorted points as a double array of (x_d, x_d), and the radial distortion
    coefficients k1, k2.  The output is a a double array of (x, y).
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
    
    ==> To Apply Radial Distorion, given coefficients k1, k2 and coordinates
      (eqn 7) r_d = r * f(r) 
      (eqn 8) using only 2 coeffs: f(r) = (1 + k1*r^2 + k2*r^4)
              This is equation #4 in Table 2 of Ma et al. 2004
    
      Ma et al. 2003 distortion model:
         x_d = x * f(r)
         y_d = y * f(r)
     
    ==> Radial Undistortion, Section 2.2 of Ma et al. 2004:
        using rf(r) and f(r) from #4 in Table 2
        f(r) = 1 + k1*r^2 + k2*r^4
        
        solving for the roots of:
        k2*r^5 + k1*r^3 + r - rd = 0.
    
        by using the companion matrix and svd
    
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
    @param xC distorted points in the camera reference frame, presumably 
    already center subtracted.  format is 3XN where N is the
    number of points.  These are (x_d, x_d) pairs in terms of Table 1 in Ma et al. 2004.
    @param k1 first radial distortion coefficient
    @param k2 second radial distortion coefficient
    @return undistorted points in the camera reference frame.  Format is 3XN where N is the
    number of points.  These are (x, y) pairs in terms of Table 1 in Ma et al. 2004.
    @throws no.uib.cipr.matrix.NotConvergedException
    */
    static double[][] removeRadialDistortion4(double[][] xC, double k1, double k2) throws NotConvergedException {
        
        if (xC.length != 3) {
            throw new IllegalArgumentException("xC.length must be 3");
        }
                
        double[][] corrected = MatrixUtil.copy(xC);
        
        //from k2*r^5 + k1*r^3 + r - rd = 0.
  
        double rd, r, fr;
        //double c2p1, fx, fy;
        double tol = 1e-5;
        int i;
        double[] roots;
        double[] coeffs = new double[]{k2, 0, k1, 0, 1, 0};
        for (i = 0; i < xC[0].length; ++i) {
            rd = Math.sqrt(corrected[0][i]*corrected[0][i] + corrected[1][i]*corrected[1][i]);
            if (Math.abs(rd) < tol) {
                continue;
            }
            coeffs[5] = -rd;
            
            roots = PolynomialRootSolver.realRoots(coeffs);
            
            if (roots == null || roots.length == 0) {
                //TODO: consider how to handle this case
                r = 0;
            } else {
                // TODO: revisit this
                r = roots[0];
            }
                                          
            // eqn (5) from Ma et al. 2004
            //u = u0 + ((ud - u0)/fr);
            //v = v0 + ((vd - v0)/fr);
            // but points are currently in image reference frame as xd, yd
            // and the radial distortion w.r.t image reference frame
            // so presumably should use eqn (4) instead
            // eqn 4): xd = x *f(r) ==> x = x_d/f(r)
            fr = 1 + k1*r*r + k2*Math.pow(r, 4.);
            corrected[0][i] /= fr;
            corrected[1][i] /= fr;
            
            // following Ma et al. 2004 Table 2,column 3 for model #3:
            // where c = y_d/x_d = y/x
            // f(x) = (1 + k1*math.sqrt(1+c^2)*x*sign(x) + k2*(1+c^2)*x^2)
            //c2p1 = Math.pow(corrected[1][i]/corrected[0][i], 2.) + 1;
            //fx = 1 + (k1*c2p1*Math.pow(corrected[0][i], 2)) + (k2*c2p1*Math.pow(corrected[0][i], 4));
            //fy = 1 + (k1*c2p1*Math.pow(corrected[1][i], 2)) + (k2*c2p1*Math.pow(corrected[1][i], 4));
            //corrected[0][i] /= fx;
            //corrected[1][i] /= fy;
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
        
        // number of features
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
     * of the next image.
               It is a 2 dimensional double array of format
               3 X (N*n) where N is the number of images.
               the first row is the x coordinates, the second row
               is the y coordinates, and the third row is "1"'s.
               The columns hold each image in order and within each image's
               columns are the features presented in the same order in each image.
               In Table 1 of Ma, Chen, & Moore 2003 "Camera Calibration"
               these are the (u_d, v_d) pairs.
     * @param u projections of the WCS feature x coordinates.  array length is n*nImages
     * @param v projections of the WCS feature y coordinates.  array length is n*nImages
     * @param cameraMatrices data structure holding the camera intrinsic parameters
     * and the extrinsic parameter matrices for each image.
     * @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
    f(r) = 1 +k1*r^2 + k2*r^4 if true,
    else use model #3 f(r) = 1 +k1*r + k2*r^2 if true.
     * @return 
     */
    static double[] solveForRadialDistortion(double[][] uvD, 
        double[] u, double[] v, 
        CameraMatrices cameraMatrices, boolean useR2R4) throws NotConvergedException {
        
        int nImages = cameraMatrices.getExtrinsics().size();
        int nFeatures = u.length/nImages;
        
        /* 
        (ud, vd) are Real observed distorted image points
        (u, v) Ideal projected undistorted image points
        [x,y,1] = A^-1 * [u, v, 1]

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
        double ui, vi, udi, vdi, xi, yi, r2, r4, r;
        double ud0 = cameraMatrices.getIntrinsics().getIntrinsic()[0][2];
        double vd0 = cameraMatrices.getIntrinsics().getIntrinsic()[1][2];
        double u0=ud0; double v0=vd0;
        double[][] xy;
        double[][] dM = new double[2*nFeatures*nImages][2];
        double[] dV = new double[2*nFeatures*nImages];
        for (i = 0; i < nImages; ++i) {
            xy = MatrixUtil.copySubMatrix(uvD, 0, 2, nFeatures*i, nFeatures*(i + 1)-1);
            xy = Camera.pixelToCameraCoordinates(xy, null, 
                cameraMatrices.getIntrinsics().getIntrinsic(), false);
            for (j = 0; j < nFeatures; ++j) {
                ui = u[nFeatures*i + j];
                vi = v[nFeatures*i + j];
                udi = uvD[0][nFeatures*i + j];
                vdi = uvD[1][nFeatures*i + j];
                xi = xy[0][j];
                yi = xy[1][j];
                r2 = xi*xi + yi*yi;
                r4 = r2*r2;
                // e.g. nFeatures=3
                //i:0 j:0          idx=0, idy=1
                //i:0 j:1          idx=2, idy=3
                //i:0 j:2          idx=4, idy=5
                //i:1 j:0  idx=6, idy=7
                //i:1 j:1  idx=8, idy=9
                if (useR2R4) {
                    dM[2*nFeatures*i + 2*j] = new double[]{(ui-u0)*r2, (ui-u0)*r4};
                    dM[2*nFeatures*i + 2*j + 1] = new double[]{(vi-v0)*r2, (vi-v0)*r4};
                } else {
                    r = Math.sqrt(r2);
                    dM[2*nFeatures*i + 2*j] = new double[]{(ui-u0)*r, (ui-u0)*r2};
                    dM[2*nFeatures*i + 2*j + 1] = new double[]{(vi-v0)*r, (vi-v0)*r2};
                }
                dV[2*nFeatures*i + 2*j] = udi -ud0 - ui;
                dV[2*nFeatures*i + 2*j + 1] = vdi -vd0 - vi;
            }
        }
           
        //k = (D^T*D)^−1*D^T*d = pseudoInv(D) * d
        double[][] dInv = MatrixUtil.pseudoinverseFullRank(dM);
        double[] k = MatrixUtil.multiplyMatrixByColumnVector(dInv, dV);
                
        return k;
    }

    static double[][] solveForHomographies(double[][] coordsI, 
        double[][] coordsW, int n, int nImages) throws NotConvergedException {
        
        double[][] h = MatrixUtil.zeros(nImages*3, 3);
        
        double[][] g, cI;
        int i, i2;
        for (i = 0; i < nImages; ++i) {
            cI = MatrixUtil.copySubMatrix(coordsI, 0, 2, n*i, n*(i + 1)-1);
            
            //3X3  and contains intrinsic camera information from the point-to-point mappings
            g = solveForHomography(cI, coordsW);
            
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
        CameraIntrinsicParameters kIntr, double[][] h, int nImages) throws NotConvergedException {
        
        List<CameraExtrinsicParameters> list = new ArrayList<CameraExtrinsicParameters>();
        //(3) for each image homography and inverse intrinsic parameter matrix,
        //    estimate the extrinsic parameters for the pose of the camera for that image.
        CameraExtrinsicParameters kExtr;
        
        int i;
        double[][] g;
        
        for (i = 0; i < nImages; ++i) {
            //                               ri,      rf,ci,cf
            g = MatrixUtil.copySubMatrix(h, 3*i, 3*i + 2, 0, 2);
            
            //NOTE: here, internal to solveForExtrinsic() would be a different place 
            // where one could remove radial distortion.
            // the method forms the image of the "absolute conic"
            kExtr = solveForExtrinsic(kIntr, g);
            
            list.add(kExtr);
        }
        
        return list;
    }
    
    /**
     * following EE 267 Virtual Reality Course Notes: 6-DOF Pose Tracking with the VRduino
       by Gordon Wetzstein
     * @param kIntr
     * @param coordsI
     * @param coordsW
     * @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
    f(r) = 1 +k1*r^2 + k2*r^4 if true,
    else use model #3 f(r) = 1 +k1*r + k2*r^2 if true.
    * note that if rCoeffs is null or empty, no radial distortion is applied.
     * @return
     * @throws NotConvergedException 
     */
    static List<CameraExtrinsicParameters> solveForExtrinsics2(
        CameraIntrinsicParameters kIntr, double[][] coordsI, double[][] coordsW,
        boolean useR2R4) 
        throws NotConvergedException {
        
        int nFeatures = coordsW[0].length;
        int nImages = coordsI[0].length/nFeatures;
        
        double[] kRadial = null;
        
        List<CameraExtrinsicParameters> list = new ArrayList<CameraExtrinsicParameters>();
        //(3) for each image homography and inverse intrinsic parameter matrix,
        //    estimate the extrinsic parameters for the pose of the camera for that image.
        CameraExtrinsicParameters kExtr;
        
        int i;
        double[][] cI;
        
        for (i = 0; i < nImages; ++i) {
            cI = MatrixUtil.copySubMatrix(coordsI, 0, 2, nFeatures*i, nFeatures*(i + 1)-1);
            
            //NOTE: here, internal to solveForExtrinsic() would be a different place 
            // where one could remove radial distortion.
            // the method forms the image of the "absolute conic"
            kExtr = solveForExtrinsic2(kIntr, kRadial, cI, coordsW, useR2R4);
            
            list.add(kExtr);
        }
        
        return list;
    }

}
