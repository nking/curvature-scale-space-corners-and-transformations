package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Camera.CameraParameters;
import algorithms.imageProcessing.transform.Camera.CameraProjection;
import algorithms.matrix.MatrixUtil;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
 * given correspondence between two images (in camera reference frame, a.k.a. calibrated coordinates)
 * and given the intrinsic and extrinsic camera
 * parameters, determine the real world position.
 * 
 * useful reading:
 * <pre>
 * http://www.cs.cmu.edu/~16385/s17/Slides/11.4_Triangulation.pdf
 * add other references here
 * </pre>
 *
 * @author nichole
 */
public class Triangulation {

    public static class WCSPt {
        /**
         * real world coordinate of the observation.  [4X1]
         */
        public double[] X;

        /**
         * a scale factor between x and X: x = alpha * P * X
         * where P = K * R * [I | -t];
         * the term 1/alpha is sometimes used as lambda, the multiple of x instead.
         */
        public double alpha;

    }

    /**
     * given the intrinsic and extrinsic camera matrices for 2 images
     * and given the matching correspondence of points between the 2 images (in camera reference frame)
     * where the correspondence are observations
     * of the same object,
     * calculate the real world coordinate of the object.
     * 
     * <pre>
     * following http://www.cs.cmu.edu/~16385/s17/Slides/11.4_Triangulation.pdf
     * </pre>
     * @param k1 intrinsic camera matrix for image 1 in units of pixels.
     * @param r1 the rotation matrix for the extrinsic camera matrix for image 1.
     * @param t1 the translation vector for the extrinsic camera matrix for image 1.
     * @param k2 intrinsic camera matrix for image 2 in units of pixels.
     * @param r2 the rotation matrix for the extrinsic camera matrix for image 2.
     * @param t2 the translation vector for the extrinsic camera matrix for image 2.
     * @param x1c the camera 1 set of correspondence points.  format is 3 x N where
     * N is the number of points.  all points are observations of real world point X
     *           and should be in camera reference frame.
     * @param x2c the camera 2 set of correspondence points.  format is 3 x N where
     *      * N is the number of points.  all points are observations of real world point X
     *      *           and should be in camera reference frame.
     * @return the point coordinates in world coordinate reference frame
     */
    public static WCSPt calculateWCSPoint(
            double[][] k1, double[][] r1, double[] t1,
            double[][] k2, double[][] r2, double[] t2,
            double[][] x1c, double[][] x2c) {
        
        if (k1.length != 3 || k1[0].length != 3) {
            throw new IllegalArgumentException("k1 must be 3 x 3");
        }
        if (k2.length != 3 || k2[0].length != 3) {
            throw new IllegalArgumentException("k2 must be 3 x 3");
        }
        if (r1.length != 3 || r1[0].length != 3) {
            throw new IllegalArgumentException("r1 must be 3 x 3");
        }
        if (r2.length != 3 || r2[0].length != 3) {
            throw new IllegalArgumentException("r2 must be 3 x 3");
        }
        if (t1.length != 3 ) {
            throw new IllegalArgumentException("t1 must be length 3");
        }
        if (t2.length != 3 ) {
            throw new IllegalArgumentException("t2 must be length 3");
        }
        if (x1c.length != 3 || x2c.length != 3) {
            throw new IllegalArgumentException("x1.length must be 3 and so must x2.length");
        }
        int n = x1c[0].length;
        if (x2c[0].length != n) {
            throw new IllegalArgumentException("x1 and x2 must be same dimensions");
        }

        // 3 x 4
        double[][] camera1 = Camera.createCamera(k1, r1, t1);
        
        double[][] camera2 = Camera.createCamera(k2, r2, t2);
        
        return calculateWCSPoint(camera1, camera2, x1c, x2c);
    }
    
    /**
     given the projection matrices and matching correspondence of points between the 2 cameras
     * (measurements of the same object, that is, as a single 3D point is returned) calculate the real world
     * coordinate of the observations.
     *
     * <pre>
     *     NOTE the projection matrix is formed using P = K * [ R | t ]
     *     and x = alpha * P * X
     *     so x1 and x2 should be in camera coordinates.
     *
     * following http://www.cs.cmu.edu/~16385/s17/Slides/11.4_Triangulation.pdf
     * </pre>
     * @param camera1 camera matrix for image 1.   the size is 3X4.
     *                the data are used to construct P = K * [ R | t ].
     * @param camera2 camera matrix for image 2. the size is 3X4.
     *                the data are used to construct P = K * [ R | t ].
     * @param x1c the set of measurements of 1 real world point X from image 1 in camera coordinate reference frame.
     * The corresponding measurements of the same point in camera 2 are in x2.
     * format is 3 x N where N is the number of measurements.
     *
     * @param x2c the set of measurements of 1 real world point X from image 2 in camera coordinate reference frame.
     *      * The corresponding measurements of the same point in camera 1 are in x2.
     *      * format is 3 x N where N is the number of measurements.
     *
     * @return the 3D coordinate of the point in world scene.  note that
     * the entire array can be normalized by the last element.
     */
    public static WCSPt calculateWCSPoint(
        CameraParameters camera1, CameraParameters camera2,
        double[][] x1c, double[][] x2c) {
        
        return calculateWCSPoint(camera1.getIntrinsicParameters().getIntrinsic(),
            camera1.getExtrinsicParameters().getRotation(),
            camera1.getExtrinsicParameters().getTranslation(),
            camera2.getIntrinsicParameters().getIntrinsic(),
            camera2.getExtrinsicParameters().getRotation(),
            camera2.getExtrinsicParameters().getTranslation(),
            x1c, x2c);
    }
    
    /**
     given the projection matrices and matching correspondence of points between the 2 images
     * (measurements of the same object, that is, as a single 3D point is returned) calculate the real world
     * coordinate of the observations.
     *
     * <pre>
     *     P = K * [ R | t ], and x1 and x2 must be in camera coordinate reference frame.
     *
     * following http://www.cs.cmu.edu/~16385/s17/Slides/11.4_Triangulation.pdf
     * </pre>
     * @param camera1 camera matrix for image 1.   the size is 3X4.
     *                P = K * [ R | t ].
     * @param camera2 camera matrix for image 2. the size is 3X4.
     *                P = K * [ R | t ].
     * @param x1c the set of measurements of 1 real world point X from camera 1 in camera coordinates.
     * The corresponding measurements of the same point in camera 2 are in x2.
     * format is 3 x N where N is the number of measurements.
     * If the data are perfect, only need 1 pair of correspondence (i.e. x1[*,0] and x2[*,0]),
     * If the data are not perfect, need more than 1 pair for best fit.
     * @param x2c the set of measurements of 1 real world point X from camera 2 in camera coordinates.
     * The corresponding measurements of the same point in camera 1 are in x1.
     * format is 3 x N where N is the number of measurements.
     * If the data are perfect, only need 1 pair of correspondence (i.e. x1[*,0] and x2[*,0]),
     * If the data are not perfect, need more than 1 pair for best fit.
     * @return the 3D coordinate of the point in world scene.  note that
     * the entire array can be normalized by the last element.
     */
    public static WCSPt calculateWCSPoint(
        CameraProjection camera1, CameraProjection camera2,
        double[][] x1c, double[][] x2c) {
        
        return calculateWCSPoint(camera1.getP(), camera2.getP(), x1c, x2c);
    }
    
     /**
     * given the projection matrices and matching correspondence of points between the 2 cameras
      * (measurements of the same object, that is, as a single 3D point is returned) calculate the real world
      * coordinate of the observations.
     *
     * <pre>
      *     P = K * [ R | t ], and x1 and x2 must be in camera coordinates(pixels).
      *
      *     NOTE: sometimes the rotation is applied before translation, then P = K * [ R | -R*t].
      *     This method uses the projection matrix rows and does not decompose it into K and R components
      *     so the assumed order of transformations does not affect this method.
      *
     * following http://www.cs.cmu.edu/~16385/s17/Slides/11.4_Triangulation.pdf
     * </pre>
     * @param camera1 camera matrix for image 1.   the size is 3X4
      *                P = K * [ R | t ].
     * @param camera2 camera matrix for image 2. the size is 3X4.
      *                 P = K * [ R | t ].
     * @param x1c the set of measurements of 1 real world point X from camera 1 in camera coordinates.
     * The corresponding measurements of the same point in camera 2 are in x2.
     * format is 3 x N where N is the number of measurements.
     * If the data are perfect, only need 1 pair of correspondence (i.e. x1[*,0] and x2[*,0]),
     * If the data are not perfect, need more than 1 pair for best fit.
     * @param x2c the set of measurements of 1 real world point X from camera 2 in camera coordinates.
     * The corresponding measurements of the same point in camera 1 are in x1.
     * format is 3 x N where N is the number of measurements.
     * If the data are perfect, only need 1 pair of correspondence (i.e. x1[*,0] and x2[*,0]),
     * If the data are not perfect, need more than 1 pair for best fit.
     * @return the 3D coordinate of the point in world scene.  note that
     * the entire array can be normalized by the last element.
     * @return the 3D coordinate of the point in world scene.  note that
     * the entire array can be normalized by the last element.
     */
    public static WCSPt calculateWCSPoint(double[][] camera1, double[][] camera2,
                                          double[][] x1c, double[][] x2c) {
        
        if (camera1.length != 3 || camera1[0].length != 4) {
            throw new IllegalArgumentException("camera1 must be 3 x 4");
        }
        if (camera2.length != 3 || camera2[0].length != 4) {
            throw new IllegalArgumentException("camera2 must be 3 x 4");
        }
        if (x1c.length != 3 || x2c.length != 3) {
            throw new IllegalArgumentException("x1.length must be 3 and so must x2.length");
        }
        int n = x1c[0].length;
        if (x2c[0].length != n) {
            throw new IllegalArgumentException("x1 and x2 must be same dimensions");
        }
      
        /*
        following CMU lectures of Kris Kitani:
        http://www.cs.cmu.edu/~16385/s17/Slides/11.4_Triangulation.pdf
        
        camera matrix P = intrinsic camera matrix times extrinsic camera matrix.
        
        note that the extrinsic matrix has a translation component, and that
        translation is not a linear transformation (see Strang chap 7), so
        the translation is kept separate in most use of it to allow further operations
        to be performed that require rotation and translation to be treated separately.
        
        K = intrinsic camera matrix.
        R = rotation matrix (see euler rotation matrix).
        t = translation vector.
        I = identity matrix.  it's size 3x3 here.
        the '|' is a separation symbol in the matrix to denote that the
            content to the right of it is concatenated to the matrix as column vectors.
        
        Note that the world coordinates can be seen to go through a translation
        then a rotation represented by the extrinsic camera matrix to result in
        homogenous coordinates in the camera reference frame.
            X_c = R * (X_wcs - t)
                = R * X_wcs - R * t
        
             4x1        4x4        4x1
           [ x_c ]                  [ x_wcs ]
           [ y_c ] = [ R  -R*t ] *  [ y_wcs ]
           [ z_c ]   [ 0   1   ]    [ z_wcs ]
           [  1  ]                  [  1    ]
        
        P = K * R * [I | -t]
        
        alternatively, write as P = K * [ R | -R*t]
        
        -----------------------------------
        since data are noisy, these equalities need to be solved as best fit:
        
        x1 = P1 * X1  and  x2 = P2 * X2
            where x1 and x2 are in homogeneous coordinates
        
        similarity relation: same direction ray, but have a scale factor alpha
        
        x = alpha * P * X 
            where alpha is a scale factor and so the projection is in the same
               direction.  it's 1./(depth of point)

       NOTE: x_c = K * x_im
             X is in WCS coordinate reference frame.
             x = alpha * P * X
                 if P is [ R | t ], then x is in image coordinates.
                 if P = K * [ R | t ], then x is in camera coordinates
                  (NOTE: sometimes the rotation is first, then P = K * [ R | -R*t])

                    [ p1  p1  p3  p4  ]   [ X ]
        x = alpha * [ p5  p6  p7  p8  ] * [ Y ]
                    [ p9  p10 p11 p12 ]   [ Z ]
                                          [ 1 ]
        
        similarity relations are solved by DLT:
        (Remove scale factor, convert to linear system and solve with SVD)
        
           let pvec_1^T = [ p1 p2 p3 p4 ], pvec_2^T = [ p5 p6 p7 p8 ] etc. 
        
                       1x4               4x1
                    [ --pvec_1^T-- ]   [ X ]
        x = alpha * [ --pvec_2^T-- ] * [ Y ]
                    [ --pvec_3^T-- ]   [ Z ]
                                       [ 1 ]
        
                    [ pvec_1^T * [X,Y,Z,1] ] <-- each row result is 1x1
        x = alpha * [ pvec_2^T * [X,Y,Z,1] ]
                    [ pvec_3^T * [X,Y,Z,1] ]
        
           let Xvec = [ X Y Z 1 ] as a row
        
                    [ pvec_1^T * Xvec ]
        x = alpha * [ pvec_2^T * Xvec ]
                    [ pvec_3^T * Xvec ]
        
        NOTE: The cross product of 2 vectors of the same direction is 0.
        
            So we have 'x cross (P * X) = 0' and alpha drops out

                        [ a2*b3 - a3*b2 ]
                a x b = [ a3*b1 - a1*b3 ]
                        [ a1*b2 - a2*b1 ]

        Can rewrite in terms of cross product:
       
        [ x ]       [ pvec_1^T * Xvec ]   [ y * pvec_3^T * Xvec - pvec_2^T * Xvec    ]   [ 0 ]
        [ y ] cross [ pvec_2^T * Xvec ] = [ pvec_1^T * Xvec - x * pvec_3^T * Xvec    ] = [ 0 ]
        [ 1 ]       [ pvec_3^T * Xvec ]   [ x * pvec_2^T * Xvec - y * pvec_1^T * Xvec]   [ 0 ]
        
        The 3rd line is a linear combination of the first and second lines. (x times the first line plus y times the second line)
        so remove it:

        [ y * pvec_3^T * Xvec - pvec_2^T * Xvec ]   [ 0 ]
        [ pvec_1^T * Xvec - x * pvec_3^T * Xvec ] = [ 0 ]
        
        This is in format A_i * Xvec = 0
        
        can concatenate the 2 rows for the 2nd image in A_i:
        
        [  y1 * p1vec_3^T - p1vec_2^T ]           [ 0 ]
        [ -x1 * p1vec_3^T + p1vec_1^T ] * Xvec  = [ 0 ]
        [  y2 * p2vec_3^T - p2vec_2^T ]           [ 0 ]
        [ -x2 * p2vec_3^T + p2vec_1^T ]           [ 0 ]
        
        solve Xvec in A * Xvec = 0 by minimizing ||A*x||^2 subject to ||x||^2 = 1
        
        Total least squares
        
        Solution is the eigenvector corresponding to smallest eigenvalue of A^T*A.
        */
                
        double u1x, u1y, u2x, u2y;
        double[] tmp;
        double[][] a = new double[4*n][4];
        int i, j;
        for (i = 0, j = 0; i < n; ++i, j+=4) {
            u1x = x1c[0][i];
            u1y = x1c[1][i];
            u2x = x2c[0][i];
            u2y = x2c[1][i];
                        
            //y1 * p1vec[2]^T - p1vec[1]^T
            tmp = Arrays.copyOf(camera1[2], 4);
            MatrixUtil.multiply(tmp, u1y);
            a[j] = MatrixUtil.subtract(tmp, camera1[1]);
            
            //p1vec[0]^T - x1 * p1vec[2]^T
            tmp = Arrays.copyOf(camera1[2], 4);
            MatrixUtil.multiply(tmp, u1x);
            a[j+1] = MatrixUtil.subtract(camera1[0], tmp);
            
            //y2 * p2vec[2]^T - p2vec[1]^T
            tmp = Arrays.copyOf(camera2[2], 4);
            MatrixUtil.multiply(tmp, u2y);
            a[j+2] = MatrixUtil.subtract(tmp, camera2[1]);
            
            //p2vec[0]^T - x2 * p2vec[2]^T
            tmp = Arrays.copyOf(camera2[2], 4);
            MatrixUtil.multiply(tmp, u2x);
            a[j+3] = MatrixUtil.subtract(camera2[0], tmp);
        }
        
        // A is  4*N X 4
        // A^T*A is 4 X 4
        double[][] aTa = MatrixUtil.createATransposedTimesA(a);
        assert(aTa.length == 4);
        assert(aTa[0].length == 4);
        
        //NOTE: SVD(A).V is the same as SVD(A^TA).V
        
        SVD svd = null;
        try {
            svd = SVD.factorize(new DenseMatrix(aTa));
        } catch (NotConvergedException ex) {
            Logger.getLogger(Triangulation.class.getName()).log(Level.SEVERE, null, ex);
            return null;
        }
        
        double[][] vT = MatrixUtil.convertToRowMajor(svd.getVt());
        assert(vT.length == 4);
        assert(vT[0].length == 4);
        
        // eigenvector corresponding to smallest eigenvector is last row in svd.V^T
        double[] X = Arrays.copyOf(vT[vT.length - 1], vT[0].length);
        
        /*
        System.out.printf("x1=\n%s\n", FormatArray.toString(x1, "%.4e"));
        System.out.printf("camera1=\n%s\n", FormatArray.toString(camera1, "%.4e"));
        System.out.printf("x2=\n%s\n", FormatArray.toString(x2, "%.4e"));
        System.out.printf("camera2=\n%s\n", FormatArray.toString(camera2, "%.4e"));
        System.out.printf("X=\n%s\n\n", FormatArray.toString(X, "%.4e"));
        */
        
        // can see that the constraint ||X||^2 = 1 is preserved
                
        double[] x1Rev = MatrixUtil.multiplyMatrixByColumnVector(camera1, X);
        double[] x2Rev = MatrixUtil.multiplyMatrixByColumnVector(camera2, X);        
        double alpha = ((1./x1Rev[2]) + (1./x2Rev[2]))/2.;
        
        MatrixUtil.multiply(x1Rev, 1./x1Rev[2]);
        MatrixUtil.multiply(x2Rev, 1./x2Rev[2]);
        
        /*
        System.out.printf("x1Rev=\n%s\n", FormatArray.toString(x1Rev, "%.4e"));
        System.out.printf("x1=\n%s\n", FormatArray.toString(MatrixUtil.extractColumn(x1, 0), "%.4e"));
        System.out.printf("x2Rev=\n%s\n", FormatArray.toString(x2Rev, "%.4e"));        
        System.out.printf("x2=\n%s\n", FormatArray.toString(MatrixUtil.extractColumn(x2, 0), "%.4e"));
        */

        //System.out.printf("alpha=\n%.3e\n", alpha);
        //MatrixUtil.multiply(X, alpha);

        WCSPt w = new WCSPt();
        w.X = X;
        w.alpha = alpha;

        return w;
    }

    /*
     * not finished
     * 
     * given a feature in 2 cameras and the rotation and translation between
     * the cameras, estimate the depths and universal scale factor to 
     * return the estimates of the 3-D position.
     <pre>
     The algorithm follows Serge Belongie lectures from Computer Vision II, CSE 252B, USSD
     who refers to Ma, Soatto, Kosecka, and Sastry 2003
     "An Invitation to 3D Vision From Images to Geometric Models" (Chap 5)
     </pre>
     * @param r rotation of camera 2 with respect to camera 1
     * @param t the translation of camera 2 with respect to camera 1
     * @param x1 coordinates of feature in image 1.   the array should be length
     * @param x2 coordinates of feature in image 1.   the array should be length
     * @return the estimated position of the 3-D point as 2 estimates which should be the same.
     * The depths and universal scale factors are returned also.
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    /*
    public static WCSResults calculateDepths(double[][] r, double[] t,
        double[] x1, double[] x2) throws NotConvergedException {
        
        // X_2 = R * X_1 + gamma * T where gamma is the universal scale factor
        //
        // lambda_2 * x_2 = lambda_1 * R * x_1 + gamma * T
        //
        // multiply both sides by skew symmetric of x_2 = [x_2]_x to use the property that 
        //     [x_2]_x * x_2 = 0 (i.e. cross product is 0).
        //
        //  [x_2]_x * lambda_2 * x_2 = [x_2]_x * lambda_1 * R * x_1 + [x_2]_x * gamma * T
        //  0 = [x_2]_x * lambda_1 * R * x_1 + [x_2]_x * gamma * T
        //
        // then [ [x_2]_x * R * x_1   [x_2]_x * T ] * [ lambda_1 ] = 0
        //                                            [   gamma  ]
        //          sizes are [ 3X2 ] * [ 2X1 ] = [ 3X1 ]
        //
        // let M = the matrix on right hand side.
        //    this assumes no noise in data
        //  then [lambda_1, gamma] = SVD(M).V^T[last row]
        //  
        //  then X_1 = lambda_1 * x_1
        // 
        //  Then back to the transformation by extrinsic parameters:
        //      lambda_2 * x_2 = lambda_1 * R * x_1 + gamma * T
        //  multiply both sides by skew symmetric of x_1 = [x_1]_x
        //      [x_1]_x * lambda_2 * x_2   =   [x_1]_x * lambda_1 * R * x_1   +   [x_1]_x * gamma * T
        //      [x_1]_x * lambda_2 * x_2   =    0  +   [x_1]_x * gamma * T
        //      [x_1]_x * lambda_2 * x_2  -  [x_1]_x * gamma * T = 0
        //     where gamma is known
        //  
        // then [ [x_1]_x * x_2   -[x_1]_x * gamma * T ] * [ lambda_2 ] = 0
        //                                                 [   1  ]
        //       sizes are [ 3X2 ] * [ 2X1 ] = [ 3X1 ]
        // then [lambda_2, 1] = SVD(M).V^T[last row]
        // X_2 = lambda_2 * x_2
        // 
        // and assert that X_2 = R * X_1 + gamma * T
        
        double[][] x1SkewSym = MatrixUtil.skewSymmetric(x1);
        double[][] x2SkewSym = MatrixUtil.skewSymmetric(x2);
        
        double[] M1Col0 = MatrixUtil.multiplyMatrixByColumnVector(
            MatrixUtil.multiply(x2SkewSym, r), x1
        );
        double[] M1Col1 = MatrixUtil.multiplyMatrixByColumnVector(x2SkewSym, t);
        double[][] M = new double[3][2];
        int i;
        for (i = 0; i < 3; ++i) {
            M[i] = new double[]{M1Col0[i], M1Col1[i]};
        }
        
        MatrixUtil.SVDProducts svd = MatrixUtil.performSVD(M);
        double lambda1 = svd.vT[svd.vT.length - 1][0];
        double gamma = svd.vT[svd.vT.length - 1][1];
        
        double[] X1 = Arrays.copyOf(x2, x2.length);
        MatrixUtil.multiply(X1, lambda1);
        
        // [ [x_1]_x * x_2   -[x_1]_x * gamma * T ]
        M1Col0 = MatrixUtil.multiplyMatrixByColumnVector(x1SkewSym, x2);
        M1Col1 = MatrixUtil.multiplyMatrixByColumnVector(x1SkewSym, t);
        MatrixUtil.multiply(M1Col1, -gamma);
        for (i = 0; i < 3; ++i) {
            M[i] = new double[]{M1Col0[i], M1Col1[i]};
        }
        svd = MatrixUtil.performSVD(M);
        double lambda2 = svd.vT[svd.vT.length - 1][0];
        double one = svd.vT[svd.vT.length - 1][1];
        assert(Math.abs(one) - 1 < 1e-2);
        
        double[] X2 = Arrays.copyOf(x2, x2.length);
        MatrixUtil.multiply(X2, lambda1);
        
        //assert that X_2 = R * X_1 + gamma * T
        double[] gt = Arrays.copyOf(t, t.length);
        MatrixUtil.multiply(gt, gamma);
        double[] checkX2 = MatrixUtil.multiplyMatrixByColumnVector(r, X1);
        checkX2 = MatrixUtil.add(checkX2, gt);
        
        for (i = 0; i < X2.length; ++i) {
            assert(Math.abs(X2[i] - checkX2[i]) < 1.e-2);
        }
        
        WCSResults w = new WCSResults();
        w.setX1(X1);
        w.setX2(X2);
        w.setDepth1(lambda1);
        w.setDepth2(lambda2);
        w.setUniversalScaleFactor(gamma);
        
        return w;
    }
    */
    
    public static class WCSResults {
        /**
         * the 3-D position of point x1 (which should be the same as X2)
         */
        private double[] X1;
        /**
         * the 3-D position of point x2 (which should be the same as X1)
         */
        private double[] X2;
        /**
         * the depth of X1
         */
        private double depth1;
        /**
         * the depth of X2
         */
        private double depth2;
        /**
         * the universal scale factor which is applied to the translation between camera positions
         */
        private double universalScaleFactor;

        /**
         * @return the X1
         */
        public double[] getX1() {
            return X1;
        }

        /**
         * @param X1 the X1 to set
         */
        public void setX1(double[] X1) {
            this.X1 = X1;
        }

        /**
         * @return the X2
         */
        public double[] getX2() {
            return X2;
        }

        /**
         * @param X2 the X2 to set
         */
        public void setX2(double[] X2) {
            this.X2 = X2;
        }

        /**
         * @return the depth1
         */
        public double getDepth1() {
            return depth1;
        }

        /**
         * @param depth1 the depth1 to set
         */
        public void setDepth1(double depth1) {
            this.depth1 = depth1;
        }

        /**
         * @return the depth2
         */
        public double getDepth2() {
            return depth2;
        }

        /**
         * @param depth2 the depth2 to set
         */
        public void setDepth2(double depth2) {
            this.depth2 = depth2;
        }

        /**
         * @return the universalScaleFactor
         */
        public double getUniversalScaleFactor() {
            return universalScaleFactor;
        }

        /**
         * @param universalScaleFactor the universalScaleFactor to set
         */
        public void setUniversalScaleFactor(double universalScaleFactor) {
            this.universalScaleFactor = universalScaleFactor;
        }
        
    }
}
