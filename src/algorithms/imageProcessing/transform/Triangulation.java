package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
 * given correspondence between two images and the camera
 * parameters as intrinsic and extrinsic parameters,
 * determine the real world position.
 * 
 * useful reading:
 * <pre>
 * http://www.cs.cmu.edu/~16385/s17/Slides/11.4_Triangulation.pdf
 * add other references here
 * </pre>
 * 
 * TODO: implement a method using minimization of the re-projection error by
 * non-linear optimization such as Levenberg-Marquardt.
 * 
 * TODO: construct a method that might be better placed in reconstruction:
 * given camera intrinsic parameters and 2-view correspondences: first calculate
 * the essential matrix, then extract the cameras extrinsic parameters from
 * them, then use triangulation to get the WCS coordinates.
 * @author nichole
 */
public class Triangulation {
    
    /**
     * given the camera matrix as intrinsic and extrinsic matrices for 2 images
     * and given the matching correspondence of points between the 2 images,
     * calculate the real world coordinate of the observations.
     * 
     * <pre>
     * following http://www.cs.cmu.edu/~16385/s17/Slides/11.4_Triangulation.pdf
     * add references here
     * </pre>
     * @param k1 intrinsic camera matrix for image 1 in units of pixels.
     * @param r1 the rotation matrix for the extrinsic camera matrix for image 1.
     * @param t1 the translation vector for the extrinsic camera matrix for image 1.
     * @param k2 intrinsic camera matrix for image 2 in units of pixels.
     * @param r2 the rotation matrix for the extrinsic camera matrix for image 2.
     * @param t2 the translation vector for the extrinsic camera matrix for image 2.
     * @param x1 the image 1 set of correspondence points.  format is 3 x N where
     * N is the number of points.  all points are observations of real world point X.
     * @param x2 the image 2 set of correspondence points.  format is 3 x N where
     * N is the number of points. .  all points are observations of real world point X.
     * @return the point coordinates in world coordinate reference frame
     */
    public static double[] calculateWCSPoint(
        double[][] k1, double[][] r1, double[] t1,
        double[][] k2, double[][] r2, double[] t2,
        double[][] x1, double[][] x2) {
        
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
        if (x1.length != 3 || x2.length != 3) {
            throw new IllegalArgumentException("x1.length must be 3 and so must x2.length");
        }
        int n = x1[0].length;
        if (x2[0].length != n) {
            throw new IllegalArgumentException("x1 and x2 must be same dimensions");
        }
        
        /*
        following CMU lectures of Kris Kitani:
        http://www.cs.cmu.edu/~16385/s17/Slides/11.4_Triangulation.pdf
        
        camera matrix P = intrinsic camera matrix times extrinsic camera matrix.
        
        note that the extrinsic matrix has a translation component, and that
        translation is not a linear transformation (see Strang chap 7), so
        the translation is kept seperate in most use of it to allow further operations
        to be performed that require rotation and translation to be treated separately.
        
        K = intrinsic camera matrix.
        R = rotation matrix (see euler roration matrix).
        t = translation vector.
        I = identity matrix.  it's size 3x3 here.
        the '|' is a seperation symbol in the matrix to denotate that the 
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
        
        or, can write as P = K * [ R | -R*t]
        
        -----------------------------------
        since data are noisy, these equalities need to be solved as best fit:
        
        x1 = P1 * X1  and  x2 = P2 * X2
            where x1 and x2 are in homogeneous coordinates
        
        x = alpha * P * X 
            where alpha is a scale factor and so the projection is in the same
               direction.
        
                    [ p1  p1  p3  p4  ]   [ X ]
        x = alpha * [ p5  p6  p7  p8  ] * [ Y ]
                    [ p9  p10 p11 p12 ]   [ Z ]
                                          [ 1 ]
        
           let pvec_1^T = [ p1 p2 p3 p4 ], etc. 
        
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
        
            So we have x cross P * X = 0 and alpha drops out

                        [ a2*b3 - a3*b2 ]
                a x b = [ a3*b1 - a1*b3 ]
                        [ a1*b2 - a2*b1 ]

        Can rewrite in terms of cross poduct:
       
        [ x ]       [ pvec_1^T * Xvec ]   [ y * pvec_3^T * Xvec - pvec_2^T * Xvec    ]   [ 0 ]
        [ y ] cross [ pvec_2^T * Xvec ] = [ pvec_1^T * Xvec - x * pvec_3^T * Xvec    ] = [ 0 ]
        [ 1 ]       [ pvec_3^T * Xvec ]   [ x * pvec_2^T * Xvec - y * pvec_1^T * Xvec]   [ 0 ]
        
        The 3rd line is a linear combination of the first and second lines. (x times the first line plus y times the second line)

        [ y * pvec_3^T * Xvec - pvec_2^T * Xvec ]   [ 0 ]
        [ pvec_1^T * Xvec - x * pvec_3^T * Xvec ] = [ 0 ]
        
        This is in format A_i * Xvec = 0
        
        can concatenate the 2 rows for the 2nd image in A_i:
        
        [ y1 * p1vec_3^T - p1vec_2^T ]           [ 0 ]
        [ p1vec_1^T - x1 * p1vec_3^T ] * Xvec  = [ 0 ]
        [ y2 * p2vec_3^T - p2vec_2^T ]           [ 0 ]
        [ p2vec_1^T - x2 * p2vec_3^T ]           [ 0 ]
        
        solve Xvec in A * Xvec = 0 by minimizing ||A*x||^2 subject to ||x||^2 = 1
        
        Solution is the eigenvector corresponding to smallest eigenvalue of A^T*A.
        
        */
        
        // 3 x 4
        double[][] camera1 = Camera.createCamera(k1, r1, t1);
        
        double[][] camera2 = Camera.createCamera(k2, r2, t2);
        
        return calculateWCSPoint(camera1, camera2, x1, x2);
    }
    
     /**
     * given the camera matrix as intrinsic and extrinsic matrices for 2 images
     * and given the matching correspondence of points between the 2 images,
     * calculate the real world coordinate of the observations.
     * 
     * <pre>
     * following http://www.cs.cmu.edu/~16385/s17/Slides/11.4_Triangulation.pdf
     * add references here
     * </pre>
     * @param camera1 camera matrix for image 1 in units of pixels.  
     * It has intrinsic and extrinsic components.
     * @param camera2 camera matrix for image 2 in units of pixels.  
     * It has intrinsic and extrinsic components.
     * @param x1 the image 1 set of measurements of real world point X.
     * The corresponding measurements of the same point in image 2 are in x2.
     * format is 3 x N where
     * N is the number of measurements.
     * @param x2 the image 2 set of measurements of real world point X.
     * The corresponding measurements of the same point in image 1 are in x1.
     * format is 3 x N where
     * N is the number of measurements.
     */
    public static double[] calculateWCSPoint(
        double[][] camera1, double[][] camera2,
        double[][] x1, double[][] x2) {
        
        if (camera1.length != 3 || camera1[0].length != 4) {
            throw new IllegalArgumentException("camera1 must be 3 x 4");
        }
        if (camera2.length != 3 || camera2[0].length != 4) {
            throw new IllegalArgumentException("camera2 must be 3 x 4");
        }
        if (x1.length != 3 || x2.length != 3) {
            throw new IllegalArgumentException("x1.length must be 3 and so must x2.length");
        }
        int n = x1[0].length;
        if (x2[0].length != n) {
            throw new IllegalArgumentException("x1 and x2 must be same dimensions");
        }
        
        /*
        following CMU lectures of Kris Kitani:
        http://www.cs.cmu.edu/~16385/s17/Slides/11.4_Triangulation.pdf
        
        camera matrix P = intrinsic camera matrix times extrinsic camera matrix.
        
        note that the extrinsic matrix has a translation component, and that
        translation is not a linear transformation (see Strang chap 7), so
        the translation is kept seperate in most use of it to allow further operations
        to be performed that require rotation and translation to be treated separately.
        
        K = intrinsic camera matrix.
        R = rotation matrix (see euler roration matrix).
        t = translation vector.
        I = identity matrix.  it's size 3x3 here.
        the '|' is a seperation symbol in the matrix to denotate that the 
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
        
        alternately, can write as P = K * [ R | -R*t]
        
        -----------------------------------
        since data are noisy, these equalities need to be solved as best fit:
        
        x1 = P1 * X1  and  x2 = P2 * X2
            where x1 and x2 are in homogeneous coordinates
        
        x = alpha * P * X 
            where alpha is a scale factor and so the projection is in the same
               direction.  it's 1./(depth of point)
        
                    [ p1  p1  p3  p4  ]   [ X ]
        x = alpha * [ p5  p6  p7  p8  ] * [ Y ]
                    [ p9  p10 p11 p12 ]   [ Z ]
                                          [ 1 ]
        
           let pvec_1^T = [ p1 p2 p3 p4 ], etc. 
        
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
        
            So we have x cross P * X = 0 and alpha drops out

                        [ a2*b3 - a3*b2 ]
                a x b = [ a3*b1 - a1*b3 ]
                        [ a1*b2 - a2*b1 ]

        Can rewrite in terms of cross poduct:
       
        [ x ]       [ pvec_1^T * Xvec ]   [ y * pvec_3^T * Xvec - pvec_2^T * Xvec    ]   [ 0 ]
        [ y ] cross [ pvec_2^T * Xvec ] = [ pvec_1^T * Xvec - x * pvec_3^T * Xvec    ] = [ 0 ]
        [ 1 ]       [ pvec_3^T * Xvec ]   [ x * pvec_2^T * Xvec - y * pvec_1^T * Xvec]   [ 0 ]
        
        The 3rd line is a linear combination of the first and second lines. (x times the first line plus y times the second line)

        [ y * pvec_3^T * Xvec - pvec_2^T * Xvec ]   [ 0 ]
        [ pvec_1^T * Xvec - x * pvec_3^T * Xvec ] = [ 0 ]
        
        This is in format A_i * Xvec = 0
        
        can concatenate the 2 rows for the 2nd image in A_i:
        
        [ y1 * p1vec_3^T - p1vec_2^T ]           [ 0 ]
        [ p1vec_1^T - x1 * p1vec_3^T ] * Xvec  = [ 0 ]
        [ y2 * p2vec_3^T - p2vec_2^T ]           [ 0 ]
        [ p2vec_1^T - x2 * p2vec_3^T ]           [ 0 ]
        
        solve Xvec in A * Xvec = 0 by minimizing ||A*x||^2 subject to ||x||^2 = 1
        
        Solution is the eigenvector corresponding to smallest eigenvalue of A^T*A.
        
        */
        
        double u1x, u1y, u2x, u2y;
        double[] tmp;
        double[][] a = new double[4*n][4];
        int i, j;
        for (i = 0, j = 0; i < n; ++i, j+=4) {
            u1x = x1[0][i];
            u1y = x1[1][i];
            u2x = x2[0][i];
            u2y = x2[1][i];
                        
            //y1 * p1vec_3^T - p1vec_2^T
            tmp = Arrays.copyOf(camera1[2], 4);
            MatrixUtil.multiply(tmp, u1y);
            a[j] = MatrixUtil.subtract(tmp, camera1[1]);
            
            //p1vec_1^T - x1 * p1vec_3^T
            tmp = Arrays.copyOf(camera1[2], 4);
            MatrixUtil.multiply(tmp, u1x);
            a[j+1] = MatrixUtil.subtract(camera1[0], tmp);
            
            //y2 * p2vec_3^T - p2vec_2^T
            tmp = Arrays.copyOf(camera2[2], 4);
            MatrixUtil.multiply(tmp, u2y);
            a[j+2] = MatrixUtil.subtract(tmp, camera2[1]);
            
            //p2vec_1^T - x2 * p2vec_3^T
            tmp = Arrays.copyOf(camera2[2], 4);
            MatrixUtil.multiply(tmp, u2x);
            a[j+3] = MatrixUtil.subtract(camera2[0], tmp);
        }
        
        // A is  4*N X 4
        // A^T*A is 4 X 4
        double[][] aTa = MatrixUtil.createATransposedTimesA(a);
        assert(aTa.length == 4);
        assert(aTa[0].length == 4);
        
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
        

        System.out.printf("x1=\n%s\n", FormatArray.toString(x1, "%.4e"));
        System.out.printf("camera1=\n%s\n", FormatArray.toString(camera1, "%.4e"));
        System.out.printf("x2=\n%s\n", FormatArray.toString(x2, "%.4e"));
        System.out.printf("camera2=\n%s\n", FormatArray.toString(camera2, "%.4e"));
        System.out.printf("X=\n%s\n\n", FormatArray.toString(X, "%.4e"));
        
        // can see that the constraint ||X||^2 = 1 is preserved
                
        double[] x1Rev = MatrixUtil.multiplyMatrixByColumnVector(camera1, X);
        double[] x2Rev = MatrixUtil.multiplyMatrixByColumnVector(camera2, X);        
        double alpha = ((1./x1Rev[2]) + (1./x2Rev[2]))/2.;
        
        MatrixUtil.multiply(x1Rev, 1./x1Rev[2]);
        MatrixUtil.multiply(x2Rev, 1./x2Rev[2]);
        
        System.out.printf("x1Rev=\n%s\n", FormatArray.toString(x1Rev, "%.4e"));
        System.out.printf("x1=\n%s\n", FormatArray.toString(MatrixUtil.extractColumn(x1, 0), "%.4e"));
        System.out.printf("x2Rev=\n%s\n", FormatArray.toString(x2Rev, "%.4e"));        
        System.out.printf("x2=\n%s\n", FormatArray.toString(MatrixUtil.extractColumn(x2, 0), "%.4e"));
        //System.out.printf("alpha=\n%.3e\n", alpha);
        //MatrixUtil.multiply(X, alpha);
        
        
        
        return X;
    }
}
