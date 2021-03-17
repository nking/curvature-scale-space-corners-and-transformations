package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
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
 * @author nichole
 */
public class Triangulation {
    
    /**
     *  create camera intrinsic matrix k with assumptions of square pixels
     * and no skew.  the focal length and optical centers should be in units of pixels.
     * NOTE that given the filed of view (FOV) and the image dimensions,
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
     * @return 
     */
    public static double[] calculateWCSPoints(
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
        int n = x1.length;
        if (x2.length != n) {
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
        double[][] camera1 = createCamera(k1, r1, t1);
        
        double[][] camera2 = createCamera(k2, r2, t2);
        
        double u1x, u1y, u2x, u2y;
        double[] tmp;
        double[][] a = new double[4*n][4];
        for (int i = 0; i < 4*n; i+=4) {
            u1x = x1[0][i];
            u1y = x1[1][i];
            u2x = x2[0][i];
            u2y = x2[1][i];
            
            //y1 * p1vec_3^T - p1vec_2^T
            tmp = Arrays.copyOf(camera1[2], 4);
            MatrixUtil.multiply(tmp, u1y);
            a[i] = MatrixUtil.subtract(tmp, camera1[1]);
            
            //p1vec_1^T - x1 * p1vec_3^T
            tmp = Arrays.copyOf(camera1[2], 4);
            MatrixUtil.multiply(tmp, u1x);
            a[i+1] = MatrixUtil.subtract(camera1[0], tmp);
            
            //y2 * p2vec_3^T - p2vec_2^T
            tmp = Arrays.copyOf(camera2[2], 4);
            MatrixUtil.multiply(tmp, u2y);
            a[i+2] = MatrixUtil.subtract(tmp, camera2[1]);
            
            //p2vec_1^T - x2 * p2vec_3^T
            tmp = Arrays.copyOf(camera2[2], 4);
            MatrixUtil.multiply(tmp, u2x);
            a[i+3] = MatrixUtil.subtract(camera2[0], tmp);
        }
        
        // A is  N X 4
        // A^T*A is 4 X 4
        //TODO: make an efficient multiplier for a^T*a in MatrixUtil:
        double[][] aTa = MatrixUtil.multiply(MatrixUtil.transpose(a), a);
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
        
        return X;
    }

    /**
     * 
     * @param k camera intrinsics matrix of size 3 x 3.
     * @param r camera extrinsics rotation matrix of size 3 x 3.
     * @param t camera extrinsics translation vector of length 3.
     * @return the camera matrix resulting from intrinsic and extrinsic parameters.
     * the size is 3 x 4.
     */
    private static double[][] createCamera(double[][] k, double[][] r, double[] t) {
        if (k.length != 3 || k[0].length != 3) {
            throw new IllegalArgumentException("k must be 3 x 3");
        }
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r must be 3 x 3");
        }
        if (t.length != 3 ) {
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
}
