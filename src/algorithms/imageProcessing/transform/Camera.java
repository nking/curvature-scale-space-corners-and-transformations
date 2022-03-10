package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import no.uib.cipr.matrix.NotConvergedException;

/**
 * utility methods for camera intrinsic and extrinsic matrices.
 * 
 * KIntrinsic is the camera intrinsic parameters
 * KExtrinsic is the camera extrinsic parameters, specifically, the camera
 * position and orientation in world coordinates.
 * 
 * to project the world coordinate system point into the camera frame
 * canonical form,
 * multiply it by the extrinsic matrix (which will translate the
 * coordinates and rotate them).
 * 
 * to project the camera frame coordinate system into the image frame,
 * multiply it by the intrinsic matrix.
 * 
 * TODO: consider implementing tangential distortion.  
 * see Heikkila, Janne, and Olli Silven 1997. "A four-step camera calibration 
 * procedure with implicit image correction." 
 * Computer Vision and Pattern Recognition, 1997. 
 * Proceedings., 1997 IEEE Computer Society Conference on. IEEE, 1997.
 * @author nichole
 */
public class Camera {
   
    public static double eps = 1e-7;
    
    /**
     *  create camera intrinsic matrix k with assumptions of square pixels
     * and no skew.  the focal length and optical centers should be in units of pixels.
     * NOTE that given the field of view (FOV) and the image dimensions,
     * one can roughly estimate the focal length as (image width/2) / tan(FOV/2).
     * @param focalLength focal length of camera in units of pixels.  NOTE that
     * the sign of the focal length remains the same as is given, so if you
     * want to use a right-hand coordinate system, you should give a negative
     * focal length.
     * @param centerX x coordinate of principal point in pixels, usually image center.
     * @param centerY y coordinate of principal point in pixels, usually image center.
     * @return intrinsic camera matrix in units of pixels.
     */
    public static double[][] createIntrinsicCameraMatrix(double focalLength,
        double centerX, double centerY) {
        
        double[][] k = new double[3][3];
        k[0] = new double[]{focalLength, 0, centerX};
        k[1] = new double[]{0, focalLength, centerY};
        k[2]= new double[]{0, 0, 1};
        
        return k;
    }
    
    /**
     *  create camera intrinsic matrix k.  the focal length and optical centers should be in units of pixels.
     * NOTE that given the field of view (FOV) and the image dimensions,
     * one can roughly estimate the focal length as (image width/2) / tan(FOV/2).
     * @param focalLengthX focal length of camera in units of pixels along x axis.
     * NOTE that
     * the sign of the focal length remains the same as is given, so if you
     * want to use a right-hand coordinate system, you should give a negative
     * focal length.
     * @param focalLengthY focal length of camera in units of pixels along y axis.
     * NOTE that
     * the sign of the focal length remains the same as is given, so if you
     * want to use a right-hand coordinate system, you should give a negative
     * focal length.
     * @param centerX x coordinate of principal point in pixels, usually image center.
     * @param centerY y coordinate of principal point in pixels, usually image center.
     * @param skew camera skew
     * @return intrinsic camera matrix in units of pixels.
     */
    public static double[][] createIntrinsicCameraMatrix(double focalLengthX,
        double focalLengthY, double centerX, double centerY, double skew) {
        
        double[][] k = new double[3][3];
        k[0] = new double[]{focalLengthX, skew, centerX};
        k[1] = new double[]{0, focalLengthY, centerY};
        k[2]= new double[]{0, 0, 1};
        
        return k;
    }
    
    /**
     *create the inverse of camera intrinsic matrix k with assumptions of square pixels
     * and no skew.  the focal length and optical centers should be in units of pixels.
     * NOTE that given the field of view (FOV) and the image dimensions,
     * one can roughly estimate the focal length as (image width/2) / tan(FOV/2).
     * @param focalLength focal length of camera in units of pixels.
     * NOTE that
     * the sign of the focal length remains the same as is given, so if you
     * want to use a right-hand coordinate system, you should give a negative
     * focal length.
     * @param centerX x coordinate of principal point in pixels, usually image center.
     * @param centerY y coordinate of principal point in pixels, usually image center.
     * @return intrinsic camera matrix in units of pixels.
     */
    public static double[][] createIntrinsicCameraMatrixInverse(double focalLength,
        double centerX, double centerY) {
        
        double[][] k = new double[3][3];
        k[0] = new double[]{1./focalLength, 0, -centerX/focalLength};
        k[1] = new double[]{0, 1./focalLength, -centerY/focalLength};
        k[2] = new double[]{0, 0, 1};
        
        return k;
    }
    
    /**
     *  create the inverse of camera intrinsic matrix k.
     * @param kIntr intrinsic camera matrix.
     * @return intrinsic camera matrix inverse.
     */
    public static double[][] createIntrinsicCameraMatrixInverse(double[][] kIntr) {
        
        /*
        translation matrix: inverse changes the signs of the translation elements, but not the diagonal.
        rotation matrix: inverse is the transpose of rotation matrix.
        scaling matrix: inverse is performed on each element, that is, the reciprocal.
        
        skew is a translation term, so will put it in its own matrix
        
             | 1  0  xc |   | fx  0   0 |   | 1  g/fx   0 | 
         K = | 0  1  yc | * | 0  fy   0 | * | 0  1      0 | 
             | 0  0   1 |   | 0   0   1 |   | 0  0      1 | 
            
             | fx  0  xc |   | 1  g/fx   0 |   | fx  g  xc |
           = | 0  fy  yc | * | 0  1      0 | = | 0  fy  yc |
             | 0  0   1  |   | 0  0      1 |   | 0  0   1  |
        
        let z = the inverse term for g/fx = -g/fx
        
               | 1  z       0 |   | 1/fx  0      0   |   | 1  0 -xc | 
        K^-1 = | 0  1       0 | * | 0     1/fy   0   | * | 0  1 -yc | 
               | 0  0       1 |   | 0     0      1/1 |   | 0  0   1 |
        
               | 1/fx   z/fy  0 |   | 1  0 -xc |   | 1/fx   z/fy    -xc/fx -yc*z/fy |
             = |  0     1/fy  0 | * | 0  1 -yc | = |  0     1/fy    -yc/fy          |
               |  0     0     1 |   | 0  0   1 |   |  0     0       1               |
        
               | 1/fx   -g/(fx*fy)   -xc/fx + g*yc/(fx*fy) |
             = | 0       1/fy        -yc/fy                |
               | 0       0           1                     |
        
        validation here:
        https://www.imatest.com/support/docs/pre-5-2/geometric-calibration-deprecated/projective-camera/
        */
        
        double[][] kInv = MatrixUtil.zeros(3, 3);
        
                
        double fx = kIntr[0][0];
        double fy = kIntr[1][1];
        double fxfy = fx*fy;
        double g = kIntr[0][1];
        double xc = kIntr[0][2];
        double yc = kIntr[1][2];
        
        kInv[0][0] = 1./fx;
        kInv[0][1] = -g/fxfy;
        kInv[0][2] = (-xc/fx) + (g*yc/fxfy);
        kInv[1][1] = 1./fy;
        kInv[1][2] = -yc/fy;
        kInv[2][2] = 1;
        
        return kInv;
    }
   
    /**
     * 
     * @param k camera intrinsics matrix of size 3 x 3.
     * @param r camera extrinsic rotation matrix of size 3 x 3.
     * @param t camera extrinsic translation vector of length 2.
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
     * @param r camera extrinsic rotation matrix of size 3 x 3.
     * @param t camera extrinsic translation vector of length 2.
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
     * The radial distortion functions to choose are from Ma et al. 2004 Table 2.
        The choices are model #4, f(r) = 1 +k1*r^2 + k2*r^4,
        or model #3, f(r) = 1 +k1*r + k2*r^2.
        Note that if rCoeffs is null or empty, no radial distortion is applied.
     * Also useful reading is NVM Tools by Alex Locher
    https://github.com/alexlocher/nvmtools.git
     * @param xC distortion-free camera centered coordinates. 
     * format is 3XN for N points.  
     * In terms of Table 1 of Ma et al. 2004, this is a double array of (x, y).
     * @param rCoeffs radial distortion vector of length 2
     * distortion vector of length 5.  can be null to skip lens distortion correction.
     * if 
     * @param focalLength focal length of camera in units of pixels.
     * @param centerX x coordinate of principal point in pixels, usually image center.
     * @param centerY y coordinate of principal point in pixels, usually image center.
     * @return pixels in the reference frame of image with distortion applied.
     * In terms of Table 1 of Ma et al. 2004, this is a double array of (u_d, v_d)
     * @param useR2R4 if true use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
        f(r) = 1 +k1*r^2 + k2*r^4,
        else if false use model #3 f(r) = 1 +k1*r + k2*r^2.
        note that if rCoeffs is null or empty, no radial distortion is applied.
     */
    public static double[][] cameraToPixelCoordinates(double[][] xC, double[] rCoeffs,
        double focalLength, double centerX, double centerY, boolean useR2R4) {
        
        // http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/parameters.html
        
        double[][] cc = MatrixUtil.copy(xC);
        for (int i = 0; i < xC[0].length; ++i) {
            // normalized pinhole projection X_c/Z_c and 
            cc[0][i] /= xC[2][i];
            cc[1][i] /= xC[2][i];
        }
        
        if (rCoeffs != null && rCoeffs.length > 0) {
            // input and output cc are in camera reference frame
            cc = CameraCalibration.applyRadialDistortion(cc, rCoeffs[0], rCoeffs[1],
                useR2R4);
        }
         
//        focalLength = Math.abs(focalLength);
        
        double[][] cameraIntr = Camera.createIntrinsicCameraMatrix(focalLength, centerX, centerY);
                       
        cc = MatrixUtil.multiply(cameraIntr, cc);
        
        return cc;
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
     * @param rCoeffs radial distortion vector of length 2
     * distortion vector of length 5.  can be null to skip lens distortion correction.
     * @param kIntr
     * @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
        f(r) = 1 +k1*r^2 + k2*r^4 if true,
        else use model #3 f(r) = 1 +k1*r + k2*r^2.
        note that if rCoeffs is null or empty, no radial distortion is applied.
     */
    public static double[][] cameraToPixelCoordinates(double[][] xC, double[] rCoeffs,
        CameraIntrinsicParameters kIntr, boolean useR2R4) {
        
        // http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/parameters.html
        
        double[][] cc = MatrixUtil.copy(xC);
        for (int i = 0; i < xC[0].length; ++i) {
            // normalized pinhole projection X_c/Z_c and 
            cc[0][i] /= xC[2][i];
            cc[1][i] /= xC[2][i];
        }
        
        if (rCoeffs != null && rCoeffs.length > 0) {
            // input and output cc are in camera reference frame
            cc = CameraCalibration.applyRadialDistortion(cc, rCoeffs[0], rCoeffs[1],
                useR2R4);
        }
                                        
        cc = MatrixUtil.multiply(kIntr.getIntrinsic(), cc);
        
        return cc;
    }
    
     /** converts pixel coordinates to normalized camera coordinates by transforming them to camera 
    reference frame then applying Lp2-normalization.
     * @param x points in the camera centered reference frame. 
     * format is 3XN for N points in row-major format. 
     * @param intrinsic 
     * @return pixels transformed to camera coordinate reerence frame then 
     * Lp2-normalized.
     * @throws no.uib.cipr.matrix.NotConvergedException
     */
    public static double[][] pixelToNormalizedCameraCoordinates(double[][] x,
        CameraIntrinsicParameters intrinsic) throws NotConvergedException {
        
        double[][] kIntrInv = Camera.createIntrinsicCameraMatrixInverse(intrinsic.getIntrinsic());
        
        // the direction of the points is calculated by K^-1 * x
        double[][] xDirection = MatrixUtil.multiply(kIntrInv, x);
        double sum;
        int i;
        for (int col = 0; col < x[0].length; ++col) {
            sum = (x[0][col]*x[0][col]) + (x[1][col]*x[1][col]) + (x[2][col]*x[2][col]);
            sum = Math.sqrt(sum);
            for (i = 0; i < 3; ++i) {
                xDirection[i][col] /= sum;
            }
        }
                
        return xDirection;
    }
    
    /**
     * transform the world coordinates xW to the camera reference frame using
     * (rot * xW + trans).
     * @param xW coordinates of objects in a world reference frame.  
     * The format is 3XN for N points.
     * @param rot the rotation matrix to apply to the translated coordinates.
     * the size is 3X3.
     * @param trans the translation to apply to coordinates xW in x, y, and z.
     * The length of the array is 3.
     * @return xW transformed to the camera reference frame.
     */
    public static double[][] worldToCameraCoordinates(double[][] xW, double[][] rot,
        double[] trans) {
        if (xW.length != 3) {
            throw new IllegalArgumentException("xW.length must be 3 (for x, y, z)");
        }
        if (trans.length != 3) {
            throw new IllegalArgumentException("trans.length must be 3 (for x, y, z)");
        }
        if (rot.length != 3 || rot[0].length != 3) {
            throw new IllegalArgumentException("rot must be 3X3");
        }
        
        int n = xW[0].length;
        double[][] xC = MatrixUtil.copy(xW);
        xC = MatrixUtil.multiply(rot, xC);
        int i, j;
        for (i = 0; i < n; ++i) {
            for (j = 0; j < 3; ++j) {
                xC[j][i] += trans[j];
            }
        }
        
        return xC;
    }
    
    /**
     * transform the world coordinate point xWPt to the camera reference frame using
     * (rot * xW + trans).
     * @param xWPt coordinates of an object in a world reference frame.  
     * The length is 3.
     * @param rot the rotation matrix to apply to the translated coordinates.
     * the size is 3X3.
     * @param trans the translation to apply to coordinates xW in x, y, and z.
     * The length of the array is 3.
     * @return the point xW transformed to the camera reference frame.
     */
    public static double[] worldToCameraCoordinates(double[] xWPt, double[][] rot,
        double[] trans) {
        if (xWPt.length != 3) {
            throw new IllegalArgumentException("xW.length must be 3 (for x, y, z)");
        }
        if (trans.length != 3) {
            throw new IllegalArgumentException("trans.length must be 3 (for x, y, z)");
        }
        if (rot.length != 3 || rot[0].length != 3) {
            throw new IllegalArgumentException("rot must be 3X3");
        }
        
        double[] xC = Arrays.copyOf(xWPt, xWPt.length);
        xC = MatrixUtil.multiplyMatrixByColumnVector(rot, xC);
        int j;
        for (j = 0; j < 3; ++j) {
            xC[j] += trans[j];
        }
        
        return xC;
    }
    
    /**
     * transform the world coordinate point xWPt to the camera reference frame using
     * (rot * xW + trans).
     * @param xWPt coordinates of an object in a world reference frame.  
     * The length is 3.
     * @param rot the rotation matrix to apply to the translated coordinates.
     * the size is 3X3.
     * @param trans the translation to apply to coordinates xW in x, y, and z.
     * The length of the array is 3.
     * @param aux an auxiliary array of length 3 used internally to hold
     * values for a calculation.  it's present to help the invoker reuse an
     * object instead of creating another.  it's values are not used.
     * @param out the point xW transformed to the camera reference frame.
     */
    public static void worldToCameraCoordinates(double[] xWPt, double[][] rot,
        double[] trans, double[] aux, double[] out) {
        if (xWPt.length != 3) {
            throw new IllegalArgumentException("xW.length must be 3 (for x, y, z)");
        }
        if (trans.length != 3) {
            throw new IllegalArgumentException("trans.length must be 3 (for x, y, z)");
        }
        if (rot.length != 3 || rot[0].length != 3) {
            throw new IllegalArgumentException("rot must be 3X3");
        }
        if (out.length != 3) {
            throw new IllegalArgumentException("out.length must be 3 (for x, y, z)");
        }
                
        //point_camera = (1/x_z) * (R * point_world + translation) is called the normalized coordinate.
        
        double[] xC = Arrays.copyOf(xWPt, xWPt.length);
      
        MatrixUtil.multiplyMatrixByColumnVector(rot, xC, aux);
        
        int j;
        for (j = 0; j < 3; ++j) {
            out[j] = aux[j] + trans[j];
        }        
    }
    
    /**
     * transform the world coordinate point xWPt to the camera reference frame using
     * (rot * xW + trans).
     * @param xWPt coordinates of an object in a world reference frame.  
     * The length is 3.
     * @param rot the rotation matrix to apply to the translated coordinates.
     * the size is 3X3.
     * @param trans the translation to apply to coordinates xW in x, y, and z.
     * The length of the array is 3.
     * @param out the point transformed to the camera reference frame.
     */
    public static void worldToCameraCoordinates(double[] xWPt, double[][] rot,
        double[] trans, double[] out) {
        if (xWPt.length != 3) {
            throw new IllegalArgumentException("xW.length must be 3 (for x, y, z)");
        }
        if (trans.length != 3) {
            throw new IllegalArgumentException("trans.length must be 3 (for x, y, z)");
        }
        if (rot.length != 3 || rot[0].length != 3) {
            throw new IllegalArgumentException("rot must be 3X3");
        }
        double[] x2 = Arrays.copyOf(xWPt, xWPt.length);
        MatrixUtil.multiplyMatrixByColumnVector(rot, x2, out);
        int j;
        for (j = 0; j < 3; ++j) {
            out[j] = x2[j] + trans[j];
        }
        
    }
    
    /** converts pixel coordinates to camera coordinates by transforming them to camera 
    reference frame then removing radial distortion.
    The input in terms of Table 1 of Ma et al. 2004 is a double array of (u_d, v_d)
    and the output is a double array of (x, y).
    Also useful reading is NVM Tools by Alex Locher
    https://github.com/alexlocher/nvmtools.git
    
     * @param x points in the camera centered reference frame. 
     * format is 3XN for N points.  
     * @param kIntr  
     * @param rCoeffs radial distortion vector of length 2 or radial and tangential
     * distortion vector of length 5.  can be null to skip lens distortion correction.
     * @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
        f(r) = 1 +k1*r^2 + k2*r^4 if true,
        else use model #3 f(r) = 1 +k1*r + k2*r^2.
        note that if rCoeffs is null or empty, no radial distortion is removed.
     * @return pixels in the reference frame of 
     * @throws no.uib.cipr.matrix.NotConvergedException 
     * @throws java.io.IOException 
     */
    public static double[][] pixelToCameraCoordinates(double[][] x, 
        CameraIntrinsicParameters kIntr, double[] rCoeffs, boolean useR2R4) 
        throws NotConvergedException, IOException {
                        
        double[][] intr = MatrixUtil.copy(kIntr.getIntrinsic());
        
        return pixelToCameraCoordinates(x, intr, rCoeffs, useR2R4);
    }
    
    /** converts pixel coordinates to camera coordinates by transforming them to camera 
    reference frame then removing radial distortion.
    The input in terms of Table 1 of Ma et al. 2004 is a double array of (u_d, v_d)
    and the output is a double array of (x, y).
    Also useful reading is NVM Tools by Alex Locher
    https://github.com/alexlocher/nvmtools.git
    
     * @param x points in the camera centered reference frame. 
     * format is 3XN for N points.  
     * @param intr intrinsic camera parameters. dimensions aer [3 X 3]
     * @param rCoeffs radial distortion vector of length 2 or radial and tangential
     * distortion vector of length 5.  can be null to skip lens distortion correction.
     * @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
        f(r) = 1 +k1*r^2 + k2*r^4 if true,
        else use model #3 f(r) = 1 +k1*r + k2*r^2.
        note that if rCoeffs is null or empty, no radial distortion is removed.
     * @return pixels in the reference frame of 
     * @throws no.uib.cipr.matrix.NotConvergedException 
     * @throws java.io.IOException 
     */
    public static double[][] pixelToCameraCoordinates(double[][] x, 
        double[][] intr, double[] rCoeffs, boolean useR2R4) 
        throws NotConvergedException, IOException {
        
        // http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/parameters.html
                        
        // use absolute value of focal lengths
//        intr[0][0] = Math.abs(intr[0][0]);
//        intr[1][1] = Math.abs(intr[1][1]);
       
        double[][] cameraIntrInv = Camera.createIntrinsicCameraMatrixInverse(
            intr);
        
        // put x into camera coordinates reference frame:
        double[][] pix = MatrixUtil.multiply(cameraIntrInv, x);
        
        if (rCoeffs != null) {
            pix = CameraCalibration.removeRadialDistortion(pix, rCoeffs[0], 
                rCoeffs[1], useR2R4);
        }
                
        return pix;
    }
    
    /**
     * transform 2D pixel measurements x_i to unit norm 3D directions
     * <pre>
     * references Szeliski 2010, eqn (6.36)
     * </pre>
     * @param x 2D pixel measurements in format 3XN where the rows are x, y, and "1"'s
     * and the columns are each point in the N points.
     * @param k camera intrinsic parameters
     * @return 
     */
    public static double[][] transformToUnitNormDirections(double[][] x, 
        CameraIntrinsicParameters k) {
        
        if (x.length != 3) {
            throw new IllegalArgumentException("x length must be 3");
        }
        
        double[][] kInv = Camera.createIntrinsicCameraMatrixInverse(k.intrinsic);
        
        //(K^-1 * x_i)/(||K^-1 * x_i||)
        double[][] norm = MatrixUtil.multiply(kInv, x);
                
        int i, j;
        double sum = 0;
        for (j = 0; j < norm.length; ++j) {
            sum = 0;
            for (i = 0; i < norm.length; ++i) {
                sum += (norm[i][j]*norm[i][j]);
            }
            sum = Math.sqrt(sum);
            for (i = 0; i < norm.length; ++i) {
                norm[i][j] /= sum;
            }
        }
        return norm;
    }
    
    public static class CameraIntrinsicParameters {
        private double[][] intrinsic;
        private double lambda;
        
        public CameraIntrinsicParameters(double[][] k) {
            this.intrinsic = k;
        }
        
        public CameraIntrinsicParameters() {
        }
        
        /**
         * @return the intrinsic parameters
         */
        public double[][] getIntrinsic() {
            return intrinsic;
        }
        /**
         * @param intrinsic the intrinsic parameters to set
         */
        public void setIntrinsic(double[][] intrinsic) {
            this.intrinsic = intrinsic;
        }

        /**
         * @return the lambda the scale factor used in projection
         */
        public double getLambda() {
            return lambda;
        }

        /**
         * @param lambda the lambda to set for scale factor of projection
         */
        public void setLambda(double lambda) {
            this.lambda = lambda;
        }
        
    }
    
    public static class CameraProjection {
        /**
         * the projection matrix of a camera which is a 3X4 matrix of
         * intrinsic times extrinsic parameter matrices
         */
        private double[][] p;
        public CameraProjection(double[][] projection) {
            this.p = projection;
        }

        /**
         * @return the p
         */
        public double[][] getP() {
            return p;
        }

        /**
         * @param p the p to set
         */
        public void setP(double[][] p) {
            this.p = p;
        }
    }
    
    public static class CameraParameters {
        private final CameraIntrinsicParameters intrinsicParameters;
        private final CameraExtrinsicParameters extrinsicParameters;

        public CameraParameters(CameraIntrinsicParameters intrinsics,
                CameraExtrinsicParameters extrinsics) {
            this.intrinsicParameters = intrinsics;
            this.extrinsicParameters = extrinsics;
        }
        
        public double[][] createProjectionMatrix() {
            
            double[][] rt = new double[3][4];
            int i, j;
            for (i = 0; i < 3; ++i) {
                rt[i] = new double[4];
                System.arraycopy(extrinsicParameters.rotation[i], 0, rt[i], 0, 3);
                rt[i][3] = extrinsicParameters.translation[i];
            }
            
            double[][] p = MatrixUtil.multiply(intrinsicParameters.getIntrinsic(), rt);
            return p;
        }
        /**
         * @return the intrinsicParameters
         */
        public CameraIntrinsicParameters getIntrinsicParameters() {
            return intrinsicParameters;
        }

        /**
         * @return the extrinsicParameters
         */
        public CameraExtrinsicParameters getExtrinsicParameters() {
            return extrinsicParameters;
        }
    }
    
    public static class CameraExtrinsicParameters {
        private double[][] rotation;
        private double[] translation;

        /**
         * @return the rotation
         */
        public double[][] getRotation() {
            return rotation;
        }

        /**
         * @param rotation the rotation to set
         */
        public void setRotation(double[][] rotation) {
            this.rotation = rotation;
        }

        /**
         * @return the translation
         */
        public double[] getTranslation() {
            return translation;
        }

        /**
         * @param translation the translation to set
         */
        public void setTranslation(double[] translation) {
            this.translation = translation;
        }
        
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("rot=\n");
            if (rotation != null) {
                sb.append(FormatArray.toString(rotation, "%.4e"));
            }
            sb.append("trans=\n");
            if (translation != null) {
                sb.append(FormatArray.toString(translation, "%.4e"));
            }
            return sb.toString();
        }
    }
    
    public static class CameraMatrices {
        private CameraIntrinsicParameters intrinsics;
        private List<CameraExtrinsicParameters> extrinsics = new ArrayList<CameraExtrinsicParameters>();
        private double[] radialDistortion;
        private boolean useR2R4;
        
        /**
         * @return the radialDistortion
         */
        public double[] getRadialDistortCoeff() {
            return radialDistortion;
        }
        
        /**
         * @return true if using radial distortion function from Ma et al. 2004 
         for model #4 in Table 2, f(r) = 1 +k1*r^2 + k2*r^4,
         else return false if using model #3 f(r) = 1 +k1*r + k2*r^2.
        */
        public boolean useR2R4() {
            return useR2R4;
        }

        /**
         * @param radialDistortion the radialDistortion to set
         * @param useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
         f(r) = 1 +k1*r^2 + k2*r^4 if true,
         else use model #3 f(r) = 1 +k1*r + k2*r^2.
         */
        public void setRadialDistortion(double[] radialDistortion, boolean useR2R4) {
            this.radialDistortion = radialDistortion;
            this.useR2R4 = useR2R4;
        }

        /**
         * @return the intrinsics
         */
        public CameraIntrinsicParameters getIntrinsics() {
            return intrinsics;
        }

        /**
         * @param intrinsics the intrinsics to set
         */
        public void setIntrinsics(CameraIntrinsicParameters intrinsics) {
            this.intrinsics = intrinsics;
        }

        /**
         * @return the extrinsics
         */
        public List<CameraExtrinsicParameters> getExtrinsics() {
            return extrinsics;
        }

        /**
         * @param extrinsics the extrinsics to set
         */
        public void addExtrinsics(CameraExtrinsicParameters extrinsics) {
            this.extrinsics.add(extrinsics);
        }
        /**
         * @param extrinsics the extrinsics to set
         */
        public void addExtrinsics(List<CameraExtrinsicParameters> extrinsics) {
            this.extrinsics.addAll(extrinsics);
        }
    }
}
