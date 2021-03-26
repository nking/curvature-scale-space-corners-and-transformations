package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import java.util.Arrays;

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
}
