package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import java.io.IOException;
import no.uib.cipr.matrix.NotConvergedException;

/**
 * images are 3042 X 3504 (previously 3042x4032)
          errors likely > 8 pixels

                img2       img1         img3
          #1 560, 532     680, 720     744, 808
          #2 2428, 512    2208, 896    2284, 528
          #3 1484, 856    1508, 1100   1416, 1000
          #4 1488, 1708   1520, 1856   1428, 1784
          #5 1228, 1896   1296, 2012   1220, 1952
          #6 2136, 1908   2028, 2000   2012, 1972
          #7 620, 2800    704, 2940    788, 2718
          #8 2364, 2844   2208, 2780   2272, 2928

          features in WCS
          #1 -11, 14, 41.5
          #2  11, 14, 41.5
          #3   0, 9.7, 41.5
          #4   0, 0, 41.5
          #5 -3.7, -3, 41.5
          #6 -8, -3, 41.5
          #7 -11, -14, 41.5
          #8  11, -14, 41.5
        
        expecting
             focalLength ~ 1604 pixels = 2.245 mm
             no skew
             xc=1521
             yc=1752
             little to no radial distortion (if was present, it is already removed)
             rotation between images = 23.4 degrees
             translation between images = 18 cm

          other information:
            pixel width = 1.4e-3mm
            FOV = 77 degrees = 1.344 radians
            
 * @author nichole
 */
public class NeurochemistryBookData {
 
    public static final int nFeatures = 8;
    public static final int mImages = 3;
    
    /**
     * get the features in world scene coordinates.
     * @return double array of size 3 X 8
     */
    public static double[][] getFeatureWCS()  {
 
        double[][] wcs = new double[nFeatures][];
        wcs[0] = new double[]{-11, 14, 41.5};
        wcs[1] = new double[]{11, 14, 41.5};
        wcs[2] = new double[]{0, 9.7, 41.5};
        wcs[3] = new double[]{0, 0, 41.5};
        wcs[4] = new double[]{-3.7, -3, 41.5};
        wcs[5] = new double[]{-8, -3, 41.5};
        wcs[6] = new double[]{-11, -14, 41.5};
        wcs[7] = new double[]{11, -14, 41.5};
        
        return MatrixUtil.transpose(wcs);
    }
     
    /**
     * get the 8 for the 3 images in pixel coordinates.  the format
     * is all features of one image followed by all features of the next image,
     * etc.
     * @return a double array of size 3 X 8*3
     */
    public static double[][] getObservedFeaturesInAllImages() {
        
        double[][] uv = MatrixUtil.zeros(3, nFeatures*mImages);
        double[][] uvi;
        int i, k;
        for (i = 0; i < mImages; ++i) {
            //[3 X (nFeatures*mImages)]
            uvi = getObservedFeaturesInImage(i);
            for (k = 0; k < 3; ++k) {
                System.arraycopy(uvi[k], 0, uv[k], i*nFeatures, nFeatures);
            }
        }
        return uv;
    }
    
    // size is 3 X (nFeatures*mImages)
    public static double[][] getObservedFeaturesInImage(int idx) {
        if (idx < 0 || idx > mImages-1) {
            throw new IllegalArgumentException("idx must be 0 through 3, inclusive");
        }
        
        double[][] uv = new double[nFeatures][];
        switch(idx) {
            case 0:
                uv[0] = new double[]{560, 532, 1};
                uv[1] = new double[]{2428, 512, 1};
                uv[2] = new double[]{1484, 856, 1};
                uv[3] = new double[]{1488, 1708, 1};
                uv[4] = new double[]{1228, 1896, 1};
                uv[5] = new double[]{2136, 1908, 1};
                uv[6] = new double[]{620, 2800, 1};
                uv[7] = new double[]{2364, 2844, 1};
                break;
            case 1:
                uv[0] = new double[]{680, 720, 1};
                uv[1] = new double[]{2208, 896, 1};
                uv[2] = new double[]{1508, 1100, 1};
                uv[3] = new double[]{1520, 1856, 1};
                uv[4] = new double[]{1296, 2012, 1};
                uv[5] = new double[]{2028, 2000, 1};
                uv[6] = new double[]{704, 2940, 1};
                uv[7] = new double[]{2208, 2780, 1};
                break;
            default:
                uv[0] = new double[]{744, 808, 1};
                uv[1] = new double[]{2284, 528, 1};
                uv[2] = new double[]{1416, 1000, 1};
                uv[3] = new double[]{1428, 1784, 1};
                uv[4] = new double[]{1220, 1952, 1};
                uv[5] = new double[]{2012, 1972, 1};
                uv[6] = new double[]{788, 2718, 1};
                uv[7] = new double[]{2272, 2928, 1};
                break;
        }        
        
        return MatrixUtil.transpose(uv);
    }
    
    public static double[][] getTransposedRotation(int imageIdx) {
        double[][] r = getRotation(imageIdx);
        return MatrixUtil.transpose(r);
    }
    
    public static double[][] getRotation(int idx) {
        if (idx < 0 || idx > mImages-1) {
            throw new IllegalArgumentException("idx must be 0 through 3, inclusive");
        }
        double[] thetas;
        switch (idx) {
            // cc along z is + for right-handed system
            case 0:
                thetas = new double[]{0, 0, -23.5*Math.PI/180.};
                break;
            case 1:
                thetas = new double[]{0, 0, 0};        
                break;
            default:
                thetas = new double[]{0, 0, 23.5*Math.PI/180.};
                break;
        }
        return Rotation.createRotationZYX(thetas);
    }
    
    /**
     * translation is in WCS units of cm.
     * @param idx
     * @return 
     */
    public static double[] getTranslation(int idx) {
        if (idx < 0 || idx > mImages-1) {
            throw new IllegalArgumentException("idx must be 0 through 3, inclusive");
        }
        switch (idx) {
            case 0:
                return new double[]{-18,0, 0};
            case 1:
                return new double[]{0, 0, 0};
            default:
                return new double[]{18, 0, 0};
        }        
    }
    
    public static double[][] getIntrinsicCameraMatrix() {
        double[][] intr = MatrixUtil.zeros(3, 3);
        intr[0][0] = 1600;
        intr[1][1] = 1600;
        intr[0][2] = 1521;
        intr[1][2] = 1752;
        intr[2][2] = 1;
        return intr;
    }
    
    /**
     * the radial distortion coefficients are for the distortion polynomial k1*r^2 + k2*r^4.
     * @return 
     */
    public static double[] getRadialDistortionR2R4() {
        return new double[]{0, 0};
    }
    
    /**
     * transform the observed to camera coordinates.
       radial distortion is removed.
     * @return
     * @throws NotConvergedException 
     */
    public static double[][] getObservedTransformedToCameraFrame() throws NotConvergedException, IOException {
        double[][] coordsI = getObservedFeaturesInAllImages();
        
        double[] kRadial = getRadialDistortionR2R4();
        
        double[][] intr = getIntrinsicCameraMatrix();
        
        boolean useR2R4 = true;
        
        double[][] coordsIC = MatrixUtil.zeros(coordsI.length, coordsI[0].length);
        double[][] x = MatrixUtil.zeros(3, nFeatures);
        double[][] xc;
        
        // transform the observed to camera coordinates:
        //   (radial distortion is removed)
        int i, j, k;
        for (j = 0; j < mImages; ++j) {
            MatrixUtil.copySubMatrix(coordsI, 0, 2, j*nFeatures, ((j+1)*nFeatures)-1, x);            
            xc = Camera.pixelToCameraCoordinates(x, intr, kRadial, useR2R4);
            for (i = 0; i < 3; ++i) {
                for (k = 0; k < nFeatures; ++k) {
                    xc[i][k] /= xc[2][k];
                }
                System.arraycopy(xc[i], 0, coordsIC[i], j*nFeatures, nFeatures);
            }
        } 
        return coordsIC;
    }
    
    /**
     * for each image, project world scene features into the camera reference
     * frame.  no distortion is added.
     * @return
     */
    public static double[][] getFeaturesProjectedToAllCameraFrames() {
        double[][] coordsW = getFeatureWCS();
             
        int i, j, k;
        
        double[] coordsWI = new double[3];
        double[] coordsWIC = new double[3];
        double[] coordsWICN = new double[3];
        double[][] r;
        double[] t;
        double[] rAux = new double[3];
        
        double[][] out = MatrixUtil.zeros(3, nFeatures*mImages);
                
        // project the world scene features to each camera reference frame
        //   (not adding distortion)
        for (i = 0; i < nFeatures; ++i) {
            //populate xWI; extract the world feature.  size [1X3]
            MatrixUtil.extractColumn(coordsW, i, coordsWI);
            for (j = 0; j < mImages; ++j) {
                r = getRotation(j);
                t = getTranslation(j);
                //transform to camera reference frame. size [1X3]
                Camera.worldToCameraCoordinates(coordsWI, r, t, rAux, coordsWIC);
                for (k = 0; k < 3; ++k) {
                    coordsWICN[k] = coordsWIC[k] / coordsWIC[2];
                    out[k][j*nFeatures + i] = coordsWICN[k];
                }
            }
        }
        return out;
    }
    
     /**
     * for each image, project world scene features into the camera reference
     * frame.  no distortion is added.
     * It uses  scale * projected = H * coordsW with H = intrinsic * | r1 r2 t|.
     * @return
     */
    public static double[][] getFeaturesProjectedToAllCameraFrames_H_LftHnd() {
        double[][] coordsW = getFeatureWCS();
             
        int i, j, k;
        
        double[][] h = MatrixUtil.zeros(3, 3);
        
        double[] coordsWI = new double[3];
        double[] coordsWIC = new double[3];
        double[] coordsWICN = new double[3];
        double[][] r;
        double[] t;
        
        double[][] out = MatrixUtil.zeros(3, nFeatures*mImages);
                
        // project the world scene features to each camera reference frame
        //   (not adding distortion)
        for (i = 0; i < nFeatures; ++i) {
            //populate xWI; extract the world feature.  size [1X3]
            MatrixUtil.extractColumn(coordsW, i, coordsWI);
            
            for (j = 0; j < mImages; ++j) {
                r = getRotation(j);
                t = getTranslation(j);
                
                h[0][0] = r[0][0];
                h[1][0] = r[1][0];
                h[2][0] = r[2][0];
                h[0][1] = r[0][1];
                h[1][1] = r[1][1];
                h[2][1] = r[2][1];
                h[0][2] = t[0];
                h[1][2] = t[1];
                h[2][2] = t[2];
                
                //transform to camera reference frame. size [1X3]
                MatrixUtil.multiplyMatrixByColumnVector(h, coordsWI, coordsWIC);
                for (k = 0; k < 3; ++k) {
                    coordsWICN[k] = coordsWIC[k] / coordsWIC[2];
                    out[k][j*nFeatures + i] = coordsWICN[k];
                }
            }
        }
        return out;
    }
    
    public static double[][] getObservedMinusProjected_Camera_Frame() throws NotConvergedException, IOException {
        
        double[][] x = getObservedTransformedToCameraFrame();
        double[][] xHat = getFeaturesProjectedToAllCameraFrames();
        
        double[][] m = MatrixUtil.elementwiseSubtract(x, xHat);
        
        return m;
    }
    
    public static double[][] getObservedMinusProjected_Camera_Frame_H_LftHnd() throws NotConvergedException, IOException {
        
        double[][] x = getObservedTransformedToCameraFrame();
        double[][] xHat = getFeaturesProjectedToAllCameraFrames_H_LftHnd();
        
        double[][] m = MatrixUtil.elementwiseSubtract(x, xHat);
        
        return m;
    }
    
    /**
     * for each image, project world scene features into the image reference
     * frame.  distortion is added internally after transformation to camera coordinates.
     * @return
     */
    public static double[][] getFeaturesProjectedToAllImageFrames() {
        double[][] coordsWCN = getFeaturesProjectedToAllCameraFrames();
        
        boolean useR2R4 = true;
        double[] rd = getRadialDistortionR2R4();
        
        double[][] intr = getIntrinsicCameraMatrix();
        
        double[] xWCNI = new double[3];
        double[] xWCNDI = new double[3];
        double[] xWDI = new double[3];
        
        int i, j, k;
                
        double[][] out = MatrixUtil.zeros(3, nFeatures*mImages);
                
        for (i = 0; i < nFeatures; ++i) {
            for (j = 0; j < mImages; ++j) {
                MatrixUtil.extractColumn(coordsWCN, j*nFeatures + i, xWCNI);                
                CameraCalibration.applyRadialDistortion(xWCNI, rd[0], rd[1], useR2R4, xWCNDI);
                
                MatrixUtil.multiplyMatrixByColumnVector(intr, xWCNDI, xWDI);
                
                for (k = 0; k < 3; ++k) {
                    out[k][j*nFeatures + i] = xWDI[k];
                }
            }
        }        
        return out;
    }
    
    /**
     * for each image, project world scene features into the image reference
     * frame.  distortion is added internally after transformation to camera coordinates.
     * @return
     */
    public static double[][] getFeaturesProjectedToAllImageFrames_H_LftHnd() {
        
        double[][] coordsWCN = getFeaturesProjectedToAllCameraFrames_H_LftHnd();
        
        boolean useR2R4 = true;
        double[] rd = getRadialDistortionR2R4();
        
        double[] xWCNI = new double[3];
        double[] xWCNDI = new double[3];
        double[] xWDI = new double[3];
        
        double[][] intr = getIntrinsicCameraMatrix();
        
        int i, j, k;
                
        double[][] out = MatrixUtil.zeros(3, nFeatures*mImages);
                
        for (i = 0; i < nFeatures; ++i) {
            for (j = 0; j < mImages; ++j) {
                MatrixUtil.extractColumn(coordsWCN, j*nFeatures + i, xWCNI);                
                CameraCalibration.applyRadialDistortion(xWCNI, rd[0], rd[1], useR2R4, xWCNDI);
                MatrixUtil.multiplyMatrixByColumnVector(intr, xWCNDI, xWDI);
                for (k = 0; k < 3; ++k) {
                    out[k][j*nFeatures + i] = xWDI[k];
                }
            }
        }
        
        // or use zhang eqn 2:
        //   H = A * |r1 r2 t|
        
        return out;
    }
    
    public static double[][] getObservedMinusProjected_Image_Frame() {
        
        double[][] x = getObservedFeaturesInAllImages();
        double[][] xHat = getFeaturesProjectedToAllImageFrames();
        
        double[][] m = MatrixUtil.elementwiseSubtract(x, xHat);
        
        return m;
    }
    
    public static double[][] getObservedMinusProjected_Image_Frame_H_LftHnd() {
        
        double[][] x = getObservedFeaturesInAllImages();
        double[][] xHat = getFeaturesProjectedToAllImageFrames_H_LftHnd();
        
        double[][] m = MatrixUtil.elementwiseSubtract(x, xHat);
        
        return m;
    }
    
    
    public static void printObservedMinusProjected_Camera_Frame() throws NotConvergedException, IOException {
        
        double[][] m = getObservedMinusProjected_Camera_Frame();
        
        System.out.printf("obs-projected in camera=\n%s\n", FormatArray.toString(m, "%.7e"));
        
        double sqsum = 0;
        
        int i, j;
        for (i = 0; i < m.length; ++i) {
            for (j = 0; j < m[i].length; ++j) {
                sqsum += (m[i][j]*m[i][j]);
            }
        }
        
        System.out.printf("sqsum=%.7e\n", sqsum);
    }
    
    public static void printObservedMinusProjected_Image_Frame() {
       
        double[][] m = getObservedMinusProjected_Image_Frame();
        
        System.out.printf("obs-projected in images=\n%s\n", FormatArray.toString(m, "%.7e"));
        
        double sqsum = 0;
        
        int i, j;
        for (i = 0; i < m.length; ++i) {
            for (j = 0; j < m[i].length; ++j) {
                sqsum += (m[i][j]*m[i][j]);
            }
        }
        
        System.out.printf("sqsum=%.7e\n", sqsum);
    }
    
    public static void printObservedMinusProjected_Camera_Frame_H_LftHnd() throws NotConvergedException, IOException {
        
        double[][] m = getObservedMinusProjected_Camera_Frame_H_LftHnd();
        
        System.out.printf("obs-projected_H_LftHnd in camera=\n%s\n", FormatArray.toString(m, "%.7e"));
        
        double sqsum = 0;
        
        int i, j;
        for (i = 0; i < m.length; ++i) {
            for (j = 0; j < m[i].length; ++j) {
                sqsum += (m[i][j]*m[i][j]);
            }
        }
        
        System.out.printf("sqsum=%.7e\n", sqsum);
    }
    
    public static void printObservedMinusProjected_Image_Frame_H_LftHnd() {
        
        double[][] m = getObservedMinusProjected_Image_Frame_H_LftHnd();
        
        System.out.printf("obs-projected_eqn2 in images=\n%s\n", FormatArray.toString(m, "%.7e"));
        
        double sqsum = 0;
        
        int i, j;
        for (i = 0; i < m.length; ++i) {
            for (j = 0; j < m[i].length; ++j) {
                sqsum += (m[i][j]*m[i][j]);
            }
        }
        
        System.out.printf("sqsum=%.7e\n", sqsum);
    }
}
