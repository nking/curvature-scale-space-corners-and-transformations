package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.statistics.Standardization;
import algorithms.util.FormatArray;
import java.io.IOException;
import java.util.Arrays;
import no.uib.cipr.matrix.NotConvergedException;

/**
 should use many more points than the 8 used here 
 and preferably more than the 3 images used here too.
  
 
 * images are 3024 X 4032
          errors likely > 8 pixels

             img2         img1            img3
        #1 678, 718     608, 530        744, 806
        #2 2210, 886    2462, 512       2286, 526
        #3 1504, 1102   1484, 858       1410, 992
        #4 1516, 1814   1486, 1678      1432, 1760
        #5 1228, 2014   1154, 1882      1164, 1940
        #6 2042, 1936   2142, 1810      2018, 1866
        #7 698, 2944    614, 2802       788, 2728
        #8 2210, 2782   2366, 2848      2276, 2930

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
             xc=1512
             yc=2016
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
        // mean=-1.4625, 0.4625, 41.5
        // stdev=8.82, 11.303, 0.
        
       // attempting to line up center of book in WCS and in the first image
       //    in order to make interpretation of calculated translation easier.
       // image at pose1 has book center to the left by about 26 pixels 
       //  (= about 18  in wcs units [cm] where have estimated the scale factor
       //   from the book width in WCS and in pixels).
       // the offset in y is about 21 in WCS units
       /*for (int i = 0; i < wcs.length;++i) {
            wcs[i][0] +=8.5; // changed this from 9 to 8.5 by trial and error until the
                             // extracted intrinsic parameters for image center were reasonable
                             // and skew near 0.
        }
       */
       
        //double[] mn = new double[3];
        //double[] sd = new double[3];
        //double[][] wcs2 = Standardization.standardUnitNormalization(
        //    wcs, mn, sd);
        // mean=-1.4625, 0.4625, 41.5
        // stdev=8.82, 11.303, 0.
        
        for (int i = 0; i < wcs.length;++i) {
            wcs[i][0] +=8.5; // changed this from 9 to 8.5 by trial and error until the
                             // extracted intrinsic parameters for image center were reasonable
                             // and skew near 0.
        }
        
        wcs = MatrixUtil.transpose(wcs);        
        return wcs;
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
        
        //TODO: need to redo these on the full images that haven't been cropped.
        
        double[][] uv = new double[nFeatures][];
        switch(idx) {
            case 0:
                uv[0] = new double[]{608, 530, 1};
                uv[1] = new double[]{2462, 512, 1};
                uv[2] = new double[]{1484, 858, 1};
                uv[3] = new double[]{1486, 1678, 1};
                uv[4] = new double[]{1154, 1882, 1};
                uv[5] = new double[]{2142, 1810, 1};
                uv[6] = new double[]{614, 2802, 1};
                uv[7] = new double[]{2366, 2848, 1};
                break;
            case 1:
                uv[0] = new double[]{678, 718, 1};
                uv[1] = new double[]{2210, 886, 1};
                uv[2] = new double[]{1504, 1102, 1};
                uv[3] = new double[]{1516, 1814, 1};
                uv[4] = new double[]{1228, 2014, 1};
                uv[5] = new double[]{2042, 1936, 1};
                uv[6] = new double[]{698, 2944, 1};
                uv[7] = new double[]{2210, 2782, 1};
                break;
            default:
                uv[0] = new double[]{744, 806, 1};
                uv[1] = new double[]{2286, 526, 1};
                uv[2] = new double[]{1410, 992, 1};
                uv[3] = new double[]{1432, 1760, 1};
                uv[4] = new double[]{1164, 1940, 1};
                uv[5] = new double[]{2018, 1866, 1};
                uv[6] = new double[]{788, 2728, 1};
                uv[7] = new double[]{2276, 2930, 1};
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
            // cc along y is + for right-handed system
            case 1:
                thetas = new double[]{0, 23.5*Math.PI/180., 0};
                break;
            case 0:
                thetas = new double[]{0, 0, 0};        
                break;
            case 2:
                thetas = new double[]{0.*Math.PI/180., -23.5*Math.PI/180., 0};
                break;
            default:
                throw new IllegalArgumentException("idx out of range");
        }
        return Rotation.createRotationZYX(thetas[0], thetas[1], thetas[2]);
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
        // offset of image center of book and center of image
        //    gives the offset in Y of the WCS reference frame origin to the image center (==camera origin).
        // that's (2016-1678)/81=4.2
        switch (idx) {
            case 1:
                //-6.636e+00, -3.037e+00, 3.537e+01
                return new double[]{-18, -4.2, 41.5};
            case 0:
                //-7.260e+00, -4.721e+00, 3.410e+01
                return new double[]{0, -4.2, 41.5};
            default:
                //-7.376e+00, -4.029e+00, 4.009e+01 
                return new double[]{18, -4.2, 41.5};
        }        
    }
    
    public static double[][] getIntrinsicCameraMatrix() {
        double[][] intr = MatrixUtil.zeros(3, 3);
        intr[0][0] = 1600; //2.189e+03
        intr[1][1] = 1600; //2.886e+03
        intr[0][2] = 1512; //1500
        intr[1][2] = 2016; //2000
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
    
    public static void printExpectedTriangulation() throws NotConvergedException {
        
        double[][] intr = new double[3][];
        intr[0] = new double[]{1600, 0, 1512};
        intr[1] = new double[]{0, 1600, 2016};
        intr[2] = new double[]{0, 0, 1};
        
        double[][] r0 = Rotation.createRotationZYX(0, 0, 0);
        
        double[][] r1 = Rotation.createRotationZYX(
            0, 23.5*Math.PI/180., 0);
        
        double[] t0 = new double[]{0, 0, 0};
        double[] t1 = new double[]{-18, 0, 0};
        
        double[][] xi0 = getObservedFeaturesInImage(0);
        double[][] xi1 = getObservedFeaturesInImage(1);
        
        double[][] x0 = MatrixUtil.zeros(3, 1);
        double[][] x1 = MatrixUtil.zeros(3, 1);
        
        int i, k;
        for (i = 0; i < xi0[0].length; ++i) {
            for (k = 0; k < 3; ++k) {
                x0[k][0] = xi0[k][i];
                x1[k][0] = xi1[k][i];
            }

            Triangulation.WCSPt wcsPt = Triangulation.calculateWCSPoint(
            intr, r0, t0,
            intr, r1, t1,
            x0, x1);

            double[] xw = wcsPt.X;
            for (k = 0; k < 4; ++k) {
                xw[k] /= xw[3];
            }
            
            System.out.printf("xw[%d]=%s\n\n", i, FormatArray.toString(xw, "%.3e"));
        }
    }
    
     /**
      * cannot use this method because translation[2] is 0 for all poses.
     * for each image, project world scene features into the camera reference
     * frame.  no distortion is added.
     * It uses  scale * projected = H * coordsW with H = intrinsic * | r1 r2 t|.
     * @return
     
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
 //cannot use this when translation[2] == 0
                //transform to camera reference frame. size [1X3]
                MatrixUtil.multiplyMatrixByColumnVector(h, coordsWI, coordsWIC);
                for (k = 0; k < 3; ++k) {
                    coordsWICN[k] = coordsWIC[k] / coordsWIC[2];
                    out[k][j*nFeatures + i] = coordsWICN[k];
                }
            }
        }
        return out;
    }*/
    
    public static double[][] getObservedMinusProjected_Camera_Frame() throws NotConvergedException, IOException {
        
        double[][] x = getObservedTransformedToCameraFrame();
        double[][] xHat = getFeaturesProjectedToAllCameraFrames();
        
        double[][] m = MatrixUtil.pointwiseSubtract(x, xHat);
        
        return m;
    }
    
    /* cannot use this method because translation[2] is 0 for all poses.
    public static double[][] getObservedMinusProjected_Camera_Frame_H_LftHnd() throws NotConvergedException, IOException {
        
        double[][] x = getObservedTransformedToCameraFrame();
        double[][] xHat = getFeaturesProjectedToAllCameraFrames_H_LftHnd();
        
        double[][] m = MatrixUtil.pointwiseSubtract(x, xHat);
        
        return m;
    }
    */
    
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
    * cannot use this method because translation[2] is 0 for all poses.
     * for each image, project world scene features into the image reference
     * frame.  distortion is added internally after transformation to camera coordinates.
     * @return
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
    */
    
    public static double[][] getObservedMinusProjected_Image_Frame() {
        
        double[][] x = getObservedFeaturesInAllImages();
        double[][] xHat = getFeaturesProjectedToAllImageFrames();
        
        double[][] m = MatrixUtil.pointwiseSubtract(x, xHat);
        
        return m;
    }
    
    /**
      * cannot use this method because translation[2] is 0 for all poses.
    public static double[][] getObservedMinusProjected_Image_Frame_H_LftHnd() {
        
        double[][] x = getObservedFeaturesInAllImages();
        double[][] xHat = getFeaturesProjectedToAllImageFrames_H_LftHnd();
        
        double[][] m = MatrixUtil.pointwiseSubtract(x, xHat);
        
        return m;
    }
    */
    
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
        
        System.out.printf("sqsum=%.7e\n\n", sqsum);
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
        
        System.out.printf("sqsum=%.7e\n\n", sqsum);
    }
    
    /*
      * cannot use this method because translation[2] is 0 for all poses.
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
        
        System.out.printf("sqsum=%.7e\n\n", sqsum);
    }
    */
    
    /*
      * cannot use this method because translation[2] is 0 for all poses.
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
        
        System.out.printf("sqsum=%.7e\n\n", sqsum);
    }
*/
}
