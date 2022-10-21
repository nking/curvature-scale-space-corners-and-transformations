package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import no.uib.cipr.matrix.NotConvergedException;

/**
 * convenience methods for camera test data Zhang 1998.
 * see testresources/zhang1998/README.txt
 * 
        each image is 640×480  [pixels^2]
        The pixel is square (aspect ratio = 1); 
        the focal length = 832.5 pixels; 
        the image center is at (303.959, 206.585); 
        there is a significant radial distortion: k1 = -0.228601, k2 = 0.190353.
 * @author nichole
 */
public class Zhang98Data {
    
    public static String ZHANGDIR = "zhang1998";
    protected final static String sep = System.getProperty("file.separator");

    /**
     * get the 256 corners of the checkerboard in world reference frame.
     * the coordinates are in inches.
     * @return double array of size 3 X 256
     * @throws IOException 
     */
    public static double[][] getFeatureWCS() throws IOException {
        
        String path = ResourceFinder.findTestResourcesDirectory() + sep + ZHANGDIR
            + sep + "Model.txt";
        
        File f = new File(path);
        if (!f.exists()) {
            throw new IOException("could not find file at " + path);
        }
        
        BufferedReader in = new BufferedReader(new FileReader(new File(path)));
        
        double[][] wcs = readFile(in);
        
        in.close();
        
        return wcs;
    }
     
    /**
     * get the 256 checkerboard corners for the 5 images in pixel coordinates.
     * @return a double array of size 3 X 256*5
     */
    public static double[][] getObservedFeaturesInAllImages() throws IOException {
        
        int mImages = 5;
        int nFeatures = 256;
        double[][] uv = MatrixUtil.zeros(3, nFeatures*mImages);
        double[][] uvi;
        int i, j;
        for (i = 0; i < mImages; ++i) {
            uvi = getObservedFeaturesInImage(i+1);
            assert(nFeatures == uvi[0].length);
            for (j = 0; j < 3; ++j) {
                System.arraycopy(uvi[j], 0, uv[j], nFeatures*i, nFeatures);
            }
        }
        return uv;
    }
    
    // data1.txt through data5.txt
    // size is 3 X 256
    public static double[][] getObservedFeaturesInImage(int idx) throws IOException {
        if (idx < 1 || idx > 5) {
            throw new IllegalArgumentException("idx must be 1 through 5, inclusive");
        }
        
        String path = ResourceFinder.findTestResourcesDirectory() + sep + ZHANGDIR
            + sep + "data" + String.valueOf(idx) + ".txt";
        
        File f = new File(path);
        if (!f.exists()) {
            throw new IOException("could not find file at " + path);
        }
        
        BufferedReader in = new BufferedReader(new FileReader(new File(path)));
        
        double[][] xy = readFile(in);
        
        in.close();
        
        return xy;
    }
    
    private static double[][] readFile(BufferedReader in) throws IOException {
        
        int nFeatures = 256;
        
        // 3 dimensions, 4*64 features 
        double[][] a = MatrixUtil.zeros(3, nFeatures);
        Arrays.fill(a[2], 1.);
        
        // reading 64 lines. each line is 4 corners of 1 square (= 8 real numbers).
        int c = 0;
        int j;
        String[] coords;
        String line = in.readLine();
        while (line != null) {
            coords = line.split("\\s+");
            assert(coords.length == 8);
            
            for (j = 0; j < 4; ++j) {
                a[0][c] = Double.parseDouble(coords[2*j]);
                a[1][c] = Double.parseDouble(coords[2*j + 1]);
                c++;
            }
            
            line = in.readLine();
        }
        in.close();
        
        return a;
    }
    
    /*
    https://www.microsoft.com/en-us/research/uploads/prod/2016/12/completecalibration.txt
    
        832.5 0.204494 832.53 303.959 206.585

        -0.228601 0.190353

        0.992759 -0.026319 0.117201
        0.0139247 0.994339 0.105341
        -0.11931 -0.102947 0.987505
        -3.84019 3.65164 12.791

        0.997397 -0.00482564 0.0719419
        0.0175608 0.983971 -0.17746
        -0.0699324 0.178262 0.981495
        -3.71693 3.76928 13.1974

        0.915213 -0.0356648 0.401389
        -0.00807547 0.994252 0.106756
        -0.402889 -0.100946 0.909665
        -2.94409 3.77653 14.2456

        0.986617 -0.0175461 -0.16211
        0.0337573 0.994634 0.0977953
        0.159524 -0.101959 0.981915
        -3.40697 3.6362 12.4551

        0.967585 -0.196899 -0.158144
        0.191542 0.980281 -0.0485827
        0.164592 0.0167167 0.98622
        -4.07238 3.21033 14.3441
    */
    
    public static double[][] getTransposedRotation(int imageIdx) {
        double[][] r = getRotation(imageIdx);
        return MatrixUtil.transpose(r);
    }
    
    public static double[][] getRotation(int idx) {
        if (idx < 1 || idx > 5) {
            throw new IllegalArgumentException("idx must be 1 through 5, inclusive");
        }
        double[][] r = new double[3][];
        switch (idx) {
            case 1:
                r[0] = new double[]{0.992759, -0.026319, 0.117201};
                r[1] = new double[]{0.0139247, 0.994339, 0.105341};
                r[2] = new double[]{-0.11931, -0.102947, 0.987505};
                break;
            case 2:
                r[0] = new double[]{0.997397, -0.00482564, 0.0719419};
                r[1] = new double[]{0.0175608, 0.983971, -0.17746};
                r[2] = new double[]{-0.0699324, 0.178262, 0.981495};        
                break;
            case 3:
                r[0] = new double[]{0.915213, -0.0356648, 0.401389};
                r[1] = new double[]{-0.00807547, 0.994252, 0.106756};
                r[2] = new double[]{-0.402889, -0.100946, 0.909665};        
                break;
            case 4:
                r[0] = new double[]{0.986617, -0.0175461, -0.16211};
                r[1] = new double[]{0.0337573, 0.994634, 0.0977953};
                r[2] = new double[]{0.159524, -0.101959, 0.981915};                
                break; 
            default:
                r[0] = new double[]{0.967585, -0.196899, -0.158144};
                r[1] = new double[]{0.191542, 0.980281, -0.0485827};
                r[2] = new double[]{0.164592, 0.0167167, 0.98622};
                break;
        }
        return r;
    }
    
    public static double[] getTranslation(int idx) {
        if (idx < 1 || idx > 5) {
            throw new IllegalArgumentException("idx must be 1 through 5, inclusive");
        }
        double[] t = null;
        switch (idx) {
            case 1:
                t = new double[]{-3.84019, 3.65164, 12.791};
                break;
            case 2:
                t = new double[]{-3.71693, 3.76928, 13.1974};
                break;
            case 3:
                t = new double[]{-2.94409, 3.77653, 14.2456};
                break;
            case 4:
                t = new double[]{-3.40697, 3.6362, 12.4551};
                break;
            default:
                t = new double[]{-4.07238, 3.21033, 14.3441};
                break;
        }        
        return t;
    }
    
    public static double[][] getIntrinsicCameraMatrix() {
        double[][] intr = MatrixUtil.zeros(3, 3);
        intr[0][0] = 832.5;
        intr[1][1] = 832.5;
        intr[0][2] = 303.959;  // dimension is 640 pixels
        intr[1][2] = 206.585;  // dimension is 480 pixels
        intr[2][2] = 1;
        return intr;
    }
    
    /**
     * the radial distortion coefficients are for the distortion polynomial k1*r^2 + k2*r^4.
     * @return 
     */
    public static double[] getRadialDistortionR2R4() {
        return new double[]{-0.228601, 0.190353};
    }
    
    /**
     * transform the observed to camera coordinates.
       radial distortion is removed.
     * @return
     * @throws IOException
     * @throws NotConvergedException 
     */
    public static double[][] getObservedTransformedToCameraFrame() throws IOException, NotConvergedException {
        double[][] coordsI = getObservedFeaturesInAllImages();
        double[][] coordsW = getFeatureWCS();
        int nFeatures = coordsW[0].length;
        int mImages = coordsI[0].length/nFeatures;
        
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
     * @throws IOException 
     */
    public static double[][] getFeaturesProjectedToAllCameraFrames() throws IOException {
        double[][] coordsI = getObservedFeaturesInAllImages();
        double[][] coordsW = getFeatureWCS();
        int nFeatures = coordsW[0].length;
        int mImages = coordsI[0].length/nFeatures;
             
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
                r = getRotation(j+1);
                t = getTranslation(j+1);
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
     * It uses Zhang equation 2 for scale * projected = H * coordsW with H = intrinsic * | r1 r2 t|.
     * @return
     * @throws IOException 
     */
    public static double[][] getFeaturesProjectedToAllCameraFramesEqn2() throws IOException {
        double[][] coordsI = getObservedFeaturesInAllImages();
        double[][] coordsW = getFeatureWCS();
        int nFeatures = coordsW[0].length;
        int mImages = coordsI[0].length/nFeatures;
             
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
                r = getRotation(j+1);
                t = getTranslation(j+1);
                
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
    
    public static double[][] getObservedMinusProjected_Camera_Frame() throws IOException, NotConvergedException {
        
        double[][] x = getObservedTransformedToCameraFrame();
        double[][] xHat = getFeaturesProjectedToAllCameraFrames();
        
        double[][] m = MatrixUtil.pointwiseSubtract(x, xHat);
        
        return m;
    }
    
    public static double[][] getObservedMinusProjected_Camera_Frame_Eqn2() throws IOException, NotConvergedException {
        
        double[][] x = getObservedTransformedToCameraFrame();
        double[][] xHat = getFeaturesProjectedToAllCameraFramesEqn2();
        
        double[][] m = MatrixUtil.pointwiseSubtract(x, xHat);
        
        return m;
    }
    
    /**
     * for each image, project world scene features into the image reference
     * frame.  distortion is added internally after transformation to camera coordinates.
     * @return
     * @throws IOException 
     */
    public static double[][] getFeaturesProjectedToAllImageFrames() throws IOException {
        double[][] coordsWCN = getFeaturesProjectedToAllCameraFrames();
        double[][] coordsW = getFeatureWCS();
        int nFeatures = coordsW[0].length;
        int mImages = coordsWCN[0].length/nFeatures;
        
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
     * @throws IOException 
     */
    public static double[][] getFeaturesProjectedToAllImageFramesEqn2() throws IOException {
        
        double[][] coordsWCN = getFeaturesProjectedToAllCameraFramesEqn2();
        double[][] coordsW = getFeatureWCS();
        
        int nFeatures = coordsW[0].length;
        int mImages = coordsWCN[0].length/nFeatures;
        
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
    
    public static double[][] getObservedMinusProjected_Image_Frame() throws IOException, NotConvergedException {
        
        double[][] x = getObservedFeaturesInAllImages();
        double[][] xHat = getFeaturesProjectedToAllImageFrames();
        
        double[][] m = MatrixUtil.pointwiseSubtract(x, xHat);
        
        return m;
    }
    
    public static double[][] getObservedMinusProjected_Image_Frame_Eqn2() throws IOException, NotConvergedException {
        
        double[][] x = getObservedFeaturesInAllImages();
        double[][] xHat = getFeaturesProjectedToAllImageFramesEqn2();
        
        double[][] m = MatrixUtil.pointwiseSubtract(x, xHat);
        
        return m;
    }
    
    
    public static void printObservedMinusProjected_Camera_Frame() throws IOException, NotConvergedException {
        
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
    
    public static void printObservedMinusProjected_Image_Frame() throws IOException, NotConvergedException {
       
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
    
    public static void printObservedMinusProjected_Camera_Frame_Eqn2() throws IOException, NotConvergedException {
        
        double[][] m = getObservedMinusProjected_Camera_Frame_Eqn2();
        
        System.out.printf("obs-projected_eqn2 in camera=\n%s\n", FormatArray.toString(m, "%.7e"));
        
        double sqsum = 0;
        
        int i, j;
        for (i = 0; i < m.length; ++i) {
            for (j = 0; j < m[i].length; ++j) {
                sqsum += (m[i][j]*m[i][j]);
            }
        }
        
        System.out.printf("sqsum=%.7e\n", sqsum);
    }
    
    public static void printObservedMinusProjected_Image_Frame_Eqn2() throws IOException, NotConvergedException {
        
        double[][] m = getObservedMinusProjected_Image_Frame_Eqn2();
        
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

    public static ImageExt getImage(int idx) throws IOException {
        if (idx < 1 || idx > 5) {
            throw new IllegalArgumentException("idx must be 1 through 5, inclusive");
        }

        String path = ResourceFinder.findTestResourcesDirectory() + sep + ZHANGDIR
                + sep + "image" + Integer.toString(idx) + ".gif";

        File f = new File(path);
        if (!f.exists()) {
            throw new IOException("could not find file at " + path);
        }

        ImageExt image = ImageIOHelper.readImageExt(path);

        return image;
    }
}
