package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Camera.CameraExtrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraMatrices;
import algorithms.matrix.MatrixUtil;
import algorithms.statistics.Standardization;
import algorithms.util.FormatArray;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class CameraCalibrationTest extends TestCase {
    
    static double eps = 1.e-5;
    
    private static final Level LEVEL = Level.INFO;
    private static final Logger log;
    static {
        log = Logger.getLogger(CameraCalibration.class.getSimpleName());
        log.setLevel(LEVEL);
    }
    
    public CameraCalibrationTest() {
    }
    
    public void estApplyRemoveRadialDistortion() 
    throws Exception {
        
        // test for model #3
        // test for model #4
        // test for k1<0, barrel distortion (happens for smaller focal lengths)
        // test for k1>0, pincushion distortion (happes for larger focal lenngths)
        
        // generate coords in radial annuli in area 100x100 then unit standard normalized
        double[][] coords = generateCoords0();
        
        log.log(LEVEL, String.format("orig coords=\n%s\n", FormatArray.toString(coords, "%.3f")));
        
        double[] k1s = new double[]{0.25, +0.25};
        double[] k2s = new double[]{0.2, +0.2};
        boolean[] useR2R4s = new boolean[]{false, true};
        int i, j;
        double k1, k2;
        double[][] xyd, xy;
        double diffDx, diffDy, diffx, diffy;
        for (boolean useR2R4 : useR2R4s) {
            for (i = 0; i < k1s.length; ++i) {
                k1 = k1s[i];
                k2 = k2s[i];
                log.log(LEVEL, "useR2R4=" + useR2R4 + "   k1="+ k1);
                xyd = CameraCalibration.applyRadialDistortion(coords, k1, k2, useR2R4);
                xy = CameraCalibration.removeRadialDistortion(xyd, k1, k2, useR2R4);
                // assert xyd != coordsI
                // assert xy == coordsI
                for (j = 0; j < coords[0].length; ++j) {
                    diffDx = Math.abs(coords[0][j] - xyd[0][j]);
                    diffDy = Math.abs(coords[1][j] - xyd[1][j]);
                    diffx = Math.abs(coords[0][j] - xy[0][j]);
                    diffy = Math.abs(coords[1][j] - xy[1][j]);
                    log.log(LEVEL, String.format("(%.4e,%.4e): diff(%.4e, %.4e)?  same(%.4e, %.4e)?\n", 
                        coords[0][j], coords[1][j], diffDx, diffDy, diffx, diffy));
                    assertTrue(diffx < eps);
                    assertTrue(diffy < eps);
                    // if barrel (k1<0), distorted is at a smaller radius
                }
            }
        }
        
        // apply radial distortion and remove radial distortion and compare
    }
    
    public void testCalibration0() throws Exception {
        // see testresources/zhang1998/README.txt
        
        // they use f(r) = 1 + k1*r + k2*r^2:
        boolean useR2R4 = false;
        
        int nFeatures = 256;
        int nImages = 5;
        
        // 3 X 256
        double[][] coordsW = Zhang98Data.getFeatureWCS();
        assertEquals(3, coordsW.length);
        assertEquals(nFeatures, coordsW[0].length);
        
        //3 X (256*5)
        double[][] coordsI = Zhang98Data.getFeaturesInAllImages();
        assertEquals(3, coordsI.length);
        assertEquals(nFeatures*nImages, coordsI[0].length);
        
        double[][] h = CameraCalibration.solveForHomographies(
            coordsI, coordsW, nFeatures, nImages);
        
        //(2) using all homographies, solve for the camera intrinsic parameters
        // this is where at least 3 points are needed per image to euqal the number of unknown intrinsic parameters.
        Camera.CameraIntrinsicParameters kIntr = CameraCalibration.solveForIntrinsic(h);
        
        double alphaE = 871.445;
        double gammaE = 0.2419;
        double u0E = 300.7676;
        double betaE = 871.1251;
        double v0E = 220.8684;
        double k1E = 0.1371;
        double k2E = -2.0101;
        
        double alpha = kIntr.getIntrinsic()[0][0];
        double gamma = kIntr.getIntrinsic()[0][1];
        double u0 = kIntr.getIntrinsic()[0][2];
        double beta = kIntr.getIntrinsic()[1][1];
        double v0 = kIntr.getIntrinsic()[1][2];
        double[] kRadial = null;
        
        assertTrue(Math.abs(alphaE - alpha) < 0.1);
        assertTrue(Math.abs(gammaE - gamma) < 0.1);
        assertTrue(Math.abs(u0E - u0) < 0.1);
        assertTrue(Math.abs(betaE - beta) < 0.1);
        assertTrue(Math.abs(v0E - v0) < 0.1);
        
        CameraMatrices cameraMatrices = new CameraMatrices();
        cameraMatrices.setIntrinsics(kIntr);
        
        List<Camera.CameraExtrinsicParameters> extrinsics = 
            CameraCalibration.solveForExtrinsics(kIntr, h, nImages);
        cameraMatrices.getExtrinsics().addAll(extrinsics);
        
        List<Camera.CameraExtrinsicParameters> extrinsics2 = CameraCalibration.solveForExtrinsics2(
            kIntr, coordsI, coordsW, useR2R4);
        
        CameraExtrinsicParameters ex1, ex2;
        for (int i = 0; i < nImages; ++i) {
            ex1 = extrinsics.get(i);
            ex2 = extrinsics2.get(i);
            log.log(LEVEL, String.format("\nimg %d:\n", i));
            log.log(LEVEL, String.format("   r1=\n%s\n", FormatArray.toString(ex1.getRotation(), "%.3e")));
            log.log(LEVEL, String.format("   r2=\n%s\n", FormatArray.toString(ex2.getRotation(), "%.3e")));
            log.log(LEVEL, String.format("   t1=\n%s\n", FormatArray.toString(ex1.getTranslation(), "%.3e")));
            log.log(LEVEL, String.format("   t2=\n%s\n", FormatArray.toString(ex2.getTranslation(), "%.3e")));
        }
        
        double[] u = new double[nFeatures*nImages];
        double[] v = new double[nFeatures*nImages];
        CameraCalibration.calculateProjected(coordsW, h, u, v);
                
        kRadial = CameraCalibration.solveForRadialDistortion(coordsI, u, v, 
            cameraMatrices, useR2R4);
        double k1 = kRadial[0];
        double k2 = kRadial[1];
        log.log(LEVEL, String.format("k=\n%s\n", FormatArray.toString(kRadial, "%.4e")));
        
        // quick test of apply and remove distortion
        // (u_d, v_d) are the distorted features in coordsI in image reference frame.
        // (x_d, y_d) are the distorted features in the camera reference frame.
        // (x, y) are the distortion-free features in the camera reference frame.
        int i, j;
        double diff, diffD;
        double[][] uvDI, xyDI, xyi, xyDUI;
        
        double[][] cIntrE = new double[3][3];
        cIntrE[0] = new double[]{alphaE, gammaE, u0E};
        cIntrE[1] = new double[]{0, betaE, v0E};
        cIntrE[2] = new double[]{0, 0, 1};
        useR2R4 = false;
        Camera.CameraIntrinsicParameters kIntrE = new Camera.CameraIntrinsicParameters(
            cIntrE);
        for (i = 0; i < nImages; ++i) {
            uvDI = MatrixUtil.copySubMatrix(coordsI, 0, 2, nFeatures*i, nFeatures*(i + 1)-1);
            
            xyDI = Camera.pixelToCameraCoordinates(uvDI, kIntrE, null, useR2R4);
            xyi = CameraCalibration.removeRadialDistortion(xyDI, k1E, k2E);
            xyDUI = CameraCalibration.applyRadialDistortion(xyi, k1E, k2E, useR2R4);
            // assert that xy != xyDI
            // assert that xyDI == xyDUI
            for (j = 0; j < nFeatures; ++j) {
                diff = Math.abs(xyi[0][j] - xyDI[0][j]);
                //assertTrue(diff > 0.1);
                diffD = Math.abs(xyDUI[0][j] - xyDI[0][j]);
                log.log(LEVEL, String.format("(%.3f,%.3f): (distorted-undistorted)=%.5e \n     (orig - removedApplied)=%.5e\n", 
                    xyDI[0][j], xyDI[1][j], diff, diffD));
                System.out.flush();
                //assertTrue(diffD < 0.1);
            }
        }
        
        //CameraMatrices cameraCalibration = CameraCalibration.estimateCamera(
        //     nFeatures, coordsI, coordsW, useR2R4);
        
    }
    
    public void estCalibration1() throws Exception {
        // see testresources/zhang1998/README.txt
        
        // they use f(r) = 1 + k1*rk2*r^2:
        boolean useR2R4 = false;
        
        int nFeatures = 256;
        int nImages = 5;
        
        // 3 X 256
        double[][] coordsW = Zhang98Data.getFeatureWCS();
        assertEquals(3, coordsW.length);
        assertEquals(nFeatures, coordsW[0].length);
        
        //3 X (256*5)
        double[][] coordsI = Zhang98Data.getFeaturesInAllImages();
        assertEquals(3, coordsI.length);
        assertEquals(nFeatures*nImages, coordsI[0].length);
        
        double alphaE = 871.445;
        double gammaE = 0.2419;
        double u0E = 300.7676;
        double betaE = 871.1251;
        double v0E = 220.8684;
        double k1E = 0.1371;
        double k2E = -2.0101;
        
        
        CameraMatrices cameraCalibration = CameraCalibration.estimateCamera(
             nFeatures, coordsI, coordsW, useR2R4);
        
    }
    
    public void estCalibration2() throws Exception {
        
        boolean useR2R4 = false;
        
        // number of features
        int n = 8;
        double[][] coordsI = new double[3][];
        coordsI[0] = new double[]{};
        coordsI[1] = new double[]{};
        coordsI[2] = new double[]{1, 1, 1, 1, 1, 1, 1, 1};
        
        double[][] coordsW = new double[3][];
        coordsW[0] = new double[]{-11, 11, 0, 0, -3.7, -8, -11, 11};
        coordsW[1] = new double[]{14,  14, 9.7, 0, -3, -3, -14, -14};
        coordsW[2] = new double[]{41.5, 41.5, 41.5, 41.5, 41.5, 41.5, 41.5, 41.5};
        
        CameraMatrices c = CameraCalibration.estimateCamera(n, coordsI, coordsW,
            useR2R4);
        
        /*
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
        */
        double eFocalLength = 1604;
        double eCenterX = 1521;
        double eCenterY = 1752;
        double[][] eR21 = Rotation.createYawRotationMatrix(23.4*(Math.PI/180.));
        double[][] eR23 = Rotation.createYawRotationMatrix(23.4*(Math.PI/180.));
        double[][] eR13 = Rotation.createYawRotationMatrix(2*23.4*(Math.PI/180.));
        
        double[] eT21 = new double[]{18, 0, 0};
        double[] eT23 = new double[]{18, 0, 0};
        double[] eT13 = new double[]{36, 0, 0};
        
        assertNotNull(c);
        
        double[][] kIntr = c.getIntrinsics().getIntrinsic();
        assertNotNull(kIntr);
        assertEquals(3, kIntr.length);
        assertEquals(3, kIntr[0].length);
        
        double focalLengthX = kIntr[0][0];
        double focalLengthY = kIntr[1][1];
        double centerX = kIntr[0][2];
        double centerY = kIntr[1][2];
        double skew = kIntr[0][1];
        /*
        double[][] collatedRotation = c.getCollatedRotation;
       
        double[][] collatedTranslation = c.getCollatedTranslation();
        
        // --------
        double[][] coordsI0 = MatrixUtil.copySubMatrix(coordsI, 0, 3, 0, n);
        double[][] h0 = Camera.solveForHomography(coordsI0, coordsW);
        assertEquals(3, h0.length);
        assertEquals(3, h0[0].length);
        //TODO: assert characteristics of h0
        */
    }

    /**
     * generate coordinates in 5 radial annuli in area 100x100 then unit standard normalized
     * @return 
     */
    private double[][] generateCoords0() {
        TDoubleList x = new TDoubleArrayList();
        TDoubleList y = new TDoubleArrayList();
        
        double rMax = 100;
        int nAnnuli = 1;
        double dR = rMax/nAnnuli;
        
        // for each annuli, want the 4 points crossing the x and y axes to test 0's.
        //     and the remaining points can be 8 points => all distributed at 2*pi/12 intervals
        double dT = 2.*Math.PI/12.;
        
        int i, j;
        double r, t;
        x.add(0);
        y.add(0);
        for (i = 1; i <= nAnnuli; ++i) {
            r = i*dR;
            for (j = 0; j < 12; ++j) {
                t = j*dT;
                //x = r cos t
                //y = r sin t
                x.add(r*Math.cos(t));
                y.add(r*Math.sin(t));
            }
        }
        
        // can use this transposed for unit standard normalization:
        //public static double[][] Standardization.standardUnitNormalization(double[][] data, 
        //    double[] outputMean, double[] outputStandardDeviation) {
        
        // write all x along column 0 instead of row 0.  and y goes in column 1
        double[][] normalized = new double[x.size()][2];
        for (i = 0; i < x.size(); ++i) {
            normalized[i] = new double[]{x.get(i), y.get(i)};
        }
        
        double[] outputMean = new double[2];
        double[] outputStDev = new double[2];
        normalized = Standardization.standardUnitNormalization(normalized,
            outputMean, outputStDev);
        
        // write for use in code above:
        //    row 0 = 's. row 1 = y's, row 2 = 1's
        double[][] coords = MatrixUtil.zeros(3, x.size());
        for (i = 0; i < x.size(); ++i) {
            coords[0][i] = normalized[i][0];
            coords[1][i] = normalized[i][1];
            coords[2][i] = 1.;
        }
        
        return coords;
    }
}
