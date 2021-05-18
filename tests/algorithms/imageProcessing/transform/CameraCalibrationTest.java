package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Camera.CameraExtrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraMatrices;
import static algorithms.imageProcessing.transform.CameraCalibration.solveForIntrinsic;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

/**
 *
 * @author nichole
 */
public class CameraCalibrationTest extends TestCase {
    
    public CameraCalibrationTest() {
    }
    
    public void testCalibration0() throws IOException, NotConvergedException {
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
            System.out.printf("\nimg %d:\n", i);
            System.out.printf("   r1=\n%s\n", FormatArray.toString(ex1.getRotation(), "%.3e"));
            System.out.printf("   r2=\n%s\n", FormatArray.toString(ex2.getRotation(), "%.3e"));
            System.out.printf("   t1=\n%s\n", FormatArray.toString(ex1.getTranslation(), "%.3e"));
            System.out.printf("   t2=\n%s\n", FormatArray.toString(ex2.getTranslation(), "%.3e"));
        }
        
        double[] u = new double[nFeatures*nImages];
        double[] v = new double[nFeatures*nImages];
        CameraCalibration.calculateProjected(coordsW, h, u, v);
                
        kRadial = CameraCalibration.solveForRadialDistortion(coordsI, u, v, 
            cameraMatrices, useR2R4);
        double k1 = kRadial[0];
        double k2 = kRadial[1];
        System.out.printf("k=\n%s\n", FormatArray.toString(kRadial, "%.4e"));
        
        
        // quick test of apply and remove distortion
        // (u_d, v_d) are the distorted features in coordsI in image reference frame.
        // (x_d, y_d) are the distorted features in the camera reference frame.
        // (x, y) are the distortion-free features in the camera reference frame.
        int i, j;
        double diff, diffD;
        double[][] uvDI, xyDI, xyi, xyDUI;
        for (i = 0; i < nImages; ++i) {
            uvDI = MatrixUtil.copySubMatrix(coordsI, 0, 2, nFeatures*i, nFeatures*(i + 1)-1);
            xyDI = Camera.pixelToCameraCoordinates(uvDI, cameraMatrices.getIntrinsics(), null, false);
            xyi = CameraCalibration.removeRadialDistortion(xyDI, k1, k2);
            xyDUI = CameraCalibration.applyRadialDistortion(xyi, k1, k2, useR2R4);
            // assert that xy != xyDI
            // assert that xyDI == xyDUI
            for (j = 0; j < nFeatures; ++j) {
                diff = Math.abs(xyi[0][j] - xyDI[0][j]);
                assertTrue(diff > 0.1);
                diffD = Math.abs(xyDUI[0][j] - xyDI[0][j]);
                System.out.printf("(distorted-undistorted)=%.5e \n(orig - removedApplied)=%.5e\n", diff, diffD);  System.out.flush();
                assertTrue(diffD < 0.1);
            }
        }
        
        //CameraMatrices cameraCalibration = CameraCalibration.estimateCamera(
        //     nFeatures, coordsI, coordsW, useR2R4);
        
    }
    
    public void estCalibration1() throws IOException, NotConvergedException {
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
    
    public void estCalibration2() throws NotConvergedException {
        
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
        double[][] eR21 = Rotation.createEulerYawRotationMatrix(23.4*(Math.PI/180.));
        double[][] eR23 = Rotation.createEulerYawRotationMatrix(23.4*(Math.PI/180.));
        double[][] eR13 = Rotation.createEulerYawRotationMatrix(2*23.4*(Math.PI/180.));
        
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
}
