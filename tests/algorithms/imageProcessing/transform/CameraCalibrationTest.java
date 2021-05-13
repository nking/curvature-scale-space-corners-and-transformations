package algorithms.imageProcessing.transform;

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
        
        double[] u = new double[nFeatures*nImages];
        double[] v = new double[nFeatures*nImages];
        CameraCalibration.calculateProjected(coordsW, h, u, v);
                
        double[] kRadial = CameraCalibration.estimateRadialDistortion(coordsI, u, v, cameraMatrices);
        double k1 = kRadial[0];
        double k2 = kRadial[1];
        System.out.printf("k=\n%s\n", FormatArray.toString(kRadial, "%.4e"));
    }
    
    public void estCalibration2() throws NotConvergedException {
        
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
        
        CameraMatrices c = CameraCalibration.estimateCamera(n, coordsI, coordsW);
        
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