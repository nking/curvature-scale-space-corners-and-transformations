package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class CameraCalibrationTest extends TestCase {
    
    public CameraCalibrationTest() {
    }
    
    public void testCalibration0() {
        use data in ~/Downloads/zhang_data/
    }
    
    public void testCalibration2() {
        
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
        
        double[][] kIntr = c.getIntrinsic();
        assertNotNull(kIntr);
        assertEquals(3, kIntr.length);
        assertEquals(3, kIntr[0].length);
        
        double focalLengthX = kIntr[0][0];
        double focalLengthY = kIntr[1][1];
        double centerX = kIntr[0][2];
        double centerY = kIntr[1][2];
        double skew = kIntr[0][1];
        
        double[][] collatedRotation = c.getCollatedRotation;
       
        double[][] collatedTranslation = c.getCollatedTranslation();
        
        // --------
        double[][] coordsI0 = MatrixUtil.copySubMatrix(coordsI, 0, 3, 0, n);
        double[][] h0 = Camera.solveForHomography(coordsI0, coordsW);
        assertEquals(3, h0.length);
        assertEquals(3, h0[0].length);
        //TODO: assert characteristics of h0
        
    }
}
