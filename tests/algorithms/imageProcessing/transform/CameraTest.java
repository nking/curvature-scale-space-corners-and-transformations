package algorithms.imageProcessing.transform;

import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

/**
 *
 * @author nichole
 */
public class CameraTest extends TestCase {
    
    public CameraTest() {
    }

    public void test0() throws NotConvergedException {
        
        double focalLength = 500;
        double centerX = 250;
        double centerY = 200;
        double tol = 1e-5;
        double diff;
        
        double[][] x = new double[3][1];
        x[0] = new double[]{10};
        x[1] = new double[]{15};
        x[2] = new double[]{1};
        
        double[][] eX = new double[3][];
        eX[0] = new double[]{(x[0][0] - centerX)/focalLength};
        eX[1] = new double[]{(x[1][0] - centerY)/focalLength};
        eX[2] = new double[]{x[2][0]};
        
        
        //double[][] cIntr = Camera.createIntrinsicCameraMatrix(focalLength,
        //    centerX, centerY);
        
        double[] rCoeffs = null;
        boolean useR2R4 = false;
        double[][] xC = Camera.pixelToCameraCoordinates(x, rCoeffs, focalLength, 
             centerX, centerY, useR2R4);
        assertEquals(eX.length, xC.length);
        assertEquals(eX[0].length, xC[0].length);
        for (int i = 0; i < eX.length; ++i) {
            for (int j = 0; j < eX[i].length; j++) {
                diff = Math.abs(eX[i][j] - xC[i][j]);
                assertTrue(diff < tol);
            }
        }
        
        
        double[][] xCP = Camera.cameraToPixelCoordinates(xC, rCoeffs, 
            focalLength, centerX, centerY, useR2R4);
        
        assertEquals(x.length, xCP.length);
        assertEquals(x[0].length, xCP[0].length);
        
        
        for (int i = 0; i < xCP.length; ++i) {
            for (int j = 0; j < xCP[i].length; j++) {
                diff = Math.abs(xCP[i][j] - x[i][j]);
                assertTrue(diff < tol);
            }
        }
    }

}
