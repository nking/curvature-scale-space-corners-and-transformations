package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class CameraPoseTest extends TestCase {
    
    public CameraPoseTest() {
    }
            
    /**
     * Test of calculateUsingEssentialMatrix method, of class CameraPose.
     */
    public void testCalculateUsingEssentialMatrix_4args() throws Exception {
        System.out.println("calculateUsingEssentialMatrix");
       
        double[][] k1 = Zhang98Data.getIntrinsicCameraMatrix();
        double[][] k2 = MatrixUtil.copy(k1);
        //x1, x2 size is 3 X 256
        double[][] x1 = Zhang98Data.getObservedFeaturesInImage(1);
        double[][] x2 = Zhang98Data.getObservedFeaturesInImage(5);
        
        Camera.CameraExtrinsicParameters[] expResult = null;
        
        Camera.CameraExtrinsicParameters[] result = 
            CameraPose.calculateUsingEssentialMatrix(k1, k2, x1, x2);
        
        assertNotNull(result);
        for (Camera.CameraExtrinsicParameters ep : result) {
            System.out.printf("\nresult:\nrot=%strans=%s\n", 
                FormatArray.toString(ep.getRotation(), "%.4e"),
                FormatArray.toString(ep.getTranslation(), "%.4e"));
        }
        
        System.out.printf("\nimg1:\nrot=%strans=%s\n", 
                FormatArray.toString(Zhang98Data.getRotation(1), "%.4e"),
                FormatArray.toString(Zhang98Data.getTranslation(1), "%.4e"));
        System.out.printf("\nimg5:\nrot=%strans=%s\n", 
                FormatArray.toString(Zhang98Data.getRotation(5), "%.4e"),
                FormatArray.toString(Zhang98Data.getTranslation(5), "%.4e"));
        
        double[][] diffRSameCenter = Rotation.procrustesAlgorithmForRotation(
            Zhang98Data.getRotation(1), Zhang98Data.getRotation(5));
        
        System.out.printf("\ndifference in rot between img1 and img5=\n%s\n", 
           FormatArray.toString(diffRSameCenter, "%.4e"));
    }

    /**
     * Test of calculateUsingEssentialMatrix method, of class CameraPose.
     */
    public void estCalculateUsingEssentialMatrix_5args() throws Exception {
        System.out.println("calculateUsingEssentialMatrix");
        double[][] k1 = null;
        double[][] k2 = null;
        double[][] x1 = null;
        double[][] x2 = null;
        double[][] outputXW = null;
        Camera.CameraExtrinsicParameters[] expResult = null;
        Camera.CameraExtrinsicParameters[] result = CameraPose.calculateUsingEssentialMatrix(k1, k2, x1, x2, outputXW);
        assertArrayEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of calculatePoseUsingDLT method, of class CameraPose.
     */
    public void estCalculatePoseUsingDLT() throws Exception {
        System.out.println("calculatePoseUsingDLT");
        double[][] x = null;
        double[][] X = null;
        Camera.CameraParameters expResult = null;
        Camera.CameraParameters result = CameraPose.calculatePoseUsingDLT(x, X);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of calculatePoseUsingCameraCalibration method, of class CameraPose.
     */
    public void estCalculatePoseUsingCameraCalibration() throws Exception {
        System.out.println("calculatePoseUsingCameraCalibration");
        Camera.CameraIntrinsicParameters intrinsics = null;
        double[][] x = null;
        double[][] X = null;
        Camera.CameraExtrinsicParameters expResult = null;
        Camera.CameraExtrinsicParameters result = CameraPose.calculatePoseUsingCameraCalibration(intrinsics, x, X);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of calculatePoseUsingPNP method, of class CameraPose.
     */
    public void estCalculatePoseUsingPNP() throws Exception {
        System.out.println("calculatePoseUsingPNP");
        Camera.CameraIntrinsicParameters intrinsics = null;
        double[][] x = null;
        double[][] X = null;
        Camera.CameraExtrinsicParameters expResult = null;
        Camera.CameraExtrinsicParameters result = CameraPose.calculatePoseUsingPNP(intrinsics, x, X);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
    
}
