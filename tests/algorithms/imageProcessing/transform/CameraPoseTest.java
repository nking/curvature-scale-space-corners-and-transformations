package algorithms.imageProcessing.transform;

import static algorithms.imageProcessing.transform.Rotation.extractThetaFromZYX;
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
    
    public void test() {}

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

}
