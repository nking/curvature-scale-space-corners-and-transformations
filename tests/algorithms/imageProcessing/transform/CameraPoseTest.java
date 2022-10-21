package algorithms.imageProcessing.transform;

import static algorithms.imageProcessing.transform.Rotation.extractThetaFromZYX;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import junit.framework.TestCase;

import java.util.Arrays;

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
    public void testCalculatePoseUsingDLT() throws Exception {

        // see testresources/zhang1998/README.txt

        // they use f(r) = 1 + k1*r + k2*r^2:
        boolean useR2R4 = false;

        int nFeatures = 256;
        int nImages = 5;

        // 3 X 256
        double[][] xW = Zhang98Data.getFeatureWCS();
        assertEquals(3, xW.length);
        assertEquals(nFeatures, xW[0].length);

        double[][] expectedKIntr = Zhang98Data.getIntrinsicCameraMatrix();
        double[] radial = Zhang98Data.getRadialDistortionR2R4();

        for (int i = 1; i <= 5; ++i) {
            System.out.println();
            double[][] x = Zhang98Data.getObservedFeaturesInImage(i);
            double[][] expectedR = Zhang98Data.getRotation(i);
            double[] expectedT = Arrays.copyOf(Zhang98Data.getTranslation(i), 3);
            MatrixUtil.multiply(expectedT, 1./expectedT[2]);

            Camera.CameraParameters result = CameraPose.calculatePoseUsingDLT(x, xW);
            Camera.CameraExtrinsicParameters extr = result.getExtrinsicParameters();
            Camera.CameraIntrinsicParameters intr = result.getIntrinsicParameters();
            double[] t = Arrays.copyOf(extr.getTranslation(), 3);
            MatrixUtil.multiply(t, 1./t[2]);

            System.out.printf("%d) r=\n%s\n", i, FormatArray.toString(extr.getRotation(), "%.3e"));
            System.out.printf("%d) t=\n%s\n", i, FormatArray.toString(t, "%.3e"));
            System.out.printf("%d) kIntr=\n%s\n", i, FormatArray.toString(intr.getIntrinsic(), "%.3e"));
            System.out.printf("    r expected=\n%s\n", FormatArray.toString(expectedR, "%.3e"));
            System.out.printf("    t expected=\n%s\n", FormatArray.toString(expectedT, "%.3e"));
            System.out.printf("    kIntr expected=\n%s\n", FormatArray.toString(expectedKIntr, "%.3e"));
        }

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
