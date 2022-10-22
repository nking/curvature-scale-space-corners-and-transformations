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

            // remove distortion from x
            double[][] xU = CameraCalibration.removeRadialDistortion(x, radial[0], radial[1], true);

            double[][] expectedR = Zhang98Data.getRotation(i);
            double[] expectedT = Arrays.copyOf(Zhang98Data.getTranslation(i), 3);
            MatrixUtil.multiply(expectedT, 1./expectedT[2]);

            double[][] xc = Camera.pixelToCameraCoordinates(x, expectedKIntr, radial, true);

            Camera.CameraPoseParameters result = CameraPose.calculatePoseAndKUsingDLT(xU, xW);

            Camera.CameraExtrinsicParameters extr2 = CameraPose.calculatePoseUsingCameraCalibration(
                    new Camera.CameraIntrinsicParameters(expectedKIntr), xU,xW);

            Camera.CameraExtrinsicParameters extr = result.getExtrinsicParameters();
            Camera.CameraIntrinsicParameters intr = result.getIntrinsicParameters();
            double[] p3 = Arrays.copyOf(result.getP3(), result.getP3().length);

            /*
            tc1 = -1*(R^-1)*p3
            tc3 = p3
            ti2 = -1 * R^-1 * K^-1 * p3
            ti4 = K^-1 * p3.   <====== best match to Zhang's translation
            */
            double[] tc1 = MatrixUtil.multiplyMatrixByColumnVector(
                    MatrixUtil.pseudoinverseFullColumnRank(expectedR), p3);
            MatrixUtil.multiply(tc1,-1);
            double[] _tc1 = MatrixUtil.multiplyMatrixByColumnVector(
                    MatrixUtil.pseudoinverseFullColumnRank(extr.getRotation()), p3);
            MatrixUtil.multiply(_tc1,-1);
            double[] tc3 = Arrays.copyOf(p3, p3.length);
            double[] ti2 = Arrays.copyOf(extr.getTranslation(), 3);
            double[] ti4 = MatrixUtil.multiplyMatrixByColumnVector(
                    MatrixUtil.pseudoinverseFullColumnRank(expectedKIntr), p3);
            MatrixUtil.multiply(ti4,1./ti4[2]);
            double[] _ti4 = MatrixUtil.multiplyMatrixByColumnVector(
                    MatrixUtil.pseudoinverseFullColumnRank(intr.getIntrinsic()), p3);
            MatrixUtil.multiply(_ti4,1./_ti4[2]);

            double[] ti2_2 = Arrays.copyOf(extr2.getTranslation(), extr2.getTranslation().length);
            MatrixUtil.multiply(ti2_2,1./ti2_2[2]);
            double[] ti4_2 = MatrixUtil.multiplyMatrixByColumnVector(
                    MatrixUtil.pseudoinverseFullColumnRank(expectedR), extr2.getTranslation());
            MatrixUtil.multiply(ti4_2,-1);
            MatrixUtil.multiply(ti4_2,1./ti4_2[2]);

            double[][] r2Orth = MatrixUtil.copy(extr2.getRotation());
            r2Orth = Rotation.orthonormalizeUsingSVD(r2Orth);

            System.out.printf("%d) r=\n%s\n", i, FormatArray.toString(extr.getRotation(), "%.3e"));
            System.out.printf("%d) r2=\n%s\n", i, FormatArray.toString(extr2.getRotation(), "%.3e"));
            System.out.printf("%d) r2Orth=\n%s\n", i, FormatArray.toString(r2Orth, "%.3e"));
            System.out.printf("    r expected=\n%s\n", FormatArray.toString(expectedR, "%.3e"));
            System.out.printf("%d) kIntr=\n%s\n", i, FormatArray.toString(intr.getIntrinsic(), "%.3e"));
            System.out.printf("    kIntr expected=\n%s\n", FormatArray.toString(expectedKIntr, "%.3e"));
            System.out.printf("%d) tc1=\n%s\n", i, FormatArray.toString(tc1, "%.3e"));
            System.out.printf("%d) _tc1=\n%s\n", i, FormatArray.toString(_tc1, "%.3e"));
            System.out.printf("%d) tc3=\n%s\n", i, FormatArray.toString(tc3, "%.3e"));
            System.out.printf("%d) ti2=\n%s\n", i, FormatArray.toString(ti2, "%.3e"));
            System.out.printf("%d) ti4=\n%s\n", i, FormatArray.toString(ti4, "%.3e"));
            System.out.printf("%d) _ti4=\n%s\n", i, FormatArray.toString(_ti4, "%.3e"));
            System.out.printf("%d) ti2_2=\n%s\n", i, FormatArray.toString(ti2_2, "%.3e"));
            System.out.printf("%d) ti4_2=\n%s\n", i, FormatArray.toString(ti4_2, "%.3e"));
            System.out.printf("    t expected=\n%s\n", FormatArray.toString(expectedT, "%.3e"));

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
