package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;

import java.io.IOException;
import java.util.Arrays;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class TriangulationTest extends TestCase {
    
    public TriangulationTest() {
    }

    /**
     * Test of calculateWCSPoint method, of class Triangulation.
     */
    public void testCalculateWCSPoint() throws IOException, NotConvergedException {
        
        //test data from:
        //http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example5.html
        // now at http://robots.stanford.edu/cs223b04/JeanYvesCalib/
        // "Fifth calibration example - Calibrating a stereo system, stereo image rectification and 3D stereo triangulation"
        // by Jean-Yves Bouguet
        // Camera Calibration Toolbox for Matlab
        //    saved as pdf in miscNotes/bouguetj_5th_calibration_example.pdf
        //
        // left camera: 
        //    focal length = 533.5, 533.5
        //    cc = 341.6, 235.2
        //    skew = 0, 0
        //    radial distortion k = -0.288, 0.097, 0.001, -0.0003, 0
        //
        // right camera: 
        //    focal length = 536.8, 536.5
        //    cc = 326.3, 250.1
        //    skew = 0, 0
        //    radial distortion k = -0.289, 0.107, 0.001, -0.0001, 0
        //
        // rotation vector om=0.00669, 0.00452, -0.0035
        // translation vector t = -99.80198, 1.12443, 0.05041
        //
        // note: the checkerboard squares are 30mm in WCS metric
        //
        // note: radial distortion should be corrected:  use on the original coordinates:
        //     x_corrected = x*(1 + k1*r^2 + k2r^4) where r is distance of point from cc.
        
        double[] k1IntrTrans = new double[]{-341.6, -234.3, 0};
        double[] k2IntrTrans = new double[]{-326.5, -249.3, 0};
        double[] k1IntrTransInv = Arrays.copyOf(k1IntrTrans, k1IntrTrans.length);
        MatrixUtil.multiply(k1IntrTransInv, -1);
        double[] k2IntrTransInv = Arrays.copyOf(k2IntrTrans, k2IntrTrans.length);
        MatrixUtil.multiply(k2IntrTransInv, -1);
        
        double[][] k1Intr = Camera.createIntrinsicCameraMatrix(533.07, k1IntrTrans[0], k1IntrTrans[1]);
        double[][] k2Intr = Camera.createIntrinsicCameraMatrix(536.7, k2IntrTrans[0], k2IntrTrans[1]);
        double[][] k2IntrInv = Camera.createIntrinsicCameraMatrixInverse(536.7, k2IntrTrans[0], k2IntrTrans[1]);
        
        double[][] k1ExtrRot = MatrixUtil.createIdentityMatrix(3);
        double[] k1ExtrTrans = new double[]{0, 0, 0};
        
        double[][] k2ExtrRot = Rotation.createRodriguesFormulaRotationMatrix(
            new double[]{0.00611, 0.00409, -0.00359});
        double[] k2ExtrTrans = new double[]{-99.85, 0.82, 0.44};
        
        System.out.printf("k1ExtrRot\n=%s\n", FormatArray.toString(k1ExtrRot, "%.3e"));
        System.out.printf("k1ExtrTrans\n=%s\n", FormatArray.toString(k1ExtrTrans, "%.3e"));
        System.out.printf("k2ExtrRot\n=%s\n", FormatArray.toString(k2ExtrRot, "%.3e"));
        System.out.printf("k2ExtrTrans\n=%s\n", FormatArray.toString(k2ExtrTrans, "%.3e"));
        
        //choosing a point that is in left01.jpg and right01.jpg
        //(307, 159)  (184, 172)
        double[][] x1 = new double[3][1];
        x1[0] = new double[]{307};
        x1[1] = new double[]{159};
        x1[2] = new double[]{1};
        
        double[][] x2 = new double[3][1];
        x2[0] = new double[]{184};
        x2[1] = new double[]{172};
        x2[2] = new double[]{1};

        Triangulation.WCSPt wcsPt = Triangulation.calculateWCSPoint(
            k1Intr, k1ExtrRot, k1ExtrTrans,
            k2Intr, k2ExtrRot, k2ExtrTrans,
            x1, x2);

        double[] xw = wcsPt.X;
        System.out.printf("\nxw=%s\n", FormatArray.toString(xw, "%.3e"));
        MatrixUtil.multiply(xw, 1./xw[xw.length - 1]);
        System.out.printf("   =%s\n", FormatArray.toString(xw, "%.3e"));
        assertTrue(Math.abs(Math.abs(xw[2]) - 425) < 150);

        double[][] camera1 = Camera.createCamera(k1Intr, k1ExtrRot, k1ExtrTrans);
        double[][] camera2 = Camera.createCamera(k2Intr, k2ExtrRot, k2ExtrTrans);
        double[][] camera1Inv = Camera.createCameraInverse(k1Intr, k1ExtrRot, k1ExtrTrans);
        double[][] camera2Inv = Camera.createCameraInverse(k2Intr, k2ExtrRot, k2ExtrTrans);
        
        System.out.printf("\nxw=%s\n", FormatArray.toString(xw, "%.3e"));

        double[] expectedx1 = MatrixUtil.multiplyMatrixByColumnVector(camera1, xw);
        double[] expectedx2 = MatrixUtil.multiplyMatrixByColumnVector(camera2, xw);
        MatrixUtil.multiply(expectedx1, 1./expectedx1[expectedx1.length - 1]);
        MatrixUtil.multiply(expectedx2, 1./expectedx2[expectedx2.length - 1]);

        System.out.printf("x1=%s\n", FormatArray.toString(MatrixUtil.transpose(x1), "%.3e"));
        System.out.printf("x2=%s\n", FormatArray.toString(MatrixUtil.transpose(x2), "%.3e"));
        System.out.printf("(camera1*xw)=%s\n", FormatArray.toString(expectedx1, "%.3e"));
        System.out.printf("(camera2*xw)=%s\n", FormatArray.toString(expectedx2, "%.3e"));

        //-----------
        
        //a correspondence closer to middle of image 1: (345, 188)  (215,238)
        x1 = new double[3][1];
        x1[0] = new double[]{345};
        x1[1] = new double[]{188};
        x1[2] = new double[]{1};
        
        x2 = new double[3][1];
        x2[0] = new double[]{215};
        x2[1] = new double[]{238};
        x2[2] = new double[]{1};

        wcsPt = Triangulation.calculateWCSPoint(
            k1Intr, k1ExtrRot, k1ExtrTrans,
            k2Intr, k2ExtrRot, k2ExtrTrans,
            x1, x2);

        xw = wcsPt.X;
        System.out.printf("\nxw=%s\n", FormatArray.toString(xw, "%.3e"));
        MatrixUtil.multiply(xw, 1./xw[xw.length - 1]);
        System.out.printf("   =%s\n", FormatArray.toString(xw, "%.3e"));

        // those were image coordinates.  use camera coordinates and extrinsic camera matrix:
        // camera2 as P = [ R | t] not P = [ R | -R*t]
        // the authors used XXc = Rc_1 * XX + Tc_1
        // (see http://robots.stanford.edu/cs223b04/JeanYvesCalib/htmls/parameters.html)

        //    radial distortion k = -0.289, 0.107, 0.001, -0.0001, 0
        double[] r1 = new double[]{-0.288, 0.097};
        double[] r2 = new double[]{-0.289, 0.107};
        double[][] p1 = Camera.createExtrinsicCameraMatrix(k1ExtrRot, k1ExtrTrans);
        double[][] p2 = Camera.createExtrinsicCameraMatrix(k2ExtrRot, k2ExtrTrans);

        double[][] x1C = Camera.pixelToCameraCoordinates(x1, k1Intr, r1, true);
        double[][] x2C = Camera.pixelToCameraCoordinates(x2, k2Intr, r2, true);

        Triangulation.WCSPt wcsPt2 = Triangulation.calculateWCSPoint(p1, p2, x1C, x2C);

        xw = wcsPt.X;
        System.out.printf("\nc: xw=%s\n", FormatArray.toString(xw, "%.3e"));
        //xw=6.670e-03, -5.187e-02, -9.986e-01, 2.102e-03
        // xw[0]/xw[3], xw[1]/xw[3], xw[2]/xw[3] = (3.173168411037107, -24.676498572787818, -475.07136060894385)

        MatrixUtil.multiply(xw, 1./xw[xw.length - 1]);
        System.out.printf("      =%s\n", FormatArray.toString(xw, "%.3e"));

        assertTrue(Math.abs(Math.abs(xw[2]) - 425) < 150);
    }
}
