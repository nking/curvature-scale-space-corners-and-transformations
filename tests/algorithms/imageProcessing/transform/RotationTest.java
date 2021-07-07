
package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import java.util.Arrays;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

/**
 *
 * @author nichole
 */
public class RotationTest extends TestCase {
    
    public RotationTest() {
    }

    public void estRodriguesFormula() {
        
        //from http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example.html
        
        double[] axis = new double[]{-1.451113, -1.827059, -0.179105};
        double[][] r = Rotation.createRodriguesFormulaRotationMatrix(axis);
        
        double[][] expected = new double[3][3];
        expected[0] = new double[]{-0.043583, 0.875946, -0.480436};
        expected[1] = new double[]{0.765974, 0.338032, 0.546825};
        expected[2] = new double[]{0.641392, -0.344170, -0.685684};
        
        //System.out.printf("r=\n%s\n", FormatArray.toString(r,"%.3e"));
        //System.out.printf("expected=\n%s\n", FormatArray.toString(expected,"%.3e"));
        
        assertEquals(expected.length, r.length);
        double tol = 1e-4;
        double diff;
        int i, j;
        for (i = 0; i < r.length; ++i) {
            for (j = 0; j < r[i].length; ++j) {
                diff = Math.abs(expected[i][j] - r[i][j]);
                assertTrue(diff < tol);
            }
        }
        
        double[] axis1 = MatrixUtil.normalizeL2(axis);
        double[] axis2 = Rotation.extractRodriguesRotationAxis(r);
        axis2 = MatrixUtil.normalizeL2(axis2);
        assertEquals(axis.length, axis2.length);
        System.out.printf("axis1=%s\n", FormatArray.toString(axis1, "%.3e"));
        System.out.printf("axis2=%s\n", FormatArray.toString(axis2, "%.3e"));
        
        
        for (i = 0; i < axis.length; ++i) {
            diff = Math.abs(axis[i] - axis2[i]);
            assertTrue(diff < tol);
        }
       
    }

    public void estProcrustesAlgorithm() throws NotConvergedException {
        double[][] a = new double[4][2];
        a[0] = new double[]{1, 2};
        a[1] = new double[]{3, 4};
        a[2] = new double[]{5, 6};
        a[3] = new double[]{7, 8};
        double[][] b = new double[4][2];
        b[0] = new double[]{1.2, 2.1};
        b[1] = new double[]{2.9, 4.3};
        b[2] = new double[]{5.2, 6.1};
        b[3] = new double[]{6.8, 8.1};
        double[][] q = new double[2][2];
        q[0] = new double[]{0.9999, -0.0126};
        q[1] = new double[]{0.0126, 0.9999};
        double[][] orthogRot = Rotation.procrustesAlgorithmForRotation(a, b);
        TestCase.assertEquals(q.length, orthogRot.length);
        double diff;
        double tol = 1.0E-4;
        for (int i = 0; i < q.length; ++i) {
            TestCase.assertEquals(q[i].length, orthogRot[i].length);
            for (int j = 0; j < q[i].length; ++j) {
                diff = Math.abs(q[i][j] - orthogRot[i][j]);
                TestCase.assertTrue(diff < tol);
            }
        }
    }

    public void estProcrustesAlgorithm2() throws NotConvergedException {
        //test data from:
        //http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example5.html
        // "Fifth calibration example - Calibrating a stereo system, stereo image rectification and 3D stereo triangulation"
        // by Jean-Yves Bouguet
        // Camera Calibration Toolbox for Matlab
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
        //check: Xc_1_right = R * Xc_1_left + T
        double[][] k1Intr = Camera.createIntrinsicCameraMatrix(533.07, 341.6, 234.3);
        double[][] k2Intr = Camera.createIntrinsicCameraMatrix(536.7, 326.5, 249.3);
        double[][] k1ExtrRot = MatrixUtil.createIdentityMatrix(3);
        double[] k1ExtrTrans = new double[]{0, 0, 0};
        double[][] k2ExtrRot = Rotation.createRodriguesFormulaRotationMatrix(new double[]{0.00611, 0.00409, -0.00359});
        double[] k2ExtrTrans = new double[]{-99.85, 0.82, 0.44};
        System.out.printf("k1ExtrRot\n=%s\n", FormatArray.toString(k1ExtrRot, "%.3e"));
        System.out.printf("k1ExtrTrans\n=%s\n\n", FormatArray.toString(k1ExtrTrans, "%.3e"));
        System.out.printf("k2ExtrRot\n=%s\n", FormatArray.toString(k2ExtrRot, "%.3e"));
        System.out.printf("k2ExtrTrans\n=%s\n\n", FormatArray.toString(k2ExtrTrans, "%.3e"));
        System.out.printf("expected k2ExtrRot from Rodrigues formula\n=%s\n", FormatArray.toString(k2ExtrRot, "%.3e"));
        /*
        [junit] =1.000e+00, 3.577e-03, 4.101e-03
        [junit] -3.602e-03, 1.000e+00, 6.103e-03
        [junit] -4.079e-03, -6.117e-03, 1.000e+00
         */
        double[][] x1 = new double[3][10];
        x1[0] = new double[]{129, 160, 140, 226, 225, 232, 341, 407, 532, 527};
        x1[1] = new double[]{145, 319, 361, 391, 61, 284, 289, 122, 48, 302};
        x1[2] = new double[]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        double[][] x2 = new double[3][10];
        x2[0] = new double[]{76, 110, 84, 164, 110, 124, 218, 275, 401, 402};
        x2[1] = new double[]{168, 331, 372, 401, 78, 295, 302, 134, 52, 318};
        x2[2] = new double[]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        int i;
        int j;
        double a;
        for (i = 0; i < x2.length; ++i) {
            for (j = 0; j < x2[i].length; ++j) {
                x2[i][j] -= k2ExtrTrans[i];
            }
        }
        double[][] orthogRot = Rotation.procrustesAlgorithmForRotation(MatrixUtil.transpose(x1), MatrixUtil.transpose(x2));
        System.out.printf("rotation from Procrustes algorithm\n=%s\n", FormatArray.toString(orthogRot, "%.3e"));
    }
    
    public void est2() {
        
        double[] theta = new double[]{25.*Math.PI/180., 35.*Math.PI/180., 55.*Math.PI/180.};
        
        //System.out.printf("original=%s\n", FormatArray.toString(theta, "%.3e"));
        
        double[][] r = Rotation.createRotationZYX(theta);
        
        double[] d = Rotation.extractThetaFromZYX(r);
        //System.out.printf("result=%s\n", FormatArray.toString(d, "%.3e"));

        assertEquals(theta.length, d.length);
        
        double diff;
        double tol = 1e-4;
        int i;
        for (i = 0; i < theta.length; ++i) {
            diff = Math.abs(theta[i] - d[i]);
            //System.out.println("diff=" + diff);
            assertTrue(diff < tol);
        }
        
        //========
        theta = new double[]{-25.*Math.PI/180., -35.*Math.PI/180., -55.*Math.PI/180.};
        
        //System.out.printf("original2=%s\n", FormatArray.toString(theta, "%.3e"));
        
        r = Rotation.createRotationZYX(theta);
        
        d = Rotation.extractThetaFromZYX(r);
        //System.out.printf("result2=%s\n", FormatArray.toString(d, "%.3e"));

        assertEquals(theta.length, d.length);
        
        for (i = 0; i < theta.length; ++i) {
            diff = Math.abs(theta[i] - d[i]);
            //System.out.println("diff=" + diff);
            assertTrue(diff < tol);
        }
    }
    
    public void test3() throws NotConvergedException {
        System.out.println("test3");
        
        /*
        at end of this method, add a test for gimbal lock.
        start with theta=[0, 90 degrees, 0]
        
        For the asymmetric matrix
            C_i_j_k = R(ax,ay,az)=R_z(az)*R_y(ay)*R_x(ax), singularity when
                ay=+-(pi/2)
        
        For the symmetrix 
            C_i_j_i = R_x(az)*R_y(ay)*R_x(ax), singularity when
                ay=+-pi or a2=0.
        */
        int i, j;
        
        double[] theta0 
            = new double[]{25.*Math.PI/180., 35.*Math.PI/180., 55.*Math.PI/180.};
        // 0.1 radians = 5.73 degrees.  1 degree = 0.0174 redians
        double[] dTheta0 
            = new double[]{2*Math.PI/180., 3*Math.PI/180., 5.5*Math.PI/180.};
        double[] theta0Up = new double[3];
        for (i = 0; i < 3; ++i) {
            theta0Up[i] = theta0[i] + dTheta0[i];
        }
        System.out.printf("theta0=\n%s\n", FormatArray.toString(theta0, "%.3e"));
        System.out.printf("dTheta0=\n%s\n", FormatArray.toString(dTheta0, "%.3e"));
        System.out.printf("theta0Up=\n%s\n", FormatArray.toString(theta0Up, "%.3e"));
        
        //System.out.printf("original=%s\n", FormatArray.toString(theta, "%.3e"));
        
        double[][] r0ZYX = Rotation.createRotationZYX(theta0);                
        double[] theta0ExZYX = Rotation.extractThetaFromZYX(r0ZYX); // recovered exactly     
        double[][] r0UpZYX = Rotation.applyRotationPerturbationZYX(theta0, dTheta0);
        double[] theta0UpExZYX = Rotation.extractThetaFromZYX(r0UpZYX);
        double[][] r0FromTheta0UpZYX = Rotation.createRotationZYX(theta0Up);
        double[] theta0UpExZYXMinusTheta0 = new double[3];
        for (i = 0; i < 3; ++i) {
            theta0UpExZYXMinusTheta0[i] = theta0UpExZYX[i] - theta0[i];
        }
        double[][] r0UpDanping = Rotation.updateRotation(dTheta0, r0ZYX);
        System.out.printf("r0ZYX=\n%s\n", FormatArray.toString(r0ZYX, "%.3e"));
        System.out.printf("r0UpZYX=\n%s\n", FormatArray.toString(r0UpZYX, "%.3e"));
        System.out.printf("r0FromTheta0UpZYX=\n%s\n", FormatArray.toString(r0FromTheta0UpZYX, "%.3e"));
        System.out.printf("r0UpDanping=\n%s\n", FormatArray.toString(r0UpDanping, "%.3e"));
        System.out.printf("theta0ExZYX=\n%s\n", FormatArray.toString(theta0ExZYX, "%.3e"));
        System.out.printf("theta0UpExZYX=\n%s\n", FormatArray.toString(theta0UpExZYX, "%.3e"));
        System.out.printf("theta0UpExZYXMinusTheta0=\n%s\n", FormatArray.toString(theta0UpExZYXMinusTheta0, "%.3e"));
                
        double[] q0Barfoot = Rotation.createQuaternionZYXFromEuler(theta0);
        double[] dPhi = Rotation.createRotationVector(theta0, dTheta0);
        double[] q0Hamilton = Rotation.createHamiltonQuaternionZYX(theta0);
        double[] q0Barfoot2 = Rotation.convertHamiltonToBarfootQuaternion(q0Hamilton);
        
        double[] a0ZYX1 = Rotation.extractRotationAxisFromZXY(r0ZYX);
        double[] a0ZYX2 = Rotation.extractRodriguesRotationAxis(r0ZYX);
        double m0ZYX1 = MatrixUtil.lPSum(a0ZYX1, 2);
        double m0ZYX2 = MatrixUtil.lPSum(a0ZYX2, 2);
        double[] q0AABarfootZYX1 = Rotation.createUnitLengthQuaternionBarfoot(a0ZYX1, m0ZYX1);
        double[] q0AABarfootZYX2 = Rotation.createUnitLengthQuaternionBarfoot(a0ZYX2, m0ZYX2);
        
        // rotate a quaternion.  
        double[] dq0Hamilton = Rotation.createHamiltonQuaternionZYX(dTheta0);
        double[] dq0Barfoot = Rotation.convertHamiltonToBarfootQuaternion(dq0Hamilton);
        double[] q0Up = Rotation.rotateVectorByQuaternion4(dq0Barfoot, q0Barfoot2);
        double[] q0SafeUp = Rotation.applySingularitySafeRotationPerturbationQuaternion(theta0, dTheta0);
        
        double[][] r0Diff = Rotation.procrustesAlgorithmForRotation(r0ZYX, r0UpZYX);
        
        double dq0 = Rotation.distanceBetweenQuaternions(q0Barfoot2, q0Up);
        double dq0Euclid = Rotation.distanceBetweenQuaternionEuclideanTransformations(q0Barfoot2, q0Up);
        
        // r0Q is ZYX
        double[][] r0Q = Rotation.createRotationMatrixFromQuaternion4(q0Barfoot2);
        double[][] r0QUpZYX = Rotation.applyRotationPerturbationZYX(theta0, dTheta0);
       
        /*
        System.out.printf("theta0ExZYX=\n%s\n", FormatArray.toString(theta0ExZYX, "%.3e"));
        System.out.printf("theta0UpExZYX=\n%s\n", FormatArray.toString(theta0UpExZYX, "%.3e"));
        System.out.printf("theta0UpExZYXMinusTheta0=\n%s\n", FormatArray.toString(theta0UpExZYXMinusTheta0, "%.3e"));
        System.out.printf("r0ZYX=\n%s\n", FormatArray.toString(r0ZYX, "%.3e"));
        System.out.printf("r0UpZYX=\n%s\n", FormatArray.toString(r0UpZYX, "%.3e"));
        System.out.printf("r0FromTheta0UpZYX=\n%s\n", FormatArray.toString(r0FromTheta0UpZYX, "%.3e"));
        System.out.printf("theta0=\n%s\n", FormatArray.toString(theta0, "%.3e"));
        System.out.printf("r0Diff=\n%s\n", FormatArray.toString(r0Diff, "%.3e"));
        */
        /*
        double[] axis0Rod = Rotation.extractRodriguesRotationAxis(r0);
        double[][] r0Rod = Rotation.createRodriguesFormulaRotationMatrix(theta0);
        System.out.printf("axis0Rod=\n%s\n", FormatArray.toString(axis0Rod, "%.3e"));
        System.out.printf("r0Rod=\n%s\n", FormatArray.toString(r0Rod, "%.3e"));
                
        double[][] r0XYZ = Rotation.createRotationXYZ(theta0);
        
        double[][] sTheta = Rotation.sTheta(theta0);
        double[][] sThetaInv = MatrixUtil.pseudoinverseRankDeficient(sTheta);
        double[] dRotVector = MatrixUtil.multiplyMatrixByColumnVector(sTheta, dTheta0);
        double[] dThetaWithPotentialSingularity = 
            MatrixUtil.multiplyMatrixByColumnVector(sThetaInv, dRotVector);
        
        System.out.printf("sTheta=\n%s\n", FormatArray.toString(sTheta, "%.3e"));
        System.out.printf("inv sTheta=\n%s\n", FormatArray.toString(sThetaInv, "%.3e"));
        
        System.out.printf("theta=\n   %s\n", FormatArray.toString(theta0, "%.3e"));
        System.out.printf("delta theta=\n   %s\n\n", FormatArray.toString(dTheta0, "%.3e"));
        System.out.printf("*theta extracted from r0 rotated=\n   %s\n", FormatArray.toString(theta0ExZYX, "%.3e"));
        System.out.printf("*theta extracted from r rotated=\n   %s\n\n", FormatArray.toString(theta0UpExZYX, "%.3e"));
        System.out.printf("delta rotation vector=\n   %s\n", FormatArray.toString(dRotVector, "%.3e"));
        System.out.printf("dThetaWithPotentialSingularity=\n   %s\n\n", 
            FormatArray.toString(dThetaWithPotentialSingularity, "%.3e"));
        
        System.out.printf("theta updated by + dT =\n   %s\n", FormatArray.toString(thetaUpdated, "%.3e"));
        System.out.printf("thetaExtrMinusThetaT =\n   %s\n", FormatArray.toString(thetaExtrMinusTheta, "%.3e"));
        
        
        System.out.printf("r0=\n%s\n", FormatArray.toString(r0ZYX, "%.3e"));
        System.out.printf("r rotated=\n%s\n", FormatArray.toString(r0UpZYX, "%.3e"));
        System.out.printf("r calculated from theta updates=\n%s\n", FormatArray.toString(rUpdated, "%.3e"));
        
        
        System.out.printf("rDiff_rodrigues=\n%s\n", FormatArray.toString(rDiff_rodrigues, "%.3e"));
        System.out.printf("rDiff_barfoot=\n%s\n", FormatArray.toString(rDiff_barfoot, "%.3e"));
        */
        
    }
}
