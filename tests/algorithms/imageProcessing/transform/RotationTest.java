
package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import java.util.Arrays;
import java.util.Random;

import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

/**
 *
 * @author nichole
 */
public class RotationTest extends TestCase {
    
    public RotationTest() {
    }

    public void testRodriguesFormula() {
        
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
        double tol = 1e-3;
        double diff;
        int i, j;
        for (i = 0; i < r.length; ++i) {
            for (j = 0; j < r[i].length; ++j) {
                diff = Math.abs(expected[i][j] - r[i][j]);
                assertTrue(diff < tol);
            }
        }
        
        double[] axis1 = MatrixUtil.normalizeL2(axis);
        double[] rotVec2 = Rotation.extractRodriguesRotationVector(r);
        double angle2 = 0;
        for (double a : rotVec2) {
            angle2 += (a*a);
        }
        angle2 = Math.sqrt(angle2);
        double[] axis2 = MatrixUtil.normalizeL2(rotVec2);
        assertEquals(axis.length, axis2.length);
        System.out.printf("\naxis1=%s\n", FormatArray.toString(axis1, "%.3e"));
        System.out.printf("axis2=%s\n", FormatArray.toString(axis2, "%.3e"));
        System.out.printf("angle2=%.3e\n", angle2);
        
        for (i = 0; i < axis1.length; ++i) {
            diff = Math.abs(axis1[i] - axis2[i]);
            assertTrue(diff < tol);
        }
       
    }

    public void testProcrustesAlgorithm() throws NotConvergedException {
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

    public void testProcrustesAlgorithm2() throws NotConvergedException {
        
        System.out.println("\ntestProcrustesAlgorithm2()");
        
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
        double[][] orthogRot = Rotation.procrustesAlgorithmForRotation(
                MatrixUtil.transpose(x1), MatrixUtil.transpose(x2));
        System.out.printf("rotation from Procrustes algorithm\n=%s\n", FormatArray.toString(orthogRot, "%.3e"));
    }
    
    public void test2() {
        
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
        System.out.printf("theta0=%s\n", FormatArray.toString(theta0, "%.3e"));
        System.out.printf("dTheta0=%s\n", FormatArray.toString(dTheta0, "%.3e"));
        System.out.printf("theta0Up=%s\n", FormatArray.toString(theta0Up, "%.3e"));
        
        //System.out.printf("original=%s\n", FormatArray.toString(theta, "%.3e"));
        
        double[][] r0ZYX = Rotation.createRotationZYX(theta0); // checked          
        double[] theta0ExZYX = Rotation.extractThetaFromZYX(r0ZYX); // ** recovered exactly, that is, same as theta0 
        
        double[][] r0UpZYX = Rotation.applySingularitySafeEulerAnglesPerturbationZYX(theta0, dTheta0);
        double[] theta0UpExZYX = Rotation.extractThetaFromZYX(r0UpZYX); // ** nearly the same as theta0Up
        double[][] r0FromTheta0UpZYX = Rotation.createRotationZYX(theta0Up);
        double[] theta0UpExZYXMinusTheta0 = new double[3];
        for (i = 0; i < 3; ++i) {
            theta0UpExZYXMinusTheta0[i] = theta0UpExZYX[i] - theta0[i];
        }
        System.out.printf("\nr0ZYX=\n%s\n", FormatArray.toString(r0ZYX, "%.3e"));
        System.out.printf("r0UpZYX=\n%s\n", FormatArray.toString(r0UpZYX, "%.3e"));
        System.out.printf("r0FromTheta0UpZYX=\n%s\n", FormatArray.toString(r0FromTheta0UpZYX, "%.3e"));
        System.out.printf("\ntheta0ExZYX=%s\n", FormatArray.toString(theta0ExZYX, "%.3e"));
        System.out.printf("theta0UpExZYX=%s\n", FormatArray.toString(theta0UpExZYX, "%.3e"));
        System.out.printf("theta0UpExZYXMinusTheta0(=dTheta?)=%s\n", FormatArray.toString(theta0UpExZYXMinusTheta0, "%.3e"));
                
        
        double[] q0Barfoot = Rotation.createQuaternionZYXFromEuler(theta0);
        double[] q0Hamilton = Rotation.createHamiltonQuaternionZYX(theta0); // checked
        double[] q0Barfoot2 = Rotation.convertHamiltonToBarfootQuaternion(q0Hamilton); // checked 
        
        double[] a0ZYX1 = Rotation.extractRotationAxisFromZXY(r0ZYX);        
        double[] rotVecZYX2 = Rotation.extractRodriguesRotationVector(r0ZYX);
        double angleZYX2 = 0;
        for (double a : rotVecZYX2) {
            angleZYX2 += (a*a);
        }
        angleZYX2 = Math.sqrt(angleZYX2);
        double[] a0ZYX2 = MatrixUtil.normalizeL2(rotVecZYX2);
        
        
        double m0ZYX1 = MatrixUtil.lPSum(a0ZYX1, 2);
        double m0ZYX2 = MatrixUtil.lPSum(a0ZYX2, 2);
        double[][] r0FromA0ZYX1 = Rotation.createRodriguesFormulaRotationMatrix(a0ZYX1);
        
        double[] q0AABarfootZYX1 = Rotation.createUnitLengthQuaternionBarfoot(a0ZYX1, m0ZYX1);
        double[] q0AABarfootZYX2 = Rotation.createUnitLengthQuaternionBarfoot(a0ZYX2, m0ZYX2);
        
        System.out.printf("\na0ZYX1=%s\n", FormatArray.toString(a0ZYX1, "%.3e"));
        System.out.printf("a0ZYX2=%s\n", FormatArray.toString(a0ZYX2, "%.3e"));
        System.out.printf("m0ZYX1=%.3e\n", m0ZYX1);
        System.out.printf("m0ZYX2=%.3e\n", m0ZYX2);
        System.out.printf("\nr0FromA0ZYX1=\n%s\n", FormatArray.toString(r0FromA0ZYX1, "%.3e"));
        
        System.out.printf("\nq0Barfoot=%s\n", FormatArray.toString(q0Barfoot, "%.3e"));
        System.out.printf("q0Hamilton=%s\n", FormatArray.toString(q0Hamilton, "%.3e"));
        System.out.printf("q0Barfoot2=%s\n", FormatArray.toString(q0Barfoot2, "%.3e"));
        
        System.out.printf("q0AABarfootZYX1=%s\n", FormatArray.toString(q0AABarfootZYX1, "%.3e"));
        System.out.printf("q0AABarfootZYX2=%s\n", FormatArray.toString(q0AABarfootZYX2, "%.3e"));
        
        // rotate a quaternion.  
        double[] dq0Hamilton = Rotation.createHamiltonQuaternionZYX(dTheta0);
        double[] dq0Barfoot = Rotation.convertHamiltonToBarfootQuaternion(dq0Hamilton);
        double[] q0Up = Rotation.rotateVectorByQuaternion4(dq0Barfoot, q0Barfoot2);
        double[] q0SafeUp = Rotation.applySingularitySafeRotationPerturbationQuaternion(theta0, dTheta0);
        double[][] r0QSafeUp = Rotation.createRotationMatrixFromQuaternion4(q0SafeUp); 
        
        double dq0Dist = Rotation.distanceBetweenQuaternions(q0Barfoot2, q0Up);
        //double dq0Euclid = Rotation.distanceBetweenQuaternionEuclideanTransformations(q0Barfoot2, q0Up);
        
        System.out.printf("\ndq0Hamilton=%s\n", FormatArray.toString(dq0Hamilton, "%.3e"));
        System.out.printf("dq0Barfoot=%s\n", FormatArray.toString(dq0Barfoot, "%.3e"));
        System.out.printf("\nq0Up=%s\n", FormatArray.toString(q0Up, "%.3e"));
        System.out.printf("q0SafeUp=%s\n", FormatArray.toString(q0SafeUp, "%.3e"));
        System.out.printf("r0QSafeUp^T=\n%s\n", FormatArray.toString(MatrixUtil.transpose(r0QSafeUp), "%.3e")); // approx equal to r0UpZYX within 0.035 radians = 2 degrees
        
        double[][] r0Q = Rotation.createRotationMatrixFromQuaternion4(q0Barfoot2); // this is transposed compared to Rotation.createRotationZYX(theta0)       
        double[][] r0Diff = Rotation.procrustesAlgorithmForRotation(r0ZYX, r0UpZYX);
        double[] thetaExR0Diff = Rotation.extractThetaFromZYX(r0Diff);

        //dPhi= S(theta) * dTheta.  length is 3.
        double[] dPhi = Rotation.createRotationVector(theta0, dTheta0);
        //sTheta is the matrix relating angular velocity to rotation angle rates.
        double[][] sTheta = Rotation.sTheta(theta0);
        double[] rodVecFromEuler = Rotation.convertEulerAnglesToRodriguesVectorForZYX(theta0);
        double[] rodPertVecFromEuler = Rotation.convertEulerAnglesToRodriguesVectorForZYX(dTheta0);
        
        System.out.printf("\nsTheta=%s\n", FormatArray.toString(sTheta, "%.3e"));
        System.out.printf("dPhi=%s\n", FormatArray.toString(dPhi, "%.3e"));
        System.out.printf("rodVecFromEuler=%s\n", FormatArray.toString(rodVecFromEuler, "%.3e"));
        System.out.printf("rodPertVecFromEuler=%s\n", FormatArray.toString(rodPertVecFromEuler, "%.3e"));
            
        System.out.printf("\ndq0Dist=%.3e\n", dq0Dist);
        //System.out.printf("dq0Euclid=%.3e\n", dq0Euclid);
        System.out.printf("r0Q=\n%s\n", FormatArray.toString(r0Q, "%.3e"));
        System.out.printf("r0QUpZYX=\n%s\n", FormatArray.toString(r0UpZYX, "%.3e"));
        System.out.printf("r0Diff=\n%s\n", FormatArray.toString(r0Diff, "%.3e"));
        System.out.printf("thetaExR0Diff(theta extracted from r)Diff)=\n%s\n", FormatArray.toString(thetaExR0Diff, "%.3e"));
                
        // XYZ
        double[][] r0XYZ = Rotation.createRotationXYZ(theta0);
        double[] theta0ExXYZ = Rotation.extractThetaFromXYZ(r0XYZ);     
        double[][] r0UpXYZ = Rotation.applySingularitySafeEulerAnglesPerturbationZYX(theta0, dTheta0);
        double[] theta0UpExXYZ = Rotation.extractThetaFromXYZ(r0UpXYZ);
        double[][] r0FromTheta0UpXYZ = Rotation.createRotationXYZ(theta0Up);
        double[] theta0UpExXYZMinusTheta0 = new double[3];
        for (i = 0; i < 3; ++i) {
            theta0UpExXYZMinusTheta0[i] = theta0UpExXYZ[i] - theta0[i];
        }
        System.out.printf("\nr0XYZ=\n%s\n", FormatArray.toString(r0XYZ, "%.3e"));
        System.out.printf("r0UpXYZ=\n%s\n", FormatArray.toString(r0UpXYZ, "%.3e"));
        System.out.printf("r0FromTheta0UpXYZ=\n%s\n", FormatArray.toString(r0FromTheta0UpXYZ, "%.3e"));
        System.out.printf("\ntheta0ExXYZ=%s\n", FormatArray.toString(theta0ExXYZ, "%.3e"));
        System.out.printf("theta0UpExXYZ=%s\n", FormatArray.toString(theta0UpExXYZ, "%.3e"));
        System.out.printf("theta0UpExXYZMinusTheta0=%s\n", FormatArray.toString(theta0UpExXYZMinusTheta0, "%.3e"));
        
        
    }
    
    public void test4() throws NotConvergedException {
         
        System.out.println("test4 for gimbal lock singularity");
        
        /*
        a test for gimbal lock.
        theta=[0, 90 degrees (==pi/2), 0]
        
        For the asymmetric matrix
            C_i_j_k = R(ax,ay,az)=R_z(az)*R_y(ay)*R_x(ax), singularity when
                ay=+-(pi/2)
        */
        int i;
        double[] theta0 
            = new double[]{0., 0.5*Math.PI/180., 0.};
        // 0.1 radians = 5.73 degrees.  1 degree = 0.0174 redians
        double[] dTheta0 
            = new double[]{2*Math.PI/180., 3*Math.PI/180., 5.5*Math.PI/180.};
        double[] theta0Up = new double[3];
        for (i = 0; i < 3; ++i) {
            theta0Up[i] = theta0[i] + dTheta0[i];
        }
        System.out.printf("theta0=%s\n", FormatArray.toString(theta0, "%.3e"));
        System.out.printf("dTheta0=%s\n", FormatArray.toString(dTheta0, "%.3e"));
        System.out.printf("theta0Up=%s\n", FormatArray.toString(theta0Up, "%.3e"));
        
        //System.out.printf("original=%s\n", FormatArray.toString(theta, "%.3e"));
        
        double[][] r0ZYX = Rotation.createRotationZYX(theta0);          
        double[] theta0ExZYX = Rotation.extractThetaFromZYX(r0ZYX); //** recovered exactly, that is, same as theta0 
        
        double[][] r0UpZYX = Rotation.applySingularitySafeEulerAnglesPerturbationZYX(theta0, dTheta0);
        double[] theta0UpExZYX = Rotation.extractThetaFromZYX(r0UpZYX); // ** nearly the same as theta0Up
        double[][] r0FromTheta0UpZYX = Rotation.createRotationZYX(theta0Up);
        double[] theta0UpExZYXMinusTheta0 = new double[3];
        for (i = 0; i < 3; ++i) {
            theta0UpExZYXMinusTheta0[i] = theta0UpExZYX[i] - theta0[i];
        }
        System.out.printf("\nr0ZYX=\n%s\n", FormatArray.toString(r0ZYX, "%.3e"));
        System.out.printf("r0UpZYX=\n%s\n", FormatArray.toString(r0UpZYX, "%.3e"));
        System.out.printf("r0FromTheta0UpZYX=\n%s\n", FormatArray.toString(r0FromTheta0UpZYX, "%.3e"));
        System.out.printf("\ntheta0ExZYX=%s\n", FormatArray.toString(theta0ExZYX, "%.3e"));
        System.out.printf("theta0UpExZYX=%s\n", FormatArray.toString(theta0UpExZYX, "%.3e"));
        System.out.printf("theta0UpExZYXMinusTheta0(=dTheta?)=%s\n", FormatArray.toString(theta0UpExZYXMinusTheta0, "%.3e"));
        
        double[] q0Barfoot = Rotation.createQuaternionZYXFromEuler(theta0);
        double[] q0Hamilton = Rotation.createHamiltonQuaternionZYX(theta0); 
        double[] q0Barfoot2 = Rotation.convertHamiltonToBarfootQuaternion(q0Hamilton); 
        
        // rotate a quaternion.  
        double[] dq0Hamilton = Rotation.createHamiltonQuaternionZYX(dTheta0);
        double[] dq0Barfoot = Rotation.convertHamiltonToBarfootQuaternion(dq0Hamilton);
        double[] q0Up = Rotation.rotateVectorByQuaternion4(dq0Barfoot, q0Barfoot2);
        double[] q0SafeUp = Rotation.applySingularitySafeRotationPerturbationQuaternion(theta0, dTheta0);
        double[][] r0QSafeUp = Rotation.createRotationMatrixFromQuaternion4(q0SafeUp); 
        
        double dq0Dist = Rotation.distanceBetweenQuaternions(q0Barfoot2, q0Up);
        //double dq0Euclid = Rotation.distanceBetweenQuaternionEuclideanTransformations(q0Barfoot2, q0Up);
        
        System.out.printf("\ndq0Hamilton=%s\n", FormatArray.toString(dq0Hamilton, "%.3e"));
        System.out.printf("dq0Barfoot=%s\n", FormatArray.toString(dq0Barfoot, "%.3e"));
        System.out.printf("\nq0Up=%s\n", FormatArray.toString(q0Up, "%.3e"));
        System.out.printf("q0SafeUp=%s\n", FormatArray.toString(q0SafeUp, "%.3e"));
        System.out.printf("r0QSafeUp^T=\n%s\n", FormatArray.toString(MatrixUtil.transpose(r0QSafeUp), "%.3e")); // approx equal to r0UpZYX within 0.035 radians = 2 degrees
        
        double[][] r0Q = Rotation.createRotationMatrixFromQuaternion4(q0Barfoot2); // this is transposed compared to Rotation.createRotationZYX(theta0)       
        double[][] r0Diff = Rotation.procrustesAlgorithmForRotation(r0ZYX, r0UpZYX); //similar to value in test3() which is good
        double[] thetaExR0Diff = Rotation.extractThetaFromZYX(r0Diff);

        //dPhi= S(theta) * dTheta.  length is 3.
        double[] dPhi = Rotation.createRotationVector(theta0, dTheta0);
        //sTheta is the matrix relating angular velocity to rotation angle rates.
        double[][] sTheta = Rotation.sTheta(theta0);
        double[] rodVecFromEuler = Rotation.convertEulerAnglesToRodriguesVectorForZYX(theta0);
        double[] rodPertVecFromEuler = Rotation.convertEulerAnglesToRodriguesVectorForZYX(dTheta0);
        
        System.out.printf("\nsTheta=%s\n", FormatArray.toString(sTheta, "%.3e"));
        System.out.printf("dPhi=%s\n", FormatArray.toString(dPhi, "%.3e"));
        System.out.printf("rodVecFromEuler=%s\n", FormatArray.toString(rodVecFromEuler, "%.3e"));
        System.out.printf("rodPertVecFromEuler=%s\n", FormatArray.toString(rodPertVecFromEuler, "%.3e"));
            
        System.out.printf("\ndq0Dist=%.3e\n", dq0Dist);
        //System.out.printf("dq0Euclid=%.3e\n", dq0Euclid);
        System.out.printf("r0Q=\n%s\n", FormatArray.toString(r0Q, "%.3e"));
        System.out.printf("r0QUpZYX=\n%s\n", FormatArray.toString(r0UpZYX, "%.3e"));
        System.out.printf("r0Diff=\n%s\n", FormatArray.toString(r0Diff, "%.3e"));
        System.out.printf("thetaExR0Diff(theta extracted from r)Diff)=\n%s\n", FormatArray.toString(thetaExR0Diff, "%.3e"));
                
    }

    public void testRodriguesRotationBouguet() throws NotConvergedException {
        // tests from the Bouguet toolbox in rodrigues.m
        System.out.println("\ntestRodriguesRotationBouguet");

        long seed = System.nanoTime();
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);

        double[] om;
        double[] dom;
        double[][] R1, dR1, R2, R2a;

        /*
        %% TEST OF dRdom:
        om = randn(3,1);
        dom = randn(3,1)/1000000;
        [R1,dR1] = rodrigues(om);
        R2 = rodrigues(om+dom);
        R2a = R1 + reshape(dR1 * dom,3,3);
        gain = norm(R2 - R1)/norm(R2 - R2a)
         */

        int i;
        om = new double[3];
        dom = new double[3];
        for (i = 0; i < 3; ++i) {
            om[i] = rand.nextDouble();
            dom[i] = rand.nextDouble()/1000000;
        }
        Rotation.RodriguesRotation rRot = Rotation.createRodriguesRotationMatrixBouguet(om);
        // [3X3]
        R1 = rRot.r;
        // [9X3]
        dR1 = rRot.dRdR;
        Rotation.RodriguesRotation rRot2 = Rotation.createRodriguesRotationMatrixBouguet(MatrixUtil.add(om, dom));
        R2 = rRot2.r;
        //R2a = R1 + reshape(dR1 * dom,3,3);
        // [9X3][3X1] = [9X1]
        double[] tmp = MatrixUtil.multiplyMatrixByColumnVector(dR1, dom);
        double[][] tmp2 = new double[3][];
        tmp2[0] = new double[]{tmp[0], tmp[3], tmp[6]};
        tmp2[1] = new double[]{tmp[1], tmp[4], tmp[7]};
        tmp2[2] = new double[]{tmp[2], tmp[5], tmp[8]};
        R2a = MatrixUtil.pointwiseAdd(R1, tmp2);
        // norm1 ~ 1E-6;  norm2 ~ 4E-7
        double norm1 = MatrixUtil.spectralNorm(MatrixUtil.pointwiseSubtract(R2, R1));
        double norm2 = MatrixUtil.spectralNorm(MatrixUtil.pointwiseSubtract(R2, R2a));
        double gain = norm1/norm2;
        System.out.printf("norm1=%.4e, norm2=%.4e, gain=%.4e\n", norm1, norm2, gain);

        double[][] Rcompare = Rotation.createRodriguesFormulaRotationMatrix(om);
        double[] omCompare = Rotation.extractRodriguesRotationVector(R1);
        double dCompare = distanceBetween(R1, Rcompare);
        double norm = MatrixUtil.spectralNorm(MatrixUtil.pointwiseSubtract(R1, Rcompare));
        System.out.printf("bouguet R=\n%s\n", FormatArray.toString(R1, "%.4e"));
        System.out.printf("existing R=\n%s\n", FormatArray.toString(Rcompare, "%.4e"));
        System.out.printf("dist between=%.4e, norm diffs=%.4e\n", dCompare, norm);

        double d21 = distanceBetween(R2, R1);
        double d2a = distanceBetween(R2, R2a);
        // d21 is similar to norm1, though maybe only for small angles
        // d2a is similar to norm2, ...

        // =========
        System.out.println("\nTEST OF dOmdR");
        double[][] R, dR, domdR;
        double[] omc, om2, om_app;
        /*
        %% TEST OF dOmdR:
        om = randn(3,1);
        R = rodrigues(om);
        dom = randn(3,1)/10000;
        dR = rodrigues(om+dom) - R;

        //[omc,domdR] = rodrigues(R);
        //[om2] = rodrigues(R+dR);
        //om_app = omc + domdR*dR(:); // [3X1] + [3X9]*[9X1]
        //gain = norm(om2 - omc)/norm(om2 - om_app)
        */

        for (i = 0; i < 3; ++i) {
            om[i] = rand.nextDouble();
            dom[i] = rand.nextDouble()/10000;
        }
        Rotation.RodriguesRotation rRot3 = Rotation.createRodriguesRotationMatrixBouguet(om);
        R = rRot3.r;

        double[] omPlusDom = MatrixUtil.add(om, dom);
        Rotation.RodriguesRotation rRot4 = Rotation.createRodriguesRotationMatrixBouguet(omPlusDom);
        dR = MatrixUtil.pointwiseSubtract(rRot4.r, R);

        double[] checkOMC = Rotation.extractRodriguesRotationVector(R);

        Rotation.RodriguesRotation rRot5 = Rotation.extractRodriguesRotationVectorBouguet(R);
        omc = rRot5.om; // om should equal omc
        domdR = rRot5.dRdR;
        Rotation.RodriguesRotation rRot6 = Rotation.extractRodriguesRotationVectorBouguet(
                MatrixUtil.pointwiseAdd(R, dR)
        );
        om2 = rRot6.om;
        om_app = MatrixUtil.add(omc,
                MatrixUtil.multiplyMatrixByColumnVector(domdR, MatrixUtil.stack(dR)));
        norm1 = MatrixUtil.lPSum(MatrixUtil.subtract(om2, omc), 2);
        norm2 = MatrixUtil.lPSum(MatrixUtil.subtract(om2, om_app), 2);
        gain = norm1/norm2;
        System.out.printf("norm1=%.4e, norm2=%.4e, gain=%.4e\n", norm1, norm2, gain);

        double[] tmp3 = Arrays.copyOf(om, 3);
        MatrixUtil.multiply(tmp3, 1./MatrixUtil.lPSum(om, 2));
        System.out.printf("om=%s norm=%.4e\n", FormatArray.toString(tmp3, "%.3e"),
                MatrixUtil.lPSum(om, 2));

        tmp3 = Arrays.copyOf(omc, 3);
        MatrixUtil.multiply(tmp3, 1./MatrixUtil.lPSum(omc, 2));
        System.out.printf("omc=%s norm=%.4e\n", FormatArray.toString(tmp3, "%.3e"),
                MatrixUtil.lPSum(omc, 2));

        tmp3 = Arrays.copyOf(om2, 3);
        MatrixUtil.multiply(tmp3, 1./MatrixUtil.lPSum(om2, 2));
        System.out.printf("om2=%s norm=%.4e\n", FormatArray.toString(tmp3, "%.3e"),
                MatrixUtil.lPSum(om2, 2));

        tmp3 = Arrays.copyOf(omPlusDom, 3);
        MatrixUtil.multiply(tmp3, 1./MatrixUtil.lPSum(omPlusDom, 2));
        System.out.printf("omPlusDom=%s norm=%.4e\n", FormatArray.toString(tmp3, "%.3e"),
                MatrixUtil.lPSum(omPlusDom, 2));

        tmp3 = Arrays.copyOf(om_app, 3);
        MatrixUtil.multiply(tmp3, 1./MatrixUtil.lPSum(om_app, 2));
        System.out.printf("omapp=%s norm=%.4e\n", FormatArray.toString(tmp3, "%.3e"),
                MatrixUtil.lPSum(om_app, 2));

        // om ~ omPlusDom as the later is a small addition
        // omc ~ omapp as the later is a small addition
        //R
        //om+dom
        //dR = om+dom - R
        System.out.printf("R(om) =\n%s\n", FormatArray.toString(R, "%.4e"));
        System.out.printf("R(om+dom)=\n%s\n", FormatArray.toString(rRot4.r, "%.4e"));
        System.out.printf("dR = R(om+dom)-R=\n%s\n", FormatArray.toString(dR, "%.4e"));

        omCompare = Rotation.extractRodriguesRotationVector(R);
        norm = MatrixUtil.lPSum(MatrixUtil.subtract(omc, omCompare), 2);
        System.out.printf("norm diffs between bouguet(omc) and existing(omc)=%.4e\n", norm);

        //==========
        System.out.println("\nOTHER BUG: (FIXED NOW!!!)");
        /*
        %% OTHER BUG: (FIXED NOW!!!)
        omu = randn(3,1);
        omu = omu/norm(omu)
        om = pi*omu;
        [R,dR]= rodrigues(om);
        [om2] = rodrigues(R);
        [om om2]*/
        double[] omu = new double[3];
        for (i = 0; i < 3; ++i) {
            omu[i] = rand.nextDouble();
        }
        om = MatrixUtil.normalizeL2(omu);
        MatrixUtil.multiply(om, Math.PI);
        double[][] _R = Rotation.createRodriguesFormulaRotationMatrix(om);
        Rotation.RodriguesRotation rRot7 = Rotation.createRodriguesRotationMatrixBouguet(om);
        R = rRot7.r;
        System.out.printf("R(om) = \n%s\n", FormatArray.toString(R, "%.4e"));
        System.out.printf("R(om) existing = \n%s\n", FormatArray.toString(_R, "%.4e"));
        dR = rRot7.dRdR;
        Rotation.RodriguesRotation rRot8 = Rotation.extractRodriguesRotationVectorBouguet(R);
        om2 = rRot8.om;
        double[] _om2 = Rotation.extractRodriguesRotationVector(R);
        // om and om2 should be  the same
        System.out.printf("rand(3)*pi = om = %s\n", FormatArray.toString(om, "%.4e"));
        System.out.printf("om2=toVec(R(om)) = %s\n", FormatArray.toString(om2, "%.4e"));
        System.out.printf("_om2=toVec(R(om)) existing = %s\n", FormatArray.toString(_om2, "%.4e"));

        System.out.printf("R(om2) = \n%s\n",
                FormatArray.toString(Rotation.createRodriguesRotationMatrixBouguet(om2).r, "%.4e"));
        System.out.printf("R(_om2) existing = \n%s\n",
                FormatArray.toString(Rotation.createRodriguesFormulaRotationMatrix(_om2), "%.4e"));

        //=======
        /*
        %% NORMAL OPERATION
        om = randn(3,1);
        [R,dR]= rodrigues(om);
        [om2] = rodrigues(R);
        [om om2]
        */
        System.out.println("\nNORMAL OPERATION");
        for (i = 0; i < 3; ++i) {
            om[i] = rand.nextDouble();
        }
        Rotation.RodriguesRotation rRot9 = Rotation.createRodriguesRotationMatrixBouguet(om);
        R = rRot9.r;
        dR = rRot9.dRdR;
        Rotation.RodriguesRotation rRot10 = Rotation.extractRodriguesRotationVectorBouguet(R);
        om2 = rRot10.om;
        System.out.printf("should be equal:\n");
        System.out.printf("om=%s\nom2=%s\n", FormatArray.toString(om, "%.4e"),
                FormatArray.toString(om2, "%.4e"));

        norm = MatrixUtil.lPSum(om, 2);
        tmp3 = Arrays.copyOf(om, 3);
        MatrixUtil.multiply(tmp3, 1./norm);
        System.out.printf("om norm=%.4e\nom normalized=%s\n", norm,
                FormatArray.toString(tmp3, "%.4e"));

        norm = MatrixUtil.lPSum(om2, 2);
        tmp3 = Arrays.copyOf(om2, 3);
        MatrixUtil.multiply(tmp3, 1./norm);
        System.out.printf("om2 norm=%.4e\nom2 normalized=%s\n", norm,
                FormatArray.toString(tmp3, "%.4e"));

        //=====================
        /*
        %% Test: norm(om) = pi
        u = randn(3,1);
        u = u / sqrt(sum(u.^2));
        om = pi*u;
        R = rodrigues(om);
        R2 = rodrigues(rodrigues(R));
        norm(R - R2)*/
        System.out.println("\nTest: norm(om) = pi");
        double[] u = new double[3];
        for (i = 0; i < 3; ++i) {
            u[i] = rand.nextDouble();
        }
        MatrixUtil.multiply(u, 1./MatrixUtil.lPSum(u, 2));
        om = Arrays.copyOf(u, u.length);
        MatrixUtil.multiply(om, Math.PI);

        Rotation.RodriguesRotation rRot11 = Rotation.createRodriguesRotationMatrixBouguet(om);
        R = rRot11.r;
        Rotation.RodriguesRotation rRot12 =
                Rotation.createRodriguesRotationMatrixBouguet(Rotation.extractRodriguesRotationVectorBouguet(R).om);
        R2 = rRot12.r;
        norm1 = MatrixUtil.spectralNorm(R);
        norm2 = MatrixUtil.spectralNorm(R2);
        norm = MatrixUtil.spectralNorm(MatrixUtil.pointwiseSubtract(R, R2));
        System.out.printf("norm(R)=%.4e, norm(R2)=%.4e, norm(R-R2)=%.4e\n", norm1, norm2, norm);

        //================
        /*
        %% Another test case where norm(om)=pi from Chen Feng (June 27th, 2014)
        R = [-0.950146567583153 -6.41765854280073e-05 0.311803617668748; ...
             -6.41765854277654e-05 -0.999999917385145 -0.000401386434914383; ...
              0.311803617668748 -0.000401386434914345 0.950146484968298];
        om = rodrigues(R)
        norm(om) - pi
        */
        System.out.println("\nAnother test case where norm(om)=pi from Chen Feng (June 27th, 2014)");
        // math.acos(1-0.950146567583153)*radToDeg = 87.142 degrees
        R[0] = new double[]{-0.950146567583153, -6.41765854280073e-05, 0.311803617668748};
        R[1] = new double[]{-6.41765854277654e-05, -0.999999917385145, -0.000401386434914383};
        R[2] = new double[]{0.311803617668748, -0.000401386434914345, 0.950146484968298};
        Rotation.RodriguesRotation rRot13 = Rotation.extractRodriguesRotationVectorBouguet(R);
        om = rRot13.om;
        System.out.printf("om=%s\nom2=%s\n", FormatArray.toString(om, "%.4e"),
                FormatArray.toString(om2, "%.4e"));
        norm = MatrixUtil.lPSum(om, 2);
        System.out.printf("norm=%.4e\n", norm);

        //================
        /*
        %% Another test case where norm(om)=pi from 余成义 (July 1st, 2014)
        R = [-0.999920129411407	-6.68593208347372e-05	-0.0126384464118876; ...
             9.53007036072085e-05	-0.999997464662094	-0.00224979713751896; ...
            -0.0126382639492467	-0.00225082189773293	0.999917600647740];
        om = rodrigues(R)
        norm(om) - pi
         */
        System.out.println("\nAnother test case where norm(om)=pi from 余成义 (July 1st, 2014)");
        R[0] = new double[]{-0.999920129411407,	-6.68593208347372e-05,	-0.0126384464118876};
        R[1] = new double[]{9.53007036072085e-05,	-0.999997464662094,	-0.00224979713751896};
        R[2] = new double[]{-0.0126382639492467,	-0.00225082189773293,	0.999917600647740};
        Rotation.RodriguesRotation rRot14 = Rotation.extractRodriguesRotationVectorBouguet(R);
        om = rRot14.om;
        norm = MatrixUtil.lPSum(om, 2);
        System.out.printf("norm=%.4e\n\n", norm);
    }

    private double distanceBetween(double[][] r0, double[][] r1) {
        return Rotation.distanceBetweenQuaternions(
                Rotation.createQuaternionZYXFromEuler(
                        Rotation.extractThetaFromZYX(r0)
                ),
                Rotation.createQuaternionZYXFromEuler(
                        Rotation.extractThetaFromZYX(r1)
                )
        );
    }
}
