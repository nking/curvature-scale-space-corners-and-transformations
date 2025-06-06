
package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath;
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

    public void testCreateRotationXYZ() throws NotConvergedException {

        double[] eulerXYZ = new double[]{-0.5, 0.12, 1.32};
        double[][] r = Rotation.createRotationXYZ(eulerXYZ[0], eulerXYZ[1], eulerXYZ[2], true);

        double c1 = Math.cos(eulerXYZ[0]);
        double c2 = Math.cos(eulerXYZ[1]);
        double c3 = Math.cos(eulerXYZ[2]);
        double s1 = Math.sin(eulerXYZ[0]);
        double s2 = Math.sin(eulerXYZ[1]);
        double s3 = Math.sin(eulerXYZ[2]);

        //from "Euler Angles, Quaternions, and Transformation Matrices", NASA Shuttle Program, July 1977
        double[][] e = new double[][]{
                {c2*c3, -c2*s3, s2},
                {s1*s2*c3 + c1*s3, -s1*s2*s3 +c1*c3, -s1*c2},
                {-c1*s2*c3 +s1*s3, c1*s2*s3 + s1*c3, c1*c2}
        };

        //System.out.printf("XYZ r=\n%s\n", FormatArray.toString(r, "%.5f"));
        //System.out.printf("e=\n%s\n", FormatArray.toString(e, "%.5f"));

        assertTrue(MiscMath.areEqual(e, r, 1E-5));
        double[] _eulerXYZ = Rotation.extractThetaFromXYZ(r, true);
        assertTrue(MiscMath.areEqual(eulerXYZ, _eulerXYZ, 1E-5));

        double[][] r2 = Rotation.createRotationXYZ(eulerXYZ[0], eulerXYZ[1], eulerXYZ[2], false);
        double[] eulerXYZ2 = Rotation.extractThetaFromXYZ(r2, false);
        assertTrue(MiscMath.areEqual(eulerXYZ, eulerXYZ2, 1E-5));

        double[] axis = new double[3];
        double angle = Rotation.createAngleAxisFromEulerAnglesXYZ(eulerXYZ, axis);
        double[] rotVec = Rotation.createRotationVectorFromAngleAxis(axis, angle);
        double[] eRV = new double[]{-0.34548961,  0.42746664,  1.25928115};
        assertTrue(MiscMath.areEqual(eRV, rotVec, 1E-5));

        double[] qB1 = Rotation.createQuaternionUnitLengthBarfootFromRotationVector(rotVec);
        double[] qB2 = Rotation.createQuaternionUnitLengthBarfoot(angle, axis);
        double[] qB3 = Rotation.createQuaternionUnitLengthBarfootFromEulerXYZ(eulerXYZ);

        double[] eQ = new double[]{-0.1594735 ,  0.19731303,  0.58126776,  0.77315171};
        assertTrue(MiscMath.areEqual(eQ, qB1, 1E-5));
        assertTrue(MiscMath.areEqual(eQ, qB2, 1E-5));
        assertTrue(MiscMath.areEqual(eQ, qB3, 1E-5));
    }

    public void testAngleAxisAndRotationVector() throws NotConvergedException {
        double[] axis = new double[]{1 / Math.sqrt(3), 1 / Math.sqrt(3), 1 / Math.sqrt(3)};//~33 degrees
        double angle = 2 * Math.PI / 3; // 120 degrees // same as 60 degrees from -1*axis

        double[] rotVec = Rotation.createRotationVectorFromAngleAxis(axis, angle);
        //[1.2091995761561454, 1.2091995761561454, 1.2091995761561454] // 69.3 degrees
        double[] _axis = new double[3];
        double _angle = Rotation.createAngleAxisFromRotationVector(rotVec, _axis);
        //_angle=, _axis=

        // since 120 degrees CCW is > |120-180| CW, will reverse the axis and sign of angle.
        // with angle2 = angle - Math.PI // -60 degrees
        // we scale the axis by (angle/angle2)
        double angle2 = angle - Math.PI;
        double[] axis2 = Arrays.copyOf(axis, axis.length);
        MatrixUtil.multiply(axis2, (angle / angle2));
        double[] rotVec2 = Rotation.createRotationVectorFromAngleAxis(axis2, angle2);
        //[1.2091995761561454, 1.2091995761561454, 1.2091995761561454]
        // same rotation vector, though the angle axis representations were different.
        assertTrue(MiscMath.areEqual(rotVec, rotVec2, 1E-5));
        double[] _axis2 = new double[3];
        double _angle2 = Rotation.createAngleAxisFromRotationVector(rotVec2, _axis2);
        // now, both rotation vectors are normalized to same standard:
        assertTrue(MiscMath.areEqual(_axis, _axis2, 1E-5));
        assertTrue(Math.abs(_angle - _angle2) < 1E-5);

        // considering inverse formula: (angle, axis) inverse = (-angle, axis) or (angle, -axis)
        double[] negRotVec = Arrays.copyOf(rotVec, rotVec.length);
        MatrixUtil.multiply(negRotVec, -1);
        double[] axis3 = Arrays.copyOf(axis, axis.length);
        double angle3 = -1 * angle;
        double[] rotVec3 = Rotation.createRotationVectorFromAngleAxis(axis3, angle3);
        assertTrue(MiscMath.areEqual(negRotVec, rotVec3, 1E-5));
        double[] axis4 = Arrays.copyOf(axis, axis.length);
        MatrixUtil.multiply(axis4, -1);
        double angle4 = angle;
        double[] rotVec4 = Rotation.createRotationVectorFromAngleAxis(axis4, angle4);
        assertTrue(MiscMath.areEqual(negRotVec, rotVec4, 1E-5));

        //R.from_rotvec([1.2091995761561454, 1.2091995761561454, 1.2091995761561454]).as_quat()
        //array([0.5, 0.5, 0.5, 0.5])
        double[] eq = new double[]{0.5, 0.5, 0.5, 0.5};
        // angle = 2.0943951023931953, axis=[0.5773502691896257, 0.5773502691896257, 0.5773502691896257]
        angle = Rotation.makeUnitLengthAngleAxis(angle, axis);
        // angle=  axis=
        double[] axis5 = Arrays.copyOf(axis, axis.length);
        //angle5=1.0471975511965974,  axis5 = [-0.5773502691896257, -0.5773502691896257, -0.5773502691896257]
        double[] qPassive = Rotation.createQuaternionUnitLengthBarfoot(angle, axis);
        assertTrue(MiscMath.areEqual(eq, qPassive, 1E-5));
        double[] qActive = Rotation.createQuaternionUnitLengthBarfoot(angle4, axis4);

        double[] p = new double[]{1, 0, 0};
        double[] p2 = Rotation.rotateVectorByQuaternion4(qActive, p);
        // for an active rotation, we expect result = j
        // for a passive rotation, we expect result = k;
        double[] ePassive = new double[]{0, 1, 0};
        double[] eActive = new double[]{0, 0, 1};
        assertTrue(MiscMath.areEqual(eActive, p2, 1E-5));

        assertTrue(MiscMath.areEqual(ePassive, Rotation.rotateVectorByQuaternion4(qPassive, p), 1E-5));

        // create rotation matrix for passive transformation and extract rotation vector
        boolean passive = true;
        double[] eulerXYZ = new double[]{0.12, -0.34, 1.34};
        double[][] _r1 = Rotation.createRotationXYZ(eulerXYZ[0], eulerXYZ[1], eulerXYZ[2], passive);
        Rotation.RodriguesRotation rr1 = Rotation.extractRotationVectorRodriguesBouguet(_r1);
        double[] rv1 = Rotation.extractRotationVectorRodrigues(_r1);
        //R.from_euler('XYZ', [0.12, -0.34, 1.34]).as_rotvec()
        //-0.12663695, -0.36569661,  1.30424048
        assertTrue(MiscMath.areEqual(new double[]{-0.12663695, -0.36569661,  1.30424048}, rv1, 1E-5));
        assertTrue(MiscMath.areEqual(new double[]{-0.12663695, -0.36569661,  1.30424048}, rr1.rotVec, 1E-5));

        // create rotation matrix for active transformation and extract rotation vector
        passive = false;
        double[][] _r2 = Rotation.createRotationXYZ(eulerXYZ[0], eulerXYZ[1], eulerXYZ[2], passive);
        Rotation.RodriguesRotation rr2 = Rotation.extractRotationVectorRodriguesBouguet(_r2);
        double[] rv2 = Rotation.extractRotationVectorRodrigues(_r2);
        assertTrue(MiscMath.areEqual(new double[]{-0.32857351,  0.20790901, -1.34495176}, rv2, 1E-5));
        assertTrue(MiscMath.areEqual(new double[]{-0.32857351,  0.20790901, -1.34495176}, rr2.rotVec, 1E-5));
    }

    public void testQuaternions() throws NotConvergedException {

        // testing for an active system
        boolean passive = false;

        // for comparison data to test Barfoot quaternions (which are active and intrinsic),
        // will use scipy Rotation for expected values
        //
        // start with euler angles as active (scipy is passive so will get -1*eulerAnglesPassive) and create an intrinsic Rotation.
        // extract the rotation vector and compare to this one.
        // extract the quaternion and compare to the Barfoot quaternion

        double[] eulerAnglesPassive = new double[]{0.12, -0.51, 1.45};
        //double[] eulerAnglesActive = Arrays.copyOf(eulerAnglesPassive, 3);
        //MatrixUtil.multiply(eulerAnglesActive, -1);

        double[][] r = Rotation.createRotationXYZ(eulerAnglesPassive[0], eulerAnglesPassive[1], eulerAnglesPassive[2], true);
        double[][] eRPassive = new double[][]{
                {0.10516813, -0.86638481, -0.48817725},
                {0.97853176,  0.17765111, -0.10447817},
                {0.17724353, -0.46670916,  0.86646828}
        };
        double[][] rActive = Rotation.createRotationXYZ(eulerAnglesPassive[0], eulerAnglesPassive[1], eulerAnglesPassive[2], false);
        double[][] eRActive = new double[][]{
                {0.10516813,  0.86638481,  0.48817725},
                {-0.99261631,  0.06162127,  0.10447817},
                {0.0604362 , -0.49556047,  0.86646828}
        };
        assertTrue(MiscMath.areEqual(eRPassive, r, 1E-5));
        assertTrue(MiscMath.areEqual(eRActive, rActive, 1E-5));

        double[] rotVec = Rotation.extractRotationVectorRodrigues(r);
        //[-0.27172186, -0.49915489,  1.38393504]
        double[] rotVecActive = Rotation.extractRotationVectorRodrigues(rActive);
        //[-0.46634452,  0.33243642, -1.4447986]

        double[] axisPassive = new double[3];
        double[] axisActive = new double[3];
        double anglePassive = Rotation.createAngleAxisFromRotationVector(rotVec, axisPassive);
        double angleActive = Rotation.createAngleAxisFromRotationVector(rotVecActive, axisActive);

        double[] qPassive = Rotation.createQuaternionHamiltonFromAngleAxis(anglePassive, axisPassive);
        double[] e = new double[]{0.73302243, -0.12354021, -0.22694421,  0.62921559};
        assertTrue(MiscMath.areEqual(e, qPassive, 1E-5));
        double[] qB = Rotation.createQuaternionBarfootFromHamilton(qPassive);
        Rotation.convertQuaternionHamiltonToBarfoot(qPassive);
        e = new double[]{-0.12354021, -0.22694421,  0.62921559,  0.73302243};
        assertTrue(MiscMath.areEqual(e, qB, 1E-5));
        assertTrue(MiscMath.areEqual(e, qPassive, 1E-5));

        // Barfoot paper uses Active and intrinsic transformations, so checking consistency of those here
        qB = Rotation.createQuaternionUnitLengthBarfoot(angleActive, axisActive);
        e = new double[]{-0.21040352,  0.14998738, -0.65185867,  0.71296173};
        assertTrue(MiscMath.areEqual(e, qB, 1E-5));

        qB = Rotation.createQuaternionUnitLengthBarfootFromRotationVector(rotVecActive);
        assertTrue(MiscMath.areEqual(e, qB, 1E-5));

        qB = Rotation.createQuaternionUnitLengthBarfootFromEuler(eulerAnglesPassive, Rotation.EulerSequence.XYZ_ACTIVE);
        assertTrue(MiscMath.areEqual(e, qB, 1E-5));

        double[][] rB4 = Rotation.createRotation4FromQuaternion(qB);
        double[][] rB3 = MatrixUtil.copySubMatrix(rB4, 0, 2, 0, 2);
        assertTrue(MiscMath.areEqual(eRActive, rB3, 1E-5));


        double[] qBInv = Rotation.inverseQuaternionBarfoot(qB);
        for (int i = 0; i < 3; ++i) {
            assertEquals(qB[i], -qBInv[i]);
        }
        assertEquals(qB[3], qBInv[3]);

        double[] ax1 =  new double[]{1./Math.sqrt(1+4+9), -2./Math.sqrt(1+4+9), 3./Math.sqrt(1+4+9)};
        double[] ax2 =  new double[]{10./Math.sqrt(100+1+16), -1./Math.sqrt(100+1+16), 4./Math.sqrt(100+1+16)};
        double ang1 = 0.5;
        double ang2 = 1.3;
        double[] q1 = Rotation.createQuaternionUnitLengthBarfoot(ang1, ax1);
        //[0.06612148940441465, -0.1322429788088293, 0.19836446821324394, 0.9689124217106447]
        double[] q2 = Rotation.createQuaternionHamiltonFromAngleAxis(ang2, ax2);
        //[0.7960837985490559, 0.5594950300243704, -0.05594950300243705, 0.2237980120097482]
        double[] q12 = Rotation.quaternionMultiply(q1, q2, false);
        //[0.8897183441754278, 0.3508917648939988, -0.15208773240925755, 0.24929011014863733]
        // agrees with numpy quaternions when account for intrinsic and hamilton
        e = new double[]{0.8897183441754278, 0.3508917648939988, -0.15208773240925755, 0.24929011014863733};
        assertTrue(MiscMath.areEqual(e, q12, 1E-5));

        q12 = Rotation.quaternionMultiply(q1, q2, true);
        //[0.8897183441754278, 0.3508917648939988, -0.15208773240925755, 0.24929011014863733]
        // agrees with numpy quaternions when account for intrinsic and hamilton
        e = new double[]{0.682548333857785, 0.67412017253864, 0.132454542793051, 0.249290110148637};
        assertTrue(MiscMath.areEqual(e, q12, 1E-5));

        double[] qDiv1 = Rotation._quaternionDivide(q1, q2, false);
        e = new double[]{-0.860122628415764, -0.410083196413325, 0.240874879688247, 0.184391237432149};
        assertTrue(MiscMath.areEqual(e, qDiv1, 1E-5));

        double[] qDiv21Passive = Rotation._quaternionDivide(q2, q1, true);
        e = new double[]{-0.860122628415764, -0.410083196413325, 0.240874879688247, 0.184391237432149};
        assertTrue(MiscMath.areEqual(e, qDiv1, 1E-5));

        double[] pt = new double[]{120, -15, 39};
        double[] pt2 = Rotation.rotateAPointByQuaternionBarfoot(q1, pt);
        double[] v2 = Rotation.rotateVectorByQuaternion4(q1, pt);
        e = new double[]{103.41623131,  23.29723934,  70.05941579};
        assertTrue(MiscMath.areEqual(e, pt2, 1E-5));
        assertTrue(MiscMath.areEqual(e, v2, 1E-5));

        double q1Mag = MatrixUtil.lPSum(q1, 2);
        double q2Mag = MatrixUtil.lPSum(q2, 2);

        // for fraction 0, should return q1
        double[] qS1 = Rotation.quaternionSlerp(q1, q2, 0., true);
        assertTrue(MiscMath.areEqual(q1, qS1, 1E-5));

        // for fraction 1, should return q2
        double[] qS2 = Rotation.quaternionSlerp(q1, q2, 1., true);
        assertTrue(MiscMath.areEqual(q2, qS2, 1E-5));

        double[] qS3 = Rotation.quaternionSlerp(q1, q2, 0.25, true);
        // expected result is from: https://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/index.htm
        // and python module quaternion
        e = new double[]{0.332932567367378,0.07726708691910845, 0.15462983810434694, 0.9269710437942627};
        assertTrue(MiscMath.areEqual(e, qS3, 1E-5));

        double[] qS4 = Rotation.quaternionSlerp(q1, q2, 0.97, true);
        e = new double[]{0.7919862213113104, 0.5490604092641526, -0.04707957064428259, 0.2628193414870057};
        assertTrue(MiscMath.areEqual(e, qS4, 1E-5));
    }

    public void testGeodesicNorm() throws NotConvergedException {
        // testing log using pyquaternion test code.
        // https://github.com/KieranWynn/pyquaternion/blob/master/pyquaternion/test/test_quaternion.py
        double[] axis = new double[]{1, 0, 0};
        double angle = Math.PI;
        double[] q = Rotation.createQuaternionUnitLengthBarfoot(angle, axis);
        double[] logQ = Rotation._log(q);
        assertTrue(MiscMath.areEqual(
                new double[]{Math.PI/2., 0, 0, 0}, logQ, 1E-5));

        q = new double[]{1,0,0, 0};
        double[] p = new double[]{0,1,0, 0};
        double d = Rotation.quaternionGeodesicNorm(q, p);
        assertTrue(Math.abs(d - Math.PI/2) < 1E-5);

        //q = Quaternion(angle=pi/2, axis=[1,0,0])
        //p = Quaternion(angle=pi/2, axis=[0,1,0])
        q = new double[]{Math.sqrt(2)/2.,0,0, Math.sqrt(2)/2.};
        p = new double[]{0, Math.sqrt(2)/2., 0, Math.sqrt(2)/2.};
        d = Rotation.quaternionGeodesicNorm(q, p);
        assertTrue(Math.abs(d - Math.PI/3) < 1E-5);

        q = new double[]{-0.5, -0.5, -0.5,  -0.5};
        p = new double[]{0.5, 0.5, 0.5,  0.5};
        d = Rotation.quaternionGeodesicNorm(q, p);
        assertTrue(Math.abs(d) < 1E-5);

        //double[] d1 = new double[]{0, 45.*(Math.PI/180), 0};
        //double[] d2 = new double[]{0, 72.*(Math.PI/180), 0};
        //double[] qd12P = Rotation.rotationBetweenTwoDirections0(d1, d2, true);
        //double[] qd12A = Rotation.rotationBetweenTwoDirections0(d1, d2, false);

        double[][] r1 = Rotation.createRotationXYZ(0.1, 0.2, 0.3);
        double[][] r2 = Rotation.createRotationXYZ(0.1, 0.2, 0.3);
        double[][] r2ToR1 = Rotation.procrustesAlgorithmForRotation(r1, r2);
        double[][] expected = MatrixUtil.createIdentityMatrix(3);
        double[][] diff = MatrixUtil.pointwiseSubtract(r2ToR1, expected);
        double diffSN = MatrixUtil.frobeniusNorm(diff);
        assertTrue(Math.abs(diffSN) < 1E-7);

        r1 = Rotation.createRotationRoll(36 * Math.PI/180);
        r2 = Rotation.createRotationRoll(-17 * Math.PI/180);
        expected = Rotation.createRotationRoll(53 * Math.PI/180.);
        r2ToR1 = Rotation.procrustesAlgorithmForRotation(r1, r2);
        diff = MatrixUtil.pointwiseSubtract(r2ToR1, expected);
        diffSN = MatrixUtil.frobeniusNorm(diff);
        assertTrue(Math.abs(diffSN) < 1E-7);

        //rotationBetweenTwoDirections0
        //rotationBetweenTwoDirections1
    }

    public void testPerturbations() throws NotConvergedException {

        double tol = 1E-5;

        boolean passive = false;

        // applying sequential perturbations to rotation matrix
        double[][] rot0 = MatrixUtil.createIdentityMatrix(3);
        double[] euler = Rotation.extractThetaFromZYX(rot0, passive);
        assertTrue(MiscMath.areEqual(new double[]{0,0,0}, euler, tol));

        double[][] r1;
        double[] perturb, theta0;
        boolean returnQuaternion;

        // small perturbations must be <= about 0.1 radians ~5.7 degrees
        // 0.5 is large, but wanting to look at how close a rotation from perturb as rotation vector is to result
        perturb = new double[]{0.1, -0.05, 0.5};
        theta0 = new double[]{0, 0, 0};
        returnQuaternion = false;

        Rotation.RotationPerturbationMatrix rP = (Rotation.RotationPerturbationMatrix)
                Rotation.applySingularitySafeRotationPerturbation(
                        theta0, perturb, Rotation.EulerSequence.ZYX_ACTIVE, returnQuaternion);
        r1 = rP.rotation;
        // for intrinsic, extracting eulerZYX=[0.09966865249116204, -0.049710870978323454, 0.4636476090008061]
        // perturb = 0.1, -0.05, 0.5
        double[] eulerR1ZYX = Rotation.extractThetaFromZYX(r1, passive);
        for (int i = 0; i < eulerR1ZYX.length; ++i) {
            double tolerance = Math.abs(0.8*perturb[i]);
            double diff = Math.abs(eulerR1ZYX[i] - perturb[i]);
            assertTrue(diff <= tolerance);
        }

        Rotation.RotationPerturbationQuaternion qP = (Rotation.RotationPerturbationQuaternion)
                Rotation.applySingularitySafeRotationPerturbation(
                        theta0, perturb, Rotation.EulerSequence.ZYX_ACTIVE, !returnQuaternion);
        double[] q1 = qP.quaternion;
        double[] q0 = Rotation.createIdentityQuaternion();
        double[] diffQ1Q0 = Rotation.quaternionMultiply(q1, Rotation.quaternionConjugateOperator(q0), passive);

        // for intrinsic, q1 = [-0.05, 0.025, -0.25, 1.0]
        // expecting intrinsic ZYX, active: [-0.04223358,  0.0365512 , -0.24580705,  0.96770824]
        double[] eQ = new double[]{-0.04223358,  0.0365512 , -0.24580705,  0.96770824};
        for (int i = 0; i < eQ.length; ++i) {
            double tolerance = 0.5*Math.abs(eQ[i]);
            double diff = Math.abs(q1[i] - eQ[i]);
            assertTrue(diff <= tolerance);
        }

        double[] perturb2 = new double[]{0.1, -0.08, 0.05};
        Rotation.RotationPerturbationMatrix rP2 = (Rotation.RotationPerturbationMatrix)
                Rotation.applySingularitySafeRotationPerturbation(rP, perturb2);
        double[][] r2 = rP2.rotation;
        double[][] diffR1ToR2 = Rotation.procrustesAlgorithmForRotation(r2, r1);
        double[] extractedP2 = Rotation.extractThetaFromZYX(diffR1ToR2, passive);
        double[][] expectedDiffR1ToR2 = Rotation.createRotationRodriguesFormula(perturb2, passive);

        // extractedP2 is on the order of perturb2 but smaller.  except last term tends towards value 1 because of the
        //    internal structure of the algorithm (col2 = i2)
        // expectedDiffR1ToR2 is on the order of diffR1ToR2 excepting [1][0] and [0][1]
        // roughly asserting here for future regression checks
        for (int i = 0; i < 2; ++i) {
            double tolerance = Math.abs(0.85*perturb2[i]);
            double diff = Math.abs(extractedP2[i] - perturb2[i]);
            assertTrue(diff <= tolerance);
        }
        for (int i = 0; i < diffR1ToR2.length; ++i) {
            for (int j = 0; j < diffR1ToR2[i].length; ++j) {
                double tolerance = 0.7 * Math.abs(expectedDiffR1ToR2[i][j]);
                if ((i==1 && j==0) || (i==0 && j==1)) continue;
                double diff = Math.abs(diffR1ToR2[i][j] - expectedDiffR1ToR2[i][j]);
                assertTrue(diff <= tolerance);
            }
        }

        Rotation.RotationPerturbationQuaternion qP2 = (Rotation.RotationPerturbationQuaternion)
                Rotation.applySingularitySafeRotationPerturbation(qP, perturb2);

        // perturb2 = 0.1, -0.08, 0.05
        double[] q2 = qP2.quaternion;

        //diff in rotation
        //q2 = q(dphi)^+ * q1
        // q(dphi)^+ = q2 * q1^-1
        double[] diffQ2Q1 = Rotation.quaternionMultiply(q2, Rotation.quaternionConjugateOperator(q1), passive);

        //R.from_euler('ZYX', [-0.05, 0.08, -0.1]).as_quat()
        double[] eQ2 = new double[]{-0.04892521,  0.04117523, -0.02294818,  0.99768948};
        for (int i = 0; i < 4; ++i) {
            if (i == 2) continue;  // skip the 3rd term. details in createQuaternionBarfootFromEuler
            double tolerance = 0.5*Math.abs(eQ2[i]);
            double diff = Math.abs(diffQ2Q1[i] - eQ2[i]);
            assertTrue(diff <= tolerance);
        }

    }

    public void testCreateRotationZYX() {
        double[] eulerXYZ = new double[]{-0.5, 0.12, 1.32};
        double[][] r = Rotation.createRotationZYX(eulerXYZ[0], eulerXYZ[1], eulerXYZ[2]);

        double c1 = Math.cos(eulerXYZ[2]);
        double c2 = Math.cos(eulerXYZ[1]);
        double c3 = Math.cos(eulerXYZ[0]);
        double s1 = Math.sin(eulerXYZ[2]);
        double s2 = Math.sin(eulerXYZ[1]);
        double s3 = Math.sin(eulerXYZ[0]);

        double[][] e = new double[][]{
                {c1*c2, c1*s2*s3 - s1*c3, c1*s2*c3 + s1*s3},
                {s1*c2, s1*s2*s3 + c1*c3, s1*s2*c3 - c1*s3},
                {-s2, c2*s3, c2*c3}
        };

        //System.out.printf("ZYX r=\n%s\n", FormatArray.toString(r, "%.5f"));
        //System.out.printf("e=\n%s\n", FormatArray.toString(e, "%.5f"));

        assertTrue(MiscMath.areEqual(e, r, 1E-5));

        double[] eulerZYX2 = Rotation.extractThetaFromZYX(r, true);
        assertTrue(MiscMath.areEqual(eulerZYX2, eulerXYZ, 1E-5));

        r = Rotation.createRotationZYX(eulerXYZ[0], eulerXYZ[1], eulerXYZ[2], false);
        double[] euler2 = Rotation.extractThetaFromZYX(r, false);
        assertTrue(MiscMath.areEqual(eulerXYZ, euler2, 1E-5));

    }

    public void testCreateRotationFromUnitLengthAngleAxis() throws NotConvergedException {

        //[0, 2*pi], [0, pi], and [0, 2*pi]
        double[] axis = new double[]{2.14, 0.13, -0.5};
        double angle = 4;
        angle = Rotation.makeUnitLengthAngleAxis(angle, axis);
        // 0.8584, [0.9720747555844925, 0.05905127019905796, -0.22712026999637674]

        // unit length has no effect on rotVec
        double[] rotVec = Rotation.createRotationVectorFromAngleAxis(axis, angle);
        //[8.56, 0.52, -2.0] if not unit length
        //[[2.452274178231329, 0.14896992671498727, -0.5729612565961049]  which is factor 3.49063742 smaller

        boolean passive = true;

        double[][] r = Rotation.createRotationRodriguesFormula(axis, angle, passive);

        Rotation.RodriguesRotation rr = Rotation.createRotationRodriguesBouguet(rotVec, passive);
        double[][] r2 = rr.r;

        double[][] e = new double[][] {
                {0.9000724 ,  0.23591439, -0.3663524},
                {-0.02759772, -0.80820664, -0.58825199},
                {-0.43486555,  0.53957987, -0.72093378}
        };
        //System.out.printf("AA r=\n%s\n", FormatArray.toString(r, "%.5f"));
        //System.out.printf("RV r=\n%s\n", FormatArray.toString(r2, "%.5f"));
        //System.out.printf("e=\n%s\n", FormatArray.toString(e, "%.5f"));

        assertTrue(MiscMath.areEqual(e, r, 1E-5));
        assertTrue(MiscMath.areEqual(e, r2, 1E-5));

        // correct:
        double[] _rv1 = Rotation.extractRotationVectorRodrigues(r2);
        double[] _rv2 = Rotation.extractRotationVectorRodriguesBouguet(r2).rotVec;
        /* from scipy Rotation
        R.from_matrix(np.array([[0.9000724 ,  0.23591439, -0.3663524], [-0.02759772, -0.80820664, -0.58825199], [-0.43486555,  0.53957987, -0.72093378]])).as_rotvec()
        array([ 2.45227417,  0.14896993, -0.57296126])
         */
        double[] e2 = new double[]{2.45227417,  0.14896993, -0.57296126};
        assertTrue(MiscMath.areEqual(e2, _rv1, 1E-5));
        assertTrue(MiscMath.areEqual(e2, _rv2, 1E-5));

        // =========================================
        // try all larger angles
        axis = new double[]{-1*(Math.PI/2 - 0.12), 0.13, -0.5};
        //[-1.4507963267948965, 0.13, -0.5]
        angle = Rotation.makeUnitLengthAngleAxis(angle, axis);

        rotVec = Rotation.createRotationVectorFromAngleAxis(axis, angle);
        //[-0.7004065598563302, 0.06276060333187992, -0.24138693589184584]
        r = Rotation.createRotationRodriguesFormula(axis, angle, passive);

        rr = Rotation.createRotationRodriguesBouguet(rotVec, passive);
        r2 = rr.r;

        e = new double[][] {
                {0.97030335,  0.19876976,  0.13784773},
                {-0.24073974,  0.73799186,  0.63040653},
                {0.02357526, -0.644871  ,  0.76392775}
        };
        //System.out.printf("AA r=\n%s\n", FormatArray.toString(r, "%.5f"));
        //System.out.printf("RV r=\n%s\n", FormatArray.toString(r2, "%.5f"));
        //System.out.printf("e=\n%s\n", FormatArray.toString(e, "%.5f"));

        assertTrue(MiscMath.areEqual(e, r, 1E-5));
        assertTrue(MiscMath.areEqual(e, r2, 1E-5));
    }

    public void _testRodriguesFormula() throws NotConvergedException {
        
        //from http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example.html

        boolean passive = false;
        double[] axis = new double[]{-1.451113, -1.827059, -0.179105};
        // in degrees: -83.1426505 , -104.68276962,  -10.26196059
        double[][] r = Rotation.createRotationRodriguesFormula(axis, passive);
        
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
        double[] rotVec2 = Rotation.extractRotationVectorRodrigues(r);
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

    public void _testProcrustesAlgorithm2() throws NotConvergedException {
        
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
        double[][] k2ExtrRot = Rotation.createRotationRodriguesFormula(
                new double[]{0.00611, 0.00409, -0.00359}, false);
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

    public void testRodriguesRotationBouguet() throws NotConvergedException {
        // tests from the Bouguet toolbox in rodrigues.m
        //System.out.println("\ntestRodriguesRotationBouguet");

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

        boolean passive = true;
        int i;
        om = new double[3];
        dom = new double[3];
        for (i = 0; i < 3; ++i) {
            om[i] = rand.nextDouble();
            dom[i] = rand.nextDouble()*1E-6;
        }
        Rotation.RodriguesRotation rRot = Rotation.createRotationRodriguesBouguet(om, passive);
        // [3X3]
        R1 = rRot.r;
        // [9X3]
        dR1 = rRot.dRdR;
        Rotation.RodriguesRotation rRot2 = Rotation.createRotationRodriguesBouguet(MatrixUtil.add(om, dom), passive);
        R2 = rRot2.r;
        //R2a = R1 + reshape(dR1 * dom,3,3);
        // [9X3][3X1] = [9X1]
        double[] tmp = MatrixUtil.multiplyMatrixByColumnVector(dR1, dom);
        double[][] tmp2 = new double[3][];
        tmp2[0] = new double[]{tmp[0], tmp[3], tmp[6]};
        tmp2[1] = new double[]{tmp[1], tmp[4], tmp[7]};
        tmp2[2] = new double[]{tmp[2], tmp[5], tmp[8]};
        R2a = MatrixUtil.pointwiseAdd(R1, tmp2);
        //gain = norm(R2 - R1)/norm(R2 - R2a)
        // norm1 ~ 1E-6;  norm2 ~ 4E-7
        double norm1 = MatrixUtil.spectralNorm(MatrixUtil.pointwiseSubtract(R2, R1));
        double norm2 = MatrixUtil.spectralNorm(MatrixUtil.pointwiseSubtract(R2, R2a));
        double gain = norm1/norm2;
        //System.out.printf("norm1=%.4e, norm2=%.4e, gain=%.4e\n", norm1, norm2, gain);

        double[][] Rcompare = Rotation.createRotationRodriguesFormula(om, passive);
        double[] omCompare = Rotation.extractRotationVectorRodrigues(R1);
        double norm = MatrixUtil.spectralNorm(MatrixUtil.pointwiseSubtract(R1, Rcompare));
        // d21 is similar to norm1, though maybe only for small angles
        // d2a is similar to norm2, ...

        // =========
        //System.out.println("\nTEST OF dOmdR");
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
            dom[i] = rand.nextDouble()*1E-4;
        }
        Rotation.RodriguesRotation rRot3 = Rotation.createRotationRodriguesBouguet(om, passive);
        R = rRot3.r;

        double[] omPlusDom = MatrixUtil.add(om, dom);
        Rotation.RodriguesRotation rRot4 = Rotation.createRotationRodriguesBouguet(omPlusDom, passive);
        dR = MatrixUtil.pointwiseSubtract(rRot4.r, R);

        double[] checkOMC = Rotation.extractRotationVectorRodrigues(R);

        Rotation.RodriguesRotation rRot5 = Rotation.extractRotationVectorRodriguesBouguet(R);
        omc = rRot5.rotVec; // om should equal omc
        domdR = rRot5.dRdR;
        Rotation.RodriguesRotation rRot6 = Rotation.extractRotationVectorRodriguesBouguet(
                MatrixUtil.pointwiseAdd(R, dR));
        om2 = rRot6.rotVec;
        om_app = MatrixUtil.add(omc,
                MatrixUtil.multiplyMatrixByColumnVector(domdR, MatrixUtil.stack(dR)));
        norm1 = MatrixUtil.lPSum(MatrixUtil.subtract(om2, omc), 2);
        norm2 = MatrixUtil.lPSum(MatrixUtil.subtract(om2, om_app), 2);
        gain = norm1/norm2;
        //System.out.printf("norm1=%.7e, norm2=%.7e, gain=%.7e\n", norm1, norm2, gain);

        double[] tmp3 = Arrays.copyOf(om, 3);
        MatrixUtil.multiply(tmp3, 1./MatrixUtil.lPSum(om, 2));
        //System.out.printf("om=%s norm=%.7e\n", FormatArray.toString(tmp3, "%.7e"),
        //        MatrixUtil.lPSum(om, 2));

        tmp3 = Arrays.copyOf(omc, 3);
        MatrixUtil.multiply(tmp3, 1./MatrixUtil.lPSum(omc, 2));
        //System.out.printf("omc=%s norm=%.7e\n", FormatArray.toString(tmp3, "%.7e"),
        //        MatrixUtil.lPSum(omc, 2));

        tmp3 = Arrays.copyOf(om2, 3);
        MatrixUtil.multiply(tmp3, 1./MatrixUtil.lPSum(om2, 2));
        //System.out.printf("om2=%s norm=%.7e\n", FormatArray.toString(tmp3, "%.7e"),
        //        MatrixUtil.lPSum(om2, 2));

        tmp3 = Arrays.copyOf(omPlusDom, 3);
        MatrixUtil.multiply(tmp3, 1./MatrixUtil.lPSum(omPlusDom, 2));
        //System.out.printf("omPlusDom=%s norm=%.7e\n", FormatArray.toString(tmp3, "%.7e"),
        //        MatrixUtil.lPSum(omPlusDom, 2));

        tmp3 = Arrays.copyOf(om_app, 3);
        MatrixUtil.multiply(tmp3, 1./MatrixUtil.lPSum(om_app, 2));
        //System.out.printf("omapp=%s norm=%.7e\n", FormatArray.toString(tmp3, "%.7e"),
        //        MatrixUtil.lPSum(om_app, 2));

        // om ~ omPlusDom as the later is a small addition
        // omc ~ omapp as the later is a small addition
        //R
        //om+dom
        //dR = om+dom - R
        //System.out.printf("R(om) =\n%s\n", FormatArray.toString(R, "%.4e"));
        //System.out.printf("R(om+dom)=\n%s\n", FormatArray.toString(rRot4.r, "%.4e"));
        //System.out.printf("dR = R(om+dom)-R=\n%s\n", FormatArray.toString(dR, "%.4e"));

        omCompare = Rotation.extractRotationVectorRodrigues(R);
        norm = MatrixUtil.lPSum(MatrixUtil.subtract(omc, omCompare), 2);
        //System.out.printf("norm diffs between bouguet(omc) and existing(omc)=%.4e\n", norm);
        assertTrue(Math.abs(norm) < 1E-7);

        //==========
        //System.out.println("\nOTHER BUG: (FIXED NOW!!!)");
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
        double[][] _R = Rotation.createRotationRodriguesFormula(om, false);
        Rotation.RodriguesRotation rRot7 = Rotation.createRotationRodriguesBouguet(om, passive);
        R = rRot7.r;
        //System.out.printf("R(om) = \n%s\n", FormatArray.toString(R, "%.4e"));
        //System.out.printf("R(om) existing = \n%s\n", FormatArray.toString(_R, "%.4e"));
        assertTrue(MiscMath.areEqual(R, _R, 1E-5));
        dR = rRot7.dRdR;
        Rotation.RodriguesRotation rRot8 = Rotation.extractRotationVectorRodriguesBouguet(R);
        om2 = rRot8.rotVec;
        // rRot8.rotVec should == rRot7.rotVec
        double[] _om2 = Rotation.extractRotationVectorRodrigues(R);
        // om and om2 should be  the same
        assertTrue(MiscMath.areEqual(rRot8.rotVec, rRot7.rotVec, 1E-5));
        assertTrue(MiscMath.areEqual(rRot8.r, rRot7.r, 1E-5));

        //=======
        /*
        %% NORMAL OPERATION
        om = randn(3,1);
        [R,dR]= rodrigues(om);
        [om2] = rodrigues(R);
        [om om2]
        */
        //System.out.println("\nNORMAL OPERATION");
        for (i = 0; i < 3; ++i) {
            om[i] = rand.nextDouble();
        }
        Rotation.RodriguesRotation rRot9 = Rotation.createRotationRodriguesBouguet(om, passive);
        R = rRot9.r;
        dR = rRot9.dRdR;
        Rotation.RodriguesRotation rRot10 = Rotation.extractRotationVectorRodriguesBouguet(R);
        om2 = rRot10.rotVec;
        assertTrue(MiscMath.areEqual(om, om2, 1E-5));

        //=====================
        /*
        %% Test: norm(om) = pi
        u = randn(3,1);
        u = u / sqrt(sum(u.^2));
        om = pi*u;
        R = rodrigues(om);
        R2 = rodrigues(rodrigues(R));
        norm(R - R2)*/
        //System.out.println("\nTest: norm(om) = pi");
        double[] u = new double[3];
        for (i = 0; i < 3; ++i) {
            u[i] = rand.nextDouble();
        }
        MatrixUtil.multiply(u, 1./MatrixUtil.lPSum(u, 2));
        om = Arrays.copyOf(u, u.length);
        MatrixUtil.multiply(om, Math.PI);

        Rotation.RodriguesRotation rRot11 = Rotation.createRotationRodriguesBouguet(om, passive);
        R = rRot11.r;
        Rotation.RodriguesRotation rRot12 =
                Rotation.createRotationRodriguesBouguet(
                        Rotation.extractRotationVectorRodriguesBouguet(R).rotVec, passive);
        R2 = rRot12.r;
        norm1 = MatrixUtil.spectralNorm(R);
        norm2 = MatrixUtil.spectralNorm(R2);
        norm = MatrixUtil.spectralNorm(MatrixUtil.pointwiseSubtract(R, R2));
        //System.out.printf("norm(R)=%.4e, norm(R2)=%.4e, norm(R-R2)=%.4e\n", norm1, norm2, norm);
        assertTrue(Math.abs(norm1 - 1) < 1E-11);
        assertTrue(Math.abs(norm2 - 1) < 1E-11);
        assertTrue(Math.abs(norm2 - norm1) < 1E-11);

        //================
        /*
        %% Another test case where norm(om)=pi from Chen Feng (June 27th, 2014)
        R = [-0.950146567583153 -6.41765854280073e-05 0.311803617668748; ...
             -6.41765854277654e-05 -0.999999917385145 -0.000401386434914383; ...
              0.311803617668748 -0.000401386434914345 0.950146484968298];
        om = rodrigues(R)
        norm(om) - pi
        */
        //System.out.println("\nAnother test case where norm(om)=pi from Chen Feng (June 27th, 2014)");
        // math.acos(1-0.950146567583153)*radToDeg = 87.142 degrees
        R[0] = new double[]{-0.950146567583153, -6.41765854280073e-05, 0.311803617668748};
        R[1] = new double[]{-6.41765854277654e-05, -0.999999917385145, -0.000401386434914383};
        R[2] = new double[]{0.311803617668748, -0.000401386434914345, 0.950146484968298};
        Rotation.RodriguesRotation rRot13 = Rotation.extractRotationVectorRodriguesBouguet(R);
        om = rRot13.rotVec;
        //System.out.printf("om=%s\nom2=%s\n", FormatArray.toString(om, "%.4e"),
        //        FormatArray.toString(om2, "%.4e"));
        norm = MatrixUtil.lPSum(om, 2);
        //System.out.printf("norm=%.4e\n", norm);
        assertTrue(Math.abs(norm - Math.PI) < 1E-7);

        //================
        /*
        %% Another test case where norm(om)=pi from 余成义 (July 1st, 2014)
        R = [-0.999920129411407	-6.68593208347372e-05	-0.0126384464118876; ...
             9.53007036072085e-05	-0.999997464662094	-0.00224979713751896; ...
            -0.0126382639492467	-0.00225082189773293	0.999917600647740];
        om = rodrigues(R)
        norm(om) - pi
         */
        //System.out.println("\nAnother test case where norm(om)=pi from 余成义 (July 1st, 2014)");
        R[0] = new double[]{-0.999920129411407,	-6.68593208347372e-05,	-0.0126384464118876};
        R[1] = new double[]{9.53007036072085e-05,	-0.999997464662094,	-0.00224979713751896};
        R[2] = new double[]{-0.0126382639492467,	-0.00225082189773293,	0.999917600647740};
        Rotation.RodriguesRotation rRot14 = Rotation.extractRotationVectorRodriguesBouguet(R);
        om = rRot14.rotVec;
        norm = MatrixUtil.lPSum(om, 2);
        //System.out.printf("norm=%.4e\n\n", norm);
        assertTrue(Math.abs(norm - Math.PI) < 1E-4);
    }

    // impl distance between from barfoot
}
