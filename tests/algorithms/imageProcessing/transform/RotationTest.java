
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

        //PAUSED HERE.  need passive, active argument:

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

        double[] p = new double[]{1, 0, 0, 0};
        double[] p2 = Rotation.rotateVectorByQuaternion4(qActive, p);
        // for an active rotation, we expect result = j
        // for a passive rotation, we expect result = k;
        double[] ePassive = new double[]{0, 1, 0, 0};
        double[] eActive = new double[]{0, 0, 1, 0};
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

        //PAUSED HERE

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

        // the operators
        //inverseQuaternionBarfoot
        // multiplyQuaternionsBarfoot
        //quaternionConjugateOperator
        //quaternionLefthandCompoundOperator
        //quaternionRighthandCompoundOperator
        //rotateAPointByQuaternionBarfoot
        //rotateVectorByQuaternion4
        //rotationBetweenTwoDirections0
        //rotationBetweenTwoDirections1
        //procrustesAlgorithmForRotation

    }


    public void estPerturbations() throws NotConvergedException {

        double tol = 1E-5;

        boolean passive = false;

        // applying sequential perturbations to rotation matrix
        double[][] rot0 = MatrixUtil.createIdentityMatrix(3);
        double[] euler = Rotation.extractThetaFromZYX(rot0, passive);
        assertTrue(MiscMath.areEqual(new double[]{0,0,0}, euler, tol));

        double[][] r1;
        double[] perturb, theta0;
        boolean returnQuaternion;

        // perturbations <= about 0.2 radians ~11 degrees
        perturb = new double[]{0.1, -0.05, 0.5};
        theta0 = new double[]{0, 0, 0};
        returnQuaternion = false;

        Rotation.RotationPerturbationMatrix rP = (Rotation.RotationPerturbationMatrix)
                Rotation.applySingularitySafeRotationPerturbation(
                        theta0, perturb, Rotation.EulerSequence.ZYX_ACTIVE, returnQuaternion);
        r1 = rP.rotation;

        // TODO: consider perspective and angles
        euler = Rotation.extractThetaFromZYX(r1, passive);

        //assertTrue(MiscMath.areEqual(perturb, euler, tol));

        int t = 1;
    }


    public void testSequences() throws NotConvergedException {

        //from active perspective now:

        double[] eulerXYZ = new double[]{-0.5, 0.12, 1.32};
        double[] eulerZYX = new double[]{1.32, 0.12, -0.5};

        double[][] rActiveXYZ = Rotation.createRotationFromEulerAngles(eulerXYZ, Rotation.EulerSequence.XYZ_ACTIVE);

        double[] rotVecXYZ = Rotation.createRotationVectorFromEulerAngles(eulerXYZ, Rotation.EulerSequence.XYZ_ACTIVE);
        //[0.5032547612511571, 0.23015806481528148, -1.3190236188603013]
        // same as R.from_euler('XYZ', [+0.5, -0.12, -1.32]).as_rotvec()
        double[] e = new double[]{0.50325476,  0.23015806, -1.31902362};
        assertTrue(MiscMath.areEqual(e, rotVecXYZ, 1E-5));

        double[] rotVecZYX = Rotation.createRotationVectorFromEulerAngles(eulerZYX, Rotation.EulerSequence.ZYX_ACTIVE);
        //[-1.2592811528236525, -0.4274666442528654, 0.3454896074433141]
        // same as R.from_euler('xyz', [-1.32, -0.12, 0.5]).as_rotvec()
        e = new double[]{-1.25928115, -0.42746664,  0.34548961};
        assertTrue(MiscMath.areEqual(e, rotVecZYX, 1E-5));

        /*
        Barfoot is using active and intrinsic.

        scipy library uses intrinsic and allows specification of active or passive, but the
        input euler or rotation must match the sequence order.

        Barfoot "1-2-3" and "alpha-beta-gamma" as names, refer to the passive Euler rotation sequence XYZ,
        then Barfoot impl. converts the sequence to active by reversing the order and making the angles negative.
        active A(a)*B(b)*C(c) == passive C(-c)*B(-b)*A(-a)

        Barfoot examples use 2 different euler sequences, active XYZ and active ZYX.
         */

        //double[] qBXYZActive = Rotation.createQuaternionBarfootFromRotationVector(rotVecXYZ,
        //        Rotation.EulerSequence.XYZ_ACTIVE);
        //[-0.2636069507418558, 0.06366132010060735, 0.5669163288230646, 0.7778589126296669]
        /* this is similar to -1 * vector portion of R.from_rotvec(rotVecXYZ).as_quat()
         scalar portions are similar but not same
         here: array([-0.26360695,  0.06366132,  0.56691633,  0.77785891])
         scipy: array([ 0.23071752,  0.10551613, -0.60470735,  0.75496013])
         where scipy: R.from_rotvec([0.50325476, 0.23015806, -1.31902362]).as_quat()
         */

        //double[] qBZYXActive = Rotation.createQuaternionBarfootFromRotationVector(rotVecZYX,
        //        Rotation.EulerSequence.ZYX_ACTIVE);
        //[0.5374232236710791, 0.0699706093309076, -0.25880741233766974, 0.7995618273829274]
        // this is similar to -1 *vector portion of R.from_rotvec(rotVecZYX).as_quat()
        // scalar portions are similar but not same

        //createQuaternionUnitLengthBarfootFromEuler(euler, seq)


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
        System.out.printf("AA r=\n%s\n", FormatArray.toString(r, "%.5f"));
        System.out.printf("RV r=\n%s\n", FormatArray.toString(r2, "%.5f"));
        System.out.printf("e=\n%s\n", FormatArray.toString(e, "%.5f"));

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
        System.out.printf("AA r=\n%s\n", FormatArray.toString(r, "%.5f"));
        System.out.printf("RV r=\n%s\n", FormatArray.toString(r2, "%.5f"));
        System.out.printf("e=\n%s\n", FormatArray.toString(e, "%.5f"));

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

    public void _testProcrustesAlgorithm() throws NotConvergedException {
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

    public void _testRodriguesRotationBouguet() throws NotConvergedException {
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

        boolean passive = true;
        int i;
        om = new double[3];
        dom = new double[3];
        for (i = 0; i < 3; ++i) {
            om[i] = rand.nextDouble();
            dom[i] = rand.nextDouble()/1000000;
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
        // norm1 ~ 1E-6;  norm2 ~ 4E-7
        double norm1 = MatrixUtil.spectralNorm(MatrixUtil.pointwiseSubtract(R2, R1));
        double norm2 = MatrixUtil.spectralNorm(MatrixUtil.pointwiseSubtract(R2, R2a));
        double gain = norm1/norm2;
        System.out.printf("norm1=%.4e, norm2=%.4e, gain=%.4e\n", norm1, norm2, gain);

        double[][] Rcompare = Rotation.createRotationRodriguesFormula(om, passive);
        double[] omCompare = Rotation.extractRotationVectorRodrigues(R1);
        double norm = MatrixUtil.spectralNorm(MatrixUtil.pointwiseSubtract(R1, Rcompare));
        System.out.printf("bouguet R=\n%s\n", FormatArray.toString(R1, "%.4e"));
        System.out.printf("existing R=\n%s\n", FormatArray.toString(Rcompare, "%.4e"));

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

        omCompare = Rotation.extractRotationVectorRodrigues(R);
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
        double[][] _R = Rotation.createRotationRodriguesFormula(om, false);
        Rotation.RodriguesRotation rRot7 = Rotation.createRotationRodriguesBouguet(om, passive);
        R = rRot7.r;
        System.out.printf("R(om) = \n%s\n", FormatArray.toString(R, "%.4e"));
        System.out.printf("R(om) existing = \n%s\n", FormatArray.toString(_R, "%.4e"));
        dR = rRot7.dRdR;
        Rotation.RodriguesRotation rRot8 = Rotation.extractRotationVectorRodriguesBouguet(R);
        om2 = rRot8.rotVec;
        double[] _om2 = Rotation.extractRotationVectorRodrigues(R);
        // om and om2 should be  the same
        System.out.printf("rand(3)*pi = om = %s\n", FormatArray.toString(om, "%.4e"));
        System.out.printf("om2=toVec(R(om)) = %s\n", FormatArray.toString(om2, "%.4e"));
        System.out.printf("_om2=toVec(R(om)) existing = %s\n", FormatArray.toString(_om2, "%.4e"));

        System.out.printf("R(om2) = \n%s\n",
                FormatArray.toString(Rotation.createRotationRodriguesBouguet(om2, passive).r, "%.4e"));
        System.out.printf("R(_om2) existing = \n%s\n",
                FormatArray.toString(Rotation.createRotationRodriguesFormula(_om2, false), "%.4e"));

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
        Rotation.RodriguesRotation rRot9 = Rotation.createRotationRodriguesBouguet(om, passive);
        R = rRot9.r;
        dR = rRot9.dRdR;
        Rotation.RodriguesRotation rRot10 = Rotation.extractRotationVectorRodriguesBouguet(R);
        om2 = rRot10.rotVec;
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

        Rotation.RodriguesRotation rRot11 = Rotation.createRotationRodriguesBouguet(om, passive);
        R = rRot11.r;
        Rotation.RodriguesRotation rRot12 =
                Rotation.createRotationRodriguesBouguet(
                        Rotation.extractRotationVectorRodriguesBouguet(R).rotVec, passive);
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
        Rotation.RodriguesRotation rRot13 = Rotation.extractRotationVectorRodriguesBouguet(R);
        om = rRot13.rotVec;
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
        Rotation.RodriguesRotation rRot14 = Rotation.extractRotationVectorRodriguesBouguet(R);
        om = rRot14.rotVec;
        norm = MatrixUtil.lPSum(om, 2);
        System.out.printf("norm=%.4e\n\n", norm);
    }

    // impl distance between from barfoot
}
