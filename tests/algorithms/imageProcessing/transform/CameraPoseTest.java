package algorithms.imageProcessing.transform;

import static algorithms.imageProcessing.transform.Rotation.extractThetaFromZYX;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath;
import algorithms.util.FormatArray;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

import java.io.IOException;
import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author nichole
 */
public class CameraPoseTest extends TestCase {
    
    public CameraPoseTest() {
    }
    
    public void estNCBook() throws IOException {
        //NeurochemistryBookData2.runCorresMaker();
    }

    public void testZhangExample() throws Exception {

        // ======== first checking calculatePoseFromP with Zhang data example ================

        // from chap 2 of Camera Calibration book, chapter by Zhang (add reference)
        //https://people.cs.rutgers.edu/~elgammal/classes/cs534/lectures/CameraCalibration-book-chapter.pdf
        double[][] p = new double[][]{
                {7.025659E-1, -2.861189E-2, -5.377696E-1, 6.24189E1},
                {2.077632E-1, 1.265804, 1.591456E-1, 1.075646E1},
                {4.634764E-4, -5.282382E-5, 4.255347E-4, 1}
        };

        boolean passive = true;//false;
        double angle = 47.7*Math.PI/180.;

        // wide angle: f <= 35 mm
        // telephoto: f > 35mm
        // K in terms of
        // focal lengths in units of mm in positions[0][0] and [1][1]
        // and the optical center in units of pixel coordinates in [0][2] and [1][2]
        // position [0][1] holds the pixel skew if any.
        //  there is skew when the camera sensor is not perpendicular to optical axis.
        double[][] K = new double[][]{
                {1380.12, 0, 246.52},
                {0, 2032.57, 243.68},
                {0, 0, 1}
        };

        // angle between 2 image axes is 90 degrees
        // aspect ratio of the pixels = alpha/beta = 0.679
        // rotation axis: (almost vertical)
        double[] axis = new double[]{-0.08573, -0.99438, 0.0621};
        double[] rotVec = Rotation.createRotationVectorFromAngleAxis(axis, angle);

        double[][] rot = Rotation.createRotationRodriguesFormula(rotVec, passive);
        double[] eulerAngles = Rotation.extractThetaFromXYZ(rot, passive);
        double[] eulerAnglesDegrees = Arrays.copyOf(eulerAngles, eulerAngles.length);
        MatrixUtil.multiply(eulerAnglesDegrees, 180./Math.PI);

        /*
        expected rotation:
        [0.6750834808562537, -0.018079003106068848, -0.7375197919525242]
        [0.07388502469134674, 0.9963303528698166, 0.043206840623514806]
        [0.7340322179586962, -0.08365989240369291, 0.6739419302889822]
         */

        // t in mm = 1.5m from camera in WCS
        double[] translation = new double[]{-211.28, -106.06, 1583.75};
        double tZ = translation[2];

        // camera projection matrix P is defined in
        //     z*x = P*X
        //         = K*[R|T]*X
        //  where z is the distance to the object from the camera
        // Zhang alters his example P to move z from the RHS to the LHS:
        //      (K/z) * [R | T]
        double[][] reproducingP =
            MatrixUtil.multiply(
                new double[][] {
                {K[0][0]/tZ, K[0][1]/tZ, K[0][2]/tZ},
                {K[1][0]/tZ, K[1][1]/tZ, K[1][2]/tZ},
                {K[2][0]/tZ, K[2][1]/tZ, K[2][2]/tZ}},

                 new double[][] {
                         {rot[0][0], rot[0][1], rot[0][2], translation[0]},
                         {rot[1][0], rot[1][1], rot[1][2], translation[1]},
                         {rot[2][0], rot[2][1], rot[2][2], translation[2]}
                 }
        );

        Camera.CameraPoseParameters c = CameraPose.calculatePoseFromP(p);
        double fs = MatrixUtil.frobeniusNorm(
                MatrixUtil.pointwiseSubtract(K, c.getIntrinsicParameters().getIntrinsic()));

        // should be identity matrix if same:
        double[][] diffR = Rotation.procrustesAlgorithmForRotation(rot, c.getExtrinsicParameters().getRotation());
        double fsR = MatrixUtil.frobeniusNorm( MatrixUtil.pointwiseSubtract(
                diffR, MatrixUtil.createIdentityMatrix(3)));

        double[] diffT = MatrixUtil.subtract(translation, c.getExtrinsicParameters().getTranslation());

        assertTrue(fs < 1);
        assertTrue(fsR < 1);

        // taking a look at P composed from K instead of K/z gives even better matches to his example params (excepting p)
        double[][] pNotDivByZ = new double[3][4];
        for (int ii = 0; ii < 3; ++ii) {
            System.arraycopy(rot[ii], 0, pNotDivByZ[ii], 0, 3);
            pNotDivByZ[ii][3] = translation[ii];
        }
        pNotDivByZ = MatrixUtil.multiply(K, pNotDivByZ);

        Camera.CameraPoseParameters _c = CameraPose.calculatePoseFromP(pNotDivByZ);
        double _fs = MatrixUtil.frobeniusNorm(
                MatrixUtil.pointwiseSubtract(K, _c.getIntrinsicParameters().getIntrinsic()));

        double[] _r = Rotation.extractRotationAxisFromZXY(_c.getExtrinsicParameters().getRotation());
        // should be identity matrix if same:
        double[][] _diffR = Rotation.procrustesAlgorithmForRotation(rot, _c.getExtrinsicParameters().getRotation());
        double _fsR = MatrixUtil.frobeniusNorm( MatrixUtil.pointwiseSubtract(
                _diffR, MatrixUtil.createIdentityMatrix(3)));

        double[] _diffT = MatrixUtil.subtract(translation, _c.getExtrinsicParameters().getTranslation());

        assertTrue(_fs < 1);
        assertTrue(_fsR < 1);
        assertTrue(MatrixUtil.lPSum(_diffT, 2) < 1);

        // generate a set of features for an image
        // create a camera intrinsic matrix K,
        // motion w.r.t. image and scene, rotation R, and translation t
        // and then projected WCS features.
        // test that can recover the projection matrix

        /*
        TODO: write tests for other K matrices:
        smartphones K intr from Wu, Chen, &Chen "Visual Positioning Indoors: Human Eyes vs SmartPhone Cameras"
        Table 1
        Model  Xiaomi 5 Huawei P9 Samsung Note5 Lenovo Tango iPhone 7P
          fx    3831.011 3096.023 4048.113 3854.211 3289.89
          fy    3832.273 3096.611 4046.466 3851.217 3289.17
          ox    1844.276 1482.911 2587.339 1492.329 1991.804
          oy    2226.916 1982.791 1556.018 2692.189 1491.939
         */

        long seed = System.currentTimeMillis();
        System.out.println("seed=" + seed);
        Random rand = new Random(seed);

        int nTests = 2;
        int zInit = 1;
        double[][] x;
        double[][] XW0, XW;

        for (int i = 0; i < nTests; ++i) {

            int nP = 100 + rand.nextInt(200);

            // image size [400 x 650]
            //XW0 = randomPoints3D(rand, new int[]{211, 106, -1584}, nP);
            XW0 = randomPoints3D(rand, new int[]{1000, 2000, -100}, nP);

            // x = P * X
            //x = MatrixUtil.multiply(p, XW0);
            x = MatrixUtil.multiply(pNotDivByZ, XW0);
            normalize(x);

            //int np2 = filter(x, XW);
            //x = MatrixUtil.copySubMatrix(x, 0, x.length-1, 0, x[0].length-1);
            //XW = MatrixUtil.copySubMatrix(XW, 0, XW.length-1, 0, XW[0].length-1);

            XW = MatrixUtil.copySubMatrix(XW0, 0, 2, 0, XW0[0].length - 1);

            Camera.CameraPoseParameters result = CameraPose.calculatePoseFromXXW(x, XW);

            double[][] resultR = result.getExtrinsicParameters().getRotation();
            double[] resultT = result.getExtrinsicParameters().getTranslation();
            double[][] resultK = result.getIntrinsicParameters().getIntrinsic();

            diffR = Rotation.procrustesAlgorithmForRotation(rot, result.getExtrinsicParameters().getRotation());
            fsR = MatrixUtil.frobeniusNorm( MatrixUtil.pointwiseSubtract(
                    diffR, MatrixUtil.createIdentityMatrix(3)));

            diffT = MatrixUtil.subtract(translation, resultT);
            double[][] diffK = MatrixUtil.pointwiseSubtract(K, resultK);
            double diffTSum = MatrixUtil.lPSum(diffT, 2);
            double diffKSum = MatrixUtil.frobeniusNorm(diffK);

            double resultLambda1 = result.getIntrinsicParameters().getLambda1();

            assertTrue(diffTSum < 1);
            assertTrue(diffKSum < 1);
            assertTrue(fsR < 1E-3);
        }

    }

    public void testZhangExample_Bougetcode() throws Exception {

        // ======== first checking calculatePoseFromP with Zhang data example ================

        // from chap 2 of Camera Calibration book, chapter by Zhang (add reference)
        //https://people.cs.rutgers.edu/~elgammal/classes/cs534/lectures/CameraCalibration-book-chapter.pdf
        double[][] p = new double[][]{
                {7.025659E-1, -2.861189E-2, -5.377696E-1, 6.24189E1},
                {2.077632E-1, 1.265804, 1.591456E-1, 1.075646E1},
                {4.634764E-4, -5.282382E-5, 4.255347E-4, 1}
        };

        boolean passive = true;//false;
        double angle = 47.7*Math.PI/180.;

        // wide angle: f <= 35 mm
        // telephoto: f > 35mm
        // K in terms of
        // focal lengths in units of mm in positions[0][0] and [1][1]
        // and the optical center in units of pixel coordinates in [0][2] and [1][2]
        // position [0][1] holds the pixel skew if any.
        //  there is skew when the camera sensor is not perpendicular to optical axis.
        double[][] K = new double[][]{
                {1380.12, 0, 246.52},
                {0, 2032.57, 243.68},
                {0, 0, 1}
        };

        // angle between 2 image axes is 90 degrees
        // aspect ratio of the pixels = alpha/beta = 0.679
        // rotation axis: (almost vertical)
        double[] axis = new double[]{-0.08573, -0.99438, 0.0621};
        double[] rotVec = Rotation.createRotationVectorFromAngleAxis(axis, angle);

        double[][] rot = Rotation.createRotationRodriguesFormula(rotVec, passive);
        double[] eulerAngles = Rotation.extractThetaFromXYZ(rot, passive);
        double[] eulerAnglesDegrees = Arrays.copyOf(eulerAngles, eulerAngles.length);
        MatrixUtil.multiply(eulerAnglesDegrees, 180./Math.PI);

        /*
        expected rotation:
        [0.6750834808562537, -0.018079003106068848, -0.7375197919525242]
        [0.07388502469134674, 0.9963303528698166, 0.043206840623514806]
        [0.7340322179586962, -0.08365989240369291, 0.6739419302889822]
         */

        // t in mm = 1.5m from camera in WCS
        double[] translation = new double[]{-211.28, -106.06, 1583.75};
        double tZ = translation[2];

        // camera projection matrix P is defined in
        //     z*x = P*X
        //         = K*[R|T]*X
        //  where z is the distance to the object from the camera
        // Zhang alters his example P to move z from the RHS to the LHS:
        //      (K/z) * [R | T]
        double[][] reproducingP =
                MatrixUtil.multiply(
                        new double[][] {
                                {K[0][0]/tZ, K[0][1]/tZ, K[0][2]/tZ},
                                {K[1][0]/tZ, K[1][1]/tZ, K[1][2]/tZ},
                                {K[2][0]/tZ, K[2][1]/tZ, K[2][2]/tZ}},

                        new double[][] {
                                {rot[0][0], rot[0][1], rot[0][2], translation[0]},
                                {rot[1][0], rot[1][1], rot[1][2], translation[1]},
                                {rot[2][0], rot[2][1], rot[2][2], translation[2]}
                        }
                );

        // taking a look at P composed from K instead of K/z gives even better matches to his example params (excepting p)
        double[][] pNotDivByZ = new double[3][4];
        for (int ii = 0; ii < 3; ++ii) {
            System.arraycopy(rot[ii], 0, pNotDivByZ[ii], 0, 3);
            pNotDivByZ[ii][3] = translation[ii];
        }
        pNotDivByZ = MatrixUtil.multiply(K, pNotDivByZ);

        // generate a set of features for an image
        // create a camera intrinsic matrix K,
        // motion w.r.t. image and scene, rotation R, and translation t
        // and then projected WCS features.
        // test that can recover the projection matrix

        long seed = System.currentTimeMillis();
        System.out.println("seed=" + seed);
        Random rand = new Random(seed);

        int nTests = 2;
        int zInit = 1;
        double[][] x;
        double[][] XW0, XW;

        Camera.CameraIntrinsicParameters intrC = new Camera.CameraIntrinsicParameters();
        intrC.setIntrinsic(K);

        for (int i = 0; i < nTests; ++i) {

            int nP = 100 + rand.nextInt(200);

            // image size [400 x 650]
            //XW0 = randomPoints3D(rand, new int[]{211, 106, -1584}, nP);
            XW0 = randomPoints3D(rand, new int[]{1000, 2000, -100}, nP);

            // x = P * X
            //x = MatrixUtil.multiply(p, XW0);
            x = MatrixUtil.multiply(pNotDivByZ, XW0);
            normalize(x);

            //int np2 = filter(x, XW);
            //x = MatrixUtil.copySubMatrix(x, 0, x.length-1, 0, x[0].length-1);
            //XW = MatrixUtil.copySubMatrix(XW, 0, XW.length-1, 0, XW[0].length-1);

            XW = MatrixUtil.copySubMatrix(XW0, 0, 2, 0, XW0[0].length - 1);

            for (int ii = 0; ii < 4; ++ii) {
                boolean refine, useBouguetForRodrigues;
                switch (ii) {
                    case 0: refine = false; useBouguetForRodrigues = false; break;
                    case 1: refine = false; useBouguetForRodrigues = true; break;
                    case 2: refine = true; useBouguetForRodrigues = false; break;
                    default: refine = true; useBouguetForRodrigues = true; break;
                }

                Camera.CameraExtrinsicParameters result = CameraPose.calculatePoseUsingBouguet(intrC, x, XW, refine,
                        useBouguetForRodrigues);

                double[][] resultR = result.getRotation();
                double[] resultT = result.getTranslation();

                double[][] diffR = Rotation.procrustesAlgorithmForRotation(rot, resultR);
                double fsR = MatrixUtil.frobeniusNorm(MatrixUtil.pointwiseSubtract(
                        diffR, MatrixUtil.createIdentityMatrix(3)));

                double[] diffT = MatrixUtil.subtract(translation, resultT);
                double diffTSum = MatrixUtil.lPSum(diffT, 2);

                assertTrue(diffTSum < 1);
                assertTrue(fsR < 1E-3);
            }
        }

    }

    private int filter(double[][] x, double[][] XW) {
        int np = x[0].length;
        int n = 0;
        int j = 0;
        for (int i = 0; i < np; ++i) {
            if (x[0][i] >= 0 && x[1][i] >= 0) {
                if (i != j) {
                    for (int k = 0; k < x.length; ++k) {
                        x[k][j] = x[k][i];
                    }
                    for (int k = 0; k < XW.length; ++k) {
                        XW[k][j] = XW[k][i];
                    }
                }
                ++n;
                ++j;
            }
        }
        return n;
    }

    private double[] randomTranslation(Random rand) {
        double[] t = new double[]{rand.nextInt(100), rand.nextInt(100), 1};
        if (rand.nextBoolean()) {
            t[0] *= -1;
        }
        if (rand.nextBoolean()) {
            t[1] *= -1;
        }
        return t;
    }

    private double[] randomTranslationSmall(Random rand) {
        double[] t = new double[]{rand.nextDouble(), rand.nextDouble(), 1};
        if (rand.nextBoolean()) {
            t[0] *= -1;
        }
        if (rand.nextBoolean()) {
            t[1] *= -1;
        }
        return t;
    }

    private double[][] randomPoints3D(Random rand, int[] xyzC, int n) {
        double[][] x = new double[4][n];
        Arrays.fill(x[3], 1);
        double r;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < 3; ++j) {
                r = rand.nextInt(Math.abs(xyzC[j]) / 10);
                if (rand.nextBoolean()) {
                    r *= -1;
                }
                x[j][i] = xyzC[j] + r;
            }
        }
        return x;
    }

    //det(R)=1 is a proper rotation matrix.  rotation angles are counterclockwise.
    private double[][] randomRotation(Random rand) {
        double x = rand.nextInt(360) * Math.PI/180.;
        double y = rand.nextInt(360) * Math.PI/180.;
        double z = rand.nextInt(360) * Math.PI/180.;

        double[][] r = Rotation.createRotationXYZ(x, y, z);
        double det = MatrixUtil.determinant(r);
        assert(Math.abs(det - 1) < 1E-7);
        return r;
    }

    protected void normalize(double[][] x) {
        for (int i = 0; i < x[0].length; ++i) {
            for (int j = 0; j < x.length; ++j) {
                x[j][i] /= x[x.length - 1][i];
            }
        }
    }
}
