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

    /**
     * Test of calculatePoseUsingDLT method, of class CameraPose.
     */
    public void estCalculatePoseUsingZhangeData() throws Exception {

        System.out.println("testCalculatePoseUsingDLT");
        // see testresources/zhang1998/README.txt

        // they use f(r) = 1 + k1*r + k2*r^2:
        boolean useR2R4 = false;

        int nFeatures = 256;
        int nImages = 5;

        // 3 X 256 size
        // range of values is [0, 17/2.54] w/ center at (0.5*17/2.54)
        double[][] xW = Zhang98Data.getFeatureWCS();
        assertEquals(3, xW.length);
        assertEquals(nFeatures, xW[0].length);

        double[][] expectedKIntr = Zhang98Data.getIntrinsicCameraMatrix();
        double[] radial = Zhang98Data.getRadialDistortionR2R4();
        Camera.CameraIntrinsicParameters expectedIntr = new Camera.CameraIntrinsicParameters(
                expectedKIntr, radial, true);

        double[][] rT0 = new double[5][];
        double[][] rT2 = new double[5][];
        double[][] rTE = new double[5][];

        for (int i = 1; i <= 5; ++i) {

            //The planar homography pose method extracts R and t much better than
            // the method using a projective matrix. see the distances beteen expected rotation
            // and estimated rotations below.

            System.out.println();
            //3 X 256
            double[][] x = Zhang98Data.getObservedFeaturesInImage(i);

            double[][] expectedR = Zhang98Data.getRotation(i);
            double detR = MatrixUtil.determinant(expectedR);
            double[] expectedT = Arrays.copyOf(Zhang98Data.getTranslation(i), 3);
            MatrixUtil.multiply(expectedT, 1./expectedT[2]);

            double[] qHamiltonE = Rotation.createHamiltonQuaternionZYX(Rotation.extractThetaFromZYX(expectedR));
            double[] qBarfootE = Rotation.convertHamiltonToBarfootQuaternion(qHamiltonE);

            // avg is -0.2946225
            double[] zMeanAndStd = null;
            double _zMean = -0.295;

            // debug, check using bouguet
            Triangulation.WCSPt[] wcs = null;

            double[] WCSZ = null;
            if (i < 5) {
                wcs = new Triangulation.WCSPt[x[0].length];
                WCSZ = new double[x[0].length];

                // changing to camera coordinates to put a distance from aperture back into
                // the coordinates for use in triangulation:

                double[][] x2 = Zhang98Data.getObservedFeaturesInImage(i+1);
                double[][] expectedR2 = Zhang98Data.getRotation(i+1);
                double[] expectedT2 = Arrays.copyOf(Zhang98Data.getTranslation(i+1), 3);
                MatrixUtil.multiply(expectedT2, 1./expectedT2[2]);

                //3 x N for a single point
                for (int j = 0; j < x[0].length; ++j) {
                    double[][] _x1 = new double[][]{MatrixUtil.extractColumn(x, j)};
                    double[][] _x2 = new double[][]{MatrixUtil.extractColumn(x2, j)};
                    _x1 = MatrixUtil.transpose(_x1);
                    _x2 = MatrixUtil.transpose(_x2);

                    // the predicted WCS are too small by a factor of 10-ish, but triangulation algorithm is
                    // correct.
                    // it looks like the image data as x, y are not homogeneous coordinates truncated to x,y
                    // so not a good idea to use for pose testing.  need one more dimension to the data.

                    wcs[j] = Triangulation.calculateWCSPoint(expectedKIntr, expectedR, expectedT,
                            expectedKIntr, expectedR2, expectedT2, _x1, _x2);
                    WCSZ[j] = wcs[j].X[2];
                }
                zMeanAndStd = MiscMath.getAvgAndStDev(WCSZ);

                // x = alpha * P * X
                // alpha is 1/depth of point
                System.out.printf("from triangulation (Z = %.5f +- %.5f):\n", zMeanAndStd[0], zMeanAndStd[1]);
            }

            boolean useBouguetForRodrigues = false;
            //double[] om = qHamiltonE;
            //double[] om = qBarfootE;
            double[] om = Rotation.extractRodriguesRotationVectorBouguet(expectedR).om;
            CameraPose.ProjectedPoints xBou = CameraPose.bouguetProjectPoints2(
                    xW, om, expectedT, expectedIntr, useBouguetForRodrigues);

            // the bouguet results are reasonable, so have an error in calculatePoseFromXXW
            boolean refine = false;
            boolean useBouguetsRodrigues = false;
            Camera.CameraExtrinsicParameters c = CameraPose.calculatePoseUsingBouguet(
                    expectedIntr, x, xW, refine, useBouguetsRodrigues);

            double[][] _xBou = new double[3][xW[0].length];
            Arrays.fill(_xBou[2], 1);
            System.arraycopy(xBou.xEst[0], 0, _xBou[0], 0, xBou.xEst[0].length);
            System.arraycopy(xBou.xEst[1], 0, _xBou[1], 0, xBou.xEst[1].length);

            //double[][] xc = Camera.pixelToCameraCoordinates(x, expectedKIntr, radial, true);
            Camera.CameraPoseParameters result = CameraPose.calculatePoseFromXXW(x, xW);
            Camera.CameraPoseParameters result11 = CameraPose.calculatePoseFromXXW(_xBou, xW);
            // looks like the results from calculatePoseFromXXW need refinement

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

            // find distances between rotations
            rT0[i-1] = Rotation.extractThetaFromZYX(extr.getRotation());
            rT0[i-1] = Arrays.copyOf(rT0[i-1], rT0[i-1].length);
            MatrixUtil.multiply(rT0[i-1], 180./Math.PI);

            rTE[i-1] = Rotation.extractThetaFromZYX(expectedR);
            rTE[i-1] = Arrays.copyOf(rTE[i-1], rTE[i-1].length);
            MatrixUtil.multiply(rTE[i-1], 180./Math.PI);

            double[] qHamilton = Rotation.createHamiltonQuaternionZYX(Rotation.extractThetaFromZYX(extr.getRotation()));
            double[] qBarfoot = Rotation.convertHamiltonToBarfootQuaternion(qHamilton);
            double distR0 = Rotation.distanceBetweenQuaternions(qBarfootE, qBarfoot);

            System.out.printf("%d) r=\n%s\n", i, FormatArray.toString(extr.getRotation(), "%.3e"));
            System.out.printf("    r expected=\n%s\n", FormatArray.toString(expectedR, "%.3e"));
            System.out.printf("dist r     =%.3e\n", distR0);
            System.out.printf("%d) kIntr=\n%s\n", i, FormatArray.toString(intr.getIntrinsic(), "%.3e"));
            System.out.printf("    kIntr expected=\n%s\n", FormatArray.toString(expectedKIntr, "%.3e"));
            System.out.printf("%d) tc1=\n%s\n", i, FormatArray.toString(tc1, "%.3e"));
            System.out.printf("%d) _tc1=\n%s\n", i, FormatArray.toString(_tc1, "%.3e"));
            System.out.printf("%d) tc3=\n%s\n", i, FormatArray.toString(tc3, "%.3e"));
            System.out.printf("%d) ti2=\n%s\n", i, FormatArray.toString(ti2, "%.3e"));
            System.out.printf("%d) ti4=\n%s\n", i, FormatArray.toString(ti4, "%.3e"));
            System.out.printf("%d) _ti4=\n%s\n", i, FormatArray.toString(_ti4, "%.3e"));
            System.out.printf("    t expected=\n%s\n", FormatArray.toString(expectedT, "%.3e"));

        }

        // looking offset in x, y, z rotations for the solutions that would put
        //    them into agreement with the expected, within reasonable stdevs of their means.

        // TODO: need to redo these for angle subtraction

        double[] mean = new double[3];
        double[] stdev = new double[3];
        double[] diff = new double[3];
        System.out.printf("rThetas1 (degrees)=\n");
        for (int i = 0; i < 5; ++i) {
            System.out.printf("  %s\n", FormatArray.toString(rT0[i], "%.3f"));
            MatrixUtil.pointwiseSubtract(rT0[i], rTE[i], diff);
            for (int j = 0; j < 3; ++j) {
                mean[j] += diff[j];
            }
        }
        for (int j = 0; j < 3; ++j) {
            mean[j] /= 5.;
        }
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 3; ++j) {
                stdev[j] += (diff[j] - mean[j])*(diff[j] - mean[j]);
            }
        }
        for (int j = 0; j < 3; ++j) {
            stdev[j] = Math.sqrt(stdev[j]/4.);
        }
        System.out.printf("rThetas1 mean diff from expected, stdev of mean diff from expected (degrees)=\n   %s\n   %s\n",
                FormatArray.toString(mean, "%.3f"), FormatArray.toString(stdev, "%.3f"));

        Arrays.fill(mean, 0);
        Arrays.fill(stdev, 0);
        System.out.printf("rThetas2 (degrees)=\n");
        for (int i = 0; i < 5; ++i) {
            System.out.printf("  %s\n", FormatArray.toString(rT2[i], "%.3f"));
            MatrixUtil.pointwiseSubtract(rT2[i], rTE[i], diff);
            for (int j = 0; j < 3; ++j) {
                mean[j] += diff[j];
            }
        }
        for (int j = 0; j < 3; ++j) {
            mean[j] /= 5.;
        }
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 3; ++j) {
                stdev[j] += (diff[j] - mean[j])*(diff[j] - mean[j]);
            }
        }
        for (int j = 0; j < 3; ++j) {
            stdev[j] = Math.sqrt(stdev[j]/4.);
        }
        System.out.printf("rThetas2 mean diff from expected, stdev of mean diff from expected (degrees)=\n   %s\n   %s\n",
                FormatArray.toString(mean, "%.3f"), FormatArray.toString(stdev, "%.3f"));

        System.out.printf("rThetasE (degrees)=\n");
        for (int i = 0; i < 5; ++i) {
            System.out.printf("  %s\n", FormatArray.toString(rTE[i], "%.3f"));
        }

    }

    public void estCalculatePoseUsingBouguet() throws IOException, NotConvergedException {

        System.out.println("testCalculatePoseUsingBouguet");
        // following unit test at bottom of
        // https://github.com/fragofer/TOOLBOX_calib/
        // compute_extrinsic.m

        int Np = 4;
        double[] s = new double[]{10, 10, 5};

        long seed = System.nanoTime();
        //seed = 691498835451920L; // this with Np=6 (reproducible for same Random.java implementations):
        //                             om= 1.9923e-01, 1.3144e-01, 3.8956e-01
        //                          omckk= 1.8306e-02, 5.2196e-01, 3.5149e+00
        //                          om-omckk= 1.8093e-01, -3.9053e-01, -3.1254e+00 <=== last is ~ 360 off, so need wrap correction for diffs
        //                              T= 0.0000e+00, 0.0000e+00, 1.0000e+02
        //                           Tckk= -6.9559e-01, -6.8940e-01, -1.0996e+02  <===-1*T
        System.out.println("seed = " + seed);
        Random random = new Random(seed);

        //om = randn(3,1);
        int i, j;
        double[] om = new double[3];
        for (i = 0; i < 3; ++i) {
            om[i] = random.nextDouble();
        }
        //T = [0;0;100];
        double[] T = new double[]{0, 0, 100}; //[3X1]
        //noise = 2/1000;
        double noise = 2./1000;
        //XX = [sx*randn(1,Np);sy*randn(1,Np);sz*randn(1,Np)];
        double[][] XX = new double[3][];

        for (i = 0; i < 3; ++i) {
            XX[i] = new double[Np];
            for (j = 0; j < Np; ++j) {
                XX[i][j] = s[i] * random.nextDouble();
            }
        }
        //xx = project_points(XX,om,T);
        double[] radial = new double[]{1e-7, 1e-7};
        Camera.CameraIntrinsicParameters intr = new Camera.CameraIntrinsicParameters(
            MatrixUtil.createIdentityMatrix(3), radial, true);
        boolean useBouguetsRodrigues = false;
        CameraPose.ProjectedPoints pp = CameraPose.bouguetProjectPoints2(XX, om, T, intr, useBouguetsRodrigues);
        double[][] xx = pp.xEst;

        System.out.printf("om=\n%s\n", FormatArray.toString(om, "%.3e"));
        System.out.printf("T=\n%s\n", FormatArray.toString(T, "%.3e"));
        //System.out.printf("XX=\n%s\n", FormatArray.toString(XX, "%.3e"));
        //System.out.printf("xEst=\n%s\n", FormatArray.toString(xx, "%.3e"));

        //xxn = xx + noise * randn(2,Np); //[2Xn]
        double[][] noiseM = new double[2][];
        for (i = 0; i < 2; ++i) {
            noiseM[i] = new double[Np];
            for (j = 0; j < Np; ++j) {
                noiseM[i][j] = noise*random.nextDouble();
            }
        }
        // and so Rtransform is w.r.t. the real world origin of XW whether or not the user calibrated that to be
        //    the center of the camera.
        //[2Xn]
        double[][] xxn = MatrixUtil.pointwiseAdd(xx, noiseM);

        //System.out.printf("xxn = xEst with noise =\n%s\n", FormatArray.toString(xxn, "%.3e"));

        //[omckk,Tckk] = compute_extrinsic(xxn,XX);
        boolean refine = true;
        Camera.CameraExtrinsicParameters c = CameraPose.calculatePoseUsingBouguet(intr, xxn, XX, refine, useBouguetsRodrigues);
        //[om omckk om-omckk]
        //[T Tckk T-Tckk]
        System.out.printf("om=\n  %s\n", FormatArray.toString(om, "%.4e"));
        System.out.printf("omckk=\n  %s\n", FormatArray.toString(c.getRodriguesVector(), "%.4e"));
        System.out.printf("om-omckk=\n  %s\n", FormatArray.toString(
                MatrixUtil.subtract(om, c.getRodriguesVector()), "%.4e"));

        System.out.printf("T=\n  %s\n", FormatArray.toString(T, "%.4e"));
        System.out.printf("Tckk=\n  %s\n", FormatArray.toString(c.getTranslation(), "%.4e"));
        System.out.printf("T-Tckk=\n  %s\n", FormatArray.toString(
                MatrixUtil.subtract(T, c.getTranslation()), "%.4e"));

        // expected errors?  can see a rough dependence of 1/(Np^2)
        // for nP=100:
        double oTol = 1E-3;
        double tTol = 5E-1;
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

        boolean passive = false;
        double rotAboutAxis= 47.7*Math.PI/180.;

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
        double[] expectedRAxis = new double[]{-0.08573, -0.99438, 0.0621};

        double[][] expectedRot = Rotation.createRotationFromUnitLengthAngleAxis(expectedRAxis,
                rotAboutAxis, passive);

        /*
        expected rotation:
        [0.6750834808562537, -0.018079003106068848, -0.7375197919525242]
        [0.07388502469134674, 0.9963303528698166, 0.043206840623514806]
        [0.7340322179586962, -0.08365989240369291, 0.6739419302889822]
         */

        // t in mm = 1.5m from camera in WCS
        double[] expectedT = new double[]{-211.28, -106.06, 1583.75};

        // camera projection matrix P is defined in
        //     z*x = P*X
        //         = K*[R|T]*X
        //  where z is the distance to the object from the camera
        // Zhang alters his example P to move z from the RHS to the LHS:
        //      (K/z) * [R | T]
        double[][] reproducingP =
            MatrixUtil.multiply(
                new double[][] {
                {K[0][0]/expectedT[2], K[0][1]/expectedT[2], K[0][2]/expectedT[2]},
                {K[1][0]/expectedT[2], K[1][1]/expectedT[2], K[1][2]/expectedT[2]},
                {K[2][0]/expectedT[2], K[2][1]/expectedT[2], K[2][2]/expectedT[2]}},

                 new double[][] {
                         {expectedRot[0][0], expectedRot[0][1], expectedRot[0][2], expectedT[0]},
                         {expectedRot[1][0], expectedRot[1][1], expectedRot[1][2], expectedT[1]},
                         {expectedRot[2][0], expectedRot[2][1], expectedRot[2][2], expectedT[2]}
                 }
        );

        Camera.CameraPoseParameters c = CameraPose.calculatePoseFromP(p);
        double fs = MatrixUtil.frobeniusNorm(
                MatrixUtil.pointwiseSubtract(K, c.getIntrinsicParameters().getIntrinsic()));

        double[] r = Rotation.extractRotationAxisFromZXY(c.getExtrinsicParameters().getRotation());
        // should be identity matrix if same:
        double[][] diffR = Rotation.procrustesAlgorithmForRotation(expectedRot, c.getExtrinsicParameters().getRotation());
        double fsR = MatrixUtil.frobeniusNorm( MatrixUtil.pointwiseSubtract(
                diffR, MatrixUtil.createIdentityMatrix(3)));

        double[] diffT = MatrixUtil.subtract(expectedT, c.getExtrinsicParameters().getTranslation());

        assertTrue(fs < 1);
        assertTrue(fsR < 1);

        // taking a look at P composed from K instead of K/z gives even better matches to his example params (excepting p)
        double[][] pNotDivByZ = new double[3][4];
        for (int ii = 0; ii < 3; ++ii) {
            System.arraycopy(expectedRot[ii], 0, pNotDivByZ[ii], 0, 3);
            pNotDivByZ[ii][3] = expectedT[ii];
        }
        pNotDivByZ = MatrixUtil.multiply(K, pNotDivByZ);

        Camera.CameraPoseParameters _c = CameraPose.calculatePoseFromP(pNotDivByZ);
        double _fs = MatrixUtil.frobeniusNorm(
                MatrixUtil.pointwiseSubtract(K, _c.getIntrinsicParameters().getIntrinsic()));

        double[] _r = Rotation.extractRotationAxisFromZXY(_c.getExtrinsicParameters().getRotation());
        // should be identity matrix if same:
        double[][] _diffR = Rotation.procrustesAlgorithmForRotation(expectedRot, _c.getExtrinsicParameters().getRotation());
        double _fsR = MatrixUtil.frobeniusNorm( MatrixUtil.pointwiseSubtract(
                _diffR, MatrixUtil.createIdentityMatrix(3)));

        double[] _diffT = MatrixUtil.subtract(expectedT, _c.getExtrinsicParameters().getTranslation());

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
            XW0 = randomPoints3D(rand, new int[]{211, 106, -1584}, nP);

            // x = P * X
            x = MatrixUtil.multiply(p, XW0);
            normalize(x);

            //int np2 = filter(x, XW);
            //x = MatrixUtil.copySubMatrix(x, 0, x.length-1, 0, x[0].length-1);
            //XW = MatrixUtil.copySubMatrix(XW, 0, XW.length-1, 0, XW[0].length-1);

            XW = MatrixUtil.copySubMatrix(XW0, 0, 2, 0, XW0[0].length - 1);

            Camera.CameraPoseParameters result = CameraPose.calculatePoseFromXXW(x, XW);

            double[][] resultR = result.getExtrinsicParameters().getRotation();
            double[] resultT = result.getExtrinsicParameters().getTranslation();
            double[][] resultK = result.getIntrinsicParameters().getIntrinsic();
            double resultLambda1 = result.getIntrinsicParameters().getLambda1();

            diffR = Rotation.procrustesAlgorithmForRotation(expectedRot, result.getExtrinsicParameters().getRotation());
            fsR = MatrixUtil.frobeniusNorm( MatrixUtil.pointwiseSubtract(
                    diffR, MatrixUtil.createIdentityMatrix(3)));

            diffT = MatrixUtil.subtract(expectedT, resultT);
            double[][] diffK = MatrixUtil.pointwiseSubtract(K, resultK);
            double diffTSum = MatrixUtil.lPSum(diffT, 2);
            double diffKSum = MatrixUtil.frobeniusNorm(diffK);

            //=============== ================ ==================== ===============
            // generate 3 sets of same camera views of WCS and format the data for planar camera method input

            /*

            estimate the angles and offsets to move a camera to see the scene

            // the rotation axis is essentially the y-axis.
            //    rotation cc about the y-axis (pitch)
            // translate camera in x and small rotation
            double[][] expectedRot2 = Rotation.createRotationFromUnitLengthAngleAxis(expectedRAxis,
                    rotAboutAxis*0.99, passive);
            double[] expectedT2 = Arrays.copyOf(expectedT, expectedT.length);
            for (int j = 0; j < expectedT2.length; ++j) {
                expectedT2[j] -= 10;
            }
            double[][] p2 = new double[3][4];
            for (int ii = 0; ii < 3; ++ii) {
                System.arraycopy(expectedRot2[ii], 0, p2[ii], 0, 3);
                p2[ii][3] = expectedT2[ii];
            }
            p2 = MatrixUtil.multiply(K, p2);

            double[][] expectedRot3 = Rotation.createRotationFromUnitLengthAngleAxis(expectedRAxis,
                    rotAboutAxis*0.98, passive);
            double[] expectedT3 = Arrays.copyOf(expectedT, expectedT.length);
            for (int j = 0; j < expectedT3.length; ++j) {
                expectedT3[j] -= 10;
            }
            double[][] p3 = new double[3][4];
            for (int ii = 0; ii < 3; ++ii) {
                System.arraycopy(expectedRot3[ii], 0, p3[ii], 0, 3);
                p3[ii][3] = expectedT3[ii];
            }
            p3 = MatrixUtil.multiply(K, p3);

            double[][] x1 = MatrixUtil.multiply(pNotDivByZ, XW0);
            normalize(x1);
            double[][] x2 = MatrixUtil.multiply(p2, XW0);
            normalize(x2);
            double[][] x3 = MatrixUtil.multiply(p3, XW0);
            normalize(x3);
            double[][] xAll = new double[3][x[0].length * 3];
            for (int ii = 0; ii < 3; ++ii) {
                System.arraycopy(x1[ii], 0, xAll[ii], 0, x1[ii].length);
                System.arraycopy(x2[ii], 0, xAll[ii], x1[ii].length, x1[ii].length);
                System.arraycopy(x3[ii], 0, xAll[ii], 2*x1[ii].length, x1[ii].length);
            }

            boolean useR24 = false;
            Camera.CameraMatrices ms = CameraCalibration.estimateCameraPlanar(x[0].length, xAll, XW, useR24);

            double[][] diffR2 = Rotation.procrustesAlgorithmForRotation(expectedRot, ms.getExtrinsics().get(0).getRotation());
            double[] diffT2 = MatrixUtil.subtract(expectedT, ms.getExtrinsics().get(0).getTranslation());
            double[][] diffK2 = MatrixUtil.pointwiseSubtract(K, ms.getIntrinsics().getIntrinsic());
            double diffTSum2 = MatrixUtil.lPSum(diffT2, 2);
            double diffKSum2 = MatrixUtil.frobeniusNorm(diffK2);

             */

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
