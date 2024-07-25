package algorithms.imageProcessing.transform;

import static algorithms.imageProcessing.transform.Rotation.extractThetaFromZYX;

import algorithms.matrix.MatrixUtil;
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
     * NOTE: not using this because need the image points to be homogeneous coordinates.
     */
    public void estCalculatePoseUsingZhangeData() throws Exception {

        System.out.println("testCalculatePoseUsingDLT");
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

            // debug, check using bouguet
            {

                Triangulation.WCSPt[] wcs = null;

                if (i < 5) {
                    wcs = new Triangulation.WCSPt[x[0].length];

                    //TODO: consider radial distortion corrections

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
                    }

                    // x = alpha * P * X
                    // alpha is 1/depth of point
                    System.out.printf("from triangulation:\n");
                    for (int j = 0; j < x[0].length; ++j) {
                        double[] _XW = MatrixUtil.extractColumn(xW, j);
                        System.out.printf("%d) alpha=%f\n", j, wcs[j].alpha);
                        System.out.printf("%s\n", FormatArray.toString(_XW, "%.5f"));
                        double[] _XWC = Arrays.copyOf(wcs[j].X, 4);
                        MatrixUtil.multiply(_XWC, 1./_XWC[3]);
                        System.out.printf("%s\n\n", FormatArray.toString(_XWC, "%.5f"));
                    }

                }

                boolean useBouguetForRodrigues = false;
                //double[] om = qHamiltonE;
                //double[] om = qBarfootE;
                double[] om = Rotation.extractRodriguesRotationVectorBouguet(expectedR).om;
                CameraPose.ProjectedPoints xBou = CameraPose.bouguetProjectPoints2(
                        xW, om, expectedT, expectedIntr, useBouguetForRodrigues);
            }
            double[][] xc = Camera.pixelToCameraCoordinates(x, expectedKIntr, radial, true);
            // remove distortion from x.  these should be in camera reference frame
            //double[][] xcU = CameraCalibration.removeRadialDistortion(xc, radial[0], radial[1], true);

            //Camera.CameraPoseParameters result = CameraPose.calculatePoseFromXXW(xcU, xW);
            Camera.CameraPoseParameters result = CameraPose.calculatePoseFromXXW(x, xW);

            Camera.CameraExtrinsicParameters extr = result.getExtrinsicParameters();
            Camera.CameraIntrinsicParameters intr = result.getIntrinsicParameters();
            double[] p3 = Arrays.copyOf(result.getP3(), result.getP3().length);

            Camera.CameraExtrinsicParameters extr2 = CameraPose.calculatePoseUsingCameraCalibration(
                    new Camera.CameraIntrinsicParameters(expectedKIntr, radial, useR2R4), x, xW);

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

            // find distances between rotations
            rT0[i-1] = Rotation.extractThetaFromZYX(extr.getRotation());
            rT0[i-1] = Arrays.copyOf(rT0[i-1], rT0[i-1].length);
            MatrixUtil.multiply(rT0[i-1], 180./Math.PI);

            rT2[i-1] = Rotation.extractThetaFromZYX(r2Orth);
            rT2[i-1] = Arrays.copyOf(rT2[i-1], rT2[i-1].length);
            MatrixUtil.multiply(rT2[i-1], 180./Math.PI);

            rTE[i-1] = Rotation.extractThetaFromZYX(expectedR);
            rTE[i-1] = Arrays.copyOf(rTE[i-1], rTE[i-1].length);
            MatrixUtil.multiply(rTE[i-1], 180./Math.PI);

            double[] qHamilton = Rotation.createHamiltonQuaternionZYX(Rotation.extractThetaFromZYX(extr.getRotation()));
            double[] qBarfoot = Rotation.convertHamiltonToBarfootQuaternion(qHamilton);
            double distR0 = Rotation.distanceBetweenQuaternions(qBarfootE, qBarfoot);

            qHamilton = Rotation.createHamiltonQuaternionZYX(Rotation.extractThetaFromZYX(r2Orth));
            qBarfoot = Rotation.convertHamiltonToBarfootQuaternion(qHamilton);
            double distR2Orth = Rotation.distanceBetweenQuaternions(qBarfootE, qBarfoot);

            System.out.printf("%d) r=\n%s\n", i, FormatArray.toString(extr.getRotation(), "%.3e"));
            System.out.printf("%d) r2=\n%s\n", i, FormatArray.toString(extr2.getRotation(), "%.3e"));
            System.out.printf("%d) r2Orth=\n%s\n", i, FormatArray.toString(r2Orth, "%.3e"));
            System.out.printf("    r expected=\n%s\n", FormatArray.toString(expectedR, "%.3e"));
            System.out.printf("dist r     =%.3e\ndist r2Orth=%.3e\n", distR0, distR2Orth);
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

    public void testCalculatePoseUsingBouguet() throws IOException, NotConvergedException {

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
}
