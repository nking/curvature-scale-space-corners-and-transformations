package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Camera.CameraExtrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraMatrices;
import algorithms.matrix.MatrixUtil;
import algorithms.statistics.Standardization;
import algorithms.util.FormatArray;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;

import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

import static junit.framework.TestCase.assertEquals;

/**
 *
 * @author nichole
 */
public class CameraCalibrationTest extends TestCase {
    
    static double eps = 1.e-5;
    
    private static final Level LEVEL = Level.FINEST;
    private static final Logger log;
    static {
        log = Logger.getLogger(CameraCalibration.class.getSimpleName());
    }
    
    public CameraCalibrationTest() {
    }
    
    public void __testApplyRemoveRadialDistortion()
    throws Exception {
        
        log.log(LEVEL, "testApplyRemoveRadialDistortion");

        // test for model #3
        // test for model #4
        // test for k1<0, barrel distortion (happens for smaller focal lengths)
        // test for k1>0, pincushion distortion (happes for larger focal lenngths)
        
        // generate coords in radial annuli in area 100x100 then unit standard normalized
        double[][] coords = generateCoords0();
        
        log.log(LEVEL, String.format("orig coords=\n%s\n", FormatArray.toString(coords, "%.3f")));
        
        double[] k1s = new double[]{0.25, +0.25};
        double[] k2s = new double[]{0.2, +0.2};
        boolean[] useR2R4s = new boolean[]{false, true};
        int i, j;
        double k1, k2;
        double[][] xyd, xy;
        double diffDx, diffDy, diffx, diffy;
        for (boolean useR2R4 : useR2R4s) {
            for (i = 0; i < k1s.length; ++i) {
                k1 = k1s[i];
                k2 = k2s[i];
                log.log(LEVEL, "useR2R4=" + useR2R4 + "   k1="+ k1);
                xyd = CameraCalibration.applyRadialDistortion(coords, k1, k2, useR2R4);
                xy = CameraCalibration.removeRadialDistortion(xyd, k1, k2, useR2R4);
                // assert xyd != coordsI
                // assert xy == coordsI
                for (j = 0; j < coords[0].length; ++j) {
                    diffDx = Math.abs(coords[0][j] - xyd[0][j]);
                    diffDy = Math.abs(coords[1][j] - xyd[1][j]);
                    diffx = Math.abs(coords[0][j] - xy[0][j]);
                    diffy = Math.abs(coords[1][j] - xy[1][j]);
                    log.log(LEVEL, String.format("(%.4e,%.4e): diff(%.4e, %.4e)?  same(%.4e, %.4e)?\n", 
                        coords[0][j], coords[1][j], diffDx, diffDy, diffx, diffy));
                    assertTrue(diffx < eps);
                    assertTrue(diffy < eps);
                    // if barrel (k1<0), distorted is at a smaller radius
                }
            }
        }
        
        // apply radial distortion and remove radial distortion and compare
    }

    public void estCalibration3DZhang() throws Exception {
        // need the 3rd dimension to use this
        //log.log(LEVEL, "testCalibration3DZhang");
        //calibrationZhang(true);
    }

    public void testCalibrationPlanarZhang() throws Exception {
        log.log(LEVEL, "testCalibrationPlanarZhang");
        calibrationZhang(false);
    }

    protected void calibrationZhang(boolean use3D) throws Exception {

        // see testresources/zhang1998/README.txt
        
        // they use f(r) = 1 + k1*r + k2*r^2:
        boolean useR2R4 = true;
        
        int nFeatures = Zhang98Data.nFeatures;
        int nImages = Zhang98Data.mImages;
        
        // 3 X 256
        double[][] coordsW = Zhang98Data.getFeatureWCS();
        assertEquals(3, coordsW.length);
        assertEquals(nFeatures, coordsW[0].length);
        
        //3 X (256*5)
        double[][] coordsI = Zhang98Data.getObservedFeaturesInAllImages();
        assertEquals(3, coordsI.length);
        assertEquals(nFeatures*nImages, coordsI[0].length);

        CameraMatrices c = use3D ? CameraCalibration.estimateCamera3D(nFeatures, coordsI, coordsW, useR2R4) :
                 CameraCalibration.estimateCameraPlanar(nFeatures, coordsI, coordsW, useR2R4);

        double alphaE = 871.445;
        double gammaE = 0.2419;
        double u0E = 300.7676;
        double betaE = 871.1251;
        double v0E = 220.8684;
        double k1E = 0.1371;
        double k2E = -2.0101;

        Camera.CameraIntrinsicParameters kIntr = c.getIntrinsics();
        
        double alpha = kIntr.getIntrinsic()[0][0];
        double gamma = kIntr.getIntrinsic()[0][1];
        double u0 = kIntr.getIntrinsic()[0][2];
        double beta = kIntr.getIntrinsic()[1][1];
        double v0 = kIntr.getIntrinsic()[1][2];
        double[] kRadial = c.getRadialDistortCoeff();

        // TODO: add a round of smaller assertions once refinements are back in
        assertTrue(Math.abs(alphaE - alpha) < 1.0);
        assertTrue(Math.abs(gammaE - gamma) < 0.1);
        assertTrue(Math.abs(u0E - u0) < 1.0);
        assertTrue(Math.abs(betaE - beta) < 1.0);
        assertTrue(Math.abs(v0E - v0) < 1.0);

        List<CameraExtrinsicParameters> extrinsics = c.getExtrinsics();
        CameraExtrinsicParameters ex1;
        for (int i = 0; i < nImages; ++i) {
            ex1 = extrinsics.get(i);
            log.log(LEVEL, String.format("\nimg %d:\n", i));
            log.log(LEVEL, String.format("   r=%s\n", FormatArray.toString(ex1.getRotation(), "%.3e")));
            log.log(LEVEL, String.format("   t=%s\n", FormatArray.toString(ex1.getTranslation(), "%.3e")));

            double[][] rExp = Zhang98Data.getRotation(i+1);
            double[] tZYX = Rotation.extractThetaFromZYX(rExp);
            double[] tXYZ = Rotation.extractThetaFromXYZ(rExp);
            System.out.printf("THETAS ZYX (%d): %s\n", i+1, FormatArray.toString(tZYX, "%.4f"));
            System.out.printf("THETAS XYZ (%d): %s\n", i+1, FormatArray.toString(tXYZ, "%.4f"));

            double[] tExp = Zhang98Data.getTranslation(i+1);
            double[][] rDiff = Rotation.procrustesAlgorithmForRotation(rExp, ex1.getRotation());
            double fsR = MatrixUtil.frobeniusNorm( MatrixUtil.pointwiseSubtract(
                    rDiff, MatrixUtil.createIdentityMatrix(3)));

            double[] tRatio = MatrixUtil.pointwiseDivision(ex1.getTranslation(), tExp);
            log.log(LEVEL, String.format("   trans result/expected==%s\n", FormatArray.toString(tRatio, "%.3e")));
            log.log(LEVEL, String.format("   rSSD=%.3f\n", fsR));

            int t = 2;
        }

        log.log(LEVEL, String.format("k=%s\n", FormatArray.toString(kRadial, "%.4e")));
        log.log(LEVEL, String.format("expected k=%.3e, %.3e\n", k1E, k2E));

        // quick test of apply and remove distortion
        // (u_d, v_d) are the distorted features in coordsI in image reference frame.
        // (x_d, y_d) are the distorted features in the camera reference frame.
        // (x, y) are the distortion-free features in the camera reference frame.
        int i, j;
        double diff, diffD;
        double[][] uvDI, xyDI, xyi, xyDUI;
        
        double[][] cIntrE = new double[3][3];
        cIntrE[0] = new double[]{alphaE, gammaE, u0E};
        cIntrE[1] = new double[]{0, betaE, v0E};
        cIntrE[2] = new double[]{0, 0, 1};
        useR2R4 = false;
        Camera.CameraIntrinsicParameters kIntrE = new Camera.CameraIntrinsicParameters(
            cIntrE, null, useR2R4);
        for (i = 0; i < nImages; ++i) {
            uvDI = MatrixUtil.copySubMatrix(coordsI, 0, 2, nFeatures*i, nFeatures*(i + 1)-1);
            
            xyDI = Camera.pixelToCameraCoordinates(uvDI, kIntrE);
            xyi = CameraCalibration.removeRadialDistortion(xyDI, k1E, k2E);
            xyDUI = CameraCalibration.applyRadialDistortion(xyi, k1E, k2E, useR2R4);
            // assert that xy != xyDI
            // assert that xyDI == xyDUI
            for (j = 0; j < nFeatures; ++j) {
                diff = Math.abs(xyi[0][j] - xyDI[0][j]);
                //assertTrue(diff > 0.1);
                diffD = Math.abs(xyDUI[0][j] - xyDI[0][j]);
                log.log(LEVEL, String.format("(%.3f,%.3f): (distorted-undistorted)=%.5e \n     (orig - removedApplied)=%.5e\n", 
                    xyDI[0][j], xyDI[1][j], diff, diffD));
                System.out.flush();
                //assertTrue(diffD < 0.1);
            }
        }

    }

    public void __testPlanarRandom() throws Exception {
        log.log(LEVEL, "testPlanarRandom");

        long seed = System.currentTimeMillis();
        System.out.println("seed=" + seed);
        Random rand = new Random(seed);

        //=============== ================ ==================== ===============
        // generate 3 sets of same camera views of WCS and format the data for planar camera method input.
        //   the images are a sequence of motions.  frame

        double[][] K = new double[][]{
                {1380.12, 0, 246.52},
                {0, 2032.57, 243.68},
                {0, 0, 1}
        };

        boolean passive = false;

        int nPoints = 100;

        // 3 X nPoints
        // each point is [rand, rand, 1]
        double[][] xW = generateRandomPlane(nPoints, rand);

        int nImages = 3;

        // to keep points within FOV, will use a central rotation and translation then random small changes around those
        double rotAboutAxis= 47.7*Math.PI/180.;
        double[] rAxis = new double[]{-0.08573, -0.99438, 0.0621};
        double[][] rot = Rotation.createRotationRodriguesFormula(rAxis, rotAboutAxis, passive);
        double[] thetaXYZ = Rotation.extractThetaFromXYZ(rot, passive);
        double[][] _rot = Rotation.createRotationXYZ(thetaXYZ[0], thetaXYZ[1], thetaXYZ[2], passive);

        // t in mm = 1.5m from camera in WCS
        double[] trans = new double[]{-211.28, -106.06, 1583.75};

        // change thetaXYZ for the 2 images by [-0.2 to +0.2] for each axis
        // change expectedT for the 2 images by [-2 to +2]

        double[][] xAll = new double[3][nImages*nPoints];

        double[][] rAll = new double[3*nImages][3];
        double[][] tAll = new double[nImages][3];

        double[] _rThetaXYZ = new double[3];
        double[][] _r;
        double[] _t = new double[3];
        double[][] _h, _x = null;
        for (int i = 0; i < nImages; ++i) {
            for (int row = 0; row < 3; ++row) {
                _rThetaXYZ[row] = rand.nextDouble() / 2.5;
                if (rand.nextBoolean()) {
                    _rThetaXYZ[row] *= -1;
                }
                _rThetaXYZ[row] += thetaXYZ[row];

                _t[row] = row < 2 ? rand.nextDouble() * 5 : rand.nextDouble() * 2;
                if (rand.nextBoolean()) {
                    _t[row] *= -1;
                }
                _t[row] += trans[row];
            }

            _r = Rotation.createRotationXYZ(_rThetaXYZ[0], _rThetaXYZ[1], _rThetaXYZ[2]);

            _x = SceneImageHelper.createImagePoints2DPlanar(xW, K, _r, _t);

            for (int j = 0; j < _r.length; ++j) {
                System.arraycopy(_r[j], 0, rAll[i*3 + j], 0, _r[j].length);
            }
            for (int j = 0; j < _x.length; ++j) {
                System.arraycopy(_x[j], 0, xAll[j], i*_x[j].length, _x[j].length);
            }
            System.arraycopy(_t, 0, tAll[i], 0, _t.length);
        }

        boolean useR24 = false;
        Camera.CameraMatrices ms = CameraCalibration.estimateCameraPlanar(nPoints, xAll, xW, useR24);

        double fsK = MatrixUtil.frobeniusNorm(
                MatrixUtil.pointwiseSubtract(K, ms.getIntrinsics().getIntrinsic()));

        System.out.printf("K frob norm of I - procrustes(result K, expectedK) = %f\n", fsK);

        for (int i = 0; i < nImages; ++i) {
            double[][] _diffR = Rotation.procrustesAlgorithmForRotation(rot, ms.getExtrinsics().get(i).getRotation());
            double fsR = MatrixUtil.frobeniusNorm( MatrixUtil.pointwiseSubtract(
                    _diffR, MatrixUtil.createIdentityMatrix(3)));

            double[] _diffT = MatrixUtil.subtract(trans, ms.getExtrinsics().get(i).getTranslation());
            double ssdT = MatrixUtil.lPSum(_diffT, 2);

            System.out.printf("image %d, fsR=%.3f, ssdT=%.3f\n", i, fsR, ssdT);
        }

    }

    // range each point is [rand, rand, 1]
    private double[][] generateRandomPlane(int nPoints, Random rand) {
        double[][] xW = new double[3][nPoints];
        Arrays.fill(xW[2], 1);
        for (int i = 0; i < nPoints; ++i) {
            for (int j = 0; j < 2; ++j) {
                xW[j][i] = rand.nextDouble();
            }
        }
        return xW;
    }

    private double[][] generateWCS(Random rand, int nPoints, double[] boundsWCS) {
        ///example bounds=-7.47e-05, -3.31e-03, 2.03e-01, 2.48e+02, 5.90e-02, 1.20e-03
        double[][] XW = new double[4][nPoints];
        for (int i = 0; i < nPoints; ++i) {
            for (int j = 0; j < 4; ++j) {
                double diff = boundsWCS[2*j+1] - boundsWCS[2*j];
                XW[j][i] = boundsWCS[2*j] + rand.nextDouble()*diff;
            }
        }
        return XW;
    }

    protected void normalize(double[][] x, int normRow) {
        for (int i = 0; i < x[0].length; ++i) {
            for (int j = 0; j < x.length; ++j) {
                x[j][i] /= x[normRow][i];
            }
        }
    }

    /**
     * generate coordinates in 5 radial annuli in area 100x100 then unit standard normalized
     * @return 
     */
    private double[][] generateCoords0() {
        TDoubleList x = new TDoubleArrayList();
        TDoubleList y = new TDoubleArrayList();
        
        double rMax = 100;
        int nAnnuli = 1;
        double dR = rMax/nAnnuli;
        
        // for each annuli, want the 4 points crossing the x and y axes to test 0's.
        //     and the remaining points can be 8 points => all distributed at 2*pi/12 intervals
        double dT = 2.*Math.PI/12.;
        
        int i, j;
        double r, t;
        x.add(0);
        y.add(0);
        for (i = 1; i <= nAnnuli; ++i) {
            r = i*dR;
            for (j = 0; j < 12; ++j) {
                t = j*dT;
                //x = r cos t
                //y = r sin t
                x.add(r*Math.cos(t));
                y.add(r*Math.sin(t));
            }
        }
        
        // can use this transposed for unit standard normalization:
        //public static double[][] Standardization.standardUnitNormalization(double[][] data, 
        //    double[] outputMean, double[] outputStandardDeviation) {
        
        // write all x along column 0 instead of row 0.  and y goes in column 1
        double[][] normalized = new double[x.size()][2];
        for (i = 0; i < x.size(); ++i) {
            normalized[i] = new double[]{x.get(i), y.get(i)};
        }
        
        double[] outputMean = new double[2];
        double[] outputStDev = new double[2];
        normalized = Standardization.standardUnitNormalization(normalized,
            outputMean, outputStDev);
        
        // write for use in code above:
        //    row 0 = 's. row 1 = y's, row 2 = 1's
        double[][] coords = MatrixUtil.zeros(3, x.size());
        for (i = 0; i < x.size(); ++i) {
            coords[0][i] = normalized[i][0];
            coords[1][i] = normalized[i][1];
            coords[2][i] = 1.;
        }
        
        return coords;
    }

    
    public void estNeurochemistryBookPoses() throws Exception {
        
        //test pose
        //test triangulation
        
        /*
        images are 3024x4032
          errors likely > 8 pixels
        
            img2         img1            img3
        #1 678, 718     608, 530        744, 806
        #2 2210, 886    2462, 512       2286, 526
        #3 1504, 1102   1484, 858       1410, 992
        #4 1516, 1814   1486, 1678      1432, 1760
        #5 1228, 2014   1154, 1882      1164, 1940
        #6 2042, 1936   2142, 1810      2018, 1866
        #7 698, 2944    614, 2802       788, 2728
        #8 2210, 2782   2366, 2848      2276, 2930

          features in WCS
          #1 -11, 14, 41.5
          #2  11, 14, 41.5
          #3   0, 9.7, 41.5
          #4   0, 0, 41.5
          #5 -3.7, -3, 41.5
          #6 -8, -3, 41.5
          #7 -11, -14, 41.5
          #8  11, -14, 41.5
        
        expecting
             focalLength ~ 1604 pixels = 2.245 mm
             no skew
             xc=1512
             yc=2016
             little to no radial distortion (if was present, it is already removed)
             rotation between images = 23.4 degrees
             translation between images = 18 cm

          other information:
            pixel width = 1.4e-3mm
            FOV = 77 degrees = 1.344 radians
        */
        
        boolean useR2R4 = false;
        
        int nFeatures = NeurochemistryBookData.nFeatures;
        int nImages = NeurochemistryBookData.mImages;
        
        // 3 X (3*8)
        double[][] coordsW = NeurochemistryBookData.getFeatureWCS();
        assertEquals(3, coordsW.length);
        assertEquals(nFeatures, coordsW[0].length);
        
        //3 X (3*8))
        double[][] coordsI = NeurochemistryBookData.getObservedFeaturesInAllImages();
        assertEquals(3, coordsI.length);
        assertEquals(nFeatures*nImages, coordsI[0].length);
        
        //log.log(LEVEL, String.format("coordsW dimensions = [%d X %d]\ncoordsI dimensions = [%d X %d]\n",
        //    coordsW.length, coordsW[0].length, coordsI.length, coordsI[0].length));
        
        
        CameraMatrices cameraMatrices = CameraCalibration.estimateCameraPlanar(
            nFeatures, coordsI, coordsW, useR2R4);
        
        Camera.CameraIntrinsicParameters kIntr = cameraMatrices.getIntrinsics();
        
        List<Camera.CameraExtrinsicParameters> extrinsics = cameraMatrices.getExtrinsics();
        
        log.log(LEVEL, String.format("intr=\n%s\n", FormatArray.toString(kIntr.getIntrinsic(), "%.3e")));
        
        Camera.CameraExtrinsicParameters ex1;
        for (int i = 0; i < nImages; ++i) {
            ex1 = extrinsics.get(i);
            log.log(LEVEL, String.format("\n"));
            log.log(LEVEL, String.format("   r%d=\n%s\n", i, FormatArray.toString(ex1.getRotation(), "%.3e")));
            log.log(LEVEL, String.format("ansR%d=\n%s\n", i, FormatArray.toString(NeurochemistryBookData.getRotation(i), "%.3e")));
            log.log(LEVEL, String.format("   t%d=\n%s\n", i,FormatArray.toString(ex1.getTranslation(), "%.3e")));
            log.log(LEVEL, String.format("ansT%d=\n%s\n", i,FormatArray.toString(NeurochemistryBookData.getTranslation(i), "%.3e")));
        }
        
        NeurochemistryBookData.printExpectedTriangulation();
        /*
        NeurochemistryBookData.printObservedMinusProjected_Camera_Frame();
        NeurochemistryBookData.printObservedMinusProjected_Image_Frame();
        */
        
        /*
        // theta_x is rotated by 180 degrees 
        double[][] q = new double[3][];
        q[0] = new double[]{9.928e-01, 7.859e-02, -9.033e-02};
        q[1] = new double[]{7.195e-02, -9.946e-01, -7.454e-02};
        q[2] = new double[]{-9.571e-02, 6.750e-02, -9.931e-01};
        double[] qt = new double[3];
        Rotation.extractThetaFromZYX(q, qt);
        for (int i = 0; i < 3; ++i) {
            System.out.printf("qt=%.2f  ", qt[i]*180./Math.PI);
        }
        System.out.println();
        */
        /*
        // theta_z is rotataed by -180
        double[][] q = new double[3][];
        q[0] = new double[]{-9.928e-01, 7.859e-02, 9.033e-02};
        q[1] = new double[]{-7.195e-02, -9.946e-01, 7.454e-02};
        q[2] = new double[]{9.571e-02, 6.750e-02, 9.931e-01};
        double[] qt = new double[3];
        Rotation.extractThetaFromZYX(q, qt);
        for (int i = 0; i < 3; ++i) {
            System.out.printf("qt=%.2f  ", qt[i]*180./Math.PI);
        }
        System.out.println();
        */
        
        double[][] q = new double[3][];
        q[0] = new double[]{9.652e-01, -4.408e-03, 2.616e-01};
        q[1] = new double[]{3.769e-02, -9.871e-01, -1.557e-01};
        q[2] = new double[]{2.589e-01, 1.601e-01, -9.525e-01};
        double[] qt = new double[3];
        Rotation.extractThetaFromZYX(q, qt);
        for (int i = 0; i < 3; ++i) {
            System.out.printf("qt=%.2f  ", qt[i]*180./Math.PI);
        }
        System.out.println();
    }

    public void testSolveForExtrinsicPlanarWetzstein() {
        //static Camera.CameraExtrinsicParameters solveForExtrinsicPlanarWetzstein(
        //        double[][] coordsC, double[][] coordsW) throws NotConvergedException

        /* generate points:
        [xc]  =  [1  0  0  tx ] * [ r11  r12  r13  0 ] * [xw] = [ r11  r12  r13  tx ] * [xw]
        [yc]     [0  1  0  ty ]   [ r21  r22  r23  0 ]   [yw]   [ r21  r22  r23  ty ]   [yw]
        [zc]     [0  0  1  tz ]   [ r31  r32  r33  0 ]   [zw]   [ r31  r32  r33  tz ]   [zw]
        [1]      [0  0  0  1  ]   [ 0    0    0    1 ]   [1 ]   [ 0    0    0    1  ]   [1 ]

         */
    }

    protected void normalize(double[][] x) {
        for (int i = 0; i < x[0].length; ++i) {
            for (int j = 0; j < x.length; ++j) {
                x[j][i] /= x[x.length - 1][i];
            }
        }
    }
}
