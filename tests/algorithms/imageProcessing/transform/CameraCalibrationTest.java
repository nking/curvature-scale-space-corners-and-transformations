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
    
    public void testApplyRemoveRadialDistortion()
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
        
        int nFeatures = 256;
        int nImages = 5;
        
        // 3 X 256
        double[][] coordsW = Zhang98Data.getFeatureWCS();
        assertEquals(3, coordsW.length);
        assertEquals(nFeatures, coordsW[0].length);
        
        //3 X (256*5)
        double[][] coordsI = Zhang98Data.getObservedFeaturesInAllImages();
        assertEquals(3, coordsI.length);
        assertEquals(nFeatures*nImages, coordsI[0].length);

        CameraMatrices c = use3D ? CameraCalibration.estimateCamera3D(nFeatures, coordsI, coordsW, useR2R4):
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

        // TODO: add assertions back once the refinements are added back to the methods
        /*
        assertTrue(Math.abs(alphaE - alpha) < 0.1);
        assertTrue(Math.abs(gammaE - gamma) < 0.1);
        assertTrue(Math.abs(u0E - u0) < 0.1);
        assertTrue(Math.abs(betaE - beta) < 0.1);
        assertTrue(Math.abs(v0E - v0) < 0.1);
         */

        List<CameraExtrinsicParameters> extrinsics = c.getExtrinsics();
        CameraExtrinsicParameters ex1;
        for (int i = 0; i < nImages; ++i) {
            ex1 = extrinsics.get(i);
            log.log(LEVEL, String.format("\nimg %d:\n", i));
            log.log(LEVEL, String.format("   r=%s\n", FormatArray.toString(ex1.getRotation(), "%.3e")));
            log.log(LEVEL, String.format("   t=%s\n", FormatArray.toString(ex1.getTranslation(), "%.3e")));

            double[][] rExp = Zhang98Data.getRotation(i+1);
            double[] tExp = Zhang98Data.getTranslation(i+1);
            double[][] rDiff = Rotation.procrustesAlgorithmForRotation(rExp, ex1.getRotation());

            double[] tDiff = MatrixUtil.subtract(tExp, ex1.getTranslation());
            double rFS = MatrixUtil.frobeniusNorm(rDiff);
            double tS = MatrixUtil.lPSum(tDiff, 2);
            log.log(LEVEL, String.format("   tSSD=%.3f\n", tS*tS));
            log.log(LEVEL, String.format("   rSSD=%.3f\n", rFS));
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

    public void testCalibration3DRandom() throws Exception {
        log.log(LEVEL, "testCalibration3DRandom");

        /*
        generate an intrinsic camera matrix.
        generate 3 or more extrinsic camera matrices (1 each for 1 image).
        generate a projection matrix P for each image.

        generate 3D points in WCS.

        project points to images

        use estimateCamera3D to recover the intrinsic and extrinsic matrices
        */

        long seed = 12345l;//System.nanoTime();
        System.out.println("seed=" + seed);
        Random rand = new Random(seed);

        double[][] kIntr = new double[][]{
                {1380, 0, 247},
                {0, 2033, 244},
                {0, 0, 1}
        };

        //TODO: change the signs in omI and tI
        //TODO: add radial distortion
        //TODO: change camera matrix and conditions for several cases: smart phone, telephoto, etc.
        int nImages = 3 + rand.nextInt(7);
        int nPoints = 200 + rand.nextInt(200);

        List<double[]> om = new ArrayList<>();
        List<double[][]> r = new ArrayList<>();
        List<double[]> t = new ArrayList<>();
        List<double[][]> p = new ArrayList<>();

        double[] boundsWCS = new double[]{
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY
        };

        double[] radial = null;

        //double[] expectedRom = new double[]{-0.08573, -0.99438, 0.0621};
        //double[][] expectedRot = Rotation.createRotationFromUnitLengthAngleAxis(expectedRom, 47.7*Math.PI/180.);
        //double[] expectedT = new double[]{-211.28, -106.06, 1583.75};

        for (int i = 0; i < nImages; ++i) {
            double e = 1E-4;

            double[] omI = new double[]{rand.nextDouble(), rand.nextDouble(), rand.nextDouble()};
            omI = MatrixUtil.normalizeLP(omI, 2);
            double[][] rI = Rotation.createRotationFromUnitLengthAngleAxis(omI, 5*rand.nextDouble() * Math.PI/180.);
            double[] tI = new double[]{-rand.nextInt((int)Math.round(kIntr[0][2])),
                    -rand.nextInt((int)Math.round(kIntr[1][2])), rand.nextInt((int)Math.round(kIntr[1][1]))};

            double[][] pI = new double[3][4];
            for (int j = 0; j < 3; ++j) {
                System.arraycopy(rI[j], 0, pI[j], 0, 3);
                pI[j][3] = tI[j];
            }
            pI = MatrixUtil.multiply(kIntr, pI);

            Camera.CameraPoseParameters pose = CameraPose.calculatePoseFromP(pI);

            // assert that we recover the camera matrices:
            double[][] diff = MatrixUtil.pointwiseSubtract(kIntr, pose.getIntrinsicParameters().getIntrinsic());
            double fs = MatrixUtil.frobeniusNorm(diff) / Math.sqrt(5);
            assertTrue(fs < 1E-7);
            diff = MatrixUtil.pointwiseSubtract(rI, pose.getExtrinsicParameters().getRotation());
            fs = MatrixUtil.frobeniusNorm(diff) / Math.sqrt(9);
            assertTrue(fs < 1E-7);
            fs = MatrixUtil.lPSum(MatrixUtil.subtract(tI, pose.getExtrinsicParameters().getTranslation()), 2);
            assertTrue(fs < 1E-7);

            // TODO: make the coord generation simpler.
            // this is a quick look at using pInv to de-project an image point into WCS.
            // looking at the bounds for the 3D points over all nImages to try to limit the range of 3D to create
            // only those within image frame of all images.
            // for the nImages random rotations and translations above, the intersection of 3D bounds over all images
            // is 0, so to use this better, would need to improve the rotation and translation combinations.
            // meanwhile, boundsWCS is still used, but not as intended.
            double[][] pIInv = MatrixUtil.pseudoinverseFullRowRank(pI);
            // 4x3 * 3xn = 4xn
            double[][] bWCS = MatrixUtil.multiply(pIInv, new double[][]{
                    {0, 2*kIntr[0][2]}, {0, 2*kIntr[1][2]}, {1, 1}
            });

            log.log(LEVEL, String.format("bWCS%d\n%s", i, FormatArray.toString(bWCS, "%10.2e")));
            //double[][] _bIm = MatrixUtil.multiply(pI, bWCS);
            //normalize(_bIm, 2);
            for (int j = 0; j < 4; ++j) {
                double min = Math.min(bWCS[j][0], bWCS[j][1]);
                double max = Math.max(bWCS[j][0], bWCS[j][1]);
                boundsWCS[2 * j] = Math.max(boundsWCS[2 * j], min);
                boundsWCS[2 * j + 1] = Math.min(boundsWCS[2 * j + 1], max);
            }

            om.add(omI);
            r.add(rI);
            t.add(tI);
            p.add(pI);
        }
        log.log(LEVEL, String.format("bounds=[%s]\n", FormatArray.toString(boundsWCS, "%.2e")));
        for (int i=0; i < boundsWCS.length/2; ++i) {
            log.log(LEVEL, String.format("boundsWCS %b\n", boundsWCS[2*i] < boundsWCS[2*i + 1]));
        }

        // roughly, find the largest lower left and smallest upper right corner in boundsWCS

        double[][] XW = generateWCS(rand, nPoints, boundsWCS);

        double[][] coordsI = new double[3][nPoints * nImages];
        List<double[][]> x = new ArrayList<>();
        for (int i = 0; i < nImages; ++i) {
            double[][] xI = MatrixUtil.multiply(p.get(i), XW);
            normalize(xI, 2);
            x.add(xI);
            for (int j = 0; j < xI[0].length; ++j) {
                if (xI[0][i] < 0 || xI[1][i] < 0) {
                    int _t=1;
                }
            }
            for (int j = 0; j < xI.length; ++j) {
                System.arraycopy(xI[j], 0, coordsI[j], i * nPoints, nPoints);
            }
        }

        normalize(XW, 3);
        double[][] coordsW = MatrixUtil.copySubMatrix(XW, 0, 2, 0, XW[0].length-1);
        //normalize(coordsW, 2); need this dimension in forming P


        boolean useR2R4 = false;
        CameraMatrices c = CameraCalibration.estimateCamera3D(nPoints, coordsI, coordsW, useR2R4);
        double[][] resultK = c.getIntrinsics().getIntrinsic();

        double[][] kDiff = MatrixUtil.pointwiseSubtract(kIntr, resultK);
        double kFS = MatrixUtil.frobeniusNorm(kDiff);
        log.log(LEVEL, FormatArray.toString(kDiff, "%10.4e"));
        log.log(LEVEL, String.format("   kSSD=%.3f\n", kFS*kFS));

        for (int i = 0; i < nImages; ++i) {
            CameraExtrinsicParameters extr = c.getExtrinsics().get(i);
            double[][] resultR = extr.getRotation();
            double[] resultT = extr.getTranslation();

            double[][] rDiff = Rotation.procrustesAlgorithmForRotation(r.get(i), resultR);
            double[] tDiff = MatrixUtil.subtract(t.get(i), resultT);
            double rFS = MatrixUtil.frobeniusNorm(rDiff);
            double tS = MatrixUtil.lPSum(tDiff, 2);
            log.log(LEVEL, String.format("   tSSD=%.3f\n", tS*tS));
            log.log(LEVEL, String.format("   rSSD=%.3f\n", rFS));
        }
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
}
