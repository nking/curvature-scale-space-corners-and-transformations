package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.features.RANSACEuclideanSolver;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.imageProcessing.matching.CorrespondenceMaker;
import algorithms.imageProcessing.matching.CorrespondenceMaker.CorrespondenceList;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.matching.ORBMatcher;
import static algorithms.imageProcessing.transform.Rotation.extractThetaFromZYX;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.*;

import java.io.IOException;
import java.text.Normalizer;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertNotNull;
import no.uib.cipr.matrix.NotConvergedException;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class ReconstructionTest extends TestCase {
    
    public ReconstructionTest() {
    }

    public void testProjectiveWithZhang() throws IOException, NotConvergedException {

        int idx1 = 1;
        int idx2 = 2;

        double[][] x1 = Zhang98Data.getObservedFeaturesInImage(idx1);
        double[][] x2 = Zhang98Data.getObservedFeaturesInImage(idx2);
        double[][] intr = Zhang98Data.getIntrinsicCameraMatrix();
        double[] radial = Zhang98Data.getRadialDistortionR2R4();
        double[][] r1 = Zhang98Data.getRotation(idx1);
        double[] t1 = Zhang98Data.getTranslation(idx1);
        double[][] r2 = Zhang98Data.getRotation(idx2);
        double[] t2 = Zhang98Data.getTranslation(idx2);
        double[][] camera1 = Camera.createCamera(intr, r1, t1);
        double[][] camera2 = Camera.createCamera(intr, r2, t2);
        double[][] x1c = Camera.pixelToCameraCoordinates(x1, intr, radial, true);
        double[][] x2c = Camera.pixelToCameraCoordinates(x2, intr, radial, true);
        double[][] XW = Zhang98Data.getFeatureWCS();

        double[] r1Vec = Rotation.extractRodriguesRotationVector(r1);
        double[] r2Vec = Rotation.extractRodriguesRotationVector(r2);
        double[][] r1Transposed = MatrixUtil.transpose(r1);
        double[][] r1I = MatrixUtil.multiply(r1Transposed, r1);
        double[][] r2RelToR1 = MatrixUtil.multiply(r1Transposed, r2);
        System.out.printf("t1=%s\n", FormatArray.toString(t1, "%.3e"));
        System.out.printf("t2=%s\n", FormatArray.toString(t2, "%.3e"));
        System.out.printf("diff t=%s\n", FormatArray.toString(MatrixUtil.subtract(t2, t1), "%.3e"));
        System.out.printf("r1=%s\n", FormatArray.toString(r1, "%.3e"));
        System.out.printf("r2=%s\n", FormatArray.toString(r2, "%.3e"));
        System.out.printf("r2RelToR1=%s\n", FormatArray.toString(r2RelToR1, "%.3e"));

        int n = x1[0].length;

        double[][] x1Pt = MatrixUtil.zeros(3, 1);
        double[][] x2Pt = MatrixUtil.zeros(3, 1);
        int ii, j, i;
        for (ii = 0; ii < n; ++ii) {
            for (j = 0; j < 3; ++j) {
                x1Pt[j][0] = x1c[j][ii];
                x2Pt[j][0] = x2c[j][ii];
            }
            // better answer produced by using camera coords xc, yc and extrinsic projection matrix P=[R|t]

            Triangulation.WCSPt wcsPt = Triangulation.calculateWCSPoint(camera1, camera2, x1Pt, x2Pt);
            //Triangulation.WCSPt wcsPt = Triangulation.calculateWCSPoint(intr, r1, t1, intr, r2, t2, x1Pt, x2Pt);
            if (wcsPt == null) {
                continue;
            }
            MatrixUtil.multiply(wcsPt.X, 1. / wcsPt.X[3]);
            System.out.printf("WCS[%d]=%s,\t alpha=%.3e\n",
                    ii, FormatArray.toString(wcsPt.X, "%.3e"),
                    wcsPt.alpha);
            System.out.printf("expected=%s\n",
                    FormatArray.toString(MatrixUtil.extractColumn(XW, ii), "%.3e"));
        }

        Reconstruction.ProjectionResults[] pr = Reconstruction.calculateProjectiveReconstruction(x1c, x2c);
        for (ii = 0; ii < pr.length; ++ii) {
            double[][] p = pr[ii].projectionMatrices;
            System.out.printf("pProj=\n%s\n", FormatArray.toString(p, "%.3e"));
            System.out.printf("r2RelToR1=%s\n", FormatArray.toString(r2RelToR1, "%.3e"));
            double[][] pXW = pr[ii].XW;
            System.out.printf("pXW[%d]=%s\n",
                    ii, FormatArray.toString(MatrixUtil.transpose(pXW), "%.3e"));
        }
    }

        public void estProjectiveWithBouguetIm2LeftRight() throws IOException, NotConvergedException {

        //test data from:
        //http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example5.html
        // now at http://robots.stanford.edu/cs223b04/JeanYvesCalib/
        // "Fifth calibration example - Calibrating a stereo system, stereo image rectification and 3D stereo triangulation"
        // by Jean-Yves Bouguet
        // Camera Calibration Toolbox for Matlab
        //    saved as pdf in miscNotes/bouguetj_5th_calibration_example.pdf
        // see testresources/bouguet_stereo
        //
        // this is left02.jpg and right02.jpg
        //
        //
        // note: the checkerboard squares are 30mm in WCS metric
        //
        // note: radial distortion should be corrected:  use on the original coordinates:
        //     x_corrected = x*(1 + k1*r^2 + k2r^4) where r is distance of point from cc.

        double[][] k1Intr = Camera.createIntrinsicCameraMatrix(533.5, 341.6, 235.2);
        double[][] k2Intr = Camera.createIntrinsicCameraMatrix(536.6, 326.3, 250.1);
        double[][] K1IntrInv = Camera.createIntrinsicCameraMatrixInverse(k1Intr);
        double[][] K2IntrInv = Camera.createIntrinsicCameraMatrixInverse(k2Intr);

        boolean useR2R4 = true;
        double[] radial1 = new double[]{-0.288, 0.097, 0.001};
        double[] radial2 = new double[]{-0.289,0.107, -0.001};

        // ===================================================================================
        /* for the 2nd image pair:
        left extrinsic matrix:
        [[-9.48402795e-02  9.75707586e-01 -1.97484246e-01 -4.67064408e+01]
         [ 7.52568521e-01  2.00130415e-01  6.27366272e-01 -8.09188900e+01]
         [ 6.51648635e-01 -8.91208345e-02 -7.53267239e-01  2.66643941e+02]
         [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  1.00000000e+00]]

         right extrinsic matrix:
         [[-9.03220351e-02  9.76271959e-01 -1.96812071e-01 -1.45613832e+02]
         [ 7.48996521e-01  1.96836697e-01  6.32660673e-01 -8.15340647e+01]
         [ 6.56388713e-01 -9.02683572e-02 -7.49002992e-01  2.66971558e+02]
         [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  1.00000000e+00]]
         */
        double[][] camera1RT = new double[3][];
        camera1RT[0] = new double[]{9.48402795e-02,  9.75707586e-01, -1.97484246e-01, -4.67064408e+01};
        camera1RT[1] = new double[]{7.52568521e-01,  2.00130415e-01,  6.27366272e-01, -8.09188900e+01};
        camera1RT[2] = new double[]{6.51648635e-01, -8.91208345e-02, -7.53267239e-01,  2.66643941e+02};
        double[][] camera2RT = new double[3][];
        camera2RT[0] = new double[]{-9.03220351e-02,  9.76271959e-01, -1.96812071e-01, -1.45613832e+02};
        camera2RT[1] = new double[]{7.48996521e-01,  1.96836697e-01,  6.32660673e-01, -8.15340647e+01};
        camera2RT[2] = new double[]{6.56388713e-01, -9.02683572e-02, -7.49002992e-01,  2.66971558e+02};

        double[][] bouguetR1 = MatrixUtil.copySubMatrix(camera1RT, 0, 2, 0,2);
        double[][] bouguetR2 = MatrixUtil.copySubMatrix(camera2RT, 0, 2, 0,2);
        double[][] bouguetR1Orth = Rotation.orthonormalizeUsingSVD(bouguetR1);
        double[][] bouguetR2Orth = Rotation.orthonormalizeUsingSVD(bouguetR2);
        double[] bougetT1 = MatrixUtil.extractColumn(camera1RT, 3);
        double[] bougetT2 = MatrixUtil.extractColumn(camera2RT, 3);
        double[][] bouguetR1Transposed = MatrixUtil.transpose(bouguetR1);
        double[][] bouguetR1I = MatrixUtil.multiply(bouguetR1Transposed, bouguetR1);
        double[][] bouguetR2RelToR1 = MatrixUtil.multiply(bouguetR1Transposed, bouguetR2);
        double[] bouguetT12 = MatrixUtil.subtract(bougetT2, bougetT1);
        double[] bouguetR1Vec = Rotation.extractRodriguesRotationVector(bouguetR1);
        double[] bouguetR2Vec = Rotation.extractRodriguesRotationVector(bouguetR2);
        double[] bouguetR2RelToR1Vec = Rotation.extractRodriguesRotationVector(bouguetR2RelToR1);
        double[] bouguetR1IVec = Rotation.extractRodriguesRotationVector(bouguetR1I);
        double[][] bouguetRot12 = Rotation.rotationBetweenTwoDirections0(bouguetR1Vec, bouguetR2Vec);
        double[] bouguetRVec12 = Rotation.extractRodriguesRotationVector(bouguetRot12);

        double[][] chk1 = MatrixUtil.multiply(bouguetRot12, bouguetR1); // same as R2?
        double[][] chk2 = MatrixUtil.multiply(bouguetRot12, bouguetR2); // same as R1?
        double[][] chk3 = MatrixUtil.multiply(MatrixUtil.transpose(bouguetRot12), bouguetR2); // same as R1?

        //expected om = 0.00611, 0.00489, -0.00359   +- 0.0027, 0.00308, 0.00029
        //expected t = -99.84929, 0.82221, 0.43647   +- 0.142, 0.11352, 0.49773

        double[][] x1 = MatrixUtil.zeros(3, 13);
        double[][] x2 = MatrixUtil.zeros(3, 13);
        Arrays.fill(x1[2], 1);
        Arrays.fill(x2[2], 1);
        x1[0][0]=5.8400e+02;  x1[1][0]=1.4800e+02;  x2[0][0]=3.8000e+02;  x2[1][0]=1.5200e+02;
        x1[0][1]=3.9600e+02;  x1[1][1]=3.9800e+02;  x2[0][1]=2.6200e+02;  x2[1][1]=4.1000e+02;
        x1[0][2]=5.9000e+02;  x1[1][2]=1.4700e+02;  x2[0][2]=3.8500e+02;  x2[1][2]=1.5400e+02;
        x1[0][3]=2.5400e+02;  x1[1][3]=3.5400e+02;  x2[0][3]=1.2400e+02;  x2[1][3]=3.6400e+02;
        x1[0][4]=5.4600e+02;  x1[1][4]=2.3200e+02;  x2[0][4]=3.6700e+02;  x2[1][4]=2.4600e+02;
        x1[0][4]=4.9000e+02;  x1[1][5]=9.6000e+01;  x2[0][5]=2.7000e+02;  x2[1][5]=1.0400e+02;
        x1[0][6]=4.8000e+02;  x1[1][6]=3.8000e+02;  x2[0][6]=3.3600e+02;  x2[1][7]=3.9600e+02;
        x1[0][7]=2.5300e+02;  x1[1][7]=7.9000e+01;  x2[0][7]=6.5000e+01;  x2[1][7]=1.0400e+02;
        x1[0][8]=3.6100e+02;  x1[1][8]=3.8300e+02;  x2[0][8]=2.2800e+02;  x2[1][8]=3.9600e+02;
        x1[0][9]=4.0200e+02;  x1[1][9]=3.8800e+02;  x2[0][9]=2.6000e+02;  x2[1][9]=4.0000e+02;
        x1[0][10]=2.0800e+02;  x1[1][10]=1.6400e+02;  x2[0][10]=4.8000e+01;  x2[1][10]=1.8400e+02;
        x1[0][11]=4.0400e+02;  x1[1][11]=3.8800e+02;  x2[0][11]=2.6400e+02;  x2[1][11]=4.0400e+02;
        x1[0][12]=4.4800e+02;  x1[1][12]=2.4800e+02;  x2[0][12]=2.6400e+02;  x2[1][12]=2.6400e+02;
        double[][] x1c = Camera.pixelToCameraCoordinates(x1, k1Intr, null /*radial1*/, useR2R4);
        double[][] x2c = Camera.pixelToCameraCoordinates(x2, k2Intr, null /*radial2*/, useR2R4);

        System.out.printf("K1Inv = \n%s\n",
                FormatArray.toString(Camera.createIntrinsicCameraMatrixInverse(k1Intr), "%.3e"));
        System.out.printf("K2Inv = \n%s\n",
                FormatArray.toString(Camera.createIntrinsicCameraMatrixInverse(k2Intr), "%.3e"));

        double[][] camera1R = MatrixUtil.copySubMatrix(camera1RT, 0, 2, 0, 2);
        double[][] camera2R = MatrixUtil.copySubMatrix(camera2RT, 0, 2, 0, 2);

        double[][] leftInv = MatrixUtil.transpose(camera1R);
        // apply inv(left) to left to see that the result is the identity matrix
        double[][] leftTleft = MatrixUtil.createATransposedTimesA(camera1R);
        // then apply inv(left) to right to be able to compare the result with the identity matrix,
        //    that is, the rotation of right with respect to left.
        double[][] leftTright = MatrixUtil.multiply(leftInv, camera2R);
        // also, the rotation that would take the right camera to the orientation of the left camera is then
        // the inverse of that = inv*(inv(left)*right) = inv(right)*left
        double[][] rightToLeft = MatrixUtil.multiply(MatrixUtil.transpose(camera2R), camera1R);
        System.out.printf("rightToLeft=\n%s\n", FormatArray.toString(rightToLeft, "%.3e"));

        double[] thetasLeftToRight = Rotation.extractThetaFromZYX(leftTright);
        double[] thetasRightToLeft = Rotation.extractThetaFromZYX(rightToLeft);
        double[] rVecRodr = Rotation.extractRodriguesRotationVector(leftTright);
        double[] rAxis = Rotation.extractRotationAxisFromZXY(leftTright);
        System.out.printf("\nthetas (R to L)=%s\n", FormatArray.toString(thetasRightToLeft, "%.3e"));
        System.out.printf("thetas (L to R)=%s\n", FormatArray.toString(thetasLeftToRight, "%.3e"));
        System.out.printf("rVecRodr=%s\n", FormatArray.toString(rVecRodr, "%.3e"));
        System.out.printf("rAxis=%s\n\n", FormatArray.toString(rAxis, "%.3e"));

        double[][] k2ExtrRot = leftTright;
        double[] k2ExtrTrans = new double[]{-4.67064408e+01- -1.45613832e+02, -8.09188900e+01 - -8.15340647e+01,
                2.66643941e+02 - 2.66971558e+02};

        double[] thetas = thetasLeftToRight;
        double[] thetaDegrees = Arrays.copyOf(thetas, thetas.length);
        MatrixUtil.multiply(thetaDegrees, 180./Math.PI);

        Camera.CameraExtrinsicParameters extrCalc = Reconstruction.calculateProjectiveMotion(x1c, x2c);

        double[] ti = Arrays.copyOf(extrCalc.getTranslation(), 3);
        //xim = K * xc.  xc = K^-1*xim
        //double[][] k2IntrInv = Camera.createIntrinsicCameraMatrixInverse(k2Intr);
        double[] ti0 = MatrixUtil.multiplyMatrixByColumnVector(k2Intr, ti);
        MatrixUtil.multiply(ti, 1./ti[2]);
        double[] ti00 = MatrixUtil.multiplyMatrixByColumnVector(k2Intr, ti);
        System.out.printf("transformed translation=\n%s\n    =%s\n", FormatArray.toString(ti0, "%.3e"),
                FormatArray.toString(ti00, "%.3e"));

        assertNotNull(extrCalc);
        assertNotNull(extrCalc.getRotation());
        double[][] extrCalcRot = extrCalc.getRotation();
        double[] extrCalcThetas = new double[3];
        Rotation.extractThetaFromZYX(extrCalcRot, extrCalcThetas);
        double[] extrCalcThetaDegrees = Arrays.copyOf(extrCalcThetas, extrCalcThetas.length);
        MatrixUtil.multiply(extrCalcThetaDegrees, 180./Math.PI);

        System.out.printf("1: derived rot=\n%s\ndegrees=\n%s\ntrans=\n%s\n",
                FormatArray.toString(extrCalcRot, "%.3e"),
                FormatArray.toString(extrCalcThetaDegrees, "%.3e"),
                FormatArray.toString(extrCalc.getTranslation(), "%.3e"));

        System.out.printf("1: expected rot=\n%s\nexpected degrees=\n%s\nexpected trans=\n%s\n",
                FormatArray.toString(k2ExtrRot, "%.3e"),
                FormatArray.toString(thetaDegrees, "%.3e"),
                FormatArray.toString(k2ExtrTrans, "%.3e"));

        System.out.printf("the derived rotation agrees with python's openCV cv2.recoverPose()\n");

        // test triangulation, given the values that Bouguet found for R and T
        System.out.printf("\ngiven Bouguet P1, P2, check triangulation:\n");
        // the triangulation results agree with the figure displaying the cameras and the checkerboards in WCS
        int n = x1[0].length;
        int j, ii;
        /*
        x_c = R * xw + t
            = [R | t] xw
        x_im = K*x_c = K*[R|t]*xw
         */
        double[][] x1Pt = MatrixUtil.zeros(3, 1);
        double[][] x2Pt = MatrixUtil.zeros(3, 1);
        for (ii = 0; ii < n; ++ii) {
            for (j = 0; j < 3; ++j) {
                x1Pt[j][0] = x1c[j][ii];
                x2Pt[j][0] = x2c[j][ii];
            }
            // better answer produced by using camera coords xc, yc and extrinsic projection matrix P=[R|t]

            Triangulation.WCSPt wcsPt = Triangulation.calculateWCSPoint(camera1RT, camera2RT, x1Pt, x2Pt);
            if (wcsPt == null) {
                continue;
            }
            MatrixUtil.multiply(wcsPt.X, 1. / wcsPt.X[3]);
            System.out.printf("xc=(%s, %s) ==> WCS[%d]=%s,\t alpha=%.3e\n",
                FormatArray.toString(MatrixUtil.transpose(x1Pt), "%.3e"),
                FormatArray.toString(MatrixUtil.transpose(x2Pt), "%.3e"),
                ii, FormatArray.toString(wcsPt.X, "%.3e"),
               wcsPt.alpha);
        }

        // the MASKS method:
        Reconstruction.ProjectionResults[] results = Reconstruction.calculateProjectiveReconstruction(x1c, x2c);
        for (Reconstruction.ProjectionResults result : results) {
            System.out.printf("projectionMatrix: \n%s\n", FormatArray.toString(result.projectionMatrices, "%.3e"));
            System.out.printf("PR XW=\n%s\n", FormatArray.toString(MatrixUtil.transpose(result.XW), "%.3e"));
        }

        // needs x1 and x2 in image coordinates:
        Reconstruction.ReconstructionResults results2 = Reconstruction.calculateUsingEssentialMatrix(
            k1Intr, k2Intr, x1, x2);
        System.out.printf("\nMethod calculateUsingEssentialMatrix(): \nRot2=\n%s\nTrans2=\n%s\nEM XW=\n%s\n",
                FormatArray.toString(results2.k2ExtrRot, "%.3e"),
                FormatArray.toString(results2.k2ExtrTrans, "%.3e"),
                FormatArray.toString(MatrixUtil.transpose(results2.XW), "%.3e"));

    }

    public void _testProjectiveWithBouguetLeftIm1Im2() throws IOException, NotConvergedException {

        System.out.printf("testProjectiveWithBouguetLeftIm1Im2\n");

        //test data from:
        //http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example5.html
        // now at http://robots.stanford.edu/cs223b04/JeanYvesCalib/
        // "Fifth calibration example - Calibrating a stereo system, stereo image rectification and 3D stereo triangulation"
        // by Jean-Yves Bouguet
        // Camera Calibration Toolbox for Matlab
        //    saved as pdf in miscNotes/bouguetj_5th_calibration_example.pdf
        // see testresources/bouguet_stereo
        //
        // this is left01.jpg and left02.jpg
        //
        // note: the checkerboard squares are 30mm in WCS metric
        //
        // note: radial distortion should be corrected:  use on the original coordinates:
        //     x_corrected = x*(1 + k1*r^2 + k2r^4) where r is distance of point from cc.

        double[][] k1Intr = Camera.createIntrinsicCameraMatrix(533.5, 341.6, 235.2);
        boolean useR2R4 = true;
        double[] radial1 = new double[]{-0.288, 0.097, 0.001};

        double[][] cameraIm2RT = new double[3][];
        cameraIm2RT[0] = new double[]{9.48402795e-02,  9.75707586e-01, -1.97484246e-01, -4.67064408e+01};
        cameraIm2RT[1] = new double[]{7.52568521e-01,  2.00130415e-01,  6.27366272e-01, -8.09188900e+01};
        cameraIm2RT[2] = new double[]{6.51648635e-01, -8.91208345e-02, -7.53267239e-01,  2.66643941e+02};
        double[][] cameraIm1RT = new double[3][];
        cameraIm1RT[0] = new double[]{9.91272809e-03,  9.62212583e-01, -2.72118876e-01, -8.96506577e+01};
        cameraIm1RT[1] = new double[]{9.86085734e-01,  3.57539594e-02,  1.62347096e-01, -1.29537145e+02};
        cameraIm1RT[2] = new double[]{1.65941746e-01, -2.69941844e-01, -9.48469682e-01,  4.77077898e+02};
        double[][] cameraIm1KRT = MatrixUtil.multiply(k1Intr, cameraIm1RT);
        double[][] cameraIm2KRT = MatrixUtil.multiply(k1Intr, cameraIm2RT);

        double[][] x1, x2;
        // =======
        x1 = MatrixUtil.zeros(3, 9);
        x2 = MatrixUtil.zeros(3, 9);
        x1[0][0]=243;  x1[1][0]=65;
        x2[0][0]=224;  x2[1][0]=350;
        x1[0][1]=512;  x1[1][1]=84;
        x2[0][1]=252;  x2[1][1]=80;
        x1[0][2]=279;  x1[1][2]=225;
        x2[0][2]=412;  x2[1][2]=366;
        x1[0][3]=476;  x1[1][3]=296;
        x2[0][3]=570;  x2[1][3]=190;
        x1[0][4]=340;  x1[1][4]=124;
        x2[0][4]=294;  x2[1][4]=288;
        x1[0][5]=444;  x1[1][5]=158;
        x2[0][5]=351;  x2[1][5]=191;
        x1[0][6]=373;  x1[1][6]=189;
        x2[0][6]=392;  x2[1][6]=280;
        x1[0][7]=408;  x1[1][7]=56;
        x2[0][7]=208;  x2[1][7]=203;
        x1[0][8]=308;  x1[1][8]=284;
        x2[0][8]=492;  x2[1][8]=360;
        double[][] x1c = Camera.pixelToCameraCoordinates(x1, k1Intr, null /*radial1*/, useR2R4);
        double[][] x2c = Camera.pixelToCameraCoordinates(x2, k1Intr, null /*radial2*/, useR2R4);

        Camera.CameraExtrinsicParameters extrCalc = Reconstruction.calculateProjectiveMotion(x1c, x2c);

        assertNotNull(extrCalc);
        assertNotNull(extrCalc.getRotation());
        double[][] extrCalcRot = extrCalc.getRotation();

        // test triangulation, given the values that Bouguet found for R and T
        System.out.printf("\ngiven Bouguet P1, P2, check triangulation:\n");
        // the triangulation results agree with the figure displaying the cameras and the checkerboards in WCS
        int n = x1[0].length;
        int j, ii;
        /*
        x_c = R * xw + t
            = [R | t] xw
        x_im = K*x_c = K*[R|t]*xw
         */
        double[][] x1Pt = MatrixUtil.zeros(3, 1);
        double[][] x2Pt = MatrixUtil.zeros(3, 1);
        for (ii = 0; ii < n; ++ii) {
            for (j = 0; j < 3; ++j) {
                x1Pt[j][0] = x1c[j][ii];
                x2Pt[j][0] = x2c[j][ii];
            }
            // this assumes the camera is stationary and the WCS object moves.
            // haven't worked through the math yet to see if can just reverse the arguments like this:
            Triangulation.WCSPt wcsPt = Triangulation.calculateWCSPoint(cameraIm1RT, cameraIm2RT, x1Pt, x2Pt);
            if (wcsPt == null) {
                continue;
            }
            MatrixUtil.multiply(wcsPt.X, 1. / wcsPt.X[3]);
            System.out.printf("WCS[%d]=%s,\t alpha=%.3e\n", ii, FormatArray.toString(wcsPt.X, "%.3e"),
                    wcsPt.alpha);
        }
    }

    public void estProjectiveWithBouguet() throws IOException, NotConvergedException {

        String sep = System.getProperty("file.separator");
        String path = ResourceFinder.findTestResourcesDirectory() + sep + "bouguet_stereo" + sep;
        String filePath1 = path + "left01_masked.png";
        String filePath2 = path + "right01_masked.png";

        int nCorners = 500;
        boolean debug = false;
        boolean useToleranceAsStatFactor = true;
        boolean recalcIterations = false;// possibly faster if set to true
        double tol = 0.5;
        ErrorType errorType = ErrorType.SAMPSONS;

        CorrespondenceList corres = CorrespondenceMaker.findUsingORB(filePath1, filePath2, nCorners, useToleranceAsStatFactor,
                recalcIterations, tol, errorType, debug);
        assertNotNull(corres);

        /*
        int i, j;
        PairIntArray x1E = new PairIntArray();
        PairIntArray x2E = new PairIntArray();
        for (i = 0; i < corres.x1[0].length; ++i) {
            x1E.add(new PairInt(corres.x1[0][i], corres.x1[1][i]));
            x2E.add(new PairInt(corres.x2[0][i], corres.x2[1][i]));
        }
        PairIntArray x1EOut = new PairIntArray();
        PairIntArray x2EOut = new PairIntArray();
        RANSACEuclideanSolver solver = new RANSACEuclideanSolver();
        EuclideanTransformationFit fitE = solver.calculateEuclideanTransformation(x1E, x2E, x1EOut, x2EOut, 1E-6);
        assertNotNull(fitE);

        // plot
        ImageExt image1 = ImageIOHelper.readImageExt(filePath1);
        ImageExt image2 = ImageIOHelper.readImageExt(filePath2);
        CorrespondencePlotter plotter = new CorrespondencePlotter(image1, image2);
        int nM = x1EOut.getN();
        for (i = 0; i < nM; ++i) {
            plotter.drawLineInAlternatingColors(x1EOut.getX(i), x1EOut.getY(i),
                    x2EOut.getX(i), x2EOut.getY(i), 1);
        }
        String ts = (new SimpleDateFormat("yyyyddMMHHmmss")).format(new Date());
        plotter.writeImage("_" + ts + "_corres");
        */

        double[][] x1 = corres.x1;
        double[][] x2 = corres.x2;

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
        double[][] left1 = new double[3][];
        left1[0] = new double[]{9.91272809e-03,  9.62212583e-01, -2.72118876e-01};
        left1[1] = new double[]{9.86085734e-01,  3.57539594e-02,  1.62347096e-01};
        left1[2] = new double[]{1.65941746e-01, -2.69941844e-01, -9.48469682e-01};
        double[][] right1 = new double[3][];
        right1[0] = new double[]{1.35496246e-02,  9.61601043e-01, -2.74116476e-01};
        right1[1] = new double[]{9.85091244e-01,  3.41816887e-02,  1.68602649e-01};
        right1[2] = new double[]{1.71498247e-01, -2.72314243e-01, -9.46801618e-01};

        // inverse of the rotation matrix is the transpose
        double[][] leftInv = MatrixUtil.transpose(left1);
        // apply inv(left) to left to see that the result is the identity matrix
        double[][] leftTleft = MatrixUtil.createATransposedTimesA(left1);
        // then apply inv(left) to right to be able to compare the result with the identity matrix,
        //    that is, the rotation of right with respect to left.
        double[][] leftTright = MatrixUtil.multiply(leftInv, right1);
        // also, the rotation that would take the right camera to the orientation of the left camera is then
        // the inverse of that = inv*(inv(left)*right) = inv(right)*left
        double[] thetasLeftToRight = Rotation.extractThetaFromZYX(leftTright);
        //thetas000  is 0.002, 0.006, 0.002.   note that (0.002-0.007) radians is 0.29 degrees.
        // expecting om=0.00669, 0.00452, -0.0035
        double[][] rightToLeft = MatrixUtil.multiply(MatrixUtil.transpose(right1), left1);
        System.out.printf("rightToLeft=\n%s\n", FormatArray.toString(rightToLeft, "%.3e"));

        double[] thetasRightToLeft = Rotation.extractThetaFromZYX(rightToLeft);
        double[] rVecRodr = Rotation.extractRodriguesRotationVector(leftTright);
        double[] rAxis = Rotation.extractRotationAxisFromZXY(leftTright);
        System.out.printf("\nthetas (R to L)=%s\n", FormatArray.toString(thetasRightToLeft, "%.3e"));
        System.out.printf("thetas (L to R)=%s\n", FormatArray.toString(thetasLeftToRight, "%.3e"));
        System.out.printf("\nthetas=%s\n", FormatArray.toString(thetasLeftToRight, "%.3e"));
        System.out.printf("rVecRodr=%s\n", FormatArray.toString(rVecRodr, "%.3e"));
        System.out.printf("rAxis=%s\n", FormatArray.toString(rAxis, "%.3e"));
        System.out.printf("expect=%s\n\n", FormatArray.toString(new double[]{0.00669, 0.00452, -0.0035}, "%.3e"));


        double[] rodriguesL1 = Rotation.extractRodriguesRotationVector(left1);
        double[] rodriguesL2 = Rotation.extractRodriguesRotationVector(right1);
        // diffRodrigues ix [-0.0014421799439685579, -0.006811388962245868, 0.003946526358523994]
        double[] diffRodrigues = MatrixUtil.subtract(rodriguesL1, rodriguesL2);
        double[][] r120 = Rotation.rotationBetweenTwoDirections0(rodriguesL1, rodriguesL2);
        double[][] r121 = Rotation.rotationBetweenTwoDirections1(rodriguesL1, rodriguesL2);

        //what rotation
        double[][] diffRSameCenter = Rotation.procrustesAlgorithmForRotation(left1, right1);
        // expecting: 0.00669, 0.00452, -0.0035.  not a bad approximation for the first image set.

        double[][] k1Intr = Camera.createIntrinsicCameraMatrix(533.5, 341.6, 235.2);
        double[][] k2Intr = Camera.createIntrinsicCameraMatrix(536.6, 326.3, 250.1);

        double[][] k1ExtrRot = MatrixUtil.createIdentityMatrix(3);
        double[] k1ExtrTrans = new double[]{0, 0, 0};
        double[][] k2ExtrRot = Rotation.createRodriguesFormulaRotationMatrix(
                new double[]{0.00611, 0.00409, -0.00359});
        double[] k2ExtrTrans = new double[]{-99.85, 0.82, 0.44};

        double[] thetas = new double[3];
        Rotation.extractThetaFromZYX(k2ExtrRot, thetas);
        double[] thetaDegrees = Arrays.copyOf(thetas, thetas.length);
        MatrixUtil.multiply(thetaDegrees, 180./Math.PI);

        boolean useR2R4 = true;
        double[] radial1 = new double[]{-0.288, 0.097, 0.001};
        double[] radial2 = new double[]{-0.289,0.107, -0.001};

        double[][] x1c = Camera.pixelToCameraCoordinates(x1, k1Intr, radial1, useR2R4);
        double[][] x2c = Camera.pixelToCameraCoordinates(x1, k1Intr, radial1, useR2R4);

        // calculate the rotation and translation from the correspondence
        Camera.CameraExtrinsicParameters extrCalc = Reconstruction.calculateProjectiveMotion(x1, x2);
        assertNotNull(extrCalc);
        assertNotNull(extrCalc.getRotation());
        double[][] extrCalcRot = extrCalc.getRotation();
        double[] extrCalcThetas = new double[3];
        Rotation.extractThetaFromZYX(extrCalcRot, extrCalcThetas);
        double[] extrCalcThetaDegrees = Arrays.copyOf(extrCalcThetas, extrCalcThetas.length);
        MatrixUtil.multiply(extrCalcThetaDegrees, 180./Math.PI);

        System.out.printf("0: derived rot=\n%s\ndegrees=\n%s\ntrans=\n%s\n",
                FormatArray.toString(extrCalcRot, "%.3e"),
                FormatArray.toString(extrCalcThetaDegrees, "%.3e"),
                FormatArray.toString(extrCalc.getTranslation(), "%.3e"));

        System.out.printf("0: expected rot=\n%s\nexpected degrees=\n%s\nexpected trans=\n%s\n",
                FormatArray.toString(k2ExtrRot, "%.3e"),
                FormatArray.toString(thetaDegrees, "%.3e"),
                FormatArray.toString(k2ExtrTrans, "%.3e"));

        // ===================================================================================
        /* for the 2nd image pair:
        left extrinsic matrix:
        [[-9.48402795e-02  9.75707586e-01 -1.97484246e-01 -4.67064408e+01]
         [ 7.52568521e-01  2.00130415e-01  6.27366272e-01 -8.09188900e+01]
         [ 6.51648635e-01 -8.91208345e-02 -7.53267239e-01  2.66643941e+02]
         [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  1.00000000e+00]]

         right extrinsic matrix:
         [[-9.03220351e-02  9.76271959e-01 -1.96812071e-01 -1.45613832e+02]
         [ 7.48996521e-01  1.96836697e-01  6.32660673e-01 -8.15340647e+01]
         [ 6.56388713e-01 -9.02683572e-02 -7.49002992e-01  2.66971558e+02]
         [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  1.00000000e+00]]
         */
        filePath1 = path + "left02_masked.png";
        filePath2 = path + "right02_masked.png";

        corres = CorrespondenceMaker.findUsingORB(filePath1, filePath2, nCorners, useToleranceAsStatFactor,
                recalcIterations, tol, errorType, debug);
        assertNotNull(corres);
        x1 = corres.x1;
        x2 = corres.x2;

        x1c = Camera.pixelToCameraCoordinates(x1, k1Intr, radial1, useR2R4);
        x2c = Camera.pixelToCameraCoordinates(x1, k1Intr, radial1, useR2R4);

        for (int i = 0; i < x1[0].length; ++i) {
            System.out.printf("x1[0][%d]=%.4e;  x1[1][%d]=%.4e;  x2[0][%d]=%.4e;  x2[1][%d]=%.4e;\n",
                    i, x1[0][i], i, x1[1][i], i, x2[0][i], i, x2[1][i]);
        }

        left1[0] = new double[]{9.48402795e-02,  9.75707586e-01, -1.97484246e-01};
        left1[1] = new double[]{7.52568521e-01,  2.00130415e-01,  6.27366272e-01};
        left1[2] = new double[]{6.51648635e-01, -8.91208345e-02, -7.53267239e-01};
        right1[0] = new double[]{-9.03220351e-02,  9.76271959e-01, -1.96812071e-01};
        right1[1] = new double[]{7.48996521e-01,  1.96836697e-01,  6.32660673e-01};
        right1[2] = new double[]{6.56388713e-01, -9.02683572e-02, -7.49002992e-01};

        leftInv = MatrixUtil.transpose(left1);
        // apply inv(left) to left to see that the result is the identity matrix
        leftTleft = MatrixUtil.createATransposedTimesA(left1);
        // then apply inv(left) to right to be able to compare the result with the identity matrix,
        //    that is, the rotation of right with respect to left.
        leftTright = MatrixUtil.multiply(leftInv, right1);
        thetasLeftToRight = Rotation.extractThetaFromZYX(leftTright);

        k2ExtrRot = leftTright;
        k2ExtrTrans = new double[]{-4.67064408e+01- -1.45613832e+02, -8.09188900e+01 - -8.15340647e+01,
                2.66643941e+02 - 2.66971558e+02};

        thetas = thetasLeftToRight;
        thetaDegrees = Arrays.copyOf(thetas, thetas.length);
        MatrixUtil.multiply(thetaDegrees, 180./Math.PI);

        extrCalc = Reconstruction.calculateProjectiveMotion(x1c, x2c);

        double[] ti = Arrays.copyOf(extrCalc.getTranslation(), 3);
        //xim = K * xc.  xc = K^-1*xim
        //double[][] k2IntrInv = Camera.createIntrinsicCameraMatrixInverse(k2Intr);
        double[] ti0 = MatrixUtil.multiplyMatrixByColumnVector(k2Intr, ti);
        MatrixUtil.multiply(ti, 1./ti[2]);
        double[] ti00 = MatrixUtil.multiplyMatrixByColumnVector(k2Intr, ti);
        System.out.printf("transformed translation=\n%s\n    =%s\n", FormatArray.toString(ti0, "%.3e"),
                FormatArray.toString(ti00, "%.3e"));

        assertNotNull(extrCalc);
        assertNotNull(extrCalc.getRotation());
        extrCalcRot = extrCalc.getRotation();
        extrCalcThetas = new double[3];
        Rotation.extractThetaFromZYX(extrCalcRot, extrCalcThetas);
        extrCalcThetaDegrees = Arrays.copyOf(extrCalcThetas, extrCalcThetas.length);
        MatrixUtil.multiply(extrCalcThetaDegrees, 180./Math.PI);

        System.out.printf("1: derived rot=\n%s\ndegrees=\n%s\ntrans=\n%s\n",
                FormatArray.toString(extrCalcRot, "%.3e"),
                FormatArray.toString(extrCalcThetaDegrees, "%.3e"),
                FormatArray.toString(extrCalc.getTranslation(), "%.3e"));

        System.out.printf("1: expected rot=\n%s\nexpected degrees=\n%s\nexpected trans=\n%s\n",
                FormatArray.toString(k2ExtrRot, "%.3e"),
                FormatArray.toString(thetaDegrees, "%.3e"),
                FormatArray.toString(k2ExtrTrans, "%.3e"));

    }

    public void estCalculateProjectiveReconstruction() throws Exception {
        System.out.println("testCalculateProjectiveReconstruction");

        double[][] k1 = Zhang98Data.getIntrinsicCameraMatrix();
        double[][] k2 = MatrixUtil.copy(k1);
        //x1, x2 size is 3 X 256
        double[][] x1 = Zhang98Data.getObservedFeaturesInImage(1);
        double[][] x2 = Zhang98Data.getObservedFeaturesInImage(5);

        double[][] x1C = Camera.pixelToCameraCoordinates(x1, k1,
                Zhang98Data.getRadialDistortionR2R4(), true);
        double[][] x2C = Camera.pixelToCameraCoordinates(x2, k2,
                Zhang98Data.getRadialDistortionR2R4(), true);

        // should be 2 solutions
        Reconstruction.ProjectionResults[] results =
                Reconstruction.calculateProjectiveReconstruction(x1C, x2C);

        Camera.CameraExtrinsicParameters kExtr2 = Reconstruction.calculateProjectiveMotion(x1C, x2C);

        System.out.printf("camera2 R from projective motion =\n%s\n",
                FormatArray.toString(kExtr2.getRotation(), "%.3e"));
        System.out.printf("camera2 T from projective motion =\n%s\n",
                FormatArray.toString(kExtr2.getTranslation(), "%.3e"));

        System.out.printf("0) projection matrices=\n%s\n",
                FormatArray.toString(results[0].projectionMatrices, "%.3e"));
        System.out.printf("0) XW=\n%s\n",
                FormatArray.toString(MatrixUtil.transpose(results[0].XW), "%.3e"));

        System.out.printf("1) projection matrices=\n%s\n",
                FormatArray.toString(results[1].projectionMatrices, "%.3e"));
        System.out.printf("1) XW=\n%s\n",
                FormatArray.toString(MatrixUtil.transpose(results[1].XW), "%.3e"));
    }


        /**
       * Test of calculateUsingEssentialMatrix method, of class CameraPose.
       */
    public void estCalculateUsingEssentialMatrix() throws Exception {
        System.out.println("calculateUsingEssentialMatrix");

        int idx1 = 1;
        int idx2 = 5;

        double[][] k1 = Zhang98Data.getIntrinsicCameraMatrix();
        double[][] k2 = MatrixUtil.copy(k1);
        //x1, x2 size is 3 X 256
        double[][] x1 = Zhang98Data.getObservedFeaturesInImage(idx1);
        double[][] x2 = Zhang98Data.getObservedFeaturesInImage(idx2);

        Reconstruction.ReconstructionResults result = 
            Reconstruction.calculateUsingEssentialMatrix(k1, k2, x1, x2);
        
        assertNotNull(result);
        
        System.out.printf("\nresult:\nrot1=%strans1=%s\n", 
            FormatArray.toString(result.k1ExtrRot, "%.4e"),
            FormatArray.toString(result.k1ExtrTrans, "%.4e"));
        System.out.printf("\nresult:\nrot2=%strans2=%s\n", 
            FormatArray.toString(result.k2ExtrRot, "%.4e"),
            FormatArray.toString(result.k2ExtrTrans, "%.4e"));
        
        System.out.printf("\nimg1:\nrot=%strans=%s\n", 
                FormatArray.toString(Zhang98Data.getRotation(idx1), "%.4e"),
                FormatArray.toString(Zhang98Data.getTranslation(idx1), "%.4e"));
        System.out.printf("\nimg5:\nrot=%strans=%s\n", 
                FormatArray.toString(Zhang98Data.getRotation(idx2), "%.4e"),
                FormatArray.toString(Zhang98Data.getTranslation(idx2), "%.4e"));
        
        double[][] diffRSameCenter = Rotation.procrustesAlgorithmForRotation(
            Zhang98Data.getRotation(1), Zhang98Data.getRotation(idx2));
        
        System.out.printf("\ndifference in rot between img1 and img5=\n%s\n", 
           FormatArray.toString(diffRSameCenter, "%.4e"));
        
        double[] out = new double[3];
        extractThetaFromZYX(diffRSameCenter, out);
        
        System.out.printf("\ndifference in rot between img1 and img5 in euler angles=\n");
        for (double a : out) {
            System.out.printf("%.2f ", 180.*a/Math.PI);
        }
        System.out.println();
        double[] diffTrans = MatrixUtil.subtract(Zhang98Data.getTranslation(1), Zhang98Data.getTranslation(5));
        System.out.printf("\ndifference in trans between img1 and img5=\n%s\n",
            FormatArray.toString(diffTrans, "%.3e"));
        System.out.flush();
    }

    /**
     * Test of calculateReconstruction method, of class Reconstruction.
     */
    public void estCalculateReconstructionWithIntrinsicCamera() throws NotConvergedException {
        
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
        // note: the checkerboard sqaures are 30mm in WCS metric
        //
        // note: radial distortion should be corrected:  use on the original coordinates:
        //     x_corrected = x*(1 + k1*r^2 + k2r^4) where r is distance of point from cc.
        
        double[][] k1Intr;
        double[][] k2Intr;
            k1Intr = Camera.createIntrinsicCameraMatrix(533.07, 341.6, 234.3);
            k2Intr = Camera.createIntrinsicCameraMatrix(536.7, 326.5, 249.3);
        
        double[][] k1ExtrRot = MatrixUtil.createIdentityMatrix(3);
        double[] k1ExtrTrans = new double[]{0, 0, 0};
        
        double[][] k2ExtrRot = Rotation.createRodriguesFormulaRotationMatrix(
            new double[]{0.00611, 0.00409, -0.00359});
        double[] k2ExtrTrans = new double[]{-99.85, 0.82, 0.44};
        //double[] k2ExtrTransRev = Arrays.copyOf(k2ExtrTrans, k2ExtrTrans.length);
        //MatrixUtil.multiply(k2ExtrTransRev, -1);
        
        System.out.printf("expected k2ExtrRot from Rodrigues formula\n=%s\n", FormatArray.toString(k2ExtrRot, "%.3e"));
        /*
        [junit] =1.000e+00, 3.577e-03, 4.101e-03 
        [junit] -3.602e-03, 1.000e+00, 6.103e-03 
        [junit] -4.079e-03, -6.117e-03, 1.000e+00 
        */
        
        double[][] x1 = new double[3][10];
        x1[0] = new double[]{129, 160, 140, 226, 225, 232, 341, 407, 532, 527};
        x1[1] = new double[]{145, 319, 361, 391, 61, 284, 289, 122, 48, 302};
        x1[2] = new double[]{1, 1, 1, 1, 1, 1, 1, 1,  1, 1};
        
        double[][] x2 = new double[3][10];
        x2[0] = new double[]{76, 110, 84, 164, 110, 124, 218, 275, 401, 402};
        x2[1] = new double[]{168, 331, 372, 401, 78, 295, 302, 134, 52, 318};
        x2[2] = new double[]{1, 1, 1, 1, 1, 1, 1, 1,  1, 1};
        
        /*
        System.out.printf("results=%s\n\n", rr.toString());
        
  //multiply by focal length?  or K?   see if triangulation in szeliski has advice
        
        System.out.println("XW points normalized by last coordinates:");
        double[] pt = new double[4];
        for (int j = 0; j < rr.XW[0].length; ++j) {    
            for (int i = 0; i < rr.XW.length; ++i) {
                pt[i] = rr.XW[i][j]/rr.XW[3][j];
            }
            System.out.printf("%s\n", FormatArray.toString(pt, "%.3e"));
        }
        */
    }
    
    public void estAffine() throws NotConvergedException, IOException {
    
        double[][] xIn = loadNCBook(false);

        Reconstruction.OrthographicProjectionResults affineProj =
            Reconstruction.calculateAffineReconstruction(xIn, 2);

    }

    public void estParaperspective() throws IOException, NotConvergedException {
        double[][] xIn = loadNCBook(false);

        Reconstruction.ParaperspectiveProjectionResults paraProj
                = Reconstruction.calculateParaperspectiveReconstruction(xIn, 2);
    }
    
    private double[][] loadNCBook(boolean rotateBy90) throws IOException, NotConvergedException {
        
        int maxDimension = 256;//512;

        ImageProcessor imageProcessor = new ImageProcessor();
        ImageSegmentation imageSegmentation = new ImageSegmentation();

        String[][] filePairs = new String[1][];
        filePairs[0] = new String[]{
            "nc_book_01.png",
            "nc_book_02.png"};
        
        boolean binImages = false;//true;
        
        int rotate = rotateBy90 ? 1 : 0;
        
        int i, ii, np;
            
        String lbl = "_";
        if (rotate == 1) {
            lbl = "rot90_";
        }

        String fileName1 = filePairs[0][0];
        String fileName2 = filePairs[0][1];

        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        int w1 = img1.getWidth();
        int h1 = img1.getHeight();
        if (binImages) {
            int binFactor1 = (int) Math.ceil(Math.max(
                    (float) w1 / maxDimension,
                    (float) h1 / maxDimension));
            img1 = imageProcessor.binImage(img1, binFactor1);
            //MiscDebug.writeImage(img1, "_"  + fileName1Root);
        }
        GreyscaleImage img1GS = img1.copyToGreyscale2();

        idx = fileName2.lastIndexOf(".");
        String fileName2Root = fileName2.substring(0, idx);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
        int w2 = img2.getWidth();
        int h2 = img2.getHeight();
        if (binImages) {
            int binFactor2 = (int) Math.ceil(Math.max(
                    (float) w2 / maxDimension,
                    (float) h2 / maxDimension));
            img2 = imageProcessor.binImage(img2, binFactor2);
        }
        GreyscaleImage img2GS = img2.copyToGreyscale2();

        Transformer tr = new Transformer();
        TransformationParameters params = new TransformationParameters();
        TransformationParameters paramsRev = new TransformationParameters();

        if (rotate == 1) {
            params.setTranslationY(img2GS.getWidth());
            params.setRotationInDegrees(90);

            img2GS = tr.applyTransformation(img2GS, params,
                    img2GS.getHeight(), img2GS.getWidth());

            paramsRev.setTranslationX(img2GS.getHeight());
            paramsRev.setRotationInDegrees(-90);
        }

        int x, y;

        if (binImages) {
            np = 300;
        } else {
            np = 600;
        }

        ORB orb1 = new ORB(img1GS, np);
        orb1.detectAndExtract();
        //orb1.overrideToAlsoCreate1stDerivKeypoints();
        //orb1.overrideToNotCreateATrousKeypoints();

        ORB orb2 = new ORB(img2GS, np);
        orb2.detectAndExtract();
        //orb2.overrideToAlsoCreate1stDerivKeypoints();
        //orb2.overrideToNotCreateATrousKeypoints();

        ORB.Descriptors d1 = orb1.getAllDescriptors();
        ORB.Descriptors d2 = orb2.getAllDescriptors();
        double[][] xKP1 = orb1.getAllKeyPointsHomogenous();
        double[][] xKP2 = orb2.getAllKeyPointsHomogenous();

        int nKP1 = xKP1[0].length;
        int nKP2 = xKP2[0].length;

        int x1, x2, y1, y2;
        Image tmp1 = img1GS.copyToColorGreyscale();
        for (i = 0; i < nKP1; ++i) {
            y = (int) xKP1[1][i];
            x = (int) xKP1[0][i];
            ImageIOHelper.addPointToImage(x, y, tmp1, 2, 255, 0, 0);
        }
        //MiscDebug.writeImage(tmp1, "_kp_gs_" + lbl + fileName1Root);

        Image tmp2 = img2GS.copyToColorGreyscale();
        for (i = 0; i < nKP2; ++i) {
            y = (int) xKP2[1][i];
            x = (int) xKP2[0][i];
            ImageIOHelper.addPointToImage(x, y, tmp2, 2, 255, 0, 0);
        }
        MiscDebug.writeImage(tmp1, "_kp_gs_" + lbl + fileName1Root);
        MiscDebug.writeImage(tmp2, "_kp_gs_" + lbl + fileName2Root);

        // kp1 and kp2 are in col major format (col, row) and we want normalized points, so converting to double[][]
        double[][] xKP1n = MatrixUtil.copy(xKP1);
        double[][] xKP2n = MatrixUtil.copy(xKP2);
        double[][] t1 = EpipolarNormalizationHelper.unitStandardNormalize(xKP1n);
        double[][] t2 = EpipolarNormalizationHelper.unitStandardNormalize(xKP2n);

        ORBMatcher.FitAndCorres fitC = null;
        if (true /*normalize points*/) {
            // col, row
            fitC = ORBMatcher.matchDescriptors(d1, d2, xKP1n, xKP2n);
        } else {
            fitC = ORBMatcher.matchDescriptors(d1, d2, xKP1, xKP2);
        }

        if (fitC.mI == null) {
            return null;
        }
        int[][] mi = fitC.mI;
        if (fitC.mIF != null) {
            mi = fitC.mIF;
        }

        int idx1, idx2;
        CorrespondencePlotter plotter = new CorrespondencePlotter(tmp1, tmp2);
        for (i = 0; i < mi.length; ++i) {
            idx1 = mi[i][0];
            idx2 = mi[i][1];
            x1 = (int) xKP1[0][idx1];
            y1 = (int) xKP1[1][idx1];
            x2 = (int) xKP2[0][idx2];
            y2 = (int) xKP2[1][idx2];
            plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 1);
        }
        plotter.writeImage("_r_corres_orb_gs_" + lbl + fileName1Root);

        /*
        x the image coordinates of feature correspondences in 2 or more
         * images.  format is 2 X (nImages * nFeatures) where row 0 holds the x-coordinates
         * and row 1 holds the y-coordinates and each image's features are given
         * before the next and the features are ordered in the same manner within
         * all images.
         * for example: row 0 = [img_0_feature_0, ... img_0_feature_n-1, ... img_m-1_feature_0,...
         *     img_m-1_feature_n-1].
         */

        double[][] xIn = MatrixUtil.zeros(2, 2*mi.length);
        for (i = 0; i < mi.length; ++i) {
            idx1 = mi[i][0];
            xIn[0][i] = (int) xKP1[0][idx1];
            xIn[1][i] = (int) xKP1[1][idx1];
        }
        for (i = 0; i < mi.length; ++i) {
            idx2 = mi[i][1];
            xIn[0][i + mi.length] = (int) xKP2[0][idx2];
            xIn[1][i + mi.length] = (int) xKP2[1][idx2];
        }

        return xIn;
    }

    /**
     * given a list of keypoints in format (col, row), write a double array of size [2 X kP.size()] where the
     * first row holds the col values, the second row holds the row values and the third row holds values "1" for
     * the z-axis representation in homogeneous coordinates.
     * @param kP input list of keypoints in format (col, row)
     * @return the data in format [2 X kP.size()]
     */
    private double[][] convertToArray(List<PairInt> kP) {
        int n = kP.size();
        double[][] x = MatrixUtil.zeros(3, n);
        for (int i = 0; i < n; ++i) {
            x[0][i] = kP.get(i).getX();
            x[1][i] = kP.get(i).getY();
        }
        Arrays.fill(x[2], 1);
        return x;
    }

}
