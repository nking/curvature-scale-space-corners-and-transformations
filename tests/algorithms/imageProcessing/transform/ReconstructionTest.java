package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.features.RANSACSolver;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.imageProcessing.matching.CorrespondenceMaker;
import algorithms.imageProcessing.matching.CorrespondenceMaker.CorrespondenceList;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.matching.ORBMatcher;
import static algorithms.imageProcessing.transform.Rotation.extractThetaFromZYX;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscDebug;
import algorithms.util.*;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertNotNull;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
 *
 * @author nichole
 */
public class ReconstructionTest extends TestCase {

    protected final static String sep = System.getProperty("file.separator");

    public ReconstructionTest() {
    }

    public void testProjectiveWithZhang() throws IOException, NotConvergedException {

        int idx1 = 1;
        int idx2 = 2;

        double[][] x1 = Zhang98Data.getObservedFeaturesInImage(idx1);
        double[][] x2 = Zhang98Data.getObservedFeaturesInImage(idx2);
        double[][] intr = Zhang98Data.getIntrinsicCameraMatrix();
        double[][] invIntr = Camera.createIntrinsicCameraMatrixInverse(intr);
        double[] radial = Zhang98Data.getRadialDistortionR2R4();
        double[][] r1 = Zhang98Data.getRotation(idx1);
        double[] t1 = Zhang98Data.getTranslation(idx1);
        double[][] r2 = Zhang98Data.getRotation(idx2);
        double[] t2 = Zhang98Data.getTranslation(idx2);
        double[][] camera1 = Camera.createCamera(intr, r1, t1); //P = intr * [R | t]
        double[][] camera2 = Camera.createCamera(intr, r2, t2); //P = intr * [R | t]
        System.out.printf("P1=K*[R1|t1]=camera1=\n%s\n", FormatArray.toString(camera1, "%.3e"));
        System.out.printf("P2=K*[R2|t2]=camera2=\n%s\n", FormatArray.toString(camera2, "%.3e"));
        //double[][] p1 = Camera.createExtrinsicCameraMatrix(r1, t1);
        //double[][] p2 = Camera.createExtrinsicCameraMatrix(r2, t2);
        double[][] x1c = Camera.pixelToCameraCoordinates(x1, intr, radial, true);
        double[][] x2c = Camera.pixelToCameraCoordinates(x2, intr, radial, true);
        double[][] XW = Zhang98Data.getFeatureWCS();

        double[] r1Vec = Rotation.extractRotationVectorRodrigues(r1);
        double[] r2Vec = Rotation.extractRotationVectorRodrigues(r2);
        double[][] r1Transposed = MatrixUtil.transpose(r1);
        double[][] r1I = MatrixUtil.multiply(r1Transposed, r1);
        double[][] r2RelToR1 = MatrixUtil.multiply(r1Transposed, r2);
        System.out.printf("t1=%s\n", FormatArray.toString(t1, "%.3e"));
        System.out.printf("t2=%s\n", FormatArray.toString(t2, "%.3e"));
        System.out.printf("diff t=%s\n", FormatArray.toString(MatrixUtil.subtract(t2, t1), "%.3e"));
        System.out.printf("r1=%s\n", FormatArray.toString(r1, "%.3e"));
        System.out.printf("r2=%s\n", FormatArray.toString(r2, "%.3e"));
        System.out.printf("r1Tr1=%s\n", FormatArray.toString(r1I, "%.3e"));
        System.out.printf("r2RelToR1=%s\n", FormatArray.toString(r2RelToR1, "%.3e"));

        int n = x1[0].length;

        double[][] x1Pt = MatrixUtil.zeros(3, 1);
        double[][] x2Pt = MatrixUtil.zeros(3, 1);
        double[][] x1iPt = MatrixUtil.zeros(3, 1);
        double[][] x2iPt = MatrixUtil.zeros(3, 1);
        int ii, j, i;

        // put x into camera coordinates reference frame:
        Triangulation.WCSPt wcsPt;
        for (ii = 0; ii < n; ++ii) {
            for (j = 0; j < 3; ++j) {
                x1Pt[j][0] = x1c[j][ii];
                x2Pt[j][0] = x2c[j][ii];
            }
            // better answer produced by using camera coords xc, yc and extrinsic projection matrix P=[R|t]
            wcsPt = Triangulation.calculateWCSPoint(camera1, camera2, x1Pt, x2Pt);
            //wcsPt = Triangulation.calculateWCSPoint(intr, r1, t1, intr, r2, t2, x1Pt, x2Pt);
            MatrixUtil.multiply(wcsPt.X, 1. / wcsPt.X[3]);

            System.out.printf("WCS[%d]=%s,\t alpha=%.3e\n",
                    ii, FormatArray.toString(wcsPt.X, "%.3e"),
                    wcsPt.alpha);
            System.out.printf("expected=%s\n",
                    FormatArray.toString(MatrixUtil.extractColumn(XW, ii), "%.3e"));
        }

        Reconstruction.ProjectionResults[] pr = Reconstruction.calculateProjectiveReconstruction(
            intr, intr, x1c, x2c);
        for (ii = 0; ii < pr.length; ++ii) {
            double[][] p = pr[ii].projectionMatrices;
            double[][] r = MatrixUtil.multiply(invIntr, MatrixUtil.copySubMatrix(p, 3, 5, 0, 2));
            double[] t = MatrixUtil.multiplyMatrixByColumnVector(invIntr,
                    MatrixUtil.extractColumn(MatrixUtil.copySubMatrix(p, 3, 5, 3, 3), 0));
            double distR = Rotation.distanceUsingRigidBodyDisplacements(r, r2RelToR1, false);
            System.out.printf("pProj=\n%s\n", FormatArray.toString(p, "%.3e"));
            System.out.printf("r=\n%s\n", FormatArray.toString(r, "%.3e"));
            System.out.printf("t=\n%s\n", FormatArray.toString(t, "%.3e"));
            System.out.printf("r2RelToR1=%s\n", FormatArray.toString(r2RelToR1, "%.3e"));
            System.out.printf("distR=%.3e\n", distR);
            double[][] pXW = pr[ii].XW;
            System.out.printf("pXW[%d]=%s\n",
                    ii, FormatArray.toString(MatrixUtil.transpose(pXW), "%.3e"));
        }

    }

    public void testProjectiveWithBouguetIm2LeftRight() throws IOException, NotConvergedException {

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
        // Page 3 of miscNotes/bouguetj_5th_calibration_example.pdf shows that the triangulation results
        // of points in the world coordinate system should have Z coordinates between 350 and 500.

        // note: the checkerboard squares are 30mm in WCS metric
        //
        // note: radial distortion should be corrected:  use on the original coordinates:
        //     x_corrected = x*(1 + k1*r^2 + k2r^4) where r is distance of point from cc.

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

        double[][] k2IntrLeft = Camera.createIntrinsicCameraMatrix(533.5, 341.6, 235.2);
        double[][] k2IntrRight = Camera.createIntrinsicCameraMatrix(536.6, 326.3, 250.1);
        double[][] K2IntrLeftInv = Camera.createIntrinsicCameraMatrixInverse(k2IntrLeft);
        double[][] K2IntrRightInv = Camera.createIntrinsicCameraMatrixInverse(k2IntrRight);

        boolean useR2R4 = true;
        double[] radial2Left = new double[]{-0.288, 0.097, 0.001};
        double[] radial2Right = new double[]{-0.289,0.107, -0.001};

        double[][] rtLeft = new double[3][];
        rtLeft[0] = new double[]{9.48402795e-02,  9.75707586e-01, -1.97484246e-01, -4.67064408e+01};
        rtLeft[1] = new double[]{7.52568521e-01,  2.00130415e-01,  6.27366272e-01, -8.09188900e+01};
        rtLeft[2] = new double[]{6.51648635e-01, -8.91208345e-02, -7.53267239e-01,  2.66643941e+02};
        double[][] rtRight = new double[3][];
        rtRight[0] = new double[]{-9.03220351e-02,  9.76271959e-01, -1.96812071e-01, -1.45613832e+02};
        rtRight[1] = new double[]{7.48996521e-01,  1.96836697e-01,  6.32660673e-01, -8.15340647e+01};
        rtRight[2] = new double[]{6.56388713e-01, -9.02683572e-02, -7.49002992e-01,  2.66971558e+02};

        double[][] r2Left = MatrixUtil.copySubMatrix(rtLeft, 0, 2, 0,2);
        double[][] r2Right = MatrixUtil.copySubMatrix(rtRight, 0, 2, 0,2);
        double[][] r2LeftOrth = Rotation.orthonormalizeUsingSVD(r2Left);
        double[][] r2RightOrth = Rotation.orthonormalizeUsingSVD(r2Right);
        double[] t2Left = MatrixUtil.extractColumn(rtLeft, 3);
        double[] t2Right = MatrixUtil.extractColumn(rtRight, 3);

        double[][] camera2Left = Camera.createCamera(k2IntrLeft, r2Left, t2Left); //P = intr * [R | t]
        double[][] camera2Right = Camera.createCamera(k2IntrRight, r2Right, t2Right); //P = intr * [R | t]
        System.out.printf("P2Left=K*[R1|t1]=camera2Left=\n%s\n", FormatArray.toString(camera2Left, "%.3e"));
        System.out.printf("P2Right=K*[R2|t2]=camera2Right=\n%s\n", FormatArray.toString(camera2Right, "%.3e"));

        double[][] r2LeftTransposed = MatrixUtil.transpose(r2Left);
        double[][] r2LeftI = MatrixUtil.multiply(r2LeftTransposed, r2Left);
        double[][] r2RightRelToLeft = MatrixUtil.multiply(r2LeftTransposed, r2Right);
        double[][] r2LeftOrthTransposed = MatrixUtil.transpose(r2LeftOrth);
        double[][] r2LeftOrthI = MatrixUtil.multiply(r2LeftOrthTransposed, r2LeftOrth);
        double[][] r2RightOrthRelToLeft = MatrixUtil.multiply(r2LeftOrthTransposed, r2RightOrth);

        double[] t2Diff = MatrixUtil.subtract(t2Right, t2Left);
        double[] r2LeftVec = Rotation.extractRotationVectorRodrigues(r2Left);
        double[] r2RightVec = Rotation.extractRotationVectorRodrigues(r2Right);
        double[] r2LeftOrthVec = Rotation.extractRotationVectorRodrigues(r2LeftOrth);
        double[] r2RightOrthVec = Rotation.extractRotationVectorRodrigues(r2RightOrth);
        double[] r2RightRelToLeftVec = Rotation.extractRotationVectorRodrigues(r2RightRelToLeft);

        double[] r2RightOrthRelToLeftVec = Rotation.extractRotationVectorRodrigues(r2RightOrthRelToLeft);
        double[] r2LeftIVec = Rotation.extractRotationVectorRodrigues(r2LeftI);
        //double[][] r2Direc = Rotation.rotationBetweenTwoDirections0(r2LeftVec, r2RightVec);
        //double[][] r2OrthDirec = Rotation.rotationBetweenTwoDirections0(r2LeftOrthVec, r2RightOrthVec);
        //double[] r2DirecVec = Rotation.extractRotationVectorRodrigues(r2Direc);
        //double[] r2OrthDirecVec = Rotation.extractRotationVectorRodrigues(r2OrthDirec);
        //expected om = 0.00611, 0.00489, -0.00359   +- 0.0027, 0.00308, 0.00029
        //expected t = -99.84929, 0.82221, 0.43647   +- 0.142, 0.11352, 0.49773

        boolean passive = false;
        double[] r2LeftToRightVec = Arrays.copyOf(r2RightRelToLeftVec, r2RightRelToLeftVec.length);
        MatrixUtil.multiply(r2LeftToRightVec, -1);
        double[][] r2LeftToRight = Rotation.createRotationRodriguesFormula(r2LeftToRightVec, passive);
        int i, j;
        {//DEBUG
            System.out.printf("r2Left=\n%s\n", FormatArray.toString(r2Left, "%.3e"));
            System.out.printf("r2Right=\n%s\n", FormatArray.toString(r2Right, "%.3e"));
            System.out.printf("r2LeftOrth=\n%s\n", FormatArray.toString(r2LeftOrth, "%.3e"));
            System.out.printf("r2RightOrth=\n%s\n", FormatArray.toString(r2RightOrth, "%.3e"));
            System.out.printf("r2LeftI=\n%s\n", FormatArray.toString(r2LeftI, "%.3e"));
            System.out.printf("r2RightRelToLeft=\n%s\n", FormatArray.toString(r2RightRelToLeft, "%.3e"));
            System.out.printf("r2LeftToRight=\n%s\n", FormatArray.toString(r2LeftToRight, "%.3e"));
            System.out.printf("r2LeftOrthI=\n%s\n", FormatArray.toString(r2LeftOrthI, "%.3e"));
            System.out.printf("r2RightOrthRelToLeft=\n%s\n", FormatArray.toString(r2RightOrthRelToLeft, "%.3e"));
            System.out.printf("\n");
            System.out.printf("r2LeftVec=\n%s\n", FormatArray.toString(r2LeftVec, "%.3e"));
            System.out.printf("r2RightVec=\n%s\n", FormatArray.toString(r2RightVec, "%.3e"));
            System.out.printf("r2LeftOrthVec=\n%s\n", FormatArray.toString(r2LeftOrthVec, "%.3e"));
            System.out.printf("r2RightOrthVec=\n%s\n", FormatArray.toString(r2RightOrthVec, "%.3e"));
            System.out.printf("r2LeftOrthVec=\n%s\n", FormatArray.toString(r2LeftOrthVec, "%.3e"));
            // r2RightRelToLeftVec is close to EXPECTED, but is of negative sign (rotation from 2 to 1)
            System.out.printf("r2RightRelToLeftVec=\n%s\n", FormatArray.toString(r2RightRelToLeftVec, "%.3e"));
            System.out.printf("r2RightOrthRelToLeftVec=\n%s\n", FormatArray.toString(r2RightOrthRelToLeftVec, "%.3e"));
            System.out.printf("r2LeftIVec=\n%s\n", FormatArray.toString(r2LeftIVec, "%.3e"));
            // r2Direc is closest to EXPECTED
            //System.out.printf("r2Direc=\n%s\n", FormatArray.toString(r2Direc, "%.3e"));
            //System.out.printf("r2OrthDirec=\n%s\n", FormatArray.toString(r2OrthDirec, "%.3e"));
            //System.out.printf("r2DirecVec=\n%s\n", FormatArray.toString(r2DirecVec, "%.3e"));
            //System.out.printf("r2OrthDirecVec=\n%s\n", FormatArray.toString(r2OrthDirecVec, "%.3e"));
            System.out.printf("EXPECTED OM=(0.00611, 0.00489, -0.00359)   +- 0.0027, 0.00308, 0.00029\n");
            System.out.printf("\n");
            System.out.printf("t2Left=\n%s\n", FormatArray.toString(t2Left, "%.3e"));
            System.out.printf("t2Right=\n%s\n", FormatArray.toString(t2Right, "%.3e"));
            System.out.printf("t2Diff=\n%s\n", FormatArray.toString(t2Diff, "%.3e"));
            //expected om = 0.00611, 0.00489, -0.00359   +- 0.0027, 0.00308, 0.00029
            //expected t = -99.84929, 0.82221, 0.43647   +- 0.142, 0.11352, 0.49773
            System.out.printf("EXPECTED T=(-99.84929, 0.82221, 0.43647)   +- 0.142, 0.11352, 0.49773\n");
            System.out.printf("\n");
        }

        System.out.printf("K1Inv = \n%s\n",
                FormatArray.toString(Camera.createIntrinsicCameraMatrixInverse(k2IntrLeft), "%.3e"));
        System.out.printf("K2Inv = \n%s\n",
                FormatArray.toString(Camera.createIntrinsicCameraMatrixInverse(k2IntrRight), "%.3e"));

        Corres corres = makeBouguetIm2LeftRightCorres();
        double[][] x1 = corres.x1;
        double[][] x2 = corres.x2;
        double[][] x1c = Camera.pixelToCameraCoordinates(x1, k2IntrLeft, radial2Left, useR2R4);
        double[][] x2c = Camera.pixelToCameraCoordinates(x2, k2IntrRight, radial2Right, useR2R4);

        double[][] x1Pt = MatrixUtil.zeros(3, 1);
        double[][] x2Pt = MatrixUtil.zeros(3, 1);
        double[][] x1iPt = MatrixUtil.zeros(3, 1);
        double[][] x2iPt = MatrixUtil.zeros(3, 1);
        int ii;
        int n = x1[0].length;

        // put x into camera coordinates reference frame.  radial distortion has been removed.
        Triangulation.WCSPt wcsPt;
        double[][] XWFromT = new double[3][n];

        for (ii = 0; ii < n; ++ii) {
            for (j = 0; j < 3; ++j) {
                x1Pt[j][0] = x1[j][ii];
                x2Pt[j][0] = x2[j][ii];
            }
            wcsPt = Triangulation.calculateWCSPoint(camera2Left, camera2Right, x1Pt, x2Pt);
            MatrixUtil.multiply(wcsPt.X, 1. / wcsPt.X[3]);

            System.out.printf("WCS[%d]=%s,\t alpha=%.3e\n",
                    ii, FormatArray.toString(wcsPt.X, "%.3e"),
                    wcsPt.alpha);

            for (j = 0; j < 3; ++j) {
                XWFromT[j][ii] = wcsPt.X[j];// last column is '1'
            }
        }

        // having x1 and XWFromT or x1 and XWFromT, can test the camera pose methods
        Camera.CameraPoseParameters cPose1 = CameraPose.calculatePoseFromXXW(x1, XWFromT);
        Camera.CameraPoseParameters cPose2 = CameraPose.calculatePoseFromXXW(x2, XWFromT);

        // differences from expected for camera pose
        //TODO: look into the intrinsic scale factors in cPose1 and cPose2
        double[][] diffIntr1 = new double[3][3];
        double[][] diffIntr2 = new double[3][3];
        double[][] diffRot1 = Rotation.procrustesAlgorithmForRotation(
                cPose1.getExtrinsicParameters().getRotation(),
                r2LeftOrth); // or rtLeft
        double[][] diffRot2 = Rotation.procrustesAlgorithmForRotation(
            cPose2.getExtrinsicParameters().getRotation(),
            r2RightOrth); // or rtRight
        double[] diffT1 = MatrixUtil.subtract(cPose1.getExtrinsicParameters().getTranslation(),
                t2Left);
        double[] diffT2 = MatrixUtil.subtract(cPose2.getExtrinsicParameters().getTranslation(),
                t2Right);

        for (int i2 = 0; i2 < 3; ++i2) {
            for (int j2 = 0; j2 < 3; ++j2) {
                diffIntr1[i2][j2] = cPose1.getIntrinsicParameters().getIntrinsic()[i2][j2] - k2IntrLeft[i2][j2];
                if (k2IntrLeft[i2][j2] != 0) {
                    diffIntr1[i2][j2] = Math.abs(diffIntr1[i2][j2] / k2IntrLeft[i2][j2]);
                }
                diffIntr2[i2][j2] = cPose2.getIntrinsicParameters().getIntrinsic()[i2][j2] - k2IntrRight[i2][j2];
                if (k2IntrRight[i2][j2] != 0) {
                    diffIntr2[i2][j2] = Math.abs(diffIntr2[i2][j2] / k2IntrRight[i2][j2]);
                }
            }
        }
        //TODO: add asserts

        boolean calibrated = true;
        double[][] x1n, x2n;
        if (calibrated) {
            x1n = MatrixUtil.copy(x1c);
            x2n = MatrixUtil.copy(x2c);
        } else {
            x1n = MatrixUtil.copy(x1);
            x2n = MatrixUtil.copy(x2);
        }
        double[][] t1 = EpipolarNormalizationHelper.unitStandardNormalize(x1n);
        double[][] t2 = EpipolarNormalizationHelper.unitStandardNormalize(x2n);
        EpipolarTransformer tr = new EpipolarTransformer();
        double[][] fmN = tr.calculateEpipolarProjection2(x1n, x2n, calibrated);
        //de-normalize to work with image coordinates
        double[][] fm = MatrixUtil.copy(fmN); // this is em when calibrated
        EpipolarNormalizationHelper.denormalizeFM(fm, t1, t2);
        double[][] e1e2 = EpipolarTransformer.calculateEpipoles(new DenseMatrix(fm));
        double[][] h = Reconstruction.calculateProjectiveHomographyWithLeastSquares(
            x1, x2, fm, e1e2[1]);
        System.out.printf("H=\n%s\n", FormatArray.toString(h, "%.3e"));
        {
            SVD svdH = SVD.factorize(new DenseMatrix(h));
            MatrixUtil.multiply(h,1./svdH.getS()[1]);
            System.out.printf("H/s[1]=\n%s\n", FormatArray.toString(h, "%.3e"));
        }
        double[][] camera2PH = MatrixUtil.multiply(k2IntrRight, h);


        // non-planar, so this is not working so well...
        Reconstruction.ProjectionResults[] pr = Reconstruction.calculateProjectiveReconstruction(
                k2IntrLeft, k2IntrRight, x1c, x2c);
        // TODO: consider refine/optimize the rotation and translation

        for (ii = 0; ii < pr.length; ++ii) {
            double[][] p = pr[ii].projectionMatrices;
            double[][] r = MatrixUtil.multiply(K2IntrRightInv, MatrixUtil.copySubMatrix(p, 3, 5, 0, 2));
            double[][] rOrth = Rotation.orthonormalizeUsingSVD(r);
            double[] t = MatrixUtil.multiplyMatrixByColumnVector(K2IntrRightInv,
                    MatrixUtil.extractColumn(MatrixUtil.copySubMatrix(p, 3, 5, 3, 3), 0));
            double distR = Rotation.distanceUsingRigidBodyDisplacements(r, r2RightRelToLeft, false);
            System.out.printf("pProj=\n%s\n", FormatArray.toString(p, "%.3e"));
            System.out.printf("r=\n%s\n", FormatArray.toString(r, "%.3e"));
            System.out.printf("rOrth=\n%s\n", FormatArray.toString(rOrth, "%.3e"));
            System.out.printf("t=\n%s\n", FormatArray.toString(t, "%.3e"));
            //System.out.printf("r2Direc=\n%s\n", FormatArray.toString(r2Direc, "%.3e"));
            System.out.printf("distR=%.3e\n", distR);

            //double[][] rDiffFromDirec = Rotation.rotationBetweenTwoDirections0(
            //        Rotation.extractRotationVectorRodrigues(r), r2DirecVec);
            //System.out.printf("rDiffFromDirec=\n%s\n", FormatArray.toString(rDiffFromDirec, "%.3e"));
            //double[] degreesDiff = Rotation.extractThetaFromZYX(rDiffFromDirec);
            //MatrixUtil.multiply(degreesDiff, 180./Math.PI);
            // can see the z-axis negative direction here:
            //System.out.printf("rDiffFromDirec degrees=%s\n", FormatArray.toString(degreesDiff, "%.3e"));

            double[][] pXW = pr[ii].XW;
            System.out.printf("pXW[%d]=%s\n",
                    ii, FormatArray.toString(MatrixUtil.transpose(pXW), "%.3e"));
        }


        Reconstruction.ReconstructionResults rr = Reconstruction.calculateUsingEssentialMatrix(
                k2IntrLeft, k2IntrRight, x1, x2);
        System.out.printf("EM: r1=\n%s\n", FormatArray.toString(rr.k1ExtrRot, "%.3e"));
        System.out.printf("EM: r2=\n%s\n", FormatArray.toString(rr.k2ExtrRot, "%.3e"));
        System.out.printf("EM: t1=\n%s\n", FormatArray.toString(rr.k1ExtrTrans, "%.3e"));
        System.out.printf("EM: t2=\n%s\n", FormatArray.toString(rr.k2ExtrTrans, "%.3e"));
        double[][] XW = rr.XW;
        for (ii = 0; ii < XW[0].length; ++ii) {
            System.out.printf("EM XW[%d]=%s\n",
                    ii, FormatArray.toString(MatrixUtil.extractColumn(XW, ii), "%.3e"));
        }

    }

    public void estProjectiveWithBouguetLeftIm1Im2() throws IOException, NotConvergedException {

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
        double[] rVecRodr = Rotation.extractRotationVectorRodrigues(leftTright);
        double[] rAxis = Rotation.extractRotationAxisFromZXY(leftTright);
        System.out.printf("\nthetas (R to L)=%s\n", FormatArray.toString(thetasRightToLeft, "%.3e"));
        System.out.printf("thetas (L to R)=%s\n", FormatArray.toString(thetasLeftToRight, "%.3e"));
        System.out.printf("\nthetas=%s\n", FormatArray.toString(thetasLeftToRight, "%.3e"));
        System.out.printf("rVecRodr=%s\n", FormatArray.toString(rVecRodr, "%.3e"));
        System.out.printf("rAxis=%s\n", FormatArray.toString(rAxis, "%.3e"));
        System.out.printf("expect=%s\n\n", FormatArray.toString(new double[]{0.00669, 0.00452, -0.0035}, "%.3e"));


        double[] rodriguesL1 = Rotation.extractRotationVectorRodrigues(left1);
        double[] rodriguesL2 = Rotation.extractRotationVectorRodrigues(right1);

        //what rotation
        double[][] diffRSameCenter = Rotation.procrustesAlgorithmForRotation(left1, right1);
        // expecting: 0.00669, 0.00452, -0.0035.  not a bad approximation for the first image set.

        double[][] k1Intr = Camera.createIntrinsicCameraMatrix(533.5, 341.6, 235.2);
        double[][] k2Intr = Camera.createIntrinsicCameraMatrix(536.6, 326.3, 250.1);

        boolean passive = false;

        double[][] k1ExtrRot = MatrixUtil.createIdentityMatrix(3);
        double[] k1ExtrTrans = new double[]{0, 0, 0};
        double[][] k2ExtrRot = Rotation.createRotationRodriguesFormula(
                new double[]{0.00611, 0.00409, -0.00359}, passive);
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

    //NOTE: should not use the Zhang data for this.  x,y are not homogenous coords
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
                Reconstruction.calculateProjectiveReconstruction(k1, k2, x1C, x2C);

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

        boolean passive = false;

        double[][] k2ExtrRot = Rotation.createRotationRodriguesFormula(
            new double[]{0.00611, 0.00409, -0.00359}, passive);
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
     * the z-axis representation in homogenous coordinates.
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

    private static class Corres {
        double[][] x1;
        double[][] x2;
    }

    public Corres makeBouguetIm2LeftRightCorres() throws IOException {

        boolean calibrated = false;
        double[][] left2 = getBouguetIm2LeftCorresDouble();
        double[][] right2 = getBouguetIm2RightCorresDouble();
        double[][] t1 = EpipolarNormalizationHelper.unitStandardNormalize(left2);
        double[][] t2 = EpipolarNormalizationHelper.unitStandardNormalize(right2);
        RANSACSolver solver = new RANSACSolver();

        boolean useToleranceAsStatFactor = true;
        boolean reCalcIterations = false;// possibly faster if set to true
        double tol = 2.;
        ErrorType errorType = ErrorType.SAMPSONS;

        EpipolarTransformationFit fit = solver.calculateEpipolarProjection(
                new DenseMatrix(left2), new DenseMatrix(right2), errorType, useToleranceAsStatFactor,
                tol, reCalcIterations, calibrated);

        System.out.printf("epipolar fit=%s\n", fit.toString());

        int i;
        String path = ResourceFinder.findTestResourcesDirectory() + sep + "bouguet_stereo"
                + sep + "left02_masked.png";
        File f = new File(path);
        if (!f.exists()) {
            throw new IOException("could not find file at " + path);
        }

        ImageExt leftImage = ImageIOHelper.readImageExt(path);
        ImageExt rightImage = ImageIOHelper.readImageExt(ResourceFinder.findTestResourcesDirectory() + sep + "bouguet_stereo"
                + sep + "right02_masked.png");

        CorrespondencePlotter plotter = new CorrespondencePlotter(leftImage, rightImage);

        left2 = getBouguetIm2LeftCorresDouble();
        right2 = getBouguetIm2RightCorresDouble();

        Corres corres = new Corres();
        corres.x1 = MatrixUtil.zeros(3, fit.getInlierIndexes().size());
        corres.x2 = MatrixUtil.zeros(3, fit.getInlierIndexes().size());
        Arrays.fill(corres.x1[2], 1);
        Arrays.fill(corres.x2[2], 1);

        int idx1, x1, y1, x2, y2, j;
        for (i = 0; i < fit.inlierIndexes.size(); ++i) {
            idx1 = fit.getInlierIndexes().get(i);
            x1 = (int)left2[0][idx1];
            y1 = (int)left2[1][idx1];
            x2 = (int)right2[0][idx1];
            y2 = (int)right2[1][idx1];
            plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 1);
            for (j = 0; j < 2; ++j) {
                corres.x1[j][i] = left2[j][idx1];
                corres.x2[j][i] = right2[j][idx1];
            }
        }
        plotter.writeImage("_bouguet_left2_right2_corres");

        return corres;
    }

    public static double[][] getBouguetIm2LeftCorresDouble() {
        double[][] x = new double[72][];
        x[0] = new double[]{206, 165, 1};
        x[1] = new double[]{213, 241, 1};
        x[2] = new double[]{218, 301, 1};
        x[3] = new double[]{221, 325, 1};
        x[4] = new double[]{224, 351, 1};
        x[5] = new double[]{253, 54, 1};
        x[6] = new double[]{252, 79, 1};
        x[7] = new double[]{249, 129, 1};
        x[8] = new double[]{252, 174, 1};
        x[9] = new double[]{252, 214, 1};
        x[10] = new double[]{252, 247, 1};
        x[11] = new double[]{254, 281, 1};
        x[12] = new double[]{253, 308, 1};
        x[13] = new double[]{256, 335, 1};
        x[14] = new double[]{254, 355, 1};
        x[15] = new double[]{256, 368, 1};
        x[16] = new double[]{306, 87, 1};
        x[17] = new double[]{304, 138, 1};
        x[18] = new double[]{302, 182, 1};
        x[19] = new double[]{298, 221, 1};
        x[20] = new double[]{297, 259, 1};
        x[21] = new double[]{294, 289, 1};
        x[22] = new double[]{295, 317, 1};
        x[23] = new double[]{292, 342, 1};
        x[24] = new double[]{290, 367, 1};
        x[25] = new double[]{292, 377, 1};
        x[26] = new double[]{369, 73, 1};
        x[27] = new double[]{366, 99, 1};
        x[28] = new double[]{358, 146, 1};
        x[29] = new double[]{353, 193, 1};
        x[30] = new double[]{348, 232, 1};
        x[31] = new double[]{343, 269, 1};
        x[32] = new double[]{339, 299, 1};
        x[33] = new double[]{336, 329, 1};
        x[34] = new double[]{330, 353, 1};
        x[35] = new double[]{326, 372, 1};
        x[36] = new double[]{326, 384, 1};
        x[37] = new double[]{428, 83, 1};
        x[38] = new double[]{424, 109, 1};
        x[39] = new double[]{414, 159, 1};
        x[40] = new double[]{405, 203, 1};
        x[41] = new double[]{396, 241, 1};
        x[42] = new double[]{389, 277, 1};
        x[43] = new double[]{382, 310, 1};
        x[44] = new double[]{373, 337, 1};
        x[45] = new double[]{370, 362, 1};
        x[46] = new double[]{363, 383, 1};
        x[47] = new double[]{362, 393, 1};
        x[48] = new double[]{485, 121, 1};
        x[49] = new double[]{469, 169, 1};
        x[50] = new double[]{457, 215, 1};
        x[51] = new double[]{443, 254, 1};
        x[52] = new double[]{435, 288, 1};
        x[53] = new double[]{427, 317, 1};
        x[54] = new double[]{416, 343, 1};
        x[55] = new double[]{406, 370, 1};
        x[56] = new double[]{540, 134, 1};
        x[57] = new double[]{524, 184, 1};
        x[58] = new double[]{510, 223, 1};
        x[59] = new double[]{493, 261, 1};
        x[60] = new double[]{478, 297, 1};
        x[61] = new double[]{468, 325, 1};
        x[62] = new double[]{454, 354, 1};
        x[63] = new double[]{445, 375, 1};
        x[64] = new double[]{605, 118, 1};
        x[65] = new double[]{591, 147, 1};
        x[66] = new double[]{553, 235, 1};
        x[67] = new double[]{539, 270, 1};
        x[68] = new double[]{521, 305, 1};
        x[69] = new double[]{509, 333, 1};
        x[70] = new double[]{486, 382, 1};
        x[71] = new double[]{468, 416, 1};
        // make it [3 X 72]
        x = MatrixUtil.transpose(x);
        return x;
    }

    public static double[][] getBouguetIm2RightCorresDouble() {
        double[][] x = new double[72][];
        x[0] = new double[]{43, 183, 1};
        x[1] = new double[]{62, 255, 1};
        x[2] = new double[]{81, 311, 1};
        x[3] = new double[]{89, 336, 1};
        x[4] = new double[]{96, 360, 1};
        x[5] = new double[]{57, 77, 1};
        x[6] = new double[]{64, 103, 1};
        x[7] = new double[]{72, 147, 1};
        x[8] = new double[]{77, 188, 1};
        x[9] = new double[]{90, 227, 1};
        x[10] = new double[]{95, 259, 1};
        x[11] = new double[]{106, 291, 1};
        x[12] = new double[]{115, 322, 1};
        x[13] = new double[]{119, 345, 1};
        x[14] = new double[]{124, 364, 1};
        x[15] = new double[]{130, 379, 1};
        x[16] = new double[]{109, 105, 1};
        x[17] = new double[]{115, 156, 1};
        x[18] = new double[]{119, 199, 1};
        x[19] = new double[]{127, 235, 1};
        x[20] = new double[]{134, 271, 1};
        x[21] = new double[]{144, 303, 1};
        x[22] = new double[]{147, 330, 1};
        x[23] = new double[]{151, 351, 1};
        x[24] = new double[]{158, 377, 1};
        x[25] = new double[]{160, 388, 1};
        x[26] = new double[]{154, 88, 1};
        x[27] = new double[]{158, 115, 1};
        x[28] = new double[]{160, 163, 1};
        x[29] = new double[]{167, 208, 1};
        x[30] = new double[]{170, 246, 1};
        x[31] = new double[]{177, 282, 1};
        x[32] = new double[]{182, 313, 1};
        x[33] = new double[]{183, 337, 1};
        x[34] = new double[]{190, 363, 1};
        x[35] = new double[]{190, 383, 1};
        x[36] = new double[]{194, 397, 1};
        x[37] = new double[]{209, 96, 1};
        x[38] = new double[]{209, 122, 1};
        x[39] = new double[]{210, 169, 1};
        x[40] = new double[]{215, 215, 1};
        x[41] = new double[]{218, 257, 1};
        x[42] = new double[]{218, 291, 1};
        x[43] = new double[]{223, 323, 1};
        x[44] = new double[]{222, 349, 1};
        x[45] = new double[]{224, 370, 1};
        x[46] = new double[]{229, 392, 1};
        x[47] = new double[]{227, 405, 1};
        x[48] = new double[]{270, 132, 1};
        x[49] = new double[]{265, 181, 1};
        x[50] = new double[]{264, 225, 1};
        x[51] = new double[]{264, 266, 1};
        x[52] = new double[]{266, 302, 1};
        x[53] = new double[]{263, 331, 1};
        x[54] = new double[]{264, 359, 1};
        x[55] = new double[]{265, 380, 1};
        x[56] = new double[]{327, 141, 1};
        x[57] = new double[]{322, 190, 1};
        x[58] = new double[]{317, 237, 1};
        x[59] = new double[]{315, 277, 1};
        x[60] = new double[]{309, 311, 1};
        x[61] = new double[]{305, 338, 1};
        x[62] = new double[]{303, 368, 1};
        x[63] = new double[]{304, 394, 1};
        x[64] = new double[]{395, 121, 1};
        x[65] = new double[]{386, 153, 1};
        x[66] = new double[]{368, 247, 1};
        x[67] = new double[]{361, 285, 1};
        x[68] = new double[]{354, 320, 1};
        x[69] = new double[]{349, 348, 1};
        x[70] = new double[]{339, 398, 1};
        x[71] = new double[]{336, 431, 1};
        // make it [3 X 72]
        x = MatrixUtil.transpose(x);
        return x;
    }

    public PairIntArray getBouguetIm2LeftCorresPairInt() {

        PairIntArray leftX = new PairIntArray();
        leftX.add(206, 165);
        leftX.add(213, 241);
        leftX.add(218, 301);
        leftX.add(221, 325);
        leftX.add(224, 351);
        leftX.add(249, 129);
        leftX.add(252, 79);
        leftX.add(252, 174);
        leftX.add(252, 214);
        leftX.add(252, 247);
        leftX.add(253, 54);
        leftX.add(253, 308);
        leftX.add(254, 281);
        leftX.add(254, 355);
        leftX.add(256, 335);
        leftX.add(256, 368);
        leftX.add(290, 367);
        leftX.add(292, 342);
        leftX.add(292, 377);
        leftX.add(294, 289);
        leftX.add(295, 317);
        leftX.add(297, 259);
        leftX.add(298, 221);
        leftX.add(302, 182);
        leftX.add(304, 138);
        leftX.add(306, 87);
        leftX.add(326, 372);
        leftX.add(326, 384);
        leftX.add(330, 353);
        leftX.add(336, 329);
        leftX.add(339, 299);
        leftX.add(343, 269);
        leftX.add(348, 232);
        leftX.add(353, 193);
        leftX.add(358, 146);
        leftX.add(362, 393);
        leftX.add(363, 383);
        leftX.add(366, 99);
        leftX.add(369, 73);
        leftX.add(370, 362);
        leftX.add(373, 337);
        leftX.add(382, 310);
        leftX.add(389, 277);
        leftX.add(396, 241);
        leftX.add(405, 203);
        leftX.add(406, 370);
        leftX.add(414, 159);
        leftX.add(416, 343);
        leftX.add(424, 109);
        leftX.add(427, 317);
        leftX.add(428, 83);
        leftX.add(435, 288);
        leftX.add(443, 254);
        leftX.add(445, 375);
        leftX.add(454, 354);
        leftX.add(457, 215);
        leftX.add(468, 325);
        leftX.add(468, 416);
        leftX.add(469, 169);
        leftX.add(478, 297);
        leftX.add(485, 121);
        leftX.add(486, 382);
        leftX.add(493, 261);
        leftX.add(509, 333);
        leftX.add(510, 223);
        leftX.add(521, 305);
        leftX.add(524, 184);
        leftX.add(539, 270);
        leftX.add(540, 134);
        leftX.add(553, 235);
        leftX.add(591, 147);
        leftX.add(605, 118);
        return leftX;
    }

    public PairIntArray getBouguetIm2RightCorresPairInt() {
        PairIntArray rightX = new PairIntArray();

        rightX.add(43, 183);
        rightX.add(57, 77);
        rightX.add(62, 255);
        rightX.add(64, 103);
        rightX.add(72, 147);
        rightX.add(77, 188);
        rightX.add(81, 311);
        rightX.add(89, 336);
        rightX.add(90, 227);
        rightX.add(95, 259);
        rightX.add(96, 360);
        rightX.add(106, 291);
        rightX.add(109, 105);
        rightX.add(115, 156);
        rightX.add(115, 322);
        rightX.add(119, 199);
        rightX.add(119, 345);
        rightX.add(124, 364);
        rightX.add(127, 235);
        rightX.add(130, 379);
        rightX.add(134, 271);
        rightX.add(144, 303);
        rightX.add(147, 330);
        rightX.add(151, 351);
        rightX.add(154, 88);
        rightX.add(158, 115);
        rightX.add(158, 377);
        rightX.add(160, 163);
        rightX.add(160, 388);
        rightX.add(167, 208);
        rightX.add(170, 246);
        rightX.add(177, 282);
        rightX.add(182, 313);
        rightX.add(183, 337);
        rightX.add(190, 363);
        rightX.add(190, 383);
        rightX.add(194, 397);
        rightX.add(209, 96);
        rightX.add(209, 122);
        rightX.add(210, 169);
        rightX.add(215, 215);
        rightX.add(218, 257);
        rightX.add(218, 291);
        rightX.add(222, 349);
        rightX.add(223, 323);
        rightX.add(224, 370);
        rightX.add(227, 405);
        rightX.add(229, 392);
        rightX.add(263, 331);
        rightX.add(264, 225);
        rightX.add(264, 266);
        rightX.add(264, 359);
        rightX.add(265, 181);
        rightX.add(265, 380);
        rightX.add(266, 302);
        rightX.add(270, 132);
        rightX.add(303, 368);
        rightX.add(304, 394);
        rightX.add(305, 338);
        rightX.add(309, 311);
        rightX.add(315, 277);
        rightX.add(317, 237);
        rightX.add(322, 190);
        rightX.add(327, 141);
        rightX.add(336, 431);
        rightX.add(339, 398);
        rightX.add(349, 348);
        rightX.add(354, 320);
        rightX.add(361, 285);
        rightX.add(368, 247);
        rightX.add(386, 153);
        rightX.add(395, 121);
        return rightX;
    }
}
