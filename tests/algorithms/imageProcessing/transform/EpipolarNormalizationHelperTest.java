package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

import java.io.IOException;
import java.util.Arrays;

public class EpipolarNormalizationHelperTest extends TestCase {

    public void test0() throws IOException, NotConvergedException {
        /*
        these numbers are from SciPy unit test
         test_fundamental_matrix_estimation() in file
         https://github.com/scikit-image/scikit-image/transform/tests/test_geometric.py

         the Scipy scikit-image github repository posts a BSD Clause-3 license
         https://github.com/scikit-image/scikit-image/blob/025757f8d5a1ece2687781000560aad79047b6c3/LICENSE.txt

         They reference the "COLMAP SfM library".
         */
        /*
        src = np.array([1.839035, 1.924743, 0.543582,  0.375221,
                0.473240, 0.142522, 0.964910,  0.598376,
                0.102388, 0.140092, 15.994343, 9.622164,
                0.285901, 0.430055, 0.091150,  0.254594]).reshape(-1, 2)
        dst = np.array([1.002114, 1.129644, 1.521742, 1.846002,
                1.084332, 0.275134, 0.293328, 0.588992,
                0.839509, 0.087290, 1.779735, 1.116857,
                0.878616, 0.602447, 0.642616, 1.028681]).reshape(-1, 2)
        */

        double[][] x1 = new double[8][];
        x1[0] = new double[]{1.839035, 1.924743, 1};
        x1[1] = new double[]{0.543582,  0.375221, 1};
        x1[2] = new double[]{0.473240, 0.142522,  1};
        x1[3] = new double[]{0.964910,  0.598376, 1};
        x1[4] = new double[]{0.102388, 0.140092, 1};
        x1[5] = new double[]{15.994343, 9.622164, 1};
        x1[6] = new double[]{0.285901, 0.430055, 1};
        x1[7] = new double[]{0.091150,  0.254594, 1};
        x1 = MatrixUtil.transpose(x1);
        // [3XN], N=8

        double[][] x2 = new double[8][];
        x2[0] = new double[]{1.002114, 1.129644, 1};
        x2[1] = new double[]{1.521742, 1.846002, 1};
        x2[2] = new double[]{1.084332, 0.275134, 1};
        x2[3] = new double[]{0.293328, 0.588992,1};
        x2[4] = new double[]{0.839509, 0.087290, 1};
        x2[5] = new double[]{1.779735, 1.116857,1};
        x2[6] = new double[]{0.878616, 0.602447, 1};
        x2[7] = new double[]{0.642616, 1.028681, 1};
        x2 = MatrixUtil.transpose(x2);
        // [3XN], N=8

        // expected fundamental matrix
        double[][] eFM = new double[3][];
        eFM[0] = new double[]{-0.217859, 0.419282, -0.0343075};
        eFM[1] = new double[]{-0.0717941, 0.0451643, 0.0216073};
        eFM[2] = new double[]{0.248062, -0.429478, 0.022101};
        eFM = MatrixUtil.transpose(eFM);
        // [3X3]

        //x2'^T * F * x1 = 0. // [N X 3]*[3 X 3]*[3 X N] = [N X N]
        double[][] check = MatrixUtil.multiply(MatrixUtil.transpose(x2),
                MatrixUtil.multiply(eFM, x1));
        int i;
        double[] ms;
        // "check"'s columns have means ~ 3E-3 and stdev ~ 0.1
        for (i = 0; i < check.length; ++i) {
            ms = MiscMath0.getAvgAndStDev(check[i]);
            System.out.printf("row %d: mean, stdev=%s\n", i, FormatArray.toString(ms, "%.3e"));
            ms = MiscMath0.getAvgAndStDev(MatrixUtil.extractColumn(check, i));
            System.out.printf("col %d: mean, stdev=%s\n", i, FormatArray.toString(ms, "%.3e"));
        }

        // expected essential matrix
        double[][] eEM = new double[3][];
        eEM[0] = new double[]{-0.0811666, 0.255449, -0.0478999};
        eEM[1] = new double[]{-0.192392, -0.0531675, 0.119547};
        eEM[2] = new double[]{0.177784, -0.22008, -0.015203};
        eEM = MatrixUtil.transpose(eEM);

        double[][] x1c = MatrixUtil.copy(x1);
        double[][] t1 = EpipolarNormalizationHelper.unitStandardNormalize(x1c);
        double tol = 1E-7;
        // assert the mean and stdev of x1c by rows
        for (i = 0; i < 2; ++i) {
            ms = MiscMath0.getAvgAndStDev(x1c[i]);
            assertTrue(Math.abs(ms[0]) < tol);
            assertTrue(Math.abs(ms[1] - 1.) < tol);
        }
        double[][] x2c = MatrixUtil.copy(x2);
        double[][] t2 = EpipolarNormalizationHelper.unitStandardNormalize(x2c);
        for (i = 0; i < 2; ++i) {
            ms = MiscMath0.getAvgAndStDev(x2c[i]);
            assertTrue(Math.abs(ms[0]) < tol);
            assertTrue(Math.abs(ms[1] - 1.) < tol);
        }

        double[][] x1cD = MatrixUtil.copy(x1c);
        double[][] x2cD = MatrixUtil.copy(x2c);
        EpipolarNormalizationHelper.denormalize(x1cD, x2cD, t1, t2);
        int j;
        for (i = 0; i < x1cD.length; ++i) {
            for (j = 0; j < x1cD[0].length; ++j) {
                assertTrue(Math.abs(x1[i][j] - x1cD[i][j]) < 1E-7);
            }
        }
        for (i = 0; i < x2cD.length; ++i) {
            for (j = 0; j < x2cD[0].length; ++j) {
                assertTrue(Math.abs(x2[i][j] - x2cD[i][j]) < 1E-7);
            }
        }

        EpipolarTransformer tr = new EpipolarTransformer();
        double[][] fm0 = tr.calculateEpipolarProjection2(x1, x2, false);
        double[][] fmN = tr.calculateEpipolarProjection2(x1c, x2c, false);
        double[][] fmD = MatrixUtil.copy(fmN);
        EpipolarNormalizationHelper.denormalizeFM(fmD, t1, t2);

        // compare the fms to eFM, but normalizes them all by element[2][2] first
        double[][] fm0c = MatrixUtil.copy(fm0);
        MatrixUtil.multiply(fm0c, 1./fm0c[1][0]);
        double[][] fmNc = MatrixUtil.copy(fmN);
        MatrixUtil.multiply(fmNc, 1./fmNc[1][0]);
        double[][] fmDc = MatrixUtil.copy(fmD);
        MatrixUtil.multiply(fmDc, 1./fmDc[1][0]);

        double[][] efmc = MatrixUtil.copy(eFM);
        MatrixUtil.multiply(efmc, 1./efmc[1][0]);

        // efmc matches the solution where we did not normalize the data x1 and x2

        double[][] diff = MatrixUtil.pointwiseSubtract(efmc, fm0c);
        double fs = MatrixUtil.frobeniusNorm(diff);
        assertTrue(fs < 0.1);

        //-----

        double[][] em0 = tr.calculateEpipolarProjection2(x1, x2, true);
        double[][] emN = tr.calculateEpipolarProjection2(x1c, x2c, true);

        double[][] emD = MatrixUtil.copy(emN);

        EpipolarNormalizationHelper.denormalizeFM(emD, t1, t2);

        // compare the fms to eFM, but normalizes them all by element[2][2] first
        double[][] em0c = MatrixUtil.copy(em0);
        MatrixUtil.multiply(em0c, 1./em0c[1][0]);
        double[][] emNc = MatrixUtil.copy(emN);
        MatrixUtil.multiply(emNc, 1./emNc[1][0]);
        double[][] emDc = MatrixUtil.copy(emD);
        MatrixUtil.multiply(emDc, 1./emDc[1][0]);

        double[][] eemc = MatrixUtil.copy(eEM);
        MatrixUtil.multiply(eemc, 1./eemc[1][0]);

        diff = MatrixUtil.pointwiseSubtract(eemc, em0c);
        fs = MatrixUtil.frobeniusNorm(diff);
        assertTrue(fs < 1);

        // in the absence of noise checkN is 0.

        //x2^T * E * x1 = 0. // [N X 3]*[3 X 3]*[3 X N] = [N X N]
        double[][] checkN = MatrixUtil.multiply(MatrixUtil.transpose(x2c),
                MatrixUtil.multiply(emN, x1c));

        for (i = 0; i < checkN.length; ++i) {
            ms = MiscMath0.getAvgAndStDev(checkN[i]);
            System.out.printf("row %d: mean, stdev = %s\n", i, FormatArray.toString(ms, "%.3e"));
            ms = MiscMath0.getAvgAndStDev(MatrixUtil.extractColumn(checkN, i));
            System.out.printf("col %d: mean, stdev= %s\n", i, FormatArray.toString(ms, "%.3e"));
        }

        /*
        System.out.printf("em0 =\n%s\n", FormatArray.toString(em0, "%.6f"));
        System.out.printf("expected EM =\n%s\n", FormatArray.toString(eEM, "%.6f"));
        System.out.printf("em0c =\n%s\n", FormatArray.toString(em0c, "%.6f"));
        */

        /*System.out.printf("denorm x1 =\n%s\n", FormatArray.toString(x1cD, "%.6f"));
        System.out.printf("x1 =\n%s\n", FormatArray.toString(x1, "%.6f"));
        System.out.printf("denorm x2 =\n%s\n", FormatArray.toString(x2cD, "%.6f"));
        System.out.printf("x2 =\n%s\n", FormatArray.toString(x2, "%.6f"));
        */

    }

    public void test2() throws IOException, NotConvergedException {

        boolean passive = true;

        // add masks reference for example test
        /*double[][] XW = new double[][]{
                {0, 1, 1, 0, 0, 1, 1, 0, 0.2, 0.8, 0.2, 0.8},
                {0, 0, 1, 1, 0, 0, 1, 1, 1.5, 1.5, 1.5, 1.5},
                {1, 1, 1, 1, 0, 0, 0, 0, 0.8, 0.8, 0.2, 0.2},
                {1, 1, 1, 1, 1, 1, 1, 1, 1,   1,   1,   1}
        };*/
        // only want 8 points
        double[][] XW = new double[][]{
                {0, 1, /*1,*/ 0, 0, /*1,*/ 1, 0, 0.2, /*0.8,*/ 0.2, 0.8},
                {0, 0, /*1,*/ 1, 0, /*0,*/ 1, 1, 1.5, /*1.5,*/ 1.5, 1.5},
                {1, 1, /*1,*/ 1, 0, /*0,*/ 0, 0, 0.8, /*0.8,*/ 0.2, 0.2},
                {1, 1, /*1,*/ 1, 1, /*1,*/ 1, 1, 1,   /*1,*/   1,   1}
        };

        /*double[][] Rinit = rot_matrix([1 1 1],0);
        Zinit = 5;
        Pinit = [ Rinit(1,:) 0 ;
                  Rinit(2,:) 0 ;
                 Rinit(3,:) Zinit;
                 0 0 0 1];
        XC(:,:,1) = Pinit*XW;
        */
        double[][] rInit = masksRotationMatrix(new double[]{1,1,1}, 0, passive);
        // same as using rotVec = new double[]{0,0,0};
        double zInit = 5;
        double[][] pInit = new double[4][4];
        for (int ii = 0; ii < 3; ++ii) {
            System.arraycopy(rInit[ii], 0, pInit[ii], 0, 3);
        }
        pInit[2][3] = zInit;
        pInit[3][3] = 1;

        double[][] xC1 = MatrixUtil.multiply(pInit, XW);

        double[][] xr1 = MatrixUtil.copySubMatrix(xC1, 0, 2, 0, xC1[0].length-1);
        int row, col;
        for (row = 0; row < xr1.length; ++row) {
            for (col = 0; col < xr1[row].length; ++col) {
                xr1[row][col] /= xr1[2][col];
            }
        }

        double[][] K = new double[][]{{600, 0, 300},
                {0, 600, 300},
                {0, 0, 1}};

        double[][] xIm1 = MatrixUtil.multiply(K, xr1);


        double[] ax = new double[]{0, 1, 0};
        double[] trans12 = new double[]{1,0,1};
        double angle = -20;
        double[] rot_axis = Arrays.copyOf(ax, ax.length);
        MatrixUtil.multiply(rot_axis, 1./MatrixUtil.lPSum(rot_axis, 2));
        double theta = (angle)*Math.PI/180;

        double[][] r12 = masksRotationMatrix(rot_axis, theta, passive);
        // same as using rotVec12 = new double[]{0,theta,0};
        //double[][] r12 = Rotation.createRotationRodriguesFormula(rot_axis, angle, passive);
        double[][] p12 = new double[4][4];
        for (int ii = 0; ii < 3; ++ii) {
            System.arraycopy(r12[ii], 0, p12[ii], 0, 3);
            p12[ii][3] = trans12[ii];
        }
        p12[3][3] = 1;

        double[][] xC2 = MatrixUtil.multiply(p12, xC1);
        double[][] xr2 = MatrixUtil.copySubMatrix(xC2, 0, 2, 0, xC2[0].length-1);
        for (row = 0; row < xr2.length; ++row) {
            for (col = 0; col < xr2[row].length; ++col) {
                xr2[row][col] /= xr2[2][col];
            }
        }

        double[][] xIm2 = MatrixUtil.multiply(K, xr2);

        //E = [T]_x * R
        double[][] skewSymT = MatrixUtil.skewSymmetric(trans12);
        double[][] eEM = MatrixUtil.multiply(skewSymT, r12);

        double[] skewTT = MatrixUtil.multiplyMatrixByColumnVector(skewSymT, trans12);
        assertTrue(MiscMath.areEqual(new double[]{0, 0, 0}, skewTT, 1E-11));

        //F = K^-T * E * K^-1
        double[][] kInv = Camera.createIntrinsicCameraMatrixInverse(K);
        double[][] eFM = MatrixUtil.multiply(
                MatrixUtil.multiply(MatrixUtil.transpose(kInv), eEM), kInv);
        MatrixUtil.multiply(eFM, 1./eFM[2][2]);

        //lambda1 * [T]_x * x2 = lambda2 * [T]_x * R * x1
        //lambda1 * [T]_x * x2 = lambda2 * E * x1
        // x2 = (lambda2/lambda1) * ([T]_x)^-1 * E * x1 = (lambda2/lambda1) * R * x1
        //double[][] _xIm2 = MatrixUtil.multiply(r12, xIm1);

        //lambda1 * [T]_x * x2 = lambda2 * [T]_x * R * x1 + [T]_x * T
        //lambda1 * [T]_x * x2 = lambda2 * E * x1 + [T]_x * T
        //double[][] lhs = MatrixUtil.multiply(skewSymT, xIm2);
        //double[][] rhs = MatrixUtil.multiply(eEM, xIm1);

        double[][] x1c = MatrixUtil.copy(xIm1);
        double[][] t1 = EpipolarNormalizationHelper.unitStandardNormalize(x1c);
        double tol = 1E-7;

        double[][] x2c = MatrixUtil.copy(xIm2);
        double[][] t2 = EpipolarNormalizationHelper.unitStandardNormalize(x2c);

        double[][] x1cD = MatrixUtil.copy(x1c);
        double[][] x2cD = MatrixUtil.copy(x2c);
        EpipolarNormalizationHelper.denormalize(x1cD, x2cD, t1, t2);
        int i,j;
        for (i = 0; i < x1cD.length; ++i) {
            for (j = 0; j < x1cD[0].length; ++j) {
                assertTrue(Math.abs(xIm1[i][j] - x1cD[i][j]) < 1E-7);
            }
        }
        for (i = 0; i < x2cD.length; ++i) {
            for (j = 0; j < x2cD[0].length; ++j) {
                assertTrue(Math.abs(xIm2[i][j] - x2cD[i][j]) < 1E-7);
            }
        }

        EpipolarTransformer tr = new EpipolarTransformer();
        double[][] fm0 = tr.calculateEpipolarProjection2(xIm1, xIm2, false);
        MatrixUtil.multiply(fm0, 1./fm0[2][2]);
        double[][] fmN = tr.calculateEpipolarProjection2(x1c, x2c, false);
        double[][] fmD = MatrixUtil.copy(fmN);
        EpipolarNormalizationHelper.denormalizeFM(fmD, t1, t2);
        MatrixUtil.multiply(fmD, 1./fmD[2][2]);

        double[][] em0 = tr.calculateEpipolarProjection2(xr1, xr2, true);
        MatrixUtil.multiply(em0, 1./em0[2][2]);
        double[][] emN = tr.calculateEpipolarProjection2(xr1, xr2, true);
        double[][] emD = MatrixUtil.copy(emN);
        EpipolarNormalizationHelper.denormalizeFM(emD, t1, t2);
        MatrixUtil.multiply(emD, 1./emD[2][2]);

        // e1 is right null space (from V)
        // e2 is left null space (from U)
        double[][] e1e2 = EpipolarTransformer.calculateEpipoles(new DenseMatrix(fmD));

        // M as a Rotation
        double[][] M = MatrixUtil.multiply(MatrixUtil.transpose(MatrixUtil.skewSymmetric(e1e2[0])), fmD);
        double[][] P1 = new double[][]{
                {1, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 1, 0}
        };
        double[][] P2 = new double[3][4];
        for (int ii = 0; ii < 3; ++ii) {
            System.arraycopy(M[ii], 0, P2[ii], 0, 3);
            P2[ii][3] = e1e2[0][ii];
        }

        EVD evd = EVD.factorize(new DenseMatrix(M));

        int t = 2;

    }
    static double[][] masksRotationMatrix(double[] omega, double theta, boolean passive) {
        double[][] omegaHat = MatrixUtil.skewSymmetric(omega);
        double normOmega = MatrixUtil.lPSum(omega, 2);
        double[][] eye = MatrixUtil.createIdentityMatrix(3);
        if (Math.abs(normOmega) > 1E-7) {

            double[][] oh1 = MatrixUtil.copy(omegaHat);
            MatrixUtil.multiply(oh1, 1./normOmega);

            double[][] oh2 = MatrixUtil.multiply(omegaHat, omegaHat);
            MatrixUtil.multiply(oh2, 1./(normOmega*normOmega));

            double s = Math.sin(normOmega*theta);
            if (!passive) {
                s *= -1;
            }
            double c = 1 - Math.cos(normOmega*theta);

            MatrixUtil.multiply(oh1, s);
            MatrixUtil.multiply(oh2, c);

            return MatrixUtil.pointwiseAdd(eye,
                    MatrixUtil.pointwiseAdd(oh1, oh2)
            );
        } else {
            return eye;
        }
    }
}
