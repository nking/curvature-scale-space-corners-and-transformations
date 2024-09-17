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

         updated location:
         https://github.com/scikit-image/scikit-image/blob/main/skimage/transform/tests/test_geometric.py
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

        // expected essential matrix
        double[][] eEM = new double[3][];
        eEM[0] = new double[]{-0.0811666, 0.255449, -0.0478999};
        eEM[1] = new double[]{-0.192392, -0.0531675, 0.119547};
        eEM[2] = new double[]{0.177784, -0.22008, -0.015203};
        eEM = MatrixUtil.transpose(eEM);

        double[][] x1c = MatrixUtil.copy(x1);
        double[][] t1 = EpipolarNormalizationHelper.unitStandardNormalize(x1c);
        int i, j;
        double[] ms;
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

        // match with tolerance in case change to use scale std of (nDim-1) instead of 1.
        double[][] eFMN = new double[][] {
                {-0.6164932107421747, -0.24156493615453212, -0.18595938005677867},
            {0.7071491331679607, 0.09063380472040895, 0.10699368212780982},
            {0.06211595906961463, -0.05257531441824141, -0.02511815364621099}
        };

        double[][] eEMN = new double[][] {
                {-0.26126038710558697, -0.4012125745099643, -0.24919857450237587},
                {0.4946767953545696, -0.18522911621742388, -0.06691716621987827},
                {0.10676310049702824, -0.22582818112148015, -0.1204101792395726}
        };

        EpipolarTransformer tr = new EpipolarTransformer();
        double[][] fm0 = tr.calculateEpipolarProjection2(x1, x2, false);
        double[][] fmN = tr.calculateEpipolarProjection2(x1c, x2c, false);
        double[][] fmD = MatrixUtil.copy(fmN);
        EpipolarNormalizationHelper.denormalizeFM(fmD, t1, t2);

        double[][] em0 = tr.calculateEpipolarProjection2(x1, x2, true);
        double[][] emN = tr.calculateEpipolarProjection2(x1c, x2c, true);
        double[][] emD = MatrixUtil.copy(emN);
        EpipolarNormalizationHelper.denormalizeFM(emD, t1, t2);

        tol = 0.5;
        for (i = 0; i < 6; ++i) {
            double[][] c;
            double[][] e;
            if (i < 3) {
                e = eFM;
            } else {
                e = eEM;
            }
            if (i == 1 || i == 4) continue;
            switch(i) {
                case 0 : c = fm0; break;
                case 1 : c = fmN; break;
                case 2 : c = fmD; break;
                case 3 : c = em0; break;
                case 4 : c = emN; break;
                default : c = emD; break;
            }
            double[][] diff = MatrixUtil.pointwiseSubtract(e, c);
            double fs = MatrixUtil.frobeniusNorm(diff);
            //System.out.printf("%d) fs=%.4e\n", i, fs);
            assertTrue(fs < tol);
        }

        // check epipolar lines

        // l2 = F * x1
        // l1 = F^T * x2
        // x2^T * F * x1 = 0
        // e2^T * F = 0
        // F * e1 = 0
        // l2 = F * x1
        // l1 = F^T * x2

        tol = 15;

        for (i = 0; i < 6; ++i) {
            if (i == 1 || i == 4) continue;

            double[][] _x1, _x2;
            double[][] _fOrEM;
            switch(i) {
                case 0 : _x1=x1; _x2=x2; _fOrEM=fm0; break;
                case 1 : _x1=x1c; _x2=x2c; _fOrEM=fmN; break;
                case 2 : _x1=x1; _x2=x2; _fOrEM=fmD; break;
                case 3 : _x1=x1; _x2=x2; _fOrEM=em0; break;
                case 4 : _x1=x1c; _x2=x2c; _fOrEM=emN; break;
                default : _x1=x1; _x2=x2; _fOrEM=emD; break;
            }

            // 3x3  3xN = 3XN
            double[][] l1 = MatrixUtil.multiply( MatrixUtil.transpose(_fOrEM), _x2);
            double[][] l2 = MatrixUtil.multiply( _fOrEM, _x1);

            // e1 = e1e2[0]
            // e2 = e1e2[1]
            double[][] e1e2 = EpipolarTransformer.calculateEpipoles(new DenseMatrix(_fOrEM));

            // x1^T * l1 = 0   which is x1^T * F^T * x2 = 0
            // x2^T * l2 = 0   which is x2^T * F * x1 = 0 and same as (x1^T * F^T * x2)^T
            // e2^T * F = 0  <=== alot faster to look at and gives same result

            // these test data points might have outliers in them.

            // the matrices which were calculated with unit standard points, have x1TL1 and e2TF closer to 0

            double[][] x1TL1 = MatrixUtil.multiply( MatrixUtil.transpose(_x1), l1);
            //double[][] x2TL2 = MatrixUtil.multiply( MatrixUtil.transpose(_x2), l2);
            double[] e2TF = MatrixUtil.multiplyRowVectorByMatrix(e1e2[0], _fOrEM);

            double fs = MatrixUtil.frobeniusNorm(x1TL1);
            //System.out.printf("%d) x1TL1 fs=%.4e\n", i, fs);

            double norm = MatrixUtil.lPSum(e2TF, 2);
            //System.out.printf("%d) e2TF norm=%.4e\n", i, fs);

            for (double[] row : x1TL1) {
                for (double element : row) {
                    assertTrue(Math.abs(element) < tol);
                }
            }

            for (double element : e2TF) {
                assertTrue(Math.abs(element) < tol);
            }
        }

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

        //F = K2^-T * E * K1^-1
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
        double[][] fmN = tr.calculateEpipolarProjection2(x1c, x2c, false);
        double[][] fmD = MatrixUtil.copy(fmN);
        EpipolarNormalizationHelper.denormalizeFM(fmD, t1, t2);

        double[][] em0 = tr.calculateEpipolarProjection2(xr1, xr2, true);
        double[][] emN = tr.calculateEpipolarProjection2(xr1, xr2, true);
        double[][] emD = MatrixUtil.copy(emN);
        EpipolarNormalizationHelper.denormalizeFM(emD, t1, t2);

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
