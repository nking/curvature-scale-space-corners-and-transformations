package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

import java.io.IOException;

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

        double[][] eFM = new double[3][];
        eFM[0] = new double[]{-0.217859, 0.419282, -0.0343075};
        eFM[1] = new double[]{-0.0717941, 0.0451643, 0.0216073};
        eFM[2] = new double[]{0.248062, -0.429478, 0.022101};
        eFM = MatrixUtil.transpose(eFM);

        //x2' * F * x1 = 0. // [N X 3]*[3 X 3]*[3 X N] = [N X N]
        double[][] check = MatrixUtil.multiply(MatrixUtil.pseudoinverseRankDeficient(x2),
                MatrixUtil.multiply(eFM, x1));
        int i;
        double[] ms;
        // "check"'s columns have means ~ 3E-3 and stdev ~ 0.1
  /*      for (i = 0; i < check.length; ++i) {
            ms = MiscMath0.getAvgAndStDev(check[i]);
            System.out.printf("r%d: %s\n", i, FormatArray.toString(ms, "%.3e"));
            ms = MiscMath0.getAvgAndStDev(MatrixUtil.extractColumn(check, i));
            System.out.printf("c%d: %s\n", i, FormatArray.toString(ms, "%.3e"));
        }
*/
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

        EpipolarTransformer tr = new EpipolarTransformer();
        double[][] fm0 = tr.calculateEpipolarProjection2(x1, x2, false);
        double[][] fmN = tr.calculateEpipolarProjection2(x1c, x2c, false);

        // fm0 resembles eFM, but scaled.
        double[][] fm0c = MatrixUtil.copy(fm0);
        double c = eFM[0][0]/fm0c[0][0];
        MatrixUtil.multiply(fm0c, c);

        //x2' * F * x1 = 0. // [N X 3]*[3 X 3]*[3 X N] = [N X N]
        double[][] checkN = MatrixUtil.multiply(MatrixUtil.pseudoinverseRankDeficient(x2c),
                MatrixUtil.multiply(fmN, x1c));
/*
        for (i = 0; i < checkN.length; ++i) {
            ms = MiscMath0.getAvgAndStDev(checkN[i]);
            System.out.printf("r%d: %s\n", i, FormatArray.toString(ms, "%.3e"));
            ms = MiscMath0.getAvgAndStDev(MatrixUtil.extractColumn(checkN, i));
            System.out.printf("c%d: %s\n", i, FormatArray.toString(ms, "%.3e"));
        }
*/
        double[][] fmD = MatrixUtil.copy(fmN);
        EpipolarNormalizationHelper.denormalizeFM(fmD, t1, t2);

        /*
        System.out.printf("fm0 =\n%s\n", FormatArray.toString(fm0, "%.6f"));
        System.out.printf("fmN =\n%s\n", FormatArray.toString(fmN, "%.6f"));
        System.out.printf("denormalized fm=\n%s\n", FormatArray.toString(fmD, "%.6f"));

        System.out.printf("%.3e * fm0 =\n%s\n", c, FormatArray.toString(fm0c, "%.6f"));
        System.out.printf("expected FM =\n%s\n", FormatArray.toString(eFM, "%.6f"));
        */

        int j;
        for (i = 0; i < eFM.length; ++i) {
            for (j = 0; j < eFM[0].length; ++j) {
                assertTrue(Math.abs(eFM[i][j] - fm0c[i][j]) < 1E-1);
            }
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
        /*System.out.printf("denorm x1 =\n%s\n", FormatArray.toString(x1cD, "%.6f"));
        System.out.printf("x1 =\n%s\n", FormatArray.toString(x1, "%.6f"));
        System.out.printf("denorm x2 =\n%s\n", FormatArray.toString(x2cD, "%.6f"));
        System.out.printf("x2 =\n%s\n", FormatArray.toString(x2, "%.6f"));
        */
    }
}
