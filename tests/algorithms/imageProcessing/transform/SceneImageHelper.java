package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;

import java.util.Arrays;

public class SceneImageHelper {

    public static double[][] createImagePoints2DPlanar(double[][] xW, double[][] k, double[][] r, double[] t) {

        if (xW.length != 3) {
            throw new IllegalArgumentException("expectign xW.length = 3");
        }
        int nFeatures = xW[0].length;

        double[][] xW4 = new double[4][xW[0].length];
        System.arraycopy(xW[0], 0, xW4[0], 0, nFeatures);
        System.arraycopy(xW[1], 0, xW4[1], 0, nFeatures);
        Arrays.fill(xW4[3], 1.);

        double[][] P4 = new double[4][4];
        for (int row = 0; row < 3; ++row) {
            System.arraycopy(r[row], 0, P4[row], 0, 3);
            P4[row][3] = t[row];
        }
        P4[3][3] = 1;
        double[][] x4 = MatrixUtil.multiply(P4, xW4);
        double[][] x3 = new double[3][nFeatures];
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < nFeatures; ++col) {
                x3[row][col] = x4[row][col]/x4[2][col];
            }
        }
        x3 = MatrixUtil.multiply(k, x3);

        return x3;
    }

    private static void scaleImagePoints2DPlanar(double[][] x, double[][] xEst) {

        int nFeatures = x[0].length;

        double[] meanScaleX = new double[x.length];
        for (int col = 0; col < nFeatures; ++col) {
            for (int row = 0; row < 3; ++row) {
                meanScaleX[row] += (x[row][col]/xEst[row][col]);
            }
        }
        for (int k = 0; k < meanScaleX.length; ++k) {
            meanScaleX[k] /= nFeatures;
        }
        double[] stdX = new double[3];
        double diff;
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < nFeatures; ++col) {
                diff = (x[row][col]/xEst[row][col]) - meanScaleX[row];
                stdX[row] += (diff*diff);
            }
        }
        for (int k = 0; k < stdX.length; ++k) {
            stdX[k] /= (nFeatures - 1.);
            stdX[k] = Math.sqrt(stdX[k]);
        }
        System.out.printf("scaleX=%s, stdX=%s\n", FormatArray.toString(meanScaleX, "%.2e"), FormatArray.toString(stdX, "%.2e"));

        // apply scale to current xEst
        for (int col = 0; col < nFeatures; ++col) {
            for (int row = 0; row < 3; ++row) {
                xEst[row][col] *= meanScaleX[0];
            }
        }
    }
}
