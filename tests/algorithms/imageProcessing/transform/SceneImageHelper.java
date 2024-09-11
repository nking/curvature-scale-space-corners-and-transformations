package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;

import java.util.Arrays;
import java.util.Random;

public class SceneImageHelper {

    public static double[][] createImagePoints2D(double[][] xW, double[][] k, double[][] r, double[] t,
                                                 boolean makeHomogenous) {
        int nFeatures = xW[0].length;
        double[][] xW4;
        if (xW.length == 2 || xW.length == 3) {
            xW4 = new double[4][xW[0].length];
            System.arraycopy(xW[0], 0, xW4[0], 0, nFeatures);
            System.arraycopy(xW[1], 0, xW4[1], 0, nFeatures);
            if (xW.length == 3) {
                //TODO: consider whether need to normalize xW[0] and xW[1] by this row
                System.arraycopy(xW[2], 0, xW4[2], 0, nFeatures);
            }
            Arrays.fill(xW4[3], 1.);
        } else if (xW.length == 4) {
            xW4 = MatrixUtil.copy(xW);
            for (int row = 0; row < 3; ++row) {
                for (int col = 0; col < nFeatures; ++col) {
                    xW4[row][col] /= xW4[row][3];
                }
            }
        } else {
            throw new IllegalArgumentException("xW.length must be 2, 3 or 4");
        }

        double[][] P4 = new double[4][4];
        for (int row = 0; row < 3; ++row) {
            System.arraycopy(r[row], 0, P4[row], 0, 3);
            P4[row][3] = t[row];
        }
        P4[3][3] = 1;
        double[][] x4 = MatrixUtil.multiply(P4, xW4);

        // last row should be all 1's
        assert(Math.abs(x4[3][0] - 1.) < 1E-7);

        // normalize if makeHomogenous
        double[][] x3 = new double[3][nFeatures];
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < nFeatures; ++col) {
                x3[row][col] = x4[row][col];
                if (makeHomogenous) {
                    x3[row][col] /= x4[2][col];
                }
            }
        }
        // convert camera coordinates to pixel coordinates
        x3 = MatrixUtil.multiply(k, x3);
        return x3;
    }
    public static double[][] createImagePoints2DPlanar(double[][] xW, double[][] k, double[][] r, double[] t) {

        if (xW.length != 3) {
            throw new IllegalArgumentException("expectign xW.length = 3");
        }
        int nFeatures = xW[0].length;

        double[][] xW4 = new double[4][xW[0].length];
        System.arraycopy(xW[0], 0, xW4[0], 0, nFeatures);
        System.arraycopy(xW[1], 0, xW4[1], 0, nFeatures);
        // row 2 is all 0s for planar solution
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

    /**
     *
     * @param rand
     * @param nPoints
     * @param xyzRange pairs of ranges for x, y, z.
     *                 in format new double[][] { {xmin, xmx}, {ymin, ymax}, {zmin, zmax}}
     * @return
     */
    public static double[][] randomPoints3D(Random rand, int nPoints, double[][] xyzRange) {
        if (xyzRange.length != 3) {
            throw new IllegalArgumentException("xyzRange.length should be 3");
        }
        if (xyzRange[0].length != 2) {
            throw new IllegalArgumentException("xyzRange[0].length should be 2");
        }
        //TODO: improve this
        double[] halfRange = new double[]{(xyzRange[0][1] - xyzRange[0][0])/2., (xyzRange[1][1] - xyzRange[1][0])/2.,
                (xyzRange[2][1] - xyzRange[2][0])/2.};
        double[][] xW = new double[4][nPoints];
        Arrays.fill(xW[3], 1);
        double sign;
        for (int i = 0; i < nPoints; ++i) {
            for (int j = 0; j < 3; ++j) {
                sign = rand.nextBoolean() ? -1 : 1;
                xW[j][i] = halfRange[j] + sign*(rand.nextDouble() * halfRange[j]);
            }
        }
        return xW;
    }
}
