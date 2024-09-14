package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import no.uib.cipr.matrix.NotConvergedException;

public class EpipolarNormalizationHelper {

    /*
    MASKS eqn (6.79)
    x1n = x coords that have been normalized to a mean of 0 and standard deviation of 1.
    x2n = y coords normalized similarly.
    let H1 and H2 be their respective normalization transformation matrices.

    x1n = H1 * x1
    x2n = H2 * x2

    x2^T * F * x1 = x2nT * H2^(-T) * F * H1^(-1) * x1n = 0

    when calculating the fundamental matrix using normalized coordinates,
    one is actually calculation Fn = H2^(-T) * F * H1^(-1)

    so to de-normalize F, we have F = H2^T * Fn * H1

    x2^T = x2n^T * H2^-T
    x1 = H1^-1 * x1n
     */

    /*
    <pre>
        Revisiting Hartley’s Normalized Eight-Point Algorithm
        Chojnacki et al. 2003

        inverse changes the order of operations,
        inverse translation matrix: inverse changes the signs of the translation elements,
                but not the diagonal.
        inverse rotation matrix: inverse is the transpose of rotation matrix.
        inverse scaling matrix: inverse is performed on each element, that is, the reciprocal.

        denormalized x2^T = transpose(normalized x2) * transpose(T2^1)
        denormalized x1 = (T1^-1) * (normalized x1)
        denormalized FM = transpose(T2) * FM_normalized * T1
         */

    /**
     * given homogenous coordinate data x in row major format [nDimension X nPoints],
     * perform unit standard normalization to modify
     * x to have mean = 0 and standard deviation = 1.
     * Note that the last row of x is assumed to be '1's and is not modified.
     * Note also that if the stdev of the zero mean subtracted data is 0 for any row, a scale
     * factor of 1 is used for it and a warning is printed to stdout.
     * the transformation matrix internally already used on x, is returned.
     * <pre>
     *     e.g. the transformation matrix for [3 X nPoint] input x:
     *     double[][] t = new double[3][];
     *         t[0] = new double[]{1./scale,       0,     -centroidX/scale};
     *         t[1] = new double[]{0,           1./scale, -centroidY/scale};
     *         t[2] = new double[]{0,           0,           1};
     *
     *     to denormalize the points and the fundamental matrix:
     *
     *     denormalized x2^T = transpose(normalized x2) * transpose(T2^1)
     *     denormalized x1 = (T1^-1) * (normalized x1)
     *     denormalized FM = transpose(T2) * FM_normalized * T1
     *
     *     where T1 is the transformation matrix from points x1
     *     and T2 is the transformation matrix from points x2.
     *
     * </pre>
     * @param x matrix of homogenous coordinate data in format [nDimension X nPoints] where row 0 has the
     *          'x-axis' coordinates, row 1 has the 'y-axis' coordinates,
     *          and if row 3 is present it is the 'z-axis' coordinates.
     *          The last row of data is assumed to be 1's and is not modified nor included in the output results.
     *
     * @return transformation matrix internally used to normalize x.
     <pre>
    e.g. for x input [3 X nPoints]:
    t[0] = new double[]{1./scale,       0,     -centroidX/scale};
    t[1] = new double[]{0,           1./scale, -centroidY/scale};
    t[2] = new double[]{0,           0,           1};
    </pre>
     */
    public static double[][] unitStandardNormalize(double[][] x) {

        int n = x[0].length;
        int nd = x.length - 1;

        double[] rowMeans = MatrixUtil.rowMeans(x);

        double[] stdevs = new double[nd];
        double diff;
        int row, col;
        for (row = 0; row < nd; ++row) {
            for (col = 0; col < n; ++col) {
                // difference from mean 0
                diff = x[row][col] - rowMeans[row];
                stdevs[row] += (diff * diff);
            }
        }
        for (row = 0; row < nd; ++row) {
            stdevs[row] = Math.sqrt(stdevs[row]/(n-1.));
            if (stdevs[row] == 0.0) {
                System.out.println("WARNING: standard deviation of row " + row
                        + " was 0, so consider using another normalization " +
                        "method like min-max instead");
            }
        }

        double[][] t = new double[x.length][x.length];
        for (row = 0; row < nd; ++row) {
            if (stdevs[row] == 0.) {
                t[row][row] = 1.;
                t[row][nd] = -rowMeans[row];
            } else {
                t[row][row] = 1. / stdevs[row];
                t[row][nd] = -rowMeans[row] / stdevs[row];
            }
        }
        t[nd][nd] = 1;

        double[][] xTrans = MatrixUtil.multiply(t, x);
        for (row = 0; row < xTrans.length; ++row) {
            System.arraycopy(xTrans[row], 0, x[row], 0, xTrans[row].length);
        }

        return t;
    }

    /**
     * denormalized x2^T = transpose(normalized x2) * transpose(T2^1)
     * denormalized x1 = (T1^-1) * (normalized x1)
     *
     * @param x1 a.k.a. xLeft are data that have been normalized using the given
     *      transformation matrix t1.
     *      x1 format is [3 X nPoints] where row 0 are the
     *      'x-axis' coordinates, row 1 are the 'y-axis' coordinates,
     *      and row 3 are the 'z-axis' coordinates which are expected
     *      to be all '1' (= homogenous coordinates).
     * @param x2 a.k.a. xRight are data that have been normalized using the given
     *     transformation matrix t1.
     *     x1 format is [3 X nPoints] where row 0 are the
     *     'x-axis' coordinates, row 1 are the 'y-axis' coordinates,
     *     and row 3 are the 'z-axis' coordinates which are expected
     *     to be all '1' (= homogenous coordinates).
     * @param t1 transformation matrix for x1. dimensions are [3 x 3]
     * @param t2 transformation matrix for x2. dimensions are [3 x 3]
     */
    public static void denormalize(double[][] x1, double[][] x2, double[][] t1,
                                   double[][] t2) {
        if (x1.length != 3) {
            throw new IllegalArgumentException("x1.length must be 3");
        }
        if (x2.length != 3) {
            throw new IllegalArgumentException("x2.length must be 3");
        }
        int n = x1[0].length;
        if (x2[0].length != n) {
            throw new IllegalArgumentException("x2.length must be the same as x1.length");
        }
        if (t1.length != 3 || t1[0].length != 3) {
            throw new IllegalArgumentException("t1 must be [3 X 3]");
        }
        if (t2.length != 3 || t2[0].length != 3) {
            throw new IllegalArgumentException("t2 must be [3 X 3]");
        }

        /*
        u1_normalized = T1 * u1
          u2_normalized = T2 * u2

          denormalized u1 = T1^-1 * u1_normalized
          denormalized u2 = T2^-1 * u2_normalized

          FM_normalized = inverse(transpose(T2)) * FM * inverse(T1)
         */
        double[][] tInv1 = inverseT(t1);
        double[][] tInv2 = inverseT(t2);

        // [NX3]*[3X3]=[NX3]
        double[][] x2D = MatrixUtil.multiply(tInv2, x2);
        double[][] x1D = MatrixUtil.multiply(tInv1, x1);

        int i;
        for (i = 0; i < 3; ++i) {
            System.arraycopy(x1D[i], 0, x1[i], 0, x1D[i].length);
        }
        for (i = 0; i < 3; ++i) {
            System.arraycopy(x2D[i], 0, x2[i], 0, x2D[i].length);
        }
    }

    /**
     * FM_normalized = inverse(transpose(T2)) * FM * inverse(T1)
     * denormalized FM = transpose(T2) * FM_normalized * T1
     *
     * @param fm the fundamental matrix calculated using normalized x1 and x2 data
     *           points, where t1 and t2 are the transformation matrices used
     *           to normalize x1 and x2, respectively.
     * @param t1 transformation matrix for x1. dimensions are [3 x 3]
     * @param t2 transformation matrix for x2. dimensions are [3 x 3]
     */
    public static void denormalizeFM(double[][] fm, double[][] t1, double[][] t2) {

        if (fm.length != 3 || fm[0].length != 3) {
            throw new IllegalArgumentException("fm must be [3 X 3]");
        }
        if (t1.length != 3 || t1[0].length != 3) {
            throw new IllegalArgumentException("t1 must be [3 X 3]");
        }
        if (t2.length != 3 || t2[0].length != 3) {
            throw new IllegalArgumentException("t2 must be [3 X 3]");
        }
        int n = fm[0].length;

        /*
        u1_normalized = T1 * u1
          u2_normalized = T2 * u2

          denormalized u1 = T1^-1 * u1_normalized
          denormalized u2 = T2^-1 * u2_normalized

          FM_normalized = inverse(transpose(T2)) * FM * inverse(T1)
            with caveat about centroid and normalization details...

         FM_normalized * T1 = inverse(transpose(T2)) * FM
         transpose(T2) * FM_normalized * T1 = FM

         denormalized FM = transpose(T2) * FM_normalized * T1

         */
        double[][] t2T = MatrixUtil.transpose(t2);
        double[][] fmD = MatrixUtil.multiply(t2T, fm);
        fmD = MatrixUtil.multiply(fmD, t1);

        for (int i = 0; i < fmD.length; ++i) {
            System.arraycopy(fmD[i], 0, fm[i], 0, fmD[i].length);
        }
    }

    static double[][] transposeInverseT(double[][] t) {
        if (t.length != 3 && t.length != 4) {
            throw new IllegalArgumentException("t.length must be 3 or 4");
        }
        int nd = t.length - 1;

        double[][] tInvT = new double[t.length][t.length];
        for (int i = 0; i < nd; ++i) {
            tInvT[i][i] = 1./t[i][i];
            //tInvT[i][t.length - 1] = -1 * t[i][t.length - 1]/t[i][i];
            // for transposed:
            tInvT[t.length - 1][i] = -1 * t[i][t.length - 1]/t[i][i];
        }
        tInvT[nd][nd] = 1;

        return tInvT;
    }

    public static double[][] inverseT(double[][] t) {

        double[][] inv = MatrixUtil.zeros(3,3);
        inv[0][0] = 1./t[0][0];
        inv[1][1] = 1./t[1][1];
        inv[2][2] = 1;
        inv[0][2] = -1.*t[0][2]/t[0][0];
        inv[1][2] = -1.*t[1][2]/t[1][1];

        /* from camera calibration inverse
           | 1/fx   0           -xc/fx  |
         = | 0       1/fy        -yc/fy  |
           | 0       0           1       |
        */

        return inv;
    }

}
