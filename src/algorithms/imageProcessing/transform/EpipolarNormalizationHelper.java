package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;

public class EpipolarNormalizationHelper {

    /*
    <pre>
        Revisiting Hartley’s Normalized Eight-Point Algorithm
        Chojnacki et al. 2003

        inverse changes the order of operations,
        inverse translation matrix: inverse changes the signs of the translation elements,
                but not the diagonal.
        inverse rotation matrix: inverse is the transpose of rotation matrix.
        inverse scaling matrix: inverse is performed on each element, that is, the reciprocal.

        denormalized x2 = transpose(normalized x2) * transpose(T2^1)
        denormalized x1 = (T1^-1) * (normalized x1)
        denormalized FM = transpose(T2) * FM_normalized * T1

                | 1/s   0  0 |   | 1  0  -xc |
            T = |  0  1/s  0 | * | 0  1  -yc |
                |  0    0  1 |   | 0  0   1  |

                   | 1  0  xc |   | s  0   0 | = | s  0  xc |
            T^-1 = | 0  1  yc | * | 0  s   0 |   | 0  s  yc |
                   | 0  0   1 |   | 0  0   1 |   | 0  0  1  |

            transposed(T^-1) = | s   0    0 |
                               | 0   s    0 |
                               | xc  yc   1 |

                          | 1  0  xc |   | s^2  0    0 |   | 1  0  0 |
            T^-1 * T^-T = | 0  1  yc | * | 0   s^2   0 | * | 0  1  0 |
                          | 0  0   1 |   | 0    0    1 |   | xc yc 1 |

                          | s^2 + xc^2   xc*yc       xc |
                        = | yc*xc        s^2 + yc^2  yc |
                          | xc           yc          1  |
         </pre>
         */

    /**
     * given data points x, perform unit standard normalization to modify
     * x to have mean = 0 and standard deviation = 1.
     * the transformation matrix is returned.
     * <pre>
     *     the transformation matrix
     *     double[][] t = new double[3][];
     *         t[0] = new double[]{1./scale,       0,     -centroidX/scale};
     *         t[1] = new double[]{0,           1./scale, -centroidY/scale};
     *         t[2] = new double[]{0,           0,           1};
     *
     *     to denormalize the points and the fundamental matrix:
     *
     *     denormalized x2 = transpose(normalized x2) * transpose(T2^1)
     *     denormalized x1 = (T1^-1) * (normalized x1)
     *     denormalized FM = transpose(T2) * FM_normalized * T1
     *
     *     where T1 is the transformation matrix from points x1
     *     and T2 is the transformation matrix from points x2.
     *
     * </pre>
     * @param x matrix of data in format [3 X nPoints] where row 0 are the
     *          'x-axis' coordinates, row 1 are the 'y-axis' coordinates,
     *          and row 3 are the 'z-axis' coordinates which are expected
     *          to be all '1' (= homogenous coordinates).
     *
     * @return transformation matrix used to normalize x.
     <pre>
    t[0] = new double[]{1./scale,       0,     -centroidX/scale};
    t[1] = new double[]{0,           1./scale, -centroidY/scale};
    t[2] = new double[]{0,           0,           1};
    </pre>
     */
    public static double[][] unitStandardNormalize(double[][] x) {
        if (x.length != 3) {
            throw new IllegalArgumentException("x.length must be 3");
        }
        int n = x[0].length;

        double[] mS = new double[4];

        int i;
        int j;
        for (i = 0; i < n; ++i) {
            for (j = 0; j < 2; ++j) {
                mS[j] += x[j][i];
            }
        }

        for (j = 0; j < 2; ++j) {
            mS[j] /= n;
        }

        double d;
        for (i = 0; i < n; ++i) {
            for (j = 0; j < 2; ++j) {
                d = (x[j][i] - mS[j]);
                mS[j + 2] += (d * d);
            }
        }
        for (j = 2; j < 4; ++j) {
            mS[j] = Math.sqrt(mS[j]/(n - 1.0));
        }

        double[][] t = new double[3][];
        t[0] = new double[]{1./mS[2],       0,     -mS[0]/mS[2]};
        t[1] = new double[]{0,           1./mS[3], -mS[1]/mS[3]};
        t[2] = new double[]{0,           0,           1};

        double[][] xyN = MatrixUtil.multiply(t, x);
        for (i = 0; i < 3; ++i) {
            System.arraycopy(xyN[i], 0, x[i], 0, xyN[i].length);
        }
        return t;
    }

    /**
     * denormalized x2 = transpose(normalized x2) * transpose(T2^1)
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
        denormalized x2 = transpose(normalized x2) * transpose(T2^1)
        denormalized x1 = (T1^-1) * (normalized x1)
        denormalized FM = transpose(T2) * FM_normalized * T1
         */
        double[][] tInv2 = transposeInverseT(t2);
        double[][] tInv1 = inverseT(t1);

        double[][] x2D = MatrixUtil.multiply(MatrixUtil.transpose(x2), tInv2);
        double[][] x1D = MatrixUtil.multiply(tInv1, x1);

        int i;
        for (i = 0; i < x1[0].length; ++i) {
            System.arraycopy(x1D[i], 0, x1, 0, x1D[i].length);
        }
        for (i = 0; i < x2[0].length; ++i) {
            System.arraycopy(x2D[i], 0, x2, 0, x2D[i].length);
        }
    }

    /**
     * denormalized FM = transpose(T2) * FM_normalized * T1
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

        //denormalized FM = transpose(T2) * FM_normalized * T1
        double[][] fmD = MatrixUtil.multiply(fm, t1);
        fmD = MatrixUtil.multiply(MatrixUtil.transpose(t2), fmD);

        for (int i = 0; i < fmD.length; ++i) {
            System.arraycopy(fmD[i], 0, fm[i], 0, fmD[i].length);
        }
    }

    private static double[][] transposeInverseT(double[][] t) {

        double[][] tinv = MatrixUtil.zeros(3, 3);
        tinv[0][0] = 1. / t[0][0];
        tinv[1][1] = 1. / t[1][1];
        tinv[2][2] = 1;
        tinv[2][0] = -1. * t[0][2];
        tinv[2][1] = -1. * t[1][2];

        return tinv;
    }

    private static double[][] inverseT(double[][] t) {

        double[][] inv = MatrixUtil.zeros(3,3);
        inv[0][0] = 1./t[0][0];
        inv[1][1] = 1./t[1][1];
        inv[2][2] = 1;
        inv[0][2] = -1.*t[0][2];
        inv[1][2] = -1.*t[1][2];

        return inv;
    }

}
