package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.ImageProcessor;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import algorithms.util.PolygonAndPointPlotter;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * various optical flow algorithms.
 *
 * TODO: consider implementing Baker-Matthews for KLT tracker
 */
public class OpticalFlow {

    private static final Level LEVEL = Level.FINE;
    private static final Logger log;
    static {
        log = Logger.getLogger(CameraCalibration.class.getSimpleName());
    }

    /**
     an algorithm that estimates the motion between 2 images, assuming constancy of brightness and short elapsed times.
     The algorithm can be used for small displacements between images.
     The algorithm is extremely sensitive to the size of the input images.
     <pre>
     to help determine the best image size to use for this algorithm, here are some empirically derived numbers from
     unit tests:
         image width 32, maxV=[32,*), alpha=0.1, maxIter=100, epsSq=1E-2, up to dx=14
         image width 64, maxV=[64,*), alpha=0.1, maxIter=100, epsSq=1E-2, up to dx=16
         image width 64, maxV=[64,*), alpha=0.1, maxIter=1000, epsSq=1E-4, up to dx=53
         image width 128, maxV=[128,*), alpha=0.1, maxIter=100, epsSq=1E-2, up to dx=9
         image width 128, maxV=[128,*), alpha=0.1, maxIter=100, epsSq=[1E-4, 1E-3], up to dx=30
         image width 128, maxV=[128,*), alpha=0.1, maxIter=1000, epsSq=1E-4, up to dx=101
         image width 128, maxV=[128,*), alpha=0.1, maxIter=[1000, 10000], epsSq=1E-3, up to dx=31
         image width 255, maxV=[128,*), alpha=0.1, maxIter=[1000, 200000] epsSq=[1E-10, 1E-3], up to dx=2

     summary of those results:
         ---------------------------------------------------
         motion    image dimension for fastest calculations
         -------  -----------------------------------------
         14 pix       32
         53 pix       64
         101 pix      128

     </pre>
     <pre>
     references for the algorithm details:

     cmu lectures for Computer Vision 16-385, lectures 14 by Kris Kitani

     Horn, B.K. and Schunck, B.G., 1980. Determining optical flow.

     boundary condition handling was adopted from:
     Meinhardt-Llopis, Pérez, and Kondermann, "Horn-Schunck Optical Flow with a Multi-Scale Strategy",
     Image Processing On Line, 3 (2013), pp. 151–172. https://doi.org/10.5201/ipol.2013.20
     </pre>
     uses default values of alpha=1, epsSq=1E-9, maxIter=1000, uInit=0, vInit=0
     * @param im1 image array.  for best results, consider the image to be integers, that is, signal range should be > 1
     *            i.e, [0, 255]
     * @param im2 image array taken shortly after im1 in time.
     * @return u and v as motion between im1 and im2 in units of pixels.  u is dx/dt and v is dy/dt.
     */
    public static List<double[][]> hornSchunck(double[][] im1, double[][] im2) {
        return hornSchunck(im1, im2, 0, 0, 0.1, 1000, 1E-4);
    }

    /**
     an algorithm that estimates the motion between 2 images, assuming constancy of brightness and short elapsed times.
     The algorithm can be used for small displacements between images that are mostly smooth gradients.
     The algorithm is extremely sensitive to the size of the input images.
     <pre>
     to help determine the best image size to use for this algorithm, here are some empirically derived numbers from
     unit tests:
         image width 32, maxV=[32,*), alpha=0.1, maxIter=100, epsSq=1E-2, up to dx=14
         image width 64, maxV=[64,*), alpha=0.1, maxIter=100, epsSq=1E-2, up to dx=16
         image width 64, maxV=[64,*), alpha=0.1, maxIter=1000, epsSq=1E-4, up to dx=53
         image width 128, maxV=[128,*), alpha=0.1, maxIter=100, epsSq=1E-2, up to dx=9
         image width 128, maxV=[128,*), alpha=0.1, maxIter=100, epsSq=[1E-4, 1E-3], up to dx=30
         image width 128, maxV=[128,*), alpha=0.1, maxIter=1000, epsSq=1E-4, up to dx=101
         image width 128, maxV=[128,*), alpha=0.1, maxIter=[1000, 10000], epsSq=1E-3, up to dx=31
         image width 255, maxV=[128,*), alpha=0.1, maxIter=[1000, 200000] epsSq=[1E-10, 1E-3], up to dx=2

     summary of those results:
         ---------------------------------------------------
         motion    image dimension for fastest calculations
         -------  -----------------------------------------
         14 pix       32
         53 pix       64
         101 pix      128

     </pre>
     <pre>
     references for the algorithm details:

     cmu lectures for Computer Vision 16-385, lectures 14 by Kris Kitani

     Horn, B.K. and Schunck, B.G., 1980. Determining optical flow.

     boundary condition handling was adopted from:
     Meinhardt-Llopis, Pérez, and Kondermann, "Horn-Schunck Optical Flow with a Multi-Scale Strategy",
     Image Processing On Line, 3 (2013), pp. 151–172. https://doi.org/10.5201/ipol.2013.20
     </pre>
     uses default values of alpha=1, epsSq=1E-9, maxIter=1000, uInit=0, vInit=0
     * @param im1 image array.  for best results, consider the image to be integers, that is, signal range should be > 1
     *            i.e, [0, 255] etc.
     * @param im2 image array taken shortly after im1 in time.
     * @param uInit initial guess for u.
     * @param vInit initial guess for v.
     * @param alpha a weight parameter used in the updates of u and v internally.  must be > 0.
     *              for best results should be in range [1,10].
     *              alpha=0.1 is recommended.
     * @param maxIter maximum number of iterations to perform.  if u,v square differences do not converge to
     *                <= epsSq, up to maxIter iterations are performed.
     * @param epsSq the condition for convergence.  sum of the square differences of u,v from previous iteration
     *              less than epsSq is convergence.
     * @return
     */
    public static List<double[][]> hornSchunck(double[][] im1, double[][] im2,
                                               final double uInit, final double vInit,
                                               final double alpha, final int maxIter, final double epsSq) {

        if (alpha <= 0) {
            throw new IllegalArgumentException("alpha must be > 0");
        }
        int h = im1.length;
        int w = im1[0].length;
        if (im2.length != h || im2[0].length != w) {
            throw new IllegalArgumentException("im2 and im2 should have same dimensions");
        }

        // image gradients
        double[][] gX = hsGradX(im1, im2);
        double[][] gY = hsGradY(im1, im2);
        // temporal differences
        double[][] gT = hsGradT(im1, im2);

        //double denom = (1./lambda) + dIx*dIx + dIy*dIy;
        //double denom = alpha*alpha + gX*gX + gY*gY;
        double[][] denom = new double[h][w];
        for (int r = 0; r < h; ++r) {
            for (int c = 0; c < w; ++c) {
                denom[r][c] = (alpha*alpha + gX[r][c]*gX[r][c] + gY[r][c]*gY[r][c]);
            }
        }

        double[][] u = new double[h][w];
        double[][] v = new double[h][w];
        if (uInit != 0) {
            for (double[] uI : u) {
                Arrays.fill(uI, uInit);
            }
        }
        if (vInit != 0) {
            for (double[] vI : v) {
                Arrays.fill(vI, vInit);
            }
        }

        double[][] uAvg = new double[h][w];
        double[][] vAvg = new double[h][w];

        int nIter = 0;

        // smaller values of lambda for smoother flow
        //double lambda = 0.5;//1./(3.*3.);//1./(15.*15.);
        // alpha squared models an addition of moise
        //Math.pow(MiscMath0.mean(MatrixUtil.rowMeans(im1)), 2);
        //System.out.println("gX mean^2=" + Math.pow(MiscMath0.mean(MatrixUtil.rowMeans(gX)), 2));
        //System.out.println("gT mean^2=" + Math.pow(MiscMath0.mean(MatrixUtil.rowMeans(gT)), 2));
        //System.out.println("im1 max^2=" + Math.pow(MiscMath0.findMax(im1), 2));

        boolean hasConverged = false;

        ImageProcessor imageProcessor = new ImageProcessor();

        double[][] prevU = null;
        double[][] prevV = null;

        double sumDiff = 0;
        while (nIter == 0 || nIter < maxIter && !hasConverged) {
            sumDiff = 0;
            ++nIter;

            hsConvolve(uAvg, u);
            hsConvolve(vAvg, v);

            for (int r = 0; r < h; ++r) {
                for (int c = 0; c < w; ++c) {
                    double numer = (gX[r][c] * uAvg[r][c]) + (gY[r][c] * vAvg[r][c]) + gT[r][c];
                    double nd = numer/denom[r][c];
                    u[r][c] = uAvg[r][c] - gX[r][c] * nd;
                    v[r][c] = vAvg[r][c] - gY[r][c] * nd;
                    if (prevU != null) {
                        sumDiff += (Math.pow(u[r][c] - prevU[r][c], 2) + Math.pow(v[r][c] - prevV[r][c], 2));
                    }
                }
            }
            if (prevU != null) {
                sumDiff /= (h * w);
                hasConverged = (sumDiff < epsSq);
            }
            prevU = MatrixUtil.copy(u);
            prevV = MatrixUtil.copy(v);
        }

        log.log(LEVEL, "nIter=" + nIter + ", sumDiff=" + sumDiff);

        List<double[][]> out = new ArrayList<>();
        out.add(u);
        out.add(v);
        return out;

        /*
        from cmu lectures for Computer Vision 16-385, lectures 14 by Kris Kitani
        1. Precompute image gradients
             StructureTensorD for each image
        2. Precompute temporal gradients
             frame differencing
        3. Initialize flow field
             u=dx/dy, v=dy/dt
             init u=0, v=0
        4. While not converged:
        Compute flow field updates for each pixel:
           u_est_{k,l} = u_{k,l}
                          - (I_x * u_{k,l} + I_y * v_{k,l} + I_t) * I_x
                          / ((1/lambda) + (I_x)^2) + (I_y)^2))
           v_est_{k,l} = v_{k,l}
                         - (I_x * u_{k,l} + I_y * v_{k,l} + I_t) * I_y
                          / ((1/lambda) + (I_x)^2) + (I_y)^2))
         */
    }

    /**
     * calculate the flow of im1 over patches of dimension pathDimension where the patch centers are given as x,y coordinates.
     * The flow u,v is returned for each patch.
     * Corners make good patch centers.
     <pre>
     reference for the algorithm details:
         cmu lectures for Computer Vision 16-385, lectures 14 by Kris Kitani
     </pre>
     <pre>
     Some unit test results on a diagonal gradient image, printed here as information about estimatable motion range
     in an image size under ideal conditions:

     using the median of the calculated u's, v's as the estimated u (=dx/dt) over all patches:
     m = image length = image width.
     max value is maximum pixel intensity in test image.

             m     patchDim   max value    allowing error m/10   nPatches  can calc up to max u
             32      5         32             3                     21         23
             64      5         64             6                     53         57
             128     5         128            12                    117        120
             256     5         256            25                    245        250   (largest error=7)
             512     5         512            51                    501        506   (largest error=7)
             1024    5         1024           102                   1013       1018  (largest error=7)
     </pre>
     * TODO: consider overloading method for features specified in some format
     * @param im1 image array.  for best results, should hold values similar to an image range [0, 255], etc.
     * @param im2 image array taken shortly after im1 in time.
     * @param patchCenters x,y coordinates of patch centers in format of size [nPatches x 2], e.g.
     *                     patchCenters[0][0] = x0, patchCenters[0][1] = y0,
     *                     patchCenters[1][0] = x1, patchCenters[1][1] = y1, ...
     *        Note, patchCenters +- pathDimension/2 should be within bounds of image.
     * @param patchDimension the width of a square patch to determine flow within. e.g. if patchDimension is 5,
     *                       the flow is determined in 5x5 patches.  should be >= 3.
     * @return list of uv pairs, one for each patchCenter.
     */
    public static List<double[]> lucasKanade(double[][] im1, double[][] im2, int[][] patchCenters, int patchDimension) throws NotConvergedException, IOException {
        int h = im1.length;
        int w = im1[0].length;
        if (im2.length != h || im2[0].length != w) {
            throw new IllegalArgumentException("im2 and im2 should have same dimensions");
        }

        if (patchDimension < 3) {
            throw new IllegalArgumentException("patchDimension should be >= 3");
        }

        List<double[]> uvList = new ArrayList<>();

        double thresh = 0.04;

        // image gradients
        double[][] gX = hsGradX(im1, im2);
        double[][] gY = hsGradY(im1, im2);
        // temporal differences
        double[][] gT = hsGradT(im1, im2);

        double[] u = new double[patchCenters.length];
        double[] v = new double[patchCenters.length];
        int uvI = 0;

        for (int[] pc : patchCenters) {
            // build A [patchDimension*patchDimension x 2]

            double[][] A = new double[2][2];
            // build b [patchDimension]
            double[] b = new double[2];

            for (int dX = -patchDimension/2; dX <= patchDimension/2; ++dX) {
                if ((patchDimension & 1) != 1 && dX == patchDimension/2) continue;
                int c = pc[0] + dX;
                if (c < 0 || c >= im2[0].length)
                    throw new IllegalArgumentException("col " + c + " out of bounds of width: " + c);
                for (int dY = -patchDimension/2; dY <= patchDimension/2; ++dY) {
                    if ((patchDimension & 1) != 1 && dY == patchDimension/2) continue;
                    int r = pc[1] + dY;
                    if (r < 0 || r >= im2.length)
                        throw new IllegalArgumentException("row " + r + " out of bounds of height: " + r);

                    A[0][0] += (gX[r][c] * gX[r][c]);
                    A[0][1] += (gX[r][c] * gY[r][c]);
                    A[1][1] += (gY[r][c] * gY[r][c]);
                    b[0] += gX[r][c] * gT[r][c];
                    b[1] += gY[r][c] * gT[r][c];
                } // end for loop dY
            }// end for loop dX
            A[1][0] = A[0][1];
            //xEst = pseudoInv(A) * b // [2x2] * [2x1] = [2X1]

            //A^T*A should be invertible
            //A^T*A should not be too small (length >= 2)
            // eigenvalue1/eigenvalue2 should not be too large

            // can use Harris corners response > 0 to ensure good patch to include
            double detA = MatrixUtil.determinant(A);
            double trA = MatrixUtil.trace(A);
            boolean cont = (detA - thresh * trA * trA) > 0;

            /* or that there is power in eigen[1], hence rank =2
            double[][] ATA = MatrixUtil.createATransposedTimesA(A);
            EVD evd = no.uib.cipr.matrix.EVD.factorize(new DenseMatrix(ATA));
            double[] eigenValues = evd.getRealEigenvalues();
            boolean cont2 = eigenValues[eigenValues.length - 1] >= thresh;
            */
            if (!cont)
                continue;

            double[] uv = MatrixUtil.multiplyMatrixByColumnVector(MatrixUtil.pseudoinverseFullRowRank(A), b);
            uvList.add(uv);

            u[uvI] = uv[0];
            v[uvI] = uv[1];
            ++uvI;
        }

        return uvList;
    }

    protected static double[][] strictColumnDiff(double[][] a) {
        double[][] b = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; ++i) {
            for (int j = 1; j < b.length; ++j) {
                b[i][j] = a[i][j] - a[i][j-1];
            }
        }
        return b;
    }

    /**
     * extract value from x w/ boundary replication when out of bounds.
     * method adopted from
     <pre>
     boundary condition handling was adopted from:
     Meinhardt-Llopis, Pérez, and Kondermann, "Horn-Schunck Optical Flow with a Multi-Scale Strategy",
     Image Processing On Line, 3 (2013), pp. 151–172. https://doi.org/10.5201/ipol.2013.20
     </pre>
     * @param x
     * @param col
     * @param row
     * @return
     */
    protected static double p(double[][] x, int row, int col) {
        int w = x[0].length;
        int h = x.length;
        if (col < 0) col = 0;
        if (row < 0) row = 0;
        if (col >= w) col = w - 1;
        if (row >= h) row = h - 1;
        return x[row][col];
    }

    /**
     create y gradient
     <pre>
     Horn & Schmunck 1980, Sect 7
     </pre>
     * @param im1
     * @param im2
     * @return
     */
    protected static double[][] hsGradY(double[][] im1, double[][] im2) {
        int h = im1.length;
        int w = im1[0].length;
        double[][] b = new double[h][w];

        for (int j = 0; j < h; ++j) {
            for (int i = 0; i < w; ++i) {
                b[j][i] = p(im1, j+1,i) - p(im1, j, i)
                        + p(im1,j+1,i+1) - p(im1,j,i+1)
                        + p(im2, j+1,i) - p(im2, j, i)
                        + p(im2,j+1,i+1) - p(im2,j,i+1);
                b[j][i] /= 4.;
            }
        }
        return b;
    }

    protected static double[][] strictRowDiff(double[][] a) {
        double[][] b = new double[a.length][a[0].length];
        for (int i =1; i < a.length; ++i) {
            for (int j = 0; j < b.length; ++j) {
                b[i][j] = a[i][j] - a[i-1][j];
            }
        }
        return b;
    }

    //Horn & Schmunck 1980, Sect 7
    protected static double[][] hsGradX(double[][] im1, double[][] im2) {
        int h = im1.length;
        int w = im1[0].length;
        double[][] b = new double[h][w];

        for (int j = 0; j < h; ++j) {
            for (int i = 0; i < w; ++i) {
                b[j][i] = p(im1, j,i+1) - p(im1, j, i)
                        + p(im1, j+1, i+1) - p(im1, j+1, i)
                        + p(im2, i+1,j) - p(im2, j, i)
                        + p(im2, j+1, i+1) - p(im2, j+1, i);
                b[j][i] /= 4.;
            }
        }
        return b;
    }

    //Horn & Schmunck 1980, Sect 7
    protected static double[][] hsGradT(double[][] im1, double[][] im2) {
        int h = im1.length;
        int w = im1[0].length;
        double[][] b = new double[h][w];

        for (int j = 0; j < h; ++j) {
            for (int i = 0; i < w; ++i) {
                b[j][i] = p(im2,j,i) - p(im1,j,i)
                        + p(im2,j,i+1) - p(im1,j,i+1)
                        + p(im2,j+1,i) - p(im1,j+1,i)
                        + p(im2,j+1,i+1) - p(im1,j+1,i+1);
                b[j][i] /= 4.;
            }
        }

        return b;
    }

    /**
     * weighted sum over the 8 neighbors.
     * Horn, Schunck 1980  section 8
     * @param out
     * @param in
     */
    private static void hsConvolve(double[][] out, double[][] in) {
        int h = in.length;
        int w = in[0].length;
        for (int j = 0; j < h; ++j) {
            for (int i = 0; i < w; ++i) {
                out[j][i] = (p(in, j-1, i) + p(in,j+1,i)
                        + p(in,j,i-1) + p(in,j,i+1))/6.;
                out[j][i] += (p(in, j-1, i-1) + p(in,j-1,i+1)
                        + p(in,j+1,i-1) + p(in,j+1,i+1))/12.;
            }
        }
    }

}
