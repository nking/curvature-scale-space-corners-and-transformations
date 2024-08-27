package algorithms.imageProcessing.features;

import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.StructureTensorD;
import algorithms.matrix.MatrixUtil;
import algorithms.sort.MiscSorter;

import java.util.Arrays;

public class ShiTomasi {

    /**
     * returns the 2nd strongest eigenvalues of structure tensor matrix:
     <pre>
       Z = [gX^2  gX*gY]
           [gXgY  gY^2 ]

       where gX and gY are the gradients.
     </pre>
     Note that row-major notation is used.
     * @param a
     * @param  outEigs is indexed by (row*a.length) + col and holds values that are the 2nd strongest eigenvalue for the
     *                 structure tensor at each point a[row][col].
     */
    static void tensorSecondEigenvalues(double[][] a, double[] outEigs) {

        StructureTensorD t = new StructureTensorD(a, 1.0f, false);
        //double[][] gTX = t.getDX();
        //double[][] gTY = t.getDY();

        ImageProcessor imageProcessor = new ImageProcessor();
        double[][] gTX = MatrixUtil.copy(a);
        double[][] gTY = MatrixUtil.copy(a);
        // {0, -1, 1} should match the paper diff in method sumSteepestDescErrImageProduct if change back to it
        imageProcessor.applyKernel1D(gTX, new double[]{0, 1, -1}, false);
        imageProcessor.applyKernel1D(gTY, new double[]{0, 1, -1}, true);

        double[][] z = new double[2][2];
        double gX, gY;
        double[] eigI = null;
        for (int row = 0, idx=0; row < a.length; ++row) {
            for (int col = 0; col < a[row].length; ++col, ++idx) {
                gX = gTX[row][col];
                gY = gTY[row][col];

                z[0][0] = gX*gX;
                z[0][1] = gX*gY;
                z[1][0] = z[0][1];
                z[1][1] = gY*gY;

                eigI = MatrixUtil.eigenvalues2X2(z);
                outEigs[idx] = eigI[1];
            }
        }
    }

    /**
     * returns the nPoints of 2nd strongest eigenvalues of the structure tensor of each point.
     * @param a
     * @param nPoints number of points to return
     * @param outMinMax is populated with the 2nd strongest eigenvalues for the first and the last points
     *                  for the coordinate list returned.
     * @return the strongest 2nd strongest eigenvalues of the structure tensor of each point if not == 0.
     */
    public static int[][] goodFeatureCoordinates(double[][] a, int nPoints, double[] outMinMax) {

        if (outMinMax.length != 2) {
            throw  new IllegalArgumentException("outMinMax length must be 2");
        }

        int w = a[0].length;
        int h = a.length;
        int wh = w*h;
        double[] vals = new double[wh];
        tensorSecondEigenvalues(a, vals);

        int[] idxs =  MiscSorter.mergeSortDecreasing(vals);

        int[][] xy = new int[nPoints][2];
        int iP = 0;
        int idx, y, x;
        for (int i = 0; i < nPoints; ++i) {
            if (vals[i] == 0.) continue;
            idx = idxs[i];
            y = idx / w;
            x = idx % w;
            xy[iP][0] = x;
            xy[iP][1] = y;
            ++iP;
        }
        xy = Arrays.copyOf(xy, iP);
        if (iP == 0) {
            Arrays.fill(outMinMax, 0);
        } else {
            outMinMax[0] = vals[0];
            if (iP > 1) {
                outMinMax[1] = vals[iP - 1];
            }
        }

        return xy;
    }
}
