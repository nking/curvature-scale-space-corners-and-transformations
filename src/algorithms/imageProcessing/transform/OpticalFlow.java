package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.StructureTensor;
import algorithms.imageProcessing.StructureTensorD;
import algorithms.matrix.MatrixUtil;
import no.uib.cipr.matrix.NotConvergedException;

import java.util.ArrayList;
import java.util.List;

public class OpticalFlow {

    /**
     * @param im1
     * @param im2
     * @return u and v for im1
     */
    public static List<double[][]> hornSchunck(double[][] im1, double[][] im2) {

        int m = im1.length;
        int n = im1[0].length;
        if (im2.length != m || im2[0].length != n) {
            throw new IllegalArgumentException("im2 and im2 should have same dimensions");
        }

        // image dradients:
        StructureTensorD s1 = new StructureTensorD(im1, 1.0f, false);
        //StructureTensorD s2 = new StructureTensorD(im1, 1.0f, false);

        // temporal differences (gradients)
        double[][] t12 = MatrixUtil.pointwiseSubtract(im1, im2);

        double[][] u = new double[m][n];
        double[][] v = new double[m][n];

        double[][] uEst;
        double[][] vEst;

        double epsSq = 1E-8;
        int maxIter = 200;//2_000;
        int nIter = 0;

        // smaller values of lambda for smoother flow
        double lambda = 1./(15.*15.);//4;// 1/L = a*a

        boolean hasConverged = false;

        while (nIter == 0 || nIter < maxIter && !hasConverged) {
            double sumDiff = 0;
            ++nIter;
            uEst = new double[m][n];
            vEst = new double[m][n];
            for (int r = 0; r < m; ++r) {
                for (int c = 0; c < n; ++c) {
                    double numer = (s1.getDX()[r][c] * u[r][c]) + (s1.getDY()[r][c] * v[r][c]) + t12[r][c];
                    double denom = (1/lambda) + s1.getDX()[r][c]*s1.getDX()[r][c] + s1.getDY()[r][c]*s1.getDY()[r][c];
                    double nd = numer/denom;
                    uEst[r][c] = u[r][c] - s1.getDX()[r][c] * nd;
                    vEst[r][c] = v[r][c] - s1.getDY()[r][c] * nd;
                    sumDiff += (Math.pow(uEst[r][c] - u[r][c], 2) + Math.pow(vEst[r][c] - v[r][c], 2));
                }
            }
            sumDiff /= (m*n);
            hasConverged = (sumDiff < epsSq);
            u = uEst;
            v = vEst;
        }

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
     * TODO: consider overloading method for features specified in some format
     * @param im1
     * @param im2
     * @param patchCenters x,y coordinates of patch centers in format of size [nPatches x 2], e.g.
     *                     patchCenters[0][0] = x0, patchCenters[0][1] = y0,
     *                     patchCenters[1][0] = x1, patchCenters[1][1] = y1, ...
     *        Note, patchCenters +- pathDimension/2 should be within bounds of image.
     * @param patchDimension the width of a square patch to determine flow within. e.g. if patchDimension is 5,
     *                       the flow is determined in 5x5 patches.
     * @return list of uv pairs, one for each patchCenter.
     */
    public static List<double[]> lucasKanade(double[][] im1, double[][] im2, int[][] patchCenters, int patchDimension) throws NotConvergedException {
        int m = im1.length;
        int n = im1[0].length;
        if (im2.length != m || im2[0].length != n) {
            throw new IllegalArgumentException("im2 and im2 should have same dimensions");
        }

        List<double[]> uvList = new ArrayList<>();

        // image gradients:
        StructureTensorD s1 = new StructureTensorD(im1, 1.0f, false);
        //StructureTensorD s2 = new StructureTensorD(im1, 1.0f, false);

        // temporal differences (gradients)
        double[][] t12 = MatrixUtil.pointwiseSubtract(im1, im2);
        for (int[] pc : patchCenters) {
            // build A [patchDimension*patchDimension x 2]
            double[][] A = new double[2][2];
            // build b [patchDimension]
            double[] b = new double[2];
            for (int dX = -patchDimension/2; dX <= patchDimension/2; ++dX) {
                if ((patchDimension & 1) != 1 && dX == patchDimension/2) continue;
                int c = pc[0] + dX;
                if (c < 0 || c >= im2[0].length) throw new IllegalArgumentException("out of bounds of width: " + c);
                for (int dY = -patchDimension/2; dY <= patchDimension/2; ++dY) {
                    if ((patchDimension & 1) != 1 && dY == patchDimension/2) continue;
                    int r = pc[1] + dY;
                    if (r < 0 || r >= im2.length) throw new IllegalArgumentException("out of bounds of height: " + r);

                    A[0][0] += (s1.getDX()[r][c] * s1.getDX()[r][c]);
                    A[0][1] += (s1.getDX()[r][c] * s1.getDY()[r][c]);
                    A[1][1] += (s1.getDY()[r][c] * s1.getDY()[r][c]);
                    A[1][0] = A[0][1];

                    b[0] += s1.getDX()[r][c] * t12[r][c];
                    b[1] += s1.getDY()[r][c] * t12[r][c];
                } // end for loop dY
            }// end for loop dX
            //xEst = pseudoInv(A) * b // [2x2] * [2x1] = [2X1]
            double[] uv = MatrixUtil.multiplyMatrixByColumnVector(MatrixUtil.pseudoinverseFullRowRank(A), b);
            uvList.add(uv);
        }
        return uvList;
    }
}
