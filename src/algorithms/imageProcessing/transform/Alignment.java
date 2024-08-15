package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.StructureTensor;
import algorithms.imageProcessing.StructureTensorD;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;

import java.io.IOException;
import java.util.Arrays;

/**
 * various 2D and 3D alignment methods
 * 
 * @author nichole
 */
public class Alignment {

    /*
    NOT READY FOR USE... 
    refine pInit so that template and image are aligned by
    computing the warp that can be applied to the template to align it with the image.
    This alignment minimizes the re-projection error between the projected image and template
    (though the error is an LP1 sum rather than LP2).
    to get a good initial pInit, one could extract features from template and image and make correspondence list
    of matching features using tools such as ransac or assumptions of translation only etc depending upon context.
    <pre>
    reference:
    Figure 4 of
    Baker & Matthews, 2016, CMU-RI-TR-02-16
    "Lucas-Kanade 20 Years On: A Unifying Framework: Part 1"
    </pre>
     @param template
     * @param image
     * @param pInit the projection matrix for a homogenous coordinate.
     pInit should be size [2 rows X 3 columns] and a good initial estimate.
     pInit is multiplied by [x, y, 1]^T to apply warp to point [x, y, 1].
     @param maxIter the maximum number of iterations to perform
     * @throws NotConvergedException
     * @throws IOException
    */
    public static Result inverseCompositional2DAffine(double[][] template, double[][] image, double[][] pInit,
                                                      int maxIter)
            throws NotConvergedException, IOException {

        if (pInit.length != 2 || pInit[0].length != 3) {
            throw new IllegalArgumentException("expecting pInit length to be 2, and pInit[0].length to be 3");
        }

        /*
          since this isn't convex, nor even quasi-convex, we need a good initial guess to get the
          translation offsets correct (p5 and p6) presumably before the remaining parameters.

          TODO: consider improving the error image construction to make it quasi-convex
              from Ke & Kanade, "Quasi-convex Optimization for Robust Geometric Reconstruction"
                  convex function of y:  any norm function g(y) = ||y||_l.
                  affine function of X: h(X) = (p_u(X), p_v(X)).
                  The composition of convex function g and affine function h is a convex function: g ◦ h
                  ==> p(X) = (g ◦ h)(X) is a convex function and p(X) >= 0
         */

        int nTR = template.length;
        int nTC = template[0].length;
        int nT = nTR*nTC;

        //======  precompute  =======

        // gradients of T:   each is [nTR X nTC]
        StructureTensorD gT = new StructureTensorD(template, 1, false);

        /*
        W(x; p) for a 2D affine warp:
                = [1+p[0]      p2  p4 ] * [x]
                  [  p[1]  1+p[3]  p5 ]   [y]
                                          [1]
        pdW/pdp = [x, 0, y, 0, 1, 0]
                  [0, x, 0, y, 0, 1]
        */

        // ==== create the steepest decent image of the template ====
        //[nTR*nTC X 6]
        // steepest descent images  gT * dWdP = [nTR X nTC] * [nTR*nTC X 6] = [nTR*nTC X 6] for x then for y
        // can format as [nTr X nTc X 6] or [nTr * nTc X 6].
        // will choose the latter to make the Hessian mult easier
        //steepest descent img = gradient * dWdP
        //                       [1x2] * [2x6]
        double[][] dTdWdp = new double[nT][];
        double[][] tmp1 = new double[1][2];
        double[][] tmp2 = new double[2][6];
        for (int y = 0, idx = 0; y < nTR; ++y) {
            for (int x = 0; x < nTC; ++x, ++idx) {

                tmp1[0][0] = gT.getDX()[y][x];
                tmp1[0][1] = gT.getDY()[y][x];

                tmp2[0][0] = x;
                tmp2[0][2] = y;
                tmp2[0][4] = 1;
                tmp2[1][1] = x;
                tmp2[1][3] = y;
                tmp2[1][5] = 1;

                // idx=x*width + y;  x=idx%width; y=idx/width
                //[1X6]
                dTdWdp[idx] = MatrixUtil.multiply(tmp1, tmp2)[0];
            }
        }

        // Hessian: [6 x nTR*nTC]*[nTR*nTC X 6] = [6x6]
        double[][] hessianTmplt = MatrixUtil.createATransposedTimesA(dTdWdp);
        assert(hessianTmplt.length == 6);
        assert(hessianTmplt[0].length == 6);

        double[][] invHessianTmplt = MatrixUtil.pseudoinverseFullRowRank(hessianTmplt);

        //======  end precompute  =======

        double[][] wX = new double[nT][3];
        double[][] wY = new double[nT][3];

        double eps = 1E-3;

        double[] p = new double[6];
        // the Baker-Matthews algorithm uses column major notation, so will set p from pInit here,
        // then set results in pInit at end of algorithm
        p[0] = pInit[0][0];
        p[1] = pInit[1][0];
        p[2] = pInit[0][1];
        p[3] = pInit[1][1];
        p[4] = pInit[0][2];
        p[5] = pInit[1][2];

        double[] deltaP = new double[6];
        double[] pSum = new double[6];

        int nY = image.length;
        int nX = image[0].length;
        boolean converged = false;
        int nIter = 0;
        double[] nPAndErr = null;
        while (nIter < maxIter && !converged) {
            // warp image I with W(x;p) to compute I(W(x;p))
            if (nIter == 0) {
                computeWarpX(wX, p, nX, nY);
                computeWarpY(wY, p, nX, nY);
            } else {
                // update the warp
                updateWarpX(wX, deltaP, nX, nY);
                updateWarpY(wY, deltaP, nX, nY);
            }

            // (1),(2),(7): warp I and subtract from T, then mult by steepest descent image, summing over all x
            nPAndErr = sumSteepestDescErrImageProduct(image, template, wX, wY, dTdWdp, pSum);

            // compute deltaP as invHessianTmplt * pSum = [6 X 6] * [6 X 1]
            MatrixUtil.multiplyMatrixByColumnVector(invHessianTmplt, pSum, deltaP);

            // norm of deltaP < eps
            converged = MatrixUtil.lPSum(deltaP, 2) < eps;

            ++nIter;
        }

        System.out.printf("nIter=%d, converged=%b\npInit=%s\n",
                nIter, converged, FormatArray.toString(pInit, "%.3f"));

        double[] pEst = extractP(wX, wY, nX, nY);

        double[][] w = condenseWarp(wX, wY, nX, nY);

        pInit[0][0] = pEst[0];
        pInit[1][0] = pEst[1];
        pInit[0][1] = pEst[2];
        pInit[1][1] = pEst[3];
        pInit[0][2] = pEst[4];
        pInit[1][2] = pEst[5];

        System.out.printf("pFinal=%s\n", FormatArray.toString(pInit, "%.3f"));

        return new Result(w, nPAndErr[1]);
    }

    private static double[] extractP(double[][] wX, double[][] wY, int nX, int nY) throws NotConvergedException {

        /* from the warp, we can extract p
         [(1+p1)   p3    p5] * [x]
         [  p2   (1+p4)  p6]   [y]
                               [1]

         nT = nX * nY

        use the discrete information in wX and wY:

             b = A * c

             if we were to solve each parameter in p, one by one, could make 6 vectors and solve
                 A^T * b = A^T*A * c
                 (1/ (A^T*A)) * A^T * b = c

                 b = [nT X 1], A = [nT X 1], c = [1X1]

                 the equations for each of the 6:
                     Wx[0] = x * (1 + p1)
                     Wx[1] = y * p3
                     Wx[2] = p5
                     Wy[0] = x * p2
                     Wy[1] = y * (1+p4)
                     Wy[2] = p6

             or we can solve all 6 in matrix format:

             [wx[0]  wx[1]  wx[2]  wy[0]  wy[1]  wy[2]] = [x  y  1  x  y  1] * diagonal matrix of [(1+p1),p3,p5,p2,(1+p4),p6)]

             [nT X 6]    =  [nT X 6] * [6 X 6]

             C = diagonal matrix
               = A^-1 * B

               = S^-1 A2 S
                 where A2 = A^-1 * B
                 where S is columns of eigenvectors of A2
         */
        int nT = nX * nY;

        /*
        double[] p = new double[6];

        double[][] a = new double[nT][6];
        double[][] b = new double[nT][6];

        for (int y = 0, idx=0; y < nY; ++y) {
            for (int x = 0; x < nX; ++x, ++idx) {
                a[idx] = new double[]{wX[idx][0], wX[idx][1], wX[idx][2], wY[idx][0], wY[idx][1], wY[idx][2]};
                b[idx] = new double[]{x, y, 1, x, y, 1};
            }
        }

        //A rank is 6
        double[][] aInv = MatrixUtil.pseudoinverseFullRowRank(a);
        double[][] aInvB = MatrixUtil.multiply(aInv, b);

        // S columns are right eigenvectors where evd.rightEigenVectors are in the columns already
        // S^-1 rows are left eigenvectors where evd.leftEigenVectors are in the columns so need to be transposed
        EVD evd = EVD.factorize(new DenseMatrix(aInvB));
        if (evd.getRealEigenvalues().length < 6) {
            // can possibly still solve parameters one by one
            throw new RuntimeException("Error in extracting p from the warp.");
        }
        double[][] s = MatrixUtil.convertToRowMajor(evd.getRightEigenvectors());
        double[][] sInv = MatrixUtil.transpose(MatrixUtil.convertToRowMajor(evd.getLeftEigenvectors()));

        //S^-1 * (A^-1*b) * S
        double[][] c = MatrixUtil.multiply(sInv, MatrixUtil.multiply(aInvB, s));

        for (int i = 0; i < 6; ++i) {
            p[i] = c[i][i];
        }
        */

        // the single parameter at a time comparison
        double[] pEst = new double[6];
        /*
        A^T * b = A^T*A * c
                 (1/ (A^T*A)) * A^T * b = c

                 b = [nT X 1], A = [nT X 1], c = [1X1]

                 the equations for each of the 6:
                     Wx[0] = x * (1 + p1)
                     Wx[1] = y * p3
                     Wx[2] = p5
                     Wy[0] = x * p2
                     Wy[1] = y * (1+p4)
                     Wy[2] = p6
         */
        double[] aa = new double[nT];
        double[] bb = new double[nT];
        for (int y = 0, idx=0; y < nY; ++y) {
            for (int x = 0; x < nX; ++x, ++idx) {
                aa[idx] = wX[idx][0];
                bb[idx] = x;
            }
        }
        pEst[0] = (1./MatrixUtil.innerProduct(aa, aa)) * MatrixUtil.innerProduct(aa, bb);
        pEst[0] -= 1;

        for (int y = 0, idx=0; y < nY; ++y) {
            for (int x = 0; x < nX; ++x, ++idx) {
                aa[idx] = wX[idx][1];
                bb[idx] = y;
            }
        }
        pEst[2] = (1./MatrixUtil.innerProduct(aa, aa)) * MatrixUtil.innerProduct(aa, bb);

        for (int y = 0, idx=0; y < nY; ++y) {
            for (int x = 0; x < nX; ++x, ++idx) {
                aa[idx] = wX[idx][2];
                bb[idx] = 1;
            }
        }
        pEst[4] = (1./MatrixUtil.innerProduct(aa, aa)) * MatrixUtil.innerProduct(aa, bb);

        for (int y = 0, idx=0; y < nY; ++y) {
            for (int x = 0; x < nX; ++x, ++idx) {
                aa[idx] = wY[idx][0];
                bb[idx] = x;
            }
        }
        pEst[1] = (1./MatrixUtil.innerProduct(aa, aa)) * MatrixUtil.innerProduct(aa, bb);

        for (int y = 0, idx=0; y < nY; ++y) {
            for (int x = 0; x < nX; ++x, ++idx) {
                aa[idx] = wY[idx][1];
                bb[idx] = y;
            }
        }
        pEst[3] = (1./MatrixUtil.innerProduct(aa, aa)) * MatrixUtil.innerProduct(aa, bb);
        pEst[3] -= 1;

        for (int y = 0, idx=0; y < nY; ++y) {
            for (int x = 0; x < nX; ++x, ++idx) {
                aa[idx] = wY[idx][2];
                bb[idx] = 1;
            }
        }
        pEst[5] = (1./MatrixUtil.innerProduct(aa, aa)) * MatrixUtil.innerProduct(aa, bb);

        return pEst;
    }

    /**
     * add the 3 components of warp together for x and y respectively.  They were kept separate for use in composition
     * in other methods.
     * @param wX
     * @param wY
     * @param nX
     * @param nY
     * @return the warped x and y coordinates in an array of size [nX * nY X 2]
     */
    private static double[][] condenseWarp(double[][] wX, double[][] wY, int nX, int nY) {

        int nT = nX * nY;
        double[][] w = new double[nT][2];
        double xSum ;
        double ySum ;
        for (int y = 0, idx = 0; y < nY; ++y) {
            for (int x = 0; x < nX; ++x, ++idx) {
                xSum = 0;
                for (int k = 0; k < 3; ++k) {
                    xSum += wX[idx][k];
                }
                ySum = 0;
                for (int k = 0; k < 3; ++k) {
                    ySum += wY[idx][k];
                }
                w[idx] = new double[] {xSum, ySum};
            }
        }
        return w;
    }

    /**
     * calculates the error image, multiplies it by the steepest descent image transposed.  the results are summed
     * over all pixels and stored in outPSum.
     * The return values are the number of points used in the calculation (might be less than all image points due to
     * warped coordinates out of bounds), and the square root of the sum of squared error image values.
     *
     * The method handles
     * steps (1),(2),(7) of Figure 4 of the "Lucas-Kanade 20 Years on A Unifying..." for inverse compositional alignment
     * algorithm.
     * @param image
     * @param t
     * @param wX
     * @param wY
     * @param steep
     * @param outPSum the output sum of p over all points
     * @return double array holding the number of points used in the calculation and the quare root of the sum of squared error image values.
     * as double[]{nP, errSSD}
     */
    private static double[] sumSteepestDescErrImageProduct(double[][] image, double[][] t,
       double[][] wX, double[][] wY, double[][] steep, double[] outPSum) {

        // since the error image is only used in the calculation of the product of the steepest descent image
        // and the error image, we will return the result instead of intermediate products

        Arrays.fill(outPSum, 0);

        ImageProcessor imageProcessor = new ImageProcessor();

        // steepX and steepY are filled by
        // x=idx%width, y=idx/width
        // idx=x*width + y

        int nP = 0;

        double x2, y2, v2, diff;
        double[] pS;
        double errSSD = 0;
        for (int y = 0, idx = 0; y < image.length; ++y) {
            for (int x = 0; x < image[0].length; ++x, ++idx) {
                // this is just in case template is smaller than image.
                //TODO: review image and template dimensions and carry any restrictions up to input arguments
                if (x < 0 || x >= t[0].length || y < 0 || y >= t.length) {
                    continue;
                }
                x2 = dotXY(x, y, wX[idx]);
                y2 = dotXY(x, y, wY[idx]);
                if (x2 < 0 || Math.ceil(x2) >= image[0].length || y2 < 0 || Math.ceil(y2) >= image.length) {
                    continue;
                }
                // method expecting col major data so reverse the coords:
                v2 = imageProcessor.biLinearInterpolation(image, y2, x2);

                // error image is I(w(x;p)) - T(x).  not squared nor abs value of
                // TODO: this could be improved by considering the fraction of the pixel represented
                diff = v2 - t[y][x];

                errSSD += diff * diff;

                pS = steep[idx];

                // add to outPSum v2*pS
                for (int k = 0; k < outPSum.length; ++k) {
                    outPSum[k] += (diff * pS[k]);
                }

                ++nP;
            }
        }
        errSSD /= nP;
        errSSD = Math.sqrt(errSSD);

        //System.out.printf("errorSSD=%.3f\n", errSSD);

        return new double[]{nP, errSSD};
    }

    private static double dotXY(int x, int y, double[] v) {
        return x*v[0] + y*v[1] + v[2];
    }

    private static void wX(double[] p, double[] outWX) {
        outWX[0] = 1 + p[0];
        outWX[1] = p[2];
        outWX[2] = p[4];
    }

    private static void wY(double[] p, double[] outWY) {
        outWY[0] = p[1];
        outWY[1] = 1 + p[3];
        outWY[2] = p[5];
    }

    private static double[] inverseWParamX(double[] deltaP) {
        double[] out = new double[3];
        double det = (1+deltaP[0])*(1+deltaP[3]) - deltaP[1]*deltaP[2];
        out[0] = (-deltaP[0] - deltaP[0]*deltaP[3] + deltaP[1]*deltaP[2])/det;
        out[1] = -deltaP[2]/det;
        out[2] = (-deltaP[4] - deltaP[3]*deltaP[4] + deltaP[2]*deltaP[5])/det;
        return out;
    }

    private static double[] inverseWParamY(double[] deltaP) {
        double[] out = new double[3];
        double det = (1+deltaP[0])*(1+deltaP[3]) - deltaP[1]*deltaP[2];
        out[0] = -deltaP[1]/det;
        out[1] = (-deltaP[3] - deltaP[0]*deltaP[3] + deltaP[1]*deltaP[2])/det;
        out[2] = (-deltaP[5] - deltaP[0]*deltaP[5] + deltaP[1]*deltaP[4])/det;
        return out;
    }

    private static void computeWarpX(double[][] wX, double[] p, int nX, int nY) {
        for (int y = 0, idx=0; y < nY; ++y) {
            for (int x = 0; x < nX; ++x, ++idx) {
                wX[idx][0] = x*(1.+p[0]);
                wX[idx][1] = y*p[2];
                wX[idx][2] = p[4];
            }
        }
    }
    private static void updateWarpX(double[][] wX, double[] deltaP, int nX, int nY) {
        double[] params = inverseWParamX(deltaP); // length 3
        for (int y = 0, idx = 0; y < nY; ++y) {
            for (int x = 0; x < nX; ++x, ++idx) {
                // w(x; p) composition w(x; deltap)^-1
                wX[idx][0] += x * params[0];
                wX[idx][1] += y * params[1];
                wX[idx][2] += params[2];
            }
        }
    }

    private static void computeWarpY(double[][] wY, double[] p, int nX, int nY) {
        for (int y = 0, idx=0; y < nY; ++y) {
            for (int x = 0; x < nX; ++x, ++idx) {
                wY[idx][0] = x*p[1];
                wY[idx][1] = y*(1. + p[3]);
                wY[idx][2] = p[5];
            }
        }
    }

    private static void updateWarpY(double[][] wY, double[] deltaP, int nX, int nY) {
        double[] params = inverseWParamY(deltaP); // length 3
        for (int y = 0, idx = 0; y < nY; ++y) {
            for (int x = 0; x < nX; ++x, ++idx) {
                wY[idx][0] += x * params[0];
                wY[idx][1] += y * params[1];
                wY[idx][2] += params[2];
            }
        }
    }

    public static class Result {
        /**
         * the warped coordinates.  the data array is size [nx*ny X 2] where each row index is a composite index of x and y
         * accessed in these ways:
         * the index = (x * imagewidth) + y.
         * conversely can extract which row and column an index represents by using x = idx % imagewidth and y = idx/imagewidth.

         the values at each index are the warped x and warped y coordinates as an array of length 2

         warp size is [nx * ny X 2].
         */
        double[][] warp;

        double errorSD;

        public Result(double[][] w, double errorSSD) {
            this.warp = w;
            this.errorSD = errorSSD;
        }
    }
}
