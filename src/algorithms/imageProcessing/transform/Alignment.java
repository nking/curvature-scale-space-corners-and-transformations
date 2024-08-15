package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.StructureTensorD;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import no.uib.cipr.matrix.NotConvergedException;

import java.io.IOException;
import java.util.Arrays;

/**
 * various 2D and 3D alignment methods
 * 
 * @author nichole
 */
public class Alignment {

    /**
     * 2-D translation is translation in x and translation in y.
     * affine 2-D matrix can contain terms for rotation in the plane, scale, shear, and translation.
     */
    enum Type{
        TRANSLATION_2D, AFFINE_2D,
        AFFINE_3D_NOT_IMPLEMENTED
    }

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
    @return returns the square root of the sum of squared error image values where error image is the difference between
    the warped image and the template.
   */
    public static double inverseCompositional2DAffine(double[][] template, double[][] image, double[][] pInit,
                                                      int maxIter)
            throws NotConvergedException, IOException {

        if (pInit.length != 2 || pInit[0].length != 3) {
            throw new IllegalArgumentException("expecting pInit length to be 2, and pInit[0].length to be 3");
        }

        return inverseCompositional(template, image, pInit, maxIter, Type.AFFINE_2D);
    }

    public static double inverseCompositional2DTranslation(double[][] template, double[][] image,
                                                           double[] xYInit,
                                                      int maxIter)
            throws NotConvergedException, IOException {
        if (xYInit.length != 2) {
            throw new IllegalArgumentException("expecting xYInit to be length 2");
        }

        double[][] pInit = new double[][] {
                {1, 0, xYInit[0]},
                {0, 1, xYInit[1]}
        };
        double errSD = inverseCompositional(template, image, pInit, maxIter, Type.TRANSLATION_2D);
        xYInit[0] = pInit[0][2];
        xYInit[1] = pInit[1][2];
        return errSD;
    }

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
     @return returns the square root of the sum of squared error image values where error image is the difference between
     the warped image and the template.
    */
    public static double inverseCompositional(double[][] template, double[][] image, double[][] pInit,
                                                      int maxIter, Type type)
            throws NotConvergedException, IOException {

        if (pInit.length != 2 || pInit[0].length != 3) {
            throw new IllegalArgumentException("expecting pInit length to be 2, and pInit[0].length to be 3");
        }
        if (type == null) {
            throw new IllegalArgumentException("type cannot be null");
        }
        if (type.equals(Type.AFFINE_3D_NOT_IMPLEMENTED)) {
            throw new IllegalArgumentException("not yet implemented");
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

        //======  precompute  =======

        // gradients of T:   each is [nTR X nTC]
        StructureTensorD gT = new StructureTensorD(template, 1, false);
        double[][] gTX = gT.getDY();
        double[][] gTY = gT.getDX();

        // ==== create the steepest decent image of the template ====
        //steepest descent img = gradient * dWdP
        double[][] dTdWdp = null;
        if (type.equals(Type.AFFINE_2D)) {
            //[ntX*nTY X 6]
            dTdWdp = createSteepestDescentImageAffine2D(gTX, gTY);
        } else if (type.equals(Type.TRANSLATION_2D)) {
            //[ntX*nTY X 2]
            dTdWdp = createSteepestDescentImageTranslation2D(gTX, gTY);
        }

        // Hessian: [6 x nTX*nTY]*[nTX*nTY X 6] = [6x6] for affine 2D;  for trans2D [2x2]
        double[][] hessianTmplt = MatrixUtil.createATransposedTimesA(dTdWdp);

        double[][] invHessianTmplt = MatrixUtil.inverse(hessianTmplt);//MatrixUtil.pseudoinverseFullRowRank(hessianTmplt);

        //======  end precompute  =======

        double[][] warp = null;

        double eps = 1E-3;

        int len = 6;
        if (type.equals(Type.TRANSLATION_2D)) {
            len = 2;
        }

        double[] deltaP = new double[len];
        double[] pSum = new double[len];

        boolean converged = false;
        int nIter = 0;
        double[] nPAndErr = null;
        while (nIter < maxIter && !converged) {

            if (nIter == 0) {
                warp = createWarp0(pInit);
            } else {
                // update the warp
                if (type.equals(Type.AFFINE_2D)) {
                    updateWarp2DAffIC(warp, deltaP);
                } else if (type.equals(Type.TRANSLATION_2D)) {
                    updateWarp2DTransIC(warp, deltaP);
                }
            }

            if (type.equals(Type.AFFINE_2D)) {
                warp[0][0] += 1;
                warp[1][1] += 1;
            }

            // (1),(2),(7): warp I and subtract from T, then mult by steepest descent image, summing over all x
            nPAndErr = sumSteepestDescErrImageProduct(image, template, warp, dTdWdp, pSum);

            if (type.equals(Type.AFFINE_2D)) {
                warp[0][0] -= 1;
                warp[1][1] -= 1;
            }

            // compute deltaP as invHessianTmplt * pSum = [6 X 6] * [6 X 1] or [2X2]*[2X1]
            MatrixUtil.multiplyMatrixByColumnVector(invHessianTmplt, pSum, deltaP);

            double norm = MatrixUtil.lPSum(deltaP, 2);
            converged = norm < eps;

            ++nIter;
        }

        System.out.printf("nIter=%d, converged=%b\npInit=%s\ndeltaP=%s\n",
                nIter, converged, FormatArray.toString(pInit, "%.3f"), FormatArray.toString(deltaP, "%.3f"));

        System.arraycopy(warp[0], 0, pInit[0], 0, pInit[0].length);
        System.arraycopy(warp[1], 0, pInit[1], 0, pInit[1].length);

        System.out.printf("pFinal=%s\n", FormatArray.toString(pInit, "%.3f"));

        return nPAndErr[1];
    }

    private static double[][] createSteepestDescentImageAffine2D(double[][] gTX, double[][] gTY) {
        //[nTR*nTC X 6]
        // steepest descent images  gT * dWdP = [nTR X nTC] * [nTR*nTC X 6] = [nTR*nTC X 6] for x then for y
        // can format as [nTr X nTc X 6] or [nTr * nTc X 6].
        // will choose the latter to make the Hessian mult easier
        //steepest descent img = gradient * dWdP
        //                       [1x2] * [2x6]

        /*
        W(x; p) for a 2D affine warp:
                = [1+p[0]      p2  p4 ] * [x]
                  [  p[1]  1+p[3]  p5 ]   [y]
                                          [1]
        pdW/pdp = [x, 0, y, 0, 1, 0]
                  [0, x, 0, y, 0, 1]
        */

        int nYT = gTX.length;
        int nXT = gTY[0].length;
        int nTT = nXT*nYT;

        double[][] dTdWdp = new double[nTT][];
        double[][] tmp1 = new double[1][2];
        double[][] tmp2 = new double[2][6];
        for (int y = 0, idx = 0; y < nYT; ++y) {
            for (int x = 0; x < nXT; ++x, ++idx) {

                tmp1[0][0] = gTX[y][x];
                tmp1[0][1] = gTY[y][x];

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
        return dTdWdp;
    }

    private static double[][] createSteepestDescentImageTranslation2D(double[][] gTX, double[][] gTY) {
        //[nTR*nTC X 2]
        // steepest descent images  gT * dWdP = [nTR X nTC] * [nTR*nTC X 6] = [nTR*nTC X 6] for x then for y
        // can format as [nTr X nTc X 6] or [nTr * nTc X 6].
        // will choose the latter to make the Hessian mult easier
        //steepest descent img = gradient * dWdP
        //                       [1x2] * [2x2]
        int nYT = gTX.length;
        int nXT = gTY[0].length;
        int nTT = nXT*nYT;

        /*
        using the column major notation of the paper
        W is [ 1  0  p1]  * [x]
             [ 0  1  p2]    [y]
                            [1]
         dW/dP = [ dWx/dp1   dWx/dp2 ]
                 [ dWy/dp1   dWy/dp2 ]
               = [  1         0 ]
                 [  0         1 ]
           which is the identity matrix

         so the steepest descent image = each pixel's gTx, gTy

         */

        double[][] dTdWdp = new double[nTT][];
        for (int y = 0, idx = 0; y < nYT; ++y) {
            for (int x = 0; x < nXT; ++x, ++idx) {
                dTdWdp[idx] = new double[]{gTX[y][x], gTY[y][x]};
            }
        }
        return dTdWdp;
    }

    private static double[][] createWarp0(double[][] pInit) {
        return new double[][]{
                {pInit[0][0], pInit[0][1], pInit[0][2]},
                {pInit[1][0], pInit[1][1], pInit[1][2]},
                {0, 0, 1}
        };
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
     * @param warp
     * @param steep
     * @param outPSum the output sum of p over all points
     * @return double array holding the number of points used in the calculation and the quare root of the sum of squared error image values.
     * as double[]{nP, errSSD}
     */
    private static double[] sumSteepestDescErrImageProduct(double[][] image, double[][] t,
       double[][] warp, double[][] steep, double[] outPSum) {

        // since the error image is only used in the calculation of the product of the steepest descent image
        // and the error image, we will return the result instead of intermediate products

        Arrays.fill(outPSum, 0);

        ImageProcessor imageProcessor = new ImageProcessor();

        // steepX and steepY are filled by
        // x=idx%width, y=idx/width
        // idx=x*width + y

        int nP = 0;

        double[] xy2;
        double[] xy = new double[]{0,0,1};
        double v2, diff;
        double[] pS;
        double errSSD = 0;
        int idx;
        for (int y = 0; y < t.length; ++y) {
            for (int x = 0; x < t[0].length; ++x) {
                xy[0] = x;
                xy[1] = y;
                xy2 = MatrixUtil.multiplyMatrixByColumnVector(warp, xy);

                if (xy2[0] < 0 || Math.ceil(xy2[0]) >= image[0].length || xy2[1] < 0 || Math.ceil(xy2[1]) >= image.length) {
                    continue;
                }
                // method expecting col major data so reverse the coords:
                v2 = imageProcessor.biLinearInterpolation(image, xy2[1], xy2[0]);

                // error image is I(w(x;p)) - T(x).  not squared nor abs value of
                // TODO: this could be improved by considering the fraction of the pixel represented
                //diff = v2 - t[y][x]; // paper has I(W(x))-T(x)
                diff = t[y][x] - v2; // needed for 2D-translation

                errSSD += diff * diff;

                // idx = r*width + c
                idx = (int)Math.round(xy2[1])*image[0].length + (int)Math.round(xy2[0]);

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

    private static double[][] inverseWParam2DTrans(double[] deltaP) throws NotConvergedException {
        //TODO: write out the determinant and cofactor terms here
        return new double[][]{
                {1, 0, -deltaP[0]},
                {0, 1, -deltaP[1]},
                {0, 0, 1}
        };
    }

    private static double[][] inverseWParam2DAff(double[] deltaP) {
        double[][] out = new double[3][3];
        double det = (1+deltaP[0])*(1+deltaP[3]) - deltaP[1]*deltaP[2];
        out[0][0] = (-deltaP[0] - deltaP[0]*deltaP[3] + deltaP[1]*deltaP[2])/det;
        out[1][0] = -deltaP[1]/det;
        out[0][1] = -deltaP[2]/det;
        out[1][1] = (-deltaP[3] - deltaP[0]*deltaP[3] + deltaP[1]*deltaP[2])/det;
        out[0][2] = (-deltaP[4] - deltaP[3]*deltaP[4] + deltaP[2]*deltaP[5])/det;
        out[1][2] = (-deltaP[5] - deltaP[0]*deltaP[5] + deltaP[1]*deltaP[4])/det;
        out[2][2] = 1;

        //3D affine needs this:
        out[0][0] += 1;
        out[1][1] += 1;
        return out;
    }

    private static void updateWarp2DTransIC(double[][] warp, double[] deltaP) throws NotConvergedException {
        double[][] params = inverseWParam2DTrans(deltaP); //[3X3]

        double[][] tmp = MatrixUtil.multiply(warp, params);
        for (int i = 0; i < warp.length; ++i) {
            System.arraycopy(tmp[i], 0, warp[i], 0, tmp[i].length);
        }

        if (Math.abs(warp[2][2] - 1.) > 1E-7) {
            int t = 2;
        }
        warp[2][0] = 0;
        warp[2][1] = 0;
        warp[2][2] = 1;
    }

    private static void updateWarp2DAffIC(double[][] warp, double[] deltaP) throws NotConvergedException {
        double[][] params = inverseWParam2DAff(deltaP); //[3X3]

        // the 3D affine needs this:
        warp[0][0] += 1;
        warp[1][1] += 1;

        double[][] tmp = MatrixUtil.multiply(warp, params);
        for (int i = 0; i < warp.length; ++i) {
            System.arraycopy(tmp[i], 0, warp[i], 0, tmp[i].length);
        }

        warp[0][0] -= 1;
        warp[1][1] -= 1;

        if (Math.abs(warp[2][2] - 1.) > 1E-7) {
            //MatrixUtil.multiply(warp, 1./warp[2][2]);
            int t = 2;
        }
        warp[2][0] = 0;
        warp[2][1] = 0;
        warp[2][2] = 1;
    }

    private static boolean areTheSame(double[][] a, double[][] b, double tol) {
        for (int i = 0; i < a.length; ++i) {
            if (!areTheSame(a[i], b[i], tol)) {
                return false;
            }
        }
        return true;
    }
    private static boolean areTheSame(double[] a, double[] b, double tol) {
        for (int i = 0; i < a.length; ++i) {
           if (Math.abs(a[i] - b[i]) > tol) {
                return false;
            }
        }
        return true;
    }

}
