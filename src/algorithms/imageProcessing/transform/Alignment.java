package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.StructureTensorD;
import algorithms.matrix.MatrixUtil;
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
        TRANSLATION_2D, AFFINE_2D
    }

    /*
   NOT READY FOR USE...
   (might be better to use Reconstruction or CameraPose or epipolar methods if need more than translation solved
   as this algorithm provides rough solutions for small displacements for the simplest cases
   of no scaling, no shear, and only translation.)
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
                                                      int maxIter, double eps)
            throws NotConvergedException, IOException {

        if (pInit.length != 2 || pInit[0].length != 3) {
            throw new IllegalArgumentException("expecting pInit length to be 2, and pInit[0].length to be 3");
        }

        return inverseCompositional(template, image, pInit, maxIter, Type.AFFINE_2D, eps);
    }

    /** find the x and y translations that best align the template to the image.
   This alignment minimizes the re-projection error between the projected image and template
   (though the error is an LP1 sum rather than LP2).
   <pre>
   reference:
   Figure 4 of
   Baker & Matthews, 2016, CMU-RI-TR-02-16
   "Lucas-Kanade 20 Years On: A Unifying Framework: Part 1"
   </pre>
    @param template
    * @param image
    * @param xYInit an array of length 2 holding the x initial guess (can be 0) and y initial guess (can be 0) for
     *               translations of a pixel in template to match a pixel in image.
     *              NOTE that the final solution is written into this array.
     *              xYInit is an in-out variable.
    @param maxIter the maximum number of iterations to perform
    * @throws NotConvergedException
    * @throws IOException
    @return returns the square root of the sum of squared error image values where error image is the difference between
    the warped image and the template.  Note that the solution x and y offsets are written into xYInit.
   */
    public static double inverseCompositional2DTranslation(double[][] template, double[][] image,
                                                           double[] xYInit,
                                                      int maxIter, double eps)
            throws NotConvergedException, IOException {
        if (xYInit.length != 2) {
            throw new IllegalArgumentException("expecting xYInit to be length 2");
        }

        double[][] pInit = new double[][] {
                {1, 0, xYInit[0]},
                {0, 1, xYInit[1]}
        };
        double errSD = inverseCompositional(template, image, pInit, maxIter, Type.TRANSLATION_2D, eps);
        xYInit[0] = pInit[0][2];
        xYInit[1] = pInit[1][2];
        return errSD;
    }

    // Gauss-Newton.
    // NOT READY FOR USE
    public static double _inverseCompositional2DTranslationGN(double[][] template, double[][] image,
                                                              double[] xYInit,
                                                              int maxIter, double eps)
            throws NotConvergedException, IOException {
        if (xYInit.length != 2) {
            throw new IllegalArgumentException("expecting xYInit to be length 2");
        }

        double[][] pInit = new double[][] {
                {1, 0, xYInit[0]},
                {0, 1, xYInit[1]}
        };
        double errSD = _inverseCompositionalGN(template, image, pInit, maxIter, Type.TRANSLATION_2D, eps);
        xYInit[0] = pInit[0][2];
        xYInit[1] = pInit[1][2];
        return errSD;
    }

    /*
    NOT READY FOR USE
    refine pInit so that template and image are aligned by
    computing the warp that can be applied to the template to align it with the image.
    This alignment minimizes the re-projection error between the projected image and template
    (though the error is an LP1 sum rather than LP2).
    to get a good initial pInit, one could extract features from template and image and make correspondence list
    of matching features using tools such as ransac or assumptions of translation only etc depending upon context.
    <pre>
    reference2:
        Figure 4 of
        Baker & Matthews, 2016, CMU-RI-TR-02-16
        "Lucas-Kanade 20 Years On: A Unifying Framework: Part 1"

    the warp update change is adapted from Bouguet, "Pyramidal Implementation of the Affine Kanade Feature Tracker
    Descriptions of the Algorithm"
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
                                                      int maxIter, Type type, double eps)
            throws NotConvergedException, IOException {

        if (pInit.length != 2 || pInit[0].length != 3) {
            throw new IllegalArgumentException("expecting pInit length to be 2, and pInit[0].length to be 3");
        }
        if (type == null) {
            throw new IllegalArgumentException("type cannot be null");
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
            dTdWdp = createSteepestDescentImagesAffine2D(gTX, gTY);
        } else if (type.equals(Type.TRANSLATION_2D)) {
            //[ntX*nTY X 2]
            dTdWdp = createSteepestDescentImageTranslation2D(gTX, gTY);
        }

        // Hessian: [6 x nTX*nTY]*[nTX*nTY X 6] = [6x6] for affine 2D;  for trans2D [2x2]
        double[][] hessianTmplt = MatrixUtil.createATransposedTimesA(dTdWdp);

        //G must be invertible.  the template image must have gradient info in x and y at each point of evaluation
        //   for the affine transformation.  for that reason, need to perform the affine algorithm over patches
        // TODO: implement a patch based version <add Bouguet reference>
        double[][] invHessianTmplt = MatrixUtil.inverse(hessianTmplt);//MatrixUtil.pseudoinverseFullRowRank(hessianTmplt);

        //======  end precompute  =======

        double[][] warp = null;

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
                if (type.equals(Type.AFFINE_2D)) {
                    warp[0][0] += 1;
                    warp[1][1] += 1;
                }
            } else {
                // update the warp
                if (type.equals(Type.AFFINE_2D)) {
                    // as suggested by Bouguet, use the forward composition warp update instead:
                    updateWarp2DAffFC(warp, deltaP);
                    //updateWarp2DAffIC(warp, deltaP);
                } else if (type.equals(Type.TRANSLATION_2D)) {
                    updateWarp2DTransIC(warp, deltaP);
                }
            }

            // (1),(2),(7): warp I and subtract from T, then mult by steepest descent image, summing over all x
            nPAndErr = sumSteepestDescErrImageProduct(image, template, warp, dTdWdp, pSum, type);

            // compute deltaP as invHessianTmplt * pSum = [6 X 6] * [6 X 1] or [2X2]*[2X1]
            MatrixUtil.multiplyMatrixByColumnVector(invHessianTmplt, pSum, deltaP);

            double norm = MatrixUtil.lPSum(deltaP, 2);
            converged = norm < eps;

            ++nIter;
        }

        System.out.printf("nIter=%d, converged=%b\n", nIter, converged);

        System.arraycopy(warp[0], 0, pInit[0], 0, pInit[0].length);
        System.arraycopy(warp[1], 0, pInit[1], 0, pInit[1].length);

        return nPAndErr[1];
    }

    /**
     *NOT READY FOR USE
     *      <pre>
     *      adapted from the single image aligment within:
     *      Bouguet, "Pyramidal Implementation of the Affine Kanade Feature Tracker
     *      Descriptions of the Algorithm"
     *      </pre>
     * @param template
     * @param image
     * @param xYInit
     * @param maxIter
     * @param eps
     * @param yXKeypoints a 2-D array of size [nKeypoints x 2] of keypoints to calculate affine warp at.
     *                  Note that the keypoints should not be closer to the image bounds than patchHalfWidth.
     * @param patchHalfWidth a box of size 2*patchHalfWidth is used for the calculations are patch centers in yXPAthes.
     * @return
     * @throws NotConvergedException
     * @throws IOException
     */
    static Warps _inverseComposition2DTranslationKeypoints(double[][] template, double[][] image,
        double[] xYInit, int maxIter, double eps, int[][] yXKeypoints, int patchHalfWidth, Type type) throws NotConvergedException {
        if (xYInit.length != 2) {
            throw new IllegalArgumentException("expecting xYInit to be length 2");
        }

        double[][] pInit = new double[][] {
                {1, 0, xYInit[0]},
                {0, 1, xYInit[1]}
        };
        return _inverseCompositionKeypoints(template, image, pInit, maxIter, eps, yXKeypoints, patchHalfWidth, type);
    }

    /**
     * NOT READY FOR USE
     <pre>
     adapted from the single image aligment within:
     Bouguet, "Pyramidal Implementation of the Affine Kanade Feature Tracker
     Descriptions of the Algorithm"
     </pre>
     * @param template
     * @param image
     * @param pInit
     * @param maxIter
     * @param eps
     * @param yXKeypoints a 2-D array of size [nKeypoints x 2] of keypoints to calculate affine warp at.
     *                  Note that the keypoints should not be closer to the image bounds than patchHalfWidth.
     * @param patchHalfWidth a box of size 2*patchHalfWidth is used for the calculations are patch centers in yXPAthes.
     * @return
     * @throws NotConvergedException
     * @throws IOException
     */
    static Warps _inverseCompositionKeypoints(double[][] template, double[][] image,
        double[][] pInit, int maxIter, double eps, int[][] yXKeypoints, int patchHalfWidth,
                                                            Type type) throws NotConvergedException {

        if (pInit.length != 2 || pInit[0].length != 3) {
            throw new IllegalArgumentException("expecting pInit length to be 2, and pInit[0].length to be 3");
        }

        //======  precompute  =======

        // gradients of T:   each is [nTR X nTC]
        ImageProcessor imageProcessor = new ImageProcessor();
        double[][] gTX = MatrixUtil.copy(template);
        double[][] gTY = MatrixUtil.copy(template);
        imageProcessor.applyKernel1D(gTX, new double[]{-1, 0, +1}, false);
        imageProcessor.applyKernel1D(gTY, new double[]{-1, 0, +1}, true);
        MatrixUtil.multiply(gTX, 0.5);
        MatrixUtil.multiply(gTY, 0.5);

        int wT = template[0].length;
        //int hT = template.length;
        int nKeypoints = yXKeypoints.length;

        // array of single coordinates representing x and y coordinates
        int[] keypoints = makeReverseLookupArray(yXKeypoints, wT);

        // ==== create the steepest decent image of the template ====
        //steepest descent img = gradient * dWdP = [keypoints.length X 6]

        double[][] dTdWdp = createSteepestDescentImagesAffine2D(gTX, gTY, keypoints, patchHalfWidth);

        // a hessian for each patch.   [nKeypoints X 6 X 6]
        double[][][] hessianTmplt = createHessians(dTdWdp);

        // [nKeypoints X 6 X 6]
        double[][][] invHessianTmplt = createInverseHessians(hessianTmplt);
        //======  end precompute  =======

        double[][][] warp = new double[nKeypoints][3][3];

        double[][] deltaP = new double[nKeypoints][6];

        boolean converged = false;
        int nIter = 0;
        Result result = null;
        while (nIter < maxIter && !converged) {

            if (nIter == 0) {
                if (type.equals(Type.AFFINE_2D)) {
                    init2DAffineWarps(warp, pInit);
                } else {
                    init2DTranslationWarps(warp, pInit);
                }
            } else {
                if (type.equals(Type.AFFINE_2D)) {
                    update2DAffineWarps(warp, deltaP);
                } else {
                    update2DTranslationWarps(warp, deltaP);
                }
            }

            // (1),(2),(7): warp I and subtract from T, then mult by steepest descent image, summing over all x in patches
            result = sumSteepestDescErrImageProducts(image, template, warp, dTdWdp, keypoints, patchHalfWidth);

            converged = true;
            // compute deltaP as invHessianTmplt * pSum = [6 X 6] * [6 X 1] for each patch
            for (int pIdx = 0; pIdx < nKeypoints; ++pIdx) {
                double[] pSum = result.pSums[pIdx];

                //MatrixUtil.multiplyMatrixByColumnVector(invHessianTmplt[pIdx], pSum, deltaP[pIdx]);
                //double norm = MatrixUtil.lPSum(deltaP[pIdx], 2);
                //converged &= (norm < eps);

                // checking eqn 34
                // hessiant * deltaP[pIdx] - pSum
                //double[] chk = MatrixUtil.subtract(MatrixUtil.multiplyMatrixByColumnVector(hessianTmplt[pIdx], deltaP[pIdx]),
                //        pSum);

                //parameter c of eqn (91) as a substitute for 1/H. when diff=I-t, result is same for perfect data over these small patches
                double numer = MatrixUtil.dot(pSum, pSum);
                double denom = MatrixUtil.dot(pSum, MatrixUtil.multiplyMatrixByColumnVector(hessianTmplt[pIdx], pSum));
                double c = numer / denom;
                if (Double.isNaN(c)) {
                    //break;
                }
                double[] _deltaP = Arrays.copyOf(pSum, pSum.length);
                MatrixUtil.multiply(_deltaP, c);
                if (true) {
                    System.arraycopy(_deltaP, 0, deltaP[pIdx], 0, _deltaP.length);
                }

                double norm = MatrixUtil.lPSum(deltaP[pIdx], 2);
                converged &= (norm < eps);
            }

            ++nIter;
        }

        System.out.printf("nIter=%d, converged=%b\n", nIter, converged);

        Warps out = new Warps();
        out.warps = warp;
        out.ssd = result.ssd;
        return out;
    }


    private static void update2DTranslationWarps(double[][][] warps, double[][] deltaPs) throws NotConvergedException {
        for (int i = 0; i < warps.length; ++i) {
            double[][] params = inverseWParam2DTrans(deltaPs[i]); //[3X3]
            double[][] tmp = MatrixUtil.multiply(warps[i], params);
            for (int i2 = 0; i2 < warps[i].length; ++i2) {
                System.arraycopy(tmp[i2], 0, warps[i][i2], 0, tmp[i2].length);
            }
        }
    }

    private static void update2DAffineWarps(double[][][] warps, double[][] deltaPs) {

        for (int i = 0; i < warps.length; ++i) {
            double[][] params = new double[][]{
                    {1 + deltaPs[i][0], deltaPs[i][2], deltaPs[i][4]},
                    {deltaPs[i][1], 1 + deltaPs[i][3], deltaPs[i][5]},
                    {0, 0, 1}
            };
            double[][] invParams = MatrixUtil.inverse(params);

            // forward composition:
            double[][] tmp1 = MatrixUtil.multiply(warps[i], params);
            // inverse composition
            double[][] tmp2 = MatrixUtil.multiply(warps[i], invParams);
            double[][] tmp = tmp1;
            for (int j = 0; j < tmp.length; ++j) {
                System.arraycopy(tmp[j], 0, warps[i][j], 0, tmp[j].length);
            }
        }
    }

    private static void init2DAffineWarps(double[][][] warps, double[][] pInit) {
        for (int i = 0; i < warps.length; ++i) {
            for (int j = 0; j < pInit.length; ++j) {
                System.arraycopy(pInit[j], 0, warps[i][j], 0, pInit[j].length);
            }
            warps[i][0][0] += 1;
            warps[i][1][1] += 1;
            warps[i][2][2] = 1;
        };
    }

    private static void init2DTranslationWarps(double[][][] warps, double[][] pInit) {
        for (int i = 0; i < warps.length; ++i) {
            for (int j = 0; j < pInit.length; ++j) {
                System.arraycopy(pInit[j], 0, warps[i][j], 0, pInit[j].length);
            }
            warps[i][2][2] = 1;
        };
    }

    protected static class Result {
        double[][] pSums;
        double ssd;
    }

    protected static class Warps {
        double[][][] warps;
        double ssd;
    }

    private static double[][][] createHessians(double[][] dTdWdp) {
        double[][][] h = new double[dTdWdp.length][dTdWdp[0].length][];
        for (int pIdx = 0; pIdx < dTdWdp.length; ++pIdx) {
            double[] steep = dTdWdp[pIdx];
            h[pIdx] = MatrixUtil.outerProduct(steep, steep);
        }
        return h;
    }

    private static double[][][] createInverseHessians(double[][][] h) throws NotConvergedException {
        double[][][] invH = new double[h.length][h[0].length][];
        for (int pIdx = 0; pIdx < h.length; ++pIdx) {
            double[][] hI = h[pIdx];
            invH[pIdx] = MatrixUtil.inverse(hI);
            if (Double.isNaN(invH[pIdx][0][0])) {
                invH[pIdx] = MatrixUtil.pseudoinverseRankDeficient(hI);
            }
        }
        return h;
    }

    private static double[][] createSteepestDescentImagesAffine2D(double[][] gTX, double[][] gTY, int[] patchIdxs,
                                                                  int patchHalfWidth) {
        int nYT = gTX.length;
        int nXT = gTY[0].length;

        //steepest descent img = gradient * dWdP = [patchIdxs.length X 6]
        // for each x sum over range -patchHalfWidth to +patchHalfWidth
        double[][] dTdWdp = new double[patchIdxs.length][6];

        double[][] tmp1 = new double[1][2];
        double[][] tmp2 = new double[2][6];
        double[][] tmp3 = new double[1][6];
        int idx, y, x, x2, y2;
        int[] yx = new int[2];
        for (int pIdx = 0; pIdx < patchIdxs.length; ++pIdx) {
            idx = patchIdxs[pIdx];
            idxToRowCol(idx, nXT, yx);
            y = yx[0];
            x = yx[1];

            // sum over box around patch center
            for (y2 = y - patchHalfWidth; y2 <= y + patchHalfWidth; ++y2) {
                for (x2 = x - patchHalfWidth; x2 <= x + patchHalfWidth; ++x2) {
                    tmp1[0][0] = getBoundaryCorrectedValue(gTX, y2, x2);
                    tmp1[0][1] = getBoundaryCorrectedValue(gTY, y2, x2);

                    TODO: here, the coordinates x2 and y2 might need to be the offsets from the keypoint
                    and might need to use a weight for the distance from the keypoint
                    tmp2[0][0] = x2;
                    tmp2[0][2] = y2;
                    tmp2[0][4] = 1;
                    tmp2[1][1] = x2;
                    tmp2[1][3] = y2;
                    tmp2[1][5] = 1;
                    /*
                    [ gx[pix]  gy[pix] ]  * [x2  0  y2  0  1  0]
                                            [ 0  x2  0  y2 0  1]
                      = [gx*x2  gy*y2  gx*y2  gy*y2  gx  gy]
                     */

                    // idx=x*width + y;  x=idx%width; y=idx/width
                    //[1X6]
                    MatrixUtil.multiply(tmp1, tmp2, tmp3);
                    for (int k = 0; k < 6; ++k) {
                        dTdWdp[pIdx][k] += tmp3[0][k];
                    }
                }
            }
        }
        return dTdWdp;
    }

    private static double getBoundaryCorrectedValue(double[][] im, int y, int x) {
        if (x < 0) {
            x = 0;
        }
        if (y < 0) {
            y = 0;
        }
        if (x >= im[0].length) {
            x = im[0].length - 1;
        }
        if (y >= im.length) {
            y = im.length - 1;
        }
        return im[y][x];
    }
    private static double getBoundaryCorrectedValue(double[][] im, double y, double x) {
        return getBoundaryCorrectedValue(im, (int)Math.round(y), (int)Math.round(x));
    }

    private static int getBoundaryCorrectedCoord(int coord, int dimLen) {
        if (coord < 0) {
            return 0;
        }
        if (coord >= dimLen) {
            return dimLen - 1;
        }
        return coord;
    }

    // NOT READY FOR USE
    public static double _inverseCompositionalGN(double[][] template, double[][] image, double[][] pInit,
                                                 int maxIter, Type type, double eps)
            throws NotConvergedException, IOException {

        if (pInit.length != 2 || pInit[0].length != 3) {
            throw new IllegalArgumentException("expecting pInit length to be 2, and pInit[0].length to be 3");
        }
        if (type == null) {
            throw new IllegalArgumentException("type cannot be null");
        }

        //======  precompute  =======

        // gradients of T:   each is [nTR X nTC]
        //StructureTensorD gT = new StructureTensorD(template, 1, false);
        //double[][] gTX = gT.getDY();
        //double[][] gTY = gT.getDX();
        ImageProcessor imageProcessor = new ImageProcessor();
        double[][] gTX = MatrixUtil.copy(template);
        double[][] gTY = MatrixUtil.copy(template);
        // {0, -1, 1} should match the paper diff in method sumSteepestDescErrImageProduct if change back to it
        imageProcessor.applyKernel1D(gTX, new double[]{0, 1, -1}, false);
        imageProcessor.applyKernel1D(gTY, new double[]{0, 1, -1}, true);

        // ==== create the steepest decent image of the template ====
        //steepest descent img = gradient * dWdP
        double[][] dTdWdp = null;
        if (type.equals(Type.AFFINE_2D)) {
            //[ntX*nTY X 6]
            dTdWdp = createSteepestDescentImagesAffine2D(gTX, gTY);
        } else if (type.equals(Type.TRANSLATION_2D)) {
            //[ntX*nTY X 2]
            dTdWdp = createSteepestDescentImageTranslation2D(gTX, gTY);
        }

        // Fig 11, step (6) of <add reference for Baker and Matthews 20 Years on..."
        // the template hessian here is summation_{all points x} pdG^2/pdp^2 =
        // Hessian: [6 x nTX*nTY]*[nTX*nTY X 6] = [6x6] for affine 2D;  for trans2D [2x2] = dTdWdp^T * dTdWdp
        double[][] hessianTmplt = MatrixUtil.createATransposedTimesA(dTdWdp);

        // used while debugging:
        double[][] _invHessianTmplt = MatrixUtil.inverse(hessianTmplt);//MatrixUtil.pseudoinverseFullRowRank(hessianTmplt);

        //======  end precompute  =======

        double[][] warp = null;

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
                if (type.equals(Type.AFFINE_2D)) {
                    warp[0][0] += 1;
                    warp[1][1] += 1;
                }
            } else {
                // update the warp
                if (type.equals(Type.AFFINE_2D)) {
                    // using the forward compositional update instead as suggested by Bouguet
                    updateWarp2DAffFC(warp, deltaP);
                    //updateWarp2DAffIC(warp, deltaP);
                } else if (type.equals(Type.TRANSLATION_2D)) {
                    updateWarp2DTransIC(warp, deltaP);
                }
            }

            // (1),(2),(7): warp I and subtract from T, then mult by steepest descent image, summing over all x
            nPAndErr = sumSteepestDescErrImageProduct(image, template, warp, dTdWdp, pSum, type);

            // compute deltaP as c * pSum = [6 X 1] or [2X1]
            System.arraycopy(pSum, 0, deltaP, 0, pSum.length);
            if (nPAndErr[0] == 0) {
                break;
            }
            //parameter c of eqn (91) as a substitute for 1/H
            double numer = MatrixUtil.dot(pSum, pSum);
            double denom = MatrixUtil.dot(pSum, MatrixUtil.multiplyMatrixByColumnVector(hessianTmplt, pSum));
            double c = numer / denom;
            if (Double.isNaN(c)) {
                break;
            }
            MatrixUtil.multiply(deltaP, c);
            // does invHess / c = I? no... diagonalization of invHess fails too

            //compare to newton update:
            double[] _deltaP = MatrixUtil.multiplyMatrixByColumnVector(_invHessianTmplt, pSum);

            double norm = MatrixUtil.lPSum(deltaP, 2);

            converged = norm < eps;

            ++nIter;
        }

        System.out.printf("nIter=%d, converged=%b\n", nIter, converged);

        System.arraycopy(warp[0], 0, pInit[0], 0, pInit[0].length);
        System.arraycopy(warp[1], 0, pInit[1], 0, pInit[1].length);

        return nPAndErr[1];
    }


    private static double[][] createSteepestDescentImagesAffine2D(double[][] gTX, double[][] gTY) {
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

    //[nx*ny X [2x3]]
    private static double[][][] createSteepestDescentImageAffine2D_23(double[][] gTX, double[][] gTY) {
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

        double[][][] dTdWdp = new double[nTT][][];
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
                //[1X6] shape into [2x3]
                //dTdWdp[idx] = MatrixUtil.multiply(tmp1, tmp2)[0];
                double[] tmp3 = MatrixUtil.multiply(tmp1, tmp2)[0];
                dTdWdp[idx] = new double[][]{
                        {tmp3[0], tmp3[2], tmp3[4]},
                        {tmp3[1], tmp3[3], tmp3[5]}
                };
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
     * steps (1),(2),(7) of Figure 4 of Baker-Matthews, 2002 or steps (1), (2), (7), (8) of Figure 11.
     * <pre>
     *   Baker-Matthews, 2002, "Lucas-Kanade 20 Years On: A Unifying Framework: Part 1"
     * </pre>
     * @param image
     * @param t
     * @param warp
     * @param steep steepest descent image, dTdWdp, i.e. gradient of T * dW/dp
     * @param outPSum the output sum of p over all points, i.e. sum over all points of (steep^T * (I(w(x;p)) - T(x))
     * @return double array holding the number of points used in the calculation,
     * the square root of the sum of squared error image values
     * as double[]{nP, errSSD}.
     */
    private static double[] sumSteepestDescErrImageProduct(double[][] image, double[][] t,
        double[][] warp, double[][] steep, double[] outPSum, Type type) {

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
        double v2, diff=0;
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
                if (type.equals(Type.AFFINE_2D)) {
                    diff = v2 - t[y][x];
                } else if (type.equals(Type.TRANSLATION_2D)) {
                    diff = t[y][x] - v2;
                }

                errSSD += diff * diff;

                // idx = r*width + c
                idx = (int)Math.round(xy2[1])*image[0].length + (int)Math.round(xy2[0]);

                pS = steep[idx];

                // add to outPSum v2*pS
                for (int k = 0; k < outPSum.length; ++k) {
                    // TODO: apply normalization
                    outPSum[k] += (diff * pS[k]);
                }

                ++nP;
            }
        }
        errSSD /= nP;
        errSSD = Math.sqrt(errSSD);


        //System.out.printf("errorSSD=%.3f\n", errSSD);
        if (nP == 0) {
            return new double[]{nP, Double.MAX_VALUE};
        }
        return new double[]{nP, errSSD};
    }

    private static Result sumSteepestDescErrImageProducts(double[][] image, double[][] template,
                                                          double[][][] warps,
                                                           double[][] dTdWdp, int[] patchIdxs, int patchHalfWidth) {

        int nPatches = patchIdxs.length;

        double[][] outPSums = new double[nPatches][dTdWdp[0].length];

        ImageProcessor imageProcessor = new ImageProcessor();

        int nYT = template.length;
        int nXT = template[0].length;

        double[] xy2;
        double[] xy = new double[]{0,0,1};
        double v2, diff=0;
        double[] pS;
        double errSSD = 0;
        int idx, y, x, x2, y2;
        int[] yx = new int[2];
        for (int pIdx = 0; pIdx < patchIdxs.length; ++pIdx) {
            idx = patchIdxs[pIdx];
            idxToRowCol(idx, nXT, yx);
            y = yx[0];
            x = yx[1];

            pS = dTdWdp[pIdx];

            // sum over box around patch center
            for (y2 = y - patchHalfWidth; y2 <= y + patchHalfWidth; ++y2) {
                for (x2 = x - patchHalfWidth; x2 <= x + patchHalfWidth; ++x2) {

                    xy[0] = getBoundaryCorrectedCoord(x2, nXT);
                    xy[1] = getBoundaryCorrectedCoord(y2, nYT);
                    xy2 = MatrixUtil.multiplyMatrixByColumnVector(warps[pIdx], xy);

                    if (xy2[0] < 0 || Math.ceil(xy2[0]) >= image[0].length || xy2[1] < 0 || Math.ceil(xy2[1]) >= image.length) {
                        //TODO: consider alternatives
                        continue;
                    }
                    // method expecting col major data so reverse the coords:
                    v2 = imageProcessor.biLinearInterpolation(image, xy2[1], xy2[0]);

                    diff = v2 - getBoundaryCorrectedValue(template, xy[1], xy[0]);

                    errSSD += diff * diff;

                    for (int k = 0; k < outPSums[pIdx].length; ++k) {
                        outPSums[pIdx][k] += (diff * pS[k]);
                    }
                }
            } // end loop over patch
        }

        errSSD /= nPatches;
        errSSD = Math.sqrt(errSSD);

        Result result = new Result();
        result.pSums = outPSums;
        result.ssd = errSSD;
        return result;
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
        // determinant of cofactor p[5]
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
    }

    private static void _updateWarp2DAffIC(double[][] warp, double[] deltaP) throws NotConvergedException {
        double[][] params = new double[][]{
                {deltaP[0] + 1, deltaP[2], deltaP[4]},
                {deltaP[1], deltaP[3] + 1, deltaP[5]},
                {0, 0, 1}
        };
        params = MatrixUtil.inverse(params);

        double[][] tmp = MatrixUtil.multiply(warp, params);
        for (int i = 0; i < warp.length; ++i) {
            System.arraycopy(tmp[i], 0, warp[i], 0, tmp[i].length);
        }
    }

    private static void updateWarp2DAffFC(double[][] warp, double[] deltaP) throws NotConvergedException {

        // for affine, Bouguet uses forward composition:
        // instead of composition w = w * w(delta)^-1 it is w = w * w(delta).

        //double[][] params = inverseWParam2DAff(deltaP); //[3X3]
        double[][] params = new double[][]{
                {1+deltaP[0], deltaP[2], deltaP[4]},
                {deltaP[1], 1+deltaP[3], deltaP[5]},
                {0, 0, 1}
        };

        double[][] tmp = MatrixUtil.multiply(warp, params);

        for (int i = 0; i < warp.length; ++i) {
            System.arraycopy(tmp[i], 0, warp[i], 0, tmp[i].length);
        }
    }
    private static void updateWarp2DAffIC(double[][] warp, double[] deltaP) throws NotConvergedException {

        //double[][] params = inverseWParam2DAff(deltaP); //[3X3]
        double[][] params = new double[][]{
                {1+deltaP[0], deltaP[2], deltaP[4]},
                {deltaP[1], 1+deltaP[3], deltaP[5]},
                {0, 0, 1}
        };
        double[][] invParams = MatrixUtil.inverse(params);

        if (Double.isNaN(invParams[0][0])) {
            invParams = MatrixUtil.pseudoinverseRankDeficient(params);
        }
        double[][] tmp = MatrixUtil.multiply(warp, invParams);

        for (int i = 0; i < warp.length; ++i) {
            System.arraycopy(tmp[i], 0, warp[i], 0, tmp[i].length);
        }
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


    /*
    representation of a 2-dimensional point as 1 index:
    idx = row*width + col
    col = idx % width
    row = idx / width
     */
    private static void idxToRowCol(int idx, int imgWidth, int[] out) {
        out[0] = idx / imgWidth;
        out[1] = idx % imgWidth;
    }

    private static int rowColToIdx(int row, int col, int imgWidth) {
        return row * imgWidth + col;
    }

    private static int[] makeReverseLookupArray(int[][] yx, int imgWidth) {
        int nP = yx.length;
        int[] out = new int[nP];
        for (int i = 0; i < nP; ++i) {
            out[i] = rowColToIdx(yx[i][0], yx[i][1], imgWidth);
        }
        return out;
    }

}
