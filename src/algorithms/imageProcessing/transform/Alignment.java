package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.Kernel1DHelper;
import algorithms.imageProcessing.StructureTensorD;
import algorithms.matrix.MatrixUtil;
import no.uib.cipr.matrix.NotConvergedException;

import java.io.IOException;
import java.util.Arrays;

/**
 * various 2D alignment methods.


     One could use a pyramidal approach (a wrapper around invocations of this method) to find the solution with
     the smallest re-projection error over all scales.

     <pre>
     for 2D affine warps:
     one could decompose a warp matrix into the 7 unknown variables:
     for 2D affine, rotation is 1 variable, shear is 2 variables, translation is 2 variables, scale is 2 variables,
     so the total is 7 unknown variables encapsulated in the 6 numbers of the affine projection matrix.

     here are the separate matrices for the 2D transformations:
     [1  0  t_x ]  for translation
     [0  1  t_y ]
     [0  0   1  ]

     [s_x  0    0 ]  for scale
     [0    s_y  0 ]
     [0    0    1 ]

     [cos(th)  -sin(th)  0 ]  for rotation
     [sin(th)  cos(th)   0 ]
     [0         0        1 ]

     [0     sh_x  0 ]  for shear
     [sh_y  0     0 ]
     [0     0     1  ]

     a[0][0] and a[1][1] terms must be <= 1 for the rotation component in 2D affine.

     For 2D affine, one could use a pyramidal approach to constrain the x and y scales within some feasible range,
     reducing the number of variables to solve for to 5 and one could isolate the upper 2x2 left of the
     warp (2D affine) matrix to these terms:

     [ sh_x * sin(th)   sh_x * cos(th) ]
     [ sh_y * cos(th)  -sh_y * sin(th)  ]
     One can form ratios of terms to solve for th, then plug to solve for sh_x and sh_y.

     tx and t_y are solved for in the last column of the warp with the caveat that the scale has been factored into them.
     </pre>

     If one needs a 2D affine solution, and this is not producing good results, consider use of
     Reconstruction or CameraPose or epipolar methods.

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

    protected static double[] kernel = new double[]{0, 1, -1};

    /*
   2D affine alignment.  Note that this method operates over the entire image, which might not be the best strategy for
   2D-affine calculations.  Instead, one might prefer to use keypoints chosen by Harris response function etc.
   and use the method inverseCompositionKeypointsCpImgs(...) with those keypoints.

   <pre>
   reference:
   Figure 4 of
   Baker & Matthews, 2016, CMU-RI-TR-02-16
   "Lucas-Kanade 20 Years On: A Unifying Framework: Part 1"
   </pre>
    @param template
    * @param image
    @param pInit the initial projection matrix of size [2 rows X 3 columns].
     It should be composed of { {1 + affine[0][0], affine[0][1], transX},
                                 {affine[1][0], 1 + affine[1][1], transY}}.
                                 e.g. {{1,0,0},{0,1,0}}
    @param maxIter the maximum number of iterations to perform
    * @throws NotConvergedException
    * @throws IOException
    @return returns the square root of the sum of squared error image values where error image is the difference between
    the warped image and the template, and returns the number of iterations
   */
    public static double[] inverseCompositional2DAffine(double[][] template, double[][] image, double[][] pInit,
                                                      int maxIter, double eps)
            throws NotConvergedException, IOException {

        if (pInit.length != 2 || pInit[0].length != 3) {
            throw new IllegalArgumentException("expecting pInit length to be 2, and pInit[0].length to be 3");
        }

        return inverseCompositional(template, image, pInit, maxIter, Type.AFFINE_2D, eps);
    }

    /** find the x and y translations that best align the template to the image.
    Alternatively, one could use keypoints chosen by Harris response function etc. along with
     inverseCompositionKeypoints(...)
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
     Also returns the number of iterations used.
   */
    public static double[] inverseCompositional2DTranslation(double[][] template, double[][] image,
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
        double[] errSD = inverseCompositional(template, image, pInit, maxIter, Type.TRANSLATION_2D, eps);
        xYInit[0] = pInit[0][2];
        xYInit[1] = pInit[1][2];
        return errSD;
    }

    /*
    2D alignment over full images.
    Note that for type 2D affine alignment, this might not be the best strategy.
    Instead, one might prefer to use keypoints chosen by Harris response function etc.
   and for 2D affine use the method inverseCompositionKeypointsCpImgs(...) with those keypoints,
   else for 2D translation use inverseCompositionKeypoints(...)

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
     @param pInit the initial projection matrix of size [2 rows X 3 columns].
     It should be composed of { {1 + affine[0][0], affine[0][1], transX},
                                 {affine[1][0], 1 + affine[1][1], transY}}.
                                 e.g. {{1,0,0},{0,1,0}}
     @param maxIter the maximum number of iterations to perform
     * @throws NotConvergedException
     * @throws IOException
     @return returns the square root of the sum of squared error image values where error image is the difference between
     the warped image and the template, and returns the number of iterations
    */
    public static double[] inverseCompositional(double[][] template, double[][] image, double[][] pInit,
                                                      int maxIter, Type type, double eps)
            throws NotConvergedException, IOException {

        if (pInit.length != 2 || pInit[0].length != 3) {
            throw new IllegalArgumentException("expecting pInit length to be 2, and pInit[0].length to be 3");
        }
        if (type == null) {
            throw new IllegalArgumentException("type cannot be null");
        }

        /*
          TODO: consider improving the error image construction to make it quasi-convex
              from Ke & Kanade, "Quasi-convex Optimization for Robust Geometric Reconstruction"
                  convex function of y:  any norm function g(y) = ||y||_l.
                  affine function of X: h(X) = (p_u(X), p_v(X)).
                  The composition of convex function g and affine function h is a convex function: g ◦ h
                  ==> p(X) = (g ◦ h)(X) is a convex function and p(X) >= 0
         */

        //======  precompute  =======

        // gradients of T:   each is [nTR X nTC]

        //StructureTensorD gT = new StructureTensorD(template, 1, false);
        //double[][] gTX = gT.getDY();
        //double[][] gTY = gT.getDX();
        //gTX = strictColumnDiff(template);
        //gTY = strictRowDiff(template);

        ImageProcessor imageProcessor = new ImageProcessor();
        double[][] gTX = MatrixUtil.copy(template);
        double[][] gTY = MatrixUtil.copy(template);
        // {0, -1, 1} should match the paper diff in method sumSteepestDescErrImageProduct if change back to it
        imageProcessor.applyKernel1D(gTX, kernel, false);
        imageProcessor.applyKernel1D(gTY, kernel, true);

        // ==== create the steepest decent image of the template ====
        //steepest descent img = gradient * dWdP
        double[][] dTdWdp = null;
        if (type.equals(Type.AFFINE_2D)) {
            //[ntX*nTY X 6]
            dTdWdp = createSteepestDescentImagesAffine2D(gTX, gTY);
            int COMPARE1 = 2;
        } else if (type.equals(Type.TRANSLATION_2D)) {
            //[ntX*nTY X 2]
            dTdWdp = createSteepestDescentImageTranslation2D(gTX, gTY);
        }

        // Hessian: [6 x nTX*nTY]*[nTX*nTY X 6] = [6x6] for affine 2D;  for trans2D [2x2]
        double[][] hessianTmplt = MatrixUtil.createATransposedTimesA(dTdWdp);

        //G must be invertible.  the template image must have gradient info in x and y at each point of evaluation
        //   for the affine transformation.  for that reason, need to perform the affine algorithm over patches
        double[][] invHessianTmplt = MatrixUtil.inverse(hessianTmplt);
        if (Double.isNaN(invHessianTmplt[0][0])) {
            invHessianTmplt = MatrixUtil.pseudoinverseRankDeficient(hessianTmplt);
        }

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
                warp = init2DTranslationWarp(pInit);
            } else {
                // update the warp
                if (type.equals(Type.AFFINE_2D)) {
                    // as suggested by Bouguet, use the forward composition warp update instead:
                    update2DAffFCWarp(warp, deltaP);
                    //updateWarp2DAffIC(warp, deltaP);
                } else if (type.equals(Type.TRANSLATION_2D)) {
                    update2DTransICWarp(warp, deltaP);
                }
            }

            // (1),(2),(7): warp I and subtract from T, then mult by steepest descent image, summing over all x
            nPAndErr = sumSteepestDescErrImageProduct(image, template, warp, dTdWdp, pSum, type);

            int COMPARE1 = 2;

            // compute deltaP as invHessianTmplt * pSum = [6 X 6] * [6 X 1] or [2X2]*[2X1]
            MatrixUtil.multiplyMatrixByColumnVector(invHessianTmplt, pSum, deltaP);

            double norm = MatrixUtil.lPSum(deltaP, 2);
            converged = norm < eps;

            ++nIter;
        }

        //System.out.printf("nIter=%d, converged=%b\n", nIter, converged);

        System.arraycopy(warp[0], 0, pInit[0], 0, pInit[0].length);
        System.arraycopy(warp[1], 0, pInit[1], 0, pInit[1].length);

        return new double[]{nPAndErr[1], nIter};
    }

    /*
   2D alignment for a single keypoint whose coordinates are x, y.   The calculations
   are performed within a window with hX half width in the x dimension and hY half width in the y dimension.
    Note that for type 2D affine alignment, this might not be the best strategy, but instead,
    one might prefer to use the method inverseCompositionKeypointsCpImgs(...) .

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
    @param x the x coordinate of the feature to determine warp alignment for
    @param y the y coordinate of the feature to determine warp alignment for
    @param hX the half width in x dimension of a window around x to use in calculations.  This half width must be large enough to
    include the possible warped point.
    @param hY the half width in y dimension of a window around y to use in calculations.  This half width must be large enough to
    include the possible warped point.
    @param pInit the initial projection matrix of size [2 rows X 3 columns].
     It should be composed of { {1 + affine[0][0], affine[0][1], transX},
                                 {affine[1][0], 1 + affine[1][1], transY}}.
                                 e.g. {{1,0,0},{0,1,0}}
    @param maxIter the maximum number of iterations to perform
    * @throws NotConvergedException
    * @throws IOException
    @return returns the square root of the sum of squared error image values where error image is the difference between
    the warped image and the template, and returns the number of iterations
   */
    public static double[] inverseCompositional(double[][] template, double[][] image,
                                                int x, int y, int hX, int hY,
                                                double[][] pInit,
                                                int maxIter, Type type, double eps)
            throws NotConvergedException, IOException {

        if (pInit.length != 2 || pInit[0].length != 3) {
            throw new IllegalArgumentException("expecting pInit length to be 2, and pInit[0].length to be 3");
        }
        if (type == null) {
            throw new IllegalArgumentException("type cannot be null");
        }

        /*
          TODO: consider improving the error image construction to make it quasi-convex
              from Ke & Kanade, "Quasi-convex Optimization for Robust Geometric Reconstruction"
                  convex function of y:  any norm function g(y) = ||y||_l.
                  affine function of X: h(X) = (p_u(X), p_v(X)).
                  The composition of convex function g and affine function h is a convex function: g ◦ h
                  ==> p(X) = (g ◦ h)(X) is a convex function and p(X) >= 0
         */

        //======  precompute  =======

        // gradients of T:   each is [nTR X nTC]
        /*
        StructureTensorD gT = new StructureTensorD(template, 1, false);
        double[][] gTX = gT.getDY();
        double[][] gTY = gT.getDX();
         */
        /*
        ImageProcessor imageProcessor = new ImageProcessor();
            double[][] gTX = MatrixUtil.copy(template);
            double[][] gTY = MatrixUtil.copy(template);
            // {0, -1, 1} should match the paper diff in method sumSteepestDescErrImageProduct if change back to it
            imageProcessor.applyKernel1D(gTX, kernel, false);
            imageProcessor.applyKernel1D(gTY, kernel, true);
         */

        // ==== create the steepest decent image of the template ====
        //steepest descent img = gradient * dWdP
        double[][] dTdWdp = null;
        if (type.equals(Type.AFFINE_2D)) {
            //[ntX*nTY X 6]
            dTdWdp = createSteepestDescentImagesAffine2D(template, x-hX, x+hX, y-hY, y+hY);
            int COMPARE2 = 2;
        } else if (type.equals(Type.TRANSLATION_2D)) {
            //[ntX*nTY X 2]
            dTdWdp = createSteepestDescentImageTranslation2D(template, x-hX, x+hX, y-hY, y+hY);
        }

        // Hessian: [6 x nTX*nTY]*[nTX*nTY X 6] = [6x6] for affine 2D;  for trans2D [2x2]
        double[][] hessianTmplt = MatrixUtil.createATransposedTimesA(dTdWdp);

        //G must be invertible.  the template image must have gradient info in x and y at each point of evaluation
        //   for the affine transformation.  for that reason, need to perform the affine algorithm over patches
        double[][] invHessianTmplt = MatrixUtil.inverse(hessianTmplt);//MatrixUtil.pseudoinverseFullRowRank(hessianTmplt);
        if (Double.isNaN(invHessianTmplt[0][0])) {
            invHessianTmplt = MatrixUtil.pseudoinverseRankDeficient(hessianTmplt);
        }

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
                warp = init2DTranslationWarp(pInit);
            } else {
                // update the warp
                if (type.equals(Type.AFFINE_2D)) {
                    // as suggested by Bouguet, use the forward composition warp update instead:
                    update2DAffFCWarp(warp, deltaP);
                    //updateWarp2DAffIC(warp, deltaP);
                } else if (type.equals(Type.TRANSLATION_2D)) {
                    update2DTransICWarp(warp, deltaP);
                }
            }

            // (1),(2),(7): warp I and subtract from T, then mult by steepest descent image, summing over all x
            nPAndErr = sumSteepestDescErrImageProduct(image, template, warp, dTdWdp, pSum, type, x-hX, x+hX, y-hY, y+hY);

            int COMPARE2 = 2;

            // compute deltaP as invHessianTmplt * pSum = [6 X 6] * [6 X 1] or [2X2]*[2X1]
            MatrixUtil.multiplyMatrixByColumnVector(invHessianTmplt, pSum, deltaP);

            double norm = MatrixUtil.lPSum(deltaP, 2);
            converged = norm < eps;

            ++nIter;
        }

        //System.out.printf("nIter=%d, converged=%b\n", nIter, converged);

        System.arraycopy(warp[0], 0, pInit[0], 0, pInit[0].length);
        System.arraycopy(warp[1], 0, pInit[1], 0, pInit[1].length);

        return new double[]{nPAndErr[1], nIter};
    }

    /**
     * 2D alignment method for each keypoint in windows of half widths hX and hY.
     * NOTE that internally, this method copies the window sections into images of those dimensions to calculate the
     * warps.  For 2D translation, the alternative method inverseComposition2DTranslationKeypoints() produces
     * slightly better results and has better performance.
     *
     * @param template
     * @param image
     * @param xYInit
     * @param maxIter
     * @param eps e.g. 0.1 or 0.5
     * @param yXKeypoints a 2-D array of size [nKeypoints x 2] of keypoints to calculate affine warp at.
     *                  Note that the keypoints should not be closer to the image bounds than patchHalfWidth.
     *                    The keypoints are given as pairs of {y, x} coordinates.
     * @param hX half width of window around keypoint in x dimension.
     * @param hY half width of window around keypoint in y dimension.
     *
     * @return
     * @throws NotConvergedException
     * @throws IOException
     */
    public static Warps inverseComposition2DTranslationKeypointsCpImgs(double[][] template, double[][] image,
                                                                       double[] xYInit, int maxIter, double eps, int[][] yXKeypoints,
                                                                       int hX, int hY, Type type) throws NotConvergedException, IOException {
        if (xYInit.length != 2) {
            throw new IllegalArgumentException("expecting xYInit to be length 2");
        }

        double[][] pInit = new double[][] {
                {1, 0, xYInit[0]},
                {0, 1, xYInit[1]}
        };
        return inverseCompositionKeypointsCpImgs(template, image, pInit, maxIter, eps, yXKeypoints, hX, hY, type);
    }

    /**
     * 2D alignment method for each keypoint in windows of half widths hX and hY.
     * NOTE that internally, the method operates on windows in the full images rather than copying the windows
     * to smaller images as is done in inverseCompositionKeypointsCpImgs.
     * For 2D translation warps, this method provides slightly better results.
     * @param template
     * @param image
     * @param xYInit
     * @param maxIter
     * @param eps e.g. 0.1 or 0.5
     * @param yXKeypoints a 2-D array of size [nKeypoints x 2] of keypoints to calculate affine warp at.
     *                  Note that the keypoints should not be closer to the image bounds than patchHalfWidth.
     * @param hX half width of window around keypoint in x dimension.
     * @param hY half width of window around keypoint in y dimension.
     *
     * @return
     * @throws NotConvergedException
     * @throws IOException
     */
    public static Warps inverseComposition2DTranslationKeypoints(double[][] template, double[][] image,
                                                                 double[] xYInit, int maxIter, double eps, int[][] yXKeypoints,
                                                                 int hX, int hY, Type type) throws NotConvergedException, IOException {
        if (xYInit.length != 2) {
            throw new IllegalArgumentException("expecting xYInit to be length 2");
        }

        double[][] pInit = new double[][] {
                {1, 0, xYInit[0]},
                {0, 1, xYInit[1]}
        };
        return inverseCompositionKeypoints(template, image, pInit, maxIter, eps, yXKeypoints, hX, hY, type);
    }

    /**
     2D alignment method for each keypoint in windows of half widths hX and hY.
     * NOTE that internally, this method copies the window sections into images of those dimensions to calculate the
     * warps.
     * For 2D translation, the alternative method inverseComposition2DTranslationKeypoints() produces
     * slightly better results and has better performance.
     * For 2D affine warps, one should prefer this method.
     *
     * @param template
     * @param image
     @param pInit the initial projection matrix of size [2 rows X 3 columns].
     It should be composed of { {1 + affine[0][0], affine[0][1], transX},
     {affine[1][0], 1 + affine[1][1], transY}}.
     e.g. {{1,0,0},{0,1,0}}
     * @param maxIter
     * @param eps
     * @param yXKeypoints a 2-D array of size [nKeypoints x 2] of keypoints to calculate affine warp at.
     *                  Note that the keypoints should not be closer to the image bounds than patchHalfWidth.
     *                    need at least 3 points for the 2D affine.
     * @param hX half width of window around keypoint in x dimension.
     * @param hY half width of window around keypoint in y dimension.
     * @return
     * @throws NotConvergedException
     * @throws IOException
     */
    public static Warps inverseCompositionKeypointsCpImgs(double[][] template, double[][] image,
                                                          double[][] pInit, int maxIter, double eps, int[][] yXKeypoints, int hX, int hY,
                                                          Type type) throws NotConvergedException, IOException {

        if (pInit.length != 2 || pInit[0].length != 3) {
            throw new IllegalArgumentException("expecting pInit length to be 2, and pInit[0].length to be 3");
        }
        if (type.equals(Type.AFFINE_2D)) {
            if (yXKeypoints.length * (2 * Math.min(hX, hY) +1) < 3) {
                throw new IllegalArgumentException("need at least 3 points to solve for affine");
            }
        }

        int nKeypoints = yXKeypoints.length;

        Warps warps = new Warps();
        warps.warps = new double[nKeypoints][][];
        double ssd = 0;
        int y, x;
        int maxOfNIters = 0;
        for (int i = 0; i < yXKeypoints.length; ++i) {
            y = yXKeypoints[i][0];
            x = yXKeypoints[i][1];

            // exclude out of bounds
            if (y - hY < 0 || x - hX < 0 || y+hY>= template.length || x+hX >= template[0].length) {
                throw new IllegalArgumentException("keypoints and their windows need to be within bounds of template image");
            }

            double[][] tImg = MatrixUtil.copySubMatrix(template, y - hY, y + hY, x - hX, x + hX);
            double[][] iImg = MatrixUtil.copySubMatrix(image, y - hY, y + hY, x - hX, x + hX);

            double[][] _pInit = MatrixUtil.copy(pInit);
            double[] errSSD = Alignment.inverseCompositional(tImg, iImg, _pInit, maxIter, type, eps);
            maxOfNIters = (int)Math.max(maxOfNIters, Math.round(errSSD[1]));

            warps.warps[i] = MatrixUtil.copy(_pInit);

            // TODO: revisit this. might need to add in quadrature
            ssd += errSSD[0];
        }

        warps.ssd = ssd;
        warps.nIterMax = maxOfNIters;
        return warps;
    }

    /**
     2D AFFINE alignment method for each keypoint in windows of half widths hX and hY.
     * NOTE that internally, this method copies the window sections into images of those dimensions to calculate the
     * warps.  Also, not that the warp is prefixed by an image-wide optical flow calculation
     * in order to initialize the projection matrix.
     *
     * @param template
     * @param image
     * @param maxIter
     * @param eps
     * @param yXKeypoints a 2-D array of size [nKeypoints x 2] of keypoints to calculate affine warp at.
     *                  Note that the keypoints should not be closer to the image bounds than patchHalfWidth.
     *                    need at least 3 points for the 2D affine.
     * @param hX half width of window around keypoint in x dimension.
     * @param hY half width of window around keypoint in y dimension.
     * @return
     * @throws NotConvergedException
     * @throws IOException
     */
     static Warps inverseCompositionKeypointsCpImgs2DAffinePre2DTrans(double[][] template, double[][] image,
                                                          int maxIter, double eps, int[][] yXKeypoints, int hX, int hY)
             throws NotConvergedException, IOException {

        if (yXKeypoints.length * (2 * Math.min(hX, hY) +1) < 3) {
            throw new IllegalArgumentException("need at least 3 points to solve for affine");
        }

        int nKeypoints = yXKeypoints.length;

        double[] xYInit = new double[2];
        double[] errSSD = Alignment.inverseCompositional2DTranslation(template, image, xYInit, maxIter, eps);

        Warps warps = new Warps();
        warps.warps = new double[nKeypoints][][];
        double ssd = 0;
        int y, x;
        int maxOfNIters = 0;
        for (int i = 0; i < yXKeypoints.length; ++i) {
            y = yXKeypoints[i][0];
            x = yXKeypoints[i][1];

            // exclude out of bounds
            if (y - hY < 0 || x - hX < 0 || y+hY>= template.length || x+hX >= template[0].length) {
                throw new IllegalArgumentException("keypoints and their windows need to be within bounds of template image");
            }

            double[][] tImg = MatrixUtil.copySubMatrix(template, y - hY, y + hY, x - hX, x + hX);
            double[][] iImg = MatrixUtil.copySubMatrix(image, y - hY, y + hY, x - hX, x + hX);

            double[][] _pInit = new double[][]{
                    {1, 0, xYInit[0]},
                    {0, 1, xYInit[1]}
            };
            errSSD = Alignment.inverseCompositional(tImg, iImg, _pInit, maxIter,
                    Type.AFFINE_2D, eps);
            maxOfNIters = (int)Math.max(maxOfNIters, Math.round(errSSD[1]));

            warps.warps[i] = MatrixUtil.copy(_pInit);

            // TODO: revisit this. might need to add in quadrature
            ssd += errSSD[0];
        }

        warps.ssd = ssd;
        warps.nIterMax = maxOfNIters;
        return warps;
    }

    /**
     * 2D alignment for each keypoint in windows of half widths hX and hY.
     * NOTE that internally, the method operates on windows in the full images rather than copying the windows
     * to smaller images as is done in inverseCompositionKeypointsCpImgs.

     For 2D translation, one would prefer this method to the alternative method inverseCompositionKeypointsCpImgs(...).
     For 2D affine warps, one should prefer inverseCompositionKeypointsCpImgs(...) to this method.
     *
     * @param template
     * @param image
     @param pInit the initial projection matrix of size [2 rows X 3 columns].
     It should be composed of { {1 + affine[0][0], affine[0][1], transX},
     {affine[1][0], 1 + affine[1][1], transY}}.
     e.g. {{1,0,0},{0,1,0}}
     * @param maxIter
     * @param eps
     * @param yXKeypoints a 2-D array of size [nKeypoints x 2] of keypoints to calculate affine warp at.
     *                  Note that the keypoints should not be closer to the image bounds than patchHalfWidth.
     *                    need at least 3 points for the 2D affine.
     * @param hX half width of window around keypoint in x dimension.
     * @param hY half width of window around keypoint in y dimension.
     * @return
     * @throws NotConvergedException
     * @throws IOException
     */
    public static Warps inverseCompositionKeypoints(double[][] template, double[][] image,
                                                    double[][] pInit, int maxIter, double eps, int[][] yXKeypoints, int hX, int hY,
                                                    Type type) throws NotConvergedException, IOException {

        if (pInit.length != 2 || pInit[0].length != 3) {
            throw new IllegalArgumentException("expecting pInit length to be 2, and pInit[0].length to be 3");
        }
        if (type.equals(Type.AFFINE_2D)) {
            if (yXKeypoints.length * (2 * Math.min(hX, hY) +1) < 3) {
                throw new IllegalArgumentException("need at least 3 points to solve for affine");
            }
        }

        int nKeypoints = yXKeypoints.length;

        Warps warps = new Warps();
        warps.warps = new double[nKeypoints][][];
        double ssd = 0;
        int y, x;
        int maxOfNIters = 0;
        for (int i = 0; i < yXKeypoints.length; ++i) {
            y = yXKeypoints[i][0];
            x = yXKeypoints[i][1];

            double[][] _pInit = MatrixUtil.copy(pInit);
            double[] errSSD = Alignment.inverseCompositional(template, image, x, y, hX, hY, _pInit, maxIter, type, eps);

            maxOfNIters = (int)Math.max(maxOfNIters, Math.round(errSSD[1]));

            warps.warps[i] = MatrixUtil.copy(_pInit);

            // TODO: revisit this. might need to add in quadrature
            ssd += errSSD[0];
        }

        warps.ssd = ssd;
        warps.nIterMax = maxOfNIters;
        return warps;
    }


    static Warps inverseCompositionKeypoints2DAffinePre2DTrans(double[][] template, double[][] image, int maxIter,
                                                               double eps, int[][] yXKeypoints, int hX, int hY) throws IOException, NotConvergedException {

        if (yXKeypoints.length * (2 * Math.min(hX, hY) +1) < 3) {
            throw new IllegalArgumentException("need at least 3 points to solve for affine");
        }
        double[] xYInit = new double[2];
        double[] errSSD = Alignment.inverseCompositional2DTranslation(template, image, xYInit, maxIter, eps);

        int nKeypoints = yXKeypoints.length;

        Warps warps = new Warps();
        warps.warps = new double[nKeypoints][][];
        double ssd = 0;
        int y, x;
        int maxOfNIters = 0;
        for (int i = 0; i < yXKeypoints.length; ++i) {
            y = yXKeypoints[i][0];
            x = yXKeypoints[i][1];

            double[][] _pInit = new double[][]{
                    {1, 0, xYInit[0]},
                    {0, 1, xYInit[1]}
            };
            errSSD = Alignment.inverseCompositional(template, image, x, y, hX, hY, _pInit, maxIter,
                    Type.AFFINE_2D, eps);

            maxOfNIters = (int)Math.max(maxOfNIters, Math.round(errSSD[1]));

            warps.warps[i] = MatrixUtil.copy(_pInit);

            // TODO: revisit this. might need to add in quadrature
            ssd += errSSD[0];
        }

        warps.ssd = ssd;
        warps.nIterMax = maxOfNIters;
        return warps;
    }


    protected static class Result {
        double[][] pSums;
        double ssd;
    }

    protected static class Warps {
        public int nIterMax;
        double[][][] warps;
        double ssd;
    }

    private static double[][] createSteepestDescentImagesAffine2D(double[][] gTX, double[][] gTY) {
        return createSteepestDescentImagesAffine2D(gTX, gTY, 0, gTX[0].length-1, 0, gTX.length - 1);
    }

    /**
     *
     * @param template
     * @param beginX first x coordinate
     * @param endX last x coordinate, inclusive
     * @param beginY first y coordinate
     * @param endY last y coordinate, inclusive
     * @return
     */
    private static double[][] createSteepestDescentImagesAffine2D(double[][] template,
                                                                  int beginX, int endX, int beginY, int endY) {
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

        int n = (endX - beginX + 1) * (endY - beginY + 1);

        double xc = (endX - beginX)/2;
        double yc = (endY - beginY)/2;

        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();

        double[][] dTdWdp = new double[n][];
        double[][] tmp1 = new double[1][2];
        double[][] tmp2 = new double[2][6];
        double gX, gY;
        for (int y = beginY, idx = 0; y <= endY; ++y) {
            for (int x = beginX; x <= endX; ++x, ++idx) {

                // method expects column major data, so swap row and col:
                gX = kernel1DHelper.convolvePointWithKernel(template, y, x,
                        beginY, endY, beginX, endX, kernel, false);

                gY = kernel1DHelper.convolvePointWithKernel(template, y, x,
                        beginY, endY, beginX, endX, kernel, true);

                tmp1[0][0] = gX;
                tmp1[0][1] = gY;

                tmp2[0][0] = x - xc;
                tmp2[0][2] = y - yc;
                tmp2[0][4] = 1;
                tmp2[1][1] = x - xc;
                tmp2[1][3] = y - yc;
                tmp2[1][5] = 1;

                //[1X6]
                dTdWdp[idx] = MatrixUtil.multiply(tmp1, tmp2)[0];
            }
        }
        return dTdWdp;
    }
    /**
     *
     * @param gTX
     * @param gTY
     * @param beginX first x coordinate
     * @param endX last x coordinate, inclusive
     * @param beginY first y coordinate
     * @param endY last y coordinate, inclusive
     * @return
     */
    private static double[][] createSteepestDescentImagesAffine2D(double[][] gTX, double[][] gTY,
                                                                  int beginX, int endX, int beginY, int endY) {
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

        int n = (endX - beginX + 1) * (endY - beginY + 1);

        double xc = (endX - beginX)/2;
        double yc = (endY - beginY)/2;

        double[][] dTdWdp = new double[n][];
        double[][] tmp1 = new double[1][2];
        double[][] tmp2 = new double[2][6];
        double gX, gY;
        for (int y = beginY, idx = 0; y <= endY; ++y) {
            for (int x = beginX; x <= endX; ++x, ++idx) {

                tmp1[0][0] = gTX[y][x];
                tmp1[0][1] = gTY[y][x];

                tmp2[0][0] = x - xc;
                tmp2[0][2] = y - yc;
                tmp2[0][4] = 1;
                tmp2[1][1] = x - xc;
                tmp2[1][3] = y - yc;
                tmp2[1][5] = 1;

                //[1X6]
                dTdWdp[idx] = MatrixUtil.multiply(tmp1, tmp2)[0];
            }
        }
        return dTdWdp;
    }

    private static double[][] createSteepestDescentImageTranslation2D(double[][] gTX, double[][] gTY) {
        return createSteepestDescentImageTranslation2D(gTX, gTY, 0, gTX[0].length-1, 0, gTX.length - 1);
    }

    private static double[][] createSteepestDescentImageTranslation2D(double[][] gTX, double[][] gTY,
                                                                      int beginX, int endX, int beginY, int endY) {
        //[nTR*nTC X 2]
        // steepest descent images  gT * dWdP = [nTR X nTC] * [nTR*nTC X 6] = [nTR*nTC X 6] for x then for y
        // can format as [nTr X nTc X 6] or [nTr * nTc X 6].
        // will choose the latter to make the Hessian mult easier
        //steepest descent img = gradient * dWdP
        //                       [1x2] * [2x2]
        int n = (endX - beginX + 1) * (endY - beginY + 1);

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

        double[][] dTdWdp = new double[n][];
        for (int y = beginY, idx = 0; y <= endY; ++y) {
            for (int x = beginX; x <= endX; ++x, ++idx) {
                dTdWdp[idx] = new double[]{gTX[y][x], gTY[y][x]};
            }
        }
        return dTdWdp;
    }

    private static double[][] createSteepestDescentImageTranslation2D(double[][] template,
                                                                      int beginX, int endX, int beginY, int endY) {
        //[nTR*nTC X 2]
        // steepest descent images  gT * dWdP = [nTR X nTC] * [nTR*nTC X 6] = [nTR*nTC X 6] for x then for y
        // can format as [nTr X nTc X 6] or [nTr * nTc X 6].
        // will choose the latter to make the Hessian mult easier
        //steepest descent img = gradient * dWdP
        //                       [1x2] * [2x2]
        int n = (endX - beginX + 1) * (endY - beginY + 1);

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

        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();

        double[][] dTdWdp = new double[n][];
        for (int y = beginY, idx = 0; y <= endY; ++y) {
            for (int x = beginX; x <= endX; ++x, ++idx) {
                dTdWdp[idx] = new double[]{
                        kernel1DHelper.convolvePointWithKernel(template, y, x,
                                beginY, endY, beginX, endX, kernel, false),
                        kernel1DHelper.convolvePointWithKernel(template, y, x,
                                beginY, endY, beginX, endX, kernel, true)
                };
            }
        }
        return dTdWdp;
    }

    private static double[][] init2DTranslationWarp(double[][] pInit) {
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

        return sumSteepestDescErrImageProduct(image, t, warp, steep, outPSum, type, 0, t[0].length-1, 0, t.length - 1);
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
        double[][] warp, double[][] steep, double[] outPSum, Type type,
                                                           int beginX, int endX, int beginY, int endY) {

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
        for (int y = beginY, idxT =0; y <= endY; ++y) {
            for (int x = beginX; x <= endX; ++x, ++idxT) {
                xy[0] = x;
                xy[1] = y;
                xy2 = MatrixUtil.multiplyMatrixByColumnVector(warp, xy);

                // roundoff
                xy2[0] = Math.round(xy2[0] * 1E2)/1E2;
                xy2[1] = Math.round(xy2[1] * 1E2)/1E2;

                if (xy2[0] < beginX || Math.ceil(xy2[0]) > endX || xy2[1] < beginY || Math.ceil(xy2[1]) > endY) {
                    continue;
                }
                // method expecting col major data so reverse the coords:
                v2 = imageProcessor.biLinearInterpolation(image, xy2[1], xy2[0]);

                // error image is I(w(x;p)) - T(x).  not squared nor abs value of
                // TODO: this could be improved by considering the fraction of the pixel represented
                //if (type.equals(Type.AFFINE_2D)) {
                //    diff = v2 - t[y][x]; // use if change from inverse warp(x; deltaP)  to forward warp(x;deltap) in update
                //} else if (type.equals(Type.TRANSLATION_2D)) {
                    diff = t[y][x] - v2;
                //}

                errSSD += diff * diff;

                pS = steep[idxT];

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
        if (nP == 0) {
            return new double[]{nP, Double.MAX_VALUE};
        }
        return new double[]{nP, errSSD};
    }

    private static double[][] inverseWParam2DTrans(double[] deltaP) throws NotConvergedException {
        return new double[][]{
                {1, 0, -deltaP[0]},
                {0, 1, -deltaP[1]},
                {0, 0, 1}
        };
    }

    private static void update2DTransICWarp(double[][] warp, double[] deltaP) throws NotConvergedException {
        double[][] params = inverseWParam2DTrans(deltaP); //[3X3]

        double[][] tmp = MatrixUtil.multiply(warp, params);
        for (int i = 0; i < warp.length; ++i) {
            System.arraycopy(tmp[i], 0, warp[i], 0, tmp[i].length);
        }
    }

    private static void update2DAffFCWarp(double[][] warp, double[] deltaP)  {

        // for affine, Bouguet uses forward composition w = w * w(delta)
        // instead of composition w = w * w(delta)^-1 from Baker-Matthews.
        // it gives similarly wrong results for test data, and similarly
        // rarely right results (possibly are same to first order?).
        // will use the Baker-Matthew's paper's inverse compositional update
        // (future proofing for when have improved gradient windows here, etc.)

        double[][] params = new double[][]{
                {1+deltaP[0], deltaP[2], deltaP[4]},
                {deltaP[1], 1+deltaP[3], deltaP[5]},
                {0, 0, 1}
        };
        double[][] invParams = MatrixUtil.inverse(params);
        /*if (Double.isNaN(invParams[0][0])) {
            try {
                invParams = MatrixUtil.pseudoinverseRankDeficient(params);
            } catch (NotConvergedException e) {
                invParams = new double[][]{
                        {0,0,0},
                        {0,0,0},
                        {0, 0, 0}
                };
            }
        }*/
        params = invParams;

        double[][] tmp = MatrixUtil.multiply(warp, params);

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

    protected static double[][] strictColumnDiff(double[][] a) {
        double[][] b = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; ++i) {
            for (int j = 1; j < b.length; ++j) {
                b[i][j] = a[i][j] - a[i][j-1];
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
}
