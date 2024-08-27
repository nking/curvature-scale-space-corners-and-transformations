package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.*;
import algorithms.imageProcessing.features.PhaseCongruencyDetector;
import algorithms.imageProcessing.features.ShiTomasi;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.misc.MiscMath0;
import algorithms.sort.MiscSorter;
import algorithms.statistics.Standardization;
import algorithms.util.FormatArray;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

import java.io.IOException;
import java.util.*;

/**
 *  TODO: add test data generation and tests for affine projected data (rotation and shear for a given
 *   pyramid scale image)
 *
 *   TODO: add test data generation and tests for multiple moving objects (specified by keypoints and found by keypoints).
 *
 * @author nichole
 */
public class AlignmentTest extends TestCase {

    boolean printDebug = false;

    public AlignmentTest() {
    }

    public void test2InverseCompositionAlignment() throws IOException, NotConvergedException {
        if (printDebug)
            System.out.printf("test2InverseCompositionAlignment\n");
        int pixMax = 32*2;
        int m = 32*2;
        //System.out.printf("pixMax=%d, m=%d\n", pixMax, m);
        int n = m;
        double[][] im1;
        double[][] im2;
        int maxIter = 100;//1000;
        double eps = 0.1;

        for (int iTest = 1; iTest < (m/2)-3; ++iTest) {

            double dX = iTest;
            double dY = iTest;

            Images images = generateImages(m,  n, dX, dY, pixMax);
            if (images == null || images.yXKeypoints.length == 0) {
                break;
            }
            im1 = images.im1;
            im2 = images.im2;

            // set all to 0 to solve for
            //MatrixUtil.fill(pInit, 0.);

            if (printDebug)
            System.out.printf("\niTest=%d\n", iTest);

            double[] errSSD;
            ///*  TEST 2-D translation
            double[] xYInit = new double[]{0, 0};
            errSSD = Alignment.inverseCompositional2DTranslation(im1, im2, xYInit, maxIter, eps);
            if (printDebug)
                System.out.printf("2D trans newton: nIter=%d, uvFinal=\n%s\n",
                    (int)Math.round(errSSD[1]), FormatArray.toString(xYInit, "%.2f"));
            assertEquals(iTest, (int)Math.round(xYInit[0]));
            assertEquals(iTest, (int)Math.round(xYInit[1]));
            //xYInit = new double[]{uInit, vInit};
            //errSSD = Alignment._inverseCompositional2DTranslationGN(im1, im2, xYInit, maxIter, eps);
            //System.out.printf("2D trans GN: uvFinal=\n%s\n", FormatArray.toString(xYInit, "%.2f"));
            //assertEquals(iTest, (int)Math.round(xYInit[0]));
            //assertEquals(iTest, (int)Math.round(xYInit[1]));
            //*/
            ///*

            double[][] pInit = new double[][]{
                    {1, 0, 0}, {0, 1, 0}
            };
            errSSD = Alignment.inverseCompositional2DAffine(im1, im2, pInit, maxIter, eps);
            if (printDebug)
                System.out.printf("2D affine Newton: nIter=%d, pFinal=\n%s\n", (int)Math.round(errSSD[1]),
                    FormatArray.toString(pInit, "%.2f"));

            //*/
        }

    }

    static class Images {
        double[][] im1;
        double[][] im2;
        double dX;
        double dY;
        int[][] yXKeypoints;
        int hX, hY;
    }
    private Images generateImages(int m, int n, double dX, double dY, double pixMax) {
        double[][] im1 = new double[m][n];
        double[][] im2 = new double[m][n];
        double deltaI = pixMax / m;

        int iP = 0;

        int hX = (int)(dX);
        int hY = (int)(dY);
        if (hX < 2) {
            hX = 2;
        }
        if (hY < 2) {
            hY = 2;
        }

        int nPatches = m*n;
        if (nPatches < 1) return null;
        int[][] yXPatches = new int[nPatches][];

        double[] i2j2 = new double[3];

        double[][] pInit = MatrixUtil.createIdentityMatrix(3);
        pInit[0][2] = dX;
        pInit[1][2] = dY;

        for (int i = 0; i < m; ++i) {//x
            double b = deltaI * i;
            b = (int) Math.round(b); // to match image processing temporarily
            for (int j = 0; j <= i; ++j) {//y
                im1[i][j] = b;
                im1[j][i] = b;

                // warp coords for im2
                MatrixUtil.multiplyMatrixByColumnVector(pInit, new double[]{i, j, 1}, i2j2);
                int i3 = (int)Math.round(i2j2[0]);
                int j3 = (int)Math.round(i2j2[1]);
                boolean b0 = (i3 >= 0 && i3 < m && j3 >= 0 && j3 < n);
                if (b0) {
                    im2[i3][j3] = b;
                }

                MatrixUtil.multiplyMatrixByColumnVector(pInit, new double[]{j, i, 1}, i2j2);
                int i4 = (int)Math.round(i2j2[0]);
                int j4 = (int)Math.round(i2j2[1]);
                boolean b1 = (i4 >= 0 && i4 < m && j4 >= 0 && j4 < n);
                if (b1) {
                    im2[i4][j4] = b;
                }
                if (!b1 || !b0) continue;

                if (i==j) {
                    if (i - hY < 0 || i + dY + hX >= m) continue;
                    if (j - hX < 0 || j + dX + hX >= n) continue;
                    if (i3 - hY < 0 || i3 + hY >= m) continue;
                    if (j3 - hX < 0 || j3 + hX >= n) continue;
                    yXPatches[iP++] = new int[]{i, i};
                }
            }
        }
        yXPatches = Arrays.copyOf(yXPatches, iP);

        // apply same min max normalization to all images, so learn it from the template:
        double[] minMax = Standardization.minMaxNormalizeImage(im1);
        Standardization.minMaxNormalizeImage(im2, minMax[0], minMax[1]);

        Images images = new Images();
        images.im1 = im1;
        images.im2 = im2;
        images.dX = dX;
        images.dY = dY;
        images.hX = hX;
        images.hY = hY;
        images.yXKeypoints = yXPatches;
        return images;
    }

    private Images generateImagesV(int m, int n, double dX, double dY, double pixMax) {
        // draws a "V" with internal angle of 90 degrees
        // (so gradients gX and gY in keypoint windows do not have 0s)
        double[][] im1 = new double[m][n];
        double[][] im2 = new double[m][n];
        double deltaI = pixMax / (m+1);

        int iP = 0;

        int hX = (int)(dX);
        int hY = (int)(dY);
        if (hX < 2) {
            hX = 2;
        }
        if (hY < 2) {
            hY = 2;
        }

        int nPatches = m*n;
        if (nPatches < 1) return null;
        int[][] yXPatches = new int[nPatches][];

        double[] x2y2 = new double[3];

        double[][] pInit = MatrixUtil.createIdentityMatrix(3);
        pInit[0][2] = dX;
        pInit[1][2] = dY;

        final int xKP = n/2;
        for (int yKP = 0; yKP < m; ++yKP) {//x
            double b = deltaI * (yKP+1);

            boolean bKP = false;

            // RHS of right angle V:
            for (int x = xKP, y = yKP; x < n && y < m; ++x, ++y) {
                im1[y][x] = b;

                MatrixUtil.multiplyMatrixByColumnVector(pInit, new double[]{x, y, 1}, x2y2);
                int x2 = (int)Math.round(x2y2[0]);
                int y2 = (int)Math.round(x2y2[1]);
                boolean b1 = (x2 >= 0 && x2 < n && y2 >= 0 && y2 < m);
                if (b1) {
                    im2[y2][x2] = b;

                    // if both windows are within image bounds, can write to keypoint list:
                    if (x == xKP && y==yKP && x-hX >= 0 && x2-hX >= 0 && x+hX < n && x2+hX<n
                    && y-hY >= 0 && y2-hY >= 0 && y+hY < m && y2+hY<m) {
                        bKP = true;
                    }
                }
            }
            // LHS of right angle V:
            for (int x = xKP-1, y=yKP+1; x >= 0 && y < m; --x, ++y) {
                im1[y][x] = b;
                MatrixUtil.multiplyMatrixByColumnVector(pInit, new double[]{x, y, 1}, x2y2);
                int x2 = (int)Math.round(x2y2[0]);
                int y2 = (int)Math.round(x2y2[1]);
                if (x2 >= 0 && x2 < n && y2 >= 0 && y2 < m) {
                    im2[y2][x2] = b;
                }
            }
            if (bKP) {
                yXPatches[iP++] = new int[]{yKP, xKP};
            }
        }
        yXPatches = Arrays.copyOf(yXPatches, iP);

        // apply same min max normalization to all images, learned from the template:
        double[] minMax = Standardization.minMaxNormalizeImage(im1);
        Standardization.minMaxNormalizeImage(im2, minMax[0], minMax[1]);

        Images images = new Images();
        images.im1 = im1;
        images.im2 = im2;
        images.dX = dX;
        images.dY = dY;
        images.hX = hX;
        images.hY = hY;
        images.yXKeypoints = yXPatches;
        return images;
    }

    private static void fillWithRandom(double[][] im, int minVal, int maxVal, Random rand) {
        for (int i = 0; i < im.length; ++i) {//x
            for (int j = 0; j < im[i].length; ++j) {//x
                im[i][j] = minVal + rand.nextInt(maxVal - minVal);
            }
        }

    }

    private Images generateImagesRandom(int m, int n, double dX, double dY, double pixMax,
                                        int type /* 0 = semi-detailed blobs, 1=big blobs */) {
        double[][] im1 = new double[m][n];
        double[][] im2 = new double[m][n];

        //the first image will be a random pattern, the 2nd will be im1 warped by projection matrix

        long seed = System.nanoTime();
        System.out.println("seed=" + seed);
        Random rand = new Random(seed);

        fillWithRandom(im1, 0, 255, rand);

        if (type == 0) {
            binaryThreshold(im1);
            blur(im1, SIGMA.ONE);
            dilateAboveBinaryThreshold(im1);
        } else {
            binaryThreshold(im1);
            blur(im1, SIGMA.THREE);
        }

        //fillWithRandom(im2, 0, 255, rand);

        int iP = 0;

        int hX = (int)(dX);
        int hY = (int)(dY);
        if (hX < 2) {
            hX = 2;
        }
        if (hY < 2) {
            hY = 2;
        }

        int nPatches = m*n;
        if (nPatches < 1) return null;
        int[][] yXPatches = new int[nPatches][];

        double[] i2j2 = new double[3];

        double[][] pInit = MatrixUtil.createIdentityMatrix(3);
        pInit[0][2] = dX;
        pInit[1][2] = dY;

        for (int i = 0; i < m; ++i) {//y
            for (int j = 0; j < n; ++j) {//x

                // warp coords for im2
                MatrixUtil.multiplyMatrixByColumnVector(pInit, new double[]{i, j, 1}, i2j2);
                int i3 = (int)Math.round(i2j2[0]);
                int j3 = (int)Math.round(i2j2[1]);
                if (i3 < 0 || i3 >= m || j3 < 0 || j3 >= n) {
                    continue;
                }
                im2[i3][j3] = im1[i][j];

                if (i == j) {
                    if (i - hY < 0 || i + dY + hX >= m) continue;
                    if (j - hX < 0 || j + dX + hX >= n) continue;
                    if (i3 - hY < 0 || i3 + hY >= m) continue;
                    if (j3 - hX < 0 || j3 + hX >= n) continue;
                    yXPatches[iP++] = new int[]{i, i};
                }
            }
        }
        yXPatches = Arrays.copyOf(yXPatches, iP);

        // apply same min max normalization to all images, so learn it from the template:
        double[] minMax = Standardization.minMaxNormalizeImage(im1);
        Standardization.minMaxNormalizeImage(im2, minMax[0], minMax[1]);

        Images images = new Images();
        images.im1 = im1;
        images.im2 = im2;
        images.dX = dX;
        images.dY = dY;
        images.hX = hX;
        images.hY = hY;
        images.yXKeypoints = yXPatches;
        return images;
    }

    public void test3InverseCompositionAlignment() throws IOException, NotConvergedException {
        if (printDebug)
            System.out.printf("test3InverseCompositionAlignment\n");

        int pixMax = 32;//32*2;
        int m = 32;//32*2;
        //System.out.printf("pixMax=%d, m=%d\n", pixMax, m);
        int n = m;
        double[][] im1;
        double[][] im2;
        int maxIter = 50;//100;
        double eps = 0.01;

        double perturb = 1E-3;

        int nPoints = m*n;

        //for (int gImg = 0; gImg < 4; ++gImg) {
        for (int iTest = 0; iTest < (m/3); ++iTest) {

            String gTitle = null;
            Images images = null;

            //for (int iTest = 1; iTest < (m/3); ++iTest) {
            for (int gImg = 0; gImg < 4; ++gImg) {

                double dX = iTest;
                double dY = iTest;

                switch (gImg) {
                    case(0) : {
                        gTitle = "bands of LLRect";
                        images = generateImages(m, n, dX, dY, pixMax);
                        break;
                    }
                    case(1) : {
                        gTitle = "V pattern";
                        images = generateImagesV(m, n, dX, dY, pixMax);
                        break;
                    }
                    case(2) : {
                        gTitle = "Random blobs with details";
                        images = generateImagesRandom(m, n, dX, dY, pixMax, 0);
                        break;
                    }
                    case(3) : {
                        gTitle = "Random big blobs";
                        images = generateImagesRandom(m, n, dX, dY, pixMax, 1);
                        break;
                    }
                }

                if (images == null) {
                    continue;
                }
                im1 = images.im1;
                im2 = images.im2;
                int hX = images.hX;
                int hY = images.hY;

                int[][] yXKeypoints0 = images.yXKeypoints;

                double[] outMinMax = new double[2];
                //shitomasi doesn't pick up the lower left corner of each rectangle band in generateImages()
                int[][] keypoints3 = ShiTomasi.goodFeatureCoordinates(im1, nPoints, outMinMax);
                int[][] yXKeypoints3 = new int[keypoints3.length][2];
                int iP = 0;
                for (int i = 0; i < keypoints3.length; ++i) {
                    yXKeypoints3[iP][0] = keypoints3[i][1];
                    yXKeypoints3[iP][1] = keypoints3[i][0];

                    if (yXKeypoints3[iP][0] - hY < 0 || yXKeypoints3[iP][1] - hX < 0
                            || yXKeypoints3[iP][0] + hY >= im1.length || yXKeypoints3[iP][1] + hX >= im1[0].length) {
                        continue;
                    }
                    ++iP;
                }
                yXKeypoints3 = Arrays.copyOf(yXKeypoints3, iP);

                GreyscaleImage gsImg = makeGSImg(im1);

                // for 2D translation, the whole image result is optical flow and is always correct with
                //   the 3 types of datasets.
                // for the other 2D translation methods, the manually created keypoints0 have better results
                //   than the orb generated keypoints.
                ORB orb1 = new ORB(gsImg, nPoints);
                orb1.overrideToNotCreateDescriptors();
                orb1.overrideToUseSingleScale();
                //orb1.overrideToCreateCurvaturePoints();
                orb1.overrideToAlsoCreate1stDerivKeypoints();
                orb1.detectAndExtract();
                List<PairInt> kp1 = orb1.getAllKeyPoints();

                int[][] yXKeypoints = new int[kp1.size()][2];
                iP = 0;
                for (int i = 0; i < kp1.size(); ++i) {
                    yXKeypoints[iP][0] = kp1.get(i).getY();
                    yXKeypoints[iP][1] = kp1.get(i).getX();
                    if (yXKeypoints[iP][0] - hY < 0 || yXKeypoints[iP][1] - hX < 0
                            || yXKeypoints[iP][0] + hY >= im1.length || yXKeypoints[iP][1] + hX >= im1[0].length) {
                        continue;
                    }
                    ++iP;
                }
                yXKeypoints = Arrays.copyOf(yXKeypoints, iP);

                /*
                Image tmp = gsImg.copyToColorGreyscale();
                for (int i = 0; i < yXKeypoints.length; ++i) {
                    ImageIOHelper.addPointToImage(yXKeypoints[i][1], yXKeypoints[i][0], tmp, 2, 255, 0, 0);
                }
                long ts = MiscDebug.getCurrentTimeFormatted();
                MiscDebug.writeImage(tmp, ts + "_orb_");
                MiscDebug.writeImage(makeGSImg(im1).copyToColorGreyscale(), ts + "_1");
                MiscDebug.writeImage(makeGSImg(im2).copyToColorGreyscale(), ts + "_2");
                */

                if (printDebug) {
                    System.out.println("-----------------------------------------");
                    System.out.printf("\niTest=%d  (nPoints=%d)\ngenerated image type=%s\n",
                            iTest, yXKeypoints.length, gTitle);
                }

                //System.out.printf("im1=\n%s", FormatArray.toString(im1, "%.2f"));
                //System.out.printf("im2=\n%s", FormatArray.toString(im2, "%.2f"));

                int chk = 0;//yXKeypoints.length / 2;

                Alignment.Warps warps;
                double[] errSSD;
                double[] xYInit;

                ///*  TEST 2-D translation
                xYInit = new double[]{perturb, perturb};

                errSSD = Alignment.inverseCompositional2DTranslation(im1, im2, xYInit, maxIter, eps);
                if (printDebug)
                    System.out.printf("2D trans (whole image):\n  nIter=%d\n  %s\n",
                        (int) Math.round(errSSD[1]), FormatArray.toString(xYInit, "%.2f"));
                //*/
                xYInit = new double[]{perturb, perturb};
                warps = Alignment.inverseComposition2DTranslationKeypointsCpImgs(im1, im2, xYInit, maxIter, eps,
                        yXKeypoints, hX, hY, Alignment.Type.TRANSLATION_2D);
                if (printDebug)
                    System.out.printf("2D trans (sub-section image copies):\n  nIter=%d,\n%s\n%s\n",
                        warps.nIterMax, FormatArray.toString(warps.warps[chk], "%.2f"),
                        FormatArray.toString(warps.warps[yXKeypoints.length/2], "%.2f"));
                //assertEquals(iTest, (int)Math.round(warps.warps[chk][0][2]));
                //assertEquals(iTest, (int)Math.round(warps.warps[chk][1][2]));

                ///*
                xYInit = new double[]{perturb, perturb};
                warps = Alignment.inverseComposition2DTranslationKeypoints(im1, im2, xYInit, maxIter, eps,
                        yXKeypoints, hX, hY, Alignment.Type.TRANSLATION_2D);

                if (printDebug)
                    System.out.printf("2D trans (image windows):\n  nIter=%d\n%s\n",
                        warps.nIterMax, FormatArray.toString(warps.warps[chk], "%.2f"));

                //assertEquals(iTest, (int)Math.round(warps.warps[chk][0][2]));
                //assertEquals(iTest, (int)Math.round(warps.warps[chk][1][2]));
                //*/

                ///* Test 2D AFFINE

                warps = Alignment.inverseCompositionKeypointsCpImgs(im1, im2,
                        new double[][]{{1, 0, perturb}, {0, 1, perturb}}, maxIter, eps,
                        yXKeypoints, hX, hY, Alignment.Type.AFFINE_2D);
                if (printDebug)
                    System.out.printf("2D AFFINE (sub-section image copies)):\n  nIter=%d\n%s\n",
                        warps.nIterMax, FormatArray.toString(warps.warps[chk], "%.2f"));

                warps = Alignment.inverseCompositionKeypointsCpImgs2DAffinePre2DTrans(
                        im1, im2, maxIter, eps, yXKeypoints, hX, hY);
                if (printDebug)
                    System.out.printf("2D trans + 2D affine (sub-section image copies):\n  nIter=%d,\n%s\n",
                        warps.nIterMax, FormatArray.toString(warps.warps[chk], "%.2f"));

                warps = Alignment.inverseCompositionKeypoints(im1, im2,
                        new double[][]{{1, 0, perturb}, {0, 1, perturb}}, maxIter, eps,
                        yXKeypoints, hX, hY, Alignment.Type.AFFINE_2D);
                if (printDebug)
                    System.out.printf("2D AFFINE (image windows)):\n  nIter=%d\n%s\n",
                        warps.nIterMax, FormatArray.toString(warps.warps[chk], "%.2f"));

                warps = Alignment.inverseCompositionKeypoints2DAffinePre2DTrans(
                        im1, im2, maxIter, eps, yXKeypoints, hX, hY);
                if (printDebug) {
                    System.out.printf("2D trans + 2D affine (image windows):\n  nIter=%d,\n%s\n",
                            warps.nIterMax, FormatArray.toString(warps.warps[chk], "%.2f"));
                    System.out.flush();
                }

                //assertEquals(iTest, (int)Math.round(pInit[0][2]));
                //assertEquals(iTest, (int)Math.round(pInit[1][2]));
                //*/

                /*
                 e.g. a first pass through
                 Alignment.inverseCompositionKeypointsCpImgs(im1, im2, pInit, maxIter, eps, yXKeypoints, patchHalfWidth, Type.AFFINE_2D);
                    resulting in
                         0.53, 0.47, -3.57
                         0.47, 0.53, -3.57
                         0.00, 0.00, 1.00
                 is not a feasible combination of th and sh_x, sh_y for scale in x and y = 1
                     [ sh_x * sin(th)   sh_x * cos(th) ]  = [ 0.53  0,47 ]
                     [ sh_y * cos(th)  -sh_y * sin(th)  ]    [0.47   0.53 ]
                 */
            }
        }

    }

    private GreyscaleImage makeGSImg(double[][] im) {
        GreyscaleImage gsImg = new GreyscaleImage(im[0].length, im.length);
        for (int x = 0; x < im[0].length; ++x) {
            for (int y = 0; y < im.length; ++y) {
                // scale v to be within 0 to 155
                int v = (int)(255 * im[y][x]);
                if (v > 255)
                    v = 255;
                gsImg.setValue(x, y, v);
            }
        }
        return gsImg;
    }

    protected void writeToPng(double[][] im1, String fileName) throws IOException {
        // scale from 0 to 255
        int h = im1.length;
        int w = im1[0].length;

        GreyscaleImage image = new GreyscaleImage(w, h);
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                int v = (int)Math.round(im1[y][x]);
                image.setValue(x, y, v);
            }
        }
        String sep = System.getProperty("file.separator");
        String filePath = ResourceFinder.findOutputTestDirectory() + sep + fileName;
        ImageIOHelper.writeOutputImage(filePath, image);
    }

    private List<PairInt> findKeyPoints(ImageExt img, int nCorners, int patchSize) {
        //double[][] b1 = luvTheta.exportRowMajor();

        /*canny = new CannyEdgeFilterAdaptive();
        canny.overrideToUseAdaptiveThreshold();
        canny.overrideToNotUseLineThinner();
        canny.applyFilter(img.copyToGreyscale2());
        EdgeFilterProducts prod = canny.getFilterProducts();
        GreyscaleImage gsImg = prod.getGradientXY();
        gsImg.normalizeToMax255();
        double[][] b1 = gsImg.exportRowMajor();*/

        // NOTE: phase congruency is useful to smooth out noisy details in trees, leaves, etc.
        // the edge map from it can be used to find good corners.
        ImageProcessor imageProcessor = new ImageProcessor();
        PhaseCongruencyDetector phaseCDetector = new PhaseCongruencyDetector();
        phaseCDetector.setToDebug();
        PhaseCongruencyDetector.PhaseCongruencyProducts products =
                phaseCDetector.phaseCongMono(img.copyToGreyscale2());

        int buffer = 5;
        int len = buffer + patchSize/2;

        int[][] thinned = products.getThinned();
        GreyscaleImage edges = img.createWithDimensions().copyToGreyscale();
        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                if (thinned[j][i] > 0) {
                    if (i - len < 0 || i+len >= img.getWidth() || j - len < 0 || j+len >= img.getHeight())
                        continue;
                    edges.setValue(i, j, thinned[j][i]);
                }
            }
        }

        // ORB uses Harris corners among many things to choose keypoints:
        ORB orb1 = new ORB(edges, nCorners);
        orb1.overrideToNotCreateDescriptors();
        orb1.overrideToUseSingleScale();
        orb1.detectAndExtract();
        List<PairInt> kp1 = orb1.getAllKeyPoints();

          // uncomment to plot keypoints over a copy of the image
        //Image tmp = edges.copyToColorGreyscale();
        Image tmp = img.copyImage();

        MiscDebug.writeImage(edges, "_thinned_" + "_" + "_phasecong_" + "_");
        for (int i = 0; i < kp1.size(); ++i) {
            int x = kp1.get(i).getX();
            int y = kp1.get(i).getY();
            if (x - len < 0 || x+len >= img.getWidth() || y - len < 0 || y+len >= img.getHeight())
                continue;
            ImageIOHelper.addPointToImage(x, y, tmp, 2, 255, 0, 0);
        }
        MiscDebug.writeImage(tmp, "_thinned_" + "_" + "_phasecong_kp_" +
                MiscDebug.getCurrentTimeFormatted());

        return kp1;
    }

    protected static void dilateAboveBinaryThreshold(double[][]im1) {
        int m = im1.length;
        int n = im1[0].length;
        OtsuThresholding oT = new OtsuThresholding();
        double thresh = oT.calculateBinaryThreshold2D(im1, 2);
        double[][] im1Copy = MatrixUtil.copy(im1);
        for (int y = 0; y < m; ++y) {//y
            for (int x = 0; x < n; ++x) {//x
                if (im1[y][x] > thresh) {
                    double v = im1Copy[y][x];
                    // dilate
                    for (int i = -1; i <= 1; ++i) {
                        int x2 = x + i;
                        if (x2 < 0 || x2 > (n - 1)) {
                            continue;
                        }
                        for (int j = -1; j <= 1; ++j) {
                            if (i == 0 && j == 0) {
                                continue;
                            }
                            int y2 = y + j;
                            if (y2 < 0 || y2 > (m - 1)) {
                                continue;
                            }
                            im1[y2][x2] = v;
                        }
                    }
                }
            }
        }
    }

    protected static void binaryThreshold(double[][]im1) {
        int m = im1.length;
        int n = im1[0].length;
        OtsuThresholding oT = new OtsuThresholding();
        double thresh = oT.calculateBinaryThreshold2D(im1, 2);
        double[][] im1Copy = MatrixUtil.copy(im1);
        for (int y = 0; y < m; ++y) {//y
            for (int x = 0; x < n; ++x) {//x
                if (im1[y][x] > thresh) {
                    im1[y][x] = 1;
                } else {
                    im1[y][x] = 0;
                }
            }
        }
    }

    protected static void blur(double[][] im1, SIGMA sigma) {

        float[] kernel = Gaussian1D.getKernel(sigma);
        double[] kernelD = new double[kernel.length];
        for (int i = 0; i < kernel.length; ++i) {
            kernelD[i] = kernel[i];
        }

        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.applyKernel1D(im1, kernelD, true);
        imageProcessor.applyKernel1D(im1, kernelD, false);
    }

    protected void check(double[][] a, double min, double max) {
        for (int i = 0; i < a.length; ++i) {
            for (int j =0; j < a[0].length; ++j) {
                if (a[i][j] < min || a[i][j] > max) {
                    int t = 2;
                }
            }
        }
    }
}
