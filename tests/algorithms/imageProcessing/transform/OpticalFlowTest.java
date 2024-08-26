package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.*;
import algorithms.imageProcessing.features.PhaseCongruencyDetector;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.*;
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
 *
 * @author nichole
 */
public class OpticalFlowTest extends TestCase {
    protected class LKResults {
        double[] mode;
        double[] mean;
        double[] median;
        double[] s;
        boolean[] modeIsUnique;
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("mode=");
            if (mode != null) {
                sb.append(String.format("%s", FormatArray.toString(mode, "%.3f")));
                sb.append(String.format(" (%s)", Arrays.toString(modeIsUnique)));
            }
            sb.append("\nmedian=");
            if (median != null) {
                sb.append(String.format("%s", FormatArray.toString(median, "%.3f")));
            }
            sb.append("\ns=");
            if (s != null) {
                sb.append(String.format("%s", FormatArray.toString(s, "%.3f")));
            }
            sb.append("\nmean=");
            if (mean != null) {
                sb.append(String.format("%s", FormatArray.toString(mean, "%.3f")));
            }
            return sb.toString();
        }
    }

    public OpticalFlowTest() {
    }
    
    public void estHornSchunck() throws Exception {

        // do not use images larger than 128 with horn schunck method

        // the tomato images are a difficult example as there is differential motion across the image
        // to (rigid body motion rather than purely translational).

        String[] fileNames = new String[]{"tomatoes_01_128.png", "tomatoes_03_128.png"};
        /*
         for the 01 and 03 of width 128 pix
        63,38   52,36     dx=11, dy=2
        76.42   65,40
        102,76   95,76    dx=7, dy=0
         */

        ImageProcessor imageProcessor = new ImageProcessor();
        SIGMA blur = null;//SIGMA.ONE;

        double uInit = 0;
        double vInit = 0;

        for (int i = 0; i < fileNames.length-1; ++i) {
            String fileName = fileNames[i];
            String filePath = ResourceFinder.findFileInTestResources(fileName);
            ImageExt img = ImageIOHelper.readImageExt(filePath);
            int w = img.getWidth();
            int h = img.getHeight();

            if (blur != null)
                imageProcessor.blur(img, blur, 0, 255);

            //double[][] b1 = img.exportHSBRowMajor().get(2);
            double[][] b1 = img.exportRowMajor().get(0);

            fileName = fileNames[i+1];
            filePath = ResourceFinder.findFileInTestResources(fileName);
            img = ImageIOHelper.readImageExt(filePath);

            if (blur != null)
                imageProcessor.blur(img, blur, 0, 255);

            //double[][] b2 = img.exportHSBRowMajor().get(2);
            double[][] b2 = img.exportRowMajor().get(0);

            List<double[][]> uvHS = OpticalFlow.hornSchunck(b1, b2);

            double[][] uH = Histogram.createHistogram(MatrixUtil.stack(uvHS.get(0)), 255);
            double[][] vH = Histogram.createHistogram(MatrixUtil.stack(uvHS.get(1)), 255);
            int uIdx = MiscMath0.findYMaxIndex(uH[1]);
            int vIdx = MiscMath0.findYMaxIndex(vH[1]);
            double uMode = uH[0][uIdx];
            double vMode = vH[0][vIdx];

            int[] uc = new int[uH[0].length];
            int[] uv = new int[uH[0].length];
            for (int j = 0; j < uc.length; ++j) {
                uc[j] = (int)Math.round(uH[1][j]);
                uv[j] = (int)Math.round(uH[0][j]);
            }
            MiscSorter.sortByDecr(uc, uv);

            /*System.out.printf("from u hist: sorted top 10\n");
            for (int j = 0; j <= Math.min(10, uc.length); ++j) {
                System.out.printf("    count=%d, u=%d\n", uc[j], uv[j]);
            }*/
            System.out.printf("uMode=%d, vMode=%d\n", (int)Math.round(uMode), (int)Math.round(vMode));
        }
    }

    public void estLucasKanade() throws Exception {

        // the tomato images are a difficult example as there is differential motion across the image
        // to (rigid body motion rather than purely translational).

        //even when smoothed by gaussian with sigma=4, and choosing best keypoints
        // and a large enough patch size to have overlapping pixels between images,
        // the method doesn't perform well for expected differences as high as 20 in x

        // each 300x225
        String[] fileNames = new String[]{"tomatoes_01.png", /*"tomatoes_02.png",*/
                "tomatoes_03.png"};
        //fileNames = new String[]{"tomatoes_01_128.png", "tomatoes_03_128.png"};

        /*
        a feature in the 3 tomato images large:
        177, 140  in 01.png
        177, 139  in 02.png
        157, 135  in 03.png  <== dx/t02 = 20, dy/t02=5

        for the 01 and 03 of width 128 pix
        63,38   52,36     dx=11, dy=2
        76.42   65,40
        102,76   95,76    dx=7, dy=0

         */
        List<double[][]> images = new ArrayList<>();

        int nCorners = 50;

        int patchSize = 5;//30;
        //patchSize = 15; // for 128 pix width image

        SIGMA blur = null;//SIGMA.FOUR;

        ImageProcessor imageProcessor = new ImageProcessor();
        CannyEdgeFilterAdaptive canny;

        for (int i = 0; i < fileNames.length-1; ++i) {
            String fileName = fileNames[i];
            String filePath = ResourceFinder.findFileInTestResources(fileName);
            ImageExt img = ImageIOHelper.readImageExt(filePath);
            int w = img.getWidth();
            int h = img.getHeight();

            if (blur != null)
                imageProcessor.blur(img, blur, 0, 255);

            //GreyscaleImage luvTheta = imageProcessor.createCIELUVTheta(img, 255);
            double[][] b1 = img.exportHSBRowMajor().get(2);

            double min = MiscMath0.findMin(b1);
            double max = MiscMath.findMax(b1);
            double range = max - min;

            GreyscaleImage gsImg = new GreyscaleImage(b1[0].length, b1.length);
            for (int x = 0; x < b1[0].length; ++x) {
                for (int y = 0; y < b1.length; ++y) {
                    // scale v to be within 0 to 155
                    int v = 255 * (int)Math.floor( (b1[y][x] - min)/range );
                    gsImg.setValue(x, y, v);
                }
            }

            int[][] patchCenters;

            // x, y
            if (true) {
                List<PairInt> kp1 = findKeyPoints(img, nCorners, patchSize);
                patchCenters = new int[kp1.size()][2];
                for (int ii = 0; ii < kp1.size(); ++ii) {
                    patchCenters[ii][0] = kp1.get(ii).getX();
                    patchCenters[ii][1] = kp1.get(ii).getY();
                }
            } else {
                // override patches manually for now for width=300 pix
                if (false) {
                    patchCenters = new int[4][2];
                    patchCenters[0] = new int[]{176, 140};
                    patchCenters[1] = new int[]{215, 146};
                    patchCenters[2] = new int[]{184, 171};
                    patchCenters[3] = new int[]{126, 121};
                } else  {
                    patchCenters = new int[3][2];
                    patchCenters[0] = new int[]{198, 103};
                    patchCenters[1] = new int[]{146, 151};
                    patchCenters[2] = new int[]{128, 100};
                }
            }

            // override for now for width = 128 pix
            /*patchCenters = new int[3][2];
            patchCenters[0] = new int[]{63,38};
            patchCenters[1] = new int[]{76,42};
            patchCenters[2] = new int[]{102,76};
           */

            fileName = fileNames[i+1];
            filePath = ResourceFinder.findFileInTestResources(fileName);
            img = ImageIOHelper.readImageExt(filePath);

            if (blur != null)
                imageProcessor.blur(img, blur, 0, 255);

            double[][] b2 = img.exportHSBRowMajor().get(2);

            // quick look at keypoints for b2 alone
            findKeyPoints(img, nCorners, patchSize);

            // expect u=-20, v=-5 for 01 to 03.   large motion, so needs a multiplscal pyrimdal approach
            List<double[]> uvLK = OpticalFlow.lucasKanade(b1, b2, patchCenters, patchSize);

            LKResults results = process(uvLK);

            System.out.printf("%s\n", results.toString());

        }
    }

    public void __test2LucasKanade() throws IOException, NotConvergedException {
        int patchDim = 5;
        int m = 32;
        int pixMax = m;
        //System.out.printf("pixMax=%d, m=%d\n", pixMax, m);
        int n = m;
        double[][] im1;
        double[][] im2;

        /*
        using the median of the calculated u's, v's and a generous error of 10% of image dimension:

            m     patchDim   max value    allowing error m/10   nPatches  can calc up to max u
            32      5         32             3                     21         23
            64      5         64             6                     53         57
            128     5         128            12                    117        120
            256     5         256            25                    245        250   (largest error=7)
            512     5         512            51                    501        506   (largest error=7)
            1024    5         1024           102                   1013       1018  (largest error=7)
        */

        for (int iTest = 1; iTest < m-9 /*m-5*/; ++iTest) {

            im1 = new double[m][n];
            im2 = new double[m][n];
            double deltaI = pixMax / m;

            for (int i = 1; i < m; ++i) {
                double b = deltaI * i;
                b = (int)Math.round(b); // to match image processing temporarily
                for (int j = 0; j <= i; ++j) {
                    im1[i][j] = b;
                    im1[j][i] = b;
                    if (i + iTest >= 0 && i + iTest < m && j + iTest >= 0 && j + iTest < n) {
                        im2[i + iTest][j + iTest] = b;
                    }
                    if (j + iTest >= 0 && j + iTest < m && i + iTest >= 0 && i + iTest < n) {
                        im2[j + iTest][i + iTest] = b;
                    }
                }
            }

            // start at 1+patchDim.  end at (m-1)-patchDim // 6 32-5
            int pIdxI = 1 + patchDim;
            int pIdxF = (m - 1) - patchDim;
            int pN = pIdxF - pIdxI + 1;

            int[][] patchCenters = new int[pN][2];
            for (int ii = 0; ii < pN; ++ii) {
                patchCenters[ii][0] = pIdxI + ii;
                patchCenters[ii][1] = pIdxI + ii;
            }

            //System.out.printf("im1:\n%s\n", FormatArray.toString(im1, "%.3f"));
            //System.out.printf("im2:\n%s\n", FormatArray.toString(im2, "%.3f"));

            List<double[]> uvLK = OpticalFlow.lucasKanade(im1, im2, patchCenters, patchDim);

            if (uvLK.isEmpty()) {
                System.out.println("Test=" + iTest + " failed to process any patches");
                fail();
            }

            LKResults results = process(uvLK);
            //System.out.println(results.toString());
            //System.out.flush();

            //TODO: calculate an error
            double tolerance = (m < 128) ? 4 : 7;

            assertTrue(Math.abs(results.median[0] - iTest) <= tolerance);
            assertTrue(Math.abs(results.median[1] - iTest) <= tolerance);

        }
    }

    protected LKResults process(List<double[]> uvLK) {

        LKResults results = new LKResults();

        // calc mode of patches
        Map<Integer, Integer> uCounts = new HashMap<>();
        Map<Integer, Integer> vCounts = new HashMap<>();

        double[] u = new double[uvLK.size()];
        double[] v = new double[uvLK.size()];
        for (int ii = 0; ii < u.length; ++ii) {
            int _u = (int)Math.round(uvLK.get(ii)[0]);
            int _v = (int)Math.round(uvLK.get(ii)[1]);
            u[ii] = _u;
            v[ii] = _v;

            //System.out.printf("(%.3f, %.3f) ", u[ii], v[ii]);
            //System.out.flush();

            uCounts.put(_u, uCounts.getOrDefault(_u, 0) + 1);
            vCounts.put(_v, vCounts.getOrDefault(_v, 0) + 1);
        }
       // System.out.printf("\n");

        int[] vals = new int[uCounts.size()];
        int[] cts = new int[uCounts.size()];
        int j = 0;
        for (Map.Entry<Integer, Integer> entry : uCounts.entrySet()) {
            vals[j] = entry.getKey();
            cts[j++] = entry.getValue();
        }
        MiscSorter.sortBy1stArg(cts, vals);
        // for larger displacements with small patches, cannot use mode:
        int uMode = vals[vals.length - 1];
        boolean uniqueU = vals.length == 1 || (cts[vals.length - 1] != cts[vals.length - 2]);
        //assertEquals(iTest, -uMode);

        vals = new int[vCounts.size()];
        cts = new int[vCounts.size()];
        j = 0;
        for (Map.Entry<Integer, Integer> entry : vCounts.entrySet()) {
            vals[j] = entry.getKey();
            cts[j++] = entry.getValue();
        }
        MiscSorter.sortBy1stArg(cts, vals);
        int vMode = vals[vals.length - 1];
        boolean uniqueV = vals.length == 1 || (cts[vals.length - 1] != cts[vals.length - 2]);
        // for larger displacements with small patches, cannot use mode:
        //assertEquals(iTest, -vals[vals.length - 1]);
        results.mode = new double[]{uMode, vMode};
        results.modeIsUnique = new boolean[]{uniqueU, uniqueV};

        results.median = new double[2];
        results.s = new double[2];
        results.mean = new double[2];

        // use the median and an error
        double[] mADMinMax = MiscMath0.calculateMedianOfAbsoluteDeviation(u);
        double[] avgStd = MiscMath.getAvgAndStDev(u);
        double s = 1.4826*mADMinMax[0];
        results.median[0] = mADMinMax[1];
        results.s[0] = s;

        mADMinMax = MiscMath0.calculateMedianOfAbsoluteDeviation(v);
        avgStd = MiscMath.getAvgAndStDev(v);
        s = 1.4826*mADMinMax[0];
        results.median[1] = mADMinMax[1];
        results.s[1] = s;
        results.mean[1] = avgStd[0];

        return results;
    }

    public void __test2HornSchunck() throws IOException {

        /*
        some empirically derived numbers from unit tests:
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
        */
        int pixMax = 32;
        int m = 32;
        //System.out.printf("pixMax=%d, m=%d\n", pixMax, m);
        int n = m;
        double[][] im1;
        double[][] im2;
        double uInit = 0;//0.5;
        double vInit = 0;//0.5;
        int maxIter = 100;
        double epsSq = 1E-2;

        double alpha = 1E-1;//1E1;

        for (int iTest = 1; iTest < 14; ++iTest) {
            //m = iTest * 10;
            //n = m;
            im1 = new double[m][n];
            im2 = new double[m][n];
            double deltaI = pixMax / m;

            //alpha = iTest;

            for (int i = 1; i < m; ++i) {
                double b = deltaI * i;
                b = (int)Math.round(b); // to match image processing temporarily
                for (int j = 0; j <= i; ++j) {
                    im1[i][j] = b;
                    im1[j][i] = b;
                    if (i + iTest >= 0 && i + iTest < m && j + iTest >= 0 && j + iTest < n) {
                        im2[i + iTest][j + iTest] = b;
                    }
                    if (j + iTest >= 0 && j + iTest < m && i + iTest >= 0 && i + iTest < n) {
                        im2[j + iTest][i + iTest] = b;
                    }
                }
            }

            //System.out.printf("im1:\n%s\n", FormatArray.toString(im1, "%.3f"));
            //System.out.printf("im2:\n%s\n", FormatArray.toString(im2, "%.3f"));

            List<double[][]> uvHS = OpticalFlow.hornSchunck(im1, im2, uInit, vInit, alpha,
                maxIter, epsSq);

            double[][] uH = Histogram.createHistogram(MatrixUtil.stack(uvHS.get(0)), 255);
            double[][] vH = Histogram.createHistogram(MatrixUtil.stack(uvHS.get(1)), 255);
            int uIdx = MiscMath0.findYMaxIndex(uH[1]);
            int vIdx = MiscMath0.findYMaxIndex(vH[1]);
            double uMode = uH[0][uIdx];
            double vMode = vH[0][vIdx];

            double avgU = MiscMath0.mean(MatrixUtil.rowMeans(
                    MatrixUtil.copySubMatrix(uvHS.get(0), 1, m - 1, 1, n - 1)
            ));
            double avgV = MiscMath0.mean(MatrixUtil.rowMeans(
                    MatrixUtil.copySubMatrix(uvHS.get(1), 1, m - 1, 1, n - 1)
            ));

            //System.out.printf("Test=%d, uEst=%.3f, vEst=%.3f (avgs=%.3f,%.3f) (uInit,vInit=%.1f,%.1f), alpha=%.4e\n",
            //        iTest, uMode, vMode, avgU, avgV, uInit, vInit, alpha);
            //System.out.printf("u:\n%s\n", FormatArray.toString(uvHS.get(0), "%.3f"));
            //System.out.printf("v:\n%s\n", FormatArray.toString(uvHS.get(1), "%.3f"));
            //System.out.flush();

            assertEquals(iTest, Math.round(uMode));
            assertEquals(iTest, Math.round(vMode));

            //writeToPng(im1, "im1_" + iTest + ".png");
            //writeToPng(im2, "im2_" + iTest + ".png");
        }

    }

    public void __test2InverseCompositionAlignment() throws IOException, NotConvergedException {
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
            if (images == null || images.yXPatches.length == 0) {
                break;
            }
            im1 = images.im1;
            im2 = images.im2;

            // set all to 0 to solve for
            //MatrixUtil.fill(pInit, 0.);

            System.out.printf("\niTest=%d\n", iTest);

            double[] errSSD;
            ///*  TEST 2-D translation
            double[] xYInit = new double[]{0, 0};
            errSSD = Alignment.inverseCompositional2DTranslation(im1, im2, xYInit, maxIter, eps);
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
        int[][] yXPatches;
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

        int nPatches = Math.min(m,n) - 2*Math.min(hX, hY)-1;
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
        images.yXPatches = yXPatches;
        return images;
    }

    private static void fillWithRandom(double[][] im, int minVal, int maxVal, Random rand) {
        for (int i = 0; i < im.length; ++i) {//x
            for (int j = 0; j < im[i].length; ++j) {//x
                im[i][j] = minVal + rand.nextInt(maxVal - minVal);
            }
        }

    }

    private Images generateImages2(int m, int n, double dX, double dY, double pixMax) {
        double[][] im1 = new double[m][n];
        double[][] im2 = new double[m][n];

        //the first image will be a random pattern, the 2nd will be im1 warped by projection matrix

        long seed = System.nanoTime();
        System.out.println("seed=" + seed);
        Random rand = new Random(seed);

        fillWithRandom(im1, 1, 255, rand);
        fillWithRandom(im2, 1, 255, rand);

        int iP = 0;

        int hX = (int)(dX);
        int hY = (int)(dY);
        if (hX < 2) {
            hX = 2;
        }
        if (hY < 2) {
            hY = 2;
        }

        int nPatches = Math.min(m,n) - 2*Math.min(hX, hY)-1;
        if (nPatches < 1) return null;
        int[][] yXPatches = new int[nPatches][];

        double[] i2j2 = new double[3];

        double[][] pInit = MatrixUtil.createIdentityMatrix(3);
        pInit[0][2] = dX;
        pInit[1][2] = dY;

        for (int i = 0; i < m; ++i) {//x
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
        images.yXPatches = yXPatches;
        return images;
    }

    public void test3InverseCompositionAlignment() throws IOException, NotConvergedException {
        System.out.printf("test3InverseCompositionAlignment\n");

        int pixMax = 32*2;
        int m = 32*2;
        //System.out.printf("pixMax=%d, m=%d\n", pixMax, m);
        int n = m;
        double[][] im1;
        double[][] im2;
        int maxIter = 100;//1000;
        double eps = 0.01;

        double perturb = 1E-3;

        for (int iTest = 1; iTest < (m/3); ++iTest) {

            double dX = iTest;
            double dY = iTest;

            Images images = generateImages(m,  n, dX, dY, pixMax);
            if (images == null || images.yXPatches.length == 0) {
                break;
            }
            im1 = images.im1;
            im2 = images.im2;
            int hX = images.hX;
            int hY = images.hY;

            int[][] yXPatches = images.yXPatches;

            System.out.printf("\niTest=%d\n", iTest);

            //System.out.printf("im1=\n%s", FormatArray.toString(im1, "%.2f"));
            //System.out.printf("im2=\n%s", FormatArray.toString(im2, "%.2f"));

            int chk = yXPatches.length/2;

            Alignment.Warps warps;
            ///*  TEST 2-D translation
            double[] xYInit = new double[]{perturb, perturb};
            warps = Alignment.inverseComposition2DTranslationKeypointsCpImgs(im1, im2, xYInit, maxIter, eps,
                    yXPatches, hX, hY, Alignment.Type.TRANSLATION_2D);
            System.out.printf("2D trans (sub-section image copies):\n  nIter=%d,\n%s\n",
                    warps.nIterMax, FormatArray.toString(warps.warps[chk], "%.2f"));
            assertEquals(iTest, (int)Math.round(warps.warps[chk][0][2]));
            assertEquals(iTest, (int)Math.round(warps.warps[chk][1][2]));

            xYInit = new double[]{perturb, perturb};
            warps = Alignment.inverseComposition2DTranslationKeypoints(im1, im2, xYInit, maxIter, eps,
                    yXPatches, hX, hY, Alignment.Type.TRANSLATION_2D);
            System.out.printf("2D trans (image windows):\n  nIter=%d\n%s\n",
                    warps.nIterMax, FormatArray.toString(warps.warps[chk], "%.2f"));
            assertEquals(iTest, (int)Math.round(warps.warps[chk][0][2]));
            assertEquals(iTest, (int)Math.round(warps.warps[chk][1][2]));

            //*/

            ///* Test 2D AFFINE

            warps = Alignment.inverseCompositionKeypointsCpImgs(im1, im2,
                    new double[][]{{1,0,perturb},{0,1,perturb}}, maxIter, eps,
                    yXPatches, hX, hY, Alignment.Type.AFFINE_2D);
            System.out.printf("2D AFFINE (sub-section image copies)):\n  nIter=%d\n%s\n",
                    warps.nIterMax, FormatArray.toString(warps.warps[chk], "%.2f"));

            warps = Alignment.inverseCompositionKeypoints(im1, im2,
                    new double[][]{{1,0,perturb},{0,1,perturb}}, maxIter, eps,
                    yXPatches, hX, hY, Alignment.Type.AFFINE_2D);
            System.out.printf("2D AFFINE (image windows)):\n  nIter=%d\n%s\n",
                    warps.nIterMax, FormatArray.toString(warps.warps[chk], "%.2f"));
            System.out.flush();

            //assertEquals(iTest, (int)Math.round(pInit[0][2]));
            //assertEquals(iTest, (int)Math.round(pInit[1][2]));
            //*/

            /*
             e.g. a first pass through
             Alignment.inverseCompositionKeypointsCpImgs(im1, im2, pInit, maxIter, eps, yXPatches, patchHalfWidth, Type.AFFINE_2D);
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

}
