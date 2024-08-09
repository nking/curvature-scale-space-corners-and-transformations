package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.*;
import algorithms.imageProcessing.features.PhaseCongruencyDetector;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.imageProcessing.segmentation.NormalizedCuts;
import algorithms.imageProcessing.segmentation.SLICSuperPixels;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.*;
import algorithms.sort.MiscSorter;
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

    public OpticalFlowTest() {
    }
    
    public void estHornSchunck() throws Exception {

        // each 300x225
        String[] fileNames = new String[]{"tomatoes_01.png", "tomatoes_02.png", /*"tomatoes_03.png"*/};
        /*
        a feature in the 3 tomatoe images:
        177, 140  in 01.png
        177, 139  in 02.png
        157, 135  in 03.png
         */

        double uInit = 0;
        double vInit = 0;

        for (int i = 0; i < fileNames.length-1; ++i) {
            String fileName = fileNames[i];
            String filePath = ResourceFinder.findFileInTestResources(fileName);
            ImageExt img = ImageIOHelper.readImageExt(filePath);
            int w = img.getWidth();
            int h = img.getHeight();

            //double[][] b1 = img.exportHSBRowMajor().get(2);
            double[][] b1 = img.exportRowMajor().get(0);

            fileName = fileNames[i+1];
            filePath = ResourceFinder.findFileInTestResources(fileName);
            img = ImageIOHelper.readImageExt(filePath);

            //double[][] b2 = img.exportHSBRowMajor().get(2);
            double[][] b2 = img.exportRowMajor().get(0);

            List<double[][]> uvHS = OpticalFlow.hornSchunck(b1, b2);

            int t = 2;
        }
    }

    public void estLucasKanade() throws Exception {

        // each 300x225
        String[] fileNames = new String[]{"tomatoes_01.png", /*"tomatoes_02.png",*/ "tomatoes_03.png"};
        /*
        a feature in the 3 tomato images:
        177, 140  in 01.png
        177, 139  in 02.png
        157, 135  in 03.png
         */
        List<double[][]> images = new ArrayList<>();

        int nCorners = 300;
        int patchSize = 30;

        ImageProcessor imageProcessor = new ImageProcessor();
        CannyEdgeFilterAdaptive canny;

        for (int i = 0; i < fileNames.length-1; ++i) {
            String fileName = fileNames[i];
            String filePath = ResourceFinder.findFileInTestResources(fileName);
            ImageExt img = ImageIOHelper.readImageExt(filePath);
            int w = img.getWidth();
            int h = img.getHeight();

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


            // x, y
            List<PairInt> kp1 = findKeyPoints(img, nCorners, patchSize);

            int[][] patchCenters = new int[kp1.size()][2];
            for (int ii = 0; ii < kp1.size(); ++ii) {
                patchCenters[ii][0] = kp1.get(ii).getX();
                patchCenters[ii][1] = kp1.get(ii).getY();
            }

            fileName = fileNames[i+1];
            filePath = ResourceFinder.findFileInTestResources(fileName);
            img = ImageIOHelper.readImageExt(filePath);

            double[][] b2 = img.exportHSBRowMajor().get(2);

            // expect u=-20, v=-5 for 01 to 03.   large motion, so needs a multiplscal pyrimdal approach
            List<double[]> uvLK = OpticalFlow.lucasKanade(b1, b2, patchCenters, patchSize);

            int t = 1;

        }
    }

    public void test2LucasKanade() throws IOException, NotConvergedException {
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
            // calc mode of patches
            Map<Integer, Integer> uCounts = new HashMap<>();
            Map<Integer, Integer> vCounts = new HashMap<>();

            double[] u = new double[uvLK.size()];
            double[] v = new double[uvLK.size()];
            //System.out.println("Test=" + iTest + ", pN=" + pN);
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
            //System.out.printf("\n");

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

            // use the median and an error
            double[] mADMinMaxU = MiscMath0.calculateMedianOfAbsoluteDeviation(u);
            double[] avgStd = MiscMath.getAvgAndStDev(u);
            double s = 1.4826*mADMinMaxU[0];

            // using 10% of image width or length as a generous error
            double error = m/10;
            // found that that largest needed was 7 for the number of patches I'm using and the gradient pattern
            error = Math.min(error, 7);
            //System.out.printf("expecting %d.  median=%.0f +- %.1f (s=%.1f, mean=%.1f, uMode=%d(%b))\n",
            //        -iTest, mADMinMaxU[1], error, s, avgStd[0], uMode, uniqueU);
            //System.out.flush();

            if (false /*uniqueU*/) {
                assertTrue(Math.abs(-uMode - iTest) <= error);
                assertTrue(Math.abs(-vMode - iTest) <= error);
            } else {
                assertTrue(Math.abs(-mADMinMaxU[1] - iTest) <= error);
                double[] mADMinMaxV = MiscMath0.calculateMedianOfAbsoluteDeviation(v);
                s = 1.4826 * mADMinMaxV[0];
                assertTrue(Math.abs(-mADMinMaxV[1] - iTest) <= error);
            }
        }
    }

    public void test2HornSchunck() throws IOException {

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


        ImageProcessor imageProcessor = new ImageProcessor();
        PhaseCongruencyDetector phaseCDetector = new PhaseCongruencyDetector();
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

        ORB orb1 = new ORB(edges, nCorners);
        orb1.overrideToNotCreateDescriptors();
        orb1.overrideToUseSingleScale();
        orb1.detectAndExtract();
        List<PairInt> kp1 = orb1.getAllKeyPoints();

        Image tmp = edges.copyToColorGreyscale();

        MiscDebug.writeImage(edges, "_thinned_" + "_" + "_phasecong_" + "_");
        for (int i = 0; i < kp1.size(); ++i) {
            int x = kp1.get(i).getX();
            int y = kp1.get(i).getY();
            if (x - len < 0 || x+len >= img.getWidth() || y - len < 0 || y+len >= img.getHeight())
                continue;
            ImageIOHelper.addPointToImage(x, y, tmp, 2, 255, 0, 0);
        }
        MiscDebug.writeImage(tmp, "_thinned_" + "_" + "_phasecong_kp_" + "_");

        return kp1;
    }

}
