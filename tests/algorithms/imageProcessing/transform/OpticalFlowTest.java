package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.*;
import algorithms.imageProcessing.features.PhaseCongruencyDetector;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.imageProcessing.segmentation.NormalizedCuts;
import algorithms.imageProcessing.segmentation.SLICSuperPixels;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.List;

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

            List<double[][]> uvHS = OpticalFlow.hornSchunck(b1, b2, uInit, vInit);

            int t = 2;
        }
    }

    public void estLucasKanade() throws Exception {

        // each 300x225
        String[] fileNames = new String[]{"tomatoes_01.png", /*"tomatoes_02.png",*/ "tomatoes_03.png"};
        /*
        a feature in the 3 tomatoe images:
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

    public void test2() {
        int m = 10, n = m;
        double[][] im1 = new double[m][n];
        double[][] im2 = new double[m][n];
        double deltaI = 1./m;
        for (int i = 1; i < m; ++i) {
            double b = deltaI*i;
            for (int j = 0; j <= i; ++j) {
                im1[i][j] = b;
                im1[j][i] = b;

                if (i+1 >= 0 && i+1 < m && j+1 >= 0 && j+1 < n) {
                    im2[i + 1][j + 1] = b;
                }
                if (j+1 >= 0 && j+1 < m && i+1 >= 0 && i+1 < n) {
                    im2[j + 1][i + 1] = b;
                }
            }
        }
        double uInit = 0;//0.5;
        double vInit = 0;//0.5;

        List<double[][]> uvHS = OpticalFlow.hornSchunck(im1, im2, uInit, vInit);

        System.out.printf("u=%s\n", FormatArray.toString(uvHS.get(0), "%.3f"));
        System.out.printf("v=%s\n", FormatArray.toString(uvHS.get(1), "%.3f"));

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
