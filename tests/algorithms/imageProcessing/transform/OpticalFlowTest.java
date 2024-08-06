package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.*;
import algorithms.imageProcessing.features.PhaseCongruencyDetector;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.imageProcessing.segmentation.NormalizedCuts;
import algorithms.imageProcessing.segmentation.SLICSuperPixels;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.misc.MiscMath0;
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
    
    public void testHornSchunck() throws Exception {

        // each 300x225
        String[] fileNames = new String[]{"tomatoes_01.png", /*"tomatoes_02.png",*/ "tomatoes_03.png"};

        List<double[][]> images = new ArrayList<>();

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

    public void testLucasKanade() throws Exception {

        // each 300x225
        String[] fileNames = new String[]{"tomatoes_01.png", /*"tomatoes_02.png",*/ "tomatoes_03.png"};

        List<double[][]> images = new ArrayList<>();

        int nCorners = 500;

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
            List<PairInt> kp1 = findKeyPoints(img, nCorners);

            int[][] patchCenters = new int[kp1.size()][2];
            for (int ii = 0; ii < kp1.size(); ++ii) {
                patchCenters[ii][0] = kp1.get(ii).getX();
                patchCenters[ii][1] = kp1.get(ii).getY();
            }

            fileName = fileNames[i+1];
            filePath = ResourceFinder.findFileInTestResources(fileName);
            img = ImageIOHelper.readImageExt(filePath);

            double[][] b2 = img.exportHSBRowMajor().get(2);
            //GreyscaleImage = luvTheta = imageProcessor.createCIELUVTheta(img, 255);
            //double[][] b2 = luvTheta.exportRowMajor();
            /*
            canny = new CannyEdgeFilterAdaptive();
            canny.overrideToUseAdaptiveThreshold();
            canny.overrideToNotUseLineThinner();
            canny.applyFilter(img.copyToGreyscale2());
            prod = canny.getFilterProducts();
            gsCanny = prod.getGradientXY();
            gsCanny.normalizeToMax255();
            double[][] b2 = gsCanny.exportRowMajor();
            */

            List<double[]> uvLK = OpticalFlow.lucasKanade(b1, b2, patchCenters, 5);

            int t = 1;

        }
    }

    private List<PairInt> findKeyPoints(ImageExt img, int nCorners) {
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

        int[][] thinned = products.getThinned();
        GreyscaleImage edges = img.createWithDimensions().copyToGreyscale();
        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                if (thinned[j][i] > 0) {
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
            ImageIOHelper.addPointToImage(kp1.get(i).getX(), kp1.get(i).getY(), tmp, 2, 255, 0, 0);
        }
        MiscDebug.writeImage(tmp, "_thinned_" + "_" + "_phasecong_kp_" + "_");

        return kp1;
    }

}
