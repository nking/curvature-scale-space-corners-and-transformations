package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.imageProcessing.features.orb.ORB.Descriptors;
import algorithms.imageProcessing.transform.EpipolarNormalizationHelper;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscDebug;
import algorithms.util.CorrespondencePlotter;

import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * a class to use harris corners filtered to local maxima and ORB descriptors to find matching
 * keypoints between 2 images.  Not that the images may need pre-processing so that at least
 * 50% of the potential features are common to both.
 * Texture matching with HOGs may return better results.
 * Also note that alternative means may be necessary for extremely repetetive patterns.
 */
public class CorrespondenceMaker {

    public static class CorrespondenceList {
        /**
         * points from image 1 that are matched to x2 from image 2.
         * points in format [3XN] where N is the number of matched points.
         * the first row is x-axis values, the second row is 'y-axis' values,
         * and the third row is z-axis values.
         * the z-axis is all '1's for homogeneous coordinates.
         */
        public double[][] x1;
        /**
         * points from image 2 that are matched to x1 from image 1.
         * points in format [3XN] where N is the number of matched points.
         * the first row is x-axis values, the second row is 'y-axis' values,
         * and the third row is z-axis values.
         * the z-axis is all '1's for homogeneous coordinates.
         */
        public double[][] x2;

        /**
         * if epipolar fit succeeded at RANSAC stage of outlier removal, the epipolar fits the fundamental matrix fm
         */
        public double[][] fm;

        /**
         * if epipolar fit succeeded at RANSAC stage of outlier removal, the errors for the fit were
         * calculated.
         */
        public double[] errors;
    }

    public static CorrespondenceList findUsingORB(String filePath1, String filePath2, int nCorners,
                  boolean debug) throws IOException {

        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        boolean useToleranceAsStatFactor = true;
        boolean recalcIterations = false;// possibly faster if set to true
        double tol = 2.;
        ErrorType errorType = ErrorType.SAMPSONS;

        return _findUsingORB(img1.copyToGreyscale2(), img2.copyToGreyscale2(), nCorners, useToleranceAsStatFactor,
                recalcIterations, tol, errorType, debug);
    }

    public static CorrespondenceList findUsingORB(GreyscaleImage image1, GreyscaleImage image2, int nCorners,
                                                  boolean debug) throws IOException {

        boolean useToleranceAsStatFactor = true;
        boolean recalcIterations = false;// possibly faster if set to true
        double tol = 2.;
        ErrorType errorType = ErrorType.SAMPSONS;

        return _findUsingORB(image1, image2, nCorners, useToleranceAsStatFactor,
                recalcIterations, tol, errorType, debug);
    }
    public static CorrespondenceList findUsingORB(String filePath1, String filePath2, int nCorners,
        boolean useToleranceAsStatFactor, boolean recalcIterations, double tol, ErrorType errorType,
        boolean debug) throws IOException {

        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        return _findUsingORB(img1.copyToGreyscale2(), img2.copyToGreyscale2(), nCorners, useToleranceAsStatFactor,
                recalcIterations, tol, errorType, debug);
    }

    public static CorrespondenceList findUsingORB(GreyscaleImage image1, GreyscaleImage image2, int nCorners,
        boolean useToleranceAsStatFactor, boolean recalcIterations, double tol, ErrorType errorType,
        boolean debug) throws IOException {

        return _findUsingORB(image1, image2, nCorners, useToleranceAsStatFactor,
                recalcIterations, tol, errorType, debug);
    }

    private static CorrespondenceList _findUsingORB(GreyscaleImage image1, GreyscaleImage image2, int nCorners,
        boolean useToleranceAsStatFactor, boolean recalcIterations, double tol, ErrorType errorType,
        boolean debug) throws IOException {

        ORB orb1 = new ORB(image1, nCorners);
        orb1.detectAndExtract();

        ORB orb2 = new ORB(image2, nCorners);
        orb2.detectAndExtract();

        ORB.Descriptors d1 = orb1.getAllDescriptors();
        Descriptors d2 = orb2.getAllDescriptors();
        //[3Xn]
        double[][] xKP1 = orb1.getAllKeyPointsHomogenous();
        double[][] xKP2 = orb2.getAllKeyPointsHomogenous();

        if (xKP1 == null || xKP1.length == 0 || xKP1[0].length == 0){
            System.err.println("could not find corners in image 1");
            return null;
        }
        if (xKP2 == null || xKP2.length == 0 || xKP2[0].length == 0){
            System.err.println("could not find corners in image 2");
            return null;
        }

        double[][] x1 = MatrixUtil.copy(xKP1);
        double[][] x2 = MatrixUtil.copy(xKP2);
        double[][] t1 = EpipolarNormalizationHelper.unitStandardNormalize(x1);
        double[][] t2 = EpipolarNormalizationHelper.unitStandardNormalize(x2);

        int i, j;
        String ts = null;

        if (debug) {
            // print the corners lists onto images
            ts = (new SimpleDateFormat("yyyyddMMHHmmss")).format(new Date());

            Image tmp1 = image1.copyToColorGreyscale();
            for (i = 0; i < xKP1[0].length; ++i) {
                ImageIOHelper.addPointToImage((int)xKP1[0][i], (int)xKP1[1][i], tmp1, 2, 255, 0, 0);
                //System.out.printf("x1[0][%d]=%d;  x1[1][%d]=%d;\n",
                //        i, (int)xKP1[0][i], i, (int)xKP1[1][i]);
            }
            MiscDebug.writeImage(tmp1, "_" + ts + "_corners_im1");
            tmp1 = image2.copyToColorGreyscale();
            for (i = 0; i < xKP2[0].length; ++i) {
                ImageIOHelper.addPointToImage((int)xKP2[0][i], (int)xKP2[1][i], tmp1, 2, 255, 0, 0);
                //System.out.printf("x2[0][%d]=%d;  x2[1][%d]=%d;\n",
                //        i, (int)xKP2[0][i], i, (int)xKP2[1][i]);
            }
            MiscDebug.writeImage(tmp1, "_" + ts + "_corners_im2");
        }

        // greedy, local max match of points:
        ORBMatcher.FitAndCorres fitAndCorres = ORBMatcher.matchDescriptors(d1, d2, x1, x2, useToleranceAsStatFactor,
                tol, errorType, recalcIterations, false);

        if (fitAndCorres == null) {
            return null;
        }

        int idx1, idx2;
        double[][] x1M, x2M;
        int nM;
        double[][] fm = null;
        double[] errors = null;
        CorrespondencePlotter plotter = debug ?
                new CorrespondencePlotter(image1.copyToColorGreyscale(), image2.copyToColorGreyscale()) : null;
        if (fitAndCorres.mIF != null) {
            System.out.printf("%s) #Matched = %d\n", ts, fitAndCorres.mI.length);
            nM = fitAndCorres.mIF.length;
            x1M = MatrixUtil.zeros(3, nM);
            x2M = MatrixUtil.zeros(3, nM);
            fm = MatrixUtil.copy(fitAndCorres.fm);
            EpipolarNormalizationHelper.denormalizeFM(fm, t1, t2);
            errors = new double[nM];
            for (i = 0; i < nM; ++i) {
                errors[i] = fitAndCorres.errors.get(i);
                idx1 = fitAndCorres.mIF[i][0];
                idx2 = fitAndCorres.mIF[i][1];
                for (j = 0; j < 3; ++j) {
                    x1M[j][i] = (int) xKP1[j][idx1];
                    x2M[j][i] = (int) xKP2[j][idx2];
                }
                if (debug) {
                    plotter.drawLineInAlternatingColors((int)x1M[0][i], (int)x1M[1][i],
                            (int)x2M[0][i], (int)x2M[1][i], 1);
                }
            }
            System.out.printf("%s) #outlierRemovedMatched = %d\n", ts, nM);
        } else {
            // fit failed, but we have the ORBMatcher greedy matching points as the 'mI' array.
            nM = fitAndCorres.mI.length;
            x1M = MatrixUtil.zeros(3, nM);
            x2M = MatrixUtil.zeros(3, nM);
            for (i = 0; i < nM; ++i) {
                idx1 = fitAndCorres.mI[i][0];
                idx2 = fitAndCorres.mI[i][1];
                for (j = 0; j < 3; ++j) {
                    x1M[j][i] = (int) xKP1[j][idx1];
                    x2M[j][i] = (int) xKP2[j][idx2];
                }
                if (debug) {
                    plotter.drawLineInAlternatingColors((int)x1M[0][i], (int)x1M[1][i],
                            (int)x2M[0][i], (int)x2M[1][i], 1);
                }
            }
        }
        if (debug) {
            plotter.writeImage("_" + ts + "_corres");
        }

        CorrespondenceList c = new CorrespondenceList();
        c.x1 = x1M;
        c.x2 = x2M;
        c.fm = fm;
        c.errors = errors;
        return c;
    }

}
