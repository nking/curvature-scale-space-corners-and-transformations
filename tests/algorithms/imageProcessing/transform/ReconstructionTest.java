package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.imageProcessing.matching.ORBMatcher;
import static algorithms.imageProcessing.transform.Rotation.extractThetaFromZYX;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.FormatArray;
import algorithms.util.PairInt;
import algorithms.util.QuadInt;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertNotNull;
import no.uib.cipr.matrix.NotConvergedException;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class ReconstructionTest extends TestCase {
    
    public ReconstructionTest() {
    }
    
      /**
     * Test of calculateUsingEssentialMatrix method, of class CameraPose.
     */
    public void testCalculateUsingEssentialMatrix() throws Exception {
        System.out.println("calculateUsingEssentialMatrix");
       
        double[][] k1 = Zhang98Data.getIntrinsicCameraMatrix();
        double[][] k2 = MatrixUtil.copy(k1);
        //x1, x2 size is 3 X 256
        double[][] x1 = Zhang98Data.getObservedFeaturesInImage(1);
        double[][] x2 = Zhang98Data.getObservedFeaturesInImage(5);
                
        Reconstruction.ReconstructionResults result = 
            Reconstruction.calculateUsingEssentialMatrix(k1, k2, x1, x2);
        
        assertNotNull(result);
        
        System.out.printf("\nresult:\nrot1=%strans1=%s\n", 
            FormatArray.toString(result.k1ExtrRot, "%.4e"),
            FormatArray.toString(result.k1ExtrTrans, "%.4e"));
        System.out.printf("\nresult:\nrot2=%strans2=%s\n", 
            FormatArray.toString(result.k2ExtrRot, "%.4e"),
            FormatArray.toString(result.k2ExtrTrans, "%.4e"));
        
        System.out.printf("\nimg1:\nrot=%strans=%s\n", 
                FormatArray.toString(Zhang98Data.getRotation(1), "%.4e"),
                FormatArray.toString(Zhang98Data.getTranslation(1), "%.4e"));
        System.out.printf("\nimg5:\nrot=%strans=%s\n", 
                FormatArray.toString(Zhang98Data.getRotation(5), "%.4e"),
                FormatArray.toString(Zhang98Data.getTranslation(5), "%.4e"));
        
        double[][] diffRSameCenter = Rotation.procrustesAlgorithmForRotation(
            Zhang98Data.getRotation(1), Zhang98Data.getRotation(5));
        
        System.out.printf("\ndifference in rot between img1 and img5=\n%s\n", 
           FormatArray.toString(diffRSameCenter, "%.4e"));
        
        double[] out = new double[3];
        extractThetaFromZYX(diffRSameCenter, out);
        
        System.out.printf("\ndifference in rot between img1 and img5 in euler angles=\n");
        for (double a : out) {
            System.out.printf("%.2f ", 180.*a/Math.PI);
        }
        System.out.println();
        double[] diffTrans = MatrixUtil.subtract(Zhang98Data.getTranslation(1), Zhang98Data.getTranslation(5));
        System.out.printf("\ndifference in trans between img1 and img5=\n%s\n",
            FormatArray.toString(diffTrans, "%.3e"));
        System.out.flush();
    }

    /**
     * Test of calculateReconstruction method, of class Reconstruction.
     */
    public void testCalculateReconstructionWithIntrinsicCamera() throws NotConvergedException {
        
        //test data from:
        //http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example5.html        
        // "Fifth calibration example - Calibrating a stereo system, stereo image rectification and 3D stereo triangulation"
        // by Jean-Yves Bouguet
        // Camera Calibration Toolbox for Matlab
        //
        // left camera: 
        //    focal length = 533.5, 533.5
        //    cc = 341.6, 235.2
        //    skew = 0, 0
        //    radial distortion k = -0.288, 0.097, 0.001, -0.0003, 0
        //
        // right camera: 
        //    focal length = 536.8, 536.5
        //    cc = 326.3, 250.1
        //    skew = 0, 0
        //    radial distortion k = -0.289, 0.107, 0.001, -0.0001, 0
        //
        // rotation vector om=0.00669, 0.00452, -0.0035
        // translation vector t = -99.80198, 1.12443, 0.05041
        //
        // note: the checkerboard sqaures are 30mm in WCS metric
        //
        // note: radial distortion should be corrected:  use on the original coordinates:
        //     x_corrected = x*(1 + k1*r^2 + k2r^4) where r is distance of point from cc.
        
        double[][] k1Intr;
        double[][] k2Intr;
            k1Intr = Camera.createIntrinsicCameraMatrix(533.07, 341.6, 234.3);
            k2Intr = Camera.createIntrinsicCameraMatrix(536.7, 326.5, 249.3);
        
        double[][] k1ExtrRot = MatrixUtil.createIdentityMatrix(3);
        double[] k1ExtrTrans = new double[]{0, 0, 0};
        
        double[][] k2ExtrRot = Rotation.createRodriguesFormulaRotationMatrix(
            new double[]{0.00611, 0.00409, -0.00359});
        double[] k2ExtrTrans = new double[]{-99.85, 0.82, 0.44};
        //double[] k2ExtrTransRev = Arrays.copyOf(k2ExtrTrans, k2ExtrTrans.length);
        //MatrixUtil.multiply(k2ExtrTransRev, -1);
        
        System.out.printf("expected k2ExtrRot from Rodrigues formula\n=%s\n", FormatArray.toString(k2ExtrRot, "%.3e"));
        /*
        [junit] =1.000e+00, 3.577e-03, 4.101e-03 
        [junit] -3.602e-03, 1.000e+00, 6.103e-03 
        [junit] -4.079e-03, -6.117e-03, 1.000e+00 
        */
        
        double[][] x1 = new double[3][10];
        x1[0] = new double[]{129, 160, 140, 226, 225, 232, 341, 407, 532, 527};
        x1[1] = new double[]{145, 319, 361, 391, 61, 284, 289, 122, 48, 302};
        x1[2] = new double[]{1, 1, 1, 1, 1, 1, 1, 1,  1, 1};
        
        double[][] x2 = new double[3][10];
        x2[0] = new double[]{76, 110, 84, 164, 110, 124, 218, 275, 401, 402};
        x2[1] = new double[]{168, 331, 372, 401, 78, 295, 302, 134, 52, 318};
        x2[2] = new double[]{1, 1, 1, 1, 1, 1, 1, 1,  1, 1};
        
        /*
        System.out.printf("results=%s\n\n", rr.toString());
        
  //multiply by focal length?  or K?   see if triangulation in szeliski has advice
        
        System.out.println("XW points normalized by last coordinates:");
        double[] pt = new double[4];
        for (int j = 0; j < rr.XW[0].length; ++j) {    
            for (int i = 0; i < rr.XW.length; ++i) {
                pt[i] = rr.XW[i][j]/rr.XW[3][j];
            }
            System.out.printf("%s\n", FormatArray.toString(pt, "%.3e"));
        }
        */
    }
    
    public void testAffine() throws NotConvergedException, IOException {
    
        double[][] xIn = loadNCBook(false);

        Reconstruction.OrthographicProjectionResults affineProj =
            Reconstruction.calculateAffineReconstruction(xIn, 2);

    }

    public void testParaperspective() throws IOException, NotConvergedException {
        double[][] xIn = loadNCBook(false);

        Reconstruction.ParaperspectiveProjectionResults paraProj
                = Reconstruction.calculateParaperspectiveReconstruction(xIn, 2);
    }
    
    private double[][] loadNCBook(boolean rotateBy90) throws IOException, NotConvergedException {
        
        int maxDimension = 256;//512;

        ImageProcessor imageProcessor = new ImageProcessor();
        ImageSegmentation imageSegmentation = new ImageSegmentation();

        String[][] filePairs = new String[1][];
        filePairs[0] = new String[]{
            "nc_book_01.png",
            "nc_book_02.png"};
        
        boolean binImages = false;//true;
        
        int rotate = rotateBy90 ? 1 : 0;
        
        int i, ii, np;
            
        String lbl = "_";
        if (rotate == 1) {
            lbl = "rot90_";
        }

        String fileName1 = filePairs[0][0];
        String fileName2 = filePairs[0][1];

        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        int w1 = img1.getWidth();
        int h1 = img1.getHeight();
        if (binImages) {
            int binFactor1 = (int) Math.ceil(Math.max(
                    (float) w1 / maxDimension,
                    (float) h1 / maxDimension));
            img1 = imageProcessor.binImage(img1, binFactor1);
            //MiscDebug.writeImage(img1, "_"  + fileName1Root);
        }
        GreyscaleImage img1GS = img1.copyToGreyscale2();

        idx = fileName2.lastIndexOf(".");
        String fileName2Root = fileName2.substring(0, idx);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
        int w2 = img2.getWidth();
        int h2 = img2.getHeight();
        if (binImages) {
            int binFactor2 = (int) Math.ceil(Math.max(
                    (float) w2 / maxDimension,
                    (float) h2 / maxDimension));
            img2 = imageProcessor.binImage(img2, binFactor2);
        }
        GreyscaleImage img2GS = img2.copyToGreyscale2();

        Transformer tr = new Transformer();
        TransformationParameters params = new TransformationParameters();
        TransformationParameters paramsRev = new TransformationParameters();

        if (rotate == 1) {
            params.setTranslationY(img2GS.getWidth());
            params.setRotationInDegrees(90);

            img2GS = tr.applyTransformation(img2GS, params,
                    img2GS.getHeight(), img2GS.getWidth());

            paramsRev.setTranslationX(img2GS.getHeight());
            paramsRev.setRotationInDegrees(-90);
        }

        int x, y;

        if (binImages) {
            np = 300;
        } else {
            np = 600;
        }

        ORB orb1 = new ORB(img1GS, np);
        orb1.detectAndExtract();
        //orb1.overrideToAlsoCreate1stDerivKeypoints();
        //orb1.overrideToNotCreateATrousKeypoints();

        ORB orb2 = new ORB(img2GS, np);
        orb2.detectAndExtract();
        //orb2.overrideToAlsoCreate1stDerivKeypoints();
        //orb2.overrideToNotCreateATrousKeypoints();

        ORB.Descriptors d1 = orb1.getAllDescriptors();
        ORB.Descriptors d2 = orb2.getAllDescriptors();
        double[][] xKP1 = orb1.getAllKeyPointsHomogenous();
        double[][] xKP2 = orb2.getAllKeyPointsHomogenous();

        int nKP1 = xKP1[0].length;
        int nKP2 = xKP2[0].length;

        int x1, x2, y1, y2;
        Image tmp1 = img1GS.copyToColorGreyscale();
        for (i = 0; i < nKP1; ++i) {
            y = (int) xKP1[1][i];
            x = (int) xKP1[0][i];
            ImageIOHelper.addPointToImage(x, y, tmp1, 2, 255, 0, 0);
        }
        //MiscDebug.writeImage(tmp1, "_kp_gs_" + lbl + fileName1Root);

        Image tmp2 = img2GS.copyToColorGreyscale();
        for (i = 0; i < nKP2; ++i) {
            y = (int) xKP2[1][i];
            x = (int) xKP2[0][i];
            ImageIOHelper.addPointToImage(x, y, tmp2, 2, 255, 0, 0);
        }
        MiscDebug.writeImage(tmp1, "_kp_gs_" + lbl + fileName1Root);
        MiscDebug.writeImage(tmp2, "_kp_gs_" + lbl + fileName2Root);

        // kp1 and kp2 are in col major format (col, row) and we want normalized points, so converting to double[][]
        double[][] xKP1n = MatrixUtil.copy(xKP1);
        double[][] xKP2n = MatrixUtil.copy(xKP2);
        double[][] t1 = EpipolarNormalizationHelper.unitStandardNormalize(xKP1n);
        double[][] t2 = EpipolarNormalizationHelper.unitStandardNormalize(xKP2n);

        int[][] matchedIdxs = null;
        if (true /*normalize points*/) {
            // col, row
            matchedIdxs = ORBMatcher.matchDescriptors(d1, d2, xKP1n, xKP2n);
            if (matchedIdxs != null) {
                System.out.printf("# matched after normalization = %d\n", matchedIdxs.length);
            }
        } else {
            matchedIdxs = ORBMatcher.matchDescriptors(d1, d2, xKP1, xKP2);
        }

        if (matchedIdxs == null) {
            return null;
        }
        int idx1, idx2;
        System.out.println(lbl + fileName1Root + " matched=" + matchedIdxs.length);
        CorrespondencePlotter plotter = new CorrespondencePlotter(tmp1, tmp2);
        for (i = 0; i < matchedIdxs.length; ++i) {
            idx1 = matchedIdxs[i][0];
            idx2 = matchedIdxs[i][1];
            x1 = (int) xKP1[0][idx1];
            y1 = (int) xKP1[1][idx1];
            x2 = (int) xKP2[0][idx2];
            y2 = (int) xKP2[1][idx2];
            plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 1);
        }
        plotter.writeImage("_r_corres_orb_gs_" + lbl + fileName1Root);

        /*
        x the image coordinates of feature correspondences in 2 or more
         * images.  format is 2 X (nImages * nFeatures) where row 0 holds the x-coordinates
         * and row 1 holds the y-coordinates and each image's features are given
         * before the next and the features are ordered in the same manner within
         * all images.
         * for example: row 0 = [img_0_feature_0, ... img_0_feature_n-1, ... img_m-1_feature_0,...
         *     img_m-1_feature_n-1].
         */
        double[][] xIn = MatrixUtil.zeros(2, 2*matchedIdxs.length);
        for (i = 0; i < matchedIdxs.length; ++i) {
            idx1 = matchedIdxs[i][0];
            xIn[0][i] = (int) xKP1[0][idx1];
            xIn[1][i] = (int) xKP1[1][idx1];
        }
        for (i = 0; i < matchedIdxs.length; ++i) {
            idx2 = matchedIdxs[i][1];
            xIn[0][i + matchedIdxs.length] = (int) xKP2[0][idx2];
            xIn[1][i + matchedIdxs.length] = (int) xKP2[1][idx2];
        }

        return xIn;
    }

    /**
     * given a list of keypoints in format (col, row), write a double array of size [2 X kP.size()] where the
     * first row holds the col values, the second row holds the row values and the third row holds values "1" for
     * the z-axis representation in homogeneous coordinates.
     * @param kP input list of keypoints in format (col, row)
     * @return the data in format [2 X kP.size()]
     */
    private double[][] convertToArray(List<PairInt> kP) {
        int n = kP.size();
        double[][] x = MatrixUtil.zeros(3, n);
        for (int i = 0; i < n; ++i) {
            x[0][i] = kP.get(i).getX();
            x[1][i] = kP.get(i).getY();
        }
        Arrays.fill(x[2], 1);
        return x;
    }

}
