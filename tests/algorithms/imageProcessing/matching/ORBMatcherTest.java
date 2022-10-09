package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.features.HOGs;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.imageProcessing.features.orb.ORB.Descriptors;
import algorithms.imageProcessing.transform.EpipolarNormalizationHelper;
import algorithms.imageProcessing.transform.Reconstruction;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.matrix.MatrixUtil;
import algorithms.matrix.MatrixUtil.SVDProducts;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.FormatArray;
import algorithms.util.PairInt;
import algorithms.util.QuadInt;
import algorithms.util.ResourceFinder;
import gnu.trove.list.TDoubleList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

/**
 *
 * @author nichole
 */
public class ORBMatcherTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public ORBMatcherTest() {
    }

    public void test0() throws Exception {

        int maxDimension = 256;//512;

        ImageProcessor imageProcessor = new ImageProcessor();
        ImageSegmentation imageSegmentation = new ImageSegmentation();

        String[][] filePairs = new String[7][];
        filePairs[0] = new String[]{
            "venturi_mountain_j6_0001.png",
            "venturi_mountain_j6_0010.png"};
        filePairs[1] = new String[]{
            "campus_010.jpg",
            "campus_011.jpg"};
        filePairs[2] = new String[]{
            "merton_college_I_001.jpg",
            "merton_college_I_002.jpg"};
         
        filePairs[3] = new String[]{
            "books_illum3_v0_695x555.png",
            "books_illum3_v6_695x555.png"};
                    
        filePairs[4] = new String[]{
            "checkerboard_01.jpg",
            "checkerboard_02.jpg"};
        
        filePairs[5] = new String[]{
            "nc_book_01.png",
            "nc_book_02.png"};

        filePairs[6] = new String[]{
                "merton_college_I_001.jpg",
                "merton_college_I_002.jpg"};
        
        boolean binImages = true;
        
        int i, ii, np;
        for (int rotate = 0; rotate < 2; ++rotate) {
        //for (int rotate = 0; rotate < 1; ++rotate) {
            
            String lbl = "_";
            if (rotate == 1) {
                lbl = "rot90_";
            }
 
            for (ii = 0; ii < filePairs.length; ++ii) {

                String fileName1 = filePairs[ii][0];
                String fileName2 = filePairs[ii][1];

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

                Descriptors d1 = orb1.getAllDescriptors();
                Descriptors d2 = orb2.getAllDescriptors();
                double[][] xKP1 = orb1.getAllKeyPointsHomogenous();
                double[][] xKP2 = orb2.getAllKeyPointsHomogenous();

                int nKP1 = xKP1[0].length;
                int nKP2 = xKP2[0].length;

                int x1, x2, y1, y2;
                Image tmp1 = img1GS.copyToColorGreyscale();
                for (i = 0; i < nKP1; ++i) {
                    y = (int)xKP1[1][i];
                    x = (int)xKP1[0][i];
                    ImageIOHelper.addPointToImage(x, y, tmp1, 2, 255, 0, 0);
                }
                //MiscDebug.writeImage(tmp1, "_kp_gs_" + lbl + fileName1Root);

                Image tmp2 = img2GS.copyToColorGreyscale();
                for (i = 0; i < nKP2; ++i) {
                    y = (int)xKP2[1][i];
                    x = (int)xKP2[0][i];
                    ImageIOHelper.addPointToImage(x, y, tmp2, 2, 255, 0, 0);
                }

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
                        System.out.printf("%d) # matched after normalization = %d\n", ii, matchedIdxs.length);
                    }
                } else {
                    matchedIdxs = ORBMatcher.matchDescriptors(d1, d2, xKP1, xKP2);
                }

                if (matchedIdxs == null) {
                    continue;
                }

                //MiscDebug.writeImage(tmp2, "_kp_gs_" + lbl + fileName2Root);
                int idx1, idx2;
                System.out.println(lbl + fileName1Root + " matched=" + matchedIdxs.length);
                CorrespondencePlotter plotter = new CorrespondencePlotter(tmp1, tmp2);
                for (i = 0; i < matchedIdxs.length; ++i) {
                    idx1 = matchedIdxs[i][0];
                    idx2 = matchedIdxs[i][1];
                    x1 = (int)xKP1[0][idx1];
                    y1 = (int)xKP1[1][idx1];
                    x2 = (int)xKP2[0][idx2];
                    y2 = (int)xKP2[1][idx2];
                    plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 1);
                }
                plotter.writeImage("_corres_orb_gs_" + lbl + fileName1Root);

                /*NOTE:
                   Further considering algorithm in paper on outlier correction:
                       https://www.researchgate.net/publication/221110532_Outlier_Correction_in_Image_Sequences_for_the_Affine_Camera
                   "Outlier Correction in Image Sequences for the Affine Camera"
                      by Huynh, Hartley, and Heydeon 2003
                      Proceedings of the Ninth IEEE International Conference on Computer Vision (ICCV’03)
                */
            }
        }
    }

}
