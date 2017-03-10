package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import gnu.trove.list.TIntList;
import java.util.logging.Logger;
import junit.framework.TestCase;

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

        String[][] filePairs = new String[5][];
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

        //for (int i = 0; i < filePairs.length; ++i) {
        for (int i = 3; i < 4; ++i) {

            String fileName1 = filePairs[i][0];

            String fileName2 = filePairs[i][1];
           
            int idx = fileName1.lastIndexOf(".");
            String fileName1Root = fileName1.substring(0, idx);
            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
            int w1 = img1.getWidth();
            int h1 = img1.getHeight();
            int binFactor1 = (int) Math.ceil(Math.max(
                (float) w1 / maxDimension,
                (float) h1 / maxDimension));
            img1 = imageProcessor.binImage(img1, binFactor1);
            //MiscDebug.writeImage(img1, "_"  + fileName1Root);
            GreyscaleImage img1GS = img1.copyToGreyscale2();
            //imageProcessor.blur(img1GS, SIGMA.ONE);
            
            idx = fileName2.lastIndexOf(".");
            String fileName2Root = fileName2.substring(0, idx);
            String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
            ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
            int w2 = img2.getWidth();
            int h2 = img2.getHeight();
            int binFactor2 = (int) Math.ceil(Math.max(
                (float) w2 / maxDimension,
                (float) h2 / maxDimension));
            img2 = imageProcessor.binImage(img2, binFactor2);
            //MiscDebug.writeImage(img2, "_"  + fileName2Root);
            GreyscaleImage img2GS = img2.copyToGreyscale2();
            //imageProcessor.blur(img2GS, SIGMA.ONE);
            
            GreyscaleImage[] lch1 = imageProcessor.createLCHForLUV(img1);
            GreyscaleImage[] lch2 = imageProcessor.createLCHForLUV(img2);

            int np = 100;       
            
            ORB orb1 = new ORB(img1GS, np);
            //orb1.overrideToAlsoCreate1stDerivKeypoints();
            //orb1.overrideToCreateCurvaturePoints();
            orb1.overrideToNotCreateATrousKeypoints();
            //orb1.overrideToUseSingleScale();
            orb1.detectAndExtract();

            ORB orb2 = new ORB(img2GS, np);
            //orb2.overrideToAlsoCreate1stDerivKeypoints();
            //orb2.overrideToCreateCurvaturePoints();
            orb2.overrideToNotCreateATrousKeypoints();
            //orb2.overrideToUseSingleScale();
            orb2.detectAndExtract();
            
            {//DEBUG
                Image tmp1 = img1.copyToGreyscale2().copyToColorGreyscale();
                for (int ii = 0; ii < 1; ++ii) {
                    TIntList pixIdxs = orb1.getKeyPointListPix(ii);
                    ImageIOHelper.addCurveToImage(pixIdxs, tmp1, 1, 255, 0, 0);
                }
                MiscDebug.writeImage(tmp1, "_kp_" + fileName1Root);
                
                Image tmp2 = img2.copyToGreyscale2().copyToColorGreyscale();
                for (int ii = 0; ii < 1; ++ii) {
                    TIntList pixIdxs = orb2.getKeyPointListPix(ii);
                    ImageIOHelper.addCurveToImage(pixIdxs, tmp2, 1, 255, 0, 0);
                }
                MiscDebug.writeImage(tmp2, "_kp_" + fileName2Root);
            }
            
            
            orb1 = new ORB(lch1[1], 500);
            //orb1.overrideToAlsoCreate1stDerivKeypoints();
            //orb1.overrideToCreateCurvaturePoints();
            orb1.overrideToUseSingleScale();
            orb1.detectAndExtract();

            orb2 = new ORB(lch2[1], 500);
            //orb2.overrideToAlsoCreate1stDerivKeypoints();
            //orb2.overrideToCreateCurvaturePoints();
            orb2.overrideToUseSingleScale();
            orb2.detectAndExtract();
            
            {//DEBUG
                Image tmp1 = lch1[1].copyToColorGreyscale();
                for (int ii = 0; ii < 1; ++ii) {
                    TIntList pixIdxs = orb1.getKeyPointListPix(ii);
                    ImageIOHelper.addCurveToImage(pixIdxs, tmp1, 1, 255, 0, 0);
                }
                MiscDebug.writeImage(tmp1, "_kp_C_" + fileName1Root);
                
                Image tmp2 = lch2[1].copyToColorGreyscale();
                for (int ii = 0; ii < 1; ++ii) {
                    TIntList pixIdxs = orb2.getKeyPointListPix(ii);
                    ImageIOHelper.addCurveToImage(pixIdxs, tmp2, 1, 255, 0, 0);
                }
                MiscDebug.writeImage(tmp2, "_kp_C_" + fileName2Root);
            }
            
        }
    }

}
