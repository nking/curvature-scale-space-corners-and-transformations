package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.*;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.imageProcessing.features.orb.ORB.Descriptors;
import algorithms.imageProcessing.features.orb.ORB2;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.PairInt;
import algorithms.util.QuadInt;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;

import java.util.List;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class ORBMatcher2Test extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public ORBMatcher2Test() {
    }

    public void test0() throws Exception {

        int maxDimension = 256;//512;

        ImageProcessor imageProcessor = new ImageProcessor();
        ImageSegmentation imageSegmentation = new ImageSegmentation();

        String[][] filePairs = new String[6][];
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
        
        boolean binImages = true;
        
        int i, ii, np;
        //for (int rotate = 0; rotate < 2; ++rotate) {
        for (int rotate = 0; rotate < 1; ++rotate) {
            
            String lbl = "_";
            if (rotate == 1) {
                lbl = "rot90_";
            }
 
            for (ii = 0; ii < filePairs.length; ++ii) {
            //for (ii = 5; ii < 6; ++ii) {

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

                ORB2 orb1 = new ORB2(img1GS, np);
                orb1.detectAndExtract();
                //orb1.overrideToAlsoCreate1stDerivKeypoints();
                //orb1.overrideToNotCreateATrousKeypoints();

                ORB2 orb2 = new ORB2(img2GS, np);
                orb2.detectAndExtract();
                //orb1.overrideToAlsoCreate1stDerivKeypoints();
                //orb2.overrideToNotCreateATrousKeypoints();

                Descriptors d1 = orb1.getAllDescriptors();
                Descriptors d2 = orb2.getAllDescriptors();
                List<PairInt> kp1 = orb1.getAllKeyPointsRC();
                List<PairInt> kp2 = orb2.getAllKeyPointsRC();
                
                QuadInt[] matched = ORBMatcher.matchDescriptorsRC(d1, d2, kp1, kp2);

                int x1, x2, y1, y2;
                int col, row;
                {//DEBUG
                    Image tmp1 = img1GS.copyToColorGreyscale();
                    for (i = 0; i < kp1.size(); ++i) {
                        col = kp1.get(i).getY();
                        row = kp1.get(i).getX();
                        ImageIOHelper.addPointToImage(col, row, tmp1, 2, 255, 0, 0);
                        /*if (x == 213 && y == 58) {
                            ImageIOHelper.addPointToImage(x, y, tmp1, 3, 255, 255, 0);
                        }*/
                    }
                    //MiscDebug.writeImage(tmp1, "_kp_gs_" + lbl + fileName1Root);

                    Image tmp2 = img2GS.copyToColorGreyscale();
                    for (i = 0; i < kp2.size(); ++i) {
                        col = kp2.get(i).getY();
                        row = kp2.get(i).getX();
                        ImageIOHelper.addPointToImage(col, row, tmp2, 2, 255, 0, 0);
                        /*if (x == 178 && y == 177) {
                            ImageIOHelper.addPointToImage(x, y, tmp2, 3, 255, 255, 0);
                        }*/
                    }
                    //MiscDebug.writeImage(tmp2, "_kp_gs_" + lbl + fileName2Root);
                    System.out.println(lbl + fileName1Root + " matched=" + matched.length);
                    CorrespondencePlotter plotter = 
                        new CorrespondencePlotter(tmp1, tmp2);
                    for (i = 0; i < matched.length; ++i) {
                        y1 = matched[i].getA();
                        x1 = matched[i].getB();
                        y2 = matched[i].getC();
                        x2 = matched[i].getD();
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
    
}