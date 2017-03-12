package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.features.HOGs;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.imageProcessing.features.orb.ORB.Descriptors;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.MiscDebug;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.PairInt;
import algorithms.util.QuadInt;
import algorithms.util.ResourceFinder;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import java.util.ArrayList;
import java.util.List;
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

        for (int i = 0; i < filePairs.length; ++i) {
        //for (int i = 2; i < 3; ++i) {

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
            GreyscaleImage img2GS = img2.copyToGreyscale2();
            //MiscDebug.writeImage(img2, "_"  + fileName2Root);
            
            
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
            
            Descriptors d1 = orb1.getDescriptorsList().get(0);
            Descriptors d2 = orb2.getDescriptorsList().get(0);
            List<PairInt> kp1 = orb1.getKeyPointListColMaj(0);
            List<PairInt> kp2 = orb2.getKeyPointListColMaj(0);
            QuadInt[] matched;
            
            if (i == 2) {
                
                HOGs hogs1 = new HOGs(img1GS, 1, 6);
                HOGs hogs2 = new HOGs(img2GS, 1, 6);
                TDoubleList orientations1 = orb1.getOrientationsList().get(0);
                TDoubleList orientations2 = orb2.getOrientationsList().get(0);
                
                matched = ORBMatcher.matchDescriptors(d1, d2, kp1, kp2,
                    hogs1, hogs2, orientations1, orientations2);
            
            } else {
            
                matched = ORBMatcher.matchDescriptors(
                    d1, d2, kp1, kp2);
            }
            
            {//DEBUG
                Image tmp1 = img1.copyToGreyscale2().copyToColorGreyscale();
                for (int ii = 0; ii < 1; ++ii) {
                    TIntList pixIdxs = orb1.getKeyPointListPix(ii);
                    ImageIOHelper.addCurveToImage(pixIdxs, tmp1, 1, 255, 0, 0);
                }
                MiscDebug.writeImage(tmp1, "_kp_" + fileName1Root);
                
                Image tmp2 = img2GS.copyToColorGreyscale();
                for (int ii = 0; ii < 1; ++ii) {
                    TIntList pixIdxs = orb2.getKeyPointListPix(ii);
                    ImageIOHelper.addCurveToImage(pixIdxs, tmp2, 1, 255, 0, 0);
                }
                MiscDebug.writeImage(tmp2, "_kp_" + fileName2Root);
                CorrespondencePlotter plotter = 
                    new CorrespondencePlotter(tmp1, tmp2);
                for (int ii = 0; ii < matched.length; ++ii) {
                    int x1 = matched[ii].getA();
                    int y1 = matched[ii].getB();
                    int x2 = matched[ii].getC();
                    int y2 = matched[ii].getD();
                    plotter.drawLineInAlternatingColors(x1, y1, 
                        x2, y2, 1);
                }
                plotter.writeImage("_corres_" + fileName1Root);
            }
            
            
            ORB orb1c = new ORB(lch1[1], np);
            //orb1c.overrideToAlsoCreate1stDerivKeypoints();
            //orb1c.overrideToCreateCurvaturePoints();
            orb1c.overrideToUseSingleScale();
            orb1c.overrideToNotCreateATrousKeypoints();
            orb1c.detectAndExtract();

            ORB orb2c = new ORB(lch2[1], np);
            //orb2c.overrideToAlsoCreate1stDerivKeypoints();
            //orb2c.overrideToCreateCurvaturePoints();
            orb2c.overrideToUseSingleScale();
            orb2c.overrideToNotCreateATrousKeypoints();
            orb2c.detectAndExtract();
            
            Descriptors d1c = orb1c.getDescriptorsList().get(0);
            Descriptors d2c = orb2c.getDescriptorsList().get(0);
            List<PairInt> kp1c = orb1c.getKeyPointListColMaj(0);
            List<PairInt> kp2c = orb2c.getKeyPointListColMaj(0);
            if (i == 2) {
                
                HOGs hogs1 = new HOGs(lch1[1], 1, 6);
                HOGs hogs2 = new HOGs(lch2[1], 1, 6);
                TDoubleList orientations1 = orb1.getOrientationsList().get(0);
                TDoubleList orientations2 = orb2.getOrientationsList().get(0);
                
                matched = ORBMatcher.matchDescriptors(d1c, d2c, kp1c, kp2c,
                    hogs1, hogs2, orientations1, orientations2);
            
            } else {
            
                matched = ORBMatcher.matchDescriptors(
                    d1c, d2c, kp1c, kp2c);
            }
            
            {//DEBUG
                Image tmp1 = lch1[1].copyToColorGreyscale();
                for (int ii = 0; ii < tmp1.getNPixels(); ++ii) {
                    tmp1.setRGB(ii, tmp1.getR(ii) + 100, 
                        tmp1.getG(ii) + 100, tmp1.getB(ii) + 100);
                }
                for (int ii = 0; ii < 1; ++ii) {
                    TIntList pixIdxs = orb1.getKeyPointListPix(ii);
                    ImageIOHelper.addCurveToImage(pixIdxs, tmp1, 1, 255, 0, 0);
                }
                MiscDebug.writeImage(tmp1, "_kp_C_" + fileName1Root);
                
                Image tmp2 = lch2[1].copyToColorGreyscale();
                for (int ii = 0; ii < tmp2.getNPixels(); ++ii) {
                    tmp2.setRGB(ii, tmp2.getR(ii) + 100, 
                        tmp2.getG(ii) + 100, tmp2.getB(ii) + 100);
                }
                for (int ii = 0; ii < 1; ++ii) {
                    TIntList pixIdxs = orb2.getKeyPointListPix(ii);
                    ImageIOHelper.addCurveToImage(pixIdxs, tmp2, 1, 255, 0, 0);
                }
                MiscDebug.writeImage(tmp2, "_kp_C_" + fileName2Root);
                CorrespondencePlotter plotter = 
                    new CorrespondencePlotter(tmp1, tmp2);
                for (int ii = 0; ii < matched.length; ++ii) {
                    int x1 = matched[ii].getA();
                    int y1 = matched[ii].getB();
                    int x2 = matched[ii].getC();
                    int y2 = matched[ii].getD();
                    plotter.drawLineInAlternatingColors(x1, y1, 
                        x2, y2, 1);
                }
                plotter.writeImage("_corres_C_" + fileName1Root);
            }
            
            // combining the descriptors
            List<Descriptors> d1CombList = new ArrayList<Descriptors>();
            d1CombList.add(d1);
            d1CombList.add(d1c);
            List<Descriptors> d2CombList = new ArrayList<Descriptors>();
            d2CombList.add(d2);
            d2CombList.add(d2c);
            
            Descriptors d1Comb = ORB.combineDescriptors(d1CombList);
            Descriptors d2Comb = ORB.combineDescriptors(d2CombList);
            kp1.addAll(kp1c);
            kp2.addAll(kp2c);
            matched = ORBMatcher.matchDescriptors(d1Comb, d2Comb, kp1, kp2);
            
            {//DEBUG
                Image tmp1 = img1GS.copyToColorGreyscale();
                ImageIOHelper.addCurveToImage(kp1, tmp1, 1, 255, 0, 0);
                MiscDebug.writeImage(tmp1, "_kp_comb_" + fileName1Root);
                
                Image tmp2 = img2GS.copyToColorGreyscale();
                ImageIOHelper.addCurveToImage(kp2, tmp2, 1, 255, 0, 0);
                MiscDebug.writeImage(tmp2, "_kp_comb_" + fileName2Root);
                CorrespondencePlotter plotter = 
                    new CorrespondencePlotter(tmp1, tmp2);
                for (int ii = 0; ii < matched.length; ++ii) {
                    int x1 = matched[ii].getA();
                    int y1 = matched[ii].getB();
                    int x2 = matched[ii].getC();
                    int y2 = matched[ii].getD();
                    plotter.drawLineInAlternatingColors(x1, y1, 
                        x2, y2, 1);
                }
                plotter.writeImage("_corres_comb_" + fileName1Root);
            }
        }
    }
    
    public void test90() throws Exception {

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

        for (int i = 0; i < filePairs.length; ++i) {
        //for (int i = 2; i < 3; ++i) {

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
            
            Transformer tr = new Transformer();
            TransformationParameters params = new TransformationParameters();
            params.setTranslationY(img2.getWidth());
            params.setRotationInDegrees(90);
            
            GreyscaleImage img2GS = img2.copyToGreyscale2();
            img2GS = tr.applyTransformation(img2GS, params, 
                img2GS.getHeight(), img2GS.getWidth());
            
            TransformationParameters paramsRev = new TransformationParameters();
            paramsRev.setTranslationX(img2GS.getHeight());
            paramsRev.setRotationInDegrees(-90);
            
            MiscDebug.writeImage(img2GS, "_img2rot90_"  + fileName2Root);
            
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
            
            Descriptors d1 = orb1.getDescriptorsList().get(0);
            Descriptors d2 = orb2.getDescriptorsList().get(0);
            List<PairInt> kp1 = orb1.getKeyPointListColMaj(0);
            List<PairInt> kp2 = orb2.getKeyPointListColMaj(0);
            
            ORB orb1c = new ORB(lch1[1], np);
            //orb1c.overrideToAlsoCreate1stDerivKeypoints();
            //orb1c.overrideToCreateCurvaturePoints();
            orb1c.overrideToUseSingleScale();
            orb1c.overrideToNotCreateATrousKeypoints();
            orb1c.detectAndExtract();

            lch2[1] = tr.applyTransformation(lch2[1], params, 
                lch2[1].getHeight(), lch2[1].getWidth());
            
            ORB orb2c = new ORB(lch2[1], np);
            //orb2c.overrideToAlsoCreate1stDerivKeypoints();
            //orb2c.overrideToCreateCurvaturePoints();
            orb2c.overrideToUseSingleScale();
            orb2c.overrideToNotCreateATrousKeypoints();
            orb2c.detectAndExtract();
            
            Descriptors d1c = orb1c.getDescriptorsList().get(0);
            Descriptors d2c = orb2c.getDescriptorsList().get(0);
            List<PairInt> kp1c = orb1c.getKeyPointListColMaj(0);
            List<PairInt> kp2c = orb2c.getKeyPointListColMaj(0);
            
            
            
            // combining the descriptors
            List<Descriptors> d1CombList = new ArrayList<Descriptors>();
            d1CombList.add(d1);
            d1CombList.add(d1c);
            List<Descriptors> d2CombList = new ArrayList<Descriptors>();
            d2CombList.add(d2);
            d2CombList.add(d2c);
            
            Descriptors d1Comb = ORB.combineDescriptors(d1CombList);
            Descriptors d2Comb = ORB.combineDescriptors(d2CombList);
            List<PairInt> kp1Comb = new ArrayList<PairInt>(kp1);
            kp1Comb.addAll(kp1c);
            List<PairInt> kp2Comb = new ArrayList<PairInt>(kp2);
            kp2Comb.addAll(kp2c);
            
            
            
            QuadInt[] matched, matchedC, matchedComb;
            
            if (true || i == 2) {
                
                HOGs hogs1 = new HOGs(img1GS, 1, 6);
                HOGs hogs2 = new HOGs(img2GS, 1, 6);
                TDoubleList orientations1 = orb1.getOrientationsList().get(0);
                TDoubleList orientations2 = orb2.getOrientationsList().get(0);
                
                matched = ORBMatcher.matchDescriptors(d1, d2, kp1, kp2,
                    hogs1, hogs2, orientations1, orientations2);
                
                hogs1 = new HOGs(lch1[1], 1, 6);
                hogs2 = new HOGs(lch2[1], 1, 6);
                orientations1 = orb1c.getOrientationsList().get(0);
                orientations2 = orb2c.getOrientationsList().get(0);
                
                matchedC = ORBMatcher.matchDescriptors(d1c, d2c, kp1c, kp2c,
                    hogs1, hogs2, orientations1, orientations2);
            
            } else {
            
                matched = ORBMatcher.matchDescriptors(
                    d1, d2, kp1, kp2);
                
                matchedC = ORBMatcher.matchDescriptors(
                    d1c, d2c, kp1c, kp2c);
            }
            
            matchedComb = ORBMatcher.matchDescriptors(d1Comb, d2Comb, kp1Comb, 
                kp2Comb);
            
            // reverse transforming the dataset2 images and points to
            // more easily see the results
            img2GS = tr.applyTransformation(img2GS, paramsRev, img2GS.getHeight(), 
                img2GS.getWidth());
            lch2[1] = tr.applyTransformation(lch2[1], paramsRev,
                lch2[1].getHeight(), lch2[1].getWidth());
            tr.applyTransformation(paramsRev, kp2);            
            tr.applyTransformation(paramsRev, kp2c);            
            tr.applyTransformation(paramsRev, kp2Comb);
               
            
            {//DEBUG
                Image tmp1 = img1.copyToGreyscale2().copyToColorGreyscale();
                ImageIOHelper.addCurveToImage(kp1, tmp1, 1, 255, 0, 0);
                
                //MiscDebug.writeImage(tmp1, "_kp_" + fileName1Root);
                
                Image tmp2 = img2GS.copyToColorGreyscale();
                ImageIOHelper.addCurveToImage(kp2, tmp2, 1, 255, 0, 0);
                
                //MiscDebug.writeImage(tmp2, "_kp_" + fileName2Root);
                
                CorrespondencePlotter plotter = 
                    new CorrespondencePlotter(tmp1, tmp2);
                for (int ii = 0; ii < matched.length; ++ii) {
                    int x1 = matched[ii].getA();
                    int y1 = matched[ii].getB();
                    
                    int x2 = matched[ii].getC();
                    int y2 = matched[ii].getD();
                    double[] xy2Tr = tr.applyTransformation(paramsRev, x2, y2);
                    
                    plotter.drawLineInAlternatingColors(x1, y1, 
                        (int)xy2Tr[0], (int)xy2Tr[1], 
                        1);
                }
                plotter.writeImage("_corres_rot90_" + fileName1Root);
            }
            
            {//DEBUG
                Image tmp1 = lch1[1].copyToColorGreyscale();
                for (int ii = 0; ii < tmp1.getNPixels(); ++ii) {
                    tmp1.setRGB(ii, tmp1.getR(ii) + 100, 
                        tmp1.getG(ii) + 100, tmp1.getB(ii) + 100);
                }
                ImageIOHelper.addCurveToImage(kp1c, tmp1, 1, 255, 0, 0);
                
                Image tmp2 = lch2[1].copyToColorGreyscale();
                for (int ii = 0; ii < tmp2.getNPixels(); ++ii) {
                    tmp2.setRGB(ii, tmp2.getR(ii) + 100, 
                        tmp2.getG(ii) + 100, tmp2.getB(ii) + 100);
                }
                
                ImageIOHelper.addCurveToImage(kp2c, tmp2, 1, 255, 0, 0);
                
                CorrespondencePlotter plotter = 
                    new CorrespondencePlotter(tmp1, tmp2);
                for (int ii = 0; ii < matched.length; ++ii) {
                    int x1 = matched[ii].getA();
                    int y1 = matched[ii].getB();
                    int x2 = matched[ii].getC();
                    int y2 = matched[ii].getD();
                    double[] xy2Tr = tr.applyTransformation(paramsRev, x2, y2);
                    
                    plotter.drawLineInAlternatingColors(x1, y1, 
                        (int)xy2Tr[0], (int)xy2Tr[1], 
                        1);
                }
                plotter.writeImage("_corres_C_rot90_" + fileName1Root);
            }
            
            
            {//DEBUG
                Image tmp1 = img1GS.copyToColorGreyscale();
                ImageIOHelper.addCurveToImage(kp1, tmp1, 1, 255, 0, 0);
                
                Image tmp2 = img2GS.copyToColorGreyscale();
                ImageIOHelper.addCurveToImage(kp2Comb, tmp2, 1, 255, 0, 0);
                CorrespondencePlotter plotter = 
                    new CorrespondencePlotter(tmp1, tmp2);
                for (int ii = 0; ii < matched.length; ++ii) {
                    int x1 = matched[ii].getA();
                    int y1 = matched[ii].getB();
                    int x2 = matched[ii].getC();
                    int y2 = matched[ii].getD();
                    double[] xy2Tr = tr.applyTransformation(paramsRev, x2, y2);
                    
                    plotter.drawLineInAlternatingColors(x1, y1, 
                        (int)xy2Tr[0], (int)xy2Tr[1], 
                        1);
                }
                plotter.writeImage("_corres_comb_rot90_" + fileName1Root);
            }
        }
    }

}
