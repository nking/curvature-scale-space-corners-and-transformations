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
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.MiscDebug;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.PairInt;
import algorithms.util.QuadInt;
import algorithms.util.ResourceFinder;
import gnu.trove.list.TDoubleList;
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
        
        int i, ii, np;
        for (int rotate = 0; rotate < 1/*2*/; ++rotate) {
            
            String lbl = "_";
            if (rotate == 1) {
                lbl = "rot90_";
            }

            //for (ii = 0; ii < filePairs.length; ++ii) {
            for (ii = 5; ii < 6; ++ii) {

                String fileName1 = filePairs[ii][0];
                String fileName2 = filePairs[ii][1];

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


                GreyscaleImage[] lch1 = imageProcessor.createLCHForLUV(img1);
                GreyscaleImage[] lch2 = imageProcessor.createLCHForLUV(img2);

                // Two ways to compensate for the discontinuity between values 
                // 255 and 0 on these images whose values are cylindrical coordinates.
                // (1) shift copies of images by values of 40-ish and if value is
                //     larger than 255, subtract 255 to wrap around the discontinuity.
                //     Then use the feature extraction algorithms on shifted and unshifted
                //     images and combine results such as keypoints.  Also need to
                //     filter out redundant results.
                // (2) find a best shift in values for an image:
                // look at the histogram stats of all pixels in an image to find the largest
                // gap of at least 20-ish in terms of the bin axis.
                // e.g. if bin size=5, the argmin bin of the histogram, or
                // if there are more than one empty bins then the center of the
                // largest consecutive number of empty bins.
                // the argmin of bins is then used as a shift to apply to the
                // image.
                //    ideally, the 2 images to compare would have the same shift
                //    applied to them, but most of the comparison methods used
                //    account for the "luminosity" difference so the different
                //    shifts usually would not affect the result, but that should
                //    be checked before using...

                // using the 40 degree shift for now

                int thetaShift = 40;
                GreyscaleImage h1Shifted = imageProcessor.copyAndShiftPolarAngleImage(
                    lch1[2], thetaShift);
                GreyscaleImage h2Shifted = imageProcessor.copyAndShiftPolarAngleImage(
                    lch2[2], thetaShift);
                
                
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

                    lch2[2] = tr.applyTransformation(lch2[2], params,
                            lch2[2].getHeight(), lch2[2].getWidth());
                }

                int x, y;
                
                np = 300;       

                ORB orb1 = new ORB(img1GS, np);
                orb1.detectAndExtract();
                orb1.overrideToAlsoCreate1stDerivKeypoints();
                orb1.overrideToNotCreateATrousKeypoints();

                ORB orb2 = new ORB(img2GS, np);
                orb2.detectAndExtract();
                orb1.overrideToAlsoCreate1stDerivKeypoints();
                orb2.overrideToNotCreateATrousKeypoints();

                Descriptors d1 = orb1.getAllDescriptors();
                Descriptors d2 = orb2.getAllDescriptors();
                List<PairInt> kp1 = orb1.getAllKeyPoints();
                List<PairInt> kp2 = orb2.getAllKeyPoints();
                QuadInt[] matched;

                if (ii == 2) {

                    HOGs hogs1 = new HOGs(img1GS, 1, 6);
                    HOGs hogs2 = new HOGs(img2GS, 1, 6);
                    TDoubleList orientations1 = orb1.getAllOrientations();
                    TDoubleList orientations2 = orb2.getAllOrientations();

                    matched = ORBMatcher.matchDescriptors(d1, d2, kp1, kp2,
                        hogs1, hogs2, orientations1, orientations2);

                } else {

                    matched = ORBMatcher.matchDescriptors(
                        d1, d2, kp1, kp2);
                }

                int x1, x2, y1, y2;
                {//DEBUG
                    Image tmp1 = img1.copyToGreyscale2().copyToColorGreyscale();
                    for (i = 0; i < kp1.size(); ++i) {
                        y = kp1.get(i).getY();
                        x = kp1.get(i).getX();
                        ImageIOHelper.addPointToImage(x, y, tmp1, 2, 255, 0, 0);
                    }
                    MiscDebug.writeImage(tmp1, "_kp_gs_" + lbl + fileName1Root);

                    Image tmp2 = img2GS.copyToColorGreyscale();
                    for (i = 0; i < kp2.size(); ++i) {
                        y = kp2.get(i).getY();
                        x = kp2.get(i).getX();
                        ImageIOHelper.addPointToImage(x, y, tmp2, 2, 255, 0, 0);
                    }
                    MiscDebug.writeImage(tmp2, "_kp_gs_" + lbl + fileName2Root);
                    System.out.println(lbl + fileName1Root + " matched=" + matched.length);
                    CorrespondencePlotter plotter = 
                        new CorrespondencePlotter(tmp1, tmp2);
                    for (i = 0; i < matched.length; ++i) {
                        x1 = matched[i].getA();
                        y1 = matched[i].getB();
                        x2 = matched[i].getC();
                        y2 = matched[i].getD();
                        plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 1);
                    }
                    plotter.writeImage("_corres_gs_" + lbl + fileName1Root);
                }


                ORB orb1c = new ORB(lch1[2], np);
                orb1c.detectAndExtract();

                ORB orb2c = new ORB(lch2[2], np);
                orb2c.detectAndExtract();

                ORB orb1csh = new ORB(h1Shifted, np);
                orb1csh.detectAndExtract();

                ORB orb2csh = new ORB(h2Shifted, np);
                orb2csh.detectAndExtract();

                Descriptors d1c = orb1c.getAllDescriptors();
                Descriptors d2c = orb2c.getAllDescriptors();
                List<PairInt> kp1c = orb1c.getAllKeyPoints();
                List<PairInt> kp2c = orb2c.getAllKeyPoints();

                Descriptors d1csh = orb1csh.getAllDescriptors();
                Descriptors d2csh = orb2csh.getAllDescriptors();
                List<PairInt> kp1csh = orb1csh.getAllKeyPoints();
                List<PairInt> kp2csh = orb2csh.getAllKeyPoints();
                TDoubleList o1csh = orb1csh.getAllOrientations();
                TDoubleList o2csh = orb2csh.getAllOrientations();
                
                Object[] combined1 = ORB.combine(d1c, kp1c, orb1c.getAllOrientations(), 
                    d1csh, kp1csh, o1csh);
                Descriptors combined1D = (Descriptors)combined1[0];
                @SuppressWarnings("unchecked")
                List<PairInt> combined1K = (List<PairInt>)combined1[1];
                TDoubleList combined1O = (TDoubleList)combined1[2];
                
                assertEquals(combined1D.descriptors.length, combined1K.size());
                assertEquals(combined1K.size(), combined1O.size());
        
                Object[] combined2 = ORB.combine(d2c, kp2c, orb2c.getAllOrientations(), 
                    d2csh, kp2csh, o2csh);
                Descriptors combined2D = (Descriptors)combined2[0];
                @SuppressWarnings("unchecked")
                List<PairInt> combined2K = (List<PairInt>)combined2[1];
                TDoubleList combined2O = (TDoubleList)combined2[2];
                assertEquals(combined2D.descriptors.length, combined2K.size());
                assertEquals(combined2O.size(), combined2K.size());

                // combining the descriptors
                List<Descriptors> d1CombList = new ArrayList<Descriptors>();
                d1CombList.add(d1);
                d1CombList.add(combined1D);
                List<Descriptors> d2CombList = new ArrayList<Descriptors>();
                d2CombList.add(d2);
                d2CombList.add(combined2D);

                Descriptors d1Comb = ORB.combineDescriptors(d1CombList);
                Descriptors d2Comb = ORB.combineDescriptors(d2CombList);
                kp1.addAll(combined1K);
                kp2.addAll(combined2K);
                matched = ORBMatcher.matchDescriptors(d1Comb, d2Comb, kp1, kp2);

                //TODO: revisit all of this and related code.
                if (matched == null) {
                    continue;
                }
                
                /*
                after results, if want to reverse transform the dataset2 images and points to
                more easily see the results
                img2GS = tr.applyTransformation(img2GS, paramsRev, img2GS.getHeight(), 
                    img2GS.getWidth());
                lch2[2] = tr.applyTransformation(lch2[2], paramsRev,
                    lch2[2].getHeight(), lch2[2].getWidth());
                tr.applyTransformation(paramsRev, kp2);            
                tr.applyTransformation(paramsRev, kp2c);            
                tr.applyTransformation(paramsRev, kp2Comb);
                */
                
                {//DEBUG
                    Image tmp1 = img1GS.copyToColorGreyscale();
                    ImageIOHelper.addCurveToImage(kp1, tmp1, 1, 255, 0, 0);
                    MiscDebug.writeImage(tmp1, "_kp_comb_" + lbl + fileName1Root);

                    Image tmp2 = img2GS.copyToColorGreyscale();
                    ImageIOHelper.addCurveToImage(kp2, tmp2, 1, 255, 0, 0);
                    MiscDebug.writeImage(tmp2, "_kp_comb_" + lbl + fileName2Root);
                    CorrespondencePlotter plotter = 
                        new CorrespondencePlotter(tmp1, tmp2);
                    for (i = 0; i < matched.length; ++i) {
                        x1 = matched[i].getA();
                        y1 = matched[i].getB();
                        x2 = matched[i].getC();
                        y2 = matched[i].getD();
                        plotter.drawLineInAlternatingColors(x1, y1, 
                            x2, y2, 1);
                    }
                    plotter.writeImage("_corres_comb_hogs_" + lbl + fileName1Root);
                }
            }
        }
    }
    
}
