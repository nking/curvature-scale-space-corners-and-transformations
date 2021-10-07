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
                
                QuadInt[] matched = ORBMatcher.matchDescriptors(d1, d2, kp1, kp2);
                
                int x1, x2, y1, y2;
                {//DEBUG
                    Image tmp1 = img1GS.copyToColorGreyscale();
                    for (i = 0; i < kp1.size(); ++i) {
                        y = kp1.get(i).getY();
                        x = kp1.get(i).getX();
                        ImageIOHelper.addPointToImage(x, y, tmp1, 2, 255, 0, 0);
                        /*if (x == 213 && y == 58) {
                            ImageIOHelper.addPointToImage(x, y, tmp1, 3, 255, 255, 0);
                        }*/
                    }
                    //MiscDebug.writeImage(tmp1, "_kp_gs_" + lbl + fileName1Root);

                    Image tmp2 = img2GS.copyToColorGreyscale();
                    for (i = 0; i < kp2.size(); ++i) {
                        y = kp2.get(i).getY();
                        x = kp2.get(i).getX();
                        ImageIOHelper.addPointToImage(x, y, tmp2, 2, 255, 0, 0);
                        /*if (x == 178 && y == 177) {
                            ImageIOHelper.addPointToImage(x, y, tmp2, 3, 255, 255, 0);
                        }*/
                    }
                    //MiscDebug.writeImage(tmp2, "_kp_gs_" + lbl + fileName2Root);
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
                    plotter.writeImage("_corres_orb_gs_" + lbl + fileName1Root);
                    
                    /*
                    Using DLT (see wikipedia):
                    
                    p2 = A * p1   where p1 is (x1, y1) and p2 is (x2, y2)
                                  A = a00 a01, a10 a11

                    A = P2 * P1^T * (P!*P1^T)^-1

                    H = 0  -1
                        1   0  anti-symmetric matrix

                    alpha * P2 = A * P1
                    multiply both sides of equation from the left by (P2)^T*H
                    (P2)^T*H * alpha * P2 = (P2)^T*H * A * P1
                    note that (P2)^T*H*P2 = 0
                    then have (P2)^T*H*A*P1 = 0

                    [x2 y2] * [0  -1 ] * [a00  a01 ] * [x1] = 0
                              [1   0 ]   [a10  a11 ]   [y1]

                    [ y2  -x2 ] * [a00  a01 ] * [x1] = 0
                                  [a10  a11 ]   [y1]

                    [ y2*a00 - x2*a10   y2*a01 - x2*a11 ] * [x1] = 0
                                                            [y1]

                      a00*y2*x1 - a10*x2*x1 + a01*y2*y1 - a11*x2*y1 = 0

                      can be written as 0 = (b)^T * a where b and a are vectors
                       A = (y2*x1)  x=(a00)
                           (y2*y1)    (a01)
                           (-x2*x1)   (a10)
                           (-x2*y1)   (a11)
                    
                  solving for x = argmin_x ( ||A*x||^2 ) subject to ||x||^2 = 1

                    SVD(A).vT[last row]
                    
                    reprojection error: 
                          y2*(a00*x1+a01*y1) + x2*(-a10*x1-a11*y1) = 0
                          y2*(a00*x1+a01*y1) = x2*(a10*x1+a11*y1) 
                          x2/y2 = ( (a00*x1+a01*y1)/(a10*x1+a11*y1) )
                    */
                    
                    SVDProducts svd = null;
                    double[][] ai = new double[1][5];
                    double[][] a = new double[matched.length][5];
                    double[] xHat;
                    double re = 0; // reprojection error for all
                    int _x2, _y2;
                    double[] res = new double[matched.length]; // reprojection error for each i
                    for (i = 0; i < matched.length; ++i) {
                        x1 = matched[i].getA();
                        y1 = matched[i].getB();
                        x2 = matched[i].getC();
                        y2 = matched[i].getD();
                        
                        /*
                        A = (y2*x1)  x=(a00)
                           (y2*y1)    (a01)
                           (-x2*x1)   (a10)
                           (-x2*y1)   (a11)
                        */
                        ai[0] = new double[]{y2*x1, y2*y1, -x2*x1, -x2*y1};
                        a[i] = Arrays.copyOf(ai[0], ai[0].length);
                       
                        svd = MatrixUtil.performSVD(ai);
                        xHat = svd.vT[svd.vT.length - 1];
                        for (int z = 0; z < xHat.length; ++z) {
                            xHat[z] /= xHat[xHat.length - 1];
                        }
                        //xHat = MatrixUtil.extractColumn(svd.u, svd.u.length - 1);
                        System.out.printf("%d (x1,y1), (x2,y2) = (%d,%d), (%d,%d)\n", i, x1, y1, x2, y2);
                        System.out.printf("xhat_%d=%s\n", i, FormatArray.toString(xHat, "%.3e"));
                        
                        /*
                          y2*(a00*x1+a01*y1) + x2*(-a10*x1-a11*y1) = 0
                          y2*(a00*x1+a01*y1) = x2*(a10*x1+a11*y1) 
                          x2/y2 = ( (a00*x1+a01*y1)/(a10*x1+a11*y1) )
                        */
                        _x2 = (int)(xHat[0]*x1+xHat[1]*y1); 
                        _y2 = (int)(xHat[2]*x1+xHat[3]*y1);
                        
                        //(u1*xhat[0]+v1*xhat[1] - u2)^2 + (u1*xhat[2]+v1*xhat[3] - v2)^2
                        res[i] = Math.pow(_x2 - x2, 2) + Math.pow(_y2 - y2, 2);
                        re += res[i];
                        System.out.printf("reprojection error=%.3f\n", res[i]);
                    }
                    svd = MatrixUtil.performSVD(a);
                    xHat = svd.vT[svd.vT.length - 1];
                    for (int z = 0; z < xHat.length; ++z) {
                        xHat[z] /= xHat[xHat.length - 1];
                    }
                    System.out.printf("xhat_all=%s\n", FormatArray.toString(xHat, "%.3e"));
                    System.out.printf("svd(a).s =%s\n", i, FormatArray.toString(svd.s, "%.3e"));
                    System.out.printf("reprojection error all/n=%.3f\n", re/(double)matched.length);
                    
                    // apply this to all image 1 points
                    tmp1 = img1GS.copyToColorGreyscale();
                    tmp2 = img2GS.copyToColorGreyscale();
                    plotter = new CorrespondencePlotter(tmp1, tmp2);
                    for (i = 0; i < kp1.size(); ++i) {
                        x1 = kp1.get(i).getX();
                        y1 = kp1.get(i).getY();
                        x2 = (int)(x1*xHat[0]+y1*xHat[1]); 
                        y2 = (int)(x1*xHat[2]+y1*xHat[3]);
                        plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 1);
                    }
                    plotter.writeImage("_corres_orb_gs_proj_" + lbl + fileName1Root);
                    
                    /*NOTE: from this next outlier removal, can see in the images
                       that the rotation of each transformation is also needed for
                       use in outlier removal.
                       
                       Rotation is obtained from the affine reconstruction:
                          OrthographicProjectionResults re = Reconstruction.calculateAffineReconstruction(
                               double[][] x, int mImages).

                           where OrthographicProjectionResults results = new OrthographicProjectionResults();
                           results.XW = s;
                           results.rotationMatrices = rotStack;
                    
                       Further considering algorithm in paper on outlier correction:
                           https://www.researchgate.net/publication/221110532_Outlier_Correction_in_Image_Sequences_for_the_Affine_Camera
                       "Outlier Correction in Image Sequences for the Affine Camera"
                          by Huynh, Hartley, and Heydeon 2003
                          Proceedings of the Ninth IEEE International Conference on Computer Vision (ICCV’03)
                    */
                    
                    tmp1 = img1GS.copyToColorGreyscale();
                    tmp2 = img2GS.copyToColorGreyscale();
                    plotter = new CorrespondencePlotter(tmp1, tmp2);
                    // should instead use  median absolute deviation and interquartile range or other robust dispersion methods
                    double[] meanStdev = MiscMath.getAvgAndStDev(res);
                    System.out.printf("reprojection error mean=%.4e stdev=%.4e\n", meanStdev[0], meanStdev[1]);
                    int count = 0;
                    re = 0;
                    for (i = 0; i < matched.length; ++i) {
                        if (Math.abs(res[i] - meanStdev[0]) <= meanStdev[1]/2.){//2.5*meanStdev[1]) {
                            x1 = matched[i].getA();
                            y1 = matched[i].getB();
                            x2 = matched[i].getC();
                            y2 = matched[i].getD();
                            re += res[i];
                            a[count] = new double[]{y2*x1, y2*y1, -x2*x1, -x2*y1};                            
                            count++;
                            plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 1);
                        }
                    }
                    plotter.writeImage("_corres_orb_gs_filtered_" + lbl + fileName1Root);
                    System.out.printf("number removed=%d\n", matched.length - count);
                    a = Arrays.copyOf(a, count);
                    svd = MatrixUtil.performSVD(a);
                    xHat = svd.vT[svd.vT.length - 1];
                    for (int z = 0; z < xHat.length; ++z) {
                        xHat[z] /= xHat[xHat.length - 1];
                    }
                    System.out.printf("xhat_filtered=%s\n", FormatArray.toString(xHat, "%.3e"));
                    System.out.printf("svd(a).s =%s\n", i, FormatArray.toString(svd.s, "%.3e"));
                    System.out.printf("reprojection error filtered/n=%.3f\n", re/(double)count);
                    
                }
            }
        }
    }
    
}
