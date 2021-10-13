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
        
        boolean binImages = true;
        
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
            orb1.overrideToAlsoCreate1stDerivKeypoints();
            orb1.overrideToNotCreateATrousKeypoints();

            ORB orb2 = new ORB(img2GS, np);
            orb2.detectAndExtract();
            orb1.overrideToAlsoCreate1stDerivKeypoints();
            orb2.overrideToNotCreateATrousKeypoints();

            ORB.Descriptors d1 = orb1.getAllDescriptors();
            ORB.Descriptors d2 = orb2.getAllDescriptors();
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

                MatrixUtil.SVDProducts svd = null;
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

                double[][] xIn = MatrixUtil.zeros(2, 2*count);
                double[][] x1In = MatrixUtil.zeros(3, count);
                double[][] x2In = MatrixUtil.zeros(3, count);
                count = 0;
                for (i = 0; i < matched.length; ++i) {
                    if (Math.abs(res[i] - meanStdev[0]) <= meanStdev[1]/2.){//2.5*meanStdev[1]) {
                        x1 = matched[i].getA();
                        y1 = matched[i].getB();
                        xIn[0][count] = x1;
                        xIn[1][count] = y1;
                        x1In[0][count] = x1;
                        x1In[1][count] = y1;
                        x1In[2][count] = 1;
                        x2In[0][count] = matched[i].getC();
                        x2In[1][count] = matched[i].getD();
                        x2In[2][count] = 1;
                        count++;
                    }
                }
                for (i = 0; i < matched.length; ++i) {
                    if (Math.abs(res[i] - meanStdev[0]) <= meanStdev[1]/2.){//2.5*meanStdev[1]) {
                        x2 = matched[i].getC();
                        y2 = matched[i].getD();
                        xIn[0][count] = x2;
                        xIn[1][count] = y2;
                        count++;
                    }
                }
                return xIn;
        }
    }
}
