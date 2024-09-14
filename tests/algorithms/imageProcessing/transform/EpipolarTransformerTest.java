package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.matching.ORBMatcher;
import algorithms.imageProcessing.transform.EpipolarTransformer.NormalizationTransformations;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath;
import algorithms.misc.MiscMath0;
import algorithms.statistics.Standardization;
import algorithms.util.*;

import java.awt.Color;
import java.io.IOException;
import java.util.*;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertNotNull;
import no.uib.cipr.matrix.DenseMatrix;

/**
 *
 * @author nichole
 */
public class EpipolarTransformerTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public EpipolarTransformerTest() {
    }
    
    public void testNormalization0() throws Exception {
        
        int m = 3;
        int n = 24;
        
        double[][] x1 = new double[m][n];
        double[][] x2 = new double[m][n];
        x1[0] = new double[]{262, 316, 260, 284, 234, 177, 216
           , 220, 248, 248, 319, 159, 176, 407, 393, 119, 117, 428, 427, 112, 109, 425, 256, 411
        };
        x1[1] = new double[]{356, 342, 305, 279, 217, 76, 63
            , 158, 27, 46, 36, 54, 76, 64, 85,       115, 141, 120, 147, 320, 375, 333, 192, 213
        };
        x1[2] = new double[n];  Arrays.fill(x1[2], 1.0);
        
        x2[0] = new double[]{156, 191, 153, 167, 135, 97, 119
            , 125, 137, 137, 183, 87, 97, 248, 238, 71, 70, 267, 267, 73, 74, 272,    147, 256
        };
        x2[1] = new double[]{308, 301, 270, 249, 202, 97, 83
            , 156, 51,  65,   49,  83, 97, 61,  81, 130, 149, 110, 134, 276, 315, 299, 180, 192
        };
        x2[2] = new double[n];  Arrays.fill(x2[2], 1.0);
        
        DenseMatrix x1M = new DenseMatrix(x1);
        DenseMatrix x2M = new DenseMatrix(x2);
        
        EpipolarTransformer.NormalizedXY normXY1 = EpipolarTransformer.normalize(x1M);
        EpipolarTransformer.NormalizedXY normXY2 = EpipolarTransformer.normalize(x2M);
        DenseMatrix leftNorm = normXY1.getXy();
        DenseMatrix rightNorm = normXY2.getXy();
        double[][] _leftNorm = MatrixUtil.convertToRowMajor(leftNorm);
        double[][] _rightNorm = MatrixUtil.convertToRowMajor(rightNorm);
         
        double tol = 1e-5;
        double tol2 = 1.e-2;
        double expected = 1;
        
        double[] meanAndStDevLeftX = MiscMath.getAvgAndStDev(_leftNorm[0]);
        double[] meanAndStDevLeftY = MiscMath.getAvgAndStDev(_leftNorm[1]);
        double[] meanAndStDevRightX = MiscMath.getAvgAndStDev(_rightNorm[0]);
        double[] meanAndStDevRightY = MiscMath.getAvgAndStDev(_rightNorm[1]);
        System.out.printf("mn&stDev left X=%s\n", FormatArray.toString(meanAndStDevLeftX, "%.3f"));
        System.out.printf("mn&stDev left Y=%s\n", FormatArray.toString(meanAndStDevLeftY, "%.3f"));
        System.out.printf("mn&stDev right X=%s\n", FormatArray.toString(meanAndStDevRightX, "%.3f"));
        System.out.printf("mn&stDev right Y=%s\n", FormatArray.toString(meanAndStDevRightY, "%.3f"));

        double rmsLeft = Math.sqrt(meanAndStDevLeftX[1] * meanAndStDevLeftX[1] +
            meanAndStDevLeftY[1]*meanAndStDevLeftY[1]);
        double rmsRight = Math.sqrt(meanAndStDevRightX[1] * meanAndStDevLeftX[1] +
            meanAndStDevLeftY[1]*meanAndStDevRightY[1]);
        System.out.printf("rmsLeft=%.4e  rmsRight=%.4e\n\n", rmsLeft, rmsRight);
        
        assertTrue(Math.abs(meanAndStDevLeftX[0]) < tol);
        assertTrue(Math.abs(meanAndStDevLeftY[0]) < tol);
        assertTrue(Math.abs(rmsLeft - expected) < tol2);
        
        assertTrue(Math.abs(meanAndStDevRightX[0]) < tol);
        assertTrue(Math.abs(meanAndStDevRightY[0]) < tol);
        assertTrue(Math.abs(rmsRight - expected) < tol2);
      
        
        NormalizationTransformations leftNT = normXY1.getNormalizationMatrices();
        NormalizationTransformations rightNT = normXY2.getNormalizationMatrices();
       
        double diffX, diffY,diffZ;
        
        double[][] x1Denorm = MatrixUtil.multiply(leftNT.tDenorm, _leftNorm);
        double[][] x2Denorm = MatrixUtil.multiply(rightNT.tDenorm, _rightNorm);
        
        for (int j = 0; j < x2[0].length; ++j) {
            diffX = Math.abs(x2[0][j]- x2Denorm[0][j]);
            diffY = Math.abs(x2[1][j]- x2Denorm[1][j]);
            diffZ = Math.abs(x2[2][j]- x2Denorm[2][j]);
            assertTrue(diffX < tol);
            assertTrue(diffY < tol);
            assertTrue(diffZ < tol);
            diffX = Math.abs(x1[0][j]- x1Denorm[0][j]);
            diffY = Math.abs(x1[1][j]- x1Denorm[1][j]);
            diffZ = Math.abs(x1[2][j]- x1Denorm[2][j]);
            assertTrue(diffX < tol);
            assertTrue(diffY < tol);
            assertTrue(diffZ < tol);
        }
        
        /*
         note:
            because translation is not a linear transformation (see Strang Chap 7)
               one has to keep it as a separate tranformation matrix when
               performing operations on a sequence of matrices such as inverse
               and transpose operations.
        
            translation matrix: inverse changes the signs of the elements.
            
            rotation matrix: inverse is the transpose of rotation matrix.

            (A*B*C)^-1 = (C^-1) * (B^-1) * (A^-1)

        also, when A * A^(-1) = I, one can use:
                        1
            A^(-1) =  ------ C^(T)  where C_ij = cofactor of a_ij
                       det A
        
        Revisiting Hartley’s Normalized Eight-Point Algorithm
          Chojnacki et al. 2003
        
          u1_normalized = T1 * u1
          u2_normalized = T2 * u2

          denormalized u1 = T1^-1 * u1_normalized 
          denormalized u2 = T2^-1 * u2_normalized

          FM_normalized = inverse(transpose(T2)) * FM * inverse(T1)
            with caveat about centroid and normalization details...

         FM_normalized * T1 = inverse(transpose(T2)) * FM
         transpose(T2) * FM_normalized * T1 = FM

      ==> denormalized FM = transpose(T2) * FM_normalized * T1
        
          u2^T * FM * u1 = u2_normalized^T * FM_normalized * u1_normalized = residual

              | 1/s   0  0 |   | 1  0  -xc |   | 1/s    0   -xc/s |
          T = |  0  1/s  0 | * | 0  1  -yc | = |   0  1/s   -yc/s |
              |  0    0  1 |   | 0  0   1  |   |   0    0      1  |

                 | 1  0  xc |   | s  0   0 |   | s   0  xc |
          T^-1 = | 0  1  yc | * | 0  s   0 | = | 0   s  yc |
                 | 0  0   1 |   | 0  0   1 |   | 0   0   1 |
   
                        | 1  0  xc |   | s^2  0    0 |   | 1  0  0 |
          T^-1 * T^-T = | 0  1  yc | * | 0   s^2   0 | * | 0  1  0 |
                        | 0  0   1 |   | 0    0    1 |   | xc yc 1 |

                        | s^2 + xc^2   xc*yc       xc |
                      = | yc*xc        s^2 + yc^2  yc |
                        | xc           yc          1  |
        
                                  | s  0   0 |   | 1  0  0 |   |  s   0   0 |
        from that can see  T^-T = | 0  s   0 | * | 0  1  0 | = |  0   s   0 |
                                  | 0  0   1 |   | xc yc 1 |   | xc  yc   1 |
        
              |  1    0    0 |   | 1/s   0  0 |   |   1/s    0    0 |
        T^T = |  0    1    0 | * |  0  1/s  0 | = |     0   1/s   0 |
              |-xc  -yc    1 |   |  0    0  1 |   | -xc/s  -yc/s  1 |
        
                   ( | 1/s   0  0 | )^-1     |  1    0    0 |   |  s  0   0 |
        (T^T)^-1 = ( |  0  1/s  0 | )     *  |  0    1    0 | = |  0  s   0 |
                   ( |  0    0  1 | )        | xc   yc    1 |   | xc  yc  1 |
        
        can see that (T^-1)^T = (T^T)^-1         
        */
       
        /*
        hartley 1997 house 
        
        temporary test until have external examples to test.
        
        checking that results of this code make sense in terms of the matrix operations
        above.
        
            FM of Hartley 1997 in his publication:
            [junit] 4.144e-06, -6.231e-05, 2.817e-02 
            [junit] 2.684e-05, -4.437e-07, -3.188e-01 
            [junit] -3.852e-02, 3.274e-01, 1.000e+00 
            [junit] 
            
            this code's normalized FM=
            [junit] -2.430e-03, -5.531e-02, 3.050e-02 
            [junit] 6.257e-02, -4.565e-03, -7.036e-01 
            [junit] -3.560e-03, 7.049e-01, -2.492e-04 
            [junit] 
            
            this code's de-normalized FM=
            [junit] 3.455e-06, 7.863e-05, -2.257e-02 
            [junit] -8.895e-05, 6.490e-06, 1.246e-01 
            [junit] 1.920e-02, -1.255e-01, 1.000e+00 
        
        [junit] T_left=
        [junit] 9.3917e-03, 0.0000e+00, -2.0586e+00 <== xc1=-219.19
        [junit] 0.0000e+00, 9.3917e-03, -2.0539e+00 <== yc=-218.69
        [junit] 0.0000e+00, 0.0000e+00, 1.0000e+00 
        [junit] 
        [junit] T_right=
        [junit] 9.3503e-03, 0.0000e+00, -2.0132e+00 <== -215.31
        [junit] 0.0000e+00, 9.3503e-03, -2.0401e+00 <== -218.19
        [junit] 0.0000e+00, 0.0000e+00, 1.0000e+00
        
        */
        
        double factor;
        
        // normalized
        double[][] fm0 = new double[3][3];
        fm0[0] = new double[]{-2.430e-03, -5.531e-02, 3.050e-02};
        fm0[1] = new double[]{6.257e-02, -4.565e-03, -7.036e-01};
        fm0[2] = new double[]{-3.560e-03, 7.049e-01, -2.492e-04};
        double[][] fm00 = MatrixUtil.copy(fm0);
        factor = 1./(fm00[2][2]);
        MatrixUtil.multiply(fm00, factor);
        
        //de-normalized
        double[][] dfm0 = new double[3][3];
        dfm0[0] = new double[]{3.455e-06, 7.863e-05, -2.257e-02};
        dfm0[1] = new double[]{-8.895e-05, 6.490e-06, 1.246e-01};
        dfm0[2] = new double[]{1.920e-02, -1.255e-01, 1.000e+00};
        
        double[][] tLeft = new double[3][3];
        tLeft[0] = new double[]{9.3917e-03, 0.0000e+00, -2.0586e+00};
        tLeft[1] = new double[]{0.0000e+00, 9.3917e-03, -2.0539e+00};
        tLeft[2] = new double[]{0.0000e+00, 0.0000e+00, 1.0000e+00};
        
        double[][] tRight = new double[3][3];
        tRight[0] = new double[]{9.3503e-03, 0.0000e+00, -2.0132e+00};
        tRight[1] = new double[]{0.0000e+00, 9.3503e-03, -2.0401e+00};
        tRight[2] = new double[]{0.0000e+00, 0.0000e+00, 1.0000e+00};
        
        double[][] tLeftInv = new double[3][3]; 
        tLeftInv[0] = new double[]{1.0648e+02, 0.0000e+00, 2.1919e+02};
        tLeftInv[1] = new double[]{0.0000e+00, 1.0648e+02, 2.1869e+02};
        tLeftInv[2] = new double[]{0.0000e+00, 0.0000e+00, 1.0000e+00};
        
        double[][] tRightInv = new double[3][3]; 
        tRightInv[0] = new double[]{1.0695e+02, 0.0000e+00, 2.1531e+02};
        tRightInv[1] = new double[]{0.0000e+00, 1.0695e+02, 2.1819e+02};
        tRightInv[2] = new double[]{0.0000e+00, 0.0000e+00, 1.0000e+00};
        double[][] tRightInvTranspose = MatrixUtil.transpose(tRightInv);
        
        tol = 1e-2;
        //denormalized FM = transpose(T2) * FM_normalized * T1
        double[][] dfm0Tst = MatrixUtil.multiply(MatrixUtil.transpose(tRight), fm0);
        dfm0Tst = MatrixUtil.multiply(dfm0Tst, tLeft);
        factor = 1./(dfm0Tst[2][2]);
        double[][] dfm0Tst0 = MatrixUtil.copy(dfm0Tst);
        MatrixUtil.multiply(dfm0Tst0, factor);
        assertTrue(assertSimilar(dfm0, dfm0Tst0, tol));
        
        //FM_normalized = inverse(transpose(T2)) * denormalized FM * inverse(T1)
        double[][] fm0Tst = MatrixUtil.multiply(tRightInvTranspose, dfm0Tst);
        fm0Tst = MatrixUtil.multiply(fm0Tst, tLeftInv);
        factor = 1./(fm0Tst[2][2]);
        MatrixUtil.multiply(fm0Tst, factor);
        // factor of 1.2 between them
        //assertTrue(assertSimilar(fm00, fm0Tst, tol));
        
        double[][] dfm0Tst2 = MatrixUtil.multiply(MatrixUtil.transpose(tRight), fm0Tst);
        dfm0Tst2 = MatrixUtil.multiply(dfm0Tst2, tLeft);
        factor = 1./(dfm0Tst2[2][2]);
        MatrixUtil.multiply(dfm0Tst2, factor);
        assertTrue(assertSimilar(dfm0, dfm0Tst2, tol));
    }
    private boolean assertSimilar(double[][] a, double[][] b, double tol) {
        double diff;
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                diff = Math.abs(a[i][j] - b[i][j]);
                if (diff > tol) {
                    return false;
                }
            }
        }
        return true;
    }

    public void estMoreThan7Points1() throws Exception {
        System.out.println("testMoreThan7Points1");

        PairIntArray leftTrueMatches = new PairIntArray();
        PairIntArray rightTrueMatches = new PairIntArray();
        getMertonCollege10TrueMatches(leftTrueMatches, rightTrueMatches);
        
        EpipolarTransformer spTransformer = new EpipolarTransformer();
        
        double[][] left = convertToDouble(leftTrueMatches);
        double[][] right = convertToDouble(leftTrueMatches);
        
        EpipolarTransformer.NormalizedXY normXY1 = EpipolarTransformer.normalize(new DenseMatrix(left));
        EpipolarTransformer.NormalizedXY normXY2 = EpipolarTransformer.normalize(new DenseMatrix(right));
        DenseMatrix leftM = normXY1.getXy();
        DenseMatrix rightM = normXY2.getXy();
        
        DenseMatrix normalizedFM = spTransformer.calculateEpipolarProjection(
            leftM, rightM, false);

        assertNotNull(normalizedFM);
        
        DenseMatrix fm = EpipolarTransformer.denormalizeTheFundamentalMatrix(
            normalizedFM, normXY1.getNormalizationMatrices(),
            normXY2.getNormalizationMatrices());
        
        boolean useToleranceAsStatFactor = true;//false;
        final double tolerance = 3.88;
        ErrorType errorType = ErrorType.DIST_TO_EPIPOLAR_LINE;
        
        Distances distances = new Distances();
        EpipolarTransformationFit fit = null;
        if (useToleranceAsStatFactor) {
            fit = distances.calculateError2(normalizedFM,
                leftM, rightM, errorType, tolerance);
        } else {
            fit = distances.calculateError(normalizedFM,
                leftM, rightM, errorType, tolerance);
        }
        
        assertNotNull(fit);
        assertNotNull(fit.getErrorType());
        System.out.println("fit=" + fit.toString());
        assertEquals(10, fit.getErrors().size());
        assertEquals(10, fit.getInlierIndexes().size());
        fit.calculateErrorStatistics();
        double mean = fit.getMeanError();
        double stdev = fit.getStDevFromMean();
        List<Double> errors = fit.getErrors();
        
        assertTrue(mean < tolerance);
        assertTrue(stdev < tolerance);
        for (Double error : errors) {
            assertTrue(error < tolerance);
        }
        
        String fileName1 = "merton_college_I_001.jpg";
        String fileName2 = "merton_college_I_002.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);

        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        int image2Width = img2.getWidth();
        int image2Height = img2.getHeight();
        
        overplotEpipolarLines(fm, 
            leftTrueMatches, rightTrueMatches,
            img1, img2, 
            image1Width, image1Height, image2Width, image2Height, 
            "e_merton_college_10pts" + Integer.valueOf(0).toString());
    }

    public void testCorners() throws IOException {
        String fileName1 = "merton_college_I_001.jpg";
        String fileName2 = "merton_college_I_002.jpg";
        //fileName1 = "brown_lowe_2003_image1.jpg";
        //fileName2 = "brown_lowe_2003_image2.jpg";
        //fileName1 = "books_illum3_v0_695x555.png";
        //fileName2 = "books_illum3_v6_695x555.png";

        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        GreyscaleImage image1 = ImageIOHelper.readImageAsGreyscaleFullRange(filePath1);
        GreyscaleImage image2 = ImageIOHelper.readImageAsGreyscaleFullRange(filePath2);
        double[][] img1 = convertToDouble(image1.toNormalizedRowMajorArray());
        double[][] img2 = convertToDouble(image2.toNormalizedRowMajorArray());

        ImageProcessor imageProcessor = new ImageProcessor();

        int minDist = 1;
        float thresholdRel = 0.01f;//0.1f;
        boolean ignore0sInThresh = true;

        // sorting by descending strength of score.
        // format is [nCorners X 2] where each row of the keypoint list is a pair of (row, col)
        int[][] keyPoints1 = ImageProcessor.calcHarrisCorners(img1, minDist, thresholdRel, ignore0sInThresh);
        int[][] keyPoints2 = ImageProcessor.calcHarrisCorners(img2, minDist, thresholdRel, ignore0sInThresh);

        //writeToCWD(keyPoints1, keyPoints2);

        overplotKeypoints(keyPoints1, keyPoints2, image1.copyToColorGreyscale(), image2.copyToColorGreyscale(),
                Integer.toString(0));

        System.out.printf("nKP=%d, %d\n", keyPoints1.length, keyPoints2.length);

        int nKP = Math.min(keyPoints1.length, keyPoints2.length);
        //nKP = 100; can reduce the number to see the matches more easily
        nKP = (int)(2.5*nKP);

        ORB orb1 = new ORB(image1, nKP);
        ORB orb2 = new ORB(image2, nKP);

        orb1.detectAndExtract();
        orb2.detectAndExtract();

        float[] scales1 = orb1.getScales();
        // [0]=24X32, [1]=40X53, [5]=160X213
        int nOct = scales1.length/3;//scales1.length - 1;
        ORB.Descriptors d1 = orb1.getDescriptorsList().get(nOct);
        ORB.Descriptors d2 = orb2.getDescriptorsList().get(nOct);
        double[][] xKP1 = orb1.getKeyPointsHomogenous(nOct);
        double[][] xKP2 = orb2.getKeyPointsHomogenous(nOct);

        /*
        d1 = orb1.getAllDescriptors();
        d2 = orb2.getAllDescriptors();
        xKP1 = orb1.getAllKeyPointsHomogenous();
        xKP2 = orb2.getAllKeyPointsHomogenous();
        */

        //writeToCWD(kp1, kp2);

        //int[][] indexesOfSeveralTrueMatches = indexesOfSeveralTrueMatchesMerton(kp1, kp2);
        //printDescriptorMatches(d1, d2, indexesOfSeveralTrueMatches);

        //System.out.printf("several true matches, indexes:\n%s\n",
        //        FormatArray.toString(indexesOfSeveralTrueMatches, "%d"));

        System.out.printf("# keypoints in 1 and 2 = %d, %d\n", xKP1[0].length, xKP2[0].length);

        String dirPath = ResourceFinder.findDirectory("bin");

        int i;
        int r1, r2, c1, c2;
        int row, col;
        Image tmp1 = image1.copyToColorGreyscale();
        for (i = 0; i < xKP1[0].length; ++i) {
            col = (int)xKP1[0][i];
            row = (int)xKP1[1][i];
            ImageIOHelper.addPointToImage(col, row, tmp1, 2, 255, 0, 0);
        }
        //MiscDebug.writeImage(tmp1, "_kp_gs_" + lbl + fileName1Root);

        Image tmp2 = image2.copyToColorGreyscale();
        for (i = 0; i < xKP2[0].length; ++i) {
            col = (int)xKP2[0][i];
            row = (int)xKP2[1][i];
            ImageIOHelper.addPointToImage(col, row, tmp2, 2, 255, 0, 0);
        }
        ImageIOHelper.writeOutputImage(
                dirPath + "/m1_harriscorners_02" + ".png", tmp1);
        ImageIOHelper.writeOutputImage(
                dirPath + "/m2_harriscorners_02" + ".png", tmp2);

        Image tmp3 = image1.copyToColorGreyscale();
        Image tmp4 = image2.copyToColorGreyscale();


        // matched are format row1, col1, row2, col2

        double[][] xKP1n = MatrixUtil.copy(xKP1);
        double[][] xKP2n = MatrixUtil.copy(xKP2);
        double[][] t1 = EpipolarNormalizationHelper.unitStandardNormalize(xKP1n);
        double[][] t2 = EpipolarNormalizationHelper.unitStandardNormalize(xKP2n);

        ORBMatcher.FitAndCorres fitC = null;
        if (true /*normalize points*/) {
            fitC = ORBMatcher.matchDescriptors(d1, d2, xKP1n, xKP2n);
        } else {
            fitC = ORBMatcher.matchDescriptors(d1, d2, xKP1, xKP2);
        }
        if (fitC.mIF == null) {
            return;
        }
        System.out.println(fileName1 + " matched=" + fitC.mIF.length);
        CorrespondencePlotter plotter =
                new CorrespondencePlotter(tmp1, tmp2);
        Color clr = null;
        int idx1, idx2;
        for (i = 0; i < fitC.mIF.length; ++i) {
            idx1 = fitC.mIF[i][0];
            idx2 = fitC.mIF[i][1];
            r1 = (int)xKP1[1][idx1];
            c1 = (int)xKP1[0][idx1];
            r2 = (int)xKP2[1][idx2];
            c2 = (int)xKP2[0][idx2];
            plotter.drawLineInAlternatingColors(c1, r1, c2, r2, 1);

            clr = getColor(clr);
            ImageIOHelper.addPointToImage(c1, r1, tmp3, 2, clr.getRed(), clr.getGreen(), clr.getBlue());
            ImageIOHelper.addPointToImage(c2, r2, tmp4, 2, clr.getRed(), clr.getGreen(), clr.getBlue());
        }
        plotter.writeImage("_corres_orb_gs_merton_college_I_001_002");

        ImageIOHelper.writeOutputImage(
                dirPath + "/matched_merton_college_I_001" + ".png", tmp3);
        ImageIOHelper.writeOutputImage(
                dirPath + "/matched_merton_college_I_002" + ".png", tmp4);

    }

    private void printDescriptorMatches(ORB.Descriptors d1, ORB.Descriptors d2,
                                        int[][] indexesOfSeveralTrueMatches) {

        int[][] cost = ORB.calcDescriptorCostMatrix(d1.descriptors, d2.descriptors);

        // statistics of cost matrix
        double[] costD = stack(cost);
        double[] mADMinMax = MiscMath0.calculateMedianOfAbsoluteDeviation(costD);
        double kMAD = 1.4826;
        double s = kMAD*mADMinMax[0];
        double r0 = mADMinMax[1] - 3 * s;
        double r1 = mADMinMax[1] + 3 * s;
        double[] medianAndIQR = MiscMath0.calcMedianAndIQR(costD);
        System.out.printf("cost matrix median=%.3e, median of abs dev=%.3e,  s=%.3e\n",
                mADMinMax[1], mADMinMax[0], s);
        System.out.printf("cost matrix median=%.3e, IQR=%.3e\n", medianAndIQR[0], medianAndIQR[1]);

        int i;
        int i1, i2;
        double c;
        for (i = 0; i < indexesOfSeveralTrueMatches.length; ++i) {
            i1 = indexesOfSeveralTrueMatches[i][0];
            i2 = indexesOfSeveralTrueMatches[i][1];
            System.out.printf("cost of idxs %d:%d = %d", i1, i2, cost[i1][i2]);
            c = (mADMinMax[1] - cost[i1][i2])/s;
            if (c > 2.5) {
                System.out.printf(" * ");
            }
            System.out.println();
        }
    }

    private double[] stack(int[][] cost) {
        int n = cost.length * cost[0].length;
        double[] s = new double[n];
        int i, j;
        int c = 0;
        for (i = 0; i < cost.length; ++i) {
            for (j = 0; j < cost[0].length; ++j) {
                s[c] = cost[i][j];
                ++c;
            }
        }
        return s;
    }

    private void writeToCWD(List<PairInt> keyPoints1, List<PairInt> keyPoints2) throws IOException {
        StringBuilder sb = new StringBuilder();
        int i, j;
        int row, col;
        sb.append("# row, col is format for compatibility with row-major data\n");
        for (i = 0; i < keyPoints1.size(); ++i) {
            col = keyPoints1.get(i).getY();
            row = keyPoints1.get(i).getX();
            sb.append(Integer.toString(row)).append(",").append(Integer.toString(col)).append("\n");
        }
        String fileName =  "merton_001_corners.txt";
        ResourceFinder.writeToCWD(sb.toString(), fileName);

        sb = new StringBuilder();
        sb.append("# row, col is format for compatibility with row-major data\n");
        for (i = 0; i < keyPoints2.size(); ++i) {
            col = keyPoints2.get(i).getY();
            row = keyPoints2.get(i).getX();
            sb.append(Integer.toString(row)).append(",").append(Integer.toString(col)).append("\n");
        }
        fileName =  "merton_002_corners.txt";
        ResourceFinder.writeToCWD(sb.toString(), fileName);
    }

    private void writeToCWD(int[][] keyPoints1, int[][] keyPoints2) throws IOException {
        StringBuilder sb = new StringBuilder();
        int i, j;
        int row, col;
        sb.append("# row, col is format for compatibility with row-major data\n");
        for (i = 0; i < keyPoints1.length; ++i) {
            col = keyPoints1[i][1];
            row = keyPoints1[i][0];
            sb.append(Integer.toString(row)).append(",").append(Integer.toString(col)).append("\n");
        }
        String fileName =  "merton_001_corners.txt";
        ResourceFinder.writeToCWD(sb.toString(), fileName);

        sb = new StringBuilder();
        sb.append("# row, col is format for compatibility with row-major data\n");
        for (i = 0; i < keyPoints2.length; ++i) {
            col = keyPoints2[i][1];
            row = keyPoints2[i][0];
            sb.append(Integer.toString(row)).append(",").append(Integer.toString(col)).append("\n");
        }
        fileName =  "merton_002_corners.txt";
        ResourceFinder.writeToCWD(sb.toString(), fileName);
    }

    public void testORB() throws IOException {
        /*
         >>> from skimage.feature import corner_harris, corner_peaks
        >>> import numpy as np
        >>> square3 = np.zeros([10, 10])
        >>> square3[2:8, 2:8] = 1
        >>> square3.astype(int)
        array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 3, 2, 1, 1, 1, 1, 0, 0],
               [0, 0, 2, 1, 1, 1, 1, 1, 0, 0],
               [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
               [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
               [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
               [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
        */
        GreyscaleImage image = new GreyscaleImage(10, 11);
        double[][] img = MatrixUtil.zeros(11, 10);
        int i;
        int j;
        for (i = 3; i < 9; ++i) {
            for (j = 2; j < 8; ++j) {
                img[i][j] = 1;
                image.setValue(j, i, 1);
            }
        }
        img[3][2] = 3;
        img[4][2] = 2;
        img[3][3] = 2;
        image.setValue(2, 3, 3);
        image.setValue(2, 4, 2);
        image.setValue(3, 3, 2);

        Transformer tr = new Transformer();
        GreyscaleImage image2 = image.copyImage();
        TransformationParameters params = new TransformationParameters();
        params.setTranslationY(image2.getWidth()-1);
        // counter-clockwise rotation by 90.
        params.setRotationInDegrees(90);
        image2 = tr.applyTransformation(image2, params, image2.getHeight(), image2.getWidth());

        int[][] expected = new int[4][2];
        expected[0] = new int[]{3, 2};
        expected[1] = new int[]{3, 7};
        expected[2] = new int[]{8, 2};
        expected[3] = new int[]{8, 7};

        int nKP = 3*expected.length;

        ORB orb1 = new ORB(image, nKP);
        ORB orb2 = new ORB(image2, nKP);

        orb1.detectAndExtract();
        orb2.detectAndExtract();

        ORB.Descriptors d1 = orb1.getAllDescriptors();
        ORB.Descriptors d2 = orb2.getAllDescriptors();
        double[][] xKP1 = orb1.getAllKeyPointsHomogenous();
        double[][] xKP2 = orb2.getAllKeyPointsHomogenous();

        double[][] xKP1n = MatrixUtil.copy(xKP1);
        double[][] xKP2n = MatrixUtil.copy(xKP2);
        double[][] t1 = EpipolarNormalizationHelper.unitStandardNormalize(xKP1n);
        double[][] t2 = EpipolarNormalizationHelper.unitStandardNormalize(xKP2n);

        ORBMatcher.FitAndCorres fitC = null;
        if (true /*normalized points*/) {
            fitC = ORBMatcher.matchDescriptors(d1, d2, xKP1n, xKP2n);
        } else {
            fitC = ORBMatcher.matchDescriptors(d1, d2, xKP1, xKP2);
        }
        /*
        test expecting in (row, col) format
                  (8, 2) <--> (7, 8)  => idx 0, 1  (bits diff=0)
                  (8, 7) <--> (2, 8)         1, 0  (bits diff=0)
                  (3, 7) <--> (2, 3)         2, 2  (bits diff=0)
                  (3, 2) <--> (7, 3)         3, 3  (bits diff=0)
         */

        // row1, col1, row2, col2
        Set<QuadInt> expectedM = new HashSet<QuadInt>();
        expectedM.add(new QuadInt(8, 2, 7, 8));
        expectedM.add(new QuadInt(8, 7, 2, 8));
        expectedM.add(new QuadInt(3, 7, 2, 3));
        expectedM.add(new QuadInt(3, 2, 7, 3));
    //    assertEquals(expectedM.size(), matchedIdxs.length);
        int idx1, idx2;
        QuadInt m;
        for (i = 0; i < fitC.mI.length; ++i) {
            idx1 = fitC.mI[i][0];
            idx2 = fitC.mI[i][1];
            m = new QuadInt((int)xKP1[1][idx1], (int)xKP1[0][idx1], (int)xKP2[1][idx2],(int)xKP2[0][idx2]);
            System.out.printf("matched=%s\n", m.toString());
            expectedM.remove(m);
        }
        assertTrue(expectedM.size() >= 2);
    }

    /*
    for more datasets:
    http://www.robots.ox.ac.uk/~vgg/data/
    */
    
    private void overplotEpipolarLines(DenseMatrix fm, PairIntArray set1,
        PairIntArray set2, 
        Image img1, Image img2, int image1Width,
        int image1Height, int image2Width, int image2Height, String outfileNumber) 
        throws IOException {

        EpipolarTransformer st = new EpipolarTransformer();
        
        DenseMatrix input1 =
            Util.rewriteInto3ColumnMatrix(set1);

        DenseMatrix input2 =
            Util.rewriteInto3ColumnMatrix(set2);

        for (int ii = 0; ii < input1.numColumns(); ii++) {
            double x = input1.get(0, ii);
            double y = input1.get(1, ii);
            ImageIOHelper.addPointToImage((float) x, (float) y, img1, 3,
                255, 0, 0);
        }
        for (int ii = 0; ii < input2.numColumns(); ii++) {
            double x2 = input2.get(0, ii);
            double y2 = input2.get(1, ii);
            ImageIOHelper.addPointToImage((float) x2, (float) y2, img2, 3,
                255, 0, 0);
        }

        EpipolarTransformer spTransformer = new EpipolarTransformer();

        Color clr = null;
        for (int ii = 0; ii < input2.numColumns(); ii++) {
            clr = getColor(clr);
            DenseMatrix epipolarLinesInLeft = 
                MatrixUtil.multiply(algorithms.matrix.MatrixUtil.transpose(fm), input2);
            PairIntArray leftLine = spTransformer.getEpipolarLine(
                epipolarLinesInLeft, image1Width, image1Height, ii);
            ImageIOHelper.addCurveToImage(leftLine, img1, 0,
                clr.getRed(), clr.getGreen(), clr.getBlue());
        }

        clr = null;
        for (int ii = 0; ii < input1.numColumns(); ii++) {
            clr = getColor(clr);
            DenseMatrix epipolarLinesInRight = MatrixUtil.multiply(fm, input1);
            PairIntArray rightLine = spTransformer.getEpipolarLine(
                epipolarLinesInRight, img2.getWidth(), img2.getHeight(), ii);
            ImageIOHelper.addCurveToImage(rightLine, img2, 0,
                clr.getRed(), clr.getGreen(), clr.getBlue());
        }

        String dirPath = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(
            dirPath + "/tmp_m_1_" + outfileNumber + ".png", img1);
        ImageIOHelper.writeOutputImage(
            dirPath + "/tmp_m_2_" + outfileNumber + ".png", img2);
    }

    private void overplotKeypoints(int[][] keypoints1, int[][] keypoints2,
            Image img1, Image img2, String outfileNumber)
            throws IOException {

        int ii;
        for (ii = 0; ii < keypoints1.length; ii++) {
            ImageIOHelper.addPointToImage(keypoints1[ii][1], keypoints1[ii][0], img1, 3,
                    255, 0, 0);
        }
        for (ii = 0; ii < keypoints2.length; ii++) {
            ImageIOHelper.addPointToImage(keypoints2[ii][1], keypoints2[ii][0], img2, 3,
                    255, 0, 0);
        }

        String dirPath = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(
                dirPath + "/m_1_harriscorners_" + outfileNumber + ".png", img1);
        ImageIOHelper.writeOutputImage(
                dirPath + "/m_2_harriscorners_" + outfileNumber + ".png", img2);
    }
    
    private Color getColor(Color clr) {
        if ((clr == null) || clr.equals(Color.MAGENTA)) {
            return Color.BLUE;
        }
        if (clr.equals(Color.BLUE)) {
            return Color.PINK;
        } else if (clr.equals(Color.PINK)) {
            return Color.GREEN;
        } else if (clr.equals(Color.GREEN)) {
            return Color.RED;
        } else if (clr.equals(Color.RED)) {
            return Color.CYAN;
        } else if (clr.equals(Color.CYAN)) {
            return Color.MAGENTA;
        } else if (clr.equals(Color.MAGENTA)) {
            return Color.LIGHT_GRAY;
        } else {
            return Color.ORANGE;
        }
    }
    
    protected void getMertonCollege10TrueMatches(PairIntArray left, 
        PairIntArray right) {
        
        /*
        58, 103   32, 100
        486, 46   474, 49
        845, 127  878, 151
        949, 430  998, 471
        541, 428  533, 460
        225, 453  213, 498
        49, 509   21, 571
        373, 239  365, 258
        737, 305  762, 335
        84, 273   60, 298
        */
        
        left.add(58, 103);  right.add(32, 100);
        left.add(486, 46);   right.add(474, 49);
        left.add(845, 127);   right.add(878, 151);
        left.add(949, 430);   right.add(998, 471);
        left.add(541, 428);   right.add(533, 460);
        left.add(225, 453);   right.add(213, 498);
        left.add(49, 509);   right.add(21, 571);
        left.add(373, 239);   right.add(365, 258);
        left.add(737, 305);   right.add(762, 335);
        left.add(84, 273);   right.add(60, 298);
    }

    protected int[][] indexesOfSeveralTrueMatchesMerton(
            List<PairInt> keypoints1, List<PairInt> keypoints2) {

        List<PairInt> m1 = new ArrayList<PairInt>();
        List<PairInt> m2 = new ArrayList<PairInt>();
        populateSeveralTrueMatchesMerton(m1, m2);
        assertEquals(m1.size(), m2.size());

        int[][] indexes = new int[m1.size()][2];

        int i;
        int idx1, idx2;
        for(i = 0; i < m1.size(); ++i) {
            idx1 = keypoints1.indexOf(m1.get(i));
            idx2 = keypoints2.indexOf(m2.get(i));

            //assert(idx1 > -1);
            //assert(idx2 > -1);
            indexes[i] = new int[]{idx1, idx2};
            System.out.printf("m1m2 [idx %d:%d], (row,col  %d,%d  %d,%d)\n",
                    idx1, idx2,
                    m1.get(i).getX(), m1.get(i).getY(),
                    m2.get(i).getX(), m2.get(i).getY());
        }
        return indexes;
    }

    protected void populateSeveralTrueMatchesMerton(
        List<PairInt> keypoints1, List<PairInt> keypoints2
    ) {
        
        /*
        extracted from list of local max filtered harris corners w/ descriptor matches bettern
        than 2.5 * st. dev of all cost combinations.

        #row, col
        58,449
        195,79
        418,430
        487,356
        436,899
        293,770
        70,851
        249,542
        */
        keypoints1.add(new PairInt(58,449));
        keypoints1.add(new PairInt(195,79));
        keypoints1.add(new PairInt(418,430));
        keypoints1.add(new PairInt(487,356));
        keypoints1.add(new PairInt(436,899));
        keypoints1.add(new PairInt(293,770));
        keypoints1.add(new PairInt(70,851));
        keypoints1.add(new PairInt(249,542));

        /*
        #row, col
        60,430
        207,55
        455,424
        532,350
        478,943
        326,806
        92,885
        270,537
         */
        keypoints2.add(new PairInt(60,430));
        keypoints2.add(new PairInt(207,55));
        keypoints2.add(new PairInt(455,424));
        keypoints2.add(new PairInt(532,350));
        keypoints2.add(new PairInt(478,943));
        keypoints2.add(new PairInt(326,806));
        keypoints2.add(new PairInt(92,885));
        keypoints2.add(new PairInt(270,537));
    }

    protected void getMertonCollege7TrueMatches(PairIntArray left,
                                                PairIntArray right) {

        /*
        58, 103   32, 100
        486, 46   474, 49
        845, 127  878, 151
        949, 430  998, 471
        541, 428  533, 460
        225, 453  213, 498
        49, 509   21, 571

        */

        left.add(58, 103);  right.add(32, 100);
        left.add(486, 46);   right.add(474, 49);
        left.add(845, 127);   right.add(878, 151);
        left.add(949, 430);   right.add(998, 471);
        left.add(541, 428);   right.add(533, 460);
        left.add(225, 453);   right.add(213, 498);
        left.add(49, 509);   right.add(21, 571);

    }

    private double[][] convertToDouble(PairIntArray a) {
        int n = a.getN();
        double[][] out = new double[3][n];
        for (int i = 0; i < 3; ++i) {
            out[i] = new double[n];
        }
        Arrays.fill(out[2], 1);
        for (int i = 0; i < n; ++i) {
            out[0][i] = a.getX(i);
            out[1][i] = a.getY(i);
        }
        return out;
    }

    private double[][] convertToDouble(float[][] a) {
        double[][] c = new double[a.length][];
        int i, j;
        for (i = 0; i < a.length; ++i) {
            c[i] = new double[a[0].length];
            for (j = 0; j < a[0].length; ++j) {
                c[i][j] = a[i][j];
            }
        }
        return c;
    }
    
    public void testIsDegenerate() {

        /* euler transformations
        
        about z-axis (yaw):           about x-axis (roll):       about the y-axis (pitch):
            | cos φ   sin φ    0 |    |     1       0       0 |  |  cos ψ  sin ψ    0 |
            |-sin φ   cos φ    0 |    |     0   cos θ   sin θ |  | -sin ψ  cos ψ    0 |
            |     0       0    1 |    |     0  -sin θ   cos θ |  |      0      0    1 |
        
        */
        /*               
             (0,15,1)
               (2.5,10,1)  (7.5,10,1)
                    (5, 5, 1)                   
        
        */
        
        double[][] x1 = new double[3][4];
        x1[0] = new double[]{5, 2.5,   7.5, 0};// x's
        x1[1] = new double[]{5,  10,   10, 15};// y's
        x1[2] = new double[]{1,   1,    1,  1};// z's
        
        double psi = 10.*Math.PI/180;
        double[][] tPitch = new double[3][3];
        tPitch[0] = new double[]{Math.cos(psi),  Math.sin(psi), 0};
        tPitch[1] = new double[]{-Math.sin(psi), Math.cos(psi), 0};
        tPitch[2] = new double[]{0,                          0, 1};
        
        double tx = 100;
        double[][] tTransX = new double[3][3];
        tTransX[0] = new double[]{1, 0, tx};
        tTransX[1] = new double[]{0, 1,  0};
        tTransX[2] = new double[]{0, 0,  1};
        
        double[][] x2 = MatrixUtil.multiply(tPitch, x1);
        x2 = MatrixUtil.multiply(tTransX, x2);        
        
        DenseMatrix x1M = new DenseMatrix(x1);
        DenseMatrix x2M = new DenseMatrix(x2);
        
        boolean d = MiscMath0.areColinear(x1, 1E-6) ||
                MiscMath0.areColinear(x2, 1E-6);
        assertTrue(d);
    
        /*               
                     (5,15,1)
               (2.5,10,1)  (7.5,10,1)
                    (5, 5, 1)                   
        
        */
        x1[0] = new double[]{5, 2.5,   7.5, 5};// x's
        x1[1] = new double[]{5,  10,   10, 15};// y's
        x1[2] = new double[]{1,   1,    1,  1};// z's

        d = MiscMath0.areColinear(x1, 1E-6);
        assertFalse(d);
    }

    public static void main(String[] args) {
        
        try {
            EpipolarTransformer test = new EpipolarTransformer();
            
            //test.testRANSAC();
                        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
    }

}
