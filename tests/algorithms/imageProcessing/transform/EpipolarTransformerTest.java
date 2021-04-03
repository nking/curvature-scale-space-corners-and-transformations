package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.transform.EpipolarTransformer.NormalizationTransformations;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath;
import algorithms.util.FormatArray;
import algorithms.util.ResourceFinder;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import java.awt.Color;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
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
        double sqrt2 = Math.sqrt(2);
        
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
        System.out.printf("rmsLeft=%.4e  rmsRight=%.4e\n sqrt(2)=%.4e\n", rmsLeft, rmsRight, sqrt2);
        
        assertTrue(Math.abs(meanAndStDevLeftX[0]) < tol);
        assertTrue(Math.abs(meanAndStDevLeftY[0]) < tol);
        assertTrue(Math.abs(rmsLeft - sqrt2) < tol2);
        
        assertTrue(Math.abs(meanAndStDevRightX[0]) < tol);
        assertTrue(Math.abs(meanAndStDevRightY[0]) < tol);
        assertTrue(Math.abs(rmsRight - sqrt2) < tol2);
      
        
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
    
     
    public void testMoreThan7Points() throws Exception {
            
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
            leftM, rightM);

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
        System.out.println("nCorr=" + fit.getInlierIndexes().size());
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
    
    public void test7Points() throws Exception {
            
        PairIntArray leftTrueMatches = new PairIntArray();
        PairIntArray rightTrueMatches = new PairIntArray();
        getMertonCollege7TrueMatches(leftTrueMatches, rightTrueMatches);
                
        double[][] left = convertToDouble(leftTrueMatches);
        double[][] right = convertToDouble(leftTrueMatches);
        DenseMatrix leftM = new DenseMatrix(left);
        DenseMatrix rightM = new DenseMatrix(right);
        
        EpipolarTransformer.NormalizedXY normXY1 = EpipolarTransformer.normalize(leftM);
        EpipolarTransformer.NormalizedXY normXY2 = EpipolarTransformer.normalize(rightM);
        
        EpipolarTransformer spTransformer = new EpipolarTransformer();

        List<DenseMatrix> normalizedFMs 
            = spTransformer.calculateEpipolarProjectionFor7Points(
            normXY1.getXy(), normXY2.getXy());

        assertNotNull(normalizedFMs);
        assertFalse(normalizedFMs.isEmpty());
                
        boolean useToleranceAsStatFactor = false;
        final double tolerance = 3;
        
        ErrorType errorType = ErrorType.SAMPSONS;
        Distances distances = new Distances();
        
        List<DenseMatrix> denormalizedFMs = new ArrayList<DenseMatrix>();
        List<EpipolarTransformationFit> fits = new ArrayList<EpipolarTransformationFit>();
        
        EpipolarTransformationFit bestFit = null;
        int bestFitIdx = -1;
        
        for (DenseMatrix nfm : normalizedFMs) {
            DenseMatrix fm = EpipolarTransformer.denormalizeTheFundamentalMatrix(
                nfm, normXY1.getNormalizationMatrices(),
                normXY2.getNormalizationMatrices());
            denormalizedFMs.add(fm);
            
            EpipolarTransformationFit fit = null;
            if (useToleranceAsStatFactor) {
                fit = distances.calculateError2(nfm,
                        leftM, rightM, errorType, tolerance);
            } else {
                fit = distances.calculateError(nfm,
                        leftM, rightM, errorType, tolerance);
            }
            
            fits.add(fit);
            
            if (fit.isBetter(bestFit)) {
                bestFit = fit;
                bestFitIdx = fits.size() - 1;
            }
        }
        
        assertNotNull(bestFit);
        assertNotNull(bestFit.getErrorType());
        assertEquals(7, bestFit.getErrors().size());
        assertEquals(7, bestFit.getInlierIndexes().size());
        bestFit.calculateErrorStatistics();
        double mean = bestFit.getMeanError();
        double stdev = bestFit.getStDevFromMean();
        List<Double> errors = bestFit.getErrors();
        
        assertTrue(mean < 0.01);
        assertTrue(stdev < 0.01);
        for (Double error : errors) {
            assertTrue(error < 0.01);
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

        overplotEpipolarLines(denormalizedFMs.get(bestFitIdx), 
            leftTrueMatches, rightTrueMatches,
            img1, img2, 
            image1Width, image1Height, image2Width, image2Height, 
            "e_merton_college_7pts" + Integer.valueOf(0).toString()); 
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
                MatrixUtil.multiply(
                algorithms.matrix.MatrixUtil.transpose(fm), input2);
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
        
        boolean d = EpipolarTransformer.isDegenerate(x1M, x2M);
        assertTrue(d);
    
        /*               
                     (5,15,1)
               (2.5,10,1)  (7.5,10,1)
                    (5, 5, 1)                   
        
        */
        x1[0] = new double[]{5, 2.5,   7.5, 5};// x's
        x1[1] = new double[]{5,  10,   10, 15};// y's
        x1[2] = new double[]{1,   1,    1,  1};// z's
        
        x2 = MatrixUtil.multiply(tPitch, x1);
        x2 = MatrixUtil.multiply(tTransX, x2);        
        
        x1M = new DenseMatrix(x1);
        x2M = new DenseMatrix(x2);
        
        d = EpipolarTransformer.isDegenerate(x1M, x2M);
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
