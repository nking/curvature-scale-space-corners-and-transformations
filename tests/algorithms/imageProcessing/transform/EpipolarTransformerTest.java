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
          u1_normalized = T1 * u1
          u2_normalized = T2 * u2

          denormalized u1 = T1^-1 * u1_normalized and similar foru2

          FM_normalized = inverse(transpose(T2)) * FM * inverse(T1)

          Revisiting Hartley’s Normalized Eight-Point Algorithm
          Chojnacki et al. 2003

        * denormalized FM = transpose(T2) * FM_normalized * T1
        
          print condition number:  largest eigenvalue/2nd smallest eigenvalue.
          when condition numbr is arge, the 2 smallest eigenvalues are close to on another
            and that makes their eigenvectors sensitive to small pertubations.

          u2^T * FM * u1 = u2_normalized^T * FM_normalized * u1_normalized = residual

              | 1/s   0  0 |   | 1  0  -xc |   | 1/s    0   -s*xc |
          T = |  0  1/s  0 | * | 0  1  -yc | = |   0  1/s   -s*yc |
              |  0    0  1 |   | 0  0   1  |   |   0    0      1  |

                 | 1  0  xc |   | s  0   0 |
          T^-1 = | 0  1  yc | * | 0  s   0 |
                 | 0  0   1 |   | 0  0   1 |
   
                        | 1  0  xc |   | s^2  0    0 |   | 1  0  0 |
          T^-1 * T^-T = | 0  1  yc | * | 0   s^2   0 | * | 0  1  0 |
                        | 0  0   1 |   | 0    0    1 |   | xc yc 1 |

                        | s^2 + xc^2   xc*yc       xc |
                      = | yc*xc        s^2 + yc^2  yc |
                        | xc           yc          1  |
        */
       
        //---- test a fundamental matrix.
        //       the values for x1, y1, x2, and y2 in the fundamental matrix can be tested
        /*
        The fundamental matrix terms solved for by the orthogonal component of best fit to A:
            A[i][0] = x1 * x2; <— FM[0][0]
            A[i][1] = x1 * y2; <— FM[1][0]
            A[i][2] = x1;      <— FM[2][0]  <======
            A[i][3] = y1 * x2; <— FM[0][1]
            A[i][4] = y1 * y2; <— FM[1][1]
            A[i][5] = y1;      <— FM[2][1]  <======
            A[i][6] = x2;      <— FM[0][2]  <======
            A[i][7] = y2;      <— FM[1][2]  <======
            A[i][8] = 1;         
        */
        
        /*
        hartley 1997 house 
        [junit] expected de-normalized FM=
    [junit] 4.144e-06, -6.231e-05, 2.817e-02 
    [junit] 2.684e-05, -4.437e-07, -3.188e-01 
    [junit] -3.852e-02, 3.274e-01, 1.000e+00
        
        before de-normalization:
        [junit] 9.773e-02, 2.826e+00, 1.605e+00
    [junit] -2.847e+00, -1.889e-01, -3.970e+00
    [junit] -1.521e+00, 3.753e+00, -6.851e-02
        */
        
    }
    
    /*
    public void _testNormalization() throws Exception {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        log.info("SEED=" + seed);
        sr.setSeed(seed);
        
        EpipolarTransformer sTransformer = new EpipolarTransformer();
        
        int w = 300;
        int h = 400;
        
        double eps = 0.00001;
        
        double sqrtTwo = Math.sqrt(2);
        
        //x_transformed = xc*s + ((x - xc)*s) + tX = x*s - xc
        //y_transformed = yc*s + ((y - yc)*s) + tY = y*s - yc
        
        int nTests = 1;
        int nPoints = 100;
        
        DenseMatrix xy = new DenseMatrix(3, nPoints);
        
        for (int i = 0; i < nTests; ++i) {
            
            for (int j = 0; j < nPoints; ++j) {
                int x = sr.nextInt(w);
                int y = sr.nextInt(h);
                xy.set(0, j, x);
                xy.set(1, j, y);
                xy.set(2, j, 1);
            }
          
            EpipolarTransformer.NormalizedXY normXY = EpipolarTransformer.normalize(xy);
            
            DenseMatrix xy2 = normXY.getXy();
            assertEquals(nPoints, xy2.numColumns());
            assertEquals(3, xy2.numRows());            
                   
            double avgDist = 0;
            for (int j = 0; j < nPoints; ++j) {
                double xt = xy2.get(0, j);
                double yt = xy2.get(1, j);
                // trnsformed center is 0,0 so their coords are their distances
                avgDist += Math.sqrt(xt*xt + yt*yt);
            }
            avgDist /= (double)nPoints;
            //System.out.println("mean dist^2=" + avgDist);
            assertTrue(avgDist <= (sqrtTwo + eps));
                       
            // 3XN      3X3        3XN
            //normXY = tMatrix dot xy
            //  normXY * inv(tMatrix) = xy
            DenseMatrix invT = algorithms.imageProcessing.util.MatrixUtil.inverse(
                normXY.getNormalizationMatrix());
            
            DenseMatrix denorm = 
                MatrixUtil.multiply(invT, normXY.getXy());
            
            for (int j = 0; j < nPoints; ++j) {
                double x = xy.get(0, j);
                double y = xy.get(1, j);
                
                double x2 = denorm.get(0, j);
                double y2 = denorm.get(1, j);
                
                assertTrue(Math.abs(x - x2) < eps);
                assertTrue(Math.abs(y - y2) < eps);
                
                //System.out.println("x=" + x + " y=" + y);
                //System.out.println("    x2=" + x2 + " y2=" + y2);
            }
        }
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
            normalizedFM, normXY1.getNormalizationMatrix(),
            normXY2.getNormalizationMatrix());
                
        
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
                nfm, normXY1.getNormalizationMatrix(),
                normXY2.getNormalizationMatrix());
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
    */
    
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
