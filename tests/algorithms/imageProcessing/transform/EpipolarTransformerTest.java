package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.matrix.MatrixUtil;
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
    
    public void testNormalization() throws Exception {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        log.info("SEED=" + seed);
        sr.setSeed(seed);
        
        EpipolarTransformer sTransformer = new EpipolarTransformer();
        
        int w = 300;
        int h = 400;
        
        double eps = 0.00001;
        
        double sqrtTwo = Math.sqrt(2);
        /*
        x_transformed = xc*s + ((x - xc)*s) + tX = x*s - xc
        y_transformed = yc*s + ((y - yc)*s) + tY = y*s - yc
        */
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
