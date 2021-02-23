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
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.SVD;
import org.junit.Test;

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
          
            EpipolarTransformer.NormalizedXY normXY = sTransformer.normalize(xy);
            
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
        Distances distances = new Distances();
        
        DenseMatrix fm = spTransformer.calculateEpipolarProjection(
            leftTrueMatches, rightTrueMatches);

        assertNotNull(fm);
        
        double tolerance = 3;
        
        DenseMatrix matchedLeftXY = Util.rewriteInto3ColumnMatrix(leftTrueMatches);
        DenseMatrix matchedRightXY = Util.rewriteInto3ColumnMatrix(rightTrueMatches);
        
        EpipolarTransformationFit fit = distances.calculateError(fm, matchedLeftXY, 
            matchedRightXY, ErrorType.DIST_TO_EPIPOLAR_LINE, tolerance);
        
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
        
    }
    
    public void test7Points() throws Exception {
            
        PairIntArray leftTrueMatches = new PairIntArray();
        PairIntArray rightTrueMatches = new PairIntArray();
        getMertonCollege7TrueMatches(leftTrueMatches, rightTrueMatches);
        
        EpipolarTransformer spTransformer = new EpipolarTransformer();
        
        Distances distances = new Distances();
        
        List<DenseMatrix> fms =
            spTransformer.calculateEpipolarProjectionFor7Points(leftTrueMatches, 
            rightTrueMatches);

        assertNotNull(fms);
        assertFalse(fms.isEmpty());
        
        double tolerance = 3;
        
        DenseMatrix matchedLeftXY = Util.rewriteInto3ColumnMatrix(leftTrueMatches);
        DenseMatrix matchedRightXY = Util.rewriteInto3ColumnMatrix(rightTrueMatches);
        EpipolarTransformationFit fit = null;
        for (DenseMatrix fm : fms) {
            EpipolarTransformationFit fitI = 
                distances.calculateError(fm, matchedLeftXY, 
                    matchedRightXY, ErrorType.DIST_TO_EPIPOLAR_LINE, tolerance);
            if (fitI.isBetter(fit)) {
                fit = fitI;
            }
        }
        assertNotNull(fit);
        assertNotNull(fit.getErrorType());
        assertEquals(7, fit.getErrors().size());
        assertEquals(7, fit.getInlierIndexes().size());
        fit.calculateErrorStatistics();
        double mean = fit.getMeanError();
        double stdev = fit.getStDevFromMean();
        List<Double> errors = fit.getErrors();
        
        assertTrue(mean < 0.01);
        assertTrue(stdev < 0.01);
        for (Double error : errors) {
            assertTrue(error < 0.01);
        }
        
    }
    
    /*
    for more datasets:
    http://www.robots.ox.ac.uk/~vgg/data/
    */
    
    private void overplotEpipolarLines(DenseMatrix fm, PairFloatArray set1,
        PairFloatArray set2, Image img1, Image img2, int image1Width,
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
