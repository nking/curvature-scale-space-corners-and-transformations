package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.util.ResourceFinder;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import java.awt.Color;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;
import org.ejml.simple.*;

/**
 *
 * @author nichole
 */
public class EpipolarTransformerTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public EpipolarTransformerTest() {
    }
    
    public void testCreateScaleTranslationMatrix() throws Exception {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        log.info("SEED=" + seed);
        sr.setSeed(seed);
        
        EpipolarTransformer sTransformer = new EpipolarTransformer();
        
        int w = 300;
        int h = 400;
        double scale = Math.sqrt(2.)/150.;
        double centroidX = (w/2);
        double centroidY = (h/2);
        
        SimpleMatrix tMatrix = sTransformer.createScaleTranslationMatrix(scale, 
            centroidX, centroidY);
        
        double sqrtTwo = Math.sqrt(2);
        /*
        x_transformed = xc*s + ((x - xc)*s) + tX = x*s - xc
        y_transformed = yc*s + ((y - yc)*s) + tY = y*s - yc
        */
        int nTests = 1;
        int nPoints = 100;
        SimpleMatrix xy = new SimpleMatrix(3, nPoints);
        double[][] xyExpectedTr = new double[2][nPoints];
        for (int j = 0; j < 2; ++j) {
            xyExpectedTr[j] = new double[nPoints];
        }
        for (int i = 0; i < nTests; ++i) {
            for (int j = 0; j < nPoints; ++j) {
                int x = sr.nextInt(w);
                int y = sr.nextInt(h);
                double xt = (x*scale) + -centroidX*scale;
                double yt = (y*scale) + -centroidY*scale;
                xy.set(0, j, x);
                xy.set(1, j, y);
                xy.set(2, j, 1);
                xyExpectedTr[0][j] = xt;
                xyExpectedTr[1][j] = yt;
                //System.out.println(
                //    String.format("(%.5f, %.5f) --> (%.5f, %.5f)", 
                //    (float)x, (float)y, (float)xt, (float)yt));
            }
            
            double[][] xyTransformed = MatrixUtil.dot(tMatrix, xy);
        
            double avgDist = 0;
            for (int j = 0; j < nPoints; ++j) {
                double xt = xyTransformed[0][j];
                double yt = xyTransformed[1][j];
                assertTrue(Math.abs(xt - xyExpectedTr[0][j]) < 0.01);
                assertTrue(Math.abs(yt - xyExpectedTr[1][j]) < 0.01);
                avgDist += Math.sqrt(xt*xt + yt*yt);
            }
            avgDist /= (double)nPoints;
            assertTrue(avgDist <= sqrtTwo);
        }
    }
    
    public void testMoreThan7Points() throws Exception {
            
        PairIntArray leftTrueMatches = new PairIntArray();
        PairIntArray rightTrueMatches = new PairIntArray();
        getMertonCollege10TrueMatches(leftTrueMatches, rightTrueMatches);
        
        EpipolarTransformer spTransformer = new EpipolarTransformer();
        
        SimpleMatrix fm = spTransformer.calculateEpipolarProjection(
            leftTrueMatches, rightTrueMatches);

        assertNotNull(fm);
        
        double tolerance = 3;
        
        SimpleMatrix matchedLeftXY = spTransformer.rewriteInto3ColumnMatrix(leftTrueMatches);
        SimpleMatrix matchedRightXY = spTransformer.rewriteInto3ColumnMatrix(rightTrueMatches);
        EpipolarTransformationFit fit = spTransformer.calculateError(fm, matchedLeftXY, 
            matchedRightXY, ErrorType.DIST_TO_EPIPOLAR_LINE, tolerance);
        
        assertNotNull(fit);
        assertNotNull(fit.getErrorType());
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
        
        EpipolarTransformationFit fit2 = spTransformer.calculateError(fm, 
            matchedLeftXY, matchedRightXY, ErrorType.SAMPSONS, tolerance);
        
        assertNotNull(fit2);
        assertNotNull(fit2.getErrorType());
        assertEquals(10, fit2.getErrors().size());
        assertEquals(10, fit2.getInlierIndexes().size());
        fit2.calculateErrorStatistics();
        double mean2 = fit2.getMeanError();
        double stdev2 = fit2.getStDevFromMean();
        List<Double> errors2 = fit2.getErrors();
        
        assertTrue(mean2 < tolerance);
        assertTrue(stdev2 < tolerance);
        for (Double error : errors2) {
            assertTrue(error < tolerance);
        }
    }
    
    public void test7Points() throws Exception {
            
        PairIntArray leftTrueMatches = new PairIntArray();
        PairIntArray rightTrueMatches = new PairIntArray();
        getMertonCollege7TrueMatches(leftTrueMatches, rightTrueMatches);
        
        EpipolarTransformer spTransformer = new EpipolarTransformer();
        
        List<SimpleMatrix> fms =
            spTransformer.calculateEpipolarProjectionFor7Points(leftTrueMatches, 
            rightTrueMatches);

        assertNotNull(fms);
        assertFalse(fms.isEmpty());
        
        double tolerance = 3;
        
        SimpleMatrix matchedLeftXY = spTransformer.rewriteInto3ColumnMatrix(leftTrueMatches);
        SimpleMatrix matchedRightXY = spTransformer.rewriteInto3ColumnMatrix(rightTrueMatches);
        EpipolarTransformationFit fit = null;
        for (SimpleMatrix fm : fms) {
            EpipolarTransformationFit fitI = 
                spTransformer.calculateError(fm, matchedLeftXY, 
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
        
        assertTrue(mean < 0.001);
        assertTrue(stdev < 0.001);
        for (Double error : errors) {
            assertTrue(error < 0.001);
        }
        
        EpipolarTransformationFit fit2 = null;
        for (SimpleMatrix fm : fms) {
            EpipolarTransformationFit fitI = 
                spTransformer.calculateError(fm, matchedLeftXY, 
                    matchedRightXY, ErrorType.SAMPSONS, tolerance);
            if (fitI.isBetter(fit2)) {
                fit2 = fitI;
            }
        }
        assertNotNull(fit2);
        assertNotNull(fit2.getErrorType());
        assertEquals(7, fit2.getErrors().size());
        assertEquals(7, fit2.getInlierIndexes().size());
        fit2.calculateErrorStatistics();
        double mean2 = fit2.getMeanError();
        double stdev2 = fit2.getStDevFromMean();
        List<Double> errors2 = fit2.getErrors();
        
        assertTrue(mean2 < 0.001);
        assertTrue(stdev2 < 0.001);
        for (Double error : errors2) {
            assertTrue(error < 0.001);
        }
    }
    
    /*
    for more datasets:
    http://www.robots.ox.ac.uk/~vgg/data/data-mview.html
    */
    
    private void overplotEpipolarLines(SimpleMatrix fm, PairFloatArray set1,
        PairFloatArray set2, Image img1, Image img2, int image1Width,
        int image1Height, int image2Width, int image2Height, String outfileNumber) 
        throws IOException {

        EpipolarTransformer st = new EpipolarTransformer();
        
        SimpleMatrix input1 =
            st.rewriteInto3ColumnMatrix(set1);

        SimpleMatrix input2 =
            st.rewriteInto3ColumnMatrix(set2);

        for (int ii = 0; ii < input1.numCols(); ii++) {
            double x = input1.get(0, ii);
            double y = input1.get(1, ii);
            ImageIOHelper.addPointToImage((float) x, (float) y, img1, 3,
                255, 0, 0);
        }
        for (int ii = 0; ii < input2.numCols(); ii++) {
            double x2 = input2.get(0, ii);
            double y2 = input2.get(1, ii);
            ImageIOHelper.addPointToImage((float) x2, (float) y2, img2, 3,
                255, 0, 0);
        }

        EpipolarTransformer spTransformer = new EpipolarTransformer();

        Color clr = null;
        for (int ii = 0; ii < input2.numCols(); ii++) {
            clr = getColor(clr);
            SimpleMatrix epipolarLinesInLeft = fm.transpose().mult(input2);
            PairIntArray leftLine = spTransformer.getEpipolarLine(
                epipolarLinesInLeft, image1Width, image1Height, ii);
            ImageIOHelper.addCurveToImage(leftLine, img1, 0,
                clr.getRed(), clr.getGreen(), clr.getBlue());
        }

        clr = null;
        for (int ii = 0; ii < input1.numCols(); ii++) {
            clr = getColor(clr);
            SimpleMatrix epipolarLinesInRight = fm.mult(input1);
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
