package algorithms.imageProcessing.features;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.EpipolarTransformer;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.awt.Color;
import java.io.IOException;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author nichole
 */
public class RANSACSolverTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public RANSACSolverTest() {
    }

    public void estRANSAC() throws Exception {
    
        //TODO: add the reference for this data here.
        
        PairIntArray leftTrueMatches = new PairIntArray();
        PairIntArray rightTrueMatches = new PairIntArray();
        getMertonCollege10TrueMatches(leftTrueMatches, rightTrueMatches);
        
        PairIntArray leftFalseMatches = new PairIntArray();
        PairIntArray rightFalseMatches = new PairIntArray();
        getMertonCollegeFalseMatch1(leftFalseMatches, rightFalseMatches);
        getMertonCollegeFalseMatch2(leftFalseMatches, rightFalseMatches);
        getMertonCollegeFalseMatch3(leftFalseMatches, rightFalseMatches);
        
        PairIntArray leftTruePlusFalse = leftTrueMatches.copy();
        PairIntArray rightTruePlusFalse = rightTrueMatches.copy();
        getMertonCollegeFalseMatch1(leftTruePlusFalse, rightTruePlusFalse);
        getMertonCollegeFalseMatch2(leftTruePlusFalse, rightTruePlusFalse);
        getMertonCollegeFalseMatch3(leftTruePlusFalse, rightTruePlusFalse);
        
        PairIntArray outputLeft = new PairIntArray(); 
        PairIntArray outputRight = new PairIntArray();
        
        RANSACSolver solver = new RANSACSolver();
        
        EpipolarTransformationFit fit = solver.calculateEpipolarProjection(
            //leftTruePlusFalse, rightTruePlusFalse, outputLeft, outputRight);
            leftTrueMatches, rightTrueMatches, outputLeft, outputRight);
        
        assertNotNull(fit);
        
        log.info("fit=" + fit.toString());
        log.info(" fm=" + fit.getFundamentalMatrix().toString());

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

        overplotEpipolarLines(fit.getFundamentalMatrix(), outputLeft, outputRight,
            img1, img2, 
            image1Width, image1Height, image2Width, image2Height, 
            Integer.valueOf(0).toString()); 
        
        int n = outputLeft.getN();
        
        log.info("leftTrueMatches=" + leftTrueMatches.getN());
        log.info("rightTrueMatches=" + rightTrueMatches.getN());
        log.info("outputRight=" + outputRight.getN());
        
        /*
        solution when all 10 are given and the first randomly chosen among
        them are used:
                  -0.000  -0.000   0.000  
          [junit]  0.000   0.000   0.036  
          [junit] -0.002  -0.045   1.000 
        */
        
        List<Double> errors = fit.getErrors();
        assertTrue(leftTrueMatches.getN() == outputLeft.getN());
        
        for (double error : errors) {
            assertTrue(error < 3);
        }
        
    }
    
    public void testSampsonErrors_stereo_01() {

        // TODO: needs corrections for normalization
        
        // 1024 X 768
        String fileName1 = "merton_college_I_001.jpg";
        String fileName2 = "merton_college_I_002.jpg";
        
        PairIntArray m1 = new PairIntArray(7);
        PairIntArray m2 = new PairIntArray(7);
        populateWithMertonMatches7(m1, m2);
        
        // in progress...
        
        PairIntArray out1 = new PairIntArray(m1.getN());
        PairIntArray out2 = new PairIntArray(m1.getN());
        RANSACSolver solver = new RANSACSolver();
        EpipolarTransformationFit fit
            = solver.calculateEpipolarProjection(m1, m2, out1, out2);

        assertEquals(m1.getN(), out1.getN());
        assertEquals(m2.getN(), out2.getN());
        
        SimpleMatrix fm = fit.getFundamentalMatrix();
        
        System.out.print("fm.errors=");
        for (Double d : fit.getErrors()) {
            System.out.println(" " + d);
        }
        System.out.println("");
        System.out.println("fm effective tolerance = " 
            + fit.getEffectiveTolerance() + " stdv=" +
            fit.getEffectiveToleranceStandardDeviation());
        
        float tolOffset = fit.getEffectiveTolerance();
        
        EpipolarTransformer eTransformer = new EpipolarTransformer();
        
        SimpleMatrix m1m = eTransformer.rewriteInto3ColumnMatrix(m1);
        SimpleMatrix m2m = eTransformer.rewriteInto3ColumnMatrix(m2);
        
        SimpleMatrix m2EpipolarLines = fm.mult(m1m);
        SimpleMatrix m1EpipolarLines = fm.transpose().mult(m2m);

        // the 2 distances:
        float[] outputDist = new float[2];
        for (int i = 0; i < m1.getN(); ++i) {
            int j = i;
            
            eTransformer.calculatePerpDistFromLines(
                m1m, m2m, 
                m2EpipolarLines, m1EpipolarLines, 
                i, j, outputDist);
                        
            System.out.println(String.format(
                "fm dists=(%.2f)  unnorm dists=(%.2f, %.2f)", 
                fit.getErrors().get(i).floatValue(), outputDist[0], outputDist[1]));
        
            assertTrue(Math.abs(fit.getErrors().get(i).doubleValue()) 
                < 0.001);
            assertTrue(Math.abs(outputDist[0]) < 0.001);
            assertTrue(Math.abs(outputDist[1]) < 0.001);
        }
                
        // some points which are matches, but aren't in the
        // 7 points list:
        /*
        58, 39    31, 26
        325, 43   281, 40
        577, 68   561, 76
        655, 426  665, 462
        */
        
        PairIntArray points1 = new PairIntArray(4);
        PairIntArray points2 = new PairIntArray(4);
        points1.add(58, 39);    points2.add(31, 26);
        points1.add(325, 43);   points2.add(281, 40);
        points1.add(577, 68);   points2.add(561, 76);
        points1.add(655, 426);  points2.add(665, 462);
        // add a few bad matches
        points1.add(58, 39);    points2.add(31, 26 + 10);
        points1.add(325, 43);   points2.add(281, 40 + 10);
        points1.add(577, 68);   points2.add(561, 76 + 10);
        points1.add(655, 426);  points2.add(665, 462 + 10);
        points1.add(58, 39);    points2.add(31, 26 + 100);
        points1.add(325, 43);   points2.add(281, 40 + 100);
        points1.add(577, 68);   points2.add(561, 76 + 100);
        points1.add(655, 426);  points2.add(665, 462 + 100);
        
        SimpleMatrix p1m = eTransformer.rewriteInto3ColumnMatrix(points1);
        SimpleMatrix p2m = eTransformer.rewriteInto3ColumnMatrix(points2);
        
        SimpleMatrix p2EpipolarLines = fm.mult(p1m);
        SimpleMatrix p1EpipolarLines = fm.transpose().mult(p2m);
        
        for (int i = 0; i < points1.getN(); ++i) {
            int j = i;
            
            eTransformer.calculatePerpDistFromLines(
                p1m, p2m, 
                p2EpipolarLines, p1EpipolarLines, 
                i, j, outputDist);
                        
            System.out.println(String.format(
                "extr: fm p unnorm dists=(%.2f, %.2f)", 
                outputDist[0], outputDist[1]));
        
            //assertTrue(Math.abs(outputDist[0]) < 0.001);
            //assertTrue(Math.abs(outputDist[1]) < 0.001);
        }
        
    }
    
    private void populateWithMertonMatches7(PairIntArray xy1,
        PairIntArray xy2) {
        /*
        104, 275   83, 300
        525, 221   516, 238
        509, 234   500, 252
        843, 80    876, 102
        870, 484   944, 536
        520, 481   533, 524
        88, 504    65, 564
        */
        xy1.add(104, 275);  xy2.add(83, 300);
        xy1.add(525, 221);  xy2.add(516, 238);
        xy1.add(509, 234);  xy2.add(500, 252);
        xy1.add(843, 80);  xy2.add(876, 102);
        xy1.add(870, 484);  xy2.add(944, 536);
        xy1.add(520, 481);  xy2.add(533, 524);
        xy1.add(88, 504);  xy2.add(65, 564);
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
    
    protected void getMertonCollege10TrueMatches(PairFloatArray left, 
        PairFloatArray right) {
        
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
    
    protected void getMertonCollege7TrueMatches(PairFloatArray left, 
        PairFloatArray right) {
        
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
    
    protected void getMertonCollegeFalseMatch1(PairIntArray left, 
        PairIntArray right) {
        //765, 487   753, 552
        left.add(765, 487);   right.add(753, 552);
    }
    protected void getMertonCollegeFalseMatch2(PairIntArray left, 
        PairIntArray right) {
        //253, 141    256, 229
        left.add(253, 141);   right.add(256, 229);
    }
    protected void getMertonCollegeFalseMatch3(PairIntArray left, 
        PairIntArray right) {
        //459, 354  432, 525
        left.add(459, 354);   right.add(432, 525);
    }
    
    private void overplotEpipolarLines(SimpleMatrix fm, PairIntArray set1,
        PairIntArray set2, Image img1, Image img2, int image1Width,
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
}
