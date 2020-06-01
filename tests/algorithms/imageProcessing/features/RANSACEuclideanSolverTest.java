package algorithms.imageProcessing.features;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.transform.Distances;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.EpipolarTransformer;
import algorithms.imageProcessing.transform.EuclideanEvaluator;
import algorithms.imageProcessing.transform.EuclideanTransformationFit;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.transform.Util;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscDebug;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.awt.Color;
import java.io.IOException;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;

/**
 *
 * @author nichole
 */
public class RANSACEuclideanSolverTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public RANSACEuclideanSolverTest() {
    }

    /**
     * this method uses 2 images
     *  merton_college_I_001.jpg and merton_college_I_002.jpg are
        from the Merton College I dataset on:
        http://www.robots.ox.ac.uk/~vgg/data/data-mview.html

     * @throws Exception 
     */
    public void testRANSAC() throws Exception {
            
        System.out.println("testRANSAC");
        
        PairIntArray leftTrueMatches = new PairIntArray();
        PairIntArray rightTrueMatches = new PairIntArray();
        getMertonCollege10TrueMatches(leftTrueMatches, rightTrueMatches);
        
        /*
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
        */
        
        PairIntArray outputLeft = new PairIntArray(); 
        PairIntArray outputRight = new PairIntArray();
        
        RANSACEuclideanSolver solver = new RANSACEuclideanSolver();
        double tolerance = 5;
        EuclideanTransformationFit fit = solver.calculateEuclideanTransformation(
            //leftTruePlusFalse, rightTruePlusFalse, outputLeft, outputRight);
            leftTrueMatches, rightTrueMatches, outputLeft, outputRight,
            tolerance);
        
        assertNotNull(fit);
        log.info("fit=" + fit.toString());
        log.info(" tr=" + fit.getTransformationParameters().toString());

        String fileName1 = "merton_college_I_001.jpg";
        String fileName2 = "merton_college_I_002.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);

        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);

        overplot(fit.getTransformationParameters(), 
            outputLeft, outputRight,
            //leftTrueMatches, rightTrueMatches,
            img1, img2, "reuclidean_" + Integer.valueOf(0).toString()); 
        
        int n = outputLeft.getN();
        
        log.info("leftTrueMatches=" + leftTrueMatches.getN());
        log.info("rightTrueMatches=" + rightTrueMatches.getN());
        log.info("outputRight=" + outputRight.getN());
        
        for (double error : fit.getErrors()) {
            assertTrue(error <= tolerance);
        }
        
    }
    
    /**
     * this method uses 2 images
     *  merton_college_I_001.jpg and merton_college_I_002.jpg are
        from the Merton College I dataset on:
        http://www.robots.ox.ac.uk/~vgg/data/data-mview.html
     * @throws IOException 
     */
    public void testErrors_stereo_01() throws IOException {
        
        // 1024 X 768
        String fileName1 = "merton_college_I_001.jpg";
        String fileName2 = "merton_college_I_002.jpg";
        
        PairIntArray outputLeft = new PairIntArray(); 
        PairIntArray outputRight = new PairIntArray();
        
        PairIntArray m1 = new PairIntArray(7);
        PairIntArray m2 = new PairIntArray(7);
        populateWithMertonMatches7(m1, m2);
        
        RANSACEuclideanSolver solver = new RANSACEuclideanSolver();
        double tolerance = 5;
        EuclideanTransformationFit fit = solver.calculateEuclideanTransformation(
            //leftTruePlusFalse, rightTruePlusFalse, outputLeft, outputRight);
            m1, m2, outputLeft, outputRight,
            tolerance);
        
        assertNotNull(fit);
        log.info("fit=" + fit.toString());
        log.info(" tr=" + fit.getTransformationParameters().toString());

        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);

        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);

        overplot(fit.getTransformationParameters(), 
            outputLeft, outputRight,
            //leftTrueMatches, rightTrueMatches,
            img1, img2, "reuclidean_" + Integer.valueOf(1).toString()); 
        
        int n = outputLeft.getN();
        
        log.info("leftTrueMatches=" + m1.getN());
        log.info("rightTrueMatches=" + m2.getN());
        log.info("outputRight=" + outputRight.getN());
        
        for (double error : fit.getErrors()) {
            assertTrue(error <= tolerance);
        }
        
    }
    
    public void testErrors_panoramic_01() throws IOException {
        
        System.out.println("testErrors_panoramic_01");

        // 617 X 874
        String fileName1 = "brown_lowe_2003_image1.jpg";
        String fileName2 = "brown_lowe_2003_image2.jpg";
        
        PairIntArray m1 = new PairIntArray(7);
        PairIntArray m2 = new PairIntArray(7);
        populateWithBrownAndLoweMatches7(m1, m2);
                
        PairIntArray out1 = new PairIntArray(m1.getN());
        PairIntArray out2 = new PairIntArray(m1.getN());
        RANSACEuclideanSolver solver = new RANSACEuclideanSolver();
        double tolerance = 5;
        EuclideanTransformationFit fit
            = solver.calculateEuclideanTransformation(m1, m2, out1, out2, tolerance);
                
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        
        overplot(fit.getTransformationParameters(), m1, m2,
            img1, img2, 
            "reuclidean_" + Integer.valueOf(2).toString());        
       
        assertTrue(out1.getN() >= (m1.getN()*0.5));
        assertEquals(m1.getN(), m2.getN());
        for (double error : fit.getErrors()) {
            assertTrue(error <= tolerance);
        }
        
        //------------
        // add a couple mroe true matches
        m1.add(377, 138);   m2.add(112, 123);
        m1.add(358, 330);   m2.add(54, 312);
        // add a few bad matches
        m1.add(377, 138);   m2.add(112, 123 + 10);
        m1.add(358, 330);   m2.add(54, 312 + 10);
        m1.add(377, 138);   m2.add(112, 123 + 100);
        m1.add(358, 330);   m2.add(54, 312 + 100);
        out1 = new PairIntArray(m1.getN());
        out2 = new PairIntArray(m1.getN());
        
        fit = solver.calculateEuclideanTransformation(m1, m2, out1, out2, 
            tolerance);

        img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        
        overplot(fit.getTransformationParameters(), m1, m2,
            img1, img2, 
            "reuclidean_b_" + Integer.valueOf(2).toString()); 
        
        for (double error : fit.getErrors()) {
            assertTrue(error <= tolerance);
        }
    }
    
    public void testErrors_moderateProjection_01() throws IOException {
        
        System.out.println("testErrors_moderateProjection_01");

        // 480 X 640
        String fileName1 = "campus_010.jpg";
        String fileName2 = "campus_011.jpg";
        
        PairIntArray m1 = new PairIntArray();
        PairIntArray m2 = new PairIntArray();
        populateWithCampusMatchesMoreThan7(m1, m2);
            
        PairIntArray out1 = new PairIntArray(m1.getN());
        PairIntArray out2 = new PairIntArray(m1.getN());
        
        RANSACEuclideanSolver solver = new RANSACEuclideanSolver();
        double tolerance = 5;
        EuclideanTransformationFit fit
            = solver.calculateEuclideanTransformation(m1, m2, out1, out2, tolerance);
        
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);        
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        
        overplot(fit.getTransformationParameters(), m1, m2,
            img1, img2, 
            "reuclidean_" + Integer.valueOf(3).toString());
                
        System.out.print("fit.errors=");
        for (Double d : fit.getErrors()) {
            System.out.println(" " + d);
        }
        System.out.println("");        
        System.out.flush();
        
        assertTrue(out1.getN() >= (m1.getN()*0.5));
        assertEquals(m1.getN(), m2.getN());
                
    }
    
    public void testErrors_moderateProjection_02() {
        
        System.out.println("testErrors_moderateProjection_02");

        // scaled versions of android 04 and 02
        
        PairIntArray m1 = new PairIntArray(7);
        PairIntArray m2 = new PairIntArray(7);
        populateWithAndroid0402SmallMatches7(m1, m2);
                
        PairIntArray out1 = new PairIntArray(m1.getN());
        PairIntArray out2 = new PairIntArray(m1.getN());
        RANSACEuclideanSolver solver = new RANSACEuclideanSolver();
        double tolerance = 5;
        EuclideanTransformationFit fit
            = solver.calculateEuclideanTransformation(m1, m2, out1, out2, tolerance);

        System.out.print("fit.errors=");
        for (Double d : fit.getErrors()) {
            System.out.println(" " + d);
        }
        System.out.println("");        
        System.out.flush();
        
        assertTrue(out1.getN() >= (m1.getN()*0.5));
        assertEquals(m1.getN(), m2.getN());
        
        // -------------------------                
        // some points which are matches, but aren't in the
        // 7 points list:
        /*
        34, 130   166, 54 
        91, 129   186, 57
        */
        
        PairIntArray points1 = new PairIntArray(4);
        PairIntArray points2 = new PairIntArray(4);
        points1.add(34, 130);   points2.add(166, 54);
        points1.add(91, 129);   points2.add(186, 57);
        // add a few bad matches
        points1.add(34, 130);   points2.add(166, 54 + 10);
        points1.add(91, 129);   points2.add(186, 57 + 10);
        points1.add(34, 130);   points2.add(166, 54 + 100);
        points1.add(91, 129);   points2.add(186, 57 + 100);
        
        out1 = new PairIntArray(m1.getN());
        out2 = new PairIntArray(m1.getN());
        
        fit = solver.calculateEuclideanTransformation(m1, m2, out1, out2, 
            tolerance);
        
        System.out.print("fit.errors=");
        for (Double d : fit.getErrors()) {
            System.out.println(" " + d);
        }
        System.out.println("");        
        System.out.flush();
        
    }
    
    public void populateWithAndroid0402SmallMatches7(PairIntArray xy1,
        PairIntArray xy2) {
        /*
        46, 165  167, 77
        84, 162  181, 81
        63, 149  173, 68
        37, 122  166, 49
        88, 123  183, 52
        54, 93   171, 30
        67, 93   175, 30
        */
        xy1.add(46, 165);  xy2.add(167, 77);
        xy1.add(84, 162);  xy2.add(181, 81);
        xy1.add(63, 149);  xy2.add(173, 68);
        xy1.add(37, 122);  xy2.add(166, 49);
        xy1.add(88, 123);  xy2.add(183, 52);
        xy1.add(54, 93);  xy2.add(171, 30);
        xy1.add(67, 93);  xy2.add(175, 30);
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
    
    private void populateWithBrownAndLoweMatches7(PairIntArray xy1,
        PairIntArray xy2) {
        /*
        295, 180  15, 149
        384, 259  97, 243
        327, 365  13, 348
        465, 154  191, 154
        494, 259  197, 254
        481, 365  168, 352
        305, 243  14, 218
        */
        xy1.add(295, 180);  xy2.add(15, 149);
        xy1.add(384, 259);  xy2.add(97, 243);
        xy1.add(327, 365);  xy2.add(13, 348);
        xy1.add(465, 154);  xy2.add(191, 154);
        xy1.add(494, 259);  xy2.add(197, 254);
        xy1.add(481, 365);  xy2.add(168, 352);
        xy1.add(305, 243);  xy2.add(14, 218);
    }
    
    private void populateWithCampusMatches7(PairIntArray xy1,
        PairIntArray xy2) {
        /*
        10, 336   272, 339
        73, 206   334, 212
        33, 126   294, 134
        192, 103  452, 101
        172, 254  433, 256
        209, 337  472, 342
        84, 272   345, 275
        */
        
        xy1.add(10, 336);  xy2.add(272, 339);
        xy1.add(73, 206);  xy2.add(334, 212);
        xy1.add(33, 126);  xy2.add(294, 134);
        xy1.add(192, 103);  xy2.add(452, 101);
        xy1.add(172, 254);  xy2.add(433, 256);
        xy1.add(209, 337);  xy2.add(472, 342);
        xy1.add(84, 272);  xy2.add(345, 275);
    }
    
    private void populateWithCampusMatchesMoreThan7(PairIntArray xy1,
        PairIntArray xy2) {
        
        xy1.add(10, 336);  xy2.add(272, 339);
        xy1.add(73, 206);  xy2.add(334, 212);
        xy1.add(33, 126);  xy2.add(294, 134);
        xy1.add(192, 103);  xy2.add(452, 101);
        xy1.add(172, 254);  xy2.add(433, 256);
        xy1.add(209, 337);  xy2.add(472, 342);
        xy1.add(84, 272);  xy2.add(345, 275);
        xy1.add(131, 153);  xy2.add(389, 156);
        xy1.add(175, 350);  xy2.add(439, 354);
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
    
    private void overplot(TransformationParameters params, PairIntArray set1,
        PairIntArray set2, Image img1, Image img2, String lbl) 
        throws IOException {
                      
        ImageIOHelper.addCurveToImage(set1, img1, 3, 255, 0, 0);
        
        ImageIOHelper.addCurveToImage(set2, img2, 3, 255, 0, 0);
        
        MiscDebug.writeImage(img1, "_kp_1_" + lbl);
        MiscDebug.writeImage(img2, "_kp_2_" + lbl);

        CorrespondencePlotter plotter = new CorrespondencePlotter(img1, img2);
        for (int ii = 0; ii < set1.getN(); ++ii) {
            int x1 = set1.getX(ii);
            int y1 = set1.getY(ii);
            int x2 = set2.getX(ii);
            int y2 = set2.getY(ii);
            plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 1);
        }
        plotter.writeImage("_corres_euclidean_" + lbl);
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
