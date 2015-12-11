package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import java.awt.Color;
import java.io.IOException;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;
import org.ejml.simple.*;

/**
 *
 * @author nichole
 */
public class StereoProjectionTransformerTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public StereoProjectionTransformerTest() {
    }
    
    public void testCalcEpipolar() throws Exception {

        String fileName1 = "checkerboard_01.jpg";
        String fileName2 = "checkerboard_02.jpg";
        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);
        
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
        
        
        FeatureMatcherSettings settings = new FeatureMatcherSettings();
        settings.setDebug(true);
        settings.setStartWithBinnedImages(true);
        settings.setUseNormalizedFeatures(false);
        
        EpipolarSolver solver = new EpipolarSolver(img1, img2, settings);
        StereoProjectionTransformerFit fit = solver.solve();
        
        int[] x1 = new int[]{248,341,154,341,341,339,341,249,341,155,};
        int[] y1 = new int[]{109,204,295,204,204,295,204,205,204,201,};
        int[] x2 = new int[]{264,169,354,169,169,170,169,264,169,353,};
        int[] y2 = new int[]{385,293,202,293,293,198,293,290,293,292,};

        StereoProjectionTransformer spTr = new StereoProjectionTransformer();

        /*
        public SimpleMatrix calculateEpipolarProjectionForPerfectlyMatched(
        PairFloatArray pointsLeftXY,  PairFloatArray pointsRightXY) {
        */

        //PairFloatArray calculateDistancesFromEpipolar(
        //SimpleMatrix fm, SimpleMatrix matchedLeftPoints, 
        //SimpleMatrix matchedRightPoints) {
        
    }
    
    /*
    for more datasets:
    http://www.robots.ox.ac.uk/~vgg/data/data-mview.html
    */
    
    public void estRANSAC() throws Exception {
    
        //TODO: add the reference for this data here.
        String fileName1 = "merton_college_I_001.jpg";
        String fileName2 = "merton_college_I_002.jpg";
            
        PairFloatArray left7 = new PairFloatArray();
        PairFloatArray right7 = new PairFloatArray();
        getMertonCollege7TrueMatches(left7, right7);
        
        StereoProjectionTransformer spTr = new StereoProjectionTransformer();
        
        SimpleMatrix left7XY = 
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(left7);
        SimpleMatrix right7XY =
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(right7);
        
        List<SimpleMatrix> fits = spTr.calculateEpipolarProjectionFor7Points(
            left7XY, right7XY);
         
        assertTrue(fits.size() == 1);
        
        SimpleMatrix fit = fits.get(0);
                        
        log.info("fit=" + fit.toString());

        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);

        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        int image2Width = img2.getWidth();
        int image2Height = img2.getHeight();

        overplotEpipolarLines(fit, left7, right7,
            img1, img2, 
            image1Width, image1Height, image2Width, image2Height, 
            Integer.valueOf(0).toString()); 
        
        /*
        solution for 7pt epipolar when points are not normalized:
        INFO: fit=Type = dense , numRows = 3 , numCols = 3
         [junit] -0.000  -0.000  -0.006  
         [junit]  0.000  -0.000   0.062  
         [junit]  0.004  -0.071   1.000 
        */
        
        assertTrue(Math.abs(fit.get(0, 0) - 0) < 0.005);
        assertTrue(Math.abs(fit.get(1, 0) - 0) < 0.005);
        assertTrue(Math.abs(fit.get(2, 0) - 0) < 0.005);
        
        assertTrue(Math.abs(fit.get(0, 1) - 0) < 0.005);
        assertTrue(Math.abs(fit.get(1, 1) - 0) < 0.005);
        assertTrue(Math.abs(fit.get(2, 1) - -0.071) < 0.005);
        
        assertTrue(Math.abs(fit.get(0, 2) - -0.006) < 0.005);
        assertTrue(Math.abs(fit.get(1, 2) - 0.062) < 0.005);
        assertTrue(Math.abs(fit.get(2, 2) - 1.0) < 0.005);
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
    
    private void overplotEpipolarLines(SimpleMatrix fm, PairFloatArray set1,
        PairFloatArray set2, Image img1, Image img2, int image1Width,
        int image1Height, int image2Width, int image2Height, String outfileNumber) 
        throws IOException {

        SimpleMatrix input1 =
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(set1);

        SimpleMatrix input2 =
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(set2);

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

        StereoProjectionTransformer spTransformer = new
            StereoProjectionTransformer();

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
    public static void main(String[] args) {
        
        try {
            StereoProjectionTransformerTest test = 
                new StereoProjectionTransformerTest();
            
            //test.testRANSAC();
                        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
    }

}
