package algorithms.imageProcessing;

import static algorithms.imageProcessing.StereoProjectionTransformer.rewriteInto3ColumnMatrix;
import algorithms.misc.MiscMath;
import algorithms.util.ResourceFinder;
import algorithms.util.LinearRegression;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import java.awt.Color;
import java.util.List;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.ejml.simple.*;

/**
 *
 * @author nichole
 */
public class PointMatcher2Test {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public PointMatcher2Test() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    /*
     * 
     * this is a scratch pad for trying a point matcher that searches the
    fundamental matrix space.
     * 
     */
    
    /*
    for more datasets:
    http://www.robots.ox.ac.uk/~vgg/data/data-mview.html
    */
    
    @Test
    public void test0() throws Exception {
                
        /*
        get FM of the merton college set.
           calc max and min of normalized values that create the FM
        
        do the same for the brown & lowe 2003 dataset which has
        a large translation.
        
        */
        
        PairIntArray matched1 = new PairIntArray();
        PairIntArray matched2 = new PairIntArray();
        //DataForTests.readMerton1Matched(matched1, matched2);
        //DataForTests.readBrownAndLoweMatches(matched1, matched2);
        DataForTests.readBrownAndLoweInflectionPointsImage1(matched1);
        DataForTests.readBrownAndLoweInflectionPointsImage2(matched2);
        
        String fileName1 = "merton_college_I_001.jpg";
        fileName1 = "brown_lowe_2003_image1.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();
        String fileName2 = "merton_college_I_002.jpg";
        fileName2 = "brown_lowe_2003_image2.jpg";
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        
        StereoProjectionTransformer spTransformer = 
            new StereoProjectionTransformer();
        
        PairIntArray outputLeft = new PairIntArray();
        PairIntArray outputRight = new PairIntArray();
        
        SimpleMatrix fm = 
            spTransformer.calculateEpipolarProjectionForUnmatched(
            matched1, matched2, image1Width, image1Height,
            outputLeft, outputRight);
        
        PairFloatArray finalMatched1 = new PairFloatArray();
        PairFloatArray finalMatched2 = new PairFloatArray();
        for (int i = 0; i < outputLeft.getN(); i++) {
            finalMatched1.add(outputLeft.getX(i), outputLeft.getY(i));
            finalMatched2.add(outputRight.getX(i), outputRight.getY(i));
        }
        
        img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        image1Width = img1.getWidth();
        image1Height = img1.getHeight();
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        
        SimpleMatrix input1 = 
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(finalMatched1);
        SimpleMatrix input2 = 
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(finalMatched2);

        for (int ii = 0; ii < input1.numCols(); ii++) {
            double x1 = input1.get(0, ii);
            double y1 = input1.get(1, ii);
            ImageIOHelper.addPointToImage((float) x1, (float) y1, img1, 3, 
                255, 0, 0);
            double x2 = input2.get(0, ii);
            double y2 = input2.get(1, ii);
            ImageIOHelper.addPointToImage((float) x2, (float) y2, img2, 3, 
                255, 0, 0);
        }
 
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
        ImageIOHelper.writeOutputImage(dirPath + "/tmp_m_1.png", img1);
        ImageIOHelper.writeOutputImage(dirPath + "/tmp_m_2.png", img2);

        int z = 1;
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
            PointMatcher2Test test = new PointMatcher2Test();
            
            test.test0();
            
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
    }

}
