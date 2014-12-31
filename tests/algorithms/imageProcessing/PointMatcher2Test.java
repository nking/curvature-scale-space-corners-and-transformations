package algorithms.imageProcessing;

import algorithms.compGeometry.PointPartitioner;
import static algorithms.imageProcessing.StereoProjectionTransformer.rewriteInto3ColumnMatrix;
import algorithms.misc.MiscMath;
import algorithms.util.ResourceFinder;
import algorithms.util.LinearRegression;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import java.awt.Color;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.ejml.simple.*;
import static org.junit.Assert.fail;

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
   
    @Test
    public void test1() throws Exception {
                
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
        
        double[][] allS = new double[matched1.getN()][2];
        for (int i = 0; i < matched1.getN(); i++) {
            allS[i] = new double[2];
            allS[i][0] = matched1.getX(i);
            allS[i][1] = matched1.getY(i);
        }
        double[][] allM = new double[matched2.getN()][2];
        for (int i = 0; i < matched2.getN(); i++) {
            allM[i] = new double[2];
            allM[i][0] = matched2.getX(i);
            allM[i][1] = matched2.getY(i);
        }
        
        int nD = 2;
        PointPartitioner partitioner = new PointPartitioner();
        PairIntArray[] partitions1 = partitioner.partition(matched1, nD);
        PairIntArray[] partitions2 = partitioner.partition(matched2, nD);
        for (int i1 = 0; i1 < partitions1.length; i1++) {
            PairIntArray set1 = partitions1[i1];
            int nPoints1 = set1.getN();
            double[][] sMatrix = new double[nPoints1][2];
            for (int i = 0; i < nPoints1; i++) {
                sMatrix[i] = new double[2];
                sMatrix[i][0] = set1.getX(i);
                sMatrix[i][1] = set1.getY(i);
            }
            
            for (int i2 = 0; i2 < partitions2.length; i2++) {
                PairIntArray set2 = partitions2[i2];
                int nPoints2 = set2.getN();
                double[][] mMatrix = new double[nPoints2][2];
                double[][] m = new double[nPoints2][2];
                for (int i = 0; i < nPoints2; i++) {
                    mMatrix[i] = new double[2];
                    mMatrix[i][0] = set2.getX(i);
                    mMatrix[i][1] = set2.getY(i);
                    m[i] = new double[2];
                    m[i][0] = set2.getX(i);
                    m[i][1] = set2.getY(i);
                }
                
                PointMatcher pointMatcher = new PointMatcher();
                float[] transXY = pointMatcher.calculateTranslation(sMatrix, 
                    mMatrix);
                
                int z = 1;
                
                if (Float.isNaN(transXY[0]) || Float.isNaN(transXY[1])) {
                    continue;
                }
                
                run(sMatrix, m, filePath1, filePath2, allS, allM);
            }
        }        
    }
    
    private double[] minSValueOfColumns = null;
    private double[] maxSValueOfColumns = null;
            
    private double[][] register(double[][] theScene, double[][] theModel,
        double[][] allScene, double[][] allModel) {
        
        double[][] p = new double[3][];
        /*p[0] = new double[]{1,    0., transX};
        p[1] = new double[]{0.0,    1, transY};
        p[2] = new double[]{0,       0, noeffect};*/
        p[0] = new double[]{1,    0., 0};
        p[1] = new double[]{0.0,   1, 0};
        p[2] = new double[]{0,     0, 1};
        
        PointMatcher pointMatcher = new PointMatcher();
        
        ProjectiveFit fit = 
            pointMatcher.calculateProjectiveTransformationUsingDownhillSimplex(
            theScene, theModel, p, allScene, allModel);
        
        double[][] params = fit.getProjection();
        
        return params;
    }
   
    private void run(double[][] sMatrix, double[][] mMatrix, String 
        filePath1, String filePath2, double[][] allS, double[][] allM) 
        throws Exception {
        
        double[][] proj = register(sMatrix, mMatrix,
            allS, allM);
        
        PointMatcher pointMatcher = new PointMatcher();
        
        double[][] transformed = pointMatcher.transformUsingProjection(
            proj, sMatrix);
                
        //Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        for (int ii = 0; ii < transformed.length; ii++) {
            double xt = transformed[ii][0];
            double yt = transformed[ii][1];
            ImageIOHelper.addPointToImage((float) xt, (float) yt, img2, 3, 
                0, 0, 255);
        }
        for (int ii = 0; ii < mMatrix.length; ii++) {
            double x2 = mMatrix[ii][0];
            double y2 = mMatrix[ii][1];
            ImageIOHelper.addPointToImage((float) x2, (float) y2, img2, 3, 
                255, 0, 0);
        }
        String dirPath = ResourceFinder.findDirectory("bin");
        //ImageIOHelper.writeOutputImage(dirPath + "/tmp_t_1.png", img1);
        ImageIOHelper.writeOutputImage(dirPath + "/tmp_t_2.png", img2);
        
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
            
            //test.test0();
            
            test.test1();
            
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
    }

}
