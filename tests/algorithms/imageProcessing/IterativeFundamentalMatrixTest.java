package algorithms.imageProcessing;

import Jama.Matrix;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.awt.Color;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class IterativeFundamentalMatrixTest {
    
    public IterativeFundamentalMatrixTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    //@Test
    public void testSolveUsingRANSAC() throws Exception {
        
        // mRows = 3;  nCols = 484
        PairFloatArray matched1 = DataForTests.readMerton1UnnormalizedXY1Data();
        
        PairFloatArray matched2 = DataForTests.readMerton1UnnormalizedXY2Data();
        
        IterativeFundamentalMatrix iterativeFM = 
            new IterativeFundamentalMatrix();
        
        StereoProjectionTransformer spTransformer = 
            iterativeFM.solveUsingRANSAC(matched1, matched2);
        
        String fileName1 = "merton_college_I_001.jpg";
        String fileName2 = "merton_college_I_002.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        
        int nLimitTo = 10;//nPoints;
        
        Color clr = null;
        PairIntArray subsetLeft = new PairIntArray();
        for (int i = 0; i < nLimitTo; i++) {
            
            clr = getColor(clr);
            
            PairIntArray leftLine = spTransformer.getEpipolarLineInLeft(
                img1.getWidth(), i);
            
            ImageIOHelper.addCurveToImage(leftLine, img1, 0, 
                clr.getRed(), clr.getGreen(), clr.getBlue());
            
            subsetLeft.add(Math.round(matched1.getX(i)), Math.round(matched1.getY(i)));
        }
        ImageIOHelper.addCurveToImage(subsetLeft, img1, 2, 255, 0, 0);
        
        clr = null;
        PairIntArray subsetRight = new PairIntArray();
        for (int i = 0; i < nLimitTo; i++) {
            clr = getColor(clr);
            
            PairIntArray rightLine = spTransformer.getEpipolarLineInRight(
                img2.getWidth(), i);
            
            ImageIOHelper.addCurveToImage(rightLine, img2, 0, 
                clr.getRed(), clr.getGreen(), clr.getBlue()); 
            
            subsetRight.add(Math.round(matched2.getX(i)), Math.round(matched2.getY(i)));
        }
        ImageIOHelper.addCurveToImage(subsetRight, img2, 2, 255, 0, 0);
    
        String dirPath = ResourceFinder.findDirectory("bin");
        
        ImageIOHelper.writeOutputImage(dirPath + "/tmp_ransac_1.png", img1);
       
        ImageIOHelper.writeOutputImage(dirPath + "/tmp_ransac_2.png", img2);
    }
    
    //@Test
    public void testSolveUsingOutlierRemoval() throws Exception {
        
        // mRows = 3;  nCols = 484
        PairFloatArray matched1 = DataForTests.readMerton1UnnormalizedXY1Data();
        
        PairFloatArray matched2 = DataForTests.readMerton1UnnormalizedXY2Data();
                
        String fileName1 = "merton_college_I_001.jpg";
        String fileName2 = "merton_college_I_002.jpg";
        
        StereoProjectionTransformer spTransformer = 
            solveUsingOutlierRemoval(matched1, matched2, fileName1, fileName2);
    
        //already contains exact matching points, so should recover all:
        
        assertTrue(matched1.getN() == spTransformer.getNumberOfMatches());
    }
    
    @Test
    public void testSolveUsingOutlierRemoval2() throws Exception {
        
        PairFloatArray matched1 = new PairFloatArray();
        PairFloatArray matched2 = new PairFloatArray();
        
        DataForTests.readBrownAndLoweMatches(matched1, matched2);
        
        String fileName1 = "brown_lowe_2003_image1.jpg";
        String fileName2 = "brown_lowe_2003_image2.jpg";
        
        //already contains exact matching points, so should recover all:
        
        StereoProjectionTransformer spTransformer = 
            solveUsingOutlierRemoval(matched1, matched2, fileName1, fileName2);
        
        assertTrue(matched1.getN() == spTransformer.getNumberOfMatches());
        
        /*
        TODO: repeat same test with corners generated w/ outdoor method
        and paired with
        
        TransformationPointFit calculateTransformation(PairIntArray set1, 
            PairIntArray set2, int image1Width, int image1Height);
        
        
        */
    
    }
    
    private StereoProjectionTransformer solveUsingOutlierRemoval(
        PairFloatArray matched1, PairFloatArray matched2, String imageFileName1, 
        String imageFileName2) throws Exception {
        
        IterativeFundamentalMatrix iterativeFM = 
            new IterativeFundamentalMatrix();
        
        StereoProjectionTransformer spTransformer = 
            iterativeFM.solveUsingOutlierRemoval(matched1, matched2);
        
        String fileName1 = imageFileName1;
        String fileName2 = imageFileName2;
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        
        int nLimitTo = 10;//spTransformer.getNumberOfMatches();
        
        matched1 = spTransformer.getLeftXYFloat();
        
        matched2 = spTransformer.getRightXYFloat();
        
        Color clr = null;
        PairIntArray subsetLeft = new PairIntArray();
        for (int i = 0; i < nLimitTo; i++) {
            
            clr = getColor(clr);
            
            PairIntArray leftLine = spTransformer.getEpipolarLineInLeft(
                img1.getWidth(), i);
            
            ImageIOHelper.addCurveToImage(leftLine, img1, 0, 
                clr.getRed(), clr.getGreen(), clr.getBlue());
            
            subsetLeft.add(Math.round(matched1.getX(i)), Math.round(matched1.getY(i)));
        }
        ImageIOHelper.addCurveToImage(subsetLeft, img1, 2, 255, 0, 0);
        
        clr = null;
        PairIntArray subsetRight = new PairIntArray();
        for (int i = 0; i < nLimitTo; i++) {
            clr = getColor(clr);
            
            PairIntArray rightLine = spTransformer.getEpipolarLineInRight(
                img2.getWidth(), i);
            
            ImageIOHelper.addCurveToImage(rightLine, img2, 0, 
                clr.getRed(), clr.getGreen(), clr.getBlue()); 
            
            subsetRight.add(Math.round(matched2.getX(i)), Math.round(matched2.getY(i)));
        }
        ImageIOHelper.addCurveToImage(subsetRight, img2, 2, 255, 0, 0);
    
        String dirPath = ResourceFinder.findDirectory("bin");
        
        ImageIOHelper.writeOutputImage(dirPath + "/tmp_outlierrm_1.png", img1);
       
        ImageIOHelper.writeOutputImage(dirPath + "/tmp_outlierrm_2.png", img2);
        
        return spTransformer;
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
