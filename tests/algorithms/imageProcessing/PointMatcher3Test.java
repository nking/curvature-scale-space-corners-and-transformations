package algorithms.imageProcessing;

import algorithms.compGeometry.PointPartitioner;
import static algorithms.imageProcessing.StereoProjectionTransformer.rewriteInto3ColumnMatrix;
import algorithms.misc.MiscMath;
import algorithms.util.ResourceFinder;
import algorithms.util.LinearRegression;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import java.awt.Color;
import java.io.IOException;
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
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.fail;

/**
 *
 * @author nichole
 */
public class PointMatcher3Test {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public PointMatcher3Test() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
   
    /*
    for more datasets:
    http://www.robots.ox.ac.uk/~vgg/data/data-mview.html
    */
    
    private double[][] createRandomPoints(SecureRandom sr, int nPoints,
        int xMin, int yMin, double xRange, double yRange) {
     
        double[][] sMatrix = new double[nPoints][2];
        for (int i = 0; i < nPoints; i++) {
            sMatrix[i] = new double[2];
            double x = xMin + (sr.nextDouble()*xRange);
            double y = yMin + (sr.nextDouble()*yRange);
            sMatrix[i][0] = x;
            sMatrix[i][1] = y;
        }
        
        return sMatrix;
    }
    
    private double[][] createCopyOfSize(double[][] m, int sizeToCreate) {
         
        int end = m.length;
        if (sizeToCreate < end) {
            end = sizeToCreate;
        }
        
        double[][] copy = new double[sizeToCreate][2];
        for (int i = 0; i < end; i++) {
            copy[i] = new double[2];
            copy[i][0] = m[i][0];
            copy[i][1] = m[i][1];
        }
        
        return copy;
    }
    
    private void populateWithRandomPoints( SecureRandom sr, double[][] m, 
        int nPoints1, int nPoints2, int xMin, int yMin,
        double xRange, double yRange) {
        
        for (int i = nPoints1; i < nPoints2; i++) {
            m[i] = new double[2];
            double x = xMin + (sr.nextDouble() * xRange);
            double y = yMin + (sr.nextDouble() * yRange);
            m[i][0] = x;
            m[i][1] = y;
        }
    }
    
    private void translateX(double[][] m, double translateX) {
        
        for (int i = 0; i < m.length; i++) {
            double x = m[i][0];
            m[i][0] = x + translateX;
        }
    }
    private void translateY(double[][] m, double translateY) {
        
        for (int i = 0; i < m.length; i++) {
            double y = m[i][1];
            m[i][1] = y + translateY;
        }
    }
   
    @Test
    public void test1() throws Exception {

        // test for dataset which already matches exactly
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nPoints1 = 70;
        
        double xRange = 400;
        double yRange = 300;
        
        int x1Min = 0;
        int y1Min = 0;
        int x2Min = 0;
        int y2Min = 0;
        
        double[][] sMatrix = createRandomPoints(sr, nPoints1,
            x1Min, y1Min, xRange, yRange);
                
        double[][] mMatrix = createCopyOfSize(sMatrix, nPoints1);
        
        PointMatcher pointMatcher = new PointMatcher();
        
        double[][] p = new double[3][];
        /*p[0] = new double[]{1,    0., transX};
        p[1] = new double[]{0.0,    1, transY};
        p[2] = new double[]{0,       0, noeffect};*/
        p[0] = new double[]{1,    0., 0};
        p[1] = new double[]{0.0,   1, 0};
        p[2] = new double[]{0,     0, 1};
        
        ProjectiveFit fit = 
            pointMatcher.calculateProjectiveTransformationUsingDownhillSimplex(
            sMatrix, mMatrix, p, sMatrix, mMatrix);
        
        double[][] params = fit.getProjection();
        
        assertNotNull(params);
        
        double[][] transformed = pointMatcher.transformUsingProjection(params, 
            sMatrix);
        overplotTransformed(transformed, mMatrix, (int)(x2Min + xRange),
           (int) (y2Min + yRange), 1);
        
        assertTrue(params[0][0] == 1);
        assertTrue(params[0][1] == 0);
        assertTrue(params[0][2] == 0);
        
        assertTrue(params[1][0] == 0);
        assertTrue(params[1][1] == 1);
        assertTrue(params[1][2] == 0);
        
        assertTrue(params[2][0] == 0);
        assertTrue(params[2][1] == 0);
        assertTrue(params[2][2] == 1);
        
    }
    
    @Test
    public void test2() throws Exception {

        // test for exact match plus noise
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nPoints1 = 70;
        
        double xRange = 400;
        double yRange = 300;
        
        int x1Min = 0;
        int y1Min = 0;
        int x2Min = 0;
        int y2Min = 0;
        
        double[][] sMatrix = createRandomPoints(sr, nPoints1,
            x1Min, y1Min, xRange, yRange);
        
        int nPoints2 = (int)(1.3f * nPoints1);
        
        double[][] mMatrix = createCopyOfSize(sMatrix, nPoints2);
        
        double translateX = 0;
        double translateY = 0;
        
        populateWithRandomPoints(sr, mMatrix, nPoints1, nPoints2, x2Min, y2Min,
            xRange, yRange);
        
        PointMatcher pointMatcher = new PointMatcher();
        
        double[][] p = new double[3][];
        /*p[0] = new double[]{1,    0., transX};
        p[1] = new double[]{0.0,    1, transY};
        p[2] = new double[]{0,       0, noeffect};*/
        p[0] = new double[]{1,    0., 0};
        p[1] = new double[]{0.0,   1, 0};
        p[2] = new double[]{0,     0, 1};
        
        ProjectiveFit fit = 
            pointMatcher.calculateProjectiveTransformationUsingDownhillSimplex(
            sMatrix, mMatrix, p, sMatrix, mMatrix);
        
        double[][] params = fit.getProjection();
        
        assertNotNull(params);
        
        double[][] transformed = pointMatcher.transformUsingProjection(params, 
            sMatrix);
        overplotTransformed(transformed, mMatrix, (int)(x2Min + xRange),
           (int) (y2Min + yRange), 1);
        
        assertTrue(params[0][0] == 1);
        assertTrue(params[0][1] == 0);
        double diff = Math.abs(params[0][2] - translateX);
        if (!(diff < 0.1)) {
            System.out.println("params[0][2]=" + params[0][2] 
                + " translateX=" + translateX);
        }
        assertTrue(diff < 0.1);
        
        assertTrue(params[1][0] == 0);
        assertTrue(params[1][1] == 1);
        diff = Math.abs(params[1][2] - translateY);
        if (!(diff < 0.1)) {
            System.out.println("params[1][2]=" + params[1][2] 
                + " translateY=" + translateY);
        }
        assertTrue(diff < 0.1);
        
        assertTrue(params[2][0] == 0);
        assertTrue(params[2][1] == 0);
        assertTrue(params[2][2] == 1);
        
    }
    
    @Test
    public void test3() throws Exception {

        // test for exact match translated in X
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1419991550083L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nPoints1 = 70;
        
        double xRange = 400;
        double yRange = 300;
        
        int x1Min = 0;
        int y1Min = 0;
        
        double[][] sMatrix = createRandomPoints(sr, nPoints1,
            x1Min, y1Min, xRange, yRange);
        
        int nPoints2 = nPoints1;
        
        double[][] mMatrix = createCopyOfSize(sMatrix, nPoints2);
        
        double translateX = 200;
        double translateY = 0;
        int x2Min = (int)translateX;
        int y2Min = (int)translateY;
        
        translateX(mMatrix, translateX);
        
        PointMatcher pointMatcher = new PointMatcher();
        
        double[][] p = new double[3][];
        /*p[0] = new double[]{1,    0., transX};
        p[1] = new double[]{0.0,    1, transY};
        p[2] = new double[]{0,       0, noeffect};*/
        p[0] = new double[]{1,    0., 0};
        p[1] = new double[]{0.0,   1, 0};
        p[2] = new double[]{0,     0, 1};
        
        ProjectiveFit fit = 
            pointMatcher.calculateProjectiveTransformationUsingDownhillSimplex(
            sMatrix, mMatrix, p, sMatrix, mMatrix);
        
        double[][] params = fit.getProjection();
        
        assertNotNull(params);
        
        int x0 = (int)(x2Min + x1Min + xRange);
        int y0 = (int)(y2Min + y1Min + yRange);
        
        double[][] transformed = pointMatcher.transformUsingProjection(params, 
            sMatrix);
        overplotTransformed(transformed, mMatrix, x0, y0, 1);
        
        assertTrue(params[0][0] == 1);
        assertTrue(params[0][1] == 0);
        assertTrue(Math.abs(params[0][2] - translateX) < 0.1);
        
        assertTrue(params[1][0] == 0);
        assertTrue(params[1][1] == 1);
        assertTrue(Math.abs(params[1][2] - translateY) < 0.1);
        
        assertTrue(params[2][0] == 0);
        assertTrue(params[2][1] == 0);
        assertTrue(params[2][2] == 1);
        
    }
    
    @Test
    public void test4() throws Exception {

        // test for exact match translated in X plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nPoints1 = 70;
        
        double xRange = 400;
        double yRange = 300;
        
        int x1Min = 0;
        int y1Min = 0;
        
        double[][] sMatrix = createRandomPoints(sr, nPoints1,
            x1Min, y1Min, xRange, yRange);
        
        int nPoints2 = (int)(1.3f * nPoints1);
        
        double[][] mMatrix = createCopyOfSize(sMatrix, nPoints2);
        
        double translateX = 200;
        double translateY = 0;
        int x2Min = (int)translateX;
        int y2Min = (int)translateY;
        
        translateX(mMatrix, translateX);
        
        populateWithRandomPoints(sr, mMatrix, nPoints1, nPoints2, x2Min, y2Min,
            xRange, yRange);
        
        PointMatcher pointMatcher = new PointMatcher();
        
        double[][] p = new double[3][];
        /*p[0] = new double[]{1,    0., transX};
        p[1] = new double[]{0.0,    1, transY};
        p[2] = new double[]{0,       0, noeffect};*/
        p[0] = new double[]{1,    0., 0};
        p[1] = new double[]{0.0,   1, 0};
        p[2] = new double[]{0,     0, 1};
        
        ProjectiveFit fit = 
            pointMatcher.calculateProjectiveTransformationUsingDownhillSimplex(
            sMatrix, mMatrix, p, sMatrix, mMatrix);
        
        double[][] params = fit.getProjection();
        
        assertNotNull(params);
        
        double[][] transformed = pointMatcher.transformUsingProjection(params, 
            sMatrix);
        overplotTransformed(transformed, mMatrix, (int)(x2Min + xRange),
           (int) (y2Min + yRange), 1);
        
        assertTrue(params[0][0] == 1);
        assertTrue(params[0][1] == 0);
        double diff = Math.abs(params[0][2] - translateX);
        if (!(diff < 0.1)) {
            System.out.println("params[0][2]=" + params[0][2] 
                + " translateX=" + translateX);
        }
        assertTrue(diff < 0.1);
        
        assertTrue(params[1][0] == 0);
        assertTrue(params[1][1] == 1);
        diff = Math.abs(params[1][2] - translateY);
        if (!(diff < 0.1)) {
            System.out.println("params[1][2]=" + params[1][2] 
                + " translateY=" + translateY);
        }
        assertTrue(diff < 0.1);
        
        assertTrue(params[2][0] == 0);
        assertTrue(params[2][1] == 0);
        assertTrue(params[2][2] == 1);
        
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
    
    private void overplotTransformed(double[][] transformed, double[][] model,
        int width, int height, int testNumber) throws IOException {
        
        Image image = new Image(width, height);
        
        for (int ii = 0; ii < transformed.length; ii++) {
            double x2 = transformed[ii][0];
            double y2 = transformed[ii][1];
            ImageIOHelper.addPointToImage((float) x2, (float) y2, image, 3, 
                0, 0, 255);
        }
        for (int ii = 0; ii < model.length; ii++) {
            double x = model[ii][0];
            double y = model[ii][1];
            ImageIOHelper.addPointToImage((float) x, (float) y, image, 2, 
                255, 0, 0);
        }
        String dirPath = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(
            dirPath + "/tmp_t2_" + testNumber + ".png", image);
        
    }
    
    public static void main(String[] args) {
        
        try {
            PointMatcher3Test test = new PointMatcher3Test();
            
            test.test1();
            test.test2();
            test.test3();
            test.test4();
            
            /*
            tests for :
            -- for same set w/ translation in y and noise
            -- for same set w/ translation in x and y and noise
            -- for same set w/ rotation and noise
            -- for same set w/ scale and noise
            -- for same set w/ rotation and translation and noise
            -- for same set w/ scale and rotation and translation and noise
            */
                        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }

}
