package algorithms.imageProcessing;

import algorithms.compGeometry.PointPartitioner;
import static algorithms.imageProcessing.StereoProjectionTransformer.rewriteInto3ColumnMatrix;
import algorithms.imageProcessing.util.MatrixUtil;
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
  
    @Test
    public void test1() throws Exception {

        // test for dataset which already matches exactly
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = nScenePoints;
        
        double scale = 1.;
        double rotation = 0.;
        int translateX = 0;
        int translateY = 0;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 1);
    }
    
    @Test
    public void test2() throws Exception {

        // test for exact match plus noise
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);

        double scale = 1.;
        double rotation = 0.;
        int translateX = 0;
        int translateY = 0;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 2);
    }
    
    @Test
    public void test3() throws Exception {

        // test for exact match translated in X
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1419991550083L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = nScenePoints;
                
        double scale = 1.;
        double rotation = 0.;
        int translateX = 200;
        int translateY = 0;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 3);
        
    }
    
    @Test
    public void test4() throws Exception {

        // test for exact match translated in X plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
       
        int nModelPoints = (int)(1.3f * nScenePoints);
                
        double scale = 1.;
        double rotation = 0.;
        int translateX = 200;
        int translateY = 0;
                
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 4);
    }
    
    @Test
    public void test5() throws Exception {

        // test for exact match translated in X
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1419991550083L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = nScenePoints;
                
        double scale = 1.;
        double rotation = 0.;
        int translateX = 200;
        int translateY = 110;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 5);
        
    }
    
    @Test
    public void test6() throws Exception {

        // test for exact match translated in X plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
                
        double scale = 1.;
        double rotation = 0.;
        int translateX = 200;
        int translateY = 110;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 6);
    }
    
    @Test
    public void test7() throws Exception {

        // test for exact match rotated by 30 degrees
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1420149550374L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = nScenePoints;
                
        double rotation = 30. * Math.PI/180.;
        double scale = 1.;
        int translateX = 0;
        int translateY = 0;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 7);
        
    }
  
    @Test
    public void test8() throws Exception {

        // test for rotation plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
                
        double rotation = 30. * Math.PI/180.;
        double scale = 1;
        int translateX = 0;
        int translateY = 0;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 8);
    }
    
    @Test
    public void test9() throws Exception {

        // test for scale plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
        
        double rotation = 0.;
        double scale = 4.;
        int translateX = 0;
        int translateY = 0;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 9);
    }
    
    @Test
    public void test10() throws Exception {

        // test for scale smaller than 1, plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1420179838620L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
        
        double rotation = 0.;
        double scale = 1./4.;
        int translateX = 0;
        int translateY = 0;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 10);
    }
    
    @Test
    public void test11() throws Exception {

        // test for rotation and translation, plus random points
       
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1420159107635L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
        
        double rotation = 30.*Math.PI/180.;
        double scale = 1.;
        int translateX = 200;
        int translateY = -10;
     
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 11);
    }

    @Test
    public void test12() throws Exception {

        // test for scale, rotation and translation, plus random points
       
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1420159107635L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
        
        double rotation = 30.*Math.PI/180.;
        double scale = 2.;
        int translateX = 200;
        int translateY = -10;
     
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 12);
    }
    
    @Test
    public void test13() throws Exception {

        // test for scale, rotation and translation, plus random points
       
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1420159107635L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
        
        double rotation = 14.*Math.PI/180.;
        double scale = 1.;
        int translateX = 280;
        int translateY = -14;

        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 13);
    }
    
    @Test
    public void test14() throws Exception {

        // test for scale close to 1
       
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1420187783647L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
        
        double rotation = 0.;
        double scale = 1.2;
        int translateX = 100;
        int translateY = -14;

        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 13);
    }
    
    private void runTest(SecureRandom sr, int nScenePoints, int nModelPoints, 
        int xRange, int yRange,
        double scale, double rotation, int translationX, int translationY,
        int testNumber) throws Exception {       
        
        int xModelMin = 0;
        int yModelMin = 0;
        
        PairIntArray model = createRandomPoints(sr, nScenePoints,
            xModelMin, yModelMin, xRange, yRange);
        
        PairIntArray scene = model.copy();
                        
        populateWithRandomPoints(sr, model, nScenePoints, nModelPoints, 
            xModelMin, yModelMin, xRange, yRange);
        
        int xModelCentroid = (xModelMin + xRange) >> 1;
        int yModelCentroid = (yModelMin + yRange) >> 1;
        
        scaleAndRotate(scene, 1./scale, -1*rotation, xModelCentroid, 
            yModelCentroid);
        int tx = (translationX == 0) ? 0 : (int)(-1*translationX/scale);
        int ty = (translationY == 0) ? 0 : (int)(-1*translationY/scale);
        translateX(scene, tx);
        translateY(scene, ty);
        
        // transform the centroid point from model to use for calculations
        TransformationParameters revParams = new TransformationParameters();
        revParams.setScale((float)(1./scale));
        revParams.setRotationInRadians(-1.f*(float)rotation);
        revParams.setTranslationX(tx);
        revParams.setTranslationY(ty);
        
        MatchedPointsTransformationCalculator tc = new 
            MatchedPointsTransformationCalculator();
        double[] xySceneCen = tc.applyTransformation(revParams, 
            xModelCentroid, yModelCentroid, xModelCentroid, yModelCentroid);
           
        int xSceneCentroid = (int)xySceneCen[0];
        int ySceneCentroid = (int)xySceneCen[1];
        
        PointMatcher pointMatcher = new PointMatcher();
        
        TransformationPointFit fit = 
            pointMatcher.calculateProjectiveTransformationWrapper(
            scene, model, xSceneCentroid, ySceneCentroid, 
            xModelCentroid, yModelCentroid);

        System.out.println("=> " + fit.toString());
        
        TransformationParameters params = fit.getParameters();

        Transformer transformer = new Transformer();

        PairIntArray transformed = transformer.applyTransformation(
            params, scene, xSceneCentroid, ySceneCentroid);

        overplotTransformed(transformed, model, xRange, yRange, testNumber);

        int count = 0;
        for (int i = 0; i < transformed.getN(); i++) {
            int x = transformed.getX(i);
            // tolerance?
            if ((x < 0) || (x > 2*xModelCentroid)) {
                continue;
            }
            int y = transformed.getY(i);
            if ((y < 0) || (y > 2*yModelCentroid)) {
                continue;
            }
            count++;
        }
        
        System.out.println("=> Number of transformed scene points within bounds = " 
            + count);
        
        int nExpected = (nScenePoints > count) ? count : nScenePoints; 
        
        assertTrue(Math.abs(nExpected - fit.getNumberOfMatchedPoints()) 
            < 0.1*nScenePoints);
        
        assertTrue(Math.abs(params.getRotationInRadians() - rotation) <= 10.0);
        assertTrue(Math.abs(params.getScale() - scale) < 1.0);
        assertTrue(Math.abs(params.getTranslationX() - translationX) <= 1.0);
        assertTrue(Math.abs(params.getTranslationY() - translationY) <= 1.0);
        
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
    
    private void populateWithRandomPoints(SecureRandom sr, PairIntArray m, 
        int nPoints1, int nPoints2, int xMin, int yMin,
        int xRange, int yRange) {
        
        for (int i = nPoints1; i < nPoints2; i++) {
            int x = xMin + sr.nextInt(xRange);
            int y = yMin + sr.nextInt(yRange);
            m.add(x, y);
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
    
    private void translateX(PairIntArray m, int translateX) {
        for (int i = 0; i < m.getN(); i++) {
            int x = m.getX(i);
            int y = m.getY(i);
            m.set(i, (x + translateX), y);
        }
    }
    
    private void translateY(PairIntArray m, int translateY) {
        for (int i = 0; i < m.getN(); i++) {
            int x = m.getX(i);
            int y = m.getY(i);
            m.set(i, x, (y + translateY));
        }
    }
    
    private void scaleAndRotate(double[][] m, 
        double scale, double rotationInRadians, 
        double centroidX, double centroidY) {
        /*
        xr[i] = centroidX1*s + ( 
                ((x - centroidX1) * scaleTimesCosine) +
                ((y - centroidY1) * scaleTimesSine));
        yr[i] = centroidY1*s + ( 
                (-(x - centroidX1) * scaleTimesSine) +
                ((y - centroidY1) * scaleTimesCosine));
        */
        double scaleTimesCosine = scale * Math.cos(rotationInRadians);
        double scaleTimesSine = scale * Math.sin(rotationInRadians);
        
        for (int i = 0; i < m.length; i++) {
            double x = m[i][0];
            double y = m[i][1];
            double rx = centroidX*scale + ( 
                ((x - centroidX) * scaleTimesCosine) +
                ((y - centroidY) * scaleTimesSine));
            double ry = centroidY*scale + ( 
                (-(x - centroidX) * scaleTimesSine) +
                ((y - centroidY) * scaleTimesCosine));
            m[i][0] = rx;
            m[i][1] = ry;
        }
    }
    
    private void scaleAndRotate(PairIntArray m, 
        double scale, double rotationInRadians, 
        double centroidX, double centroidY) {
        
        //rotationInRadians *= -1;
        
        /*
        xr[i] = centroidX1*s + ( 
                ((x - centroidX1) * scaleTimesCosine) +
                ((y - centroidY1) * scaleTimesSine));
        yr[i] = centroidY1*s + ( 
                (-(x - centroidX1) * scaleTimesSine) +
                ((y - centroidY1) * scaleTimesCosine));
        */
        double scaleTimesCosine = scale * Math.cos(rotationInRadians);
        double scaleTimesSine = scale * Math.sin(rotationInRadians);
        
        for (int i = 0; i < m.getN(); i++) {
            double x = m.getX(i);
            double y = m.getY(i);
            double rx = centroidX*scale + ( 
                ((x - centroidX) * scaleTimesCosine) +
                ((y - centroidY) * scaleTimesSine));
            double ry = centroidY*scale + ( 
                (-(x - centroidX) * scaleTimesSine) +
                ((y - centroidY) * scaleTimesCosine));
            
            m.set(i, (int)Math.round(rx), (int)Math.round(ry));
        }
    }
      
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
    
    private PairIntArray createRandomPoints(SecureRandom sr, int nPoints,
        int xMin, int yMin, int xRange, int yRange) {
     
        PairIntArray output = new PairIntArray(nPoints);
        for (int i = 0; i < nPoints; i++) {
            int x = xMin + sr.nextInt(xRange);
            int y = yMin + sr.nextInt(yRange);
            output.add(x, y);
        }
        
        return output;
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
    
    public static void main(String[] args) {
        
        try {
            PointMatcher3Test test = new PointMatcher3Test();
            
            test.test1();
            test.test2();
            test.test3();
            test.test4();
            test.test5();
            test.test6();
            test.test7();      
            test.test8();
            test.test9();
            test.test10();          
            test.test11();
            test.test12();
            test.test13();
            test.test14();
            
            /*
            tests for :
            -- for same set w/ projection
            -- for same set w/ projection and noise
                        
            tests for scales which are close to 1 and less than 2
            */
        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }

    private void overplotTransformed(PairIntArray transformed, 
        PairIntArray model, int width, int height, int testNumber) 
        throws IOException {
        
        Image image = new Image(width, height);
        
        for (int ii = 0; ii < transformed.getN(); ii++) {
            double x2 = transformed.getX(ii);
            double y2 = transformed.getY(ii);
            ImageIOHelper.addPointToImage((float) x2, (float) y2, image, 3, 
                0, 0, 255);
        }
        for (int ii = 0; ii < model.getN(); ii++) {
            double x = model.getX(ii);
            double y = model.getY(ii);
            ImageIOHelper.addPointToImage((float) x, (float) y, image, 2, 
                255, 0, 0);
        }
        String dirPath = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(
            dirPath + "/tmp_t2_" + testNumber + ".png", image);
    }

}
