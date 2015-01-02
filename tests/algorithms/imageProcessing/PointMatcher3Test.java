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
        
        int nPoints1 = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int x1Min = 0;
        int y1Min = 0;
        int x2Min = 0;
        int y2Min = 0;
        
        PairIntArray scene = createRandomPoints(sr, nPoints1, x1Min, y1Min, 
            xRange, yRange);
        
        PairIntArray model = scene.copy();
        
        PointMatcher pointMatcher = new PointMatcher();
        
        TransformationPointFit fit = 
            pointMatcher.calculateProjectiveTransformation(
            scene, model, xRange, yRange);
        
        Transformer transformer = new Transformer();
        PairIntArray transformed = transformer.applyTransformation(
            fit.getParameters(), scene, xRange, yRange);
        
        overplotTransformed(transformed, model, (x2Min + xRange),
           (y2Min + yRange), 1);
        
        assertTrue(Math.abs(fit.getScale() - 1) < 0.1);
        assertTrue(Math.abs(fit.getRotationInRadians() - 0) < 0.1);
        assertTrue(Math.abs(fit.getTranslationX() - 0) < 1.0);
        assertTrue(Math.abs(fit.getTranslationY() - 0) < 1.0);
        
    }
    
    @Test
    public void test2() throws Exception {

        // test for exact match plus noise
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nPoints1 = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int x1Min = 0;
        int y1Min = 0;
        int x2Min = 0;
        int y2Min = 0;
        
        PairIntArray scene = createRandomPoints(sr, nPoints1,
            x1Min, y1Min, xRange, yRange);
        
        int nPoints2 = (int)(1.3f * nPoints1);
        
        PairIntArray model = scene.copy();
        
        double translateX = 0;
        double translateY = 0;
        
        populateWithRandomPoints(sr, model, nPoints1, nPoints2, x2Min, y2Min,
            xRange, yRange);
        
        PointMatcher pointMatcher = new PointMatcher();
        
        TransformationPointFit fit = 
            pointMatcher.calculateProjectiveTransformation(
            scene, model, xRange, yRange);
        
        Transformer transformer = new Transformer();
        PairIntArray transformed = transformer.applyTransformation(
            fit.getParameters(), scene, xRange, yRange);
        
        overplotTransformed(transformed, model, (x2Min + xRange),
           (y2Min + yRange), 2);
        
        assertTrue(Math.abs(fit.getScale() - 1) < 0.1);
        assertTrue(Math.abs(fit.getRotationInRadians() - 0) < 0.1);
        assertTrue(Math.abs(fit.getTranslationX() - 0) < 1.0);
        assertTrue(Math.abs(fit.getTranslationY() - 0) < 1.0);
        
    }
    
    @Test
    public void test3() throws Exception {

        // test for exact match translated in X
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1419991550083L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nPoints1 = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int x1Min = 0;
        int y1Min = 0;
        
        PairIntArray scene = createRandomPoints(sr, nPoints1,
            x1Min, y1Min, xRange, yRange);
        
        int nPoints2 = nPoints1;
        
        PairIntArray model = scene.copy();
        
        int translateX = 200;
        int translateY = 0;
        int x2Min = translateX;
        int y2Min = translateY;
        
        translateX(model, translateX);

        PointMatcher pointMatcher = new PointMatcher();
        
        TransformationPointFit fit = 
            pointMatcher.calculateProjectiveTransformation(
            scene, model, xRange, yRange);
        
        Transformer transformer = new Transformer();
        PairIntArray transformed = transformer.applyTransformation(
            fit.getParameters(), scene, xRange, yRange);
        
        overplotTransformed(transformed, model, (x2Min + xRange),
           (y2Min + yRange), 3);
        
        assertTrue(Math.abs(fit.getScale() - 1) < 0.1);
        assertTrue(Math.abs(fit.getRotationInRadians() - 0) < 0.1);
        assertTrue(Math.abs(fit.getTranslationX() - translateX) < 1.0);
        assertTrue(Math.abs(fit.getTranslationY() - translateY) < 1.0);
        
    }
    
    @Test
    public void test4() throws Exception {

        // test for exact match translated in X plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nPoints1 = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int x1Min = 0;
        int y1Min = 0;
        
        PairIntArray scene = createRandomPoints(sr, nPoints1,
            x1Min, y1Min, xRange, yRange);
        
        int nPoints2 = (int)(1.3f * nPoints1);
        
        PairIntArray model = scene.copy();
        
        int translateX = 200;
        int translateY = 0;
        int x2Min = translateX;
        int y2Min = translateY;
        
        translateX(model, translateX);
        
        populateWithRandomPoints(sr, model, nPoints1, nPoints2, x2Min, y2Min,
            xRange, yRange);
        
        PointMatcher pointMatcher = new PointMatcher();
                
        TransformationPointFit fit = 
            pointMatcher.calculateProjectiveTransformation(
            scene, model, xRange, yRange);
        
        System.out.println("=> " + fit.toString());
        
        Transformer transformer = new Transformer();
        PairIntArray transformed = transformer.applyTransformation(
            fit.getParameters(), scene, xRange, yRange);
        
        overplotTransformed(transformed, model, (x2Min + xRange),
           (y2Min + yRange), 4);
        
        assertTrue(Math.abs(fit.getScale() - 1) < 0.1);
        assertTrue(Math.abs(fit.getRotationInRadians() - 0) <= 1.0);
        assertTrue(Math.abs(fit.getTranslationX() - translateX) <= 1.0);
        assertTrue(Math.abs(fit.getTranslationY() - translateY) <= 1.0);
        
    }
    
    @Test
    public void test5() throws Exception {

        // test for exact match translated in X
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1419991550083L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nPoints1 = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int x1Min = 0;
        int y1Min = 0;
        
        PairIntArray scene = createRandomPoints(sr, nPoints1,
            x1Min, y1Min, xRange, yRange);
        
        int nPoints2 = nPoints1;
        
        PairIntArray model = scene.copy();
        
        int translateX = 200;
        int translateY = 110;
        int x2Min = translateX;
        int y2Min = translateY;
        
        translateX(model, translateX);
        translateY(model, translateY);
        
        PointMatcher pointMatcher = new PointMatcher();
        
        TransformationPointFit fit = 
            pointMatcher.calculateProjectiveTransformation(
            scene, model, xRange, yRange);
        
        System.out.println("=> " + fit.toString());
        
        Transformer transformer = new Transformer();
        PairIntArray transformed = transformer.applyTransformation(
            fit.getParameters(), scene, xRange, yRange);
        
        overplotTransformed(transformed, model, (x2Min + xRange),
           (y2Min + yRange), 5);
        
        assertTrue(Math.abs(fit.getScale() - 1) < 0.1);
        assertTrue(Math.abs(fit.getRotationInRadians() - 0) <= 1.0);
        assertTrue(Math.abs(fit.getTranslationX() - translateX) <= 1.0);
        assertTrue(Math.abs(fit.getTranslationY() - translateY) <= 1.0);
        
    }
    
    @Test
    public void test6() throws Exception {

        // test for exact match translated in X plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nPoints1 = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int x1Min = 0;
        int y1Min = 0;
        
        PairIntArray scene = createRandomPoints(sr, nPoints1,
            x1Min, y1Min, xRange, yRange);
        
        int nPoints2 = (int)(1.3f * nPoints1);
        
        PairIntArray model = scene.copy();
        
        int translateX = 200;
        int translateY = 110;
        int x2Min = translateX;
        int y2Min = translateY;
        
        translateX(model, translateX);
        translateY(model, translateY);
        
        populateWithRandomPoints(sr, model, nPoints1, nPoints2, x2Min, y2Min,
            xRange, yRange);
        
        PointMatcher pointMatcher = new PointMatcher();
        
        TransformationPointFit fit = 
            pointMatcher.calculateProjectiveTransformation(
            scene, model, xRange, yRange);
        
        System.out.println("=> " + fit.toString());
        
        Transformer transformer = new Transformer();
        PairIntArray transformed = transformer.applyTransformation(
            fit.getParameters(), scene, xRange, yRange);
        
        overplotTransformed(transformed, model, (x2Min + xRange),
           (y2Min + yRange), 6);
        
        assertTrue(Math.abs(fit.getScale() - 1) < 0.1);
        assertTrue(Math.abs(fit.getRotationInRadians() - 0) <= 1.0);
        assertTrue(Math.abs(fit.getTranslationX() - translateX) <= 1.0);
        assertTrue(Math.abs(fit.getTranslationY() - translateY) <= 1.0);
        
    }
    
    @Test
    public void test7() throws Exception {

        // test for exact match rotated by 30 degrees
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1420149550374L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nPoints1 = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int x1Min = 0;
        int y1Min = 0;
        
        PairIntArray scene = createRandomPoints(sr, nPoints1,
            x1Min, y1Min, xRange, yRange);
        
        int nPoints2 = nPoints1;
        
        PairIntArray model = scene.copy();
        
        double rotation = 30. * Math.PI/180.;
        double scale = 1.;
        double translateX = 0;
        double translateY = 0;
        int x2Min = (int)translateX;
        int y2Min = (int)translateY;
        
        int centroidX1 = (x1Min + xRange) >> 1;
        int centroidY1 = (y1Min + yRange) >> 1;
        int centroidX2 = (x2Min + xRange) >> 1;
        int centroidY2 = (y2Min + yRange) >> 1;
        
        scaleAndRotate(scene, scale, -1*rotation, centroidX1, centroidY1);
        
        PointMatcher pointMatcher = new PointMatcher();
        
        TransformationPointFit fit = 
            pointMatcher.calculateProjectiveTransformation(
            scene, model, xRange, yRange);
                
        System.out.println("=> " + fit.toString());
        
        Transformer transformer = new Transformer();
        
        PairIntArray transformed = transformer.applyTransformation(
            fit.getParameters(), scene, centroidX1, centroidY1);
        
        overplotTransformed(transformed, model, xRange, yRange, 7);
        
        assertTrue(Math.abs(fit.getRotationInRadians() - rotation)
            < 1.0);
        assertTrue(Math.abs(fit.getScale() - fit.getScale()) < 1.0);
        assertTrue(Math.abs(fit.getTranslationX() - translateX) < 1.0);
        assertTrue(Math.abs(fit.getTranslationY() - translateY) < 1.0);
        assertTrue(fit.getNumberOfMatchedPoints() == 70);
        
    }
  
    @Test
    public void test8() throws Exception {

        // test for rotation plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nPoints1 = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int x1Min = 0;
        int y1Min = 0;
        
        PairIntArray scene = createRandomPoints(sr, nPoints1,
            x1Min, y1Min, xRange, yRange);
        
        int nPoints2 = (int)(1.3f * nPoints1);
        
        PairIntArray model = scene.copy();
        
        double rotation = 30. * Math.PI/180.;
        double scale = 1;
        int translateX = 0;
        int translateY = 0;
        int x2Min = translateX;
        int y2Min = translateY;
        
        int centroidX1 = (x1Min + xRange) >> 1;
        int centroidY1 = (y1Min + yRange) >> 1;
        int centroidX2 = (x2Min + xRange) >> 1;
        int centroidY2 = (y2Min + yRange) >> 1;
        
        populateWithRandomPoints(sr, model, nPoints1, nPoints2, x2Min, y2Min,
            xRange, yRange);
        
        scaleAndRotate(scene, scale, -1*rotation, centroidX1, centroidY1);
        
        PointMatcher pointMatcher = new PointMatcher();
        
        TransformationPointFit fit = 
            pointMatcher.calculateProjectiveTransformation(
            scene, model, xRange, yRange);
                
        System.out.println("=> " + fit.toString());
        
        Transformer transformer = new Transformer();
        
        PairIntArray transformed = transformer.applyTransformation(
            fit.getParameters(), scene, centroidX1, centroidY1);
        
        overplotTransformed(transformed, model, xRange, yRange, 8);
        
        assertTrue(Math.abs(fit.getRotationInRadians() - rotation)
            < 1.0);
        assertTrue(Math.abs(fit.getScale() - fit.getScale()) < 1.0);
        assertTrue(Math.abs(fit.getTranslationX() - translateX) < 1.0);
        assertTrue(Math.abs(fit.getTranslationY() - translateY) < 1.0);
        assertTrue(fit.getNumberOfMatchedPoints() == nPoints1);
    }
    
    @Test
    public void test9() throws Exception {

        // test for scale plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nPoints1 = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int x1Min = 0;
        int y1Min = 0;
        
        PairIntArray scene = createRandomPoints(sr, nPoints1,
            x1Min, y1Min, xRange, yRange);
        
        int nPoints2 = (int)(1.3f * nPoints1);
        
        PairIntArray model = scene.copy();
        
        double rotation = 0.;
        double scale = 4.;
        int translateX = 0;
        int translateY = 0;
        int x2Min = translateX;
        int y2Min = translateY;
        
        int centroidX1 = (x1Min + xRange) >> 1;
        int centroidY1 = (y1Min + yRange) >> 1;
        int centroidX2 = (x2Min + xRange) >> 1;
        int centroidY2 = (y2Min + yRange) >> 1;
        
        populateWithRandomPoints(sr, model, nPoints1, nPoints2, x2Min, y2Min,
            xRange, yRange);
        
        scaleAndRotate(scene, 1./scale, -1*rotation, centroidX1, centroidY1);
        
        PointMatcher pointMatcher = new PointMatcher();
        
        TransformationPointFit fit = 
            pointMatcher.calculateProjectiveTransformation(
            scene, model, xRange, yRange);

        System.out.println("=> " + fit.toString());
        
        Transformer transformer = new Transformer();
        
        PairIntArray transformed = transformer.applyTransformation(
            fit.getParameters(), scene, centroidX1, centroidY1);
        
        overplotTransformed(transformed, model, xRange, yRange, 9);
        
        assertTrue(Math.abs(fit.getRotationInRadians() - rotation) <= 1.0);
        assertTrue(Math.abs(fit.getScale() - fit.getScale()) < 1.0);
        assertTrue(Math.abs(fit.getTranslationX() - translateX) <= 1.0);
        assertTrue(Math.abs(fit.getTranslationY() - translateY) <= 1.0);
        assertTrue(fit.getNumberOfMatchedPoints() == nPoints1);
    }
    
    @Test
    public void test10() throws Exception {

        // test for scale smaller than 1, plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nPoints1 = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int x1Min = 0;
        int y1Min = 0;
        
        PairIntArray scene = createRandomPoints(sr, nPoints1,
            x1Min, y1Min, xRange, yRange);
        
        int nPoints2 = (int)(1.3f * nPoints1);
        
        PairIntArray model = scene.copy();
        
        double rotation = 0.;
        double scale = 1./4.;
        int translateX = 0;
        int translateY = 0;
        int x2Min = translateX;
        int y2Min = translateY;
        
        int centroidX1 = (x1Min + xRange) >> 1;
        int centroidY1 = (y1Min + yRange) >> 1;
        int centroidX2 = (x2Min + xRange) >> 1;
        int centroidY2 = (y2Min + yRange) >> 1;
        
        populateWithRandomPoints(sr, model, nPoints1, nPoints2, x2Min, y2Min,
            xRange, yRange);
        
        scaleAndRotate(scene, 1./scale, -1*rotation, centroidX1, centroidY1);
        
        PointMatcher pointMatcher = new PointMatcher();
        
        TransformationPointFit fit = 
            pointMatcher.calculateProjectiveTransformationWrapper(
            scene, model, xRange, yRange, xRange, yRange);

        System.out.println("=> " + fit.toString());
        
        Transformer transformer = new Transformer();
        
        PairIntArray transformed = transformer.applyTransformation(
            fit.getParameters(), scene, centroidX1, centroidY1);
        
        overplotTransformed(transformed, model, xRange, yRange, 9);
        
        assertTrue(Math.abs(fit.getRotationInRadians() - rotation) <= 1.0);
        assertTrue(Math.abs(fit.getScale() - fit.getScale()) < 1.0);
        assertTrue(Math.abs(fit.getTranslationX() - translateX) <= 1.0);
        assertTrue(Math.abs(fit.getTranslationY() - translateY) <= 1.0);
        assertTrue(fit.getNumberOfMatchedPoints() == nPoints1);
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
            
            /*
            tests for :
            -- for same set w/ rotation and translation and noise
            -- for same set w/ scale and rotation and translation and noise
            -- for same set w/ projection
            -- for same set w/ projection and noise
            
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
