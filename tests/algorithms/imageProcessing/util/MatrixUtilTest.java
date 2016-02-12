package algorithms.imageProcessing.util;

import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Logger;
import junit.framework.TestCase;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;
import org.ejml.simple.*;
    
/**
 *
 * @author nichole
 */
public class MatrixUtilTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public MatrixUtilTest() {
    }
    
    public void testDot() throws Exception {
       
        double[][] m1 = new double[2][3];
        m1[0] = new double[]{0, 1, 0};
        m1[1] = new double[]{1000, 100, 10};
        
        double[][] m2 = new double[3][2];
        m2[0] = new double[]{2, 1};
        m2[1] = new double[]{3, 0};
        m2[2] = new double[]{4, 0};
         
        /*
        0     1     0     2  1
        1000  100  10     3  0
                          4  0
        
        0*2    + 1*3   + 0*0     0*1    +  1*0   +  0*0
        1000*2 + 100*3 + 10*4    1000*1 +  100*0 + 10*0
        */
       
        double[][] m = MatrixUtil.dot(new SimpleMatrix(m1), 
            new SimpleMatrix(m2));
        
        assertTrue(m.length == 2);
        assertTrue(m[0].length == 2);
        assertTrue(m[1].length == 2);
        
        assertTrue(m[0][0] == 3);
        assertTrue(m[1][0] == 2340);
        assertTrue(m[0][1] == 0);
        assertTrue(m[1][1] == 1000);
    }
    
    public void testMultiply() throws Exception {

        double[][] a = new double[2][3];
        a[0] = new double[]{1, 2, 3};
        a[1] = new double[]{2, 3, 4};

        double[] b = new double[]{4, 3};

        double[] m = MatrixUtil.multiply(a, b);

        assertTrue( m[0] ==  10 );
        assertTrue( m[1] ==  17 );
        assertTrue( m[2] ==  24 );
     
        double[][] c = new double[3][];
        for (int i = 0; i < 3; i++) {
            c[i] = new double[3];
            for (int j = 0; j < 3; j++) {
                c[i][j] = 1;
            }
        }
        c[1][2] = 0;
        
           
        /*
        a b c      p0 p1 p2
        d e f      p3 p4 p5
                   p6 p7 p8        
        a*p0+... a*p a*p
        d*p0+... d*p d*p
        
        1 2 3    1 1 1
        2 3 4    1 1 0
                 1 1 1
        
        (1*1 + 2*1 + 3*1)  (1*1 + 2*1 + 3*1)  (1*1 + 0 + 3*1)
        (2*1 + 3*1 + 4*1)  (2*1 + 3*1 + 4*1)  (2*1 + 0 + 4*1)
        
        6  6 4
        9  9 6
        */
        
        double[][] d = MatrixUtil.multiply(a, c);

        assertTrue(d[0][0] == 6);
        assertTrue(d[0][1] == 6);
        assertTrue(d[0][2] == 4);
        assertTrue(d[1][0] == 9);
        assertTrue(d[1][1] == 9);
        assertTrue(d[1][2] == 6);
        
        /*
        example:  m is 1 2 3
                       4 5 6
                       7 8 9
        
                  n is 100  101  1
                       200  201  1
        
        multiply m by transpose of n:
        
        1 2 3     100  200
        4 5 6     101  201
        7 8 9      1   1
        
        (1*100 + 2*101 + 3*1)   (1*200 + 2*201 + 3*1)     305   605
        (4*100 + 5*101 + 6*1)   (4*200 + 5*201 + 6*1)  =  911  1811
        (7*100 + 8*101 + 9*1)   (7*200 + 8*201 + 9*1)    1517  3017
        */
        double[][] aa = new double[3][];
        aa[0] = new double[]{1, 2, 3};
        aa[1] = new double[]{4, 5, 6};
        aa[2] = new double[]{7, 8, 9};
        
        double[][] bb = new double[2][];
        bb[0] = new double[]{100, 101, 1};
        bb[1] = new double[]{200, 201, 1};
        
        double[][] cc = MatrixUtil.multiplyByTranspose(aa, bb);
        
        assertTrue(cc[0][0] == 305);
        assertTrue(cc[0][1] == 605);
        assertTrue(cc[1][0] == 911);
        assertTrue(cc[1][1] == 1811);
        assertTrue(cc[2][0] == 1517);
        assertTrue(cc[2][1] == 3017);
        
        int[] z = new int[]{0, 1, 2, 3};
        int factor = 2;
        int[] expectedZ = new int[]{0, 2, 4, 6};
        MatrixUtil.multiply(z, factor);
        for (int i = 0; i < z.length; i++) {
            assertTrue(z[i] == expectedZ[i]);
        }
    }

    public void testAdd() throws Exception {

        double[] a = new double[]{1, 2, 3, 4};
        double[] b = new double[]{100, 100, 100, 100};

        double[] expected = new double[]{101, 102, 103, 104};
        
        double[] c = MatrixUtil.add(a, b);
        
        assertTrue(Arrays.equals(expected, c));
    }
    
    public void testAdd2() throws Exception {

        float[] a = new float[]{1, 2, 3, 4};
        float[] b = new float[]{100, 100, 100, 100};

        float[] expected = new float[]{101, 102, 103, 104};
        
        float[] c = MatrixUtil.add(a, b);
        
        assertTrue(Arrays.equals(expected, c));
    }
    
    public void testAdd3() throws Exception {

        int[] a = new int[]{1, 2, 3, 4};
        int add = -1;
        
        int[] expected = new int[]{0, 1, 2, 3};
        
        MatrixUtil.add(a, add);
        
        assertTrue(Arrays.equals(expected, a));
    }
    
    public void testSubtract() throws Exception {

        float[] a = new float[]{100, 100, 100, 100};
        float[] b = new float[]{1, 2, 3, 4};

        float[] expected = new float[]{99, 98, 97, 96};
        
        float[] c = MatrixUtil.subtract(a, b);
        
        assertTrue(Arrays.equals(expected, c));
    }
    
    public void testTranspose() throws Exception {
        
        /*
        100  101  1    100  200
        200  201  1    101  201
                        1    1
        */
        float[][] bb = new float[2][];
        bb[0] = new float[]{100, 101, 1};
        bb[1] = new float[]{200, 201, 1};
        
        float[][] expected = new float[3][];
        expected[0] = new float[]{100, 200};
        expected[1] = new float[]{101, 201};
        expected[2] = new float[]{1, 1};
        
        float[][] cc = MatrixUtil.transpose(bb);
        
        for (int i = 0; i < cc.length; i++) {
            for (int j = 0; j < cc[i].length; j++) {
                assertTrue(expected[i][j] == cc[i][j]);
            }
        }
        
        float[][] dd = MatrixUtil.transpose(cc);
        
        for (int i = 0; i < dd.length; i++) {
            for (int j = 0; j < dd[i].length; j++) {
                assertTrue(bb[i][j] == dd[i][j]);
            }
        }
        
    }
    
    public void testScaleToUnitVariance() throws Exception {
        
        SimpleMatrix[] dataAndClasses = readIrisDataset();
        
        double v0 = dataAndClasses[0].get(0, 0);
        double v1 = dataAndClasses[0].get(1, 0);
        double v2 = dataAndClasses[0].get(2, 0);
        double v3 = dataAndClasses[0].get(3, 0);
        double[] expected = new double[]{5.1,3.5,1.4,0.2};
        assertTrue(Math.abs(v0 - expected[0]) <  0.05*Math.abs(expected[0]));
        assertTrue(Math.abs(v1 - expected[1]) <  0.05*Math.abs(expected[1]));
        assertTrue(Math.abs(v2 - expected[2]) <  0.05*Math.abs(expected[2]));
        assertTrue(Math.abs(v3 - expected[3]) <  0.05*Math.abs(expected[3]));
        
        SimpleMatrix normData = MatrixUtil.scaleToUnitStandardDeviation(
            dataAndClasses[0]);
        
        /*
        assert first 4
         [ -9.0068e-01   1.0321e+00  -1.3413e+00  -1.3130e+00]
         [ -1.1430e+00  -1.2496e-01  -1.3413e+00  -1.3130e+00]
         [ -1.3854e+00   3.3785e-01  -1.3981e+00  -1.3130e+00]
         [ -1.5065e+00   1.0645e-01  -1.2844e+00  -1.3130e+00]
         [ -1.0218e+00   1.2635e+00  -1.3413e+00  -1.3130e+00]        
        */
        v0 = normData.get(0, 0);
        v1 = normData.get(1, 0);
        v2 = normData.get(2, 0);
        v3 = normData.get(3, 0);
        expected = new double[]{-9.0068e-01, 1.0321e+00, -1.3413e+00, -1.3130e+00};
        assertTrue(Math.abs(v0 - expected[0]) <  0.05*Math.abs(expected[0]));
        assertTrue(Math.abs(v1 - expected[1]) <  0.05*Math.abs(expected[1]));
        assertTrue(Math.abs(v2 - expected[2]) <  0.05*Math.abs(expected[2]));
        assertTrue(Math.abs(v3 - expected[3]) <  0.05*Math.abs(expected[3]));
        
        // assert mean = 0
        // assert var = 1
        
        int n = normData.numCols();
        int nRows = normData.numRows();
                
        double[] mean = new double[nRows];
        for (int i = 0; i < n; ++i) {            
            for (int j = 0; j < nRows; ++j) {
                mean[j] += normData.get(j, i);
            }
        }
        for (int j = 0; j < nRows; ++j) {
            mean[j] /= (double)n;
            assertTrue(Math.abs(mean[j]) < 0.1);
        }
        
        double[] stdev = new double[nRows];
        for (int i = 0; i < n; ++i) {            
            for (int j = 0; j < nRows; ++j) {
                double d = normData.get(j, i) - mean[j];
                stdev[j] += (d * d);
            }
        }
        for (int j = 0; j < nRows; ++j) {
            stdev[j] = Math.sqrt(stdev[j]/(double)(n - 1));
            assertTrue(Math.abs(stdev[j] - 1.) < 0.1);
        }
        
    }
    
    public void testCreateLDATrasformation() throws Exception {
        
        SimpleMatrix[] dataAndClasses = readIrisDataset();
        SimpleMatrix classes = dataAndClasses[1].copy();
        
        SimpleMatrix w = MatrixUtil.createLDATransformation(dataAndClasses[0], 
            dataAndClasses[1]);
        
        assertEquals(2, w.numRows());
        assertEquals(4, w.numCols());
        
        assertTrue(Math.abs(w.get(0, 0) - 0.15) < 0.01);
        assertTrue(Math.abs(w.get(0, 1) - 0.148) < 0.01);
        assertTrue(Math.abs(w.get(0, 2) - -0.851) < 0.01);
        assertTrue(Math.abs(w.get(0, 3) - -0.481) < 0.01);
        
        assertTrue(Math.abs(w.get(1, 0) - 0.010) < 0.01);
        assertTrue(Math.abs(w.get(1, 1) - 0.327) < 0.01);
        assertTrue(Math.abs(w.get(1, 2) - -0.575) < 0.01);
        assertTrue(Math.abs(w.get(1, 3) - 0.750) < 0.01);
        
        
        SimpleMatrix normData = MatrixUtil.scaleToUnitStandardDeviation(dataAndClasses[0]);
                               
        // transforms from integer classes to zero based counting with delta of 1
        // for example:  [1, 2, 5, ...] becomes [0, 1, 2, ...]
        int nClasses = MatrixUtil.transformToZeroBasedClasses(classes);
        
        SimpleMatrix w2 = MatrixUtil.createLDATransformation2(normData, classes, nClasses);
        
        assertEquals(2, w2.numRows());
        assertEquals(4, w2.numCols());
        
        assertTrue(Math.abs(w2.get(0, 0) - 0.15) < 0.01);
        assertTrue(Math.abs(w2.get(0, 1) - 0.148) < 0.01);
        assertTrue(Math.abs(w2.get(0, 2) - -0.851) < 0.01);
        assertTrue(Math.abs(w2.get(0, 3) - -0.481) < 0.01);
        
        assertTrue(Math.abs(w2.get(1, 0) - 0.010) < 0.01);
        assertTrue(Math.abs(w2.get(1, 1) - 0.327) < 0.01);
        assertTrue(Math.abs(w2.get(1, 2) - -0.575) < 0.01);
        assertTrue(Math.abs(w2.get(1, 3) - 0.750) < 0.01);
        
        
        int nr = w2.numRows();
        int nc = w2.numCols();
        int nr2 = normData.numRows();
        int nc2 = normData.numCols();        
        // 2 X 150
        SimpleMatrix dataTransformed = new SimpleMatrix(MatrixUtil.dot(w, normData));
        
        float minX = Float.MAX_VALUE;
        float maxX = Float.MIN_VALUE;
        float minY = Float.MAX_VALUE;
        float maxY = Float.MIN_VALUE;
        int[] countClasses = new int[nClasses];
        for (int col = 0; col < dataTransformed.numCols(); ++col) {
            int k = (int)Math.round(classes.get(0, col));
            countClasses[k]++;
            float x = (float)dataTransformed.get(0, col);
            float y = (float)dataTransformed.get(1, col);
            if (x < minX) {
                minX = x;
            }
            if (y < minY) {
                minY = y;
            }
            if (x > maxX) {
                maxX = x;
            }
            if (y > maxY) {
                maxY = y;
            }
        }
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(minX, maxX, 
            minY, maxY);
        
        for (int k = 0; k < nClasses; ++k) {
            float[] xPoint = new float[countClasses[k]];
            float[] yPoint = new float[countClasses[k]];
            int count = 0;
            for (int col = 0; col < dataTransformed.numCols(); ++col) {
                if ((int)Math.round(classes.get(0, col)) != k) {
                    continue;
                }
                xPoint[count] = (float)dataTransformed.get(0, col);
                yPoint[count] = (float)dataTransformed.get(1, col);
                count++;
            }
            float[] xPoly = null;
            float[] yPoly = null;
            plotter.addPlot(xPoint, yPoint, xPoly, yPoly, "class " + k);
        }
        String file1 = plotter.writeFile();
        
        /*System.out.println(String.format("(%.3f, %.3f)", 
                (float)dataTransformed.get(0, 0), 
                (float)dataTransformed.get(1, 0)));
        System.out.println(String.format("(%.3f, %.3f)", 
                (float)dataTransformed.get(0, 1), 
                (float)dataTransformed.get(1, 1)));
        System.out.println(String.format("(%.3f, %.3f)", 
                (float)dataTransformed.get(0, 2), 
                (float)dataTransformed.get(1, 2)));*/
        
        assertTrue(Math.abs(dataTransformed.get(0, 0) - 1.791) < 0.01);
        assertTrue(Math.abs(dataTransformed.get(1, 0) - 0.115) < 0.01);
        
        assertTrue(Math.abs(dataTransformed.get(0, 1) - 1.583) < 0.01);
        assertTrue(Math.abs(dataTransformed.get(1, 1) - -0.265) < 0.01);
        
        assertTrue(Math.abs(dataTransformed.get(0, 2) - 1.664) < 0.01);
        assertTrue(Math.abs(dataTransformed.get(1, 2) - -0.084) < 0.01);

        // --- to make a transformation usable on features not normalized:
        SimpleMatrix w3 = MatrixUtil.createLDATransformation(
            dataAndClasses[0], dataAndClasses[1]);
        SimpleMatrix dataTransformed3 = new SimpleMatrix(MatrixUtil.dot(w3, 
            dataAndClasses[0]));
        
        minX = Float.MAX_VALUE;
        maxX = Float.MIN_VALUE;
        minY = Float.MAX_VALUE;
        maxY = Float.MIN_VALUE;
        countClasses = new int[nClasses];
        for (int col = 0; col < dataTransformed3.numCols(); ++col) {
            int k = (int)Math.round(classes.get(0, col));
            countClasses[k]++;
            float x = (float)dataTransformed3.get(0, col);
            float y = (float)dataTransformed3.get(1, col);
            if (x < minX) {
                minX = x;
            }
            if (y < minY) {
                minY = y;
            }
            if (x > maxX) {
                maxX = x;
            }
            if (y > maxY) {
                maxY = y;
            }
        }
        
        PolygonAndPointPlotter plotter2 = new PolygonAndPointPlotter(minX, maxX, 
            minY, maxY);
        
        for (int k = 0; k < nClasses; ++k) {
            float[] xPoint = new float[countClasses[k]];
            float[] yPoint = new float[countClasses[k]];
            int count = 0;
            for (int col = 0; col < dataTransformed3.numCols(); ++col) {
                if ((int)Math.round(classes.get(0, col)) != k) {
                    continue;
                }
                xPoint[count] = (float)dataTransformed3.get(0, col);
                yPoint[count] = (float)dataTransformed3.get(1, col);
                count++;
            }
            float[] xPoly = null;
            float[] yPoly = null;
            plotter2.addPlot(xPoint, yPoint, xPoly, yPoly, "class " + k);
        }
        String file2 = plotter2.writeFile2();
        
       
    }
    
    public void testCreatePCATrasformation() throws Exception {
        
        SimpleMatrix[] dataAndClasses = readHalfMoonDataset();
        SimpleMatrix classes = dataAndClasses[1].copy();
        
        SimpleMatrix w = MatrixUtil.createRBFKernelPCATransformation2(
            dataAndClasses[0]);
        
        float minX = Float.MAX_VALUE;
        float maxX = Float.MIN_VALUE;
        float minY = Float.MAX_VALUE;
        float maxY = Float.MIN_VALUE;
        
        int n = w.numCols();
        
        float[] xPoint = new float[n];
        float[] yPoint = new float[n];
        int count = 0;
        for (int col = 0; col < n; ++col) {
            float x = (float)w.get(0, col);
            float y = (float)w.get(1, col);
            xPoint[count] = x;
            yPoint[count] = y;
            count++;
            if (x < minX) {
                minX = x;
            }
            if (y < minY) {
                minY = y;
            }
            if (x > maxX) {
                maxX = x;
            }
            if (y > maxY) {
                maxY = y;
            }
        }
        
        PolygonAndPointPlotter plotter2 = new PolygonAndPointPlotter(minX, maxX, 
            minY, maxY);
        
        float[] xPoly = null;
        float[] yPoly = null;
        plotter2.addPlot(xPoint, yPoint, xPoly, yPoly, "rbf");

        String file = plotter2.writeFile();
        int z = 1;
    }
    
    public void testLDATrasformation() throws Exception {
        
        SimpleMatrix[] dataAndClasses = readSegmentationDataset();
        SimpleMatrix classes = dataAndClasses[1].copy();
        
        int nClasses = 2;
              
        float minX = Float.MAX_VALUE;
        float maxX = Float.MIN_VALUE;
        float minY = Float.MAX_VALUE;
        float maxY = Float.MIN_VALUE;
        int[] countClasses = new int[nClasses];
        
        // --- to make a transformation usable on features not normalized:
        SimpleMatrix w3 = MatrixUtil.createLDATransformation(
            dataAndClasses[0], dataAndClasses[1]);
        SimpleMatrix dataTransformed3 = new SimpleMatrix(MatrixUtil.dot(w3, 
            dataAndClasses[0]));
        
        for (int col = 0; col < dataTransformed3.numCols(); ++col) {
            int k = (int)Math.round(classes.get(0, col));
            countClasses[k]++;
            float x = (float)dataTransformed3.get(0, col);
            float y = (float)dataTransformed3.get(1, col);
            if (x < minX) {
                minX = x;
            }
            if (y < minY) {
                minY = y;
            }
            if (x > maxX) {
                maxX = x;
            }
            if (y > maxY) {
                maxY = y;
            }
        }
        
        PolygonAndPointPlotter plotter2 = new PolygonAndPointPlotter(minX, maxX, 
            minY, maxY);
        
        for (int k = 0; k < nClasses; ++k) {
            float[] xPoint = new float[countClasses[k]];
            float[] yPoint = new float[countClasses[k]];
            int count = 0;
            for (int col = 0; col < dataTransformed3.numCols(); ++col) {
                if ((int)Math.round(classes.get(0, col)) != k) {
                    continue;
                }
                xPoint[count] = (float)dataTransformed3.get(0, col);
                yPoint[count] = (float)dataTransformed3.get(1, col);
                count++;
            }
            float[] xPoly = null;
            float[] yPoly = null;
            plotter2.addPlot(xPoint, yPoint, xPoly, yPoly, "class " + k);
        }
        String file2 = plotter2.writeFile2();
        
        log.info("transfromation matrix = \n" + w3.toString());
        
        int z = 1;
    }
    
    private SimpleMatrix[] readSegmentationDataset() throws Exception {
        
        BufferedReader bReader = null;
        FileReader reader = null;
        
        String filePath = ResourceFinder.findFileInTestResources("segmentation_colors1.csv");
        
        try {
            reader = new FileReader(new File(filePath));
            
            bReader = new BufferedReader(reader);
            
            SimpleMatrix data = new SimpleMatrix(4, 28);
            SimpleMatrix classes = new SimpleMatrix(1, 28);
            
            String line = bReader.readLine();
            line = bReader.readLine();
            
            int count = 0;
            
            while (line != null) {
                
                String[] items = line.split(",");
                if (items.length != 5) {
                    throw new IllegalStateException("expecting 5 items in a line");
                }
 
                for (int j = 0; j < 4; ++j) {
                    data.set(j, count, Double.valueOf(items[j]));
                }
                
                double classValue = 0;
                if (items[4].equals("1")) {
                    classValue = 0;
                } else if (items[4].equals("2")) {
                    classValue = 1;
                }
                
                classes.set(count, classValue);
                
                line = bReader.readLine();
                
                count++;
            }
            
            return new SimpleMatrix[]{data, classes};
            
        } catch (IOException e) {
            log.severe(e.getMessage());
        } finally {
            if (reader == null) {
                reader.close();
            }
            if (bReader == null) {
                bReader.close();
            }
        }
        
        return null;
    }
    
    private SimpleMatrix[] readIrisDataset() throws Exception {
        
        BufferedReader bReader = null;
        FileReader reader = null;
        
        String filePath = ResourceFinder.findFileInTestResources("iris.data");
        
        try {
            reader = new FileReader(new File(filePath));
            
            bReader = new BufferedReader(reader);
            
            SimpleMatrix data = new SimpleMatrix(4, 150);
            SimpleMatrix classes = new SimpleMatrix(1, 150);
            
            String line = bReader.readLine();
            
            int count = 0;
            
            while (line != null) {
                
                String[] items = line.split(",");
                if (items.length != 5) {
                    throw new IllegalStateException("expecting 5 items in a line");
                }
 
                for (int j = 0; j < 4; ++j) {
                    data.set(j, count, Double.valueOf(items[j]));
                }
                
                double classValue = 4;
                if (items[4].equals("Iris-setosa")) {
                    classValue = 1;
                } else if (items[4].equals("Iris-versicolor")) {
                    classValue = 2;
                }
                //Iris-virginica
                
                classes.set(count, classValue);
                
                line = bReader.readLine();
                
                count++;
            }
            
            return new SimpleMatrix[]{data, classes};
            
        } catch (IOException e) {
            log.severe(e.getMessage());
        } finally {
            if (reader == null) {
                reader.close();
            }
            if (bReader == null) {
                bReader.close();
            }
        }
        
        return null;
    }
    
    private SimpleMatrix[] readHalfMoonDataset() throws Exception {
        
        SimpleMatrix data = new SimpleMatrix(2, 100);
        SimpleMatrix classes = new SimpleMatrix(1, 100);
        
        data.set( 0 , 0 , 0.871318704123 );
        data.set( 1 , 0 , 0.490717552004 );
        data.set( 0 , 1 , 0.715472413369 );
        data.set( 1 , 1 , -0.458667853037 );
        data.set( 0 , 2 , 1.46253829024 );
        data.set( 1 , 2 , -0.386599306373 );
        data.set( 0 , 3 , -0.222520933956 );
        data.set( 1 , 3 , 0.974927912182 );
        data.set( 0 , 4 , 0.327699109739 );
        data.set( 1 , 4 , -0.240277997075 );
        data.set( 0 , 5 , 1.0 );
        data.set( 1 , 5 , 0.0 );
        data.set( 0 , 6 , 0.949055747011 );
        data.set( 1 , 6 , 0.315108218024 );
        data.set( 0 , 7 , 0.0 );
        data.set( 1 , 7 , 0.5 );
        data.set( 0 , 8 , 1.40478334312 );
        data.set( 1 , 8 , -0.414412623016 );
        data.set( 0 , 9 , 0.967294863039 );
        data.set( 1 , 9 , 0.25365458391 );
        data.set( 0 , 10 , 0.0960230259077 );
        data.set( 1 , 10 , 0.995379112949 );
        data.set( 0 , 11 , 0.427883339878 );
        data.set( 1 , 11 , -0.320172254597 );
        data.set( 0 , 12 , 1.09602302591 );
        data.set( 1 , 12 , -0.495379112949 );
        data.set( 0 , 13 , 0.198586378132 );
        data.set( 1 , 13 , -0.0981105304912 );
        data.set( 0 , 14 , 0.0320515775717 );
        data.set( 1 , 14 , 0.999486216201 );
        data.set( 0 , 15 , -0.900968867902 );
        data.set( 1 , 15 , 0.433883739118 );
        data.set( 0 , 16 , 1.15959989503 );
        data.set( 1 , 16 , -0.487181783414 );
        data.set( 0 , 17 , -0.761445958369 );
        data.set( 1 , 17 , 0.648228395308 );
        data.set( 0 , 18 , 0.073083242654 );
        data.set( 1 , 18 , 0.124732995121 );
        data.set( 0 , 19 , 1.03205157757 );
        data.set( 1 , 19 , -0.499486216201 );
        data.set( 0 , 20 , -0.623489801859 );
        data.set( 1 , 20 , 0.781831482468 );
        data.set( 0 , 21 , 1.76144595837 );
        data.set( 1 , 21 , -0.148228395308 );
        data.set( 0 , 22 , 0.345365054421 );
        data.set( 1 , 22 , 0.93846842205 );
        data.set( 0 , 23 , -0.284527586631 );
        data.set( 1 , 23 , 0.958667853037 );
        data.set( 0 , 24 , -0.404783343122 );
        data.set( 1 , 24 , 0.914412623016 );
        data.set( 0 , 25 , 1.87131870412 );
        data.set( 1 , 25 , 0.00928244799606 );
        data.set( 0 , 26 , 1.62348980186 );
        data.set( 1 , 26 , -0.281831482468 );
        data.set( 0 , 27 , 0.838088104892 );
        data.set( 1 , 27 , 0.545534901211 );
        data.set( 0 , 28 , 0.0184408430089 );
        data.set( 1 , 28 , 0.308841371299 );
        data.set( 0 , 29 , -0.871318704123 );
        data.set( 1 , 29 , 0.490717552004 );
        data.set( 0 , 30 , 0.222520933956 );
        data.set( 1 , 30 , 0.974927912182 );
        data.set( 0 , 31 , 1.83808810489 );
        data.set( 1 , 31 , -0.0455349012105 );
        data.set( 0 , 32 , -0.518392568311 );
        data.set( 1 , 32 , 0.855142763005 );
        data.set( 0 , 33 , 0.654634945579 );
        data.set( 1 , 33 , -0.43846842205 );
        data.set( 0 , 34 , 1.57211666012 );
        data.set( 1 , 34 , -0.320172254597 );
        data.set( 0 , 35 , 1.7183493501 );
        data.set( 1 , 35 , -0.195682550603 );
        data.set( 0 , 36 , 1.96729486304 );
        data.set( 1 , 36 , 0.24634541609 );
        data.set( 0 , 37 , 1.99179001382 );
        data.set( 1 , 37 , 0.372122838315 );
        data.set( 0 , 38 , 0.281650649902 );
        data.set( 1 , 38 , -0.195682550603 );
        data.set( 0 , 39 , 0.718349350098 );
        data.set( 1 , 39 , 0.695682550603 );
        data.set( 0 , 40 , 0.284527586631 );
        data.set( 1 , 40 , 0.958667853037 );
        data.set( 0 , 41 , 1.80141362187 );
        data.set( 1 , 41 , -0.0981105304912 );
        data.set( 0 , 42 , -0.718349350098 );
        data.set( 1 , 42 , 0.695682550603 );
        data.set( 0 , 43 , 0.161911895108 );
        data.set( 1 , 43 , -0.0455349012105 );
        data.set( 0 , 44 , 0.99794539275 );
        data.set( 1 , 44 , 0.0640702199807 );
        data.set( 0 , 45 , 0.967948422428 );
        data.set( 1 , 45 , -0.499486216201 );
        data.set( 0 , 46 , 0.761445958369 );
        data.set( 1 , 46 , 0.648228395308 );
        data.set( 0 , 47 , 1.28452758663 );
        data.set( 1 , 47 , -0.458667853037 );
        data.set( 0 , 48 , 0.623489801859 );
        data.set( 1 , 48 , 0.781831482468 );
        data.set( 0 , 49 , 0.032705136961 );
        data.set( 1 , 49 , 0.24634541609 );
        data.set( 0 , 50 , 0.518392568311 );
        data.set( 1 , 50 , 0.855142763005 );
        data.set( 0 , 51 , -0.0960230259077 );
        data.set( 1 , 51 , 0.995379112949 );
        data.set( 0 , 52 , 0.00205460724966 );
        data.set( 1 , 52 , 0.435929780019 );
        data.set( 0 , 53 , -0.967294863039 );
        data.set( 1 , 53 , 0.25365458391 );
        data.set( 0 , 54 , 0.926916757346 );
        data.set( 1 , 54 , 0.375267004879 );
        data.set( 0 , 55 , 1.99794539275 );
        data.set( 1 , 55 , 0.435929780019 );
        data.set( 0 , 56 , -0.345365054421 );
        data.set( 1 , 56 , 0.93846842205 );
        data.set( 0 , 57 , -0.949055747011 );
        data.set( 1 , 57 , 0.315108218024 );
        data.set( 0 , 58 , 0.840400104967 );
        data.set( 1 , 58 , -0.487181783414 );
        data.set( 0 , 59 , -0.926916757346 );
        data.set( 1 , 59 , 0.375267004879 );
        data.set( 0 , 60 , 0.572116660122 );
        data.set( 1 , 60 , 0.820172254597 );
        data.set( 0 , 61 , 1.94905574701 );
        data.set( 1 , 61 , 0.184891781976 );
        data.set( 0 , 62 , 0.404783343122 );
        data.set( 1 , 62 , 0.914412623016 );
        data.set( 0 , 63 , 0.672300890261 );
        data.set( 1 , 63 , 0.740277997075 );
        data.set( 0 , 64 , 0.159599895033 );
        data.set( 1 , 64 , 0.987181783414 );
        data.set( 0 , 65 , 0.801413621868 );
        data.set( 1 , 65 , 0.598110530491 );
        data.set( 0 , 66 , 0.128681295877 );
        data.set( 1 , 66 , 0.00928244799606 );
        data.set( 0 , 67 , 0.777479066044 );
        data.set( 1 , 67 , -0.474927912182 );
        data.set( 0 , 68 , 0.376510198141 );
        data.set( 1 , 68 , -0.281831482468 );
        data.set( 0 , 69 , 0.981559156991 );
        data.set( 1 , 69 , 0.191158628701 );
        data.set( 0 , 70 , -0.838088104892 );
        data.set( 1 , 70 , 0.545534901211 );
        data.set( 0 , 71 , -0.572116660122 );
        data.set( 1 , 71 , 0.820172254597 );
        data.set( 0 , 72 , -0.159599895033 );
        data.set( 1 , 72 , 0.987181783414 );
        data.set( 0 , 73 , 0.00820998617675 );
        data.set( 1 , 73 , 0.372122838315 );
        data.set( 0 , 74 , 0.900968867902 );
        data.set( 1 , 74 , 0.433883739118 );
        data.set( 0 , 75 , -0.99794539275 );
        data.set( 1 , 75 , 0.0640702199807 );
        data.set( 0 , 76 , 0.238554041631 );
        data.set( 1 , 76 , -0.148228395308 );
        data.set( 0 , 77 , 1.92691675735 );
        data.set( 1 , 77 , 0.124732995121 );
        data.set( 0 , 78 , 2.0 );
        data.set( 1 , 78 , 0.5 );
        data.set( 0 , 79 , -0.801413621868 );
        data.set( 1 , 79 , 0.598110530491 );
        data.set( 0 , 80 , 0.991790013823 );
        data.set( 1 , 80 , 0.127877161685 );
        data.set( 0 , 81 , 0.537461709759 );
        data.set( 1 , 81 , -0.386599306373 );
        data.set( 0 , 82 , 0.0509442529893 );
        data.set( 1 , 82 , 0.184891781976 );
        data.set( 0 , 83 , -1.0 );
        data.set( 1 , 83 , 1.22464679915e-16 );
        data.set( 0 , 84 , 0.595216656878 );
        data.set( 1 , 84 , -0.414412623016 );
        data.set( 0 , 85 , 1.34536505442 );
        data.set( 1 , 85 , -0.43846842205 );
        data.set( 0 , 86 , -0.672300890261 );
        data.set( 1 , 86 , 0.740277997075 );
        data.set( 0 , 87 , 1.22252093396 );
        data.set( 1 , 87 , -0.474927912182 );
        data.set( 0 , 88 , 1.98155915699 );
        data.set( 1 , 88 , 0.308841371299 );
        data.set( 0 , 89 , -0.0320515775717 );
        data.set( 1 , 89 , 0.999486216201 );
        data.set( 0 , 90 , -0.981559156991 );
        data.set( 1 , 90 , 0.191158628701 );
        data.set( 0 , 91 , -0.462538290241 );
        data.set( 1 , 91 , 0.886599306373 );
        data.set( 0 , 92 , 0.903976974092 );
        data.set( 1 , 92 , -0.495379112949 );
        data.set( 0 , 93 , -0.991790013823 );
        data.set( 1 , 93 , 0.127877161685 );
        data.set( 0 , 94 , 1.67230089026 );
        data.set( 1 , 94 , -0.240277997075 );
        data.set( 0 , 95 , 0.0990311320976 );
        data.set( 1 , 95 , 0.0661162608824 );
        data.set( 0 , 96 , 1.51839256831 );
        data.set( 1 , 96 , -0.355142763005 );
        data.set( 0 , 97 , 0.462538290241 );
        data.set( 1 , 97 , 0.886599306373 );
        data.set( 0 , 98 , 1.9009688679 );
        data.set( 1 , 98 , 0.0661162608824 );
        data.set( 0 , 99 , 0.481607431689 );
        data.set( 1 , 99 , -0.355142763005 );
            
        classes.set( 0 , 0 , 0 );
        classes.set( 0 , 1 , 1 );
        classes.set( 0 , 2 , 1 );
        classes.set( 0 , 3 , 0 );
        classes.set( 0 , 4 , 1 );
        classes.set( 0 , 5 , 0 );
        classes.set( 0 , 6 , 0 );
        classes.set( 0 , 7 , 1 );
        classes.set( 0 , 8 , 1 );
        classes.set( 0 , 9 , 0 );
        classes.set( 0 , 10 , 0 );
        classes.set( 0 , 11 , 1 );
        classes.set( 0 , 12 , 1 );
        classes.set( 0 , 13 , 1 );
        classes.set( 0 , 14 , 0 );
        classes.set( 0 , 15 , 0 );
        classes.set( 0 , 16 , 1 );
        classes.set( 0 , 17 , 0 );
        classes.set( 0 , 18 , 1 );
        classes.set( 0 , 19 , 1 );
        classes.set( 0 , 20 , 0 );
        classes.set( 0 , 21 , 1 );
        classes.set( 0 , 22 , 0 );
        classes.set( 0 , 23 , 0 );
        classes.set( 0 , 24 , 0 );
        classes.set( 0 , 25 , 1 );
        classes.set( 0 , 26 , 1 );
        classes.set( 0 , 27 , 0 );
        classes.set( 0 , 28 , 1 );
        classes.set( 0 , 29 , 0 );
        classes.set( 0 , 30 , 0 );
        classes.set( 0 , 31 , 1 );
        classes.set( 0 , 32 , 0 );
        classes.set( 0 , 33 , 1 );
        classes.set( 0 , 34 , 1 );
        classes.set( 0 , 35 , 1 );
        classes.set( 0 , 36 , 1 );
        classes.set( 0 , 37 , 1 );
        classes.set( 0 , 38 , 1 );
        classes.set( 0 , 39 , 0 );
        classes.set( 0 , 40 , 0 );
        classes.set( 0 , 41 , 1 );
        classes.set( 0 , 42 , 0 );
        classes.set( 0 , 43 , 1 );
        classes.set( 0 , 44 , 0 );
        classes.set( 0 , 45 , 1 );
        classes.set( 0 , 46 , 0 );
        classes.set( 0 , 47 , 1 );
        classes.set( 0 , 48 , 0 );
        classes.set( 0 , 49 , 1 );
        classes.set( 0 , 50 , 0 );
        classes.set( 0 , 51 , 0 );
        classes.set( 0 , 52 , 1 );
        classes.set( 0 , 53 , 0 );
        classes.set( 0 , 54 , 0 );
        classes.set( 0 , 55 , 1 );
        classes.set( 0 , 56 , 0 );
        classes.set( 0 , 57 , 0 );
        classes.set( 0 , 58 , 1 );
        classes.set( 0 , 59 , 0 );
        classes.set( 0 , 60 , 0 );
        classes.set( 0 , 61 , 1 );
        classes.set( 0 , 62 , 0 );
        classes.set( 0 , 63 , 0 );
        classes.set( 0 , 64 , 0 );
        classes.set( 0 , 65 , 0 );
        classes.set( 0 , 66 , 1 );
        classes.set( 0 , 67 , 1 );
        classes.set( 0 , 68 , 1 );
        classes.set( 0 , 69 , 0 );
        classes.set( 0 , 70 , 0 );
        classes.set( 0 , 71 , 0 );
        classes.set( 0 , 72 , 0 );
        classes.set( 0 , 73 , 1 );
        classes.set( 0 , 74 , 0 );
        classes.set( 0 , 75 , 0 );
        classes.set( 0 , 76 , 1 );
        classes.set( 0 , 77 , 1 );
        classes.set( 0 , 78 , 1 );
        classes.set( 0 , 79 , 0 );
        classes.set( 0 , 80 , 0 );
        classes.set( 0 , 81 , 1 );
        classes.set( 0 , 82 , 1 );
        classes.set( 0 , 83 , 0 );
        classes.set( 0 , 84 , 1 );
        classes.set( 0 , 85 , 1 );
        classes.set( 0 , 86 , 0 );
        classes.set( 0 , 87 , 1 );
        classes.set( 0 , 88 , 1 );
        classes.set( 0 , 89 , 0 );
        classes.set( 0 , 90 , 0 );
        classes.set( 0 , 91 , 0 );
        classes.set( 0 , 92 , 1 );
        classes.set( 0 , 93 , 0 );
        classes.set( 0 , 94 , 1 );
        classes.set( 0 , 95 , 1 );
        classes.set( 0 , 96 , 1 );
        classes.set( 0 , 97 , 0 );
        classes.set( 0 , 98 , 1 );
        classes.set( 0 , 99 , 1 );
        return new SimpleMatrix[]{data, classes};
    }
}
