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
            if (x < maxX) {
                maxX = x;
            }
            if (y > maxY) {
                maxY = y;
            }
        }
        maxX = 3;
        maxY = 3;
        
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
            if (x < maxX) {
                maxX = x;
            }
            if (y > maxY) {
                maxY = y;
            }
        }
        maxX = -1*minX;
        maxY = -1*minY;
        
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
}
