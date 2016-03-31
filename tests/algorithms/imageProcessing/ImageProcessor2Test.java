package algorithms.imageProcessing;

import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.Complex;
import algorithms.misc.MedianSmooth;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.Arrays;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;
import thirdparty.ca.uol.aig.fftpack.Complex1D;
import thirdparty.ca.uol.aig.fftpack.ComplexDoubleFFT;

/**
 *
 * @author nichole
 */
public class ImageProcessor2Test extends TestCase {

    public ImageProcessor2Test(String testName) {
        super(testName);
    }
    
    public void test0_2DFFT_compare() throws Exception {
        
        /*
        usng a 2-d signum function
        */
        
        int nRows = 4;
        int nCols = 4;
        
        Complex[][] a = new Complex[nRows][];
        for (int i = 0; i < nRows; ++i) {
            a[i] = new Complex[nCols];
            a[i][0] = new Complex(1, 0);
            a[i][1] = new Complex(1, 0);
            a[i][2] = new Complex(-1, 0);
            a[i][3] = new Complex(-1, 0);
        }
        
        Complex1D[] a2 = new Complex1D[nRows];
        for (int i = 0; i < nRows; ++i) {
            a2[i] = new Complex1D();
            a2[i].x = new double[nCols];
            a2[i].y = new double[nCols];
            a2[i].x[0] = 1;
            a2[i].y[0] = 0;
            a2[i].x[1] = 1;
            a2[i].y[1] = 0;
            a2[i].x[2] = -1;
            a2[i].y[2] = 0;
            a2[i].x[3] = -1;
            a2[i].y[3] = 0;
        }
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        // --- perform unnormalized fft transforms -----
        Complex[][] ccTr1 = imageProcessor.create2DFFT(a, false, true);
        Complex1D[] ccTr2 = imageProcessor.create2DFFT2(a2, false, true);  
        assertEquals(ccTr1.length, ccTr2.length);
        assertEquals(ccTr1[0].length, ccTr2[0].x.length);
        for (int row = 0; row < ccTr1.length; ++row) {
            for (int col = 0; col < ccTr1[0].length; ++col) {
                double r1 = ccTr1[row][col].re();
                double im1 = ccTr1[row][col].im();
                
                double r2 = ccTr2[row].x[col];
                double im2 = ccTr2[row].y[col];
                
                System.out.println(String.format("** %f %f", (float)r1, (float)r2));
                
                assertTrue(Math.abs(r2 - r1) < 0.01);
                assertTrue(Math.abs(im2 - im1) < 0.01);
            }
        } 
        
        double norm = a.length * a[0].length;
        
        // --- perform unnormalized inverse fft transforms -----
        Complex[][] ccTrTr1 = imageProcessor.create2DFFT(ccTr1, false, false);
        Complex1D[] ccTrTr2 = imageProcessor.create2DFFT2(ccTr2, false, false);
        for (int row = 0; row < ccTr1.length; ++row) {
            for (int col = 0; col < ccTr1[0].length; ++col) {
                double r1 = ccTrTr1[row][col].re();
                double im1 = ccTrTr1[row][col].im();
                
                double r2 = ccTrTr2[row].x[col];
                double im2 = ccTrTr2[row].y[col];
                
                double r0 = a[row][col].re();
                
                System.out.println(String.format("*** %f %f  %f", 
                    (float)r1/norm, (float)r2/norm, (float)r0));
                
                assertTrue(Math.abs(r2 - r1) < 0.01);
                assertTrue(Math.abs(im2 - im1) < 0.01);
                
                
                //assertTrue(Math.abs(r0 - r1) < 0.01);
            }
        }
    }
    
    /*public void test0() throws Exception {
        
        String fileName = "lab.gif";
        String filePath = ResourceFinder.findFileInTestResources(fileName);     
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
   
        ImageProcessor imageProcessor = new ImageProcessor();
        Complex[][] cc = imageProcessor.convertImage(img);
        Complex[][] ccOut = imageProcessor.create2DFFT(cc, true, true);
        Complex[][] ccOut2 = imageProcessor.create2DFFT(ccOut, true, false);
           
        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                double vTr = ccOut[i][j].re();
                double vTrTr = ccOut2[i][j].re();
                int v = img.getValue(i, j);
                double diff = v - vTrTr;
                double r = vTrTr/v;
                System.out.println("img=" + v + " FFT(v)=" + vTr 
                    + " iFFT(FFT(v)=" + vTrTr + " diff = " + diff
                    + " ratio=" + r  + " nPix=" + img.getNPixels());
                //assertTrue(Math.abs(diff) < 1);
            }
        }
    }*/
    
    public void est() throws Exception {
        
        String fileName = "lab.gif";
        String filePath = ResourceFinder.findFileInTestResources(fileName);     
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
   
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Complex[][] ccOut = imageProcessor.create2DFFTWithSwapMajor(img, true);
        ccOut = imageProcessor.create2DFFT(ccOut, false);
   
        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                double vTrTr = ccOut[j][i].re();
                int v = img.getValue(i, j);
                assertTrue(Math.abs(v - vTrTr) < 1);
            }
        }
    }
    
}
