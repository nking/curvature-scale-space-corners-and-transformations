package algorithms.imageProcessing;

import algorithms.imageProcessing.GreyscaleImage.Type;
import algorithms.misc.Complex;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PeriodicFFTTest extends TestCase {

    public PeriodicFFTTest(String testName) {
        super(testName);
    }
    
    public void estFFT2D() throws Exception {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        String fileName = "lab.gif";
        //String fileName = "lena.jpg";
        //String fileName = "tmp.png";
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
        
        int nCols = img.getWidth();
        int nRows = img.getHeight();
        
        double maxValue = Double.MIN_VALUE;
        for (int row = 0; row < nRows; row++) {
            for (int col = 0; col < nCols; col++) {
                int re = img.getValue(col, row);
                if (re > maxValue) {
                    maxValue = re;
                }
            }
        }
        
        ImageDisplayer.displayImage(fileName, img);
                        
        Complex[][] ccOut = imageProcessor.create2DFFTWithSwapMajor(img, true);
        maxValue = Double.MIN_VALUE;
        for (int row = 0; row < nRows; row++) {
            for (int col = 0; col < nCols; col++) {
                int re = img.getValue(col, row);
                if (re > maxValue) {
                    maxValue = re;
                }
            }
        }
        
        GreyscaleImage img1 = img.createFullRangeIntWithDimensions();
        for (int row = 0; row < nRows; row++) {
            for (int col = 0; col < nCols; col++) {
                int re = (int)ccOut[row][col].re();
                img1.setValue(col, row, re);
            }
        }                           
        ImageDisplayer.displayImage("FFT of " + fileName, img1);

        // ---- inverse FFT -------
        ccOut = imageProcessor.create2DFFT(ccOut, false);
        maxValue = Double.MIN_VALUE;
        for (int row = 0; row < nRows; row++) {
            for (int col = 0; col < nCols; col++) {
                int re = (int)ccOut[row][col].re();
                if (re > maxValue) {
                    maxValue = re;
                }
            }
        }
        System.out.println("maxValue=" + maxValue);
        for (int row = 0; row < nRows; row++) {
            for (int col = 0; col < nCols; col++) {
                int re = (int)ccOut[row][col].re();
                img1.setValue(col, row, re);
            }
        }
        ImageDisplayer.displayImage("FFT^-1 of FFT " + fileName, img1);
        
        int maxDiff = Integer.MIN_VALUE;
        for (int i = 0; i < img1.getNPixels(); ++i) {
            int v = Math.abs(img.getValue(i) - img1.getValue(i));
            if (v > maxDiff) {
                maxDiff = v;
            }
        }
        System.out.println("maxDiff=" + maxDiff);
        assertTrue(maxDiff < 5);  
        
        // --- the periodic fft -----
        img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
        PeriodicFFT pfft = new PeriodicFFT();
        
        // results use notation a[row][col]
        // fftSmooth, fftPeriodic, smooth, periodic
        Complex[][][] results = pfft.perfft2(img, true);

        // comparison by eye to paper's Figure 3:
        // http://www.mi.parisdescartes.fr/~moisan/papers/2009-11r.pdf
        
        for (int i = 0; i < 4; ++i) {
            String label;
            Complex[][] t = results[i];
            if (i == 0) {
                label = "fftSmooth";
            } else if (i == 1) {
                label = "fftPeriodic";
            } else if (i == 2) {
                label = "smooth";
            } else {
                label = "periodic";
            }
            
            GreyscaleImage imgT = img.createFullRangeIntWithDimensions();
            imageProcessor.writeToImageWithSwapMajor(imgT, t);
                       
            ImageDisplayer.displayImage(label, imgT);
            
            double min = Double.MAX_VALUE;
            for (int row = 0; row < nRows; ++row) {
                for (int col = 0; col < nCols; ++col) {
                    double re = t[row][col].re();
                    if (re < min) {
                        min = re;
                    }
                }
            }
            if (min <= 0) {
                min += 0.001;
            } else {
                min = 0;
            }
            double[] logValues = new double[t.length * t[0].length];
            for (int row = 0; row < nRows; ++row) {
                for (int col = 0; col < nCols; ++col) {
                    int idx = (row * nCols) + col;
                    double re = t[row][col].re();
                    logValues[idx] = Math.log(re - min);
                }
            }
            
            int[] scaled = MiscMath.rescale(logValues, 0, 255);
            GreyscaleImage img2_ = img.createFullRangeIntWithDimensions();
            for (int row = 0; row < nRows; ++row) {
                for (int col = 0; col < nCols; ++col) {
                    int idx = (row * nCols) + col;
                    img2_.setValue(col, row, scaled[idx]);
                }
            }
            MiscDebug.writeImage(img2_, label + "_log_" + fileName);

            GreyscaleImage img3_ = img.createFullRangeIntWithDimensions();
            for (int row = 0; row < nRows; ++row) {
                for (int col = 0; col < nCols; ++col) {
                    double re = t[row][col].re();
                    img3_.setValue(col, row, (int) re);
                }
            }
            HistogramEqualization hEq = new HistogramEqualization(img3_);
            hEq.applyFilter();

            MiscDebug.writeImage(img3_, label + "_" + fileName);
        }
        int z = 1;
    }
    
    public void test0() throws Exception {
                
        //String fileName = "lena.jpg";
        //String fileName = "tmp.png";
        String fileName = "lab.gif";
        //String fileName = "house.gif";
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
        
        PeriodicFFT pfft = new PeriodicFFT();
        Complex[][][] pC = pfft.perfft2(img, true);
        Complex[][] periodicComponent = pC[0];
        
    }
    
    public void estPerfft2() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources("checkerboard_01.jpg");  
        GreyscaleImage img0 = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
        GreyscaleImage img = new GreyscaleImage(img0.getWidth(), img0.getHeight(),
            Type.Bits32FullRangeInt);
        for (int col = 0; col < img0.getWidth(); col++) {
            for (int row = 0; row < img0.getHeight(); row++) {
                int v = img0.getValue(col, row);
                img.setValue(col, row, v);
            }
        }
        
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.apply2DFFT(img, true);
        
        ImageDisplayer.displayImage("FFT of checkerboard", img);
        
        
        img = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
        PeriodicFFT pfft = new PeriodicFFT();
        // output of perfft2 uses notation a[row][col]
        Complex[][][] pC = pfft.perfft2(img, false);
        Complex[][] periodicComponent = pC[0];
        
        int nRows = pC.length;
        int nCols = pC[0].length;
        
        double[] logValues = new double[nRows * nCols];
        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {
                int idx = (row * nCols) + col;
                double re = periodicComponent[row][col].re();
                logValues[idx] = Math.log(re);
            }
        }
        
        int[] scaled = MiscMath.rescale(logValues, 0, 255);
        GreyscaleImage img2 = new GreyscaleImage(img.getWidth(), img.getHeight(),
            Type.Bits32FullRangeInt);
        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {
                int idx = (row * nCols) + col;
                img2.setValue(col, row, scaled[idx]);
            }
        }
        ImageDisplayer.displayImage("log of Periodic FFT of checkerboard", img2);
        
        GreyscaleImage img3 = new GreyscaleImage(img.getWidth(), img.getHeight(),
            Type.Bits32FullRangeInt);
        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {
                double re = periodicComponent[row][col].re();
                img3.setValue(col, row, (int)re);
            }
        }
        HistogramEqualization hEq = new HistogramEqualization(img3);
        hEq.applyFilter();
        
        ImageDisplayer.displayImage("Periodic FFT of checkerboard", img3);
        
    }
    
}
