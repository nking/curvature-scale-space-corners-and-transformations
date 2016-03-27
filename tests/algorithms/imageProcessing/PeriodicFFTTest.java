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
    
    public void testFFT2D() throws Exception {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        String filePath = ResourceFinder.findFileInTestResources("lena.jpg");        
        GreyscaleImage img2 = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
        
        ImageDisplayer.displayImage("lena", img2);
        
        Complex[][] cc = imageProcessor.convertImage(img2);
                
        Complex[][] ccOut = imageProcessor.create2DFFT(cc, true);
        
        GreyscaleImage img3 = img2.createFullRangeIntWithDimensions();
        for (int col = 0; col < img3.getWidth(); col++) {
            for (int row = 0; row < img3.getHeight(); row++) {
                double re = ccOut[col][row].re();
                img3.setValue(col, row, (int)re);
            }
        }
                           
        ImageDisplayer.displayImage("FFT of lena", img3);
        
        ccOut = imageProcessor.create2DFFT(ccOut, false);
        for (int col = 0; col < img3.getWidth(); col++) {
            for (int row = 0; row < img3.getHeight(); row++) {
                double re = ccOut[col][row].re();
                img3.setValue(col, row, (int)re);
            }
        }
        
        ImageDisplayer.displayImage("FFT^-1 of FFT lena", img3);
        
        int maxDiff = Integer.MIN_VALUE;
        for (int i = 0; i < img2.getNPixels(); ++i) {
            int v = Math.abs(img2.getValue(i) - img3.getValue(i));
            if (v > maxDiff) {
                maxDiff = v;
            }
        }
        assertTrue(maxDiff < 5);  
        
        // --- the periodic fft -----
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
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
            
            GreyscaleImage imgT = img2.createFullRangeIntWithDimensions();
            imageProcessor.writeToImageWithSwapMajor(imgT, t);
            ImageDisplayer.displayImage(label, imgT);
            
            int nRows = t.length;
            int nCols = t[0].length;
            
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
            GreyscaleImage img2_ = new GreyscaleImage(img.getWidth(), img.getHeight(),
                Type.Bits32FullRangeInt);
            for (int row = 0; row < nRows; ++row) {
                for (int col = 0; col < nCols; ++col) {
                    int idx = (row * nCols) + col;
                    img2_.setValue(col, row, scaled[idx]);
                }
            }
            MiscDebug.writeImage(img2_, label + "_log_lena");

            GreyscaleImage img3_ = new GreyscaleImage(img.getWidth(), img.getHeight(),
                Type.Bits32FullRangeInt);
            for (int row = 0; row < nRows; ++row) {
                for (int col = 0; col < nCols; ++col) {
                    double re = t[row][col].re();
                    img3_.setValue(col, row, (int) re);
                }
            }
            HistogramEqualization hEq = new HistogramEqualization(img3_);
            hEq.applyFilter();

            MiscDebug.writeImage(img3_, label + "_lena");
        }
        
    }
    
    public void testPerfft2() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources("checkerboard_01.jpg");  
        GreyscaleImage img0 = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
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
