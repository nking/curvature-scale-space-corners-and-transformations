package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.features.ORB;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.PairInt;
import algorithms.misc.Complex;
import algorithms.misc.ComplexModifiable;
import algorithms.misc.Histogram;
import algorithms.misc.MedianSmooth;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.StatsInSlidingWindow;
import algorithms.util.ResourceFinder;
import algorithms.util.TwoDFloatArray;
import algorithms.util.VeryLongBitString;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Logger;
import thirdparty.ca.uol.aig.fftpack.Complex1D;
import thirdparty.ca.uol.aig.fftpack.ComplexDoubleFFT;

/**
 *
 * @author nichole
 */
public class ImageProcessor {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    /**
     * use a sobel (first deriv gaussian sigma=0.5, binomial [-1, 0,1]
     * and return gradients in X and y. note the image may contain
     * negative values.
     * @param input
     * @return 
     */
    public GreyscaleImage[] createSobelGradients(GreyscaleImage input) {
    
        float[] kernel = Gaussian1DFirstDeriv.getBinomialKernelSigmaZeroPointFive();
        
        GreyscaleImage gX = input.copyToFullRangeIntImage();
        GreyscaleImage gY = input.copyToFullRangeIntImage();
        applyKernel1D(gX, kernel, true);
        applyKernel1D(gY, kernel, false);
        
        return new GreyscaleImage[]{gX, gY};
    }

    public void applySobelKernel(GreyscaleImage input) {

        float[] kernel = Gaussian1DFirstDeriv.getBinomialKernelSigmaZeroPointFive();
        
        GreyscaleImage gX = input.copyToFullRangeIntImage();
        GreyscaleImage gY = input.copyToFullRangeIntImage();
        applyKernel1D(gX, kernel, true);
        applyKernel1D(gY, kernel, false);
        
        GreyscaleImage img2 = combineConvolvedImages(gX, gY);

        input.resetTo(img2);
    }
    
    public void applySobelX(float[][] input) {

        float[] kernel = Gaussian1DFirstDeriv.getBinomialKernelSigmaZeroPointFive();
        
        applyKernel1D(input, kernel, true);
    }
    
    public void applySobelY(float[][] input) {

        float[] kernel = Gaussian1DFirstDeriv.getBinomialKernelSigmaZeroPointFive();
        
        applyKernel1D(input, kernel, false);
    }
    
    public Map<PairInt, Integer> applySobelKernel(GreyscaleImage input, 
        Set<PairInt> points) {
       
        float[] kernel = Gaussian1DFirstDeriv.getBinomialKernel(
            SIGMA.ZEROPOINTSEVENONE);

        return applyKernel(input, points, kernel);
    }
        
    /**
     * calculate the sobel gradient of the color image using CIELAB DeltaE 2000
     * and return gX, gY, and gXY with array indices being pixel
     * indexes of the image.
     * @param img
     * @return float[][]{gX, gY, gXY}
     */
    public float[][] calculateGradientUsingDeltaE2000(ImageExt img) {
        
        int n = img.getNPixels();
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        float jnd = 2.3f;
        
        // using 1D sobel kernel -1,0,1, calculating deltaE
        // between the pixels to either side of center pixel
        
        int x1, y1, x2, y2;
        
        float[] outX = new float[n];                    
        for (int i = 0; i < w; i++) {
            
            x1 = i - 1;
            if (x1 < 0) {
                x1 = 0;
            }
            x2 = i + 1;
            if (x2 > (w - 1)) {
                x2 = w - 1;
            }
            
            for (int j = 0; j < h; j++) {
                
                float[] lab1 = img.getCIELAB(x1, j);
                float[] lab2 = img.getCIELAB(x2, j);
                
                double deltaE = cieC.calcDeltaECIE2000(
                    lab1, lab2);
                
                outX[img.getInternalIndex(i, j)] = (float)deltaE;
            }                    
        }
        
        float[] outY = new float[n];                    
        for (int i = 0; i < w; i++) {
           
            for (int j = 0; j < h; j++) {
                
                y1 = j - 1;
                if (y1 < 0) {
                    y1 = 0;
                }
                y2 = j + 1;
                if (y2 > (h - 1)) {
                    y2 = h - 1;
                }
                
                float[] lab1 = img.getCIELAB(i, y1);
                float[] lab2 = img.getCIELAB(i, y2);
                
                double deltaE = cieC.calcDeltaECIE2000(
                    lab1, lab2);
                
                outY[img.getInternalIndex(i, j)] = (float)deltaE;
            }                    
        }
        
        // make a combined array
        float[] outXY = new float[outX.length];
        for (int i = 0; i < outX.length; ++i) {
            double gXY = Math.sqrt(outX[i] * outX[i] + outY[i] * outY[i]);
            outXY[i] = (float)gXY;
        }
      
        return new float[][]{outX, outY, outXY};
    }
    
    /**
     * apply a sobel kernel (gaussian first derivative, binomial approx 
     * for sigma=sqrt(2)/2) to the points in the region bounded by
     * (xLL, yLL) to (xUR, yUR), inclusive and return the results as a map.
     * @param pointValues values for points in the bounding region
     * @param xLL lower left corner x coordinate of bounding box.  xLL is less than xUR.
     * @param yLL lower left corner y coordinate of bounding box.  yLL is less than yUR.
     * @param xUR upper right corner x coordinate of bounding box.  xUR is larger than xLL.
     * @param yUR upper right corner y coordinate of bounding box.  yUR is larger than yLL.
     * @param outputGradientValues
     */
    public void applySobelKernel(Map<PairInt, Integer> pointValues,
        int xLL, int yLL, int xUR, int yUR,
        Map<PairInt, Integer> outputGradientValues) {
       
        float[] kernel = Gaussian1DFirstDeriv.getBinomialKernel(
            SIGMA.ZEROPOINTSEVENONE);

        applyKernel(pointValues, xLL, yLL, xUR, yUR, kernel, 
            outputGradientValues);
    }

    public void applyLaplacianKernel(GreyscaleImage input) {

        IKernel kernel = new Laplacian();
        Kernel kernelXY = kernel.getKernel();

        float norm = kernel.getNormalizationFactor();

        applyKernel(input, kernelXY, norm);
    }

    /**
     * apply the kernels to the input.  Note that the current image format
     * only accepts value between 0 and 255, inclusive.
     * @param input
     * @param kernelX
     * @param kernelY
     * @param normFactorX
     * @param normFactorY
     */
    protected void applyKernels(Image input, Kernel kernelX, Kernel kernelY,
        float normFactorX, float normFactorY) {

        /*
        assumes that kernelX is applied to a copy of the img
        and kernelY is applied to a separate copy of the img and
        then they are added in quadrature for the final result.
        */

        Image imgX = input.copyImage();

        Image imgY = input.copyImage();

        applyKernel(imgX, kernelX, normFactorX);

        applyKernel(imgY, kernelY, normFactorY);

        Image img2 = combineConvolvedImages(imgX, imgY);

        input.resetTo(img2);
    }

    protected void applyKernels(GreyscaleImage input, Kernel kernelX, Kernel kernelY,
        float normFactorX, float normFactorY) {

        GreyscaleImage[] gXgY = convolveWithKernels(input, kernelX, kernelY, normFactorX,
            normFactorY);

        GreyscaleImage img2 = combineConvolvedImages(gXgY[0], gXgY[1]);

        input.resetTo(img2);
    }
    
    protected GreyscaleImage[] convolveWithKernels(GreyscaleImage input, Kernel kernelX, Kernel kernelY,
        float normFactorX, float normFactorY) {

        /*
        assumes that kernelX is applied to a copy of the img
        and kernelY is applied to a separate copy of the img and
        then they are added in quadrature for the final result.
        */

        GreyscaleImage imgX = input.copyImage();

        GreyscaleImage imgY = input.copyImage();

        applyKernel(imgX, kernelX, normFactorX);

        applyKernel(imgY, kernelY, normFactorY);

        return new GreyscaleImage[]{imgX, imgY};
    }

    protected Map<PairInt, Integer> applyKernel(GreyscaleImage input, 
        Set<PairInt> points, float[] kernel) {

        /*
        assumes that kernelX is applied to a copy of the img
        and kernelY is applied to a separate copy of the img and
        then they are added in quadrature for the final result.
        */

        Map<PairInt, Integer> convX = applyKernel(input, points, kernel, true);

        Map<PairInt, Integer> convY = applyKernel(input, points, kernel, false);
        
        Map<PairInt, Integer> output = new HashMap<PairInt, Integer>();
        
        for (PairInt p : points) {
            
            int vX = convX.get(p).intValue();
            
            int vY = convY.get(p).intValue();
            
            int v = (int)Math.round(Math.sqrt(vX * vX + vY * vY));
            
            output.put(p, Integer.valueOf(v));
        }
        
        return output;
    }

    protected void applyKernel(Map<PairInt, Integer> pointValues, 
        int xLL, int yLL, int xUR, int yUR, float[] kernel,
        Map<PairInt, Integer> outputGradientValues) {

        /*
        assumes that kernelX is applied to a copy of the img
        and kernelY is applied to a separate copy of the img and
        then they are added in quadrature for the final result.
        */

        Map<PairInt, Integer> convX = applyKernel(pointValues, xLL, yLL, xUR, yUR,
            kernel, true);

        Map<PairInt, Integer> convY = applyKernel(pointValues, xLL, yLL, xUR, yUR,
            kernel, false);
                
        for (int xp = xLL; xp <= xUR; ++xp) {
            
            for (int yp = yLL; yp <= yUR; ++yp) {
                
                PairInt p = new PairInt(xp, yp);
                
                if (!convX.containsKey(p) || !convY.containsKey(p)) {
                    continue;
                }
                
                int vX = convX.get(p).intValue();

                int vY = convY.get(p).intValue();
                
                int v = (int) Math.round(Math.sqrt(vX * vX + vY * vY));

                if (v != 0) {
                    outputGradientValues.put(p, Integer.valueOf(v));
                }
            }
        }        
    }

    public Image combineConvolvedImages(Image imageX, Image imageY) {

        Image img2 = new Image(imageX.getWidth(), imageX.getHeight());

        for (int i = 0; i < imageX.getWidth(); i++) {
            for (int j = 0; j < imageX.getHeight(); j++) {

                int rX = imageX.getR(i, j);
                int gX = imageX.getG(i, j);
                int bX = imageX.getB(i, j);

                int rY = imageY.getR(i, j);
                int gY = imageY.getG(i, j);
                int bY = imageY.getB(i, j);

                double r = Math.sqrt(rX*rX + rY*rY);
                double g = Math.sqrt(gX*gX + gY*gY);
                double b = Math.sqrt(bX*bX + bY*bY);

                r = (r > 255) ? 255 : r;
                g = (g > 255) ? 255 : g;
                b = (b > 255) ? 255 : b;

                //int rgb = (int)(((rSum & 0x0ff) << 16)
                //    | ((gSum & 0x0ff) << 8) | (bSum & 0x0ff));

                img2.setRGB(i, j, (int)r, (int)g, (int)b);
            }
        }

        return img2;
    }

    /**
     * @param imageX
     * @param imageY
     * @return
     */
    public GreyscaleImage combineConvolvedImages(final GreyscaleImage imageX,
        final GreyscaleImage imageY) {

        GreyscaleImage img2 = imageX.createWithDimensions();

        for (int i = 0; i < imageX.getWidth(); i++) {
            for (int j = 0; j < imageX.getHeight(); j++) {

                int gX = imageX.getValue(i, j);

                int gY = imageY.getValue(i, j);

                double g = Math.sqrt(gX*gX + gY*gY);

                if (g > 255) {
                    g = 255;
                }
                
                img2.setValue(i, j, (int)g);
            }
        }

        return img2;
    }

    public void capToRange(GreyscaleImage img, int minV, int maxV) {
        int w = img.getWidth();
        int h = img.getHeight();
        for (int col = 0; col < w; ++col) {
            for (int row = 0; row < h; ++row) {
                int v = img.getValue(col, row);
                if (v < minV) {
                    img.setValue(col, row, minV);
                } else if (v > maxV) {
                    img.setValue(col, row, maxV);
                }
            }
        }
    }
    
    public int[][] applyKernel(int[][] a, int[][] b) {
        
        final int na0 = a.length;
        final int na1 = a[0].length;
        final int nb0 = b.length;
        final int nb1 = b[0].length;
        
        if (nb1 > nb0) {
            throw new IllegalArgumentException("expecting first argument to be " 
            + " data and the 2nd to be the kernel to convolve the data with");
        }

        int nbmid0 = (nb0 - 1)/2;
        int nbmid1 = (nb1 - 1)/2;
        
        int[][] output = new int[na0][na1];
        
        for (int i = 0; i < na0; i++) {

            output[i] = new int[na1];

            for (int j = 0; j < na1; j++) {

                int sum = 0;

                for (int col = 0; col < nb0; col++) {

                    int x = col - nbmid0;

                    int imgX = i + x;

                    // edge corrections.  use replication
                    if (imgX < 0) {
                        imgX = -1 * imgX - 1;
                    } else if (imgX >= na0) {
                        int diff = imgX - na0;
                        imgX = na0 - diff - 1;
                    }

                    for (int row = 0; row < nb1; row++) {

                        int y = row - nbmid1;

                        int imgY = j + y;

                        // edge corrections.  use replication
                        if (imgY < 0) {
                            imgY = -1 * imgY - 1;
                        } else if (imgY >= na1) {
                            int diff = imgY - na1;
                            imgY = na1 - diff - 1;
                        }

                        int k = b[col][row];

                        sum += k * a[imgX][imgY];
                    }
                }
if (sum > 511) {
    int z = 1;//  59, 125
}
                
                output[i][j] = sum;
            }
        }

        return output;
    }

    /**
     * apply kernel to input. NOTE, that because the image is composed of
     * vectors that should have values between 0 and 255, inclusive, if the
     * kernel application results in a value outside of that range, the value
     * is reset to 0 or 255.
     * @param input
     * @param kernel
     * @param normFactor
     */
    protected void applyKernel(Image input, Kernel kernel, float normFactor) {

        int h = (kernel.getWidth() - 1) >> 1;

        Image output = new Image(input.getWidth(), input.getHeight());

        for (int i = 0; i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {

                long rValue = 0;
                long gValue = 0;
                long bValue = 0;

                // apply the kernel to pixels centered in (i, j)

                for (int col = 0; col < kernel.getWidth(); col++) {

                    int x = col - h;

                    int imgX = i + x;

                    // edge corrections.  use replication
                    if (imgX < 0) {
                        imgX = -1 * imgX - 1;
                    } else if (imgX >= input.getWidth()) {
                        int diff = imgX - input.getWidth();
                        imgX = input.getWidth() - diff - 1;
                    }

                    for (int row = 0; row < kernel.getHeight(); row++) {

                        int y = row - h;

                        int imgY = j + y;

                        // edge corrections.  use replication
                        if (imgY < 0) {
                            imgY = -1 * imgY - 1;
                        } else if (imgY >= input.getHeight()) {
                            int diff = imgY - input.getHeight();
                            imgY = input.getHeight() - diff - 1;
                        }

                        int rPixel = input.getR(imgX, imgY);
                        int gPixel = input.getG(imgX, imgY);
                        int bPixel = input.getB(imgX, imgY);

                        int k = kernel.getValue(col, row);

                        rValue += k * rPixel;
                        gValue += k * gPixel;
                        bValue += k * bPixel;
                    }
                }

                rValue *= normFactor;
                gValue *= normFactor;
                bValue *= normFactor;

                /*
                if ((rValue > 255) || (rValue < 0)) {
                    throw new IllegalStateException("rValue is " + rValue);
                }
                if ((gValue > 255) || (gValue < 0)) {
                    throw new IllegalStateException("gValue is " + gValue);
                }
                if ((bValue > 255) || (bValue < 0)) {
                    throw new IllegalStateException("bValue is " + bValue);
                }*/

                if (rValue < 0) {
                    rValue = 0;
                }
                if (rValue > 255) {
                    rValue = 255;
                }
                if (gValue < 0) {
                    gValue = 0;
                }
                if (gValue > 255) {
                    gValue = 255;
                }
                if (bValue < 0) {
                    bValue = 0;
                }
                if (bValue > 255) {
                    bValue = 255;
                }

                output.setRGB(i, j, (int)rValue, (int)gValue, (int)bValue);
            }
        }

        input.resetTo(output);
    }
    
    /**
     * apply kernel to input. NOTE, that because the image is composed of
     * vectors that should have values between 0 and 255, inclusive, if the
     * kernel application results in a value outside of that range, the value
     * is reset to 0 or 255.
     * @param input
     * @param kernel
     * @param normFactor
     */
    protected void applyKernel(float[][] input, Kernel kernel, float normFactor) {

        int h = (kernel.getWidth() - 1) >> 1;

        int width = input.length;
        int height = input[0].length;
        
        float[][] output = new float[width][];
        for (int i = 0; i < width; ++i) {
            output[i] = new float[height];
        }
        
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {

                double value = 0;

                // apply the kernel to pixels centered in (i, j)

                for (int col = 0; col < kernel.getWidth(); col++) {

                    int x = col - h;

                    int imgX = i + x;

                    // edge corrections.  use replication
                    if (imgX < 0) {
                        imgX = -1 * imgX - 1;
                    } else if (imgX >= width) {
                        int diff = imgX - width;
                        imgX = width - diff - 1;
                    }

                    for (int row = 0; row < kernel.getHeight(); row++) {

                        int y = row - h;

                        int imgY = j + y;

                        // edge corrections.  use replication
                        if (imgY < 0) {
                            imgY = -1 * imgY - 1;
                        } else if (imgY >= height) {
                            int diff = imgY - height;
                            imgY = height - diff - 1;
                        }

                        float v = input[imgX][imgY];

                        int k = kernel.getValue(col, row);

                        value += k * v;
                    }
                }

                value *= normFactor;

                output[i][j] = (float)value;
            }
        }
        
        for (int i = 0; i < width; ++i) {
            System.arraycopy(output[i], 0, input[i], 0, height);
        }
    }
    
    /**
     * apply kernel to input for pixels in points.
     * @param input
     * @param points
     * @param kernel
     * @param calcForX
     * @return 
     */
    protected Map<PairInt, Integer> applyKernel(GreyscaleImage input, 
        Set<PairInt> points, float[] kernel, boolean calcForX) {

        int h = (kernel.length - 1) >> 1;

        Map<PairInt, Integer> output = new HashMap<PairInt, Integer>();

        for (PairInt p : points) {
            
            int i = p.getX();
            int j = p.getY();

            float sum = 0;

            // apply the kernel to pixels centered in (i, j)
            for (int g = 0; g < kernel.length; g++) {
                float gg = kernel[g];
                if (gg == 0) {
                    continue;
                }
                
                int x2, y2;
                if (calcForX) {
                    int delta = g - h;
                    x2 = i + delta;
                    y2 = j;
                    // edge corrections.  use replication
                    if (x2 < 0) {
                        x2 = -1 * x2 - 1;
                    } else if (x2 >= input.getWidth()) {
                        int diff = x2 - input.getWidth();
                        x2 = input.getWidth() - diff - 1;
                    }
                } else {
                    int delta = g - h;
                    y2 = j + delta;
                    x2 = i;
                    // edge corrections.  use replication
                    if (y2 < 0) {
                        y2 = -1 * y2 - 1;
                    } else if (y2 >= input.getHeight()) {
                        int diff = y2 - input.getHeight();
                        y2 = input.getHeight() - diff - 1;
                    }
                }
                
                int v = input.getValue(x2, y2);

                sum += gg * v;
            }
            
            output.put(p, Integer.valueOf((int) sum));
        }

        return output;
    }

    protected Map<PairInt, Integer> applyKernel(
        Map<PairInt, Integer> pointValues, int xLL, int yLL, int xUR, int yUR,
        float[] kernel, boolean calcForX) {

        int h = (kernel.length - 1) >> 1;

        Map<PairInt, Integer> output = new HashMap<PairInt, Integer>();

        for (int xp = xLL; xp <= xUR; ++xp) {
            
            for (int yp = yLL; yp <= yUR; ++yp) {
                
                PairInt p = new PairInt(xp, yp);
                
                float sum = 0;

                // apply the kernel to pixels centered in (i, j)
                for (int g = 0; g < kernel.length; g++) {
                    float gg = kernel[g];
                    if (gg == 0) {
                        continue;
                    }

                    int x2, y2;
                    if (calcForX) {
                        int delta = g - h;
                        x2 = xp + delta;
                        y2 = yp;
                        // edge corrections.  use replication
                        if (x2 < 0) {
                            x2 = -1 * x2 - 1;
                        } else if (x2 > xUR) {
                            if (!pointValues.containsKey(new PairInt(x2, y2))) {
                                x2 = xUR;
                            }
                        }
                    } else {
                        int delta = g - h;
                        y2 = yp + delta;
                        x2 = xp;
                        // edge corrections.  use replication
                        if (y2 < 0) {
                            y2 = -1 * y2 - 1;
                        } else if (y2 > yUR) {
                            if (!pointValues.containsKey(new PairInt(x2, y2))) {
                                y2 = yUR;
                            }
                        }
                    }
                    
                    // TODO: revisit this for normalization
                    PairInt p2 = new PairInt(x2, y2);
                    if (!pointValues.containsKey(p2)) {
                        continue;
                    }
                    
                    int v = pointValues.get(p2).intValue();

                    sum += gg * v;
                    
                } // end sum over kernl for a pixel
                
                int v = Math.round(sum);
                if (v != 0) {
                    output.put(p, Integer.valueOf(v));
                }
            }
        }

        return output;
    }

    /**
     * apply kernel to input. NOTE, that because the image is composed of vectors
     * that should have values between 0 and 255, inclusive, if the kernel application
     * results in a value outside of that range, the value is reset to 0 or
     * 255.
     * @param input
     * @param kernel
     * @param normFactor
     */
    protected void applyKernel(GreyscaleImage input, Kernel kernel, float normFactor) {

        int h = (kernel.getWidth() - 1) >> 1;

        GreyscaleImage output = input.createWithDimensions();

        //TODO: consider changing normalization to be similar to Kernel1DHelper

        for (int i = 0; i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {

                double value = 0;

                // apply the kernel to pixels centered in (i, j)

                for (int col = 0; col < kernel.getWidth(); col++) {

                    int x = col - h;

                    int imgX = i + x;

                    // edge corrections.  use replication
                    if (imgX < 0) {
                        imgX = -1 * imgX - 1;
                    } else if (imgX >= input.getWidth()) {
                        int diff = imgX - input.getWidth();
                        imgX = input.getWidth() - diff - 1;
                    }

                    for (int row = 0; row < kernel.getHeight(); row++) {

                        int y = row - h;

                        int imgY = j + y;

                        // edge corrections.  use replication
                        if (imgY < 0) {
                            imgY = -1 * imgY - 1;
                        } else if (imgY >= input.getHeight()) {
                            int diff = imgY - input.getHeight();
                            imgY = input.getHeight() - diff - 1;
                        }

                        int pixel = input.getValue(imgX, imgY);
                        int k = kernel.getValue(col, row);

                        value += k * pixel;
                    }
                }
                value *= normFactor;
                                
                int v = (int)value;

                if (v < 0) {
                    v = 0;
                }
                if (v > 255) {
                    v = 255;
                }
                output.setValue(i, j, v);
            }
        }

        input.resetTo(output);
    }

    /**
     * calculate theta from the gradient x and y images.
     * <pre>
     * The results are given as quadrants of values from 0 to 90.
     * 
     *     -45    90    45          y/x
                -  |  +
            0 -----|----- 0
                +  |  -
            45    90    -45

           when X is 0: if Y > 0, theta is 90
           when Y is 0: if X >= 0, theta is 0
     * </pre>
     * @param convolvedX
     * @param convolvedY
     * @return 
     */
    public GreyscaleImage computeTheta(final GreyscaleImage convolvedX,
        final GreyscaleImage convolvedY) {

        GreyscaleImage output = convolvedX.createSignedWithDimensions();

        for (int i = 0; i < convolvedX.getWidth(); i++) {
            for (int j = 0; j < convolvedX.getHeight(); j++) {

                double gX = convolvedX.getValue(i, j);

                double gY = convolvedY.getValue(i, j);

                int thetaG = calculateTheta(gX, gY);

                output.setValue(i, j, thetaG);

            }
        }

        return output;
    }

    /**
     * calculate theta from the gradient x and y images and transform to
     * range 0 to 180.
     * <pre>
     * 
     *           90    45          y/x
                -  |  +
          180 -----|----- 0
                +  |  -
                         
     * </pre>
     * 
     * @param convolvedX
     * @param convolvedY
     * @return 
     */
    public GreyscaleImage computeTheta180(final GreyscaleImage convolvedX,
        final GreyscaleImage convolvedY) {

        GreyscaleImage output = convolvedX.createWithDimensions();

        for (int i = 0; i < convolvedX.getWidth(); i++) {
            for (int j = 0; j < convolvedX.getHeight(); j++) {

                double gX = convolvedX.getValue(i, j);

                double gY = convolvedY.getValue(i, j);
                
                if (gY < 0) {
                    gX *= -1;
                    gY *= -1;
                }

                double radians = Math.atan2(gY, gX);
             
                int theta = (int)(radians * 180./Math.PI);
                if (theta == 180) {
                    theta = 0;
                }

                output.setValue(i, j, theta);

            }
        }

        return output;
    }

    /**
     * calculate theta, transforming values from -pi to pi to range 0 to 360.
     * Note that the value of theta for pixel with gradientX or gradientY
     * smaller than the resolution of the data (FWHM of PSF) is actually
     * orthogonal to it's real value, so those need to be interpreted
     * differently.
     * To calculate the minimum resolvable angle, combine the sigmas in
     * quadrature for all stages used to create the gradients.
     * For example, a blur using sigma=1 followed by a gradient using
     * sigma sqrt(1)/2 results in a sigma of
     * sqrt( (1*1) + (sqrt(1)/2)*(sqrt(1)/2) ).  The FWHM of the combined
     * operations is then approx 2.35 * sigma, so that's 3 pixels.
     * The minimum resolvable angle is then math.atan2(1, 3)*180./Math.PI is
     * 18.4 degrees for this example, so any theta within 19 degrees of
     * horizontal or vertical where there is signal in the image, has to be
     * corrected by 90 degrees.
     * Such correction isn't made here to allow the method to be used without
     * such knowledge.
     * <pre>
     *           Y
     *          90
     *     135   |    +45
     *           |
     *   180------------ 0   X
     *           |
     *    225    |   315
     *          270
     * </pre>
     * @param gradientX
     * @param gradientY
     * @return
     */
    public GreyscaleImage computeTheta360_0(final GreyscaleImage gradientX,
        final GreyscaleImage gradientY) {

        GreyscaleImage output = gradientX.createFullRangeIntWithDimensions();

        for (int i = 0; i < gradientX.getWidth(); i++) {
            for (int j = 0; j < gradientX.getHeight(); j++) {

                double gX = gradientX.getValue(i, j);

                double gY = gradientY.getValue(i, j);

                // -pi to pi radians
                double theta = Math.atan2(gY, gX);

                // transform to 0 to 2*pi radians
                if (theta < 0) {
                    theta += 2. * Math.PI;
                }

                int d = (int)Math.round(theta * 180./Math.PI);

                output.setValue(i, j, d);
            }
        }

        return output;
    }

    /**
     * compute theta as a polar angle that increases in a clockwise manner
     * and has a range from 0 to 359, inclusive.
     *
     * @param convolvedX
     * @param convolvedY
     * @return
     */
    public GreyscaleImage computeTheta360(final GreyscaleImage convolvedX,
        final GreyscaleImage convolvedY) {

        GreyscaleImage output = convolvedX.createFullRangeIntWithDimensions();

        for (int i = 0; i < convolvedX.getWidth(); i++) {
            for (int j = 0; j < convolvedX.getHeight(); j++) {

                double gX = convolvedX.getValue(i, j);

                double gY = convolvedY.getValue(i, j);

                double thetaR = (2. * Math.PI) - AngleUtil.polarAngleCCW(gX, gY);

                int thetaG = (int)Math.round(thetaR * 180./Math.PI);

                if (thetaG > 359) {
                    thetaG = thetaG - 360;
                }

                output.setValue(i, j, thetaG);

            }
        }

        return output;
    }

    public GreyscaleImage subtractImages(final GreyscaleImage image,
        final GreyscaleImage subtrImage) {

        if (image.getWidth() != subtrImage.getWidth()) {
            throw new IllegalArgumentException("image widths must be the same");
        }
        if (image.getHeight() != subtrImage.getHeight()) {
            throw new IllegalArgumentException("image heights must be the same");
        }

        GreyscaleImage output = image.createFullRangeIntWithDimensions();

        for (int i = 0; i < image.getWidth(); i++) {
            for (int j = 0; j < image.getHeight(); j++) {

                int diff = image.getValue(i, j) - subtrImage.getValue(i, j);

                output.setValue(i, j, diff);
            }
        }

        return output;
    }

    protected int calculateTheta(double gradientX, double gradientY) {

        /*  -45    90    45          y/x
                -  |  +
            0 -----|----- 0
                +  |  -
            45    90    -45

           when X is 0: if Y > 0, theta is 90
           when Y is 0: if X >= 0, theta is 0
        */

        if (gradientX == 0 && (gradientY != 0)) {
            return 90;
        }

        if (gradientY == 0) {
            return 0;
        }

        double theta = Math.atan(gradientY/gradientX)*180./Math.PI;

        int angle = (int)theta;

        // +x, +y -> +
        // -x, +y -> -
        // -x, -y -> +
        // +x, -y -> -

        if (!(gradientX < 0) && !(gradientY < 0)) {
            if (angle < 0) {
                // make it positive if negative
                angle *= -1;
            }
        } else if ((gradientX < 0) && !(gradientY < 0)) {
            if (!(angle < 0)) {
                // make it negative if it's not
                angle *= -1;
            }
        } else if ((gradientX < 0) && (gradientY < 0)) {
            if (angle < 0) {
                // make it positive if negative
                angle *= -1;
            }
        } else if (!(gradientX < 0) && (gradientY < 0)) {
            if (!(angle < 0)) {
                // make it negative if it's not
                angle *= -1;
            }
        }

        return angle;
    }

    /**
     * images bounded by zero's have to be shrunk to the columns and rows
     * of the first non-zeroes in order to keep the lines that should be
     * attached to the image edges from eroding completely.
     * Note, this expects image has values only non-negative numbers.
     * @param input
     * @return
     */
    public int[] shrinkImageToFirstNonZeros(final GreyscaleImage input) {

        int w = input.getWidth();
        int h = input.getHeight();

        int xNZFirst = -1;
        int xNZLast = -1;
        int yNZFirst = -1;
        int yNZLast = -1;

        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                int pixValue = input.getValue(i, j);
                if (pixValue > 0) {
                    xNZFirst = i;
                    break;
                }
            }
            if (xNZFirst > -1) {
                break;
            }
        }
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                int pixValue = input.getValue(i, j);
                if (pixValue > 0) {
                    yNZFirst = j;
                    break;
                }
            }
            if (yNZFirst > -1) {
                break;
            }
        }

        for (int i = (w - 1); i > -1; i--) {
            for (int j = (h - 1); j > -1; j--) {
                if (input.getValue(i, j) > 0) {
                    xNZLast = i;
                    break;
                }
            }
            if (xNZLast > -1) {
                break;
            }
        }

        for (int j = (h - 1); j > -1; j--) {
            for (int i = (w - 1); i > -1; i--) {
                int pixValue = input.getValue(i, j);
                if (pixValue > 0) {
                    yNZLast = j;
                    break;
                }
            }
            if (yNZLast > -1) {
                break;
            }
        }

        if ((xNZFirst > 0) || (xNZLast < (w - 1)) || (yNZFirst > 0) ||
            (yNZLast < (h - 1))) {

            //add a 2 pix border
            xNZFirst -= 2;
            yNZFirst -= 2;
            if (xNZFirst < 0) {
                xNZFirst = 0;
            }
            if (yNZFirst < 0) {
                yNZFirst = 0;
            }
            if (xNZLast == -1) {
                xNZLast = input.getWidth() - 1;
            } else if (xNZLast < (input.getWidth() - 2)) {
                // add a 1 pix border
                xNZLast += 2;
            } else if (xNZLast < (input.getWidth() - 1)) {
                // add a 1 pix border
                xNZLast++;
            }
            if (yNZLast == -1) {
                yNZLast = input.getHeight() - 1;
            } else if (yNZLast < (input.getHeight() - 2)) {
                // add a 1 pix border
                yNZLast += 2;
            } else if (yNZLast < (input.getHeight() - 1)) {
                // add a 1 pix border
                yNZLast++;
            }

            int xLen = xNZLast - xNZFirst + 1;

            int yLen = yNZLast - yNZFirst + 1;

            GreyscaleImage output = new GreyscaleImage(xLen, yLen, input.getType());
            output.setXRelativeOffset(xNZFirst);
            output.setYRelativeOffset(yNZFirst);

            for (int i = xNZFirst; i <= xNZLast; i++) {

                int iIdx = i - xNZFirst;

                for (int j = yNZFirst; j <= yNZLast; j++) {

                    int jIdx = j - yNZFirst;

                    output.setValue(iIdx, jIdx, input.getValue(i, j));
                }
            }

            input.resetTo(output);

            return new int[]{xNZFirst, yNZFirst};
        }

        return new int[]{0, 0};
    }

    public void shrinkImage(final GreyscaleImage input,
        int[] offsetsAndDimensions) {

        int w2 = offsetsAndDimensions[2];
        int h2 = offsetsAndDimensions[3];

        int offset1X = offsetsAndDimensions[0];
        int offset1Y = offsetsAndDimensions[1];

        GreyscaleImage output = new GreyscaleImage(w2, h2, input.getType());
        output.setXRelativeOffset(offset1X);
        output.setYRelativeOffset(offset1Y);

        int x = 0;

        int endCol = (offset1X + w2);
        if (endCol > input.getWidth()) {
            endCol = input.getWidth();
        }
        int endRow =  (offset1Y + h2);
        if (endRow > input.getHeight()) {
            endRow = input.getHeight();
        }

        for (int col = offset1X; col < endCol; col++) {

            int y = 0;

            for (int row = offset1Y; row < endRow; row++) {

                int v = input.getValue(col, row);

                output.setValue(x, y, v);

                y++;
            }

            x++;
        }

        input.resetTo(output);
    }

    /**
     * change coordinates of the input as if they were cropped to the given
     * offset and dimensions.
     * @param input
     * @param offsetsAndDimensions int[]{xOffset, yOffset, finalWidth, finalHeight}
     */
    public void shrinkImage(final Set<PairInt> input,
        int[] offsetsAndDimensions) {

        /*
          -------        ____
          |     |  ==>  |    |
          |     |       |____|
          -------
        */
        //xOffset, yOffset, width, height
        // subtract xOffset from x in input and yOffset from y in input

        //TODO: remove points out of bounds of final image

        for (PairInt p : input) {
            p.setX(p.getX() - offsetsAndDimensions[0]);
            p.setY(p.getY() - offsetsAndDimensions[1]);
        }

    }

    /**
     * change coordinates of the input as if they were cropped to the given
     * offset and dimensions.
     * @param input
     * @param offsetsAndDimensions int[]{xOffset, yOffset, finalWidth, finalHeight}
     */
    public void shrinkImage(final PairIntArray input,
        int[] offsetsAndDimensions) {

        /*
          -------        ____
          |     |  ==>  |    |
          |     |       |____|
          -------
        */
        //xOffset, yOffset, width, height
        // subtract xOffset from x in input and yOffset from y in input

        //TODO: remove points out of bounds of final image

        for (int i = 0; i < input.getN(); ++i) {
            int x = input.getX(i) - offsetsAndDimensions[0];
            int y = input.getY(i) - offsetsAndDimensions[1];
            input.set(i, x, y);
        }

    }

    public void convertToBinaryImage(GreyscaleImage input) {
        for (int col = 0; col < input.getWidth(); col++) {
            for (int row = 0; row < input.getHeight(); row++) {
                int v = input.getValue(col, row);
                if (v != 0) {
                    input.setValue(col, row, 1);
                }
            }
        }
    }

    /**
     * NOT READY FOR USE YET.
     *
     * @param theta
     * @return
     * @throws java.io.IOException
     * @throws java.security.NoSuchAlgorithmException
     */
    public GreyscaleImage createRoughSkyMask(GreyscaleImage theta) throws
        IOException, NoSuchAlgorithmException {

        if (theta == null) {
            throw new IllegalArgumentException("theta cannot be null");
        }

        theta = theta.copyImage();

        ImageSegmentation imageSegmentation = new ImageSegmentation();
        imageSegmentation.applyUsingKMPP(theta, 2);

        subtractMinimum(theta);

        convertToBinaryImage(theta);

        removeSpurs(theta);

        throw new UnsupportedOperationException("not ready for use yet");
        //return theta;
    }

    public void subtractMinimum(GreyscaleImage input) {

        int min = input.getMin();

        for (int col = 0; col < input.getWidth(); col++) {

            for (int row = 0; row < input.getHeight(); row++) {

                int v = input.getValue(col, row);

                int f = v - min;

                input.setValue(col, row, f);
            }
        }
    }

    /**
     * multiply these images, that is pixel by pixel multiplication.
     * input2 is assumed to be 0 or 1
     * @param input1
     * @param input2 the mask of 0's and 1's to apply to input1
     */
    public void multiplyBinary(Image input1, GreyscaleImage input2)  {

        if (input1 == null) {
            throw new IllegalArgumentException("input1 cannot be null");
        }
        if (input2 == null) {
            throw new IllegalArgumentException("input2 cannot be null");
        }
        if (input1.getWidth() != input2.getWidth()) {
            throw new IllegalArgumentException(
            "input1 and input2 must have same widths");
        }
        if (input1.getHeight()!= input2.getHeight()) {
            throw new IllegalArgumentException(
            "input1 and input2 must have same heights");
        }

        for (int col = 0; col < input1.getWidth(); col++) {

            for (int row = 0; row < input1.getHeight(); row++) {

                int m = input2.getValue(col, row);

                if (m == 0) {

                    input1.setRGB(col, row, 0, 0, 0);
                }
            }
        }
    }

    /**
     * compare each pixel and set output to 0 if both inputs are 0, else set
     * output to 1.
     * @param input1
     * @param input2
     * @return
     */
    public GreyscaleImage binaryOr(GreyscaleImage input1, GreyscaleImage input2)  {

        if (input1 == null) {
            throw new IllegalArgumentException("input1 cannot be null");
        }
        if (input2 == null) {
            throw new IllegalArgumentException("input2 cannot be null");
        }
        if (input1.getWidth() != input2.getWidth()) {
            throw new IllegalArgumentException(
            "input1 and input2 must have same widths");
        }
        if (input1.getHeight()!= input2.getHeight()) {
            throw new IllegalArgumentException(
            "input1 and input2 must have same heights");
        }

        GreyscaleImage output = input1.createWithDimensions();

        for (int col = 0; col < input1.getWidth(); col++) {

            for (int row = 0; row < input1.getHeight(); row++) {

                int v1 = input1.getValue(col, row);

                int v2 = input2.getValue(col, row);

                if ((v1 != 0) || (v2 != 0)) {
                    output.setValue(col, row, 1);
                }
            }
        }

        return output;
    }

    protected void blur(GreyscaleImage input, float[] kernel) {

        applyKernel1D(input, kernel, true);

        applyKernel1D(input, kernel, false);
    }
    
    protected void applyKernelTwo1Ds(int[][] input, float[] kernel) {

        applyKernel1D(input, kernel, true);

        applyKernel1D(input, kernel, false);
    }
    
    public void applyKernelTwo1Ds(float[][] input, float[] kernel) {

        applyKernel1D(input, kernel, true);

        applyKernel1D(input, kernel, false);
        
    }

    /**
     * apply a sigma=0.5 first derivative of Gaussian ([0.5, 0, -0.5], a.k.a. Sobel)
     * @param input
     */
    public GreyscaleImage createSmallFirstDerivGaussian(GreyscaleImage input) {

        float[] kernel = Gaussian1DFirstDeriv.getKernel(
            SIGMA.ZEROPOINTFIVE);

        GreyscaleImage gX = input.copyImage();
        GreyscaleImage gY = input.copyImage();
        applyKernel1D(gX, kernel, true);
        applyKernel1D(gY, kernel, false);

        return combineConvolvedImages(gX, gY);
    }

    protected void blur(GreyscaleImage input, float[] kernel, int minValue, int maxValue) {

        applyKernel1D(input, kernel, true, minValue, maxValue);

        applyKernel1D(input, kernel, false, minValue, maxValue);
    }

    public void blur(GreyscaleImage input, float sigma) {

        float[] kernel = Gaussian1D.getKernel(sigma);

        blur(input, kernel);
    }
    
    /**
     * in order to make a smoother "blur" operation, the full gaussian
     * profile down to 0.001 * HWZI is created and convolved with image input
     * at a smaller sigma, recursively, until the resulting sigma equals 
     * the sigma given to the method.
     * gaussian profiles add in quadrature, 
     * sigma_tot^2 = sigma_1^2 + sigma_2^2.
     * 
     * @param input
     * @param sigma should be a quadrature factor of sqrt(2)/2, that is 
     *     sqrt(2)/2, 1, sqrt(1.5), sqrt(2.0), sqrt(2.5), sqrt(3.0), 
     *     sqrt(3.5), sqrt(4.0), sqrt(4.5)...
     *     in other words, the square of sigma should be 0.5 or a positive non zero 
     *     integer or a positive non zero integer plus 0.5.
     */
    public void recursiveBlur(float[][] input, SIGMA sigma) {
        
        float sigmaF = SIGMA.getValue(sigma);
        float sigmaFSQ = sigmaF * sigmaF;
        int sigmaISQ = (int)Math.floor(sigmaF);
        float remainder = sigmaFSQ - sigmaISQ;
        if (!(
            ((sigmaISQ >= 0) && (Math.abs(remainder - 0.5f) < 0.01)) 
            || ((sigmaISQ > 0) && (remainder < 0.01))
            )) {
            throw new IllegalArgumentException("sigma has to be a value that"
                + "is the result of recursive combinations of sqrt(2)/2");
        }
        
        float sigma0Sq = (float)Math.sqrt(2.)/2.f;
        sigma0Sq *= sigma0Sq;
        
        float currentSigmaSq = sigma0Sq;
        
        float finalSigmaSq = SIGMA.getValue(sigma);
        finalSigmaSq *= finalSigmaSq;
        
        //float[] kernel = Gaussian1D.getKernel((float)Math.sqrt(2.)/2.f, 0);
        //for sigma=sqrt(2)/2 kernel
        //   =[6.962646E-5, 0.010333488, 0.20755373, 0.5641896, 
        //     0.20755373, 0.010333488, 6.962646E-5]
        // which is essentially 1, 20, 56, 20, 1
        final float[] kernel = new float[]{0.010333488f, 0.20755373f, 0.5641896f, 
             0.20755373f, 0.010333488f};
        
        do {
            applyKernel1D(input, kernel, true);
            applyKernel1D(input, kernel, false);
            
            currentSigmaSq = currentSigmaSq + sigma0Sq;
            
        } while ((finalSigmaSq > currentSigmaSq) && Math.abs(finalSigmaSq - currentSigmaSq) > 0.01);
        
    }

    public void blur(GreyscaleImage input, SIGMA sigma) {

        float[] kernel = Gaussian1D.getKernel(sigma);

        blur(input, kernel);
    }

    public void blur(GreyscaleImage input, SIGMA sigma, int minValue, int maxValue) {
        float[] kernel = Gaussian1D.getKernel(sigma);
        blur(input, kernel, minValue, maxValue);
    }

    /**
     * blur the r, g, b vectors of image input by sigma.
     * @param input
     * @param sigma
     */
    public void blur(Image input, float sigma) {

        float[] kernel = Gaussian1D.getKernel(sigma);

        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();

        int w = input.getWidth();
        int h = input.getHeight();
        Image output = input.copyImage();

        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                double[] conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, true);
                output.setRGB(i, j, (int)conv[0], (int)conv[1], (int)conv[2]);
            }
        }

        input.resetTo(output);

        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                double[] conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, false);
                output.setRGB(i, j, (int)conv[0], (int)conv[1], (int)conv[2]);
            }
        }

        input.resetTo(output);
    }

    /**
     * blur the r, g, b vectors of image input by sigma.
     * @param input
     * @param sigma
     */
    public void blur(Image input, float sigma, int minValue, int maxValue) {

        float[] kernel = Gaussian1D.getKernel(sigma);

        blur(input, kernel, minValue, maxValue);
    }

    /**
     * blur the r, g, b vectors of image input by sigma.
     * @param input
     * @param sigma
     */
    public void blur(Image input, SIGMA sigma, int minValue, int maxValue) {

        float[] kernel = Gaussian1D.getKernel(sigma);

        blur(input, kernel, minValue, maxValue);
    }
    
    /**
     * blur the r, g, b vectors of image input by sigma.
     * @param input
     */
    public void blur(Image input, float[] kernel, int minValue, int maxValue) {

        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();

        int w = input.getWidth();
        int h = input.getHeight();
        Image output = input.copyImage();

        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                double[] conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, true);
                int r = (int)Math.round(conv[0]);
                int g = (int)Math.round(conv[1]);
                int b = (int)Math.round(conv[2]);
                if (r < minValue) {
                    r = minValue;
                } else if (r > maxValue) {
                    r = maxValue;
                }
                if (g < minValue) {
                    g = minValue;
                } else if (g > maxValue) {
                    g = maxValue;
                }
                if (b < minValue) {
                    b = minValue;
                } else if (b > maxValue) {
                    b = maxValue;
                }
                output.setRGB(i, j, r, g, b);
            }
        }

        input.resetTo(output);

        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                double[] conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, false);
                int r = (int)Math.round(conv[0]);
                int g = (int)Math.round(conv[1]);
                int b = (int)Math.round(conv[2]);
                if (r < minValue) {
                    r = minValue;
                } else if (r > maxValue) {
                    r = maxValue;
                }
                if (g < minValue) {
                    g = minValue;
                } else if (g > maxValue) {
                    g = maxValue;
                }
                if (b < minValue) {
                    b = minValue;
                } else if (b > maxValue) {
                    b = maxValue;
                }
                output.setRGB(i, j, r, g, b);
            }
        }

        input.resetTo(output);
    }

    public void applyFirstDerivGaussian(GreyscaleImage input, SIGMA sigma,
        int minValueRange, int maxValueRange) {

        float[] kernel = Gaussian1DFirstDeriv.getBinomialKernel(sigma);

        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();

        int w = input.getWidth();
        int h = input.getHeight();
        GreyscaleImage output = input.copyImage();

        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {

                double conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, true);

                int v = (int)Math.round(conv);

                if (v < minValueRange) {
                    v = minValueRange;
                } else if (v > maxValueRange) {
                    v = maxValueRange;
                }

                output.setValue(i, j, v);
            }
        }

        input.resetTo(output);

        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {

                double conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, false);

                int v = (int)Math.round(conv);

                if (v < minValueRange) {
                    v = minValueRange;
                } else if (v > maxValueRange) {
                    v = maxValueRange;
                }

                output.setValue(i, j, v);
            }
        }

        input.resetTo(output);
    }

    public void applySecondDerivGaussian(GreyscaleImage input, SIGMA sigma,
        int minValueRange, int maxValueRange) {

        float[] kernel = Gaussian1DSecondDeriv.getBinomialKernel(sigma);

        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();

        int w = input.getWidth();
        int h = input.getHeight();
        GreyscaleImage output = input.copyImage();

        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {

                double conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, true);

                int v = (int)Math.round(conv);

                if (v < minValueRange) {
                    v = minValueRange;
                } else if (v > maxValueRange) {
                    v = maxValueRange;
                }

                output.setValue(i, j, v);
            }
        }

        input.resetTo(output);

        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {

                double conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, false);

                int v = (int)Math.round(conv);

                if (v < minValueRange) {
                    v = minValueRange;
                } else if (v > maxValueRange) {
                    v = maxValueRange;
                }

                output.setValue(i, j, v);
            }
        }

        input.resetTo(output);
    }
    
    public int[] performSecondDerivGaussian(GreyscaleImage input, 
        SIGMA sigma) {

        float[] kernel = Gaussian1DSecondDeriv.getBinomialKernel(sigma);

        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();

        int w = input.getWidth();
        int h = input.getHeight();
        
        int[] output = new int[w * h];

        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {

                double conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, true);

                output[input.getInternalIndex(i, j)] = 
                    (int)Math.round(conv);
            }
        }
        
        int[] output2 = Arrays.copyOf(output, 
            output.length);

        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {

                double conv = kernel1DHelper.convolvePointWithKernel(
                    output2, w, h, i, j, kernel, false);

                output[input.getInternalIndex(i, j)] = 
                    (int)Math.round(conv);
            }
        }

        return output;
    }

    public void applyKernel1D(GreyscaleImage input, float[] kernel,
        boolean calcForX) {

        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();

        GreyscaleImage output = input.createWithDimensions();

        for (int i = 0; i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {
                double conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, calcForX);
                int g = (int)Math.round(conv);
                output.setValue(i, j, g);
            }
        }

        input.resetTo(output);
    }
    
    public void applyKernel1D(int[][] input, float[] kernel,
        boolean calcForX) {

        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();

        int w = input.length;
        int h = input[0].length;
        
        int[][] output = new int[w][];

        for (int i = 0; i < w; i++) {
            output[i] = new int[h];
            for (int j = 0; j < h; j++) {
                double conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, calcForX);
                int g = (int)Math.round(conv);
                output[i][j] = g;
            }
        }

        for (int i = 0; i < w; i++) {
            System.arraycopy(output[i], 0, input[i], 0, h);
        }
    }
    
    public void applyKernel1D(float[][] input, float[] kernel,
        boolean calcForX) {

        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();

        int w = input.length;
        int h = input[0].length;
        
        float[][] output = new float[w][];

        for (int i = 0; i < w; i++) {
            output[i] = new float[h];
            for (int j = 0; j < h; j++) {
                double conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, calcForX);
                float g = (float)conv;
                output[i][j] = g;
            }
        }

        for (int i = 0; i < w; i++) {
            System.arraycopy(output[i], 0, input[i], 0, h);
        }
    }

    public void applyKernel1D(GreyscaleImage input, float[] kernel,
        boolean calcForX, int minValue, int maxValue) {

        GreyscaleImage output = convolveWithKernel1D(input, kernel, calcForX, 
            minValue, maxValue);

        input.resetTo(output);
    }
    
    public GreyscaleImage convolveWithKernel1D(GreyscaleImage input, float[] kernel,
        boolean calcForX, int minValue, int maxValue) {

        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();

        GreyscaleImage output = input.createWithDimensions();

        for (int i = 0; i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {
                double conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, calcForX);
                int g = (int)Math.round(conv);
                if (g < minValue) {
                    g = minValue;
                } else if (g > maxValue) {
                    g = maxValue;
                }
                output.setValue(i, j, g);
            }
        }

        return output;
    }

    /**
     * blur the r, g, b vectors of image input by sigma.
     * @param input
     * @param kernel iD kernel
     */
    protected void blur(Image input, float[] kernel) {

        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();

        int w = input.getWidth();
        int h = input.getHeight();
        Image output = (ImageExt)input.copyImage();

        for (int i = 0; i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {
                double[] conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, true);
                output.setRGB(i, j, (int)conv[0], (int)conv[1], (int)conv[2]);
            }
        }

        input.resetTo(output);

        for (int i = 0; i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {
                double[] conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, false);
                output.setRGB(i, j, (int)conv[0], (int)conv[1], (int)conv[2]);
            }
        }

        input.resetTo(output);
    }

    public void divideByBlurredSelf(GreyscaleImage input, float sigma) {

        GreyscaleImage input2 = input.copyImage();

        blur(input, sigma);

        for (int i = 0; i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {
                int v = input.getValue(i, j);
                int vorig = input2.getValue(i, j);
                if (v != 0) {
                    float r = (float)vorig/(float)v;
                    input.setValue(i, j, (int)(100*r));
                }
            }
        }
    }

    /**
     * make a binary mask with the given zeroCoords as a group of starter points
     * for the mask and also set to '0' any points within zeroCoords' bounds.
     *
     * @param theta
     * @param zeroCoords
     * @return
     */
    public GreyscaleImage createMask(GreyscaleImage theta, PairIntArray zeroCoords) {

        GreyscaleImage out = theta.createWithDimensions();

        out.fill(1);

        for (int pIdx = 0; pIdx < zeroCoords.getN(); pIdx++) {

            int x = zeroCoords.getX(pIdx);
            int y = zeroCoords.getY(pIdx);
            out.setValue(x, y, 0);
        }

        return out;
    }

    /**
     * make a binary mask with the given zeroCoords as a group of starter points
     * for the mask and also set to '0' any points within zeroCoords' bounds.
     *
     * @param theta
     * @param nonzeroCoords
     * @return
     */
    public GreyscaleImage createInvMask(GreyscaleImage theta,
        PairIntArray nonzeroCoords) {

        GreyscaleImage out = theta.createWithDimensions();

        for (int pIdx = 0; pIdx < nonzeroCoords.getN(); pIdx++) {

            int x = nonzeroCoords.getX(pIdx);
            int y = nonzeroCoords.getY(pIdx);
            out.setValue(x, y, 1);
        }

        return out;
    }

    /**
     * this is meant to operate on an image with only 0's and 1's
     * @param input
     */
    public void removeSpurs(GreyscaleImage input) {

        int width = input.getWidth();
        int height = input.getHeight();

        int nIterMax = 1000;
        int nIter = 0;
        int numRemoved = 1;

        while ((nIter < nIterMax) && (numRemoved > 0)) {

            numRemoved = 0;

            for (int col = 0; col < input.getWidth(); col++) {

                if ((col < 2) || (col > (width - 3))) {
                    continue;
                }

                for (int row = 0; row < input.getHeight(); row++) {

                    if ((row < 2) || (row > (height - 3))) {
                       continue;
                    }

                    int v = input.getValue(col, row);

                    if (v == 0) {
                        continue;
                    }

                    // looking for pixels having only one neighbor who subsequently
                    // has only 1 or 2 neighbors
                    // as long as neither are connected to image boundaries

                    int neighborIdx = getIndexIfOnlyOneNeighbor(input, col, row);

                    if (neighborIdx > -1) {
                        int neighborX = input.getCol(neighborIdx);
                        int neighborY = input.getRow(neighborIdx);

                        int nn = count8RegionNeighbors(input, neighborX, neighborY);

                        if (nn <= 2) {
                            input.setValue(col, row, 0);
                            numRemoved++;
                        }
                    } else {
                        int n = count8RegionNeighbors(input, col, row);
                        if (n == 0) {
                            input.setValue(col, row, 0);
                            numRemoved++;
                        }
                    }
                }
            }

            log.fine("numRemoved=" + numRemoved + " nIter=" + nIter);

            nIter++;
        }

    }

    public void removeSpurs(Set<PairInt> points, int width, int height) {

        int nIterMax = 1000;
        int nIter = 0;
        int numRemoved = 1;

        while ((nIter < nIterMax) && (numRemoved > 0)) {

            numRemoved = 0;

            Set<PairInt> rm = new HashSet<PairInt>();

            for (PairInt p : points) {

                // looking for pixels having only one neighbor who subsequently
                // has only 1 or 2 neighbors
                // as long as neither are connected to image boundaries

                PairInt neighbor = getIndexIfOnlyOneNeighbor(points, p,
                    width, height);

                if (neighbor != null) {

                    int nn = count8RegionNeighbors(points, neighbor, width,
                        height);

                    if (nn <= 2) {
                        rm.add(p);
                        numRemoved++;
                    }
                } else {
                    int n = count8RegionNeighbors(points, p, width, height);
                    if (n == 0) {
                        rm.add(p);
                        numRemoved++;
                    }
                }
            }

            for (PairInt p : rm) {
                points.remove(p);
            }

            log.fine("numRemoved=" + numRemoved + " nIter=" + nIter);

            nIter++;
        }

    }

    protected int count8RegionNeighbors(GreyscaleImage input, int x, int y) {

        int width = input.getWidth();
        int height = input.getHeight();

        int count = 0;
        for (int c = (x - 1); c <= (x + 1); c++) {
            if ((c < 0) || (c > (width - 1))) {
                continue;
            }
            for (int r = (y - 1); r <= (y + 1); r++) {
                if ((r < 0) || (r > (height - 1))) {
                    continue;
                }
                if ((c == x) && (r == y)) {
                    continue;
                }
                int v = input.getValue(c, r);
                if (v > 0) {
                    count++;
                }
            }
        }

        return count;
    }

    protected int count8RegionNeighbors(GreyscaleImage input, int x, int y,
        int edgeValue) {

        int width = input.getWidth();
        int height = input.getHeight();

        int count = 0;
        for (int c = (x - 1); c <= (x + 1); c++) {
            if ((c < 0) || (c > (width - 1))) {
                continue;
            }
            for (int r = (y - 1); r <= (y + 1); r++) {
                if ((r < 0) || (r > (height - 1))) {
                    continue;
                }
                if ((c == x) && (r == y)) {
                    continue;
                }
                int v = input.getValue(c, r);
                if (v == edgeValue) {
                    count++;
                }
            }
        }

        return count;
    }

    protected int count8RegionNeighbors(Set<PairInt> points, PairInt point,
        int width, int height) {

        int x = point.getX();
        int y = point.getY();

        int count = 0;

        for (int c = (x - 1); c <= (x + 1); c++) {
            if ((c < 0) || (c > (width - 1))) {
                continue;
            }
            for (int r = (y - 1); r <= (y + 1); r++) {
                if ((r < 0) || (r > (height - 1))) {
                    continue;
                }
                if ((c == x) && (r == y)) {
                    continue;
                }
                PairInt tmp = new PairInt(c, r);
                if (points.contains(tmp)) {
                    count++;
                }
            }
        }

        return count;
    }

    protected int getIndexIfOnlyOneNeighbor(GreyscaleImage input, int x, int y) {

        int width = input.getWidth();
        int height = input.getHeight();

        int count = 0;
        int xNeighbor = -1;
        int yNeighbor = -1;

        for (int c = (x - 1); c <= (x + 1); c++) {
            if ((c < 0) || (c > (width - 1))) {
                continue;
            }
            for (int r = (y - 1); r <= (y + 1); r++) {
                if ((r < 0) || (r > (height - 1))) {
                    continue;
                }
                if ((c == x) && (r == y)) {
                    continue;
                }
                int v = input.getValue(c, r);
                if (v > 0) {
                    if (count > 0) {
                        return -1;
                    }
                    xNeighbor = c;
                    yNeighbor = r;
                    count++;
                }
            }
        }

        if (count == 0) {
            return -1;
        }

        int index = input.getIndex(xNeighbor, yNeighbor);

        return index;
    }

    protected int getIndexIfOnlyOneNeighbor(GreyscaleImage input, int x, int y,
        int edgeValue) {

        int width = input.getWidth();
        int height = input.getHeight();

        int count = 0;
        int xNeighbor = -1;
        int yNeighbor = -1;

        for (int c = (x - 1); c <= (x + 1); c++) {
            if ((c < 0) || (c > (width - 1))) {
                continue;
            }
            for (int r = (y - 1); r <= (y + 1); r++) {
                if ((r < 0) || (r > (height - 1))) {
                    continue;
                }
                if ((c == x) && (r == y)) {
                    continue;
                }
                int v = input.getValue(c, r);
                if (v == edgeValue) {
                    if (count > 0) {
                        return -1;
                    }
                    xNeighbor = c;
                    yNeighbor = r;
                    count++;
                }
            }
        }

        if (count == 0) {
            return -1;
        }

        int index = input.getIndex(xNeighbor, yNeighbor);

        return index;
    }

    protected PairInt getIndexIfOnlyOneNeighbor(Set<PairInt> points,
        PairInt point, int width, int height) {

        int x = point.getX();
        int y = point.getY();

        int count = 0;
        PairInt neighbor = null;

        for (int c = (x - 1); c <= (x + 1); c++) {
            if ((c < 0) || (c > (width - 1))) {
                continue;
            }
            for (int r = (y - 1); r <= (y + 1); r++) {
                if ((r < 0) || (r > (height - 1))) {
                    continue;
                }
                if ((c == x) && (r == y)) {
                    continue;
                }
                PairInt tmp = new PairInt(c, r);
                if (points.contains(tmp)) {
                    if (count > 0) {
                        return null;
                    }
                    neighbor = tmp;
                    count++;
                }
            }
        }

        if (count == 0) {
            return null;
        }

        return neighbor;
    }

    /**
     * returns avg r, avg g, avg b
     * @param points
     * @param theta
     * @param originalImage
     * @param addAlongX
     * @param addAmount
     * @return
     */
    private int[] getAvgMinMaxColor(PairIntArray points, GreyscaleImage theta,
        Image originalImage, boolean addAlongX, int addAmount) {

        int xOffset = theta.getXRelativeOffset();
        int yOffset = theta.getYRelativeOffset();

        double rSum = 0;
        double gSum = 0;
        double bSum = 0;

        int count = 0;

        for (int pIdx = 0; pIdx < points.getN(); pIdx++) {

            int x = points.getX(pIdx);
            int y = points.getY(pIdx);

            int ox = x + xOffset;
            int oy = y + yOffset;

            //TODO: this may need corrections for other orientations
            if (addAlongX) {
                ox += addAmount;
            } else {
                oy += addAmount;
            }

            if ((ox < 0) || (ox > (originalImage.getWidth() - 1))) {
                continue;
            }
            if ((oy < 0) || (oy > (originalImage.getHeight() - 1))) {
                continue;
            }

            int rgb = originalImage.getRGB(ox, oy);
            int r = (rgb >> 16) & 0xFF;
            int g = (rgb >> 8) & 0xFF;
            int b = rgb & 0xFF;

            rSum += r;
            gSum += g;
            bSum += b;

            count++;
        }

        if (count == 0) {
            return new int[]{0, 0, 0};
        }

        rSum /= (double)count;
        gSum /= (double)count;
        bSum /= (double)count;

        return new int[]{(int)Math.round(rSum), (int)Math.round(gSum),
            (int)Math.round(bSum)};
    }

    public void applyErosionFilter(GreyscaleImage img) {
        ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
        lt.applyFilter(img);
    }

    public void applyErosionFilterOnZeroes(GreyscaleImage img) {
        ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
        lt.applyFilterOnZeros(img);
    }

    public void applyInvert255(GreyscaleImage img) {
        // assumption that pixels lie in range 0 to 255

        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                int v = img.getValue(i, j);
                int vInv = 255 - v;
                img.setValue(i, j, vInv);
            }
        }
    }

    public void applyInvert255(Image img) {
        // assumption that pixels lie in range 0 to 255

        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                img.setRGB(i, j,
                    255 - img.getR(i, j),
                    255 - img.getG(i, j),
                    255 - img.getB(i, j));
            }
        }
    }

    public GreyscaleImage binImageToKeepZeros(GreyscaleImage img,
        int binFactor) {

        if (img == null) {
            throw new IllegalArgumentException("img cannot be null");
        }

        int w0 = img.getWidth();
        int h0 = img.getHeight();

        int w1 = w0/binFactor;
        int h1 = h0/binFactor;

        GreyscaleImage out = new GreyscaleImage(w1, h1, img.getType());
        out.setXRelativeOffset(Math.round(img.getXRelativeOffset()/binFactor));
        out.setYRelativeOffset(Math.round(img.getYRelativeOffset()/binFactor));

        for (int i = 0; i < w1; i++) {

            for (int j = 0; j < h1; j++) {

                int vSum = 0;
                int count = 0;
                boolean isZero = false;

                // if there's a zero in the binFactor x binFactor block,
                // v is set to 0

                for (int ii = (i*binFactor); ii < ((i + 1)*binFactor); ii++) {
                    for (int jj = (j*binFactor); jj < ((j + 1)*binFactor); jj++) {

                        if ((ii < 0) || (ii > (w0 - 1))) {
                            continue;
                        }
                        if ((jj < 0) || (jj > (h0 - 1))) {
                            continue;
                        }

                        int v = img.getValue(ii, jj);

                        if (v == 0) {
                            isZero = true;
                            vSum = 0;
                            break;
                        }

                        vSum += v;
                        count++;
                    }
                    if (isZero) {
                        break;
                    }
                }

                if (vSum > 0) {
                    float v = (float)vSum/(float)count;
                    vSum = Math.round(v);
                }

                out.setValue(i, j, vSum);
            }
        }

        return out;
    }

    public GreyscaleImage binImage(GreyscaleImage img, int binFactor) {

        if (img == null) {
            throw new IllegalArgumentException("img cannot be null");
        }

        int w0 = img.getWidth();
        int h0 = img.getHeight();

        int w1 = w0/binFactor;
        int h1 = h0/binFactor;

        GreyscaleImage out = new GreyscaleImage(w1, h1, img.getType());
        out.setXRelativeOffset(Math.round(img.getXRelativeOffset()/binFactor));
        out.setYRelativeOffset(Math.round(img.getYRelativeOffset()/binFactor));

        for (int i = 0; i < w1; i++) {

            for (int j = 0; j < h1; j++) {

                int vSum = 0;
                int count = 0;

                for (int ii = (i*binFactor); ii < ((i + 1)*binFactor); ii++) {
                    for (int jj = (j*binFactor); jj < ((j + 1)*binFactor); jj++) {

                        if ((ii < 0) || (ii > (w0 - 1))) {
                            continue;
                        }
                        if ((jj < 0) || (jj > (h0 - 1))) {
                            continue;
                        }

                        int v = img.getValue(ii, jj);

                        vSum += v;
                        count++;
                    }
                }

                if (count > 0) {
                    float v = (float)vSum/(float)count;
                    vSum = Math.round(v);
                }

                out.setValue(i, j, vSum);
            }
        }

        return out;
    }
    
    public Image binImage(Image img,  int binFactor) {

        if (img == null) {
            throw new IllegalArgumentException("img cannot be null");
        }

        int w0 = img.getWidth();
        int h0 = img.getHeight();

        int w1 = w0/binFactor;
        int h1 = h0/binFactor;

        Image out = new Image(w1, h1, !img.is64Bit);

        binImage(img, binFactor, out);

        return out;
    }

    public ImageExt binImage(ImageExt img,  int binFactor) {

        if (img == null) {
            throw new IllegalArgumentException("img cannot be null");
        }

        int w0 = img.getWidth();
        int h0 = img.getHeight();

        int w1 = w0/binFactor;
        int h1 = h0/binFactor;

        ImageExt out = new ImageExt(w1, h1, !img.is64Bit);

        binImage(img, binFactor, out);

        return out;
    }
    
    public int[] binArray(int[] a, Image img, int binFactor) {
        
        if (img == null) {
            throw new IllegalArgumentException("img cannot be null");
        }

        int w0 = img.getWidth();
        int h0 = img.getHeight();

        int w1 = w0/binFactor;
        int h1 = h0/binFactor;

        int[] output = new int[w1 * h1];

        for (int i = 0; i < w1; i++) {

            for (int j = 0; j < h1; j++) {

                int aSum = 0;
                int count = 0;

                for (int ii = (i*binFactor); ii < ((i + 1)*binFactor); ii++) {
                    for (int jj = (j*binFactor); jj < ((j + 1)*binFactor); jj++) {

                        if ((ii < 0) || (ii > (w0 - 1))) {
                            continue;
                        }
                        if ((jj < 0) || (jj > (h0 - 1))) {
                            continue;
                        }

                        int pixIdx2 = img.getInternalIndex(ii, jj);
                       
                        aSum += a[pixIdx2];
                        
                        count++;
                    }
                }

                if (count > 0) {
                    aSum = Math.round((float)aSum/(float)count);
                }
                
                int pixIdx = (j * w1) + i;
                
                output[pixIdx] = aSum;
            }
        }

        return output;
    }
    
    /**
     * given an array with indexes of pixels in reference
     * frame of img, expand the array to the size of
     * an image with (resultWidth, resultHeight) which
     * is roughly a factor of binFactor (number resolution loss
     * means need to pass in the result width and height
     * to calculate the output pixel indexes).
     * @param input
     * @param img
     * @param binFactor
     * @param resultWidth
     * @param resultHeight
     * @return 
     */
    public int[] unbinArray(int[] input, 
        Image img, int binFactor, int resultWidth,
        int resultHeight) {

        if (input == null) {
            throw new IllegalArgumentException("input cannot be null");
        }

        int w0 = img.getWidth();
        int h0 = img.getHeight();

        int w1 = resultWidth;
        int h1 = resultHeight;
        
        int[] output = new int[w1 * h1];
        
        for (int i = 0; i < w0; i++) {
            for (int j = 0; j < h0; j++) {

                int pixIdx = img.getInternalIndex(i, j);
                int v = input[pixIdx];
                
                int stop1 = ((i + 1)*binFactor);
                if (stop1 > (w1 - 1)) {
                    stop1 = w1;
                }
                for (int ii = (i*binFactor); ii < stop1; ii++) {
                    int stop2 = ((j + 1)*binFactor);
                    if (stop2 > (h1 - 1)) {
                        stop2 = h1;
                    }
                    for (int jj = (j*binFactor); jj < stop2; jj++) {
                        int pixIdx2 = (jj * w1) + ii;
                        output[pixIdx2] = v;
                    }
                    if (j == (h0 - 1)) {
                        // just in case excess unset past binFactor
                        for (int jj = stop2; jj < h1; jj++) {
                            int pixIdx2 = (jj * w1) + ii;
                            output[pixIdx2] = v;
                        }
                    }
                }
                if (i == (w0 - 1)) {
                    // just in case excess unset past binFastor
                    for (int ii = stop1; ii < w1; ii++) {
                        int stop2 = ((j + 1)*binFactor);
                        if (stop2 > (h1 - 1)) {
                            stop2 = h1;
                        }
                        for (int jj = (j*binFactor); jj < stop2; jj++) {
                            int pixIdx2 = (jj * w1) + ii;
                            output[pixIdx2] = v;
                        }
                        if (j == (h0 - 1)) {
                            // just in case excess unset
                            for (int jj = stop2; jj < h1; jj++) {
                                int pixIdx2 = (jj * w1) + ii;
                                output[pixIdx2] = v;
                            }
                        }
                    }
                }
            }
        }

        return output;
    }

    /**
     * given an array with indexes of pixels in reference
     * frame of img, expand the array to the size of
     * an image with (resultWidth, resultHeight) which
     * is roughly a factor of binFactor (number resolution loss
     * means need to pass in the result width and height
     * to calculate the output pixel indexes).
     * @param input
     * @param img
     * @param binFactor
     * @param resultWidth
     * @param resultHeight
     * @return 
     */
    public List<Set<PairInt>> unbinSets(List<Set<PairInt>> input, 
        Image img, int binFactor, int resultWidth,
        int resultHeight) {

        if (input == null) {
            throw new IllegalArgumentException("input cannot be null");
        }

        int w0 = img.getWidth();
        int h0 = img.getHeight();

        int w1 = resultWidth;
        int h1 = resultHeight;
        
        List<Set<PairInt>> output = new ArrayList<Set<PairInt>>();

        for (Set<PairInt> a : input) {
            
            Set<PairInt> aOut = new HashSet<PairInt>(a.size() * binFactor);
            output.add(aOut);
            
            for (PairInt p : a) {
                int i = p.getX();
                int j = p.getY();
                int pixIdx = img.getInternalIndex(i, j);
                int stop1 = ((i + 1)*binFactor);
                if (stop1 > (w1 - 1)) {
                    stop1 = w1;
                }
                for (int ii = (i*binFactor); ii < stop1; ii++) {
                    int stop2 = ((j + 1)*binFactor);
                    if (stop2 > (h1 - 1)) {
                        stop2 = h1;
                    }
                    for (int jj = (j*binFactor); jj < stop2; jj++) {
                        aOut.add(new PairInt(ii, jj));
                    }
                    if (j == (h0 - 1)) {
                        // just in case excess unset past binFactor
                        for (int jj = stop2; jj < h1; jj++) {
                            aOut.add(new PairInt(ii, jj));
                        }
                    }
                }
                if (i == (w0 - 1)) {
                    // just in case excess unset past binFastor
                    for (int ii = stop1; ii < w1; ii++) {
                        int stop2 = ((j + 1)*binFactor);
                        if (stop2 > (h1 - 1)) {
                            stop2 = h1;
                        }
                        for (int jj = (j*binFactor); jj < stop2; jj++) {
                            aOut.add(new PairInt(ii, jj));
                        }
                        if (j == (h0 - 1)) {
                            // just in case excess unset
                            for (int jj = stop2; jj < h1; jj++) {
                                aOut.add(new PairInt(ii, jj));
                            }
                        }
                    }
                }
            }
        }
        
        return output;
    }

    private void binImage(Image inputImg,  int binFactor, Image outputImg) {

        if (inputImg == null) {
            throw new IllegalArgumentException("img cannot be null");
        }

        int w0 = inputImg.getWidth();
        int h0 = inputImg.getHeight();

        int w1 = outputImg.getWidth();
        int h1 = outputImg.getHeight();

        for (int i = 0; i < w1; i++) {

            for (int j = 0; j < h1; j++) {

                long rSum = 0;
                long gSum = 0;
                long bSum = 0;

                int count = 0;

                for (int ii = (i*binFactor); ii < ((i + 1)*binFactor); ii++) {
                    for (int jj = (j*binFactor); jj < ((j + 1)*binFactor); jj++) {

                        if ((ii < 0) || (ii > (w0 - 1))) {
                            continue;
                        }
                        if ((jj < 0) || (jj > (h0 - 1))) {
                            continue;
                        }

                        int rgb = inputImg.getRGB(ii, jj);

                        int r = (rgb >> 16) & 0xFF;
                        int g = (rgb >> 8) & 0xFF;
                        int b = rgb & 0xFF;

                        rSum += r;
                        gSum += g;
                        bSum += b;

                        count++;
                    }
                }

                if (count > 0) {
                    rSum = Math.round((float)rSum/(float)count);
                    gSum = Math.round((float)gSum/(float)count);
                    bSum = Math.round((float)bSum/(float)count);
                }

                outputImg.setRGB(i, j, (int)rSum, (int)gSum, (int)bSum);
            }
        }
    }

    public GreyscaleImage unbinMask(GreyscaleImage mask, int binFactor,
        GreyscaleImage originalTheta) {

        if (mask == null) {
            throw new IllegalArgumentException("mask cannot be null");
        }

        if (originalTheta == null) {
            throw new IllegalArgumentException("originalTheta cannot be null");
        }

        GreyscaleImage out = originalTheta.createWithDimensions();

        int w0 = mask.getWidth();
        int h0 = mask.getHeight();

        int w1 = out.getWidth();
        int h1 = out.getHeight();

        for (int i = 0; i < w0; i++) {
            for (int j = 0; j < h0; j++) {

                int v = mask.getValue(i, j);

                int stop1 = ((i + 1)*binFactor);
                if (stop1 > (w1 - 1)) {
                    stop1 = w1;
                }
                for (int ii = (i*binFactor); ii < stop1; ii++) {
                    int stop2 = ((j + 1)*binFactor);
                    if (stop2 > (h1 - 1)) {
                        stop2 = h1;
                    }
                    for (int jj = (j*binFactor); jj < stop2; jj++) {
                        out.setValue(ii, jj, v);
                    }
                    for (int jj = stop2; jj < h1; jj++) {
                        out.setValue(ii, jj, v);
                    }
                }
                for (int ii = stop1; ii < w1; ii++) {
                    int stop2 = ((j + 1)*binFactor);
                    if (stop2 > (h1 - 1)) {
                        stop2 = h1;
                    }
                    for (int jj = (j*binFactor); jj < stop2; jj++) {
                        out.setValue(ii, jj, v);
                    }
                    for (int jj = stop2; jj < h1; jj++) {
                        out.setValue(ii, jj, v);
                    }
                }
            }
        }

        // TODO: consider correction for oversampling at location of skyline
        // using originalTheta

        return out;
    }

    public GreyscaleImage expandBy2UsingBilinearInterp(GreyscaleImage input) {

        if (input == null) {
            throw new IllegalArgumentException("input cannot be null");
        }

        int w1 = 2 * input.getWidth();
        int h1 = 2 * input.getHeight();

        return expandBy2UsingBilinearInterp(input, w1, h1);
    }

    /**
     * expand image to final size by a factor of 2, and use the given output
     * widths and heights which are expected to be either twice the input
     * or twice plus 1.
     * @param input
     * @param outWidth
     * @param outHeight
     * @return
     */
    public GreyscaleImage expandBy2UsingBilinearInterp(GreyscaleImage input,
        int outWidth, int outHeight) {

        if (input == null) {
            throw new IllegalArgumentException("input cannot be null");
        }

        int w0 = input.getWidth();
        int h0 = input.getHeight();

        if ((2*w0 != outWidth) && ((2*w0 + 1) != outWidth)) {
            throw new IllegalArgumentException(
            "outWidth should be 2 * input.getWidth() or (2 * input.getWidth()) + 1");
        }
        if ((2*h0 != outHeight) && ((2*h0 + 1) != outHeight)) {
            throw new IllegalArgumentException(
            "outHeight should be 2 * input.getHeight() or (2 * input.getHeight()) + 1");
        }

        GreyscaleImage out = input.createWithDimensions(outWidth, outHeight);

        for (int i = 0; i < outWidth; ++i) {
            for (int j = 0; j < outHeight; ++j) {
                int v = upsampleBy2UsingBilinearInterp(input, i, j);
                out.setValue(i, j, v);
            }
        }

        return out;
    }
    
    public int upsampleBy2UsingBilinearInterp(GreyscaleImage input,
        int x, int y) {
        
        int w0 = input.getWidth();
        int h0 = input.getHeight();

        if (((x & 1) != 1) && ((y & 1) != 1)) {
            int x0 = x/2;
            int y0 = y/2;
            if ((x0 < w0) && (y0 < h0)) {
                return input.getValue(x0, y0);
            }
        }

        float x0 = (float)x/2.f;
        float y0 = (float)y/2.f;

        if (x0 > (w0 - 1)) {
            x0 = w0 - 1;
        }
        if (y0 > (h0 - 1)) {
            y0 = h0 - 1;
        }

        double v2 = biLinearInterpolation(input, x0, y0);

        return (int)Math.round(v2);
    }
    
    // 
    public double upsampleBy2UsingBilinearInterp(double[][] input,
        int x, int y) {
        
        int w0 = input.length;
        int h0 = input[0].length;

        if (((x & 1) != 1) && ((y & 1) != 1)) {
            int x0 = x/2;
            int y0 = y/2;
            if ((x0 < w0) && (y0 < h0)) {
                return input[x0][y0];
            }
        }

        float x0 = (float)x/2.f;
        float y0 = (float)y/2.f;

        if (x0 > (w0 - 1)) {
            x0 = w0 - 1;
        }
        if (y0 > (h0 - 1)) {
            y0 = h0 - 1;
        }

        double v2 = biLinearInterpolation(input, x0, y0);

        return v2;
    }
    
    public double upsampleBy2UsingBilinearInterp(Complex[][] input,
        int x, int y, boolean calcForReal) {
        
        int w0 = input.length;
        int h0 = input[0].length;

        if (((x & 1) != 1) && ((y & 1) != 1)) {
            int x0 = x/2;
            int y0 = y/2;
            if ((x0 < w0) && (y0 < h0)) {
                if (calcForReal) {
                    return input[x0][y0].re();
                } else {
                    return input[x0][y0].im();
                }
            }
        }

        float x0 = (float)x/2.f;
        float y0 = (float)y/2.f;

        if (x0 > (w0 - 1)) {
            x0 = w0 - 1;
        }
        if (y0 > (h0 - 1)) {
            y0 = h0 - 1;
        }

        double v2 = biLinearInterpolation(input, x0, y0, calcForReal);

        return v2;
    }
    
    public double upsampleUsingBilinearInterp(Complex[][] input,
        int x, int y, boolean calcForReal, int factor) {
        
        int w0 = input.length;
        int h0 = input[0].length;

        if (((x & 1) != 1) && ((y & 1) != 1)) {
            int x0 = x/factor;
            int y0 = y/factor;
            if ((x0 < w0) && (y0 < h0)) {
                if (calcForReal) {
                    return input[x0][y0].re();
                } else {
                    return input[x0][y0].im();
                }
            }
        }

        float x0 = (float)x/(float)factor;
        float y0 = (float)y/(float)factor;

        if (x0 > (w0 - 1)) {
            x0 = w0 - 1;
        }
        if (y0 > (h0 - 1)) {
            y0 = h0 - 1;
        }

        double v2 = biLinearInterpolation(input, x0, y0, calcForReal);

        return v2;
    }
    
    public GreyscaleImage unbinImage(GreyscaleImage input, int binFactor) {

        if (input == null) {
            throw new IllegalArgumentException("input cannot be null");
        }

        int w0 = input.getWidth();
        int h0 = input.getHeight();

        GreyscaleImage out = input.createWithDimensions(
            binFactor* w0, binFactor * h0);

        int w1 = out.getWidth();
        int h1 = out.getHeight();

        for (int i = 0; i < w0; i++) {
            for (int j = 0; j < h0; j++) {

                int v = input.getValue(i, j);

                for (int ii = (i*binFactor); ii < ((i + 1)*binFactor); ii++) {
                    for (int jj = (j*binFactor); jj < ((j + 1)*binFactor); jj++) {
                        out.setValue(ii, jj, v);
                    }
                    for (int jj = ((j + 1)*binFactor); jj < h1; jj++) {
                        out.setValue(ii, jj, v);
                    }
                }
                for (int ii = ((i + 1)*binFactor); ii < w1; ii++) {
                    for (int jj = (j*binFactor); jj < ((j + 1)*binFactor); jj++) {
                        out.setValue(ii, jj, v);
                    }
                    for (int jj = ((j + 1)*binFactor); jj < h1; jj++) {
                        out.setValue(ii, jj, v);
                    }
                }
            }
        }

        return out;
    }

    public List<PairIntArray> unbinZeroPointLists(List<PairIntArray> zeroPointLists,
        int binFactor) {

        if (zeroPointLists == null) {
            throw new IllegalArgumentException("mask cannot be null");
        }

        List<PairIntArray> output = new ArrayList<PairIntArray>();

        for (PairIntArray zeroPointList : zeroPointLists) {

            PairIntArray transformed = new PairIntArray(zeroPointList.getN() *
                binFactor);

            for (int i = 0; i < zeroPointList.getN(); i++) {

                int x = zeroPointList.getX(i);
                int y = zeroPointList.getY(i);

                for (int ii = (x*binFactor); ii < ((x + 1)*binFactor); ii++) {
                    for (int jj = (y*binFactor); jj < ((y + 1)*binFactor); jj++) {
                        transformed.add(ii, jj);
                    }
                }
            }

            output.add(transformed);

        }

        return output;
    }

    public Set<PairInt> unbinZeroPointLists(Set<PairInt> zeroPoints,
        int binFactor) {

        if (zeroPoints == null) {
            throw new IllegalArgumentException("zeroPoints cannot be null");
        }

        Set<PairInt> output = new HashSet<PairInt>();

        for (PairInt zeroPoint : zeroPoints) {

            int x = zeroPoint.getX();
            int y = zeroPoint.getY();

            for (int ii = (x*binFactor); ii < ((x + 1)*binFactor); ii++) {
                for (int jj = (y*binFactor); jj < ((y + 1)*binFactor); jj++) {

                    PairInt p = new PairInt(ii, jj);

                    output.add(p);
                }
            }
        }

        return output;
    }

    public void printImageColorContrastStats(Image image, int rgbSkyAvg,
        int plotNumber) throws IOException {

        /*
        http://dilnxsrv.king.ac.uk/papers/wses2001.pdf
           Y   | 16  |   | 0.256  0.504  0.098 | |R|
           U = | 128 | + |-0.148 -0.291  0.439 | |G|
           V   | 128 |   | 0.439 -0.368 -0.072 | |B|
        */
        double[][] m = new double[3][];
        m[0] = new double[]{0.256, 0.504, 0.098};
        m[1] = new double[]{-0.148, -0.291, 0.439};
        m[2] = new double[]{0.439, -0.368, -0.072};

        int rSky = (rgbSkyAvg >> 16) & 0xFF;
        int gSky = (rgbSkyAvg >> 8) & 0xFF;
        int bSky = rgbSkyAvg & 0xFF;
        double[] yuvSky = MatrixUtil.multiply(m, new double[]{rSky, gSky, bSky});

        double t313 = Math.pow(3, (1./3.));

        int w = image.getWidth();
        int h = image.getHeight();
        int slice = 1;//10;

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        for (int i = 0; i < 6; i++) {

            int startCol = -1;
            int stopCol = -1;
            int startRow = -1;
            int stopRow = -1;
            boolean plotAlongRows = true;
            String labelSuffix = null;

            switch(i) {
                case 0:
                    //horizontal at low y
                    startCol = 0;
                    stopCol = w - 1;
                    startRow = slice;
                    stopRow = startRow + slice;
                    plotAlongRows = false;
                    labelSuffix = "horizontal stripe at low y";
                    break;
                case 1:
                    //horizontal at mid y
                    startCol = 0;
                    stopCol = w - 1;
                    startRow = (h - slice)/2 ;
                    stopRow = startRow + slice;
                    plotAlongRows = false;
                    labelSuffix = "horizontal stripe at mid y";
                    break;
                case 2:
                    //horizontal at high y
                    startCol = 0;
                    stopCol = w - 1;
                    startRow = (h - 2*slice) - 1;
                    stopRow = startRow + slice;
                    plotAlongRows = false;
                    labelSuffix = "horizontal stripe at high y";
                    break;
                case 3:
                    //vertical at low x
                    startCol = slice;
                    stopCol = startCol + slice;
                    startRow = 0;
                    stopRow = h - 1;
                    plotAlongRows = true;
                    labelSuffix = "vertical stripe at low x";
                    break;
                case 4:
                    //vertical at mid x
                    startCol = (w - slice)/2;
                    stopCol = startCol + slice;
                    startRow = 0;
                    stopRow = h - 1;
                    plotAlongRows = true;
                    labelSuffix = "vertical stripe at mid x";
                    break;
                default:
                    //vertical at high x
                    startCol = (w - 2*slice) - 1;
                    stopCol = startCol + slice;
                    startRow = 0;
                    stopRow = h - 1;
                    plotAlongRows = true;
                    labelSuffix = "vertical stripe at high x";
                    break;
            }

            // contrast as y
            // hue
            // blue
            // red
            float[] contrast = null;
            float[] hue = null;
            float[] red = null;
            float[] blue = null;
            float[] white = null;
            float[] axis = null;

            if (!plotAlongRows) {

                // plot along columns

                contrast = new float[w];
                hue = new float[w];
                red = new float[w];
                blue = new float[w];
                white = new float[w];
                axis = new float[w];

                for (int col = startCol; col <= stopCol; col++) {

                    int row = startRow;

                    int r = image.getR(col, row);
                    int g = image.getG(col, row);
                    int b = image.getB(col, row);
                    double[] rgb = new double[]{r, g, b};

                    double[] yuv = MatrixUtil.multiply(m, rgb);
                    yuv = MatrixUtil.add(yuv, new double[]{16, 128, 128});
                    double hueValue = Math.atan2(t313 * (g - b), ((2 * r) - g - b));

                    double contrastValue = (yuvSky[0] - yuv[0])/yuv[0];

                    double whiteValue = (r + g + b)/3.;

                    contrast[col] = (float)contrastValue;
                    hue[col] = (float)hueValue;
                    blue[col] = (float)b;
                    red[col] = (float)r;
                    white[col] = (float)whiteValue;

                    axis[col] = col;
                }

            } else {
                // plot along rows
                contrast = new float[h];
                hue = new float[h];
                red = new float[h];
                blue = new float[h];
                white = new float[h];
                axis = new float[h];

                for (int row = startRow; row <= stopRow; row++) {

                    int col = startCol;

                    int r = image.getR(col, row);
                    int g = image.getG(col, row);
                    int b = image.getB(col, row);
                    double[] rgb = new double[]{r, g, b};

                    double[] yuv = MatrixUtil.multiply(m, rgb);
                    yuv = MatrixUtil.add(yuv, new double[]{16, 128, 128});
                    double hueValue = Math.atan2(t313 * (g - b), ((2 * r) - g - b));

                    double contrastValue = (yuvSky[0] - yuv[0])/yuv[0];

                    double whiteValue = (r + g + b)/3.;

                    contrast[row] = (float)contrastValue;
                    hue[row] = (float)hueValue;
                    blue[row] = (float)b;
                    red[row] = (float)r;
                    white[row] = (float)whiteValue;

                    axis[row] = row;
                }

            }

            float xmn = MiscMath.findMin(axis);
            float xmx = MiscMath.findMax(axis);

            float ymn = MiscMath.findMin(contrast);
            float ymx = 1.1f * MiscMath.findMax(contrast);
            plotter.addPlot(xmn, xmx, ymn, ymx,
                axis, contrast, null, null, null, null,
                "contrast " + labelSuffix);

            ymn = MiscMath.findMin(hue);
            ymx = 1.1f * MiscMath.findMax(hue);
            plotter.addPlot(xmn, xmx, ymn, ymx,
                axis, hue, null, null, null, null, "hue " + labelSuffix);

            ymn = MiscMath.findMin(blue);
            ymx = 1.1f * MiscMath.findMax(blue);
            plotter.addPlot(xmn, xmx, ymn, ymx,
                axis, blue, null, null, null, null, "blue " + labelSuffix);

            ymn = MiscMath.findMin(red);
            ymx = 1.1f * MiscMath.findMax(red);
            plotter.addPlot(xmn, xmx, ymn, ymx,
                axis, red, null, null, null, null, "red " + labelSuffix);

            ymn = MiscMath.findMin(white);
            ymx = 1.1f * MiscMath.findMax(white);
            plotter.addPlot(xmn, xmx, ymn, ymx,
                axis, white, null, null, null, null, "white " + labelSuffix);

            plotter.writeFile(plotNumber);
        }
    }

    public double[] calculateYRGB(PairIntArray points, Image originalColorImage,
        int xOffset, int yOffset, boolean addAlongX, int addAmount) {

        double[][] m = new double[3][];
        m[0] = new double[]{0.256, 0.504, 0.098};
        m[1] = new double[]{-0.148, -0.291, 0.439};
        m[2] = new double[]{0.439, -0.368, -0.072};

        double avgY = 0;
        double avgR = 0;
        double avgG = 0;
        double avgB = 0;

        for (int i = 0; i < points.getN(); i++) {

            int x = points.getX(i);
            int y = points.getY(i);

            int ox = x + xOffset;
            int oy = y + yOffset;

            if (addAlongX) {
                ox += addAmount;
            } else {
                oy += addAmount;
            }
            if ((ox < 0) || (ox > (originalColorImage.getWidth() - 1))) {
                continue;
            }
            if ((oy < 0) || (oy > (originalColorImage.getHeight() - 1))) {
                continue;
            }

            int r = originalColorImage.getR(x, y);
            int g = originalColorImage.getG(x, y);
            int b = originalColorImage.getB(x, y);
            double[] rgb = new double[]{r, g, b};
            double[] yuv = MatrixUtil.multiply(m, rgb);

            avgY += yuv[0];

            avgR += r;
            avgG += g;
            avgB += b;
        }

        avgY /= (double)points.getN();
        avgR /= (double)points.getN();
        avgG /= (double)points.getN();
        avgB /= (double)points.getN();

        return new double[]{avgY, avgR, avgG, avgB};
    }

    private GreyscaleImage padUpToPowerOfTwo(GreyscaleImage input) {

        int w0 = input.getWidth();
        int h0 = input.getHeight();

        int w = 1 << (int)(Math.ceil(Math.log(w0)/Math.log(2)));
        int h = 1 << (int)(Math.ceil(Math.log(h0)/Math.log(2)));

        int xOffset = w - w0;
        int yOffset = h - h0;

        if (xOffset == 0 && yOffset == 0) {
            return input.copyImage();
        }

        int xOffsetOrig = input.getXRelativeOffset();
        int yOffsetOrig = input.getYRelativeOffset();

        GreyscaleImage output = new GreyscaleImage(w, h, input.getType());
        output.setXRelativeOffset(xOffset + xOffsetOrig);
        output.setYRelativeOffset(yOffset + yOffsetOrig);

        for (int i = 0; i < w0; ++i) {
            for (int j = 0; j < h0; ++j) {
                int v = input.getValue(i, j);
                output.setValue(i + xOffset, j + yOffset, v);
            }
        }

        return output;
    }

    private Complex[][] padUpToPowerOfTwoComplex(GreyscaleImage input) {

        int w0 = input.getWidth();
        int h0 = input.getHeight();

        int w = 1 << (int)(Math.ceil(Math.log(w0)/Math.log(2)));
        int h = 1 << (int)(Math.ceil(Math.log(h0)/Math.log(2)));

        int xOffset = w - w0;
        int yOffset = h - h0;

        Complex[][] output = new Complex[w][];
        for (int i = 0; i < w; ++i) {
            output[i] = new Complex[h];
        }

        if (xOffset == 0 && yOffset == 0) {
            for (int i = 0; i < w; ++i) {
                output[i] = new Complex[h];
                for (int j = 0; j < h; ++j) {
                    output[i][j] = new Complex(input.getValue(i, j), 0);
                }
            }
            return output;
        }

        for (int i = 0; i < w; ++i) {
            output[i] = new Complex[h];
            for (int j = 0; j < h; ++j) {
                if ((i < xOffset) || (j < yOffset)) {
                    output[i][j] = new Complex(0, 0);
                } else {
                    output[i][j] = new Complex(input.getValue(i - xOffset, j - yOffset), 0);
                }
            }
        }

        return output;
    }

    /**
     *
     * @param input a double array without complex interleaving items in columns
     * @return
     */
    private Complex[][] padUpToPowerOfTwoComplex(double[][] input) {

        int w0 = input.length;
        int h0 = input[0].length;

        int w = 1 << (int)(Math.ceil(Math.log(w0)/Math.log(2)));
        int h = 1 << (int)(Math.ceil(Math.log(h0)/Math.log(2)));

        int xOffset = w - w0;
        int yOffset = h - h0;

        Complex[][] output = new Complex[w][];
        for (int i = 0; i < w; ++i) {
            output[i] = new Complex[h];
        }

        if (xOffset == 0 && yOffset == 0) {

            for (int i = 0; i < w; ++i) {
                for (int j = 0; j < h; ++j) {
                    output[i][j] = new Complex(input[i][j], 0);
                }
            }

            return output;
        }

        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                if ((i < xOffset) || (j < yOffset)) {
                    output[i][j] = new Complex(0, 0);
                } else {
                    output[i][j] = new Complex(input[i - xOffset][j - yOffset], 0);
                }
            }
        }

        return output;
    }

    /**
     * add zeros to beginning of arrays to pad up to a size of the power of 2.
     * 
     * @param input
     * @return 
     */
    public Complex[][] padUpToPowerOfTwo(final Complex[][] input) {

        int n0 = input.length;
        int n1 = input[0].length;

        int nn0 = 1 << (int)(Math.ceil(Math.log(n0)/Math.log(2)));
        int nn1 = 1 << (int)(Math.ceil(Math.log(n1)/Math.log(2)));

        int offset0 = nn0 - n0;
        int offset1 = nn1 - n1;

        if (offset0 == 0 && offset1 == 0) {
            Complex[][] output = new Complex[n0][];
            for (int i0 = 0; i0 < n0; ++i0) {
                output[i0] = Arrays.copyOf(input[i0], n1);
            }
            return output;
        }

        Complex[][] output = new Complex[nn0][];
        for (int i0 = 0; i0 < nn0; ++i0) {
            output[i0] = new Complex[nn1];
        }

         for (int i0 = 0; i0 < nn0; ++i0) {
            for (int i1 = 0; i1 < nn1; ++i1) {
                if ((i0 < offset0) || (i1 < offset1)) {
                    output[i0][i1] = new Complex(0, 0);
                } else {
                    output[i0][i1] = input[i0 - offset0][i1 - offset1].copy();
                }
            }
        }

        return output;
    }

    /**
     * apply 2D FFT transform using the JFFTPack.
     *
     * @param input input image, which should probably be type full range int
     * @param forward if true, apply FFT transform, else inverse FFT transform
     */
    public void apply2DFFT2(GreyscaleImage input, boolean forward) {

         Complex1D[] ccOut = create2DFFT2WithSwapMajor(input, forward);

         assert(ccOut.length == input.getHeight());
         assert(ccOut[0].x.length == input.getWidth());

         for (int i0 = 0; i0 < ccOut.length; ++i0) {
             for (int i1 = 0; i1 < ccOut[i0].x.length; ++i1) {
                 double re = ccOut[i0].x[i1];
                 input.setValue(i1, i0, (int)re);
             }
         }
    }

    /**
     * apply 2D FFT transform using the efficient iterative power of 2 method
     * that uses the butterfly operation if image dimensions are a power of
     * 2, else uses an alternative.
     *
     * @param input
     * @param forward if true, apply FFT transform, else inverse FFT transform
     */
    public void apply2DFFT(GreyscaleImage input, boolean forward) {

        int n0 = input.getWidth();
        int n1 = input.getHeight();

        int nn0 = 1 << (int)(Math.ceil(Math.log(n0)/Math.log(2)));
        int nn1 = 1 << (int)(Math.ceil(Math.log(n1)/Math.log(2)));

        if (nn0 > n0 || nn1 > n1) {
            apply2DFFT2(input, forward);
            return;
        }

        // initialize matrix of complex numbers as real numbers from image (imaginary are 0's)
        Complex[][] cc = convertImage(input);

        Complex[][] ccOut = create2DFFT(cc, forward);

        input.fill(0);
        for (int col = 0; col < input.getWidth(); col++) {
            for (int row = 0; row < input.getHeight(); row++) {
                double re = ccOut[col][row].re();
                input.setValue(col, row, (int)re);
            }
        }

    }
    
    public Complex[][] create2DFFT(double[][] input, boolean forward) {

        // performs normalization by default
        return create2DFFT(input, true, forward);
    }
    
    /**
     * perform fft on input.  
     * @param input
     * @param doNormalize
     * @param forward
     * @return
     */
    public Complex[][] create2DFFT(final double[][] input, boolean doNormalize,
        boolean forward) {
        
        Complex[][] input2 = new Complex[input.length][];
        for (int i = 0; i < input.length; ++i) {
            input2[i] = new Complex[input[0].length];
            for (int j = 0; j < input[0].length; ++j) {
                input2[i][j] = new Complex(input[i][j], 0);
            }
        }
        
        return create2DFFT(input2, doNormalize, forward);
    }

    public Complex[][] create2DFFT(Complex[][] input, boolean forward) {

        // performs normalization by default
        return create2DFFT(input, true, forward);
    }

    /**
     * runtime complexity: is O(N*lg_2(N)) for N not power of 2,
     * else is 
     * 
     * perform fft on input.
     * @param input
     * @param doNormalize
     * @param forward
     * @return
     */
    public Complex[][] create2DFFT(final Complex[][] input, boolean doNormalize,
        boolean forward) {

        final int n0 = input.length;
        final int n1 = input[0].length;
        
        int nn0 = 1 << (int)(Math.ceil(Math.log(n0)/Math.log(2)));
        int nn1 = 1 << (int)(Math.ceil(Math.log(n1)/Math.log(2)));
        
        if (nn0 > n0 || nn1 > n1) {
            Complex1D[] input2 = convertToComplex1D(input);
            Complex1D[] output = create2DFFT2(input2, doNormalize, forward);
            Complex[][] output2 = convertToComplex(output);
            return output2;
        }
            
        Complex[][] output = copy(input);
        
        // padding is at front of cols and rows
    
        FFT fft = new FFT();
        if (!doNormalize) {
            fft.setToNotNormalize();
        }

        // ----- perform FFT by dimension 0 -----
        for (int i0 = 0; i0 < nn0; i0++) {
            if (forward) {
                output[i0] = fft.fft(output[i0]);
            } else {
                output[i0] = fft.ifft(output[i0]);
            }
        }
        
        // re-use array for the FFT by dimension 1
        Complex[] tmp = new Complex[nn0];

        /*
        nn0
         |
        \|/
        [0]  ..........nn1-1
        [1]  ..........nn1-1
        */
       
        // ----- perform the FFT on dimension 1 ------
        for (int i1 = 0; i1 < nn1; ++i1) {

            // store each column in tmp array and perform fft on it then
            // recopy values back into columns
            for (int i0 = 0; i0 < nn0; ++i0) {
                tmp[i0] = output[i0][i1];
            }

            if (forward) {
                tmp = fft.fft(tmp);
            } else {
                tmp = fft.ifft(tmp);
            }

            for (int i0 = 0; i0 < nn0; ++i0) {
                output[i0][i1] = tmp[i0];
            }
        }
        
        return output;        
    }

    /**
     * perform a 2-dimension FFT using the JFFTPack library using the input
     * img and return the results as a complex two dimensional array
     * which uses the format a[row][col].
     *
     * @param img
     * @param forward
     * @return
     */
    public Complex1D[] create2DFFT2WithSwapMajor(GreyscaleImage img, boolean forward) {

        // swap major axes for input to FFT 2D algorithm 2
        Complex1D[] cInput = new Complex1D[img.getWidth()];
        for (int i = 0; i < img.getWidth(); ++i) {
            cInput[i] = new Complex1D();
            cInput[i].x = new double[img.getHeight()];
            cInput[i].y = new double[img.getHeight()];
            for (int j = 0; j < img.getHeight(); ++j) {
                cInput[i].x[j] = img.getValue(i, j);
            }
        }

        Complex1D[] output2 = create2DFFT2(cInput, forward);

        return output2;
    }

    /**
     * apply 2D FFT transform using the efficient iterative power of 2 method
     * that uses the butterfly operation.
     *
     * @param img
     * @param forward if true, apply FFT transform, else inverse FFT transform
     * @return 2d fft results in format a[row][col]
     */
    public Complex[][] create2DFFTWithSwapMajor(GreyscaleImage img, boolean forward) {

        // normalize by default
        return create2DFFTWithSwapMajor(img, true, forward);
    }

    /**
     * apply 2D FFT transform .
     *
     * @param input
     * @param doNormalize perform FFT normalization if true
     * @param forward if true, apply FFT transform, else inverse FFT transform
     * @return 2d fft results in format a[row][col]
     */
    public Complex[][] create2DFFTWithSwapMajor(GreyscaleImage input,
        boolean doNormalize, boolean forward) {

        // initialize matrix of complex numbers as real numbers from image (imaginary are 0's)
        Complex[][] cc = convertImageWithSwapMajor(input);

        Complex[][] ccFFT = create2DFFT(cc, doNormalize, forward);

        return ccFFT;
    }

    /**
     * perform a 2-dimension FFT using the JFFTPack library.
     *
     * @param input double array of complex data in format double[nRows][2*nColumns]
     * where the column elements are alternately the complex real number and the
     * complex imaginary number.
     * @param forward
     * @return two dimensional complex array of size Complex[nRows][input.nCols/2)
     */
    public Complex1D[] create2DFFT2(Complex1D[] input, boolean forward) {

        // perform normalization by default
        return create2DFFT2(input, true, forward);
    }

    /**
     * perform a 2-dimension FFT using the JFFTPack library.
     * 
     * runtime complexity: is O(N*lg_2(N)) for N not power of 2.
     * 
     * @param input double array of complex data in format double[nRows][2*nColumns]
     * where the column elements are alternately the complex real number and the
     * complex imaginary number.
     * @param forward
     * @return two dimensional complex array of size Complex[nRows][input.nCols/2)
     */
    public Complex1D[] create2DFFT2(Complex1D[] input, boolean performNormalization,
        boolean forward) {

        final int n0 = input.length;
        final int n1 = input[0].x.length;

        Complex1D[] output = Arrays.copyOf(input, input.length);

        ComplexDoubleFFT fft1 = new ComplexDoubleFFT(n1);

        final double norm1 = 1./Math.sqrt(n1);

        // ----- perform FFT by dimension 0 -----
        for (int i0 = 0; i0 < n0; i0++) {

            if (forward) {
                fft1.ft(output[i0]);
            } else {
                fft1.bt(output[i0]);
            }

            // normalize the data
            if (performNormalization) {
                Complex1D a = output[i0];
                for (int idx = 0; idx < a.x.length; ++idx) {
                    a.x[idx] *= norm1;
                    a.y[idx] *= norm1;
                }
            }
        }
        
        // re-use array for the FFT by dimension 1 (across rows)
        Complex1D tmp = new Complex1D();
        tmp.x = new double[n0];
        tmp.y = new double[n0];

        ComplexDoubleFFT fft0 = new ComplexDoubleFFT(n0);

        final double norm0 = performNormalization ? (1./Math.sqrt(n0)) : 1.;
        
        // ----- perform the FFT on dimension 1 ------
        for (int i1 = 0; i1 < n1; ++i1) {

            // store each column in tmp array and perform fft on it then
            // recopy values back into columns
            for (int i0 = 0; i0 < n0; ++i0) {
                tmp.x[i0] = output[i0].x[i1];
                tmp.y[i0] = output[i0].y[i1];
            }

            if (forward) {
                fft0.ft(tmp);
            } else {
                fft0.bt(tmp);
            }

            for (int i0 = 0; i0 < n0; ++i0) {
                output[i0].x[i1] = tmp.x[i0] * norm0;
                output[i0].y[i1] = tmp.y[i0] * norm0;
            }
        }

        return output;
    }

    public void writeToImageWithSwapMajor(GreyscaleImage img, Complex[][] cc) {

        img.fill(0);

        // write back to original image
        for (int col = 0; col < img.getWidth(); col++) {
            for (int row = 0; row < img.getHeight(); row++) {
                double re = cc[row][col].re();
                img.setValue(col, row, (int)re);
            }
        }

    }

    public void writePositiveRealToImage(GreyscaleImage img, Complex[][] cc) {

        img.fill(0);

        // write back to original image
        for (int col = 0; col < img.getWidth(); col++) {
            for (int row = 0; row < img.getHeight(); row++) {
                double re = cc[col][row].re();
                double a = cc[col][row].abs();
                double p = cc[col][row].phase();
                if (re > 0) {
                    img.setValue(col, row, (int)re);
                }
            }
        }

    }

    /**
     * create a complex double array with format a[col][row]
     * @param input
     * @return
     */
    protected Complex[][] convertImage(GreyscaleImage input) {

        // initialize matrix of complex numbers as real numbers from image
        Complex[][] cc = new Complex[input.getWidth()][];

        for (int col = 0; col < input.getWidth(); col++) {

            cc[col] = new Complex[input.getHeight()];

            for (int row = 0; row < input.getHeight(); row++) {
                cc[col][row] = new Complex(input.getValue(col, row), 0);
            }
        }

        return cc;
    }

    /**
     * create an array of size double[nRows][*nCols] where the column elements
     *    are alternately the complex real and complex imaginary numbers
     *    (and the imaginary are 0 for this being real input).
     * @param input
     * @return
     */
    protected double[][] createInterleavedComplexSwapMajor(GreyscaleImage input) {

        int nCols = input.getWidth();
        int nRows = input.getHeight();

         // initialize matrix of complex numbers as real numbers from image
        double[][] d = new double[nRows][];

        for (int row = 0; row < nRows; row++) {

            d[row] = new double[2 * nCols];

            for (int col = 0; col < nCols; ++col) {
                d[row][2*col] = input.getValue(col, row);
                d[row][(2*col) + 1] = 0;
            }
        }

        return d;
    }

    /**
     * create a complex double array with format a[row][col]
     * @param input
     * @return
     */
    protected Complex[][] convertImageWithSwapMajor(GreyscaleImage input) {

        int nCols = input.getWidth();
        int nRows = input.getHeight();

        // initialize matrix of complex numbers as real numbers from image
        Complex[][] cc = new Complex[nRows][];

        for (int row = 0; row < nRows; row++) {

            cc[row] = new Complex[nCols];

            for (int col = 0; col < nCols; col++) {
                cc[row][col] = new Complex(input.getValue(col, row), 0);
            }
        }

        return cc;
    }

    protected Complex[][] convertImage(double[][] input) {

        int w = input.length;
        int h = input[0].length;

        // initialize matrix of complex numbers as real numbers from image
        Complex[][] cc = new Complex[w][];

        for (int col = 0; col < w; col++) {

            cc[col] = new Complex[h];

            for (int row = 0; row < h; row++) {
                cc[col][row] = new Complex(input[col][row], 0);
            }
        }

        return cc;
    }

    /**
     * NOT READY FOR USE YET
     *
     * @param input
     */
    public void applyDeconvolution(GreyscaleImage input) throws IOException {

        //TODO NOT READY FOR USE YET...

        applyWienerFilter(input);

    }

    /**
     * NOT READY FOR USE YET
     *
     * @param input
     */
    public void applyWienerFilter(GreyscaleImage input) throws IOException {

        //TODO NOT READY FOR USE YET...

        CannyEdgeFilterAdaptive cef = new CannyEdgeFilterAdaptive();

        // note, this is not scaled for total sum = 1 yet
        GreyscaleImage psf = cef.createGradientPSFForTesting();
        double sum = 0;
        for (int col = 0; col < psf.getWidth(); col++) {
            for (int row = 0; row < psf.getHeight(); row++) {
                int v = psf.getValue(col, row);
                sum += v;
            }
        }
        psf = padToNearestPowerOf2Dimensions(psf);
        Complex[][] psfNorm = new Complex[psf.getWidth()][];
        for (int col = 0; col < psf.getWidth(); col++) {
            psfNorm[col] = new Complex[psf.getHeight()];
            for (int row = 0; row < psf.getHeight(); row++) {
                int v = psf.getValue(col, row);
                double vn = v / sum;
                psfNorm[col][row] = new Complex(vn, 0);
            }
        }
        psfNorm = create2DFFT(psfNorm, true);

        // filter out low values?
        for (int i = 0; i < psfNorm.length; i++) {
            for (int j = 0; j < psfNorm[0].length; j++) {
                double r = psfNorm[i][j].re();
                if (r < 0.1) {
                    psfNorm[i][j] = new Complex(0, 0);
                }
            }
        }

        GreyscaleImage img0 = padToNearestPowerOf2Dimensions(input);

        ImageDisplayer.displayImage("before deconv", img0);

        Complex[][] imgCC = convertImage(img0);

        Complex[][] imgFFT = create2DFFT(imgCC, true);

        /*
        complex division:
           a times reciprocal of b

        reciprocal:
            double scale = re*re + im*im;
            r = Complex(re / scale, -im / scale);

        times:
            real = a.re * b.re - a.im * b.im;
            imag = a.re * b.im + a.im * b.re;
        */

        Complex[][] ccDeconv = new Complex[imgFFT.length][];
        int pXH =  psfNorm.length >> 1;
        int pYH =  psfNorm[0].length >> 1;
        for (int col = 0; col < imgFFT.length; col++) {

            ccDeconv[col] = new Complex[imgFFT[0].length];

            for (int row = 0; row < imgFFT[0].length; row++) {

                Complex v = imgFFT[col][row];

                // for convolution, each element of kernel and neighboring
                // pixel (including center pixel) were multiplied and result
                // is given to center pixel.

                // for deconvolution, the sums of the division are calculated

                Complex pixSum = new Complex(v.re(), v.im());

                for (int pXIdx = 0; pXIdx < psfNorm.length; pXIdx++) {
                    int pixXIdx = col + (pXIdx - pXH);

                    // correct for out of bounds of image
                    if (pixXIdx < 0) {
                        // replicate
                        pixXIdx = -1*pixXIdx - 1;
                        if (pixXIdx > (img0.getWidth() - 1)) {
                            pixXIdx = img0.getWidth() - 1;
                        }
                    } else if (pixXIdx >= img0.getWidth()) {
                        int diff = pixXIdx - img0.getWidth();
                        pixXIdx = img0.getWidth() - diff - 1;
                        if (pixXIdx < 0) {
                            pixXIdx = 0;
                        }
                    }

                    for (int pYIdx = 0; pYIdx < psfNorm.length; pYIdx++) {

                        if (psfNorm[pXIdx][pYIdx].abs() == 0) {
                            continue;
                        }

                        int pixYIdx = row + (pYIdx - pYH);

                        // correct for out of bounds of image
                        if (pixYIdx < 0) {
                            // replicate
                            pixYIdx = -1*pixYIdx - 1;
                            if (pixYIdx > (img0.getHeight() - 1)) {
                                pixYIdx = img0.getHeight() - 1;
                            }
                        } else if (pixYIdx >= img0.getHeight()) {
                            int diff = pixYIdx - img0.getHeight();
                            pixYIdx = img0.getHeight() - diff - 1;
                            if (pixYIdx < 0) {
                                pixYIdx = 0;
                            }
                        }

                        Complex vk = imgFFT[pixXIdx][pixYIdx];

                        if (vk.abs() == 0) {
                            continue;
                        }

                        Complex kRecip = psfNorm[pXIdx][pYIdx].reciprocal();

                        Complex vDivPSF = vk.times(kRecip);

                        pixSum = pixSum.plus(vDivPSF);
                    }
                }


                ccDeconv[col][row] = pixSum;
            }
        }

        GreyscaleImage img2 = img0.createFullRangeIntWithDimensions();

        writePositiveRealToImage(img2, ccDeconv);

        ImageDisplayer.displayImage("FFT(img0)/FFT(PSF)", img2);


        Complex[][] inverse = create2DFFT(ccDeconv, false);

        GreyscaleImage img4 = img0.createFullRangeIntWithDimensions();

        writePositiveRealToImage(img4, inverse);

        ImageDisplayer.displayImage("ifft of FFT(img0)/FFT(PSF)", img4);


        for (int i = 0; i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {
                double f = inverse[i][j].re();
                int v = input.getValue(i, j);
                if (v > 0 && f > 0) {
                    // apply it to the original image?  f*v or v or f?
                    input.setValue(i, j, v);
                } else {
                    input.setValue(i, j, 0);
                }
            }
        }
    }

    public GreyscaleImage padToNearestPowerOf2Dimensions(GreyscaleImage img) {

        int w = img.getWidth();
        int h = img.getHeight();

        boolean xIsPowerOf2 = MiscMath.isAPowerOf2(w);
        boolean yIsPowerOf2 = MiscMath.isAPowerOf2(h);
        if (xIsPowerOf2 && yIsPowerOf2) {
            return img.copyImage();
        }

        int w2 = w;
        int h2 = h;
        if (!xIsPowerOf2) {
            double p2X = Math.ceil(Math.log(w)/Math.log(2));
            w2 = (1 << (int)p2X);
        }
        if (!yIsPowerOf2) {
            double p2Y = Math.ceil(Math.log(h)/Math.log(2));
            h2 = (1 << (int)p2Y);
        }

        GreyscaleImage img2 = new GreyscaleImage(w2, h2, img.getType());

        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {
                int v = img.getValue(col, row);
                img2.setValue(col, row, v);
            }
        }

        return img2;
    }

    /**
     * read the image and store the non-zero pixels in a set.  note that negative
     * values will also be stored in the output set.
     * @param img
     * @return
     */
    public Set<PairInt> readNonZeroPixels(GreyscaleImage img) {

        Set<PairInt> set = new HashSet<PairInt>();

        for (int col = 0; col < img.getWidth(); col++) {
            for (int row = 0; row < img.getHeight(); row++) {
                int v = img.getValue(col, row);
                if (v != 0) {
                    set.add(new PairInt(col, row));
                }
            }
        }

        return set;
    }

    public void writeAsBinaryToImage(GreyscaleImage img, Set<PairInt>
        nonZeroPoints) {

        img.fill(0);

        for (PairInt p : nonZeroPoints) {
            int x = p.getX();
            int y = p.getY();
            img.setValue(x, y, 1);
        }

    }

    /**
     * find contiguous zeros in image and if the number of pixels in a groups
     * is less than contiguousZerosLimit, fill in the pixels with the
     * value of the neighboring pixels.
     * NOTE: this is set to use the 4-neighbor region, but can be set to use
     * 8-neighbors if needed.
     */
    public void fillInPixels(GreyscaleImage img, final int valueToFill,
        final int contiguousZerosLimit) {

        DFSContiguousValueFinder finder = new DFSContiguousValueFinder(img);
        finder.setMinimumNumberInCluster(1);
        finder.findGroups(valueToFill);

        int nGroups = finder.getNumberOfGroups();

        for (int i = 0; i < nGroups; ++i) {

            int n = finder.getNumberofGroupMembers(i);

            if (n <= contiguousZerosLimit) {

                PairIntArray group = finder.getXY(i);

                // find the adjacent non-zero pixels to these
                Set<PairInt> neighbors = new HashSet<PairInt>();
                for (int j = 0; j < group.getN(); ++j)  {
                    getNeighborsNotThisValue(img, group.getX(j), group.getY(j),
                        valueToFill, neighbors);
                }

                // get thier average intensities
                double avgV = 0;
                for (PairInt p : neighbors) {
                    int v = img.getValue(p.getX(), p.getY());
                    avgV += v;
                }
                avgV /= (double)neighbors.size();
                int vRepl = Math.round((float)avgV);
                for (int j = 0; j < group.getN(); ++j)  {
                    int x = group.getX(j);
                    int y = group.getY(j);
                    img.setValue(x, y, vRepl);
                }
            }
        }

    }

    public void getNeighborsNotThisValue(GreyscaleImage input, int x, int y,
        final int value, Set<PairInt> outputNeighbors) {

        int width = input.getWidth();
        int height = input.getHeight();

        for (int c = (x - 1); c <= (x + 1); c++) {
            if ((c < 0) || (c > (width - 1))) {
                continue;
            }
            for (int r = (y - 1); r <= (y + 1); r++) {
                if ((r < 0) || (r > (height - 1))) {
                    continue;
                }
                if ((c == x) && (r == y)) {
                    continue;
                }
                int v = input.getValue(c, r);
                if (v != value) {
                    PairInt p = new PairInt(c, r);
                    outputNeighbors.add(p);
                }
            }
        }
    }

    public Set<PairInt> getNeighbors(Image input, PairInt p) {

        int w = input.getWidth();
        int h = input.getHeight();

        int x = p.getX();
        int y = p.getY();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        Set<PairInt> nbrs = new HashSet<PairInt>();
        
        for (int k = 0; k < dxs.length; ++k) {
            int x2 = x + dxs[k];
            int y2 = y + dys[k];
            if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                continue;
            }
            nbrs.add(new PairInt(x2, y2));
        }
        
        return nbrs;
    }

    /**
     * NOT YET TESTED
     *
     http://en.wikipedia.org/wiki/Bilinear_interpolation
     http://en.wikipedia.org/wiki/Bilinear_interpolation#/media/File:Bilinear_interpolation_visualisation.svg
     * @param x
     * @param y
     * @return
     */
    public double biLinearInterpolation(GreyscaleImage gsImg, float x, float y) {

        double x1 = Math.floor(x);

        double x2 = Math.ceil(x);

        double y1 = Math.floor(y);

        double y2 = Math.ceil(y);

        double v1, v2;

        if (x1 == x2) {

            v1 = gsImg.getValue((int)x1, (int)y1);

            if (y1 == y2) {
                return v1;
            }

            v2 = gsImg.getValue((int)x1, (int)y2);

        } else {

            // interpolate over row y1
            v1 = ((x2 - x)/(x2 - x1)) * gsImg.getValue((int)x1, (int)y1) +
                ((x - x1)/(x2 - x1)) * gsImg.getValue((int)x2, (int)y1);

            if (y1 == y2) {
                return v1;
            }

            // interpolate over row y2
            v2 = ((x2 - x)/(x2 - x1)) * gsImg.getValue((int)x1, (int)y2) +
                ((x - x1)/(x2 - x1)) * gsImg.getValue((int)x2, (int)y2);
        }

        // interpolate the fraction of v1 and v2 over rows
        double v = ((y2 - y)/(y2 - y1)) * v1 + ((y - y1)/(y2 - y1)) * v2;

        return v;
    }
    
    /**
     * NOT YET TESTED
     *
     http://en.wikipedia.org/wiki/Bilinear_interpolation
     http://en.wikipedia.org/wiki/Bilinear_interpolation#/media/File:Bilinear_interpolation_visualisation.svg
     * @param x
     * @param y
     * @return
     */
    public double biLinearInterpolation(Complex[][] img, float x, float y,
        boolean calcForReal) {

        double x1 = Math.floor(x);

        double x2 = Math.ceil(x);

        double y1 = Math.floor(y);

        double y2 = Math.ceil(y);

        double v1, v2;

        if (x1 == x2) {

            if (calcForReal) {
                v1 = img[(int)x1][(int)y1].re();
            } else {
                v1 = img[(int)x1][(int)y1].im();
            }

            if (y1 == y2) {
                return v1;
            }

            if (calcForReal) {
                v2 = img[(int)x1][(int)y].re();
            } else {
                v2 = img[(int)x1][(int)y].im();
            }
        } else {

            double a, b;
            if (calcForReal) {
                a = img[(int)x1][(int)y1].re();
                b = img[(int)x2][(int)y1].re();
            } else {
                a = img[(int)x1][(int)y1].im();
                b = img[(int)x2][(int)y1].im();
            }
            
            // interpolate over row y1
            v1 = ((x2 - x)/(x2 - x1)) * a + ((x - x1)/(x2 - x1)) * b;

            if (y1 == y2) {
                return v1;
            }

            double c, d;
            if (calcForReal) {
                c = img[(int)x1][(int)y2].re();
                d = img[(int)x2][(int)y2].re();
            } else {
                c = img[(int)x1][(int)y2].im();
                d = img[(int)x2][(int)y2].im();
            }
            
            // interpolate over row y2
            v2 = ((x2 - x)/(x2 - x1)) * c + ((x - x1)/(x2 - x1)) * d;
        }

        // interpolate the fraction of v1 and v2 over rows
        double v = ((y2 - y)/(y2 - y1)) * v1 + ((y - y1)/(y2 - y1)) * v2;

        return v;
    }
    
    /**
     * NOT YET TESTED
     *
     http://en.wikipedia.org/wiki/Bilinear_interpolation
     http://en.wikipedia.org/wiki/Bilinear_interpolation#/media/File:Bilinear_interpolation_visualisation.svg
     * @param x
     * @param y
     * @return
     */
    public double biLinearInterpolation(double[][] gsImg, float x, float y) {

        double x1 = Math.floor(x);

        double x2 = Math.ceil(x);

        double y1 = Math.floor(y);

        double y2 = Math.ceil(y);

        double v1, v2;

        if (x1 == x2) {

            v1 = gsImg[(int)x1][(int)y1];

            if (y1 == y2) {
                return v1;
            }

            v2 = gsImg[(int)x1][(int)y];

        } else {

            // interpolate over row y1
            v1 = ((x2 - x)/(x2 - x1)) * gsImg[(int)x1][(int)y1] +
                ((x - x1)/(x2 - x1)) * gsImg[(int)x2][(int)y1];

            if (y1 == y2) {
                return v1;
            }

            // interpolate over row y2
            v2 = ((x2 - x)/(x2 - x1)) * gsImg[(int)x1][(int)y2] +
                ((x - x1)/(x2 - x1)) * gsImg[(int)x2][(int)y2];
        }

        // interpolate the fraction of v1 and v2 over rows
        double v = ((y2 - y)/(y2 - y1)) * v1 + ((y - y1)/(y2 - y1)) * v2;

        return v;
    }

    /**
    NOT YET TESTED
     http://en.wikipedia.org/wiki/Bilinear_interpolation
     http://en.wikipedia.org/wiki/Bilinear_interpolation#/media/File:Bilinear_interpolation_visualisation.svg
     *
     * @param x
     * @param y
     * @return
     */
    public double[] biLinearInterpolation(Image clrImg, float x, float y) {

        double x1 = Math.floor(x);

        double x2 = Math.ceil(x);

        double y1 = Math.floor(y);

        double y2 = Math.ceil(y);

        double r1, r2, g1, g2, b1, b2;

        if (x1 == x2) {

            r1 = clrImg.getR((int)x1, (int)y1);
            g1 = clrImg.getG((int)x1, (int)y1);
            b1 = clrImg.getB((int)x1, (int)y1);

            if (y1 == y2) {
                return new double[]{r1, g1, b1};
            }

            r2 = clrImg.getR((int)x1, (int)y2);
            g2 = clrImg.getG((int)x1, (int)y2);
            b2 = clrImg.getB((int)x1, (int)y2);

        } else {

            double v1X2Frac = ((x2 - x)/(x2 - x1));
            double v1X1Frac = ((x - x1)/(x2 - x1));

            // interpolate over row y1
            r1 = v1X2Frac * clrImg.getR((int)x1, (int)y1) +
                v1X1Frac * clrImg.getR((int)x2, (int)y1);

            g1 = v1X2Frac * clrImg.getG((int)x1, (int)y1) +
                v1X1Frac * clrImg.getG((int)x2, (int)y1);

            b1 = v1X2Frac * clrImg.getB((int)x1, (int)y1) +
                v1X1Frac * clrImg.getB((int)x2, (int)y1);

            if (y1 == y2) {
                return new double[]{r1, g1, b1};
            }

            // interpolate over row y2
            r2 = v1X2Frac * clrImg.getR((int)x1, (int)y2) +
                v1X1Frac * clrImg.getR((int)x2, (int)y2);

            g2 = v1X2Frac * clrImg.getG((int)x1, (int)y2) +
                v1X1Frac * clrImg.getG((int)x2, (int)y2);

            b2 = v1X2Frac * clrImg.getB((int)x1, (int)y2) +
                v1X1Frac * clrImg.getB((int)x2, (int)y2);
        }

        double v1Y2Frac = ((y2 - y)/(y2 - y1));
        double v1Y1Frac = ((y - y1)/(y2 - y1));

        // interpolate the fraction of v1 and v2 over rows
        double r = v1Y2Frac * r1 + v1Y1Frac * r2;

        // interpolate the fraction of v1 and v2 over rows
        double g = v1Y2Frac * g1 + v1Y1Frac * g2;

        // interpolate the fraction of v1 and v2 over rows
        double b = v1Y2Frac * b1 + v1Y1Frac * b2;

        return new double[]{r, g, b};
    }

    public void applyAdaptiveMeanThresholding(GreyscaleImage img) {

        applyAdaptiveMeanThresholding(img, 3);
    }

    public void applyAdaptiveMeanThresholding(GreyscaleImage img,
        int halfDimension) {

        GreyscaleImage imgM = img.copyImage();

        applyCenteredMean2(imgM, halfDimension);

        int c = 7;

        int foreground = 255;//1;
        int background = 0;

        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                int v = img.getValue(i, j);
                int m = imgM.getValue(i, j);
                int t = m - c;
                if (v > t) {
                    img.setValue(i, j, foreground);
                } else {
                    img.setValue(i, j, background);
                }
            }
        }

        imgM = null;

        //System.gc();
    }

    /**
     * @param input
     * @return
     */
    public GreyscaleImage makeWatershedFromAdaptiveMedian(GreyscaleImage input) {

        int w = input.getWidth();
        int h = input.getHeight();

        int[] dxs0, dys0;
        dxs0 = Misc.dx8;
        dys0 = Misc.dy8;
        GreyscaleImage tmpImg2 = input.copyImage();
        // fill in gaps of size 1 flooded the whole image. invert afterwards had same result.
        // increase the 0's by 1 pixel then invert however, is interesting.
        // where there is a '0', make all neighbors a '0':
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int v = input.getValue(i, j);
                if (v != 0) {
                    continue;
                }
                for (int k = 0; k < dxs0.length; ++k) {
                    int x1 = i + dxs0[k];
                    int y1 = j + dys0[k];
                    if (x1 < 0 || (x1 > (w - 1)) || y1 < 0 || (y1 > (h - 1))) {
                        continue;
                    }
                    tmpImg2.setValue(x1, y1, 0);
                }
            }
        }

        // invert image
        for (int i = 0; i < tmpImg2.getNPixels(); ++i) {
            int v = tmpImg2.getValue(i);
            tmpImg2.setValue(i, 255 - v);
        }

        WaterShed ws = new WaterShed();
        int[][] labelled = ws.createLabelledImage(tmpImg2.copyImage());
        GreyscaleImage wsImg = tmpImg2.createFullRangeIntWithDimensions();
        if (labelled == null) {
            return wsImg;
        }
        for (int j = 0; j < h; ++j) {
            for (int i = 0; i < w; ++i) {
                int v = labelled[i][j];
                wsImg.setValue(i, j, v);
            }
        }

        return wsImg;
    }

    /**
     * an algorithm that takes as input an image array of values and for
     * each unique values larger than 0, searches for connected components with a lowerLimitSize
     * and returns those as lists of pixels.
     * The runtime complexity is due to a DFS traversal so is dependent upon
     * the connectivity, that is O(|V| + |E|).
     * @param input
     * @param lowerLimitSize, the minimum length of a connected component to
     * keep and return as an item in the results list.
     * @return
     */
    public List<Set<PairInt>> extractConnectedComponents(int[][] input,
        int lowerLimitSize) {

        Map<Integer, Set<PairInt>> valuePixelsMap = new HashMap<Integer, Set<PairInt>>();

        int w = input.length;
        int h = Integer.MIN_VALUE;
        for (int i = 0; i < input.length; ++i) {
            for (int j = 0; j < input[i].length; ++j) {
                int v = input[i][j];
                if (v < 1) {
                    continue;
                }
                Integer key = Integer.valueOf(v);
                Set<PairInt> set = valuePixelsMap.get(key);
                if (set == null) {
                    set = new HashSet<PairInt>();
                    valuePixelsMap.put(key, set);
                }
                set.add(new PairInt(i, j));
                if (input[i].length > h) {
                    h = input[i].length;
                }
            }
        }

        List<Set<PairInt>> outputLists = new ArrayList<Set<PairInt>>();

        for (Entry<Integer, Set<PairInt>> entry : valuePixelsMap.entrySet()) {

            Set<PairInt> set = entry.getValue();

            DFSConnectedGroupsFinder finder = new DFSConnectedGroupsFinder();
            finder.setMinimumNumberInCluster(lowerLimitSize);
            finder.findConnectedPointGroups(set);

            int nGroups = finder.getNumberOfGroups();

            for (int i = 0; i < nGroups; ++i) {
                Set<PairInt> group = finder.getXY(i);
                Set<PairInt> set2 = new HashSet<PairInt>(group);
                outputLists.add(set2);
            }
        }

        return outputLists;
    }

    /**
     * an algorithm to operate on the results of adaptive mean algorithm,
     * that is, expecting the input is all 0 or positive numbers, and that
     * the 0's are the segments that are searched to filter out the connected
     * segments shorter than lowerLimitSize.
     * @param img
     * @param lowerLimitSize, the minimum length of a connected component to
     * keep and return as an item in the results list.
     * @param mask
     * @return
     */
    public List<Set<PairInt>> extractConnectedComponents(GreyscaleImage img,
        int lowerLimitSize, Set<PairInt> mask, int edgeValue) {

        if (img == null) {
            throw new IllegalArgumentException("img canot be null");
        }
        if (mask == null) {
            throw new IllegalArgumentException(
            "mask can be empty, but cannot be null");
        }

        int w = img.getWidth();
        int h = img.getHeight();

        // for input being adaptive mean, most pixels are 255, and the edges are '0'
        // so we are looking for the edges not in the mask.
        Set<PairInt> pixels = new HashSet<PairInt>();
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int v = img.getValue(i, j);
                if (v != edgeValue) {
                    continue;
                }
                PairInt p = new PairInt(i, j);
                if (!mask.contains(p)) {
                    pixels.add(p);
                }
            }
        }

        DFSConnectedGroupsFinder finder = new DFSConnectedGroupsFinder();
        finder.setMinimumNumberInCluster(lowerLimitSize);
        finder.findConnectedPointGroups(pixels);

        List<Set<PairInt>> outputLists = new ArrayList<Set<PairInt>>();

        int nGroups = finder.getNumberOfGroups();

        for (int i = 0; i < nGroups; ++i) {
            Set<PairInt> group = finder.getXY(i);
            Set<PairInt> set = new HashSet<PairInt>(group);
            outputLists.add(set);
        }

        return outputLists;
    }

    /**
     * create an image of the mean of the surrounding dimension x dimension
     * pixels for each pixel.  The calculation starts at 0 and the end
     * dimension pixels are averaged using the decreasing number of pixels.
     * <pre>
     * for example, image:
     * [10] [12] [12]
     * [10] [12] [12]
     *
     * for dimension = 2 becomes:
     * [11] [12] [12]
     * [11] [12] [12]
     * </pre>
     * runtime complexity is O(N_pixels)
     * This can be used as part of adaptive mean thresholding.
     *
     * @param img
     * @param dimension
     */
    public void applyBoxcarMean(GreyscaleImage img, int dimension) {

        if ((img.getWidth() < dimension) || (img.getHeight() < dimension)) {
            throw new IllegalArgumentException("dimension is larger than image"
                + " dimensions.  method not yet handling that.");
        }

        /*
        becomes efficient when dimension > 3

        sum along columns first using dynamic programming:
        sumCol[j=0] = sum_j=0_to_dim of row[i]
        sumCol[j=1] = sumCol[0] - row[j-1] + row[dim + j - 1]
        sumCol[j=2] = sumCol[1] - row[j-1] + row[dim + j - 1]
        */

        int w = img.getWidth();
        int h = img.getHeight();

        int[] mean = new int[img.getNPixels()];

        // sum along rows
        for (int i = 0; i < w; ++i) {
            int sum0 = 0;
            for (int j = 0; j < dimension; ++j) {
                sum0 += img.getValue(i, j);
            }
            mean[img.getInternalIndex(i, 0)] = sum0;
            for (int j = 1; j <= (h - dimension); ++j) {
                int vp = img.getValue(i, j - 1);
                int vl =  img.getValue(i, dimension + j - 1);
                sum0 = sum0 - vp + vl;
                mean[img.getInternalIndex(i, j)] = sum0;
            }
            // last dimension - 1 rows: sum along them, divide by count then mult by dimension
            for (int j = (h - dimension + 1); j < h; ++j) {
                float count = h - j;
                float sum = 0;
                for (int k = j; k < h; ++k) {
                    sum += img.getValue(i, k);
                }
                sum /= count;
                sum *= dimension;
                mean[img.getInternalIndex(i, j)] = Math.round(sum);
            }
        }

        int[] mean2 = new int[img.getNPixels()];

        // sum along columns
        for (int j = 0; j < h; ++j) {
            int sum0 = 0;
            for (int i = 0; i < dimension; ++i) {
                sum0 += mean[img.getInternalIndex(i, j)];
            }
            mean2[img.getInternalIndex(0, j)] = sum0;
            for (int i = 1; i <= (w - dimension); ++i) {
                int vp = mean[img.getInternalIndex(i - 1, j)];
                int vl = mean[img.getInternalIndex(dimension + i - 1, j)];
                sum0 = sum0 - vp + vl;
                mean2[img.getInternalIndex(i, j)] = sum0;
            }

            // last dimension - 1 cols: sum along them, divide by count then mult by dimension
            for (int i = (w - dimension + 1); i < w; ++i) {
                float count = h - i;
                float sum = 0;
                for (int k = i; k < w; ++k) {
                    sum += mean[img.getInternalIndex(k, j)];
                }
                sum /= count;
                sum *= dimension;
                mean2[img.getInternalIndex(i, j)] = Math.round(sum);
            }
        }

        // divide each value by dimension * dimension
        float dsq = dimension * dimension;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int v = mean2[img.getInternalIndex(i, j)];
                v = Math.round((float)v/dsq);
                img.setValue(i, j, v);
            }
        }
    }
    
    /**
     * create an image of the mean of the surrounding dimension x dimension
     * pixels for each pixel centered on each pixel.  For the starting
     * and ending (dimension/2) pixels, the average uses a decreasing
     * number of pixels.
     * <pre>
     * for example, image:
     * [10] [12] [12]
     * [10] [12] [12]
     *
     * for halfDimension = 1 becomes:
     * [11] [11] [12]
     * [11] [11] [12]
     * </pre>
     * runtime complexity is O(N_pixels) and never more than
     * 4 times O(N_pixels).
     * @param img
     * @param halfDimension the pixel center + and - this value in x and y
     * are averaged
     */
    public void applyCenteredMean2(GreyscaleImage img, int halfDimension) {

        SummedAreaTable sumTable = new SummedAreaTable();
        
        GreyscaleImage imgS = sumTable.createAbsoluteSummedAreaTable(img);
        
        imgS = sumTable.applyMeanOfWindowFromSummedAreaTable(imgS, 
            2*halfDimension + 1);
        
        img.resetTo(imgS);
    }

    public double[][] createUnitStandardDeviation(GreyscaleImage img, int halfDimension) {

        SummedAreaTable sumTable = new SummedAreaTable();
        
        GreyscaleImage imgM = sumTable.createAbsoluteSummedAreaTable(img);
        imgM = sumTable.applyMeanOfWindowFromSummedAreaTable(imgM, 
            2*halfDimension + 1);
         
        int w = img.getWidth();
        int h = img.getHeight();
        double[][] out = new double[w][];
        for (int i = 0; i < w; ++i) {
            out[i] = new double[h];
            for (int j = 0; j < h; ++j) {
                int m = imgM.getValue(i, j);
                double v = img.getValue(i, j) - m;
                if (m == 0) {
                    v = 0;
                } else {
                    v = v / (Math.sqrt(2)/m);
                }
                out[i][j] = v;
            }
        }
       
        return out;
    }
    
    /**
     * create an image of the mean of the surrounding dimension x dimension
     * pixels for each pixel centered on each pixel.  For the starting
     * and ending (dimension/2) pixels, the average uses a decreasing
     * number of pixels.
     * <pre>
     * for example, image:
     * [10] [12] [12]
     * [10] [12] [12]
     *
     * for halfDimension = 1 becomes:
     * [11] [11] [12]
     * [11] [11] [12]
     * </pre>
     * runtime complexity is O(N_pixels) and never more than
     * 4 times O(N_pixels).
     * @param img
     * @param halfDimension the pixel center + and - this value in x and y
     * are averaged
     */
    public void applyCenteredMean2(double[][] img, int halfDimension) {

        SummedAreaTable sumTable = new SummedAreaTable();
        
        double[][] imgS = sumTable.createAbsoluteSummedAreaTable(img);
        
        imgS = sumTable.applyMeanOfWindowFromSummedAreaTable(imgS, 
            2*halfDimension + 1);
        
        for (int i = 0; i < img.length; ++i) {
            System.arraycopy(imgS[i], 0, img[i], 0, imgS[i].length);
        }
    }
    
    /**
     * create an image of the mean of the surrounding dimension x dimension
     * pixels for each pixel centered on each pixel.  For the starting
     * and ending (dimension/2) pixels, the average uses a decreasing
     * number of pixels.
     * <pre>
     * for example, image:
     * [10] [12] [12]
     * [10] [12] [12]
     *
     * for halfDimension = 1 becomes:
     * [11] [11] [12]
     * [11] [11] [12]
     * </pre>
     * runtime complexity is O(N_pixels), but is also dependent upon 
     * halfDimension.  Prefer to use applyCenteredMean2 which is always
     * less than 4 times O(N) in runtime complexity.
     * @param img
     * @param halfDimension the pixel center + and - this value in x and y
     * are averaged
     */
    public void applyCenteredMean(GreyscaleImage img, int halfDimension) {

        if ((img.getWidth() < 2*halfDimension) ||
            (img.getHeight() < 2*halfDimension)) {
            throw new IllegalArgumentException("dimension is larger than image"
                + " dimensions.  method not yet handling that.");
        }
        
        /*
        becomes efficient when halfDimension > 1

        sum along columns first using dynamic programming, then rows
        */

        int dimension = 2*halfDimension + 1;

        int w = img.getWidth();
        int h = img.getHeight();

        int[] mean = new int[img.getNPixels()];

        int[] imgValues = img.getValues();

        // sum along rows
        for (int i = 0; i < w; ++i) {

            /* pixels before halfDimension
            halfDimension = 2
            0 1 2 3 4 5 6
                <
            */
            for (int j = 0; j < halfDimension; ++j) {
                float count = halfDimension - j;
                float sum = 0;
                for (int k = j; k < halfDimension; ++k) {
                    int pixIdx = img.getIndex(i, k);
                    sum += imgValues[pixIdx];
                }
                sum /= count;
                sum *= dimension;
                int pixIdx = img.getIndex(i, j);
                mean[pixIdx] = Math.round(sum);
            }

            /* pixels between halfDimension and j-halfDimension
            halfDimension = 2
            0 1 2 3 4 5
            |   *   |  sum from idx - halfDimension to idx + halfDimension, incl
            but store in idx
            */
            int sum0 = 0;
            for (int j = 0; j <= 2*halfDimension; ++j) {
                int pixIdx = img.getIndex(i, j);
                sum0 += imgValues[pixIdx];
            }
            int pixIdx = img.getIndex(i, halfDimension);
            mean[pixIdx] = sum0;
            /*
            halfDimension = 2
            0 1 2 3 4 5 6
              |   *   |
                |   *   |

            */
            for (int j = halfDimension + 1; j < (h - halfDimension); ++j) {
                pixIdx = img.getIndex(i, j - halfDimension - 1);
                int vp = imgValues[pixIdx];
                pixIdx = img.getIndex(i, j + halfDimension);
                int vl =  imgValues[pixIdx];
                sum0 = sum0 - vp + vl;
                pixIdx = img.getIndex(i, j);
                mean[pixIdx] = sum0;
            }
            /* last halfDimension pixels
            0 1 2 3   n=4, halfDimension = 2
                >
            */
            for (int j = (h - halfDimension); j < h; ++j) {
                float count = h - j + 1;
                float sum = 0;
                for (int k = (j - 1); k < h; ++k) {
                    pixIdx = img.getIndex(i, k);
                    sum +=  imgValues[pixIdx];
                }
                sum /= count;
                sum *= dimension;
                pixIdx = img.getIndex(i, j);
                mean[pixIdx] = Math.round(sum);
            }
        }

        // sum along columns
        for (int j = 0; j < h; ++j) {

            /* pixels before halfDimension
            halfDimension = 2
            0 1 2 3 4 5 6
                <
            */
            for (int i = 0; i < halfDimension; ++i) {
                float count = halfDimension - i;
                float sum = 0;
                for (int k = i; k < halfDimension; ++k) {
                    int pixIdx = img.getIndex(k, j);
                    sum += mean[pixIdx];
                }
                sum /= count;
                sum *= dimension;
                int pixIdx = img.getIndex(i, j);
                imgValues[pixIdx] =  Math.round(sum);
            }

            /* pixels between halfDimension and j-halfDimension
            halfDimension = 2
            0 1 2 3 4 5
            |   *   |  sum from idx - halfDimension to idx + halfDimension, incl
            but store in idx
            */
            int sum0 = 0;
            for (int i = 0; i <= 2*halfDimension; ++i) {
                int pixIdx = img.getIndex(i, j);
                sum0 += mean[pixIdx];
            }
            int pixIdx = img.getIndex(halfDimension, j);
            imgValues[pixIdx] = sum0;
            /*
            halfDimension = 2
            0 1 2 3 4 5 6
              |   *   |
                |   *   |

            */
            for (int i = halfDimension + 1; i < (w - halfDimension); ++i) {
                pixIdx = img.getIndex(i - halfDimension - 1, j);
                int vp = mean[pixIdx];
                pixIdx = img.getIndex(i + halfDimension, j);
                int vl =  mean[pixIdx];
                sum0 = sum0 - vp + vl;
                pixIdx = img.getIndex(i, j);
                imgValues[pixIdx] = sum0;
            }
            /* last halfDimension pixels
            0 1 2 3   n=4, halfDimension = 2
                >
            */
            for (int i = (w - halfDimension); i < w; ++i) {
                float count = w - i + 1;
                float sum = 0;
                for (int k = (i - 1); k < w; ++k) {
                    pixIdx = img.getIndex(k, j);
                    sum += mean[pixIdx];
                }
                sum /= count;
                sum *= dimension;
                pixIdx = img.getIndex(i, j);
                imgValues[pixIdx] = Math.round(sum);
            }
        }

        float dsq = dimension * dimension;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int pixIdx = img.getIndex(i, j);
                int v = imgValues[pixIdx];
                v = Math.round((float)v/dsq);
                imgValues[pixIdx] = v;
            }
        }
        for (int i = 0; i < img.getNPixels(); ++i) {
            img.setValue(i, imgValues[i]);
        }
    }

    public GreyscaleImage createSmallImage(int bufferSize, Set<PairInt> points,
        int pointValue) {

        //minMaxXY int[]{xMin, xMax, yMin, yMax}
        int[] minMaxXY = MiscMath.findMinMaxXY(points);

        int xOffset = minMaxXY[0] - bufferSize;
        int yOffset = minMaxXY[2] - bufferSize;

        int width = (minMaxXY[1] - minMaxXY[0]) + (2 * bufferSize);
        int height = (minMaxXY[3] - minMaxXY[2]) + (2 * bufferSize);

        GreyscaleImage img = new GreyscaleImage(width, height);
        img.setXRelativeOffset(xOffset);
        img.setYRelativeOffset(yOffset);
        for (PairInt p : points) {
            int x = p.getX() - xOffset;
            int y = p.getY() - yOffset;
            img.setValue(x, y, pointValue);
        }

        return img;
    }
    
    public Set<PairInt> extract2ndDerivPoints(GreyscaleImage img,
        Set<PairInt> filterToPoints, int nApprox) {
        
        Set<PairInt> set = extract2ndDerivPoints(img, nApprox, true);
        
        Set<PairInt> output = new HashSet<PairInt>();
        
        for (PairInt p : filterToPoints) {
            if (set.contains(p)) {
                output.add(p);
            }
        }
        
        return output;
    }

    public Set<PairInt> extract2ndDerivPoints(GreyscaleImage img) {

        GreyscaleImage gsImg = img.copyImage();

        applySecondDerivGaussian(gsImg, SIGMA.ONE, 0, 255);

        PairIntArray valueCounts = Histogram.createADescendingSortbyFrequencyArray(gsImg);
        int v0 = 0;
        int c0 = 0;
        int v1 = 0;
        int c1 = 0;
        for (int i = (valueCounts.getN() - 1); i > -1; --i) {
            int v = valueCounts.getX(i);
            int c = valueCounts.getY(i);
            if (v0 == 0) {
                if (c > 12) {
                    v0 = v;
                    c0 = c;
                }
            } else if (c < (2.5 * c0)) {
                v1 = v;
                c1 = c;
            } else {
                break;
            }
        }

        Set<PairInt> pixels = new HashSet<PairInt>();
        for (int i = 0; i < gsImg.getNPixels(); ++i) {
            int v = gsImg.getValue(i);
            int x = gsImg.getCol(i);
            int y = gsImg.getRow(i);
            if (v >= v1) {
                pixels.add(new PairInt(x, y));
            }
        }

        log.info("before nPoints=" + pixels.size());

        reduceTo4NeighborCentroids(pixels);

        log.info("after nPoints=" + pixels.size());

        return pixels;
    }

    /**
     * NOT READY FOR USE YET
     * extract the high value points in the second derivative gaussian of
     * img to a number of points less than or equal to maxNPoints and
     * if the variable reduceForNoise is true, then look for patterns
     * of noise and reduce the maximum value extracted from the 2nd deriv
     * points until no noise patterns are seen.
     * @param img
     * @param maxNPoints
     * @param reduceForNoise
     * @return
     */
    public Set<PairInt> extract2ndDerivPoints(GreyscaleImage img, int maxNPoints,
        boolean reduceForNoise) {

        GreyscaleImage gsImg = img.copyImage();

        applySecondDerivGaussian(gsImg, SIGMA.ONE, 0, 255);

        PairIntArray valueCounts = 
            Histogram.createADescendingSortByKeyArray(gsImg);
        int nTot = 0;
        int v1 = -1;
        for (int i = 0; i < valueCounts.getN(); ++i) {
            int c = valueCounts.getY(i);
            int nTmp = nTot + c;
            if (nTmp < maxNPoints) {
                nTot += c;
                v1 = valueCounts.getX(i);
            } else {
                if (v1 == -1) {
                    v1 = valueCounts.getX(i);
                }
                break;
            }
        }

        int w = gsImg.getWidth();
        int h = gsImg.getHeight();

        Set<PairInt> pixels = new HashSet<PairInt>();
        for (int i = 0; i < gsImg.getNPixels(); ++i) {
            int v = gsImg.getValue(i);
            if (v >= v1) {
                int x = gsImg.getCol(i);
                int y = gsImg.getRow(i);

                // avoid points on image boundaries
                if (x == 0 || y == 0 || (x > (w - 1)) || (y > (h - 1))) {
                    continue;
                }
                pixels.add(new PairInt(x, y));
            }
        }

        log.info("before nPoints=" + pixels.size());

        reduceTo4NeighborCentroids(pixels);

        log.info("after nPoints=" + pixels.size());

        if (reduceForNoise) {
            // look for patterns of noise and reduce v1 until not present
        }

        return pixels;
    }
    
    /**
     * NOT READY FOR USE YET
     * extract the high value points in the second derivative gaussian of
     * img to a number of points less than or equal to maxNPoints and
     * if the variable reduceForNoise is true, then look for patterns
     * of noise and reduce the maximum value extracted from the 2nd deriv
     * points until no noise patterns are seen.
     * @param img
     * @param labels labels for segmented regions for which
     * the maximum number of points are applied to limit the
     * total number returned.
     * @param maxNPoints
     * @param reduceForNoise
     * @return
     */
    public List<Set<PairInt>> extract2ndDerivPoints(
        GreyscaleImage img, int[] labels, int maxNPoints,
        boolean reduceForNoise) {

        int[] secondDerivs = 
            performSecondDerivGaussian(img, SIGMA.ONE);
        for (int idx = 0; idx < secondDerivs.length; ++idx) {
            if (secondDerivs[idx] < 0) {
                secondDerivs[idx] *= -1;
            }
        }
        
        TIntObjectMap<TIntSet> labelMap 
            = new TIntObjectHashMap<TIntSet>();
        for (int pixIdx = 0; pixIdx < labels.length; ++pixIdx) {
            int label = labels[pixIdx];
            TIntSet set = labelMap.get(label);
            if (set == null) {
                set = new TIntHashSet();
                labelMap.put(label, set);
            }
            set.add(pixIdx);
        }
        
        List<Set<PairInt>> output = new ArrayList<Set<PairInt>> ();
        
        TIntObjectIterator<TIntSet> iter = labelMap.iterator();
        
        for (int ii = labelMap.size(); ii-- > 0;) {

            iter.advance();
            
            int label = iter.key();
            
            TIntSet indexes = iter.value();
            PairIntArray valueCounts = 
                Histogram.createADescendingSortByKeyArray(
                indexes, secondDerivs);
                    
            int nTot = 0;
            int v1 = -1;
            for (int i = 0; i < valueCounts.getN(); ++i) {
                int c = valueCounts.getY(i);
                int nTmp = nTot + c;
                if (nTmp < maxNPoints) {
                    nTot += c;
                    v1 = valueCounts.getX(i);
                } else {
                    if (v1 == -1) {
                        v1 = valueCounts.getX(i);
                    }
                    break;
                }
            }

            int w = img.getWidth();
            int h = img.getHeight();

            Set<PairInt> pixels = new HashSet<PairInt>();
            TIntIterator iter2 = indexes.iterator();
            while (iter2.hasNext()) {
                
                int idx = iter2.next();
                
                int v = secondDerivs[idx];
                
                int x = img.getCol(idx);
                int y = img.getRow(idx);
                PairInt p = new PairInt(x, y);
                
                if (v >= v1) {
                    // avoid points on image boundaries
                    if (x == 0 || y == 0 || (x > (w - 1)) || (y > (h - 1))) {
                        continue;
                    }
                    pixels.add(p);
                }
            }

            log.info("before nPoints=" + pixels.size());

            reduceTo4NeighborCentroids(pixels);

            log.info("after nPoints=" + pixels.size());
            
            output.add(pixels);
        }

        if (reduceForNoise) {
            // look for patterns of noise and reduce v1 until not present
        }

        return output;
    }

    private void reduceTo4NeighborCentroids(Set<PairInt> pixels) {

        Set<PairInt> processed = new HashSet<PairInt>();

        Set<PairInt> output = new HashSet<PairInt>();

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;

        Set<PairInt> neighbors = new HashSet<PairInt>();

        for (PairInt p : pixels) {

            if (processed.contains(p)) {
                continue;
            }

            curveHelper.findNeighbors(p.getX(), p.getY(), pixels, processed,
                dxs, dys, neighbors);

            processed.add(p);
            processed.addAll(neighbors);

            if (neighbors.size() == 0) {
                output.add(p);
            } else {
                double[] xyCen = curveHelper.calculateXYCentroids(neighbors);
                int x = (int)Math.round(xyCen[0]);
                int y = (int)Math.round(xyCen[1]);
                assert(Math.abs(x - p.getX()) <= 2);
                assert(Math.abs(y - p.getY()) <= 2);
                output.add(new PairInt(x, y));
            }
        }

        pixels.clear();
        pixels.addAll(output);

    }

    public double determineLowerThreshold(GreyscaleImage input,
        double lowThresholdFractionOfTotal) {

        int n = input.getNPixels();

        // value, count
        PairIntArray sortedFreq = Histogram.createADescendingSortbyFrequencyArray(input);
        double sum = 0;
        for (int i = 0; i < sortedFreq.getN(); ++i) {
            sum += sortedFreq.getY(i);
        }

        int thresh = (int)Math.round(lowThresholdFractionOfTotal * sum);

        if (thresh == 0) {
            return 0;
        }

        double critValue = -1;
        double sum2 = 0;
        for (int i = (sortedFreq.getN() - 1); i > -1; --i) {
            sum2 += sortedFreq.getY(i);
            if (sum2 > thresh) {
                return sortedFreq.getX(i);
            }
        }

        return 0;
    }

    public double highPassIntensityFilter(GreyscaleImage input,
        double lowThresholdFractionOfTotal) {

        int n = input.getNPixels();

        double critValue = determineLowerThreshold(input,
            lowThresholdFractionOfTotal);

        for (int i = 0; i < n; ++i) {
            int v = input.getValue(i);
            if (v < critValue) {
                input.setValue(i, 0);
            }
        }

        return critValue;
    }

    public void twoLayerIntensityFilter(GreyscaleImage input,
        double highThresholdFractionOfTotal) {

        int n = input.getNPixels();

        // value, count
        PairIntArray sortedFreq = Histogram.createADescendingSortbyFrequencyArray(input);
        double sum = 0;
        for (int i = 0; i < sortedFreq.getN(); ++i) {
            sum += sortedFreq.getY(i);
        }

        int thresh = (int)Math.round(highThresholdFractionOfTotal * sum);

        if (thresh == 0) {
            return;
        }

        double critValue = -1;
        double sum2 = 0;
        for (int i = (sortedFreq.getN() - 1); i > -1; --i) {
            sum2 += sortedFreq.getY(i);
            if (sum2 > thresh) {
                critValue = sortedFreq.getX(i);
                break;
            }
        }

        Stack<Integer> stack = new Stack<Integer>();
        Set<Integer> visited = new HashSet<Integer>();
        for (int i = 0; i < n; ++i) {
            int v = input.getValue(i);
            if (v >= critValue) {
                stack.add(Integer.valueOf(i));
            }
        }

        int critValue2 = (int)Math.round(0.5 * critValue);

        int w = input.getWidth();
        int h = input.getHeight();
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        GreyscaleImage tmp = input.createWithDimensions();

        while (!stack.isEmpty()) {
            Integer pixIndex = stack.pop();
            if (visited.contains(pixIndex)) {
                continue;
            }

            int x = input.getCol(pixIndex.intValue());
            int y = input.getRow(pixIndex.intValue());

            tmp.setValue(pixIndex.intValue(), 255);

            for (int i = 0; i < dxs.length; ++i) {
                int x2 = x + dxs[i];
                int y2 = y + dys[i];
                if (x2 < 0 || (x2 > (w - 1)) || (y2 < 0) || (y2 > (h - 1))) {
                    continue;
                }
                int v = input.getValue(x2, y2);
                if (v > critValue2) {
                    int pixIdx2 = input.getInternalIndex(x2, y2);
                    tmp.setValue(pixIdx2, 255);
                    stack.add(Integer.valueOf(pixIdx2));
                }
            }
            visited.add(pixIndex);
        }

        input.resetTo(tmp);
    }

    public Complex1D[] convertToComplex1D(Complex[][] input) {
        
        int n0 = input.length;
        int n1 = input[0].length;
        
        Complex1D[] output = new Complex1D[n0];
        for (int i = 0; i < n0; ++i) {
            output[i] = new Complex1D();
            output[i].x = new double[n1];
            output[i].y = new double[n1];
            for (int j = 0; j < n1; ++j) {
                output[i].x[j] = input[i][j].re();
                output[i].y[j] = input[i][j].im();
            }
        }
        
        return output;
    }

    public Complex[][] copy(Complex[][] input) {
        
        int n0 = input.length;
        
        Complex[][] output = new Complex[n0][];
        for (int i = 0; i < n0; ++i) {
            output[i] = Arrays.copyOf(input[i], input[i].length);
        }
        
        return output;
    }

    public double[][] copy(double[][] input) {
        
        int n0 = input.length;
        
        double[][] output = new double[n0][];
        for (int i = 0; i < n0; ++i) {
            output[i] = Arrays.copyOf(input[i], input[i].length);
        }
        
        return output;
    }
    
    /**
     * output is column major format
     * @param input
     * @return 
     */
    public double[][] copy(GreyscaleImage input) {
        
        int n0 = input.getNPixels();
        int w = input.getWidth();
        int h = input.getHeight();
        
        double[][] output = new double[w][h];
        for (int i = 0; i < w; ++i) {
            output[i] = new double[h];
            for (int j = 0; j < h; ++j) {
                output[i][j] = input.getValue(i, j);
            }
        }
        
        return output;
    }
    
    /**
     * output is row major format
     * @param input
     * @return 
     */
    public float[][] copyToRowMajor(GreyscaleImage input) {
        
        int n0 = input.getNPixels();
        int h = input.getWidth();
        int w = input.getHeight();
        
        float[][] output = new float[w][h];
        for (int i = 0; i < w; ++i) {
            output[i] = new float[h];
            for (int j = 0; j < h; ++j) {
                output[i][j] = input.getValue(j, i);
            }
        }
        
        return output;
    }
    
    public Complex[][] convertToComplex(Complex1D[] input) {
        
        int n0 = input.length;
        int n1 = input[0].x.length;
            
        Complex[][] output = new Complex[n0][];
        for (int i = 0; i < n0; ++i) {
            output[i] = new Complex[n1];
            for (int j = 0; j < n1; ++j) {
                output[i][j] = new Complex(input[i].x[j], input[i].y[j]);
            }
        }
        
        return output;
    }

    /**
     * convolve (x, y) with 1D kernel using sparseValueMap, but user must assert
     * that (x,y) += half kernel size is all contained within the sparseValueMap
     * because no approximations are made when a value is not in the map
     * and the results will not be normalized correctly.
     * 
     * @param sparseValueMap
     * @param x
     * @param y
     * @param kernel
     * @param calcX if true, convolve for x, else for y
     * @return 
     */
    public float convolve1D(Map<PairInt, ? extends Number> sparseValueMap, int x, int y, 
        float[] kernel, boolean calcX) {
        
        int h = (kernel.length - 1) >> 1;
        
        float sum = 0;
        for (int g = 0; g < kernel.length; g++) {
            float gg = kernel[g];
            if (gg == 0) {
                continue;
            }
            int x2, y2;
            if (calcX) {
                x2 = x + g - h;
                y2 = y;
            } else {
                y2 = y + g - h;
                x2 = x;
            }
            
            Number v = sparseValueMap.get(new PairInt(x2, y2));
            
            if (v == null) {
                throw new IllegalArgumentException(
                    "x,y += half kernel length is not in map");
            }
              
            sum += (gg * v.floatValue());
        }
        
        return sum;
    }
    
    public TIntObjectMap<TIntSet> unbin(
        TIntObjectMap<TIntSet> binnedMap,
        int binFactor, int binnedWidth, int binnedHeight,
        int resultWidth, int resultHeight) {
        
        if (binnedMap == null) {
            throw new IllegalArgumentException(
            "binnedMap cannot be null");
        }

        TIntObjectMap<TIntSet> output 
            = new TIntObjectHashMap<TIntSet>(
            binnedMap.size() * binFactor);

        int w0 = binnedWidth;
        int h0 = binnedHeight;
        
        int w1 = resultWidth;
        int h1 = resultHeight;
        
        TIntObjectIterator<TIntSet> iter =
            binnedMap.iterator();
        
        for (int lIdx = 0; lIdx < binnedMap.size(); ++lIdx) {
            
            iter.advance();
            
            TIntSet p0 = iter.value();
            TIntSet pOut = new TIntHashSet(p0.size() * binFactor);
            output.put(iter.key(), pOut);
            
            TIntIterator iter2 = p0.iterator();
            while (iter2.hasNext()) {
                int pixIdx = iter2.next();
                int j = pixIdx/w0;
                int i = pixIdx - (j * w0);
                int stop1 = ((i + 1)*binFactor);
                if (stop1 > (w1 - 1)) {
                    stop1 = w1;
                }
                for (int ii = (i*binFactor); ii < stop1; ii++) {
                    int stop2 = ((j + 1)*binFactor);
                    if (stop2 > (h1 - 1)) {
                        stop2 = h1;
                    }
                    for (int jj = (j*binFactor); jj < stop2; jj++) {
                        int pixIdx2 = (jj * w1) + ii;
                        pOut.add(pixIdx2);
                    }
                    if (j == (h0 - 1)) {
                        // just in case excess unset past binFactor
                        for (int jj = stop2; jj < h1; jj++) {
                            int pixIdx2 = (jj * w1) + ii;
                            pOut.add(pixIdx2);
                        }
                    }
                }
                if (i == (w0 - 1)) {
                    // just in case excess unset past binFastor
                    for (int ii = stop1; ii < w1; ii++) {
                        int stop2 = ((j + 1)*binFactor);
                        if (stop2 > (h1 - 1)) {
                            stop2 = h1;
                        }
                        for (int jj = (j*binFactor); jj < stop2; jj++) {
                            int pixIdx2 = (jj * w1) + ii;
                            pOut.add(pixIdx2);
                        }
                        if (j == (h0 - 1)) {
                            // just in case excess unset
                            for (int jj = stop2; jj < h1; jj++) {
                                int pixIdx2 = (jj * w1) + ii;
                                pOut.add(pixIdx2);
                            }
                        }
                    }
                }
            }
        }
        
        return output;
    }
    
    public List<PairIntArray> unbinLists(
        List<PairIntArray> contigBinnedList, 
        int binFactor, 
        int binnedWidth, int binnedHeight,
        int resultWidth, int resultHeight) {
        
        if (contigBinnedList == null) {
            throw new IllegalArgumentException(
            "contigBinnedList cannot be null");
        }

        List<PairIntArray> output 
            = new ArrayList<PairIntArray>(contigBinnedList.size());

        int w0 = binnedWidth;
        int h0 = binnedHeight;
        
        int w1 = resultWidth;
        int h1 = resultHeight;
        
        for (int lIdx = 0; lIdx < contigBinnedList.size(); ++lIdx) {
            
            PairIntArray p0 = contigBinnedList.get(lIdx);
            PairIntArray pOut = new PairIntArray(p0.getN() * binFactor);
            output.add(pOut);
            
            for (int idx = 0; idx < p0.getN(); ++idx) {
                int i = p0.getX(idx);
                int j = p0.getY(idx);
                int stop1 = ((i + 1)*binFactor);
                if (stop1 > (w1 - 1)) {
                    stop1 = w1;
                }
                for (int ii = (i*binFactor); ii < stop1; ii++) {
                    int stop2 = ((j + 1)*binFactor);
                    if (stop2 > (h1 - 1)) {
                        stop2 = h1;
                    }
                    for (int jj = (j*binFactor); jj < stop2; jj++) {
                        pOut.add(ii, jj);
                    }
                    if (j == (h0 - 1)) {
                        // just in case excess unset past binFactor
                        for (int jj = stop2; jj < h1; jj++) {
                            pOut.add(ii, jj);
                        }
                    }
                }
                if (i == (w0 - 1)) {
                    // just in case excess unset past binFastor
                    for (int ii = stop1; ii < w1; ii++) {
                        int stop2 = ((j + 1)*binFactor);
                        if (stop2 > (h1 - 1)) {
                            stop2 = h1;
                        }
                        for (int jj = (j*binFactor); jj < stop2; jj++) {
                            pOut.add(ii, jj);
                        }
                        if (j == (h0 - 1)) {
                            // just in case excess unset
                            for (int jj = stop2; jj < h1; jj++) {
                                pOut.add(ii, jj);
                            }
                        }
                    }
                }
            }
        }
        
        return output;
    }

    /**
     * create a two-dimensional float array of the img multiplied by
     * factor, but returned in row-major format [row][col].
     * @param img
     * @param factor
     * @return 
     */
    public float[][] multiply(GreyscaleImage img, float factor) {
     
        int nRows = img.getHeight();
        int nCols = img.getWidth();
        float[][] out = new float[nRows][nCols];
        for (int j = 0; j < nRows; ++j) {
            out[j] = new float[nCols];
        }

        for (int j = 0; j < nRows; ++j) {
            for (int i = 0; i < nCols; ++i) {
                out[j][i] = (float)img.getValue(i, j) * factor;
            }
        }

        return out;
    }

    public void createFirstDerivKeyPoints(float[][] image, 
        float sigma, TIntList outKeypoints0, TIntList outKeypoints1) {

        boolean createCurvatureComponents = false;
        
        StructureTensor tensorComponents = new StructureTensor(image, 
            sigma, createCurvatureComponents);
       
        float hLimit = 0.09f;//0.05f;
        
        createFirstDerivKeyPoints(tensorComponents, outKeypoints0, 
            outKeypoints1, hLimit);
    }
    
    public void createFirstDerivKeyPoints(
        StructureTensor tensorComponents, TIntList outKeypoints0, 
        TIntList outKeypoints1, float hLimit) {
       
        TIntList kp0 = new TIntArrayList();
        TIntList kp1 = new TIntArrayList();

        // square of 1st deriv:
        float[][] firstDeriv = add(tensorComponents.getDXSquared(), 
            tensorComponents.getDYSquared());

        peakLocalMax(firstDeriv, 1, 0.1f, kp0, kp1);

        float[][] detA = tensorComponents.getDeterminant();
        float[][] traceA = tensorComponents.getTrace();
        
        //float min = MiscMath.findMin(secondDeriv);
        //float max = MiscMath.findMax(secondDeriv);
        //System.out.println("min=" + min + " max=" + max);
        //System.out.println("nRows=" + nRows + " nCols=" + nCols);
        for (int i = 0; i < kp0.size(); ++i) {
            int x = kp0.get(i);
            int y = kp1.get(i);
            // harmonic mean, Brown, Szeliski, and Winder (2005),
            float hMean = (detA[x][y]/traceA[x][y]);
            if (hMean > hLimit) {
                //System.out.println(String.format("(%d,%d) detA/tr=%.4f",
                //    x, y, hMean));
                outKeypoints0.add(x);
                outKeypoints1.add(y);
            }
        }
        
        /*{// DEBUG
            float[][] a = copy(firstDeriv);
            MiscMath.applyRescale(a, 0, 255);
            MiscDebug.writeImage(a, "_fitsr_deriv_" 
                + MiscDebug.getCurrentTimeFormatted());
        }*/
        
    }

    private TFloatList calculateImageScales(List<GreyscaleImage> images) {

        TFloatList scales = new TFloatArrayList();
        scales.add(1);
        
        float prevScl = 1;
        
        for (int i = 1; i < images.size(); ++i) {
            float w = images.get(i).getWidth();
            float h = images.get(i).getHeight();
        
            float x = (float)images.get(i - 1).getWidth()/w;
            float y = (float)images.get(i - 1).getHeight()/h;
            
            float scale = prevScl * (x + y)/2.f;
            
            scales.add(scale);
            
            prevScl = scale;
        }
        
        return scales;
    }

    public Set<PairInt> binPoints(Set<PairInt> points, int binFactor) {
    
        Set<PairInt> out = new HashSet<PairInt>();
        for (PairInt p : points) {
            int x = Math.round((float)p.getX()/(float)binFactor);
            int y = Math.round((float)p.getY()/(float)binFactor);
            out.add(new PairInt(x, y));
        }
        
        return out;
    }

    public static class Colors {
        private final float[] colors;
        public Colors(float[] theColors) {
            colors = theColors;
        }
        public float[] getColors() {
            return colors;
        }
    }

    public Colors calculateAverageLAB(ImageExt input, Set<PairInt> points) {

        double labA = 0;
        double labB = 0;
        double labL = 0;

        for (PairInt p : points) {
            float[] lab = input.getCIELAB(p.getX(), p.getY());
            labL += lab[0];
            labA += lab[1];
            labB += lab[2];
        }
        labL /= (double)points.size();
        labA /= (double)points.size();
        labB /= (double)points.size();

        float[] labAvg = new float[]{(float)labL, (float)labA, (float)labB};

        Colors c = new Colors(labAvg);

        return c;
    }
    
    /**
     * calculate CIE LAB mean and standard deviation of
     * the mean.
     * @param input
     * @param points
     * @return 
     */
    public double[] calculateCIELABStats(ImageExt input, 
        Set<PairInt> points) {

        int n = points.size();
        double[] labA = new double[n];
        double[] labB = new double[n];
        double[] labL = new double[n];

        int count = 0;
        for (PairInt p : points) {
            float[] lab = input.getCIELAB(p.getX(), p.getY());
            labL[count] = lab[0];
            labA[count] = lab[1];
            labB[count] = lab[2];
            count++;
        }
        
        double[] lAS = MiscMath.getAvgAndStDev(labL);
        double[] aAS = MiscMath.getAvgAndStDev(labA);
        double[] bAS = MiscMath.getAvgAndStDev(labB);        

        double[] out = new double[6];
        out[0] = lAS[0];
        out[1] = lAS[1];
        out[2] = aAS[0];
        out[3] = aAS[1];
        out[4] = bAS[0];
        out[5] = bAS[1];
        
        return out;
    }

    public int calculateAverageHueAngle(ImageExt input, Set<PairInt> points) {

        double hueAngle = 0;

        for (PairInt p : points) {

            float[] lab = input.getCIELAB(p.getX(), p.getY());

            double ha;
            if (lab[1] == 0) {
                ha = 0;
            } else {
                ha = (Math.atan2(lab[2], lab[1]) * 180. / Math.PI);
                if (ha < 0) {
                    ha += 360.;
                }
            }
            hueAngle += ha;
        }
        hueAngle /= (double)points.size();

        return (int)Math.round(hueAngle);
    }

    public Colors calculateAverageRGB(ImageExt input, Set<PairInt> points) {

        float r = 0;
        float g = 0;
        float b = 0;

        for (PairInt p : points) {
            int idx = input.getInternalIndex(p.getX(), p.getY());
            r += input.getR(idx);
            g += input.getG(idx);
            b += input.getB(idx);
        }
        r /= (float)points.size();
        g /= (float)points.size();
        b /= (float)points.size();

        float[] rgbAvg = new float[]{r, g, b};

        Colors c = new Colors(rgbAvg);

        return c;
    }

    /**
     * opposite to shift zero-frequency component to the center of the spectrum
     * in that it shifts the zero-frequency component to the smallest indexes
     * in the arrays.
     *
     * adapted from
     * https://github.com/numpy/numpy/blob/master/LICENSE.txt
     * which has copyright
     *
     * Copyright (c) 2005-2016, NumPy Developers.
        All rights reserved.

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are
        met:

            * Redistributions of source code must retain the above copyright
               notice, this list of conditions and the following disclaimer.

            * Redistributions in binary form must reproduce the above
               copyright notice, this list of conditions and the following
               disclaimer in the documentation and/or other materials provided
               with the distribution.

            * Neither the name of the NumPy Developers nor the names of any
               contributors may be used to endorse or promote products derived
               from this software without specific prior written permission.

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
        "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
        LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
        A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
        OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
        SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
        LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
        DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
        THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
        (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
        OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
     * @param a
     * @return
     */
    public double[][] ifftShift(double[][] a) {

        double[][] b = new double[a.length][];

        // ---- reorder columns ----
        int nc = a.length;
        int p2 = nc - ((nc + 1)/2);
        List<Integer> range = new ArrayList<Integer>();
        for (int i = p2; i < nc; ++i) {
            range.add(Integer.valueOf(i));
        }
        for (int i = 0; i < p2; ++i) {
            range.add(Integer.valueOf(i));
        }
        int count = 0;
        for (Integer index : range) {
            int idx = index.intValue();
            b[count] = Arrays.copyOf(a[idx], a[idx].length);
            count++;
        }

        // ---- reorder rows ------
        nc = a[0].length;
        p2 = nc - ((nc + 1)/2);
        range = new ArrayList<Integer>();
        for (int i = p2; i < nc; ++i) {
            range.add(Integer.valueOf(i));
        }
        for (int i = 0; i < p2; ++i) {
            range.add(Integer.valueOf(i));
        }

        double[][] c = new double[a.length][];
        for (int i = 0; i < c.length; ++i) {
            c[i] = new double[a[0].length];
        }
        count = 0;
        for (Integer index : range) {
            int j = index.intValue();
            for (int col = 0; col < a.length; ++col) {
                c[col][count] = b[col][j];
            }
            count++;
        }

        return c;
    }

    /**
     * opposite to shift zero-frequency component to the center of the spectrum
     * in that it shifts the zero-frequency component to the smallest indexes
     * in the arrays.
     *
     * adapted from
     * https://github.com/numpy/numpy/blob/master/LICENSE.txt
     * which has copyright
     *
     * Copyright (c) 2005-2016, NumPy Developers.
        All rights reserved.

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are
        met:

            * Redistributions of source code must retain the above copyright
               notice, this list of conditions and the following disclaimer.

            * Redistributions in binary form must reproduce the above
               copyright notice, this list of conditions and the following
               disclaimer in the documentation and/or other materials provided
               with the distribution.

            * Neither the name of the NumPy Developers nor the names of any
               contributors may be used to endorse or promote products derived
               from this software without specific prior written permission.

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
        "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
        LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
        A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
        OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
        SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
        LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
        DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
        THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
        (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
        OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
     * @param a
     * @return
     */
    public Complex[][] ifftShift(Complex[][] a) {

        Complex[][] b = new Complex[a.length][];

        // ---- reorder dimension 0 ----
        int n0 = a.length;
        int p2 = n0 - ((n0 + 1)/2);
        List<Integer> range = new ArrayList<Integer>();
        for (int i = p2; i < n0; ++i) {
            range.add(Integer.valueOf(i));
        }
        for (int i = 0; i < p2; ++i) {
            range.add(Integer.valueOf(i));
        }
        int count = 0;
        for (Integer index : range) {
            int idx = index.intValue();
            b[count] = Arrays.copyOf(a[idx], a[idx].length);
            count++;
        }

        // ---- reorder dimension 1 ------
        int n1 = a[0].length;
        p2 = n1 - ((n1 + 1)/2);
        range = new ArrayList<Integer>();
        for (int i = p2; i < n1; ++i) {
            range.add(Integer.valueOf(i));
        }
        for (int i = 0; i < p2; ++i) {
            range.add(Integer.valueOf(i));
        }

        Complex[][] c = new Complex[a.length][];
        for (int i = 0; i < c.length; ++i) {
            c[i] = new Complex[a[0].length];
        }
        count = 0;
        for (Integer index : range) {
            int j = index.intValue();
            for (int col = 0; col < a.length; ++col) {
                c[col][count] = b[col][j];
            }
            count++;
        }

        return c;
    }
    
    public double[][] threshold(double[][] a, double threshold) {
        
        double[][] output = new double[a.length][];
        for (int i = 0; i < output.length; ++i) {
            output[i] = Arrays.copyOf(a[i], a[i].length);
            for (int j = 0; j < output[i].length; ++j) {
                double v = output[i][j];
                if (v < threshold) {
                    output[i][j] = 0;
                }
            }
        }
        
        return output;
    }
    
    /**
     * 
     * @param img
     * @return 
     */
    public GreyscaleImage[] createLabAandB(ImageExt img) {
        
        int w = img.getWidth();
        int h = img.getHeight();

        long t0 = System.currentTimeMillis();

        float[] labA = new float[w * h];
        float[] labB = new float[w * h];
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            float[] lab = img.getCIELAB(i);
            labA[i] = lab[1];
            labB[i] = lab[2];
        }
        
        labA = MiscMath.rescale(labA, 0, 255);
        labB = MiscMath.rescale(labB, 0, 255);
        
        GreyscaleImage labAImg = new GreyscaleImage(w, h,
            GreyscaleImage.Type.Bits32FullRangeInt);
        
        GreyscaleImage labBImg = new GreyscaleImage(w, h,
            GreyscaleImage.Type.Bits32FullRangeInt);
        
        for (int i = 0; i < labA.length; ++i) {
            labAImg.setValue(i, Math.round(labA[i]));
            labBImg.setValue(i, Math.round(labB[i]));
        }

        return new GreyscaleImage[]{labAImg, labBImg};
    }
    
    /**
     * create a color difference image equivalent to grey - o1,
     * r + g + b - (r - g) = 2*g + b
     * @param img
     * @return 
     */
    public GreyscaleImage createGPlusB(ImageExt img) {
        
        // consider scaling them by a set range of minimum and maximum possible
        // instead of minimum and maximum present
        // r + g + b - (r - g)--> 2*g + b
        int w = img.getWidth();
        int h = img.getHeight();

        float[] values = new float[w * h];
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            int g = img.getG(i);
            int b = img.getB(i);
            values[i] = 2*g + b;
        }
        
        GreyscaleImage out = new GreyscaleImage(w, h,
            GreyscaleImage.Type.Bits32FullRangeInt);
        
        values = MiscMath.rescale(values, 0, 255);
            
        for (int i = 0; i < values.length; ++i) {
            out.setValue(i, Math.round(values[i]));
        }

        return out;
    }
    
    /**
     * 
     * @param img
     * @return 
     */
    public GreyscaleImage createO1(ImageExt img) {
        
        int w = img.getWidth();
        int h = img.getHeight();

        float[] o1 = new float[w * h];
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            int r = img.getR(i);
            int g = img.getG(i);
            o1[i] = (r - g);
        }
        
        GreyscaleImage o1Img = new GreyscaleImage(w, h,
            GreyscaleImage.Type.Bits32FullRangeInt);
        
        o1 = MiscMath.rescale(o1, 0, 255);
            
        for (int i = 0; i < o1.length; ++i) {
            o1Img.setValue(i, Math.round(o1[i]));
        }

        return o1Img;
    }
    
    /**
     * 
     * @param img
     * @return 
     */
    public GreyscaleImage createGMinusB(ImageExt img) {
        
        int w = img.getWidth();
        int h = img.getHeight();

        float[] gb = new float[w * h];
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            int b = img.getB(i);
            int g = img.getG(i);
            gb[i] = (g - b);
        }
        
        GreyscaleImage gbImg = new GreyscaleImage(w, h,
            GreyscaleImage.Type.Bits32FullRangeInt);
        
        gb = MiscMath.rescale(gb, 0, 255);
            
        for (int i = 0; i < gb.length; ++i) {
            gbImg.setValue(i, Math.round(gb[i]));
        }

        return gbImg;
    }
    
    public void apply2LayerFilterOtsu(GreyscaleImage input) {
     
        apply2LayerFilterOtsu(input, 0.75f, 2.f);
    }
    
    public void apply2LayerFilterOtsu(GreyscaleImage input, float otsuFactor,
        float lowToHighFactor) {
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        if (w < 3 || h < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }
    
        OtsuThresholding ot = new OtsuThresholding();
            
        float tHigh = tHigh = otsuFactor * ot.calculateBinaryThreshold256(input);
        float tLow = tHigh/lowToHighFactor;
            
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        int n = input.getNPixels();
        
        GreyscaleImage img2 = input.createWithDimensions();
        
        for (int i = 0; i < img2.getNPixels(); ++i) {
            
            int v = input.getValue(i);
            
            if (v < tLow) {
                continue;
            } else if (v > tHigh) {
                img2.setValue(i, v);
                continue;
            }
            
            int x = input.getCol(i);
            int y = input.getRow(i);
            
            boolean foundHigh = false;
            boolean foundMid = false;
            
            for (int k = 0; k < dxs.length; ++k) {                
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                    continue;
                }
                int v2 = input.getValue(x2, y2);
                if (v2 > tHigh) {
                    foundHigh = true;
                    break;
                } else if (v2 > tLow) {
                    foundMid = true;
                }
            }
            if (foundHigh) {
                img2.setValue(i, v);
                continue;
            }
            if (!foundMid) {
                continue;
            }
            // search the 5 by 5 region for a "sure edge" pixel
            for (int dx = -2; dx <= 2; ++dx) {
                int x2 = x + dx;
                if ((x2 < 0) || (x2 > (w - 1))) {
                    continue;
                }
                for (int dy = -2; dy <= 2; ++dy) {
                    int y2 = y + dy;
                    if ((y2 < 0) || (y2 > (h - 1))) {
                        continue;
                    }
                    if (x2 == x && y2 == y) {
                        continue;
                    }
                    int v2 = input.getValue(x2, y2);
                    if (v2 > tHigh) {
                        img2.setValue(i, v);
                        foundHigh = true;
                        break;
                    }
                }
                if (foundHigh) {
                    break;
                }
            }
        }
        
        input.resetTo(img2);
        
        // apply post thinning corrections?
    }

    public int[] getAverageRGB(Image img, PairIntArray pArr) {
    
        if (pArr.getN() == 0) {
            return null;
        }
        
        int rSum = 0;
        int gSum = 0;
        int bSum = 0;
        for (int i = 0; i < pArr.getN(); ++i) {
            int x = pArr.getX(i);
            int y = pArr.getY(i);
            rSum += img.getR(x, y);
            gSum += img.getG(x, y);
            bSum += img.getB(x, y);
        }
        rSum /= pArr.getN();
        gSum /= pArr.getN();
        bSum /= pArr.getN();
        
        return new int[]{rSum, gSum, bSum};
    }
    
    public int[] getAverageRGB(Image img, Collection<PairInt> pArr) {
    
        if (pArr.isEmpty()) {
            return null;
        }
        
        int rSum = 0;
        int gSum = 0;
        int bSum = 0;
        for (PairInt p : pArr) {
            int x = p.getX();
            int y = p.getY();
            rSum += img.getR(x, y);
            gSum += img.getG(x, y);
            bSum += img.getB(x, y);
        }
        rSum /= pArr.size();
        gSum /= pArr.size();
        bSum /= pArr.size();
        
        return new int[]{rSum, gSum, bSum};
    }
    
    /**
     * NOTE: needs testing...invoker should trim for image bounds
     * where needed.
     * @param points
     * @param sigma
     */
    public void blur(Set<PairInt> points, SIGMA sigma) {
        
        // gaussian smoothing by sigma
        float[] kernel = Gaussian1D.getKernel(sigma);
        
        applyKernel(points, kernel);
    }
    
    /**
     * NOTE: needs testing...invoker should trim for image bounds
     * where needed.
     * @param points
     * @param sigma
     */
    public void applyKernel(Set<PairInt> points, float[] kernel) {
        
        int[] minMaxXY = MiscMath.findMinMaxXY(points);
        
        int xLL = minMaxXY[0] - 10;
        if (xLL < 0) {
            xLL = 0;
        }
        int yLL = minMaxXY[2] - 10;
        if (yLL < 0) {
            yLL = 0;
        }
        
        int xUR = minMaxXY[1] + 10;
        int yUR = minMaxXY[3] + 10;
        
        int w = xUR - xLL + 1;
        int h = yUR - yLL + 1;
        
        int[][] img = new int[w][];
        for (int i = 0; i < w; ++i) {
            img[i] = new int[h];
        }
        
        for (PairInt p : points) {
            int x = p.getX() - xLL;
            int y = p.getY() - yLL;
            img[x][y] = 126;
        }
        
        applyKernelTwo1Ds(img, kernel);
        
        points.clear();
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                if (img[i][j] <= 0) {
                    continue;
                }
                int x = i + xLL;
                int y = j + yLL;
                PairInt p = new PairInt(x, y);
                points.add(p);
            }
        }
    }
    
    /**
     * blur the points by sigma and trim any extending beyond image bounds.
     * @param points
     * @param sigma
     * @param imgWidth
     * @param imgHeight 
     */
    public void blurAndTrim(Set<PairInt> points, SIGMA sigma, int imgWidth, 
        int imgHeight) {
        
        blur(points, sigma);
        
        // trim any points extending beyond image bounds
        Set<PairInt> rm = new HashSet<PairInt>();
        for (PairInt p : points) {
            if ((p.getX() > (imgWidth - 1)) || (p.getY() > (imgHeight - 1)) ||
                (p.getX() < 0) || (p.getY() < 0)) {
                rm.add(p);
            }
        }
        
        points.removeAll(rm);
    }
    
    /**
     * NOTE: modifies input by the blur step.
     * @param contiguousPoints
     * @param sigma
     * @param imgWidth
     * @param imgHeight
     * @return 
     */
    public PairIntArray extractSmoothedOrderedBoundary(
        Set<PairInt> contiguousPoints, SIGMA sigma, int imgWidth, int imgHeight) {
                
        blurAndTrim(contiguousPoints, sigma, imgWidth, imgHeight);
        
        PerimeterFinder2 finder = new PerimeterFinder2();
        PairIntArray ordered = finder.extractOrderedBorder(
            contiguousPoints);
    
        return ordered;
    }
    
    /**
     * remove links from the adjacency map in the colorspace for pairs with
     * smaller histogram similarity than threshold.
     * current impl is using hsv color histograms.
     * @param img
     * @param listOfSets
     * @param adjMap
     * @param threshold 
     */
    public void filterAdjacencyMap(ImageExt img, List<Set<PairInt>> listOfSets,
        TIntObjectMap<TIntSet> adjMap, float threshold) {
        
        int n = listOfSets.size();
        
        ColorHistogram clrHist = new ColorHistogram();
        
        int[][][] hsvH = new int[n][][];
        //int[][][] cielabH = new int[n][][];
        
        List<PairInt> centroids = new ArrayList<PairInt>();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        for (int i = 0; i < n; ++i) {
            Set<PairInt> set = listOfSets.get(i);
            hsvH[i] = clrHist.histogramHSV(img, set);
            //cielabH[i] = clrHist.histogramCIELAB(img, set);
            centroids.add(curveHelper.calculateXYCentroids2(set));
        }
        
        Set<PairInt> rm = new HashSet<PairInt>();
        
        Set<PairInt> visited = new HashSet<PairInt>();
        
        TIntObjectIterator<TIntSet> iter = adjMap.iterator();
        for (int i = 0; i < adjMap.size(); ++i) {
            
            iter.advance();
            
            TIntSet set = iter.value();
            int idx1 = iter.key();
                        
            TIntIterator iter2 = set.iterator();
            while (iter2.hasNext()) {
                int idx2 = iter2.next();
                PairInt p = null;
                if (idx1 < idx2) {
                    p = new PairInt(idx1, idx2);
                } else {
                    p = new PairInt(idx2, idx1);
                }
                if (rm.contains(p) || visited.contains(p)) {
                    continue;
                }
                visited.add(p);
                
                int[][] hsv1 = hsvH[p.getX()];
                int[][] hsv2 = hsvH[p.getY()];
                
                //int[][] cie1 = cielabH[p.getX()];
                //int[][] cie2 = cielabH[p.getY()];
                
                float hsvInter = clrHist.intersection(hsv1, hsv2);
                
                //float cieInter = clrHist.intersection(cie1, cie2);
                
                System.out.println(
                    "cen1=" + centroids.get(idx1)
                    + " cen2=" + centroids.get(idx2)
                    + " hsvInter=" + hsvInter
                    //+ " cieInter=" + cieInter
                );
                
                if (hsvInter < threshold) {
                    rm.add(p);
                }
            }
        }

        for (PairInt r : rm) {
            TIntSet set1 = adjMap.get(r.getX());
            set1.remove(r.getY());

            set1 = adjMap.get(r.getY());
            set1.remove(r.getX());
        }
        
        System.out.println("adjacecny map filter removed " + rm.size());
    }
    
    /**
     * apply a dilate operator of size 3 x 3 to any pixel in image with value
     * greater than 0.
     * Note that if the img is not binary, the result may not be ideal because
     * no attempt has been made to account for existing pixel value when 
     * overwritten during dilation of adjacent pixel.
     * 
     * @param img 
     * @return  
     */
    public GreyscaleImage dilate(GreyscaleImage img) {
        return dilateOrErode(img, true);
    }
    
    /**
     * apply an erosion operator of size 3 x 3 to any pixel in image with value
     * 0.
     * 
     * @param img 
     * @return  
     */
    public GreyscaleImage erode(GreyscaleImage img) {
       return dilateOrErode(img, false);
    }
    
    /**
     * 
     * @param img 
     * @return  
     */
    private GreyscaleImage dilateOrErode(GreyscaleImage img, boolean dilate) {
       
        int w = img.getWidth();
        int h = img.getHeight();
        int n = img.getNPixels();
        
        GreyscaleImage out = img.copyImage();
        
        for (int pixIdx = 0; pixIdx < n; ++pixIdx) {
            int v = img.getValue(pixIdx);
            if (dilate && v < 1) {
                continue;
            } else if (!dilate && v != 0) {
                continue;
            }
            int x = img.getCol(pixIdx);
            int y = img.getRow(pixIdx);
            for (int i = -1; i <= 1; ++i) {
                int x2 = x + i;
                if (x2 < 0 || x2 > (w - 1)) {
                    continue;
                }
                for (int j = -1; j <= 1; ++j) {
                    if (i == 0 && j == 0) {
                        continue;
                    }
                    int y2 = y + j;
                    if (y2 < 0 || y2 > (h - 1)) {
                        continue;
                    }
                    out.setValue(x2, y2, v);
                }
            }
        }
        
        return out;
    }
    
    /**
     * 
     * @param img 
     * @return  
     */
    private Set<PairInt> dilate(Set<PairInt> set, 
        int imageWidth, int imageHeight, boolean dilate) {
       
        int w = imageWidth;
        int h = imageHeight;
        
        Set<PairInt> out = new HashSet<PairInt>();
        
        for (PairInt p : set) {
            int x = p.getX();
            int y = p.getY();
            for (int i = -1; i <= 1; ++i) {
                int x2 = x + i;
                if (x2 < 0 || x2 > (w - 1)) {
                    continue;
                }
                for (int j = -1; j <= 1; ++j) {
                    if (i == 0 && j == 0) {
                        continue;
                    }
                    int y2 = y + j;
                    if (y2 < 0 || y2 > (h - 1)) {
                        continue;
                    }
                    out.add(new PairInt(x2, y2));
                }
            }
        }
        
        return out;
    }
    
    /**
     * apply erode then dilate
     * 
     * @param img
     * @return 
     */
    public GreyscaleImage opening(GreyscaleImage img) {
        GreyscaleImage erode = erode(img);
        return dilate(erode);
    }
    
    /**
     * apply dilate then erode
     * 
     * @param img
     * @return 
     */
    public GreyscaleImage closing(GreyscaleImage img) {
        GreyscaleImage dilate = dilate(img);
        return erode(dilate);
    }
    
    /**
     * apply morphological thinning
     * @param img
     */
    public void applyThinning2(GreyscaleImage img) {

        int n0 = img.getWidth();
        int n1 = img.getHeight();
        
        int[][] morphInput = new int[n0][];
        for (int i = 0; i < n0; ++i) {
            morphInput[i] = new int[n1];
        }
        for (int i = 0; i < n0; ++i) {
            for (int j = 0; j < n1; ++j) {
                if (img.getValue(i, j) > 0) {                    
                    morphInput[i][j] = 1;
                } else {
                    morphInput[i][j] = 0;
                }
            }
        }        
        
        MorphologicalFilter mFilter = new MorphologicalFilter();
        int[][] skel = mFilter.bwMorphThin(morphInput, Integer.MAX_VALUE);

        Set<PairInt> points = new HashSet<PairInt>();
        
        for (int i = 0; i < n0; ++i) {
            for (int j = 0; j < n1; ++j) {
                int m = skel[i][j];
                int v = img.getValue(i, j) * m;
                img.setValue(i, j, v);
                if (v > 0) {
                    points.add(new PairInt(i, j));
                }
            }
        }
        
        PostLineThinnerCorrections pLTC = new PostLineThinnerCorrections();
        pLTC._correctForArtifacts(points, n0, n1);
        for (int i = 0; i < n0; ++i) {
            for (int j = 0; j < n1; ++j) {
                PairInt p = new PairInt(i, j);
                if (!points.contains(p)) {
                    img.setValue(i, j, 0);
                }
            }
        }
    }
    
    /**
     * apply 8 hit or miss filters iteratively until convergence to thin the
     * image.  the operation is performed on all pixels with value > 0.
     * @param img
     */
    public void applyThinning(GreyscaleImage img) {

        //from https://en.wikipedia.org/wiki/Hit-or-miss_transform
        // and thinning
        GreyscaleImage out = img.copyImage();

        // x,y pairs are sequential in these
        int[] c1 = new int[]{0, 0, -1, -1, 0, -1, 1, -1};
        int[] d1 = new int[]{-1, 1, 0, 1, 1, 1};
        int[] c2 = new int[]{-1, 0, 0, 0, -1, -1, 0, -1};
        int[] d2 = new int[]{0, 1, 1, 1, 1, 0};

        /*
        
            - - -        - -
              +        + + -
            + + +      + +
        
        */
        PairInt[][] neighborCoordOffsets
            = AbstractLineThinner.createCoordinatePointsForEightNeighbors(
            0, 0);

        int nEdited = 0;
        int nIter = 0;
        do {
            nEdited = 0;

            // test c1, d1 and it rotated by 90 3 times
            // test c2, d2 and it rotated by 90 3 times
            // need to alternate direction of approach
            for (int t = 0; t < 2; ++t) {
                int[] tmpC;
                int[] tmpD;
                if (t == 0) {
                    tmpC = Arrays.copyOf(c1, c1.length);
                    tmpD = Arrays.copyOf(d1, d1.length);
                } else {
                    tmpC = Arrays.copyOf(c2, c2.length);
                    tmpD = Arrays.copyOf(d2, d2.length);
                }
                for (int r = 0; r < 4; ++r) {
                    if (r > 0) {
                        rotatePairsBy90(tmpC);
                        rotatePairsBy90(tmpD);
                    }
                    
                    for (int x = 1; x < (img.getWidth() - 1); ++x) {
                        for (int y = 1; y < (img.getHeight() - 1); ++y) {
                            int v = img.getValue(x, y);
                            if (v == 0) {
                                continue;
                            }
                            if (allArePresent(img, x, y, tmpC)
                                && allAreNotPresent(img, x, y, tmpD)) {
                                if (!ImageSegmentation.doesDisconnect(out,
                                    neighborCoordOffsets, x, y)) {
                                    out.setValue(x, y, 0);
                                    nEdited++;
                                }
                            }
                        }
                    }
                    //MiscDebug.writeImage(out, "_thin_");
                }
                
                img.resetTo(out);
            }
            nIter++;
        } while (nEdited > 0);
                
        img.resetTo(out);
        
        Set<PairInt> points = readNonZeroPixels(img);
        
        PostLineThinnerCorrections pLTC = new PostLineThinnerCorrections();
        pLTC._correctForArtifacts(points, img.getWidth(), 
            img.getHeight());
        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                PairInt p = new PairInt(i, j);
                if (!points.contains(p)) {
                    img.setValue(i, j, 0);
                }
            }
        }
    }
    
    /**
     * apply 8 hit or miss filters iteratively until convergence to thin the
     * image.  the operation is performed on all pixels with value > 0.
     */
    public void applyThinning(Set<PairInt> points, int imageWidth, int imageHeight,
        boolean usePostCorrections) {

        //from https://en.wikipedia.org/wiki/Hit-or-miss_transform
        // and thinning
        Set<PairInt> out = new HashSet<PairInt>(points);

        // x,y pairs are sequential in these
        int[] c1 = new int[]{0, 0, -1, -1, 0, -1, 1, -1};
        int[] d1 = new int[]{-1, 1, 0, 1, 1, 1};
        int[] c2 = new int[]{-1, 0, 0, 0, -1, -1, 0, -1};
        int[] d2 = new int[]{0, 1, 1, 1, 1, 0};

        /*
        
            - - -        - -
              +        + + -
            + + +      + +
        
        */
        PairInt[][] neighborCoordOffsets
            = AbstractLineThinner.createCoordinatePointsForEightNeighbors(
            0, 0);

        int nEdited = 0;
        int nIter = 0;
        do {
            nEdited = 0;

            // test c1, d1 and it rotated by 90 3 times
            // test c2, d2 and it rotated by 90 3 times
            // need to alternate direction of approach
            for (int t = 0; t < 2; ++t) {
                int[] tmpC;
                int[] tmpD;
                if (t == 0) {
                    tmpC = Arrays.copyOf(c1, c1.length);
                    tmpD = Arrays.copyOf(d1, d1.length);
                } else {
                    tmpC = Arrays.copyOf(c2, c2.length);
                    tmpD = Arrays.copyOf(d2, d2.length);
                }
                for (int r = 0; r < 4; ++r) {
                    if (r > 0) {
                        rotatePairsBy90(tmpC);
                        rotatePairsBy90(tmpD);
                    }
                    
                    for (PairInt p : points) {
                        int x = p.getX();
                        int y = p.getY();
                        if (allArePresent(points, x, y, tmpC)
                            && allAreNotPresent(points, x, y, tmpD)) {
                            if (!ImageSegmentation.doesDisconnect(out,
                                neighborCoordOffsets, x, y, imageWidth, 
                                imageHeight)) {
                                out.remove(p);
                                nEdited++;
                            }
                        }
                    }
                    
                    //MiscDebug.writeImage(out, "_thin_");
                }
                
                points.clear();
                points.addAll(out);
            }
            nIter++;
        } while (nEdited > 0);
                
        points.clear();
        points.addAll(out);
     
        if (usePostCorrections) {
            PostLineThinnerCorrections pLTC = new PostLineThinnerCorrections();
            pLTC._correctForArtifacts(points, imageWidth, imageHeight);
        }
    }
    
    private void rotatePairsBy90(int[] xy) {
         
        /*
        int cos90 = 0;
        int sin90 = 1;
        scale = 1
        xc = yc = 0
        rotX = - (xc*scale + (((x0-xc)*scale*math.cos(theta))
                 + ((y0-yc)*scale*math.sin(theta)))
               = - (0 + (y0-0))
               = -y0

        rotY = - (yc*scale + ((-(x0-xc)*scale*math.sin(theta))
                 + ((y0-yc)*scale*math.cos(theta)))
               = - (0 + ((-(x0-0)))
               = +x0
        */
        
        for (int i = 0; i < xy.length; i += 2) {
            int x = xy[i];
            int y = xy[i + 1];
            xy[i] = -y;
            xy[i + 1] = x;
        }
    }
    
    private boolean allArePresent(GreyscaleImage img, int x, int y, int[] xy) {
        
        for (int k = 0; k < xy.length; k += 2) {
            int tx = x + xy[k];
            int ty = y + xy[k + 1];
            if (img.getValue(tx, ty) == 0) {
                return false;
            }
        }
        
        return true;
    }
    
    private boolean allArePresent(Set<PairInt> points, int x, int y, int[] xy) {
        
        for (int k = 0; k < xy.length; k += 2) {
            int tx = x + xy[k];
            int ty = y + xy[k + 1];
            if (!points.contains(new PairInt(tx, ty))) {
                return false;
            }
        }
        
        return true;
    }
    
    private boolean allAreNotPresent(GreyscaleImage img, int x, int y, int[] xy) {
        
        for (int k = 0; k < xy.length; k += 2) {
            int tx = x + xy[k];
            int ty = y + xy[k + 1];
            if (img.getValue(tx, ty) != 0) {
                return false;
            }
        }
        
        return true;
    }
    
    private boolean allAreNotPresent(Set<PairInt> points, int x, int y, int[] xy) {
        
        for (int k = 0; k < xy.length; k += 2) {
            int tx = x + xy[k];
            int ty = y + xy[k + 1];
            if (points.contains(new PairInt(tx, ty))) {
                return false;
            }
        }
        
        return true;
    }
    
    /**
     * create a two dimensional row-major format array of curvature of the
     * image img.  Note that img is expected to have all values >= 0.
     * Also note that sigma should be equal to or greater than
     * sqrt(2)/2.
     * @param img
     * @param sigma
     * @return 
     */
    public float[][] createCurvatureImage(GreyscaleImage img, float sigma) {
        
        /* curvature:
               dot is the degree of derivative...see ScaleSpaceCurvature
                      X_dot(t,o~) * Y_dot_dot(t,o~) - Y_dot(t,o~) * X_dot_dot(t,o~) 
            k(t,o~) = -------------------------------------------------------------
                                   (X_dot^2(t,o~) + Y_dot^2(t,o~))^1.5
        */
        
        float max = img.max();
        
        if (max <= 0) {
            throw new IllegalArgumentException("img values must be"
                + " >= 0 and maximum must be > 0");
        }
        
        // -- switch to row-major ----
        float[][] image = multiply(img, 1.f/max);
        
        assert(image.length == img.getHeight());
        assert(image[0].length == img.getWidth());
        
        return createCurvatureImage(image, sigma);
    }
    
    /**
     * create a two dimensional row-major format array of curvature of the
     * image img.  Note that img is expected to have all values >= 0.
     * Also note that sigma should be equal to or greater than
     * sqrt(2)/2.
     * @param image row-major formatted image
     * @param sigma
     * @return 
     */
    public float[][] createCurvatureImage(float[][] image, float sigma) {
        
        /* curvature:
               dot is the degree of derivative...see ScaleSpaceCurvature
                      X_dot(t,o~) * Y_dot_dot(t,o~) - Y_dot(t,o~) * X_dot_dot(t,o~) 
            k(t,o~) = -------------------------------------------------------------
                                   (X_dot^2(t,o~) + Y_dot^2(t,o~))^1.5
        */
        
        // --- create Sobel derivatives (gaussian 1st deriv sqrt(2)/2 = 0.707)----
        
        // switch X and Y sobel operations for row major

        TwoDFloatArray[] components =
            createCurvatureComponents(image, sigma);
        
        float[][] gX = components[0].a;
        float[][] gY = components[1].a;
        float[][] gX2 = components[2].a;
        float[][] gY2 = components[3].a;
        
        int nRows = gX.length;
        int nCols = gX[0].length;
        
        float[][] curvature = copy(gX);
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                
                float gX2gX2 = gX2[i][j] * gX2[i][j];
                float gY2gY2 = gY2[i][j] * gY2[i][j];
                if (gX2gX2 == 0 && gY2gY2 == 0) {
                    curvature[i][j] = Float.MAX_VALUE;
                    continue;
                }
                //(dx * dy(dy) - dy * dx(dx)) / (dx(dx)*dx(dx) + dy(dy)*dy(dy))
                curvature[i][j] = (float)(
                    (gX[i][j] * gY2[i][j] - gY[i][j] * gX2[i][j])
                    / Math.pow((gX2gX2 + gY2gY2), 1.5));
            }
        }
        
        return curvature;
    }
    
    
    /**
     * create an image segmented by curvature zero-crossings.
     * calculates the curvature in O(N_pixels) but using transcendental
     * operations, and then sets the output to 255 where curvature is
     * greater than 0 and sets output to 0 where curvature is less than 0.
     * the "zeroHandling" parameter determines how curvature exact values
     * of 0 are handled.
     * for zeroHandling = 0, curvature equal to zero sets output to 1;
     * for zeroHandling = 1, curvature equal to zero sets output to 0;
     * for zeroHandling = 2, curvature equal to zero sets output to 255.
     * 
     * This method is used in MedialAxis1.java.
     * 
     * @param image image in row-major format.
     * @param sigma (note, the internal first derivative, first step is
     * sigma=sqrt(2)/2 so this given sigma should be that or larger.
     * @param zeroHandling allowed values are 0, 1, or 2.
     * 0 = set curvature exact zero values to value 1,
     * 1 = set curvature exact zero values to value 0,
     * 2 = set curvature exact zero values to value 255.
     * @return 
     */
    public GreyscaleImage createZeroCrossingsCurvature(float[][] image,
        float sigma, int zeroHandling) {
                
        if ((zeroHandling < 0) || (zeroHandling > 2)) {
            throw new IllegalArgumentException("zeroHandling must be "
                + "between 0 and 2, inclusive.");
        }
        
        // -- switch to row-major until output ----
        float[][] curvature = createCurvatureImage(image, sigma);
        
        return createZeroCrossingsCurvatureImage(curvature, zeroHandling);
    }
    
    /**
     * create an image segmented by curvature zero-crossings.
     * calculates the curvature in O(N_pixels) but using transcendental
     * operations, and then sets the output to 255 where curvature is
     * greater than 0 and sets output to 0 where curvature is less than 0.
     * the "zeroHandling" parameter determines how curvature exact values
     * of 0 are handled.
     * for zeroHandling = 0, curvature equal to zero sets output to 1;
     * for zeroHandling = 1, curvature equal to zero sets output to 0;
     * for zeroHandling = 2, curvature equal to zero sets output to 255.
     * 
     * This method is used in MedialAxis1.java.
     * 
     * @param img
     * @param sigma (note, the internal first derivative, first step is
     * sigma=sqrt(2)/2 so this given sigma should be that or larger.
     * @param zeroHandling allowed values are 0, 1, or 2.
     * 0 = set curvature exact zero values to value 1,
     * 1 = set curvature exact zero values to value 0,
     * 2 = set curvature exact zero values to value 255.
     * @return 
     */
    public GreyscaleImage createZeroCrossingsCurvature(GreyscaleImage img,
        float sigma, int zeroHandling) {
                
        if ((zeroHandling < 0) || (zeroHandling > 2)) {
            throw new IllegalArgumentException("zeroHandling must be "
                + "between 0 and 2, inclusive.");
        }
        
        // -- switch to row-major until output ----
        float[][] curvature = createCurvatureImage(img, sigma);
        
        return createZeroCrossingsCurvatureImage(curvature, zeroHandling);
    }
    
    /**
     * create an image segmented by curvature zero-crossings.
     * calculates the curvature in O(N_pixels) but using transcendental
     * operations, and then sets the output to 255 where curvature is
     * greater than 0 and sets output to 0 where curvature is less than 0.
     * the "zeroHandling" parameter determines how curvature exact values
     * of 0 are handled.
     * for zeroHandling = 0, curvature equal to zero sets output to 1;
     * for zeroHandling = 1, curvature equal to zero sets output to 0;
     * for zeroHandling = 2, curvature equal to zero sets output to 255.
     * 
     * This method is used in MedialAxis1.java.
     * 
     * @param curvature
     * @param zeroHandling allowed values are 0, 1, or 2.
     * 0 = set curvature exact zero values to value 0,
     * 1 = set curvature exact zero values to value 1,
     * 2 = set curvature exact zero values to value 255.
     * @return 
     */
    private GreyscaleImage createZeroCrossingsCurvatureImage(float[][] curvature,
        int zeroHandling) {
                
        if ((zeroHandling < 0) || (zeroHandling > 2)) {
            throw new IllegalArgumentException("zeroHandling must be "
                + "between 0 and 2, inclusive.");
        }
        
        int nRows = curvature.length;
        int nCols = curvature[0].length;
        
        GreyscaleImage out = new GreyscaleImage(nCols, nRows);
        
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                float v = curvature[i][j];
                if (v > 0) {
                    out.setValue(j, i, 255);
                } else if (v == 0) {
                    switch(zeroHandling) {
                        case 0:
                            out.setValue(j, i, 0);
                            break;
                        case 1:
                            out.setValue(j, i, 1);
                            break;
                        default:
                            out.setValue(j, i, 255);
                            break;
                    }
                } else {
                    out.setValue(j, i, 0);
                }
            }
        }
        
        return out;
    }
    
    /**
     * create texture transforms from 
     * "Textured Image Segmentation" by Laws, 1980.
     * 
     * The transforms are combinations of filters based on
     *  L5 level = [1 4 6 4 1]     
             gaussian, binomial for sigma=1 
             B3 spline function, used in ATrous wavelet
        E5 edge  = [-1 -2 0 2 1]   
             1st deriv of gaussian, binomial for sigma=1
        S5 spot =   [-1 0 2 0 -1]  
             -1 times 2nd deriv binomial for sigma=sqrt(2)/2,... LOG
        R5 ripple = [1 -4 6 -4 1]  
              3rd deriv gaussian, ...Gabor
     
     NOTE: bright clumps in R5 R5 look useful for finding vegetation.
        It finds the bounds of the vegetation... places where the 
        the change of the change of the gradient is large (and dense).
        can apply adaptive means to the feature image to find the 
        brightest of these.
        L5 S5 looks useful for finding horizontal lines such as edge segments
        of windows.
    
     * @param img 
       @param state 0=do not process derivatives further,
       1=subtract mean, 2=subtract mean and square to make variance,
       3=make zero mean, unit standard derivative, but multiplied by
       255 to put into integer range for result.
     * @return textureTransforms 
       GreyscaleImage[]{
       L5E5/E5L5, L5S5/S5L5, L5R5/R5L5, E5E5.
       E5S5/S5E5, E5R5/R5E5, S5S5, S5R5/R5S5,
       R5R5}
     */
    public Map<String, GreyscaleImage> createTextureTransforms(
        GreyscaleImage img, int state) {
        
        if (state < 0 || state > 3) {
            throw new IllegalArgumentException("state must be between"
                + "0 and 3, inclusive");
        }
        
        /*
        adapted from a cs lecture on texture filters from uw
        (https://courses.cs.washington.edu/courses/cse455/09wi/Lects/lect12.pdf
        which possibly uses:
        "Statistical Texture Analysis" by Srinivasan and Shobha 2008)
        
        Both contain content from "Textured Image Segmentation" by Laws, 1980.
        
        filters:
            L5 level = [1 4 6 4 1]     
                 gaussian, binomial for sigma=1 
                 B3 spline function, used in ATrous wavelet
            E5 edge  = [-1 -2 0 2 1]   
                 1st deriv of gaussian, binomial for sigma=1
            S5 spot =   [-1 0 2 0 -1]  
                 -1 times 2nd deriv binomial for sigma=sqrt(2)/2; a.k.a. LOG
            R5 ripple = [1 -4 6 -4 1]  
                  3rd deriv gaussian, a.k.a. Gabor
            W5 waves = [-1, 2, 0, -2, -1]

        - the 2D masks are created by multiplying the 1D masks to make a 5x5 matrix
              E5 X L5 = -1 -4 -6 -4 -1
                        -2 -8 -12 -8 -1
                         0  0  0  0   0
                         2  8  12  8  1
                         1  4  6   4  1
        - there are 9 feature vectors one could make.  can see that some compose
          tensors, and different keypoint algorithms.
        
        there are 9 feature vectors one could make.
          created by subtracting the mean neighborhood intensity from pixel
             filter the neighborhood with the 16 5 x 5 masks
             then energy at each pixel is summing abs value of filter output
                across neighbor region and storing result for the center pixel.
          The 9 features made from those 16 combinations of 4 filters are:
              L5L5, L5E5/E5L5, L5S5/S5L5, L5R5/R5L5, E5E5.
              E5S5/S5E5, E5R5/R5E5, S5S5, S5R5/R5S5,
              R5R5
        */
        
        float[] kernelL5 = new float[]{ 1,  4, 6, 4, 1};
        float[] kernelE5 = new float[]{-1, -2, 0, 2, 1};
        float[] kernelS5 = new float[]{-1,  0, 2, 0, -1};
        float[] kernelR5 = new float[]{ 1, -4, 6, -4, 1};
        float[][] kernels = new float[4][];
        kernels[0] = kernelL5;
        kernels[1] = kernelE5;
        kernels[2] = kernelS5;
        kernels[3] = kernelR5;
        String[] labels = new String[]{"L5", "E5", "S5", "R5"};
                
        Map<String, GreyscaleImage> transformed = new
            HashMap<String, GreyscaleImage>();

        for (int dir = 0; dir < 2; ++dir) {
            for (int l0 = 0; l0 < labels.length; ++l0) {
                for (int l1 = l0; l1 < labels.length; ++l1) {
                    int i = (dir == 0) ? l0 : l1;
                    int j = (dir == 0) ? l1 : l0;
                    float[] filter1 = kernels[i];
                    if ((dir == 1) && (i == j)) {
                        continue;
                    }
                    float[] filter2 = kernels[j];
                    GreyscaleImage img2 = img.copyToFullRangeIntImage();
                    applyKernel1D(img2, filter1, true);
                    applyKernel1D(img2, filter2, false);

                    if (i != j) {
                        GreyscaleImage img3 = img.copyToFullRangeIntImage();
                        applyKernel1D(img3, filter2, true);
                        applyKernel1D(img3, filter1, false);
                        img2 = divide(img2, img3);
                    }

                    GreyscaleImage imgM = null;

                    /*
                    0=do not process derivatives further,
                    1=subtract mean, 
                    2=subtract mean and square to make variance,
                    3=make zero mean, unit standard derivative, 
                      but multiplied by 255 to put into integer range for result.
                    */

                    if (state > 0) {
                        imgM = img2.copyToFullRangeIntImage();
                        applyCenteredMean2(imgM, 2);
                        img2 = subtractImages(img2, imgM);
                    }

                    if (state == 3) {
                        // make unit standard deviation image, but mult by 255
                        // for storage in integer.
                        // NOTE: considering change of output to float array for
                        // comparison to databases
                        for (int ii = 0; ii < img2.getNPixels(); ++ii) {
                            double m = imgM.getValue(ii);
                            double v = 255.*img2.getValue(ii)/(Math.sqrt(2)/m);
                            img2.setValue(ii, (int)Math.round(v));
                        }
                    } else if (state == 2) {
                        // square img2 to result in variance
                        for (int ii = 0; ii < img2.getNPixels(); ++ii) {
                            int v = img2.getValue(ii);
                            v *= v;
                            img2.setValue(ii, v);
                        }
                    }

                    String label = labels[i] + labels[j];

                    transformed.put(label, img2);

                    /*{
                        GreyscaleImage img3 = img2.copyImage();
                        MiscMath.rescale(img3, 0, 255);
                        MiscDebug.writeImage(img3, "_" + labels[i] + labels[j] + "_feature_");
                        applyAdaptiveMeanThresholding(img3, 2);
                        MiscDebug.writeImage(img3, "_" + labels[i] + labels[j] + "_feature_adap_means_");
                    }*/
                }
            }
        }
        
        return transformed;
    }
    
    public void exploreTextures() throws IOException {
        
        /*
        textures in frequency space:
        need to simplify the number of points contributing to the
        frequency domain pattern,
        so will calculate key points as sparse representation.
        */
        
        int maxDimension = 256;//512;
        
        String fileName1 = "android_statues_02.jpg";
        fileName1 = "merton_college_I_001.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath1);
        
        long ts = MiscDebug.getCurrentTimeFormatted();

        int w1 = img.getWidth();
        int h1 = img.getHeight();

        int binFactor1 = (int) Math.ceil(Math.max(
            (float) w1 / maxDimension,
            (float) h1 / maxDimension));

        img = binImage(img, binFactor1);
        
        float max = img.max();
        
        if (max <= 0) {
            throw new IllegalArgumentException("img values must be"
                + " >= 0 and maximum must be > 0");
        }
        
        int nRows = img.getHeight();
        int nCols = img.getWidth();
        
        // axis 0 coordinates
        TIntList keypoints0 = new TIntArrayList();
        
        // axis 1 coordinates
        TIntList keypoints1 = new TIntArrayList();
        
        /*ORB orb = new ORB(10000);
        orb.overrideToNotCreateDescriptors();
        orb.overrideToAlsoCreate1stDerivKeypoints();
        //orb.overrideToCreateCurvaturePoints();
        orb.detectAndExtract(img.copyToColorGreyscale());
        keypoints0.addAll(orb.getAllKeyPoints0());
        keypoints1.addAll(orb.getAllKeyPoints1());
        */
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        EdgeFilterProducts products = imageSegmentation.createGradient(
            img.copyToColorGreyscaleExt(), 2, ts);
                
        GreyscaleImage gradient = products.getGradientXY();
        
        float sigma = SIGMA.getValue(SIGMA.ZEROPOINTSEVENONE);
        
        // -- switch to row-major ----
        float[][] image = multiply(img, 1.f/max);
                
        //  strong high density responses in r5r5 for edges of vegetation.
        //  textures such as bricks or roof tiles are present in r5r5
        //  but so are strong edges, so possibly need
        //  to use the gradient edges here to distinguish between
        //  corner and the numerous points that are not good matching
        //  points.
        
        // thresh is usually 0.01f
        //createCurvatureKeyPoints(image, sigma, keypoints0, keypoints1, 
        //   0.001f);
        
        createR5R5KeyPoints(image, keypoints0, keypoints1);
        //createE5E5KeyPoints(image, keypoints0, keypoints1);
        //createL5E5KeyPoints(image, keypoints0, keypoints1);
        //createS5S5KeyPoints(image, keypoints0, keypoints1);        
        //createFirstDerivKeyPoints(image, sigma, keypoints0, keypoints1);
        
        //Image kpImg = img.copyToColorGreyscale();
        Image kpImg = gradient.copyToColorGreyscale();
        for (int i = 0; i < keypoints0.size(); ++i) {
            int x = keypoints1.get(i);
            int y = keypoints0.get(i);
            kpImg.setRGB(x, y, 255, 0, 0);
        }
        MiscDebug.writeImage(kpImg, "_keypoints_1_");
        
        /*
        Complex1D[] ccOut = create2DFFT2WithSwapMajor(kpImg, true);

        assert(nRows == ccOut[0].x.length);
        assert(nCols == ccOut.length);
        
        double[][] kpFreqR = new double[nCols][];
        double[][] kpFreqI = new double[nCols][];
        for (int i0 = 0; i0 < ccOut.length; ++i0) {
            kpFreqR[i0] = new double[nRows];
            kpFreqI[i0] = new double[nRows];
            for (int i1 = 0; i1 < nRows; ++i1) {
                kpFreqR[i0][i1] = ccOut[i0].x[i1];
                kpFreqI[i0][i1] = ccOut[i0].y[i1];
            }
        }
        
        TIntList plotRows = new TIntArrayList();
        TIntList plotCols = new TIntArrayList();
        plotRows.add(10);
        plotRows.add(50);
        plotRows.add(100);
        plotRows.add(110);
        
        plotCols.add(10);
        plotCols.add(50);
        plotCols.add(100);
        plotCols.add(150);
        plotCols.add(200);
        plotCols.add(210);
        
        String lbl = "keypoints_freq";
        
        MiscDebug.plot(kpFreqR, plotCols, plotRows, lbl);
       
        MiscMath.applyRescale(kpFreqR, 0, 255);
        Image kpFreqRImg = img.copyToColorGreyscale();
        for (int i = 0; i < nCols; ++i) {
            for (int j = 0; j < nRows; ++j) {
                kpFreqRImg.setRGB(i, j, (int)Math.round(kpFreqR[i][j]), 0, 0);
            }
        }
        MiscDebug.writeImage(kpFreqRImg, "_keypoints_freq_");
        */
        
   /*
        // ---- edited _keypoints_1_ image to keep only a characteristic section ---
        String filePath = ResourceFinder.findFileInTestResources(
            "vegetation_peak_texture.png");
        
        GreyscaleImage imgPattern = ImageIOHelper.readImage(
            filePath).copyToGreyscale();
        
        Complex[][] fftPattern = PhaseCongruencyDetector
            .createLowPassFreqDomainFilter(imgPattern);
       
        PeriodicFFT perfft2 = new PeriodicFFT();
        Complex[][][] perfResults = perfft2.perfft2(img, false);
        Complex[][] fftImage = perfResults[1];
        
        // --- image in the frequency domain convolved with texture patch ----
        Complex[][] freqDomainImageTimesPattern = 
            convolveWithKernel(fftImage, fftPattern);
        
        // ----- transform that to spatial domain ----
        Complex[][] fComplex = create2DFFT(freqDomainImageTimesPattern, false, false);    
        double[][] transformedReal = new double[nCols][];
        for (int i0 = 0; i0 < nCols; ++i0) {
            transformedReal[i0] = new double[nRows];
            for (int i1 = 0; i1 < nRows; ++i1) {
                transformedReal[i0][i1] = fComplex[i1][i0].abs();
            }
        }
        
        MiscMath.applyRescale(transformedReal, 0, 255);
        GreyscaleImage kpFreqR2Img = new GreyscaleImage(nCols, nRows);
        for (int i = 0; i < nCols; ++i) {
            for (int j = 0; j < nRows; ++j) {
                kpFreqR2Img.setValue(i, j, 
                    (int)Math.round(transformedReal[i][j]));
            }
        }
        MiscDebug.writeImage(kpFreqR2Img, "_keypoints_freq2_spatial_");
    */
    }
    
    /**
     * 
     * @param image
     * @param sigma
     * @param outputKeypoints0
     * @param outputKeypoints1
     * @param thresholdRel default is 0.01f
     */
    public void createCurvatureKeyPoints(float[][] image, float sigma,
        TIntList outputKeypoints0, TIntList outputKeypoints1,
        float thresholdRel) {
        
        boolean doCreateCurvatureKeyPoints = true;
        
        StructureTensor tensorComponents = new StructureTensor(image, 
            sigma, doCreateCurvatureKeyPoints);
        
        createCurvatureKeyPoints(tensorComponents, outputKeypoints0,
            outputKeypoints1, thresholdRel);
    }
    
    /**
     * 
     * @param tensorComponents
     * @param outputKeypoints0
     * @param outputKeypoints1
     * @param thresholdRel default is 0.01f
     */
    public void createCurvatureKeyPoints(StructureTensor tensorComponents,
        TIntList outputKeypoints0, TIntList outputKeypoints1,
        float thresholdRel) {
        
        // square of 1st deriv:
        //float[][] firstDeriv = add(tensorComponents.getDXSquared(), 
        //    tensorComponents.getDYSquared());
        //float max1stDeriv = MiscMath.findMax(firstDeriv);
        //float f = max1stDeriv / 10;

        float[][] dx = tensorComponents.getDX();
        float[][] dx2 = tensorComponents.getDDX();
        float[][] dy = tensorComponents.getDY();
        float[][] dy2 = tensorComponents.getDDY();
        float[][] curvature = copy(dx);
        
        int nRows = dx.length;
        int nCols = dx[0].length;
        
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                float dx2dx2 = dx2[i][j] * dx2[i][j];
                float dy2dy2 = dy2[i][j] * dy2[i][j];
                if (dx2dx2 == 0 && dy2dy2 == 0) {
                    curvature[i][j] = Float.MAX_VALUE;
                    continue;
                }
                //(dx * dy(dy) - dy * dx(dx)) / (dx(dx)*dx(dx) + dy(dy)*dy(dy))^1.5
                curvature[i][j] = (float)(
                    (dx[i][j] * dy2[i][j] - dy[i][j] * dx2[i][j])
                    / (dx2dx2 + dy2dy2));
                    /// Math.pow((dx2dx2 + dy2dy2), 1.5));
            }
        }

        peakLocalMax(curvature, 1, thresholdRel, outputKeypoints0, 
            outputKeypoints1);
        
        /*{// DEBUG
            float[][] a = copy(curvature);
            System.out.println("min=" + MiscMath.findMin(a) + " max=" + MiscMath.findMax(a));
            MiscMath.applyAbsoluteValue(a);
            MiscMath.applyRescale(a, 0, 255);
            MiscDebug.writeImage(a, "_curvature_" 
                + MiscDebug.getCurrentTimeFormatted());
        }*/
        
        /*MiscMath.applyRescale(curvature, 0, 255);
        GreyscaleImage kpImg = new GreyscaleImage(nCols, nRows);
        for (int i = 0; i < outputKeypoints0.size(); ++i) {
            int x = outputKeypoints1.get(i);
            int y = outputKeypoints0.get(i);
            kpImg.setValue(x, y, Math.round(curvature[y][x]));
        }
        MiscDebug.writeImage(kpImg, "_curvature_");
        */
    }
    
    public void createR5R5KeyPoints(float[][] image,
        TIntList outputKeypoints0, TIntList outputKeypoints1) {
        
        //3rd deriv gaussian, a.k.a. Gabor
        float[] kernel = new float[]{ 1, -4, 6, -4, 1};
        
        createLawKeyPoints(image, kernel, kernel, outputKeypoints0, 
            outputKeypoints1);
    }
    
    public void createS5S5KeyPoints(float[][] image,
        TIntList outputKeypoints0, TIntList outputKeypoints1) {
        
        //-1 times 2nd deriv binomial for sigma=sqrt(2)/2; a.k.a. LOG
        float[] kernel = new float[]{-1, 0, 2, 0, -1};
        
        createLawKeyPoints(image, kernel, kernel, outputKeypoints0, 
            outputKeypoints1);
    }
    
    public void createE5E5KeyPoints(float[][] image,
        TIntList outputKeypoints0, TIntList outputKeypoints1) {
        
        //1st deriv of gaussian, binomial for sigma=1
        float[] kernel = new float[]{-1, -2, 0, 2, 1};
        
        createLawKeyPoints(image, kernel, kernel, outputKeypoints0, 
            outputKeypoints1);
    }
    
    public void createL5E5KeyPoints(float[][] image,
        TIntList outputKeypoints0, TIntList outputKeypoints1) {
        
        //1st deriv of gaussian, binomial for sigma=1
        float[] kernelE5 = new float[]{-1, -2, 0, 2, 1};
        
        // gaussian, binomial for sigma=1; B3 spline used in ATrous wavelet
        float[] kernelL5 = new float[]{1, 4, 6, 4, 1};
        
        createLawKeyPoints(image, kernelL5, kernelE5, outputKeypoints0, 
            outputKeypoints1);
    }
    
    private void createLawKeyPoints(float[][] image, float[] kernel1,
        float[] kernel2,
        TIntList outputKeypoints0, TIntList outputKeypoints1) {
        
        float[][] image2 = copy(image);
        
        // row major, so need to use y operations for x and vice versa
        applyKernel1D(image2, kernel1, true);
        applyKernel1D(image2, kernel2, false);
        
        /*if (!Arrays.equals(kernel1, kernel2)) {
            float[][] image3 = copy(image2);
            applyKernel1D(image3, kernel1, false);
            applyKernel1D(image3, kernel2, true);
            image2 = divide(image2, image3);
        }*/
  
        // put float back into integer scale, 0 to 255
        MiscMath.applyRescale2(image2, 0, 255);
        
        int nCols = image2.length;
        int nRows = image2[0].length;
        
        GreyscaleImage imageM = new GreyscaleImage(nCols, nRows);
        for (int i = 0; i < nCols; ++i) {
            for (int j = 0; j < nRows; ++j) {
                imageM.setValue(i, j, Math.round(image2[i][j]));
            }
        }
        
        applyCenteredMean2(imageM, 2);
       
        // subtract mean
        for (int i = 0; i < nCols; ++i) {
            for (int j = 0; j < nRows; ++j) {
                image2[i][j] -= imageM.getValue(i, j);
            }
        }
                   
        // NOTE: square image2 to result in variance if preferred.
        for (int i = 0; i < image2.length; ++i) {
            for (int j = 0; j < image2[i].length; ++j) {
                float v = image2[i][j];
                //v *= v;
                image2[i][j] = v;
            }
        }
       
        GreyscaleImage img2 = imageM.createFullRangeIntWithDimensions();
        for (int i = 0; i < image2.length; ++i) {
            for (int j = 0; j < image2[i].length; ++j) {
                float v = image2[i][j];
                img2.setValue(i, j, Math.round(v));
            }
        }
            
        // use adaptive means to extract centers
        applyAdaptiveMeanThresholding(img2, 1);
        for (int i = 3; i < img2.getWidth(); ++i) {
            for (int j = 3; j < img2.getHeight(); ++j) {
                if (img2.getValue(i, j) == 0) {
                    outputKeypoints0.add(i);
                    outputKeypoints1.add(j);
                }
            }
        }

        // or, extract many points all over the image:
        //peakLocalMax(image2, 1, 0.01f, outputKeypoints0, 
        //    outputKeypoints1);
    }
    
    /**
     * return first derivatives from sigma=sqrt(2)/2 and
     * then second derivatives with sigma convolved from the
     * first derivatives.
     * @param image
     * @param sigma
     * @return [d/dx, d/dy, d/dx(dx), d/dy(dy)]
     */
    protected TwoDFloatArray[] createCurvatureComponents(
        float[][] image, float sigma) {

        // --- create Sobel derivatives (gaussian 1st deriv sqrt(2)/2 = 0.707)----
        
        // switch X and Y sobel operations for row major

        float[][] gX = copy(image);
        applySobelY(gX);

        float[][] gY = copy(image);
        applySobelX(gY);

        float[] kernel = Gaussian1D.getKernel(sigma);
       
        // for curvature, need d/dy(dy) and d/dx(dx)
        
        float[][] gX2 = copy(gX);
        float[][] gY2 = copy(gY);
        
        // row major, so need to use y operations for x and vice versa
        applyKernel1D(gX2, kernel, false);
        applyKernel1D(gY2, kernel, true);
        
        TwoDFloatArray[] components = new TwoDFloatArray[4];
        components[0] = new TwoDFloatArray(gX);
        components[1] = new TwoDFloatArray(gY);
        components[2] = new TwoDFloatArray(gX2);
        components[3] = new TwoDFloatArray(gY2);
        
        return components;
    }
    
    /**
     * get centers of high spatial and intensity variability using the
     * gaussian third derivative and adaptive means.
     *
     * @param img 
       
     * @return points of high density spatial variability and intensity.
     */
    public Set<PairInt> findHighDensityHighVariabilityPoints(GreyscaleImage img) {
        
        float[] kernelR5 = new float[]{ 1, -4, 6, -4, 1};
        
        GreyscaleImage img2 = img.copyToFullRangeIntImage();
        applyKernel1D(img2, kernelR5, true);
        applyKernel1D(img2, kernelR5, false);
                    
        GreyscaleImage imgM = img2.copyToFullRangeIntImage();
        applyCenteredMean2(imgM, 2);
        img2 = subtractImages(img2, imgM);
                    
        for (int ii = 0; ii < img2.getNPixels(); ++ii) {
            int v = img2.getValue(ii);
            v *= v;
            img2.setValue(ii, v);
        }
        
        applyAdaptiveMeanThresholding(img, 1);
        
        Set<PairInt> points = new HashSet<PairInt>();
        for (int ii = 0; ii < img2.getNPixels(); ++ii) {
            int v = img2.getValue(ii);
            if (v == 0) {
                points.add(new PairInt(img2.getCol(ii), img2.getRow(ii)));
            }
        }            
        
        return points;
    }
    
    public GreyscaleImage divide(GreyscaleImage img1, GreyscaleImage img2) {
        
        if (img1.getNPixels() != img2.getNPixels() || img1.getWidth() != 
            img2.getWidth() || img1.getHeight() != img2.getHeight()) {
            throw new IllegalArgumentException("img1 and img2 must be same size");
        }
        
        int n = img1.getNPixels();
        
        GreyscaleImage out = img1.copyToFullRangeIntImage();

        for (int i = 0; i < n; ++i) {
            int v = out.getValue(i);
            int v2 = img2.getValue(i);
            if (v2 > 0) {
                v /= v2;                
            }
            out.setValue(i, v);
        }
        
        return out;
    }
    
    public float[][] copy(float[][] a) {

        int n1 = a.length;
        int n2 = a[0].length;

        float[][] out = new float[n1][n2];
        for (int i = 0; i < n1; ++i) {
            out[i] = Arrays.copyOf(a[i], a[i].length);
        }

        return out;
    }
    
    public double[][] copyToDouble(float[][] a) {

        int n1 = a.length;
        int n2 = a[0].length;

        double[][] out = new double[n1][n2];
        for (int i = 0; i < n1; ++i) {
            out[i] = new double[n2];
            for (int j = 0; j < n2; ++j) {
                out[i][j] = a[i][j];
            }
        }

        return out;
    }

    public float[][] multiply(float[][] a, float[][] b) {

        float[][] c = copy(a);

        for (int i = 0; i < c.length; ++i) {
            for (int j = 0; j < c[0].length; ++j) {
                c[i][j] *= b[i][j];
            }
        }

        return c;
    }

    public float[][] divide(float[][] a, float[][] b) {

        float[][] c = copy(a);

        for (int i = 0; i < c.length; ++i) {
            for (int j = 0; j < c[0].length; ++j) {
                float v = b[i][j];
                if (v == 0) {
                    c[i][j] = Float.MAX_VALUE;
                } else {
                    c[i][j] /= v;
                }
            }
        }

        return c;
    }

    public float[][] add(float[][] a, float[][] b) {

        float[][] c = copy(a);

        for (int i = 0; i < c.length; ++i) {
            for (int j = 0; j < c[0].length; ++j) {
                c[i][j] += b[i][j];
            }
        }

        return c;
    }

    public float[][] subtract(float[][] a, float[][] b) {

        float[][] c = copy(a);

        for (int i = 0; i < c.length; ++i) {
            for (int j = 0; j < c[0].length; ++j) {
                c[i][j] -= b[i][j];
            }
        }

        return c;
    }
    
    /**
     * if image dimensions are a power of 2, returns null, else, increases the
     * dimensions to the power of 2 larger than each dimension.
     * the padding zeros will be at the end of the axes.
     * @param img
     * @return 
     */
    public GreyscaleImage increaseToPowerOf2(GreyscaleImage img) {
        
        int n0 = img.getWidth();
        int n1 = img.getHeight();

        int nn0 = 1 << (int)(Math.ceil(Math.log(n0)/Math.log(2)));
        int nn1 = 1 << (int)(Math.ceil(Math.log(n1)/Math.log(2)));

        if (nn0 == n0 && nn1 == n1) {
            return null;
        }
        
        GreyscaleImage img2 = new GreyscaleImage(nn0, nn1, 
            img.getType());
        for (int i = 0; i < n0; ++i) {
            for (int j = 0; j < n1; ++j) {
                img2.setValue(i, j, img.getValue(i, j));
            }
        }
        
        return img2;
    }
    
    /**
     * if image dimensions are a power of 2, returns null, else, increases the
     * dimensions to the power of 2 larger than each dimension.
     * the padding zeros will be at the end of the axes.
     * @param img
     * @return 
     */
    public Image increaseToPowerOf2(Image img) {
        
        int n0 = img.getWidth();
        int n1 = img.getHeight();

        int nn0 = 1 << (int)(Math.ceil(Math.log(n0)/Math.log(2)));
        int nn1 = 1 << (int)(Math.ceil(Math.log(n1)/Math.log(2)));

        if (nn0 == n0 && nn1 == n1) {
            return null;
        }
        
        Image img2 = img.createWithDimensions(nn0, nn1);
        for (int i = 0; i < n0; ++i) {
            for (int j = 0; j < n1; ++j) {
                img2.setRGB(i, j, img.getRGB(i, j));
            }
        }
        
        return img2;
    }
    
    public List<Set<PairInt>> findConnectedSameValueGroups(GreyscaleImage 
        img) {
        
        TIntObjectMap<Set<PairInt>> valuePointMap 
            = new TIntObjectHashMap<Set<PairInt>>();
        
        int w = img.getWidth();
        int h = img.getHeight();
        // O(N)
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int v = img.getValue(i, j);
                PairInt p = new PairInt(i, j);
                Set<PairInt> set = valuePointMap.get(v);
                if (set == null) {
                    set = new HashSet<PairInt>();
                    valuePointMap.put(v, set);
                }
                set.add(p);
            }
        }
        
        // TODO: consider replacing with a small tolerance around 
        //    index values, that is a clustering step by small amount in
        //    intensity along with search for contiguous
        
        List<Set<PairInt>> contigSets = new ArrayList<Set<PairInt>>();

        // DFS search is O(N*lg2(N)).  these are partitioned into 
        // indiv sets of DFS but unique sets so the overall runtime should
        // be similar, espec if improve the finder to create adjacency
        // maps instead of loop over all possible neighbor offsets.        
        TIntObjectIterator<Set<PairInt>> iter = valuePointMap.iterator();
        for (int i = 0; i < valuePointMap.size(); ++i) {
            
            iter.advance();
            
            // find the contiguous among set.
            Set<PairInt> set = iter.value();
            
            DFSConnectedGroupsFinder finder = new DFSConnectedGroupsFinder();
            finder.setToUse8Neighbors();
            finder.setMinimumNumberInCluster(1);
            finder.findConnectedPointGroups(set);
            int ns = finder.getNumberOfGroups();
            for (int j = 0; j < ns; ++j) {
                contigSets.add(finder.getXY(j));
            }
        }
        
        return contigSets;
    }
    
    /**
     Find peaks in an image as coordinate list 
     Peaks are the local maxima in a region of `2 * min_distance + 1`
     (i.e. peaks are separated by at least `min_distance`).
     If peaks are flat (i.e. multiple adjacent pixels have identical
     intensities), the coordinates of all such pixels are returned.
     If both `threshold_abs` and `threshold_rel` are provided, the maximum
     of the two is chosen as the minimum intensity threshold of peaks.

      adapted from
     https://github.com/scikit-image/scikit-image/blob/92a38515ac7222aab5e606f9de46caf5f503a7bd/skimage/feature/peak.py

     The implementation below is adapted from the scipy implementation which has
     * the following copyright:

     https://github.com/scikit-image/scikit-image/blob/master/LICENSE.txt

    -- begin scipy, skimage copyright ---
    Unless otherwise specified by LICENSE.txt files in individual
    directories, all code is

    Copyright (C) 2011, the scikit-image team
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in
        the documentation and/or other materials provided with the
        distribution.
     3. Neither the name of skimage nor the names of its contributors may be
        used to endorse or promote products derived from this software without
        specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
    IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
    INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
    IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
    -- end scipy, skimage copyright ---

    Note, in some places, scipy functions have been
    replaced with existing functions in this project in this implementation below.


     * @param img
     * @param minDistance
        Minimum number of pixels separating peaks in a region of `2 *
        min_distance + 1` (i.e. peaks are separated by at least
        `min_distance`).
        To find the maximum number of peaks, use `min_distance=1`.
      @param outputKeypoints0 the output row coordinates of keypoints
     * @param outputKeypoints1 the output col coordinates of keypoints
     */
    public void peakLocalMax(float[][] img, int minDistance,
        float thresholdRel,
        TIntList outputKeypoints0, TIntList outputKeypoints1) {

        int excludeBorder = minDistance;
        int numPeaks = Integer.MAX_VALUE;
        //int numPeaksPerLabel = Integer.MAX_VALUE;

        /*
        The peak local maximum function returns the coordinates of local peaks
        (maxima) in an image. A maximum filter is used for finding local maxima.
        This operation dilates the original image. After comparison of the dilated
        and original image, this function returns the coordinates or a mask of the
        peaks where the dilated image equals the original image.
        */

        int nRows = img.length;
        int nCols = img[0].length;

        //# Non maximum filter
        int size = 2 * minDistance + 1;
        float[][] imageMax = maximumFilter(img, size);
        assert(nRows == imageMax.length);
        assert(nCols == imageMax[0].length);
        //mask = image == image_max

        //debugPrint("before shift imageMax=", imageMax);

        // a fudge to match results of scipy which must store same windows at
        // locations shifted by minDistance or so in x and y from the
        // beginning of the sliding window
        applyShift(imageMax, minDistance, nRows, nCols);

        /*
        {//DEBUG
            float min = MiscMath.findMin(img);
            float max = MiscMath.findMax(img);
            System.out.println("min=" + min + " max=" + max);
            float[][] img2 = copy(img);
            MiscMath.applyRescale(img2, 0, 255);
            GreyscaleImage gsImg = new GreyscaleImage(nRows, nCols);
            for (int i = 0; i < nRows; ++i) {
                for (int j = 0; j < nCols; ++j) {
                    int v = Math.round(img2[i][j]);
                    if (v > 255) {
                        v = 255;
                    }
                    gsImg.setValue(i, j, v);
                }
            }
            MiscDebug.writeImage(gsImg, "_CURVATURE_");
        }*/

        //TODO: should be able to simplify the mask here

        // 1's where same, else 0's
        int[][] mask = new int[nRows][nCols];
        for (int i = 0; i < nRows; ++i) {
            mask[i] = new int[nCols];
            for (int j = 0; j < nCols; ++j) {
                if (img[i][j] == imageMax[i][j]) {
                    mask[i][j] = 1;
                }
            }
        }

        //debugPrint("0 mask=", mask);


        // exclude border
        for (int i = 0; i < nRows; ++i) {
            if ((i < excludeBorder) || (i > (nRows - 1 - excludeBorder))){
                Arrays.fill(mask[i], 0);
            } else {
                Arrays.fill(mask[i], 0, excludeBorder, 0);
                Arrays.fill(mask[i], nCols - excludeBorder, nCols, 0);
            }
        }


        // find top peak candidates above a threshold.
        // TODO: should this use mask so excluding borders?
        float thresholdAbs = MiscMath.findMin(img);
        float thresholdMax = thresholdRel * MiscMath.findMax(img);
        thresholdAbs = Math.max(thresholdAbs, thresholdMax);

        // mask &= image > 0.1
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                if (imageMax[i][j] > thresholdAbs) {
                    mask[i][j] &= 1;
                } else {
                    mask[i][j] = 0;
                }
            }
        }

        /*
        {//DEBUG
            try{
            int min = MiscMath.findMin(mask);
            int max = MiscMath.findMax(mask);
            System.out.println("min=" + min + " max=" + max);
            float factor = 255.f;
            GreyscaleImage gsImg = new GreyscaleImage(nRows, nCols);
            for (int i = 0; i < nRows; ++i) {
                for (int j = 0; j < nCols; ++j) {
                    int v = Math.round(factor * mask[i][j]);
                    if (v > 255) {
                        v = 255;
                    }
                    gsImg.setValue(i, j, v);
                }
            }
            ImageDisplayer.displayImage("mask", gsImg);
            int z = 1;
            } catch(Exception e) {}
        }
        */

        //debugPrint("mask &= image > " + thresholdAbs, mask);

        // Select highest intensities (num_peaks)
        // expected output is [row index, col index, ...]

        //TODO: should num_peaks be this.nKeypoints?  re-read paper...
        if (numPeaks == Integer.MAX_VALUE) {
            // find non-zero pixels in mask
            float[] values = new float[nRows * nCols];
            int[] pixIdxs = new int[values.length];
            int count = 0;
            for (int i = 0; i < mask.length; ++i) {
                for (int j = 0; j < mask[i].length; ++j) {
                    if (mask[i][j] > 0.f) {
                        values[count] = img[i][j];
                        //(row * width) + col
                        pixIdxs[count] = (j * nRows) + i;
                        count++;
                    }
                }
            }
            values = Arrays.copyOf(values, count);
            pixIdxs = Arrays.copyOf(pixIdxs, count);
            MultiArrayMergeSort.sortByDecr(values, pixIdxs);
            for (int i = 0; i < values.length; ++i) {
                int pixIdx = pixIdxs[i];
                int jj = pixIdx/nRows;
                int ii = pixIdx - (jj * nRows);
                outputKeypoints0.add(ii);
                outputKeypoints1.add(jj);
            }
        } else {
            //need to sort to keep top numPeaks
            FixedSizeSortedVector<ORB.Pix> vec = new
                FixedSizeSortedVector<ORB.Pix>(numPeaks, ORB.Pix.class);
            for (int i = 0; i < mask.length; ++i) {
                for (int j = 0; j < mask[i].length; ++j) {
                    if (mask[i][j] > 0.f) {
                        ORB.Pix pix = new ORB.Pix(i, j, Float.valueOf(img[i][j]));
                        vec.add(pix);
                    }
                }
            }
            for (int i = 0; i < vec.getNumberOfItems(); ++i) {
                ORB.Pix pix = vec.getArray()[i];
                outputKeypoints0.add(pix.i);
                outputKeypoints1.add(pix.j);
            }
        }
    }

    private void applyShift(float[][] imageMax, int minDistance, int nRows,
        int nCols) {

        for (int i = 0; i < nRows; ++i) {
            System.arraycopy(imageMax[i], 0, imageMax[i], minDistance,
                nCols - minDistance);
            for (int j = 0; j < minDistance; ++j) {
                imageMax[i][j] = 0;
            }
            for (int j = (nCols - minDistance); j < nCols; ++j) {
                imageMax[i][j] = 0;
            }
        }
        //debugPrint("shift 0 imageMax=", imageMax);
        for (int j = 0; j < nCols; ++j) {
            for (int i = (nRows - minDistance); i >= minDistance; --i) {
                imageMax[i][j] = imageMax[i - minDistance][j];
            }
            for (int i = 0; i < minDistance; ++i) {
                imageMax[i][j] = 0;
            }
            for (int i = (nRows - minDistance); i < nRows; ++i) {
                imageMax[i][j] = 0;
            }
        }
        //debugPrint("shifted imageMax=", imageMax);
    }
    
    /**
     * @author nichole
     * @param img
     * @param size
     * @return
     */
    public float[][] maximumFilter(float[][] img, int size) {

        int nRows = img.length;
        int nCols = img[0].length;

        // return_value = out
        float[][] out = new float[nRows][nCols];
        for (int i = 0; i < nRows; ++i) {
            out[i] = new float[nCols];
        }

        // have adapted median window algorithm for this:
        StatsInSlidingWindow maxWindow = new StatsInSlidingWindow();
        maxWindow.calculateMaximum(img, out, size, size);

        return out;
    }
    
    /**
     * 
     * @param img
     * @param g the kernel to convolve img with at point (x,y).
     * Note that it's assumed the kernel is already normalized to sum 0.
     * @return 
     */
    public Complex[][] convolveWithKernel(final Complex[][] img, 
        Complex[][] g) {
        
        Complex[][] c = new Complex[img.length][];
        
        for (int x = 0; x < img.length; ++x) {
            c[x] = new Complex[img[0].length];
            for (int y = 0; y < img[0].length; ++y) {
                ComplexModifiable sum = 
                    convolvePointWithKernel(img, x, y, g);
                c[x][y] = new Complex(sum.re(), sum.im());
            }
        }
       
        return c;
    }
    
     /**
     * 
     * @param img
     * @param x
     * @param y
     * @param g the kernel to convolve img with at point (x,y).
     * Note that it's assumed the kernel is already normalized to sum 0.
     * @return sum of convolution at point x,y
     */
    public ComplexModifiable convolvePointWithKernel(final Complex[][] img, int x, 
        int y, Complex[][] g) {
        
        int n0 = g.length;
        int n1 = g[0].length;
        
        int h0 = (n0 >> 1);
        int h1 = (n1 >> 1);
        
        int w = img.length;
        int h = img[0].length;
        
        ComplexModifiable sum = new ComplexModifiable(0, 0);

        for (int i = -h0; i < (n0 - h0); i++) {
            int x2 = x + i;
            if (x2 < 0) {
               // replicate
               x2 = -1*x2 - 1;
               if (x2 > (w - 1)) {
                   x2 = w - 1;
               }
            } else if (x2 > (w - 1)) {
                int diff = x2 - w;
                x2 = w - diff - 1;
                if (x2 < 0) {
                    x2 = 0;
                }
            }
            for (int j = -h1; j < (n1 - h1); j++) {
                int y2 = y + j;
                if (y2 < 0 || (y2 > (h - 1))) {
                    continue;
                }
                if (y2 < 0) {
                   // replicate
                   y2 = -1*y2 - 1;
                   if (y2 > (h - 1)) {
                       y2 = h - 1;
                   }
                } else if (y2 > (h - 1)) {
                    int diff = y2 - h;
                    y2 = h - diff - 1;
                    if (y2 < 0) {
                        y2 = 0;
                    }
                }
            
                Complex gg = g[i + h0][j + h1];

                Complex point = img[x2][y2];
                
                Complex m = point.times(gg);

                sum.plus(m);                
            }
        }
        
        return sum;
    } 
    
    public Set<PairInt> createWindowOfPoints(int x, int y, 
        int windowHalfWidth, int imageWidth, int imageHeight) {
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        for (int i = (x - windowHalfWidth); i <= (x + windowHalfWidth); ++i) {
            if (i < 0 || i > (imageWidth - 1)) {
                continue;
            }
            for (int j = (y - windowHalfWidth); j <= (y + windowHalfWidth); ++j) {
                if (j < 0 || j > (imageHeight - 1)) {
                    continue;
                }
                points.add(new PairInt(i, j));
            }
        }
        
        return points;
    }
    
    /**
     * given an input image, creates a decimation pyramid with 
     * median smoothing followed by either integer or bilinear
     * interpolation down-sampling.
     * This method returns images down-sampled by scale sizes that
     * are a factor of 2 from each until image size
     * is smaller than decimationLimit.  Then the
     * method creates discrete scales in between the factors of 2 
     * pyramid, such as 1.25, 1.5, 1.75, then bisecting scales
     * with 3, 5, 12, etc.
     * @param input
     * @param decimationLimit limit to smallest image dimension size returned
     * @return 
     */
    public List<GreyscaleImage> buildPyramid2(GreyscaleImage input,
        int decimationLimit) {
    
        int nBetween = 3;
        
        return buildPyramid2(input, decimationLimit, nBetween);
    }
    
    public List<List<GreyscaleImage>> buildColorPyramid(Image img,
        boolean buildLarger) {
        
        List<List<GreyscaleImage>> pyr = new ArrayList<List<GreyscaleImage>>();
        
        GreyscaleImage r = img.copyRedToGreyscale();
        GreyscaleImage g = img.copyGreenToGreyscale();
        GreyscaleImage b = img.copyBlueToGreyscale();
        List<GreyscaleImage> gsR;
        List<GreyscaleImage> gsG;
        List<GreyscaleImage> gsB;
        
        if (buildLarger) {
            
            ImageProcessor imageProcessor = new ImageProcessor();
            
            gsR = imageProcessor.buildPyramid2(r, 32);
            gsG = imageProcessor.buildPyramid2(g, 32);
            gsB = imageProcessor.buildPyramid2(b, 32);
            
        } else {
            
            MedianTransform mt = new MedianTransform();
            
            gsR = new ArrayList<GreyscaleImage>();
            gsG = new ArrayList<GreyscaleImage>();
            gsB = new ArrayList<GreyscaleImage>();
            
            mt.multiscalePyramidalMedianTransform2(r, gsR, 32);
            mt.multiscalePyramidalMedianTransform2(g, gsG, 32);
            mt.multiscalePyramidalMedianTransform2(b, gsB, 32);
        }
        
        assert(gsR.size() == gsG.size());
        assert(gsR.size() == gsB.size());

        for (int i = 0; i < gsR.size(); ++i) {
            GreyscaleImage r2 = gsR.get(i);
            GreyscaleImage g2 = gsG.get(i);
            GreyscaleImage b2 = gsB.get(i);

            List<GreyscaleImage> rgb = new ArrayList<GreyscaleImage>();
            rgb.add(r2);
            rgb.add(g2);
            rgb.add(b2);
            
            pyr.add(rgb);
        }
        
        return pyr;
    }
    
    public List<GreyscaleImage> buildPyramid(GreyscaleImage img,
        boolean buildLarger) {
                
        GreyscaleImage cp = img.copyImage();
        
        List<GreyscaleImage> out;
        
        if (buildLarger) {
            
            ImageProcessor imageProcessor = new ImageProcessor();
            
            out = imageProcessor.buildPyramid2(cp, 32);
            
        } else {
            
            MedianTransform mt = new MedianTransform();
            
            out = new ArrayList<GreyscaleImage>();
            
            mt.multiscalePyramidalMedianTransform2(cp, out, 32);
        }
        
        return out;
    }
    
    /**
     * given an input image, creates a decimation pyramid with 
     * median smoothing followed by either integer or bilinear
     * interpolation down-sampling.
     * This method returns images down-sampled by scale sizes that
     * are a factor of 2 from each until image size
     * is smaller than decimationLimit.  Then the
     * method creates discrete scales in between the factors of 2 
     * pyramid, such as 1.25, 1.5, 1.75, then bisecting scales
     * with 3, 5, 12, etc.
     * @param input
     * @param decimationLimit limit to smallest image dimension size returned
     * @param nBetween the number of images to create in between the powers
     * of 2.
     * @return 
     */
    public List<GreyscaleImage> buildPyramid2(GreyscaleImage input,
        int decimationLimit, int nBetween) {
        
        List<GreyscaleImage> output = new ArrayList<GreyscaleImage>();

        MedianTransform mt = new MedianTransform();
        mt.multiscalePyramidalMedianTransform2(input, output, decimationLimit);

        if (output.size() == 1) {
            return output;
        }

        List<GreyscaleImage> output2 = new ArrayList<GreyscaleImage>();
        
        // add an image in between each after output[2]
        for (int i = 0; i < output.size() - 1; ++i) {
            output2.add(mt.decimateImage(output.get(i), 1.5f, 0, 255));
        }
       
        output.addAll(output2);
       
        Collections.sort(output, new DecreasingSizeComparator());

        return output;
        
        /*
        a gaussian pyramid based upon a kernel of sigma 0.707
        can be built recursively, but the number of iterations
        to reach scale factors larger than 3 is increasingly
        very large number of recursions.
        
        FWHM = 2.35 * sigma
        for recursive gaussian and kernel0 w/ sigma=sqrt(2)/2 = 0.707
                                                     sigma   FWHM
        s1^2 = (sqrt(2)/2)^2 = 0.5                    0.7     1.67
        s2^2 = 0.5 + (sqrt(2)/2)^2 = 1                1       2.35                 
        s3^2 = 1 + (sqrt(2)/2)^2                      1.22    2.9
        s4^2 = 1.5 + (sqrt(2)/2)^2                    1.4     3.3
        s5^2 = 2 + (sqrt(2)/2)^2                      1.58    3.72                      
        s6^2 = 2.5 + (sqrt(2)/2)^2                    1.7     4.07
        s7^2 = 3.0 + (sqrt(2)/2)^2                    1.87    4.4  
        ...
        s17^2 = 8.0 + (sqrt(2)/2)^2                   2.83    6.65  <-- factor of 4 from s1
        
        Alternatively, making a kernel of size sigma>2 needs larger
        kernels too so the O(N) has a large constant factor in front of it.
        
        A hybrid solution would be to use the current
        pyramidal median transform which returns images 
        blurred and downsampled by factors of 2.
        Then add to that, discrete samplings of scale sizes
        that improve the pyramid.
        1   2   4   8   16
        for example, calculate for scale factors 1.25, 1.5, 1.75
        then bisecting the existing pyramidal scales to make
        3, 6, and 12, etc.
        Might consider a version of this method which provides more
        scales at higher factors...
        */
    }

    /**
     * use bilinear interpolation to downsample a two dimensional array.
     * 
     * adapted from pseudocode at
     * http://tech-algorithm.com/articles/bilinear-image-scaling/
    
     * @param pixels
     * @param w2 output length of first dimension of array
     * @param h2 output length of second dimension of array
     * @return 
     */
    public float[][] bilinearDownSampling(float[][] pixels,  
        int w2, int h2) {
        
        float[][] out = new float[w2][h2];
        for (int i = 0; i < w2; ++i) {
            out[i] = new float[h2];  
        }
        
        int w = pixels.length; 
        int h = pixels[0].length;
        
        int x, y;
        float A, B, C, D, gray;
        
        float xRatio = (float)w/(float)w2;
        float yRatio = (float)h/(float)h2;
        float xDiff, yDiff;
        int offset = 0 ;
        for (int i = 0; i < w2; i++) {
            for (int j = 0; j < h2; j++) {
                x = (int)(xRatio * i);
                y = (int)(yRatio * j);
                xDiff = (xRatio * i) - x;
                yDiff = (yRatio * j) - y;

                A = pixels[x][y];
                B = pixels[x][y+1];
                C = pixels[x+1][y];
                D = pixels[x+1][y+1];

                // Y = A(1-w)(1-h) + B(w)(1-h) + C(h)(1-w) + Dwh
                gray = 
                    A * (1 - xDiff) * (1 - yDiff) + 
                    C * (xDiff) * (1 - yDiff) +
                    B * (yDiff) * (1 - xDiff) +  
                    D * (xDiff * yDiff);

                out[i][j] = gray;                                   
            }
        }
        
        return out;
    }
    
    /**
     * use bilinear interpolation to downsample a two dimensional array.
     * Assumes that valid pixel values are 0 to 255 and clamps if outside range.
     * 
     * adapted from pseudocode at
     * http://tech-algorithm.com/articles/bilinear-image-scaling/
    
     * @param input
     * @param w2 output length of first dimension of array
     * @param h2 output length of second dimension of array
     * @return 
     */
    public GreyscaleImage bilinearDownSampling(GreyscaleImage input,
        int w2, int h2) {
   
        return bilinearDownSampling(input, w2, h2, 0, 255);
    }
   
    /**
     * use bilinear interpolation to downsample a two dimensional array.
     * 
     * adapted from pseudocode at
     * http://tech-algorithm.com/articles/bilinear-image-scaling/
    
     * @param input
     * @param w2 output length of first dimension of array
     * @param h2 output length of second dimension of array
     * @param minValue value to clamp results to if less than this
     * @param maxValue value to clamp results to if larger than this
     * @return 
     */
    public GreyscaleImage bilinearDownSampling(GreyscaleImage input,
        int w2, int h2, int minValue, int maxValue) {
        
        GreyscaleImage out = input.createWithDimensions(w2, h2);
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        int x, y;
        float A, B, C, D, gray;
        
        float xRatio = (float)w/(float)w2;
        float yRatio = (float)h/(float)h2;
        float xDiff, yDiff;
        int offset = 0 ;
        for (int i = 0; i < w2; i++) {
            for (int j = 0; j < h2; j++) {
                x = (int)(xRatio * i);
                y = (int)(yRatio * j);
                xDiff = (xRatio * i) - x;
                yDiff = (yRatio * j) - y;

                A = input.getValue(x, y);
                B = input.getValue(x, y+1);
                C = input.getValue(x+1, y);
                D = input.getValue(x+1, y+1);

                // Y = A(1-w)(1-h) + B(w)(1-h) + C(h)(1-w) + Dwh
                gray = 
                    A * (1 - xDiff) * (1 - yDiff) + 
                    C * (xDiff) * (1 - yDiff) +
                    B * (yDiff) * (1 - xDiff) +  
                    D * (xDiff * yDiff);

                int v = Math.round(gray);
                if (v < minValue) {
                    v = minValue;
                } else if (v > maxValue) {
                    v = maxValue;
                }
                
                out.setValue(i, j, v);                              
            }
        }
        
        return out;
    }
    
    /**
     * use bilinear interpolation to downsample a two dimensional array.
     * 
     * adapted from pseudocode at
     * http://tech-algorithm.com/articles/bilinear-image-scaling/
    
     * @param input
     * @param w2 output length of first dimension of array
     * @param h2 output length of second dimension of array
     * @param minValue value to clamp results to if less than this
     * @param maxValue value to clamp results to if larger than this
     * @return 
     */
    public Image bilinearDownSampling(Image input,
        int w2, int h2, int minValue, int maxValue) {
        
        Image out = input.createWithDimensions(w2, h2);
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        int x, y, rI, gI, bI;
        float A, B, C, D, red, green, blue;
        
        float xRatio = (float)w/(float)w2;
        float yRatio = (float)h/(float)h2;
        float xDiff, yDiff;
        int offset = 0 ;
        for (int i = 0; i < w2; i++) {
            for (int j = 0; j < h2; j++) {
                x = (int)(xRatio * i);
                y = (int)(yRatio * j);
                xDiff = (xRatio * i) - x;
                yDiff = (yRatio * j) - y;

                A = input.getR(x, y);
                B = input.getR(x, y+1);
                C = input.getR(x+1, y);
                D = input.getR(x+1, y+1);

                // Y = A(1-w)(1-h) + B(w)(1-h) + C(h)(1-w) + Dwh
                red = A * (1 - xDiff) * (1 - yDiff) + 
                    C * (xDiff) * (1 - yDiff) +
                    B * (yDiff) * (1 - xDiff) +  
                    D * (xDiff * yDiff);

                rI = Math.round(red);
                if (rI < minValue) {
                    rI = minValue;
                } else if (rI > maxValue) {
                    rI = maxValue;
                }
                
                A = input.getG(x, y);
                B = input.getG(x, y+1);
                C = input.getG(x+1, y);
                D = input.getG(x+1, y+1);
                green = A * (1 - xDiff) * (1 - yDiff) + 
                    C * (xDiff) * (1 - yDiff) +
                    B * (yDiff) * (1 - xDiff) +  
                    D * (xDiff * yDiff);
                gI = Math.round(green);
                if (gI < minValue) {
                    gI = minValue;
                } else if (gI > maxValue) {
                    gI = maxValue;
                }
                
                A = input.getB(x, y);
                B = input.getB(x, y+1);
                C = input.getB(x+1, y);
                D = input.getB(x+1, y+1);
                blue = A * (1 - xDiff) * (1 - yDiff) + 
                    C * (xDiff) * (1 - yDiff) +
                    B * (yDiff) * (1 - xDiff) +  
                    D * (xDiff * yDiff);
                bI = Math.round(blue);
                if (bI < minValue) {
                    bI = minValue;
                } else if (bI > maxValue) {
                    bI = maxValue;
                }
                
                out.setRGB(i, j, rI, gI, bI);                              
            }
        }
        
        return out;
    }
    
    /**
     * comparator that assumes can compare by widths along,
     * unless there is a tie, then it uses height.
     */
    private class DecreasingSizeComparator implements 
        Comparator<GreyscaleImage> {

        @Override
        public int compare(GreyscaleImage o1, 
            GreyscaleImage o2) {
        
            int w1 = o1.getWidth();
            int w2 = o2.getWidth();
            
            if (w1 > w2) {
                return -1;
            } else if (w1 < w2) {
                return 1;
            }
        
            int h1 = o1.getHeight();
            int h2 = o2.getHeight();
            
            if (h1 > h2) {
                return -1;
            } else if (h1 < h2) {
                return 1;
            }
            
            return 0;
        }        
    }
    
    public GreyscaleImage[] createCIELABImages(ImageExt img) {
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        GreyscaleImage ells = new GreyscaleImage(w, h);
        GreyscaleImage as = new GreyscaleImage(w, h);
        GreyscaleImage bs = new GreyscaleImage(w, h);
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        float[] mins = new float[]{0,  -190, -113};
        float[] maxs = new float[]{105, 103, 99};
        float[] scales = new float[3];
        
        for (int i = 0; i < scales.length; ++i) {
            float r = maxs[i] - mins[i];
            scales[i] = 255.f/r;
        }
        
        int n = img.getNPixels();
        int v;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                
                int r = img.getR(i, j);
                int g = img.getG(i, j);
                int b = img.getB(i, j);
                
                float[] lab = cieC.rgbToCIELAB1931(r, g, b);
             
                v = (int)((lab[0] - mins[0])*scales[0]);
                ells.setValue(i, j, v);
                
                v = (int)((lab[1] - mins[1])*scales[1]);
                as.setValue(i, j, v);
                
                v = (int)((lab[2] - mins[2])*scales[2]);
                bs.setValue(i, j, v);
            }
        }
        
        return new GreyscaleImage[]{ells, as, bs};
    }
    
    /**
     * convert the image to cie l*a*b* and then use a and b
     * to calculate polar angle around 0 in degrees.
     * If maxV of 360, returns full value image, 
     * else if is 255, scales the values to max value of 255, etc.
     * @param img
     * @param maxV
     * @return 
     */
    public GreyscaleImage createCIELAB1931Theta(Image img, int maxV) {
        
        return createCIELAB1931Theta(img, maxV, 0);
    }
    
    /**
     * convert the image to cie l*a*b* and then use a and b
     * to calculate polar angle around 0 in degrees.
     * If maxV of 360, returns full value image, 
     * else if is 255, scales the values to max value of 255, etc.
     * @param img
     * @param maxV
     * @param offsetTheta an offset in degrees from 0 to shift
     *    values by.  This is useful for a quick look between 
     * results of offset=0 and offset=10, for example, to look
     * at the wrap around point for 360 to 0
     * @return 
     */
    public GreyscaleImage createCIELAB1931Theta(Image img, int maxV, 
        int offset) {
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        GreyscaleImage theta = null;
        if (maxV < 256) {
            theta = new GreyscaleImage(w, h);
        } else {
            theta = new GreyscaleImage(w, h, 
                GreyscaleImage.Type.Bits32FullRangeInt);
        }
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        float[] mins = new float[]{0,  -190, -113};
        float[] maxs = new float[]{105, 103, 99};
        float[] scales = new float[3];
        
        for (int i = 0; i < scales.length; ++i) {
            float r = maxs[i] - mins[i];
            scales[i] = 255.f/r;
        }
        
        double ts = (double)maxV/(double)359;
        
        int n = img.getNPixels();
        int v;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                
                int r = img.getR(i, j);
                int g = img.getG(i, j);
                int b = img.getB(i, j);
                
                float[] lab = cieC.rgbToCIELAB1931(r, g, b);
             
                //float v1 = (lab[1] - mins[1])*scales[1];
                //float v2 = (lab[2] - mins[2])*scales[2];
                float v1 = lab[1];
                float v2 = lab[2];
               
                double t = Math.atan2(v2, v1);
                t *= (180./Math.PI);
                t += offset;
                if (t < 0) {
                    t += 360;
                } else if (t > 359) {
                    t -= 360;
                }
                t *= ts;
                v = (int)t;
                
                theta.setValue(i, j, v);
            }
        }
        
        return theta;
    }
    
    /**
     * convert the image to cie l*a*b* and then use a and b
     * to calculate polar angle around 0 in degrees.
     * If maxV of 360, returns full value image, 
     * else if is 255, scales the values to max value of 255, etc.
     * @param img
     * @param maxV
     * @return 
     */
    public int calculateCIELABTheta(int red, int green, int blue, int maxV) {
        
        CIEChromaticity cieC = new CIEChromaticity();
       
        double ts = (double)maxV/(double)359;
        
        float[] lab = cieC.rgbToCIELAB1931(red, green, blue);
        
        float v1 = lab[1];
        float v2 = lab[2];

        double t = Math.atan2(v2, v1);
        t *= (180./Math.PI);
        if (t < 0) {
            t += 360;
        } else if (t > 359) {
            t -= 360;
        }
        t *= ts;
        
        return (int)t;
    }
    
    /**
     * create 3 images of the LCG color space where LCH is the 
     * luminosity, magnitude, and polar angle of CIE LUV color space.
     * @param img
     * @return 
     */
    public GreyscaleImage[] createLCHForLUV(Image img) {
        
        int w = img.getWidth();
        int h = img.getHeight();
        int n = img.getNPixels();
        
        /*
        range of CIE LUV using default standard illumination of
        D65 daylight is:
        L       0 to 104.5
        u   -86.9 to 183.8
        v  -141.4 to 112.3
        luminosity L*  0 to 104.5
        magnitude, m:  sqrt(2) * 183.8 = 260
        angle,     a:  0 to 359
        */
        
        float[] factors = new float[]{255.f/104.5f, 255.f/260.f, 255.f/359.f};
        
        GreyscaleImage[] output = new GreyscaleImage[3];
        for (int i = 0; i < output.length; ++i) {
            output[i] = new GreyscaleImage(w, h);
        }
        
        CIEChromaticity cieC = new CIEChromaticity();
                
        int v;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                
                int r = img.getR(i, j);
                int g = img.getG(i, j);
                int b = img.getB(i, j);
                
                float[] lma = cieC.rgbToPolarCIELUV(r, g, b);
             
                for (int k = 0; k < output.length; ++k) {
                    if (k == 2 && lma[k] > 359.f) {
                        lma[k] = 360 - lma[k];
                    }
                    v = Math.round(factors[k] * lma[k]);
                    if (v > 255) {
                        v = 255;
                    } else if (v < 0) {
                        v = 0;
                    }
                    output[k].setValue(i, j, v);
                }
            }
        }
        
        return output;
    }
   
    /**
     * convert the image to cie luv and then calculate polar angle of u and v
     * around 0 in degrees.
     * If maxV of 360, returns full value image, 
     * else if is 255, scales the values to max value of 255, etc.
     * @param img
     * @param maxV
     * @return 
     */
    public GreyscaleImage createCIELUVTheta(Image img, int maxV) {
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        GreyscaleImage theta = null;
        if (maxV < 256) {
            theta = new GreyscaleImage(w, h);
        } else {
            theta = new GreyscaleImage(w, h, 
                GreyscaleImage.Type.Bits32FullRangeInt);
        }
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        double ts = (double)maxV/(double)359;
        
        int n = img.getNPixels();
        int v;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                
                int r = img.getR(i, j);
                int g = img.getG(i, j);
                int b = img.getB(i, j);
                
                float[] lma = cieC.rgbToPolarCIELUV(r, g, b);
             
                double t = lma[2];
                t *= ts;
                v = (int)t;
                
                theta.setValue(i, j, v);
            }
        }
        
        return theta;
    }
    
    /**
     * convert the image to cie lab and then calculate polar angle of u and v
     * around 0 in degrees.
     * If maxV of 360, returns full value image, 
     * else if is 255, scales the values to max value of 255, etc.
     * @param img
     * @param maxV
     * @return 
     */
    public GreyscaleImage createCIELABTheta(Image img, int maxV) {
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        GreyscaleImage theta = null;
        if (maxV < 256) {
            theta = new GreyscaleImage(w, h);
        } else {
            theta = new GreyscaleImage(w, h, 
                GreyscaleImage.Type.Bits32FullRangeInt);
        }
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        double ts = (double)maxV/(double)359;
        
        int n = img.getNPixels();
        int v;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                
                int r = img.getR(i, j);
                int g = img.getG(i, j);
                int b = img.getB(i, j);
                
                float[] lab = cieC.rgbToCIELAB(r, g, b);
             
                double t = Math.atan2(lab[2], lab[1]);
                t *= (180. / Math.PI);
                if (t < 0) {
                    t += 360;
                } else if (t > 359) {
                    t -= 360;
                }

                t *= ts;
                v = (int)t;
                
                theta.setValue(i, j, v);
            }
        }
        
        return theta;
    }
    
    public TIntObjectMap<VeryLongBitString> createAdjacencyMap(
        TObjectIntMap<PairInt> pointIndexMap, List<Set<PairInt>> labeledPoints) {

        int n = labeledPoints.size();

        TIntObjectMap<VeryLongBitString> output
            = new TIntObjectHashMap<VeryLongBitString>();

        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;

        for (int label = 0; label < n; ++label) {
            Set<PairInt> set = labeledPoints.get(label);
            VeryLongBitString nbrs = new VeryLongBitString(n);
            
            boolean aloSet = false;
            
            for (PairInt p : set) {
                int x = p.getX();
                int y = p.getY();
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    PairInt p2 = new PairInt(x2, y2);
                    if (pointIndexMap.containsKey(p2)) {
                        int label2 = pointIndexMap.get(p2);
                        if (label2 != label) {
                            nbrs.setBit(label2);
                            aloSet = true;
                        }
                    }
                }
            }
            if (aloSet) {
                output.put(label, nbrs);
            }
        }

        return output;
    }
    
    public void applyUnsharpMask(GreyscaleImage img) {
        
        //NOTE: if make a color version of this, should use a color
        // space such as CIE LAB or LUV and operate on the L only to
        // preserve color (or B of HSB).
        
        // NOTE: useful in looking at the general concept of
        //  original, blurred and threshold comparison was code copyrighted by
        //  Romain Guy available here:
        //  http://www.java2s.com/Code/Java/Advanced-Graphics/UnsharpMaskDemo.htm
        //  might need to place the copytight here for derivative of.
        
        float amount = 0.2f;
        int radius = 50;
        float threshold = 0;
        
        applyUnsharpMask(img, amount, radius, threshold);
        
    }
    
    /**
     * if a single pixel differs from its neighbors by more than 
     * 2o or 2 sigma, the value gets set to the neighbors
     * @param img 
     */
    public void singlePixelFilter(GreyscaleImage img) {
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        int n = img.getNPixels();
        int w = img.getWidth();
        int h = img.getHeight();
        
        for (int i = 0; i < n; ++i) {
            
            int x = img.getCol(i);
            int y = img.getRow(i);
            
            float diff;
            float avg = 0;
            int ns = 0;
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < 0 || y2 < 0 || x2 >= w || y2 >= h) {
                    continue;
                }
                avg += img.getValue(x2, y2);
                ns++;
            }
            avg /= (float)ns;
            if (ns == 0) {
                continue;
            }
            
            float stdv = 0;
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < 0 || y2 < 0 || x2 >= w || y2 >= h) {
                    continue;
                }
                diff = img.getValue(x2, y2) - avg;
                stdv += (diff * diff);
            }
            stdv = (float)Math.sqrt(stdv/((float)ns - 1.0f));
        
            diff = Math.abs(img.getValue(x, y) - avg);
            if (diff > (2.*stdv)) {
                img.setValue(x, y, (int)avg);
            }
        }
    }
    
    public void applyUnsharpMask(GreyscaleImage img, float percentage,
        int radius, float threshold) {
   
        MedianSmooth ms = new MedianSmooth();
        GreyscaleImage blurred = ms.calculate(img, 3, 3);
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            int v = img.getValue(i);
            int blur = blurred.getValue(i);
            
            int diff = v - blur;
            if (Math.abs(diff) >= threshold) {
                v = (int) (percentage * diff + v);
                if (v < 0) {
                    v = 0;
                } else if (v > 255) {
                    v = 255;
                }
                img.setValue(i, v);
            }
        }
    }
    
    /**
     * use an adaptive threshold with window size to add the threshold
     * to the image where binarization would have created a "1".
     * This enhances shadows and dark edges to an extreme that might not
     * be desirable for all uses.
     * 
     * @param img 
     */
    public void enhanceContrast(GreyscaleImage img, int window) {
       
        GreyscaleImage inImg = img;
        double[][] outImg;
        
        int w = inImg.getWidth();
        int h = inImg.getHeight();
        
        double[][] sTable = new double[w][];
        for (int i = 0; i < w; ++i) {
            sTable[i] = new double[h];
            for (int j = 0; j < h; ++j) {
                sTable[i][j] = inImg.getValue(i, j);
            }
        }
        
        AdaptiveThresholding th =
            new AdaptiveThresholding();
        outImg = th.createAdaptiveThresholdImage(sTable, window, 0.2);
        
        double v;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                double t = outImg[i][j];
                v = inImg.getValue(i, j);
                if (v > t) {
                    // adding threshold back to emphasize contrast
                    // instead of binarization:
                    v += t;                    
                }
                outImg[i][j] = v;
            }
        }
        
        MiscMath.applyRescale(outImg, 0, 255);
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                inImg.setValue(i, j, (int)outImg[i][j]);
            }
        }
    }

    // TODO: implement the methods in 
    // http://www.merl.com/publications/docs/TR2008-030.pdf
    // for an O(n) filter.
    // "Constant Time O(1) Bilateral Filtering" by Porikli
    //public void applyBiLateralFilter(Image img) {
    //}
    
    // and trilateral filter by Tumblin et al. 2003
    
}
