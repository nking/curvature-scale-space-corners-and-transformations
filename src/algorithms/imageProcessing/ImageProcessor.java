package algorithms.imageProcessing;

import algorithms.FixedSizeSortedVector;
import algorithms.bipartite.MinHeapForRT2012;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.util.AngleUtil;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.*;
import algorithms.util.*;
import algorithms.VeryLongBitString;
import algorithms.heapsAndPQs.HeapNode;
import algorithms.sort.MultiArrayMergeSort;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import no.uib.cipr.matrix.Matrix;
import thirdparty.ca.uol.aig.fftpack.Complex1D;

/**
 *
 * @author nichole
 */
public class ImageProcessor {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    /**
     * use a sobel 
     * (1D Gaussian w/ sigma=0.71 convolution then
     * 1D first deriv gaussian sigma=0.5 convolution)
     * and return gradients in X and y. note the image may contain
     * negative values.
     * @param input
     * @return gX, gY
     */
    public GreyscaleImage[] createSobelGradients(GreyscaleImage input) {

        float[] kernel = Gaussian1DFirstDeriv.getBinomialKernelSigmaZeroPointFive();

        GreyscaleImage gX = input.copyToFullRangeIntImage();
        /*
         1
         2  * [1 0 -1]
         1
        */
        applyKernel1D(gX, new float[]{0.5f, 1.f, 0.5f}, false);
        applyKernel1D(gX, kernel, true);
        
        GreyscaleImage gY = input.copyToFullRangeIntImage();
        applyKernel1D(gY, kernel, false);
        applyKernel1D(gY, new float[]{0.5f, 1.f, 0.5f}, true);

        return new GreyscaleImage[]{gX, gY};
    }
    
    /**
     * apply a 1D first deriv gaussian sigma=0.5 convolution)
     * and return gradients in X and y. note the image may contain
     * negative values.
     * @param input
     * @return
     */
    public GreyscaleImage[] createGradients(GreyscaleImage input) {

        //[1 0 -1]
        float[] kernel = Gaussian1DFirstDeriv.getBinomialKernelSigmaZeroPointFive();

        GreyscaleImage gX = input.copyToFullRangeIntImage();
        
        applyKernel1D(gX, kernel, true);
        
        GreyscaleImage gY = input.copyToFullRangeIntImage();
        applyKernel1D(gY, kernel, false);

        return new GreyscaleImage[]{gX, gY};
    }
    
    /**
     * apply a 1D first deriv gaussian sigma=0.5 convolution)
     * and return gradients in X and y. note the image may contain
     * negative values.
     * @param input
     * @return
     */
    public GreyscaleImage[] createCentralDifferenceGradients(GreyscaleImage input) {

        if (input.getWidth() < 3 || input.getHeight() < 3) {
            throw new IllegalStateException("image dimensions must be at least 3 pixels each");
        }
        //[1 0 -1]
        GreyscaleImage gX = input.copyToFullRangeIntImage();
        
        GreyscaleImage gY = input.copyToFullRangeIntImage();
        
        int w = input.getWidth();
        int h = input.getHeight();
     
        int v0, v1;
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {                
                if (i > 0 && (i + 1) < w) {
                    //gx[i][j] = img[i-1][j] - img[i+1][j]
                    v0 = input.getValue(i - 1, j);
                    v1 = input.getValue(i + 1, j);
                    gX.setValue(i, j, (v1 - v0));
                } else if (i == 0) {
                    // handle boundaries
                    v0 = input.getValue(0, j);
                    v1 = input.getValue(1, j);
                    gX.setValue(0, j, (v1 - v0));
                } else {
                    // handle boundaries
                    v0 = input.getValue(i - 1, j);
                    v1 = input.getValue(w - 1, j);
                    gX.setValue(i, j, (v1 - v0));
                }
                if (j > 0 && (j + 1) < h) {
                    //gy[i][j] = img[i][j-1] - img[i][j+1]
                    v0 = input.getValue(i, j - 1);
                    v1 = input.getValue(i, j + 1);
                    gY.setValue(i, j, (v1 - v0));
                } else if (j == 0) {
                    // handle boundaries
                    v0 = input.getValue(i, 0);
                    v1 = input.getValue(i, 1);
                    gY.setValue(i, j, (v1 - v0));
                } else {
                    // handle boundaries
                    v0 = input.getValue(i, j - 1);
                    v1 = input.getValue(i, h - 1);
                    gY.setValue(i, j, (v1 - v0));
                }
            }
        }
        
        return new GreyscaleImage[]{gX, gY};
    }

    public void applySobelKernel(GreyscaleImage input) {

        float[] kernel = Gaussian1DFirstDeriv.getBinomialKernelSigmaZeroPointFive();

        GreyscaleImage[] gXY = createSobelGradients(input);
       
        GreyscaleImage img2 = combineConvolvedImages(gXY[0], gXY[1]);

        input.resetTo(img2);
    }

    /**
     * given a color image array with first dimension being color index
     * and the second dimension being the image pixel index,
     * apply the kernel [1,0,-1] to each pixel and combine the results
     * as SSD.
     * @param ptImg polar theta image of a color space such as
     * H of LCH that contains values between 0 and 255.
     * @param lowerDiff value in degrees for which a difference in
     * pixels results in a final value of "1".  For example,
     * 20 degrees.
     */
    public GreyscaleImage createBinary1stDerivForPolarTheta(
        GreyscaleImage ptImg, int lowerDiff) {

        int nPix = ptImg.getNPixels();
        int w = ptImg.getWidth();
        int h = ptImg.getHeight();

        GreyscaleImage out = ptImg.createWithDimensions();

        // kernel is .5, 0, -.5 so looking for difference in pixels on either
        //   side being .lte. lowerDiff
        int[] diffs = new int[4];
        int offset;
        int above;
        for (int i = 1; i < w - 1; ++i) {
            for (int j = 1; j < h - 1; ++j) {

                diffs[0] = ptImg.getValue(i - 1, j);
                diffs[1] = ptImg.getValue(i + 1, j);
                diffs[2] = ptImg.getValue(i, j - 1);
                diffs[3] = ptImg.getValue(i, j + 1);
                offset = 0;
                above = 0;
                for (int k = 0; k < 2; ++k) {
                    // in case there is wrap around, test adding a phase
                    //   and take the smaller of the results for each diff.
                    if (diffs[offset] > diffs[offset + 1]) {
                        // add a phase to next value if it's closer to current with addition
                        if ((diffs[offset] - diffs[offset + 1]) >
                            (diffs[offset + 1] + 255) - diffs[offset]) {
                            diffs[offset + 1] += 255;
                        }
                    } else if (diffs[offset + 1] > diffs[offset]) {
                        // add a phase to next value if it's closer to current with addition
                        if ((diffs[offset + 1] - diffs[offset]) >
                            (diffs[offset] + 255) - diffs[offset + 1]) {
                            diffs[offset] += 255;
                        }
                    }
                    int d = diffs[offset] - diffs[offset + 1];
                    if (Math.abs(d) >= lowerDiff) {
                        above = 1;
                        break;
                    }
                    offset += 2;
                }

                if (above == 1) {
                    out.setValue(i, j, 1);
                }
            }
        }

        return out;
    }

    /**
     * create  a float array from the image (the image is not scaled).
     * @param img
     * @return
     */
    public float[] copyToFloat(GreyscaleImage img) {

        float[] a = new float[img.getNPixels()];

        for (int i = 0; i < a.length; ++i) {
            a[i] = img.getValue(i);
        }

        return a;
    }

    /**
     * create  a int array from the image (the image is not scaled).
     * @param img
     * @return
     */
    public int[] copyToInt(GreyscaleImage img) {

        int[] a = new int[img.getNPixels()];

        for (int i = 0; i < a.length; ++i) {
            a[i] = img.getValue(i);
        }

        return a;
    }

    /**
     * using the binary results from createBinary1stDerivForPolarTheta
     * and the greyscale results from sobel operator,
     * scale the greyscale result so that the maximum value is 1.f,
     * add both and divide by 2.
     * NOTE: for other uses, may want to make a method which does not
     * scale the greyscale results or uses a different weighting
     * in the addition.
     *
     * @param gsImg
     * @param ptImg
     * @param lowerDiff
     * @return
     */
    public float[] createSobelColorScores(GreyscaleImage gsImg,
        GreyscaleImage ptImg, int lowerDiff) {

        int nPix = gsImg.getNPixels();
        int w = gsImg.getWidth();
        int h = gsImg.getHeight();

        if (ptImg.getWidth() != w || ptImg.getHeight() != h) {
            throw new IllegalArgumentException("images must be same size");
        }

        GreyscaleImage ptGrad = createBinary1stDerivForPolarTheta(
            ptImg, lowerDiff);

        float[] out = new float[nPix];
        for (int i = 0; i < gsImg.getNPixels(); ++i) {
            out[i] = gsImg.getValue(i);
        }

        out = createSobelConvolution(out, w, h);
        float maxV = MiscMath.findMax(out);
        float factor = 0.5f/maxV;
        int pixIdx;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                pixIdx = (j * w) + i;

                int v0 = ptGrad.getValue(pixIdx);
                float v = out[pixIdx] * factor;
                if (v0 == 1) {
                    v += 0.5f;
                }
                out[pixIdx] = v;
            }
        }

        return out;
    }

    /**
     * create a greyscale adaptive threshold gradient with canny algorithm
     * and then a color contrast gradient with "H" of LCH, and 1st deriv with
     * a threshold of 20 degrees for binarization, scale them to
     * 127 and add them.
     * The color binary pixels are scaled to 1/4th the maximum
     * of the greyscale gradient.
     *
     * The results could be improved in various ways, but for now
     * is a quick way to look at completing greyscale intensity
     * gradient contours with the color contrast gradient.
     *
     * @param img
     * @return
     */
    public GreyscaleImage createGradientWithColorAndGreyscale(Image img) {

        int w = img.getWidth();
        int h = img.getHeight();

        GreyscaleImage gsImg = img.copyToGreyscale2();

        CannyEdgeFilterAdaptive canny = new CannyEdgeFilterAdaptive();
        canny.overrideToUseAdaptiveThreshold();
        canny.overrideToNotUseLineThinner();
        canny.applyFilter(gsImg);
        EdgeFilterProducts prod = canny.getFilterProducts();

        float[] gsCanny = copyToFloat(prod.getGradientXY());

        GreyscaleImage scaled = MiscMath.rescaleAndCreateImage(gsCanny, w, h);

        //TODO: could consider using the sobel polar theta in canny edges
        //   as additional cues for strong edges in the 2-layer filter
        //   and then only use that result here.

        GreyscaleImage ptImg = createCIELUVTheta(img, 255);
        GreyscaleImage ptGrad =
            //createBinary2ndDerivForPolarTheta(ptImg, 20);
createBinary1stDerivForPolarTheta(ptImg, 20);

        /*
        ptGrad.multiply(255);
        applyAdaptiveMeanThresholding(ptGrad, 1);
        for (int j = 0; j < ptGrad.getNPixels(); ++j) {
            ptGrad.setValue(j, 255 - ptGrad.getValue(j));
        }*/

        float[] ptGradFloat = copyToFloat(ptGrad);

        for (int j = 0; j < ptGradFloat.length; ++j) {
            ptGradFloat[j] *= 63;
            ptGradFloat[j] += scaled.getValue(j);
        }

        scaled = MiscMath.rescaleAndCreateImage(ptGradFloat, w, h);

        return scaled;
    }

    /**
     * given a color image array with first dimension being color index
     * and the second dimension being the image pixel index,
     * apply the sobel kernel to each pixel and combine the results
     * as SSD.
     * @param colorInput with first dimension being color index
     * and the second dimension being the image pixel index
     * @param imgWidth
     * @return 
     */
    public float[] createSobelConvolution(float[][] colorInput, int imgWidth,
        int imgHeight) {

        int nClrs = colorInput.length;
        int nPix = colorInput[0].length;

        if (nPix != (imgWidth * imgHeight)) {
            throw new IllegalArgumentException("image width X height must equal "
                + " colorInput[0].length");
        }
        
        Kernel1DHelper kernelHelper = new Kernel1DHelper();
        float[] kernelG = new float[]{0.5f, 1.f, 0.5f};
        float[] kernel = new float[]{0.5f, 0, -0.5f};
                
        // apply kernelG to gX then kernel to gX
        // apply kernel to gY then kernelG
        // then combine both

        float[] outX = createConvolution(colorInput, imgWidth, imgHeight,
            kernelG, false);
        outX = createConvolution(outX, imgWidth, imgHeight,
            kernel, true);
        
        float[] outY = createConvolution(colorInput, imgWidth, imgHeight,
            kernel, true);
        outY = createConvolution(outY, imgWidth, imgHeight,
            kernelG, false);
         
        float[] out = new float[nPix];
        for (int i = 0; i < nPix; ++i) {
            float v = outX[i] * outX[i] + outY[i] * outY[i];
            v = (float)Math.sqrt(v/2.f);
            out[i] = v;
        }
        
        return out;
    }
    
    private float[] createConvolution(float[][] colorInput, int imgWidth,
        int imgHeight, float[] kernel, boolean calcForX) {

        int nClrs = colorInput.length;
        int nPix = colorInput[0].length;

        if (nPix != (imgWidth * imgHeight)) {
            throw new IllegalArgumentException("image width X height must equal "
                + " colorInput[0].length");
        }
        
        Kernel1DHelper kernelHelper = new Kernel1DHelper();
        
        float[] out = new float[nPix];

        for (int i = 0; i < imgWidth; ++i) {
            for (int j = 0; j < imgHeight; ++j) {

                double sum = 0;

                for (int c = 0; c < nClrs; ++c) {

                    if (calcForX) {
                        sum += kernelHelper.convolvePointWithKernel(
                            colorInput[c], i, j, kernel, true, imgWidth, imgHeight);
                    } else {
                        sum += kernelHelper.convolvePointWithKernel(
                            colorInput[c], i, j, kernel, false, imgWidth, imgHeight);
                    }
                }

                int pixIdx = (j * imgWidth) + i;

                out[pixIdx] = (float)(sum/(float)nClrs);
            }
        }

        return out;
    }
    
    private float[] createConvolution(float[] input, int imgWidth,
        int imgHeight, float[] kernel, boolean calcForX) {

        int nPix = input.length;

        if (nPix != (imgWidth * imgHeight)) {
            throw new IllegalArgumentException("image width X height must equal "
                + " input.length");
        }
        
        Kernel1DHelper kernelHelper = new Kernel1DHelper();
        
        float[] out = new float[nPix];

        for (int i = 0; i < imgWidth; ++i) {
            for (int j = 0; j < imgHeight; ++j) {

                float sum = 0;

                if (calcForX) {
                    sum += kernelHelper.convolvePointWithKernel(
                        input, i, j, kernel, true, imgWidth, imgHeight);
                } else {
                    sum += kernelHelper.convolvePointWithKernel(
                        input, i, j, kernel, false, imgWidth, imgHeight);
                }

                int pixIdx = (j * imgWidth) + i;

                out[pixIdx] = sum;
            }
        }

        return out;
    }

    /**
     * given a greyscale image
     * apply the kernel to each pixel and combine the results
     * as SSD.
     * @param greyscaleInput with index being the image pixel index
     */
    public float[] createSobelConvolution(float[] greyscaleInput, int imgWidth,
        int imgHeight) {

        int nPix = greyscaleInput.length;

        if (nPix != (imgWidth * imgHeight)) {
            throw new IllegalArgumentException("image width X height must equal "
                + " colorInput[0].length");
        }

        Kernel1DHelper kernelHelper = new Kernel1DHelper();
        float[] kernelG = new float[]{0.5f, 1.f, 0.5f};
        float[] kernel = new float[]{0.5f, 0, -0.5f};
         
        // apply kernelG to gX then kernel to gX
        // apply kernel to gY then kernelG
        // then combine both

        float[] outX = createConvolution(greyscaleInput, imgWidth, imgHeight,
            kernelG, false);
        outX = createConvolution(outX, imgWidth, imgHeight,
            kernel, true);
        
        float[] outY = createConvolution(greyscaleInput, imgWidth, imgHeight,
            kernel, true);
        outY = createConvolution(outY, imgWidth, imgHeight,
            kernelG, false);
         
        float[] out = new float[nPix];
        for (int i = 0; i < nPix; ++i) {
            float v = outX[i] * outX[i] + outY[i] * outY[i];
            v = (float)Math.sqrt(v/2.f);
            out[i] = v;
        }
        
        return out;
    }

    public void applySobelX(float[][] input) {

        float[] kernel = Gaussian1DFirstDeriv.getBinomialKernelSigmaZeroPointFive();

        /*
         1
         2  * [1 0 -1]
         1
        */
        applyKernel1D(input, new float[]{0.5f, 1.f, 0.5f}, false);
        
        applyKernel1D(input, kernel, true);
    }
    public void applySobelX(double[][] input) {

        double[] kernel = MiscMath0.convertFloatToDouble(Gaussian1DFirstDeriv.getBinomialKernelSigmaZeroPointFive());

        /*
         1
         2  * [1 0 -1]
         1
        */
        applyKernel1D(input, new double[]{0.5, 1., 0.5}, false);

        applyKernel1D(input, kernel, true);
    }

    public void applySobelY(float[][] input) {

        /*
         1
         0  * [1 2 1]
         -1
        */

        //0.5f, 0.0f, -0.5f
        float[] kernel = Gaussian1DFirstDeriv.getBinomialKernelSigmaZeroPointFive();

        applyKernel1D(input, kernel, false);
        
        applyKernel1D(input, new float[]{0.5f, 1.f, 0.5f}, true);
    }

    public void applySobelY(double[][] input) {

        /*
         1
         0  * [1 2 1]
         -1
        */

        //0.5f, 0.0f, -0.5f
        double[] kernel = MiscMath0.convertFloatToDouble(Gaussian1DFirstDeriv.getBinomialKernelSigmaZeroPointFive());

        applyKernel1D(input, kernel, false);

        applyKernel1D(input, new double[]{0.5, 1., 0.5}, true);
    }
    
    /**
     * calculate the 1st deriv gradient of the color image using CIELAB DeltaE 2000
     * and return gX, gY, and gXY with array indices being pixel
     * indexes of the image.
     * @param img
     * @return float[][]{gX, gY, gXY}
     */
    public float[][] calculateGradientUsingDeltaE2000(ImageExt img) {

        //TODO: consider making this sobel by applying the gaussian kernel too
        
        int n = img.getNPixels();

        CIEChromaticity cieC = new CIEChromaticity();

        int w = img.getWidth();
        int h = img.getHeight();

        float jnd = 2.3f;

        // using 1D 1st deriv kernel -1,0,1, calculating deltaE
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
            double gXY = Math.sqrt((outX[i] * outX[i] + outY[i] * outY[i])/2.f);
            outXY[i] = (float)gXY;
        }

        return new float[][]{outX, outY, outXY};
    }

    public GreyscaleImage applyLaplacianKernel(GreyscaleImage input,
        int minValue, int maxValue) {

        IKernel kernel = new Laplacian();
        Kernel kernelXY = kernel.getKernel();

        float norm = kernel.getNormalizationFactor();

        return applyKernel(input, kernelXY, norm, minValue, maxValue);
    }

    public Image combineConvolvedImages(Image imageX, Image imageY) {

        Image img2 = imageX.createWithDimensions();

        for (int i = 0; i < imageX.getWidth(); i++) {
            for (int j = 0; j < imageX.getHeight(); j++) {

                int rX = imageX.getR(i, j);
                int gX = imageX.getG(i, j);
                int bX = imageX.getB(i, j);

                int rY = imageY.getR(i, j);
                int gY = imageY.getG(i, j);
                int bY = imageY.getB(i, j);

                double r = Math.sqrt((rX*rX + rY*rY)/2.f);
                double g = Math.sqrt((gX*gX + gY*gY)/2.f);
                double b = Math.sqrt((bX*bX + bY*bY)/2.f);

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

                double g = Math.sqrt((gX*gX + gY*gY)/2.f);

                if (g > 255) {
                    g = 255;
                }

                img2.setValue(i, j, (int)g);
            }
        }

        return img2;
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
     * apply kernel to input. NOTE, that because the image is composed of vectors
     * that should have values between 0 and 255, inclusive, if the kernel application
     * results in a value outside of that range, the value is reset to 0 or
     * 255.
     * @param input
     * @param kernel
     * @param normFactor
     * @return the convolved image
     */
    protected GreyscaleImage applyKernel(GreyscaleImage input, Kernel kernel, 
        float normFactor, int minValue, int maxValue) {

        int h = (kernel.getWidth() - 1) >> 1;

        GreyscaleImage output = input.createFullRangeIntWithDimensions();

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
                
                if (v < minValue) {
                    v = minValue;
                } else if (v > maxValue) {
                    v = maxValue;
                }

                output.setValue(i, j, v);
            }
        }

        return output;
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
     * @param gradientX
     * @param gradientY
     * @return
     */
    public GreyscaleImage computeTheta180(final GreyscaleImage gradientX,
        final GreyscaleImage gradientY) {

        GreyscaleImage output = gradientX.createWithDimensions();

        for (int i = 0; i < gradientX.getWidth(); i++) {
            for (int j = 0; j < gradientX.getHeight(); j++) {

                double gX = gradientX.getValue(i, j);

                double gY = gradientY.getValue(i, j);

                if (gY < 0) {
                    gX *= -1;
                    gY *= -1;
                }

                double radians = Math.atan2(gY, gX);
                assert(radians >= 0. && radians <= 2.*Math.PI);
                
                int theta = (int)(radians * 180./Math.PI);
                if (theta == 180) {
                    theta = 0;
                }
                assert(theta >= 0 && theta < 180);
                
                output.setValue(i, j, theta);
            }
        }

        return output;
    }
    
    /**
     * calculate theta from the gradient x and y images and transform to
     * range 0 to 360.
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
    public GreyscaleImage computeTheta360_2(final GreyscaleImage convolvedX,
        final GreyscaleImage convolvedY) {

        GreyscaleImage output = convolvedX.createFullRangeIntWithDimensions();

        for (int i = 0; i < convolvedX.getWidth(); i++) {
            for (int j = 0; j < convolvedX.getHeight(); j++) {

                double gX = convolvedX.getValue(i, j);

                double gY = convolvedY.getValue(i, j);

                // -PI to PI then add PI
                double radians = Math.atan2(gY, gX) + Math.PI;

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

    public void applyKernelTwo1Ds(double[][] input, double[] kernel) {

        applyKernel1D(input, kernel, true);

        applyKernel1D(input, kernel, false);
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

    /**
     *
     * @param input
     * @param kernel
     * @param calcForX convolve along columns if true, else rows
     */
    public void applyKernel1D(float[][] input, float[] kernel,
        boolean calcForX) {

        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();

        int w = input.length;
        int h = input[0].length;

        float[][] output = new float[w][];

        for (int i = 0; i < w; i++) { // rows
            output[i] = new float[h];
            for (int j = 0; j < h; j++) { // cols
                //if calcX is true: convolution is over the column at the fixed row.
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

    /**
     *
     * @param input
     * @param kernel
     * @param calcForX convolve along columns if true, else rows
     */
    public void applyKernel1D(double[][] input, double[] kernel, boolean calcForX) {

        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();

        int w = input.length;
        int h = input[0].length;

        double[][] output = new double[w][];

        int i;
        int j;
        for (i = 0; i < w; i++) { // rows
            output[i] = new double[h];
            for (j = 0; j < h; j++) { // cols
                //if calcX is true: convolution is over the column at the fixed row.
                output[i][j] = kernel1DHelper.convolvePointWithKernel(input, i, j, kernel, calcForX);
            }
        }

        for (i = 0; i < w; i++) {
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

    /**
     * expand image to final size using the given output
     * widths and heights.
     * @param input
     * @param outWidth
     * @param outHeight
     * @param minValue minimum value to clamp results to
     * @param maxValue maximum value to clamp results to
     * @return
     */
    public GreyscaleImage upsampleUsingBilinear(GreyscaleImage input,
        int outWidth, int outHeight, int minValue, int maxValue) {

        if (input == null) {
            throw new IllegalArgumentException("input cannot be null");
        }

        int w0 = input.getWidth();
        int h0 = input.getHeight();
        int w2 = outWidth;
        int h2 = outHeight;

        if (w2 < w0 || h2 < h0) {
            throw new IllegalArgumentException("output dimensions cannot be"
                + " less than input dimensions for upsample");
        }

        /*
        example, 1D:
        output scale = 3.333

         0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
                |         |        |

          |  |  |
         0  1  2  3  4  5  6  7

        can solve using 2 different patterns:
        (1) integer upscale to Math.ceil(factor) == 4 in example
            then downsample to the final output size for easier math.
        or,
        (2) most pixels are integer copies of the current replicated pixel
        from the input, but the pixels which are integer muliples of the
        factor are composed of replicated current and next pixel
        as a fraction of sums to be computed for each.

        the first would be easier to maintain, but the later would be more
        efficient.

        will implement (2)
        */

        GreyscaleImage output = null;
        if (minValue >= 0 && maxValue <= 255) {
            output = input.createWithDimensions(w2, h2);
        } else {
            output = new GreyscaleImage(w2, h2,
                GreyscaleImage.Type.Bits32FullRangeInt);
        }

        float xFactor = (float)w2/(float)w0;
        float yFactor = (float)h2/(float)h0;

        // init vars
        int i2 = 0;
        int i2End = 0;

        for (int i0 = 0; i0 < w0; ++i0) {

            if (i0 > 0) {
                i2 = i2End + 1;
            }
            i2End = (int)Math.floor(i2 + xFactor);

            /* example, factor = 3.3
             0 : 2, 3 is combination
             4 : 6, 7 is combination
             7 : 9
            */

            float fractionX = (xFactor * (i0 + 1.f)) - i2End;
            if (fractionX < 0) {
                fractionX += 1.f;
            }

            int _i2End = (i2End < w2) ? i2End : w2;
            
            //System.out.format(
            //    "i2=%s i2End=%d,%d  frcX=%.3f\n", 
            //    i2, i2End, _i2End, fractionX);
            
            // init vars
            int j2 = 0;
            int j2End = 0;

            for (int j0 = 0; j0 < h0; ++j0) {

                // replication of integer pixels in this range:
                if (j0 > 0) {
                    j2 = j2End + 1;
                }
                j2End = (int)Math.floor(j2 + yFactor);

                int _j2End = (j2End < h2) ? j2End : h2;
                
                float fractionY = (yFactor * (j0 + 1.f)) - j2End;
                if (fractionY < 0) {
                    fractionY += 1.f;
                }
                
                //System.out.format(
                //    "  j2=%s j2End=%d,%d  frcY=%.3f\n", 
                //    j2, j2End, _j2End, fractionY);
                
                int v = input.getValue(i0, j0);
                v = (v < minValue) ? minValue : ((v > maxValue) ? maxValue : v);

                for (int ii = i2; ii < _i2End; ++ii) {
                    for (int jj = j2; jj < _j2End; ++jj) {
                
                        //System.out.format(
                        //"    set(%d,%d)=%d\n", ii, jj, v);
                
                        output.setValue(ii, jj, v);
                    }
                }

                // example, factor = 2.5
                //   i2 = 0, i2End = 2, (i1+1) * xFactor=2.5  -> fractional
                //           so pixel i2End gets contributions from i1=0 and i1=2
                //   i2 = 3, i2End = 5, (i1+1) * xFactor=5 -> integer, so skip

                // --- handle the pixels which get fractional contributions
                //     from 2 input pixels
                float vI2End = 0;
                if (fractionX > 0 && (i2End < w2)) {
                    vI2End = (fractionX * (float)v);
                    if ((i0 + 1) < w0) {
                        vI2End += ((1.f - fractionX) * 
                            (float)input.getValue(i0 + 1, j0));
                    }
                    int vI2 = Math.round(vI2End);
                    vI2 = (vI2 < minValue) ? minValue :
                        ((vI2 > maxValue) ? maxValue : vI2);
                    
                    for (int jj = j2; jj < _j2End; ++jj) {
                        
                        //System.out.format(
                        //"     *set(%d,%d)=%d\n", i2End, jj, vI2);
                        
                        output.setValue(i2End, jj, vI2);
                    }
                    
                }

                // last pixel is either a fractional sum of 2 input pixels or
                //  is start of next interval
                
                float vJ2End = 0;
                if (fractionY > 0 && (j2End < h2)) {
                    vJ2End = (fractionY * v);
                    if ((j0 + 1) < h0) {
                        vJ2End += ((1.f - fractionY) * input.getValue(i0, j0 + 1));
                    }
                    int vJ2 = Math.round(vJ2End);
                    vJ2 = (vJ2 < minValue) ? minValue :
                        ((vJ2 > maxValue) ? maxValue : vJ2);
                    
                    //System.out.format(
                    //"     **set(%d,%d)=%d\n", i2, j2End, vJ2);
                    
                    for (int iii = i2; iii < _i2End; ++iii) {
                        
                        //System.out.format(
                        // "     *set(%d,%d)=%d\n", iii, j2End, vJ2);
                        
                        output.setValue(iii, j2End, vJ2);
                    }                    
                }
                
                if (vI2End > 0 && vJ2End > 0) {
                    
                    int avg = Math.round((vI2End + vJ2End)/2.f);
                    
                    avg = (avg < minValue) ? minValue :
                        ((avg > maxValue) ? maxValue : avg);
                    
                    //System.out.format("     **set(%d,%d)=%d\n", i2End, j2End, avg);
                    
                    output.setValue(i2End, j2End, avg);
                }

                if (fractionY == 0.f) {
                    // subtract so next j2 starts at this value
                    j2End--;
                }
            }

            if (fractionX == 0.f) {
                // subtract so next i2 starts at this value
                i2End--;
            }
        }

        return output;
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
        
        FFTUtil fftUtil = new FFTUtil();

        // initialize matrix of complex numbers as real numbers from image (imaginary are 0's)
        Complex[][] cc = copyToComplex2D(input);

        Complex[][] ccOut = fftUtil.create2DFFT(cc, forward);

        input.fill(0);
        for (int col = 0; col < input.getWidth(); col++) {
            for (int row = 0; row < input.getHeight(); row++) {
                double re = ccOut[col][row].re();
                input.setValue(col, row, (int)re);
            }
        }

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
        
        FFTUtil fftUtil = new FFTUtil();

        Complex1D[] output2 = fftUtil.create2DFFT2(cInput, forward);

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
        Complex[][] cc = copyToComplexWithSwapMajor(input);

        FFTUtil fftUtil = new FFTUtil();
        
        Complex[][] ccFFT = fftUtil.create2DFFT(cc, doNormalize, forward);

        return ccFFT;
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
    protected Complex[][] copyToComplex2D(GreyscaleImage input) {

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
     * create a complex double array with format a[row][col]
     * @param input
     * @return
     */
    protected Complex[][] copyToComplexWithSwapMajor(GreyscaleImage input) {

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

    protected Complex[][] copyToComplex2D(double[][] input) {

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
     * NOT YET TESTED.  this assumes x and y should be interpreted as integers that
     * will be indexes for gsImg.
     *
     http://en.wikipedia.org/wiki/Bilinear_interpolation
     http://en.wikipedia.org/wiki/Bilinear_interpolation#/media/File:Bilinear_interpolation_visualisation.svg
     * @param x
     * @param y
     * @return
     */
    public double biLinearInterpolation(double[][] gsImg, double x, double y) {

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
     http://en.wikipedia.org/wiki/Bilinear_interpolation
     http://en.wikipedia.org/wiki/Bilinear_interpolation#/media/File:Bilinear_interpolation_visualisation.svg
     *
     * @param clrImg
     * @param x
     * @param y
     * @param outRGB r, g, b interpolated values.
     * @param auxArr1 auxiliary array of length 2 to use internally.  contents
     *    are written over and do not hold results for this method.
     * @param auxArr2 auxiliary array of size 2X2 to use internally.  contents
     *    are written over and do not hold results for this method.
     * @param auxArr3 auxiliary array of length 2 to use internally.  contents
     *    are written over and do not hold results for this method.
     */
    public void biLinearInterpolation(Image clrImg, double x, double y, double[] outRGB,
        double[] auxArr1, double[][] auxArr2, double[] auxArr3) {
        
        int x1 = (int)Math.floor(x);
        int x2 = (int)Math.ceil(x);
        int y1 = (int)Math.floor(y);
        int y2 = (int)Math.ceil(y);
        
        if ((x1 == x2) || (y1 == y2)) {
            biLinearInterpolation(clrImg, (float)x, (float)y, outRGB);
            return;
        }
                
        assert(x2 != x1);
        assert(y2 != y1);

        double d;
        auxArr1[0] = x2 - x;
        auxArr1[1] = x - x1;
        auxArr3[0] = y2 - y;
        auxArr3[1] = y - y1;
        d = 1./((x2-x1)*(y2-y1));
                        
        auxArr2[0][0] = clrImg.getR(x1, y1);
        auxArr2[0][1] = clrImg.getR(x1, y2);
        auxArr2[1][0] = clrImg.getR(x2, y1);
        auxArr2[1][1] = clrImg.getR(x2, y2);
        outRGB[0] = MatrixUtil.innerProduct(
            MatrixUtil.multiplyRowVectorByMatrix(auxArr1, auxArr2), auxArr3);
        outRGB[0] *= d;
        
        auxArr2[0][0] = clrImg.getG(x1, y1);
        auxArr2[0][1] = clrImg.getG(x1, y2);
        auxArr2[1][0] = clrImg.getG(x2, y1);
        auxArr2[1][1] = clrImg.getG(x2, y2);
        outRGB[1] = MatrixUtil.innerProduct(
            MatrixUtil.multiplyRowVectorByMatrix(auxArr1, auxArr2), auxArr3);
        outRGB[1] *= d;
        
        auxArr2[0][0] = clrImg.getB(x1, y1);
        auxArr2[0][1] = clrImg.getB(x1, y2);
        auxArr2[1][0] = clrImg.getB(x2, y1);
        auxArr2[1][1] = clrImg.getB(x2, y2);
        outRGB[2] = MatrixUtil.innerProduct(
            MatrixUtil.multiplyRowVectorByMatrix(auxArr1, auxArr2), auxArr3);
        outRGB[2] *= d;
    }
    
    /**
     http://en.wikipedia.org/wiki/Bilinear_interpolation
     http://en.wikipedia.org/wiki/Bilinear_interpolation#/media/File:Bilinear_interpolation_visualisation.svg
     *
     * @param x
     * @param y
     * @return
     */
    public double[] biLinearInterpolation(Image clrImg, float x, float y) {
        double[] rgb = new double[3];
        biLinearInterpolation(clrImg, x, y, rgb);
        return rgb;
    }

    /**
     http://en.wikipedia.org/wiki/Bilinear_interpolation
     http://en.wikipedia.org/wiki/Bilinear_interpolation#/media/File:Bilinear_interpolation_visualisation.svg
     *
     * @param x
     * @param y
     * @param outRGB
     */
    public void biLinearInterpolation(Image clrImg, float x, float y, double[] outRGB) {

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
                outRGB[0] = r1;
                outRGB[1] = g1;
                outRGB[2] = b1;
                return;
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
                outRGB[0] = r1;
                outRGB[1] = g1;
                outRGB[2] = b1;
                return;
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

        outRGB[0] = r;
        outRGB[1] = g;
        outRGB[2] = b;
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
     * @param m
     * @return
     */
    public double[][] copyToDouble2D(Matrix m) {

        double[][] output = new double[m.numRows()][];
        for (int i = 0; i < m.numRows(); ++i) {
            output[i] = new double[m.numColumns()];
            for (int j = 0; j < m.numColumns(); ++j) {
                output[i][j] = m.get(i, j);
            }
        }

        return output;
    }

    /**
     * @param input
     * @return
     */
    public float[][] copyToFloat2D(GreyscaleImage input) {

        int w = input.getWidth();
        int h = input.getHeight();
        
        float[][] output = new float[w][];
        for (int i = 0; i < w; ++i) {
            output[i] = new float[h];
            for (int j = 0; j < h; ++j) {
                output[i][j] = input.getValue(i, j);
            }
        }

        return output;
    }
    
    /**
     * @param input
     * @return
     */
    public int[][] copyToInt2D(GreyscaleImage input) {

        int w = input.getWidth();
        int h = input.getHeight();
        
        int[][] output = new int[w][];
        for (int i = 0; i < w; ++i) {
            output[i] = new int[h];
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

    /**
     *
     * @param image 2 dimensional array of data, assumed row-major format
     * @param sigma
     * @param outKeypoints0
     * @param outKeypoints1
     */
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
        float[][] firstDeriv = 
            MatrixUtil.add(tensorComponents.getDXSquared(),
                tensorComponents.getDYSquared());

        peakLocalMax(firstDeriv, 1, 0.1f, true, kp0, kp1);

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

    public Set<PairInt> binPoints(Set<PairInt> points, int binFactor) {

        Set<PairInt> out = new HashSet<PairInt>();
        for (PairInt p : points) {
            int x = Math.round((float)p.getX()/(float)binFactor);
            int y = Math.round((float)p.getY()/(float)binFactor);
            out.add(new PairInt(x, y));
        }

        return out;
    }

    /**
     * downsample the image to size w2,h2 with the given clamp values.
     * Note that the interpolation is currently a bilinear filter,
     * but in the future, a more accurate algorithm may be present.
     * possibly a windowed sinc pre-filter).
     * @param input
     * @param w2 the width of the output image
     * @param h2 the height of the output image
     * @param minValue minimum value to clamp results to
     * @param maxValue maximum value to clamp results to
     * @return
     */
    public GreyscaleImage downSample(GreyscaleImage input,
        int w2, int h2, int minValue, int maxValue) {

        // uses bilinear interpolation over pixels of input
        // contributing to down sampled output

        GreyscaleImage output = null;
        if (minValue >= 0 && maxValue <= 255) {
            output = input.createWithDimensions(w2, h2);
        } else {
            output = new GreyscaleImage(w2, h2,
                GreyscaleImage.Type.Bits32FullRangeInt);
        }

        int w0 = input.getWidth();
        int h0 = input.getHeight();

        float rW = (float)w0/(float)w2;
        float rH = (float)h0/(float)h2;

        int cX = Math.round(rW);
        int cY = Math.round(rH);

        //System.out.println("rX=" + rW + " rY=" + rH);

        for (int i = 0; i < w2; ++i) {

            float i2f = rW * i;

            for (int j = 0; j < h2; ++j) {

                double sum = 0;
                int np = 0;

                // integrate the points in input for offsets up to cX
                for (float ii = i2f; ii < (i2f + cX); ++ii) {
                    if (ii < 0 || Math.ceil(ii) >= w0) {
                        continue;
                    }

                    float j2f = rH * j;

                    // integrate the points in input for offsets up to cY
                    for (float jj = j2f; jj < (j2f + cY); ++jj) {
                        if (jj < 0 || Math.ceil(jj) >= h0) {
                            continue;
                        }

                        double v = biLinearInterpolation(input, ii, jj);

                        sum += v;

                        np++;
                    }
                }

                int v2 = (np > 0) ? (int)Math.round(sum/(float)np) : 0;

                if (v2 < 0) {
                    v2 = 0;
                } else if (v2 > 255) {
                    v2 = 255;
                }

                //System.out.format("(%d,%d) v==>%d\n", i, j, v2);

                output.setValue(i, j, v2);
            }
        }

        return output;
    }

    public int countNonZeroes(GreyscaleImage in) {
        int n = 0;
        for (int pixIdx = 0; pixIdx < in.getNPixels(); ++pixIdx) {
            if (in.getValue(pixIdx) > 0) {
                n++;
            }
        }
        return n;
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

    public int[] getAverageRGB(GreyscaleImage rImg, GreyscaleImage gImg,
        GreyscaleImage bImg, Collection<PairInt> points) {

        if (points.isEmpty()) {
            return null;
        }

        int rSum = 0;
        int gSum = 0;
        int bSum = 0;
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            rSum += rImg.getValue(x, y);
            gSum += gImg.getValue(x, y);
            bSum += bImg.getValue(x, y);
        }
        rSum /= points.size();
        gSum /= points.size();
        bSum /= points.size();

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
     * @param kernel
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
                    /*if (x == 0 && x2 == 0) {
                        continue;
                    }
                    if (y == 0 && y2 == 0) {
                        continue;
                    }
                    if (x == (w - 1) && x2 == (w - 1)) {
                        continue;
                    }
                    if (y == (h - 1) && y2 == (h - 1)) {
                        continue;
                    }*/

                    out.setValue(x2, y2, v);
                }
            }
        }

        return out;
    }

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
     * apply morphological thinning.
     * prefer this line thinner over applyThinning()
     * @param img
     * @return number of points > 0 after thinning
     */
    public int applyThinning2(GreyscaleImage img) {

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

        int n1s = 0;

        for (int i = 0; i < n0; ++i) {
            for (int j = 0; j < n1; ++j) {
                int m = skel[i][j];
                int v = img.getValue(i, j) * m;
                img.setValue(i, j, v);
                if (v > 0) {
                    n1s++;
                }
            }
        }

        return n1s;
    }

    /**
     * apply 8 hit or miss filters iteratively until convergence to thin the
     * image.  the operation is performed on all pixels with
     * value > 0.  prefer applyThinning2() to this method.
     * @param img
     */
    public void applyThinning(GreyscaleImage img) {

        //from https://en.wikipedia.org/wiki/Hit-or-miss_transform
        // and thinning

        // x,y pairs are sequential in these
        int[] c1 = new int[]{0,  0, -1, -1, 0,  -1, 1, -1};
        int[] d1 = new int[]{-1, 1, 0,   1, 1,  1};
        int[] c2 = new int[]{-1, 0, 0,   0, -1, -1, 0, -1};
        int[] d2 = new int[]{0,  1, 1,   1, 1,  0};

        int maxValue = img.getMax();
        int nBits = 1 + (int)Math.ceil(Math.log(maxValue)/Math.log(2));
        if (nBits > 31) {
            nBits = 31;
        }
        
        /*
            - - -        - -
              +        + + -
            + + +      + +
        */
        PairInt[][] neighborCoordOffsets
            = AbstractLineThinner.createCoordinatePointsForEightNeighbors(0, 0);

        int w = img.getWidth();
        int h = img.getHeight();
        int n = img.getNPixels();
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        int nEdited = 0;
        int nIter = 0;
        do {
            nEdited = 0;

            //GreyscaleImage tmp = out.copyImage();
            //tmp.multiply(255.f);
            //MiscDebug.writeImage(tmp, "_editing_" + MiscDebug.getCurrentTimeFormatted());

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

                    // to try to make it more symmetric, collecting all
                    // nullable pixels and counting the set neighbors,
                    // then revisiting by order of fewest set neighbors
                    MinHeapForRT2012 heap = new MinHeapForRT2012(9, n, nBits);

                    for (int x = 1; x < (img.getWidth() - 1); ++x) {
                        for (int y = 1; y < (img.getHeight() - 1); ++y) {
                            int v = img.getValue(x, y);

                            if (v == 0) {
                                continue;
                            }
                            if (allArePresent(img, x, y, tmpC)
                                && allAreNotPresent(img, x, y, tmpD)) {
                                if (!ImageSegmentation.doesDisconnect(img,
                                    neighborCoordOffsets, x, y)) {

                                    // number of neighbors that are not '1s
                                    int nn = 0;
                                    for (int k = 0; k < dxs.length; ++k) {
                                        int x2 = x + dxs[k];
                                        int y2 = y + dys[k];
                                        if (x2 < 0 || y2 < 0 || x2 >= w || y2 >= h) {
                                            continue;
                                        }
                                        if (img.getValue(x, y) > 0) {
                                            nn++;
                                        }
                                    }

                                    //long key = 8 - nn;
                                    long key = nn;
                                    HeapNode node = new HeapNode(key);
                                    int pixIdx = (y * w) + x;
                                    node.setData(Integer.valueOf(pixIdx));
                                    heap.insert(node);
                                }
                            }
                        }
                    }

                    while (heap.getNumberOfNodes() > 0) {

                        HeapNode node = heap.extractMin();

                        assert(node != null);

                        int pixIdx = ((Integer)node.getData()).intValue();
                        int y = pixIdx/w;
                        int x = pixIdx - (y * w);

                        int v = img.getValue(x, y);

                        if (v == 0) {
                            continue;
                        }
                        if (allArePresent(img, x, y, tmpC)
                            && allAreNotPresent(img, x, y, tmpD)) {
                            if (!ImageSegmentation.doesDisconnect(img,
                                neighborCoordOffsets, x, y)) {

                                img.setValue(x, y, 0);
                                nEdited++;
                            }
                        }
                    }
                }
            }
            nIter++;
        } while (nEdited > 0);
    }

    /**
     * apply 8 hit or miss filters iteratively until convergence to thin the
     * image.  the operation is performed on all pixels with value > 0.
     */
    public void applyThinning(Set<PairInt> points, int imageWidth, int imageHeight) {

        //from https://en.wikipedia.org/wiki/Hit-or-miss_transform
        // and thinning

        // x,y pairs are sequential in these
        int[] c1 = new int[]{0, 0, -1, -1, 0, -1, 1, -1};
        int[] d1 = new int[]{-1, 1, 0, 1, 1, 1};
        int[] c2 = new int[]{-1, 0, 0, 0, -1, -1, 0, -1};
        int[] d2 = new int[]{0, 1, 1, 1, 1, 0};

        int maxValue = Math.max(imageWidth, imageHeight);
        int nBits = 1 + (int)Math.ceil(Math.log(maxValue)/Math.log(2));
        if (nBits > 31) {
            nBits = 31;
        }
        
        /*
            - - -        - -
              +        + + -
            + + +      + +
        */
        PairInt[][] neighborCoordOffsets
            = AbstractLineThinner.createCoordinatePointsForEightNeighbors(
            0, 0);

        int w = imageWidth;
        int h = imageHeight;
        int n = points.size();
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        int nEdited = 0;
        int nIter = 0;
        do {
            nEdited = 0;

            //GreyscaleImage tmp = out.copyImage();
            //tmp.multiply(255.f);
            //MiscDebug.writeImage(tmp, "_editing_" + MiscDebug.getCurrentTimeFormatted());

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

                    // to try to make it more symmetric, collecting all
                    // nullable pixels and counting the set neighbors,
                    // then revisiting by order of fewest set neighbors
                    MinHeapForRT2012 heap = new MinHeapForRT2012(9, n, nBits);

                    for (PairInt p : points) {
                        int x = p.getX();
                        int y = p.getY();
                        if (allArePresent(points, x, y, tmpC)
                            && allAreNotPresent(points, x, y, tmpD)) {
                            if (!ImageSegmentation.doesDisconnect(points,
                                neighborCoordOffsets, x, y, imageWidth,
                                imageHeight)) {

                                // number of neighbors that are not '1s
                                int nn = 0;
                                for (int k = 0; k < dxs.length; ++k) {
                                    int x2 = x + dxs[k];
                                    int y2 = y + dys[k];
                                    if (x2 < 0 || y2 < 0 || x2 >= w || y2 >= h) {
                                        continue;
                                    }
                                    if (points.contains(new PairInt(x2, y2))) {
                                        nn++;
                                    }
                                }

                                //long key = 8 - nn;
                                long key = nn;
                                HeapNode node = new HeapNode(key);
                                int pixIdx = (y * w) + x;
                                node.setData(Integer.valueOf(pixIdx));
                                heap.insert(node);
                            }
                        }
                    }

                    while (heap.getNumberOfNodes() > 0) {

                        HeapNode node = heap.extractMin();

                        assert(node != null);

                        int pixIdx = ((Integer)node.getData()).intValue();
                        int y = pixIdx/w;
                        int x = pixIdx - (y * w);

                        if (allArePresent(points, x, y, tmpC)
                            && allAreNotPresent(points, x, y, tmpD)) {
                            if (!ImageSegmentation.doesDisconnect(points,
                                neighborCoordOffsets, x, y, imageWidth,
                                imageHeight)) {

                                points.remove(new PairInt(x, y));
                                nEdited++;
                            }
                        }
                    }
                }
            }
            nIter++;
        } while (nEdited > 0);
    }

    /**
     * apply 8 hit or miss filters iteratively until convergence to thin the
     * image.  the operation is performed on all pixels with value > 0.
     */
    public void applyThinning(TIntSet pixIdxs, int imageWidth, int imageHeight) {

        //from https://en.wikipedia.org/wiki/Hit-or-miss_transform
        // and thinning

        // x,y pairs are sequential in these
        int[] c1 = new int[]{0, 0, -1, -1, 0, -1, 1, -1};
        int[] d1 = new int[]{-1, 1, 0, 1, 1, 1};
        int[] c2 = new int[]{-1, 0, 0, 0, -1, -1, 0, -1};
        int[] d2 = new int[]{0, 1, 1, 1, 1, 0};

        int maxValue = Math.max(imageWidth, imageHeight);
        int nBits = 1 + (int)Math.ceil(Math.log(maxValue)/Math.log(2));
        if (nBits > 31) {
            nBits = 31;
        }
        
        /*
            - - -        - -
              +        + + -
            + + +      + +
        */
        PairInt[][] neighborCoordOffsets
            = AbstractLineThinner.createCoordinatePointsForEightNeighbors(
            0, 0);

        int w = imageWidth;
        int h = imageHeight;
        int n = pixIdxs.size();
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

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

                    // to try to make it more symmetric, collecting all
                    // nullable pixels and counting the set neighbors,
                    // then revisiting by order of fewest set neighbors
                    MinHeapForRT2012 heap = new MinHeapForRT2012(9, n, nBits);

                    TIntIterator iter = pixIdxs.iterator();
                    while (iter.hasNext()) {
                        int pixIdx = iter.next();
                        int y = pixIdx/w;
                        int x = pixIdx - (y * w);
                        if (allArePresent(pixIdxs, x, y, tmpC, w, h)
                            && allAreNotPresent(pixIdxs, x, y, tmpD, w, h)) {
                            if (!ImageSegmentation.doesDisconnect(pixIdxs,
                                neighborCoordOffsets, x, y, w, h)) {

                                // number of neighbors that are not '1s
                                int nn = 0;
                                for (int k = 0; k < dxs.length; ++k) {
                                    int x2 = x + dxs[k];
                                    int y2 = y + dys[k];
                                    if (x2 < 0 || y2 < 0 || x2 >= w || y2 >= h) {
                                        continue;
                                    }
                                    int pixIdx2 = (y2 * w) + x2;
                                    if (pixIdxs.contains(pixIdx2)) {
                                        nn++;
                                    }
                                }

                                //long key = 8 - nn;
                                long key = nn;
                                HeapNode node = new HeapNode(key);
                                node.setData(Integer.valueOf(pixIdx));
                                heap.insert(node);
                            }
                        }
                    }

                    while (heap.getNumberOfNodes() > 0) {

                        HeapNode node = heap.extractMin();

                        assert(node != null);

                        int pixIdx = ((Integer)node.getData()).intValue();
                        int y = pixIdx/w;
                        int x = pixIdx - (y * w);

                        if (allArePresent(pixIdxs, x, y, tmpC, w, h)
                            && allAreNotPresent(pixIdxs, x, y, tmpD, w, h)) {
                            if (!ImageSegmentation.doesDisconnect(pixIdxs,
                                neighborCoordOffsets, x, y, w, h)) {

                                pixIdxs.remove(pixIdx);
                                nEdited++;
                            }
                        }
                    }
                }
            }
            nIter++;
        } while (nEdited > 0);
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

    private boolean allArePresent(TIntSet pixIdxs, int x, int y, int[] xy,
        int imgWidth, int imgHeight) {

        for (int k = 0; k < xy.length; k += 2) {
            int tx = x + xy[k];
            int ty = y + xy[k + 1];
            if (tx < 0 || ty < 0 || (tx > (imgWidth - 1)) || (ty > (imgHeight - 1))) {
                continue;
            }
            int pixIdx2 = (ty * imgWidth) + tx;
            if (!pixIdxs.contains(pixIdx2)) {
                return false;
            }
        }

        return true;
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

    private boolean allAreNotPresent(TIntSet pixIdxs, int x, int y, int[] xy,
        int imgWidth, int imgHeight) {

        for (int k = 0; k < xy.length; k += 2) {
            int tx = x + xy[k];
            int ty = y + xy[k + 1];
            if (tx < 0 || ty < 0 || (tx > (imgWidth - 1)) || (ty > (imgHeight - 1))) {
                continue;
            }
            int pixIdx2 = (ty * imgWidth) + tx;
            if (pixIdxs.contains(pixIdx2)) {
                return false;
            }
        }

        return true;
    }

    /**
     * NOTE: this is not the same as the scale space image
     * or inflection points in the scaleSpace package.
     * It is a very quick look at the position derivatives.
     *
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
     * NOTE: this is not the same as the scale space image
     * or inflection points in the scaleSpace package.
     * It is a very quick look at the position derivatives.
     *
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

        // --- create 1st derivatives (gaussian 1st deriv sqrt(2)/2 = 0.707)----

        // switch X and Y 1st deriv operations for row major

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
     * * NOTE: this is not the same as the scale space image
     * or inflection points in the scaleSpace package.
     * It is a very quick look at the position derivative
     * zero crossings.
     *
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
        
        blur(img, SIGMA.ONE);

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
                curvature[i][j] = (
                    (dx[i][j] * dy2[i][j] - dy[i][j] * dx2[i][j])
                    / (dx2dx2 + dy2dy2));
                    /// Math.pow((dx2dx2 + dy2dy2), 1.5));
            }
        }

        peakLocalMax(curvature, 1, thresholdRel, true,
            outputKeypoints0,
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

        // --- create 1st derivatives (gaussian 1st deriv sqrt(2)/2 = 0.707)----

        // switch X and Y 1st deriv operations for row major

        float[] kernel = Gaussian1DFirstDeriv.getBinomialKernelSigmaZeroPointFive();        
        
        float[][] gX = copy(image);
        applyKernel1D(gX, kernel, true);

        float[][] gY = copy(image);
        applyKernel1D(gY, kernel, false);

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

    /**
     * calculate Harris Corner responses for each pixel in the image then use
     * local maximum thresholding with a minimum separation between maxima, to find the peaks.
     *
     * @param img
     * @return image pairs in row major format. For example:
     *   keypoint[0] = int[]{row_i0, col_i0};
     *   keypoint[1] = int[]{row_i1, col_i1};
     *   ...
     */
    public int[][] calcHarrisCorners(double[][] img) {
        int minDist = 1;
        float thresholdRel = 0.1f;
        boolean ignore0sInThresh = true;
        return calcHarrisCorners(img, minDist, thresholdRel, ignore0sInThresh);
    }

    /**
     * calculate the Harris Corners with local max threshholding and return a stacked array of pairs
     * of (row, col) for each corner.
     * @param img
     * @param minDist
     * @param thresholdRel
     * @param ignore0sInThresh
     * @return
     */
    public static int[][] calcHarrisCorners(double[][] img, int minDist, float thresholdRel,
                                            boolean ignore0sInThresh) {

        float sigmaBlur = 1;
        StructureTensorD tensorComponents = new
                StructureTensorD(img, sigmaBlur, false);

        double[][] detA = tensorComponents.getDeterminant();

        double[][] traceA = tensorComponents.getTrace();

        double[][] hr = ORB.cornerHarris(img);

        //System.out.printf("harris response=\n%s\n", FormatArray.toString(hr, "%.1f"));

        ImageProcessor imageProcessor = new ImageProcessor();
        // keypoints0 are the row numbers
        TIntList keypoints0 = new TIntArrayList();
        // keypoints1 are the column numbers
        TIntList keypoints1 = new TIntArrayList();

        imageProcessor.peakLocalMax(MatrixUtil.convertToFloat(hr),
                minDist, thresholdRel,
                ignore0sInThresh, keypoints0, keypoints1);

        int[][] keypoints = new int[keypoints0.size()][];
        for (int i = 0; i < keypoints0.size(); ++i) {
            keypoints[i] = new int[]{keypoints0.get(i), keypoints1.get(i)};
        }

        return keypoints;
    }
    
    public TIntSet convertPointsToIndexes(Set<PairInt> points, int width) {
        TIntSet set = new TIntHashSet(points.size());
        for (PairInt p : points) {
            int pixIdx = (p.getY() * width) + p.getX();
            set.add(pixIdx);
        }
        return set;
    }

    public Set<PairInt> convertIndexesToPoints(TIntSet pixIdxs, int width) {
        PixelHelper ph = new PixelHelper();
        int[] xyout = new int[2];      
        Set<PairInt> set = new HashSet<PairInt>();
        TIntIterator iter = pixIdxs.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            ph.toPixelCoords(pixIdx, width, xyout);
            int y = xyout[1];
            int x = xyout[0];
            set.add(new PairInt(x, y));
        }
        return set;
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

    NOTE also that the method has been adapted to include a flag to treat 0's as
    masked out numbers and have made changes to handle negative numbers in the 
    image. Previously, the scipy method ignored all negative numbers so images
    of only negative previously returned no peaks.
    
    @param img output image 
    @param minDistance
        Minimum number of pixels separating peaks in a region of `2 *
        min_distance + 1` (i.e. peaks are separated by at least
        `min_distance`).
        To find the maximum number of peaks, use `min_distance=1`.
     * @param thresholdRel a positive real number that becomes a factor to
     * apply to the maximum value.
     * @param ignore0sInThreshold if true, method handles negative numbers,
     * else expects only the non-negative numbers should be used.
      @param outputKeypoints0 the output row coordinates of keypoints
     * @param outputKeypoints1 the output col coordinates of keypoints
     */
    public void peakLocalMax(float[][] img, int minDistance,
        float thresholdRel, boolean ignore0sInThreshold,
        TIntList outputKeypoints0, TIntList outputKeypoints1) {
        
        int excludeBorder = minDistance;
        int numPeaks = Integer.MAX_VALUE;
        //int numPeaksPerLabel = Integer.MAX_VALUE;

        float sentinel = -10.f * Float.MIN_VALUE;
        
        /*
        The peak local maximum function returns the coordinates of local peaks
        (maxima) in an image. A maximum filter is used for finding local maxima.
        This operation dilates the original image. After comparison of the dilated
        and original image, this function returns the coordinates or a mask of the
        peaks where the dilated image equals the original image.
        */

        int nRows = img.length;
        int nCols = img[0].length;
        
        // if "ignore zeroes" is set, need to make those negative numbers
        // so that the maximum filter won't find those as maximum if all
        // significant pixels are the nonzeroes and are negative numbers.
        float[][] origImg = null;
        if (ignore0sInThreshold) {
            origImg = copy(img);
            for (int i = 0; i < img.length; ++i) {
                for (int j = 0; j < img[0].length; ++j) {
                    if (img[i][j] == 0.f) {
                        img[i][j] = sentinel;
                    }
                }
            }
        }

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

        if (ignore0sInThreshold) {
            for (int i = 0; i < imageMax.length; ++i) {
                for (int j = 0; j < imageMax[0].length; ++j) {
                    if (img[i][j] == sentinel) {
                        imageMax[i][j] = sentinel;
                    }
                }
            }
        }
        
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
        
        // for the ORB images such as harris corner response images,
        // the zeroes are where there are no values sometimes, so
        // if user has set flag ignore0sInThreshold, ignore them
        // in determining the max
        float thresholdMax;
        if (ignore0sInThreshold) {
            thresholdMax = sentinel;
            for (int i = 0; i < img.length; ++i) {
                for (int j = 0; j < img[0].length; ++j) {
                    float v = img[i][j];
                    if (v > thresholdMax) {
                        thresholdMax = v;
                    }
                }
            }
           
            if (thresholdMax <= 0.f) {
                float delta = thresholdMax * thresholdRel;
                thresholdMax += delta; 
            } else {
                thresholdMax *= thresholdRel;
            }
        } else {
            float mx = MiscMath.findMax(img);
            if (mx <= 0.f) {
                float delta = mx * thresholdRel;
                thresholdMax = mx + delta;
            } else {
                thresholdMax = thresholdRel * mx;
            }
        }
        
        if (ignore0sInThreshold) {
            if (thresholdAbs == 0.0f) {
                thresholdAbs = thresholdMax;
            } else {
                thresholdAbs = Math.max(thresholdAbs, thresholdMax);
            }
        } else {
            thresholdAbs = Math.max(thresholdAbs, thresholdMax);
        }
        
        // mask &= image > 0.1
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                float v = imageMax[i][j];
                if (v >= thresholdAbs) {
                    mask[i][j] &= 1;
                } else {
                    mask[i][j] = 0;
                }
            }
        }

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
        
        if (ignore0sInThreshold) {
            for (int i = 0; i < img.length; ++i) {
                for (int j = 0; j < img[0].length; ++j) {
                    if (img[i][j] == sentinel && origImg[i][j] == 0.f) {
                        img[i][j] = 0.0f;
                    }
                }
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
        return Filters.maximumFilter(img, size);
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

            gsR = buildPyramid2(r, 32);
            gsG = buildPyramid2(g, 32);
            gsB = buildPyramid2(b, 32);

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

            out = buildPyramid2(cp, 32);

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
     * @return
     */
    public List<GreyscaleImage> buildPyramid2(GreyscaleImage input,
        int decimationLimit) {

        List<GreyscaleImage> output = new ArrayList<GreyscaleImage>();

        MedianTransform mt = new MedianTransform();
        mt.multiscalePyramidalMedianTransform2(input, output, decimationLimit);

        if (output.size() == 1) {
            return output;
        }
        
        float factor = 1.2f;//1.5f;

        List<GreyscaleImage> output2 = new ArrayList<GreyscaleImage>();

        // add an image in between each after output[2]
        for (int i = 0; i < output.size() - 1; ++i) {
            output2.add(mt.decimateImage(output.get(i), factor, 0, 255));
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

    /**
     * convert the image to cie l*a*b* and then use a and b
     * to calculate polar angle around 0 in degrees.
     * If maxV of 360, returns full value image,
     * else if is 255, scales the values to max value of 255, etc.
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
     * create sobel gradient images for the given images specified by idxs
     * @param images
     * @param idxs
     * @return
     */
    public GreyscaleImage[] createSobels(GreyscaleImage[] images, int[] idxs) {

        GreyscaleImage[] sobels = new GreyscaleImage[2];

        int count = 0;
        for (int idx : idxs) {
            GreyscaleImage img2 = images[idx];
            sobels[count] = img2.copyImage();
            applySobelKernel(sobels[count]);
            count++;
        }
        return sobels;
    }

    public GreyscaleImage[] createSobelLCForLUV(Image img) {

        GreyscaleImage[] lch = createLCHForLUV(img);

        GreyscaleImage[] sobels = createSobels(lch, new int[]{0, 1});

        return sobels;
    }

    /**
     * create 3 images of LCH where L, C, and H are the
     * luminosity, magnitude, and polar angle of CIE LUV color space, respectively.
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
     * around 0 in degrees (a.k.a. the "H" of LCH color space, but with
     * the 1976 CIE LAB which is LUV).
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
     * convert the image to cie luv and then calculate polar angle of u and v
     * around 0 in degrees (a.k.a. the "H" of LCH color space, but with
     * the 1976 CIE LAB which is LUV).
     * If maxV of 360, returns full value image,
     * else if is 255, scales the values to max value of 255, etc.
     * @param img
     * @param maxV
     * @return
     */
    public GreyscaleImage createCIELUVThetaBilinear(Image img, int maxV) {

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
     * uses standard illuminant "wide range lightness".
     * Wide-range Lightness data were generated by Fairchild et al., 
     * who conducted two different experiments to scale lightness above and 
     * below diffuse white (CIE * L = 100 ). 
     * In the Scaling Lightness Experiment 1 (SL1) they used a luminance
       range from 156 to 2 3692cd m with 2 842 Y cd m n = 
     (Yn represents the luminance of reference white) 
     * whereas in the Scaling Lightness Experiment 2 (SL2) the
      luminance range was extended from 0 to 2 7432cd m with 2 997 Y cd m n = . 
     The SL2 data set was used to drive the adapted lightness ( z J ) 
     formula of the proposed color space (see later) while the SL1 data set was 
     used as a test data set. Each of the sets includes 19 samples. 

     * https://www.osapublishing.org/DirectPDFAccess/B810E9AE-C7C9-E594-C72DC7FBE1424F0A_368272/oe-25-13-15131.pdf?da=1&id=368272&seq=0&mobile=no
     * 
     * SL2 Training D65/2° (x,y,z)=968.08 997 883.51 
     * L_a=199 
     * C=0.69. N_C=1, F=1 
     * 
     * convert the image to cie luv and then calculate polar angle of u and v
     * around 0 in degrees (a.k.a. the "H" of LCH color space, but with
     * the 1976 CIE LAB which is LUV).
     * If maxV of 360, returns full value image,
     * else if is 255, scales the values to max value of 255, etc.
     * @param img
     * @param maxV
     * @return
     */
    public GreyscaleImage createCIELUVTheta_WideRangeLightness(Image img, int maxV) {

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

        float xd = 96.808f;
        float yd = 99.7f;
        float zd = 88.351f;
        int n = img.getNPixels();
        int v;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {

                int r = img.getR(i, j);
                int g = img.getG(i, j);
                int b = img.getB(i, j);

                float[] lma = cieC.rgbToPolarCIELUV(r, g, b, xd, yd, zd);

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

    /**
     * 
     * @param pointIndexMap map with key=pixelIndex, value=label.
     * @param labeledPoints list of point sets that are the edges of blobs.  NOTE
     * that the indexes of the list are the labels and they are the same labels
     * as present in pointIndexMap.  NOTE that pointIndexMap can have more 
     * points than are present in labeledPoints, but the one context that this
     * method was created for has same points in both data structures.
     * @param imgWidth
     * @param imgHeight
     * @return a map with key=label, value = bitstring marking the neighboring labels.
     */
    public TIntObjectMap<VeryLongBitString> createAdjacencyMap(
        TIntIntMap pointIndexMap, List<TIntSet> labeledPoints,
        int imgWidth, int imgHeight) {

        int n = labeledPoints.size();

        TIntObjectMap<VeryLongBitString> output
            = new TIntObjectHashMap<VeryLongBitString>();

        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;

        // making an adjacency list for the labels
        for (int label = 0; label < n; ++label) {

            TIntSet pixIdxs = labeledPoints.get(label);

            VeryLongBitString nbrs = new VeryLongBitString(n);
            
            TIntIterator iter2 = pixIdxs.iterator();
            while (iter2.hasNext()) {
                int pixIdx = iter2.next();
                int y = pixIdx/imgWidth;
                int x = pixIdx - (y * imgWidth);
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    if (x2 < 0 || y2 < 0 || x2 >= imgWidth || y2 >= imgHeight) {
                        continue;
                    }
                    int pixIdx2 = (y2 * imgWidth) + x2;
                    if (pointIndexMap.containsKey(pixIdx2)) {
                        int label2 = pointIndexMap.get(pixIdx2);
                        if (label2 != label) {
                            nbrs.setBit(label2);
                        }
                    }
                }
            }
            output.put(label, nbrs);
        }
        
        assert(output.size() == labeledPoints.size());

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

    /**
     * a utility method for polar theta images that have had their range
     * scaled from 0 to 255.  This method copies the image and applies
     * the valueShift amount to each pixel value and if result is larger than
     * 255, subtracts 255 to wrap around to 0 again in the polar coordinates.
     * NOTE: the algorithms that use this to compensate for the discontinuity
     * from 255 to 0 can probably be changed to use a difference in angles.
     * @param polarImage
     * @param valueShift
     * @return 
     */
    public GreyscaleImage copyAndShiftPolarAngleImage(GreyscaleImage polarImage, int valueShift) {

        GreyscaleImage shiftedImg = polarImage.createWithDimensions();
        for (int i = 0; i < polarImage.getNPixels(); ++i) {
            int v = polarImage.getValue(i);
            v += valueShift;
            if (v < 0) {
                v += 255;
            } else if (v > 255) {
                v -= 255;
            }
            shiftedImg.setValue(i, v);
        }

        return shiftedImg;
    }
    
    // TODO: implement the methods in
    // http://www.merl.com/publications/docs/TR2008-030.pdf
    // for an O(n) filter.
    // "Constant Time O(1) Bilateral Filtering" by Porikli
    //public void applyBiLateralFilter(Image img) {
    //}

    // and trilateral filter by Tumblin et al. 2003

}
