package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.EllipseHelper;
import algorithms.compGeometry.PerimeterFinder;
import algorithms.compGeometry.clustering.KMeansPlusPlus;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.Errors;
import algorithms.util.PairFloat;
import algorithms.util.PairInt;
import algorithms.util.PolynomialFitter;
import java.awt.Color;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class ImageProcessor {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());

    public void applySobelKernel(Image input) {
        
        IKernel kernel = new SobelX();
        Kernel kernelX = kernel.getKernel();
        
        float normX = kernel.getNormalizationFactor();
        
        kernel = new SobelY();
        Kernel kernelY = kernel.getKernel();
        
        float normY = kernel.getNormalizationFactor();
       
        applyKernels(input, kernelX, kernelY, normX, normY);
    }
    
    
    public void applySobelKernel(GreyscaleImage input) {
        
        IKernel kernel = new SobelX();
        Kernel kernelX = kernel.getKernel();
        
        float normX = kernel.getNormalizationFactor();
        
        kernel = new SobelY();
        Kernel kernelY = kernel.getKernel();
        
        float normY = kernel.getNormalizationFactor();
       
        applyKernels(input, kernelX, kernelY, normX, normY);
    }
    
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
                
        /*
        assumes that kernelX is applied to a copy of the img
        and kernelY is applied to a separate copy of the img and
        then they are added in quadrature for the final result.
        */
        
        GreyscaleImage imgX = input.copyImage();
        
        GreyscaleImage imgY = input.copyImage();
        
        applyKernel(imgX, kernelX, normFactorX);
        
        applyKernel(imgY, kernelY, normFactorY);
        
        GreyscaleImage img2 = combineConvolvedImages(imgX, imgY);
        
        input.resetTo(img2);
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
     * process only the green channel and set red and blue to zero
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
                
                //double g = Math.sqrt(0.5*(gX*gX + gY*gY));
                
                //g = (g > 255) ? 255 : g;
                
                double g = Math.sqrt(gX*gX + gY*gY);
                
                if (g > 255) {
                    g = 255;
                }
                                    
                //int rgb = (int)(((rSum & 0x0ff) << 16) 
                //    | ((gSum & 0x0ff) << 8) | (bSum & 0x0ff));
                 
                img2.setValue(i, j, (int)g);
            }
        }
        
        return img2;
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
                
                long value = 0;
                
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
               
                if (value < 0) {
                    value = 0;
                }
                if (value > 255) {
                    value = 255;
                }
                output.setValue(i, j, (int)value);
            }
        }
        
        input.resetTo(output);
    }
 
    public Image computeTheta(Image convolvedX, Image convolvedY) {

        Image output = new Image(convolvedX.getWidth(), convolvedX.getHeight());
        
        for (int i = 0; i < convolvedX.getWidth(); i++) {
            for (int j = 0; j < convolvedX.getHeight(); j++) {
                
                double rX = convolvedX.getR(i, j);
                double gX = convolvedX.getG(i, j);
                double bX = convolvedX.getB(i, j);
                
                double rY = convolvedY.getR(i, j);
                double gY = convolvedY.getG(i, j);
                double bY = convolvedY.getB(i, j);
                
                int thetaR = calculateTheta(rX, rY);
                int thetaG = calculateTheta(gX, gY);
                int thetaB = calculateTheta(bX, bY);
                
                output.setRGB(i, j, thetaR, thetaG, thetaB);
            }
        }
        
        return output;
    }
    
    public GreyscaleImage computeTheta(final GreyscaleImage convolvedX, 
        final GreyscaleImage convolvedY) {

        GreyscaleImage output = convolvedX.createWithDimensions();
        
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
    
    public GreyscaleImage subtractImages(final GreyscaleImage image, 
        final GreyscaleImage subtrImage) {
        
        if (image.getWidth() != subtrImage.getWidth()) {
            throw new IllegalArgumentException("image widths must be the same");
        }
        if (image.getHeight() != subtrImage.getHeight()) {
            throw new IllegalArgumentException("image heights must be the same");
        }
        
        GreyscaleImage output = image.createWithDimensions();
        
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
        
        double div = gradientY/gradientX;
                
        double theta = Math.atan(div)*180./Math.PI;
        
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
     * 
     * @param input
     * @return 
     */
    public int[] shrinkImageToFirstNonZeros(final GreyscaleImage input) {
        
        int xNZFirst = -1;
        int xNZLast = -1;
        int yNZFirst = -1;
        int yNZLast = -1;
        
        for (int i = 0; i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {
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
        for (int j = 0; j < input.getHeight(); j++) {
            for (int i = 0; i < input.getWidth(); i++) {
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
        
        for (int i = (input.getWidth() - 1); i > -1; i--) {
            for (int j = (input.getHeight() - 1); j > -1; j--) {
                if (input.getValue(i, j) > 0) {
                    xNZLast = i;                    
                    break;
                }
            }
            if (xNZLast > -1) {
                break;
            }
        }
        
        for (int j = (input.getHeight() - 1); j > -1; j--) {
            for (int i = (input.getWidth() - 1); i > -1; i--) {
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
        
        if ((xNZFirst > 0) || (xNZLast < (input.getWidth() - 1)) 
            || (yNZFirst > 0) || (yNZLast < (input.getHeight() - 1))) {
            
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
            
            GreyscaleImage output = new GreyscaleImage(xLen, yLen);
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
    
        //xOffset, yOffset, width, height
        GreyscaleImage output = new GreyscaleImage(offsetsAndDimensions[2], 
            offsetsAndDimensions[3]);
        output.setXRelativeOffset(offsetsAndDimensions[0]);
        output.setYRelativeOffset(offsetsAndDimensions[1]);
        
        int x = 0;
        
        for (int col = offsetsAndDimensions[0]; col < offsetsAndDimensions[2];
            col++) {
            
            int y = 0;
            for (int row = offsetsAndDimensions[1]; row < offsetsAndDimensions[3];
                row++) {
                
                int v = input.getValue(col, row);
                
                output.setValue(x, y, v);
                
                y++;
            }
            
            x++;
        }
    }        

    public void applyImageSegmentation(GreyscaleImage input, int kBands) 
        throws IOException, NoSuchAlgorithmException {
                
        KMeansPlusPlus instance = new KMeansPlusPlus();
        instance.computeMeans(kBands, input);
        
        int[] binCenters = instance.getCenters();
        
        for (int col = 0; col < input.getWidth(); col++) {
            
            for (int row = 0; row < input.getHeight(); row++) {
                
                int v = input.getValue(col, row);
                                
                for (int i = 0; i < binCenters.length; i++) {
                    
                    int vc = binCenters[i];
                  
                    int bisectorBelow = ((i - 1) > -1) ? 
                        ((binCenters[i - 1] + vc) / 2) : 0;

                    int bisectorAbove = ((i + 1) > (binCenters.length - 1)) ? 
                        255 : ((binCenters[i + 1] + vc) / 2);
                    
                    if ((v >= bisectorBelow) && (v <= bisectorAbove)) {
                        
                        input.setValue(col, row, vc);
                        
                        break;
                    }
                }
            }
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
     * using the gradient's theta image, find the sky as the largest set of
     * contiguous 0 values and apply the edge filter to it to reduce the
     * boundary to a single pixel curve.  
     * 
     * NOTE that the theta image has a boundary that has been increased by 
     * the original image blur and then the difference of gaussians to make 
     * the gradient, so the distance of the skyline from the real image horizon 
     * is several pixels.
     * For example, the canny edge filter used in "outdoorMode" results in
     * gaussian kernels applied twice to give an effective sigma of 
     * sqrt(2*2 + 0.5*0.5) = 2.1.  The FWHM of such a spread is then 
     * 2.355*2.1 = 6 pixels.   The theta image skyline is probably blurred to a 
     * width larger than the combined FWHM, however, making it 7 or 8 pixels.  
     * Therefore, it's recommended that the image returned from this be followed 
     * with: edge extraction; then fit the edges to the intermediate canny edge 
     * filter product (the output of the 2 layer filter) by making a translation
     * of the extracted skyline edge until the correlation with the filter2
     * image is highest (should be within 10 pixel shift).  Because this
     * method does not assume orientation of the image, the invoker needs
     * to also retrieve the centroid of the sky, so that is also returned 
     * in an output variable given in the arguments.
     * 
     * @param theta
     * @param gradientXY
     * @param originalImage
     * @param outputSkyCentroid container to hold the output centroid of 
     * the sky.
     * @param edgeSettings
     * @return
     * @throws IOException
     * @throws NoSuchAlgorithmException 
     */
    public GreyscaleImage createSkyline(GreyscaleImage theta, 
        GreyscaleImage gradientXY, Image originalImage,
        CannyEdgeFilterSettings edgeSettings, PairIntArray outputSkyCentroid) 
        throws IOException, NoSuchAlgorithmException {        
      
        GreyscaleImage mask = createBestSkyMask(theta, gradientXY, originalImage, 
            edgeSettings, outputSkyCentroid);
        
        if (mask != null) {
            
            multiply(mask, 255);
            
            CannyEdgeFilter filter = new CannyEdgeFilter();
            
            filter.setFilterImageTrim(theta.getXRelativeOffset(), 
                theta.getYRelativeOffset(), theta.getWidth(), 
                theta.getHeight());
            
            filter.applyFilter(mask);
            
            return mask;
        }
        
        return null;
    }
    
    /**
     * NOT READY FOR USE YET.
     * 
     * @param theta
     * @return 
     */
    public GreyscaleImage createRoughSkyMask(GreyscaleImage theta) throws 
        IOException, NoSuchAlgorithmException {
        
        if (theta == null) {
            throw new IllegalArgumentException("theta cannot be null");
        }
        
        theta = theta.copyImage();
        
        applyImageSegmentation(theta, 2);
        
        subtractMinimum(theta);
        
        convertToBinaryImage(theta);
        
        removeSpurs(theta);
        
        throw new UnsupportedOperationException("not ready for use yet");
        //return theta;
    }

    protected int determineBinFactorForSkyMask(int numberOfThetaPixels) {
        
        //TODO: this can be adjusted by the jvm settings for stack size
        int defaultLimit = 87000;
        
        if (numberOfThetaPixels <= defaultLimit) {
            return 1;
        }
                    
        double a = (double)numberOfThetaPixels/87000.;
        // rounds down
        int f2 = (int)a/2;
        int binFactor = f2 * 2;
        if ((a - binFactor) > 0) {
            binFactor += 2;
        }
            
        return binFactor;
    }
    
    /**
     * NOT READY FOR USE
     * 
     * create a mask for what is interpreted as sky in the image and return
     * a mask with 0's for sky and 1's for non-sky.
     * 
     * Internally, the method looks at contiguous regions of zero value pixels 
     * in the theta image and it looks at color in the original image.  
     * The camera image plane can have a rotation such that the 
     * horizon might not be along rows in the image, that is the method
     * looks for sky that is not necessarily at the top of the image, but 
     * should be on the boundary of the image
     * (note that reflection of sky can be found the same way, but this
     * method does not try to find sky that is not on the boundary of the
     * image).
     * (the "boundary" logic may change, in progress...)
     * 
     * Note that if the image contains a silhouette of featureless
     * foreground and a sky full of clouds, the method will interpret the
     * foreground as sky so it is up to the invoker to invert the mask.
     * 
     * NOTE: the cloud finding logic will currently fail if originalColorImage
     * is black and white.
     * 
     * @param theta
     * @param gradientXY
     * @param originalColorImage
     * @param edgeSettings
     * @param outputSkyCentroid container to hold the output centroid of 
     * the sky.
     * @return 
     * @throws java.io.IOException 
     * @throws java.security.NoSuchAlgorithmException 
     */
    public GreyscaleImage createBestSkyMask(final GreyscaleImage theta,
        GreyscaleImage gradientXY,
        Image originalColorImage, CannyEdgeFilterSettings edgeSettings,
        PairIntArray outputSkyCentroid) throws 
        IOException, NoSuchAlgorithmException {
        
        if (theta == null) {
            throw new IllegalArgumentException("theta cannot be null");
        }
                
        Image colorImg = originalColorImage;
        GreyscaleImage thetaImg = theta;
        GreyscaleImage gXYImg = gradientXY;
        
        int binFactor = determineBinFactorForSkyMask(theta.getNPixels());

        log.info("binFactor=" + binFactor);

        if (binFactor > 1) {
            thetaImg = binImage(theta, binFactor);
            colorImg = binImage(originalColorImage, binFactor);
            gXYImg = binImage(gradientXY, binFactor);
        }

        List<PairIntArray> zeroPointLists = getSortedContiguousZeros(thetaImg);
        
        if (zeroPointLists.isEmpty()) {
            
            GreyscaleImage mask = theta.createWithDimensions();
               
            // return an image of all 1's
            mask.fill(1);
            
            return mask;
        }
      
        // now the coordinates in zeroPointLists are w.r.t. thetaImg

        removeSetsWithNonCloudColors(zeroPointLists, colorImg, thetaImg,
            true, 0);
        
        if (zeroPointLists.isEmpty()) {
            
            GreyscaleImage mask = theta.createWithDimensions();
               
            // return an image of all 1's
            mask.fill(1);
            
            return mask;
        }
        
        //TODO:  assumes that largest smooth component of image is sky.
        // if sky is small and a foreground object is large and featureless
        // and not found as dark, this will fail. 
        // will adjust for that one day, possibly with color validation
        reduceToLargest(zeroPointLists);
        
        if (zeroPointLists.isEmpty()) {
            
            GreyscaleImage mask = theta.createWithDimensions();
               
            // return an image of all 1's
            mask.fill(1);
            
            return mask;
        }
        
        double[] avgYRGB = calculateYRGB(zeroPointLists.get(0), colorImg,
             thetaImg.getXRelativeOffset(), thetaImg.getYRelativeOffset(), 
             true, 0);
        
        int nBeforeHighContrastRemoval = count(zeroPointLists);
        
        Set<PairInt> highContrastRemoved = new HashSet<PairInt>();
        removeHighContrastPoints(zeroPointLists, colorImg, 
            thetaImg, avgYRGB[0], true, 0, highContrastRemoved);
        
        if (zeroPointLists.isEmpty()) {
            GreyscaleImage mask = theta.createWithDimensions();
            // return an image of all 1's
            mask.fill(1);
            return mask;
        }
        
        int nAfterHighContrastRemoval = 0;
        for (int ii = 0; ii < zeroPointLists.size(); ii++) {
            nAfterHighContrastRemoval += zeroPointLists.get(ii).getN();
        }
        
        removeReflectedSun(zeroPointLists, colorImg, thetaImg);
        
        if (zeroPointLists.isEmpty()) {
            GreyscaleImage mask = theta.createWithDimensions();
            mask.fill(1);
            return mask;
        }
        
        Set<PairInt> points = combine(zeroPointLists);
       
        int valueToSubtract = extractSkyFromGradientXY(gXYImg, points,
            highContrastRemoved);
        
        if (binFactor > 1) {
            
            thetaImg = theta;
            colorImg = originalColorImage;
            gXYImg = gradientXY;
            
            points = unbinZeroPointLists(points, binFactor);
        }
        
        PerimeterFinder perimeterFinder = new PerimeterFinder();
        int[] skyRowMinMax = new int[2];
        Map<Integer, PairInt> skyRowColRange = perimeterFinder.find(points, 
            skyRowMinMax);
        //TODO:  revisit this... shouldn't be necessary when grow to sky is impl 
        fillInRightBoundarySkyPoints(points, binFactor, 
            skyRowColRange, skyRowMinMax, colorImg,
            thetaImg.getXRelativeOffset(), thetaImg.getYRelativeOffset());
        
        GreyscaleImage mask = gradientXY.createWithDimensions();
        mask.fill(1);
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();            
            mask.setValue(x, y, 0);
        }
        
        int nSkyPointsBeforeFindClouds = points.size();
        
        Map<Integer, PixelColors> pixelColorsMap = new 
            HashMap<Integer, PixelColors>();
        
        Map<PairInt, Set<PixelColors> > skyColorsMap = new HashMap<PairInt, 
            Set<PixelColors> >();
        
        populatePixelColorMaps(points, colorImg, mask, pixelColorsMap, 
            skyColorsMap);
        
        
        findClouds(points, colorImg, mask, pixelColorsMap, skyColorsMap);
        
        
        Set<PairInt> sunPoints = findSunConnectedToSkyPoints(points, 
            colorImg, thetaImg.getXRelativeOffset(), 
            thetaImg.getYRelativeOffset());
        
        Set<PairInt> rainbowPoints = findRainbowConnectedToSkyPoints(
            sunPoints, points, colorImg, thetaImg.getXRelativeOffset(), 
            thetaImg.getYRelativeOffset());
        
        // find embedded non-sky pixels:
        PerimeterFinder finder = new PerimeterFinder();
        int[] rowMinMax = new int[2];
        Map<Integer, PairInt> rowColRange = finder.find(points, rowMinMax);

        DFSContiguousValueFinder contiguousNonZeroFinder = new DFSContiguousValueFinder(mask);
        contiguousNonZeroFinder.setMinimumNumberInCluster(1);
        contiguousNonZeroFinder.findEmbeddedGroupsNotThisValue(0,
            rowColRange, rowMinMax);

        Set<PairInt> embeddedPoints = new HashSet<PairInt>();
        for (int i = 0; i < contiguousNonZeroFinder.getNumberOfGroups(); i++) {
            contiguousNonZeroFinder.getXY(i, embeddedPoints);
        }
        
        if (embeddedPoints.size() > 0) {

            HistogramHolder[] brightnessHistogram = new HistogramHolder[1];

            // brightest sky is in bin [2], and dimmest is in [0]
            GroupPixelColors[] skyPartitionedByBrightness = 
                partitionInto3ByBrightness(points, originalColorImage, 
                thetaImg.getXRelativeOffset(), 
                thetaImg.getYRelativeOffset(), brightnessHistogram);

            findEmbeddedClouds(embeddedPoints, points, highContrastRemoved,
                colorImg, mask,
                pixelColorsMap, skyColorsMap,
                brightnessHistogram, skyPartitionedByBrightness);

            if ((nBeforeHighContrastRemoval - nAfterHighContrastRemoval) < 
                (int)(((float)nBeforeHighContrastRemoval)*0.1f)) {
            
                // also, do not perform this if we can tell that the sky is mostly
                // found except the nBinFactor pixels near the border.
                // (wanting to avoid overrunning low contrast skylines such as
                // snow covered peaks under a hazy sky, etc)
                if (embeddedPoints.size() > 0) {
               //   findSeparatedClouds(sunPoints, points, highContrastRemoved,
               //       colorImg, mask,
               //       pixelColorsMap, skyColorsMap,
               //       brightnessHistogram, skyPartitionedByBrightness);
                }
            }
        }
         
        //TODO: refine this number comparison
        float f = (float)points.size()/(float)nSkyPointsBeforeFindClouds;
        if (f > 2.0) {
        //    findMoreClouds(points, colorImg, mask);
        }
        
        //TODO: one more round of remove unconnected? or rather, only
        // keep the largest group of connected zeros?
        // see the stone henge test
        
        for (PairInt p : sunPoints) {
            int x = p.getX();
            int y = p.getY();            
            mask.setValue(x, y, 0);
        }
        
        points.addAll(sunPoints);
        
        /*
        grow the pixels towards the sky boundary.
        Looks like contrast and color should be used here as a safer boundary
        that works for low contrast regions too (else, could have just used
        the difference of gaussians alone).
        */
        //growToSkyline(points, skyRowColRange, skyRowMinMax, colorImg, 
        //    thetaImg.getXRelativeOffset(), thetaImg.getYRelativeOffset());
            
        /*
        avgYRGB = calculateYRGB(zeroPointLists.get(0), colorImg,
             theta.getXRelativeOffset(), theta.getYRelativeOffset(), 
             makeCorrectionsAlongX, convDispl);
        */
        
        //Set<PairInt> added = growPointsToSkyline(points, colorImg, 
        //    mask, avgYRGB);

        /*
        // TODO: if sky is red, use blue segmentation.
        GreyscaleImage tmpRed = colorImg.copyRedToImage();
        GreyscaleImage tmpBlue = colorImg.copyBlueToImage();
        MatrixUtil.multiply(tmpRed.getValues(), mask.getValues());
        MatrixUtil.multiply(tmpBlue.getValues(), mask.getValues());
        applyImageSegmentation(tmpRed, 3);        
        applyImageSegmentation(tmpBlue, 3);
        debugPlot(new HashSet<PairInt>(), tmpRed.copyImageToGreen(), 0, 0, 
            "red_seg_k3");
        debugPlot(new HashSet<PairInt>(), tmpBlue.copyImageToGreen(), 0, 0, 
            "blue_seg_k3");
       */
        
        for (PairInt p : rainbowPoints) {
            int x = p.getX();
            int y = p.getY();            
            mask.setValue(x, y, 0);
        }
        points.addAll(rainbowPoints);
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        double[] xycen = curveHelper.calculateXYCentroids(points);
 
        outputSkyCentroid.add((int)Math.round(xycen[0]), (int)Math.round(xycen[1]));

        if (mask != null) {
        //    removeSpurs(mask);
        }         
   
        return mask;
    }
    
    public Set<PairInt> combine(List<PairIntArray> points) {
        Set<PairInt> set = new HashSet<PairInt>();
        for (PairIntArray p : points) {
            for (int i = 0; i < p.getN(); i++) {
                int x = p.getX(i);
                int y = p.getY(i);
                PairInt pi = new PairInt(x, y);
                set.add(pi);
            }
        }
        return set;
    }
    
    public List<PairIntArray> getLargestSortedContiguousZeros(GreyscaleImage theta) {
                 
        return getSortedContiguousValues(theta, 0, false, true);
    }
    
    public List<PairIntArray> getSortedContiguousZeros(GreyscaleImage theta) {
                 
        return getSortedContiguousValues(theta, 0, false, false);
    }
    
    public List<PairIntArray> getLargestSortedContiguousNonZeros(GreyscaleImage theta) {
                 
        return getSortedContiguousValues(theta, 0, true, true);
    }
    
    private List<PairIntArray> getSortedContiguousValues(GreyscaleImage theta,
        int value, boolean excludeValue, boolean limitToLargest) {
        
        DFSContiguousValueFinder zerosFinder = new DFSContiguousValueFinder(theta);
        
        if (excludeValue) {
            zerosFinder.findGroupsNotThisValue(value);
        } else {
            zerosFinder.findGroups(value);
        }
        
        int nGroups = zerosFinder.getNumberOfGroups();
        
        if (nGroups == 0) {
            return new ArrayList<PairIntArray>();
        }
        // ====== find the group(s) with the largest number of zero pixels =====
        
        int[] groupIndexes = new int[nGroups];
        int[] groupN = new int[nGroups];
        for (int gId = 0; gId < nGroups; gId++) {
            int n = zerosFinder.getNumberofGroupMembers(gId);
            groupIndexes[gId] = gId;
            groupN[gId] = n;
        }
        
        MultiArrayMergeSort.sortByDecr(groupN, groupIndexes);
        
        List<Integer> groupIds = new ArrayList<Integer>();
        groupIds.add(Integer.valueOf(groupIndexes[0]));
        
        if (nGroups > 1) {
                        
            float n0 = (float)groupN[0];
            
            for (int i = 1; i < groupN.length; i++) {
                
                if (limitToLargest) {
                    float number = groupN[i];

                    float frac = number/n0;
    System.out.println(number);    
                    //TODO: this should be adjusted by some metric.
                    //      a histogram?
                    // since most images should have been binned to <= 300 x 300 pix,
                    // making an assumption about a group >= 100 pixels 
                    if ((1 - frac) < 0.4) {
                    //if (number > 100) {
                        groupIds.add(Integer.valueOf(groupIndexes[i]));
                    } else {
                        break;
                    }
                } else {
                    groupIds.add(Integer.valueOf(groupIndexes[i]));
                }
            }
        }

        List<PairIntArray> list = new ArrayList<PairIntArray>();
        
        for (Integer gIndex : groupIds) {
            
            int gIdx = gIndex.intValue();
            
            PairIntArray points = zerosFinder.getXY(gIdx);
            
            list.add(points);
        }
        
        return list;
    }
    
    public void multiply(GreyscaleImage input, float m) {
        
        for (int col = 0; col < input.getWidth(); col++) {
            
            for (int row = 0; row < input.getHeight(); row++) {
                
                int v = input.getValue(col, row);
                
                int f = (int)(m * v);
                
                input.setValue(col, row, f);
            }
        }
    }
    
    public void subtractMinimum(GreyscaleImage input) {
        
        int min = MiscMath.findMin(input.getValues());
        
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
     * No corrections are made for integer overflow.
     * @param input1
     * @param input2
     * @return 
     */
    public GreyscaleImage multiply(GreyscaleImage input1, GreyscaleImage input2)  {
        
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
                
                int v = input1.getValue(col, row) * input2.getValue(col, row);
                                
                output.setValue(col, row, v);
            }
        }
        
        return output;
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
    
    public void blur(GreyscaleImage input, float sigma) {
                        
        float[] kernel = Gaussian1D.getKernel(sigma);
        
        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();
        
        GreyscaleImage output = input.createWithDimensions();
        
        for (int i = 0; i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {
                double conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, true);                
                int g = (int) conv;
                output.setValue(i, j, g);
            }
        }
        
        input.resetTo(output);
        
        for (int i = 0; i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {
                double conv = kernel1DHelper.convolvePointWithKernel(
                    input, i, j, kernel, false);                
                int g = (int) conv;
                output.setValue(i, j, g);
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
                    if ((i==250) && (j >= 50) && (j <= 150)) {
                        log.info(Float.toString(r));
                    }
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
            
            log.info("numRemoved=" + numRemoved + " nIter=" + nIter);
            
            nIter++;
        }
        
        if (nIter  > 30) {
            try {
                multiply(input, 255);
                ImageDisplayer.displayImage("segmented for sky", input);
                int z = 1;
            } catch (IOException ex) {
                Logger.getLogger(ImageProcessor.class.getName()).log(Level.SEVERE, null, ex);
            }
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

    public GreyscaleImage binImageToKeepZeros(GreyscaleImage img, 
        int binFactor) {
        
        if (img == null) {
            throw new IllegalArgumentException("img cannot be null");
        }
        
        int w0 = img.getWidth();
        int h0 = img.getHeight();
        
        int w1 = w0/binFactor;
        int h1 = h0/binFactor;
      
        GreyscaleImage out = new GreyscaleImage(w1, h1);
        out.setXRelativeOffset(Math.round(img.getXRelativeOffset()/2.f));
        out.setYRelativeOffset(Math.round(img.getYRelativeOffset()/2.f));
        
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
    
    public GreyscaleImage binImage(GreyscaleImage img, 
        int binFactor) {
        
        if (img == null) {
            throw new IllegalArgumentException("img cannot be null");
        }
        
        int w0 = img.getWidth();
        int h0 = img.getHeight();
        
        int w1 = w0/binFactor;
        int h1 = h0/binFactor;
      
        GreyscaleImage out = new GreyscaleImage(w1, h1);
        out.setXRelativeOffset(Math.round(img.getXRelativeOffset()/2.f));
        out.setYRelativeOffset(Math.round(img.getYRelativeOffset()/2.f));
        
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
      
        Image out = new Image(w1, h1);
        
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
                        
                        int rgb = img.getRGB(ii, jj);
                        
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
                
                out.setRGB(i, j, (int)rSum, (int)gSum, (int)bSum);
            }
        }
        
        return out;
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
                
                for (int ii = (i*binFactor); ii < ((i + 1)*binFactor); ii++) {
                    for (int jj = (j*binFactor); jj < ((j + 1)*binFactor); jj++) {
                        out.setValue(ii, jj, v);
                    }
                }
            }
        }
        
        if ((originalTheta.getWidth() & 1) == 1) {
            // copy next to last column into last column
            int i = originalTheta.getWidth() - 2;
            for (int j = 0; j < h1; j++) {
                int v = out.getValue(i, j);
                out.setValue(i + 1, j, v);
            }
        }
        if ((originalTheta.getHeight() & 1) == 1) {
            // copy next to last row into last row
            int j = originalTheta.getHeight() - 2;
            for (int i = 0; i < w1; i++) {
                int v = out.getValue(i, j);
                out.setValue(i, j + 1, v);
            }
        }
        
        // TODO: consider correction for oversampling at location of skyline
        // using originalTheta
        
        return out;
    }

    private List<PairIntArray> unbinZeroPointLists(List<PairIntArray> zeroPointLists, 
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

    private Set<PairInt> unbinZeroPointLists(Set<PairInt> zeroPoints, 
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
    
    private void removeSetsWithNonCloudColors(List<PairIntArray> 
        zeroPointLists, Image originalColorImage, GreyscaleImage theta,
        boolean addAlongX, int addAmount) {
        
        int colorLimit = 100;
        
        List<Integer> remove = new ArrayList<Integer>();
        
        int xOffset = theta.getXRelativeOffset();
        int yOffset = theta.getYRelativeOffset();
        
        for (int gId = 0; gId < zeroPointLists.size(); gId++) {
            
            PairIntArray points  = zeroPointLists.get(gId);
            
            int nBelowLimit = 0;
            int nBelowLimit2 = 0;
            int nBelowLimit3 = 0;
            
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
                
                int r = originalColorImage.getR(ox, oy);
                int g = originalColorImage.getG(ox, oy);
                int b = originalColorImage.getB(ox, oy);
                
                if ((r < colorLimit) && (b < colorLimit) && (g < colorLimit)) {
                    nBelowLimit++;
                }
                // tracking color of sand
                if ((Math.abs((r/255.) - 0.35) < .10) 
                    && (Math.abs((g/255.) - 0.35) < .10) 
                    && (Math.abs((b/255.) - 0.35) < .10)
                    && (b < g)) {
                    nBelowLimit2++;
                }
                // TODO: there has to be a better way to filter green and cyan
                // green colors
                if ((Math.abs((r/255.) - 0.55) < .05) 
                    && (Math.abs((g/255.) - 0.62) < .05) 
                    && (Math.abs((b/255.) - 0.25) < .05)
                    && (b < g)) {
                    nBelowLimit3++;
                } else if (((g - r) >= 20) && ((g - b) >= 20)) {
                    nBelowLimit3++;
                } else if ((Math.abs((r/255.) - 0.25) < .03) 
                    && (Math.abs((g/255.) - 0.38) < .05) 
                    && (Math.abs((b/255.) - 0.40) < .05)
                    && (g > 80)) {
                    nBelowLimit3++;
                } else if ((Math.abs((r/255.) - 0.41) < .03) 
                    && (Math.abs((g/255.) - 0.48) < .03) 
                    && (Math.abs((b/255.) - 0.49) < .03)
                    && (g > 80)) {
                    nBelowLimit3++;
                }
            }
            
            log.fine(gId + ") nBelowLimit=" + nBelowLimit
                + " (" + ((double)nBelowLimit/(double)points.getN()) + ")");
         
            if (((double)nBelowLimit/(double)points.getN()) > 0.5) {
                
                remove.add(Integer.valueOf(gId));
                
            } else if (((double)nBelowLimit2/(double)points.getN()) > 0.5) {
                
                remove.add(Integer.valueOf(gId));
                
            } else if (((double)nBelowLimit3/(double)points.getN()) > 0.5) {
                
                remove.add(Integer.valueOf(gId));
            }
        }
        
        if (!remove.isEmpty()) {
            for (int i = (remove.size() - 1); i > -1; i--) {
                zeroPointLists.remove(remove.get(i).intValue());
            }
        }
        
    }

    private void reduceToLargest(List<PairIntArray> zeroPointLists) {
        
        int rmIdx = -1;
        
        if (zeroPointLists.size() > 1) {
                        
            float n0 = (float)zeroPointLists.get(0).getN();
            
            for (int i = 1; i < zeroPointLists.size(); i++) {
                
                float number = zeroPointLists.get(i).getN();

                float frac = number/n0;
    System.out.println(number + " n0=" + n0);    
                //TODO: this should be adjusted by some metric.
                //      a histogram?
                // since most images should have been binned to <= 300 x 300 pix,
                // making an assumption about a group >= 100 pixels 
                if (frac < 0.1) {
                    rmIdx = i;
                    break;
                }
            }
        }
        
        if (rmIdx > -1) {
            List<PairIntArray> out = new ArrayList<PairIntArray>();
            for (int i = 0; i < rmIdx; i++) {
                out.add(zeroPointLists.get(i));
            }
            zeroPointLists.clear();
            zeroPointLists.addAll(out);
        }
    }
    
    /**
     * remove high contrast points from the sky points.  this helps to remove
     * points that are present due to "blind spots" in gradientXY on the scale
     * of the combined convolution of gaussians that created the gradientXY.
     * For example, repetitive structure like skyscraper windows are objects
     * in the color image which may be missing an outline in the gradientXY.
     * These features have higher contrast in the color image than the normal
     * sky because they are objects, not sky, so this method tries to find
     * those and remove them from the sky points.  Note that the method prefers
     * to err on the side of over subtracting because later steps can find 
     * any removed connected sky as long as the majority of sky points remain.
     * 
     * @param zeroPointLists
     * @param originalColorImage
     * @param theta
     * @param avgY
     * @param addAlongX
     * @param addAmount 
     * @param outputRemovedPoints is populated with removed points if 
     * the object is not null
     */
    private void removeHighContrastPoints(List<PairIntArray> 
        zeroPointLists, Image originalColorImage, GreyscaleImage theta,
        double avgY, boolean addAlongX, int addAmount, 
        Set<PairInt> outputRemovedPoints) {
        
        int xOffset = theta.getXRelativeOffset();
        int yOffset = theta.getYRelativeOffset();
        
        // remove points that have contrast larger than tail of histogram
        HistogramHolder h = createContrastHistogram(avgY, zeroPointLists, 
            originalColorImage, xOffset, yOffset, addAlongX, 
            addAmount);
        
        if (h == null) {
            return;
        }

try {
    h.plotHistogram("contrast", 1);
    // printed as bin/classes/points_and_polygon1.html
} catch (IOException e) {
    log.severe(e.getMessage());
}
        List<Integer> strongPeaks = MiscMath.findStrongPeakIndexes(h, 0.1f);
         
        if (strongPeaks == null || strongPeaks.isEmpty()) {
            return;
        }
       
        int lastPeakIdx = strongPeaks.get(strongPeaks.size() - 1).intValue();
        
        //TODO: this is sensitive to the histogram formation so needs a wide
        // variety of data for testing to make sure it finds the right characteristic
        
        int yPeakIdx = lastPeakIdx;
        int tailXIdx = h.getXHist().length - 1;
        if (tailXIdx > yPeakIdx) {
            float yPeak =  h.getYHist()[yPeakIdx];
            float crit = 0.03f;
            float dy = Float.MIN_VALUE;
            for (int i = (yPeakIdx + 1); i < h.getYHist().length; i++) {
                
                float f = (float)h.getYHist()[i]/yPeak;
                dy = Math.abs(h.getYHist()[i] - h.getYHist()[i - 1]);
                
                System.out.println("x=" + h.getXHist()[i] + " f=" + f + " dy=" + dy);
                
                if (f < crit) {
                    tailXIdx = i;
                    break;
                }
            }
        }
        
        double[][] m = new double[3][];
        m[0] = new double[]{0.256, 0.504, 0.098};
        m[1] = new double[]{-0.148, -0.291, 0.439};
        m[2] = new double[]{0.439, -0.368, -0.072};
        
        //remove points w/ contrast higher than the tail of the histogram
        double critContrast = h.getXHist()[tailXIdx];
             
        for (int gId = 0; gId < zeroPointLists.size(); gId++) {
                        
            PairIntArray points  = zeroPointLists.get(gId);
                        
            Set<PairInt> pointsSet = new HashSet<PairInt>();
            for (int i = 0; i < points.getN(); i++) {
                int x = points.getX(i);
                int y = points.getY(i);
                PairInt pi = new PairInt(x, y);
                pointsSet.add(pi);
            }

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
                yuv = MatrixUtil.add(yuv, new double[]{16, 128, 128});

                float contrast = (float)((avgY - yuv[0]) / yuv[0]);
                
                if (contrast > critContrast) {
                    
                    PairInt pi0 = new PairInt(x, y);
                    
                    pointsSet.remove(pi0);
                    
                    if (outputRemovedPoints != null) {
                        outputRemovedPoints.add(pi0);
                    }

                    for (int xx = (x - 1); xx <= (x + 1); xx++) {
                        if ((xx < 0) || (xx > (theta.getWidth() - 1))) {
                            continue;
                        }
                        for (int yy = (y - 1); yy <= (y + 1); yy++) {
                            if ((yy < 0) || (yy > (theta.getHeight() - 1))) {
                                continue;
                            }
                            if ((xx == x) && (yy == y)) {
                                continue;
                            }
                            
                            PairInt pi1 = new PairInt(xx, yy);
                            
                            if (pointsSet.contains(pi1)) {
                                pointsSet.remove(pi1);
                                if (outputRemovedPoints != null) {
                                    outputRemovedPoints.add(pi1);
                                }
                            }
                        }
                    }
                }
            }

            if (pointsSet.size() != points.getN()) {
                PairIntArray points2 = new PairIntArray();
                for (PairInt pi : pointsSet) {
                    points2.add(pi.getX(), pi.getY());
                }
                points.swapContents(points2);
            }
        }
        
        // remove empty sets
        for (int i = (zeroPointLists.size() - 1); i > -1; i--) {
            PairIntArray point = zeroPointLists.get(i);
            if (point.getN() == 0) {
                zeroPointLists.remove(i);
            }
        }
        
debugPlot(outputRemovedPoints, originalColorImage, xOffset, yOffset, "filtered_out_high_contrast");

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
    
    /**
     * if the range of contrast is large, return a contrast histogram,
     * else return null.  
     * TODO: refactor to move the logic to return null to the invoker. 
     * For now, the only use of this method is simpler if it does not
     * return a histogram when it won't be needed.  definitely should
     * be refactored...
     * 
     * @param avgY
     * @param zeroPointLists
     * @param originalColorImage
     * @param xOffset
     * @param yOffset
     * @param addAlongX
     * @param addAmount
     * @return 
     */
    private HistogramHolder createContrastHistogram(double avgY,
        List<PairIntArray> zeroPointLists, Image originalColorImage,
        int xOffset, int yOffset, boolean addAlongX, int addAmount) {
        
        if (zeroPointLists.isEmpty()) {
            return null;
        }
        
        double[][] m = new double[3][];
        m[0] = new double[]{0.256, 0.504, 0.098};
        m[1] = new double[]{-0.148, -0.291, 0.439};
        m[2] = new double[]{0.439, -0.368, -0.072};
                
        int nPoints = 0;
        for (int gId = 0; gId < zeroPointLists.size(); gId++) {
            nPoints += zeroPointLists.get(gId).getN();
        }
        
        float[] yValues = new float[nPoints];
        
        int count = 0;
        for (int gId = 0; gId < zeroPointLists.size(); gId++) {
            
            PairIntArray points  = zeroPointLists.get(gId);
                        
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
                
                int r = originalColorImage.getR(ox, oy);
                int g = originalColorImage.getG(ox, oy);
                int b = originalColorImage.getB(ox, oy);
                
                double[] rgb = new double[]{r, g, b};
                        
                double[] yuv = MatrixUtil.multiply(m, rgb);
                yuv = MatrixUtil.add(yuv, new double[]{16, 128, 128});

                float contrastValue = (float)((avgY - yuv[0]) / yuv[0]);
                
                yValues[count] = contrastValue;
                
                count++;
            }
        }
        
        float[] yErr = Errors.populateYErrorsBySqrt(yValues);
        HistogramHolder h = Histogram.createSimpleHistogram(yValues, 
            yErr);

        if (h != null) {
            
            int lastZeroIdx = MiscMath.findLastZeroIndex(h);
            
            int nBins = h.getXHist().length;
            
            int nLastZeros = nBins - lastZeroIdx;
            
            if ((nLastZeros > 4) && (lastZeroIdx > -1)) {
                
                float halfBinWidth = (h.getXHist()[1] - h.getXHist()[1]) / 2.f;

                nBins = lastZeroIdx - 1;
                float xMin = h.getXHist()[0] - halfBinWidth;
                float xMax = h.getXHist()[lastZeroIdx - 1] + halfBinWidth;

                float contrastRange = h.getXHist()[lastZeroIdx - 1]
                    - h.getXHist()[0];

                if (contrastRange > 1.5) {

                    h = Histogram.calculateSturgesHistogramRemoveZeroTail(yValues, yErr);
                    
                    if ((h != null) && (h.getXHist().length == 1)) {
                        h = Histogram.calculateSturgesHistogram(xMin, xMax, yValues, yErr);
                    }
                    
                    return h;

                }
            }
        }
        
        return null;
    }
    
    private void transformPointsToOriginalReferenceFrame(Set<PairInt> points,
        GreyscaleImage theta, boolean makeCorrectionsAlongX, int addAmount) {
        
         // transform points to original color image frame
        int totalXOffset = theta.getXRelativeOffset();
        int totalYOffset = theta.getYRelativeOffset();

        if (makeCorrectionsAlongX) {
            totalXOffset += addAmount;
        } else {
            totalYOffset += addAmount;
        }
        
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            x += totalXOffset;
            y += totalYOffset;
            p.setX(x);
            p.setY(y);
        }
    }

    /**
     * NOT READY FOR USE YET.  THIS should be run on the binned images.
     * The images which are more than half sky should especially be
     * processed at binned size.
     * 
     * @param points
     * @param originalColorImage
     * @param theta
     * @param avgYRGB
     * @return 
     */
    private Set<PairInt> growPointsToSkyline(Set<PairInt> points, 
        Image originalColorImage, GreyscaleImage mask, double[] avgYRGB) {
        
        // mask needs to be in the frame of the theta image
        
        // transform points to original color image frame
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();

        //transformPointsToOriginalReferenceFrame(points, theta, 
        //    makeCorrectionsAlongX, addAmount);
        
        double yAvg = avgYRGB[0];
        double rColor = avgYRGB[1];
        double gColor = avgYRGB[2];
        double bColor = avgYRGB[3];
        
        // (r-b)/255 < 0.2 
        boolean useBlue = (((rColor - bColor)/255.) < 0.2);
        if (useBlue && (((rColor/bColor) > 1.0)) && (bColor < 128)) {
            useBlue = false;
        }
        
        // determine contrast and blue or red for each point in points
        Map<PairInt, PairFloat> contrastAndColorMap = calculateContrastAndBOrR(
            points, useBlue, originalColorImage, avgYRGB,
            xOffset, yOffset);

        // determine differences in contrast and blue or red for each point 
        // from their neighbors (if they are in the map) and return
        // the avg and st.dev. of 
        // {avg dContrast, stdDev dContrast, avg dBlueOrRed, stDev dBlueOrRed}
        float[] diffsAvgAndStDev = calculateAvgAndStDevOfDiffs(
            contrastAndColorMap);

        log.fine("diffsAvgAndStDev=" + Arrays.toString(diffsAvgAndStDev));
        
        float diffContrastAvg = diffsAvgAndStDev[0];
        float diffContrastStDev = diffsAvgAndStDev[1];
        float diffBlueOrRedAvg = diffsAvgAndStDev[2];
        float diffBlueOrRedStDev = diffsAvgAndStDev[3];

        // TODO: improve this setting.  configurable by stats?
        /*
        could start at reasonable value and if too high as seen by the entire
        image found as sky, reduce the values and iterate...
        */
        /*
        for halfdome test image:
            float contrastFactor = 3.5f;
            float colorFactor = -2.5f;
        for NM test image:
            
        */
        float contrastFactor = 3.5f;
        float colorFactor = -2.368f; 

        Set<PairInt> added = new HashSet<PairInt>();
        
        /*
        given map of points and their contrasts and colors and the avg changes
        in those and the std dev of that,
        use dfs to look for neighbors of the points that fall within the critical
        limits.  when the point is within limit, it gets marked as a '0'
        in the output image (which is otherwise '1')
        */
        
        int width = mask.getWidth();
        int height = mask.getHeight();
        
        double[][] m = new double[3][];
        m[0] = new double[]{0.256, 0.504, 0.098};
        m[1] = new double[]{-0.148, -0.291, 0.439};
        m[2] = new double[]{0.439, -0.368, -0.072};
        
        java.util.Stack<PairInt> stack = new java.util.Stack<PairInt>();
        
        //O(N_sky)
        for (PairInt p : points) {
            stack.add(p);
        }
       
        // null = unvisited, presence = visited
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(stack.peek());
        
        while (!stack.isEmpty()) {

            PairInt uPoint = stack.pop();
            
            int uX = uPoint.getX();
            int uY = uPoint.getY();

            //(1 + frac)*O(N) where frac is the fraction added back to stack
            
            for (int vX = (uX - 1); vX <= (uX + 1); vX++) {
                
                if ((vX < 0) || (vX > (width - 1))) {
                    continue;
                }
                
                for (int vY = (uY - 1); vY <= (uY + 1); vY++) {
                    if ((vY < 0) || (vY > (height - 1))) {
                        continue;
                    }
                    
                    PairInt vPoint = new PairInt(vX, vY);
                    
                    if (vPoint.equals(uPoint)) {
                        continue;
                    }
                                    
                    if (visited.contains(vPoint)) {
                        continue;
                    }
                    
                    if (contrastAndColorMap.containsKey(vPoint)) {
                        continue;
                    }
                    
                    int ox = vX + xOffset;
                    int oy = vY + yOffset;
                    
                    visited.add(vPoint);
                    
                    if ((ox < 0) || (ox > (originalColorImage.getWidth() - 1))) {
                        continue;
                    }
                    if ((oy < 0) || (oy > (originalColorImage.getHeight() - 1))) {
                        continue;
                    }
                    
                    int rV = originalColorImage.getR(ox, oy);
                    int gV = originalColorImage.getG(ox, oy);
                    int bV = originalColorImage.getB(ox, oy);
                    
                    double[] rgbV = new double[]{rV, gV, bV};
                    double[] yuv = MatrixUtil.multiply(m, rgbV);
                    yuv = MatrixUtil.add(yuv, new double[]{16, 128, 128});

                    float contrastV = (float)((yAvg - yuv[0]) / yuv[0]);
                    
                    // if delta constrast and delta blue or red are within
                    // limits, add to stack and set in mask
                    
                    PairFloat uContrastAndColor = contrastAndColorMap.get(uPoint);
                    
                    float vMinusUContrast = contrastV - uContrastAndColor.getX();
                    
                    boolean withinLimits = true;
                    
// contrast:  ucontrast >= (vContrast + 16*diffContrastStDev), v:42,  s:26, sh:246, norw:22, hd:5,   nm:64,  az: 1
// color:                               -18                  , v:-48, s:-12, sh:-80, norw:-15, hd:-16, nm:-79, az:+15 ((r-b)=132)
// 
if ((ox > 560) && (ox <= 570) && (oy > 425) && (oy <= 431) && (contrastV > uContrastAndColor.getX())) {
    log.info("\ny = " + oy + " diffsAvgAndStDev=" + Arrays.toString(diffsAvgAndStDev));
    int color = (useBlue ? bV : rV);
    
    float vColor = useBlue ? bV : rV;

    float VMinusUColor = vColor - uContrastAndColor.getY();

    String str = String.format("useBlue=%b u(c, c)=(%f,%f)  v(c,c)=(%f,%d)  dContrast=%f  dColor=%f  (vMinusUContrast/diffContrastAvg)=%f  (VMinusUColor/diffBlueOrRedAvg)=%f",
        useBlue, uContrastAndColor.getX(), uContrastAndColor.getY(),
        contrastV, color, (contrastV - uContrastAndColor.getX()),
        (color - uContrastAndColor.getY()), (vMinusUContrast/diffContrastAvg),
        (VMinusUColor/diffBlueOrRedAvg));
    log.info(str);
    int z = 1;
}

                    if ((vMinusUContrast/diffContrastAvg) >= contrastFactor) {
                        
                        // see if color has decreased
                        float vColor = useBlue ? bV : rV;
                        
                        float VMinusUColor = vColor - uContrastAndColor.getY();
                        
                        //TODO: consider either only -10 or abs
                        if ((VMinusUColor/diffBlueOrRedAvg) < colorFactor) {
                        //if (Math.abs(VMinusUColor/diffBlueOrRedAvg) > Math.abs(colorFactor)) {
                            withinLimits = false;
                        }
                    }
                    
                    if (withinLimits) {
                        
                        stack.add(vPoint);
                        
                        added.add(vPoint);
                        
                        float vColor = useBlue ? bV : rV;
                        
                        PairFloat vCC = new PairFloat(contrastV, vColor);
                        
                        contrastAndColorMap.put(vPoint, vCC);
                    }
                }
            }
        }
       
        points.addAll(added);
        
        for (PairInt p : added) {
            int x = p.getX();
            int y = p.getY();            
            mask.setValue(x, y, 0);
        }
        
        return added;
    }

    private Map<PairInt, PairFloat> calculateContrastAndBOrR(
        Set<PairInt> points, boolean useBlue, Image originalColorImage, 
        double[] avgYRGB, int totalXOffset, int totalYOffset) {
        
        double[][] m = new double[3][];
        m[0] = new double[]{0.256, 0.504, 0.098};
        m[1] = new double[]{-0.148, -0.291, 0.439};
        m[2] = new double[]{0.439, -0.368, -0.072};
        
        double yColor = avgYRGB[0];
        
        Map<PairInt, PairFloat> map = new HashMap<PairInt, PairFloat>();
        
        for (PairInt p : points) {
            
            int x = p.getX();
            int y = p.getY();
            
            x += totalXOffset;
            y += totalYOffset;
            
            int r = originalColorImage.getR(x, y);
            int g = originalColorImage.getG(x, y);
            int b = originalColorImage.getB(x, y);

            double[] rgb = new double[]{r, g, b};
            double[] yuv = MatrixUtil.multiply(m, rgb);
            yuv = MatrixUtil.add(yuv, new double[]{16, 128, 128});

            float contrast= (float) ((yColor - yuv[0]) / yuv[0]);
            
            PairFloat crb = new PairFloat();
            crb.setX(contrast);
            if (useBlue) {
                crb.setY(b);
            } else {
                crb.setY(r);
            }
            
            map.put(p, crb);
        }
        
        return map;
    }

    private float[] calculateAvgAndStDevOfDiffs(Map<PairInt, PairFloat> 
        contrastAndColorMap) {
        
        // average difference from neighbors
        double avgContrast = 0;
        
        double avgColor = 0;
        
        int count = 0;
        
        Iterator<Entry<PairInt, PairFloat> > iter = 
            contrastAndColorMap.entrySet().iterator();
        
        while (iter.hasNext()) {
            
            Entry<PairInt, PairFloat> entry = iter.next();
            
            PairInt i = entry.getKey();
            int x = i.getX();
            int y = i.getY();
            
            for (int xx = (x - 1); xx <= (x + 1); xx++) {
                for (int yy = (y - 1); yy <= (y + 1); yy++) {
                    
                    PairInt j = new PairInt(xx, yy);
                    
                    if (contrastAndColorMap.containsKey(j)) {
                        
                        PairFloat iCC = entry.getValue();
                        
                        PairFloat jCC = contrastAndColorMap.get(j);
                        
                        float diffContrast = Math.abs(iCC.getX() - jCC.getX());
                        
                        float diffColor = Math.abs(iCC.getY() - jCC.getY());
                        
                        avgContrast += diffContrast;
                        
                        avgColor += diffColor;
                        
                        count++;
                    }
                }
            }            
        }
        
        avgContrast /= (double)count;
        
        avgColor /= (double)count;
        
        // standard deviation of avg difference from neighbors
        double stDevContrast = 0;
        
        double stDevColor = 0;
        
        iter = contrastAndColorMap.entrySet().iterator();
        
        while (iter.hasNext()) {
            
            Entry<PairInt, PairFloat> entry = iter.next();
            
            PairInt i = entry.getKey();
            int x = i.getX();
            int y = i.getY();
            
            for (int xx = (x - 1); xx <= (x + 1); xx++) {
                for (int yy = (y - 1); yy <= (y + 1); yy++) {
                    
                    PairInt j = new PairInt(xx, yy);
                    
                    if (contrastAndColorMap.containsKey(j)) {
                        
                        PairFloat iCC = entry.getValue();
                        
                        PairFloat jCC = contrastAndColorMap.get(j);
                        
                        float diffContrast = Math.abs(iCC.getX() - jCC.getX());
                        diffContrast -= avgContrast;
                        
                        float diffColor = Math.abs(iCC.getY() - jCC.getY());
                        diffColor -= avgColor;
                        
                        stDevContrast += (diffContrast * diffContrast);
                        
                        stDevColor += (diffColor * diffColor);
                    }
                }
            } 
        }
        
        stDevContrast = Math.sqrt(stDevContrast/((double)count - 1));
        
        stDevColor = Math.sqrt(stDevColor/((double)count - 1));
        
        return new float[]{(float)avgContrast, (float)stDevContrast, 
            (float)avgColor, (float)stDevColor};
    }

    /**
     * using adaptive "thresholding" to subtract intensity levels from
     * gradientXY, find the contiguous zero values connected to skyPoints
     * and add them to skyPoints.
     * @param gradientXY
     * @param skyPoints
     * @param excludeThesePoints
     * @return 
     */
    public int extractSkyFromGradientXY(GreyscaleImage gradientXY,
        Set<PairInt> skyPoints, Set<PairInt> excludeThesePoints) {
                
        GreyscaleImage gXY2 = gradientXY.copyImage();
                
        // x is pixelValue , y is number of pixels holding that value
        // last in array is for the smallest pixelValue
        PairIntArray gXYValues = Histogram.createADescendingSortByKeyArray(gXY2);

        float sumF = 0;
        int valueFor0Point9 = 0;
        for (int i = (gXYValues.getN() - 1); i > -1; i--) {
            float f = (float)gXYValues.getY(i)/(float)gradientXY.getNPixels();
            sumF += f;
            if (sumF > 0.9) {
                break;
            }
            valueFor0Point9 = gXYValues.getX(i);
        }
/*        
float sumFrac = 0;
StringBuilder sb = new StringBuilder("gXY:\n");
for (int i = (gXYValues.getN() - 1); i > -1; i--) {
int sumToHighValues = 0;
for (int ii = i; ii > -1; ii--) {
    sumToHighValues += gXYValues.getY(ii);
}
float frac = (float)gXYValues.getY(i)/(float)gradientXY.getNPixels();
sumFrac += frac;
sb.append(String.format(" value=%d count=%d  f=%f  sumToEnd=%d", 
gXYValues.getX(i), gXYValues.getY(i), frac, sumToHighValues));
sb.append("sumF=").append(Float.toString(sumFrac));
sb.append("\n");
}
log.info(sb.toString());
*/
        int subtract = 0;
        int lastHistIdx = gXYValues.getN();
        
        float c0 = (float)gXYValues.getY(gXYValues.getN() - 1)/(float)gradientXY.getNPixels();
        float c1 = (float)gXYValues.getY(gXYValues.getN() - 2)/(float)gradientXY.getNPixels();
        //float c01 = c1 + c0;
        float cm01 = c0 - c1;
        if ((cm01 < 0.0) && (valueFor0Point9 <= 2)) {
            subtract = 0;
            lastHistIdx = gXYValues.getN();
        } else if (cm01 < 0.3) {
            if (valueFor0Point9 <= 2) {
                subtract = 2;
                lastHistIdx = gXYValues.getN() - 2;
            } else {
                subtract = 1;
                lastHistIdx = gXYValues.getN() - 1;
            }
        } else if ((cm01 >= 0.3) && (cm01 <= 0.45) && (valueFor0Point9 <= 2)) {
            subtract = 1;
            lastHistIdx = gXYValues.getN() - 1;
        }
        
        int nIter = 0;
        
        float originalMaxValue = gXYValues.getY(gXYValues.getN() - 1);
                                
        while ((subtract < originalMaxValue) && (nIter == 0)) {
            
            if (nIter > 0) {
                gXY2 = gradientXY.copyImage();
            }
            
            lastHistIdx--;
            if (lastHistIdx < 1) {
                break;
            }
            subtract = gXYValues.getX(lastHistIdx);

            subtractWithCorrectForNegative(gXY2, subtract);
           
            // ==== find contiguous zeros =====  
            
            growZeroValuePoints(skyPoints, excludeThesePoints, gXY2);
            
            //TODO: remove this section or move it into an aspect
                
            // === count number of embedded groups of non-zeros in skyPoints ===
            
            Set<PairInt> embeddedPoints = new HashSet<PairInt>();
            
            int nCorrectedEmbeddedGroups = extractEmbeddedGroupPoints(
                skyPoints, gXY2, embeddedPoints);
            
            log.fine("nIter=" + nIter + ")" 
                + " nCorrectedEmbeddedGroups=" + nCorrectedEmbeddedGroups
                + " nEmbeddedPixels=" + embeddedPoints.size()
                + " out of " + skyPoints.size()
                + " (level=" + ((float)gXYValues.getY(lastHistIdx)/originalMaxValue) 
                + " subtract=" + subtract + " out of max=" + gXYValues.getX(0)
                + ")"
            );
                        
            nIter++;
        }
/*
try {
    Image img1 = gradientXY.copyImageToGreen();
    ImageIOHelper.addToImage(skyPoints, 0, 0, img1);
    ImageDisplayer.displayImage("sky points subtract=" + subtract, img1);
} catch (IOException ex) {
    log.severe(ex.getMessage());
}
*/         
        return subtract;
    }

    private void growZeroValuePoints(Set<PairInt> points, 
        Set<PairInt> excludeThesePoints, GreyscaleImage gradientXY) {
  
        java.util.Stack<PairInt> stack = new java.util.Stack<PairInt>();
        
        //O(N_sky)
        for (PairInt p : points) {
            if (!excludeThesePoints.contains(p)) {
                stack.add(p);
            }
        }
        
        // null = unvisited, presence = visited
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(stack.peek());
        
        int width = gradientXY.getWidth();
        int height = gradientXY.getHeight();
        
        while (!stack.isEmpty()) {

            PairInt uPoint = stack.pop();
            
            int uX = uPoint.getX();
            int uY = uPoint.getY();

            //(1 + frac)*O(N) where frac is the fraction added back to stack
            for (int vX = (uX - 1); vX <= (uX + 1); vX++) {
                
                if ((vX < 0) || (vX > (width - 1))) {
                    continue;
                }
                
                for (int vY = (uY - 1); vY <= (uY + 1); vY++) {
                    
                    if ((vY < 0) || (vY > (height - 1))) {
                        continue;
                    }
                    
                    PairInt vPoint = new PairInt(vX, vY);
                    
                    if (vPoint.equals(uPoint)) {
                        continue;
                    }
                
                    if (excludeThesePoints.contains(vPoint)) {
                        continue;
                    }
                    
                    if (visited.contains(vPoint)) {
                        continue;
                    }
                                      
                    if ((vX < 0) || (vX > (width - 1))) {
                        continue;
                    }
                    if ((vY < 0) || (vY > (height - 1))) {
                        continue;
                    }
                    
                    visited.add(vPoint);
                    
                    int v = gradientXY.getValue(vX, vY);
                                        
                    if (v == 0) {
                        
                        stack.add(vPoint);
                        
                        if (!points.contains(vPoint)) {
                            points.add(vPoint);
                        }
                    }
                }
            }
        }
    }
    
    boolean isPerimeterUnbound(Map<Integer, PairInt> gRowColRange, 
        int[] gRowMinMax, Map<Integer, PairInt> boundingRowColRange, 
        int[] boundingRowMinMax,
        int xMinImage, int xMaxImage, int yMinImage, int yMaxImage) {
        
        // check top and bottom rows of group are within bounds
        if (gRowMinMax[0] == boundingRowMinMax[0]) {
            if (gRowMinMax[0] != yMinImage) {
                return true;
            }
        } else if (gRowMinMax[0] < boundingRowMinMax[0]) {
            return true;
        }
        if (gRowMinMax[1] == boundingRowMinMax[1]) {
            if (gRowMinMax[1] != yMaxImage) {
                return true;
            }
        } else if (gRowMinMax[1] > boundingRowMinMax[1]) {
            return true;
        }
        
        for (int r = gRowMinMax[0]; r <= gRowMinMax[1]; r++) {
                    
            PairInt cRange = gRowColRange.get(Integer.valueOf(r));
            
            // see if each point in cRange is on the boundary of rowColRange
            //   or within it
           
            // check left half            
            int x = cRange.getX();
            int y = r;
            
            PairInt cRange2 = boundingRowColRange.get(Integer.valueOf(y));
            if (cRange2 == null) {
                return true;
            }
            int x2 = cRange2.getX();
            
            if (x == x2) {
                if (x > xMinImage) {
                    return true;
                }
            } else if (x < x2) {
                // x is outside of the larger region defined by gRowColRange
                return true;
            }
            
            // check right half
            x2 = cRange2.getY();
            
            if (x == x2) {
                if (x < xMaxImage) {
                    return true;
                }
            } else if (x > x2) {
                // x is outside of the larger region defined by gRowColRange
                return true;
            }
        }
        
        return false;
    }
    
    private boolean isPerimeterUnbound(Map<Integer, PairInt> gRowColRange, 
        int[] gRowMinMax, Set<PairInt> skyPoints, double[] groupXYCen,
        int imageWidth, int imageHeight) {
        
        boolean unbounded = false;
        
        for (int r = gRowMinMax[0]; r <= gRowMinMax[1]; r++) {
                    
            PairInt cRange = gRowColRange.get(Integer.valueOf(r));
            
            for (int k = 0; k < 2; k++) {
                int c;
                switch(k) {
                    case 0:
                        c = cRange.getX();
                        break;
                    default:
                        c = cRange.getY();
                        break;
                }
                
                if (c < groupXYCen[0]) {
                
                    // look for points to left
                    int xt = c - 1;
                    if (xt < 0) {
                        // bounded by edge of image
                        continue;
                    }
                
                    if (r < groupXYCen[1]) {
                
                        //look for points to left and top (=lower y)                            
                        int yt = r;
                        PairInt p = new PairInt(xt, yt);
                        if (!skyPoints.contains(p)) {
                            // not bounded on left
                            unbounded = true;
                            break;
                        }
                        //found a sky point to the left
                        yt--;
                        if (yt < 0) {
                            // bounded by edge of image
                            continue;
                        } else {
                            p = new PairInt(xt, yt);
                            if (!skyPoints.contains(p)) {
                                // not bounded on left
                                unbounded = true;
                                break;
                            }
                        }
                        
                    } else {
                        
                        //look for bounding points to left, bottom (=higher y)
                        int yt = r;
                        PairInt p = new PairInt(xt, yt);
                        if (!skyPoints.contains(p)) {
                            // not bounded on left
                            unbounded = true;
                            break;
                        }
                        yt++;
                        if (yt > (imageHeight - 1)) {
                            // bounded by edge of image
                            continue;
                        } else {
                            p = new PairInt(xt, yt);
                            if (!skyPoints.contains(p)) {
                                // not bounded on left
                                unbounded = true;
                                break;
                            }
                        }
                    }
                
                } else {

                    // look for points to the right
                    int xt = c + 1;
                    if (xt > (imageWidth - 1)) {
                        // bounded by edge of image
                        continue;
                    }
                
                    if (r < groupXYCen[1]) {

                        //look for bounding points to right, top (=lower y),

                        int yt = r;
                        PairInt p = new PairInt(xt, yt);
                        if (!skyPoints.contains(p)) {
                            // not bounded on left
                            unbounded = true;
                            break;
                        }
                        yt--;
                        if (yt < 0) {
                            // bounded by edge of image
                            continue;
                        } else {
                            p = new PairInt(xt, yt);
                            if (!skyPoints.contains(p)) {
                                // not bounded on left
                                unbounded = true;
                                break;
                            }
                        }

                    } else {
                    
                        //look for bounding points to right, bottom (=higher y)

                        int yt = r;
                        PairInt p = new PairInt(xt, yt);
                        if (!skyPoints.contains(p)) {
                            // not bounded on left
                            unbounded = true;
                            break;
                        }
                        yt++;
                        if (yt > (imageHeight - 1)) {
                            // bounded by edge of image
                            continue;
                        } else {
                            p = new PairInt(xt, yt);
                            if (!skyPoints.contains(p)) {
                                // not bounded on left
                                unbounded = true;
                                break;
                            }
                        }
                    }
                }
            }
            
            if (unbounded) {
                break;
            }
        }
        
        return unbounded;
    }

    /**
     * attempt to find within pixels connected to skyPoints, pixels that
     * look like sun pixels by color (hsb) and whose x,y distribution
     * resemble and ellipse (circle w/ possible occlusion).  
     * those sun points are then added to the skyPoints.
     * Note that if the sun is present in sky and in reflection, such as
     * water, their location in x,y must be fittable by an ellipse, else they 
     * may not be found as sun points.
     * @param skyPoints
     * @param clr
     * @param xOffset
     * @param yOffset 
     * @return the extracted sun points that are connected to skyPoints
     */
    protected Set<PairInt> findSunConnectedToSkyPoints(Set<PairInt> skyPoints, 
        Image clr, int xOffset, int yOffset) {
        
        Set<PairInt> yellowPoints = new HashSet<PairInt>();

        java.util.Stack<PairInt> yellowStack = new java.util.Stack<PairInt>();
        
        int width = clr.getWidth();
        int height = clr.getHeight();
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        for (PairInt p : skyPoints) {

            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            
            for (int k = 0; k < dxs.length; k++) {
            
                int xx = x + dxs[k];
                int yy = y + dys[k];
                
                if ((xx < 0) || (xx > (width - 1)) || (yy < 0) || 
                    (yy > (height - 1))) {
                    continue;
                }
            
                PairInt p2 = new PairInt(xx - xOffset, yy - yOffset);

                if (yellowPoints.contains(p2)) {
                    continue;
                }

                int r = clr.getR(xx, yy);
                int g = clr.getG(xx, yy);
                int b = clr.getB(xx, yy);

                // all normalized from 0 to 1
                float[] hsb = new float[3];
                Color.RGBtoHSB(r, g, b, hsb);

                float h2 = hsb[0] * 360.f;
                float s2 = hsb[1] * 100.f;

                // increasing s2 to < 60 finds larger diameter sun and scattered light
                if ((r >= 240/*251*/) && (hsb[2] >= 0.97/*0.98*/) && (s2 < 30)) {

                    yellowPoints.add(p2);

                    yellowStack.add(p2);
                }
            }
        }
        
        if (yellowStack.isEmpty()) {
            return new HashSet<PairInt>();
        }
                
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(yellowStack.peek());
       
        while (!yellowStack.isEmpty()) {

            PairInt uPoint = yellowStack.pop();
            
            int uX = uPoint.getX() + xOffset;
            int uY = uPoint.getY() + yOffset;

            //(1 + frac)*O(N) where frac is the fraction added back to stack
            
            for (int k = 0; k < dxs.length; k++) {
                    
                int dx = dxs[k];
                int dy = dys[k];
                
                int vX = uX + dx;
                int vY = uY + dy;
                
                if ((vX < 0) || (vX > (width - 1)) || (vY < 0) || 
                    (vY > (height - 1))) {
                    continue;
                }
            
                PairInt vPoint = new PairInt(vX - xOffset, vY - yOffset);

                // skypoints has already been check for same color criteria
                // so no need to check again here
                if (visited.contains(vPoint) || skyPoints.contains(vPoint)
                    || yellowPoints.contains(vPoint)) {
                    continue;
                }

                visited.add(vPoint);

                int r = clr.getR(vX, vY);
                int g = clr.getG(vX, vY);
                int b = clr.getB(vX, vY);

                // all normalized from 0 to 1
                float[] hsb = new float[3];
                Color.RGBtoHSB(r, g, b, hsb);

                float h2 = hsb[0] * 360.f;
                float s2 = hsb[1] * 100.f;

                // increasing s2 to < 60 finds larger diameter sun and scattered light
                if ((r >= 240/*251*/) && (hsb[2] >= 0.97/*0.98*/) && (s2 < 30)) {

                    yellowPoints.add(vPoint);

                    yellowStack.add(vPoint);
                }
            }
        }
        
        if (yellowPoints.size() < 6) {
            return new HashSet<PairInt>();
        }
        
        //fit ellipse to yellowPoints.  ellipse because of possible occlusion.
        EllipseHelper ellipseHelper = new EllipseHelper();
        double[] params = ellipseHelper.fitEllipseToPoints(yellowPoints);
        
        if (params == null) {
            // not close to an ellipse
            return new HashSet<PairInt>();
        }
        
        float xc = (float)params[0];
        float yc = (float)params[1];
        float a = (float)params[2];
        float b = (float)params[3];
        float alpha = (float)params[4];
                
        //double[] stats = ellipseHelper.calculateEllipseResidualStats(
        //    yellowPoints.getX(), yellowPoints.getY(), xc, yc, a, b, alpha);
        
        return yellowPoints;
    }
    
    /**
     * this assumes findSunConnectedToSkyPoints has been run to populate
     * sunPoints and that those have not been filtered for a non-circular or
     * elliptical shape.  If there are no points in sunPoints (which are bright
     * yellow), then this method will not search for a rainbow.
     * @param sunPoints
     * @param points
     * @param colorImg
     * @param xRelativeOffset
     * @param yRelativeOffset
     * @return 
     */
    private Set<PairInt> findRainbowConnectedToSkyPoints(Set<PairInt> sunPoints, 
        Set<PairInt> points, Image colorImg, int xOffset, int yOffset) {
        
        if (sunPoints.isEmpty() || (sunPoints.size() < 12)) {
            return new HashSet<PairInt>();
        }
        
        /*
        looking for connected ribbons, first connected to the yellow sun points
        which may actually be rainbow points.
        
        */
        
        Set<PairInt> rainbowPoints = new HashSet<PairInt>();
        rainbowPoints.addAll(sunPoints);
                
        int width = colorImg.getWidth();
        int height = colorImg.getHeight();
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        java.util.Stack<PairInt> rainbowStack = new java.util.Stack<PairInt>();
        
        int nLast = rainbowPoints.size();
        
        for (int nBands = 0; nBands < 10; nBands++) {
            
            if ((nBands > 0) && (nLast == rainbowPoints.size())) {
                break;
            }
            
            nLast = rainbowPoints.size();
            
            int hueLo = 45;
            int hueHi = 65;
            int saturationLo = 0;
            int saturationHi = 35;
            int brightnessLo = 40;
            int brightnessHi = 100;
            
            if (nBands == 1) {
                hueLo = 5;
                saturationHi = 45;
            } else if (nBands == 2) {
                hueLo = 65;
                hueHi = 140;
                saturationLo = 0;
                saturationHi = 15;
            } else if ((nBands == 3) || (nBands == 4)) {
                hueHi = 200;
                saturationLo = 0;
            } else if (nBands >= 5) {
                hueLo = 200;
                hueHi = 360;
                saturationHi = 10;
            }
            
            rainbowStack.addAll(rainbowPoints);
            
            Set<PairInt> visited = new HashSet<PairInt>();
            visited.add(rainbowStack.peek());
       
            while (!rainbowStack.isEmpty()) {

                PairInt uPoint = rainbowStack.pop();
            
                int uX = uPoint.getX() + xOffset;
                int uY = uPoint.getY() + yOffset;

                //(1 + frac)*O(N) where frac is the fraction added back to stack

                for (int k = 0; k < dxs.length; k++) {
                    
                    int dx = dxs[k];
                    int dy = dys[k];

                    int vX = uX + dx;
                    int vY = uY + dy;

                    if ((vX < 0) || (vX > (width - 1)) || (vY < 0) || 
                        (vY > (height - 1))) {
                        continue;
                    }
            
                    PairInt vPoint = new PairInt(vX - xOffset, vY - yOffset);

                    // skypoints has already been check for same color criteria
                    // so no need to check again here
                    if (visited.contains(vPoint) || rainbowStack.contains(vPoint)
                        ) {
                        continue;
                    }

                    visited.add(vPoint);

                    int r = colorImg.getR(vX, vY);
                    int g = colorImg.getG(vX, vY);
                    int b = colorImg.getB(vX, vY);

                    // all normalized from 0 to 1
                    float[] hsb = new float[3];
                    Color.RGBtoHSB(r, g, b, hsb);

                    float hue = hsb[0] * 360.f;
                    float saturation = hsb[1] * 100.f;
                    float brightness = hsb[2] * 100.f;

                    // increasing s2 to < 60 finds larger diameter sun and scattered light
                    if ((hue >= hueLo) && (hue <= hueHi) && 
                        (saturation >= saturationLo) && (saturation <= saturationHi)
                        && (brightness >= brightnessLo) && (brightness <= brightnessHi)
                        ) {

                        if (nBands < 2) {
                            rainbowStack.add(vPoint);
                        }

                        rainbowPoints.add(vPoint);
                    }
                }
            }
        }
        
        Set<PairInt> tmp = new HashSet<PairInt>(rainbowPoints);
        tmp.removeAll(sunPoints);
        
        if (tmp.isEmpty()) {
            return tmp;
        }
        
        if (tmp.size() < 12) {
            return new HashSet<PairInt>();
        }
        
        int[] minMaxXY = MiscMath.findMinMaxXY(tmp);
        log.info("rainbow range in x: " + minMaxXY[0] + " to " + minMaxXY[1]);
        
        // TODO: consider filtering here by the size the rainbow is in the image.
        
        log.info("rainbow - sun = " + tmp.size() + " points");
        
        PolynomialFitter polyFitter = new PolynomialFitter();
        //y = c0*1 + c1*x[i] + c2*x[i]*x[i]
        float[] coef = polyFitter.solveAfterRandomSampling(rainbowPoints);
        
        if (coef == null) {
            return new HashSet<PairInt>();
        }
        
        log.info("rainbow polynomial coefficients = " + Arrays.toString(coef));
        log.info("image dimensions are " + colorImg.getWidth() + " X " + 
            colorImg.getHeight() + " pixels^2");
        
        // expecting the polynomial to have coef[1] negative
        // and close to a 45 degree angle (dx/dy about 0.5)
        // but that's for an image with 0 degrees of camera rotation about z axis.
        //TODO: revisit this and add tests with rainbows and rotated camera images
        
        if (Math.abs((Math.abs(coef[1]) - 0.5)) > 0.25) {
            return new HashSet<PairInt>();
        }
        if (coef[1] > 0) {
            return new HashSet<PairInt>();
        }
        if (coef[2] > 0.01) {
            return new HashSet<PairInt>();
        }        
        
        polyFitter.plotFit(coef, rainbowPoints, 234, colorImg.getWidth(),
            colorImg.getHeight(), "rainbow points");
        
        return rainbowPoints;
    }

    private void addBackMissingZeros(Set<PairInt> zeroPoints,
        GreyscaleImage gXYImg, int binFactor, int valueToSubtract) {
        
        int width = gXYImg.getWidth();
        int height = gXYImg.getHeight();
        
        GreyscaleImage img = gXYImg.copyImage();
        MatrixUtil.add(img.getValues(), -1*valueToSubtract);
        
        Set<PairInt> addPoints = new HashSet<PairInt>();
        
        for (PairInt p : zeroPoints) {

            int x = p.getX();
            int y = p.getY();

            for (int c = (x - binFactor); c <= (x + binFactor); c++) {
                if ((c < 0) || (c > (width - 1))) {
                    continue;
                }
                for (int r = (y - binFactor); r <= (y + binFactor); r++) {
                    if ((r < 0) || (r > (height - 1))) {
                        continue;
                    }
                    if ((c == x) && (r == y)) {
                        continue;
                    }

                    int neighborIdx = img.getIndex(c, r);

                    Integer index = Integer.valueOf(neighborIdx);

                    if (addPoints.contains(index)) {
                        continue;
                    }
                    if (zeroPoints.contains(index)) {
                        continue;
                    }

                    int v = img.getValue(c, r);
                    if (v == 0) {
                        addPoints.add(new PairInt(c, r));
                    }
                }
            }
        }
        
        zeroPoints.addAll(addPoints);
    }

    private void subtractWithCorrectForNegative(GreyscaleImage gXY2, int subtract) {

        int nz = 0;
        
        if (subtract > 0) {
            
            for (int i = 0; i < gXY2.getNPixels(); i++) {
                int v = gXY2.getValue(i);
                v -= subtract;
                if (v < 0) {
                    v = 0;
                }
                gXY2.setValue(i, v);
                if (v == 0) {
                    nz++;
                }
            }
        }

        log.info("number of set 0's=" + nz);

    }

    public int extractEmbeddedGroupPoints(
        List<Set<PairInt>> embeddedGroups, Map<Integer, PairInt> rowColRange, 
        int[] rowMinMax, Set<PairInt> outputEmbeddedPoints,
        int xMinImage, int xMaxImage, int yMinImage, int yMaxImage) {
        
        int nCorrectedEmbeddedGroups = 0;
        
        PerimeterFinder finder = new PerimeterFinder();
        
        for (int gId = 0; gId < embeddedGroups.size(); gId++) {

            Set<PairInt> groupPoints = embeddedGroups.get(gId);

            int[] gRowMinMax = new int[2];
            Map<Integer, PairInt> gRowColRange = finder.find(groupPoints, 
                gRowMinMax);

            boolean unbounded = isPerimeterUnbound(gRowColRange, gRowMinMax,
                rowColRange, rowMinMax, xMinImage, xMaxImage, yMinImage, yMaxImage);

            if (!unbounded) {
                
                nCorrectedEmbeddedGroups++;
                
                outputEmbeddedPoints.addAll(groupPoints);
            }
        }
        
        return nCorrectedEmbeddedGroups;
    }
    
    private int extractEmbeddedGroupPoints(Set<PairInt> skyPoints, 
        GreyscaleImage gXY2, Set<PairInt> outputEmbeddedPoints) {
        
        PerimeterFinder finder = new PerimeterFinder();
        int[] rowMinMax = new int[2];
        Map<Integer, PairInt> rowColRange = finder.find(skyPoints, rowMinMax);
        DFSContiguousValueFinder contiguousNonZeroFinder = 
            new DFSContiguousValueFinder(gXY2);
        contiguousNonZeroFinder.findEmbeddedGroupsNotThisValue(0, rowColRange,
            rowMinMax);
        int nEmbeddedGroups = contiguousNonZeroFinder.getNumberOfGroups();

        List<Set<PairInt>> embeddedGroups = new ArrayList<Set<PairInt>>();
        
        for (int gId = 0; gId < nEmbeddedGroups; gId++) {
            
            Set<PairInt> groupPoints = new HashSet<PairInt>();

            contiguousNonZeroFinder.getXY(gId, groupPoints);
            
            embeddedGroups.add(groupPoints);
            
            outputEmbeddedPoints.addAll(groupPoints);
        }
                
        return nEmbeddedGroups;
    }

    /**
     * fill in the rightmost sky points that would be missing due to down sizing
     * the image to find sky points then up sizing the image to current scale.
     * 
     * @param skyPoints
     * @param binFactor the size of the former down sizing to find sky points.
     * @param skyRowColRange the minimum row number of sky points with respect
     * to the canny edge filter product images (theta, gradientXY).
     * @param skyRowMinMax the maximum row number of sky points with respect
     * to the canny edge filter product images (theta, gradientXY).
     * @param originalColorImage the original color image
     * @param xRelativeOffset the offset in x of the canny edge filter intermediate
     * product images from the reference frame of the originalColorImage.
     * @param yRelativeOffset the offset in y of the canny edge filter intermediate
     * product images from the reference frame of the originalColorImage.
     */
    private void fillInRightBoundarySkyPoints(Set<PairInt> skyPoints, int binFactor, 
        Map<Integer, PairInt> skyRowColRange, int[] skyRowMinMax,
        Image originalColorImage, int xRelativeOffset, int yRelativeOffset) {
        
        if (binFactor == 1) {
            return;
        }
        
        int lastCol = (originalColorImage.getWidth() - 1) - xRelativeOffset;
        
        for (int r = skyRowMinMax[0]; r <= skyRowMinMax[1]; r++) {
            
            final int row = r;
            final Integer rowIndex = Integer.valueOf(row);
            
            PairInt cRange = skyRowColRange.get(rowIndex);
            
            int rightCol = cRange.getY();
            
            if (rightCol < (lastCol - binFactor + 1)) {
                continue;
            }
            
            for (int i = 1; i <= binFactor; i++) {
                
                int x = rightCol + i;
                
                if ((x + xRelativeOffset) > (originalColorImage.getWidth() - 1)) {
                    break;
                }
                
                PairInt rightPoint = new PairInt(x, row);
                
                skyPoints.add(rightPoint);
                
                skyRowColRange.put(rowIndex, new PairInt(cRange.getX(), x));
            }
        }        
    }

    private void addPointsAndUpdateRowColRange(Set<PairInt> skyPoints, 
        Set<PairInt> sunPoints, Map<Integer, PairInt> skyRowColRange, 
        int[] skyRowMinMax) {
        
        if (skyPoints == null) {
            throw new IllegalArgumentException("skyPoints cannot be null");
        }
        if (sunPoints == null) {
            throw new IllegalArgumentException("sunPoints cannot be null");
        }
        if (skyRowColRange == null) {
            throw new IllegalArgumentException("skyRowColRange cannot be null");
        }
        if (skyRowMinMax == null) {
            throw new IllegalArgumentException("skyRowMinMax cannot be null");
        }
        
        if (sunPoints.isEmpty()) {
            return;
        }
        
        skyPoints.addAll(sunPoints);
        
        PerimeterFinder finder = new PerimeterFinder();
        Map<Integer, PairInt> skyRowColRange2 = finder.find(skyPoints, 
            skyRowMinMax);
        
        skyRowColRange.clear();
        
        skyRowColRange.putAll(skyRowColRange2);
        
    }

    /**
     * given seed skyPoints to start from, use conservative limits on contrast
     * and color difference to add neighbors to skyPoints.  The conservative
     * limits are meant to help avoid overrunning the skyline for low
     * contrast such as a hazy sky and snow covered peaks, for example.
     * The conservative limits do not necessarily find all sky points.
     * 
     * @param skyPoints
     * @param originalColorImage
     * @param mask
     * @param pixelColorsMap
     * @param skyColorsMap 
     */
    private void findClouds(Set<PairInt> skyPoints,
        Image originalColorImage, GreyscaleImage mask,
        Map<Integer, PixelColors> pixelColorsMap,
        Map<PairInt, Set<PixelColors> > skyColorsMap
        ) {
        
        int maskWidth = mask.getWidth();
        int maskHeight = mask.getHeight();
        
        java.util.Stack<PairInt> cloudStack = new java.util.Stack<PairInt>();
        
        Set<PairInt> candidateCloudPoints = new HashSet<PairInt>();
       
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();
              
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        // TODO: some adjustment should be made for the perimeter of the
        // sun points to not erode the skyline.  
                
        for (PairInt skyPoint : skyPoints) {
            
            int x = skyPoint.getX();
            int y = skyPoint.getY();
            Integer index = Integer.valueOf(mask.getIndex(x, y));
        
            for (int k = 0; k < dxs.length; k++) {
            
                int xx = x + dxs[k];
                int yy = y + dys[k];
                
                if ((xx < 0) || (xx > (maskWidth - 1)) || (yy < 0) || 
                    (yy > (maskHeight - 1))) {
                    continue;
                }
                
                PairInt p = new PairInt(xx, yy);
                
                if (skyPoints.contains(p)) {
                    continue;
                }
                
                //add colors for index to the local sky points entry for p
                populatePixelColorMaps(p, index, originalColorImage, 
                    xOffset, yOffset, pixelColorsMap, skyColorsMap);
        
                if (candidateCloudPoints.contains(p)) {
                    continue;
                }
                
                // including perimeter points too (they're missing from  
                // skyPoints due to convolution widening of gradientXY
                // features).  see caveat above sunPoints.
                // so it's safe to add these without a color or contrast check:
                candidateCloudPoints.add(p);
                
                cloudStack.add(p);
            }
        }
    
        Set<PairInt> r0 = new HashSet<PairInt>();
        Set<PairInt> r1 = new HashSet<PairInt>();
        Set<PairInt> r2 = new HashSet<PairInt>();
        Set<PairInt> r3 = new HashSet<PairInt>();
        
        GroupPixelColors allSkyColor = new GroupPixelColors(skyPoints,
            originalColorImage, xOffset, yOffset);
        
        double rDivB = allSkyColor.getAvgRed() / allSkyColor.getAvgBlue();
        boolean skyIsRed = (rDivB > 1);
       
        log.info("==> r/b=" + rDivB
            + " redStdev=" + allSkyColor.getStdDevRed()
            + " blueStDev=" + allSkyColor.getStdDevBlue());
        
        Map<PairInt, GroupPixelColors> localSkyColors = new
            HashMap<PairInt, GroupPixelColors>();
        
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(cloudStack.peek());
       
        dxs = new int[]{-1,  0, 1, 0};
        dys = new int[]{ 0, -1, 0, 1};
        
        while (!cloudStack.isEmpty()) {

            PairInt uPoint = cloudStack.pop();
            
            int uX = uPoint.getX();
            int uY = uPoint.getY();

            //(1 + frac)*O(N) where frac is the fraction added back to stack
                      
            for (int k = 0; k < dxs.length; k++) {
                            
                int vX = uX + dxs[k];
                int vY = uY + dys[k];
                
                if ((vX < 0) || (vX > (maskWidth - 1)) || (vY < 0) || 
                    (vY > (maskHeight - 1))) {
                    continue;
                }
            
                PairInt vPoint = new PairInt(vX, vY);

                if (visited.contains(vPoint) || skyPoints.contains(vPoint)
                    || candidateCloudPoints.contains(vPoint)) {
                    continue;
                }

                visited.add(vPoint);
                    
                GroupPixelColors localSky = localSkyColors.get(uPoint);
                if (localSky == null) {
                    // this should never be null:
                    Set<PixelColors> skyColors = skyColorsMap.get(uPoint);
                    localSky = new GroupPixelColors(skyColors);
                    localSkyColors.put(uPoint, localSky);
                }

                int rV = originalColorImage.getR(vX + xOffset, vY + yOffset);
                int gV = originalColorImage.getG(vX + xOffset, vY + yOffset);
                int bV = originalColorImage.getB(vX + xOffset, vY + yOffset);

                float totalRGBV = rV + gV + bV;
                
                /*is this within color and contrast limits?  if yes,
                   add to candidateCloudPoints 
                   add to cloudStack
                   add skyColors to candidateSkyColorsMap for vPoint
                   add localSky to localSkyColors for vPoint
                */

                double contrastV = localSky.calcContrastToOther(rV, gV, bV);

                double colorDiffV = localSky.calcColorDiffToOther(rV, gV, bV);

                double skyStDevContrast = localSky.getStdDevContrast();

                double skyStDevColorDiff = localSky.getStdDevColorDiff();
                    
                /*log.info(String.format(
                    "(%d, %d) contrast=%f clrDiff=%f contrastStDev=%f clrStDev=%f",
                    vX + xOffset, vY + yOffset,
                    contrastV, colorDiffV, 
                    skyStDevContrast,
                    localSky.getStandardDeviationColorDifference())
                );*/

                boolean doNotAddToStack = false;
                 
 //TODO: this needs adjustments...
                float rPercentV = (float)rV/totalRGBV;
                float gPercentV = (float)gV/totalRGBV;
                float bPercentV = (float)bV/totalRGBV;
                
                boolean isBrown = (Math.abs(rPercentV - 0.5) < 0.4)
                    && (Math.abs(gPercentV - 0.32) < 0.1)
                    && (Math.abs(bPercentV - 0.17) < 0.1);
/*            
if (((vX + xOffset) > 347) && ((vX + xOffset) < 351) && ((vY + yOffset) > 289) && ((vY + yOffset) < 291)) {
    log.info(String.format(
        "(%d,%d) k=%d : rgb=(%d,%d,%d)   contrastV=%f  colorDiffV=%f stdevContrast=%f stdevColorDiff=%f  percent rgb=(%f,%f,%f)  isBrown=%b",
        vX + xOffset, vY + yOffset, k,
        rV, gV, bV,
        contrastV, colorDiffV, skyStDevContrast, skyStDevColorDiff,
        rPercentV, gPercentV, bPercentV, isBrown
        ));
}*/               
                if (isBrown) {
                    
                    // trying to skip over foreground such as land or sunset + water
                    
                    float[] hsb = new float[3];
                    Color.RGBtoHSB(rV, gV, bV, hsb);
                    
                    if ((colorDiffV > 15*skyStDevColorDiff) && (hsb[2] < 0.5)) {
                        
                        if (hsb[2] > 0.4) {
                            
                            r0.add(vPoint);
                            
                            continue;
                            
                        } else if ((colorDiffV > 50*skyStDevColorDiff) || 
                            (Math.abs(contrastV) > 10.*Math.abs(skyStDevContrast))
                            ) {
                            
                            r1.add(vPoint);
                            
                            continue;
                        }
                    }
                }
                
                if (skyIsRed) {
                    
                    // if contrast is '+' and large, may be the skyline boundary
                    if (
                        Double.isInfinite(skyStDevContrast)
                        || (
                            !Double.isInfinite(skyStDevContrast)
                            && (skyStDevContrast != 0.)
                            &&
                            (
                                ((contrastV > 0.1) && (Math.abs(contrastV) > 3.*skyStDevContrast))
                                ||
                                ((contrastV > 0) && (colorDiffV > 15*skyStDevColorDiff))
                            )
                        )
                        ) {

                        r2.add(vPoint);

                        continue;
                        
                    } else if (skyStDevContrast == 0.) {

                        if (contrastV >= 0.) {
                            doNotAddToStack = true;
                        }

                    } else {

                        //TODO:  if there are sun points, need a zone of
                        // avoidance to not erode the foreground 
                    }

                } else {

                    //TODO: blue filters

                    continue;
                }
                
                candidateCloudPoints.add(vPoint);

                if (!doNotAddToStack) {
                    cloudStack.add(vPoint);
                }

                localSkyColors.put(vPoint, localSky);
                skyColorsMap.put(vPoint, skyColorsMap.get(uPoint));

            }
        }
        
debugPlot(r0, r1, r2, r3, originalColorImage, xOffset, yOffset, "filtered_out_of_clouds");        
        
        skyPoints.addAll(candidateCloudPoints);
        
        for (PairInt p : skyPoints) {
            int x = p.getX();
            int y = p.getY();            
            mask.setValue(x, y, 0);
        }
    }

    private void findMoreClouds(Set<PairInt> skyPoints, Image originalColorImage, 
        GreyscaleImage mask) {
        
        /* for cloudy skies, findClouds can find the large majority of sky,
        but sometimes there are dark cloud bands separating brighter
        clouds or sky from the majority of connected sky pixels.
        findMoreClouds attempts to connect the nearby, but unconnected sky
        with the large majority of connected sky by looking at the color
        properties of skyPoints as 3 separate groups.
        */
        
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();
        
        java.util.Stack<PairInt> stack = new java.util.Stack<PairInt>();
        
        //O(N_sky)
        for (PairInt p : skyPoints) {
            stack.add(p);
        }
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(stack.peek());
        
        int[] counts = new int[3];
        GroupPixelColors[] allSkyColors = partitionInto3ByColorDifference(skyPoints,
            originalColorImage, xOffset, yOffset, counts);
        int maxPartionIdx = MiscMath.findYMaxIndex(counts);
        
        int maskWidth = mask.getWidth();
        int maskHeight = mask.getHeight();
        
        //int[] dxs = new int[]{-1,  0, 1, 0};
        //int[] dys = new int[]{ 0, -1, 0, 1};
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        while (!stack.isEmpty()) {
            
            PairInt skyPoint = stack.pop();
            
            int x = skyPoint.getX();
            int y = skyPoint.getY();
        
            for (int k = 0; k < dxs.length; k++) {
                    
                int dx = dxs[k];
                int dy = dys[k];
                
                int xx = x + dx;
                int yy = y + dy;
                
                if ((xx < 0) || (xx > (maskWidth - 1)) || (yy < 0) || 
                    (yy > (maskHeight - 1))) {
                    continue;
                }
                               
                PairInt p = new PairInt(xx, yy);
                
                if (skyPoints.contains(p) || visited.contains(p)) {
                    continue;
                }
               
                visited.add(p);
                
                xx += xOffset;
                yy += yOffset;
                
                int rV = originalColorImage.getR(xx, yy);
                int gV = originalColorImage.getG(xx, yy);
                int bV = originalColorImage.getB(xx, yy);
  
                int ii = maxPartionIdx;
                float clrDiff = allSkyColors[ii].calcColorDiffToOther(rV, gV, bV);
                float avgClrDiff = allSkyColors[ii].getAvgColorDiff();
                double stDev = allSkyColors[ii].getStdDevColorDiff();
                // between 10.0 and 15.0
                if (Math.abs(clrDiff - avgClrDiff) <= 15.0*stDev) {
                    float contrast = allSkyColors[ii].calcContrastToOther(rV, gV, bV);
                    float avgContrast = allSkyColors[ii].getAvgContrast();
                    stDev = allSkyColors[ii].getStdDevContrast();
                    if (Math.abs(contrast - avgContrast) <= 15.0*stDev) {

                        skyPoints.add(p);

                        mask.setValue(xx, yy, 0);

                        stack.add(p);
                        
                        break;
                    }
                }
         
                if (false)
                log.info(String.format("(%d, %d) dClr=%f,%f,%f   dContrast=%f,%f,%f   (%f,%f,%f) (%f,%f,%f)",
                    xx, yy, 
                    allSkyColors[0].calcColorDiffToOther(rV, gV, bV),
                    allSkyColors[1].calcColorDiffToOther(rV, gV, bV),
                    allSkyColors[2].calcColorDiffToOther(rV, gV, bV),
                    
                    allSkyColors[0].calcContrastToOther(rV, gV, bV),
                    allSkyColors[1].calcContrastToOther(rV, gV, bV),
                    allSkyColors[2].calcContrastToOther(rV, gV, bV),
                    
                    (allSkyColors[0].calcColorDiffToOther(rV, gV, bV) 
                        - allSkyColors[0].getAvgColorDiff() ) /
                        allSkyColors[0].getStdDevColorDiff(),
                    (allSkyColors[1].calcColorDiffToOther(rV, gV, bV) 
                        - allSkyColors[1].getAvgColorDiff() ) /
                        allSkyColors[1].getStdDevColorDiff(),
                    (allSkyColors[2].calcColorDiffToOther(rV, gV, bV) 
                        - allSkyColors[2].getAvgColorDiff() ) /
                        allSkyColors[2].getStdDevColorDiff(),
                    
                    (allSkyColors[0].calcContrastToOther(rV, gV, bV)
                        - allSkyColors[0].getAvgContrast()) /
                        allSkyColors[0].getStdDevContrast(),
                    (allSkyColors[1].calcContrastToOther(rV, gV, bV)
                        - allSkyColors[1].getAvgContrast()) /
                        allSkyColors[1].getStdDevContrast(),
                    (allSkyColors[2].calcContrastToOther(rV, gV, bV)
                        - allSkyColors[2].getAvgContrast()) /
                        allSkyColors[2].getStdDevContrast()
                    )
                );
            }
        }
        
    }

    private GroupPixelColors[] partitionInto3ByColorDifference(Set<PairInt> skyPoints, 
        Image originalColorImage, int xOffset, int yOffset, int[] outputCounts) {
         
        GroupPixelColors allSkyColor = new GroupPixelColors(skyPoints,
            originalColorImage, xOffset, yOffset);
        
        int i = 0;
        float[] colorDiffs = new float[skyPoints.size()];
        PairInt[] ps = new PairInt[colorDiffs.length];
        for (PairInt p : skyPoints) {
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            int r = originalColorImage.getR(x, y);
            int g = originalColorImage.getG(x, y);
            int b = originalColorImage.getB(x, y);
            colorDiffs[i] = allSkyColor.calcColorDiffToOther(r, g, b);
            ps[i] = p;
            i++;
        }
        
        float minColorDiff = MiscMath.findMin(colorDiffs);
        float maxColorDiff = MiscMath.findMax(colorDiffs);
        float[] yErr = Errors.populateYErrorsBySqrt(colorDiffs);
        HistogramHolder h = Histogram.createSimpleHistogram(minColorDiff,
            maxColorDiff, 3, colorDiffs, yErr);
        
        for (i = 0; i < 3; i++) {
            outputCounts[i] = h.getYHist()[i];
        }

        Set<PairInt> set0 = new HashSet<PairInt>();
        Set<PairInt> set1 = new HashSet<PairInt>();
        Set<PairInt> set2 = new HashSet<PairInt>();
        
        float binSize = h.getXHist()[1] - h.getXHist()[0];
        
        for (i = 0; i < colorDiffs.length; i++) {
            float cd = colorDiffs[i];
            int binN = (int)((cd - minColorDiff)/binSize);
            switch(binN) {
                case 0:
                    set0.add(ps[i]);
                    break;
                case 1:
                    set1.add(ps[i]);
                    break;
                default:
                    set2.add(ps[i]);
                    break;
            }
        }
        
        GroupPixelColors[] sets = new GroupPixelColors[3];
        sets[0] = new GroupPixelColors(set0, originalColorImage, xOffset, yOffset);
        sets[1] = new GroupPixelColors(set1, originalColorImage, xOffset, yOffset);
        sets[2] = new GroupPixelColors(set2, originalColorImage, xOffset, yOffset);
        
        return sets;
    }
    
    private GroupPixelColors[] partitionInto3ByBrightness(
        Set<PairInt> skyPoints, 
        Image originalColorImage, int xOffset, int yOffset, HistogramHolder[] h) {
         
        int i = 0;
        float[] brightness = new float[skyPoints.size()];
        PairInt[] ps = new PairInt[brightness.length];
        for (PairInt p : skyPoints) {
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            int r = originalColorImage.getR(x, y);
            int g = originalColorImage.getG(x, y);
            int b = originalColorImage.getB(x, y);
            
            float[] hsb = new float[3];
            Color.RGBtoHSB(r, g, b, hsb);
        
            brightness[i] = hsb[2];
            ps[i] = p;
            i++;
        }
        
        float min = MiscMath.findMin(brightness);
        float max = MiscMath.findMax(brightness);
        float[] yErr = Errors.populateYErrorsBySqrt(brightness);
        h[0] = Histogram.createSimpleHistogram(min, max, 3, 
            brightness, yErr);
        
try {
    h[0].plotHistogram("sky brightness", 235);
} catch(Exception e) {
    
}     

        Set<PairInt> set0 = new HashSet<PairInt>();
        Set<PairInt> set1 = new HashSet<PairInt>();
        Set<PairInt> set2 = new HashSet<PairInt>();
        
        float binSize = h[0].getXHist()[1] - h[0].getXHist()[0];
        
        for (i = 0; i < brightness.length; i++) {
            float cd = brightness[i];
            int binN = (int)((cd - min)/binSize);
            switch(binN) {
                case 0:
                    set0.add(ps[i]);
                    break;
                case 1:
                    set1.add(ps[i]);
                    break;
                default:
                    set2.add(ps[i]);
                    break;
            }
        }
        
        GroupPixelColors[] sets = new GroupPixelColors[3];
        sets[0] = new GroupPixelColors(set0, originalColorImage, xOffset, yOffset);
        sets[1] = new GroupPixelColors(set1, originalColorImage, xOffset, yOffset);
        sets[2] = new GroupPixelColors(set2, originalColorImage, xOffset, yOffset);
        
        return sets;
    }
    
    private void findEmbeddedClouds(Set<PairInt> embeddedPoints, 
        Set<PairInt> skyPoints, Set<PairInt> excludeThesePoints,
        Image originalColorImage, GreyscaleImage mask,
        Map<Integer, PixelColors> pixelColorsMap,
        Map<PairInt, Set<PixelColors> > skyColorsMap,
        HistogramHolder[] brightnessHistogram,
        GroupPixelColors[] skyPartitionedByBrightness) {
        
        findClouds2(embeddedPoints, skyPoints, excludeThesePoints,
            originalColorImage, mask,
            pixelColorsMap, skyColorsMap,
            brightnessHistogram, skyPartitionedByBrightness, true, true);
    }
    
    private void findSeparatedClouds(Set<PairInt> sunPoints, 
        Set<PairInt> skyPoints, Set<PairInt> excludeThesePoints,
        Image originalColorImage, GreyscaleImage mask,
        Map<Integer, PixelColors> pixelColorsMap,
        Map<PairInt, Set<PixelColors> > skyColorsMap,
        HistogramHolder[] brightnessHistogram,
        GroupPixelColors[] skyPartitionedByBrightness) {
        
        findClouds2(sunPoints, skyPoints, excludeThesePoints,
            originalColorImage, mask,
            pixelColorsMap, skyColorsMap,
            brightnessHistogram, skyPartitionedByBrightness, false, false);
    }

    /**
     * Find clouds using given data and one of 2 methods:
     * (1) if "findOnlyConnected" is true, pointsToConnectOrExclude should hold 
     * only the points one wants to find if they pass color and contrast filters
     * and from there, any adjacent neighbors passing the filters will be
     * found and also added to the skyPoints;
     * or (2) if "findOnlyConnected" is false, pointsToConnectOrExclude should
     * hold any points that should be excluded if found as a candidate 
     * neighbor while the search is over each pixel in originalColorImage.
     * 
     * Note that in case (2), skyPoints and pointsToConnectOrExclude might
     * be the same point set.
     * @param pointsToConnectOrExclude
     * @param skyPoints
     * @param excludeThesePoints
     * @param originalColorImage
     * @param mask
     * @param brightnessHistogram
     * @param skyPartitionedByBrightness
     * @param pointsAreEmbeddedInSky 
     * @param pointIsEmbeddedInSky if true, the pixel filter will add a pass
     * for bluish or pinkish white clouds, else false.  it should be safe to
     * use this as true when searching for clouds embedded already in sky that
     * could not be confused with snowy mountain peaks or foreground.
     */
    private void findClouds2(Set<PairInt> pointsToConnectOrExclude, 
        final Set<PairInt> skyPoints, Set<PairInt> excludeThesePoints,
        Image originalColorImage, GreyscaleImage mask,
        Map<Integer, PixelColors> pixelColorsMap,
        Map<PairInt, Set<PixelColors> > skyColorsMap,
        HistogramHolder[] brightnessHistogram,
        GroupPixelColors[] skyPartitionedByBrightness,
        boolean findOnlyConnected, boolean pointsAreEmbeddedInSky) {
  
        if (pointsToConnectOrExclude == null) {
            if (!findOnlyConnected) {
                pointsToConnectOrExclude = new HashSet<PairInt>();
            } else {
                throw new IllegalArgumentException(
                    "if findOnlyConnected is true, then " + 
                    "pointsToConnectOrExclude has to contain points that " +
                    "hold the points to connect to");
            }
        }
        
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();
     
        /*
        scattering:
            for atmospheric particles, their size is much smaller than
            lambda, the wavelength of optical light:  
                Rayleigh scattering is prop to lambda^-4.
                this is responsible for blue or red skies compared to 
                the yellow color of the sun 
                (photosphere peak is near 5500 Angstroms).
        
            for clouds, the water droplets are the scatters and their
            size is comparable to lambda for optical light:
                Mie scattering affects optical colors roughly equally,
                so leads to less light without a color change.
        
        One could use the position of the sun to approximate the airmass for
        different depths of features and calculate the transmission of the solar 
        spectrum through atmosphere with absorption by water, scattering and 
        absorption by aerosols and Rayleigh scattering in addition to Mie 
        scattering by water droplets in the clouds
        combined with a model for water vapor in the atmosphere and a range 
        of cloud optical depths and locations and a camera response function 
        to estimate the colors expected in images...
        
        The source function upon the clouds could be modeled with
        http://rredc.nrel.gov/solar/pubs/spectral/model/section5.html
        
        Since that isn't feasible for now, could take skyPoints colors and colors
        from the outer patch of sun and look for some variations of those in the
        image close to, but unconnected to existing skyPoints.
        */
     
        // find all image pixels that are not in skyPoints or sunPoints
        // that are within the color range of brightest and next brightest sky
        
        Set<PairInt> cloudPoints = new HashSet<PairInt>();
        
        double rDivB = skyPartitionedByBrightness[2].getAvgRed() / 
            skyPartitionedByBrightness[2].getAvgBlue();
        
        boolean skyIsRed = (rDivB > 1);
        
        boolean skyIsPurple = skyIsRed && (rDivB < 1.5) &&
            (Math.abs(
                skyPartitionedByBrightness[2].getAvgGreen() - 
                skyPartitionedByBrightness[2].getAvgBlue()) 
                < 
                0.1 * skyPartitionedByBrightness[2].getAvgBlue()
            );
        
        boolean hasDarkGreyClouds = false;
        
        //only perform k=0 if the sky has a narrow dark section
        boolean useKEqualsZero = 
            (skyPartitionedByBrightness[0].getStdDevGreen() < 10) && 
            (skyPartitionedByBrightness[0].getStdDevBlue() < 10) && 
            (skyPartitionedByBrightness[0].getStdDevContrast() < 10)
            && (skyPartitionedByBrightness[0].getStdDevColorDiff()
            < 10);
        if (useKEqualsZero) {
            if (skyIsRed) {
                useKEqualsZero = useKEqualsZero && 
                (skyPartitionedByBrightness[0].getStdDevRed() < 20);
            } else {
                useKEqualsZero = useKEqualsZero && 
                (skyPartitionedByBrightness[0].getStdDevRed() < 10);
            }
        } else {
            if (!skyIsRed) {
                // if the sky is grey, may have dark clouds
                //TODO: may need to edit this as foreground objects can be sim in color too
                float gDivR = skyPartitionedByBrightness[0].getAvgGreen()/
                    skyPartitionedByBrightness[0].getAvgRed();
                float gDivB = skyPartitionedByBrightness[0].getAvgGreen()/
                    skyPartitionedByBrightness[0].getAvgBlue();
                if ((gDivR >= 1.0) && (gDivR <= 1.25) && (gDivB >= 0.85)
                    && (gDivB <= 1.00) 
                    && (skyPartitionedByBrightness[0].getAvgGreen() > 100)) {
                    
                    useKEqualsZero = true;
                    hasDarkGreyClouds = true;
                }
            }
        }
       
        java.util.Stack<PairInt> stack = new java.util.Stack<PairInt>();
        
        if (!findOnlyConnected) {
        
            for (int col = 0; col < originalColorImage.getWidth(); col++) {
                for (int row = 0; row < originalColorImage.getHeight(); row++) {

                    PairInt p = new PairInt(col - xOffset, row - yOffset);
                    
                    if (excludeThesePoints.contains(p)) {
                        continue;
                    }
                    
                    filterToAddCloudPixel(pointsToConnectOrExclude, skyPoints, 
                        originalColorImage, xOffset, yOffset, p,
                        brightnessHistogram, skyPartitionedByBrightness,
                        skyIsRed, skyIsPurple, hasDarkGreyClouds, useKEqualsZero,
                        pointsAreEmbeddedInSky,
                        cloudPoints, stack
                    );                    
                }
            }
            
            if (stack.isEmpty()) {
                return;
            }

            // this block does not use the stack further because it has already
            // searched each pixel of the image
            
            // find the number of groups of connected pixels within skyPoints and
            // extSkyPoints

debugPlot(cloudPoints, originalColorImage, xOffset, yOffset, "exterior_clouds");
        
            skyPoints.addAll(cloudPoints);

            for (PairInt p : cloudPoints) {
                int x = p.getX();
                int y = p.getY();            
                mask.setValue(x, y, 0);
            }

            int n = cloudPoints.size();

            log.info("NUMBER of separated cloud points = " + n);

            DFSConnectedGroupsFinder groupsFinder = new DFSConnectedGroupsFinder();
            groupsFinder.setMinimumNumberInCluster(100);
            groupsFinder.findConnectedPointGroups(skyPoints, mask.getWidth(), 
                mask.getHeight());

            log.info("NUMBER of groups of connected skyPoints=" 
                + groupsFinder.getNumberOfGroups());

            if (groupsFinder.getNumberOfGroups() > 2 /*1*/) {
                
                findClouds(skyPoints, originalColorImage, mask, 
                    pixelColorsMap, skyColorsMap);
                
                for (PairInt p : skyPoints) {
                    int x = p.getX();
                    int y = p.getY();
                    mask.setValue(x, y, 0);
                }
            }
        
        } else {
        
            // == case (2)  find only points in or connected to embedded points
            
            // find the points in and connected to pointsToConnectOrExclude
            Set<PairInt> excludePoints = new HashSet<PairInt>();
            
            for (PairInt uPoint : pointsToConnectOrExclude) {
                
                // populate stack and cloud points with the initial points
                // that pass the color and contrast filter:
                
log.info("uPoint=(" + uPoint.getX() + ", " + uPoint.getY() + ")");
 
                if (excludePoints.contains(uPoint)) {
                    continue;
                }
                
                filterToAddCloudPixel(excludePoints, skyPoints, 
                    originalColorImage, xOffset, yOffset, uPoint,
                    brightnessHistogram, skyPartitionedByBrightness,
                    skyIsRed, skyIsPurple, hasDarkGreyClouds, useKEqualsZero,
                    pointsAreEmbeddedInSky,
                    cloudPoints, stack
                );
            }
            
            if (stack.isEmpty()) {
                return;
            }
            
            Set<PairInt> visited = new HashSet<PairInt>();
            visited.add(stack.peek());
            
            // find the neighbors of the points in the stack
            
            while (!stack.isEmpty()) {
                
                PairInt uPoint = stack.pop();
                
                int uX = uPoint.getX();
                int uY = uPoint.getY();

                for (int vX = (uX - 1); vX <= (uX + 1); vX++) {

                    for (int vY = (uY - 1); vY <= (uY + 1); vY++) {
                    
                        PairInt vPoint = new PairInt(vX, vY);
                        
                        if (vPoint.equals(uPoint)) {
                            continue;
                        }
                        
                        if (excludePoints.contains(vPoint)) {
                            continue;
                        }
                        
                        if (visited.contains(vPoint)) {
                            continue;
                        }

                        filterToAddCloudPixel(excludePoints, skyPoints, 
                            originalColorImage, xOffset, yOffset, vPoint,
                            brightnessHistogram, skyPartitionedByBrightness,
                            skyIsRed, skyIsPurple, hasDarkGreyClouds, 
                            useKEqualsZero, pointsAreEmbeddedInSky,
                            cloudPoints, stack
                        );

                        visited.add(vPoint);
                    }
                }
            }
            
debugPlot(cloudPoints, originalColorImage, xOffset, yOffset, "embedded_clouds");
        
            skyPoints.addAll(cloudPoints);
            
debugPlot(skyPoints, originalColorImage, xOffset, yOffset, "sky_added_embedded_clouds");

            for (PairInt p : cloudPoints) {
                int x = p.getX();
                int y = p.getY();            
                mask.setValue(x, y, 0);
            }
            
        } // end case (2)
        
    }

    private GroupPixelColors calculateColorOfOuterEdgeOfSun(
        Set<PairInt> sunPoints, Image originalColorImage, int xOffset, int yOffset) {
       
        EllipseHelper ellipseHelper = new EllipseHelper();
        
        double[] params = ellipseHelper.fitEllipseToPoints(sunPoints);
        
        if (params == null) {
            // not close to an ellipse
            return null;
        }
        
        float xc = (float)params[0];
        float yc = (float)params[1];
        float a = (float)params[2];
        float b = (float)params[3];
        float alpha = (float)params[4];
        
        int xMin = Integer.MAX_VALUE;
        int xMax = Integer.MIN_VALUE;
        int yMin = Integer.MAX_VALUE;
        int yMax = Integer.MIN_VALUE;
       
        for (PairInt p : sunPoints) {
            int x = p.getX();
            int y = p.getY();
            if (x < xMin) {
                xMin = x;
            }
            if (y < yMin) {
                yMin = y;
            }
            if (x > xMax) {
                xMax = x;
            }
            if (y > yMax) {
                yMax = y;
            }
        }
        
        double radius = (alpha < 0.5) ? 0.5*(xMax - xMin) : Math.sqrt(a*a - b*b);
        
        double limit = radius - 1;
        
        Set<PairInt> outer = new HashSet<PairInt>();
        for (PairInt p : sunPoints) {
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            float dx = xc - x;
            float dy = yc - y;
                
            double dist = Math.sqrt(dx*dx + dy*dy);
            
            if (dist >= limit) {
                outer.add(p);
            }
        }
        
        return new GroupPixelColors(outer, originalColorImage, xOffset, yOffset);
    }

    private int count(List<PairIntArray> zeroPointLists) {
        
        int n = 0;
        for (PairIntArray p : zeroPointLists) {
            n += p.getN();
        }
        
        return n;
    }

    private void debugPlot(Set<PairInt> extSkyPoints, Image originalColorImage, 
        int xOffset, int yOffset, String outputPrefixForFileName) {
        
        //plot is made in aspects
        
    }
    
    private void debugPlot(Set<PairInt> r0, Set<PairInt> r1, Set<PairInt> r2, 
        Set<PairInt> r3, Image originalColorImage, int xOffset, int yOffset, 
        String filtered_out_of_clouds) {
        
        //plot is made in aspects
    }

    private void filterToAddCloudPixel(Set<PairInt> excludeFromThesePoints, 
        Set<PairInt> skyPoints, Image origColorImg, 
        int xOffset, int yOffset, PairInt p,
        HistogramHolder[] brightnessHistogram, 
        GroupPixelColors[] skyBinsByBrightness,
        boolean skyIsRed, boolean skyIsPurple, boolean hasDarkGreyClouds,
        boolean useKEqualsZero, boolean pointIsEmbeddedInSky,
        Set<PairInt> outputCloudPoints, Stack<PairInt> outputStack) {
        
        if (skyPoints.contains(p) || excludeFromThesePoints.contains(p)) {
            return;
        }
        
        int col = p.getX() + xOffset;
        int row = p.getY() + yOffset;
        
        if ((col < 0) || (col > (origColorImg.getWidth() - 1)) ||
            (row < 0) || (row > (origColorImg.getHeight() - 1))) {
            return;
        }

        int r = origColorImg.getR(col, row);
        int g = origColorImg.getG(col, row);
        int b = origColorImg.getB(col, row);
        
        // because of Rayleigh scattered illumination, need large tolerance
        int limit = 55;
        boolean isGrey = (Math.abs(r - g) < limit) && (Math.abs(g - b) < limit) 
            && (Math.abs(r - b) < 1.8*limit);

        float[] hsb = new float[3];
        Color.RGBtoHSB(r, g, b, hsb);

        for (int k = 2; k >= 0; k--) {                    

            if ((k == 0) && !useKEqualsZero) {
                continue;
            } 

            double contrast = skyBinsByBrightness[k].calcContrastToOther(r, g, b);

            double diffContrast = contrast - 
                skyBinsByBrightness[k].getAvgContrast();
            
            if (skyIsRed && (diffContrast > 0) && !pointIsEmbeddedInSky) {
                continue;
            }
            
            double diffR = r - skyBinsByBrightness[k].getAvgRed();
            if (diffR < 0) {
                diffR *= -1;
            }
            double diffG = g - skyBinsByBrightness[k].getAvgGreen();
            if (diffG < 0) {
                diffG *= -1;
            }
            double diffB = b - skyBinsByBrightness[k].getAvgBlue();
            if (diffB < 0) {
                diffB *= -1;
            }
          
if ((col > 95) && (col < 97) && (row > 285) && (row < 287)) {
log.info(String.format(
"(%d,%d) k=%d : rgb=(%d,%d,%d) skyRGB=(%f,%f,%f) diffRGB=(%f,%f,%f) stdevRed=%f stdevGreen=%f stdevBlue=%f  contrast=%f  diffContrast=%f",
col, row, k,
r, g, b,
diffR, diffG, diffB,
skyBinsByBrightness[k].getAvgRed(),
skyBinsByBrightness[k].getAvgGreen(),
skyBinsByBrightness[k].getAvgBlue(),
skyBinsByBrightness[k].getStdDevRed(),
skyBinsByBrightness[k].getStdDevGreen(),
skyBinsByBrightness[k].getStdDevBlue(),
contrast, diffContrast
));
int z = 1;
}
log.info("xy=" + col + " , " + row);

            if (
                (pointIsEmbeddedInSky /*&& isGrey*/) ||
                (
                (
                    (diffR <= skyBinsByBrightness[k].getStdDevRed())
                    ||
                    (skyIsRed && (r > skyBinsByBrightness[k].getAvgRed()))
                    ||
                    (skyIsRed && (r > 155) && 
                        (diffR <= 3.5*skyBinsByBrightness[k].getStdDevRed())
                    )
                    ||
                    (!skyIsRed && hasDarkGreyClouds && 
                        (diffR <= 2.0*skyBinsByBrightness[k].getStdDevRed())
                    )
                )
                && (
                    (diffG <= 1.5*skyBinsByBrightness[k].getStdDevGreen())
                    ||
                    //consider (contrast < 0.) && (contrast > -0.05)
                    (skyIsRed && (contrast < 0.) && (g > 130) && 
                        (diffG <= 2.0*skyBinsByBrightness[k].getStdDevGreen())
                    )
                    ||
                    ((skyIsPurple || skyIsRed) && (contrast > 0.) && (g < 130) && 
                        (diffG <= 2.0*skyBinsByBrightness[k].getStdDevGreen())
                    )
                    ||
                    (!skyIsRed && hasDarkGreyClouds && 
                        (diffG <= 2.0*skyBinsByBrightness[k].getStdDevGreen())
                    )
                )
                && (
                    (diffB <= 1.2*skyBinsByBrightness[k].getStdDevBlue())
                    ||
                    (skyIsPurple && (contrast > 0.) && (b < 130) && 
                        (diffB <= 2.5*skyBinsByBrightness[k].getStdDevBlue())
                    )
                    ||
                    (!skyIsRed && hasDarkGreyClouds && 
                        (diffB <= 1.5*skyBinsByBrightness[k].getStdDevBlue())
                    )
                )
                )
                ) {

                outputCloudPoints.add(p);

                outputStack.add(p);

                break;
            }
        }
    }

    private void populatePixelColorMaps(Set<PairInt> points, 
        Image originalColorImage, GreyscaleImage mask, 
        Map<Integer, PixelColors> pixelColorsMap, 
        Map<PairInt, Set<PixelColors>> skyColorsMap) {
        
        int xOffset = mask.getXRelativeOffset(); 
        int yOffset = mask.getYRelativeOffset();
        
        for (PairInt p : points) {
            
            int index = mask.getIndex(p.getX(), p.getY());
            
            populatePixelColorMaps(p, index, 
                originalColorImage, xOffset, yOffset, 
                pixelColorsMap, skyColorsMap);
        }
    }

    /**
     * populate pixelColorsMap and skyColorsMap by adding the colors from
     * pixelIndex to skyColorsMap for key p while also adding the colors
     * from pixelIndex to pixelColorsMap.
     * @param p
     * @param pixelIndex
     * @param originalColorImage
     * @param xOffset
     * @param yOffset
     * @param pixelColorsMap
     * @param skyColorsMap 
     */
    private void populatePixelColorMaps(PairInt p, int pixelIndex,
        Image originalColorImage, int xOffset, int yOffset, 
        Map<Integer, PixelColors> pixelColorsMap, 
        Map<PairInt, Set<PixelColors>> skyColorsMap) {
        
        int x = p.getX();
        int y = p.getY();
        
        Set<PixelColors> skyColors = skyColorsMap.get(p);
        if (skyColors == null) {
            skyColors = new HashSet<PixelColors>();
            skyColorsMap.put(p, skyColors);
        }
        PixelColors skyPC = pixelColorsMap.get(pixelIndex);
        if (skyPC == null) {
            int r = originalColorImage.getR(x + xOffset, y + yOffset);
            int g = originalColorImage.getG(x + xOffset, y + yOffset);
            int b = originalColorImage.getB(x + xOffset, y + yOffset);
            skyPC = new PixelColors(r, g, b);
            pixelColorsMap.put(pixelIndex, skyPC);
        }
        skyColors.add(skyPC);
    }

    private void removeReflectedSun(List<PairIntArray> zeroPointLists, 
        Image colorImg, GreyscaleImage thetaImg) {
        
        List<Integer> remove = new ArrayList<Integer>();
        
        for (int i = 0; i < zeroPointLists.size(); i++) {
            
            PairIntArray zeroPointList = zeroPointLists.get(i);
            
            Set<PairInt> pointSet = new HashSet<PairInt>();
            
            for (int ii = 0; ii < zeroPointList.getN(); ii++) {
                int x = zeroPointList.getX(ii);
                int y = zeroPointList.getY(ii);
                PairInt p = new PairInt(x, y);
                pointSet.add(p);
            }
            
            Set<PairInt> sunPoints = findSunConnectedToSkyPoints(pointSet,
                colorImg, thetaImg.getXRelativeOffset(),
                thetaImg.getYRelativeOffset());
            
            float fracSun = (float)sunPoints.size()/(float)pointSet.size();
            
            if (fracSun >= 0.8) {
                remove.add(Integer.valueOf(i));
            }
        }
        
        if (!remove.isEmpty()) {
            for (int i = (remove.size() - 1); i > -1; i--) {
                int idx = remove.get(i).intValue();
                zeroPointLists.remove(idx);
            }
        }
    }

}
