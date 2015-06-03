package algorithms.imageProcessing;

import algorithms.compGeometry.clustering.KMeansPlusPlus;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.PairInt;
import algorithms.misc.Complex;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.util.Errors;
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
    
    public void applyLaplacianKernel(GreyscaleImage input) {
        
        IKernel kernel = new Laplacian();
        Kernel kernelXY = kernel.getKernel();
        
        float norm = kernel.getNormalizationFactor();
        
        applyKernel(input, kernelXY, norm);
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
        
        applyImageSegmentation(theta, 2);
        
        subtractMinimum(theta);
        
        convertToBinaryImage(theta);
        
        removeSpurs(theta);
        
        throw new UnsupportedOperationException("not ready for use yet");
        //return theta;
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
      
        ImageExt out = new ImageExt(w1, h1);
        
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
    
    /**
     * 
     * @param input
     * @param forward if true, apply FFT transform, else inverse FFT transform
     */
    public void apply2DFFT(GreyscaleImage input, boolean forward) {
        
        //TODO: apply padding to nearest power of 2
        
        // initialize matrix of complex numbers as real numbers from image
        Complex[][] cc = convertImage(input);
                
        Complex[][] ccOut = apply2DFFT(cc, forward);
        
        writeToImage(input, ccOut);
    }
    
    protected Complex[][] apply2DFFT(Complex[][] cc, boolean forward) {
              
        // perform FFT by column
        for (int col = 0; col < cc.length; col++) {
            if (forward) {
                cc[col] = FFT.fft(cc[col]);
            } else {
                cc[col] = FFT.ifft(cc[col]);
            }
        }
        
        
        //transpose the matrix
        cc = MatrixUtil.transpose(cc);
                
        // perform FFT by column (originally rows)
        for (int col = 0; col < cc.length; col++) {
            if (forward) {
                cc[col] = FFT.fft(cc[col]);
            } else {
                cc[col] = FFT.ifft(cc[col]);
            }
        }
        
        //transpose the matrix
        cc = MatrixUtil.transpose(cc);
        
        
        return cc;
    }
    
    /**
     * 
     * @param cc
     */
    public void writeToImage(GreyscaleImage img, Complex[][] cc) {

        img.fill(0);
        
        // write back to original image
        for (int col = 0; col < img.getWidth(); col++) {
            for (int row = 0; row < img.getHeight(); row++) {
                double re = cc[col][row].re();
                double a = cc[col][row].abs();
                double p = cc[col][row].phase();
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
     * write a 2-D complex array from the image
     * 
     * @param input
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
        
        CannyEdgeFilter cef = new CannyEdgeFilter();
        
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
        psfNorm = apply2DFFT(psfNorm, true);
        
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
                
        Complex[][] imgFFT = apply2DFFT(imgCC, true);
        
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
        
        GreyscaleImage img2 = img0.createWithDimensions();
        
        writePositiveRealToImage(img2, ccDeconv);
        
        ImageDisplayer.displayImage("FFT(img0)/FFT(PSF)", img2);
        
        
        Complex[][] inverse = apply2DFFT(ccDeconv, false);
        
        GreyscaleImage img4 = img0.createWithDimensions();
        
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
            return img;
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
        
        GreyscaleImage img2 = new GreyscaleImage(w2, h2);
            
        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {
                int v = img.getValue(col, row);
                img2.setValue(col, row, v);
            }
            if (h2 > h) {
                for (int row = h; row < h2; row++) {
                    img2.setValue(col, row, 0);
                }
            }
        }
                    
        if (!xIsPowerOf2) {
            for (int col = w; col < w2; col++) {
                for (int row = 0; row < h2; row++) {
                    img2.setValue(col, row, 0);
                }
            }
        }
            
        return img2;
    }
    
    /**
     * convert the color image into an image scaled into values 0 to 255
     * by the polar theta angle in CIE XY space.  
     * 
     * runtime complexity is O(N) + O(N*lg_2(N))
     * 
     * @param input
     * @return 
     */
    public GreyscaleImage convertToCIEXYPolarTheta(ImageExt input) {
       
        int w = input.getWidth();
        int h = input.getHeight();
        
        input.setRadiusForPopulateOnDemand(0);
        
        GreyscaleImage output = new GreyscaleImage(w, h);
        
        Map<PairInt, Float> pixThetaMap = new HashMap<PairInt, Float>();
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        float[] thetaValues = new float[input.getNPixels()];
        int thetaCount = 0;
        
        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {
                
                PairInt p = new PairInt(col, row);
                
                int r = input.getR(col, row);
                int g = input.getG(col, row);
                int b = input.getB(col, row);
                
                if ((r == 0) && (g == 0) && (b == 0)) {
                    continue;
                } else if ((r == 255) && (g == 255) && (b == 255)) {
                    output.setValue(col, row, 255);
                    continue;
                }
                
                float cieX = input.getCIEX(col, row);
                float cieY = input.getCIEY(col, row);
                
                if (cieC.isWhite(cieX, cieY)) {
                    output.setValue(col, row, 255);
                } else {
                    
                    double theta = cieC.calculateXYTheta(cieX, cieY);
                    
System.out.println(String.format("r,g,b=(%d,%d,%d) cieX,cieY=(%f,%f) theta=%f  pix=(%d,%d)", 
r, g, b, cieX, cieY, (float)theta, col, row));

                    thetaValues[thetaCount] = (float)theta;
                    
                    pixThetaMap.put(p, Float.valueOf((float)theta));
                    
                    thetaCount++;
                }
            }
        }
        
        thetaValues = Arrays.copyOf(thetaValues, thetaCount);
        
        float minValue = MiscMath.findMin(thetaValues);
        float maxValue = MiscMath.findMax(thetaValues);
                
        HistogramHolder hist = Histogram.createSimpleHistogram(minValue,
            maxValue, 254, thetaValues, 
            Errors.populateYErrorsBySqrt(thetaValues));
        
        try {
            hist.plotHistogram("cie XY theta histogram", 76543);
        } catch (Exception e) {}
        
        int nonZeroCount = 0;
        for (int i = 0; i < hist.getXHist().length; i++) {
            int c = hist.getYHist()[i];
            if (c > 0) {
                nonZeroCount++;
            }
        }
        
        float[] startBins = new float[nonZeroCount];
        
        float halfBinWidth = (hist.getXHist()[1] - hist.getXHist()[0])/2.f;
        
        nonZeroCount = 0;
        for (int i = 0; i < hist.getXHist().length; i++) {
            int c = hist.getYHist()[i];
            if (c > 0) {
                startBins[nonZeroCount] = hist.getXHist()[i] - halfBinWidth;
                nonZeroCount++;
            }
        }
        
        Iterator<Entry<PairInt, Float> > iter = pixThetaMap.entrySet().iterator();
        
        // O(N * lg_2(N))
        while (iter.hasNext()) {
            
            Entry<PairInt, Float> entry = iter.next();
            
            PairInt p = entry.getKey();
            
            float theta = entry.getValue().floatValue();
            
            int idx = Arrays.binarySearch(startBins, theta);
            
            // if it's negative, (-(insertion point) - 1)
            if (idx < 0) {
                // idx = -*idx2 - 1
                idx = -1*(idx + 1);                
            }
            
            int mappedValue = 254 - startBins.length + idx + 2;
            
            output.setValue(p.getX(), p.getY(), mappedValue);
        }
        
        return output;
    }
}
