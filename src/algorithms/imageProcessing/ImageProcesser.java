package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.clustering.KMeansPlusPlus;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author nichole
 */
public class ImageProcesser {
    
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
     * @param originalImage
     * @param outputSkyCentroid container to hold the output centroid of 
     * the sky.
     * @param edgeSettings
     * @return
     * @throws IOException
     * @throws NoSuchAlgorithmException 
     */
    public GreyscaleImage createSkyline(GreyscaleImage theta, 
        Image originalImage,
        CannyEdgeFilterSettings edgeSettings, PairIntArray outputSkyCentroid) 
        throws IOException, NoSuchAlgorithmException {        
      
        GreyscaleImage mask = createBestSkyMask(theta, originalImage, 
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
     * @param originalColorImage
     * @param edgeSettings
     * @param outputSkyCentroid container to hold the output centroid of 
     * the sky.
     * @return 
     * @throws java.io.IOException 
     * @throws java.security.NoSuchAlgorithmException 
     */
    public GreyscaleImage createBestSkyMask(final GreyscaleImage theta, 
        Image originalColorImage, CannyEdgeFilterSettings edgeSettings,
        PairIntArray outputSkyCentroid) throws 
        IOException, NoSuchAlgorithmException {
        
        if (theta == null) {
            throw new IllegalArgumentException("theta cannot be null");
        }
        
        GreyscaleImage thetaImg = theta;
        
        int binFactor = determineBinFactorForSkyMask(theta.getNPixels());
        
        if (binFactor > 1) {

            thetaImg = binImage(theta, binFactor);
        }

        List<PairIntArray> zeroPointLists = getSortedContiguousZeros(thetaImg);
        
        if (zeroPointLists.isEmpty()) {
            
            GreyscaleImage mask = thetaImg.createWithDimensions();
               
            // return an image of all 1's
            mask.fill(1);
            
            return mask;
        }
        
        if (binFactor > 1) {
                        
            thetaImg = theta;
        
            zeroPointLists = unbinZeroPointLists(zeroPointLists, binFactor);
            
        }
        
        //TODO:  this needs to be improved as soon as the rest of the
        //       algorithm is finished.
        //  currently not a robust way to determine that sky is horizontal
        //  or vertical and might not be correct for an angle not 90 degrees.
        //  probably needs equiv tranform for deconvolution of the image...
        
        int xMin = MiscMath.findMin(zeroPointLists.get(0).getX());
        int xMax = MiscMath.findMax(zeroPointLists.get(0).getX());
        int yMin = MiscMath.findMin(zeroPointLists.get(0).getY());
        int yMax = MiscMath.findMax(zeroPointLists.get(0).getY());
        
        double xLen = (double)(xMax - xMin)/(double)theta.getWidth();
        
        double yLen = (double)(yMax - yMin)/(double)theta.getHeight();
        
        boolean makeCorrectionsAlongX = (xLen < yLen) ? true : false;
        
        int convDispl = 6;
        
        // now the coordinates in zeroPointLists are w.r.t. thetaImg

        removeSetsThatAreDark(zeroPointLists, originalColorImage, thetaImg,
            makeCorrectionsAlongX, convDispl);
        
        reduceToLargest(zeroPointLists);
        
        removeHighContrastPoints(zeroPointLists, originalColorImage, 
            thetaImg, 
            makeCorrectionsAlongX, convDispl);
        
        if (binFactor > 1) {
        
Image img1 = theta.createWithDimensions().copyImageToGreen();
ImageIOHelper.addAlternatingColorCurvesToImage(zeroPointLists, img1);
ImageDisplayer.displayImage("before oversampling corrections", img1);

            int topToCorrect = zeroPointLists.size();

            //make corrections for resolution:
            addBackMissingZeros(zeroPointLists, theta, binFactor, topToCorrect);
/*
img1 = theta.createWithDimensions().copyImageToGreen();
ImageIOHelper.addAlternatingColorCurvesToImage(zeroPointLists, img1);
ImageDisplayer.displayImage("added missing zeros", img1);
*/
        }
           
        int[] rgb = getAvgMinMaxColor(zeroPointLists.get(0), thetaImg, 
            originalColorImage, makeCorrectionsAlongX, convDispl);

        Set<PairInt> points = combine(zeroPointLists);

        GreyscaleImage mask = growPointsToSkyline(points, originalColorImage, 
            theta, rgb, makeCorrectionsAlongX, convDispl);
       
log.info("SKY avg: " + rgb[0] + " min=" + rgb[1] + " max=" + rgb[2]);
ImageProcesser imageProcesser = new ImageProcesser();
imageProcesser.printImageColorContrastStats(originalColorImage, 161, 501);

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        double[] xycen = curveHelper.calculateXYCentroids(points);
 
        outputSkyCentroid.add((int)Math.round(xycen[0]), (int)Math.round(xycen[1]));

        if (mask != null) {
            removeSpurs(mask);
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
                Logger.getLogger(ImageProcesser.class.getName()).log(Level.SEVERE, null, ex);
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
    
    private int[] getAvgMinMaxColor(PairIntArray points, GreyscaleImage theta, 
        Image originalImage, boolean addAlongX, int addAmount) {
        
        int xOffset = theta.getXRelativeOffset();
        int yOffset = theta.getYRelativeOffset();
        
        double rSum = 0;
        double gSum = 0;
        double bSum = 0;
        
        double rgbMinSum = Double.MAX_VALUE;        
        double rgbMaxSum = Double.MIN_VALUE;
        
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
            
            if (rgb < rgbMinSum) {
                rgbMinSum = rgb;
            }
            if (rgb > rgbMaxSum) {
                rgbMaxSum = rgb;
            }
            
            count++;
        }
                
        if (count == 0) {
            return new int[]{0, 0, 0};
        }
        
        rSum /= (double)count;
        gSum /= (double)count;
        bSum /= (double)count;
        
        rgbMinSum /= (double)count;
        rgbMaxSum /= (double)count;
        
        double rgbSum = (rSum + gSum + bSum)/3.;        
        
        return new int[]{(int)Math.round(rgbSum), (int)Math.round(rgbMinSum),
            (int)Math.round(rgbMaxSum)};
    }

    /**
     * look for gaps in zeroValuePoints and if the color of the points within
     * the originalImage appears to be whiter than sky, consider the points
     * to be clouds.  cloud points are then added to zeroValuePoints because 
     * it's known that zeroValuePoints is subsequently used to mask out pixels
     * that should not be used for a skyline edge (and hence corners).
     * 
     * NOTE: this makes a correction for gaussian blurring, that is subtracting
     * the 6 pixels from convolving w/ gaussian of sigma=2, then sigma=0.5
     * which is what happens when "outdoorMode" is used to create the theta
     * image.
     * 
     * Note, this method won't find the clouds which are touching the horizon.
     * The Brown & Lowe 2003 panoramic images of a snowy mountain shows 
     * 2 examples of complications with cloud removal.
     * To the far left on the mountain horizon, one can see that clouds do
     * obscure the mountain, and in the middle of the horizon of image, 
     * one can see that clouds and the mountain peak are surrounded by 
     * sky to left and right, though not completely below.
     * In those cases where the clouds are touching the horizon, one probably 
     * wants to wait until have registered more than one image to understand 
     * motion and hence how to identify the clouds.
     * 
     * @param points
     * @param theta
     * @param originalImage
     * @return 
     */
    private void removeClouds(PairIntArray zeroValuePoints, GreyscaleImage theta, 
        Image originalImage, boolean addAlongX, int addAmount, 
        int rgbSkyAvg, int rgbSkyMin, int rgbSkyMax) {
        
        /*
        find the gaps in the set zeroValuePoints and determine if their
        color is whiteish compared to the background sky or looks like the
        background sky too.
        
        easiest way, though maybe not fastest:
           -- determine min and max of x and y of zeroValuePoints.
           -- scan along each row to find the start column of the row and the
              end column of the row.
           -- repeat the scan of the row only within column boundaries and 
              search for each point in zeroValuePoints
              -- if it is not present in zeroValuePoints,
                 check the originalImage value.  if it is white with respect
                 to rgbSky, add it to zeroValuePoints
           Can see that snowy mountain tops may be removed, so need to add to
           the block within the row scan, a scan above and below the possible
           cloud pixel to make sure that it is enclosed by sky pixels above
           and below too.
        
        would like to make a fast search of zeroValuePoints so using pixel index
        as main data and putting all into a hash set
        */
        
        if (zeroValuePoints.getN() == 0) {
            return;
        }
        
        Set<Integer> zpSet = new HashSet<Integer>();
        for (int pIdx = 0; pIdx < zeroValuePoints.getN(); pIdx++) {
            int x = zeroValuePoints.getX(pIdx);
            int y = zeroValuePoints.getY(pIdx);
            int idx = theta.getIndex(x, y);
            zpSet.add(Integer.valueOf(idx));
        }
        
        int xOffset = theta.getXRelativeOffset();
        int yOffset = theta.getYRelativeOffset();
          
        int rSky = (rgbSkyAvg >> 16) & 0xFF;
        int gSky = (rgbSkyAvg >> 8) & 0xFF;
        int bSky = rgbSkyAvg & 0xFF;
        
        int rMinSky = (rgbSkyMin >> 16) & 0xFF;
        int gMinSky = (rgbSkyMin >> 8) & 0xFF;
        int bMinSky = rgbSkyMin & 0xFF;
        
        int rMaxSky = (rgbSkyMax >> 16) & 0xFF;
        int gMaxSky = (rgbSkyMax >> 8) & 0xFF;
        int bMaxSky = rgbSkyMax & 0xFF;
        
        /*
        blue sky:  143, 243, 253
        white clouds on blue sky:  192, 242, 248  (increased red, roughly same blue)
        
        red sky:
        white clouds on red sky: (increased blue and green, roughly same red)
        */
        
        boolean skyIsBlue = (bSky > rSky);
        
        int xMin = MiscMath.findMin(zeroValuePoints.getX());
        int xMax = MiscMath.findMax(zeroValuePoints.getX());
        int yMin = MiscMath.findMin(zeroValuePoints.getY());
        int yMax = MiscMath.findMax(zeroValuePoints.getY());
        
        for (int row = yMin; row <= yMax; row++) {
            int start = -1;
            int stop = 0;
            for (int col = xMin; col <= xMax; col++) {
                int idx = theta.getIndex(col, row);
                if (zpSet.contains(Integer.valueOf(idx))) {
                    stop = col;
                    if (start == -1) {
                        start = col;
                    }
                }
            }
            
            if (start == -1) {
                continue;
            }
            
            // any pixels not in set are potentially cloud pixels
            for (int col = start; col <= stop; col++) {
                int idx = theta.getIndex(col, row);
                if (!zpSet.contains(Integer.valueOf(idx))) {
                    int x = theta.getCol(idx);
                    int y = theta.getRow(idx);
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
                    
                    int red = originalImage.getR(ox, oy);
                    int green = originalImage.getG(ox, oy);
                    int blue = originalImage.getB(ox, oy);

                    // is it white or near the color?
                    boolean looksLikeACloudPixel = false;
                    boolean looksLikeASkyPixel = false;
                    
                    looksLikeACloudPixel = ((blue >= bSky) && (red >= rSky)
                        && (green >= gSky) &&
                        ((skyIsBlue && (blue > bSky) && (green > gSky))
                        ||
                        (!skyIsBlue && (red > rSky))));
                    
                    if (skyIsBlue && !looksLikeACloudPixel) {
                        // expect blue ~ b, but red > r
                        int db = Math.abs(blue - bSky);
                        int dr = red - rSky;
                        // if dr is < 0, it's 'bluer', hence similar to sky
                        if ((db < 10) && (dr > 25)) {
                            // could be dark part of cloud, but could also be
                            // a rock formation for example
                            /*if ((((double)blue/(double)green) < 0.1) &&
                                (((double)blue/(double)red) > 2.)) {
                                looksLikeACloudPixel = true;
                            }*/
                            
                        } else if ((dr < 0) && (db >= 0)) {
                            if ((green - gSky) > 0) {
                                looksLikeASkyPixel = true;
                            }
                        }
                    } else if (!looksLikeACloudPixel) {
                        // expect red ~ r, but blue > b and green > g
                        int dr = Math.abs(red - rSky);
                        int db = blue - bSky;
                        int dg = green - gSky;
                        if ((dr < 10) && (db > 25) && (dg > 25)) {
                            looksLikeACloudPixel = true;
                        } else if ((db < 0) && (dr >= 0)) {
                            if ((green - gSky) < 0) {
                                looksLikeASkyPixel = true;
                            }
                        }
                    }
                      
                    if (looksLikeASkyPixel) {
                        zeroValuePoints.add(col, row);
                        zpSet.add(Integer.valueOf(idx));
                        if (col < xMin) {
                            xMin = col;
                        } 
                        if (col > xMax) {
                            xMax = col;
                        }
                        if (row < yMin) {
                            yMin = row;
                        } 
                        if (row > yMax) {
                            yMax = row;
                        }
                        continue;
                    }
                    
                    if (looksLikeACloudPixel) {
                        
                        //further scan up and down to make sure enclosed by sky
                        
                        boolean foundSkyPixel = false;
                        
                        // search for sky pixel above current (that is, lower y)
                        if ((row - 1) == -1) {
                            foundSkyPixel = true;
                        }
                        if (!foundSkyPixel) {
                            for (int rowI = (row - 1); rowI >= yMin; rowI--) {
                                int idxI = theta.getIndex(col, rowI);
                                if (zpSet.contains(Integer.valueOf(idxI))) {
                                    foundSkyPixel = true;
                                    break;
                                }
                            }
                        }
                        
                        if (!foundSkyPixel) {
                            continue;
                        }
                        
                        // search below for sky pixels, that is search higher y
                        
                        foundSkyPixel = false;
                        
                        if ((row + 1) == theta.getHeight()) {
                            foundSkyPixel = true;
                        }
                        
                        if (!foundSkyPixel) {
                            for (int rowI = (row + 1); rowI <= yMax; rowI++) {
                                int idxI = theta.getIndex(col, rowI);
                                if (zpSet.contains(Integer.valueOf(idxI))) {
                                    foundSkyPixel = true;
                                    break;
                                }
                            }
                        }
                        
                        if (!foundSkyPixel) {
                            // might be the peak of a snowy mountain
                            continue;
                        }
                        
                        // this looks like a cloud pixel, so region should be
                        // considered 'sky'
                        
                        zeroValuePoints.add(col, row);
                    }
                }
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

    /**
     * iterate over each point in zeroPointLists and visit its 8 neighbors
     * looking for those not in it's list.  if not in list and is in the
     * image as a zero value pixel, place it in the list.  note that this is a method to use
     * for correcting the zero points lists after down sampling to make the
     * list and then up sampling to use it.
     * 
     * @param zeroPointLists
     * @param theta 
     */
    private void addBackMissingZeros(List<PairIntArray> zeroPointLists, 
        GreyscaleImage theta, int binFactor, int topNumberToCorrect) {
                
        int width = theta.getWidth();
        int height = theta.getHeight();
        
        int end = topNumberToCorrect;
        if (zeroPointLists.size() < end) {
            end = zeroPointLists.size();
        }
        
        List<Set<Integer> > sets = new ArrayList<Set<Integer> >();
        for (int ii = 0; ii < end; ii++) {
            
            PairIntArray zeroValuePoints = zeroPointLists.get(ii);
            
            Set<Integer> zpSet = new HashSet<Integer>();
            for (int pIdx = 0; pIdx < zeroValuePoints.getN(); pIdx++) {
                int x = zeroValuePoints.getX(pIdx);
                int y = zeroValuePoints.getY(pIdx);
                int idx = theta.getIndex(x, y);
                zpSet.add(Integer.valueOf(idx));
            }
            
            sets.add(zpSet);
        }
        
        for (int ii = 0; ii < end; ii++) {

            PairIntArray zeroValuePoints = zeroPointLists.get(ii);

            Set<Integer> zpSet = sets.get(ii);

            Set<Integer> add = new HashSet<Integer>();
            
            for (int pIdx = 0; pIdx < zeroValuePoints.getN(); pIdx++) {

                int x = zeroValuePoints.getX(pIdx);
                int y = zeroValuePoints.getY(pIdx);

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

                        int neighborIdx = theta.getIndex(c, r);

                        Integer index = Integer.valueOf(neighborIdx);

                        if (zpSet.contains(index)) {
                            continue;
                        }

                        int v = theta.getValue(c, r);
                        if (v == 0) {
                            add.add(index);
                        }
                    }
                }
            }
            if (!add.isEmpty()) {
                for (Integer a : add) {
                    int c = theta.getCol(a.intValue());
                    int r = theta.getRow(a.intValue());
                    zeroValuePoints.add(c, r);
                    zpSet.add(a);
                }
            }
        }
       
    }

    /**
     * remove points in zeropoint lists where there is a non-zero pixel in the 
     * theta image.
     * 
     * @param zeroPointLists
     * @param theta 
     */
    private void removeOverSampledZeros(List<PairIntArray> zeroPointLists, 
        GreyscaleImage theta) {
        
        for (PairIntArray zeroValuePoints : zeroPointLists) {
            
            List<Integer> remove = new ArrayList<Integer>();
            
            for (int pIdx = 0; pIdx < zeroValuePoints.getN(); pIdx++) {
                
                int x = zeroValuePoints.getX(pIdx);
                int y = zeroValuePoints.getY(pIdx);
                
                int v = theta.getValue(x, y);
                
                if (v > 0) {
                    remove.add(Integer.valueOf(pIdx));
                }
            }
            
            if (!remove.isEmpty()) {
                for (int i = (remove.size() - 1); i > -1; i--) {
                    int idx = remove.get(i).intValue();
                    zeroValuePoints.removeRange(idx, idx);
                }
            }
        }
    }
    
    public void printImageColorContrastStats(Image image, List<PairIntArray> 
        zeroPoints) {
        
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
    
    private void removeSetsThatAreDark(List<PairIntArray> 
        zeroPointLists, Image originalColorImage, GreyscaleImage theta,
        boolean addAlongX, int addAmount) {
        
        int colorLimit = 100;
        
        List<Integer> remove = new ArrayList<Integer>();
        
        int xOffset = theta.getXRelativeOffset();
        int yOffset = theta.getYRelativeOffset();
        
        for (int gId = 0; gId < zeroPointLists.size(); gId++) {
            
            PairIntArray points  = zeroPointLists.get(gId);
            
            int nBelowLimit = 0;
            
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
                                
                if ((r < colorLimit) && (b < colorLimit) && (g < colorLimit)) {
                    nBelowLimit++;
                }                
            }
                        
            log.fine(gId + ") nBelowLimit=" + nBelowLimit
                + " (" + ((double)nBelowLimit/(double)points.getN()) + ")");
            
            if (((double)nBelowLimit/(double)points.getN()) > 0.5) {
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
    
    private void removeHighContrastPoints(List<PairIntArray> 
        zeroPointLists, Image originalColorImage, GreyscaleImage theta,
        boolean addAlongX, int addAmount) {
        
        int xOffset = theta.getXRelativeOffset();
        int yOffset = theta.getYRelativeOffset();
        
        double avgY = calculateY(zeroPointLists.get(0), originalColorImage,
             xOffset, yOffset, addAlongX, addAmount);
        
        // remove points that have contrast larger than tail of histogram
        HistogramHolder h = createContrastHistogram(avgY, zeroPointLists, 
            originalColorImage, xOffset, yOffset, addAlongX, 
            addAmount);
        
        if (h == null) {
            return;
        }
/*        
try {
    h.plotHistogram("contrast", 1);
} catch (IOException e) {
    log.severe(e.getMessage());
}
*/        
        int yPeakIdx = MiscMath.findYMaxIndex(h.getYHist());
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
    }

    public double calculateY(PairIntArray points, Image originalColorImage,
        int xOffset, int yOffset, boolean addAlongX, int addAmount) {
        
        double[][] m = new double[3][];
        m[0] = new double[]{0.256, 0.504, 0.098};
        m[1] = new double[]{-0.148, -0.291, 0.439};
        m[2] = new double[]{0.439, -0.368, -0.072};
        
        double avgY = 0;
        
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
        }
        avgY /= (double)points.getN();
        
        return avgY;
    }
    
    private HistogramHolder createContrastHistogram(List<PairIntArray> 
        zeroPointLists, Image originalColorImage, int xRelativeOffset,
        int yRelativeOffset, boolean addAlongX, int addAmount) {
                
        if (zeroPointLists.isEmpty()) {
            return null;
        }
        
        double avgY = calculateY(zeroPointLists.get(0), originalColorImage,
            xRelativeOffset, yRelativeOffset, addAlongX, addAmount);
        
        return createContrastHistogram(avgY, zeroPointLists, 
            originalColorImage, xRelativeOffset, yRelativeOffset, 
            addAlongX, addAmount);
    }
    
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
                
                int r = originalColorImage.getR(x, y);
                int g = originalColorImage.getG(x, y);
                int b = originalColorImage.getB(x, y);
                
                double[] rgb = new double[]{r, g, b};
                        
                double[] yuv = MatrixUtil.multiply(m, rgb);
                yuv = MatrixUtil.add(yuv, new double[]{16, 128, 128});

                float contrastValue = (float)((avgY - yuv[0]) / yuv[0]);
                
                yValues[count] = contrastValue;
                
                count++;
            }
        }
        
        HistogramHolder h = Histogram.calculateSturgesHistogramRemoveZeroTail(
            yValues, Errors.populateYErrorsBySqrt(yValues));

        return h;
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

    private GreyscaleImage growPointsToSkyline(Set<PairInt> points, 
        Image originalColorImage, GreyscaleImage theta, int[] rgb,
        boolean makeCorrectionsAlongX, int addAmount) {
        
        transformPointsToOriginalReferenceFrame(points, theta, 
            makeCorrectionsAlongX, addAmount);
        
        // dfs w/ boundary found by contrast and color
        
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
