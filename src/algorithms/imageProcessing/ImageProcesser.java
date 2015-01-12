package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.PointInPolygon;
import algorithms.compGeometry.clustering.KMeansPlusPlus;
import algorithms.compGeometry.convexHull.GrahamScanInt;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

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
          
        GreyscaleImage img2 = new GreyscaleImage(imageX.getWidth(), imageX.getHeight());
        
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
        
        GreyscaleImage output = new GreyscaleImage(input.getWidth(), input.getHeight());
                 
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

        GreyscaleImage output = new GreyscaleImage(convolvedX.getWidth(), 
            convolvedX.getHeight());
        
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
        
        GreyscaleImage output = new GreyscaleImage(image.getWidth(), image.getHeight());
        
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
     * @param outputSkyCentroid container to hold the output centroid of 
     * the sky.
     * @return
     * @throws IOException
     * @throws NoSuchAlgorithmException 
     */
    public GreyscaleImage createSkyline(GreyscaleImage theta, 
        PairIntArray outputSkyCentroid) throws IOException, 
        NoSuchAlgorithmException {        
      
        GreyscaleImage mask = createBestSkyMask(theta, outputSkyCentroid);
        
        if (mask != null) {
            
            multiply(mask, 255);
            
            CannyEdgeFilter filter = new CannyEdgeFilter();
            
            filter.applyFilter(mask);
            
            return mask;
        }
        
        return null;
    }
    
    /**
     * NOT READY FOR USE YET. needs correction for gaussian convolution radius,
     * (@see createSkyline comments)
     * @param input
     * @param theta
     * @throws IOException
     * @throws NoSuchAlgorithmException 
     */
    public void applyImageSegmentationForSky(GreyscaleImage input, 
        GreyscaleImage theta) throws IOException, NoSuchAlgorithmException {        
      
        PairIntArray outputSkyCentroid = new PairIntArray();
        
        GreyscaleImage mask = createBestSkyMask(theta, outputSkyCentroid);
        
        //GreyscaleImage mask = createRoughSkyMask(theta);
        
        if (mask != null) {
            
            GreyscaleImage out = multiply(input, mask);
            
            input.resetTo(out);
        }
    }
    
    /**
     * 
     * NOT READY FOR USE YET
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
    
    /**
     * create a mask for the largest contiguous zero value area in theta
     * or for the top similar largest contiguous zero value areas if there
     * are more than one.  the method doesn't make an assumption about the
     * location of 'sky'.
     * 
     * @param theta
     * @param outputSkyCentroid container to hold the output centroid of 
     * the sky.
     * @return 
     */
    public GreyscaleImage createBestSkyMask(final GreyscaleImage theta, 
        PairIntArray outputSkyCentroid) throws 
        IOException, NoSuchAlgorithmException {
        
        if (theta == null) {
            throw new IllegalArgumentException("theta cannot be null");
        }
        
        //theta = theta.copyImage();
        
        //applyImageSegmentation(theta, 2);
        
        //subtractMinimum(theta);
                
        List<PairIntArray> zeroPointLists = getLargestContiguousZeros(theta);
       
        if (zeroPointLists.isEmpty()) {
            
            GreyscaleImage mask = new GreyscaleImage(theta.getWidth(), 
                theta.getHeight());
            
            // return an image of all 1's
            mask.fill(1);
            
            return mask;
        }
        
        GreyscaleImage mask = null;
                
        for (PairIntArray points : zeroPointLists) {
            
            GreyscaleImage output = createNoiseRemovingMask(theta, points);
            
            if (mask == null) {
                mask = output;
            } else {
                mask = binaryOr(mask, output);
            }
        }
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        double[] xycen = curveHelper.calculateXYCentroids(zeroPointLists);
        
        outputSkyCentroid.add((int)Math.round(xycen[0]), (int)Math.round(xycen[1]));
        
        if (mask != null) {
            removeSpurs(mask);
        }
              
        return mask;
    }
    
    public List<PairIntArray> getLargestContiguousZeros(GreyscaleImage theta) {
                 
        return getLargestContiguousValue(theta, 0, false);
    }
    
    public List<PairIntArray> getLargestContiguousNonZeros(GreyscaleImage theta) {
                 
        return getLargestContiguousValue(theta, 0, true);
    }
    
    private List<PairIntArray> getLargestContiguousValue(GreyscaleImage theta,
        int value, boolean excludeValue) {
        
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
                
                float number = groupN[i];
                
                float frac = number/n0;
System.out.println(number);    
                //TODO: this should be adjusted by some metric.
                //      a histogram?
                if ((1 - frac) < 0.4) {
                //if (number > 100) {
                    groupIds.add(Integer.valueOf(groupIndexes[i]));
                } else {
                    break;
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
        
        GreyscaleImage output = new GreyscaleImage(input1.getWidth(), 
            input1.getHeight());
        
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
        
        GreyscaleImage output = new GreyscaleImage(input1.getWidth(), 
            input1.getHeight());
        
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
        
        GreyscaleImage output = new GreyscaleImage(input.getWidth(), 
            input.getHeight());
        
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
        
        int width = theta.getWidth();
        
        int height = theta.getHeight();
        
        GreyscaleImage out = new GreyscaleImage(width, height);
           
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
     * @param zeroCoords
     * @return 
     */
    public GreyscaleImage createNoiseRemovingMask(GreyscaleImage theta, 
        PairIntArray zeroCoords) {
        
        int width = theta.getWidth();
        
        int height = theta.getHeight();
        
        GreyscaleImage out = new GreyscaleImage(width, height);
           
        out.fill(1);
                
        float critFrac = 0.8f;
        
        int hw = 10;
                    
        for (int pIdx = 0; pIdx < zeroCoords.getN(); pIdx++) {

            int x = zeroCoords.getX(pIdx);
            int y = zeroCoords.getY(pIdx);
            out.setValue(x, y, 0);
            
            // if a direct neighbor is not 0, look at the fraction of zeros
            // surounding that neighbor and set it to 0 if larger than crit
            for (int col = (x - 1); col <= (x + 1); col++) {
                if ((col < 0) || (col > (width - 1))) {
                    continue;
                }
                for (int row = (y - 1); row <= (y + 1); row++) {
                    if ((row < 0) || (row > (height - 1))) {
                        continue;
                    }
                    int v = theta.getValue(col, row);
                    if ((v == 0) || (out.getValue(col, row) == 0)) {
                        continue;
                    }

                    // see if this is noise in the sky, hence sky
                    int zCount = 0;
                    int aCount = 0;
                    //see if surrounding area is mostly zeros to null the pixel
                    for (int c = (col - hw); c <= (col + hw); c++) {
                        if ((c < 0) || (c > (width - 1))) {
                            continue;
                        }
                        for (int r = (row - hw); r <= (row + hw); r++) {
                            if ((r < 0) || (r > (height - 1))) {
                                continue;
                            }
                            aCount++;
                            if (theta.getValue(c, r) == 0) {
                                zCount++;
                            }
                        }
                    }
                    if (aCount == 0) {
                        continue;
                    }
                    float fracZ = (float) zCount / (float) aCount;
                    if (fracZ > critFrac) {
                        out.setValue(col, row, 0);
                    }
                }
            }
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
        
        int width = theta.getWidth();
        
        int height = theta.getHeight();
        
        GreyscaleImage out = new GreyscaleImage(width, height);
                                                       
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
}
