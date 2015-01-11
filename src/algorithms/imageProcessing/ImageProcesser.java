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
     * apply kernel to input. NOTE, that because the image is composed of vectors
     * that should have values between 0 and 255, inclusive, if the kernel application
     * results in a value outside of that range, the value is reset to 0 or
     * 255.
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
    
    public void applyImageSegmentationForSky(GreyscaleImage input, 
        GreyscaleImage theta) throws IOException, NoSuchAlgorithmException {
                
        applyImageSegmentation(theta, 2);
        
        subtractMinimum(theta);
      
        GreyscaleImage mask = createSkyMask(theta);
        
        if (mask != null) {
            
            GreyscaleImage out = multiply(input, mask);
            
            input.resetTo(out);
        }
        
    }
    
    /**
     * create a mask for the largest contiguous zero value area in theta
     * or for the top similar largest contiguous zero value areas if there
     * are more than one.  the method doesn't make an assumption about the
     * location of 'sky'.
     * 
     * @param theta
     * @return 
     */
    public GreyscaleImage createSkyMask(GreyscaleImage theta) {
        
        if (theta == null) {
            throw new IllegalArgumentException("theta cannot be null");
        }
        
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
            
            GreyscaleImage output = createMask(theta, points);
            
            if (mask == null) {
                mask = output;
            } else {
                mask = binaryOr(mask, output);
            }
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
}
