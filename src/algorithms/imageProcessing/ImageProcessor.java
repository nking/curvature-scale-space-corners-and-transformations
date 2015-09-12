package algorithms.imageProcessing;

import algorithms.compGeometry.clustering.KMeansPlusPlus;
import algorithms.compGeometry.clustering.KMeansPlusPlusFloat;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.PairInt;
import algorithms.misc.Complex;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscDebug;
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

    /**
     * <pre>
     * [1, 0, -1]
     * this is the n=2 binomial filter for a Gaussian first derivative,
       that is sigma = sqrt(2)/2 = 0.707
       </pre>
     * @param input 
     */
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

        GreyscaleImage output = convolvedX.createWithDimensions();

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

        int w2 = offsetsAndDimensions[2];
        int h2 = offsetsAndDimensions[3];

        int offset1X = offsetsAndDimensions[0];
        int offset1Y = offsetsAndDimensions[1];

        GreyscaleImage output = new GreyscaleImage(w2, h2);
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

    /**
     * applies KMeansPlusPlus algorithm to the values in input to create
     * kBands segmented image (operates on input).
     * @param input
     * @param kBands
     * @throws IOException
     * @throws NoSuchAlgorithmException 
     */
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
     * multiply these images, that is pixel by pixel multiplication.
     * input2 is assumed to be 0 or 1
     * @param input1
     * @param input2 the mask of 0's and 1's to apply to input1
     * @return
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

    public void blur(GreyscaleImage input, float sigma) {

        float[] kernel = Gaussian1D.getKernel(sigma);

        blur(input, kernel);
    }

    public void blur(GreyscaleImage input, SIGMA sigma) {

        float[] kernel = Gaussian1D.getKernel(sigma);

        blur(input, kernel);
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

        GreyscaleImage out = new GreyscaleImage(w1, h1);
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

        GreyscaleImage out = new GreyscaleImage(w1, h1);
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
    
    public GreyscaleImage padUpToPowerOfTwo(GreyscaleImage input) {
        
        int w0 = input.getWidth();
        int h0 = input.getHeight();
        
        int w = 1 << (int)(Math.ceil(Math.log(w0)/Math.log(2)));
        int h = 1 << (int)(Math.ceil(Math.log(h0)/Math.log(2)));

        int xOffset = w - w0;
        int yOffset = h - h0;
        
        if (xOffset == 0 && yOffset == 0) {
            return input;
        }
        
        int xOffsetOrig = input.getXRelativeOffset();
        int yOffsetOrig = input.getYRelativeOffset();
        
        GreyscaleImage output = new GreyscaleImage(w, h);
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

    /**
     *
     * @param input
     * @param forward if true, apply FFT transform, else inverse FFT transform
     */
    public void apply2DFFT(GreyscaleImage input, boolean forward) {

        int xOffsetOrig = input.getXRelativeOffset();
        int yOffsetOrig = input.getYRelativeOffset();
            
        //TODO remove the other power of 2 padding method
        GreyscaleImage tmp = padUpToPowerOfTwo(input);

        // initialize matrix of complex numbers as real numbers from image
        Complex[][] cc = convertImage(tmp);

        Complex[][] ccOut = apply2DFFT(cc, forward);

        writeToImage(tmp, ccOut);
        
        if (tmp.getNPixels() > input.getNPixels()) {
                        
            int xOffset = tmp.getXRelativeOffset();
            int yOffset = tmp.getYRelativeOffset();
            
            // padding is at front of cols and rows
            int x = 0;
            for (int col = xOffset; col < tmp.getWidth(); col++) {
                int y = 0;
                for (int row = yOffset; row < tmp.getHeight(); row++) {
                    int v = tmp.getValue(col, row);
                    input.setValue(x, y, v);
                    y++;
                }
                x++;
             }
            input.setXRelativeOffset(xOffsetOrig);
            input.setYRelativeOffset(yOffsetOrig);
        }
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
     * create an image segmented by color...useful for experimenting with corners
     * due to color differences rather than original intensity differences.
     * The default provided here uses 8 color segmentation.
     * @param input
     * @return
     */
    public GreyscaleImage createGreyscaleFromColorSegmentation(ImageExt input) {

        int kColors = 8;

        return createGreyscaleFromColorSegmentation(input, kColors);
    }
    
    /**
     * create an image segmented by color...useful for experimenting with corners
     * due to color differences rather than original intensity differences.
     * The default provided here uses 8 color segmentation.
     * @param input
     * @return
     */
    public GreyscaleImage createGreyscaleFromColorSegmentationKMPP(ImageExt 
        input) {

        int kColors = 8;

        return createGreyscaleFromColorSegmentationKMPP(input, kColors);
    }

    /**
     * create an image segmented by color...useful for experimenting with corners
     * due to color differences rather than original intensity differences.
     * Internally, is using a blur of '1' on the image before segmentation.
     * @param input
     * @param kColors the number of colors to bin the image by.  max allowed value
     * is 253.
     * @return
     */
    public GreyscaleImage createGreyscaleFromColorSegmentation(ImageExt input,
        int kColors) {

        boolean useBlur = true;

        return createGreyscaleFromColorSegmentation(input, kColors, useBlur);
    }
    
    /**
     * create an image segmented by color...useful for experimenting with corners
     * due to color differences rather than original intensity differences.
     * Internally, is using a blur of '1' on the image before segmentation.
     * @param input
     * @param kColors the number of colors to bin the image by.  max allowed value
     * is 253.
     * @return
     */
    public GreyscaleImage createGreyscaleFromColorSegmentationKMPP(ImageExt input,
        int kColors) {

        boolean useBlur = true;

        return createGreyscaleFromColorSegmentationKMPP(input, kColors, useBlur);
    }

    /**
     * create an image segmented by color...useful for experimenting with corners
     * due to color differences rather than original intensity differences.
     * Internally, is using a blur of '1' on the image before segmentation.
     * @param input
     * @param kColors the number of colors to bin the image by.  max allowed value
     * is 253.
     * @param useBlur
     * @return
     */
    public GreyscaleImage createGreyscaleFromColorSegmentation(ImageExt input,
        int kColors, boolean useBlur) {

        if (kColors > 253) {
            throw new IllegalArgumentException("kColors must be <= 253");
        }
        if (kColors < 2) {
            throw new IllegalArgumentException("kColors must be >= 2");
        }

        GreyscaleImage img;

        int minNeighborLimit;

        if (useBlur) {

            Image input2 = input.copyImage();

            blur(input2, 1/*(float)Math.sqrt(2)/2.f*/);

            img = convertToCIEXYPolarTheta(input2, kColors);

            minNeighborLimit = 6;

        } else {

            img = convertToCIEXYPolarTheta(input, kColors);

            minNeighborLimit = 5;
        }

        int w = img.getWidth();
        int h = img.getHeight();

        // ----replace pixel, if 7 or more neighbors have same color -----
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};

        Map<Integer, Integer> freqMap = new HashMap<Integer, Integer>();

        int nChanged = 1;
        int nIterMax = 100;
        int nIter = 0;

        while (!useBlur && (nIter < nIterMax) && (nChanged > 0)) {

            log.fine("nIter=" + nIter + " nChanged=" + nChanged);

            nChanged = 0;

            for (int col = 0; col < w; col++) {
                for (int row = 0; row < h; row++) {

                    freqMap.clear();

                    Integer maxCountValue = null;
                    int maxCount = Integer.MIN_VALUE;

                    for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                        int x = dxs[nIdx] + col;
                        int y = dys[nIdx] + row;

                        if ((x < 0) || (x > (w - 1)) || (y < 0) || (y > (h - 1))) {
                            break;
                        }

                        Integer v = Integer.valueOf(img.getValue(x, y));

                        Integer c = freqMap.get(v);
                        if (c == null) {
                            c = Integer.valueOf(1);
                        } else {
                            c = Integer.valueOf(c.intValue() + 1);
                        }
                        freqMap.put(v, c);

                        if (c.intValue() > maxCount) {
                            maxCount = c.intValue();
                            maxCountValue = v;
                        }
                    }

                    if ((maxCount >= minNeighborLimit) &&
                        (img.getValue(col, row) != maxCountValue.intValue())) {

                        img.setValue(col, row, maxCountValue.intValue());
                        nChanged++;
                    }
                }
            }

            nIter++;
        }

        // rescale the image
        HistogramEqualization hEq = new HistogramEqualization(img);
        hEq.applyFilter();

        return img;
    }

    /**
     * create an image segmented by color...useful for experimenting with corners
     * due to color differences rather than original intensity differences.
     * Internally, is using a blur of '1' on the image before segmentation.
     * @param input
     * @param kColors the number of colors to bin the image by.  max allowed value
     * is 253.
     * @param useBlur
     * @return
     */
    public GreyscaleImage createGreyscaleFromColorSegmentationKMPP(ImageExt input,
        int kColors, boolean useBlur) {

        if (kColors > 253) {
            throw new IllegalArgumentException("kColors must be <= 253");
        }
        if (kColors < 2) {
            throw new IllegalArgumentException("kColors must be >= 2");
        }

        GreyscaleImage img;

        int minNeighborLimit;

        if (useBlur) {

            Image input2 = input.copyImage();

            blur(input2, 1/*(float)Math.sqrt(2)/2.f*/);

            img = convertToCIEXYPolarThetaKMPP(input2, kColors);

            minNeighborLimit = 6;

        } else {

            img = convertToCIEXYPolarThetaKMPP(input, kColors);

            minNeighborLimit = 5;
        }

        int w = img.getWidth();
        int h = img.getHeight();

        // ----replace pixel, if 7 or more neighbors have same color -----
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};

        Map<Integer, Integer> freqMap = new HashMap<Integer, Integer>();

        int nChanged = 1;
        int nIterMax = 100;
        int nIter = 0;

        while (!useBlur && (nIter < nIterMax) && (nChanged > 0)) {

            log.fine("nIter=" + nIter + " nChanged=" + nChanged);

            nChanged = 0;

            for (int col = 0; col < w; col++) {
                for (int row = 0; row < h; row++) {

                    freqMap.clear();

                    Integer maxCountValue = null;
                    int maxCount = Integer.MIN_VALUE;

                    for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                        int x = dxs[nIdx] + col;
                        int y = dys[nIdx] + row;

                        if ((x < 0) || (x > (w - 1)) || (y < 0) || (y > (h - 1))) {
                            break;
                        }

                        Integer v = Integer.valueOf(img.getValue(x, y));

                        Integer c = freqMap.get(v);
                        if (c == null) {
                            c = Integer.valueOf(1);
                        } else {
                            c = Integer.valueOf(c.intValue() + 1);
                        }
                        freqMap.put(v, c);

                        if (c.intValue() > maxCount) {
                            maxCount = c.intValue();
                            maxCountValue = v;
                        }
                    }

                    if ((maxCount >= minNeighborLimit) &&
                        (img.getValue(col, row) != maxCountValue.intValue())) {

                        img.setValue(col, row, maxCountValue.intValue());
                        nChanged++;
                    }
                }
            }

            nIter++;
        }

        // rescale the image
        HistogramEqualization hEq = new HistogramEqualization(img);
        hEq.applyFilter();

        return img;
    }
    
    /**
     * convert the color image into an image scaled into values 0 to 255
     * by the polar theta angle in CIE XY color space.  The colors are
     * divided into 254 bins plus black and white.
     *
     * runtime complexity is O(N) + O(N*lg_2(N))
     *
     * @param input
     * @return
     */
    public GreyscaleImage convertToCIEXYPolarTheta(Image input) {

        return convertToCIEXYPolarTheta(input, 254);
    }

    /**
     * convert the color image into an image scaled into values 0 to 255
     * by the polar theta angle in CIE XY color space.  The colors are divided
     * into kColors number of bins.
     *
     * runtime complexity is O(N) + O(N*lg_2(N))
     *
     * @param input
     * @param kColors the number of color bins to use for the image segmentation.
     * The minimum allowed value is 2 and the maximum allowed value is 253.
     * @return
     */
    public GreyscaleImage convertToCIEXYPolarTheta(Image input, int kColors) {

        if (kColors > 253) {
            throw new IllegalArgumentException("kColors must be <= 253");
        }
        if (kColors < 2) {
            throw new IllegalArgumentException("kColors must be >= 2");
        }

        /*
        TODO: needs some improvements in color mapping.
           -- for regions that are very spotty, might consider using the
              intensity image to help define the region and take the highest
              number density color in that region and assign it to all.
        */

        int w = input.getWidth();
        int h = input.getHeight();

        float[] tmpColorBuffer = new float[2];

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

                if ((r < 25) && (g < 25) && (b < 25)) {
                    continue;
                } else if ((r > 230) && (g > 230) && (b > 230)) {//might need to use 195 as lower limit
                    output.setValue(col, row, 255);
                    continue;
                }

                float[] cieXY = tmpColorBuffer;
                cieC.rgbToXYChromaticity(r, g, b, cieXY);

                if (cieC.isWhite(cieXY[0], cieXY[1])) {
                    output.setValue(col, row, 255);
                } else {

                    double theta = cieC.calculateXYTheta(cieXY[0], cieXY[1]);

                    thetaValues[thetaCount] = (float)theta;

                    pixThetaMap.put(p, Float.valueOf((float)theta));

                    thetaCount++;
                }
            }
        }

        thetaValues = Arrays.copyOf(thetaValues, thetaCount);

        createAndApplyHistMapping(output, pixThetaMap, thetaValues, kColors);

        return output;
    }

    /**
     * convert the color image into an image scaled into values 0 to 255
     * by the polar theta angle in CIE XY color space.  The colors are divided
     * into kColors number of bins.
     *
     * runtime complexity is O(N) + O(N*lg_2(N))
     *
     * @param input
     * @param kColors the number of color bins to use for the image segmentation.
     * The minimum allowed value is 2 and the maximum allowed value is 253.
     * @return
     */
    public GreyscaleImage convertToCIEXYPolarThetaKMPP(Image input, int kColors) {

        if (kColors > 253) {
            throw new IllegalArgumentException("kColors must be <= 253");
        }
        if (kColors < 2) {
            throw new IllegalArgumentException("kColors must be >= 2");
        }

        /*
        TODO: needs some improvements in color mapping.
           -- for regions that are very spotty, might consider using the
              intensity image to help define the region and take the highest
              number density color in that region and assign it to all.
        */

        int w = input.getWidth();
        int h = input.getHeight();

        float[] tmpColorBuffer = new float[2];

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

                if ((r < 25) && (g < 25) && (b < 25)) {
                    continue;
                } else if ((r > 230) && (g > 230) && (b > 230)) {//might need to use 195 as lower limit
                    output.setValue(col, row, 255);
                    continue;
                }

                float[] cieXY = tmpColorBuffer;
                cieC.rgbToXYChromaticity(r, g, b, cieXY);

                if (cieC.isWhite(cieXY[0], cieXY[1])) {
                    output.setValue(col, row, 255);
                } else {

                    double theta = cieC.calculateXYTheta(cieXY[0], cieXY[1]);

                    thetaValues[thetaCount] = (float)theta;

                    pixThetaMap.put(p, Float.valueOf((float)theta));

                    thetaCount++;
                }
            }
        }

        thetaValues = Arrays.copyOf(thetaValues, thetaCount);

        createAndApplyKMPPMapping(output, pixThetaMap, thetaValues, kColors);

        return output;
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

    private void createAndApplyHistMapping(GreyscaleImage output,
        Map<PairInt, Float> pixThetaMap, float[] thetaValues,
        final int kColors) {

        float minValue = MiscMath.findMin(thetaValues);
        float maxValue = MiscMath.findMax(thetaValues);

        log.fine("minTheta=" + (minValue * 180./Math.PI) +
            " maxTheta=" + (maxValue * 180./Math.PI));

        int nReserved = 254 - kColors;

        HistogramHolder hist = Histogram.createSimpleHistogram(minValue,
            maxValue, (256 - nReserved - 1), thetaValues,
            Errors.populateYErrorsBySqrt(thetaValues));

        try {
            hist.plotHistogram("cie XY theta histogram", "cieXY_hist_" 
                + MiscDebug.getCurrentTimeFormatted());
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

            int mappedValue = 255 - startBins.length + idx;

            output.setValue(p.getX(), p.getY(), mappedValue);
        }
    }

    private void createAndApplyKMPPMapping(GreyscaleImage output,
        Map<PairInt, Float> pixThetaMap, float[] thetaValues,
        final int kColors) {

        //TODO: assert kColors.  The invoker is reserving 2 bands for 
        // B & W, so nBins should probably be (kColors - 2)...
        // correct this for the invoker when testing
        int nBins = kColors;
        
        KMeansPlusPlusFloat kmpp = new KMeansPlusPlusFloat();
        kmpp.computeMeans(nBins, thetaValues);
        
        float minValue = kmpp.getMinValue();
        float maxValue = kmpp.getMaxValue();

        float[] binCenters = kmpp.getCenters();
        
        Iterator<Entry<PairInt, Float> > iter = pixThetaMap.entrySet().iterator();

        while (iter.hasNext()) {

            Entry<PairInt, Float> entry = iter.next();

            PairInt p = entry.getKey();

            float theta = entry.getValue().floatValue();

            for (int i = 0; i < binCenters.length; i++) {

                float vc = binCenters[i];

                float bisectorBelow = ((i - 1) > -1) ?
                    ((binCenters[i - 1] + vc) / 2) : minValue;

                float bisectorAbove = ((i + 1) > (binCenters.length - 1)) ?
                    maxValue : ((binCenters[i + 1] + vc) / 2);

                if ((theta >= bisectorBelow) && (theta <= bisectorAbove)) {

                    //TODO: check this
                    int mappedValue = 255 - nBins + i;
                    
                    output.setValue(p.getX(), p.getY(), mappedValue);

                    break;
                }
            }
            
            /* 
            // if binCenters is ordered, use binary search for faster results
            int idx = Arrays.binarySearch(startBins, theta);

            // if it's negative, (-(insertion point) - 1)
            if (idx < 0) {
                // idx = -*idx2 - 1
                idx = -1*(idx + 1);
            }
            int mappedValue = 255 - startBins.length + idx;

            output.setValue(p.getX(), p.getY(), mappedValue);
            */
        }
    }
    public void applyAdaptiveMeanThresholding(GreyscaleImage img) {
        
        applyAdaptiveMeanThresholding(img, 3);
    }
    
    public void applyAdaptiveMeanThresholding(GreyscaleImage img, 
        int halfDimension) {
        
        GreyscaleImage imgM = img.copyImage();

        /*
        7 x 7 averaging
        */
        applyCenteredMean(imgM, halfDimension);

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
        
        GreyscaleImage mean = img.createWithDimensions();
        
        // sum along rows
        for (int i = 0; i < w; ++i) {
            int sum0 = 0;
            for (int j = 0; j < dimension; ++j) {
                sum0 += img.getValue(i, j);
            }
            mean.setValue(i, 0, sum0);
            for (int j = 1; j <= (h - dimension); ++j) {                
                int vp = img.getValue(i, j - 1);
                int vl =  img.getValue(i, dimension + j - 1);
                sum0 = sum0 - vp + vl;
                mean.setValue(i, j, sum0);
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
                mean.setValue(i, j, Math.round(sum));
            }
        }
        
        // sum along columns
        for (int j = 0; j < h; ++j) {
            int sum0 = 0;
            for (int i = 0; i < dimension; ++i) {
                sum0 += mean.getValue(i, j);
            }
            img.setValue(0, j, sum0);
            for (int i = 1; i <= (w - dimension); ++i) {
                int vp = mean.getValue(i - 1, j);
                int vl = mean.getValue(dimension + i - 1, j);
                sum0 = sum0 - vp + vl;
                img.setValue(i, j, sum0);
            }
            
            // last dimension - 1 cols: sum along them, divide by count then mult by dimension
            for (int i = (w - dimension + 1); i < w; ++i) {
                float count = h - i;
                float sum = 0;
                for (int k = i; k < w; ++k) {
                    sum += mean.getValue(k, j);
                }
                sum /= count;
                sum *= dimension;
                img.setValue(i, j, Math.round(sum));
            }            
        }

        // divide each value by dimension * dimension
        float dsq = dimension * dimension;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int v = img.getValue(i, j);
                v = Math.round((float)v/dsq);
                img.setValue(i, j, v);
            }
        } 
    }
    
    /**
     * create an image of the mean of the surrounding dimension x dimension 
     * pixels for each pixel centered on each pixel.  For the starting
     * and ending (dimension/2) pixels, the average uses a descreasing
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
     * runtime complexity is O(N_pixels)
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
        
        GreyscaleImage mean = img.createWithDimensions();
        
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
                    sum += img.getValue(i, k);
                }
                sum /= count;
                sum *= dimension;
                mean.setValue(i, j, Math.round(sum));
            }
            
            /* pixels between halfDimension and j-halfDimension
            halfDimension = 2            
            0 1 2 3 4 5
            |   *   |  sum from idx - halfDimension to idx + halfDimension, incl
            but store in idx
            */
            int sum0 = 0;
            for (int j = 0; j <= 2*halfDimension; ++j) {
                sum0 += img.getValue(i, j);
            }
            mean.setValue(i, halfDimension, sum0);
            /*
            halfDimension = 2            
            0 1 2 3 4 5 6
              |   *   | 
                |   *   |
                    
            */
            for (int j = halfDimension + 1; j < (h - halfDimension); ++j) {                
                int vp = img.getValue(i, j - halfDimension - 1);
                int vl =  img.getValue(i, j + halfDimension);
                sum0 = sum0 - vp + vl;
                mean.setValue(i, j, sum0);
            }
            /* last halfDimension pixels
            0 1 2 3   n=4, halfDimension = 2 
                >
            */
            for (int j = (h - halfDimension); j < h; ++j) {
                float count = h - j + 1;
                float sum = 0;
                for (int k = (j - 1); k < h; ++k) {
                    sum += img.getValue(i, k);
                }
                sum /= count;
                sum *= dimension;
                mean.setValue(i, j, Math.round(sum));
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
                    sum += mean.getValue(k, j);
                }
                sum /= count;
                sum *= dimension;
                img.setValue(i, j, Math.round(sum));
            }
            
            /* pixels between halfDimension and j-halfDimension
            halfDimension = 2            
            0 1 2 3 4 5
            |   *   |  sum from idx - halfDimension to idx + halfDimension, incl
            but store in idx
            */
            int sum0 = 0;
            for (int i = 0; i <= 2*halfDimension; ++i) {
                sum0 += mean.getValue(i, j);
            }
            img.setValue(halfDimension, j, sum0);
            /*
            halfDimension = 2            
            0 1 2 3 4 5 6
              |   *   | 
                |   *   |
                    
            */
            for (int i = halfDimension + 1; i < (w - halfDimension); ++i) {                
                int vp = mean.getValue(i - halfDimension - 1, j);
                int vl =  mean.getValue(i + halfDimension, j);
                sum0 = sum0 - vp + vl;
                img.setValue(i, j, sum0);
            }
            /* last halfDimension pixels
            0 1 2 3   n=4, halfDimension = 2 
                >
            */
            for (int i = (w - halfDimension); i < w; ++i) {
                float count = w - i + 1;
                float sum = 0;
                for (int k = (i - 1); k < w; ++k) {
                    sum += mean.getValue(k, j);
                }
                sum /= count;
                sum *= dimension;
                img.setValue(i, j, Math.round(sum));
            }
        }
        
        // divide each value by dimension * dimension
        float dsq = dimension * dimension;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int v = img.getValue(i, j);
                v = Math.round((float)v/dsq);
                img.setValue(i, j, v);
            }
        } 
    }
}
