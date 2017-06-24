package algorithms.imageProcessing;

import java.util.List;

/**
 *
 * @author nichole
 */
public class ATrousWaveletTransform {

    /**
    Uses the first step of an A Trous wavelet transform with B3Spline scaling.
    * @param input
     * @return 
    */
    public GreyscaleImage smoothToLevel1B3Spline(GreyscaleImage input) {
        
        B3SplineFunction scalingFunction = new B3SplineFunction();
        
        GreyscaleImage smoothed = scalingFunction.calculate(input);
        
        return smoothed;
    }
    
    /**
    Uses the first step of an A Trous wavelet transform using triangle scaling.
     * @param input
     * @return 
    */
    public GreyscaleImage smoothFirstLevelTriangle(GreyscaleImage input) {
        
        TriangleFunction scalingFunction = new TriangleFunction();
        
        GreyscaleImage smoothed = scalingFunction.calculateNextLevel(input, 0);
        
        return smoothed;
    }

    /**
    Uses the first step of an A Trous wavelet transform, which is two 1D 
    convolutions of a binomial kernel for sigma=0.707 (=sqrt(2)/2).
     * @param input
     * @return 
    */
    public Image smoothToSigmaZeroPointSevenOne(Image input) {
        
        TriangleFunction scalingFunction = new TriangleFunction();
             
        GreyscaleImage[] smoothed = new GreyscaleImage[3];
        smoothed[0] = input.copyRedToGreyscale();
        smoothed[1] = input.copyGreenToGreyscale();
        smoothed[2] = input.copyRedToGreyscale();
        
        for (int i = 0; i < 3; ++i) {
            smoothed[i] = scalingFunction.calculateNextLevel(smoothed[i], 0);
        }
        
        Image out = input.createWithDimensions();
        for (int i = 0; i < input.getNPixels(); ++i) {
            out.setRGB(i, smoothed[0].getValue(i), 
                smoothed[1].getValue(i), smoothed[2].getValue(i));
        }
        
        return out;        
    }
    
    /**
     * The a trous algorithm is a fast implementation of a wavelet transform 
     * with no downsampling.   It is non-orthogonal, semi-linear runtime
     * complexity, is invariant under translation, and the transform is 
     * isotropic.
     * Implemented from pseudocode in http://www.multiresolution.com/svbook.pdf
       The scaling function used is the lower resolution choice, the triangle
       * function.
       * <pre>
       * The method uses recursive convolution operations, including previous
       * result to make next.
       * Each convolution uses two passes of one dimensional binomial kernels,
       * starting with the equivalent of sigma=sqrt(2)/2 = 0.707.
       * For each step, the equivalent resulting sigma is from 
       * sigma^2 = sigma_1^2 + sigma_2^2.
       * 
       * outputTransformed[1] = sigma=sqrt(2)/2 convolution
       * outputTransformed[2] = sigma=1 convolution
       * outputTransformed[3] = sqrt( (1)^2 + 1/2) = sqrt(1 + 1/2) convolution
       * outputTransformed[4] = sqrt( 1 + 1/2 + 1/2) = sqrt(2)
       * outputTransformed[5] = sqrt( 2 + 1/2)       = sqrt(2.5)
       * outputTransformed[6] = sqrt( 2 + 1/2 + 1/2) = sqrt(3)
       *  ...
       * </pre>
     * @param input
     * @param outputTransformed
     * @param outputCoeff 
     */
    public void calculateWithTriangleScalingFunction(GreyscaleImage input,
        List<GreyscaleImage> outputTransformed, List<GreyscaleImage> outputCoeff) {

        int imgDimen = Math.min(input.getWidth(), input.getHeight());

        int nr = (int)(Math.log(imgDimen)/Math.log(2));

        TriangleFunction scalingFunction = new TriangleFunction();
        
        outputTransformed.add(input.copyToSignedImage());
        
        outputCoeff.add(input.createSignedWithDimensions());

        for (int j = 0; j < nr; ++j) {
            
            GreyscaleImage cJ = outputTransformed.get(j);
 
            GreyscaleImage cJPlus1 = scalingFunction.calculateNextLevel(cJ, j);
           
            outputTransformed.add(cJPlus1);
            
            // c_(j,k) − c_(j+1,k)
            GreyscaleImage wJPlus1 = cJ.subtract(cJPlus1);//scalingFunction.subtractLevels(cJ, j);
          
            outputCoeff.add(wJPlus1);
        }
    }

    /**
     * The a trous algorithm is a fast implementation of a wavelet transform 
     * with no downsampling.   It is non-orthogonal, has semi-linear runtime
     * complexity, is invariant under translation, and the transform is 
     * isotropic.
     * Implemented from pseudocode in http://www.multiresolution.com/svbook.pdf
     * The scaling function used is the higher resolution choice, the 3rd order 
     * B Spline function.
     * <pre>
     * The method uses recursive convolution operations, including previous
       * result to make next.
       * Each convolution uses two passes of one dimensional binomial kernels,
       * starting with the equivalent of sigma=1.
       * For each step, the equivalent resulting sigma is from 
       * sigma^2 = sigma_1^2 + sigma_2^2.
       * 
       * outputTransformed[1] = sigma = 1 convolution
       * outputTransformed[2] = sqrt( (1)^2 + (1)^2 ) = sqrt(2) convolution
       * outputTransformed[3] = sqrt( 2 + 1 ) = sqrt(3) convolution
       * outputTransformed[4] = sqrt( 3 + 1 ) = sqrt(4) = 2 convolution
       * outputTransformed[5] = sqrt( 4 + 1 ) = sqrt(5) convolution
       * outputTransformed[6] = sqrt( 5 + 1 ) = sqrt(6) convolution
       * ...
       * outputTransformed[8] = sqrt( 8 + 1 ) = sqrt(9) = 3 convolution
       * </pre>
     * @param input
     * @param outputTransformed
     * @param outputCoeff 
     */
    public void calculateWithB3SplineScalingFunction(GreyscaleImage input,
        List<GreyscaleImage> outputTransformed, List<GreyscaleImage> outputCoeff) {

        int imgDimen = Math.min(input.getWidth(), input.getHeight());

        int nr = (int)(Math.log(imgDimen)/Math.log(2));

        B3SplineFunction scalingFunction = new B3SplineFunction();
        
        outputTransformed.add(input.copyToSignedImage());
        
        outputCoeff.add(input.createSignedWithDimensions());

        for (int j = 0; j < nr; ++j) {
            
            GreyscaleImage cJ = outputTransformed.get(j);
 
            GreyscaleImage cJPlus1 = scalingFunction.calculate(cJ);
           
            outputTransformed.add(cJPlus1);
            
            // c_(j,k) − c_(j+1,k)
            GreyscaleImage wJPlus1 = cJ.subtract(cJPlus1);
            
            outputCoeff.add(wJPlus1);
        }
        
    }
    
    /**
     * The a trous algorithm is a fast implementation of a wavelet transform 
     * with no downsampling.   It is non-orthogonal, has semi-linear runtime
     * complexity, is invariant under translation, and the transform is 
     * isotropic.
     * Implemented from pseudocode in http://www.multiresolution.com/svbook.pdf
     * The scaling function used is the higher resolution choice, the 3rd order 
     * B Spline function.
     * <pre>
     * The method uses recursive convolution operations, including previous
       * result to make next.
       * 
       * This method takes an argument stopIter to indicate that only stopIter
       * number of iterations are needed.  For instance, to retrieve only the
       * first populated wavelet, use stopIter = 2 (the first is initialization,
       * and the second is the calculation).
       * </pre>
     * @param input
     * @param outputTransformed
     * @param outputCoeff
     * @param stopIter
     */
    public void calculateWithB3SplineScalingFunction(GreyscaleImage input,
        List<GreyscaleImage> outputTransformed, List<GreyscaleImage> outputCoeff,
        int stopIter) {

        int imgDimen = Math.min(input.getWidth(), input.getHeight());

        int nr = (int)(Math.log(imgDimen)/Math.log(2));

        if (nr > stopIter) {
            nr = stopIter;
        }
        
        B3SplineFunction scalingFunction = new B3SplineFunction();
        
        outputTransformed.add(input.copyToSignedImage());
        
        outputCoeff.add(input.createSignedWithDimensions());

        for (int j = 0; j < nr; ++j) {
            
            GreyscaleImage cJ = outputTransformed.get(j);
 
            GreyscaleImage cJPlus1 = scalingFunction.calculate(cJ);
           
            outputTransformed.add(cJPlus1);
            
            // c_(j,k) − c_(j+1,k)
            GreyscaleImage wJPlus1 = cJ.subtract(cJPlus1);
            
            outputCoeff.add(wJPlus1);
        }
        
    }
   
    /**
     * same as calculateWithB3SplineScalingFunction except that it uses a 2-D
     * scaling function and the runtime is 2.5 times longer.
     * 
     * @param input
     * @param outputTransformed
     * @param outputCoeff 
     */
    void calculateWithB3SplineScalingFunction2(GreyscaleImage input,
        List<GreyscaleImage> outputTransformed, List<GreyscaleImage> outputCoeff) {

        int imgDimen = Math.min(input.getWidth(), input.getHeight());

        int nr = (int)(Math.log(imgDimen)/Math.log(2));

        B3SplineFunction scalingFunction = new B3SplineFunction();
        
        outputTransformed.add(input.copyToSignedImage());
        
        outputCoeff.add(input.createSignedWithDimensions());

        for (int j = 0; j < nr; ++j) {
            
            GreyscaleImage cJ = outputTransformed.get(j);
 
            GreyscaleImage cJPlus1 = scalingFunction.calculate2D(cJ);
           
            outputTransformed.add(cJPlus1);
            
            // c_(j,k) − c_(j+1,k)
            GreyscaleImage wJPlus1 = cJ.subtract(cJPlus1);
            
            outputCoeff.add(wJPlus1);
        }
        
    }

    public GreyscaleImage reconstruct(GreyscaleImage c0, 
        List<GreyscaleImage> mmCoeff) {

        int nr = mmCoeff.size();

        GreyscaleImage output = c0.copyToFullRangeIntImage();

        for (int j = 0; j < nr; ++j) {
           
            output = output.add(mmCoeff.get(j));
        }

        return output;
    }
    
    /**
     * Following
     * "Edge-Optimized À-Trous Wavelets for Local Contrast Enhancement with 
     * Robust Denoising" by Johannes Hanika, Holger Dammertz, and Hendrik Lensch
     * https://jo.dreggn.org/home/2011_atrous.pdf
     * to estimate del c_i_jj as part of creating an error image.
     * The authors calculate gradient c_i_jj using Cranley Patterson rotation 
       sampling within the A Trous B3Spline window (which is 25 pixels).
       
       This looks a little like calculating auto-correlation, except not wanting 
       the center pixel as the fixed first pixel of the difference.
       
       If gradient c_i_jj is meant to be a measure of local image noise, would 
       presumably want to select only differences between adjacent pixel pairs.
       So the use of Cranley Patterson rotation must be in selecting the second
       point using an offset chosen from the vector U of values.
       That offset is applied uniformly to the set to help choose the 2nd point.
       The universe of offsets U can only be the offsets to result in the 8
       neighbor region.
        
       Not sure, but I think that is what the authors implemented.
        
       Given to this method are the center pixel index for the A Trous window
       and the offsets as dx and dy chosen from the universe U of 8 neigbhor
       offsets.
        
       For each pixel in the window, will determine its intensity difference 
       from the pixel and the pixel that is it's coordinates plus the offsets.
       The result returned will be the average of those.
       
       Note that another paper 
       ("Efficient Multidimensional Sampling" by Kollig and Keller, 
       http://www.uni-kl.de/AG-Heinrich/EMS.pdf)
       suggests different sampling methods, so may change this in the future.
        
     * @return 
     */
    private double estimateLocalNoise(GreyscaleImage img, int pixIdx, int xOffset, 
        int yOffset) {
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        int x = img.getCol(pixIdx);
        int y = img.getRow(pixIdx);
        
        int count = 0;
        double diff = 0;
        // iterate within window to find first pixel
        for (int dx = -2; dx <= 2; ++dx) {
            int x1 = x + dx;
            if (x1 < 0 || x1 > (w - 1)) {
                continue;
            }
            for (int dy = -2; dy <= 2; ++dy) {
                int y1 = y + dy;
                if (y1 < 0 || y1 > (h - 1)) {
                    continue;
                }
                                
                int x2 = x1 + xOffset;
                int y2 = y1 + yOffset;
                if ((x2 < 0) || (x2 > (w - 1)) || (y2 < 0) ||
                    (y2 > (h - 1))) {
                    continue;
                }
                                
                diff += Math.abs(img.getValue(x1, y1) - img.getValue(x2, y2));
                
                count++;
            }
        }
        assert(count > 0);
        
        diff /= (double)count;
        
        return diff;
    }
}
