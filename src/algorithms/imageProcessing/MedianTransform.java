package algorithms.imageProcessing;

import algorithms.misc.MedianSmooth;
import java.util.List;

/**
 *
 * @author nichole
 */
public class MedianTransform {

    /**
     * A computationally expensive multiscale median transform.
     * see, the pyramidal mean transform version instead.
     * This method has a runtime complexity of n_iter * O(N_pixels * lg2(windowArea))
     * where windowArea grows from 1 to 2*2*lg2(imageDimension) + 1 
     * and nIter = lg2(imageDimension).
     * @param input
     * @param outputTransformed
     * @param outputCoeff
     */
    public void multiscaleMedianTransform(GreyscaleImage input,
        List<GreyscaleImage> outputTransformed, List<GreyscaleImage> outputCoeff) {

        int imgDimen = Math.min(input.getWidth(), input.getHeight());

        GreyscaleImage img0 = input.copyImage();

        int nr = (int)(Math.ceil(Math.log(imgDimen)/Math.log(2)));
        int s = 1;

        outputTransformed.add(img0.copyToSignedImage());
        
        outputCoeff.add(img0.createSignedWithDimensions());

        for (int j = 0; j < (nr - 1); ++j) {

            int winL = 2*s + 1;
            
            MedianSmooth med = new MedianSmooth();
            
            outputTransformed.add(med.calculate(outputTransformed.get(j), winL, winL));
            
            outputCoeff.add(outputTransformed.get(j).subtract(outputTransformed.get(j + 1)));

            s = 2*s;
        }
        
    }

    public GreyscaleImage reconstructMultiscaleMedianTransform(GreyscaleImage
        c0, List<GreyscaleImage> mmCoeff) {

        int nr = mmCoeff.size();

        GreyscaleImage output = c0.copyToSignedImage();

        for (int j = 0; j < nr; ++j) {
            output = output.add(mmCoeff.get(j));
        }

        return output;
    }

    /**
     * pyramidal median transform (faster than multiscalePyramidalMedianTransform
     * but reconstruction from coefficients is not exact, so prefer
     * multiscalePyramidalMedianTransform if exact is needed);
     * following pseudocode in http://www.multiresolution.com/svbook.pdf
     * 
     * @param input
     * @param outputTransformed
     * @param outputCoeff 
     */
    public void multiscalePyramidalMedianTransform2(GreyscaleImage input,
        List<GreyscaleImage> outputTransformed, List<GreyscaleImage> outputCoeff) {

        int imgDimen = Math.min(input.getWidth(), input.getHeight());

        GreyscaleImage img0 = input.copyImage();

        int nr = (int)(Math.ceil(Math.log(imgDimen)/Math.log(2)));
        int s = 1;
        int winL = 2*s + 1;
        
        ImageProcessor imageProcessor = new ImageProcessor();

        outputTransformed.add(img0.copyToSignedImage());
        
        outputCoeff.add(img0.createSignedWithDimensions());

        for (int j = 0; j < (nr - 1); ++j) {
            
            MedianSmooth med = new MedianSmooth();
            
            GreyscaleImage cJ = outputTransformed.get(j);
            
            GreyscaleImage cJPlus1Ast = med.calculate(cJ, winL, winL);   
            
            // decimation:
            GreyscaleImage cJPlus1 = imageProcessor.binImage(cJPlus1Ast, 2);
            
            GreyscaleImage wJPlus1 = cJ.subtract(cJPlus1Ast);
            
            outputTransformed.add(cJPlus1);
            
            outputCoeff.add(wJPlus1);
            
            assert(cJ.getWidth() == wJPlus1.getWidth());
        }
        
        // empty full size image
        outputCoeff.remove(0);
    }
    
    /**
     * pyramidal median transform for exact reconstruction.
     * following pseudocode in http://www.multiresolution.com/svbook.pdf
     * 
     * This method has a runtime complexity of 
     * n_iter * (O(N_pixels * 1.6) + 5*O(N_pixels))
     * where nIter = lg2(imageDimension) - 1 and N_pixels is decreasing
     * in size by a factor of 4 for each iteration.
     * @param input
     * @param outputTransformed
     * @param outputCoeff 
     */
    public void multiscalePyramidalMedianTransform(GreyscaleImage input,
        List<GreyscaleImage> outputTransformed, List<GreyscaleImage> outputCoeff) {

        int imgDimen = Math.min(input.getWidth(), input.getHeight());

        GreyscaleImage img0 = input.copyImage();

        int nr = (int)(Math.ceil(Math.log(imgDimen)/Math.log(2)));
        int s = 1;
        int winL = 2*s + 1;
        
        ImageProcessor imageProcessor = new ImageProcessor();

        outputTransformed.add(img0.copyToSignedImage());
        
        outputCoeff.add(img0.createSignedWithDimensions());

        for (int j = 0; j < (nr - 1); ++j) {
            
            MedianSmooth med = new MedianSmooth();
            
            GreyscaleImage cJ = outputTransformed.get(j);
                        
            // median filter and decimation:
            GreyscaleImage cJPlus1 = imageProcessor.binImage(med.calculate(cJ, winL, winL), 2);
            
            //interpolation of cJPlus1 to size cJ
            GreyscaleImage cJPlus1Ast = imageProcessor.expandBy2UsingBilinearInterp(cJPlus1);
            
            GreyscaleImage wJPlus1 = cJ.subtract(cJPlus1Ast);
            
            outputTransformed.add(cJPlus1);
            
            outputCoeff.add(wJPlus1);
            
            assert(cJ.getWidth() == wJPlus1.getWidth());
        }
        
        // empty full size image
        outputCoeff.remove(0);
    }

    /**
     * reconstruct image from products of pyramidal median transform.
     * following pseudocode in http://www.multiresolution.com/svbook.pdf
     * 
     * @param c0
     * @param mmCoeff
     * @return 
     */
    public GreyscaleImage reconstructPyramidalMultiscaleMedianTransform(
        GreyscaleImage c0, List<GreyscaleImage> mmCoeff) {

        int nr = mmCoeff.size();

        ImageProcessor imageProcessor = new ImageProcessor();

        GreyscaleImage output = c0.copyToSignedImage();

        for (int j = (nr - 1); j > -1; --j) {

            // expand by factor of 2.
            GreyscaleImage cJPrime = imageProcessor.expandBy2UsingBilinearInterp(output);

            GreyscaleImage wJ = mmCoeff.get(j);
            
            output = cJPrime.add(wJ);
        }

        return output;
    }

    public void multiscaleMedianWaveletTransform(GreyscaleImage input) {

        if (true) {
            throw new UnsupportedOperationException("not yet implemented");
        }
        
        /*
        from: Sparse Image and Signal Processing, Second Edition,
        by Starck, Murtagh, and Fadili

        estimate st dev using Donoho and Johnstone (1994) based on wavelet
        coeff of noisy data Y at the finest resolution level.
        The wavelet coeff of Y at finest scale tend to be mostly noise, while
        wavelet coeff of X at same scale can be viewed as outliers.
        sigma = MAD(w_1)/0.6745 = median(|w_1 - median(w_1)|)/0.6745

        where MAD stands for the median absolute deviation
        w_1 are the orthogonal wavelet coefficients of Y at the finest scale.
        For 2-D images, the above estimator is to be applied with the
        diagonal subband of the 2-D separable orthogonal wavelet transform.

We now turn to the estimation of . As the noise is additive, we have , and it is easy to see that


If the atoms in the dictionary  all have equal unit -norm, then obviously .
        This formula is also easily implemented if the -norms were known
        analytically, as is the case for the curvelet tight frame
        (see Section 5.4.2.2). But if these norms are not known in closed
        form, they can be estimated in practice by taking the transform
        of a Dirac, and then computing the -norm of each subband.
        */
        /*
         detect in w_(j+1) the significant coefficients:
              |w_(j+1)| > tau * MAD(w_(j+1))/0.6745
              where MAD stands for the median absolute deviation used as an estimator
              of the noise standard deviation.  see eqn (6.9) and tau a threshold
              chosen large enough to avoid false detections, for instance tau=5.
          set to zero all significant coefficients in w_(j+1).
          compute c_prime_j = w_(j+1) + c_(j+1).  hence c_prime_j is a version of c_j,
              but without the detected significant structures.
          compute the 2D starlet transform of c_prime_j with j+1 scales.
              we get w={w_prime_j,...w_prime_(j+1), c_prime_(j+1)}
          set c_(j+1) = c_prime_(j+1). therefore, c_(j+1) is smoothed with wavelets,
              but strong features have been extracted with median
          compute the median-wavelet coefficients: w_(j+1) = c_j - c_(j+1).
          s = 2*s

          output: w={w_1,...w_j, c_j}
        */
        
    }

}
