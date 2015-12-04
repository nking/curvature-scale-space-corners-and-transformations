package algorithms.imageProcessing;

import java.util.List;

/**
 *
 * @author nichole
 */
public class ATrousWaveletTransform {

    /*
    implemented from pseudocode in http://www.multiresolution.com/svbook.pdf
    
    One assumes that the sampled data {c0,l} are the scalar products, 
    at pixels l, 
    of the function f(x) 
    with a scaling function φ(x) which corresponds to a low-pass filter.

    The wavelet function ψ(x) obeys the dilation equation:
         (1/2) * ψ(x/2) = summation over k( g(k) * φ(x - k) )
    The scalar projucts are calculated as 
          (1/(2^j))<f(x), ψ((x-l)/(2^j))>
    that is
         w_((j+1),l) = summation over k( g(k)*c_(j,(l+(2^j)*k))
    j=1 corresponds to finest scale (== high frequencies)
         w_((j+1),l) = c_(j,l) - c_(j+1,l)
    
         c_((j+1),l) = summation over k( h(k)*c_(j,(l+(2^j)*k))
    
    algorithm:
        1. We initialize j to 0 and we start with the data c_(j,k) where k
           is iteration over each pixel.
        2. We carry out a discrete convolution of the data c_(j,k) 
           using the filter h.
           The distance between the central pixel and the adjacent ones is 2^j.
    
           using the triangle function as the scaling function:
               c_(j + 1,k) = (1/4)*c_(j,k-(2^j)) + (1/2)*c_(j,k) + (1/4)*c_(j,k+(2^j))
    
            Note, could use another scaling function such as B_3 spline.
            (see page 272 of http://www.multiresolution.com/svbook.pdf)
        
        3. After this smoothing, we obtain the discrete wavelet transform from 
           the difference c_(j,k) − c_(j+1,k).
    
           using the triangle function as the scaling function:
               w_(j+1,k) = (-1/4)*c_(j,k-(2^j)) + (1/2)*c_(j,k) - (1/4)*c_(j,k+(2^j))
    
        4. If j is less than the number J of resolutions we want to compute, we
           increment j and then go to step 2.
        5. The set W = {w1 , ..., wJ , cJ } represents the wavelet transform 
           of the data.
    
    reconstruction algorithm:
        The algorithm allowing us to rebuild the data-frame is immediate: the
        last smoothed array c_J is added to all the differences, wj.
    
        c_(0,l) = c_(J,l) + summation over j(w_(j,l))
    
    */
    
    /**
     * The a trous algorithm is a fast implementation of a wavelet transform 
     * with no downsampling.   It is non-orthogonal, has reasonable runtime
     * complexity, is invariant under translation, and the transform is 
     * isotropic.
     * Implemented from pseudocode in http://www.multiresolution.com/svbook.pdf
       The scaling function used is the lower resolution choice, the triangle
       * function.
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
     * with no downsampling.   It is non-orthogonal, has reasonable runtime
     * complexity, is invariant under translation, and the transform is 
     * isotropic.
     * Implemented from pseudocode in http://www.multiresolution.com/svbook.pdf
     * The scaling function used is the higher resolution choice, the 3rd order 
     * B Spline function.
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

    public GreyscaleImage reconstruct(GreyscaleImage c0, 
        List<GreyscaleImage> mmCoeff) {

        int nr = mmCoeff.size();

        GreyscaleImage output = c0.copyToFullRangeIntImage();

        for (int j = 0; j < nr; ++j) {
           
            output = output.add(mmCoeff.get(j));
        }

        return output;
    }
}
