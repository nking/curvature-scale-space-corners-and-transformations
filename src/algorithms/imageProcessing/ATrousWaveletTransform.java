package algorithms.imageProcessing;

import java.util.List;

/**
 *
 * @author nichole
 */
public class ATrousWaveletTransform {

    /**
     * The a trous algorithm is a fast implementation of a wavelet transform 
     * with no downsampling.   It is non-orthogonal, semi-linear runtime
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
     * with no downsampling.   It is non-orthogonal, has semi-linear runtime
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
    
    /*
    looking at
    A Novel Edge-Aware A-Trous Filter For Single Image Dehazing `
        Baojun Qi, Tao Wu, and Hangen He
    2012 IEEE International Conference on Information Science and Technology
        Wuhan, Hubei, China; March 23-25, 2012
    
    the input transformed images are found above as c_(j,k), c_(j+1,k), etc
                                    (   | (c(j, k) - c(j + 1, k) | )
        edge stopping function = exp( - -------------------------- )
                                    (         (sigma_r)^2          )
    
        so this is using what the above refers to as the coefficients (retained
        for use in reconstruction).
    
    the transformed image is multiplied by the edge stopping function divided
        by the normalization.
    
    the normalization is needed too.
    
    looking for edge optimization which avoids exponentiation if possible.
    https://jo.dreggn.org/home/2011_atrous.pdf
    similar transform as implemented in above class plus
    locally adaptive (per-pixel) edge-weights, improved
    denoising using BayesShrink, and a fast and simple to implement
    algorithm
    
    see pg 37 of svn book
    */
}
