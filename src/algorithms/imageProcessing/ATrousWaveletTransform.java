package algorithms.imageProcessing;

import java.util.ArrayList;
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
     * The a trous algorithm is a fast implementation of a wavelet transform 
     * with no downsampling.   It is non-orthogonal, has semi-linear runtime
     * complexity, is invariant under translation, and the transform is 
     * isotropic.
     * Implemented from pseudocode in http://www.multiresolution.com/svbook.pdf
     * The scaling function used is the higher resolution choice, the 3rd order 
     * B Spline function.
     * Edge optimized factors have been included following
     * "Edge-Optimized À-Trous Wavelets for Local Contrast Enhancement with 
     * Robust Denoising" by Johannes Hanika, Holger Dammertz, and Hendrik Lensch
     * https://jo.dreggn.org/home/2011_atrous.pdf
     * 
     * @param input
     * @param outputTransformed
     * @param outputCoeff 
     */
    private void calculateForEdgeOptimization(GreyscaleImage input,
        List<GreyscaleImage> outputTransformed, List<GreyscaleImage> outputCoeff) {

        int imgDimen = Math.min(input.getWidth(), input.getHeight());

        int nr = (int)(Math.log(imgDimen)/Math.log(2));

        B3SplineFunction scalingFunction = new B3SplineFunction();
        
        outputTransformed.add(input.copyToSignedImage());
        
        outputCoeff.add(input.createSignedWithDimensions());
        
        int jjMax = 4;
        double lambda = 0.003;
        
        // ----- decomposition -----
        
        for (int j = 0; j < nr; ++j) {
            
            GreyscaleImage cJ = outputTransformed.get(j);
 
            GreyscaleImage cJPlus1 = scalingFunction.calculate(cJ);
            
            outputTransformed.add(cJPlus1);
            
            // c_(j,k) − c_(j+1,k)
            GreyscaleImage wJPlus1 = cJ.subtract(cJPlus1);
            
            outputCoeff.add(wJPlus1);
            
            /*
            double[][] eI = new double[jjMax*(j + 1)][];
            for (int jj = 0; jj < jjMax*(j + 1); ++jj) {
                float sigmaRJJ = 0.5f + jj;
                double sumW = 0;
                double[] ws = new double[wJPlus1.getNPixels()];
                for (int p = 0; p < wJPlus1.getNPixels(); ++p) {
                    assert(wJPlus1.getValue(p) >= 0);
                    double exp = Math.exp(
                        (wJPlus1.getValue(p) * wJPlus1.getValue(p))/sigmaRJJ);
                    ws[jj] = exp;
                    sumW += exp;
                }
                eI[jj] = new double[cJ.getNPixels()];
                for (int p = 0; p < cJ.getNPixels(); ++p) {
                    double cIJJ = (ws[p]/sumW) * cJ.getValue(p);
                    double dIJJ = cIJJ - cJPlus1.getValue(p);
                    // del c_i_jj calculated w/ Cranley Patterson rotation
                    double delC =
                    eI[jj][p] = (dIJJ * dIJJ) + (lambda * delC);
                }
                // evaluate error image:  e_jj = d_j_jj squared?  + λ · || ∇c_j_jj ||
            }
            */
            /*
            estimate sigma per pixel:
                compute multiple decompositions with different parameters and 
                then choose optimal parameters for each pixel. 
            The parameters to iterate over are sigma_r_j for jj=0 to jjmax.

            For each jj we separate the current image into coarse c_j_jj and
            detail using the current σ_r_jj as global edge weight. 
            The resulting decomposition is used to evaluate an error image
              e_jj = d_j_jj squared?  + λ · || ∇c_j_jj ||
            
            The authors use a fixed  λ = 0.003 and jmax=4.
            The error measure is evaluated using same window as used in b3 spline
            function, using 25 samples on a regular grid with 
            Cranley Patterson rotation [CP76] 
            to decorrelate adjacent pixels.
            
            After all jj have been processed at scale j, we search for the
            minimum error e_jj per pixel 
            and determine the optimal edge weight σ_r_k with k = argmin_jj{e_jj}. 
            In order to avoid too drastic changes in the edge weights, 
            the per-pixel weights are smoothed once into a buffer z
               z = σ_r_k ∗ h_0   (h_0 represents the b3 spline window function)
            
            Finally, the actual wavelet decomposition step takes place
            as in Section 2.1, but taking the locally optimal parameters
            σ_r from the blended buffer z at each pixel’s position.
            We found jmax = 4 sufficient for all images in this paper,
            and chose a simple linear sampling of the interval
            [0,4∗(j+1)). 
           
            Synthesis is recombining the coarse base layer with the coefficient
            (detail) layers that are potentially scaled and soft thresholded.
            
            BayesShrink [CYV00] is a soft thresholding method designed for a 
            wide range of input signal priors, such as Laplacian and 
            Gaussian distributions assuming additive iid Gaussian noise.
            As a consistent and robust estimator for the noise standard
            deviation, we use the median absolute deviation (MAD) over
            all pixels at the finest level, and obtain
                σ_n = median(|d0|) / 0.6745
            Now we iterate for each wavelet scale, this time starting at
            the coarsest level i = imax...0. Given the signal variance (σ_y_i)^2
            the soft shrinkage threshold T is calculated as in [CYV00]
            to minimize the risk of destroying the signal, by
                T = ( (σ_n_i)^2 ) / sqrt( max{0, (σ_y_i)^2 - (σ_n_i)^2} )

                with 
                   (σ_y_i)^2 = (1/N) * summation over pixels of (d_i)^2 
                and
                    (σ_n_i)^2 = σ_n,i = σ_n / (2^i)
            
                With this threshold, shrinkage and contrast boost is applied
                to the detail coefficients at the same time, by setting
                    modified d_i = max{0,|di| −T} ·sign(di)
                      and ci−1 = ci +β · modified d_i
                where β is the contrast boost parameter. In case denoising is
                not wanted, only the last step is relevant with modified d = d.

            Our edge-optimized wavelet scheme tends to move as
            much edge information to the coarser levels as possible, so
            the detail coefficients even correspond more to the Gaussian
            distribution than they would with usual wavelet transformation.
            Also, the edge data will not be affected by thresholding
            this way, resulting in a better PSNR of the denoised results.
            */
           
        }
        
        /*
        for (int i = i; i < nr; ++i) {
            
            GreyscaleImage w = ws.get(i);
            
            GreyscaleImage cJ = outputTransformed.get(i - 1);
            
            multiplyAndDivide(cJ, w, wSum);
        }
        */
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
