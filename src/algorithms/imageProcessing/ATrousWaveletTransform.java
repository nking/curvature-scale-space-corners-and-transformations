package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author nichole
 */
public class ATrousWaveletTransform {

    /**
    Uses the first step of an A Trous wavelet transform, which is two 1D 
    convolutions of a binomial kernel for sigma=1.
    * @param input
     * @return 
    */
    public GreyscaleImage smoothToSigmaOne(GreyscaleImage input) {
        
        B3SplineFunction scalingFunction = new B3SplineFunction();
        
        GreyscaleImage smoothed = scalingFunction.calculate(input);
        
        return smoothed;
    }
    
    /**
    Uses the first step of an A Trous wavelet transform, which is two 1D 
    convolutions of a binomial kernel for sigma=0.707 (=sqrt(2)/2).
     * @param input
     * @return 
    */
    public GreyscaleImage smoothToSigmaZeroPointSevenOne(GreyscaleImage input) {
        
        TriangleFunction scalingFunction = new TriangleFunction();
        
        GreyscaleImage smoothed = scalingFunction.calculateNextLevel(input, 0);
        
        return smoothed;
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
     * Edge optimized factors have been included following
     * "Edge-Optimized À-Trous Wavelets for Local Contrast Enhancement with 
     * Robust Denoising" by Johannes Hanika, Holger Dammertz, and Hendrik Lensch
     * https://jo.dreggn.org/home/2011_atrous.pdf
     * 
     *  Not Ready For Use.
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
        
        int lastIdx = -1;
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        // ----- decomposition -----
        
        for (int j = 0; j < nr; ++j) {
            
            GreyscaleImage cJ = outputTransformed.get(j);
 
            GreyscaleImage cJPlus1 = scalingFunction.calculate(cJ);
                        
            // c_(j,k) − c_(j+1,k)
            GreyscaleImage wJPlus1 = cJ.subtract(cJPlus1);
            
//NOTE: the recombined (synthesized) image looks correct, but the
// detail images (in outputCoeff) do not resemble the figure in 
// the authors' paper (unless the authors are plotting the unmodified 
// coefficient image).   
// the negative exponential diminishes the edges
// but a positive exponential, even when sigma is scaled for the range of
// data in wJPlus1 so that it doesn't overflow the range, only enhances the
// pixels that are already bright pixels in the coefficient image, so the
// step seems expensive (exponential transcendantal function) compared to 
// thresholding those same pixels.
//
// the authors' goals may already be met by the use of the A Trous wavelet and
// the b3 spline scaling function already.
            
            // -- alter cJPlus1 and wJPlus1 with edge optimization weights -----
          
            int w = cJ.getWidth();
            int h = cJ.getHeight();
            
            double[][] eI = new double[jjMax*(j + 1)][];
            for (int jj = 0; jj < jjMax*(j + 1); ++jj) {
                double sigmaRJJ = 0.5 + jj;
                double sum0 = 0;
                double sumW = 0;
                double[] ws = new double[wJPlus1.getNPixels()]; 
                for (int p = 0; p < wJPlus1.getNPixels(); ++p) {
                    double d = wJPlus1.getValue(p);
                    double m = (d*d)/sigmaRJJ;
                    //range of exp is 0 to 1 if use -m instead
                    double exp = Math.exp(-m);
                    ws[p] = exp;
                    sum0 += cJPlus1.getValue(p);
                    sumW += (exp * cJPlus1.getValue(p));
                }
               
                {
double[] cp = Arrays.copyOf(ws, ws.length);
Arrays.sort(cp);
System.out.println("j=" + j + " *minW=" + cp[0] + " maxW=" + cp[cp.length - 1]
+ " 1Qw=" +  cp[(int)(1.*(cp.length - 1)/4.)] + " 2Qw=" +  cp[(int)(3.*(cp.length - 1)/4.)]
+ " 3QW=" + cp[(int)(3.*(cp.length - 1)/4.)]);
                } 
                
                // normalize ws
                double f0 = (sum0/sumW);
                for (int p = 0; p < cJ.getNPixels(); ++p) {
                    ws[p] *= f0;
                }
                
                eI[jj] = new double[cJ.getNPixels()];
 
                {
double[] cp = Arrays.copyOf(ws, ws.length);
Arrays.sort(cp);
System.out.println("j=" + j + " minW=" + cp[0] + " maxW=" + cp[cp.length - 1]
+ " 1Qw=" +  cp[(int)(1.*(cp.length - 1)/4.)] + " 2Qw=" +  cp[(int)(3.*(cp.length - 1)/4.)]
+ " 3QW=" + cp[(int)(3.*(cp.length - 1)/4.)]);
                }
 
                for (int p = 0; p < cJ.getNPixels(); ++p) {
                                        
                    double dIJJ = cJ.getValue(p) - (ws[p] * cJPlus1.getValue(p));
                    
                    ++lastIdx;
                    if (lastIdx > (dxs.length - 1)) {
                        lastIdx = 0;
                    }

                    // NOTE: assuming authors intend gC to describe the unmodified cJ
                    double gC = estimateLocalNoise(cJ, p, dxs[lastIdx], dys[lastIdx]);
                    
                    eI[jj][p] = (dIJJ * dIJJ) + (lambda * gC);
                }
            }
            double[] sigmas = new double[cJ.getNPixels()];
            
            // for each pixel p, find the minimum in eI and store it's sigma.
            // then use a B3Spline on the sigma image to result in the sigmas
            // to be used per pixel on cJ and then calculate the detail image again
            for (int p = 0; p < cJ.getNPixels(); ++p) {
                double minE = Double.MAX_VALUE;
                int minEJJIdx = -1;
                for (int jj = 0; jj < jjMax*(j + 1); ++jj) {
                    double e = eI[jj][p];
                    if (e < minE) {
                        minE = e;
                        minEJJIdx = jj;
                    }
                }
                sigmas[p] =  0.5 + minEJJIdx;
            }
            sigmas = scalingFunction.calculate(sigmas, w, h);
            
            // use the sigmas and weighting function to alter cJPlus1 and recalc wJPlus1
            double sumW = 0;
            double sum0 = 0;
            double[] ws = new double[wJPlus1.getNPixels()];
            for (int p = 0; p < wJPlus1.getNPixels(); ++p) {
                double sigma = sigmas[p];
                assert(sigma > 0.0);
                double d = wJPlus1.getValue(p);
                double m = (d*d)/sigma;
                double exp = Math.exp(-m);
                ws[p] = exp;
                sum0 += cJPlus1.getValue(p);
                sumW += (exp * cJPlus1.getValue(p));
            }
            double f0 = (sum0/sumW);
            for (int p = 0; p < cJ.getNPixels(); ++p) {
                double f = ws[p] * f0;
                int cIJJPlus1 = (int)Math.round(f * cJPlus1.getValue(p));
                if (cIJJPlus1 > 255) {
                    System.err.println("larger than 255.  factor=" + f);
                    cIJJPlus1 = 255;
                }
                
                cJPlus1.setValue(p, cIJJPlus1);
            }
            
            wJPlus1 = cJ.subtract(cJPlus1);
                             
            outputTransformed.add(cJPlus1);
                        
            outputCoeff.add(wJPlus1);
            
            /*
            estimate sigma per pixel:
                compute multiple decompositions with different parameters and 
                then choose optimal parameters for each pixel. 
            The parameters to iterate over are sigma_r_j for jj=0 to jjmax.

            For each jj we separate the current image into coarse c_j_jj and
            detail using the current σ_r_jj as global edge weight. 
            The resulting decomposition is used to evaluate an error image
              e_jj = d_j_jj squared  + λ · || ∇c_j_jj ||
            
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
