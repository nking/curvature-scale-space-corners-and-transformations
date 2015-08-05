package algorithms.imageProcessing;

/**
                   1             ( -(x - mu)^2 )
     f = ------------------ * exp( ----------- )
         sigma * sqrt(2*pi)      (    2o~^2    )
       
 * @author nichole
 */
public class Gaussian1D {
   
    protected static float estimateHWZI(SIGMA sigma, float fractionMax) {
      
        float s = SIGMA.getValue(sigma);

        return estimateHWZI(s, fractionMax);
    }
    
    protected static float estimateHWZI(float sigma, float fractionMax) {
        
        /*
           exp( (-(x0 - mu)^2)/2o~^2) / exp( (-(xcenter - mu)^2)/2o~^2) = fractionMax
        
           where xcenter = mu is the center of the gaussian
        
           exp( (-(x0 - mu)^2)/2o~^2) / 1 = fractionMax
        
           (-(x0 - mu)^2)/2o~^2) = ln(fractionMax)
        
           -(x0 - mu)^2 = -1*(2o~^2 * ln(fractionMax))
        
           x0 - mu = math.sqrt(-1*(2o~^2 * ln(fractionMax)))
        
           x0 = mu + o~ * math.sqrt(-1*(2 * ln(fractionMax)))
        
        */

        float x0 = (float)(sigma * Math.sqrt(-2. * Math.log(fractionMax)));
                
        return x0;
    }
  
    /**
     * get the kernel for sigma and deltaX.  
     * 
     * To create the x-axis for the kernel:
       <pre>
       int halfWidthInPixels = kernel.length >> 1; 
       int[] x = new int[kernel.length];
       for (int i = 0; i < x.length; i++) {
            x[i] = i + (int) (mu - halfWidthInPixels);
       }
       </pre>
     * 
     * 
     * @param sigma
     * @return 
     */
    public static float[] getKernel(float sigma) {
        
        if (Math.abs(sigma - 0.4247f) < 0.01) {
            //0.42466090014400953f
            return getBinomialKernelSigmaZeroPointFourTwoSeven();
        }
        
        return getKernel(sigma, 0);
    }
    
    public static float[] getKernel(float sigma, float mu) {

        float normalization = (float)(sigma * Math.sqrt(2.f * Math.PI));
        
        float hwi = estimateHWZI(sigma, 0.001f);
        if (hwi < 0) {
            hwi *= -1.f;
        }
        int halfWidthInPixels = (int)Math.ceil(hwi);
        
        int start = -1*halfWidthInPixels;
        int stopExcl = halfWidthInPixels + 1;
        
        float d, dsq;
       
        int nPoints = stopExcl - start;
        
        float[] yPoints = new float[nPoints];
        int count = 0;
        for (int i = start; i < stopExcl; i++) {
            
            float x = i;
            
            d = (x - mu);
            dsq = d*d;
            
            float y = (float) Math.exp(-1.f * dsq/(2.f * sigma * sigma));
            
            yPoints[count] = y/normalization;
                        
            count++;
        }
                                
        return yPoints;
    }
    
    /**
     * get the kernel for sigma.  
     * 
     * To create the x-axis for the kernel:
       <pre>
       int halfWidthInPixels = kernel.length >> 1; 
       int[] x = new int[kernel.length];
       for (int i = 0; i < x.length; i++) {
            x[i] = i + (int) (mu - halfWidthInPixels);
       }
       </pre>
     * 
     * 
     * @param sigma
     * @return 
     */
    public static float[] getKernel(SIGMA sigma) {
        
        if (sigma.ordinal() == SIGMA.ZEROPOINTFIVE.ordinal()) {
            return getBinomialKernelSigmaZeroPointFive();
        } else if (sigma.ordinal() == SIGMA.ZEROPOINTSEVENONE.ordinal()) {
            return getBinomialKernelSigmaZeroPointSevenOne();
        } else if (sigma.ordinal() == SIGMA.ONE.ordinal()) {
            return getBinomialKernelSigmaOne();
        } else if (sigma.ordinal() == SIGMA.ONESQRT2.ordinal()) {
            return getBinomialKernelSigmaOneSQRT2();
        } else if (sigma.ordinal() == SIGMA.ONEPOINTFIVE.ordinal()) {
            return getBinomialKernelSigmaOnePointFive();
        } else if (sigma.ordinal() == SIGMA.TWO.ordinal()) {
            return getKernelSigmaTwo();
        } else if (sigma.ordinal() == SIGMA.TWOSQRT2.ordinal()) {
            return getKernelSigmaTwoSQRT2();
        } else if (sigma.ordinal() == SIGMA.THREE.ordinal()) {
            return getKernelSigmaThree();
        } else if (sigma.ordinal() == SIGMA.FOUR.ordinal()) {
            return getKernelSigmaFour();
        } else if (sigma.ordinal() == SIGMA.FOURSQRT2.ordinal()) {
            return getKernelSigmaFourSQRT2();
        } else if (sigma.ordinal() == SIGMA.EIGHT.ordinal()) {
            return getKernelSigmaEight();
        } else {
            return getKernel(SIGMA.getValue(sigma));
        }
    }
    
    protected static float[] getBinomialKernel(SIGMA sigma) {
        
        if (sigma.ordinal() == SIGMA.ZEROPOINTFIVE.ordinal()) {
            return getBinomialKernelSigmaZeroPointFive();
        } else if (sigma.ordinal() == SIGMA.ZEROPOINTSEVENONE.ordinal()) {
            return getBinomialKernelSigmaZeroPointSevenOne();
        } else if (sigma.ordinal() == SIGMA.ONE.ordinal()) {
            return getBinomialKernelSigmaOne();
        } else if (sigma.ordinal() == SIGMA.ONESQRT2.ordinal()) {
            return getBinomialKernelSigmaOneSQRT2();
        } else if (sigma.ordinal() == SIGMA.ONEPOINTFIVE.ordinal()) {
            return getBinomialKernelSigmaOnePointFive();
        } else {
            return null;
        }
    }
    
    /**
     * for sigma = 0.42466090014400953, this results in multiplication by 1.
     */ 
    protected static float[] getBinomialKernelSigmaZeroPointFourTwoSeven() {
        return new float[]{1};
    }
    
    protected static float[] getBinomialKernelSigmaZeroPointFive() {
        // 1 1, normalization = 1/2^1 = 0.5 ... replacing w/ 1  8  1
        float a = 1.f/9.f;
        return new float[]{a, 8*a, a};
    }
    
    protected static float[] getBinomialKernelSigmaZeroPointSevenOne() {
        // 1 2 1,  norm = 1./2^2 = 0.25... actually was 0.33333
        return new float[]{0.25f, 0.5f, 0.25f};
        //return new float[]{0.3333333333333333f, 0.6666666666666666f,
        //    0.3333333333333333f};
    }
    
    protected static float[] getBinomialKernelSigmaOne() {
        // 1 4 6 4 1,
        return new float[]{0.0625f, 0.25f, 0.375f, 0.25f, 0.0625f};
    }
    
    protected static float[] getBinomialKernelSigmaOneSQRT2() {
        // 1     8    28       56      70        56      28      8     1
        // norm = 1./2^8 = 0.00390625
        float a = 1.f/256.f;
        return new float[]{a*1, a*8, a*28, a*56, 
            a*70, a*56, a*28, a*8, 
            a*1};
    }
    
    protected static float[] getBinomialKernelSigmaOnePointFive() {
        // 1  9    36     84      126            126     84      36      9     1
        // edited to:
        // 1  5    22     66      126   <160>    126     66      22      5     1
        float a = 1.f/600.f;
        return new float[]{a, 5*a, 22*a, 66*a, 126*a, 160*a, 126*a, 66*a, 22*a,
            5*a, a};
    }
  
    protected static float[] getKernelSigmaTwo() {
        return new float[]{(float)6.691511E-5, (float)4.3634136E-4, 
            0.0022159242f, 0.008764151f, 0.026995484f, 0.0647588f, 0.12098537f, 
            0.17603266f, 0.19947115f, 0.17603266f, 0.12098537f, 0.0647588f, 
            0.026995484f, 0.008764151f, 0.0022159242f, (float)4.3634136E-4, 
            (float)6.691511E-5};
    }
    
    protected static float[] getKernelSigmaTwoSQRT2() {
        return new float[]{(float)7.3284624E-5, (float)2.722854E-4, 
            (float)8.927895E-4, 0.002583372f, 0.006596873f, 0.014866283f, 
            0.029565137f, 0.051888432f, 0.08036638f, 0.10984783f, 0.13250177f, 
            0.1410474f, 0.13250177f, 0.10984783f, 0.08036638f, 0.051888432f,
            0.029565137f, 0.014866283f, 0.006596873f, 0.002583372f, 
            (float)8.927895E-4, (float)2.722854E-4, (float)7.3284624E-5};
    }
    
    protected static float[] getKernelSigmaThree() {
        return new float[]{(float)4.4610075E-5, (float)1.6009022E-4, 
            (float)5.1409315E-4, 0.0014772828f, 0.003798662f, 0.008740629f, 
            0.01799699f, 0.033159047f, 0.054670025f, 0.08065691f, 0.10648267f, 
            0.12579441f, 0.13298076f, 0.12579441f, 0.10648267f, 0.08065691f, 
            0.054670025f, 0.033159047f, 0.01799699f, 0.008740629f, 0.003798662f, 
            0.0014772828f, (float)5.1409315E-4, (float)1.6009022E-4, 
            (float)4.4610075E-5};
    }
    
    protected static float[] getKernelSigmaFour() {
        return new float[]{(float)8.814892E-5, (float)2.1817068E-4, 
            (float)5.07262E-4, 0.0011079621f, 0.0022733908f, 0.0043820753f, 
            0.007934913f, 0.013497742f, 0.02156933f, 0.0323794f, 0.04566227f, 
            0.060492683f, 0.075284354f, 0.08801633f, 0.09666703f, 0.09973557f, 
            0.09666703f, 0.08801633f, 0.075284354f, 0.060492683f, 0.04566227f, 
            0.0323794f, 0.02156933f, 0.013497742f, 0.007934913f, 0.0043820753f, 
            0.0022733908f, 0.0011079621f, (float)5.07262E-4, 
            (float)2.1817068E-4, (float)8.814892E-5};
    }
    
    protected static float[] getKernelSigmaFourSQRT2() {
        return new float[]{(float)3.6642312E-5, (float)7.1742164E-5, 
            (float)1.361427E-4, (float)2.5040476E-4, (float)4.4639476E-4, 
            (float)7.713009E-4, 0.001291686f, 0.0020966139f, 0.0032984365f, 
            0.0050295154f, 0.0074331416f, 0.010647485f, 0.014782568f, 
            0.01989212f, 0.025944216f, 0.032796565f, 0.04018319f, 0.047718722f, 
            0.054923914f, 0.061272055f, 0.06625088f, 0.06943033f, 0.0705237f, 
            0.06943033f, 0.06625088f, 0.061272055f, 0.054923914f, 0.047718722f, 
            0.04018319f, 0.032796565f, 0.025944216f, 0.01989212f, 0.014782568f, 
            0.010647485f, 0.0074331416f, 0.0050295154f, 0.0032984365f, 
            0.0020966139f, 0.001291686f, (float)7.713009E-4, 
            (float)4.4639476E-4, (float)2.5040476E-4, (float)1.361427E-4, 
            (float)7.1742164E-5, (float)3.6642312E-5};
    }
    
    protected static float[] getKernelSigmaEight() {
        return new float[]{(float)4.407446E-5, (float)6.988269E-5, 
            (float)1.0908534E-4, (float)1.6763987E-4, (float)2.53631E-4, 
            (float)3.7778224E-4, (float)5.5398105E-4, (float)7.9976505E-4, 
            0.0011366954f, 0.0015905226f, 0.0021910376f, 0.0029714876f, 
            0.0039674565f, 0.0052151233f, 0.006748871f, 0.008598284f, 
            0.010784665f, 0.013317284f, 0.0161897f, 0.019376533f, 0.022831134f, 
            0.02648458f, 0.030246342f, 0.034006875f, 0.037642177f, 0.04102012f, 
            0.044008166f, 0.04648189f, 0.048333514f, 0.04947971f, 0.049867786f, 
            0.04947971f, 0.048333514f, 0.04648189f, 0.044008166f, 0.04102012f, 
            0.037642177f, 0.034006875f, 0.030246342f, 0.02648458f, 
            0.022831134f, 0.019376533f, 0.0161897f, 0.013317284f, 
            0.010784665f, 0.008598284f, 0.006748871f, 0.0052151233f, 
            0.0039674565f, 0.0029714876f, 0.0021910376f, 0.0015905226f, 
            0.0011366954f, (float)7.9976505E-4, (float)5.5398105E-4, (float)3.7778224E-4, 
            (float)2.53631E-4, (float)1.6763987E-4, (float)1.0908534E-4, (float)6.988269E-5, 
            (float)4.407446E-5};
    }
    
    protected static double[] getHalfKernelUsingBinomialFilterSigmaOne() {
        // 1, 4, 6, 4, 1  with normalization = 16
        return new double[]{0.0625f, 0.25f, 0.375f};
    }
    
    /*<pre>
                                       binomial filter
    
    http://venus.inrialpes.fr/jlc/papers/Crowley-ScaleSpace03.pdf
    
    binomial filters are obtained w/ a cascaded convolution of a kernel composed of [1,1].

    the coefficients for the nth filter in the series, b_n(m) are:
        b_n(m) = [1,1]^(*n)  where (*n) is n auto-convolutions.

    This series provide the best (least sum of squares error) approximation to a
    Gaussian function by an integer coefficient sequence of finite duration

    given the nth binomial, there are n coefficients whose sum is 2^n.
    the mid point is the coefficients at m=n/2 and the variance is sigma^2=n/4

    b_2(m) = [1,2,1]  is a monotonic low-pass filter with no ripples in the stop band

    b_4(m) = [1, 4, 6, 4, 1]
    
    The filters b_2(m) and b_4(m) have variances of 0.5 and 1, respectively.
    The filter b_4(m) is equivalent to b_2(m) * b_2(m).
    Thus, a s=1 Gaussian filter can be computed by two convolutions with the
    kernel [1, 2, 1] at a cost of two multiplications and 4 additions per pixel.

    [1,1] * [1,1] -> [1,2,1]        
    [1,1] * [1,2,1] -> [1,3,3,1]    
    [1,1] * [1,3,3,1] -> [1,4,6,4,1]
    
    can form a triangle where
    
    n is the number of the 1-D elements in the filter - 1

    r = position of elemen in the filter kernel (0, 1, 2, ...)
               n!       ( n )
    a_n_r = --------- = ( r )
            r!*(n-r)!

    midpoint of a filter is n/2.  sigma=sqrt(n)/2
    0.1  0.8  0.1
  n=0                                          1                          
  n=1 replacing with 0.1  0.8  0.1?       1         1                     sigma=sqrt(1)/2 = 0.5
  n=2                                 1        2         1                sigma=sqrt(2)/2 = 0.707
  n=3                             1      3         3        1            sigma=sqrt(3)/2 = 0.866
  n=4                          1      4        6         4       1   <=== sigma = 1 
  n=5                      1       5      10        10       5       1
  n=6                   1      6      15       20        15      6       1
  n=7                1     7      21      35       35       21       7     1
  n=8             1     8    28       56       70        56      28      8     1      <=== sigma=sqrt(2)
  n=9          1     9    36     84      126       126      84      36      9     1   <=== sigma=3/2
     
  recursive convolution operations to reuse kernel and save computation time
  for cases when the entire pyramid is needed such as when making scale space
  images.
    sigma=sqrt(1)/2            has n=1.   if convolved w/ self => sigma=sqrt(2)
    sigma=sqrt(2)/2 = 0.707    has n=2.   if convolved w/ self => sigma=1
    sigma=1                    has n=4.   if convolved w/ self => sigma=sqrt(2)
    sigma=sqrt(2)              has n=8.   if convolved w/ self => sigma=2
    sigma=2                    has n=16.  if convolved w/ self => sigma=2*sqrt(2)
    sigma=2*sqrt(2)            has n=32.  if convolved w/ self => sigma=4
    sigma=4                    has n=64.  if convolved w/ self => sigma=4*sqrt(2)
    sigma=4*sqrt(2)            has n=128. if convolved w/ self => sigma=8
    sigma=8                    has n=256. if convolved w/ self => sigma=8*sqrt(2)
    sigma=8*sqrt               has n=512. if convolved w/ self => sigma=16
 
  also, the sigma=3 chain, used in ECSS:
    sigma=3/2           has n=9.   if convolved w/ self => 3*sqrt(2)/2
    sigma=3*sqrt(2)/2   has n=18.  if convolved w/ self => sigma=3
    sigma=3             has n=36.
    
    NOTE: if compute gaussian and find where curve falls to 0.01 of peak 
       we can derive the number of points in a kernal as:
           sigma=2, numberOfPoints = 12.14  compared to n=16 for binomial
           sigma=4, numberOfPoints = 24.28              n=64 for binomial
           sigma=8, numberOfPoints = 48.56              n=128 for binomial
    
    Some of the derivatives worked out for sigma=1:
         b_4(m) = [1, 4, 6, 4, 1]  
         first_deriv_4(m) = [1, 3, 2, -2, -3, -1]
         second_deriv_4(m) = [1, 2, -1, -4, -1, 2, 1]
    </pre>
    */
    
}
