package algorithms.imageProcessing;

import static algorithms.imageProcessing.Gaussian1D.estimateHWZI;

/**
 * class to retrieve a kernel for the first derivative of a gaussian of
 * a given sigma.
 * 
 * For the smaller sigma, curvature tests show that a binomial filter should
 * be used instead of the detailed calculation of the formula due to errors
 * in normalization due to the tails of the function.
 * Therefore, the general methods are returning binomial filter curves
 * when possible.
 * Sometimes those have been edited from the binomial pattern
 * due to centering.
 * 
 * @author nichole
 */
public class Gaussian1DFirstDeriv {
        
    //               1             ( -(x - mu)^2 )
    // f = ------------------ * exp( ----------- )
    //     sigma * sqrt(2*pi)      (    2o~^2    )
    //
    //           -(x - mu)              ( -(x - mu)^2 )
    // dfdx = -------------------- * exp( ----------- )
    //        sigma^3 * sqrt(2*pi)      (    2o~^2    )
    //
        
    /**
     * get the kernel for the given sigma and mu.  nPoints is the number of
     * points in the Gaussian kernel (determined by the approximate FWZI), so
     * mu - half of nPoints determines the starting x in the kernel.
     * 
     * the function is:
     * <pre>
     *           -(x - mu)              ( -(x - mu)^2 )
     * dfdx = -------------------- * exp( ----------- )
     *        sigma^3 * sqrt(2*pi)      (    2o~^2    )
     * 
     * </pre>
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
     * @param sigma
     * @param mu
     * @param nPoints
     * @return 
     */
    public static float[] getKernel(SIGMA sigma, float mu, int nPoints) {
        
        float s = SIGMA.getValue(sigma);
        
        return getKernel(s, mu, nPoints);
    }
    
    /**
     * get the kernel for the given sigma and mu.  nPoints is the number of
     * points in the Gaussian kernel (determined by the approximate FWZI), so
     * mu - half of nPoints determines the starting x in the kernel.
     * 
     * the function is:
     * <pre>
     *           -(x - mu)              ( -(x - mu)^2 )
     * dfdx = -------------------- * exp( ----------- )
     *        sigma^3 * sqrt(2*pi)      (    2o~^2    )
     * 
     * </pre>
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
     * @param sigma
     * @param mu
     * @param nPoints
     * @return 
     */
    public static float[] getKernel(float sigma, float mu, int nPoints) {
                
        float normalization = (float)(Math.pow(sigma, 3) * 
            Math.sqrt(2.f * Math.PI));
        
        int halfWidthInPixels = nPoints >> 1;
        
        int start = (int)(mu - halfWidthInPixels);
        int stopExcl = (int)mu + halfWidthInPixels + 1;
            
        float d;
        float dsq;
                        
        int count = 0;
        float[] yPoints = new float[nPoints];
        int i;
        for (i = start; i < stopExcl; i++) {
             
            float x = i;
            
            d = (x - mu);
            dsq = d*d;
            
            float y = (float)((-1.f * d) *
                Math.exp(-1.f * dsq/(2.f * sigma * sigma)));
                        
            yPoints[count] = y/normalization;
                        
            count++;
        }
      
        return yPoints;
    }
    
    /**
     * get the kernel for the given sigma and mu. 
     * 
     * the function is:
     * <pre>
     *           -(x - mu)              ( -(x - mu)^2 )
     * dfdx = -------------------- * exp( ----------- )
     *        sigma^3 * sqrt(2*pi)      (    2o~^2    )
     * 
     * </pre>
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
     * @param sigma
     * @param mu
     * @return 
     */
    public static float[] getKernel(float sigma, float mu) {

        float hwi = estimateHWZI(sigma, 0.001f);
        if (hwi < 0) {
            hwi *= -1.f;
        }
        int halfWidthInPixels = (int)Math.ceil(hwi);
        
        int nPoints = 2*halfWidthInPixels + 1;
        
        return getKernel(sigma, mu, nPoints);
    }
    
    /**
     * get the kernel for the given sigma. a mu of 0 is used. 
     * 
     * the function is:
     * <pre>
     *           -(x - mu)              ( -(x - mu)^2 )
     * dfdx = -------------------- * exp( ----------- )
     *        sigma^3 * sqrt(2*pi)      (    2o~^2    )
     * 
     * </pre>
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
     * @param sigma
     * @return 
     */
    public static float[] getKernel(float sigma) {
        
        return getKernel(sigma, 0);
    }
 
    /**
     * get the kernel for the given sigma and deltaX. NOTE: the factor -(x-mu)
     * has not been applied yet to allow re-use of the kernel, so when using this
     * kernel, be sure to multiply each element by -(x-mu);
     * 
     * the function is:
       <pre>
                 -(x - mu)              ( -(x - mu)^2 )
       dfdx = -------------------- * exp( ----------- )
              sigma^3 * sqrt(2*pi)      (    2o~^2    )
       
       </pre>
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
     * @param sigma
     * @return 
     */
    public static float[] getKernel(SIGMA sigma) {
      
        //using the binomial filter removes artifact sometimes introduced in
        // convolution due to nearly zero elements in the kernel and subsequent
        // renormalization by a very small number, so using binomial when
        // possible
        
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
            return getBinomialKernelSigmaTwo();
        } else if (sigma.ordinal() == SIGMA.ONEPOINTFIVESQRT2.ordinal()) {
            return getBinomialKernelSigmaOnePointFiveSQRT2();
        } else if (sigma.ordinal() == SIGMA.TWOSQRT2.ordinal()) {
            return getBinomialKernelSigmaTwoSQRT2();
        } else if (sigma.ordinal() == SIGMA.THREE.ordinal()) {
            return getBinomialKernelSigmaThree();
        } else if (sigma.ordinal() == SIGMA.FOUR.ordinal()) {
            return getKernelSigmaFour();
        } else if (sigma.ordinal() == SIGMA.FOURSQRT2.ordinal()) {
            return getKernel(SIGMA.getValue(sigma));
        } else if (sigma.ordinal() == SIGMA.EIGHT.ordinal()) {
            return getKernel(SIGMA.getValue(sigma));
        } else if (sigma.ordinal() == SIGMA.EIGHTSQRT2.ordinal()) {
            return getKernel(sigma, 0, 87);
        } else if (sigma.ordinal() == SIGMA.SIXTEEN.ordinal()) {
            return getKernel(sigma, 0, 121);
        } else if (sigma.ordinal() == SIGMA.SIXTEENSQRT2.ordinal()) {
            return getKernel(sigma, 0, 171);
        } else if (sigma.ordinal() == SIGMA.THIRTYTWO.ordinal()) {
            return getKernel(sigma, 0, 239);
        } else if (sigma.ordinal() == SIGMA.THIRTYTWOSQRT2.ordinal()) {
            return getKernel(sigma, 0, 339);
        } else if (sigma.ordinal() == SIGMA.SIXTYFOUR.ordinal()) {
            return getKernel(sigma, 0, 477);
        } else if (sigma.ordinal() == SIGMA.SIXTYFOURSQRT2.ordinal()) {
            return getKernel(sigma, 0, 675);
        } else if (sigma.ordinal() == SIGMA.ONEHUNDREDANDTWENTYEIGHT.ordinal()) {
            return getKernel(sigma, 0, 953);
        } else if (sigma.ordinal() == SIGMA.ONEHUNDREDANDTWENTYEIGHTSQRT2.ordinal()) {
            return getKernel(sigma, 0, 1347);
        } else if (sigma.ordinal() == SIGMA.TWOHUNDREDANDFIFTYSIX.ordinal()) {
            return getKernel(sigma, 0, 1905);
        } else if (sigma.ordinal() == SIGMA.TWOHUNDREDANDFIFTYSIX.ordinal()) {
            return getKernel(SIGMA.getValue(sigma));
        }
        
        throw new IllegalArgumentException("haven't implemented a method for " 
        + " sigma=" + sigma);
    }
    
    protected static float[] getBinomialKernel(SIGMA sigma) {
      
        //using the binomial filter removes artifact sometimes introduced in
        // convolution due to nearly zero elements in the kernel and subsequent
        // renormalization by a very small number, so using binomial when
        // possible
        
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
            return getBinomialKernelSigmaTwo();
        } else if (sigma.ordinal() == SIGMA.ONEPOINTFIVESQRT2.ordinal()) {
            return getBinomialKernelSigmaOnePointFiveSQRT2();
        } else if (sigma.ordinal() == SIGMA.TWOSQRT2.ordinal()) {
            return getBinomialKernelSigmaTwoSQRT2();
        } else if (sigma.ordinal() == SIGMA.THREE.ordinal()) {
            return getBinomialKernelSigmaThree();
        } else if (sigma.ordinal() == SIGMA.FOUR.ordinal()) {
            return getBinomialKernelSigmaFour();
        } else {
            return null;
        }
        
    }
    
    /*<pre>
    binomial filter
   
      n=                                                                  sigma=sqrt(1)/2 = 0.5
      n=2                                   1      0   -1                 sigma=sqrt(2)/2 = 0.707
      n=3                                 1    1     -1  -1               
      n=4                               1    2     0   -2   -1            <=== sigma = 1   sigma=sqrt(n)/2
      n=5                             1    3    2    -2   -3   -1         
      n=6                           1    4    5    0    -5   -4   -1 
      n=7                         1    5    9   5    -5   -9    -5   -1   
      n=8                       1   6    14   14   0    -14  -14    -6   -1     <=== sigma=sqrt(2)
      n=9                     1   7   20   28   14   -14   -28   -20   -7   -1  <=== sigma=3/2
      n=10                  1   8   27   48   42   0    -42   -48    -27   -8  -1  
     ...
    </pre>
    */
  
    protected static float[] getBinomialKernelSigmaZeroPointFive() {
        //1 -1 norm=0.5
        return new float[]{0.5f, 0.0f, -0.5f};
    }
   
    protected static float[] getBinomialKernelSigmaZeroPointSevenOne() {
        // 1    1   -1  -1, n=2, sigma=sqrt(2)/2 = 0.707
        // edit to:
        float a = 1.f/200.f;
        return new float[]{0, 9*a, 83*a, 0, -83*a, -9*a, 0};
    }
   
    protected static float[] getBinomialKernelSigmaOne() {
        // 1       3     2        -2    -3     -1  norm=1/12 = 0.08333333333333333
        // editing to:
        float a = 1.f/16.f;
        return new float[]{0, 2*a, 4*a, 0, -4*a, -2*a, 0};
    }
 
    protected static float[] getKernelSigmaOneSQRT2() {
        return new float[]{0, 0.0013614271f, 0.010333488f, 
            0.04459885f, 0.103776865f, 0.10984783f, 
            -0.0f, -0.10984783f, -0.103776865f, -0.04459885f, -0.010333488f, 
            -0.0013614271f, 0};
    }
    
    protected static float[] getBinomialKernelSigmaOneSQRT2() {
        //norm=256, n=8
        /*
        n=7         1   6    14   14   0    -14  -14    -6   -1
        n=8       1   7   20   28   14   -14   -28   -20   -7   -1    <=== sigma=sqrt(2)
        n=9     1   8   27   48   42   0    -42   -48    -27   -8  -1  <=== sigma=3/2
        edit to:
        */
        float a = 1.f/256.f;
        return new float[]{0, a*3, a*12, a*27, a*28,
            0, -28*a, -27*a, -12*a, -3*a, 0};
    }
    
    protected static float[] getBinomialKernelSigmaOnePointFive() {
        float a = 1.f/512.f;
        return new float[]{a*1, a*8, a*27, a*48, 
            a*42, 
            0, -1*a*42, -1*a*48, -1*a*27, 
            -1*a*8, -1*a*1};
    }

    protected static float[] getKernelSigmaTwo() {
        return new float[]{0, 0, 
            0.0033238865f, 0.010955188f, 0.026995484f, 0.048569098f, 
            0.060492683f, 0.044008166f, -0.0f, -0.044008166f, -0.060492683f, 
            -0.048569098f, -0.026995484f, -0.010955188f, -0.0033238865f, 
            0, 0};
    }
    
    protected static float[] getBinomialKernelSigmaTwo() {
        //sigma=2 n=16   appears to be too wide, so the real 'sigma' for this is > 2
        double a = 1.f/33300l;
        return new float[]{
            (float) (1l * a), (float) (14l * a), (float) (90l * a),
            (float) (350l * a), (float) (910l * a), (float) (1638l * a),
            (float) (2002l * a), (float) (1430l * a), (float) (0l * a),
            (float) (-1430l * a), (float) (-2002l * a), (float) (-1638l * a),
            (float) (-910l * a), (float) (-350l * a), (float) (-90l * a),
            (float) (-14l * a), (float) (-1l * a)
        };
    }
    
    protected static float[] getBinomialKernelSigmaOnePointFiveSQRT2() {
        //sigma=3*sqrt(2)/2   n=18
        double a = 1./133000l;
        return new float[]{
            (float) (1l * a), (float) (16l * a), (float) (119l * a),
            (float) (544l * a), (float) (1700l * a), (float) (3808l * a),
            (float) (6188l * a), (float) (7072l * a), (float) (4862l * a),
            (float) (0l * a), (float) (-4862l * a), (float) (-7072l * a),
            (float) (-6188l * a), (float) (-3808l * a), (float) (-1700l * a),
            (float) (-544l * a), (float) (-119l * a), (float) (-16l * a),
            (float) (-1l * a)
        };
    }
   
    protected static float[] getKernelSigmaTwoSQRT2() {
        return new float[]{0, 0, 
            0.0010043882f, 0.002583372f, 0.005772264f, 0.011149713f, 
            0.018478211f, 0.025944216f, 0.030137392f, 0.027461957f, 
            0.01656272f, -0.0f, -0.01656272f, -0.027461957f, -0.030137392f, 
            -0.025944216f, -0.018478211f, -0.011149713f, -0.005772264f, 
            -0.002583372f, -0.0010043882f, 0, 0};
    }
    
    protected static float[] getBinomialKernelSigmaTwoSQRT2() {
        // sigma=2*sqrt(2)  n=32
        double a = 1./2153976500l;
        return new float[]{
            (float) (1l * a), (float) (30l * a), (float) (434l * a),
            (float) (4030l * a), (float) (26970l * a), (float) (138446l * a),
            (float) (566370l * a), (float) (1893294l * a), (float) (5259150l * a),
            (float) (12271350l * a), (float) (24192090l * a), (float) (40320150l * a),
            (float) (56448210l * a), (float) (65132550l * a), (float) (58929450l * a),
            (float) (35357670l * a), (float) (0l * a), (float) (-35357670l * a),
            (float) (-58929450l * a), (float) (-65132550l * a), (float) (-56448210l * a),
            (float) (-40320150l * a), (float) (-24192090l * a), (float) (-12271350l * a),
            (float) (-5259150l * a), (float) (-1893294l * a), (float) (-566370l * a),
            (float) (-138446l * a), (float) (-26970l * a), (float) (-4030l * a),
            (float) (-434l * a), (float) (-30l * a), (float) (-1l * a)
        };
    }
   
    protected static float[] getKernelSigmaThree() {
        return new float[]{0, 0, 0, 0.0014772827f, 0.0033765885f, 0.006798267f, 
            0.011997993f, 0.018421693f, 0.024297789f, 0.026885636f, 
            0.023662815f, 0.013977156f, -0.0f, -0.013977156f, -0.023662815f, 
            -0.026885636f, -0.024297789f, -0.018421693f, -0.011997993f, 
            -0.006798267f, -0.0033765885f, -0.0014772827f, 0, 0, 0};
    }
    
    protected static float[] getBinomialKernelSigmaThree() {
        //sigma=3, n=36
        double a = 1./34500000000l;
        return new float[]{
            (float) (1l * a), (float) (34l * a), (float) (560l * a),
            (float) (5950l * a), (float) (45815l * a), (float) (272272l * a),
            (float) (1298528l * a), (float) (5101360l * a), (float) (16811300l * a),
            (float) (47071640l * a), (float) (112971936l * a), (float) (233646504l * a),
            (float) (417225900l * a), (float) (641886000l * a), (float) (843621600l * a),
            (float) (927983760l * a), (float) (811985790l * a), (float) (477638700l * a),
            (float) (0l * a), (float) (-477638700l * a), (float) (-811985790l * a),
            (float) (-927983760l * a), (float) (-843621600l * a), (float) (-641886000l * a),
            (float) (-417225900l * a), (float) (-233646504l * a), (float) (-112971936l * a),
            (float) (-47071640l * a), (float) (-16811300l * a), (float) (-5101360l * a),
            (float) (-1298528l * a), (float) (-272272l * a), (float) (-45815l * a),
            (float) (-5950l * a), (float) (-560l * a), (float) (-34l * a),
            (float) (-1l * a)
        };
    }
    
    protected static float[] getKernelSigmaFour() {
        return new float[]{0, 0, 0, 0, 0.0015629561f, 
            0.002738797f, 0.0044633886f, 0.006748871f, 0.009436582f, 
            0.0121422745f, 0.014269461f, 0.015123171f, 0.014115817f, 
            0.011002041f, 0.0060416893f, -0.0f, -0.0060416893f, -0.011002041f, 
            -0.014115817f, -0.015123171f, -0.014269461f, -0.0121422745f, 
            -0.009436582f, -0.006748871f, -0.0044633886f, -0.002738797f, 
            -0.0015629561f, 0, 0, 0, 0};
    }
    
    /**
     * curvature tests for the non-binomial kernel show better results so
     * prefer that for sigma=4.
     * 
     * @Deprecated 
     * @return 
     */
    protected static float[] getBinomialKernelSigmaFour() {
        // sigma=4  n=64
        // NOTE: the curvature tests for the non-binomial kernal sigma four 
        // are better, so use those.
        float a = 1.f/82.0f;
        return new float[]{
            (float) (1.4677788298671463E-5) * a, (float) (5.577559553495156E-5) * a,
            (float) (1.924887094025772E-4) * a, (float) (6.059829740451504E-4) * a,
            (float) (0.0017466568075419043) * a, (float) (0.004623503314081511) * a,
            (float) (0.011267648817280127) * a, (float) (0.02533103004787036) * a,
            (float) (0.05261060086865382) * a, (float) (0.10104575722392242) * a,
            (float) (0.17954411407556464) * a, (float) (0.2950768483502758) * a,
            (float) (0.4480796586059744) * a, (float) (0.627311522048364) * a,
            (float) (0.8065433854907539) * a, (float) (0.9459459459459458) * a,
            (float) (1.0) * a, (float) (0.9310344827586206) * a, (float) (0.7241379310344828) * a,
            (float) (0.3971078976640712) * a, (float) (0.0) * a, (float) (-0.3971078976640712) * a,
            (float) (-0.7241379310344828) * a, (float) (-0.9310344827586206) * a,
            (float) (-1.0) * a, (float) (-0.9459459459459458) * a, (float) (-0.8065433854907539) * a,
            (float) (-0.627311522048364) * a, (float) (-0.4480796586059744) * a,
            (float) (-0.2950768483502758) * a, (float) (-0.17954411407556464) * a,
            (float) (-0.10104575722392242) * a, (float) (-0.05261060086865382) * a,
            (float) (-0.02533103004787036) * a, (float) (-0.011267648817280127) * a,
            (float) (-0.004623503314081511) * a, (float) (-0.0017466568075419043) * a,
            (float) (-6.059829740451504E-4) * a, (float) (-1.924887094025772E-4) * a,
            (float) (-5.577559553495156E-5) * a, (float) (-1.4677788298671463E-5) * a
        };
    }
    
}
