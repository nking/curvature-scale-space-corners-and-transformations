package algorithms.imageProcessing;

import static algorithms.imageProcessing.Gaussian1D.estimateHWZI;
import algorithms.misc.MiscMath;

/**
 * class to retrieve a kernel for the first derivative of a gaussian of
 * a given sigma.
 * 
 * For the smaller sigma, curvature tests show that a binomial filter should
 * be used instead of the detailed calculation of the formula due to errors
 * in normalization due to the tails of the function.
 * Therefore, the general methods are returning binomial filter curves
 * when possible.  Sometimes those have been edited from the binomial pattern
 * due to centering.
 * 
 * @author nichole
 */
public class Gaussian1DSecondDeriv {
        
    // mu is t, the parameterization scale axis
    //
    //               1             ( -(x - mu)^2 )
    // f = ------------------ * exp( ----------- )
    //     sigma * sqrt(2*pi)      (    2o~^2    )
    //
    //           -(x - mu)              ( -(x - mu)^2 )
    // dfdx = -------------------- * exp( ----------- )
    //        sigma^3 * sqrt(2*pi)      (    2o~^2    )
    
    //
    //
    // Lowe's paper has this for first factors:  
    //      (x - mu)^2 - sigma^2
    //    ----------------------    
    //     sigma^5 * sqrt(2*pi)
    //
    //            (x - mu - sigma) * (x - mu + sigma)      ( -(x - mu)^2 )
    // d^2fdx^2 = ----------------------------------- * exp( ----------- )
    //                   sigma^5 * sqrt(2*pi)              (    2o~^2    )
    //
  
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
       <pre>
                  (x - mu - sigma) * (x - mu + sigma)      ( -(x - mu)^2 )
       d^2fdx^2 = ----------------------------------- * exp( ----------- )
                         sigma^5 * sqrt(2*pi)              (    2o~^2    )
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
     * @param mu
     * @param nPoints
     * @return 
     */
    public static float[] getKernel(float sigma, float mu, int nPoints) {

        float normalization = (float)(Math.pow(sigma, 5) * 
            Math.sqrt(2.f * Math.PI));
        
        if (sigma < 1) {
            // hack to correct for the power introduced by (x-sigma)*(x+sigma)
            int power = MiscMath.findPowerOf10(sigma*sigma);
            normalization *= (1./(sigma*sigma))*Math.pow(10, -1*power);
        }
        
        int halfWidthInPixels = nPoints >> 1;
        
        int start = (int) (mu - halfWidthInPixels);
        int stopExcl = (int)mu + halfWidthInPixels + 1;
            
        float d, dsq;
                
        int count = 0;
        float[] yPoints = new float[nPoints];
        for (int i = start; i < stopExcl; i++) {
                         
            if (count >= nPoints) {
                break;
            }
            
            float x = i;
            
            d = (x - mu);
            dsq = d*d;
            //(x - mu)^2 - sigma^2
            // OR (d - sigma) * (d + sigma)
            float y = (dsq - sigma*sigma) * (float)(
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
       <pre>
                  (x - mu - sigma) * (x - mu + sigma)      ( -(x - mu)^2 )
       d^2fdx^2 = ----------------------------------- * exp( ----------- )
                         sigma^5 * sqrt(2*pi)              (    2o~^2    )
       </pre>

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
    public static float[] getKernel(SIGMA sigma, float mu) {
        
        return getKernel(SIGMA.getValue(sigma), mu);
    }
    
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
     * get the kernel for the given sigma. mu = 0 is used.
     * 
     * the function is:
       <pre>
                  (x - mu - sigma) * (x - mu + sigma)      ( -(x - mu)^2 )
       d^2fdx^2 = ----------------------------------- * exp( ----------- )
                         sigma^5 * sqrt(2*pi)              (    2o~^2    )
       </pre>

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
     * get the kernel for the given sigma and deltaX. 
     * NOTE: working out whether an offset can be applied to results...
     * has not been applied yet to allow re-use of the kernel, so when using this
     * kernel, be sure to multiply each element by -(x-mu);
     * 
     * the function is:
     <pre>
                  (x - mu - sigma) * (x - mu + sigma)      ( -(x - mu)^2 )
       d^2fdx^2 = ----------------------------------- * exp( ----------- )
                         sigma^5 * sqrt(2*pi)              (    2o~^2    )
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
            return getKernel(SIGMA.getValue(sigma));
        } else if (sigma.ordinal() == SIGMA.FOURSQRT2.ordinal()) {
            return getKernel(SIGMA.getValue(sigma));
        } else if (sigma.ordinal() == SIGMA.EIGHT.ordinal()) {
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
    binomial filter for the second derivative:
   
      n=0? extrapolated            1 -2 1                                    
      n=1? extrapolated          1  -1 -1  1                sigma=sqrt(1)/2 = 0.5
      n=2? extrapolated       1   0  -2  0    1             sigma=sqrt(2)/2 = 0.707
      n=3                   1   1   -2  -2   1   1          sigma=sqrt(3)/2 = 0.866
      n=4                 1   2   -1  -4  -1   2   1        <=== sigma = 1   sigma=sqrt(n)/2
      n=5               1   3   1   -5  -5   1   3    1 
      n=6             1   4   4   -4  -10  -4   4   4    1
      n=7 extr      1   5   8   0  -14   -14    0   8    5    1
      n=8 extr    1   6  13   8   -14  -28   -14   8   13   6   1     <=== sigma=sqrt(2)
      n=9 extr  1  7   19  21   -6  -42  -42    -6   21   19   7    1  <=== sigma=3/2
    
      n=10                 1   8  26  40   15  -48  -84  -48  15  40  26   8   1
      n=11               1   9  34  66   55  -33 -132  -132 -33  55  66  34   9   1
      n=12             1  10  43  100 121  22  -165  -264  -165  22  121  100   43   10   1
      n=13           1  11  53  143 221 143 -143  -429 -429 -143  143  221  143   53   11   1
      n=14         1 12  64  196  364 364  0   -572 -858  -572  0   364 364  196   64  12  1
      n=15       1 13  76  260  560 728  364 -572 -1430  -1430 -572  364  728  560 260 76 13 1
      n=16   * 1 14  89  336 820 1288 1092 -208 -2002 -2860 -2002 -208 1092 1288 820 336 89 14 1 <== sigma=2
    </pre>
    */
    
    //TODO:  These should be changed to include zeros on the ends
    //        to complete the 2n deriv profile
  
    protected static float[] getBinomialKernelSigmaZeroPointFive() {
        // 1 -1  -1  1
        // edited:   1 -2  1  which is condensing n=4...
        return new float[]{0.037f, -0.074f, 0.037f};
    }
    
    protected static float[] getKernelSigmaZeroPointSevenOne() {
        return new float[]{0, 0.0072334423f, 0.020755375f, 
            -0.05641896f, 0.020755375f, 0.0072334423f, 0};
    }
    
    protected static float[] getBinomialKernelSigmaZeroPointSevenOne() {
        // extrapolated isn't correct: 1   0  -2  0    1
        // but condensing n=5 looks correct:  1   3  -8   3    1
        float a = 1.f/141.8f;
        return new float[]{a, 3*a, -8*a, 3*a, a};
    }
    
    protected static float[] getKernelSigmaOne() {
        return new float[]{0.0020074535f, 0.035454787f, 0.16197291f, -0.0f, 
            -0.3989423f, 0.0f, 0.16197291f, 0.035454787f, 0.0020074535f};
    }
       
    protected static float[] getBinomialKernelSigmaOne() {
        // 1   2   -1  -4  -1   2   1 , norm=1./9.
        // editing to:
        //float a = 1.f/30.f;
        //return new float[]{a, 5*a, 0, -12*a, 0, 5*a, a};
        return new float[]{1, 2, -1, -4, -1, 2, 1};
    }
    
    protected static float[] getBinomialKernelSigmaOneSQRT2() {
        // norm=1./?, n=8, sigma=sqrt(2)
        // editing to:
        float a = 1.f/256.f;
        return new float[]{1*a, 5*a, 13*a, 13*a, -14*a, 
            -36*a, -14*a, 13*a, 13*a, 5*a, 1*a
        };
    }
    
    protected static float[] getBinomialKernelSigmaOnePointFive() {
          
        // best results are for close to n=8 rather than n=9 extrapolation!!
        
        float a = 1.f/256.f;
        return new float[]{1*a, 5*a, 12*a, 10*a, -13*a, -30*a,
            -13*a, 10*a, 12*a, 5*a, 1*a
        };
    }
    
    protected static float[] getBinomialKernelSigmaTwo() {
        //sigma=2 for n=16
        float a = 1.f/64000.f;
        return new float[]{1*a, 14*a, 89*a, 336*a, 820*a, 1288*a, 1092*a,
            -208*a, -2002*a, 
            -2860*a, -2002*a, -208*a, 1092*a, a*1288, a*820, a*336, a*89, 
            a*14, a};
    }
    
    protected static float[] getBinomialKernelSigmaOnePointFiveSQRT2() {
        //sigma=3*sqrt(2)/2    n=18
        double a = 1./245000l;
        return new float[]{
            (float) (1l * a), (float) (16l * a), (float) (118l * a),
            (float) (528l * a), (float) (1581l * a), (float) (3264l * a),
            (float) (4488l * a), (float) (3264l * a), (float) (-1326l * a),
            (float) (-7072l * a), (float) (-9724l * a), (float) (-7072l * a),
            (float) (-1326l * a), (float) (3264l * a), (float) (4488l * a),
            (float) (3264l * a), (float) (1581l * a), (float) (528l * a),
            (float) (118l * a), (float) (16l * a), (float) (1l * a)
        };
    }
  
    protected static float[] getBinomialKernelSigmaTwoSQRT2() {
        //sigma=2*sqrt(2)    n=32
        double a = 1./4100000000l;
        return new float[]{
            (float) (1l * a), (float) (30l * a), (float) (433l * a),
            (float) (4000l * a), (float) (26536l * a), (float) (134416l * a),
            (float) (539400l * a), (float) (1754848l * a), (float) (4692780l * a),
            (float) (10378056l * a), (float) (18932940l * a), (float) (28048800l * a),
            (float) (32256120l * a), (float) (24812400l * a), (float) (2481240l * a),
            (float) (-29774880l * a), (float) (-58929450l * a), (float) (-70715340l * a),
            (float) (-58929450l * a), (float) (-29774880l * a), (float) (2481240l * a),
            (float) (24812400l * a), (float) (32256120l * a), (float) (28048800l * a),
            (float) (18932940l * a), (float) (10378056l * a), (float) (4692780l * a),
            (float) (1754848l * a), (float) (539400l * a), (float) (134416l * a),
            (float) (26536l * a), (float) (4000l * a), (float) (433l * a),
            (float) (30l * a), (float) (1l * a)
        };
    }
   
    protected static float[] getBinomialKernelSigmaThree() {
        //sigma=3    n=36
        double a = 1./65000000000l;
        return new float[]{
            (float) (1l * a), (float) (34l * a), (float) (559l * a),
            (float) (5916l * a), (float) (45255l * a), (float) (266322l * a),
            (float) (1252713l * a), (float) (4829088l * a), (float) (15512772l * a),
            (float) (41970280l * a), (float) (96160636l * a), (float) (186574864l * a),
            (float) (304253964l * a), (float) (408239496l * a), (float) (426395700l * a),
            (float) (286097760l * a), (float) (-31635810l * a), (float) (-450345060l * a),
            (float) (-811985790l * a), (float) (-955277400l * a), (float) (-811985790l * a),
            (float) (-450345060l * a), (float) (-31635810l * a), (float) (286097760l * a),
            (float) (426395700l * a), (float) (408239496l * a), (float) (304253964l * a),
            (float) (186574864l * a), (float) (96160636l * a), (float) (41970280l * a),
            (float) (15512772l * a), (float) (4829088l * a), (float) (1252713l * a),
            (float) (266322l * a), (float) (45255l * a), (float) (5916l * a),
            (float) (559l * a), (float) (34l * a), (float) (1l * a)
        };
    }
    
    protected static float[] getBinomialKernelSigmaFour() {
        //sigma=4    n=64
        //1./300000000000
        float a = 1.f/315.f;
        return new float[]{
             (float)(0.0013496157933984045)*a, (float)(0.003812253148120798)*a,
             (float)(0.00985466410118652)*a, (float)(0.023354250937084688)*a,
             (float)(0.050793931240985246)*a, (float)(0.10141062660346886)*a,
             (float)(0.18571762697885447)*a, (float)(0.31132423871622195)*a,
             (float)(0.4757716476259789)*a, (float)(0.6580431670807847)*a,
             (float)(0.8130047466820091)*a, (float)(0.8746801742701372)*a,
             (float)(0.7727267283900054)*a, (float)(0.4606218146694771)*a,
             (float)(0.05400657042753325)*a,                       (float)(-0.6908945681321693)*a,
             (float)(-1.2995036536356718)*a, (float)(-1.682898921828288)*a,
             
             (float)(-1.9582123932012366)*a, (float)(-1.682898921828288)*a,
             (float)(-1.2995036536356718)*a, (float)(-0.6908945681321693)*a,
             (float)(0.05400657042753325)*a, (float)(0.4606218146694771)*a,
             (float)(0.7727267283900054)*a, (float)(0.8746801742701372)*a, (float)(0.8130047466820091)*a,
             (float)(0.6580431670807847)*a, (float)(0.4757716476259789)*a,
             (float)(0.31132423871622195)*a, (float)(0.18577391691321074)*a,
             (float)(0.10141612056970953)*a, (float)(0.050794238557124116)*a,
             (float)(0.023354255369729467)*a, (float)(0.00985466410118652)*a,
             (float)(0.003812253148120798)*a, (float)(0.0013496157933984045)*a
        };
    }
    
}
