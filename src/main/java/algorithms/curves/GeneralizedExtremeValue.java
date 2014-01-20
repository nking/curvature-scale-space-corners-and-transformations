package algorithms.curves;

import algorithms.misc.MiscMath;

/**
  <pre>
  Generate curves following the Generalized Extreme Value probability density
  function.
 
                           (   (      ( x-mu))-(1/k))
                           (-1*(1 + k*(-----))      )
                  1        (   (      (sigma))      )   (      ( x-mu))(-1-(1/k))
  y = y_const * ----- * exp                           * (1 + k*(-----))
                sigma                                   (      (sigma))
 
  mu is  the location parameter
  sigma is the scale parameter and is > 0
  k is the shape parameter
 
 
  if k != 0,
      1 + (k * (x-mu)/sigma) > 0
 
 
  Let z = (x-mu)/sigma
 
  Extreme Value Type I:
                      1
      y = y_const * ----- * exp( -z -exp(-z))
                    sigma
      
      sigma > 0
      k = 0
 
  Extreme Value Type II:
                      k     (sigma)^(k+1)      (   (sigma)^k)
      y = y_const * ----- * (-----)       * exp( - (-----)  )
                    sigma   (  x  )            (   (  x  )  )
 
      k  > 0
      sigma > 0
      x > 0
 
  Extreme Value Type III:
                      k     (x - mu)^(k+1)      (   (x - mu)^k)
      y = y_const * ----- * (------)       * exp( - (------)  )
                    sigma   (sigma )            (   (sigma )  )
 
                      k
        = y_const * ----- * (z)^(k+1) * exp( -1*(z)^k )
                    sigma
 
      k < 0
      sigma > 0
      x > 0
 
 </pre>
 * @author nichole
 */
public class GeneralizedExtremeValue implements ICurveGenerator {

    protected final float[] x;
    protected final float[] y;

    protected final float[] dx;
    protected final float[] dy;

    public GeneralizedExtremeValue(float[] xPoints, float[] yPoints, float[] dXPoints, float[] dYPoints) {
        this.x = xPoints;
        this.y = yPoints;
        this.dx = dXPoints;
        this.dy = dYPoints;
    }

    public float[] generateNormalizedCurve(float[] parameters, float yConst) {
        if (parameters == null) {
            throw new IllegalArgumentException("parameters cannot be null");
        }
        if (parameters.length != 3) {
            throw new IllegalArgumentException("parameters must hold k, sigma, and mu");
        }

        return generateNormalizedCurve(parameters[0], parameters[1], parameters[2], yConst);
    }

    public float[] generateNormalizedCurve(float[] parameters) {
        if (parameters == null) {
            throw new IllegalArgumentException("parameters cannot be null");
        }
        if (parameters.length != 3) {
            throw new IllegalArgumentException("parameters must hold k, sigma, and mu");
        }

        return generateNormalizedCurve(parameters[0], parameters[1], parameters[2]);
    }

    public float[] generateNormalizedCurve(float k, float sigma, float mu) {

        float[] yGEV = generateCurve(x, k, sigma, mu);
        if (yGEV == null) {
            return null;
        }
        float yConst = determineYConstant(yGEV, mu);

        return generateNormalizedCurve(k, sigma, mu, yConst);
    }

    /**
     * generate a GEV curve w/ parameters k, sigma, and mu and then multiply it
     * by yConst so that the highest peak has value yConst
     *
     * @param k
     * @param sigma
     * @param mu
     * @param yConst
     * @return
     */
    public float[] generateNormalizedCurve(float k, float sigma, float mu, float yConst) {

        float[] yGEV = generateCurve(x, k, sigma, mu);
        if (yGEV == null) {
            return null;
        }

        float yMax = MiscMath.findMax(yGEV);

        for (int i = 0; i < yGEV.length; i++) {
            yGEV[i] *= yConst/yMax;
        }
        return yGEV;
    }

    /* <pre>
     *                          (   (      ( x-mu))-(1/k))
     *                          (-1*(1 + k*(-----))      )
     *                 1        (   (      (sigma))      )   (      ( x-mu))(-1-(1/k))
     * y = y_const * ----- * exp                           * (1 + k*(-----))
     *               sigma                                   (      (sigma))
     *
     * mu is  the location parameter
     * sigma is the scale parameter and is > 0
     * k is the shape parameter
     *
     *
     * Components needed in the derivatives:
     *
     *   Let z = (1 + k*( (x-mu)/sigma )
     *
     *   f1 = exp( -1. * ( z^(-1./k) ) )
     *
     *   f2 = z^(-1. - (1/k))
     *
     *   then y = (yconst/sigma) * f1 * f2
     *</pre>
     */
    public static Double generateYGEV(float xPoint, float k, float sigma, float mu) {

        if (sigma == 0) {
            //throw new IllegalArgumentException("sigma must be > 0");
            return null;
        }

        // if k =~ 0, use Gumbel which is Type I
        if (k == 0 || Float.isInfinite(1.f/k)) {
            return generateYEVTypeI(xPoint, sigma, mu);
        }        
        // if k > 0, use Frechet which is Type II
        // if k < 0, use Weibull which is Type III.  Not using k < 0 in this project

        float z = 1.f + k * ((xPoint - mu)/sigma);
        float a,b;

        boolean zIsNegative = (z < 0);
        // When z is negative, need to use alternative methods for exponentiation:
        // 
        // For z^(g/h)
        //
        // For the continuous real exponentiation operator, negative base isn't allowed
        // For the discrete real exponentiation operator,
        //    fractional exponents with odd denominators are allowed.
        //    (-z)^(g/h) = ((-z)^g)^(1/h) = ((-1)^g)(z^(g/h))
        // For the complex exponentiation operator, that is complex bases, 
        //    the results are not continuous and are infinite.
        //    the principal value is
        //        (-1)^(g/h) * ((1./z)^(g/h))

        if (zIsNegative) {
            float invNegZ = -1.0f*(1.0f/z);
            float neg1Pow = -1.0f; // TODO:  revisit this
            a = -1.f * neg1Pow * (float) Math.pow(invNegZ, (-1.f/k));
            b = neg1Pow * (float) Math.pow(invNegZ, (-1.f - (1.f/k)));
        } else {
            a = -1.f * (float) Math.pow(z, (-1.f/k));
            b = (float) Math.pow(z, (-1.f - (1.f/k)));
        }
        
        if (Float.isNaN(a) || Float.isNaN(b)) {
            // or return 0?
            return null;
        } else {
            float t = (float) ((1.f/sigma) * Math.exp(a) * b);
            return Double.valueOf(t);
        }
    }

    public static Double generateYEVTypeI(float xPoint, float sigma, float mu) {

        if (sigma == 0) {
            throw new IllegalArgumentException("sigma must be > 0");
        }

        double z = (xPoint - mu)/sigma;

        double a = (float) (-1.f*z - Math.exp(-1.0f*z));

        double yGEV = (float) ((1.f/sigma) * Math.exp(a));

        return yGEV;
    }

    public float[] generateCurve(float[] x1, float k, float sigma, float mu) {

        if (sigma == 0) {
            return null;
        }
        // if k =~ 0, use Gumbel which is Type I
        if (k == 0 || Float.isInfinite(1.f/k)) {
            return generateEVTypeICurve(x1, sigma, mu);
        }
        
        float[] yGEV = new float[x1.length];

        for (int i = 0; i < x1.length; i++) {

            float z = 1.f + k*((x1[i] - mu)/sigma);
            float a,b;

            boolean zIsNegative = (z < 0);
            // When z is negative, need to use alternative methods for exponentiation:
            //
            // For z^(g/h)
            //
            // For the continuous real exponentiation operator, negative base isn't allowed
            // For the discrete real exponentiation operator,
            //    fractional exponents with odd denominators are allowed.
            //    (-z)^(g/h) = ((-z)^g)^(1/h) = ((-1)^g)(z^(g/h))
            // For the complex exponentiation operator, that is complex bases,
            //    the results are not continuous and are infinite.
            //    the principal value is
            //        (-1)^(g/h) * ((1./z)^(g/h))
            
            if (zIsNegative) {
                float invNegZ = -1.0f*(1.0f/z);
                float neg1Pow = -1.0f; // TODO:  revisit this
                a = -1.f * neg1Pow * (float) Math.pow(invNegZ, (-1.f/k));
                b = neg1Pow * (float) Math.pow(invNegZ, (-1.f - (1.f/k)));
            } else {
                a = -1.f * (float) Math.pow(z, (-1.f/k));
                b = (float) Math.pow(z, (-1.f - (1.f/k)));
            }

            if (Float.isNaN(a) || Float.isNaN(b)) {
                yGEV[i] = 0;
            } else {
                float t = (float) ((1.f/sigma) * Math.exp(a) * b);
                if (t < 0 || Float.isNaN(t)) {
                    yGEV[i] = 0;
                } else {
                    yGEV[i] = t;
                }
            }
        }

        return yGEV;
    }

    public static float[] generateEVTypeICurve(float[] x1, float sigma, float mu) {

        if (sigma == 0) {
            throw new IllegalArgumentException("sigma must be > 0");
        }

        float[] yGEV = new float[x1.length];

        for (int i = 0; i < x1.length; i++) {

            float z = (x1[i] - mu)/sigma;

            float a = (float) (-1.f*z - Math.exp(-1.0f*z));

            yGEV[i] = (float) ((1.f/sigma) * Math.exp(a));
        }

        return yGEV;
    }

    public float determineYConstant(float[] yGEV, float mu) {

        int index = MiscMath.findYMaxIndex(y);
        if (index == -1) {
            return 0;
        }

        return determineYConstant(yGEV, mu, index);
    }

    public float determineYConstant(float[] yGEV, float mu, int index) {

        float yConst = y[index]/yGEV[index];

        return yConst;
    }

    /**
     * Errors in fitting the curve can be calculated via chain rule of derivatives of the fit:
     *
     *
     * @param y1 the y data array
     * @param yGEV the generated GEV y model array
     * @return the mean error of the fit
     */
    public static float calculateChiSq(float[] y1, float[] yGEV) {

        float yNorm = MiscMath.findMax(y1);

        float chiSum = 0.f;

        for (int i = 0; i < yGEV.length; i++) {

            float z = (yGEV[i] - y1[i])/yNorm;
            z *= z;

            chiSum += z;
        }

        return chiSum;
    }

    /**
     <pre>
     return the error in fitting the GEV curve

                                 | dy_fit |^2            | dy_fit |^2            | dy_fit|^2            |dy_fit|^2
       (err_y_fit)^2 =  (err_x)^2|--------|   + (err_k)^2|--------|   + (err_s)^2|-------|   + (err_m)^2|------|
                                 |   dx   |              |   dk   |              |  ds   |              |  dm  |
     </pre>
     * @param bestFit2
     * @return
     */
    public static double calculateFittingErrorSquared(GEVYFit yFit, float xPoint) {

        if (yFit == null) {
            return Float.POSITIVE_INFINITY;
        }

        int yMaxIdx = MiscMath.findYMaxIndex(yFit.getYFit());

        if (yMaxIdx == -1) {
            return Float.POSITIVE_INFINITY;
        }

        float xPeak = yFit.getX()[yMaxIdx];

        // since we don't have an error in the parameters, we'll make a rough
        //   guess with the parameters in bestFit compared to prior bestFit.
        //  the error in the parameters, k, sigma, and mu are assumed to be smaller
        //  than those deltas in order for the best fit to have a better fit.
        //  we can use the full delta to overestimate the error, or assume that
        //  the error has to be smaller than half of the delta otherwise the prior Fit
        //  value would have been kept.  very very rough approx for parameter errors...

        float kDelta = yFit.getKResolution()/2.f;
        float sigmaDelta = yFit.getSigmaResolution()/2.f;
        float muDelta = yFit.getMuSolutionResolution()/2.f;

        float yConst = yFit.getYScale();
        float mu = yFit.getMu();
        float k = yFit.getK();
        float sigma = yFit.getSigma();

        float xDelta = yFit.getX()[1] - yFit.getX()[0];

        double dydx = DerivGEV.derivWRTX(yConst, mu, k, sigma, xPoint);

        double dydk = DerivGEV.derivWRTK(yConst, mu, k, sigma, xPoint);

        double dyds = DerivGEV.derivWRTSigma(yConst, mu, k, sigma, xPoint);

        double dydm = DerivGEV.derivWRTMu(yConst, mu, k, sigma, xPoint);

        // since x is always given as a number, exclude it from propagation
        double err = /*Math.pow(xDelta*dydx, 2)*/ + Math.pow(kDelta*dydk, 2) + Math.pow(sigmaDelta*dyds, 2) + Math.pow(muDelta*dydm, 2);

//System.out.println("error in x alone: " + Math.sqrt( Math.pow(xDelta*dydx, 2) )/yFit.getYScale());

        return (float)err/(yFit.getYScale() * yFit.getYScale());
    }

    /**
      <pre>
      calculate the error in an area / height calculation for y=0 to y > yLimitFraction where y is
      yLimitFraction to to the right of the peak.  This is useful for determining errors for things
      like FWHM for example.
     
      For FWHM we have sum of f = sum(X_i*Y_i)_(i < yLimit)/ Y_i

          err^2 = xError^2*(Y_i/Y_i) = xError^2

          it reduces to the sum of the errors in x.  no pde's...
     </pre>
     
     * @param yFit
     * @param yLimitFraction
     * @return
     */
    public static double calculateWidthFittingError(GEVYFit yFit, float yMaxFactor) {

        if (yFit == null) {
            return Float.POSITIVE_INFINITY;
        }

        int yPeakIdx = MiscMath.findYMaxIndex(yFit.getOriginalScaleX());
        float yLimit = yMaxFactor * yFit.getOriginalScaleYFit()[yPeakIdx];
        int yLimitIdx = -1;
        for (int i = 0; i < yFit.getOriginalScaleX().length; i++) {
            if (i > yPeakIdx) {
                if (yFit.getOriginalScaleYFit()[i] > yLimit) {
                    yLimitIdx = i;
                } else {
                    break;
                }
            } else {
                yLimitIdx = i;
            }
        }

        // xError[i] should be formally calculated, but the approximation that it can be determined no
        //   better than  the bin center +- binwidth/2  is a minimum error.
        //   a safe addition to that (added in quadrature) would be an error in x derived from chi square
        //   but that is not done here

        float xDelta = yFit.getX()[1] - yFit.getX()[0];
        float xErrorSq = (xDelta*xDelta/4.f);

        float sum = 0.0f;
        for (int i = 0; i <= yLimitIdx; i++) {
            sum += xErrorSq;
        }

        sum = (float) Math.sqrt(sum);

        return sum;
    }

    public static float[] generateNormalizedCurve(float[] x1, float k, float sigma, float mu) {

        if (sigma == 0) {
            //throw new IllegalArgumentException("sigma must be > 0");
            return null;
        }

        float[] yGEV = genCurve(x1, k, sigma, mu);

        float yMax = MiscMath.findMax(yGEV);

        for (int i = 0; i < yGEV.length; i++) {
            yGEV[i] /= yMax;
        }

        return yGEV;
    }

    public static float[] genCurve(float[] x1, float k, float sigma, float mu) {

        if (sigma == 0) {
            return null;
        }
        // if k =~ 0, use Gumbel which is Type I
        if (k == 0 || Float.isInfinite(1.f/k)) {
            return generateEVTypeICurve(x1, sigma, mu);
        }

        float[] yGEV = new float[x1.length];

        for (int i = 0; i < x1.length; i++) {

            float z = 1.f + k * ((x1[i] - mu)/sigma);
            float a,b;

            boolean zIsNegative = (z < 0);
            // When z is negative, need to use alternative methods for exponentiation:
            //
            // For z^(g/h)
            //
            // For the continuous real exponentiation operator, negative base isn't allowed
            // For the discrete real exponentiation operator,
            //    fractional exponents with odd denominators are allowed.
            //    (-z)^(g/h) = ((-z)^g)^(1/h) = ((-1)^g)(z^(g/h))
            // For the complex exponentiation operator, that is complex bases,
            //    the results are not continuous and are infinite.
            //    the principal value is
            //        (-1)^(g/h) * ((1./z)^(g/h))

            if (zIsNegative) {
                float invNegZ = -1.0f*(1.0f/z);
                float neg1Pow = -1.0f; // TODO:  revisit this
                a = -1.f * neg1Pow * (float) Math.pow(invNegZ, (-1.f/k));
                b = neg1Pow * (float) Math.pow(invNegZ, (-1.f - (1.f/k)));
            } else {
                a = -1.f * (float) Math.pow(z, (-1.f/k));
                b = (float) Math.pow(z, (-1.f - (1.f/k)));
            }

            if (Float.isNaN(a) || Float.isNaN(b)) {
                yGEV[i] = 0;
            } else {
                float t = (float) ((1.f/sigma) * Math.exp(a) * b);
                if (t < 0 || Float.isNaN(t)) {
                    yGEV[i] = 0;
                } else {
                    yGEV[i] = t;
                }
            }
        }

        return yGEV;
    }
}
