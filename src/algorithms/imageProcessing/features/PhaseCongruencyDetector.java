package algorithms.imageProcessing.features;

import algorithms.QuickSort;
import algorithms.imageProcessing.AdaptiveThresholding;
import algorithms.imageProcessing.DistanceTransform;
import algorithms.imageProcessing.FilterGrid;
import algorithms.imageProcessing.FilterGrid.FilterGridProducts;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.LowPassFilter;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.MorphologicalFilter;
import algorithms.imageProcessing.NonMaximumSuppression;
import algorithms.imageProcessing.PeriodicFFT;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.Complex;
import algorithms.misc.ComplexModifiable;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import com.climbwithyourfeet.clustering.DTClusterFinder;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 An edge detector that uses principles of phase congruency to create an edge
 * map and orientation and phase angle images.
 * Phase congruency operates in the frequency domain of fourier transforms and
 * with the inverse FFT produces an image that is summed over scales and
 * cleaned of some of the noise.
 * The phase congruency refers to the overlapping of sine waves at same phases
 * in the frequency domain produced from feature edges in the spatial domain.
 * <pre>
 * The phase congruency method is thought to be better able to find edges
 * under varying illumination conditions.  It also has the characteristic of
 * producing a single response to an edge (in contrast to many spatial gradient
 * methods which for blurry edges, especially, produce a double response in
 * the gradient image).
 * </pre>
 * 
 * For further reading other than the references below, a summary of band-pass
 * quadrature filters is in
 * https://www.utc.fr/~dboukerr/Papers/Qf_JMIV_2004.pdf
 *
 * Listings of copyrights for the original source codes in languages Matlab and
 * python follow:
 *
 adapted from
  http://www.peterkovesi.com/matlabfns/PhaseCongruency/phasecongmono.m
  which has copyright:
  -------- begin Kovesi copyright --------
  Copyright (c) 1996-2013 Peter Kovesi
  Centre for Exploration Targeting
  The University of Western Australia
  peter.kovesi at uwa edu au

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  The Software is provided "as is", without warranty of any kind.
  -------- end Kovesi copyright --------
 
 * useful also was the python phasepack port by Alistair Muldal
 *  http://pydoc.net/Python/phasepack/1.4/phasepack.phasecongmono/
 * which has the following copyright:
 * # MIT License:
-------- begin phasepack copyright --------
# Permission is hereby  granted, free of charge, to any  person obtaining a
# copy of this software and associated  documentation files (the "Software"),
# to deal in the Software without restriction, subject to the following
# conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# The software is provided "as is", without warranty of any kind.
-------- end phasepack copyright --------
* 
% =========================================
*
% Notes on specifying parameters:
%
% The convolutions are done via the FFT.  Many of the parameters relate to the
% specification of the filters in the frequency plane.  The values do not seem
% to be very critical and the defaults are usually fine.  You may want to
% experiment with the values of 'nscales' and 'k', the noise compensation
% factor.
*
% Notes on filter settings to obtain even coverage of the spectrum
% sigmaOnf       .85   mult 1.3
% sigmaOnf       .75   mult 1.6     (filter bandwidth ~1 octave)
% sigmaOnf       .65   mult 2.1
% sigmaOnf       .55   mult 3       (filter bandwidth ~2 octaves)
%
% Note that better results are achieved using the large bandwidth filters.
% I generally use a sigmaOnf value of 0.55 or even smaller.
%
% References:
%
%     Peter Kovesi, "Image Features From Phase Congruency". Videre: A
%     Journal of Computer Vision Research. MIT Press. Volume 1, Number 3,
%     Summer 1999 http://mitpress.mit.edu/e-journals/Videre/001/v13.html
*
*     Michael Felsberg and Gerald Sommer, "A New Extension of Linear Signal
%     Processing for Estimating Local Properties and Detecting Features". DAGM
%     Symposium 2000, Kiel
%
%     Michael Felsberg and Gerald Sommer. "The Monogenic Signal" IEEE
%     Transactions on Signal Processing, 49(12):3136-3144, December 2001
%
%     Peter Kovesi, "Phase Congruency Detects Corners and Edges". Proceedings
%     DICTA 2003, Sydney Dec 10-12
* 
   Additions by Nichole King:
    -- use of an adaptive binary threshold in a 2-layer filter.
       That algorithm is in class AdaptiveThresholding.java and it is from
       "Efficient Implementation of Local Adaptive Thresholding
       Techniques Using Integral Images"
           by Shafaita, Keysersa, and Breuelb
    -- a flag 'extract noise' to extract the noise as a difference of the
       result of k factors and to extract a subset of the noise as candidates
       for finding textures.  It uses a  
       density based clustering algorithm
       http://nking.github.io/two-point-correlation/
       which has an MIT license
      ---- begin nking copyright ----
      The MIT License (MIT)
      Copyright (c) 2013-* Nichole King
      http://nking.github.io/two-point-correlation/

        Permission is hereby granted, free of charge, to any person obtaining 
        a copy of this software and associated documentation files 
        (the "Software"), to deal in the Software without restriction, 
        including without limitation the rights to use, copy, modify, merge, 
        publish, distribute, sublicense, and/or sell copies of the Software, 
        and to permit persons to whom the Software is furnished to do so, 
        subject to the following conditions:

        The above copyright notice and this permission notice shall be included 
        in all copies or substantial portions of the Software.
        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
        OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
        MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
        IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
        CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
        TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
        SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
     ---- end nking copyright ---- 
    
     runtime complexity of majority of algorithm is a constant factor, that
     depends upon parameters chosen, times O(N*lg_2(N)),
     but each operation uses a transcendental function.

     the advantage to using phase congruency is that it's possible to produce
     consistent clean edges with less noise and representing a wide range in
     intensity contrast.

 */
public class PhaseCongruencyDetector {

    final private static double epsilon = 1E-4;

    private boolean doPlot = false;

    private boolean extractNoise = false;

    /**
     * number of wavelet scales.  a lower value reveals more fine
     * scale features.  The default is 5.
     */
    private int nScale = 5;
    /**
     * wavelength of smallest scale filter.  The default is 3.
     */
    private int minWavelength = 3;
    /**
     * scaling factor between successive filters.  The default is 2.1.
     */
    private float mult = 2.1f;
    /**
     * ratio of standard deviation of Gaussian describing the
     * log Gabor's filter's transfer function in the frequency domain to the
     * filter center frequency.  The default is 0.55f.
     */
    private float sigmaOnf = 0.55f;
    /**
     * number of standard deviations of the noise energy beyond the
     * mean at which we set the noise threshold point.  You may want to vary this
       up to a value of 10 or 20 for noisy images.
       The default is 5.
     */
    private int k = 5;//2;
    /**
     * The fractional measure of frequency spread below which phase
     * congruency values get penalized.  The default is 0.5f.
     */
    private float cutOff = 0.5f;
    /**
     * Controls the sharpness of the transition in the sigmoid function
     * used to weight phase congruency for frequency spread.  The default is 10.
     */
    private float g = 10;
    /**
     * factor to apply to the calculated phase
     * deviation result.  Increasing this sharpens the edge responses, but can
     * also attenuate their magnitude if the gain is too large.  Sensible values
     * to use lie in the range 1-2.  The default is 1.5f.
     */
    private float deviationGain = 1.5f;
    /**
     * Parameter specifies method used to determine noise
     * statistics: -1 use median of smallest scale filter responses;
     * -2 use mode of smallest scale filter responses;
     * 0 turns off all noise compensation; and
     * > 0 use noiseMethod value as the fixed noise threshold.
     * The default is -1.
     */
    private int noiseMethod = -1;
    /**
     * the low threshold fraction of 1
     */
    private double tLow = 0.1;
    /**
     * the high threshold fraction of 1.  The default is 0.3.
     * Note that this is ignored if useAdaptiveThreshold is true, which it is
     * by default.
     */
    private double tHigh = 0.3;
    /**
     * a flag indicating use of adaptive thresholds in the 2 layer filter
     */
    private boolean useAdaptiveThreshold = true;

    public PhaseCongruencyDetector() {

    }

    public void setToDebug() {
        doPlot = true;
    }

    /**
     * noise is extracted and textures are examined.  this
     * feature name may change depending upon results in progress.
     * Note that k must be the default of 5 or less for this feature.
     */
    public void setToExtractNoise() {
        if (k > 5 || noiseMethod != -1) {
            throw new IllegalStateException("currently, the algorihm needs"
                + " k <= 5 and noiseMethod == -1 for the differencing of pixels.  "
                + " k=" + k + " and noiseMethod=" + noiseMethod);
        }
        this.extractNoise = true;
    }

    /**
     * set the number of wavelet scales.
     * @param n number of wavelet scales.  a lower value reveals more fine
     * scale features.  The default is 5.
     */
    public void setNScales(int n) {
        nScale = n;
    }

    /**
     * set the scale number of the smallest scale filter.
     * @param m wavelength of smallest scale filter.  The default is 3.
     */
    public void setMinWavelength(int m) {
        minWavelength = m;
    }

    /**
     *
     * @param f scaling factor between successive filters.  The default is 2.1.
     */
    public void setMult(float f) {
        mult = f;
    }

    /**
     *
     * @param s ratio of standard deviation of Gaussian describing the
     * log Gabor's filter's transfer function in the frequency domain to the
     * filter center frequency.  The default is 0.55f.
     */
    public void setSigmaOnf(float s) {
        sigmaOnf = s;
    }

    /**
     *
     * @param k number of standard deviations of the noise energy beyond the
     * mean at which we set the noise threshold point.  You may want to vary
     * this up to a value of 10 or 20 for noisy images.
       The default is 5.
     */
    public void setK(int k) {
        if (extractNoise && k > 5) {
            throw new IllegalStateException("currently, the algorihm needs"
                + " k <= 5 for the differencing of pixels requested by"
                + " 'extractNoise'.");
        }
        this.k = k;
    }

    /**
     *
     * @param c The fractional measure of frequency spread below which phase
     * congruency values get penalized.  The default is 0.5f.
     */
    public void setCutOff(float c) {
        cutOff = c;
    }

    /**
     *
     * @param g Controls the sharpness of the transition in the sigmoid function
     * used to weight phase congruency for frequency spread.  The default is 10.
     */
    public void setG(int g) {
        this.g = g;
    }

    /**
     *
     * @param d Amplification to apply to the calculated phase
     * deviation result.  Increasing this sharpens the edge responses, but can
     * also attenuate their magnitude if the gain is too large.  Sensible values
     * to use lie in the range 1-2.  The default is 1.5f.
     */
    public void setDevitionGain(float d) {
        deviationGain = d;
    }

    /**
     *
     * @param m Parameter specifies method used to determine noise
     * statistics: -1 use median of smallest scale filter responses;
     * -2 use mode of smallest scale filter responses;
     * 0 turns off all noise compensation; and
     * > 0 use noiseMethod value as the fixed noise threshold.
     * The default is -1.
     */
    public void setNoiseMethod(int m) {
        if (extractNoise && (m != -1)) {
            throw new IllegalStateException("the 'extract noise' feature is set"
                + " and this requires noiseMethod==-1.");
        }
        this.noiseMethod = m;
    }

    /**
     *
     * @param tLow the low threshold fraction of 1.
     * The default is 0.1.
     */
    public void setLowThreshold(double tLow) {
        this.tLow = tLow;
    }

    /**
     *
     * @param t the high threshold fraction of 1.  The default is 0.3.
     * Note that this is ignored if useAdaptiveThreshold is true, which it is
     * by default.
     */
    public void setHighThreshold(double t) {
        this.tHigh = t;
    }

    /**
     * turn off the use of adaptive thresholds in the 2 layer filter, and
     * instead use tLow and tHigh only.
     */
    public void setToNotUseAdaptiveThreshold() {
        this.useAdaptiveThreshold = false;
    }

    /**
     * use the phase congruency method of transformations to create
     * edge maps, orientation, phase angle and a suggested threshold.
     * @param img
     * @return
       <pre>
       NOTE: the return products use row major and column major notation, so
       read the method products getters carefully.

       Returned values:
         phaseCongruency  - Phase congruency indicating edge significance
                            (values are in range 0 to 1.)
         orientation      - Orientation image in integer degrees 0-180,
                            positive anticlockwise.
                            0 corresponds to a vertical edge, 90 is horizontal.
         phaseAngle       - Local weighted mean phase angle at every point in
                            the image. A value of
                                pi/2 corresponds to a bright line,
                                0 corresponds to a step and
                                -pi/2 is a dark line.
         threshold        - Calculated noise threshold (can be useful for
                            diagnosing noise characteristics of images).
                            Once you know this you can then specify fixed
                            thresholds and save some computation time.
      </pre>
     */
    public PhaseCongruencyProducts phaseCongMono(GreyscaleImage img) {

        long t0 = System.currentTimeMillis();

        int nCols = img.getWidth();
        int nRows = img.getHeight();

        //Periodic Fourier transform of image, using default normalization
        // perfft2 results use notation a[row][col]
        PeriodicFFT perfft2 = new PeriodicFFT();
        //IM = perfft2(im);                   %
        //S, P, s, p where S = FFT of smooth, P = FFT of periodic, s=spatial smooth, p = spatial p
        Complex[][][] perfResults = perfft2.perfft2(img, false);
        Complex[][] capIm = perfResults[1];

        /*
        sumAn  = zeros(rows,cols);          % Matrix for accumulating filter response
                                            % amplitude values.
        sumf   = zeros(rows,cols);          % ft is phase angle
        sumh1  = zeros(rows,cols);
        sumh2  = zeros(rows,cols);
        */
        double[][] sumAn = new double[nRows][];
        double[][] sumF = new double[sumAn.length][];
        double[][] sumH1 = new double[sumAn.length][];
        double[][] sumH2 = new double[sumAn.length][];
        for (int row = 0; row < nRows; ++row) {
            sumAn[row] = new double[nCols];
            sumF[row] = new double[nCols];
            sumH1[row] = new double[nCols];
            sumH2[row] = new double[nCols];
        }

        /*
        Generate grid data for constructing filters in the frequency domain
        [radius, u1, u2] = filtergrid(rows, cols);
        */
        // results use notation a[row][col]
        FilterGrid fg = new FilterGrid();
        FilterGridProducts fgProducts = fg.filtergrid(nRows, nCols);

        /*
        Get rid of the 0 radius value in the middle (at top left corner after
        fftshifting) so that taking the log of the radius, or dividing by the
        radius, will not cause trouble.
         radius(1,1) = 1;
        */
        fgProducts.getRadius()[0][0] = 1;

        /*
         % Construct the monogenic filters in the frequency domain.  The two
         % filters would normally be constructed as follows
         %    H1 = i*u1./radius;
         %    H2 = i*u2./radius;
         % However the two filters can be packed together as a complex valued
         % matrix, one in the real part and one in the imaginary part.  Do this by
         % multiplying H2 by i and then adding it to H1 (note the subtraction
         % because i*i = -1).  When the convolution is performed via the fft the
         % real part of the result will correspond to the convolution with H1 and
         % the imaginary part with H2.  This allows the two convolutions to be
         % done as one in the frequency domain, saving time and memory.
         H = (1i*u1 - u2)./radius;
        */
        // results use notation a[row][col]
        double[][] u1 = fgProducts.getU1();
        double[][] u2 = fgProducts.getU2();
        double[][] radius = fgProducts.getRadius();
        Complex[][] capH = new Complex[nRows][];
        for (int row = 0; row < nRows; ++row) {
            capH[row] = new Complex[nCols];
            for (int col = 0; col < nCols; ++col) {
                double re = -u2[row][col]/radius[row][col];
                double im = u1[row][col]/radius[row][col];
                capH[row][col] = new Complex(re, im);
            }
        }

        /*
         % First construct a low-pass filter that is as large as possible, yet falls
         % away to zero at the boundaries.  All filters are multiplied by
         % this to ensure no extra frequencies at the 'corners' of the FFT are
         % incorporated as this can upset the normalisation process when
         % calculating phase congruency
         lp = lowpassfilter([rows,cols],.45,15);    % Radius .4, 'sharpness' 15
        */
        // results use notation a[row][col]
        LowPassFilter lpFilter = new LowPassFilter();
        double[][] lp = lpFilter.lowpassfilter(nRows, nCols, 0.45f, 15);

        ImageProcessor imageProcessor = new ImageProcessor();

        double[][] maxAN = null;

        double tau = noiseMethod;
        // keeping taus in case need to increase noise estimate
        double sqml4 = Math.sqrt(Math.log(4));
        double logGaborDenom = 2. * Math.pow(Math.log(sigmaOnf), 2);

        double[][] width = new double[nRows][];
        double[][] weight = new double[nRows][];
        for (int row = 0; row < nRows; ++row) {
            width[row] = new double[nCols];
            weight[row] = new double[nCols];
        }

        for (int s = 0; s < nScale; ++s) {

            // Centre frequency of filter.
            double wavelength = minWavelength * Math.pow(mult, s);

            double fo = 1.0/wavelength;

            // use notation a[row][col]
            double[][] logGabor = new double[nRows][];
            for (int row = 0; row < nRows; ++row) {
                logGabor[row] = new double[nCols];
                for (int col = 0; col < nCols; ++col) {
                    double v = Math.log(radius[row][col]/fo);
                    v *= v;
                    v = Math.exp(-v/logGaborDenom);
                    //logGabor = logGabor.*lp;
                    logGabor[row][col] = lp[row][col] * v;
                }
            }
            logGabor[0][0] = 0;

            //DEBUG(logGabor, "s=" + s + " logGabor*low pass filter", s);

            // uses notation a[row][col]
            //Bandpassed image in the frequency domain
            Complex[][] capIMF = new Complex[nRows][];
            for (int row = 0; row < nRows; ++row) {
                capIMF[row] = new Complex[nCols];
                for (int col = 0; col < nCols; ++col) {
                   capIMF[row][col] = capIm[row][col].times(logGabor[row][col]);
                }
            }

            // uses notation a[row][col]
            //  Bandpassed image in spatial domain.
            //  f = real(ifft2(IMF));
            // the functions used in other code are not normalized on fft,
            // but are by inverse fft so need a combined division here by nomr=nRows*nCols
            Complex[][] fComplex = imageProcessor.create2DFFT(capIMF, false, false);

            double norm = nRows * nCols;

            //h = ifft2(IMF.*H);
            Complex[][] capIMFH = new Complex[nRows][];
            for (int row = 0; row < nRows; ++row) {
                capIMFH[row] = new Complex[nCols];
                for (int col = 0; col < nCols; ++col) {
                    capIMFH[row][col] = capIMF[row][col].times(capH[row][col]);
                }
            }
            // result needs to be divided by norm=nRows*nCols
            Complex[][] h = imageProcessor.create2DFFT(capIMFH, false, false);

            /*
            h1 = real(h);
            h2 = imag(h);
            An = sqrt(f.^2 + h1.^2 + h2.^2); % Amplitude of this scale component.
            sumAn = sumAn + An;              % Sum of component amplitudes over scale.
            sumf  = sumf  + f;
            sumh1 = sumh1 + h1;
            sumh2 = sumh2 + h2;
            */
            // uses notation a[row][col]
            double[][] aN = new double[nRows][];
            for (int row = 0; row < nRows; ++row) {
                aN[row] = new double[nCols];
                for (int col = 0; col < nCols; ++col) {
                    // results of inverse transforms need normalization
                    double f0 = fComplex[row][col].re()/norm;
                    double h1 = h[row][col].re()/norm;
                    double h2 = h[row][col].im()/norm;
                    aN[row][col] = Math.sqrt(f0*f0 + h1*h1 + h2*h2);
                    sumAn[row][col] += aN[row][col];
                    sumF[row][col] += f0;
                    sumH1[row][col] += h1;
                    sumH2[row][col] += h2;
                }
            }

            /*
            At the smallest scale estimate noise characteristics from the
            distribution of the filter amplitude responses stored in sumAn.
            tau is the Rayleigh parameter that is used to describe the
            distribution.
            */
            if (s == 0) {
                if (noiseMethod == -1) {
                    //Use median to estimate noise statistics
                    //tau = median(sumAn(:))/sqrt(log(4));
                    //TODO: could improve this below O(N*lg2(N)) with histograms
                    //      and assumptions of bin sizes...when have an
                    //      implementation of Multi-Level-Buckets, revisit this
                    double median = MiscMath.findMedian(sumAn);
                    tau = median/sqml4;
                } else if (noiseMethod == -2) {
                    //Use mode to estimate noise statistics
                    //tau = rayleighmode(sumAn(:));
                    tau = rayleighMode(sumAn);
                }
                maxAN = aN;
            } else {
                // Record maximum amplitude of components across scales.  This is needed
                // to determine the frequency spread weighting.
                //maxAN = max(maxAN, An);
                // uses notation a[row][col]
                for (int row = 0; row < nRows; ++row) {
                    for (int col = 0; col < nCols; ++col) {
                        maxAN[row][col] = Math.max(maxAN[row][col], aN[row][col]);
                    }
                }
            }

            /*
            Form weighting that penalizes frequency distributions that are
            particularly narrow.  Calculate fractional 'width' of the frequencies
            present by taking the sum of the filter response amplitudes and dividing
            by the maximum component amplitude at each point on the image.  If
            there is only one non-zero component width takes on a value of 0, if
            all components are equal width is 1.
            width = (sumAn./(maxAn + epsilon) - 1) / (nscale-1);

            Now calculate the sigmoidal weighting function.
            weight = 1.0 ./ (1 + exp( (cutOff - width)*g));
            */
            // uses notation a[row][col]

            double dn = (double)nScale - 1.;

            if (dn > 0) {
                for (int row = 0; row < nRows; ++row) {
                    for (int col = 0; col < nCols; ++col) {
                        double a = sumAn[row][col]/(maxAN[row][col] + epsilon);
                        width[row][col] = (a - 1.)/dn;
                        double v = Math.exp(g*(cutOff - width[row][col]));
                        weight[row][col] = 1./(1. + v);
                    }
                }
            }
        } // end for each scale

        /*
        Automatically determine noise threshold

        Assuming the noise is Gaussian the response of the filters to noise will
        form Rayleigh distribution.  We use the filter responses at the smallest
        scale as a guide to the underlying noise level because the smallest scale
        filters spend most of their time responding to noise, and only
        occasionally responding to features. Either the median, or the mode, of
        the distribution of filter responses can be used as a robust statistic to
        estimate the distribution mean and standard deviation as these are related
        to the median or mode by fixed constants.  The response of the larger
        scale filters to noise can then be estimated from the smallest scale
        filter response according to their relative bandwidths.

        This code assumes that the expected reponse to noise on the phase
        congruency calculation is simply the sum of the expected noise responses
        of each of the filters.  This is a simplistic overestimate, however these
        two quantities should be related by some constant that will depend on the
        filter bank being used.  Appropriate tuning of the parameter 'k' will
        allow you to produce the desired output. (though the value of k seems to
        be not at all critical)
        */

        // uses notation a[row][col]
        // ft is phase angle
        double[][] orientation = new double[nRows][];
        double[][] ft = new double[nRows][];
        double[][] energy = new double[nRows][];
        double[][] pc = new double[nRows][];
        for (int row = 0; row < nRows; ++row) {
            orientation[row] = new double[nCols];
            ft[row] = new double[nCols];
            energy[row] = new double[nCols];
            pc[row] = new double[nCols];
        }

        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {
                double gY = -sumH2[row][col];
                double gX = sumH1[row][col];
                if (gY < 0) {
                    gX *= -1;
                    gY *= -1;
                }
                orientation[row][col] = Math.atan2(gY, gX);
                if (orientation[row][col] == Math.PI) {
                    orientation[row][col] = 0;
                }
                // orientation values now range 0 - pi
                // Quantize to 0 - 180 degrees (for NONMAXSUP)
                orientation[row][col] = (int)(orientation[row][col]*180./Math.PI);

                double h1Sq = sumH1[row][col];
                h1Sq *= h1Sq;
                double h2Sq = sumH2[row][col];
                h2Sq *= h2Sq;

                //Feature type - a phase angle -pi/2 to pi/2.
                ft[row][col] = Math.atan2(sumF[row][col], Math.sqrt(h1Sq + h2Sq));

                //overall energy
                double v0 = sumF[row][col];
                v0 *= v0;
                energy[row][col] = Math.sqrt(v0 + h1Sq + h2Sq);
            }
        }

        /*
        % Compute phase congruency.  The original measure,
        % PC = energy/sumAn
        % is proportional to the weighted cos(phasedeviation).  This is not very
        % localised so this was modified to
        % PC = cos(phasedeviation) - |sin(phasedeviation)|
        % (Note this was actually calculated via dot and cross products.)  This measure
        % approximates
        % PC = 1 - phasedeviation.

        % However, rather than use dot and cross products it is simpler and more
        % efficient to simply use acos(energy/sumAn) to obtain the weighted phase
        % deviation directly.  Note, in the expression below the noise threshold is
        % not subtracted from energy immediately as this would interfere with the
        % phase deviation computation.  Instead it is applied as a weighting as a
        % fraction by which energy exceeds the noise threshold.  This weighting is
        % applied in addition to the weighting for frequency spread.  Note also the
        % phase deviation gain factor which acts to sharpen up the edge response. A
        % value of 1.5 seems to work well.  Sensible values are from 1 to about 2.

        PC = weight.*max(1 - deviationGain*acos(energy./(sumAn + epsilon)),0) ...
              .* max(energy-T,0)./(energy+epsilon);
        */

        double threshold;
        double thresholdLowNz = -1;
        if (noiseMethod >= 0) {
            //fixed noise threshold
            threshold = noiseMethod;
        } else {
            //Estimate the effect of noise on the sum of the filter responses as
            //the sum of estimated individual responses (this is a simplistic
            //overestimate). As the estimated noise response at succesive scales
            //is scaled inversely proportional to bandwidth we have a simple
            //geometric sum.
            //totalTau = tau * (1 - (1/mult)^nscale)/(1-(1/mult));
            double totalTau = tau * (1. - Math.pow((1./mult), nScale))/(1. - (1./mult));

            // Calculate mean and std dev from tau using fixed relationship
            // between these parameters and tau. See
            // http://mathworld.wolfram.com/RayleighDistribution.html
            double EstNoiseEnergyMean = totalTau * Math.sqrt(Math.PI/2.);
            double EstNoiseEnergySigma = totalTau * Math.sqrt((4. - Math.PI)/2.);

            threshold = Math.max(EstNoiseEnergyMean
                + ((float)k) * EstNoiseEnergySigma, epsilon);

            if (extractNoise) {
                // empirically, k=10 is chosen
                thresholdLowNz = Math.max(EstNoiseEnergyMean
                    + 7.f * EstNoiseEnergySigma, epsilon);
            }
        }

        double[][] lowNzFactors = null;
        if (extractNoise) {
            lowNzFactors = new double[nRows][];
            for (int row = 0; row < nRows; ++row) {
                lowNzFactors[row] = new double[nCols];
            }
        }

        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {

                double eDiv = Math.acos(energy[row][col]/(sumAn[row][col] + epsilon));

                double a = weight[row][col]
                    * Math.max(1. - deviationGain * eDiv, 0)
                    / (energy[row][col] + epsilon);

                double eMax = Math.max(energy[row][col] - threshold, 0);

                pc[row][col] = a * eMax;

                if (extractNoise) {
                    lowNzFactors[row][col] = Math.max(energy[row][col] -
                        thresholdLowNz, 0);
                    lowNzFactors[row][col] /= eMax;
                }
            }
        }

        PhaseCongruencyProducts products = new PhaseCongruencyProducts(pc,
            orientation, ft, threshold);

        NonMaximumSuppression ns = new NonMaximumSuppression();

        double[][] thinnedPC = ns.nonmaxsup(products.getPhaseCongruency(),
            products.getOrientation(), 1.2, new HashSet<PairInt>());

        products.setParameters(nScale, minWavelength, mult, sigmaOnf, k, cutOff,
            g, deviationGain, noiseMethod, tLow, tHigh);

        int[][] thinned = createEdges(products.getPhaseCongruency(), thinnedPC,
            products.getPhaseAngle(), tLow, tHigh);

        products.setThinnedImage(thinned);

        if (extractNoise) {

            int distMinimum = 70;
            
            /*
            using a higher k=7, subtracting the thinned results from the default
            thinned image.
            */

            double[][] pcLowNz = copy(products.getPhaseCongruency());
            for (int row = 0; row < nRows; ++row) {
                for (int col = 0; col < nCols; ++col) {
                    if (thinned[row][col] == 0) {
                        pcLowNz[row][col] = 0;
                    } else {
                        pcLowNz[row][col] *= lowNzFactors[row][col];
                    }
                }
            }

            int[][] thinnedLowNz = applyHysThresh(pcLowNz,
                products.getPhaseAngle(), tLow, tHigh, true);

            Set<PairInt> noisePoints = new HashSet<PairInt>();

            products.setNoiseyPixels(new HashSet<PairInt>());
            
            int dtN = 0;
            int[][] dt = new int[nRows][];
            for (int row = 0; row < nRows; ++row) {
                dt[row] = new int[nCols];
                for (int col = 0; col < nCols; ++col) {
                    // where thinned is > 0 and thinnedLowNz == 0 is a noise pixel
                    int v = thinned[row][col] - thinnedLowNz[row][col];
                    if (v > 0) {
                        // noise = points in thinned, but not in low noise thinned
                        noisePoints.add(new PairInt(row, col));
                        products.getNoiseyPixels().add(new PairInt(col, row));
                    } else {
                        // points that are signal in both
                        if (thinnedLowNz[row][col] > 0) {
                            dt[row][col] = 1;
                            dtN++;
                        }
                    }
                }
            }

            products.setThinnedLowNoise(thinnedLowNz);            

            if (doPlot) {
                int count = 0;
                Image dbg0 = new Image(nCols, nRows);
                for (int row = 0; row < nRows; ++row) {
                    for (int col = 0; col < nCols; ++col) {
                        int v = thinnedLowNz[row][col];
                        if (v > 0) {
                            dbg0.setRGB(col, row, v, v, v);
                        }
                    }
                }
                for (PairInt p : noisePoints) {
                    int[] clr = ImageIOHelper.getNextRGB(count);
                    ImageIOHelper.addPointToImage(p.getY(), p.getX(), dbg0,
                        1, clr[0], clr[1], clr[2]);
                    count++;
                }
                MiscDebug.writeImage(dbg0, "_a_noise_");
            }

            /*
            -- finding distance to non-noise points of the thinned default image
               from the noise points using a distance transform
            -- sorting those points by decreasing distance.
            */

            DistanceTransform dTrans = new DistanceTransform();
            dt = dTrans.applyMeijsterEtAl(dt);

            Set<PairIntWithIndex> points2
                = new HashSet<PairIntWithIndex>();
            for (PairInt p : noisePoints) {
                points2.add(new PairIntWithIndex(p.getX(), p.getY(),
                    points2.size()));
            }
            DTClusterFinder<PairIntWithIndex> cFinder
                = new DTClusterFinder<PairIntWithIndex>(points2,
                nRows + 1, nCols + 1);
            //cFinder.setToDebug();
            cFinder.setMinimumNumberInCluster(1);
            cFinder.calculateCriticalDensity();
            cFinder.findClusters();
            final int n = cFinder.getNumberOfClusters();

            float[] clusterSizes = new float[n];
            Image dbg = new Image(nCols, nRows);
            for (int i = 0; i < n; ++i) {
                Set<PairIntWithIndex> set = cFinder.getCluster(i);
                clusterSizes[i] = set.size();
                int[] clr = ImageIOHelper.getNextRGB(i);
                for (PairIntWithIndex p : set) {
                    ImageIOHelper.addPointToImage(p.getY(), p.getX(), dbg,
                        1, clr[0], clr[1], clr[2]);
                }
            }
            MiscDebug.writeImage(dbg, "_a_noise_clusters_");
            
            HistogramHolder hist = Histogram.createSimpleHistogram(
                0, 40, 20, clusterSizes, 
                Errors.populateYErrorsBySqrt(clusterSizes));
            int peakIdx = MiscMath.findYMaxIndex(hist.getYHist());
            assert(peakIdx != -1); 
            float sizeLimit = hist.getXHist()[peakIdx];
            if (peakIdx < (hist.getXHist().length - 2)) {
                sizeLimit = hist.getXHist()[peakIdx + 2];
            }
            
            System.out.println("nClusters=" + n 
                + " peakHist y=" + hist.getYHist()[peakIdx] 
                + " peakHist x=" + hist.getXHist()[peakIdx] 
                + " nPix=" + img.getNPixels() + 
                " peakX/nPix=" + 
                hist.getXHist()[peakIdx]/(float)img.getNPixels() +
                " peakY/nPix=" + 
                (float)hist.getYHist()[peakIdx]/(float)img.getNPixels()
                + " peakY/peakX=" +
                (float)hist.getYHist()[peakIdx]/hist.getXHist()[peakIdx]
            );
            
            try {
                hist.plotHistogram("cluster sizes", "_cluster_sizes_");
            } catch (IOException ex) {
                Logger.getLogger(PhaseCongruencyDetector.class.getName()).log(Level.SEVERE, null, ex);
            }

            int[] dist = new int[n];           
            int[] indexes = new int[n];
            int count = 0;
            for (int i = 0; i < n; ++i) {
                Set<PairIntWithIndex> set = cFinder.getCluster(i);
                if (set.size() > sizeLimit) {
                    continue;
                }
                int sumD = 0;
                int countD = 0;
                for (PairIntWithIndex p : set) {
                    if (noisePoints.contains(new PairInt(p.getX(), p.getY())) && dt[p.getX()][p.getY()] > 0) {
                        sumD += dt[p.getX()][p.getY()];
                        countD++;
                    }
                }
                if (countD > 0) {
                    indexes[count] = i;
                    sumD /= countD;
                    dist[count] = sumD;
                    count++;
                }
            }
            indexes = Arrays.copyOf(indexes, count);
            dist = Arrays.copyOf(dist, count);
            
            // sort by decr dist
            QuickSort.descendingSort(dist, indexes);
            
            if (doPlot) {
                //histogram of point distances
                float[] values = new float[n];
                for (int i = 0; i < dist.length; ++i) {
                    values[i] = dist[i];
                }
                //TODO: consider xmax based on image size
                hist = Histogram.createSimpleHistogram(
                    0, 400, 20, values,
                    Errors.populateYErrorsBySqrt(values));

                try {
                    hist.plotHistogram("dist from edges", "_noise_distances_");
                } catch (IOException ex) {
                    Logger.getLogger(PhaseCongruencyDetector.class.getName()).log(Level.SEVERE, null, ex);
                }
            }

            int end = 100;
            if (dist.length < end) {
                end = dist.length;
            }
            System.out.println("dist[0]=" + dist[0] + " dist[end-1]=" + dist[end-1]);

            int np = 0;

            List<Set<PairInt>> subsetNoise = new ArrayList<Set<PairInt>>();
            Image dbg0 = new Image(nCols, nRows);
            for (int i = 0; i < end; ++i) {
                if (dist[i] < distMinimum) {
                    break;
                }
                int idx = indexes[i];
                int[] clr = ImageIOHelper.getNextRGB(i);
                Set<PairIntWithIndex> set = cFinder.getCluster(idx);
                Set<PairInt> set2 = new HashSet<PairInt>();
                for (PairIntWithIndex p : set) {
                    ImageIOHelper.addPointToImage(p.getY(), p.getX(), dbg0,
                        1, clr[0], clr[1], clr[2]);
                    set2.add(new PairInt(p.getY(), p.getX()));
                }
                subsetNoise.add(set2);
                np += set2.size();
            }
            for (int row = 0; row < nRows; ++row) {
                for (int col = 0; col < nCols; ++col) {
                    int v = thinnedLowNz[row][col];
                    if (v > 0) {
                        dbg0.setRGB(col, row, 255, 255, 255);
                    }
                }
            }
            MiscDebug.writeImage(dbg0, "_a_texture_candidates_");
         
            System.out.println("np=" + np);
            
            products.setSubsetNoise(subsetNoise);
        }

        long t1 = System.currentTimeMillis();

        System.out.println(((t1 - t0)*1E-3) + " seconds for phasecongmono");

        return products;
    }

    private double[][] copy(double[][] a) {
        double[][] cp = new double[a.length][];
        for (int i = 0; i < a.length; ++i) {
            cp[i] = Arrays.copyOf(a[i], a[i].length);
        }
        return cp;
    }

    /**
     * creating edges using thinned phase angle image.
     * @param pc
     * @param thinnedPC
     * @param phaseAngle
     * @param tLow
     * @param tHigh
     * @return
     */
    public int[][] createEdges(double[][] pc, double[][] thinnedPC,
        double[][] phaseAngle, double tLow, double tHigh) {

        int nRows = phaseAngle.length;
        int nCols = phaseAngle[0].length;

        Set<PairInt> brightLinePoints = new HashSet<PairInt>();
        Set<PairInt> darkLinePoints = new HashSet<PairInt>();
        Set<PairInt> stepPoints = new HashSet<PairInt>();

        double piDiv4 = Math.PI/4.;
        double tolerance = piDiv4/2.;

        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {

                if (thinnedPC[row][col] < tLow) {
                    continue;
                }
                // placing values closer to -pi/2 than 0 are darklines
                //         values closer to +p1/2 than 0 are brightlines
                //         else they are closer to 0 are steps
                double v = phaseAngle[row][col];
                if (Math.abs(-piDiv4 - v) < tolerance) {
                    darkLinePoints.add(new PairInt(row, col));
                } else if (Math.abs(piDiv4 - v) < tolerance) {
                    brightLinePoints.add(new PairInt(row, col));
                } else {
                    stepPoints.add(new PairInt(row, col));
                }
            }
        }

        //DEBUG
        if (doPlot) {
            Image paImage = new Image(nCols, nRows);
            int nExtraForDot = 0;
            for (PairInt p : brightLinePoints) {
                int x = p.getY();
                int y = p.getX();
                ImageIOHelper.addPointToImage(x, y, paImage, nExtraForDot,
                    0, 255, 0);
            }
            for (PairInt p : stepPoints) {
                int x = p.getY();
                int y = p.getX();
                ImageIOHelper.addPointToImage(x, y, paImage, nExtraForDot,
                    255, 0, 0);
            }
            for (PairInt p : darkLinePoints) {
                int x = p.getY();
                int y = p.getX();
                ImageIOHelper.addPointToImage(x, y, paImage, nExtraForDot,
                    127, 0, 255);
            }
            MiscDebug.writeImage(paImage, "_a_0_pa_components_" +
                MiscDebug.getCurrentTimeFormatted());
        }

        stepPoints.addAll(brightLinePoints);
        stepPoints.addAll(darkLinePoints);

        double[][] thinned1 = new double[nRows][nCols];
        for (int i = 0; i < nRows; ++i) {
            thinned1[i] = new double[nCols];
        }
        for (PairInt p : stepPoints) {
            thinned1[p.getX()][p.getY()] = pc[p.getX()][p.getY()];
        }
        stepPoints = null;

        int[][] thinned = applyHysThresh(thinned1, phaseAngle, tLow, tHigh, true);

        if (doPlot) {
            Image pImage = new Image(nCols, nRows);
            Image tImage = new Image(nCols, nRows);
            int nExtraForDot = 0;
            for (int row = 0; row < nRows; ++row) {
                for (int col = 0; col < nCols; ++col) {
                    int v = (int)Math.round(255 * pc[row][col]);
                    if (thinned[row][col] > 0) {
                        ImageIOHelper.addPointToImage(col, row, pImage,
                            nExtraForDot, 255, 255, 255);
                        ImageIOHelper.addPointToImage(col, row, tImage,
                            nExtraForDot, v, v, v);
                    } else {
                        ImageIOHelper.addPointToImage(col, row, tImage,
                            nExtraForDot, v, 0, 0);
                    }
                }
            }
            MiscDebug.writeImage(pImage, "_a_1_thinned_1_" +
                MiscDebug.getCurrentTimeFormatted());
            MiscDebug.writeImage(tImage, "_a_2_pc_1_" +
                MiscDebug.getCurrentTimeFormatted());
        }

        int[][] thinnedbw = new int[nRows][nCols];
        for (int i = 0; i < nRows; ++i) {
            thinnedbw[i] = new int[nCols];
            for (int j = 0; j < nCols; ++j) {
                if (thinned[i][j] > 0) {
                    thinnedbw[i][j] = 1;
                }
            }
        }

        MorphologicalFilter mFilter = new MorphologicalFilter();
        int[][] skel = mFilter.bwMorphThin(thinnedbw, Integer.MAX_VALUE);

        Set<PairInt> points = new HashSet<PairInt>();
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                int m = skel[i][j];
                thinned[i][j] *= m;
                if (thinned[i][j] > 0) {
                    points.add(new PairInt(i, j));
                }
            }
        }

        //DEBUG
        if (doPlot) {
            Image paImage = new Image(nCols, nRows);
            int nExtraForDot = 0;
            for (PairInt p : points) {
                int x = p.getY();
                int y = p.getX();
                ImageIOHelper.addPointToImage(x, y, paImage, nExtraForDot,
                    255, 255, 255);
            }
            MiscDebug.writeImage(paImage, "_a_3_edge_from_phase_angle_thin_0_" +
                MiscDebug.getCurrentTimeFormatted());
        }

        //ImageProcessor imageProcessor = new ImageProcessor();
        //imageProcessor.applyThinning(points,nRows, nCols);

        //pltc.correctForLineHatHoriz(points, nRows, nCols);
        //pltc.correctForLineHatVert(points, nRows, nCols);

        // there are a very small number of clumps thicker than 1 pixel.
        //ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
        //lt.applyLineThinner(points, 0, nRows - 1, 0, nCols - 1);

        /*
        PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
        pltc.correctForIsolatedPixels(stepPoints);

        for (int i = 0; i < nRows; ++i) {
            Arrays.fill(thinned[i], 0);
        }
        for (PairInt p : points) {
            thinned[p.getX()][p.getY()] = 255;
        }

        //DEBUG
        if (doPlot) {
            Image paImage = new Image(nCols, nRows);
            int nExtraForDot = 0;
            for (PairInt p : points) {
                int x = p.getY();
                int y = p.getX();
                ImageIOHelper.addPointToImage(x, y, paImage, nExtraForDot,
                    255, 0, 0);
            }
            MiscDebug.writeImage(paImage, "_a_3_edge_from_phase_angle_thin_1_" +
                MiscDebug.getCurrentTimeFormatted());
        }
        */

        return thinned;
    }

    /**
     * creating edges using thinned phase angle image.
     * @param products
     * @param tLow
     * @param tHigh
     */
    private int countEdgePoints(PhaseCongruencyProducts products,
        double[][] thinnedPC, double tLow, double tHigh) {

        double[][] phaseAngle = products.getPhaseAngle();

        int nRows = phaseAngle.length;
        int nCols = phaseAngle[0].length;

        int count = 0;

        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {
                if (thinnedPC[row][col] < tLow) {
                    continue;
                }
                count++;
            }
        }

        return count;
    }

    public class PhaseCongruencyProducts {

        /**
         * indicates edge significance
         */
        private final double[][] phaseCongruency;

        /**
         * Orientation image in integer degrees 0-180 with positive anticlockwise.
         */
        private final double[][] orientation;

        /**
         * Local weighted mean phase angle at every point in the image.
         * A value of pi/2 corresponds to a bright line, 0 corresponds to a
         * step and -pi/2 is a dark line.
         */
        private final double[][] phaseAngle;

        /**
         * Calculated noise threshold (can be useful for diagnosing noise
         * characteristics of images).  Once you know this you can then specify
         * fixed thresholds and save some computation time.
         */
        private final double threshold;

        private int[][] thinned = null;
        
        private int[][] thinnedLowNoise = null;

        /**
         * noisey pixels are found during the 2 level threshold hysteresis
         * phase and they sometimes contain textures such as vegetation or
         * bricks and sometimes contain edge points that should be restored
         * before this stage.
         */
        private Set<PairInt> noiseyPixels = null;

        /**
         * these are a subset of noiseyPixels that are the smallest
           clusters of points and furthest away from edges
           in the lower noise image
         */
        private List<Set<PairInt>> subsetNoise = null;
        
        private PhaseCongruencyParameters parameters = null;

        public PhaseCongruencyProducts(double[][] pc, double[][] or,
            double[][] ft, double thr) {
            this.phaseCongruency = copy(pc);
            this.orientation = copy(or);
            this.phaseAngle = copy(ft);
            this.threshold = thr;
        }

        /**
         * set the thinned phase congruence image, a.k.a. the edge image.
         * the image should be in row-major notation.
         * @param thImg
         */
        public void setThinnedImage(int[][] thImg) {
            thinned = copy(thImg);
        }
        
        /**
         * set the thinned low noise image (created if the user specified
         * 'extract noise').
         * the image should be in row-major notation.
         * @param thinnedLowNz 
         */
        private void setThinnedLowNoise(int[][] thinnedLowNz) {
            thinnedLowNoise = thinnedLowNz;
        }

        /**
         * get the thinned phase congruence image, a.k.a. the edge image.
         * Note that the array is accessed as a[row][column], row major notation.
         * @@return edgeImg
         */
        public int[][] getThinned() {
            return thinned;
        }
        
        /**
         * get the thinned phase congruence image that was built if the
         * user specified 'extract noise'.
         * Note that the array is accessed as a[row][column], row major notation.
         * @return 
         */
        public int[][] getThinnedLowNoise() {
            return thinnedLowNoise;
        }

        /**
         * @return noisey pixels are found by subtracting the resulting pc image
         * made with the user given k from a default k which allows
         * more noise into the image.  The resulting "noisey" pixels
         * may be useful for designing texture filters to remove
         * such pixels from better keypoints.
         * Note that the coordinates are in the reference frame
         * of the GreyscaleImage instance, col major notation.
         */
        public Set<PairInt> getNoiseyPixels(){
            return noiseyPixels;
        }
        
        /**
         * @return the subset of noiseyPixels that are the smallest
           clusters of points and furthest away from edges
           in the lower noise image.
           Note that the coordinates are in the reference frame
         * of the GreyscaleImage instance, col major notation.
         */
        public List<Set<PairInt>> getSubsetNoise() {
            return subsetNoise;
        }

        /**
         * set pixels derived from the 'extract noise' feature of the phase
         * congruency detector.
         * Note that the
         * points should be using coordinates that are in the reference frame
         * of the GreyscaleImage instance, col major notation.
         * @param pixels 
         */
        public void setNoiseyPixels(Set<PairInt> pixels){
            noiseyPixels = pixels;
        }

        /**
         * these are a subset of noiseyPixels that are the smallest
           clusters of points and furthest away from edges
           in the lower noise image.  Note that the
         * points should be using coordinates that are in the reference frame
         * of the GreyscaleImage instance, col major notation.
         * @param subsetNoise 
         */
        private void setSubsetNoise(List<Set<PairInt>> subsetNoise) {
            this.subsetNoise = subsetNoise;
        }

        /**
         * return the gradient image produced by phase congruency as a double
         * array of values in range 0 to 1.0.
         * Note that the array is accessed as a[row][column], row major notation.
         * @return the phaseCongruency
         */
        public double[][] getPhaseCongruency() {
            return phaseCongruency;
        }

        /**
         * return the orientation image.
         * Note that the array is accessed as a[row][column], row major notation.
         * @return the orientation
         */
        public double[][] getOrientation() {
            return orientation;
        }

        /**
         * return the phase angle image.
         * <pre>
         * Local weighted mean phase angle at every point in the image.
         * A value of
                pi/2 corresponds to a bright line,
                0 corresponds to a step and
                -pi/2 is a dark line.
         * </pre>
         * Note that the array is accessed as a[row][column], row major notation
         * @return the phaseAngle
         */
        public double[][] getPhaseAngle() {
            return phaseAngle;
        }

        /**
         *
         * @return the threshold
         */
        public double getThreshold() {
            return threshold;
        }

        public PhaseCongruencyParameters getParameters() {
            return parameters;
        }

        private void setParameters(int nScale, int minWavelength, float mult,
            float sigmaOnf, int k, float cutOff, float g, float deviationGain,
            int noiseMethod, double tLow, double tHigh) {

            this.parameters = new PhaseCongruencyParameters();

            parameters.setParameters(nScale, minWavelength, mult, sigmaOnf, k,
                cutOff, g, deviationGain, noiseMethod, tLow, tHigh);
        }

        private void setParameters(PhaseCongruencyParameters params) {
            this.parameters = params;
        }

    }

    private Complex[][] copy(Complex[][] a) {

        Complex[][] b = new Complex[a.length][];
        for (int i = 0; i < b.length; ++i) {
            b[i] = Arrays.copyOf(a[i], a[i].length);
        }

        return b;
    }
    private int[][] copy(int[][] a) {

        int[][] b = new int[a.length][];
        for (int i = 0; i < b.length; ++i) {
            b[i] = Arrays.copyOf(a[i], a[i].length);
        }

        return b;
    }

    /**
    adapted from Kovesis phasecongmono.m as documented in class comments above.

      Mode is computed by forming a histogram of the data over 50 bins and then
      finding the maximum value in the histogram.  Mean and standard deviation
      can then be calculated from the mode as they are related by fixed
      constants.

      mean = mode * sqrt(pi/2)
      std dev = mode * sqrt((4-pi)/2)

      See
      http://mathworld.wolfram.com/RayleighDistribution.html
      http://en.wikipedia.org/wiki/Rayleigh_distribution

     * @param data data assumed to come from a Rayleigh distribution
     * @return
     */
    private double rayleighMode(double[][] data) {

        int nBins = 50;

        float[] values = new float[data.length];
        float max = Float.MIN_VALUE;
        int count = 0;
        for (int j = 0; j < data.length; ++j) {
            for (int i = 0; i < data[j].length; ++i) {
                values[count] = (float)data[j][i];
                if (values[count] > max) {
                    max = values[i];
                }
                count++;
            }
        }
        float[] errs = Errors.populateYErrorsBySqrt(values);

        HistogramHolder hist = Histogram.createSimpleHistogram(0, max, nBins,
            values, errs);

        int yMaxIdx = MiscMath.findYMaxIndex(hist.getYHist());

        if (yMaxIdx == -1) {
            //should not happen
            throw new IllegalStateException("Error in algorithm for extreme data case. "
                + " yMaxIdx=" + yMaxIdx);
        } else if (yMaxIdx == (hist.getXHist().length - 1)) {
            return hist.getXHist()[yMaxIdx];
        }

        double rMode = (hist.getXHist()[yMaxIdx] + hist.getXHist()[yMaxIdx + 1])/2.;

        return rMode;
    }

    /**
     * apply a 2 level threshold hysteresis filter to the image and use an
     * association radius of 2 and value > t1 from any pixel within a
     * radius of 2 of a pixel with value > t2.
     *
     * @param img the phase congruence image
     * @param t1 low threshold
     * @param t2 high threshold
     * @param restore
     * @return
     */
    int[][] applyHysThresh(double[][] img, double[][] pa, double t1, double t2,
        boolean restore) {

        // note that the kovesi code uses the octave bwselect and bwfill,
        // which results in points > t2 which is thresholding just for
        // the high value.
        // here instead will use an adaptive 2-layer threshold for values > t1.

        int w = img.length;
        int h = img[0].length;
        int n = w * h;

        if (w < 3 || h < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        double[][] threshImg = null;
        if (useAdaptiveThreshold) {
            AdaptiveThresholding th = new AdaptiveThresholding();
            threshImg = th.createAdaptiveThresholdImage(img, 15, 0.2);
            if (doPlot) {//DEBUG
                double[][] imgCp = imageProcessor.copy(img);
                double[][] imgCp2 = imageProcessor.copy(threshImg);
                for (int i = 0; i < w; ++i) {
                    for (int j = 0; j < h; ++j) {
                        double t = threshImg[i][j];
                        if (imgCp[i][j] > t && imgCp[i][j] > t1) {
                            imgCp[i][j] = 255.;
                        } else {
                            imgCp[i][j] = 0;
                        }
                        imgCp2[i][j] *= 255.;
                    }
                }
                MiscDebug.writeImage(imgCp, "img_a_thresholded_.png");
                MiscDebug.writeImage(imgCp2, "img_a_adaptive_threshold_.png");
            }
        }

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        double tHigh = t2;
        double tLow = t1;

        System.out.println("tHigh=" + tHigh +
            " replace with adaptive threshold image = " + (threshImg != null));

        int[][] img2 = new int[w][];
        for (int i = 0; i < w; ++i) {
            img2[i] = new int[h];
        }

        // store pixels w/ v > tHigh
        // and store any of it's neighbors w/ v > tLow

        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {

                double v = img[x][y];

                double tHigh0, tLow0;
                if (threshImg != null) {
                    tHigh0 = threshImg[x][y];
                    tLow0 = tHigh0/2;
                } else {
                    tHigh0 = tHigh;
                    tLow0 = tLow;
                }

                if (v < tHigh0 || v < tLow) {
                    continue;
                }

                img2[x][y] = 255;

                // store any adjacent w/ v > tLow
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }
                    double v2 = img[x2][y2];
                    if (v2 > tLow0) {
                        img2[x2][y2] = 255;
                    }
                }
            }
        }

        return img2;
    }

    private void DEBUG(Complex[][] tmp, String label) {

        try {
            algorithms.util.PolygonAndPointPlotter plotter
                = new algorithms.util.PolygonAndPointPlotter();

            int nc = tmp[0].length;
            int nr = tmp.length;

            float[] x = new float[nc];
            for (int ii = 0; ii < nc; ++ii) {
                x[ii] = ii;
            }
            float[] y = new float[nc];
            float[] xPolygon = null;
            float[] yPolygon = null;

            // plot rows 0.25*nRows, 0.5*nRows, and 0.75*nRows
            for (int nf = 1; nf < 4; nf++) {
                int rowNumber = (int) (((float) nf) * 0.25f * nr);
                for (int ii = 0; ii < nc; ++ii) {
                    y[ii] = (float) tmp[rowNumber][ii].re();
                }
                float minY = MiscMath.findMin(y);
                float maxY = MiscMath.findMax(y);
                plotter.addPlot(-1, nc + 1, minY, maxY, x, y, xPolygon,
                    yPolygon, label + " row=" + rowNumber + " REAL");
            }

            // do same for complex
            for (int nf = 1; nf < 4; nf++) {
                int rowNumber = (int) (((float) nf) * 0.25f * nr);
                for (int ii = 0; ii < nc; ++ii) {
                    y[ii] = (float) tmp[rowNumber][ii].im();
                }
                float minY = MiscMath.findMin(y);
                float maxY = MiscMath.findMax(y);
                plotter.addPlot(-1, nc + 1, minY, maxY, x, y, xPolygon,
                    yPolygon, label + " row=" + rowNumber + " IMAGINARY");
            }

            plotter.writeFile();
        } catch (Exception e) {
        }
        int z = 1;
    }

    private void DEBUG(double[][] tmp, String label, long fileNumber) {

        try {
            algorithms.util.PolygonAndPointPlotter plotter
                = new algorithms.util.PolygonAndPointPlotter();

            int nc = tmp[0].length;
            int nr = tmp.length;

            float[] x = new float[nc];
            for (int ii = 0; ii < nc; ++ii) {
                x[ii] = ii;
            }
            float[] y = new float[nc];
            float[] xPolygon = null;
            float[] yPolygon = null;

            float minY = (float)MiscMath.findMin(tmp);
            float maxY = (float)MiscMath.findMax(tmp);

            // plot rows 0.25*nRows, 0.5*nRows, and 0.75*nRows
            for (int nf = 1; nf < 4; nf++) {
                int rowNumber = (int) (((float) nf) * 0.25f * nr);
                for (int ii = 0; ii < nc; ++ii) {
                    y[ii] = (float) tmp[rowNumber][ii];
                }
                plotter.addPlot(-1, nc + 1, minY, maxY, x, y, xPolygon,
                    yPolygon, label + " row=" + rowNumber);
            }
            plotter.writeFile(fileNumber);
            MiscDebug.writeImage(tmp, label.replaceAll("\\s", "") + ".png");
        } catch (Exception e) {
        }
        int z = 1;
    }

    /**
     * for single purpose use to examine the monogenic transform
     * in domain space before inverse transform to frequency space.
     * @param img
     * @param minWavelength
     * @param mult
     * @param sigmaOnf
     * @return
     */
    public static double[][] logGaborFreqDomainFilter(GreyscaleImage img,
        final int minWavelength, final float mult, final float sigmaOnf,
        int s) {

        int noiseMethod = -1;

        int nCols = img.getWidth();
        int nRows = img.getHeight();

        FilterGrid fg = new FilterGrid();
        FilterGrid.FilterGridProducts fgProducts = fg.filtergrid(nRows, nCols);
        fgProducts.getRadius()[0][0] = 1;
        double[][] radius = fgProducts.getRadius();

        LowPassFilter lpFilter = new LowPassFilter();
        double[][] lp = lpFilter.lowpassfilter(nRows, nCols, 0.45f, 15);

        // keeping taus in case need to increase noise estimate
        double logGaborDenom = 2. * Math.pow(Math.log(sigmaOnf), 2);

        // Centre frequency of filter.
        double wavelength = minWavelength * Math.pow(mult, s);

        double fo = 1.0/wavelength;

        // use notation a[row][col]
        double[][] logGabor = new double[nRows][];
        for (int row = 0; row < nRows; ++row) {
            logGabor[row] = new double[nCols];
            for (int col = 0; col < nCols; ++col) {
                double v = Math.log(radius[row][col]/fo);
                v *= v;
                v = Math.exp(-v/logGaborDenom);
                //logGabor = logGabor.*lp;
                logGabor[row][col] = lp[row][col] * v;
            }
        }
        logGabor[0][0] = 0;

        return logGabor;
    }

    /**
     * for single purpose use to examine creating a filter that can be
     * applied to other images in frequency domain space before
     * transforming back to spatial domain.
     * @param img
     * @return row major result
     */
    public static Complex[][] createLowPassFreqDomainFilter(GreyscaleImage img) {

        int nCols = img.getWidth();
        int nRows = img.getHeight();

        FilterGrid fg = new FilterGrid();
        FilterGrid.FilterGridProducts fgProducts = fg.filtergrid(nRows, nCols);
        fgProducts.getRadius()[0][0] = 1;
        double[][] radius = fgProducts.getRadius();

        LowPassFilter lpFilter = new LowPassFilter();
        double[][] lp = lpFilter.lowpassfilter(nRows, nCols, 0.45f, 15);


        //Periodic Fourier transform of image, using default normalization
        // perfft2 results use notation a[row][col]
        PeriodicFFT perfft2 = new PeriodicFFT();
        //IM = perfft2(im);                   %
        //S, P, s, p where S = FFT of smooth, P = FFT of periodic, s=spatial smooth, p = spatial p
        Complex[][][] perfResults = perfft2.perfft2(img, false);
        Complex[][] capIm = perfResults[1];

        ComplexModifiable sumAbs = new ComplexModifiable(0, 0);
        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {
                capIm[row][col] = capIm[row][col].times(lp[row][col]);
                sumAbs.plus(capIm[row][col]);
            }
        }

        //normalize this to sum of zero and unit standard deviation
        sumAbs.times(1./((double)nCols * nRows));
        Complex mean = new Complex(sumAbs.re(), sumAbs.im());
        Complex factor = new Complex(sumAbs.re(), sumAbs.im());
        factor = factor.times(1./Math.sqrt(2));

        if (sumAbs.abs() != 0) {
            for (int row = 0; row < nRows; ++row) {
                for (int col = 0; col < nCols; ++col) {
                    Complex v = capIm[row][col].minus(mean);
                    capIm[row][col] = v.times(factor);
                }
            }
        }

        return capIm;
    }
}
