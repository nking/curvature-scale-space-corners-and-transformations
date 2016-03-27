package algorithms.imageProcessing.features;

import algorithms.imageProcessing.FilterGrid;
import algorithms.imageProcessing.FilterGrid.FilterGridProducts;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.LowPassFilter;
import algorithms.imageProcessing.PeriodicFFT;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.misc.Complex;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import java.util.Arrays;

/**
 Not yet tested or ready for use.
 
 adapted from 
  http://www.peterkovesi.com/matlabfns/PhaseCongruency/phasecongmono.m
  which has copyright:
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
 
 * useful also in looking at the python phasepack port by Alistair Muldal
 *  http://pydoc.net/Python/phasepack/1.4/phasepack.phasecongmono/
 * which has the following copyright:
 * # MIT License:

# Permission is hereby  granted, free of charge, to any  person obtaining a
# copy of this software and associated  documentation files (the "Software"),
# to deal in the Software without restriction, subject to the following
# conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# The software is provided "as is", without warranty of any kind.
# 
* 
 There are potentially many arguments, here is the full usage:
%
%   [PC or ft T] =  ...
%                phasecongmono(im, nscale, minWaveLength, mult, ...
%                         sigmaOnf, k, cutOff, g, deviationGain, noiseMethod)
%
% However, apart from the image, all parameters have defaults and the
% usage can be as simple as:
%
%    phaseCong = phasecongmono(im);
% 
% Arguments:
%              Default values      Description
%
%    nscale           4    - Number of wavelet scales, try values 3-6
%                            A lower value will reveal more fine scale
%                            features. A larger value will highlight 'major'
%                            features.
%    minWaveLength    3    - Wavelength of smallest scale filter.
%    mult             2.1  - Scaling factor between successive filters.
%    sigmaOnf         0.55 - Ratio of the standard deviation of the Gaussian 
%                            describing the log Gabor filter's transfer function 
%                            in the frequency domain to the filter center frequency.
%    k                3.0  - No of standard deviations of the noise energy beyond
%                            the mean at which we set the noise threshold point.
%                            You may want to vary this up to a value of 10 or
%                            20 for noisy images 
%    cutOff           0.5  - The fractional measure of frequency spread
%                            below which phase congruency values get penalized.
*    g                10   - Controls the sharpness of the transition in
%                            the sigmoid function used to weight phase
%                            congruency for frequency spread.                        
%    deviationGain    1.5  - Amplification to apply to the calculated phase
%                            deviation result. Increasing this sharpens the
%                            edge responses, but can also attenuate their
%                            magnitude if the gain is too large.  Sensible
%                            values to use lie in the range 1-2.
%    noiseMethod      -1   - Parameter specifies method used to determine
%                            noise statistics. 
%                              -1 use median of smallest scale filter responses
%                              -2 use mode of smallest scale filter responses
%                               0+ use noiseMethod value as the fixed noise threshold 
%                            A value of 0 will turn off all noise compensation.
%
% Returned values:
*    PC         - Phase congruency indicating edge significance (values are in range 0 to 1.)
%    or         - Orientation image in integer degrees 0-180,
%                 positive anticlockwise.
%                 0 corresponds to a vertical edge, 90 is horizontal.
%    ft         - Local weighted mean phase angle at every point in the
%                 image.  A value of pi/2 corresponds to a bright line, 0
%                 corresponds to a step and -pi/2 is a dark line.
%    T          - Calculated noise threshold (can be useful for
%                 diagnosing noise characteristics of images).  Once you know
%                 this you can then specify fixed thresholds and save some
%                 computation time.
%
% Notes on specifying parameters:  
%
% The parameters can be specified as a full list eg.
%  >> PC = phasecongmono(im, 5, 3, 2.5, 0.55, 2.0);
%
% or as a partial list with unspecified parameters taking on default values
%  >> PC = phasecongmono(im, 5, 3);
%
% or as a partial list of parameters followed by some parameters specified via a
% keyword-value pair, remaining parameters are set to defaults, for example:
%  >> PC = phasecongmono(im, 5, 3, 'k', 2.5);
% 
% The convolutions are done via the FFT.  Many of the parameters relate to the
% specification of the filters in the frequency plane.  The values do not seem
% to be very critical and the defaults are usually fine.  You may want to
% experiment with the values of 'nscales' and 'k', the noise compensation
% factor.
* 
* Typical sequence of operations to obtain an edge image:
%
%  >> [PC, or] = phasecongmono(imread('lena.tif'));
%  >> nm = nonmaxsup(PC, or, 1.5);   % nonmaxima suppression
%  >> bw = hysthresh(nm, 0.1, 0.3);  % hysteresis thresholding 0.1 - 0.3
%  >> show(bw)
%
% Notes on filter settings to obtain even coverage of the spectrum
% sigmaOnf       .85   mult 1.3
% sigmaOnf       .75   mult 1.6     (filter bandwidth ~1 octave)
% sigmaOnf       .65   mult 2.1  
% sigmaOnf       .55   mult 3       (filter bandwidth ~2 octaves)
%
% Note that better results are achieved using the large bandwidth filters.  
% I generally use a sigmaOnf value of 0.55 or even smaller.
%
% See Also:  PHASECONG, PHASECONG3, PHASESYMMONO, GABORCONVOLVE,
% PLOTGABORFILTERS, FILTERGRID

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
 */
public class PhaseCongruencyDetector {
        
    /**
     * @param img
     * @return NOTE: the return products use notation a[row][col]
     */
    public PhaseCongruencyProducts phaseCongMono(GreyscaleImage img) {
    
        // number of wavelet scales.  a lower value reveals more fine scale features.
        int nScale = 5;
        
        // wavelength of smallest scale filter
        int minWavelength = 3;
        
        // scaling factor between successive filters
        float mult = 2.1f;
        
        // ratio of standard deviation of Gaussian describing the log Gabor's filter's
        // transfer function in the frequency domain to the filter center frquency
        float sigmaOnf = 0.55f;
        
        //number of standard deviations of the noise energy beyond the mean
        // at which we set the noise threshold point.  You may want to vary this
        // up to a value of 10 or 20 for noisy images.
        float k = 2.0f;
        
        // The fractional measure of frequency spread below which phase 
        // congruency values get penalized.
        float cutOff = 0.5f;
        
        //Controls the sharpness of the transition in the sigmoid function used 
        // to weight phase congruency for frequency spread. 
        float g = 10;
        
        //Amplification to apply to the calculated phase deviation result. 
        // Increasing this sharpens the edge responses, but can also attenuate 
        // their magnitude if the gain is too large.  Sensible values to use 
        //  lie in the range 1-2.
        float deviationGain = 1.5f; 
        
        //Parameter specifies method used to determine noise statistics. 
        //     -1 use median of smallest scale filter responses
        //     -2 use mode of smallest scale filter responses
        //      0+ use noiseMethod value as the fixed noise threshold 
        // A value of 0 will turn off all noise compensation.
        int noiseMethod = -1;
        
        final double epsilon = 1E-4;
        
        int nCols = img.getWidth();
        int nRows = img.getHeight();
                
        //Periodic Fourier transform of image
        // perfft2 results use notation a[row][col]
        PeriodicFFT perfft2 = new PeriodicFFT();
        perfft2.setToNotNormalize();
        //IM = perfft2(im);                   % 
        //S, P, s, p where S = FFT of smooth, P = FFT of periodic, s=spatial smooth, p = spatial p
        Complex[][][] perfResults = perfft2.perfft2(img, true);
        Complex[][] capIm = perfResults[1];
                      
        /*
        sumAn  = zeros(rows,cols);          % Matrix for accumulating filter response
                                            % amplitude values.
        sumf   = zeros(rows,cols);                                  
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
        
        double tau = Double.NaN;
        double sqml4 = Math.sqrt(Math.log(4));
        double logGaborDenom = 2. * Math.pow(Math.log(sigmaOnf), 2);
        double threshold = epsilon;
        
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
            // but are by inverse fft so need a combined division here by nRows*nCols
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
            // result needs to be divided by nRows*nCols
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

            for (int row = 0; row < nRows; ++row) {
                for (int col = 0; col < nCols; ++col) {
                    double a = sumAn[row][col]/(maxAN[row][col] + epsilon);
                    width[row][col] = (a - 1.)/dn;
                    double v = Math.exp(g*(cutOff - width[row][col]));
                    weight[row][col] = 1./(1. + v);
                }
            }
        
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

                threshold = Math.max(EstNoiseEnergyMean + k*EstNoiseEnergySigma, epsilon);
            }
            
        } // end for each scale
        
        // uses notation a[row][col]
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
        
        // uses notation a[row][col]
        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {
                orientation[row][col] = Math.atan(-sumH2[row][col]/sumH1[row][col]);
                if (orientation[row][col] < 0) {
                    orientation[row][col] += Math.PI;
                }
                // orientation values now range 0 - pi
                // Quantize to 0 - 180 degrees (for NONMAXSUP)
                orientation[row][col] = Math.floor(orientation[row][col]*180./Math.PI);
                
                //Feature type - a phase angle -pi/2 to pi/2.
                double v1 = sumH1[row][col];
                v1 *= v1;
                double v2 = sumH2[row][col];
                v2 *= v2;
                ft[row][col] = Math.atan2(sumF[row][col], Math.sqrt(v1 + v2));
                
                //overall energy
                double v0 = sumF[row][col];
                v0 *= v0;
                energy[row][col] = Math.sqrt(v0 + v1 + v2);
                
                double eDiv = Math.acos(energy[row][col]/(sumAn[row][col] + epsilon));
                
                pc[row][col] = weight[row][col] 
                    * Math.max(1. - deviationGain * eDiv, 0)
                    * Math.max(energy[row][col] - threshold, 0)
                    / (energy[row][col] + epsilon);
                
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
        
        PhaseCongruencyProducts products = new PhaseCongruencyProducts(pc, 
            orientation, ft, threshold);
        
        return products;
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
        
        public PhaseCongruencyProducts(double[][] pc, double[][] or, 
            double[][] ft, double thr) {
            this.phaseCongruency = pc;
            this.orientation = or;
            this.phaseAngle = ft;
            this.threshold = thr;
        }

        /**
         * @return the phaseCongruency
         */
        public double[][] getPhaseCongruency() {
            return phaseCongruency;
        }

        /**
         * @return the orientation
         */
        public double[][] getOrientation() {
            return orientation;
        }

        /**
         * @return the phaseAngle
         */
        public double[][] getPhaseAngle() {
            return phaseAngle;
        }

        /**
         * @return the threshold
         */
        public double getThreshold() {
            return threshold;
        }
        
    }

    private Complex[][] copy(Complex[][] a) {
        
        Complex[][] b = new Complex[a.length][];
        for (int i = 0; i < b.length; ++i) {
            b[i] = Arrays.copyOf(a[i], a[i].length);
        }
        
        return b;
    }
    
    /**
     * adapted from Kovesis phasecongmono.m as documented in class comments above.
       
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
     * apply a 2 level threshold hysteresis filter to the image.  
     * 
     * 
     * 
     * @param img
     * @param t1
     * @param t2
     * @return 
     */
    int[][] applyHysThresh(double[][] img, double t1, double t2) {
        
        // note that the kovesi code uses the octave bwselect and bwfill,
        // which results in points > t2.
        // so will instead adapt my canny edge detector 2-layer threshold
        // here.
        // if t1 and t2 suggested settings from his usage are not quite
        // right for this, will use the otsu value and a scale factor.
        
        // copying apply2LayerFilter from the adaptive canny directly here,
        // but if t1 and t2 are not producing best results, will edit
        // the other apply2LayerFilter to accept the img here as int[][]
        // instead
        
        int w = img.length;
        int h = img[0].length;
        int n = w * h;
        
        if (w < 3 || h < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }
    
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
            
        double tHigh = t2;
        double tLow = t1;
        
        int[][] img2 = new int[w][];
        for (int i = 0; i < w; ++i) {
            img2[i] = new int[h];
        }
                
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
            
                double v = img[x][y];
            
                if (v < tLow) {
                    continue;
                } else if (v > tHigh) {
                    img2[x][y] = 255;
                    continue;
                }
                
                boolean foundHigh = false;
                boolean foundMid = false;
            
                for (int k = 0; k < dxs.length; ++k) {                
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }
                    double v2 = img[x2][y2];
                    if (v2 > tHigh) {
                        foundHigh = true;
                        break;
                    } else if (v2 > tLow) {
                        foundMid = true;
                    }
                }
                if (foundHigh) {
                    img2[x][y] = 255;
                    continue;
                }
                if (!foundMid) {
                    continue;
                }
                // search the 5 by 5 region for a "sure edge" pixel
                for (int dx = -2; dx <= 2; ++dx) {
                    int x2 = x + dx;
                    if ((x2 < 0) || (x2 > (w - 1))) {
                        continue;
                    }
                    for (int dy = -2; dy <= 2; ++dy) {
                        int y2 = y + dy;
                        if ((y2 < 0) || (y2 > (h - 1))) {
                            continue;
                        }
                        if (x2 == x && y2 == y) {
                            continue;
                        }
                        double v2 = img[x2][y2];
                        if (v2 > tHigh) {
                            img2[x][y] = 255;
                            foundHigh = true;
                            break;
                        }
                    }
                    if (foundHigh) {
                        break;
                    }
                }
            }
        }
        
        return img2;
    }
    
}
