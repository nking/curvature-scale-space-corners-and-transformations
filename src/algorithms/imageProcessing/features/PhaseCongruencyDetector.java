package algorithms.imageProcessing.features;

import algorithms.compGeometry.HoughTransform;
import algorithms.imageProcessing.EdgeExtractorSimple;
import algorithms.imageProcessing.FilterGrid;
import algorithms.imageProcessing.FilterGrid.FilterGridProducts;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.LowPassFilter;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.MorphologicalFilter;
import algorithms.imageProcessing.NonMaximumSuppression;
import algorithms.imageProcessing.PeriodicFFT;
import algorithms.imageProcessing.PostLineThinnerCorrections;
import algorithms.imageProcessing.ZhangSuenLineThinner;
import algorithms.imageProcessing.scaleSpace.CSSCornerMaker;
import algorithms.misc.Complex;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.CornerArray;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

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
 * Currently, the corners are produced using curvature scale curvature 
 * calculations for the extracted edge points.
 * Alternative corner methods may be offered in the future.
 * Note that edge detector follows the codes referenced below.
 * Also note that the Peter Kovesi implementation of Phase Congruency detector 
 * before this version which uses monogenic filters, produced a corner map using
 * minimum moments from the phase congruency energy over 6 angles.  That method
 * may be implemented separately from this class in the future because
 * the corners produced from it look very good (but are at the expense of a
 * longer runtime).
 * Note, line drawings and images which are solid shapes are probably best
 * handled by the CannyEdgeDetectorFilterAdpative with setToUseLineDrawingMode().
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
 */
public class PhaseCongruencyDetector {
    
    final private static double epsilon = 1E-4;
    
    private boolean determineCorners = false;
    
    public void setToCreateCorners() {
        this.determineCorners = true;
    }
    
    /**
     * <pre>
     * use the phase congruency method of transformations to create
     * edge maps, orientation, phase angle and a suggested threshold.
     * The default values are:
        int nScale = 5;        
        int minWavelength = 3;        
        float mult = 2.1f;
        float sigmaOnf = 0.55f;
        int k = 2;
        float cutOff = 0.5f; 
        float g = 10;
        float deviationGain = 1.5f;
        int noiseMethod = -1;
     * </pre>
     * @param img
     * @return 
       <pre>
       NOTE: the return products use notation a[row][col]
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
    
        int nScale = 5;        
        int minWavelength = 3;        
        float mult = 2.1f;
        float sigmaOnf = 0.55f;
        int k = 2;
        float cutOff = 0.5f; 
        float g = 10;
        float deviationGain = 1.5f;
        int noiseMethod = -1;
        double tLow = 0.1;
        double tHigh = 0.3;
        boolean increaseKIfNeeded = true;
        
        return phaseCongMono(img, nScale, minWavelength, mult, sigmaOnf, k, 
            increaseKIfNeeded, cutOff, g, deviationGain, noiseMethod, 
            tLow, tHigh);
    }
    
    /**
     * 
     * @param img
     * @param nScale number of wavelet scales.  a lower value reveals more fine 
     * scale features.
     * @param k number of standard deviations of the noise energy beyond the 
     * mean at which we set the noise threshold point.  You may want to vary this
       up to a value of 10 or 20 for noisy images.
       @param increaseKIfNeeded if the number of points in the thinned pc
       * is very high, this allows k to be increased, pc to be recalculated,
       * and a result of smaller number of points in the thinned image.
     * @return 
       <pre>
       NOTE: the return products use notation a[row][col]
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
    public PhaseCongruencyProducts phaseCongMono(GreyscaleImage img,
        final int nScale, int k, boolean increaseKIfNeeded) {
        
        int minWavelength = 3;        
        float mult = 2.1f;
        float sigmaOnf = 0.55f;
        float cutOff = 0.5f; 
        float g = 10;
        float deviationGain = 1.5f;
        int noiseMethod = -1;
        double tLow = 0.1;
        double tHigh = 0.3;
        
        return phaseCongMono(img, nScale, minWavelength, mult, sigmaOnf, k, 
            increaseKIfNeeded, cutOff, g, deviationGain, noiseMethod, tLow, tHigh);
    }
    
    /**
     * an edge detector based upon the principal moments of phase congruency 
     * to create an edge operator that is highly localized and has responses 
     * that are invariant to image contrast. 
     * @param img
     * @param nScale number of wavelet scales.  a lower value reveals more fine 
     * scale features.
     * @param minWavelength wavelength of smallest scale filter
     * @param mult scaling factor between successive filters
     * @param sigmaOnf ratio of standard deviation of Gaussian describing the 
     * log Gabor's filter's transfer function in the frequency domain to the 
     * filter center frequency
     * @param k number of standard deviations of the noise energy beyond the 
     * mean at which we set the noise threshold point.  You may want to vary this
       up to a value of 10 or 20 for noisy images.
       @param increaseKIfNeeded if the number of points in the thinned pc
       * is very high, this allows k to be increased, pc to be recalculated,
       * and a result of smaller number of points in the thinned image.
     * @param cutOff The fractional measure of frequency spread below which phase 
     * congruency values get penalized
     * @param g Controls the sharpness of the transition in the sigmoid function 
     * used to weight phase congruency for frequency spread. 
     * @param deviationGain Amplification to apply to the calculated phase 
     * deviation result.  Increasing this sharpens the edge responses, but can 
     * also attenuate their magnitude if the gain is too large.  Sensible values 
     * to use lie in the range 1-2.
     * @param noiseMethod Parameter specifies method used to determine noise 
     * statistics: -1 use median of smallest scale filter responses; 
     * -2 use mode of smallest scale filter responses;
     * 0 turns off all noise compensation; and
     * > 0 use noiseMethod value as the fixed noise threshold;
     * @param tLow the low threshold fraction of 1.  this is usually 0.1.
     * @param tHigh the high threshold fraction of 1.  this is usually 0.3.
     * 
     * @return 
       <pre>
       NOTE: the return products use notation a[row][col]
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
    public PhaseCongruencyProducts phaseCongMono(GreyscaleImage img,
        final int nScale, final int minWavelength, final float mult,
        final float sigmaOnf, int k, final boolean increaseKIfNeeded,
        final float cutOff,
        final float g, final float deviationGain, final int noiseMethod,
        final double tLow, final double tHigh) {
              
        if (increaseKIfNeeded && (noiseMethod >= 0)) {
            throw new IllegalArgumentException(
                "if noiseMethod > -1, there is no dependency on k");
        }
        
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
        }  // end for each scale
        
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
        // fit is phase angle
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
                orientation[row][col] = Math.atan2(-sumH2[row][col], sumH1[row][col]);
                if (orientation[row][col] < 0) {
                    orientation[row][col] += Math.PI;
                }
                // orientation values now range 0 - pi
                // Quantize to 0 - 180 degrees (for NONMAXSUP)
                orientation[row][col] = (int)(orientation[row][col]*180./Math.PI);
                
                //Feature type - a phase angle -pi/2 to pi/2.
                double h1Sq = sumH1[row][col];
                h1Sq *= h1Sq;
                double h2Sq = sumH2[row][col];
                h2Sq *= h2Sq;
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
        } 
        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {
                
                double eDiv = Math.acos(energy[row][col]/(sumAn[row][col] + epsilon));
                
                pc[row][col] = weight[row][col] 
                    * Math.max(1. - deviationGain * eDiv, 0)
                    * Math.max(energy[row][col] - threshold, 0)
                    / (energy[row][col] + epsilon);
            }
        }
                
        PhaseCongruencyProducts products = new PhaseCongruencyProducts(pc, 
            orientation, ft, threshold);      
        
        NonMaximumSuppression ns = new NonMaximumSuppression();
        
        double[][] thinnedPC = ns.nonmaxsup(products.getPhaseCongruency(), 
            products.getOrientation(), 1.2, new HashSet<PairInt>());  
        
{
GreyscaleImage tmp = new GreyscaleImage(nCols, nRows);
for (int i = 0; i < nCols; ++i) {
    for (int j = 0; j < nRows; ++j) {
        tmp.setValue(i, j, (int)(255.*thinnedPC[j][i]));
    }
}
MiscDebug.writeImage(tmp, "_THINNED_PC_");
}        
        
        // NOTE: limit does not scale with resolution, so user may want to
        // pre-process images to a common resolution or size depending upon goal
        if (increaseKIfNeeded) {
            
            // count number edge points
            int nEdgePoints = countEdgePoints(products, thinnedPC, tLow, tHigh);
            
            System.out.println("nEdgePoints=" + nEdgePoints);
        
            int limit = 15000;
            int lastK = k;
            while (nEdgePoints > limit) {
                
                // k can be as high as 20
                int deltaK = Math.round(0.333f * (20 - lastK));
                if (deltaK == 0) {
                    deltaK = 1;
                }
                k = lastK + deltaK;
                
                if (k >= 20) {
                    break;
                }
                
                lastK = k;
                
                if (noiseMethod < 0) {
                    
                    double totalTau = tau * (1. - Math.pow((1./mult), nScale))/(1. - (1./mult));
                    double EstNoiseEnergyMean = totalTau * Math.sqrt(Math.PI/2.);
                    double EstNoiseEnergySigma = totalTau * Math.sqrt((4. - Math.PI)/2.);
                    threshold = Math.max(EstNoiseEnergyMean 
                        + ((float)k) * EstNoiseEnergySigma, epsilon);
                }
                
                for (int row = 0; row < nRows; ++row) {
                    for (int col = 0; col < nCols; ++col) {
                        double eDiv = Math.acos(energy[row][col] / (sumAn[row][col] + epsilon));
                        pc[row][col] = weight[row][col]
                            * Math.max(1. - deviationGain * eDiv, 0)
                            * Math.max(energy[row][col] - threshold, 0)
                            / (energy[row][col] + epsilon);
                    }
                }

                products = new PhaseCongruencyProducts(pc, orientation, ft, 
                    threshold);

                thinnedPC = ns.nonmaxsup(products.getPhaseCongruency(), 
                    products.getOrientation(), 1.2,  new HashSet<PairInt>());

                nEdgePoints = countEdgePoints(products, thinnedPC, tLow, tHigh);
                
                System.out.println("nEdgePoints=" + nEdgePoints);
            }
        }
        
        createEdges(products, thinnedPC, tLow, tHigh);
                
{        
GreyscaleImage tmp = new GreyscaleImage(nCols, nRows);
for (int i = 0; i < nRows; ++i) {
    for (int j = 0; j < nCols; ++j) {
        if (products.getThinned()[i][j] > 0) {
            tmp.setValue(j, i, 255);
        }
    }
}
MiscDebug.writeImage(tmp, "_EDGES_");
}

        if (determineCorners) {
            
            calculateHoughTransforms(products);
                                                
            calculateCorners(products, img);
        }
        
        return products;
    }
    
    protected void calculateCorners(PhaseCongruencyProducts products,
        GreyscaleImage img) {
        
        // only computing curvature for edge points
        int[][] edgeImage = products.getThinned();
        
        int nRows = edgeImage.length;
        int nCols = edgeImage[0].length;        
        
        EdgeExtractorSimple edgeExtractor = new EdgeExtractorSimple(edgeImage);
        edgeExtractor.extractEdges();
        
        Set<PairInt> junctions = edgeExtractor.getJunctions();
        
        List<PairIntArray> theEdges = edgeExtractor.getEdges();
        
        // ---- placing these in coordinate reference frame of images -------
        Set<PairInt> tmp = new HashSet<PairInt>();
        for (PairInt p : junctions) {
            tmp.add(new PairInt(p.getY(), p.getX()));
        }
        junctions = tmp;
        
        for (PairIntArray edge : theEdges) {
            for (int i = 0; i < edge.getN(); ++i) {
                int x = edge.getX(i);
                int y = edge.getY(i);
                edge.set(i, y, x);
            }
        }
        
        CSSCornerMaker cornerMaker = new CSSCornerMaker(img.getWidth(), img.getHeight());
        cornerMaker.doNotStoreCornerRegions();
        List<CornerArray> cornerList =
            cornerMaker.findCornersInScaleSpaceMaps(theEdges);
        
        /*
        filter out small curvature:
        see line 64 of comments in CornerRegion.java.
        for curvature smaller than 0.2 won't see changes in slope in the
             neighboring 2 points on either side.
        */
        Map<PairInt, Float> cornerMap = new HashMap<PairInt, Float>();
        for (int i = 0; i < cornerList.size(); ++i) {
            CornerArray ca = cornerList.get(i);
            for (int idx = 0; idx < ca.getN(); ++idx) {
                float curvature = ca.getCurvature(idx);
                
                if (Math.abs(curvature) > 0.05) {// should not set above 0.07...
                                        
                    cornerMap.put(
                        new PairInt(Math.round(ca.getX(idx)),
                        Math.round(ca.getY(idx))), 
                        Float.valueOf(ca.getCurvature(idx)));
                }
            }
        }        
        
        // theEdges, corners, and junctions are now in reference frame of
        // GreyscaleImage instance
        cornerMaker.useHoughTransformationToFilterCornersForOrdered(theEdges, 
            cornerMap, junctions, products.getHoughLines(),
            nCols, nRows);        
        
        Set<PairInt> outputCorners = new HashSet<PairInt>(cornerMap.keySet());
        outputCorners.addAll(junctions);
        
        // this sets the edges, junctions and corners in the GreyscaleImage reference frame
        products.setEdges(theEdges);
        products.setJunctions(junctions);
        products.setCorners(outputCorners);
        
        /*        
        Image imgCp = img.copyToColorGreyscale();
        ImageIOHelper.addAlternatingColorCurvesToImage(theEdges, imgCp, 0, 0, 0);
        ImageIOHelper.addCurveToImage(outputCorners, imgCp, 3, 255, 0, 0);
        MiscDebug.writeImage(imgCp, "_CORNERS_0");
        */
    }

    private double[][] copy(double[][] a) {
        double[][] cp = new double[a.length][];
        for (int i = 0; i < a.length; ++i) {
            cp[i] = Arrays.copyOf(a[i], a[i].length);
        }
        return cp;
    }

    private Set<PairInt> extractNonZeroPoints(int[][] binaryImage) {
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        for (int i0 = 0; i0 < binaryImage.length; ++i0) {
            for (int i1 = 0; i1 < binaryImage[0].length; ++i1) {
                if (binaryImage[i0][i1] > 0) {
                    points.add(new PairInt(i0, i1));
                }
            }
        }
        
        return points;
    }

    /**
     * creating edges using thinned phase angle image.
     * @param products
     * @param tLow
     * @param tHigh 
     */
    private void createEdges(PhaseCongruencyProducts products, 
        double[][] thinnedPC, double tLow, double tHigh) {
        
        double[][] phaseAngle = products.getPhaseAngle();
        
        int nRows = phaseAngle.length;
        int nCols = phaseAngle[0].length;
 
        Set<PairInt> brightLinePoints = new HashSet<PairInt>();
        Set<PairInt> darkLinePoints = new HashSet<PairInt>();
        Set<PairInt> stepPoints = new HashSet<PairInt>();
                    
        double piDiv2 = Math.PI/2.;
        double piDiv4 = piDiv2/2.;
        
        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {
                
                if (thinnedPC[row][col] < tLow) {
                    continue;
                }
                // placing values closer to -pi/2 than 0 into img0,
                //         values closer to +p1/2 than 0 into img2,
                //         else they are closer to 0 and in img1
                // TODO: when see the results, might consider a smaller tolerance
                //        than piDiv4 for steps
                double v = phaseAngle[row][col];
                if (v < 0) {
                    if (Math.abs(v - -piDiv2) < piDiv4) {
                        brightLinePoints.add(new PairInt(row, col));
                    } else {
                        stepPoints.add(new PairInt(row, col));
                    }
                } else if (v == 0) {
                    stepPoints.add(new PairInt(row, col));
                } else {
                    if (Math.abs(v - piDiv2) < piDiv4) {
                        stepPoints.add(new PairInt(row, col));
                    } else {
                        stepPoints.add(new PairInt(row, col));
                    }
                }
            }
        }
        
        PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
        pltc.correctForIsolatedPixels(brightLinePoints, nRows, nCols);
        pltc.correctForIsolatedPixels(darkLinePoints, nRows, nCols);
        
        stepPoints.addAll(brightLinePoints);
        stepPoints.addAll(darkLinePoints);
        
        pltc.correctForIsolatedPixels(stepPoints, nRows, nCols);
        pltc.correctForLine2SpurHoriz(stepPoints, nRows, nCols);
        pltc.correctForLine2SpurVert(stepPoints, nRows, nCols);
        pltc.correctForLineHatHoriz(stepPoints, nRows, nCols);
        pltc.correctForLineHatVert(stepPoints, nRows, nCols);
        
        
        // ---- complete the edges where possible ----
        double[][] pc = products.getPhaseCongruency();
        Set<PairInt> gaps = pltc.findGapsOf1(stepPoints, nRows, nCols);
        for (PairInt p : gaps) {
            int x = p.getX();
            int y = p.getY();
            if (pc[x][y] > 0) {
                stepPoints.add(new PairInt(x, y));
            }
        }
        
        int[][] thinned = new int[nRows][nCols];
        for (int i = 0; i < nRows; ++i) {
            thinned[i] = new int[nCols];
        }
        for (PairInt p : stepPoints) {
            thinned[p.getX()][p.getY()] = 1;
        }
        stepPoints = null;
        gaps = null;
        
        MorphologicalFilter mFilter = new MorphologicalFilter();
        int[][] skel = mFilter.bwMorphThin(thinned, Integer.MAX_VALUE);

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
        pltc.correctForLineHatHoriz(points, nRows, nCols);
        pltc.correctForLineHatVert(points, nRows, nCols);
        
        // there are a very small number of clumps thicker than 1 pixel.
        ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
        lt.applyLineThinner(points, 0, nRows - 1, 0, nCols - 1);
        
        for (int i = 0; i < nRows; ++i) {
            Arrays.fill(thinned[i], 0);
        }
        for (PairInt p : points) {
            thinned[p.getX()][p.getY()] = 255;
        }
        
        products.setThinnedImage(thinned);
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

    /**
     * calculate the hough transform lines, and use them to thin the 
     * staircases within inclined lines.
     * @param products 
     */
    private void calculateHoughTransforms(PhaseCongruencyProducts products) {

        int[][] thinned = products.getThinned();

        Set<PairInt> points = new HashSet<PairInt>();
        for (int i = 0; i < thinned.length; ++i) {
            for (int j = 0; j < thinned[i].length; ++j) {
                if (thinned[i][j] > 0) {
                    points.add(new PairInt(i, j));
                }
            }
        }
        
        // ----- find lines w/ hough transform, then thin line staircase with it ---
        Set<PairInt> pointCp = new HashSet<PairInt>(points);
        HoughTransform ht = new HoughTransform();
        Map<Set<PairInt>, PairInt> lines = ht.findContiguousLines(points, 3);

        PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
        pltc.thinLineStaircases(lines, points, thinned.length, thinned[0].length);
        
        pointCp.removeAll(points);        
        
        for (PairInt p : pointCp) {
            Set<PairInt> line = null;
            for (Set<PairInt> hLine : lines.keySet()) {
                if (hLine.contains(p)) {
                    line = hLine;
                    break;
                }
            }
            if (line != null) {
                line.remove(p);
            }
        }
        
        // ------ put hough lines into reference frame of GreyscaleImage instance
        Map<Set<PairInt>, PairInt> transformedLines = new HashMap<Set<PairInt>, PairInt>();
        for (Entry<Set<PairInt>, PairInt> entry : lines.entrySet()) {
            Set<PairInt> line = entry.getKey();
            Set<PairInt> transformedLine = new HashSet<PairInt>();
            for (PairInt p : line) {
                int x = p.getX();
                int y = p.getY();
                transformedLine.add(new PairInt(y, x));
            }
            transformedLines.put(transformedLine, entry.getValue());
        }
        
        products.setHoughLines(transformedLines);        
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
                
        private List<PairIntArray> edgeList = null;
        
        private Set<PairInt> junctions = null;
        
        private Set<PairInt> corners = null;
        
        /**
         * a map with key being a hough line set of points and value being
         * the theta and radius for the hough line.  Note that the coordinates
         * in the key are using the same notation as the thinnedImage, that is,
         * a[row][col], but pairint.x is row, and pairint.y is col.
         */
        private Map<Set<PairInt>, PairInt> houghLines = null;
        
        public PhaseCongruencyProducts(double[][] pc, double[][] or, 
            double[][] ft, double thr) {
            this.phaseCongruency = copy(pc);
            this.orientation = copy(or);
            this.phaseAngle = copy(ft);
            this.threshold = thr;
        }
        
        /**
         * set the thinned phase congruence image, a.k.a. the edge image.
         * @param thImg 
         */
        public void setThinnedImage(int[][] thImg) {
            thinned = copy(thImg);
        }
        
        /**
         * get the thinned phase congruence image, a.k.a. the edge image.
         * Note that the array is accessed as a[row][column].
         * @@return edgeImg 
         */
        public int[][] getThinned() {
            return thinned;
        }

        /**
         * return the gradient image produced by phase congruency as a double
         * array of values in range 0 to 1.0.
         * Note that the array is accessed as a[row][column].
         * @return the phaseCongruency
         */
        public double[][] getPhaseCongruency() {
            return phaseCongruency;
        }

        /**
         * return the orientation image.
         * Note that the array is accessed as a[row][column].
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
         * Note that the array is accessed as a[row][column].
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

        /**
         * set the edge list extracted from the thinned image using coordinates
         * that are in the reference frame of the GreyscaleIamge instance.
         * Note that the closed curves are converted to instances of
         * PairIntArrayWithColor.
         * @param theEdgeList 
         */
        private void setEdges(List<PairIntArray> theEdgeList) {
            
            this.edgeList = new ArrayList<PairIntArray>();
            
            MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
            
            for (PairIntArray p : theEdgeList) {
                
                PairIntArray cp = p.copy();
                if (curveHelper.isAdjacent(p, 0, p.getN() - 1, 1)) {
                    cp = new PairIntArrayWithColor(cp);
                }
                this.edgeList.add(cp);
            }
        }

        /**
         * set the junction points found in the edge list.  Note that the
         * points should be using coordinates that are in the reference frame 
         * of the GreyscaleIamge instance.
         * @param theJunctions 
         */
        private void setJunctions(Set<PairInt> theJunctions) {
            this.junctions = new HashSet<PairInt>(theJunctions);
        }

        /**
         * set the corners found in the edge list.  Note that the points should 
         * be using coordinates that are in the reference frame of the 
         * GreyscaleIamge instance.
         * @param theCorners 
         */
        private void setCorners(Set<PairInt> theCorners) {
            this.corners = new HashSet<PairInt>(theCorners);
        }

        /**
         * get the list of extracted edges as points that are in the reference 
         * frame of the GreyscaleImage instance.
         * @return the edgeList
         */
        public List<PairIntArray> getEdgeList() {
            return edgeList;
        }

        /**
         * get the set of junction points found within the edges in the reference
         * frame of the GreyscaleIamge instance.
         * @return the junctions
         */
        public Set<PairInt> getJunctions() {
            return junctions;
        }

        /**
         * get the corners found in the edges as set of points that are in the
         * reference frame of the GreyscaleImage instance.
         * @return the corners
         */
        public Set<PairInt> getCorners() {
            return corners;
        }

        /**
           a map with key being a hough line set of points and value being
         * the theta and radius for the hough line.  Note that the
         * setter has placed the coordinates into the coordinate reference
         * frame of the GreyscaleImage instance.
         * @param lines 
         */
        public void setHoughLines(Map<Set<PairInt>, PairInt> lines) {
            
             this.houghLines = new HashMap<Set<PairInt>, PairInt>(lines);
        }
        
        /**
         * a map with key being a hough line set of points and value being
         * the theta and radius for the hough line.  Note that the coordinates
         * in the key are using the same notation GreyscaleImage instance, 
         * that x and y are the same as in the image
         * @return 
         */
        public Map<Set<PairInt>, PairInt> getHoughLines() {
            return houghLines;
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
     * apply a 2 level threshold hysteresis filter to the image and use an
     * association radius of 2 and value > t1 from any pixel within a
     * radius of 2 of a pixel with value > t2.
     * 
     * @param img the phase congruence image
     * @param t1 low threshold
     * @param t2 high threshold
     * @param extendAssocRadiusTo2 if true, searches the neighbors of
     *    neighbors for a strong point
     * @return 
     */
    int[][] applyHysThresh(double[][] img, double t1, double t2, 
        boolean extendAssocRadiusTo2) {
        
        // note that the kovesi code uses the octave bwselect and bwfill,
        // which results in points > t2 which is thresholding just for 
        // the high value.
        
        // so will instead adapt my canny edge detector 2-layer threshold
        // here.
        
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
                
                if (extendAssocRadiusTo2) {
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
        }
        
        return img2;
    }
    
    /*
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
    
    private void DEBUG(double[][] tmp, String label) {
       
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
                    y[ii] = (float) tmp[rowNumber][ii];
                }
                float minY = MiscMath.findMin(y);
                float maxY = MiscMath.findMax(y);
                plotter.addPlot(-1, nc + 1, minY, maxY, x, y, xPolygon,
                    yPolygon, label + " row=" + rowNumber);
            }
            plotter.writeFile();
        } catch (Exception e) {
        }
        int z = 1;
    }
    */
    
}
