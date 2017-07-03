package algorithms.imageProcessing.features;

import algorithms.imageProcessing.FFTUtil;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.LowPassFilter;
import algorithms.imageProcessing.NonMaximumSuppression;
import algorithms.imageProcessing.PolarFilterGrid;
import algorithms.misc.Complex;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import java.util.Arrays;
import java.util.HashSet;

/**
 * Non-monogenic filter phase congruence edge detector.
 * 
 * Listings of copyrights for the original source codes in languages Matlab and 
 * python follow:
 * 
 adapted from 
  http://www.peterkovesi.com/matlabfns/PhaseCongruency/phasecong2.m
  which has copyright:
  Copyright (c) 1996-2009 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
  
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
*/
public class PhaseCongruencyDetector2 {
    
    final private static double epsilon = 1E-4;

    /**
     * an edge detector based upon the principal moments of phase congruency 
     * to create an edge operator that is highly localized and has responses 
     * that are invariant to image contrast. 
     * @param img
     * @param nScale number of wavelet scales.  a lower value reveals more fine 
     * scale features.
     * @param nOrient = 6.  Number of filter orientations.
     * @param minWavelength wavelength of smallest scale filter
     * @param mult scaling factor between successive filters
     * @param sigmaOnf ratio of standard deviation of Gaussian describing the 
     * log Gabor's filter's transfer function in the frequency domain to the 
     * filter center frequency
     * @param k number of standard deviations of the noise energy beyond the 
     * mean at which we set the noise threshold point.  You may want to vary this
       up to a value of 10 or 20 for noisy images.
      
     * @param cutOff The fractional measure of frequency spread below which phase 
     * congruency values get penalized
     * @param g Controls the sharpness of the transition in the sigmoid function 
     * used to weight phase congruency for frequency spread. 
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
    public PhaseCongruencyProducts phaseCong(GreyscaleImage img,
        final int nScale, final int nOrient, final int minWavelength, final float mult,
        final float sigmaOnf,
        int k, 
        final float cutOff,
        final float g, int noiseMethod, double tLow, double tHigh,
        boolean doStoreConvolution) {
        
        long t0 = System.currentTimeMillis();
              
        int nCols = img.getWidth();
        int nRows = img.getHeight();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        // result is in format a[row][col]
        Complex[][] IMComplex = imageProcessor.create2DFFTWithSwapMajor(
            img, false, true); 
        
        // 3 dimensional sum of energy
        double[][][] EnergyV = new double[3][][];
        double[][] pcSum = new double[nRows][];
        double[][] covx2 = new double[nRows][];
        double[][] covy2 = new double[nRows][];
        double[][] covxy = new double[nRows][];
        for (int row = 0; row < nRows; ++row) {
            pcSum[row] = new double[nCols];
            covx2[row] = new double[nCols];
            covy2[row] = new double[nCols];
            covxy[row] = new double[nCols];
            for (int d = 0; d < 3; ++d) {
                EnergyV[d] = new double[nRows][nCols];
                EnergyV[d][row] = new double[nCols];
            }
        }
        
        // lists to contain convolution results and phase congruency images
        // [nOrientation][nScale][nRows][nCols]
        Complex[][][][] EO = new Complex[nOrient][][][];
        double[][][] PC = new double[nOrient][][];
        for (int or = 9; or < nOrient; ++or) {
            EO[or] = new Complex[nScale][][];
        }
        
        double[] threshold = new double[nOrient];
        
        double tau = noiseMethod;
        // keeping taus in case need to increase noise estimate
        double sqml4 = Math.sqrt(Math.log(4));
        
        // Construct a bank of log-Gabor filters at different spatial scales
        // Filters are constructed in terms of two components.
        // 1) The radial component, which controls the frequency band that the
        //    filter responds to
        // 2) The angular component, which controls the orientation that the filter
        //    responds to.
        // The two components are multiplied together to construct the overall
        // filter.
        // Construct the radial filter components... First construct a low-pass
        // filter that is as large as possible, yet falls away to zero at the
        // boundaries. All log Gabor filters are multiplied by this to ensure no
        // extra frequencies at the 'corners' of the FFT are incorporated as this
        // seems to upset the normalisation process when calculating phase
        // congruency.
        
        // results use notation a[row][col]
        PolarFilterGrid fg = new PolarFilterGrid();
        PolarFilterGrid.FilterGridProducts fgProducts = fg.filtergrid(nRows, nCols);    
        
        double[][] sinTheta = fgProducts.getSinTheta();
        double[][] cosTheta = fgProducts.getCosTheta();
        double[][] radius = fgProducts.getRadius();
        
        // construct a low-pass filter that is as large as possible, yet falls
        // away to zero at the boundaries.  All filters are multiplied by
        // this to ensure no extra frequencies at the 'corners' of the FFT are
        // incorporated as this can upset the normalisation process when
        // calculating phase congruency
        // lp = lowpassfilter([rows,cols],.45,15);    % Radius .4, 'sharpness' 15
        // results use notation a[row][col]
        LowPassFilter lpFilter = new LowPassFilter();
        double[][] lp = lpFilter.lowpassfilter(nRows, nCols, 0.45f, 15);
        
        double[][][] logGabors = new double[nScale][][];
        
        double logGaborDenom = 2. * Math.pow(Math.log(sigmaOnf), 2);
        
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
            
            logGabors[s] = logGabor;
        }
        
        for (int or = 0; or < nOrient; ++or) {
                                
            double angl = or * (Math.PI / nOrient);
            double cosAngl = Math.cos(angl);
            double sinAngl = Math.sin(angl);
            
            //TODO: calc of dTheta can be simplified
            double[][] dTheta = new double[nRows][];
            double[][] spread = new double[nRows][];
            for (int row = 0; row < nRows; ++row) {
                dTheta[row] = new double[nCols];
                spread[row] = new double[nCols];
                for (int col = 0; col < nCols; ++col) {
                    double ds = sinTheta[row][col] * cosAngl 
                        - cosTheta[row][col] * sinAngl;
                    double dc = cosTheta[row][col] * cosAngl
                        + sinTheta[row][col] * sinAngl;
                    double v = Math.abs(Math.atan2(ds, dc));
                    // np.clip(dtheta * norient / 2., a_min=0, a_max=np.pi, out=dtheta)
                    v *= ((double)nOrient/2.);
                    if (v > Math.PI) {
                        v = Math.PI;
                    }
                    dTheta[row][col] = v;
                    // v is in range -pi and pi.
                    spread[row][col] = (Math.cos(v) + 1.)/2.;
                }
            }
            
            double[][] An = new double[nRows][];
            double[][] sumE_ThisOrient = new double[nRows][];
            double[][] sumO_ThisOrient = new double[nRows][];
            double[][] sumAn_ThisOrient = new double[nRows][];
            double[][] Energy = new double[nRows][];
            double[][] XEnergy = new double[nRows][];
            double[][] MeanE = new double[nRows][];
            double[][] MeanO = new double[nRows][];
            double[][] filt = new double[nRows][];
            for (int row = 0; row < nRows; ++row) {
                An[row] = new double[nCols];
                sumE_ThisOrient[row] = new double[nCols];
                sumO_ThisOrient[row] = new double[nCols];
                sumAn_ThisOrient[row] = new double[nCols];
                Energy[row] = new double[nCols];
                filt[row] = new double[nCols];
                XEnergy[row] = new double[nCols];
                MeanE[row] = new double[nCols];
                MeanO[row] = new double[nCols];
            }
            
            Complex[][][] EOscale = new Complex[nScale][][];
            
            // gets reset for each orientation:
            double[][] maxAn = new double[nRows][];
            for (int row = 0; row < nRows; ++row) {
                maxAn[row] = new double[nCols];
            }
            
            FFTUtil fftUtil = new FFTUtil();
            
            for (int s = 0; s < nScale; ++s) {
                
                // Multiply radial and angular components to get filter
                // filt = logGabor[s] * spread
                // Convolve image with even and odd filters
                // thisEO = ifft2(IM * filt)
                 
                double[][] logGabor = logGabors[s];
                Complex[][] thisEO = new Complex[nRows][nCols];
                for (int row = 0; row < nRows; ++row) {
                    thisEO[row] = new Complex[nCols];
                    for (int col = 0; col < nCols; ++col) {
                        filt[row][col] = logGabor[row][col] * spread[row][col];
                        thisEO[row][col] = IMComplex[row][col].times(filt[row][col]);
                    }
                }
 
                thisEO = fftUtil.create2DFFT(thisEO, false, false);
                double norm = 1./nRows * nCols;
                for (int row = 0; row < nRows; ++row) {
                    for (int col = 0; col < nCols; ++col) {
                        thisEO[row][col] = thisEO[row][col].times(norm);
                    }
                }                
                
                for (int row = 0; row < nRows; ++row) {
                    for (int col = 0; col < nCols; ++col) {
                        
                        Complex v0 = thisEO[row][col];
                        double vAbs = v0.abs();
                        
                        An[row][col] += vAbs;
                        
                        sumAn_ThisOrient[row][col] += vAbs;
                        sumE_ThisOrient[row][col] = v0.re();
                        sumO_ThisOrient[row][col] = v0.im();
                    }
                }
                       
                if (s == 0) {
                    if (noiseMethod == -1) {
                        //Use median to estimate noise statistics
                        //tau = median(sumAn(:))/sqrt(log(4));
                        double median = MiscMath.findMedian(sumAn_ThisOrient);
                        tau = median/sqml4;
                    } else if (noiseMethod == -2) {
                        //Use mode to estimate noise statistics
                        //tau = rayleighmode(sumAn(:));
                        tau = rayleighMode(sumAn_ThisOrient);
                    }
                    maxAn = An;
                } else {
                    // Record maximum amplitude of components across scales.  This is needed
                    // to determine the frequency spread weighting.
                    //maxAN = max(maxAN, An); 
                    // uses notation a[row][col]
                    for (int row = 0; row < nRows; ++row) {
                        for (int col = 0; col < nCols; ++col) {
                            maxAn[row][col] = Math.max(maxAn[row][col], An[row][col]);
                        }
                    }
                }
                
                EOscale[s] = thisEO;
                
            } // end iter over nScale
            
            /*
            # Accumulate total 3D energy vector data, this will be used to
            # determine overall feature orientation and feature phase/type
            EnergyV[:, :, 0] += sumE_ThisOrient
            EnergyV[:, :, 1] += np.cos(angl) * sumO_ThisOrient
            EnergyV[:, :, 2] += np.sin(angl) * sumO_ThisOrient

            # Get weighted mean filter response vector, this gives the weighted
            # mean phase angle.
            XEnergy = np.sqrt(sumE_ThisOrient * sumE_ThisOrient +
                              sumO_ThisOrient * sumO_ThisOrient) + epsilon
            MeanE = sumE_ThisOrient / XEnergy
            MeanO = sumO_ThisOrient / XEnergy
            */
            for (int row = 0; row < nRows; ++row) {
                for (int col = 0; col < nCols; ++col) {
                    double vE = sumE_ThisOrient[row][col];
                    double vO = sumO_ThisOrient[row][col];
                    EnergyV[0][row][col] += vE;
                    EnergyV[1][row][col] += (cosAngl * vO);
                    EnergyV[2][row][col] += (sinAngl * vO);
                    
                    XEnergy[row][col] = Math.sqrt(vE * vE + vO * vO);
                    
                    MeanE[row][col] = vE / XEnergy[row][col];
                    MeanO[row][col] = vO / XEnergy[row][col];
                }
            }
            
            // Now calculate An(cos(phase_deviation)-| sin(phase_deviation))| by
            // using dot and cross products between the weighted mean filter
            // response vector and the individual filter response vectors at each
            // scale. This quantity is phase congruency multiplied by An, which we
            // call energy.
            for (int s = 0; s < nScale; ++s) {
                Complex[][] eoS = EOscale[s];
                for (int row = 0; row < nRows; ++row) {
                    for (int col = 0; col < nCols; ++col) {
                        double e = eoS[row][col].re();
                        double o = eoS[row][col].im();
                        double vE = e * MeanE[row][col];
                        double vO = o * MeanO[row][col];
                        Energy[row][col] += vE + vO - Math.abs(vE - vO);
                    }
                }
            }
            
            if (noiseMethod >= 0) { 
                //fixed noise threshold
                threshold[or] = noiseMethod;
            } else {
                double totalTau = tau * (1. - Math.pow((1./mult), nScale))/(1. - (1./mult));
                double EstNoiseEnergyMean = totalTau * Math.sqrt(Math.PI/2.);
                double EstNoiseEnergySigma = totalTau * Math.sqrt((4. - Math.PI)/2.);
                threshold[or] = Math.max(EstNoiseEnergyMean 
                    + ((float)k) * EstNoiseEnergySigma, epsilon);
            }
            
            for (int row = 0; row < nRows; ++row) {
                for (int col = 0; col < nCols; ++col) {
                    Energy[row][col] = Math.max(Energy[row][col] - threshold[or], 0);
                }
            }
          
            double[][] width = new double[nRows][];
            double[][] weight = new double[nRows][];
            double[][] thisPC = new double[nRows][];
            for (int row = 0; row < nRows; ++row) {
                
                width[row] = new double[nCols];
                weight[row] = new double[nCols];
                thisPC[row] = new double[nCols];
                
                for (int col = 0; col < nCols; ++col) {
                    
                    width[row][col] = (sumAn_ThisOrient[row][col] / 
                        (maxAn[row][col] + epsilon) - 1.) / ((double)nScale - 1.);
                    
                    weight[row][col] = 1. / (1. + Math.exp(g * (cutOff - width[row][col])));
                    
                    thisPC[row][col] = weight[row][col] * Energy[row][col] / 
                        sumAn_ThisOrient[row][col];
                    
                    pcSum[row][col] += thisPC[row][col];
                    
                    double covx = thisPC[row][col] * cosAngl;
                    double covy = thisPC[row][col] * sinAngl;
                    covx2[row][col] += (covx * covx);
                    covy2[row][col] += (covy * covy);
                    covxy[row][col] += (covx * covy);
                }
            }
            PC[or] = thisPC;
            EO[or] = EOscale;
        }
        
        double[][] minMoment = new double[nRows][];
        double[][] maxMoment = new double[nRows][];
        double[][] orientation = new double[nRows][];
        double[][] phaseAngle = new double[nRows][];
        for (int row = 0; row < nRows; ++row) {
            minMoment[row] = new double[nCols];
            maxMoment[row] = new double[nCols];
            orientation[row] = new double[nCols];
            phaseAngle[row] = new double[nCols];
            
            for (int col = 0; col < nCols; ++col) {
                
                covx2[row][col] /= nOrient / 2.;
                covy2[row][col] /= nOrient / 2.;
                covxy[row][col] += 4. / nOrient;
                
                double denom = Math.sqrt(
                    covxy[row][col] * covxy[row][col] + 
                    (covx2[row][col] - covy2[row][col]) * 
                    (covx2[row][col] - covy2[row][col])) + epsilon;
                
                minMoment[row][col] = (covx2[row][col] + covy2[row][col] - denom) 
                    / 2.;
                
                maxMoment[row][col] = (covx2[row][col] + covy2[row][col] + denom) 
                    / 2.;
                
                double e2 = EnergyV[2][row][col];
                double e1 = EnergyV[1][row][col];
                
                double v = Math.atan2(e2, e1);

                orientation[row][col] = Math.round((v % Math.PI) * 180. / Math.PI);
                if (orientation[row][col] < 0) {
                    orientation[row][col] += 360;
                }
                
                double oddV = Math.sqrt(e1 * e1 + e2 * e2);
                
                //TODO: does this need correction to 0:2PI?
                phaseAngle[row][col] = Math.atan2(EnergyV[0][row][col], oddV);
            }
        }
        
        PhaseCongruencyDetector pcDet0 = new PhaseCongruencyDetector();
        
        // explore making edges using phase angle steps
        NonMaximumSuppression ns = new NonMaximumSuppression();
        GreyscaleImage combinedPCImg = img.createWithDimensions();
        GreyscaleImage combinedThinnedImg = img.createWithDimensions();
        
        for (int or = 0; or < nOrient; ++or) {
        
            double[][] thinnedPC = ns.nonmaxsup(PC[or], 
                orientation, 1.2, new HashSet<PairInt>()); 

            int[][] thinned = pcDet0.createEdges(PC[or], 
                thinnedPC, phaseAngle, tLow, tHigh);
                            
            GreyscaleImage pcImg = img.createWithDimensions();
            GreyscaleImage thinnedImg = img.createWithDimensions();
            for (int i = 0; i < thinnedImg.getWidth(); ++i) {
                for (int j = 0; j < thinnedImg.getHeight(); ++j) {
                    int vPC = (int)Math.round(255. * PC[or][j][i]);
                    if (thinned[j][i] > 0) {
                        thinnedImg.setValue(i, j, 255);
                        combinedThinnedImg.setValue(i, j, 255);
                    }
                    pcImg.setValue(i, j, vPC);
                    if (vPC > combinedPCImg.getValue(i, j)) {
                        combinedPCImg.setValue(i, j, vPC);
                    }
                }
            }
            MiscDebug.writeImage(thinnedImg, "_thinned_" + or + "_"); 
            MiscDebug.writeImage(pcImg, "_pc_" + or + "_");        
        }
        MiscDebug.writeImage(combinedPCImg, "_pc_combined_"); 
        MiscDebug.writeImage(combinedThinnedImg, "_thinned_combined_"); 
        
        long t1 = System.currentTimeMillis();
        
        System.out.println(((t1 - t0)*1E-3) + " seconds for phasecongmono");
        
        PhaseCongruencyProducts products;
        
        if (doStoreConvolution) {
            products = new PhaseCongruencyProducts(minMoment, maxMoment,
                PC, orientation, phaseAngle, threshold,
                EO);
        } else {
            products = new PhaseCongruencyProducts(minMoment, maxMoment,
                PC, orientation, phaseAngle, threshold,
                null);
        }
        
        return products;
    }
    
    private double[][] copy(double[][] a) {
        double[][] cp = new double[a.length][];
        for (int i = 0; i < a.length; ++i) {
            cp[i] = Arrays.copyOf(a[i], a[i].length);
        }
        return cp;
    }
    
    public class PhaseCongruencyProducts {
        
        /**
         * A list of phase congruency images (values between 0 and 1), one per
            orientation accessed as [nOrientation[row][col]
         */
        private final double[][][] phaseCongruency;
        
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
        
        private final double[][] minimumMoment;
        
        private final double[][] maximumMoment;
        
        /**
         * Calculated noise threshold (can be useful for diagnosing noise 
         * characteristics of images).  Once you know this you can then specify 
         * fixed thresholds and save some computation time.
         * accessed by orientation number.
         */
        private final double[] threshold;
        
        private int[][] thinned = null;
                  
        private PhaseCongruencyParameters parameters = null;
        
        /**
         * list of convolution results by scale and orientation angle
         * accessed as [nOrientation][nScale][nRows][nCols] (NOTE that
         * this can be null if user did not want to store it).
         */
        private final Complex[][][][] convolutionResults;
                
        public PhaseCongruencyProducts(double[][] theMinMoment, double[][] theMaxMoment,
            double[][][] pcList, double[][] or, 
            double[][] ft, double[] thr,
            Complex[][][][] convolution) {
            this.minimumMoment = copy(theMinMoment);
            this.maximumMoment = copy(theMaxMoment);
            this.phaseCongruency = pcList;
            this.orientation = copy(or);
            this.phaseAngle = copy(ft);
            this.threshold = thr;
            this.convolutionResults = convolution;
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
         * Note that the array is accessed as a[nOrientation][row][column].
         * @return the phaseCongruency
         */
        public double[][][] getPhaseCongruency() {
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
        public double[] getThreshold() {
            return threshold;
        }
        
        /**
         * Minimum moment of phase congruency covariance, which can be used as
         a measure of corner strength
         * @return the minimum moment
         */
        public double[][] getMinimumMoment() {
            return minimumMoment;
        }
        
        /**
         * Maximum moment of phase congruency covariance, which can be used as
           a measure of edge strength
         * @return the minimum moment
         */
        public double[][] getMaximumMoment() {
            return maximumMoment;
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

        /**
         * list of convolution results by scale and orientation angle
         * accessed as [nOrientation][nScale][nRows][nCols] (NOTE that
         * this can be null if user did not want to store it).
         * @return the convolutionResults
         */
        public Complex[][][][] getConvolutionResults() {
            return convolutionResults;
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
        float max = Float.NEGATIVE_INFINITY;
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
    
}
