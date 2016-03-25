package algorithms.imageProcessing;

import algorithms.misc.Complex;
import java.util.Arrays;

/**
 * adapted from 
 * http://www.peterkovesi.com/matlabfns/FrequencyFilt/perfft2.m
 * 
 * The following comments are directly from the kovesi perfft2.m
 * which was written for matlab.
 * 
 * <pre>
 * PERFFT2  2D Fourier transform of Moisan's periodic image component
 
  Usage: [P, S, p, s] = perfft2(im)
 
  Argument:  im - Image to be transformed
  Returns:    P - 2D fft of periodic image component
              S - 2D fft of smooth component
              p - Periodic component (spatial domain)
              s - Smooth component (spatial domain)
 
  Moisan's "Periodic plus Smooth Image Decomposition" decomposes an image 
  into two components
         im = p + s
  where s is the 'smooth' component with mean 0 and p is the 'periodic'
  component which has no sharp discontinuities when one moves cyclically across
  the image boundaries.
  
  This wonderful decomposition is very useful when one wants to obtain an FFT of
  an image with minimal artifacts introduced from the boundary discontinuities.
  The image p gathers most of the image information but avoids periodization
  artifacts.
 
  The typical use of this function is to obtain a 'periodic only' fft of an
  image 
    >>  P = perfft2(im);
 
  Displaying the amplitude spectrum of P will yield a clean spectrum without the
  typical vertical-horizontal 'cross' arising from the image boundaries that you
  would normally see.
 
  Note if you are using the function to perform filtering in the frequency
  domain you may want to retain s (the smooth component in the spatial domain)
  and add it back to the filtered result at the end.  
 
  The computational cost of obtaining the 'periodic only' FFT involves taking an
  additional FFT.
 
   Reference: 
  This code is adapted from Lionel Moisan's Scilab function 'perdecomp.sci' 
  "Periodic plus Smooth Image Decomposition" 07/2012 available at
 
    http://www.mi.parisdescartes.fr/~moisan/p+s
 
  Paper:
  L. Moisan, "Periodic plus Smooth Image Decomposition", Journal of
  Mathematical Imaging and Vision, vol 39:2, pp. 161-179, 2011.
 </pre>
 * @author nichole
 */
public class PeriodicFFT {
    
    /**
     * calculates the 2D fft of periodic image component and returns that and
     * the 2D fft of smooth component.
     * 
     * perfft2 is ported from the python phasepack tools code which has the 
     * following copyright:
     * <pre>
     * Authorship
        ----------
        These functions were originally written for MATLAB by Peter Kovesi, and
        were ported to Python by Alistair Muldal. The original MATLAB code, as well as
        further explanatory information and references are available from `Peter
        Kovesi's website
        <http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/index.html#phasecong>`_.

        MIT License
        -----------
        Permission is hereby  granted, free of charge, to any  person obtaining a copy
        of this software and associated  documentation files (the "Software"), to deal
        in the Software without restriction, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.

        The software is provided "as is", without warranty of any kind.
     * </pre>
     * Moisan's Periodic plus Smooth Image Decomposition. The image is
     * decomposed into two parts:
     *
     * im = s + p
     *
     * where 's' is the 'smooth' component with mean 0, and 'p' is the
     * 'periodic' component which has no sharp discontinuities when one moves
     * cyclically across the image boundaries.
     *
     * useage: P, S = perfft2(im)
     * 
     * @param img 
     * @param computeSpatial if set to true, returns new Complex[][][]{S, P, s, p}
     * @return new Complex[][][]{S, P}
     *    where S is the FFT of the smooth component
     *          P is the FFT of the periodic component
     *          s smooth component in the spatial domain
     *          p periodic component in the spatial domain
     */
    public Complex[][][] perfft2(GreyscaleImage img, boolean computeSpatial) {
        
        int nRows = img.getHeight();
        int nCols = img.getWidth();
        
        /*
        Compute the boundary image which is equal to the image discontinuity
        values across the boundaries at the edges and is 0 elsewhere
        s = np.zeros_like(im)
        s[0, :] = im[0, :] - im[-1, :]
        s[-1, :] = -s[0, :]
        s[:, 0] = s[:, 0] + im[:, 0] - im[:, -1]
        s[:, -1] = s[:, -1] - im[:, 0] + im[:, -1]
        */
        double[][] s = new double[nCols][];
        for (int i = 0; i < s.length; ++i) {
            s[i] = new double[nRows];
        }
        
        for (int col = 0; col < nCols; ++col) {
            //s[row=0, col=:] = im[row=0, col=:] - im[row=-1, col=:]
            int v = img.getValue(col, 0) - img.getValue(col, nRows - 1);
            s[col][0] = v;
        }
        for (int col = 0; col < nCols; ++col) {
            //s[row=-1, col=:] = -s[row=0, col=:]
            double v = s[col][0];
            s[col][nRows - 1] = -v;
        }
        for (int row = 0; row < nRows; ++row) {
            //s[row=:, col=0] = s[row=:, col=0] + im[row=:, col=0] - im[rod=:, col=-1]
            double v = s[0][row] + img.getValue(0, row) - img.getValue(nCols - 1, row);
            s[0][row] = v;
        }
        for (int row = 0; row < nRows; ++row) {
            //s[row=:, col=-1] = s[row=:, col=-1] - im[row=:, col=0] + im[row=:, col=-1]
            double v = s[nCols - 1][row] - img.getValue(0, row) + img.getValue(nCols - 1, row);
            s[nCols - 1][row] = v;
        }
                
        /*
        Generate grid upon which to compute the filter for the boundary image in
        the frequency domain.  Note that cos() is cyclic hence the grid values can
        range from 0 .. 2*pi rather than 0 .. pi and then pi .. 0
        
        x, y = (2 * np.pi * np.arange(0, v) / float(v) for v in (cols, rows))
        cx, cy = np.meshgrid(x, y)
        */
        double[][] cx = new double[nCols][];
        for (int i = 0; i < nCols; ++i) {
            cx[i] = new double[nRows];
            double v = (2. * Math.PI * i)/(double)nCols;
            Arrays.fill(cx[i], v);
        }
        double[][] cy = new double[nCols][];
        for (int i = 0; i < nCols; ++i) {
            cy[i] = new double[nRows];
        }
        for (int row = 0; row < nRows; ++row) {
            double v = (2. * Math.PI * row)/(double)nRows;
            for (int col = 0; col < nCols; ++col) {
                cy[col][row] = v;
            }
        }
        
        /*
        denom = (2. * (2. - np.cos(cx) - np.cos(cy)))
        denom[0, 0] = 1.     # avoid / 0
        
        Generate FFT of smooth component
        S = fft2(s)./denom
        */
        
        // the fft algorithm using powers of 2 sizes is faster, so pad the data
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Complex[][] capS = imageProcessor.create2DFFT(s, true);
        for (int col = 0; col < nCols; ++col) {
            for (int row = 0; row < nRows; ++row) {
                double denom;
                if (col == 0 && row == 0) {
                    denom = 1;
                } else {
                    denom = (2. * (2. - Math.cos(cx[col][row]) - Math.cos(cy[col][row])));
                }
                capS[col][row] = capS[col][row].divided(new Complex(denom, 0));
            }
        }
        
        /*
        % The (1,1) element of the filter will be 0 so S(1,1) may be Inf or NaN
        S(1,1) = 0;          % Enforce 0 mean 
        */
        capS[0][0] = new Complex(0, 0);
        
        /*
        P = fft2(im) - S;    % FFT of periodic component
        */
        // initialize matrix of complex numbers as real numbers from image (imaginary are 0's)
        Complex[][] capP = imageProcessor.create2DFFT(imageProcessor.convertImage(img), true);
             
        for (int i = 0; i < nCols; ++i) {
            for (int j = 0; j < nRows; ++j) {
                Complex v0 = capP[i][j];
                Complex s0 = capS[i][j];
                capP[i][j] = v0.minus(s0);
            }
        }
        
            
        //DEBUG
        /*{
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < capP.length; ++i) {
            for (int j = 0; j < capP[0].length; ++j) {
                Complex v = capP[i][j];
                sb.append(String.format("S (%d,%d) %f + i %f\n", i, j, (float)v.re(),
                    (float)v.im()));
            }
        }
        System.out.println(sb.toString());
        }*/
        
        
        if (computeSpatial) {
            //s = real(ifft2(S)); 
            //p = im - s;
            
            Complex[][] lowerCaseS = imageProcessor.ifftShift(capS);
            for (int i = 0; i < lowerCaseS.length; ++i) {
                for (int j = 0; j < lowerCaseS[0].length; ++j) {
                    lowerCaseS[i][j] = new Complex(lowerCaseS[i][j].re(), 0);
                }
            }
            
            Complex[][] lowerCaseP = new Complex[lowerCaseS.length][];
            for (int i = 0; i < lowerCaseS.length; ++i) {
                lowerCaseP[i] = new Complex[lowerCaseS[i].length];
                for (int j = 0; j < lowerCaseS[i].length; ++j) {
                    lowerCaseP[i][j] = new Complex(
                        img.getValue(i, j) - lowerCaseS[i][j].re(), 0);
                }
            }
            
            return new Complex[][][]{capS, capP, lowerCaseS, lowerCaseP};
        }
        /*
        if nargout > 2       % Generate spatial domain results 
            s = real(ifft2(S)); 
            p = im - s;         
        end
        */
        
        return new Complex[][][]{capS, capP};
    }
    
}
