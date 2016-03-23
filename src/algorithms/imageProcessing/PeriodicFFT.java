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
     * @param img
     * @return 
     */
    public Complex[][][] perfft2(GreyscaleImage img) {
        
        int nRows = img.getHeight();
        int nCols = img.getWidth();
        
        /*
        Compute the boundary image which is equal to the image discontinuity
        values across the boundaries at the edges and is 0 elsewhere
                
         s = zeros(size(im));
         s(1,:)   = im(1,:) - im(end,:);
         s(end,:) = -s(1,:);
         s(:,1)   = s(:,1)   + im(:,1) - im(:,end);
         s(:,end) = s(:,end) - im(:,1) + im(:,end);
         */
        double[][] s = new double[nCols][];
        for (int i = 0; i < s.length; ++i) {
            s[i] = new double[nRows];
        }
        
        for (int row = 0; row < nRows; ++row) {
            int v = img.getValue(0, row) - img.getValue(nCols - 1, row);
            s[0][row] = v;
        }
        for (int row = 0; row < nRows; ++row) {
            double v = s[0][row];
            s[nCols - 1][row] = -v;
        }
        for (int col = 0; col < nCols; ++col) {
            double v = s[col][0] + img.getValue(col, 0) -
                img.getValue(col, nRows - 1);
            s[col][0] = v;
        }
        for (int col = 0; col < nCols; ++col) {
            double v = s[col][nRows - 1] - img.getValue(col, 0) +
                img.getValue(col, nRows - 1);
            s[col][nRows - 1] = v;
        }
                
        /*
        Generate grid upon which to compute the filter for the boundary image in
        the frequency domain.  Note that cos() is cyclic hence the grid values can
        range from 0 .. 2*pi rather than 0 .. pi and then pi .. 0
        [cx, cy] = meshgrid(
            2*pi*[0:cols-1]/cols, 
            2*pi*[0:rows-1]/rows); 
        
        matlab vector operation is applies to all the members:
            [0, 1, 2, 3, 4, ... nCols - 1]
        */
        double[][] cx = new double[nCols][];
        for (int i = 0; i < nCols; ++i) {
            cx[i] = new double[nRows];
            double v = (2. * Math.PI * i)/nCols;
            Arrays.fill(cx[i], v);
        }
        double[][] cy = new double[nCols][];
        for (int i = 0; i < nCols; ++i) {
            cy[i] = new double[nRows];
        }
        for (int row = 0; row < nRows; ++row) {
            double v = (2. * Math.PI * row)/nRows;
            for (int col = 0; col < nCols; ++col) {
                cy[col][row] = v;
            }
        }
        
        /*
        Generate FFT of smooth component
        S = fft2(s)./(2*(2 - cos(cx) - cos(cy)));
        */
        ImageProcessor imageProcessor = new ImageProcessor();
        s = imageProcessor.padUpToPowerOfTwo2(s);
        // initialize matrix of complex numbers as real numbers from image (imaginary are 0's)
        Complex[][] cc = imageProcessor.convertImage(s);
        Complex[][] capS = imageProcessor.apply2DFFT(cc, true);
        for (int i = 0; i < capS.length; ++i) {
            if (i > (nCols - 1)) {
                break;
            }
            for (int j = 0; j < capS[0].length; ++j) {
                if (j > (nRows - 1)) {
                    break;
                }
                //(2*(2 - cos(cx) - cos(cy)))
                double f = 2. * (2. - Math.cos(cx[i][j]) - Math.cos(cy[i][j]));
                double v = capS[i][j].re() / f;
                capS[i][j] = new Complex(v, 0);
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
        GreyscaleImage p = imageProcessor.padUpToPowerOfTwo2(img);
        // initialize matrix of complex numbers as real numbers from image (imaginary are 0's)
        cc = imageProcessor.convertImage(p);
        Complex[][] capP = imageProcessor.apply2DFFT(cc, true);
        for (int i = 0; i < capS.length; ++i) {
            for (int j = 0; j < capS[0].length; ++j) {
                Complex v0 = capP[i][j];
                Complex s0 = capS[i][j];
                capP[i][j] = v0.minus(s0);
            }
        }
        
        /*
        if nargout > 2       % Generate spatial domain results 
            s = real(ifft2(S)); 
            p = im - s;         
        end
        */
        
        // trim datasets back down
        Complex[][] capP2 = new Complex[nCols][];
        for (int col = 0; col < nCols; ++col) {
            capP2[col] = new Complex[nRows];
        }
        for (int col = 0; col < nCols; ++col) {
            for (int row = 0; row < nRows; ++row) {
                capP2[col][row] = capP[col][row];
            }
        }
        
        Complex[][] capS2 = new Complex[nCols][];
        for (int col = 0; col < nCols; ++col) {
            capS2[col] = new Complex[nRows];
        }
        for (int col = 0; col < nCols; ++col) {
            for (int row = 0; row < nRows; ++row) {
                capS2[col][row] = capS[col][row];
            }
        }
        
        //[P, S, p, s]
        return new Complex[][][]{capP2, capS2};
    }
    
}
