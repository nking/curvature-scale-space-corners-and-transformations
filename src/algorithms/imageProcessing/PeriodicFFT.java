package algorithms.imageProcessing;

import algorithms.misc.Complex;
import algorithms.misc.MiscMath;
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
     *    *NOTE that the output data use the notation a[row]col]
     */
    public Complex[][][] perfft2(GreyscaleImage img, boolean computeSpatial) {
                
        // using notation in PhaseCongruencyDetector,
        // that is, arrays are a[row][col]            
        int nCols = img.getWidth();
        int nRows = img.getHeight();
        
        /*
        Compute the boundary image which is equal to the image discontinuity
        values across the boundaries at the edges and is 0 elsewhere
        s = np.zeros_like(im)
        s[0, :] = im[0, :] - im[-1, :]
        s[-1, :] = -s[0, :]
        s[:, 0] = s[:, 0] + im[:, 0] - im[:, -1]
        s[:, -1] = s[:, -1] - im[:, 0] + im[:, -1]
        */
        // using notation a[row][col]
        double[][] s = new double[nRows][];
        for (int row = 0; row < nRows; ++row) {
            s[row] = new double[nCols];
        }
        
        for (int col = 0; col < nCols; ++col) {
            //s[row=0, col=:] = im[row=0, col=:] - im[row=-1, col=:]
            int v = img.getValue(col, 0) - img.getValue(col, nRows - 1);
            s[0][col] = v;
        }
        for (int col = 0; col < nCols; ++col) {
            //s[row=-1, col=:] = -s[row=0, col=:]
            double v = s[0][col];
            s[nRows - 1][col] = -v;
        }
        for (int row = 0; row < nRows; ++row) {
            //s[:, 0] = s[:, 0]  + im[:, 0] - im[:, -1]
            double v = s[row][0] + img.getValue(0, row) - img.getValue(nCols - 1, row);
            s[row][0] = v;
        }
        for (int row = 0; row < nRows; ++row) {
            //s[row=:, col=-1] = s[row=:, col=-1] - im[row=:, col=0] + im[row=:, col=-1]
            double v = s[row][nCols - 1] - img.getValue(0, row) + img.getValue(nCols - 1, row);
            s[row][nCols - 1] = v;
        }
                
        /*
        Generate grid upon which to compute the filter for the boundary image in
        the frequency domain.  Note that cos() is cyclic hence the grid values can
        range from 0 .. 2*pi rather than 0 .. pi and then pi .. 0
        
        x, y = (2 * np.pi * np.arange(0, v) / float(v) for v in (cols, rows))
        cx, cy = np.meshgrid(x, y)
        */
        // using notation a[row][col]
        double[][] cy = new double[nRows][];
        for (int row = 0; row < nRows; ++row) {
            cy[row] = new double[nCols];
            double v = (2. * Math.PI * row)/(double)nRows;
            Arrays.fill(cy[row], v);
        }
        double[][] cx = new double[nRows][];
        for (int row = 0; row < nRows; ++row) {
            cx[row] = new double[nCols];
        }
        for (int col = 0; col < nCols; ++col) {
            double v = (2. * Math.PI * col)/(double)nCols;
            for (int row = 0; row < nRows; ++row) {
                cx[row][col] = v;
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

        // using notation a[row][col]
        
        //double normCofactor = 1./(2. * Math.PI);
        
        Complex[][] capS = imageProcessor.create2DFFT(s, false, true);
        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {
                double denom;
                if (col == 0 && row == 0) {
                    denom = 1;
                } else {
                    denom = 1./(2. * (2. - Math.cos(cx[row][col]) - Math.cos(cy[row][col])));
                }
                capS[row][col] = capS[row][col].times(denom);
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
        // using notation a[row][col]
        Complex[][] capP = imageProcessor.create2DFFTWithSwapMajor(img, false, true);
        assert(capP.length == capS.length);
        assert(capP[0].length == capS[0].length);
        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {
                Complex v0 = capP[row][col];
                Complex s0 = capS[row][col];
                capP[row][col] = v0.minus(s0);
            }
        }
        
        if (computeSpatial) {
            //s = real(ifft2(S)); 
            //p = im - s;
            
            double postNorm = nRows * nCols;
            
            Complex[][] lowerCaseS = imageProcessor.create2DFFT(capS, false, false);
            for (int row = 0; row < nRows; ++row) {
                for (int col = 0; col < nCols; ++col) {
                    lowerCaseS[row][col] 
                        = new Complex(lowerCaseS[row][col].re()/postNorm, 0);
                }
            }

            Complex[][] lowerCaseP = new Complex[lowerCaseS.length][];
            for (int row = 0; row < nRows; ++row) {
                lowerCaseP[row] = new Complex[nCols];
                for (int col = 0; col < nCols; ++col) {
                    lowerCaseP[row][col] = new Complex(
                        img.getValue(col, row) - lowerCaseS[row][col].re(), 0);
                }
            }
            
            return new Complex[][][]{capS, capP, lowerCaseS, lowerCaseP};
        }
        
        return new Complex[][][]{capS, capP};
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
            float[] x2 = new float[nr];
            for (int ii = 0; ii < nr; ++ii) {
                x2[ii] = ii;
            }
            float[] y2 = new float[nr];
            // plot cols 0.25*nCols, 0.5*nCols, and 0.75*nCols
            for (int nf = 1; nf < 4; nf++) {
                int colNumber = (int) (((float) nf) * 0.25f * nc);
                for (int ii = 0; ii < nr; ++ii) {
                    y2[ii] = (float) tmp[ii][colNumber].re();
                }
                float minY = MiscMath.findMin(y2);
                float maxY = MiscMath.findMax(y2);
                plotter.addPlot(-1, nr + 1, minY, maxY, x2, y2, xPolygon,
                    yPolygon, label + " col=" + colNumber + " REAL");
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
            int z = 1;
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

}
