package algorithms.imageProcessing;

import algorithms.misc.Complex;
import java.util.Arrays;

/**
 * adapted from 
 * http://www.peterkovesi.com/matlabfns/FrequencyFilt/lowpassfilter.m
 * which has copyright:
 * 
 * Copyright (c) 1999 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% October 1999
% August  2005 - Fixed up frequency ranges for odd and even sized filters
%                (previous code was a bit approximate)
* 
 * @author nichole
 */
public class LowPassFilter {
 
    /**
     * Constructs a low-pass butterworth filter
     *                 1
     *   f = -------------------- 
     *         1.0 + (w/cutoff)^(2n)
     * 
     * The frequency origin of the returned filter is at the corners.
     * @param nRows number of rows to use in filter
     * @param nCols number of columns to use in filter
     * @param cutoff cutoff frequency of the filter, 0 to 0.5
     * @param n is the order of the filter, the higher n is the sharper the 
     * transition is. (n must be an integer >= 1).  Note that n is doubled so 
     * that it is always an even integer.
     * @return 
     */
    public double[][] lowpassfilter(int nRows, int nCols, float cutoff, int n) {
        
        if ((cutoff < 0) || (cutoff > 0.5)) {
            throw new IllegalArgumentException(
                "cutoff frequency must be between 0 and 0.5");
        }
        
        if (n < 1) {
            throw new IllegalArgumentException( "n must be >= 1");
        }
        
        /*
        Set up X and Y matrices with ranges normalised to +/- 0.5
        The following code adjusts things appropriately for odd and even values
        of rows and columns.
        if mod(cols,2)
               xrange = [-(cols-1)/2:(cols-1)/2]/(cols-1);
        else
               xrange = [-cols/2:(cols/2-1)]/cols;     
        end

         if mod(rows,2)
               yrange = [-(rows-1)/2:(rows-1)/2]/(rows-1);
         else
               yrange = [-rows/2:(rows/2-1)]/rows;     
         end
    
         [x,y] = meshgrid(xrange, yrange);
         radius = sqrt(x.^2 + y.^2);        % A matrix with every pixel = radius relative to centre.
         f = ifftshift( 1.0 ./ (1.0 + (radius ./ cutoff).^(2*n)) );   % The filter
        */
        double[] xRange = new double[nCols];
        if ((nCols & 1) == 1) {
            //u1range = [-(cols-1)/2:(cols-1)/2]/(cols-1);
            //if nCols=3, this becomes [-1, 0, 1] --> [-0.5, 0, 0.5]
            xRange[0] = -(nCols-1)/2.;
            for (int i = 1; i < nCols; ++i) {
                xRange[i] = xRange[0] + i;
            }
            for (int i = 0; i < nCols; ++i) {
                xRange[i] /= (nCols - 1.);
            }            
        } else {
            //u1range = [-cols/2:(cols/2-1)]/cols; 
            //if nCols=4, this becomes [-2, -1, 0, 1] --> [-0.5, -0.25, 0, 0.25]
            xRange[0] = -nCols/2.;
            for (int i = 1; i < nCols; ++i) {
                xRange[i] = xRange[0] + i;
            }
            for (int i = 0; i < nCols; ++i) {
                xRange[i] /= (double)(nCols);
            }            
        }
        
        double[] yRange = new double[nRows];
        if ((nRows & 1) == 1) {
            //u2range = [-(rows-1)/2:(rows-1)/2]/(rows-1);
            //if nRows=3, this becomes [-1, 0, 1] --> [-0.5, 0, 0.5]
            yRange[0] = -(nRows-1)/2.;
            for (int i = 1; i < nRows; ++i) {
                yRange[i] = yRange[0] + i;
            }
            for (int i = 0; i < nRows; ++i) {
                yRange[i] /= (nRows - 1.);
            }            
        } else {
            //u2range = [-rows/2:(rows/2-1)]/rows; 
            //if nRows=4, this becomes [-2, -1, 0, 1] --> [-0.5, -0.25, 0, 0.25]
            yRange[0] = -nRows/2.;
            for (int i = 1; i < nRows; ++i) {
                yRange[i] = yRange[0] + i;
            }
            for (int i = 0; i < nRows; ++i) {
                yRange[i] /= (double)(nRows);
            }            
        }
        
        // nRows X nCols
        //[x,y] = meshgrid(xrange, yrange);
        double[][] x = new double[nRows][];
        double[][] y = new double[nRows][];
        for (int i = 0; i < nRows; ++i) {
            x[i] = new double[nCols];
            y[i] = new double[nCols];            
            System.arraycopy(xRange, 0, x[i], 0, nCols);
        }
        for (int i = 0; i < nRows; ++i) {
            double v = yRange[i];
            for (int j = 0; j < nCols; ++j) {
                y[i][j] = v;
            }
        }
        
        //% A matrix with every pixel = radius relative to centre.
        //radius = sqrt(x.^2 + y.^2);
        double[][] radius = new double[x.length][];
        for (int i = 0; i < radius.length; ++i) {
            radius[i] = new double[x[i].length];
            for (int j = 0; j < radius[i].length; ++j) {
                double x0 = x[i][j];
                double y0 = y[i][j];
                radius[i][j] = Math.sqrt(x0 * x0 + y0 * y0);
            }
        }
        
        //% The filter
        //f = ifftshift( 1.0 ./ (1.0 + (radius ./ cutoff).^(2*n)) );
        for (int i = 0; i < radius.length; ++i) {
            for (int j = 0; j < radius[i].length; ++j) {
                double v = radius[i][j];
                v /= cutoff;
                v = Math.pow(v, 2*n) + 1;
                radius[i][j] = 1./v;
            }
        }
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        double[][] f = imageProcessor.ifftShift(radius);

        return f;
    }
}
