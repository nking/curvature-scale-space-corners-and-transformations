package algorithms.imageProcessing;

import algorithms.misc.Complex;
import java.util.Arrays;

/**
 * adapted from
 * http://pydoc.net/Python/phasepack/1.4/phasepack.phasecong/
 * which has copyright:
 * # MIT License:

# Permission is hereby  granted, free of charge, to any  person obtaining a
# copy of this software and associated  documentation files (the "Software"),
# to deal in the Software without restriction, subject to the following
# conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# The software is provided "as is", without warranty of any kind.

# Original MATLAB version by Peter Kovesi
# <http://www.csse.uwa.edu.au/~pk/research/matlabfns/PhaseCongruency/phasecong3.m>

#Python translation by Alistair Muldal
# <alistair muldal@pharm ox ac uk>
* 
* Generates grid for constructing frequency domain filters using radius and theta
%
% Usage:  [radius, cosTheta, sinTheta] = polarFiltergrid(rows, cols)
%
% Arguments:  rows, cols - Size of image/filter
%
% Returns:        
 * 
 */
public class PolarFilterGrid {
    
    /**
     * usage: filtergrid(nRows, nCols)
     * 
     * @param nRows
     * @param nCols 
     * @returns [radius, u1, u2] where 
          radius - Grid of size [nRows nCols] containing normalised
                   radius values from 0 to 0.5.  Grid is quadrant
                   shifted so that 0 frequency is at radius(1,1)
          cosTheta, sinTheta - Grids containing normalised frequency values
                   ranging from -0.5 to 0.5 in x and y directions
                   respectively. u1 and u2 are quadrant shifted.
          NOTE: the returned results use notation a[row][col]
    */
    public FilterGridProducts filtergrid(int nRows, int nCols) {
                
        /*
        Set up X and Y spatial frequency matrices, u1 and u2, with ranges
        normalised to +/- 0.5 The following code adjusts things appropriately for
        odd and even values of rows and columns so that the 0 frequency point is
        placed appropriately.
        
        if mod(cols,2)
            u1range = [-(cols-1)/2:(cols-1)/2]/(cols-1);
        else
            u1range = [-cols/2:(cols/2-1)]/cols; 
        end

        if mod(rows,2)
            u2range = [-(rows-1)/2:(rows-1)/2]/(rows-1);
        else
            u2range = [-rows/2:(rows/2-1)]/rows; 
        end

        [u1,u2] = meshgrid(u1range, u2range);
        */
        
        double[] u1Range = new double[nCols];
        double div;
        if ((nCols & 1) == 0) {
            //u1range = [-(cols-1)/2:(cols-1)/2]/(cols-1);
            //if nCols=3, this becomes [-1, 0, 1] --> [-0.5, 0, 0.5]
            u1Range[0] = -(nCols-1)/2.;
            div = nCols - 1;            
        } else {
            //u1range = [-cols/2:(cols/2-1)]/cols; 
            //if nCols=4, this becomes [-2, -1, 0, 1] --> [-0.5, -0.25, 0, 0.25]
            u1Range[0] = -nCols/2;
            div = nCols;
        }
        for (int i = 1; i < u1Range.length; ++i) {
            u1Range[i] = u1Range[i - 1] + 1.;
        }
        for (int i = 0; i < u1Range.length; ++i) {
            u1Range[i] /= div;
        }
        
        double[] u2Range = new double[nRows];
        if ((nRows & 1) == 0) {
            u2Range[0] = -(nRows-1)/2.;
            div = nRows - 1;            
        } else {
            u2Range[0] = -nRows/2.;
            div = nRows;
        }
        for (int i = 1; i < u2Range.length; ++i) {
            u2Range[i] = u2Range[i - 1] + 1.;
        }
        for (int i = 0; i < u2Range.length; ++i) {
            u2Range[i] /= div;
        }

        // nRows X nCols
        //[u1,u2] = meshgrid(u1range, u2range);
        double[][] u1 = new double[nRows][];
        double[][] u2 = new double[nRows][];
        for (int i = 0; i < nRows; ++i) {
            u1[i] = new double[nCols];
            u2[i] = new double[nCols];            
            System.arraycopy(u1Range, 0, u1[i], 0, nCols);
        }
        for (int i = 0; i < nRows; ++i) {
            double v = u2Range[i];
            for (int j = 0; j < nCols; ++j) {
                u2[i][j] = v;
            }
        }
        
        int len0 = u1.length;
        int len1 = u1[0].length;
        double[][] theta0 = new double[len0][];
        for (int i = 0; i < len0; ++i) {
            theta0[i] = new double[len1];
            for (int j = 0; j < len1; ++j) {
                double v = u1[i][j] * u1[i][j] + u2[i][j] * u2[i][j];
                theta0[i][j] = Math.atan2(-u2[i][j], u1[i][j]);
            }
        }
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        u1 = imageProcessor.ifftShift(u1);
        u2 = imageProcessor.ifftShift(u2);
       
        theta0 = imageProcessor.ifftShift(theta0);
        
        double[][] radius = new double[len0][];
        double[][] sinTheta = new double[len0][];
        double[][] cosTheta = new double[len0][];
        for (int i = 0; i < len0; ++i) {
            radius[i] = new double[len1];
            sinTheta[i] = new double[len1];
            cosTheta[i] = new double[len1];
            for (int j = 0; j < len1; ++j) {
                double v = u1[i][j] * u1[i][j] + u2[i][j] * u2[i][j];
                radius[i][j] = Math.sqrt(v);
                sinTheta[i][j] = Math.sin(theta0[i][j]);
                cosTheta[i][j] = Math.cos(theta0[i][j]);
            }
        }

        radius[0][0] = 1;
        
        FilterGridProducts products = new FilterGridProducts(radius, 
            cosTheta, sinTheta);
        
        return products;
    }
 
    public class FilterGridProducts {
        
        private final double[][] radius;
        private final double[][] cosTheta;
        private final double[][] sinTheta;
        
        public FilterGridProducts(double[][] theRadius, double[][] theCosTheta,
            double[][] theSinTheta) {
            radius = theRadius;
            cosTheta = theCosTheta;
            sinTheta = theSinTheta;
        }

        /**
         * @return the radius
         */
        public double[][] getRadius() {
            return radius;
        }

        /**
         * @return the cosTheta
         */
        public double[][] getCosTheta() {
            return cosTheta;
        }

        /**
         * @return the sinTheta
         */
        public double[][] getSinTheta() {
            return sinTheta;
        }
        
    }
}
