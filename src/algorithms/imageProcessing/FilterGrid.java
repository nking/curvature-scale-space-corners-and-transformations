package algorithms.imageProcessing;

/**
 * adapted from
 * http://www.peterkovesi.com/matlabfns/FrequencyFilt/filtergrid.m
 * which has copyright:
 * Copyright (c) 1996-2013 Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% peter.kovesi at uwa edu au
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%
% May 2013
* 
* Note that http://pydoc.net/Python/phasepack/1.4/phasepack.filtergrid/
* was also useful in finding ifftshift.
* 
* Generates grid for constructing frequency domain filters
%
% Usage:  [radius, u1, u2] = filtergrid(rows, cols)
%         [radius, u1, u2] = filtergrid([rows, cols])
%
% Arguments:  rows, cols - Size of image/filter
%
% Returns:        
 * 
 * @author nichole
 */
public class FilterGrid {
    
    /**
     * usage: filtergrid(nRows, nCols)
     * 
     * @param nRows
     * @param nCols 
     * @returns [radius, u1, u2] where 
          radius - Grid of size [nRows nCols] containing normalised
                   radius values from 0 to 0.5.  Grid is quadrant
                   shifted so that 0 frequency is at radius(1,1)
          u1, u2 - Grids containing normalised frequency values
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
        if ((nCols & 1) == 1) {
            //u1range = [-(cols-1)/2:(cols-1)/2]/(cols-1);
            //if nCols=3, this becomes [-1, 0, 1] --> [-0.5, 0, 0.5]
            u1Range[0] = -(nCols-1)/2.;
            for (int i = 1; i < nCols; ++i) {
                u1Range[i] = u1Range[0] + i;
            }
            for (int i = 0; i < nCols; ++i) {
                u1Range[i] /= (nCols - 1.);
            }            
        } else {
            //u1range = [-cols/2:(cols/2-1)]/cols; 
            //if nCols=4, this becomes [-2, -1, 0, 1] --> [-0.5, -0.25, 0, 0.25]
            u1Range[0] = -nCols/2.;
            for (int i = 1; i < nCols; ++i) {
                u1Range[i] = u1Range[0] + i;
            }
            for (int i = 0; i < nCols; ++i) {
                u1Range[i] /= (double)(nCols);
            }            
        }
        
        double[] u2Range = new double[nRows];
        if ((nRows & 1) == 1) {
            //u2range = [-(rows-1)/2:(rows-1)/2]/(rows-1);
            //if nRows=3, this becomes [-1, 0, 1] --> [-0.5, 0, 0.5]
            u2Range[0] = -(nRows-1)/2.;
            for (int i = 1; i < nRows; ++i) {
                u2Range[i] = u2Range[0] + i;
            }
            for (int i = 0; i < nRows; ++i) {
                u2Range[i] /= (nRows - 1.);
            }            
        } else {
            //u2range = [-rows/2:(rows/2-1)]/rows; 
            //if nRows=4, this becomes [-2, -1, 0, 1] --> [-0.5, -0.25, 0, 0.25]
            u2Range[0] = -nRows/2.;
            for (int i = 1; i < nRows; ++i) {
                u2Range[i] = u2Range[0] + i;
            }
            for (int i = 0; i < nRows; ++i) {
                u2Range[i] /= (double)(nRows);
            }            
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
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        u1 = imageProcessor.ifftShift(u1);
        u2 = imageProcessor.ifftShift(u2);
       
        int len0 = u1.length;
        int len1 = u1[0].length;
        double[][] radius = new double[len0][];
        for (int i = 0; i < len0; ++i) {
            radius[i] = new double[len1];
            for (int j = 0; j < len1; ++j) {
                double v = u1[i][j] * u1[i][j] + u2[i][j] * u2[i][j];
                radius[i][j] = Math.sqrt(v);
            }
        }
        
        FilterGridProducts products = new FilterGridProducts(radius, u1, u2);
        
        return products;
    }
 
    public class FilterGridProducts {
        
        private final double[][] radius;
        private final double[][] u1;
        private final double[][] u2;
        
        public FilterGridProducts(double[][] theRadius, double[][] theU1,
            double[][] theU2) {
            radius = theRadius;
            u1 = theU1;
            u2 = theU2;
        }

        /**
         * @return the radius
         */
        public double[][] getRadius() {
            return radius;
        }

        /**
         * @return the u1
         */
        public double[][] getU1() {
            return u1;
        }

        /**
         * @return the u2
         */
        public double[][] getU2() {
            return u2;
        }
        
    }
}
