package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;

/**
 * a version of non-maximum suppression compatible with the PhaseConguencyDetector.
 * 
 * adapted from
 * http://www.peterkovesi.com/matlabfns/Spatial/nonmaxsup.m
 * which has copyright:
 * 
 * Usage:
%          [im,location] = nonmaxsup(inimage, orient, radius);
%
% Function for performing non-maxima suppression on an image using an
% orientation image.  It is assumed that the orientation image gives 
% feature normal orientation angles in degrees (0-180).
%
% Input:
%   inimage - Image to be non-maxima suppressed.
% 
%   orient  - Image containing feature normal orientation angles in degrees
%             (0-180), angles positive anti-clockwise.
% 
%   radius  - Distance in pixel units to be looked at on each side of each
%             pixel when determining whether it is a local maxima or not.
%             This value cannot be less than 1.
%             (Suggested value about 1.2 - 1.5).
*             Caveat: found a value of 1.0 was better for the house.gif test.

* Returns:
%   im        - Non maximally suppressed image.
%   location  - Complex valued image holding subpixel locations of edge
%               points. For any pixel the real part holds the subpixel row
%               coordinate of that edge point and the imaginary part holds
%               the column coordinate.  (If a pixel value is 0+0i then it
%               is not an edgepoint.)
%               (Note that if this function is called without 'location'
%               being specified as an output argument is not computed)
%
% Notes:
%
% The suggested radius value is 1.2 - 1.5 for the following reason. If the
% radius parameter is set to 1 there is a chance that a maxima will not be
% identified on a broad peak where adjacent pixels have the same value.  To
% overcome this one typically uses a radius value of 1.2 to 1.5.  However
% under these conditions there will be cases where two adjacent pixels will
% both be marked as maxima.  Accordingly there is a final morphological
% thinning step to correct this.
%
% This function is slow.  It uses bilinear interpolation to estimate
% intensity values at ideal, real-valued pixel locations on each side of
% pixels to determine if they are local maxima.

% Copyright (c) 1996-2013 Peter Kovesi
* Centre for Exploration Targeting
% The University of Western Australia
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
 * 
 * adapted by nichole to check that gaps are not added to previously connected
 * curves.
 */
public class NonMaximumSuppression {
    
    /**
     * 
     * @param img
     * @param orientation
     * @param radius
     * @param useLowerThreshold if true, uses a lower threshold when checking
     * whether to keep points - the lower threshold is useful for example, when
     * thinning the steps from the phase angle image.
     * @return 
     */
    public double[][] nonmaxsup(double[][] img, double[][] orientation, 
        double radius, boolean useLowerThreshold) {
        
        if (img.length != orientation.length || img[0].length != orientation[0].length) {
            throw new IllegalArgumentException("img and orientation must be same size");
        }
        
        if (radius < 1) {
            throw new IllegalArgumentException("radius must be >= 1");
        }

        int nCols = img.length;
        int nRows = img[0].length;
        
        double[][] output = new double[nCols][];
        for (int i = 0; i < nCols; ++i) {
            output[i] = new double[nRows];
        }
        
        int iRadius = (int)Math.ceil(radius);
        
        // for matlab: Orientations start at 0 degrees but arrays start with index 1
        //orient = fix(orient)+1;
        
        double dToR = Math.PI/180.;
        double[] angle = new double[180];
        double[] xOff = new double[180];
        double[] yOff = new double[180];
        double[] hFrac = new double[180];
        double[] vFrac = new double[180];
        for (int i = 0; i < 180; ++i) {
            angle[i] = i * dToR;
            xOff[i] = radius * Math.cos(angle[i]);
            yOff[i] = radius * Math.sin(angle[i]);
            hFrac[i] = xOff[i] - Math.floor(xOff[i]);
            vFrac[i] = yOff[i] - Math.floor(yOff[i]);
        }
         
        Set<PairInt> checkForDisconnects = new HashSet<PairInt>();
        
        // run through the image interpolating grey values on each side
        // of the centre pixel to be used for the non-maximal suppression
        for (int row = (iRadius + 1); row < (nRows - iRadius); ++row) {
            for (int col = (iRadius+1); col < (nCols - iRadius) ; ++col) {
                
                int or = (int)orientation[col][row];
                double x = col + xOff[or];
                double y = row - yOff[or];
                
                //Get integer pixel locations that surround location x,y
                int fx = (int)Math.floor(x);
                int cx = (int)Math.ceil(x);
                int fy = (int)Math.floor(y);
                int cy = (int)Math.ceil(y);
                
                // top left
                double tl = img[fx][fy];
                // top right
                double tr = img[cx][fy];
                // bottom left
                double bl = img[fx][cy];
                // bottom right
                double br = img[cx][cy];
                
                // use bilinear interpolation to estimate value at x,y
                double upperAvg = tl + hFrac[or] * (tr - tl);
                double lowerAvg = bl + hFrac[or] * (br - bl);
                double v1 = upperAvg + vFrac[or] * (lowerAvg - upperAvg);
if (((col==145 || col==146) && row==83)) {
    int z = 1;//v1=0.106  img=0.059    for col=146, img=0.2437  v1=0.319
}
                if ((img[col][row] > v1) ||
                    (useLowerThreshold && (img[col][row] > 0.95*v1))) {
                    //x, y location on the `other side' of the point in question
                    x = col - xOff[or];
                    y = row + yOff[or];
                
                    fx = (int)Math.floor(x);
                    cx = (int)Math.ceil(x);
                    fy = (int)Math.floor(y);
                    cy = (int)Math.ceil(y);
                
                    // top left
                    tl = img[fx][fy];
                    // top right
                    tr = img[cx][fy];
                    // bottom left
                    bl = img[fx][cy];
                    // bottom right
                    br = img[cx][cy];
                
                    // use bilinear interpolation to estimate value at x,y
                    upperAvg = tl + hFrac[or] * (tr - tl);
                    lowerAvg = bl + hFrac[or] * (br - bl);
                    double v2 = upperAvg + vFrac[or] * (lowerAvg - upperAvg);
                     
                    if (useLowerThreshold && (img[col][row] <= v2)) {
                        
                        checkForDisconnects.add(new PairInt(col, row));
                    
                    } else if (img[col][row] > v2) {
                        
                        // this is the local maximum
                        output[col][row] = img[col][row];
                        
                        /* if wanted to create the localization image:
                            // Solve for coefficients of parabola that passes through 
                            // [-1, v1]  [0, inimage] and [1, v2]. 
                            // v = a*r^2 + b*r + c
                            c = img[col][row];
                            a = (v1 + v2)/2 - c;
                            b = a + c - v1;
        
                            // location where maxima of fitted parabola occurs
                            r = -b/(2*a);
                            location[col][row] = new Complex(col - r * xOff[or],
                                row + r * yOff[or]);
                        */
                    }
                    
                } else if (useLowerThreshold && (img[col][row] > 0.8*v1)) {
                    checkForDisconnects.add(new PairInt(col, row));
                }
            }
        }
               
        PairInt[][] neighborCoordOffsets
            = AbstractLineThinner.createCoordinatePointsForEightNeighbors(
                0, 0);
        
        for (PairInt p : checkForDisconnects) {
            int i = p.getX();
            int j = p.getY();
            
            boolean disconnects = useLowerThreshold ? 
                ImageSegmentation.doesDisconnect(output, 
                    neighborCoordOffsets, i, j) : false;
                       
            if (disconnects) {
                output[i][j] = img[i][j];
            }
        }
        
        /*
        // Finally thin the 'nonmaximally suppressed' image by pointwise
        // multiplying itself with a morphological skeletonization of itself.

        // I know it is oxymoronic to thin a nonmaximally supressed image but 
        // fixes the multiple adjacent peaks that can arise from using a radius
        // value > 1.
        */
        
        // operating on output so that can add a check for whether a new 0
        // disconnects a line
        
        // make binary image for bwmorph input
        int[][] morphInput = new int[output.length][];
        for (int i = 0; i < output.length; ++i) {
            morphInput[i] = new int[output[0].length];
        }
        for (int i = 0; i < output.length; ++i) {
            for (int j = 0; j < output[i].length; ++j) {
                                
                if (output[i][j] > 0) {
                    
                    morphInput[i][j] = 1;
                    
                } else {

                    // make sure does not add a gap into a curve
                    //boolean disconnects = ImageSegmentation.doesDisconnect(output,
                    //    neighborCoordOffsets, i, j);
                    //if (!disconnects) {
                        morphInput[i][j] = 0;
                        output[i][j] = 0;
                    //}
                }
            }
        }
        
        MorphologicalFilter mFilter = new MorphologicalFilter();
        int[][] skel = mFilter.bwMorphThin(morphInput, Integer.MAX_VALUE);
        
        for (int i = 0; i < output.length; ++i) {
            for (int j = 0; j < output[i].length; ++j) {
                int m = skel[i][j];
if ((i==145 && j == 152) || (i==146 && j==152)) {
     int z = 1;
}                
                output[i][j] *= m;
            }
        }
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        curveHelper.additionalThinning45DegreeEdges(orientation, output);
        
        return output;
    }
}
