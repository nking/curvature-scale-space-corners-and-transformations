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
     * a non-maximal implementation that expects orientation image to have
     * values in range 0 to 180.  img can have any set of values.
     * 
     * @param img
     * @param orientation
     * @param radius
     * @param useLowerThreshold if true, uses a lower threshold when checking
     * whether to keep points - the lower threshold is useful for example, when
     * thinning the steps from the phase angle image.
     * @param outputCandidateJunctionsToRestore the points removed within
     * threshold range are put into this set which can be tested after
     * 2-layer filter to see if restoring the pixels would restore their
     * values.
     * @return 
     */
    public double[][] nonmaxsup(double[][] img, double[][] orientation, 
        double radius, boolean useLowerThreshold, 
        Set<PairInt> outputCandidateJunctionsToRestore) {
        
        if (img.length != orientation.length || img[0].length != orientation[0].length) {
            throw new IllegalArgumentException("img and orientation must be same size");
        }
                
        if (radius < 1) {
            throw new IllegalArgumentException("radius must be >= 1");
        }

        int n0 = img.length;
        int n1 = img[0].length;
        
        double[][] output = new double[n0][];
        for (int i = 0; i < n0; ++i) {
            output[i] = new double[n1];
        }
        
        int iRadius = (int)Math.ceil(radius);
        
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
                 
        // run through the image interpolating grey values on each side
        // of the centre pixel to be used for the non-maximal suppression
        for (int i1 = (iRadius + 1); i1 < (n1 - iRadius); ++i1) {
            for (int i0 = (iRadius+1); i0 < (n0 - iRadius) ; ++i0) {
                
                int or = (int)orientation[i0][i1];
                if (or > 179) {
                    or -= 180;
                }
                double ii0 = i0 + xOff[or];
                double ii1 = i1 - yOff[or];
                
                //Get integer pixel locations that surround location x,y
                int f0 = (int)Math.floor(ii0);
                int c0 = (int)Math.ceil(ii0);
                int f1 = (int)Math.floor(ii1);
                int c1 = (int)Math.ceil(ii1);
                
                // top left
                double tl = img[f0][f1];
                // top right
                double tr = img[c0][f1];
                // bottom left
                double bl = img[f0][c1];
                // bottom right
                double br = img[c0][c1];
                
                // use bilinear interpolation to estimate value at x,y
                double upperAvg = tl + hFrac[or] * (tr - tl);
                double lowerAvg = bl + hFrac[or] * (br - bl);
                double v1 = upperAvg + vFrac[or] * (lowerAvg - upperAvg);
                
                if ((img[i0][i1] > v1) ||
                    (useLowerThreshold && (img[i0][i1] > 0.95*v1))) {
                    //x, y location on the `other side' of the point in question
                    ii0 = i0 - xOff[or];
                    ii1 = i1 + yOff[or];
                
                    f0 = (int)Math.floor(ii0);
                    c0 = (int)Math.ceil(ii0);
                    f1 = (int)Math.floor(ii1);
                    c1 = (int)Math.ceil(ii1);
                
                    // top left
                    tl = img[f0][f1];
                    // top right
                    tr = img[c0][f1];
                    // bottom left
                    bl = img[f0][c1];
                    // bottom right
                    br = img[c0][c1];
                
                    // use bilinear interpolation to estimate value at x,y
                    upperAvg = tl + hFrac[or] * (tr - tl);
                    lowerAvg = bl + hFrac[or] * (br - bl);
                    double v2 = upperAvg + vFrac[or] * (lowerAvg - upperAvg);
                    
                    if (useLowerThreshold && (img[i0][i1] <= v2)) {
                        
                        outputCandidateJunctionsToRestore.add(new PairInt(i0, i1));
                    
                    } else if (img[i0][i1] > v2) {
                        
                        // this is the local maximum
                        output[i0][i1] = img[i0][i1];
                        
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
                    } else {
                        outputCandidateJunctionsToRestore.add(new PairInt(i0, i1));
                    }
                    
                } else if (useLowerThreshold && (img[i0][i1] > 0.8*v1)) {
                    outputCandidateJunctionsToRestore.add(new PairInt(i0, i1));
                }
            }
        }
               
        PairInt[][] neighborCoordOffsets
            = AbstractLineThinner.createCoordinatePointsForEightNeighbors(
                0, 0);
        
        Set<PairInt> addedBack = new HashSet<PairInt>();
        for (PairInt p : outputCandidateJunctionsToRestore) {
            int i0 = p.getX();
            int i1 = p.getY();
             
            boolean disconnects = useLowerThreshold ? 
                ImageSegmentation.doesDisconnect(output, 
                    neighborCoordOffsets, i0, i1) : false;
                       
            if (disconnects) {
                output[i0][i1] = img[i0][i1];
                addedBack.add(p);
            }
        }
        
        outputCandidateJunctionsToRestore.removeAll(addedBack);
        
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
        int[][] morphInput = new int[n0][];
        for (int i = 0; i < n0; ++i) {
            morphInput[i] = new int[n1];
        }
        for (int i = 0; i < n0; ++i) {
            for (int j = 0; j < n1; ++j) {
                                
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
        
        for (int i = 0; i < n0; ++i) {
            for (int j = 0; j < n1; ++j) {
                int m = skel[i][j];                 
                output[i][j] *= m;
            }
        }
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        curveHelper.additionalThinning45DegreeEdges2(orientation, output);
        
        return output;
    }
    
    /**
     * a non-maximal implementation that expects orientation image to have
     * values in range 0 to 180.  img can have any set of values.
     * 
     * @param img
     * @param orientation
     * @param radius
     * @param useLowerThreshold if true, uses a lower threshold when checking
     * whether to keep points - the lower threshold is useful for example, when
     * thinning the steps from the phase angle image.
     * @param outputCandidateJunctionsRemoved
     */
    public void nonmaxsup(GreyscaleImage img, GreyscaleImage orientation, 
        double radius, boolean useLowerThreshold, 
        Set<PairInt> outputCandidateJunctionsRemoved) {
        
        if (img.getWidth() != orientation.getWidth() || 
            img.getHeight() != orientation.getHeight()) {
            throw new IllegalArgumentException("img and orientation must be same size");
        }
        
        if (radius < 1) {
            throw new IllegalArgumentException("radius must be >= 1");
        }
        
        int n0 = img.getWidth();
        int n1 = img.getHeight();

        double[][] a = new double[n0][];
        double[][] or = new double[n0][];
        for (int i = 0; i < n0; ++i) {
            a[i] = new double[n1];
            or[i] = new double[n1];
            for (int j = 0; j < n1; ++j) {
                a[i][j] = img.getValue(i, j);
                or[i][j] = orientation.getValue(i, j);
            }
        }
                
        double[][] thinned = nonmaxsup(a, or, radius, useLowerThreshold, 
            outputCandidateJunctionsRemoved);
        
        // apply thinning to the image
        for (int i = 0; i < n0; ++i) {
            for (int j = 0; j < n1; ++j) {
                if (thinned[i][j] == 0) {
                    img.setValue(i, j, 0);
                }
            }
        }
    }
    
}
