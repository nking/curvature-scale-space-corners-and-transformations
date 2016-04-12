package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;

/**
 * a version of non-maximum suppression compatible with the 
 * PhaseConguencyDetector, that is accepts double arrays or integers 
 * and expects for the double arrays a convention used in
 * the PhaseConguencyDetector of format a[row][col].
% 
% @author nichole
 */
public class NonMaximumSuppression {
    
    /**
     * a non-maximal suppression algorithm that keeps the inspected pixel
     * if it is the maximum within the radius and direction given by
     * the orientation image.
     * Note that orientation and img arrays are using notation
     * [row][col] and that orientation values should be between 0 and 180.
     * 
     * @param img double array using format a[row][col] of an image to be
     * thinned using non-maximum suppression.
     * @param orientation double array using format a[row][col] of an image
     * holding values between 0 and 180 in a counter clockwise reference
     * frame.
     * @param radius suggested value is 1.2 to 1.5
     * @param outputCandidateJunctionsToRestore the points removed within
     * threshold range are put into this set which can be tested after
     * 2-layer filter to see if restoring the pixels would restore their
     * values.
     * @return 
     */
    public double[][] nonmaxsup(double[][] img, double[][] orientation, 
        double radius, Set<PairInt> outputCandidateJunctionsToRestore) {
        
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
        double[] i0Off = new double[180];
        double[] i1Off = new double[180];
        for (int i = 0; i < 180; ++i) {
            angle[i] = i * dToR;
            i1Off[i] = radius * Math.abs(Math.cos(angle[i]));
            i0Off[i] = radius * Math.abs(Math.sin(angle[i]));
        }
        
        /*
        Where orientation is closer to 90, it is a horizontal line,
            want to compare pixels above and below
            (dx=r*cos(90 in radians), dy=r*sin(90 in radians))
        
        When orientation is closer to 180 or 0, it is a vertical line,
            want to compare pixels to the either side
            (dx=r*cos(90 in radians), dy=r*sin(90 in radians))
        
        This diagonal is within 22.5 degrees of 45:
                  #
               #
            #
            (dx=r*cos(45 in radians), mult by -1* dy=-r*sin(45 in radians))
            where +y is up and +y is to right in sketch
        
        This diagonal is within 22.5 degrees of 135:
            #
               #
                  #
            (dx= r*cos(45 in radians), dy=r*sin(45 in radians))
            (mult dx by +1)
            where +y is up and +y is to right in sketch
        */
        
        // run through the image interpolating grey values on each side
        // of the centre pixel to be used for the non-maximal suppression
        for (int i1 = (iRadius + 1); i1 < (n1 - iRadius); ++i1) {
            for (int i0 = (iRadius+1); i0 < (n0 - iRadius) ; ++i0) {
                
                double v = img[i0][i1];
                
                if (v == 0) {
                    continue;
                }
                
                int or = (int)orientation[i0][i1];
                if (or > 179) {
                    or -= 180;
                }
                                
                int di0 = (int)Math.round(i0Off[or]);
                int di1 = (int)Math.round(i1Off[or]);
                int i0start = i0 - di0;
                int i1start = i1 - di1;
                int i0end = i0 + di0;
                int i1end = i1 + di1;
                boolean is45 = (Math.abs(or - 45) < 22.5);
                if (is45) {
                    int swap = i0start;
                    i0start = i0end;
                    i0end = swap;
                    di0 *= -1;
                }
                
                int ii0 = i0start;
                int ii1 = i1start;
                boolean isMax = true;
                while ((ii1 <= i1end) && (
                    (!is45 && (ii0 <= i0end)) || (is45 && (ii0 >= i0end)))){
                    
                    if ((i0 == ii0 && i1 == ii1) || (ii1 < 0) || (ii1 > (n1 - 1)) 
                        || (ii0 < 0) || (ii0 > (n0 - 1))) {
                        ii0 += di0;
                        ii1 += di1;
                        continue;
                    }
                    if (img[ii0][ii1] > v) {
                        isMax = false;
                        break;
                    }
                    ii0 += di0;
                    ii1 += di1;
                }
                if (isMax) {
                    output[i0][i1] = img[i0][i1];
                }
            }
        }
        
        int[][] morphInput = new int[n0][];
        for (int i = 0; i < n0; ++i) {
            morphInput[i] = new int[n1];
        }
        for (int i = 0; i < n0; ++i) {
            for (int j = 0; j < n1; ++j) {
                if (output[i][j] > 0) {                    
                    morphInput[i][j] = 1;
                } else {
                    morphInput[i][j] = 0;
                    output[i][j] = 0;
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
     * values in range 0 to 180.
     * 
     * @param img
     * @param orientation
     * @param radius
     * @param outputCandidateJunctionsRemoved
     */
    public void nonmaxsup(GreyscaleImage img, GreyscaleImage orientation, 
        double radius, Set<PairInt> outputCandidateJunctionsRemoved) {
        
        if (img.getWidth() != orientation.getWidth() || 
            img.getHeight() != orientation.getHeight()) {
            throw new IllegalArgumentException("img and orientation must be same size");
        }
        
        if (radius < 1) {
            throw new IllegalArgumentException("radius must be >= 1");
        }
        
        int n0 = img.getWidth();
        int n1 = img.getHeight();

        double[][] a = new double[n1][];
        double[][] or = new double[n1][];
        for (int i = 0; i < n1; ++i) {
            a[i] = new double[n0];
            or[i] = new double[n0];
            for (int j = 0; j < n0; ++j) {
                a[i][j] = img.getValue(j, i);
                or[i][j] = orientation.getValue(j, i);
            }
        }
                
        double[][] thinned = nonmaxsup(a, or, radius, 
            outputCandidateJunctionsRemoved);
        
        // apply thinning to the image
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n0; ++j) {
                if (thinned[i][j] == 0) {
                    img.setValue(j, i, 0);
                }
            }
        }
    }
    
}
