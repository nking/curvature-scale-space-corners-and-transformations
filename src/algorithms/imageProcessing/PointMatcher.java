package algorithms.imageProcessing;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.compGeometry.PointPartitioner;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.ArrayPair;
import algorithms.util.Errors;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Logger;
import thirdparty.HungarianAlgorithm;
import org.ejml.simple.*;

/**
 * class to match the points extracted from two images.
 * 
 * the transformation parameters of translation, rotation and scale are 
 * found given the two sets of points.
 * 
 <pre>
   Details of determining the transformation from points.
   
   After points are matched:
       Estimating scale:
           Take 2 pairs of points in both datasets and compute the distance 
           between them and then take the ratio:

           scale = (distance between pair in set 1) 
                   / (distance between pair in set 2)
      
           Note: scale is roughly determined from contour matching too in the
              inflection matcher.
      
       Estimating rotation:
           Take the same 2 pairs and determine the difference in their angles:
               tan(theta) = delta y / delta x
      
           rotation = atan((delta y between pair in set 1)
                          /(delta x between pair in set 1)) 
                      -
                      atan((delta y between pair in set 2)
                          /(delta x between pair in set 2))
      
       Estimate translation:
           Performed on one point in set 1 with its candidate match in set 2:
           From the full transformation equation, we can rewrite:
               transX = xt0 - xc*scale - 
                   (((x0-xc)*scale*math.cos(theta)) + ((y0-yc)*scale*math.sin(theta)))
       
               transY = yt0 - yc*scale - 
                    ((-(x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))
       
               where (xc, yc) is the center of the first image
 
 Matching the points:
 
   For the contour matching, we were able to choose subsets of the entire
   sets of contours for better matching characteristics (subsets that had larger
   sigma values for their peaks could be used).
   
   For the points given here, there isn't a clear indicator for a subset that 
   could be used preferably so that all nodes might not need to be visited.
   
   The need to use pairs of points in set1 matched against pairs in set 2
   means that if one tries every combination of pairs, the runtime complexity 
   is exponential.
   
       number of ways to make pairs in set 1 times the number of ways to make
       pairs in set 2 = 
             n_1!            n_2!
         ------------  X  ------------
         2*(n_1 - 2)!     2*(n_2 - 2)!
    
   This is not feasible as soon as the sets get larger than a dozen points.
    
   Alternately, one could try a search algorithm to visit the space of all 
   possible combinations of parameters instead of all combinations of subsets 
   of pairs.
       
        max translation x possible = width of image 2
        max translation y possible = height of image 2
        max scale is set to some reasonable number like 10?
        max rotation is 1 degree changes = 360 (or quadrants?)
    
        then total number of permutations 
            = maxTransX * maxTransY * maxScale * 360
            = 1000 * 500 * 10 * 360 = 1.8 billion for pure grid search over 
              parameter space
  
        then trying the parameters on all points puts another factor into there
        of nPoints set 1 * nPoints set 2.
        
  Note, because the cosine and sine terms in the transformation equation due to
  rotation work against one another and don't proceed in a total combined
  single direction with theta, we can't use a gradient descent solver.
 
       (note, this hasn't been edited for positive Y up yet)
       positive Y is down 
       positive X is right
       positive theta starts from Y=0, X>=0 and proceeds CW
                270
           QIII  |  QIV   
                 |
          180--------- 0   +X
                 |   
           QII   |   QI
                 90
                 +Y
           
  theta        cos(theta)  sin(theta)
  0     0        1.0         0    -----
  30    0.5236   0.87        0.5       |
  45    0.785    0.707       0.707     | QI
  60    1.047    0.5         0.866     |
  90    1.57     0.0         1.0  -----
  120   2.09    -0.5         0.866     | QII
  150   2.618   -0.866       0.5       |
  180   3.1416  -1.0         0.0  -----
  210   3.6652  -0.866      -0.5       | QIII
  240                                  |
  270   4.7124   0.0        -1.0  -----
  300   5.236    0.5        -0.866     | QIV
  330                                  |
                                  -----
  So then, it looks like a grid search over rotation and scale intervals
  followed by fitting for translation to get into the local neighborhood
  of the transformation solution, followed by the use of the 
  Nelder-Mead Downhill Simplex to refine the transformation is a better 
  solution.
  
  runtime complexity of:
  
      grid search: nRotationIntervals * nScaleIntervals * nSet1 * nSet2 
  
      downhill search: not polynomial nor deterministic.  varies by dataset.
 
     In the searches, if rotation and scale are fixed, and transX and transY 
     are to be found, one can use either:

         (1) compare the offsets implied by each combination of points in set 1
             and set 2.
             for brown_lowe_2003 image1 and image2, the number of corners is
             n1 = 78
             n2 = 88
             n1 * n2 = 5616

             How does one determine the best solution among those translations?  
             One needs an assumed tolerance for the translation, and then to 
             count the number of matches to a point in the 2nd set.
             The best solution has the highest number of matches with the same
             rules for an assumed tolerance and the smallest avg difference
             with predicted positions for matches in set 2.
              
             For the tolerance, the reasons that translation might be different
             as a function of position in the image might be:
               -- due to rounding to a pixel.  these are very small errors.
               -- due to projection effects for different epipolar geometry 
                  (camera nadir perpendicular to different feature, for example).  
                  these are potentially very large and are not
                  solved in the transformation with this point matcher though a
                  tolerance is made for a small amount of it.  
                  (NOTE however, that once the points are matched, the true 
                  epipolar geometry can be calculated with the StereoProjection 
                  classes).
               -- due to errors in the location of the corner.  this can be due
                  to the edge detector.  these are small errors.
               -- due to image warping such as distortions from the shape of the
                  lens.  one would need a point matcher tailored for the 
                  specific geometric projection.
               -- due to a camera not completely perpendicularly aligned with
                  the optical axis. presumably, this is an extreme case and
                  you'd want better data...
             
             For the Brown & Lowe 2003 points, 
                 transX=293.1 (stdDev=10.3)
                 transY=14.3 (stdDev=5.9)
             the large standard deviation appears to be due to projection 
             effects.  the skyline is rotated about 13 degrees w.r.t. skyline 
             in second image while the features at the bottom remain horizontal 
             in both.
             
             The spread in standard deviation appears to be correlated w/
             the image dimensions, as would be expected with projection being 
             the largest reason for a spread in translation.
             That is, the X axis is twice the size of the Y and so are their
             respective standard deviations.
             
             Could make a generous allowance for projection effects by assuming 
             a maximum present such as that in the Brown & Lowe images, that is 
             image size times an error due to a differential spread 
             over those pixels for a maximum difference in rotation such as 
             20 degrees or something.  For Brown & Lowe images, the tolerance 
             would be number of pixels for dimension times 0.02.
             
             If there is no reasonable solution using only scale, rotation, and 
             translation, then a more computationally expensive point matcher 
             that solves for epipolar geometry too or a differential rotation
             and associated variance in other parameters
             is needed (and hasn't been implemented 
             here yet).  **OR, even better, an approach using contours and an
             understanding of occlusion (the later possibly requires shape
             identification) can be made with the contour matcher in this 
             project.**
             The contour matcher approach is currently not **commonly**
             possible with this project, because most edges are not closed 
             curves.

    OR

         (2) try every possible combination of translation for X and Y which
             would be width of image 2 in pixels times the height of image 2
             in pixels.
             for brown_lowe_2003 image1 and image2,
             image2 width = 517, height = 374
             w * h = 193358
             so this method is 34 times more

             The best solution among those translations is found in the same
             way as (1) above.
 </pre>

 * @author nichole
 */
public final class PointMatcher {
     
    private final Logger log = Logger.getLogger(this.getClass().getName());
 
    /**
     * 
     * calculate the rotation, scale, and translation that can be applied
     * to set1 to match points to set2 where the matches are not already known.  
     * The image1Width and imageHeight are used to create a tolerance in 
     * translation for matches.
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * @param set1 set of points such as corners from an image
     * @param set2 set of points such as corners from another image.
     * @param image1Width
     * @param image1Height
     * @return 
     */
    public TransformationPointFit calculateTransformation(PairIntArray set1, 
        PairIntArray set2, int image1Width, int image1Height) {
        
        boolean setsAreMatched = false;
        
        int rotDelta = 90;
        
        TransformationPointFit fit = calculateTransformationWithGridSearch(
            set1, set2, image1Width, image1Height,
            0, 360, rotDelta,
            1, 10, 1, setsAreMatched, set1, set2
            );
        
        if (fit == null) {
            return null;
        }
       
        int scale = (int)Math.round(fit.getScale());
        int rot = (int)Math.round(fit.getRotationInRadians());
        
        int rotStart = rot - (rotDelta/2);
        if (rotStart < 0) {
            rotStart += 360;
        }
        int rotStop = rot + (rotDelta/2);
        if (rotStop > 360) {
            rotStop -= 360;
        }
        if (rotStop < rotStart) {
            rotStop += 360;
        }
        
        int scaleStart = scale - 1;
        if (scaleStart < 1) {
            scaleStart = 1;
        }

        fit = calculateTransformationWithGridSearch(
            set1, set2, image1Width, image1Height,
            rotStart, rotStop, 10,
            scaleStart, scale + 1, 1, 
            setsAreMatched, set1, set2
            );
        
        if (fit == null) {
            return null;
        }
                  
        int convergence = (set1.getN() < set2.getN()) ? set1.getN() 
            : set2.getN();
        
        if ((fit.getNumberOfMatchedPoints() == convergence)
            && (fit.getMeanDistFromModel() == 0)) {
            return fit;
        }

        TransformationPointFit fit2 = refineTransformationWithDownhillSimplex(
            fit.getParameters(),
            set1, set2, image1Width, image1Height, setsAreMatched);
        
        if (fitIsBetter(fit, fit2)) {
            return fit2;
        }
        
        return fit;
    }
    
    /**
     * 
     * calculate the rotation, scale, and translation that can be applied
     * to set1 to match points to set2 for unmatched points.  The method
     * is useful for points which may be part of stereo projections or other
     * projections more complex than euclidean for which euclidean can
     * help with rough registration of data as a first step, but is not
     * the precise solution.
     * 
     * The image1Width and imageHeight are used to create a tolerance in 
     * translation for matches.
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * @param set1 set of points such as corners from an image
     * @param set2 set of points such as corners from another image.
     * @param image1Width
     * @param image1Height
     * @return 
     */
    public TransformationPointFit calculateRoughTransformationForUnmatched(
        PairIntArray set1, PairIntArray set2, int image1Width, 
        int image1Height) {
        
        /*
        divides the points into partitions for each set.
        solves each combination of the partitions quadrants.
        applies each of the transformations to all points and takes
        the best fitting transformation.
        */
        
        PairIntArray input1 = set1.copy();
        PairIntArray input2 = set2.copy();
        
        int rotDelta = 10;
        int scaleStop = 6;
        boolean setsAreMatched = false; 
        
        PointPartitioner partitioner = new PointPartitioner();
        
        TransformationPointFit bestFit = null;
        TransformationPointFit lastBestFit = null;
        int bestP1 = -1;
        int bestP2 = -1;
        
        int nD = 2;
        int nIter = 0;
        int nIterMax = 1;//(input1.getN()/20)/(nD * nD);
        
        while (nIter < nIterMax) {
            
            PairIntArray[] partitions1 = partitioner.partition(input1, nD);
        
            PairIntArray[] partitions2 = partitioner.partition(input2, nD);
        
            TransformationPointFit[] fits 
                = new TransformationPointFit[
                    partitions1.length*partitions2.length];
            
            int count = 0;

            for (int i = 0; i < partitions1.length; i++) {

                PairIntArray edge1 = partitions1[i];
                
                for (int j = 0; j < partitions2.length; j++) {
//if (!(((i==1)&&(j==0)) || ((i==3)&&(j==2)))){ continue;}                
                    PairIntArray edge2 = partitions2[j];

                    //TODO: need to alter the fit for partition solutions
                    // because sets have different number of points
                    
                    fits[count] = calculateTransformationWithGridSearch(
                        edge1, edge2, image1Width, image1Height,
                        0, 360, rotDelta,
                        1, scaleStop, 1, setsAreMatched, set1, set2
                    );

                    if (fitIsBetter(bestFit, fits[count])) {
                        bestFit = fits[count];
                        bestP1 = i;
                        bestP2 = j;
                    }

                    count++;
                }                
            }
            
            sortByDescendingMatches(fits, 0, (fits.length - 1));
                
            int z = 1;
                
            if (bestFit == null) {
                return null;
            }
            
            if (lastBestFit == bestFit) {
                break;
            }
            
            log.info("bestFit=" + bestFit.toString());
                        
            if ((partitions1[bestP1].getN() < 20) || (partitions2[bestP2].getN() < 20)) {
                break;
            }
            
            input1 = partitions1[bestP1];
            input2 = partitions2[bestP2];
            
            lastBestFit = bestFit;
            
            nIter++;
        }
        
        /*
        keep track of the quadrant that the best fit occurred in.
        
        -- one more divide by 4 pass if enough points?
        -- 
        */
                
        input1 = set1.copy();
        input2 = set2.copy();
        
        if (false) {
            
            // ==== match points and exclude the largest outliers ====
            
            double tolerance = Double.MAX_VALUE;
            
            float[][] matchIndexesAndDiffs = calculateMatchUsingOptimal(
                input1, input2, bestFit, tolerance,
                image1Width >> 1, image1Height >> 1);

            tolerance = findToleranceForOutliers(
                bestFit, matchIndexesAndDiffs,
                image1Width >> 1, image1Height >> 1);

            PairIntArray outputMatched1 = new PairIntArray();
            PairIntArray outputMatched2 = new PairIntArray();
            
            matchPoints(input1, input2, tolerance, matchIndexesAndDiffs,
                outputMatched1, outputMatched2);
            
            input1 = outputMatched1.copy();
            input2 = outputMatched2.copy();
            
            setsAreMatched = true;
        }
                       
        //TODO:  run tests to see whether this last is helpful:
        
        int convergence = (input1.getN() < input2.getN()) 
            ? input1.getN() 
            : input2.getN();

        if ((bestFit.getNumberOfMatchedPoints() == convergence)
            && (bestFit.getMeanDistFromModel() == 0)) {
            return bestFit;
        }

        bestFit = refineTransformationWithDownhillSimplex(
            bestFit.getParameters(), input1, input2, image1Width, image1Height, 
            setsAreMatched);
        
        return bestFit;
    }
    
    /**
     * 
     * calculate the rotation, scale, and translation that can be applied
     * to set1 to match points to set2 where the matches are not already known.  
     * The image1Width and imageHeight are used to create a tolerance in 
     * translation for matches.
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * @param set1 set of points such as corners from an image
     * @param set2 set of points such as corners from another image.
     * @param image1Width
     * @param image1Height
     * @param outputMatched1
     * @param outputMatched2
     * @return 
     */
    public TransformationPointFit calculateTransformationAndMatch(
        PairIntArray set1, PairIntArray set2, int image1Width, int image1Height,
        PairIntArray outputMatched1, PairIntArray outputMatched2) {
        
        /*
        attempts rough euclidean transformation and uses that to match
        points in set1 and set2, then removes outliers from the matched set.
        */
        
        boolean setsAreMatched = false;
        
        int rotDelta = 90;
        
        TransformationPointFit fit = calculateTransformationWithGridSearch(
            set1, set2, image1Width, image1Height,
            0, 360, rotDelta,
            1, 11, 2, setsAreMatched, set1, set2
            );
        
        if (fit == null) {
            return null;
        }
        
        int nIter = 0;
        int nMaxIter = 2;
        
        boolean repeat = true;
        
        PairIntArray input1 = set1;
        PairIntArray input2 = set2;
            
        while ((nIter < nMaxIter) && repeat) {
            
            int scale = (int)Math.round(fit.getScale());
            int rot = (int)Math.round(fit.getRotationInRadians());

            // if scale ~ 1 and rotation ~ 0, use an early euclidean match 

            float rotComp = (rot == 0) ? 0 :
                (fit.getParameters().getRotationInDegrees() <= 180)
                ? fit.getParameters().getRotationInDegrees() : 
                fit.getParameters().getRotationInDegrees() - 360.f;
        
            if ((scale == 1) && (Math.abs(rotComp) < 5) && !setsAreMatched) {
            
                // ==== match points and exclude the largest outliers ====
                
                double tolerance = Double.MAX_VALUE;
                
                float[][] matchIndexesAndDiffs = calculateMatchUsingOptimal(
                    input1, input2, fit, tolerance,
                    image1Width >> 1, image1Height >> 1);
                    
                tolerance = findToleranceForOutliers(
                    fit, matchIndexesAndDiffs,
                    image1Width >> 1, image1Height >> 1);
                
                matchPoints(input1, input2, tolerance, matchIndexesAndDiffs,
                    outputMatched1, outputMatched2);
                
                input1 = outputMatched1.copy();
                input2 = outputMatched2.copy();
                                
                // rematch remaining points
                float[][] matchIndexesAndDiffs2 = calculateMatchUsingOptimal(
                    input1, input2, fit, tolerance,
                    image1Width >> 1, image1Height >> 1);
                
                // clear out outputMatched1 and 2
                PairIntArray tmp = new PairIntArray();
                outputMatched1.swapContents(tmp);
                tmp = new PairIntArray();
                outputMatched2.swapContents(tmp);
                
                matchPoints(input1, input2, tolerance, matchIndexesAndDiffs2,
                    outputMatched1, outputMatched2);
                
                input1 = outputMatched1;
                input2 = outputMatched2;
                
                setsAreMatched = true;
                
                repeat = false;
            }

            int rotStart = 0;
            int rotStop = 360;

            int scaleStart = scale - 1;
            if (scaleStart < 1) {
                scaleStart = 1;
            }
        
            fit = calculateTransformationWithGridSearch(
                input1, input2, image1Width, image1Height,
                rotStart, rotStop, 10,
                scaleStart, scale + 1, 1, 
                setsAreMatched, input1, input2
                );

            if (fit == null) {
                return null;
            }
            
            int convergence = (input1.getN() < input2.getN()) 
                ? input1.getN() 
                : input2.getN();

            if ((fit.getNumberOfMatchedPoints() == convergence)
                && (fit.getMeanDistFromModel() == 0)) {
                return fit;
            }
        
            fit = refineTransformationWithDownhillSimplex(fit.getParameters(),
                input1, input2, image1Width, image1Height, 
                setsAreMatched);
        
            if (!setsAreMatched) {
                
                // ==== match points and exclude the largest outliers ====
                
                int nBefore = input1.getN();
                
                double tolerance = Double.MAX_VALUE;
                
                float[][] matchIndexesAndDiffs = calculateMatchUsingOptimal(
                    input1, input2, fit, tolerance,
                    image1Width >> 1, image1Height >> 1);
                    
                tolerance = findToleranceForOutliers(
                    fit, matchIndexesAndDiffs,
                    image1Width >> 1, image1Height >> 1);
                
                matchPoints(input1, input2, tolerance, matchIndexesAndDiffs,
                    outputMatched1, outputMatched2);
                
                input1 = outputMatched1;
                input2 = outputMatched2;
                
                setsAreMatched = true;
                
                if ((outputMatched1.getN() - nBefore) == 0) {
                    repeat = false;
                }
            }
            
            nIter++;
        }
        
        return fit;
    }

    public boolean fitIsBetter(TransformationPointFit bestFit, 
        TransformationPointFit compareFit) {
        
        if (compareFit == null) {
            return false;
        }
        if (bestFit == null) {
            return true;
        }
        
        /*
        double d = compareFit.getMeanDistFromModel();
                
        if (d < bestFit.getMeanDistFromModel()) {
            return true;
        } else if (d == bestFit.getMeanDistFromModel()) {
            if (compareFit.getNumberOfMatchedPoints()
                > bestFit.getNumberOfMatchedPoints()) {
                return true;
            } else if (compareFit.getNumberOfMatchedPoints()
                == bestFit.getNumberOfMatchedPoints()) {
                if (compareFit.getStDevFromMean() < bestFit.getStDevFromMean()) {
                    return true;
                }
            }
        }
        return false;
        */
        int nMatches = compareFit.getNumberOfMatchedPoints();
        
        if (nMatches > bestFit.getNumberOfMatchedPoints()) {
            return true;
        } else if (nMatches == bestFit.getNumberOfMatchedPoints()) {
            if (!Double.isNaN(compareFit.getMeanDistFromModel()) && (
                compareFit.getMeanDistFromModel()
                < bestFit.getMeanDistFromModel())) {
                return true;
            } else if (compareFit.getMeanDistFromModel()
                == bestFit.getMeanDistFromModel()) {
                if (compareFit.getStDevFromMean() < bestFit.getStDevFromMean()) {
                    return true;
                }
            }
        }
        return false;
    }
    
    /**
     * find the Euclidean transformation using a grid search to find the best
     * match between set1 and set2, but evaluate the fit by applying the
     * transformation to allPoints1 and comparing to allPoints2.
     * 
     * @param set1
     * @param set2
     * @param image1Width
     * @param image1Height
     * @param rotStart start of rotation search in degrees
     * @param rotStop stop (exclusive) or rotations search in degrees
     * @param rotDelta change in rotation to add to reach next step in rotation search in degrees
     * @param scaleStart
     * @param scaleStop
     * @param scaleDelta
     * @param setsAreMatched
     * @return 
     */
    private TransformationPointFit calculateTransformationWithGridSearch(
        PairIntArray set1, PairIntArray set2, 
        int image1Width, int image1Height,
        int rotStart, int rotStop, int rotDelta,
        int scaleStart, int scaleStop, int scaleDelta,
        boolean setsAreMatched, PairIntArray allPoints1, PairIntArray allPoints2
        ) {
        
        int centroidX1 = image1Width >> 1;
        int centroidY1 = image1Height >> 1;
        
        double tolTransX = 2.0f * image1Width * 0.02f;
        double tolTransY = 2.0f * image1Height * 0.02f;
        
        Transformer transformer = new Transformer();
        
        int convergence = (set1.getN() < set2.getN()) ? set1.getN() 
            : set2.getN();
        
        TransformationPointFit bestFit = null;
                
        for (int scale = scaleStart; scale < scaleStop; scale += scaleDelta) {
            for (int rot = rotStart; rot < rotStop; rot += rotDelta) {
            
                TransformationParameters params;
                
                if (setsAreMatched) {
                    params = calculateTranslationForMatched(set1, set2, 
                        rot*Math.PI/180., scale, centroidX1, centroidY1);                    
                } else {
                    params = calculateTranslationForUnmatched(set1, set2, 
                        rot*Math.PI/180., scale, centroidX1, centroidY1);
                }
                
                PairFloatArray allPoints1Tr = transformer.applyTransformation(
                    params, centroidX1, centroidY1, allPoints1);
                
                TransformationPointFit fit;
                
                if (setsAreMatched) {
                    fit = evaluateFitForMatchedTransformed(params, 
                        allPoints1Tr, allPoints2);                    
                } else {
                    fit = evaluateFitForUnmatchedTransformed(params, 
                        allPoints1Tr, allPoints2, tolTransX, tolTransY);
                }
        
                // correction for intersection area?
  
                if (fitIsBetter(bestFit, fit)) {
                    
log.info("==> " + " tx=" + fit.getTranslationX() + " ty=" + fit.getTranslationY()
+ " rot=" + rot + " scale=" + scale + " fit=" + fit.toString());

                    bestFit = fit;
                    
                    if ((bestFit.getNumberOfMatchedPoints() == convergence) 
                        && (bestFit.getMeanDistFromModel() == 0)) {
                        if (bestFit.getParameters().getRotationInRadians() > 2 * Math.PI) {
                            float rot2 = bestFit.getParameters().getRotationInRadians();
                            while (rot2 >= 2 * Math.PI) {
                                rot2 -= 2 * Math.PI;
                            }
                            bestFit.getParameters().setRotationInRadians(rot2);
                        }
                        return bestFit;
                    }
                }
            }
        }
        
        if ((bestFit != null) && (bestFit.getParameters().getRotationInRadians() 
            > 2.*Math.PI)) {
            float rot = bestFit.getParameters().getRotationInRadians();
            while (rot >= 2*Math.PI) {
                rot -= 2*Math.PI;
            }
            bestFit.getParameters().setRotationInRadians(rot);
        }
        
        return bestFit;
    }
    
    /**
     * apply the parameters to set1 and find the matches to points in set2
     * within the given tolerance for translations.
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * 
     * ALSO NOTE: if you know a better solution exists for translation 
     * parameters that matches fewer points, but has a small avg dist from
     * model and smaller standard deviation from the avg dist from model,
     * then transXTol and transYTol should be smaller.
     * @param edges1
     * @param edges2
     * @param params
     * @param tolTransX
     * @param tolTransY
     * @param centroidX1
     * @param centroidY1
     * @return 
     */
    TransformationPointFit transform(PairIntArray[] edges1, 
        PairIntArray[] edges2, TransformationParameters params, 
        double tolTransX, double tolTransY, int centroidX1, int centroidY1) {
        
        if (edges1 == null || edges2 == null) {
            throw new IllegalArgumentException(
            "neither edges1 nor edges2 can be null");
        }
        if (edges1.length != edges2.length) {
            throw new IllegalArgumentException(
            "edges1 must be the same length as edges2");
        }
                
        double scale = params.getScale();
        
        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");
            
            log.severe("scale cannot be smaller than 1");
            
            return null;
        }        
        
        int nMaxMatchable = 0;
        for (int i = 0; i < edges1.length; i++) {
            int n = edges1[i].getN();
            if (edges2[i].getN() > n) {
                n = edges2[i].getN();
            }
            nMaxMatchable += n;
        }
        
        double[] diffFromModel = new double[nMaxMatchable];
        double[] avgDiffModel = new double[1];
        
        int nMatched = 0;
        int avgDiffModelSum = 0;
        for (int i = 0; i < edges1.length; i++) {
            
            PairIntArray set1 = edges1[i];
            PairIntArray set2 = edges2[i];
            
            int n = populateDiffFromModel(set1, set2, params, 
                tolTransX, tolTransY, centroidX1, centroidY1,
                diffFromModel, nMatched, avgDiffModel);
                        
            avgDiffModelSum += (avgDiffModel[0] * n);
            
            nMatched += n;
            
        }       
        
        double finalAvgDiffModel = avgDiffModelSum / (double)nMatched;
        
        double stDevDiffModel = 0;
        for (int i = 0; i < nMatched; i++) {
            double d = diffFromModel[i] - avgDiffModel[0];
            stDevDiffModel += (d * d);
        }
        
        stDevDiffModel = (Math.sqrt(stDevDiffModel/(nMatched - 1.0f)));
        
        TransformationPointFit fit = new TransformationPointFit(params, 
            nMatched, finalAvgDiffModel, stDevDiffModel,
            Math.sqrt(tolTransX*tolTransX + tolTransY*tolTransY)
        );
        
        return fit;
    }

     /**
     * apply the parameters to set1 and find the matches to points in set2
     * within the given tolerance for translations.
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * 
     * ALSO NOTE: if you know a better solution exists for translation 
     * parameters that matches fewer points, but has a small avg dist from
     * model and smaller standard deviation from the avg dist from model,
     * then transXTol and transYTol should be smaller.
     * @param set1 points matched to set 2
     * @param set2 points matched to set 1
     * @param transXTol
     * @param transYTol
     * @param centroidX1
     * @param centroidY1
     * @return 
     */
    TransformationPointFit transform(PairIntArray set1, PairIntArray set2, 
        double scale, double rotationInRadians, double translationX,
        double translationY, double transXTol, double transYTol, 
        int centroidX1, int centroidY1) {
    
        TransformationParameters params = new TransformationParameters();
        params.setScale((float)scale);
        params.setRotationInRadians((float)rotationInRadians);
        params.setTranslationX((float)translationX);
        params.setTranslationY((float)translationY);
        
        return transform(set1, set2, params, transXTol, transYTol, centroidX1,
            centroidY1);
    }
    
    /**
     * apply the parameters to set1 and find the matches to points in set2
     * within the given tolerance for translations.
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * 
     * ALSO NOTE: if you know a better solution exists for translation 
     * parameters that matches fewer points, but has a small avg dist from
     * model and smaller standard deviation from the avg dist from model,
     * then transXTol and transYTol should be smaller.
     * @param set1 points matched to set 2
     * @param set2 points matched to set 1
     * @param params
     * @param tolTransX
     * @param tolTransY
     * @param centroidX1
     * @param centroidY1
     * @return 
     */
    TransformationPointFit transform(PairIntArray set1, PairIntArray set2, 
        TransformationParameters params, double tolTransX, double tolTransY, 
        int centroidX1, int centroidY1) {
        
        if (set1 == null || set2 == null) {
            throw new IllegalArgumentException(
            "neither set1 nor set2 can be null");
        }
                
        double scale = params.getScale();
        
        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");
            
            log.severe("scale cannot be smaller than 1");
            
            return null;
        }
        
        int nMaxMatchable = set1.getN();
        if (set2.getN() > nMaxMatchable) {
            nMaxMatchable = set2.getN();
        }
        
        double[] diffFromModel = new double[nMaxMatchable];
        double[] avgDiffModel = new double[1];
        
        int nMatched = populateDiffFromModel(set1, set2, 
            params, tolTransX, tolTransY, centroidX1, centroidY1,
            diffFromModel, 0, avgDiffModel);
        
        double stDevDiffModel = 0;
        for (int i = 0; i < nMatched; i++) {
            double d = diffFromModel[i] - avgDiffModel[0];
            stDevDiffModel += (d * d);
        }
        
        stDevDiffModel = (Math.sqrt(stDevDiffModel/(nMatched - 1.0f)));
        
        TransformationPointFit fit = new TransformationPointFit(params, 
            nMatched, avgDiffModel[0], stDevDiffModel,
            Math.sqrt(tolTransX*tolTransX + tolTransY*tolTransY)
        );
        
        return fit;
    }
    
    /**
     * apply the parameters to matched1 and return the fit
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * 
     * @param matched1 points matched to set 2
     * @param matched2 points matched to set 1
     * @param params
     * @param transXTol
     * @param transYTol
     * @param centroidX1
     * @param centroidY1
     * @return 
     */
    TransformationPointFit transformMatched(PairIntArray matched1, 
        PairIntArray matched2, 
        TransformationParameters params, 
        int centroidX1, int centroidY1) {
        
        if (matched1 == null || matched2 == null) {
            throw new IllegalArgumentException(
            "neither matched1 nor matched2 can be null");
        }
        if (matched1.getN() != matched2.getN()) {
            throw new IllegalArgumentException(
            "matched1 and matched2 have to be the same length");
        }
                
        double scale = params.getScale();
        
        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");
            
            log.severe("scale cannot be smaller than 1");
            
            return null;
        }
       
        double[] diffFromModel = new double[matched1.getN()];
        double[] avgDiffModel = new double[1];
        
        populateDiffFromModelWithMatchedPoints(matched1, matched2, 
            params, centroidX1, centroidY1, diffFromModel, avgDiffModel);
        
        double stDevDiffModel = 0;
        for (int i = 0; i < matched1.getN(); i++) {
            double d = diffFromModel[i] - avgDiffModel[0];
            stDevDiffModel += (d * d);
        }
        
        stDevDiffModel = (Math.sqrt(stDevDiffModel/(matched1.getN() - 1.0f)));
        
        TransformationPointFit fit = new TransformationPointFit(params, 
            matched1.getN(), avgDiffModel[0], stDevDiffModel,
            Double.MAX_VALUE
        );
        
        return fit;
    }
     
    /**
     * given the scale, rotation and set 1's reference frame centroids,
     * calculate the translation between set1 and set2 assuming that not all
     * points will match.  transXTol and transYTol allow a tolerance when
     * matching the predicted position of a point in set2.
     * 
     * It's expected that the invoker of this method is trying to solve for
     * translation for sets of points like corners in images.  This assumption
     * means that the number of point pair combinations is always far less
     * than the pixel combinations of translations over x and y.
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * 
     * ALSO NOTE: if you know a better solution exists for translation 
     * parameters that matches fewer points, but has a small avg dist from
     * model and smaller standard deviation from the avg dist from model,
     * then transXTol and transYTol should be set to a smaller value and passed
     * to this method.
     * 
     * @param set1 set of points from image 1 to match to image2.
     * @param set2 set of points from image 2 to be matched with image 1
     * @param rotation given in radians with value between 0 and 2*pi, exclusive
     * @param scale
     * @param centroidX1 the x coordinate of the center of image 1 from which
     * set 1 point are from.
     * @param centroidY1 the y coordinate of the center of image 1 from which
     * set 1 point are from.
     * @param setsAreMatched indicates whether set1 and set2 are already matched
     * @return 
     */
    public TransformationPointFit calculateTranslation(PairIntArray set1, 
        PairIntArray set2, double rotation, 
        double scale, int centroidX1, int centroidY1, boolean setsAreMatched) {
        
        if (set1 == null) {
            throw new IllegalArgumentException("set1 cannot be null");
        }
        if (set2 == null) {
            throw new IllegalArgumentException("set2 cannot be null");
        }
        if (set1.getN() < 2) {
            return null;
        }
        if (set2.getN() < 2) {
            return null;
        }
        
        if (setsAreMatched) {
            
            TransformationParameters params = 
                calculateTranslationForMatched(set1, set2, rotation, scale, 
                centroidX1, centroidY1);
            
            return transformMatched(set1, set2, params, centroidX1, centroidY1);
        }
        
        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");
            
            log.severe("scale cannot be smaller than 1");
            
            return null;
        }
        
        //float tolTransX = 4.f * centroidX1 * 0.02f;
        //float tolTransY = 4.f * centroidY1 * 0.02f;
        
        float tolTransX = Float.MAX_VALUE;
        float tolTransY = Float.MAX_VALUE;
        if (tolTransX < 1) {
            tolTransX = 1;
        }
        if (tolTransY < 1) {
            tolTransY = 1;
        }
        
        float s = (float)scale;
        float scaleTimesCosine = (float)(s * Math.cos(rotation));
        float scaleTimesSine = (float)(s * Math.sin(rotation));
        
        int nTrans = set1.getN() * set2.getN();
        int count = 0;
        float[] transX = new float[nTrans];
        float[] transY = new float[nTrans];
        
        // store xr and yr for evaluation of fits
        float[] xr = new float[set1.getN()];
        float[] yr = new float[xr.length];
        
        /*
        Since the points are unmatched, determine the translation for each 
        possible point pair and keep the most frequent answer as the
        estimation for translation in X and in Y.
        */
                
        for (int i = 0; i < set1.getN(); i++) {
            
            int x = set1.getX(i);
            int y = set1.getY(i);
            
            xr[i] = centroidX1*s + ( 
                ((x - centroidX1) * scaleTimesCosine) +
                ((y - centroidY1) * scaleTimesSine));
            
            yr[i] = centroidY1*s + ( 
                (-(x - centroidX1) * scaleTimesSine) +
                ((y - centroidY1) * scaleTimesCosine));
                            
            for (int j = 0; j < set2.getN(); j++) {
                
                int x2 = set2.getX(j);
                int y2 = set2.getY(j);
                
                transX[count] = x2 - xr[i];
                transY[count] = y2 - yr[i];
                
                count++;
            }
        }
        
        float peakTransX;
        float peakTransY;
        
        // when there aren't enough points for useful histogram,
        // will make a frequency map of round to integer,
        // and take the peak if its larger than next peak,
        // else, take the average of largest frequencies.
        
        if ((set1.getN() > 5) && (set2.getN() > 5)) {
            
            int nBins = 15;
            
            HistogramHolder hX = Histogram
                .createSimpleHistogram(nBins,
                //.defaultHistogramCreator(
                transX, Errors.populateYErrorsBySqrt(transX));
        
            HistogramHolder hY = Histogram
                .createSimpleHistogram(nBins,
                //.defaultHistogramCreator(
                transY, Errors.populateYErrorsBySqrt(transY));

            try {
                hX.plotHistogram("transX", 1);
                hY.plotHistogram("transY", 2);
            } catch (IOException e) {
                log.severe(e.getMessage());
            }

            int idxX =  MiscMath.findYMaxIndex(hX.getYHist());
            int idxY =  MiscMath.findYMaxIndex(hY.getYHist());
            peakTransX = hX.getXHist()[idxX];
            peakTransY = hY.getXHist()[idxY];
            
            /*
            float frac = 0.5f;
            
            ArrayPair xy = LinesAndAngles.createPolygonOfTopFWFractionMax(
                hX.getXHist(), hX.getYHistFloat(), null, null, frac);
            
            //area, x, y
            float[] areaAndCentroid = 
                LinesAndAngles.calcAreaAndCentroidOfSimplePolygon(
                    xy.getX(), xy.getY());
            
            peakTransX = areaAndCentroid[1];
            
            xy = LinesAndAngles.createPolygonOfTopFWFractionMax(
                hY.getXHist(), hY.getYHistFloat(), null, null, frac);
            
            areaAndCentroid = 
                LinesAndAngles.calcAreaAndCentroidOfSimplePolygon(
                    xy.getX(), xy.getY());
            
            peakTransY = areaAndCentroid[1];
            */
        } else {
            // keep unique x, y pairs
            List<Integer> freq = new ArrayList<Integer>();
            List<String> xy = new ArrayList<String>();
            for (int i = 0; i < transX.length; i++) {
                int tx = Math.round(transX[i]);
                int ty = Math.round(transY[i]);
                String key = String.format("%d:%d", tx, ty);
                int idx = xy.indexOf(key);
                if (idx > -1) {
                    Integer c = freq.get(idx);
                    freq.set(idx, Integer.valueOf(c.intValue() + 1));
                } else {
                    Integer c = Integer.valueOf(1);
                    freq.add(c);
                    xy.add(key);
                }
            }
            int n0Idx = -1;
            int maxN0 = Integer.MIN_VALUE;
            for (int i = 0; i < freq.size(); i++) {
                int v = freq.get(i).intValue();
                if (v > maxN0) {
                    maxN0 = v;
                    n0Idx = i;
                }
            }
            int n1Idx = -1;
            int maxN1 = Integer.MIN_VALUE;
            for (int i = 0; i < freq.size(); i++) {
                int v = freq.get(i).intValue();
                if (v < maxN0) {
                    if (v > maxN1) {
                        maxN1 = v;
                        n1Idx = i;
                    }
                }
            }
            if ((n1Idx > -1) && ((maxN0 - maxN1) > 1)) {
                //TODO: consider whether (maxN0 - maxN1) should == set1.getN()
                String[] xyItems = xy.get(n0Idx).split(":");
                peakTransX = Float.valueOf(xyItems[0]);
                peakTransY = Float.valueOf(xyItems[1]);
            } else {
                //take average of all points
                peakTransX = 0;
                peakTransY = 0;
                for (int i = 0; i < transX.length; i++) {
                    peakTransX += transX[i];
                    peakTransY += transY[i];
                }
                peakTransX /= (float)transX.length;
                peakTransY /= (float)transY.length;
            }
        }
                
        log.fine("peakTransX=" + peakTransX + "  peakTransY=" + peakTransY);
        
        // grid search or simplex around these 2 peaks to refine the answer?
        
        /*
        to refine the solution from peakTransX and peakTransY
        could either:
        (1) iterate over []transX, []transY and calculate the fit only for the
            values that are within +-20 of peakTransX and peakTransY and keep 
            the best.
            cons: might be lot of redundant calculations
        OR
        (2) create a simplex search for about 5 starter points from -20 to -1
            for translation in X and in Y
            cons: making a fine solution w/ simplex requires more starter
                  points to replace a grid search without using as many
                  tries as a grid search.
        OR
        (3) perform a grid search over a 20pix x 20 pix region surrounding
            peakTransX and peakTransY
            cons: that's 400 fits
        
        ==> choosing (2)
        */
     
        /*
        // for refinement, need a smaller than infinite translation tolerance
        tolTransX = 2.f * centroidX1 * 0.02f;
        tolTransY = 2.f * centroidY1 * 0.02f;
        
        TransformationPointFit bestFit = refineTranslationWithDownhillSimplex(
            xr, yr, peakTransX, peakTransY, tolTransX, tolTransY,
            100.f, 100.f, set2, s, (float)rotation, setsAreMatched
        );
        */
        
        // instead of refinement, evaluate the fit:
         
        
        TransformationPointFit bestFit = evaluateFit(
            xr, yr, peakTransX, peakTransY, tolTransX, tolTransY,
            set2, s, (float)rotation, setsAreMatched);
        
        
        /*
        TransformationPointFit bestFit = evaluateFitForUnMatchedOptimal( 
            xr, yr, peakTransX, peakTransY, tolTransX, tolTransY,
            set2, s, (float)rotation);
        */ 
        
        return bestFit;
    }
    
    /**
     * given the scale, rotation and set 1's reference frame centroids,
     * calculate the translation between set1 and set2 assuming that not all
     * points will match.  transXTol and transYTol allow a tolerance when
     * matching the predicted position of a point in set2.
     * 
     * It's expected that the invoker of this method is trying to solve for
     * translation for sets of points like corners in images.  This assumption
     * means that the number of point pair combinations is always far less
     * than the pixel combinations of translations over x and y.
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * 
     * ALSO NOTE: if you know a better solution exists for translation 
     * parameters that matches fewer points, but has a small avg dist from
     * model and smaller standard deviation from the avg dist from model,
     * then transXTol and transYTol should be set to a smaller value and passed
     * to this method.
     * 
     * @param set1 set of points from image 1 to match to image2.
     * @param set2 set of points from image 2 to be matched with image 1
     * @param rotation given in radians with value between 0 and 2*pi, exclusive
     * @param scale
     * @param centroidX1 the x coordinate of the center of image 1 from which
     * set 1 point are from.
     * @param centroidY1 the y coordinate of the center of image 1 from which
     * set 1 point are from.
     * @return 
     */
    public TransformationParameters calculateTranslationForUnmatched(
        PairIntArray set1, PairIntArray set2, double rotation, double scale, 
        int centroidX1, int centroidY1) {
        
        if (set1 == null) {
            throw new IllegalArgumentException("set1 cannot be null");
        }
        if (set2 == null) {
            throw new IllegalArgumentException("set2 cannot be null");
        }
        if (set1.getN() < 2) {
            return null;
        }
        if (set2.getN() < 2) {
            return null;
        }
        
        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");
            
            log.severe("scale cannot be smaller than 1");
            
            return null;
        }
        
        float s = (float)scale;
        float scaleTimesCosine = (float)(s * Math.cos(rotation));
        float scaleTimesSine = (float)(s * Math.sin(rotation));
        
        int nTrans = set1.getN() * set2.getN();
        int count = 0;
        float[] transX = new float[nTrans];
        float[] transY = new float[nTrans];
        
        // store xr and yr for evaluation of fits
        float[] xr = new float[set1.getN()];
        float[] yr = new float[xr.length];
        
        /*
        Since the points are unmatched, determine the translation for each 
        possible point pair and keep the most frequent answer as the
        estimation for translation in X and in Y.
        */
                
        for (int i = 0; i < set1.getN(); i++) {
            
            int x = set1.getX(i);
            int y = set1.getY(i);
            
            xr[i] = centroidX1*s + ( 
                ((x - centroidX1) * scaleTimesCosine) +
                ((y - centroidY1) * scaleTimesSine));
            
            yr[i] = centroidY1*s + ( 
                (-(x - centroidX1) * scaleTimesSine) +
                ((y - centroidY1) * scaleTimesCosine));
                            
            for (int j = 0; j < set2.getN(); j++) {
                
                int x2 = set2.getX(j);
                int y2 = set2.getY(j);
                
                transX[count] = x2 - xr[i];
                transY[count] = y2 - yr[i];
                
                count++;
            }
        }
        
        float peakTransX;
        float peakTransY;
        
        // when there aren't enough points for useful histogram,
        // will make a frequency map of round to integer,
        // and take the peak if its larger than next peak,
        // else, take the average of largest frequencies.
        
        if ((set1.getN() > 5) && (set2.getN() > 5)) {
            
            int nBins = 15;
            
            HistogramHolder hX = Histogram
                .createSimpleHistogram(nBins,
                //.defaultHistogramCreator(
                transX, Errors.populateYErrorsBySqrt(transX));
        
            HistogramHolder hY = Histogram
                .createSimpleHistogram(nBins,
                //.defaultHistogramCreator(
                transY, Errors.populateYErrorsBySqrt(transY));

            try {
                hX.plotHistogram("transX", 1);
                hY.plotHistogram("transY", 2);
            } catch (IOException e) {
                log.severe(e.getMessage());
            }
            
            int idxX =  MiscMath.findYMaxIndex(hX.getYHist());
            int idxY =  MiscMath.findYMaxIndex(hY.getYHist());
            peakTransX = hX.getXHist()[idxX];
            peakTransY = hY.getXHist()[idxY];
            
            /*
            float frac = 0.5f;

            ArrayPair xy = LinesAndAngles.createPolygonOfTopFWFractionMax(
                hX.getXHist(), hX.getYHistFloat(), null, null, frac);
            
            //area, x, y
            float[] areaAndCentroid = 
                LinesAndAngles.calcAreaAndCentroidOfSimplePolygon(
                    xy.getX(), xy.getY());
            
            peakTransX = areaAndCentroid[1];
            
            xy = LinesAndAngles.createPolygonOfTopFWFractionMax(
                hY.getXHist(), hY.getYHistFloat(), null, null, frac);
            
            areaAndCentroid = 
                LinesAndAngles.calcAreaAndCentroidOfSimplePolygon(
                    xy.getX(), xy.getY());
            
            peakTransY = areaAndCentroid[1];
            */
            
        } else {
            // keep unique x, y pairs
            List<Integer> freq = new ArrayList<Integer>();
            List<String> xy = new ArrayList<String>();
            for (int i = 0; i < transX.length; i++) {
                int tx = Math.round(transX[i]);
                int ty = Math.round(transY[i]);
                String key = String.format("%d:%d", tx, ty);
                int idx = xy.indexOf(key);
                if (idx > -1) {
                    Integer c = freq.get(idx);
                    freq.set(idx, Integer.valueOf(c.intValue() + 1));
                } else {
                    Integer c = Integer.valueOf(1);
                    freq.add(c);
                    xy.add(key);
                }
            }
            int n0Idx = -1;
            int maxN0 = Integer.MIN_VALUE;
            for (int i = 0; i < freq.size(); i++) {
                int v = freq.get(i).intValue();
                if (v > maxN0) {
                    maxN0 = v;
                    n0Idx = i;
                }
            }
            int n1Idx = -1;
            int maxN1 = Integer.MIN_VALUE;
            for (int i = 0; i < freq.size(); i++) {
                int v = freq.get(i).intValue();
                if (v < maxN0) {
                    if (v > maxN1) {
                        maxN1 = v;
                        n1Idx = i;
                    }
                }
            }
            if ((n1Idx > -1) && ((maxN0 - maxN1) > 1)) {
                //TODO: consider whether (maxN0 - maxN1) should == set1.getN()
                String[] xyItems = xy.get(n0Idx).split(":");
                peakTransX = Float.valueOf(xyItems[0]);
                peakTransY = Float.valueOf(xyItems[1]);
            } else {
                //take average of all points
                peakTransX = 0;
                peakTransY = 0;
                for (int i = 0; i < transX.length; i++) {
                    peakTransX += transX[i];
                    peakTransY += transY[i];
                }
                peakTransX /= (float)transX.length;
                peakTransY /= (float)transY.length;
            }
        }
                
        log.fine("peakTransX=" + peakTransX + "  peakTransY=" + peakTransY);
      
        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians((float)rotation);
        params.setScale(s);
        params.setTranslationX(peakTransX);
        params.setTranslationY(peakTransY);
        
        return params;
    }
    
    /**
     * given the scale, rotation and set 1's reference frame centroids,
     * calculate the translation between set1 and set2 assuming that not all
     * points will match.  transXTol and transYTol allow a tolerance when
     * matching the predicted position of a point in set2.
     * 
     * It's expected that the invoker of this method is trying to solve for
     * translation for sets of points like corners in images.  This assumption
     * means that the number of point pair combinations is always far less
     * than the pixel combinations of translations over x and y.
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * 
     * ALSO NOTE: if you know a better solution exists for translation 
     * parameters that matches fewer points, but has a small avg dist from
     * model and smaller standard deviation from the avg dist from model,
     * then transXTol and transYTol should be set to a smaller value and passed
     * to this method.
     * 
     * @param matched1 set of points from image 1 to match to image2.
     * @param matched2 set of points from image 2 to be matched with image 1
     * @param rotation given in radians with value between 0 and 2*pi, exclusive
     * @param scale
     * @param centroidX1 the x coordinate of the center of image 1 from which
     * set 1 point are from.
     * @param centroidY1 the y coordinate of the center of image 1 from which
     * set 1 point are from.
     * @return 
     */
    public TransformationParameters calculateTranslationForMatched(
        PairIntArray matched1, 
        PairIntArray matched2, double rotation, 
        double scale, int centroidX1, int centroidY1) {
        
        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");
            
            log.severe("scale cannot be smaller than 1");
            
            return null;
        }
                        
        double scaleTimesCosine = scale * Math.cos(rotation);
        double scaleTimesSine = scale * Math.sin(rotation);
        
        double avgTransX = 0;
        double avgTransY = 0;
        
        for (int i = 0; i < matched1.getN(); i++) {
            
            int x = matched1.getX(i);
            int y = matched1.getY(i);
            
            double xr = centroidX1*scale + ( 
                ((x - centroidX1) * scaleTimesCosine) +
                ((y - centroidY1) * scaleTimesSine));
            
            double yr = centroidY1*scale + ( 
                (-(x - centroidX1) * scaleTimesSine) +
                ((y - centroidY1) * scaleTimesCosine));
                            
            int x2 = matched2.getX(i);
                
            int y2 = matched2.getY(i);
                
            avgTransX += (int)Math.round(x2 - xr);
                
            avgTransY += (int)Math.round(y2 - yr);
        }
        
        avgTransX /= (double)matched1.getN();
        
        avgTransY /= (double)matched1.getN();
                                
        TransformationParameters params = 
            new TransformationParameters();
        params.setRotationInRadians((float) rotation);
        params.setScale((float) scale);
        params.setTranslationX((float)avgTransX);
        params.setTranslationY((float)avgTransY);

        return params;
    }
    
    public TransformationPointFit calculateTranslation(PairIntArray set1, 
        PairIntArray set2, 
        double prevNearTransX, double prevNearTransY,
        double transXTol, double transYTol, double rotation, 
        double scale, int centroidX1, int centroidY1) {
        
        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");
            
            log.severe("scale cannot be smaller than 1");
            
            return null;
        }
        
        int nMax = set1.getN() * set2.getN();
        
        int maxFound = Integer.MIN_VALUE;
        TransformationPointFit bestFit = null;
        
        double scaleTimesCosine = scale * Math.cos(rotation);
        double scaleTimesSine = scale * Math.sin(rotation);

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        double[] xyCentroids1 = curveHelper.calculateXYCentroids(set1);
        double[] xyCentroids2 = curveHelper.calculateXYCentroids(set2);
        
        double xr = centroidX1*scale + ( 
            ((xyCentroids1[0] - centroidX1) * scaleTimesCosine)
            + ((xyCentroids1[1] - centroidY1) * scaleTimesSine));
        
        double yr = centroidY1 * scale 
            + ((-(xyCentroids1[0] - centroidX1) * scaleTimesSine)
            + ((xyCentroids1[1] - centroidY1) * scaleTimesCosine));
        
        int transXC = (int)Math.round(xyCentroids2[0] - xr);
        int transYC = (int)Math.round(xyCentroids2[1] - yr);
        
        //TODO: need tests for this.  the assumed change is small
        double deltaTransX = prevNearTransX - transXC;
        if (deltaTransX < 0) {
            deltaTransX *= -1;
        }
        if (deltaTransX < 1) {
            deltaTransX = 7;
        }
        double deltaTransY = prevNearTransY - transYC;
        if (deltaTransY < 0) {
            deltaTransY *= -1;
        }
        if (deltaTransY < 1) {
            deltaTransY = 7;
        }

        int tx0 = (int)Math.round
            (((prevNearTransX + transXC)/2.) - 1.2*deltaTransX);
        int ty0 = (int)Math.round
            (((prevNearTransY + transYC)/2.) - 1.2*deltaTransY);

        int tx1 = (int)Math.round
            (((prevNearTransX + transXC)/2.) + 1.2*deltaTransX);
        int ty1 = (int)Math.round
            (((prevNearTransY + transYC)/2.) + 1.2*deltaTransY);

        if (ty0 > -2) {
            ty0 = -2;
        }
        if (ty1 < 2) {
            ty1 = 2;
        }
        
        for (int transX = tx0; transX < tx1; transX++) {
            for (int transY = ty0; transY < ty1; transY++) {
                
                TransformationParameters params = 
                    new TransformationParameters();
                params.setRotationInRadians((float) rotation);
                params.setScale((float) scale);
                params.setTranslationX(transX);
                params.setTranslationY(transY);

                TransformationPointFit fit = transform(set1, set2, params, 
                    transXTol, transYTol, centroidX1, centroidY1);

                if ((fit != null) && fitIsBetter(bestFit, fit)) {
                    
                    bestFit = fit;
                
                    if ((fit.getMeanDistFromModel() == 0) &&
                        (fit.getNumberOfMatchedPoints() > (0.3*set1.getN())) 
                        ) {
                        break;
                    }
                }
            }
        }
        
        return (bestFit != null) ? bestFit : null;
    }
     
    /**
     * refine the transformation params for set1 to better match set2.
     * @param params
     * @param set1
     * @param set2
     * @param image1Width
     * @param image1Height
     * @param setsAreMatched
     * @return 
     */
    public TransformationPointFit refineTransformationWithDownhillSimplex(
        TransformationParameters params,
        PairIntArray set1, PairIntArray set2, 
        int image1Width, int image1Height, boolean setsAreMatched) {
        
        int centroidX1 = image1Width >> 1;
        
        int centroidY1 = image1Height >> 1;
        
        // projection effects from different camera positions or orientations
        // are not calculated in this point matcher, 
        // but one can set a
        // tolerance in translation to try to allow for a small amount of it.
        // 
        // as an example: rotation of 13 degrees at top of image, 
        // and 0 at bottom from turning the camera slightly during panoramic
        // photos can result in delta translation of 0.02 per pix in X and 
        // 0.016 per pix in Y
        
        //TODO: consider allowing these to be larger
        double transXTol = image1Width * 0.02;
        
        double transYTol = image1Height * 0.02;
        
        double r = params.getRotationInRadians();
        double s = params.getScale();

        double[] drs = new double[] {            
            -5.0 * Math.PI/180., 
            -2.5 * Math.PI/180.,
            -1.0 * Math.PI/180.,
            1.0 * Math.PI/180.,
            2.5 * Math.PI/180.,
            5.0 * Math.PI/180.
        };
        if (r == 0) {
             drs = new double[]{0};
        }

        double rMin = r - (10 * Math.PI/180);
        double rMax = r + (10 * Math.PI/180);
        double sMin = s - 0.5;
        double sMax = s + 0.5;
        
        double[] dss = new double[] {
            0.1
            //-1.0, -0.1, -0.05, 0.05
        };
        if (s == 1) {
            dss = new double[]{0};
            sMin = 1;
        }
        if (rMin < 0) {
            rMin = 0;
        }
        
        int n = (1 + dss.length) * (1 + drs.length);
        
        // start with simplex for at least 3 points (fitting 2 parameters)
        TransformationPointFit[] fits = new TransformationPointFit[n];
        
        int count = 0;
        for (int i = 0; i <= dss.length; i++) {
            
            double scale = (i == 0) ? s : s + dss[i - 1];
            
            for (int j = 0; j <= drs.length; j++) {
                
                double rotation = (j == 0) ? r : r + drs[j - 1];
                
                fits[count] = calculateTranslation(set1, set2, 
                    rotation, scale, centroidX1, centroidY1, setsAreMatched);
                
                if (fits[count] != null) {
                    count++;
                }
            }
        }
        
        if (count < n) {
            fits = Arrays.copyOf(fits, count);
        }
        
        float reflectionCoeff = 1;   // > 0
        float expansionCoeff = 2;   // > 1
        float contractionCoeff = -0.5f;
        float reductionCoeff = 0.5f;

        TransformationPointFit bestFit = fitWithDownhillSimplex(
            set1, set2, 
            fits, centroidX1, centroidY1, transXTol, transYTol,
            r, s, rMin, rMax, sMin, sMax,
            reflectionCoeff, expansionCoeff, contractionCoeff, reductionCoeff,
            setsAreMatched
        );
        
        if ((bestFit != null) && 
            (bestFit.getParameters().getRotationInRadians() > 2.*Math.PI)) {
            float rot = bestFit.getParameters().getRotationInRadians();
            while (rot >= 2*Math.PI) {
                rot -= 2*Math.PI;
            }
            bestFit.getParameters().setRotationInRadians(rot);
        }
        
        return bestFit;
    }
  
    /**
     * 
     * @param set1
     * @param set2
     * @param fits
     * @param centroidX1
     * @param centroidY1
     * @param transXTol
     * @param transYTol
     * @param r
     * @param s
     * @param rMin
     * @param rMax
     * @param sMin
     * @param sMax
     * @param alpha parameter used to adjust size of changes for "reflection"
     * @param gamma parameter used to adjust size of changes for "expansion"
     * @param beta parameter used to adjust size of changes for "contraction"
     * @param tau parameter used to adjust size of changes for "reduction"
     * 
     * @return 
     */
    private TransformationPointFit fitWithDownhillSimplex(
        PairIntArray set1, PairIntArray set2, 
        TransformationPointFit[] fits, int centroidX1, int centroidY1,
        double transXTol, double transYTol,
        double r, double s,
        double rMin, double rMax, double sMin, double sMax,
        float alpha, float gamma, float beta, float tau, boolean setsAreMatched
        ) {
        
        if (alpha < 0) {
            throw new IllegalArgumentException("alpha must be > 0");
        }
        if (gamma < 1) {
            throw new IllegalArgumentException("gamma must be > 1");
        }
        
        int convergence = (set1.getN() < set2.getN()) ? set1.getN() 
            : set2.getN();
      
        if (s == 1) {
            sMin = 1;
        }
        if (rMin < 0) {
            rMin = 0;
        }

        boolean go = true;

        int nMaxIter = 100;
        int nIter = 0;
        
        int bestFitIdx = 0;
        int secondWorstFitIdx = fits.length - 2;
        int worstFitIdx = fits.length - 1;

        int lastNMatches = Integer.MIN_VALUE;
        double lastAvgDistModel = Double.MAX_VALUE;
        int nIterSameMin = 0;
        
        while (go && (nIter < nMaxIter)) {

            // this has changed to allow smaller avg +- stdv to be better
            // than more matches when close
            sortByDescendingMatches(fits, 0, (fits.length - 1));
            
            if ((lastNMatches == fits[0].getNumberOfMatchedPoints()) &&
                (Math.abs(lastAvgDistModel - fits[0].getMeanDistFromModel())
                < 0.01)) {
                nIterSameMin++;
                if (nIterSameMin >= 5) {
                    break;
                }
            } else {
                nIterSameMin = 0;
            }
            lastNMatches = fits[0].getNumberOfMatchedPoints();
            lastAvgDistModel = fits[0].getMeanDistFromModel();
            

            // determine center for all points excepting the worse fit
            double rSum = 0.0;
            double sSum = 0.0;
            for (int i = 0; i < (fits.length - 1); i++) {
                rSum += fits[i].getRotationInRadians();
                sSum += fits[i].getScale();
            }
            r = rSum / (double)(fits.length - 1);
            s = sSum / (double)(fits.length - 1);

            // "Reflection"
            double rReflect = r + (alpha * 
                (r - fits[worstFitIdx].getRotationInRadians()));
            double sReflect = s +
                (s - alpha * (fits[worstFitIdx].getScale()));
         
            if (sReflect < 1) {
                sReflect = 1;
            }

            TransformationPointFit fitReflected = 
                calculateTranslation(set1, set2, 
                    rReflect, sReflect,
                    centroidX1, centroidY1, setsAreMatched);
            
            boolean relectIsWithinBounds = 
                (rReflect >= rMin) && (rReflect <= rMax) 
                && (sReflect >= sMin) && (sReflect <= sMax);

            if (fitIsBetter(fits[secondWorstFitIdx], fitReflected)
                && relectIsWithinBounds
                && !fitIsBetter(fits[bestFitIdx], fitReflected) ) {
                
                fits[worstFitIdx] = fitReflected;
                
            } else {
                
                if (fitIsBetter(fits[bestFitIdx], fitReflected)
                    && relectIsWithinBounds) {
                    
                    // "Expansion"
                    double rExpansion = r + (gamma * 
                        (r - fits[worstFitIdx].getRotationInRadians()));
                    double sExpansion = s + (gamma * 
                        (s - fits[worstFitIdx].getScale()));
                    
                    if (sExpansion < 1) {
                        sExpansion = 1;
                    }
                    
                    TransformationPointFit fitExpansion = 
                        calculateTranslation(set1, set2, 
                            rExpansion, sExpansion,
                            centroidX1, centroidY1, setsAreMatched);
                    
                    if (fitIsBetter(fitReflected, fitExpansion)
                        && ((rExpansion >= rMin) && (rExpansion <= rMax)
                        && (sExpansion >= sMin) && (sExpansion <= sMax))) {

                        fits[worstFitIdx] = fitExpansion;
                        
                    } else {
                        
                        fits[worstFitIdx] = fitReflected;
                    }
                    
                } else {
                
                    // we know that the reflection fit is worse than the 2nd worse

                    // "Contraction"
                    double rContraction = r + (beta * 
                        (r - fits[worstFitIdx].getRotationInRadians()));
                    double sContraction = s + (beta * 
                        (s - fits[worstFitIdx].getScale()));
                    
                    TransformationPointFit fitContraction = 
                        calculateTranslation(set1, set2, 
                            rContraction, sContraction,
                            centroidX1, centroidY1, setsAreMatched);
                
                    if (sContraction < 1) {
                        sContraction = 1;
                    }

                    if (fitIsBetter(fits[worstFitIdx], fitContraction)
                        && 
                        (rContraction >= rMin) && (rContraction <= rMax) 
                        && (sContraction >= sMin) && (sContraction <= sMax)
                    ) {

                        fits[worstFitIdx] = fitContraction;
                        
                    } else {
                                                
                        // "Reduction"
                        for (int i = 1; i < fits.length; i++) {
                            
                            double rReduction = 
                                fits[bestFitIdx].getRotationInRadians()
                                + (tau * 
                                (fits[i].getRotationInRadians()
                                - fits[bestFitIdx].getRotationInRadians()));
                            double sReduction = 
                                fits[bestFitIdx].getScale()
                                + (tau * 
                                (fits[i].getScale()
                                - fits[bestFitIdx].getScale()));
                            
                            if (sReduction < 1) {
                                sReduction = 1;
                            }
                            
                            //NOTE: there's a possibility of a null fit.
                            //  instead of re-writing the fits array, will 
                            //  assign a fake infinitely bad fit which will 
                            //  fall to the bottom of the list after the next 
                            //  sort.
                            TransformationPointFit fit = calculateTranslation(
                                set1, set2, 
                                rReduction, sReduction, centroidX1, centroidY1,
                                setsAreMatched);
                            
                            if (fit != null) {
                                fits[i] = fit;
                            } else {
                                fits[i] = new TransformationPointFit(
                                    new TransformationParameters(),
                                    0, Double.MAX_VALUE, Double.MAX_VALUE,
                                    Double.MAX_VALUE);
                            }
                        }
                    }
                }
            }

            log.finest("best fit so far: " + 
                fits[bestFitIdx].getNumberOfMatchedPoints());
            
            nIter++;

            if ((fits[bestFitIdx].getNumberOfMatchedPoints() == convergence) 
                && (fits[bestFitIdx].getMeanDistFromModel() == 0)) {
                go = false;
            } else if ((r > rMax) || (r < rMin)) {
                go = false;
            } else if ((s > sMax) || (s < sMin)) {
                go = false;
            }
        }
        
        // if rotation > 2PI, subtract 2PI
        if ((fits[bestFitIdx] != null) && 
            (fits[bestFitIdx].getParameters().getRotationInRadians() 
            > 2.*Math.PI)) {
            float rot = fits[bestFitIdx].getParameters().getRotationInRadians();
            while (rot >= 2*Math.PI) {
                rot -= 2*Math.PI;
            }
            fits[bestFitIdx].getParameters().setRotationInRadians(rot);
        }
        
        return fits[bestFitIdx];
    }
   
    /**
     * refine the transformation params to make a better match of edges1 to
     * edges2 where the points within edges in both sets are not necessarily
     * 1 to 1 matches (that is, the input is not expected to be matched 
     * already).
     * 
     * TODO: improve transformEdges to find translation for all edges
     * via a search method rather than trying all pairs of points.
     * 
     * @param edges1
     * @param edges2
     * @param params
     * @param centroidX1
     * @param centroidY1
     * @param centroidX2
     * @param centroidY2
     * @return 
     */
    public TransformationParameters refineTransformation(PairIntArray[] edges1, 
        PairIntArray[] edges2, final TransformationParameters params,
        final int centroidX1, final int centroidY1,
        final int centroidX2, final int centroidY2) {
                    
        //TODO: set this empirically from tests
        double convergence = 0;
                
        double r = params.getRotationInRadians();
        double s = params.getScale();
        
        double rMin = r - (10 * Math.PI/180);
        double rMax = r + (10 * Math.PI/180);
        double sMin = s - 1.5;
        double sMax = s + 1.5;
        
        // the positive offsets can be found w/ reflection?
        // TODO: needs testing for starter points.  these are supplying the
        // "grid search" portion of exploring more than local space
        double[] drs = new double[] {
            -5.0 * Math.PI/180., 
            -2.5 * Math.PI/180.,
            -1.0 * Math.PI/180.,
            1.0 * Math.PI/180.,
            2.5 * Math.PI/180.,
            5.0 * Math.PI/180.
        };
        double[] dss = new double[] {
            -1.0, -0.1, -0.05, 0.05 /*, 0.05, 0.1, 1.0*/
        };
        if (r == 0) {
             drs = new double[]{0};
        }
        if (s == 1) {
            dss = new double[]{0};
            sMin = 1;
        }
        if (sMin < 1) {
            sMin = 1;
        }
        if (rMin < 0) {
            rMin = 0;
        }
        
        int n = (1 + dss.length) * (1 + drs.length);
        
        TransformationPointFit[] fits = new TransformationPointFit[n];
        
        int count = 0;
        for (int i = 0; i <= dss.length; i++) {
            
            double scale = (i == 0) ? s : s + dss[i - 1];
                        
            for (int j = 0; j <= drs.length; j++) {
                
                double rotation = (j == 0) ? r : r + drs[j - 1];
                                
                fits[count] = calculateTranslationAndTransformEdges(
                    rotation, scale, edges1, edges2, 
                    params.getTranslationX(), params.getTranslationY(),
                    centroidX1, centroidY1, centroidX2, centroidY2);
                
                if (fits[count] != null) {
                    count++;
                }
            }
        }
        
        if (count < n) {
            fits = Arrays.copyOf(fits, count);
        }
        
        float alpha = 1;   // > 0
        float gamma = 2;   // > 1
        float beta = -0.5f; 
        float tau = 0.5f;
            
        boolean go = true;

        int nMaxIter = 100;
        int nIter = 0;
        
        int bestFitIdx = 0;
        int secondWorstFitIdx = fits.length - 2;
        int worstFitIdx = fits.length - 1;

        int lastNMatches = Integer.MIN_VALUE;
        double lastAvgDistModel = Double.MAX_VALUE;
        int nIterSameMin = 0;
        
        while (go && (nIter < nMaxIter)) {

            if (fits.length == 0) {
                break;
            }
            
            sortByDescendingMatches(fits, 0, (fits.length - 1));

/*if (fits.length > 0) {
    log.info("best fit:  n=" + fits[bestFitIdx].getNumberOfMatchedPoints()
    + " dm=" + fits[bestFitIdx].getMeanDistFromModel()
    + " params:\n" + fits[bestFitIdx].getParameters().toString());
}*/          

            if ((lastNMatches == fits[bestFitIdx].getNumberOfMatchedPoints()) &&
                (Math.abs(lastAvgDistModel - 
                fits[bestFitIdx].getMeanDistFromModel()) < 0.01)) {
                
                nIterSameMin++;
                if (nIterSameMin >= 5) {
                    break;
                }
                
            } else {
                nIterSameMin = 0;
            }
            
            lastNMatches = fits[bestFitIdx].getNumberOfMatchedPoints();
            lastAvgDistModel = fits[bestFitIdx].getMeanDistFromModel();
            
            // determine center for all points excepting the worse fit
            double rSum = 0.0;
            double sSum = 0.0;
            for (int i = 0; i < (fits.length - 1); i++) {
                rSum += fits[i].getRotationInRadians();
                sSum += fits[i].getScale();
            }
            r = rSum / (double)(fits.length - 1);
            s = sSum / (double)(fits.length - 1);

            // "Reflection"
            double rReflect = r + (alpha * 
                (r - fits[worstFitIdx].getRotationInRadians()));
            double sReflect = s +
                (s - alpha * (fits[worstFitIdx].getScale()));
            
            TransformationPointFit fitReflected = 
                calculateTranslationAndTransformEdges(
                rReflect, sReflect, 
                edges1, edges2, 
                params.getTranslationX(), params.getTranslationY(),
                centroidX1, centroidY1, centroidX2, centroidY2);
            
            boolean relectIsWithinBounds = 
                /*(rReflect >= rMin) && (rReflect <= rMax) 
                &&*/ (sReflect >= sMin) && (sReflect <= sMax);

            if (fitIsBetter(fits[secondWorstFitIdx], fitReflected)
                && relectIsWithinBounds
                && !fitIsBetter(fits[bestFitIdx], fitReflected) ) {
                                
                fits[worstFitIdx] = fitReflected;
                
            } else {
                
                if (fitIsBetter(fits[bestFitIdx], fitReflected)
                    && relectIsWithinBounds) {
                    
                    // "Expansion"
                    double rExpansion = r + (gamma * 
                        (r - fits[worstFitIdx].getRotationInRadians()));
                    double sExpansion = s + (gamma * 
                        (s - fits[worstFitIdx].getScale()));
                    
                    TransformationPointFit fitExpansion = 
                        calculateTranslationAndTransformEdges(
                        rExpansion, sExpansion, edges1, edges2,
                        params.getTranslationX(), params.getTranslationY(),
                        centroidX1, centroidY1, centroidX2, centroidY2);
                    
                    if (fitIsBetter(fitReflected, fitExpansion)
                        && (/*(rExpansion >= rMin) && (rExpansion <= rMax)
                        &&*/ (sExpansion >= sMin) && (sExpansion <= sMax))) {

                        fits[worstFitIdx] = fitExpansion;
                        
                    } else {
                        
                        fits[worstFitIdx] = fitReflected;
                    }
                    
                } else {
                
                    // the reflection fit is worse than the 2nd worse

                    // "Contraction"
                    double rContraction = r + (beta * 
                        (r - fits[worstFitIdx].getRotationInRadians()));
                    double sContraction = s + (beta * 
                        (s - fits[worstFitIdx].getScale()));
                    
                    TransformationPointFit fitContraction = 
                        calculateTranslationAndTransformEdges(
                        rContraction, sContraction, edges1, edges2,
                        params.getTranslationX(), params.getTranslationY(),
                        centroidX1, centroidY1, centroidX2, centroidY2);
                    
                    if (fitIsBetter(fits[worstFitIdx], fitContraction)
                        /*&& 
                        (rContraction >= rMin) && (rContraction <= rMax)*/
                        && (sContraction >= sMin) && (sContraction <= sMax)
                    ) {

                        fits[worstFitIdx] = fitContraction;
                        
                    } else {
                        
                        // "Reduction"
                        for (int i = 1; i < fits.length; i++) {
                            
                            double rReduction = 
                                fits[bestFitIdx].getRotationInRadians()
                                + (tau * 
                                (fits[i].getRotationInRadians()
                                - fits[bestFitIdx].getRotationInRadians()));
                            
                            double sReduction = 
                                fits[bestFitIdx].getScale()
                                + (tau * 
                                (fits[i].getScale()
                                - fits[bestFitIdx].getScale()));
                                                        
                            //NOTE: there's a possibility of a null fit.
                            //  instead of re-writing the fits array, will 
                            //  assign a fake infinitely bad fit which will 
                            //  fall to the bottom of the list after the next 
                            //  sort.
                            TransformationPointFit fit = 
                                calculateTranslationAndTransformEdges(
                                rReduction, sReduction, 
                                edges1, edges2, 
                                params.getTranslationX(), 
                                params.getTranslationY(),
                                centroidX1, centroidY1,
                                centroidX2, centroidY2);
                            
                            if (fit != null) {
                                fits[i] = fit;
                            } else {
                                fits[i] = new TransformationPointFit(
                                    new TransformationParameters(),
                                    0, Double.MAX_VALUE, Double.MAX_VALUE,
                                    Double.MAX_VALUE);
                            }
                        }                        
                    }
                }
            }

            log.finest("best fit so far: nMatches=" 
                + fits[bestFitIdx].getNumberOfMatchedPoints() + 
                " diff from model=" + fits[bestFitIdx].getMeanDistFromModel()
                );
            
            nIter++;

            if ((fits[bestFitIdx].getNumberOfMatchedPoints() == convergence) 
                && (fits[bestFitIdx].getMeanDistFromModel() == 0)) {
                go = false;
            /*} else if ((r > rMax) || (r < rMin)) {
                go = false;*/
            } else if ((s > sMax) || (s < sMin)) {
                go = false;
            }
        }
        
        // if rotation > 2PI, subtract 2PI
        if ((fits[bestFitIdx] != null) && 
            (fits[bestFitIdx].getParameters().getRotationInRadians() 
            > 2.*Math.PI)) {
            
            float rot = fits[bestFitIdx].getParameters().getRotationInRadians();
            while (rot >= 2*Math.PI) {
                rot -= 2*Math.PI;
            }
            fits[bestFitIdx].getParameters().setRotationInRadians(rot);
        }
        
        return fits[bestFitIdx].getParameters();
    }
    
    /**
     * TODO: improve transformEdges to find translation for all edges
     * via a search method rather than trying all pairs of points.
     * 
     * Given edges1 and edges2 which we already know are matched edges due to
     * contour matching or other means, and given the rotation and scale,
     * determine the translation between the edges and return the fit.
     * 
     * @param rotInRad
     * @param scl
     * @param edges1
     * @param edges2
     * @param centroidX1
     * @param centroidY1
     * @return 
     */
    private TransformationPointFit calculateTranslationAndTransformEdges(
        double rotInRad, double scl, 
        PairIntArray[] edges1, PairIntArray[] edges2, 
        double prevNearTransX, double prevNearTransY,
        int centroidX1, int centroidY1, int centroidX2, int centroidY2) {
        
        if (edges1 == null) {
            throw new IllegalArgumentException("edges1 cannot be null");
        }
        if (edges2 == null) {
            throw new IllegalArgumentException("edges2 cannot be null");
        }
        if ((edges1.length != edges2.length)) {
            throw new IllegalArgumentException(
            "edges1 and edges2 must be the same length");
        }
    
        /*
        edges1 and edges2 are matched edges, but should be considered clouds
        of points rather than point to point matches within the edges.
        
        For each paired edge in edges1 and edges2, determine the implied
        translation by their centroids.
        These are combined in weighted average where the
        weight is the length of the edge.
        
        Then a small area surrounding the averaged translation is searched
        to find the best fit.
        */
        
        float s = (float)scl;
        float scaleTimesCosine = (float)(s * Math.cos(rotInRad));
        float scaleTimesSine = (float)(s * Math.sin(rotInRad));
        
        //TODO: revisit this:
        float tolTransX = 2.f * centroidX1 * 0.02f;
        float tolTransY = 2.f * centroidY1 * 0.02f;
       
        int nTotal = 0;
        float[] weights = new float[edges1.length];
        for (int i = 0; i < edges1.length; i++) {
            weights[i] = edges1[i].getN();
            nTotal += edges1[i].getN();
        }
        for (int i = 0; i < edges1.length; i++) {
            weights[i] /= (float)nTotal;
        }
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        float translationX = 0;
        float translationY = 0;
        for (int i = 0; i < edges1.length; i++) {
            
            PairIntArray edge1 = edges1[i];
            PairIntArray edge2 = edges2[i];
            
            double[] xycen1 = curveHelper.calculateXYCentroids(edge1);
            
            double[] xycen2 = curveHelper.calculateXYCentroids(edge2);
            
            double srX1 = centroidX1*s + ( 
                ((xycen1[0] - centroidX1) * scaleTimesCosine) +
                ((xycen1[1] - centroidY1) * scaleTimesSine));
            
            double srY1 = centroidY1*s + ( 
                (-(xycen1[0] - centroidX1) * scaleTimesSine) +
                ((xycen1[1] - centroidY1) * scaleTimesCosine));
            
            double tx = xycen2[0] - srX1;
            
            double ty = xycen2[1] - srY1;
            
            translationX += weights[i] * tx;
            
            translationY += weights[i] * ty;
        }
                 
        TransformationPointFit bestFit = refineTranslationWithDownhillSimplex(
            edges1, edges2, translationX, translationY, tolTransX, tolTransY,
            20.f, 20.f, s, (float)rotInRad, centroidX1, centroidY1);
        
        return bestFit;
    }
    
    /**
     * apply the parameters to set1 and find the matches to points in set2
     * within the given tolerance for translations.
     * 
     * runtime complexity is O(n_set1 * n_set2)
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * 
     * ALSO NOTE: if you know a better solution exists for translation 
     * parameters that matches fewer points, but has a small avg dist from
     * model and smaller standard deviation from the avg dist from model,
     * then transXTol and transYTol should be smaller.
     * @param set1
     * @param set2
     * @param params
     * @param transXTol
     * @param transYTol
     * @param centroidX1
     * @param centroidY1
     * @param diffFromModel output variable holding difference of matched
     *    points from model
     * @param diffFromModelOffset input variable holding the offset in 
     * array diffFromModel for which points are added from this method.
     * @param avgDiffModel output variable holding the average difference from
     * the model.
     * @return number of matched points.  returns -1 if the scale is less than 1
     */
    int populateDiffFromModel(PairIntArray set1, PairIntArray set2, 
        TransformationParameters params, double transXTol, double transYTol, 
        int centroidX1, int centroidY1,
        double[] diffFromModel, final int diffFromModelOffset,
        double[] avgDiffModel) {
                
        if (diffFromModel == null) {
            throw new IllegalArgumentException("diffFromModel cannot be null");
        }
        if (diffFromModelOffset < -1) {
            throw new IllegalArgumentException(
            "diffFromModelOffset cannot be < 0");
        }
        if ((diffFromModelOffset + set1.getN()) > diffFromModel.length) {
            throw new IllegalArgumentException(
            "diffFromModelOffset + potentially all matched points is larger"
            + " than array size.");
        }
            
        double scale = params.getScale();
        double rotation = params.getRotationInRadians();
        
        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");
            
            log.severe("scale cannot be smaller than 1");
            
            return -1;
        }
        
        int nMatched = 0;
        double scaleTimesCosine = scale * Math.cos(rotation);
        double scaleTimesSine = scale * Math.sin(rotation);
        
        avgDiffModel[0] = 0;
        
        for (int i = 0; i < set1.getN(); i++) {
        
            int x = set1.getX(i);
            int y = set1.getY(i);
            
            double xmcx1 = x - centroidX1;
            double ymcy1 = y - centroidY1;
           
            double xr = centroidX1 * scale
                + ((xmcx1 * scaleTimesCosine) + (ymcy1 * scaleTimesSine));

            double yr = centroidY1 * scale
                + ((-xmcx1 * scaleTimesSine) + (ymcy1 * scaleTimesCosine));

            int xt = (int)Math.round(xr + params.getTranslationX());
            int yt = (int)Math.round(yr + params.getTranslationY());
            
            int lowerX = xt - (int)transXTol;
            int higherX = xt + (int)transXTol;
            int lowerY = yt - (int)transYTol;
            int higherY = yt + (int)transYTol;
            
            for (int j = 0; j < set2.getN(); j++) {
                
                int x2 = set2.getX(j);
                int y2 = set2.getY(j);
                
                double d = Math.sqrt(Math.pow(xt - x2, 2) + Math.pow(yt - y2, 2)
                );
                
                if ((x2 < lowerX) || (x2 > higherX) || (y2 < lowerY) 
                    || (y2 > higherY)) {
                    
                    continue;
                }
                
                diffFromModel[diffFromModelOffset + nMatched] = d;
                
                avgDiffModel[0] += d;
                
                nMatched++;
                
                break;
            }
        }
        
        avgDiffModel[0] /= (double)nMatched;
        
        return nMatched;
    }
    
    /**
     * apply the parameters to matched1 and return the residuals between
     * items in matched2 and the model from matched1.
     * 
     * runtime complexity is O(N)
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * 
     * @param matched1
     * @param matched2
     * @param params
     * @param centroidX1
     * @param centroidY1
     * @param diffFromModel output variable holding difference of matched
     *    points from model
     * @param avgDiffModel output variable holding the average difference from
     * the model.
     */
    public void populateDiffFromModelWithMatchedPoints(PairIntArray matched1, 
        PairIntArray matched2, 
        TransformationParameters params, 
        int centroidX1, int centroidY1,
        double[] diffFromModel, double[] avgDiffModel) {
                
        if (diffFromModel == null) {
            throw new IllegalArgumentException("diffFromModel cannot be null");
        }
        if (avgDiffModel == null) {
            throw new IllegalArgumentException("avgDiffModel cannot be null");
        }
        if (matched1 == null) {
            throw new IllegalArgumentException("matched1 cannot be null");
        }
        if (matched2 == null) {
            throw new IllegalArgumentException("matched2 cannot be null");
        }
        if (matched1.getN() != matched2.getN()) {
            throw new IllegalArgumentException(
            "matched1 and matched2 should be the same size");
        }
            
        double scale = params.getScale();
        double rotation = params.getRotationInRadians();
        
        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");
            
            log.severe("scale cannot be smaller than 1");
            
            return;
        }
        
        double scaleTimesCosine = scale * Math.cos(rotation);
        double scaleTimesSine = scale * Math.sin(rotation);
        
        avgDiffModel[0] = 0;
        
        for (int i = 0; i < matched1.getN(); i++) {
        
            int x = matched1.getX(i);
            int y = matched1.getY(i);
            
            double xmcx1 = x - centroidX1;
            double ymcy1 = y - centroidY1;
           
            double xr = centroidX1 * scale
                + ((xmcx1 * scaleTimesCosine) + (ymcy1 * scaleTimesSine));

            double yr = centroidY1 * scale
                + ((-xmcx1 * scaleTimesSine) + (ymcy1 * scaleTimesCosine));

            int xt = (int)Math.round(xr + params.getTranslationX());
            int yt = (int)Math.round(yr + params.getTranslationY());
            
            int x2 = matched2.getX(i);
            int y2 = matched2.getY(i);
               
            double d = Math.sqrt(Math.pow(xt - x2, 2) + Math.pow(yt - y2, 2));
                
            diffFromModel[i] = d;
                
            avgDiffModel[0] += d;
        }
        
        avgDiffModel[0] /= (double)matched1.getN();        
    }

    /**
     * given unordered unmatched points set1 and set2 from image1 and image2
     * respectively, apply the transformation to set 1 and then find the
     * nearest matching in set2 within tolerance and return the 
     * matched information as a two dimensional array of 
     * {index from set1, index from set2, diff of point in set2 from
     * model generation by point in set1}
     * 
     * @param set1
     * @param set2
     * @param fit
     * @param tolerance
     * @param centroidX1
     * @param centroidY1
     * @return a two dimensional array holding the matched indexes and
     * the distances between the model and the point for that pair.
     * each row holds float[]{idx1, idx2, diff}
     */
    public float[][] calculateMatchUsingOptimal(
        PairIntArray set1, PairIntArray set2, TransformationPointFit fit,
        double tolerance, int centroidX1, int centroidY1) {
        
        double scale = fit.getParameters().getScale();
        double rotation = fit.getParameters().getRotationInRadians();

        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");

            log.severe("scale should not be smaller than 1");
        }

        double scaleTimesCosine = scale * Math.cos(rotation);
        double scaleTimesSine = scale * Math.sin(rotation);
        float transX = fit.getParameters().getTranslationX();
        float transY = fit.getParameters().getTranslationY();

        int nPoints1 = set1.getN();
        int nPoints2 = set2.getN();
        
        float[][] diffsAsCost = new float[nPoints1][nPoints2];

        // the algorithm modifies diffsAsCost, so make a copy
        float[][] diffsAsCostCopy = new float[nPoints1][nPoints2];
        
        for (int i = 0; i < set1.getN(); i++) {

            diffsAsCost[i] = new float[nPoints2];
            diffsAsCostCopy[i] = new float[nPoints2];
            
            int x = set1.getX(i);
            int y = set1.getY(i);

            double xmcx1 = x - centroidX1;
            double ymcy1 = y - centroidY1;

            double xr = centroidX1 * scale
                + ((xmcx1 * scaleTimesCosine) + (ymcy1 * scaleTimesSine));

            double yr = centroidY1 * scale
                + ((-xmcx1 * scaleTimesSine) + (ymcy1 * scaleTimesCosine));

            int xt = (int) Math.round(xr + transX);
            int yt = (int) Math.round(yr + transY);

            for (int j = 0; j < set2.getN(); j++) {

                int x2 = set2.getX(j);
                int y2 = set2.getY(j);

                double dist = Math.sqrt(Math.pow(xt - x2, 2) 
                    + Math.pow(yt - y2, 2));

                diffsAsCost[i][j] = (float)dist;
                    
                diffsAsCostCopy[i][j] = (float)dist;
            }         
        }
        
        boolean transposed = false;
        if (nPoints1 > nPoints2) {
            diffsAsCostCopy = MatrixUtil.transpose(diffsAsCostCopy);
            transposed = true;
        }
        
        HungarianAlgorithm b = new HungarianAlgorithm();
        int[][] match = b.computeAssignments(diffsAsCostCopy);
        
        // count the number of matches
        int count = 0;
        for (int i = 0; i < match.length; i++) {
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            count++;
        }
        
        float[][] output = new float[count][];
        
        count = 0;
        
        for (int i = 0; i < match.length; i++) {
           
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            
            double diff;
            
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }
            
            diff = diffsAsCost[idx1][idx2];

            output[count] = new float[3];
            output[count][0] = idx1;
            output[count][1] = idx2;
            output[count][2] = (float)diff;
            
            count++;
        }

        return output;
    }
    
    /**
     * given unordered unmatched points set1 and set2 from image1 and image2
     * respectively, apply the transformation to set 1 and then find the
     * nearest matching in set2 within transXTol and transYTol.
     * Note that the output matches if any will be in the argument variables
     * outputMatched1 and outputMatched2;
     * 
     * @param set1
     * @param set2
     * @param tolerance
     * @param matchedIndexesAndDiffs two dimensional array holding the matched 
     * indexes and the distances between the model and the point for that pair.
     * each row holds {idx1, idx2, diff}
     * @param outputMatched1 the container to hold the output matching points
     * for image 1 that are paired with outputMatched2 as a result of running
     * this method.
     * @param outputMatched2 the container to hold the output matching points
     * for image 2 that are paired with outputMatched1 as a result of running
     * this method.
     */
    public void matchPoints(
        PairIntArray set1, PairIntArray set2,
        double tolerance,
        float[][] matchedIndexesAndDiffs,
        PairIntArray outputMatched1, PairIntArray outputMatched2) {
                    
        for (int i = 0; i < matchedIndexesAndDiffs.length; i++) {
            
            int idx1 = (int)matchedIndexesAndDiffs[i][0];
            int idx2 = (int)matchedIndexesAndDiffs[i][1];
            float diff = matchedIndexesAndDiffs[i][2];
            
            if (diff < tolerance) {
                outputMatched1.add(set1.getX(idx1), set1.getY(idx1));
                outputMatched2.add(set2.getX(idx2), set2.getY(idx2));
            }
        }
    }
    
    /**
     * given unordered unmatched points set1 and set2 from image1 and image2
     * respectively, apply the transformation to set 1 and then find the
     * nearest matching in set2 within transXTol and transYTol.
     * Note that the output matches if any will be in the argument variables
     * outputMatched1 and outputMatched2;
     * 
     * @param set1
     * @param set2
     * @param fit
     * @param transXTol
     * @param transYTol
     * @param centroidX1
     * @param centroidY1
     * @param outputMatched1 the container to hold the output matching points
     * for image 1 that are paired with outputMatched2 as a result of running
     * this method.
     * @param outputMatched2 the container to hold the output matching points
     * for image 2 that are paired with outputMatched1 as a result of running
     * this method.
     */
    public void matchPointsUsingModelWithGreedyMatch(
        PairIntArray set1, PairIntArray set2, TransformationPointFit fit,
        double transXTol, double transYTol, int centroidX1, int centroidY1,
        PairIntArray outputMatched1, PairIntArray outputMatched2) {

        /*
         implements a linear search for the nearest neighbors of each point
         within translation tolerances.        
         */
        double scale = fit.getParameters().getScale();
        double rotation = fit.getParameters().getRotationInRadians();

        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");

            log.severe("scale should not be smaller than 1");
        }

        double scaleTimesCosine = scale * Math.cos(rotation);
        double scaleTimesSine = scale * Math.sin(rotation);
        float transX = fit.getParameters().getTranslationX();
        float transY = fit.getParameters().getTranslationY();       

        Set<Integer> chosen = new HashSet<Integer>();
        
        for (int i = 0; i < set1.getN(); i++) {

            int x = set1.getX(i);
            int y = set1.getY(i);

            double xmcx1 = x - centroidX1;
            double ymcy1 = y - centroidY1;

            double xr = centroidX1 * scale
                + ((xmcx1 * scaleTimesCosine) + (ymcy1 * scaleTimesSine));

            double yr = centroidY1 * scale
                + ((-xmcx1 * scaleTimesSine) + (ymcy1 * scaleTimesCosine));

            int xt = (int) Math.round(xr + transX);
            int yt = (int) Math.round(yr + transY);

            int lowerX = xt - (int) transXTol;
            int higherX = xt + (int) transXTol;
            int lowerY = yt - (int) transYTol;
            int higherY = yt + (int) transYTol;

            int bestX2 = -1;
            int bestY2 = -1;
            double bestDiff = Double.MAX_VALUE;
            int best2Idx = -1;

            for (int j = 0; j < set2.getN(); j++) {

                if (chosen.contains(Integer.valueOf(j))) {
                    continue;
                }
                
                int x2 = set2.getX(j);
                int y2 = set2.getY(j);

                if ((x2 < lowerX) || (x2 > higherX) || (y2 < lowerY)
                    || (y2 > higherY)) {

                    continue;
                }

                double d = Math.sqrt(Math.pow(xt - x2, 2) + Math.pow(yt - y2, 2)
                );

                if (d < bestDiff) {
                    bestDiff = d;
                    bestX2 = x2;
                    bestY2 = y2;
                    best2Idx = j;
                }
            }
            
            if (bestDiff < Double.MAX_VALUE) {

                outputMatched1.add(x, y);
                outputMatched2.add(bestX2, bestY2);
                
                chosen.add(Integer.valueOf(best2Idx));
            }            
        }
    }
    
    /**
     * sort the fits by descending number of matches. 
     * @param fits
     * @param idxLo
     * @param idxHi 
     */
    void sortByDescendingMatches(TransformationPointFit[] fits, int idxLo, 
        int idxHi) {
        
        if (idxLo < idxHi) {

            int idxMid = partition(fits, idxLo, idxHi);

            sortByDescendingMatches(fits, idxLo, idxMid - 1);

            sortByDescendingMatches(fits, idxMid + 1, idxHi);
        }
    }

    private int partition(TransformationPointFit[] fits, int idxLo, int idxHi) {
        
        TransformationPointFit x = fits[idxHi];
        
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            if (fitIsBetter(x, fits[i])) {
                store++;
                TransformationPointFit swap = fits[store];
                fits[store] = fits[i];
                fits[i] = swap;
            }
        }
        
        store++;
        TransformationPointFit swap = fits[store];
        fits[store] = fits[idxHi];
        fits[idxHi] = swap;
        
        return store;
    }

    /**
     * create an array of values starting with start and ending before stop
     * with separations of size delta.
     * @param start start of intervals, inclusive
     * @param stop end of intervals, exclusive
     * @param delta difference between returned sequential values
     * @return 
     */
    double[] createIntervals(double start, double stop, double delta) {
        
        int n = (int)Math.ceil((stop - start)/delta);
        
        double[] values = new double[n];
        
        int count = 0;
        for (double v = start; v < stop; v += delta) { 
            values[count] = v;
            count ++;
        }
        
        return values;
    }

    TransformationPointFit evaluateFit(
        float[] scaledRotatedX, float[] scaledRotatedY, float transX, 
        float transY, float tolTransX, float tolTransY,
        PairIntArray set2, 
        final float scale, final float rotationRadians,
        boolean setsAreMatched) {
        
        if (set2 == null) {
            throw new IllegalArgumentException(
            "set2 cannot be null");
        }
        if (scaledRotatedX == null || scaledRotatedY == null) {
            throw new IllegalArgumentException(
            "neither scaledRotatedX nor scaledRotatedY cannot be null");
        }
        
        if (setsAreMatched) {
            return evaluateFitForMatched(scaledRotatedX, scaledRotatedY, 
                transX, transY, set2, scale, rotationRadians);
        }
      
        return evaluateFitForUnMatched(scaledRotatedX, scaledRotatedY, 
            transX, transY, tolTransX, tolTransY,
            set2, scale, rotationRadians);
    }
    
    /**
     * given the model x y that have already been scaled and rotated, add the
     * transX and transY, respectively and calculated the average residual
     * between that and set2 and the standard deviation from the average.
     * Note that set2 and (scaledRotatedX, scaledRotatedY) are known to be
     * matched points so index 1 in set2 and index1 in scaledRotated represent
     * matching points.
     * @param set2
     * @param scaledRotatedX the model x points scaled and rotated
     * @param scaledRotatedY the model y points scaled and rotated
     * @param transX the x translation to apply to the model points
     * @param transY the y translation to apply to the model points
     * @return 
     */
    private TransformationPointFit evaluateFitForMatched( 
        float[] scaledRotatedX, float[] scaledRotatedY, 
        float transX, float transY, PairIntArray set2, 
        final float scale, final float rotationRadians) {
        
        if (set2 == null) {
            throw new IllegalArgumentException(
            "set2 cannot be null");
        }
        if (scaledRotatedX == null || scaledRotatedY == null) {
            throw new IllegalArgumentException(
            "neither scaledRotatedX nor scaledRotatedY cannot be null");
        }
        if (scaledRotatedX.length != scaledRotatedY.length) {
            throw new IllegalArgumentException(
            "scaledRotated X and Y must be the same length");
        }
        if (set2.getN() != scaledRotatedX.length) {
            throw new IllegalArgumentException(
            "for matched sets, set2 must be the same length as scaledRotated" 
            + " X and Y");
        }
        
        double sum = 0;
        double[] diff = new double[set2.getN()];
        
        for (int i = 0; i < set2.getN(); i++) {
            
            float modelX = scaledRotatedX[i] + transX;
            float modelY = scaledRotatedY[i] + transY;
            
            double dx = set2.getX(i) - modelX;
            double dy = set2.getY(i) - modelY;
            
            diff[i] = Math.sqrt(dx*dx + dy+dy);
            
            sum += diff[i];
        }
        
        double avgDiff = sum/(double)set2.getN();
        
        sum = 0;
        for (int i = 0; i < set2.getN(); i++) {
            sum += (diff[i] * diff[i]);
        }
        double stDev = Math.sqrt(sum/((double)set2.getN() - 1));
        
        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians(rotationRadians);
        params.setScale(scale);
        params.setTranslationX(transX);
        params.setTranslationY(transY);
        
        TransformationPointFit fit = new TransformationPointFit(params, 
            set2.getN(), avgDiff, stDev, Double.MAX_VALUE);
        
        return fit;
    }

    /**
     * given the model x y that have already been scaled and rotated, add the
     * transX and transY, respectively and calculated the average residual
     * between that and set2 and the standard deviation from the average.
     * Note that set2 and (scaledRotatedX, scaledRotatedY) are NOT known to be
     * matched points so the residuals are minimized for each point in
     * the model to find the matching in set2 before computing the
     * average and standard deviation.
     * 
     * @param set2
     * @param scaledRotatedX the model x points scaled and rotated
     * @param scaledRotatedY the model y points scaled and rotated
     * @param transX the x translation to apply to the model points
     * @param transY the y translation to apply to the model points
     * @return 
     */
    private TransformationPointFit evaluateFitForUnMatched( 
        float[] scaledRotatedX, float[] scaledRotatedY, 
        float transX, float transY, 
        float tolTransX, float tolTransY,
        PairIntArray set2, 
        final float scale, final float rotationRadians) {
        
        //return evaluateFitForUnMatchedOptimal(scaledRotatedX,
        //    scaledRotatedY, transX, transY, tolTransX, tolTransY, 
        //    set2, scale, rotationRadians);
    
        return evaluateFitForUnMatchedGreedy(scaledRotatedX,
            scaledRotatedY, transX, transY, 
            tolTransX, tolTransY,
            set2, scale, rotationRadians);
    }
    
    private TransformationPointFit evaluateFitForUnmatchedTransformed(
        TransformationParameters params, PairFloatArray unmatched1Transformed, 
        PairIntArray unmatched2, double tolTransX, double tolTransY) {
        
        
        return evaluateFitForUnMatchedTransformedOptimal(params, 
            unmatched1Transformed, unmatched2, tolTransX, tolTransY);
        
        
        /*
        return evaluateFitForUnMatchedTransformedGreedy(params, 
            unmatched1Transformed, unmatched2, tolTransX, tolTransY);
        */
    }
    
    private TransformationPointFit evaluateFitForUnMatchedTransformedOptimal(
        TransformationParameters params, PairFloatArray unmatched1Transformed, 
        PairIntArray unmatched2, double tolTransX, double tolTransY) {
        
        if (unmatched1Transformed == null) {
            throw new IllegalArgumentException(
            "unmatched1Transformed cannot be null");
        }
        if (unmatched2 == null) {
            throw new IllegalArgumentException(
            "unmatched2 cannot be null");
        }
        
        int n1 = unmatched1Transformed.getN();
        int n2 = unmatched2.getN();
        
        float[][] diffsAsCost = new float[n1][];
        float[][] diffsAsCostCopy = new float[n1][];
        
        for (int i = 0; i < n1; i++) {
            
            float modelX = unmatched1Transformed.getX(i);
            float modelY = unmatched1Transformed.getY(i);
            
            diffsAsCost[i] = new float[n2];
            diffsAsCostCopy[i] = new float[n2];
            
            for (int j = 0; j < n2; j++) {
                
                float dx = unmatched2.getX(j) - modelX;
                float dy = unmatched2.getY(j) - modelY;
                
                float diff = (float)Math.sqrt(dx*dx + dy*dy);
                
                diffsAsCost[i][j] = diff;
                diffsAsCostCopy[i][j] = diff;
            }
        }
        
        double tolerance = Math.sqrt(tolTransX*tolTransX + tolTransY*tolTransY);

        boolean transposed = false;
        if (diffsAsCostCopy.length > diffsAsCostCopy[0].length) {
            diffsAsCostCopy = MatrixUtil.transpose(diffsAsCostCopy);
            transposed = true;
        }
        
        HungarianAlgorithm b = new HungarianAlgorithm();
        
        // modifies diffsAsCostCopy:
        int[][] match = b.computeAssignments(diffsAsCostCopy);
        
        double sum = 0;
        int nMatched = 0;
        for (int i = 0; i < match.length; i++) {
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }
            
            double diff = diffsAsCost[idx1][idx2];
            
            if (diff > tolerance) {
                continue;
            }
            
            sum += diff;
            nMatched++;
        }
        
        double avgDiff = sum/(double)nMatched;
        
        sum = 0;
        for (int i = 0; i < match.length; i++) {
            
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }
            
            double dist = diffsAsCost[idx1][idx2];
            
            if (dist > tolerance) {
                continue;
            }
            
            double diff = dist - avgDiff;
            
            sum += (diff * diff);
        }
        
        double stDev = Math.sqrt(sum/((double)nMatched - 1.));
        
        TransformationPointFit fit = new TransformationPointFit(params, nMatched, 
            avgDiff, stDev, tolerance);
        
        return fit;
    }
    
    /**
     * given the model x y that have already been scaled and rotated, add the
     * transX and transY, respectively and calculated the average residual
     * between that and set2 and the standard deviation from the average.
     * Note that set2 and (scaledRotatedX, scaledRotatedY) are NOT known to be
     * matched points so the residuals are minimized for each point in
     * the model to find the matching in set2 before computing the
     * average and standard deviation.
     * 
     * @param set2
     * @param scaledRotatedX the model x points scaled and rotated
     * @param scaledRotatedY the model y points scaled and rotated
     * @param transX the x translation to apply to the model points
     * @param transY the y translation to apply to the model points
     * @return 
     */
    private TransformationPointFit evaluateFitForUnMatchedOptimal( 
        float[] scaledRotatedX, float[] scaledRotatedY, 
        float transX, float transY, float tolTransX, float tolTransY,
        PairIntArray set2, 
        final float scale, final float rotationRadians) {
        
        if (set2 == null) {
            throw new IllegalArgumentException(
            "set2 cannot be null");
        }
        if (scaledRotatedX == null || scaledRotatedY == null) {
            throw new IllegalArgumentException(
            "neither scaledRotatedX nor scaledRotatedY cannot be null");
        }
        if (scaledRotatedX.length != scaledRotatedY.length) {
            throw new IllegalArgumentException(
            "scaledRotated X and Y must be the same length");
        }
        
        float[][] diffsAsCost = new float[scaledRotatedX.length][];
        float[][] diffsAsCostCopy = new float[scaledRotatedX.length][];
        
        for (int i = 0; i < scaledRotatedX.length; i++) {
            
            float modelX = scaledRotatedX[i] + transX;
            float modelY = scaledRotatedY[i] + transY;
            
            diffsAsCost[i] = new float[set2.getN()];
            diffsAsCostCopy[i] = new float[set2.getN()];
            
            for (int j = 0; j < set2.getN(); j++) {
                
                float dx = set2.getX(j) - modelX;
                float dy = set2.getY(j) - modelY;
                
                float diff = (float)Math.sqrt(dx*dx + dy*dy);
                
                diffsAsCost[i][j] = diff;
                diffsAsCostCopy[i][j] = diff;
            }
        }
        
        double tolerance = Math.sqrt(tolTransX*tolTransX + tolTransY*tolTransY);

        boolean transposed = false;
        if (diffsAsCostCopy.length > diffsAsCostCopy[0].length) {
            diffsAsCostCopy = MatrixUtil.transpose(diffsAsCostCopy);
            transposed = true;
        }
        
        HungarianAlgorithm b = new HungarianAlgorithm();
        
        // modifies diffsAsCostCopy:
        int[][] match = b.computeAssignments(diffsAsCostCopy);
        
        double sum = 0;
        int nMatched = 0;
        for (int i = 0; i < match.length; i++) {
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }
            
            double diff = diffsAsCost[idx1][idx2];
            
            if (diff > tolerance) {
                
                continue;
            }
            
            sum += diff;
            nMatched++;
        }
        
        double avgDiff = sum/(double)nMatched;
        
        sum = 0;
        for (int i = 0; i < match.length; i++) {
            
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }
            
            double dist = diffsAsCost[idx1][idx2];
            
            if (dist > tolerance) {
                continue;
            }
            
            double diff = dist - avgDiff;
            
            sum += (diff * diff);
        }
        
        double stDev = Math.sqrt(sum/((double)nMatched - 1));
        
        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians(rotationRadians);
        params.setScale(scale);
        params.setTranslationX(transX);
        params.setTranslationY(transY);
        
        TransformationPointFit fit = new TransformationPointFit(params, nMatched, 
            avgDiff, stDev, tolerance);
        
        return fit;
    }
    
    private TransformationPointFit evaluateFitForUnMatchedTransformedGreedy(
        TransformationParameters params, PairFloatArray unmatched1Transformed, 
        PairIntArray unmatched2, double tolTransX, double tolTransY) {
        
        if (unmatched1Transformed == null) {
            throw new IllegalArgumentException(
            "unmatched1Transformed cannot be null");
        }
        if (unmatched2 == null) {
            throw new IllegalArgumentException(
            "unmatched2 cannot be null");
        }
        
        int n = unmatched1Transformed.getN();
        
        Set<Integer> chosen = new HashSet<Integer>();
        
        double[] diffs = new double[n];
        int nMatched = 0;
        double avg = 0;
        
        for (int i = 0; i < n; i++) {
            
            float modelX = unmatched1Transformed.getX(i);
            float modelY = unmatched1Transformed.getY(i);
            
            double minDiff = Double.MAX_VALUE;
            int min2Idx = -1;

            for (int j = 0; j < unmatched2.getN(); j++) {

                if (chosen.contains(Integer.valueOf(j))) {
                    continue;
                }
                
                float dx = unmatched2.getX(j) - modelX;
                float dy = unmatched2.getY(j) - modelY;
                
                float diff = (float)Math.sqrt(dx*dx + dy*dy);
                
                if ((Math.abs(dx) > tolTransX) || (Math.abs(dy) > tolTransY)) {
                    
                    continue;
                } 
                
                if (diff < minDiff) {
                    minDiff = diff;
                    min2Idx = j;
                }
            }
            
            if (minDiff < Double.MAX_VALUE) {
                diffs[nMatched] = minDiff;
                nMatched++;
                chosen.add(Integer.valueOf(min2Idx));
                avg += minDiff;
            }
        }
        
        avg = (nMatched == 0) ? Double.MAX_VALUE :
            avg / (double)nMatched;
        
        double stDev = 0;
        for (int i = 0; i < nMatched; i++) {
            double d = diffs[i] - avg;
            stDev += (d * d);
        }
        
        stDev = (nMatched == 0) ? Double.MAX_VALUE :
            Math.sqrt(stDev/((double)nMatched - 1.));
        
        TransformationPointFit fit = new TransformationPointFit(params, nMatched, 
            avg, stDev, Math.sqrt(tolTransX*tolTransX + tolTransY*tolTransY)
        );
        
        return fit;
    }
    
    /**
     * given the model x y that have already been scaled and rotated, add the
     * transX and transY, respectively and calculated the average residual
     * between that and set2 and the standard deviation from the average.
     * Note that set2 and (scaledRotatedX, scaledRotatedY) are NOT known to be
     * matched points so the residuals are minimized for each point in
     * the model to find the matching in set2 before computing the
     * average and standard deviation.
     * 
     * @param set2
     * @param scaledRotatedX the model x points scaled and rotated
     * @param scaledRotatedY the model y points scaled and rotated
     * @param transX the x translation to apply to the model points
     * @param transY the y translation to apply to the model points
     * @return 
     */
    private TransformationPointFit evaluateFitForUnMatchedGreedy( 
        float[] scaledRotatedX, float[] scaledRotatedY, 
        float transX, float transY, 
        float tolTransX, float tolTransY,
        PairIntArray set2, 
        final float scale, final float rotationRadians) {
        
        if (set2 == null) {
            throw new IllegalArgumentException(
            "set2 cannot be null");
        }
        if (scaledRotatedX == null || scaledRotatedY == null) {
            throw new IllegalArgumentException(
            "neither scaledRotatedX nor scaledRotatedY cannot be null");
        }
        if (scaledRotatedX.length != scaledRotatedY.length) {
            throw new IllegalArgumentException(
            "scaledRotated X and Y must be the same length");
        }
        
        Set<Integer> chosen = new HashSet<Integer>();
        
        double[] diffs = new double[scaledRotatedX.length];
        int nMatched = 0;
        double avg = 0;
        
        for (int i = 0; i < scaledRotatedX.length; i++) {
            
            float modelX = scaledRotatedX[i] + transX;
            float modelY = scaledRotatedY[i] + transY;
            
            double minDiff = Double.MAX_VALUE;
            int min2Idx = -1;

            for (int j = 0; j < set2.getN(); j++) {

                if (chosen.contains(Integer.valueOf(j))) {
                    continue;
                }
                
                float dx = set2.getX(j) - modelX;
                float dy = set2.getY(j) - modelY;
                
                float diff = (float)Math.sqrt(dx*dx + dy*dy);
                
                if ((Math.abs(dx) > tolTransX) || (Math.abs(dy) > tolTransY)) {
                    
                    continue;
                } 
                
                if (diff < minDiff) {
                    minDiff = diff;
                    min2Idx = j;
                }
            }
            
            if (minDiff < Double.MAX_VALUE) {
                diffs[nMatched] = minDiff;
                nMatched++;
                chosen.add(Integer.valueOf(min2Idx));
                avg += minDiff;
            }
        }
        
        avg = (nMatched == 0) ? Double.MAX_VALUE :
            avg / (double)nMatched;
        
        double stDev = 0;
        for (int i = 0; i < nMatched; i++) {
            double d = diffs[i] - avg;
            stDev += (d * d);
        }
        
        stDev = (nMatched == 0) ? Double.MAX_VALUE :
            Math.sqrt(stDev/((double)nMatched - 1.));
        
        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians(rotationRadians);
        params.setScale(scale);
        params.setTranslationX(transX);
        params.setTranslationY(transY);
       
        TransformationPointFit fit = new TransformationPointFit(params, nMatched, 
            avg, stDev, Math.sqrt(tolTransX*tolTransX + tolTransY*tolTransY)
        );
        
        return fit;
    }
    
    /**
     * given the scale, rotation and set 1's reference frame centroids,
     * calculate the translation between set1 and set2 assuming that not all
     * points will match.  transXTol and transYTol allow a tolerance when
     * matching the predicted position of a point in set2.
     * 
     * It's expected that the invoker of this method is trying to solve for
     * translation for sets of points like corners in images.  This assumption
     * means that the number of point pair combinations is always far less
     * than the pixel combinations of translations over x and y.
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * 
     * ALSO NOTE: if you know a better solution exists for translation 
     * parameters that matches fewer points, but has a small avg dist from
     * model and smaller standard deviation from the avg dist from model,
     * then transXTol and transYTol should be set to a smaller value and passed
     * to this method.
     * @param params transformation parameters to apply to matched1
     * @param matched1Transformed
     * @param matched2 set of points from image 2 that are matched to points in
     * matched1 with same indexes
     * @param centroidX1 the x coordinate of the center of image 1
     * @param centroidY1 the y coordinate of the center of image 1
     * @return 
     */
    public TransformationPointFit evaluateFitForMatchedTransformed(
        TransformationParameters params, PairFloatArray matched1Transformed, 
        PairIntArray matched2) {
        
        if (matched1Transformed == null) {
            throw new IllegalArgumentException(
                "matched1Transformed cannot be null");
        }
        if (matched2 == null) {
            throw new IllegalArgumentException("matched2 cannot be null");
        }
        if (matched1Transformed.getN() != matched2.getN()) {
            throw new IllegalArgumentException(
                "the 2 point sets must have the same length");
        }
        
        double[] diff = new double[matched1Transformed.getN()];
        
        double avg = 0;
                
        for (int i = 0; i < matched1Transformed.getN(); i++) {
            
            float modelX = matched1Transformed.getX(i);
            float modelY = matched1Transformed.getY(i);
            
            int x2 = matched2.getX(i);
            int y2 = matched2.getY(i);
            
            double dx = x2 - modelX;
            double dy = y2 - modelY;
                            
            diff[i] = Math.sqrt(dx*dx + dy*dy);
            
            avg += diff[i];
        }
        
        avg /= (double)matched2.getN();
        
        double stdDev = 0;
        
        for (int i = 0; i < matched2.getN(); i++) {
            
            double d = diff[i] - avg;
            
            stdDev += (d * d);
        }
        
        stdDev = Math.sqrt(stdDev/((double)matched2.getN() - 1));
        
        TransformationPointFit fit = new TransformationPointFit(params,
            matched2.getN(), avg, stdDev, Double.POSITIVE_INFINITY);
                
        return fit;
    }

    private TransformationPointFit evaluateFit(
        PairIntArray[] edges1, PairIntArray[] edges2, 
        float translationX, float translationY, 
        float tolTransX, float tolTransY, float scale, 
        float rotationRadians, int centroidX1, int centroidY1) {
        
        List<Double> residuals = new ArrayList<Double>();
        
        int nTotal = 0;
        
        for (int i = 0; i < edges1.length; i++) {
            
            PairIntArray edge1 = edges1[i];
            
            PairIntArray edge2 = edges2[i];
        
            nTotal += edge1.getN();
            
            calculateTranslationResidualsForUnmatchedMultiplicity( 
                edge1, edge2, translationX, translationY, tolTransX, tolTransY,
                scale, rotationRadians,
                centroidX1, centroidY1, residuals);
        }
        
        double avg = 0;
        for (Double diff : residuals) {
            avg += diff.doubleValue();
        }
        avg /= (double)residuals.size();
        
        double stdDev = 0;
        for (Double diff : residuals) {
            double d = diff.doubleValue() - avg;
            stdDev += (d * d);
        }
        stdDev = Math.sqrt(stdDev/((double)residuals.size() - 1));
        
        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians(rotationRadians);
        params.setScale(scale);
        params.setTranslationX(translationX);
        params.setTranslationY(translationY);
       
        TransformationPointFit fit = new TransformationPointFit(params, 
            residuals.size(), avg, stdDev, 
            Math.sqrt(tolTransX*tolTransX + tolTransY*tolTransY)
        );
        
        return fit;
    }
    
    /**
     * calculate the residuals between edge1 points and edge2 points with the
     * assumption that there may not be one to one point matches, but instead
     * will be general point matches.  For that reason, a greedy search for
     * nearest neighbor with the possibility of multiple matches to a neighbor 
     * is used.
     * 
     * @param edge1
     * @param edge2
     * @param transX
     * @param transY
     * @param tolTransX
     * @param tolTransY
     * @param scale
     * @param rotationRadians
     * @param outputResiduals
     */
    private void calculateTranslationResidualsForUnmatchedMultiplicity( 
        PairIntArray edge1, PairIntArray edge2,
        float transX, float transY, float tolTransX, float tolTransY,
        final float scale, final float rotationRadians,
        int centroidX1, int centroidY1, List<Double> outputResiduals) {
        
        if (edge1 == null) {
            throw new IllegalArgumentException(
            "edge1 cannot be null");
        }
        if (edge2 == null) {
            throw new IllegalArgumentException(
            "edge2 cannot be null");
        }
        
        float scaleTimesCosine = (float)(scale * Math.cos(rotationRadians));
        float scaleTimesSine = (float)(scale * Math.sin(rotationRadians));
                       
        for (int i = 0; i < edge1.getN(); i++) {
            
            int x1 = edge1.getX(i);
            int y1 = edge1.getY(i);
            
            float modelX = (centroidX1*scale + ( 
                ((x1 - centroidX1) * scaleTimesCosine) +
                ((y1 - centroidY1) * scaleTimesSine))) + transX;
            
            float modelY = (centroidY1*scale + ( 
                (-(x1 - centroidX1) * scaleTimesSine) +
                ((y1 - centroidY1) * scaleTimesCosine))) + transY;
            
            double minDiff = Double.MAX_VALUE;

            for (int j = 0; j < edge2.getN(); j++) {

                int x2 = edge2.getX(j);
                int y2 = edge2.getY(j);
                
                float dx = x2 - modelX;
                float dy = y2 - modelY;
                
                float diff = (float)Math.sqrt(dx*dx + dy*dy);
                
                if ((Math.abs(dx) > tolTransX) || (Math.abs(dy) > tolTransY)) {
                    
                    continue;
                } 
                
                if (diff < minDiff) {
                    minDiff = diff;
                }
            }
            
            if (minDiff < Double.MAX_VALUE) {
                outputResiduals.add(Double.valueOf(minDiff));
            }
        }       
    }

    private TransformationPointFit refineTranslationWithDownhillSimplex(
        float[] scaledRotatedX, float[] scaledRotatedY, 
        float transX, float transY, float tolTransX, float tolTransY,
        float plusMinusTransX, float plusMinusTransY,
        PairIntArray set2, final float scale, final float rotationRadians,
        boolean setsAreMatched) {
                    
        //TODO: set this empirically from tests
        double convergence = 0;

        float txMin = transX - plusMinusTransX;
        float txMax = transX + plusMinusTransX;
        float tyMin = transY - plusMinusTransY;
        float tyMax = transY + plusMinusTransY;

        // the positive offsets can be found w/ reflection
        float[] dtx = new float[]{
            -plusMinusTransX, -0.5f*plusMinusTransX, 
            -0.25f*plusMinusTransX,
            -0.125f*plusMinusTransX, -0.0625f*plusMinusTransX,
            -0.03125f*plusMinusTransX
        };
        float[] dty = new float[]{
            -plusMinusTransY, -0.5f*plusMinusTransY, 
            -0.25f*plusMinusTransY,
            -0.125f*plusMinusTransY, -0.0625f*plusMinusTransY,
            -0.03125f*plusMinusTransY
        };
        int n = (1 + dtx.length) * (1 + dty.length);

        TransformationPointFit[] fits = new TransformationPointFit[n];

        int count = 0;
        for (int i = 0; i <= dtx.length; i++) {

            float tx = (i == 0) ? transX : (transX + dtx[i - 1]);

            for (int j = 0; j <= dty.length; j++) {

                float ty = (i == 0) ? transY : (transY + dty[i - 1]);

                fits[count] = evaluateFit(scaledRotatedX, scaledRotatedY, 
                    tx, ty, tolTransX, tolTransY,
                    set2, scale, rotationRadians, setsAreMatched);
                   
                if (fits[count] != null) {
                    count++;
                }
            }
        }

        if (count < n) {
            fits = Arrays.copyOf(fits, count);
        }

        float alpha = 1;   // > 0
        float gamma = 2;   // > 1
        float beta = -0.5f;
        float tau = 0.5f;

        boolean go = true;

        int nMaxIter = 100;
        int nIter = 0;

        int bestFitIdx = 0;
        int secondWorstFitIdx = fits.length - 2;
        int worstFitIdx = fits.length - 1;

        int lastNMatches = Integer.MIN_VALUE;
        double lastAvgDistModel = Double.MAX_VALUE;
        int nIterSameMin = 0;

        while (go && (nIter < nMaxIter)) {

            if (fits.length == 0) {
                break;
            }

            sortByDescendingMatches(fits, 0, (fits.length - 1));

            /*if (fits.length > 0) {
                log.info("best fit:  n=" + 
                    fits[bestFitIdx].getNumberOfMatchedPoints()
                    + " dm=" + fits[bestFitIdx].getMeanDistFromModel()
                    + " params:\n" 
                    + fits[bestFitIdx].getParameters().toString());
            }*/

            if ((lastNMatches == fits[bestFitIdx].getNumberOfMatchedPoints())
                && (Math.abs(lastAvgDistModel
                    - fits[bestFitIdx].getMeanDistFromModel()) < 0.01)) {

                nIterSameMin++;
                if (nIterSameMin >= 5) {
                    break;
                }

            } else {
                nIterSameMin = 0;
            }

            lastNMatches = fits[bestFitIdx].getNumberOfMatchedPoints();
            lastAvgDistModel = fits[bestFitIdx].getMeanDistFromModel();

            // determine center for all points excepting the worse fit
            float txSum = 0.0f;
            float tySum = 0.0f;
            for (int i = 0; i < (fits.length - 1); i++) {
                txSum += fits[i].getTranslationX();
                tySum += fits[i].getTranslationY();
            }
            transX = txSum/(float)(fits.length - 1);
            transY = tySum/(float)(fits.length - 1);

            // "Reflection"
            float txReflect = transX + (alpha
                * (transX - (float)fits[worstFitIdx].getTranslationX()));
            float tyReflect = transY + (alpha
                * (transY - (float)fits[worstFitIdx].getTranslationY()));

            TransformationPointFit fitReflected
                = evaluateFit(scaledRotatedX, scaledRotatedY, 
                    txReflect, tyReflect, 
                    tolTransX, tolTransY,
                    set2, 
                    scale, rotationRadians, setsAreMatched);
            
            boolean relectIsWithinBounds
                = (txReflect >= txMin) && (txReflect <= txMax)
                && (tyReflect >= tyMin) && (tyReflect <= tyMax);

            if (fitIsBetter(fits[secondWorstFitIdx], fitReflected)
                && relectIsWithinBounds
                && !fitIsBetter(fits[bestFitIdx], fitReflected)) {

                fits[worstFitIdx] = fitReflected;

            } else {

                if (fitIsBetter(fits[bestFitIdx], fitReflected)
                    && relectIsWithinBounds) {

                    // "Expansion"
                    float txExpansion = transX + (gamma
                        * (transX - (float)fits[worstFitIdx].getTranslationX()));
                    float tyExpansion = transY + (gamma
                        * (transY - (float)fits[worstFitIdx].getTranslationY()));
                    
                    TransformationPointFit fitExpansion
                        = evaluateFit(scaledRotatedX, scaledRotatedY, 
                            txExpansion, tyExpansion, 
                            tolTransX, tolTransY,
                            set2, 
                            scale, rotationRadians, setsAreMatched);

                    if (fitIsBetter(fitReflected, fitExpansion)
                        && (txExpansion >= txMin) && (txExpansion <= txMax)
                        && (tyExpansion >= tyMin) && (tyExpansion <= tyMax)) {

                        fits[worstFitIdx] = fitExpansion;

                    } else {

                        fits[worstFitIdx] = fitReflected;
                    }

                } else {

                    // the reflection fit is worse than the 2nd worse
                    // "Contraction"
                    float txContraction = transX + (beta
                        * (transX - (float)fits[worstFitIdx].getTranslationX()));
                    float tyContraction = transY + (beta
                        * (transY - (float)fits[worstFitIdx].getTranslationY()));

                    TransformationPointFit fitContraction
                        = evaluateFit(scaledRotatedX, scaledRotatedY, 
                            txContraction, tyContraction, 
                            tolTransX, tolTransY,
                            set2, 
                            scale, rotationRadians, setsAreMatched);

                    if (fitIsBetter(fits[worstFitIdx], fitContraction)
                        && (txContraction >= txMin) && (txContraction <= txMax)
                        && (tyContraction >= tyMin) && (tyContraction <= tyMax)) {

                        fits[worstFitIdx] = fitContraction;

                    } else {

                        // "Reduction"
                        for (int i = 1; i < fits.length; i++) {

                            float txReduction
                                = (float)(fits[bestFitIdx].getTranslationX()
                                + (tau
                                * (fits[i].getTranslationX()
                                - fits[bestFitIdx].getTranslationX())));

                            float tyReduction
                                = (float)(fits[bestFitIdx].getTranslationY()
                                + (tau
                                * (fits[i].getTranslationY()
                                - fits[bestFitIdx].getTranslationY())));

                            //NOTE: there's a possibility of a null fit.
                            //  instead of re-writing the fits array, will 
                            //  assign a fake infinitely bad fit which will 
                            //  fall to the bottom of the list after the next 
                            //  sort.
                            TransformationPointFit fit
                                = evaluateFit(scaledRotatedX, scaledRotatedY, 
                                  txReduction, tyReduction, 
                                  tolTransX, tolTransY,
                                  set2, 
                                  scale, rotationRadians, setsAreMatched);

                            if (fit != null) {
                                fits[i] = fit;
                            } else {
                                fits[i] = new TransformationPointFit(
                                    new TransformationParameters(),
                                    0, Double.MAX_VALUE, Double.MAX_VALUE,
                                    Double.MAX_VALUE
                                );
                            }
                        }
                    }
                }
            }

            log.finest("best fit so far: nMatches="
                + fits[bestFitIdx].getNumberOfMatchedPoints()
                + " diff from model=" + fits[bestFitIdx].getMeanDistFromModel()
            );

            nIter++;

            if ((fits[bestFitIdx].getNumberOfMatchedPoints() == convergence)
                && (fits[bestFitIdx].getMeanDistFromModel() == 0)) {
                go = false;
            } else if ((transX > txMax) || (transX < txMin)) {
                go = false;
            } else if ((transY > tyMax) || (transY < tyMin)) {
                go = false;
            }
        }

        return fits[bestFitIdx];
    }

    double findToleranceForOutliers(TransformationPointFit fit, 
        float[][] matchIndexesAndDiffs, int centroidX1, int centroidY1) {
        
        float[] diffs = new float[matchIndexesAndDiffs.length];
        for (int i = 0; i < matchIndexesAndDiffs.length; i++) {
            diffs[i] = matchIndexesAndDiffs[i][2];
        }
        
        HistogramHolder h = Histogram.createSimpleHistogram(diffs, 
            Errors.populateYErrorsBySqrt(diffs));
        /*
        try {
            h.plotHistogram("outliers", 4);
        } catch (IOException e) {
            e.printStackTrace();
        }
        */
        int idx = MiscMath.findYMaxIndex(h.getYHistFloat());

        //TODO: this may need to be improved
        
        float yLimit = 0.3f * h.getYHistFloat()[idx];
        
        int idxLimit = h.getYHist().length - 1;
        
        for (int i = (idx + 1); i < h.getYHist().length; i++) {
            if (h.getYHist()[i] < yLimit) {
                idxLimit = i;
                break;
            }
        }
        
        int idxLimit0 = idxLimit - 1;
        if (idxLimit0 < 0) {
            idxLimit0 = 0;
        }
        
        /*
        int idxLimit0 = idxLimit + 1;
        if (idxLimit0 > (h.getYHist().length - 1)) {
            idxLimit0 = h.getYHist().length - 1;
        }*/
        
        //double tolerance = (h.getXHist()[idxLimit] + h.getXHist()[idxLimit0])/2.;
        double tolerance = h.getXHist()[idxLimit];
        
        log.info("tolerance=" + tolerance);
        
        return tolerance;
    }

    private TransformationPointFit refineTranslationWithDownhillSimplex(
        final PairIntArray[] edges1, final PairIntArray[] edges2,
        float transX, float transY, float tolTransX, float tolTransY,
        float plusMinusTransX, float plusMinusTransY,
        final float scale, final float rotationRadians,
        int centroidX1, int centroidY1) {
                    
        double convergence = 0;

        float txMin = transX - plusMinusTransX;
        float txMax = transX + plusMinusTransX;
        float tyMin = transY - plusMinusTransY;
        float tyMax = transY + plusMinusTransY;

        // the positive offsets can be found w/ reflection
        float[] dtx = new float[]{
            -plusMinusTransX, -0.5f*plusMinusTransX, 
            -0.25f*plusMinusTransX,
            -0.125f*plusMinusTransX, -0.0625f*plusMinusTransX
        };
        float[] dty = new float[]{
            -plusMinusTransY, -0.5f*plusMinusTransY, 
            -0.25f*plusMinusTransY,
            -0.125f*plusMinusTransY, -0.0625f*plusMinusTransY
        };
        int n = (1 + dtx.length) * (1 + dty.length);

        TransformationPointFit[] fits = new TransformationPointFit[n];

        int count = 0;
        for (int i = 0; i <= dtx.length; i++) {

            float tx = (i == 0) ? transX : (transX + dtx[i - 1]);

            for (int j = 0; j <= dty.length; j++) {

                float ty = (i == 0) ? transY : (transY + dty[i - 1]);

                fits[count] = evaluateFit(edges1, edges2,
                    tx, ty, tolTransX, tolTransY,
                    scale, rotationRadians, centroidX1, centroidY1);
                   
                if (fits[count] != null) {
                    count++;
                }
            }
        }

        if (count < n) {
            fits = Arrays.copyOf(fits, count);
        }

        float alpha = 1;   // > 0
        float gamma = 2;   // > 1
        float beta = -0.5f;
        float tau = 0.5f;

        boolean go = true;

        int nMaxIter = 100;
        int nIter = 0;

        int bestFitIdx = 0;
        int secondWorstFitIdx = fits.length - 2;
        int worstFitIdx = fits.length - 1;

        int lastNMatches = Integer.MIN_VALUE;
        double lastAvgDistModel = Double.MAX_VALUE;
        int nIterSameMin = 0;

        while (go && (nIter < nMaxIter)) {

            if (fits.length == 0) {
                break;
            }

            sortByDescendingMatches(fits, 0, (fits.length - 1));

            /*if (fits.length > 0) {
                log.info("best fit:  n=" + 
                    fits[bestFitIdx].getNumberOfMatchedPoints()
                    + " dm=" + fits[bestFitIdx].getMeanDistFromModel()
                    + " params:\n" 
                    + fits[bestFitIdx].getParameters().toString());
            }*/

            if ((lastNMatches == fits[bestFitIdx].getNumberOfMatchedPoints())
                && (Math.abs(lastAvgDistModel
                    - fits[bestFitIdx].getMeanDistFromModel()) < 0.01)) {

                nIterSameMin++;
                if (nIterSameMin >= 5) {
                    break;
                }

            } else {
                nIterSameMin = 0;
            }

            lastNMatches = fits[bestFitIdx].getNumberOfMatchedPoints();
            lastAvgDistModel = fits[bestFitIdx].getMeanDistFromModel();

            // determine center for all points excepting the worse fit
            float txSum = 0.0f;
            float tySum = 0.0f;
            for (int i = 0; i < (fits.length - 1); i++) {
                txSum += fits[i].getTranslationX();
                tySum += fits[i].getTranslationY();
            }
            transX = txSum/(float)(fits.length - 1);
            transY = tySum/(float)(fits.length - 1);

            // "Reflection"
            float txReflect = transX + (alpha
                * (transX - (float)fits[worstFitIdx].getTranslationX()));
            float tyReflect = transY + (alpha
                * (transY - (float)fits[worstFitIdx].getTranslationY()));

            TransformationPointFit fitReflected
                = evaluateFit(edges1, edges2,
                    txReflect, tyReflect,  tolTransX, tolTransY,
                    scale, rotationRadians, centroidX1, centroidY1);
                    
            boolean relectIsWithinBounds
                = (txReflect >= txMin) && (txReflect <= txMax)
                && (tyReflect >= tyMin) && (tyReflect <= tyMax);

            if (fitIsBetter(fits[secondWorstFitIdx], fitReflected)
                && relectIsWithinBounds
                && !fitIsBetter(fits[bestFitIdx], fitReflected)) {

                fits[worstFitIdx] = fitReflected;

            } else {

                if (fitIsBetter(fits[bestFitIdx], fitReflected)
                    && relectIsWithinBounds) {

                    // "Expansion"
                    float txExpansion = transX + (gamma
                        * (transX - (float)fits[worstFitIdx].getTranslationX()));
                    float tyExpansion = transY + (gamma
                        * (transY - (float)fits[worstFitIdx].getTranslationY()));
                    
                    TransformationPointFit fitExpansion
                        = evaluateFit(edges1, edges2,
                            txExpansion, tyExpansion,  tolTransX, tolTransY,
                            scale, rotationRadians, centroidX1, centroidY1);

                    if (fitIsBetter(fitReflected, fitExpansion)
                        && (txExpansion >= txMin) && (txExpansion <= txMax)
                        && (tyExpansion >= tyMin) && (tyExpansion <= tyMax)) {

                        fits[worstFitIdx] = fitExpansion;

                    } else {

                        fits[worstFitIdx] = fitReflected;
                    }

                } else {

                    // the reflection fit is worse than the 2nd worse
                    // "Contraction"
                    float txContraction = transX + (beta
                        * (transX - (float)fits[worstFitIdx].getTranslationX()));
                    float tyContraction = transY + (beta
                        * (transY - (float)fits[worstFitIdx].getTranslationY()));

                    TransformationPointFit fitContraction
                        = evaluateFit(edges1, edges2,
                            txContraction, tyContraction, tolTransX, tolTransY,
                            scale, rotationRadians, centroidX1, centroidY1);

                    if (fitIsBetter(fits[worstFitIdx], fitContraction)
                        && (txContraction >= txMin) && (txContraction <= txMax)
                        && (tyContraction >= tyMin) && (tyContraction <= tyMax)) {

                        fits[worstFitIdx] = fitContraction;

                    } else {

                        // "Reduction"
                        for (int i = 1; i < fits.length; i++) {

                            float txReduction
                                = (float)(fits[bestFitIdx].getTranslationX()
                                + (tau
                                * (fits[i].getTranslationX()
                                - fits[bestFitIdx].getTranslationX())));

                            float tyReduction
                                = (float)(fits[bestFitIdx].getTranslationY()
                                + (tau
                                * (fits[i].getTranslationY()
                                - fits[bestFitIdx].getTranslationY())));

                            //NOTE: there's a possibility of a null fit.
                            //  instead of re-writing the fits array, will 
                            //  assign a fake infinitely bad fit which will 
                            //  fall to the bottom of the list after the next 
                            //  sort.
                            TransformationPointFit fit
                                = evaluateFit(edges1, edges2,
                                    txReduction, tyReduction, 
                                    tolTransX, tolTransY, 
                                    scale, rotationRadians, 
                                    centroidX1, centroidY1);
                                    
                            if (fit != null) {
                                fits[i] = fit;
                            } else {
                                fits[i] = new TransformationPointFit(
                                    new TransformationParameters(),
                                    0, Double.MAX_VALUE, Double.MAX_VALUE,
                                    Double.MAX_VALUE
                                );
                            }
                        }
                    }
                }
            }

            log.finest("best fit so far: nMatches="
                + fits[bestFitIdx].getNumberOfMatchedPoints()
                + " diff from model=" + fits[bestFitIdx].getMeanDistFromModel()
            );

            nIter++;

            if ((fits[bestFitIdx].getNumberOfMatchedPoints() == convergence)
                && (fits[bestFitIdx].getMeanDistFromModel() == 0)) {
                go = false;
            } else if ((transX > txMax) || (transX < txMin)) {
                go = false;
            } else if ((transY > tyMax) || (transY < tyMin)) {
                go = false;
            }
        }

        return fits[bestFitIdx];
    }

}
