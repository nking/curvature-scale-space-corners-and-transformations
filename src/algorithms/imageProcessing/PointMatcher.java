package algorithms.imageProcessing;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.compGeometry.PointPartitioner;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.ArrayPair;
import algorithms.util.Errors;
import algorithms.util.LinearRegression;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
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
 
    protected static int minTolerance = 5;
    
    // experimenting w/ use of model diff only for projective solutions.
    // problem w/ this is that a repeated pattern in an image can give
    //  multiple possible solutions
    private boolean costIsDiff = false;
    
    private boolean costIsNumAndDiff = false;
    
    //TODO: this has to be a high number for sets with projection.
    // the solution is sensitive to this value.
    private float generalTolerance = 8;
    
    public void setCostToDiffFromModel() {
        costIsDiff = true;
    }
    
    public void setCostToNumMatchedAndDiffFromModel() {
        costIsNumAndDiff = true;
    }
    
    /**
     * NOT READY FOR USE
     * 
     * @param unmatchedLeftXY
     * @param unmatchedRightXY
     * @param image1CentroidX
     * @param image1CentroidY
     * @param image2CentroidX
     * @param image2CentroidY
     * @param outputMatchedLeftXY
     * @param outputMatchedRightXY
     */
    public TransformationPointFit performPartitionedMatching0(
        PairIntArray unmatchedLeftXY, PairIntArray unmatchedRightXY,
        int image1CentroidX, int image1CentroidY,
        int image2CentroidX, int image2CentroidY,
        PairIntArray outputMatchedLeftXY, PairIntArray outputMatchedRightXY) {
        
        /*
        -- perform whole sets matching
        -- perform vertical partition matching
        -- perform horizontal partition matching
        -- if the vert or horiz results are better that whole sets, consider
           partition into 4 
        
        TODO: methods are using a special fitting function specifically for skyline
        matches, but the function may not be ideal for whole image corner
        matching, so need to allow the choice of fitness function
        to be passed as an argument.
        */
        
        // ====== whole sets match =======
        PairIntArray allPointsLeftMatched = new PairIntArray();
        PairIntArray allPointsRightMatched = new PairIntArray();
        TransformationPointFit allPointsFit = performMatching(
            unmatchedLeftXY, unmatchedRightXY,
            image1CentroidX, image1CentroidY, 
            image2CentroidX, image2CentroidY,
            allPointsLeftMatched, allPointsRightMatched, 1.0f); 
        float allPointsNStat = (float)allPointsFit.getNumberOfMatchedPoints()/
            (float)allPointsFit.getNMaxMatchable();
                
        log.info("all points set nStat=" + allPointsNStat + " Euclidean fit=" + 
            allPointsFit.toString());
        
        
        // ====== vertical partitioned matches =======
        TransformationPointFit verticalPartitionedFit = 
            performVerticalPartitionedMatching(2,
            unmatchedLeftXY, unmatchedRightXY,
            image1CentroidX, image1CentroidY,
            image2CentroidX, image2CentroidY);
        
        /*
        // ====== horizontal partitioned matches =======
        TransformationPointFit horizontalPartitionedFit = 
            performHorizontalPartitionedMatching(2,
            unmatchedLeftXY, unmatchedRightXY,
            image1CentroidX, image1CentroidY,
            image2CentroidX, image2CentroidY);
        */
        
        TransformationPointFit bestFit = allPointsFit;
        
        if (fitIsBetter2(bestFit, verticalPartitionedFit)) {
            //fitIsBetterNStat(bestFit, verticalPartitionedFit)) {
            
            bestFit = verticalPartitionedFit;
        }
        
        if (bestFit == null) {
            return null;
        }
        
        // TODO: compare to horizontal fits when implemented
        
        
        Transformer transformer = new Transformer();
        
        PairFloatArray transformedLeft = transformer.applyTransformation2(
            bestFit.getParameters(), unmatchedLeftXY, image1CentroidX, 
            image1CentroidY);
        
        double tolerance = generalTolerance;

        float[][] matchIndexesAndDiffs = calculateMatchUsingOptimal(
            transformedLeft, unmatchedRightXY, tolerance);

        matchPoints(unmatchedLeftXY, unmatchedRightXY, tolerance, 
            matchIndexesAndDiffs, outputMatchedLeftXY, outputMatchedRightXY);
        
        return bestFit;
    }
    
    /**
     * NOT READY FOR USE
     * 
     * @param numberOfPartitions the number of vertical partitions to make.  
     * the maximum value accepted is 3 and minimum is 1.
     * @param unmatchedLeftXY
     * @param unmatchedRightXY
     * @param image1CentroidX
     * @param image1CentroidY
     * @param image2CentroidX
     * @param image2CentroidY
     * @return best fitting transformation between unmatched left and right 
     */
    public TransformationPointFit performVerticalPartitionedMatching(
        final int numberOfPartitions,
        PairIntArray unmatchedLeftXY, PairIntArray unmatchedRightXY,
        int image1CentroidX, int image1CentroidY,
        int image2CentroidX, int image2CentroidY) {
       
        if (numberOfPartitions > 3) {
            throw new IllegalArgumentException("numberOfPartitions max value is 3");
        }
        if (numberOfPartitions < 1) {
            throw new IllegalArgumentException("numberOfPartitions min value is 1");
        }
        
        int nXDiv = numberOfPartitions;
        int nYDiv = 1;
        float setsFractionOfImage = 0.5f;
        int n = (int)Math.pow(nXDiv * nYDiv, 2);
        
        TransformationPointFit[] vertPartitionedFits = new TransformationPointFit[n];
        float[] nStat = new float[n];
                            
        PointPartitioner partitioner = new PointPartitioner();
        PairIntArray[] vertPartitionedLeft = partitioner.partitionVerticalOnly(
            unmatchedLeftXY, nXDiv);
        PairIntArray[] vertPartitionedRight = partitioner.partitionVerticalOnly(
            unmatchedRightXY, nXDiv);
        
        int bestFitIdx = -1;
        
        int count = 0;
        for (int p1 = 0; p1 < vertPartitionedLeft.length; p1++) {
            for (int p2 = 0; p2 < vertPartitionedRight.length; p2++) {
 //if (!(p1 == 1 && p2 == 0)) { count++; continue;}

                // determine fit only with partitioned points

                PairIntArray part1 = vertPartitionedLeft[p1];
                PairIntArray part2 = vertPartitionedRight[p2];

                PairIntArray part1Matched = new PairIntArray();
                PairIntArray part2Matched = new PairIntArray();
                
                TransformationPointFit fit = performMatching(part1, part2,
                    image1CentroidX, image1CentroidY, 
                    image2CentroidX, image2CentroidY,
                    part1Matched, part2Matched, setsFractionOfImage);
                
                if (fit == null) {
                    count++;
                    continue;
                }
                
                nStat[count] = (float)fit.getNumberOfMatchedPoints()/
                    (float)fit.getNMaxMatchable();
                
                log.info(Integer.toString(count) 
                    + ", p1=" + p1 + " p2=" + p2 + ") nStat=" + nStat[count] +
                    " Euclidean fit=" + fit.toString());
                
                vertPartitionedFits[count] = fit;
                                
                if (bestFitIdx == -1) {
                    bestFitIdx = count;
                } else {
                    // if nStat != inf, best has smallest st dev from mean
                    if (fitIsBetter2(vertPartitionedFits[bestFitIdx], 
                        vertPartitionedFits[count])) {
                        bestFitIdx = count;
                    }
                }
                
                count++;
            }
        }
        
        if (bestFitIdx == -1) {
            return null;
        }
        
        /*
        nDiv=2
        0: [0][1]       1: [0][1]
           [0]             [1]      
        
        2: [0][1]       3: [0][1]     
              [0]             [1]
        
        nDiv=3
        0: [0][1][2]    1: [0][1][2]    2: [0][1][2]
           [0]             [1]             [2]
        
        3: [0][1][2]    4: [0][1][2]    5: [0][1][2]
              [0]             [1]             [2]
        
        6: [0][1][2]    7: [0][1][2]    8: [0][1][2]
                 [0]             [1]             [2]
        
        Consistent solutions for nDiv=2
            would have similar solutions for comparisons 0: and 3:
            (which is index 0 and n+1)
        
        Consistent solutions for nDiv=3
            would have similar solutions for comparisons 0: and 4: and 8:
            (which is index 0 and (n+1) and 2*(n+1))
          OR
            3: and 7:
            (which is index n and (2*n)+1)
        */
        
        TransformationPointFit bestFit = vertPartitionedFits[bestFitIdx];
        
        float scaleTolerance = 0.1f; 
        float rotInDegreesTolerance = 20;
        float translationTolerance = 20;
        
        int idx2 = -1;
        int idx3 = -1;
            
        if (numberOfPartitions == 2) {
            
            if (bestFitIdx == 0) {
                
                // is this similar to vertPartitionedFits[3] ?
                idx2 = 3;
                
            } else if (bestFitIdx == 3) {
                
                // is this similar to vertPartitionedFits[0] ?
                idx2 = 0;
            }
            
        } else if (numberOfPartitions == 3) {
            
            if (bestFitIdx == 0) {
                
                // is this similar to vertPartitionedFits[4] and vertPartitionedFits[8] ?
                idx2 = 4;
                idx3 = 8;
                
            } else if (bestFitIdx == 4) {
                
                // is this similar to vertPartitionedFits[0] and vertPartitionedFits[8] ?
                idx2 = 0;
                idx3 = 8;
                
            } else if (bestFitIdx == 8) {
                
                // is this similar to vertPartitionedFits[0] and vertPartitionedFits[4] ?
                idx2 = 0;
                idx3 = 4;
                
            } else if (bestFitIdx == 3) {
                
                // is this similar to vertPartitionedFits[7] ?
                idx2 = 7;
                idx3 = -1;
                
            } else if (bestFitIdx == 7) {
                
                // is this similar to vertPartitionedFits[3] ?
                idx2 = 3;
                idx3 = -1;
            }
            
            if (idx2 > -1) {
                
                boolean fitIsSimilar = fitIsSimilar(bestFit, 
                    vertPartitionedFits[idx2], scaleTolerance, 
                    rotInDegreesTolerance, translationTolerance);

                if (fitIsSimilar) {

                    log.info("similar solutions for idx=" 
                        + Integer.toString(bestFitIdx) + " and " 
                        + Integer.toString(idx2) + ":" +
                    " Euclidean fit" + Integer.toString(bestFitIdx) 
                        + "=" + bestFit.toString() + 
                    " Euclidean fit" + Integer.toString(idx2) + "=" 
                        + vertPartitionedFits[idx2].toString());

                    //TODO: combine these?
                }

                if (idx3 > -1) {

                    fitIsSimilar = fitIsSimilar(bestFit, 
                        vertPartitionedFits[idx3], scaleTolerance, 
                        rotInDegreesTolerance, translationTolerance);

                    if (fitIsSimilar) {

                        log.info("similar solutions for idx=" 
                            + Integer.toString(bestFitIdx) + " and " 
                            + Integer.toString(idx3) + ":" +
                        " Euclidean fit" + Integer.toString(bestFitIdx) 
                            + "=" + bestFit.toString() + 
                        " Euclidean fit" + Integer.toString(idx3) + "=" 
                            + vertPartitionedFits[idx3].toString());

                        //TODO: combine these?
                    }
                }
            }
        }
        
        return bestFit;
    }
    
    /**
     * NOT READY FOR USE
     * 
     * @param numberOfPartitions the number of vertical partitions to make.  
     * the maximum value accepted is 3.
     * @param unmatchedLeftXY
     * @param unmatchedRightXY
     * @param image1CentroidX
     * @param image1CentroidY
     * @param image2CentroidX
     * @param image2CentroidY
     * @param outputMatchedLeftXY
     * @param outputMatchedRightXY
     * @return best fitting transformation between unmatched points sets
     * left and right
     */
    public TransformationPointFit performVerticalPartitionedMatching(
        final int numberOfPartitions,
        PairIntArray unmatchedLeftXY, PairIntArray unmatchedRightXY,
        int image1CentroidX, int image1CentroidY,
        int image2CentroidX, int image2CentroidY,
        PairIntArray outputMatchedLeftXY, PairIntArray outputMatchedRightXY) {
       
        if (numberOfPartitions > 3) {
            throw new IllegalArgumentException("numberOfPartitions max value is 3");
        }
        if (numberOfPartitions < 1) {
            throw new IllegalArgumentException("numberOfPartitions min value is 1");
        }
        
        TransformationPointFit bestFit = performVerticalPartitionedMatching(
            numberOfPartitions, unmatchedLeftXY, unmatchedRightXY,
            image1CentroidX, image1CentroidY, image2CentroidX, image2CentroidY);
        
        if (bestFit == null) {
            return null;
        }
        
        Transformer transformer = new Transformer();
        
        PairFloatArray transformedLeft = transformer.applyTransformation2(
            bestFit.getParameters(), unmatchedLeftXY, image1CentroidX, 
            image1CentroidY);
        
        double tolerance = generalTolerance;

        float[][] matchIndexesAndDiffs = calculateMatchUsingOptimal(
            transformedLeft, unmatchedRightXY, tolerance);

        matchPoints(unmatchedLeftXY, unmatchedRightXY, tolerance, 
            matchIndexesAndDiffs, outputMatchedLeftXY, outputMatchedRightXY);
        
        return bestFit;
    }
    
    
    /**
     * NOT READY FOR USE
     * 
     * @param unmatchedLeftXY
     * @param unmatchedRightXY
     * @param image1CentroidX
     * @param image1CentroidY
     * @param image2CentroidX
     * @param image2CentroidY
     * @param outputMatchedLeftXY
     * @param outputMatchedRightXY
     */
    public void performPartitionedMatching(
        PairIntArray unmatchedLeftXY, PairIntArray unmatchedRightXY,
        int image1CentroidX, int image1CentroidY,
        int image2CentroidX, int image2CentroidY,
        PairIntArray outputMatchedLeftXY, PairIntArray outputMatchedRightXY) {
        
        /*
        need to match the points, but the sets may require a projection
        transformation rather than a simpler euclidean scale, rotation and
        translation transformation.
        
        trying every combination of point pairs isn't feasible for n > 12:
            n_1!/(2*(n_1 - 2)!) X n_2!/(2*(n_2 - 2)!)
        
        and it ignores the location of points among other points and pixel
        features.        
        
        searching the fundamental matrix after getting an initial calculation
        isn't feasible either.  For example, if used a randomly paired set of 
        points to provide the first estimate of a fundamental matrix only for 
        matching points, and then used
        intervals of factor of 10 from 20 to 1./200, then number of
        permutations is 22^8 = 54875873536.0 and thats
        not including negative factors yet.
        
        Partitioning the image into a few sections and solving the combinations 
        of partitions is an improvement to using the whole image if projection
        is present, and that does not add too much more to the runtime.

        The number of partitions should be kept low, so the partitions formed 
        will be quadrants for each image.  
        That might be changed in the future.  
        
        For each of the partition combinations:
            — finding the best fitting euclidean transformation w/ a cost 
              function that prefers fits with larger number of matches.
              (Note: might make other cost functions available in the future).
        
        To learn which partitions are the correct matches,
        
        */

        int nDiv = 2;
        
        int n = 16;
        
        TransformationPointFit[] bestFits = new TransformationPointFit[n];
        PairIntArray[] matchedPart1 = new PairIntArray[n];
        PairIntArray[] matchedPart2 = new PairIntArray[n];
        
        int count = 0;
                    
        PointPartitioner partitioner = new PointPartitioner();
        PairIntArray[] parts1 = partitioner.partition(unmatchedLeftXY, nDiv);
        PairIntArray[] parts2 = partitioner.partition(unmatchedRightXY, nDiv);

        float setsFractionOfImage = 1.f/(float)(nDiv*nDiv);
        
        for (int p1 = 0; p1 < parts1.length; p1++) {
            for (int p2 = 0; p2 < parts2.length; p2++) {
//if (!((p1 == 1 && p2 == 0) || (p1 == 3 && p2 ==2))) { continue;}

                // determine fit only with partitioned points

                PairIntArray part1 = parts1[p1];
                PairIntArray part2 = parts2[p2];

                PairIntArray part1Matched = new PairIntArray();
                PairIntArray part2Matched = new PairIntArray();
                
                TransformationPointFit  fit = performMatching(part1, part2,
                    image1CentroidX, image1CentroidY, 
                    image2CentroidX, image2CentroidY,
                    part1Matched, part2Matched, setsFractionOfImage);
                
                if (fit == null) {
                    count++;
                    continue;
                }
                
                float nStat = (float)fit.getNumberOfMatchedPoints()/
                    (float)fit.getNMaxMatchable();
                
                log.info(Integer.toString(count) 
                    + ", p1=" + p1 + " p2=" + p2 + ") nStat=" + nStat +
                    " Euclidean fit=" + fit.toString());
                
                bestFits[count] = fit;
                matchedPart1[count] = part1Matched;
                matchedPart2[count] = part2Matched;
                count++;
            }
        }
        
        // ====== examine partition best fits and transformations ====
        
        log.info("examine partition transformations");
        
        //TODO:  this may need alot of changes
        
        // find which quadrant matches represent a consistent intersection:
        // quick look at any having similar transformation parameters:
        double rotTolDegrees = 20;
        double transTol = 20;
        double scaleTol = 0.1;
        
        Set<List<Integer> > similarTransformations = new 
            HashSet<List<Integer> >();
        
        for (int i = 0; i < count; i++) {
            
            if ((bestFits[i] == null) || (bestFits[i].getNMaxMatchable() == 0)) {
                continue;
            }
        
            log.info(String.format(
                "%d) %6f %6.2f %7.1f %7.1f   [%d %d]  %6.2f  %6.2f ", 
                i,
                bestFits[i].getParameters().getScale(),
                bestFits[i].getParameters().getRotationInDegrees(),
                bestFits[i].getParameters().getTranslationX(),
                bestFits[i].getParameters().getTranslationY(),
                bestFits[i].getNumberOfMatchedPoints(),
                bestFits[i].getNMaxMatchable(),
                bestFits[i].getMeanDistFromModel(),
                bestFits[i].getStDevFromMean()
            ));
            
            TransformationParameters params1 = bestFits[i].getParameters();
            
            List<Integer> similar = new ArrayList<Integer>();
            
            for (int j = (i + 1); j < count; j++) {
            
                if ((bestFits[j] == null) || (bestFits[j].getNMaxMatchable() == 0)) {
                    continue;
                }
                
                TransformationParameters params2 = bestFits[j].getParameters();
                
                if ((Math.abs(params1.getScale() - params2.getScale()) <= scaleTol)
                    && (Math.abs(params1.getRotationInDegrees() - params2.getRotationInDegrees()) <= rotTolDegrees)
                    && (Math.abs(params1.getTranslationX() - params2.getTranslationX()) <= transTol)
                    && (Math.abs(params1.getTranslationY() - params2.getTranslationY()) <= transTol)
                    ) {
                    
                    similar.add(Integer.valueOf(j));
                }
            }
            
            if (!similar.isEmpty()) {
                similar.add(Integer.valueOf(i));
                similarTransformations.add(similar);
            }
        }
        
        log.info("similarTransformations.size()=" + similarTransformations.size());
        
        List<Integer> keepIndexes = null;
        if (!similarTransformations.isEmpty()) {
            
            double minAvg = Double.MAX_VALUE;
            List<Integer> minInd = null;
            Iterator<List<Integer> > iter = similarTransformations.iterator();
            while (iter.hasNext()) {
                List<Integer> countIndexes = iter.next();
                double sumAvg = 0;
                for (Integer countIndex : countIndexes) {
                    int cIdx = countIndex.intValue();
                    sumAvg += bestFits[cIdx].getMeanDistFromModel();
                }
                sumAvg /= (double)countIndexes.size();
                if (sumAvg < minAvg) {
                    sumAvg = minAvg;
                    minInd = countIndexes;
                }
            }
            keepIndexes = minInd;
            
        } else if (similarTransformations.size() == 1) {
            
            keepIndexes = similarTransformations.iterator().next();
           
        } else {
        
            /*
            since only one of the bestFits is correct,
            find it as the largest nMatched/nMaxMatchable.
            if there's more than one with same value and it's > 0,
            then 
            (best meanDistFromModel - best stDevMean) > (next best meanDistFromModel)
            else those are considered same 
            and should be re-evaluated against all points.
            In the later case, the new matches should be made for the
            final solution too.
            */
            
            float bestNStatistic = Float.MIN_VALUE;
            
            List<Integer> indexesWithBestNStat = new ArrayList<Integer>();
            
            for (int ii = 0; ii < bestFits.length; ii++) {
                
                if ((bestFits[ii] == null) || 
                    (bestFits[ii].getNMaxMatchable() == 0) ||
                    (bestFits[ii].getNumberOfMatchedPoints() == 0)) {
                    continue;
                }
        
                float stat = (float)bestFits[ii].getNumberOfMatchedPoints()/
                    (float)bestFits[ii].getNMaxMatchable();
                
                if (stat > bestNStatistic) {
                    bestNStatistic = stat;
                    if (!indexesWithBestNStat.isEmpty()) {
                        indexesWithBestNStat.clear();
                    }
                    indexesWithBestNStat.add(Integer.valueOf(ii));
                } else if (Math.abs(stat - bestNStatistic) < 0.001f) {
                    indexesWithBestNStat.add(Integer.valueOf(ii));
                }
            }
            
            if (indexesWithBestNStat.size() == 1) {
                
                // TODO: this is a small number of matched points.  consider
                // re-evaluating the fit for this solution but on all points,
                // and then re-matching the points from it.
                
                keepIndexes = new ArrayList<Integer>(indexesWithBestNStat);
                
            } else {
                
                // compare mean-stdev > next best mean
                
                int bestIdx = -1;
                double bestMD = Double.MIN_VALUE;
                
                for (int ii = 0; ii < indexesWithBestNStat.size(); ii++) {
                    
                    int idx0 = indexesWithBestNStat.get(ii);
                    double md0 = bestFits[idx0].getMeanDistFromModel();
                    double sd0 = bestFits[idx0].getStDevFromMean();
                
                    boolean best = true;
                    for (int j = 0; j < indexesWithBestNStat.size(); j++) {
                        if (ii == j) {
                            continue;
                        }
                        int idx1 = indexesWithBestNStat.get(j);
                        double md1 = bestFits[idx1].getMeanDistFromModel();
                        if ((md0 + sd0) >= md1) {
                            best = false;
                            break;
                        }
                    }
                    if (best) {
                        if ((bestIdx == -1) || ((md0 + sd0) < bestMD)) {
                            bestIdx = idx0;
                            bestMD = md0;
                        }
                    }                                   
                }
                
                if (bestIdx > -1) {
                    
                    // TODO: this is a small number of matched points.  consider
                    // re-evaluating the fit for this solution but on all points,
                    // and then re-matching the points from it.
                    
                    keepIndexes = new ArrayList<Integer>();
                    keepIndexes.add(Integer.valueOf(bestIdx));
                    
                } else {
                    
                    // there was no best.  either evaluate all best against
                    // all points, or consider this no solution.
                    
                    Transformer transformer = new Transformer();
                    
                    TransformationPointFit bestFit2 = null;
                    PairFloatArray bestFit2Transformed = null;
                    
                    for (int ii = 0; ii < bestFits.length; ii++) {
                        if ((bestFits[ii] == null) || 
                            (bestFits[ii].getNMaxMatchable() == 0)) {
                            continue;
                        }                        

                        PairFloatArray transformed = 
                            transformer.applyTransformation2(
                            bestFits[ii].getParameters(), unmatchedLeftXY, 
                            image1CentroidX, image1CentroidY);
                                
                        TransformationPointFit fit =
                            evaluateFitForUnMatchedTransformedOptimal(
                            bestFits[ii].getParameters(), 
                            transformed, unmatchedRightXY, 
                            generalTolerance, generalTolerance);
                        
                        if (fitIsBetter(bestFit2, fit)) {
                            bestFit2 = fit;
                            bestFit2Transformed = transformed;
                        }
                    }
                    
                    if (bestFit2 != null) {
                    
                        double tolerance = generalTolerance;
                        
                        float[][] matchIndexesAndDiffs = 
                            calculateMatchUsingOptimal(bestFit2Transformed, 
                            unmatchedRightXY, tolerance);
                
                        matchPoints(unmatchedLeftXY, unmatchedRightXY, 
                            tolerance, matchIndexesAndDiffs, 
                            outputMatchedLeftXY, outputMatchedRightXY);
                    
                        log.info("final best fit=" + bestFit2Transformed.toString());
                        
                    } else {
                        keepIndexes = null;
                    }
                }
            }
        }
        
        if (keepIndexes != null) {
             // put the matched points into the output list
            for (Integer countIndex : keepIndexes) {
 
                log.info("--> using solution from count=" + countIndex.toString()
                    + " :" + bestFits[countIndex.intValue()]);               
                
                int cIdx = countIndex.intValue();
                PairIntArray set1 = matchedPart1[cIdx];
                PairIntArray set2 = matchedPart2[cIdx];
                for (int ii = 0; ii < set1.getN(); ii++) {
                    outputMatchedLeftXY.add(set1.getX(ii), set1.getY(ii));
                }
                for (int ii = 0; ii < set2.getN(); ii++) {
                    outputMatchedRightXY.add(set2.getX(ii), set2.getY(ii));
                }
            }
        }
    }
    
    /**
     * NOT READY FOR USE
     * 
     * @param unmatchedLeftXY
     * @param unmatchedRightXY
     * @param image1CentroidX
     * @param image1CentroidY
     * @param image2CentroidX
     * @param image2CentroidY
     * @param outputMatchedLeftXY
     * @param outputMatchedRightXY
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * 
     * @return 
     */
    public TransformationPointFit performMatching(
        PairIntArray unmatchedLeftXY, PairIntArray unmatchedRightXY,
        int image1CentroidX, int image1CentroidY,
        int image2CentroidX, int image2CentroidY,
        PairIntArray outputMatchedLeftXY, PairIntArray outputMatchedRightXY,
        float setsFractionOfImage) {
        
        Transformer transformer = new Transformer();
        
        PairIntArray part1 = unmatchedLeftXY;
        PairIntArray part2 = unmatchedRightXY;

        TransformationPointFit transFit = 
            calculateProjectiveTransformationWrapper(
            part1, part2, 
            image1CentroidX, image1CentroidY,
            image2CentroidX, image2CentroidY, setsFractionOfImage);

        if (transFit == null) {
            return null;
        }
              
        // -- filter part1 and part2 to keep only intersection region

        PairIntArray filtered1 = new PairIntArray();
        PairIntArray filtered2 = new PairIntArray();
        populateIntersectionOfRegions(transFit.getParameters(),
            part1, part2, image1CentroidX, image1CentroidY,
            image2CentroidX, image2CentroidY,
            filtered1, filtered2);

        int nMaxMatchable = (filtered1.getN() < filtered2.getN()) ?
            filtered1.getN() : filtered2.getN();
                
        transFit.setMaximumNumberMatchable(nMaxMatchable);
        
        if (nMaxMatchable == 0) {
            return transFit;
        }
                
        // -- transform filtered1 for matching and evaluation

        PairFloatArray transformedFiltered1 = 
            transformer.applyTransformation2(transFit.getParameters(), 
            filtered1, image1CentroidX, image1CentroidY);

        double tolerance = generalTolerance;

        //evaluate the fit and store a statistical var: nmatched/nmaxmatchable

        float[][] matchIndexesAndDiffs = calculateMatchUsingOptimal(
            transformedFiltered1, filtered2, tolerance);

        PairFloatArray part1MatchedTransformed = new PairFloatArray();
        PairIntArray part2Matched = new PairIntArray();
                
        matchPoints(transformedFiltered1, filtered2, tolerance, 
            matchIndexesAndDiffs, part1MatchedTransformed, part2Matched);

        PairIntArray part1Matched = new PairIntArray();
        part2Matched = new PairIntArray();

        matchPoints(filtered1, filtered2, tolerance, matchIndexesAndDiffs, 
            part1Matched, part2Matched);

        TransformationPointFit fit2 = evaluateFitForMatchedTransformed(
            transFit.getParameters(), part1MatchedTransformed, part2Matched);
        
        fit2.setTolerance(tolerance);
                          
        fit2.setMaximumNumberMatchable(nMaxMatchable);

        float nStat = (nMaxMatchable > 0) ? 
            (float)part1Matched.getN()/(float)nMaxMatchable : 0;

        log.info("nStat=" + nStat + " nMaxMatchable=" + nMaxMatchable +
            " Euclidean fit=" + fit2.toString()
            + " original fit=" + transFit.toString());
           
        outputMatchedLeftXY.swapContents(part1Matched);
        
        outputMatchedRightXY.swapContents(part2Matched);        
        
        return fit2;
    }
 
    /**
     * 
     * calculate a projective transformation that results in the best match
     * between the unmatched points in set1 and set2.
     * Note that output1 and output2 have mismatches in them due to projection
     * effects, so one should probably follow this with the use of RANSAC
     * to further remove outliers.
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * 
     * NOTE: this method is not ready for use yet.
     * 
     * @param set1 set of points such as corners from an image
     * @param set2 set of points such as corners from another image.
     * @param image1CentroidX
     * @param image1CentroidY
     * @param output1
     * @param output2
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * 
     * @return 
     */
    public StereoProjectionTransformerFit 
        matchPointsUsingProjectiveTransformation(PairIntArray set1, 
        PairIntArray set2, int image1CentroidX, int image1CentroidY,
        PairIntArray output1, PairIntArray output2, 
        float setsFractionOfImage) {
                
        //TODO:  this method is changing.  needs alot of improvement!
            
        int rotDelta = 10;
        int rotStart = 0;
        int rotStop = 360; 
        int scaleStart = 1;
        int scaleStop = 10;
        int scaleDelta = 1;
        
        double tolTransX = image1CentroidX * 0.02f;
        double tolTransY = image1CentroidY * 0.02f;
        if (tolTransX < minTolerance) {
            tolTransX = minTolerance;
        }
        if (tolTransY < minTolerance) {
            tolTransY = minTolerance;
        }
        
        //TODO: improve this estimate because the result is sensitive to it
        double tolerance = Math.sqrt(tolTransX*tolTransX + tolTransY*tolTransY);
        tolerance *= 2;
        
        int convergence = (set1.getN() < set2.getN()) ? set1.getN() 
            : set2.getN();
      
        StereoProjectionTransformerFit bestFit = null;
        PairIntArray bestFitMatches1 = new PairIntArray();
        PairIntArray bestFitMatches2 = new PairIntArray();
        TransformationParameters bestParams = null;
        
        //TODO: once this is correct, find a faster approach to the local search area
        // for example, if bestFit for a scale worsens as scale increases,
        // one should be able to break from the scale
        
        StereoProjectionTransformerFit bestFitForScale = null;
        
        for (int scale = scaleStart; scale < scaleStop; scale += scaleDelta) {
                        
            for (int rot = rotStart; rot < rotStop; rot += rotDelta) {
            
                // requires estimation of translation x and y
                
                TransformationParameters params = 
                    calculateTranslationForUnmatched(set1, set2, 
                        rot*Math.PI/180., scale, 
                        image1CentroidX, image1CentroidY, setsFractionOfImage);
               
                // the optimal match is necessary
                float[][] matchIndexesAndDiffs = calculateMatchUsingOptimal(
                    set1, set2, 
                    params.getScale(), params.getRotationInRadians(), 
                    params.getTranslationX(), params.getTranslationY(),
                    tolerance, image1CentroidX, image1CentroidY);
                
                /*float[][] matchIndexesAndDiffs = calculateMatchUsingGreedy(
                    set1, set2, 
                    params.getScale(), params.getRotationInRadians(), 
                    params.getTranslationX(), params.getTranslationY(), 
                    tolerance, image1Width >> 1, image1Height >> 1);
                */
                
                if ((matchIndexesAndDiffs == null) || 
                    (matchIndexesAndDiffs.length < 8)) {
                    continue;
                }
                
                PairFloatArray subset1 = new PairFloatArray();
                PairFloatArray subset2 = new PairFloatArray();
                LinkedHashSet<Integer> matchedSet1 = new LinkedHashSet<Integer>();
                LinkedHashSet<Integer> matchedSet2 = new LinkedHashSet<Integer>();
                for (int ii = 0; ii < matchIndexesAndDiffs.length; ii++) {
                    int idx1 = (int)matchIndexesAndDiffs[ii][0];
                    int idx2 = (int)matchIndexesAndDiffs[ii][1];
                    subset1.add(set1.getX(idx1), set1.getY(idx1));
                    subset2.add(set2.getX(idx2), set2.getY(idx2));
                    matchedSet1.add(Integer.valueOf(idx1));
                    matchedSet2.add(Integer.valueOf(idx2));
                }
               
                PairFloatArray unmatched1 = new PairFloatArray();
                for (int ii = 0; ii < set1.getN(); ii++) {
                    if (matchedSet1.contains(Integer.valueOf(ii))) {
                        continue;
                    }
                    unmatched1.add(set1.getX(ii), set1.getY(ii));
                }
                PairFloatArray unmatched2 = new PairFloatArray();
                for (int ii = 0; ii < set2.getN(); ii++) {
                    if (matchedSet2.contains(Integer.valueOf(ii))) {
                        continue;
                    }
                    unmatched2.add(set2.getX(ii), set2.getY(ii));
                }
                
                // evaluate the fit of the epipolar solution with all points
                // and keep the best solution.
                StereoProjectionTransformer sTransformer = 
                    new StereoProjectionTransformer();
                
                SimpleMatrix fm = 
                    sTransformer.calculateEpipolarProjectionForPerfectlyMatched(
                    subset1, subset2);
                  
                SimpleMatrix leftPoints = 
                    StereoProjectionTransformer.rewriteInto3ColumnMatrix(unmatched1);
                
                SimpleMatrix rightPoints = 
                    StereoProjectionTransformer.rewriteInto3ColumnMatrix(unmatched2);
        
                // x is index from leftPoints.  y is index from rightPoints.
                PairIntArray matchedIndexes = new PairIntArray();
            
//TODO: problem here w/ eval using only perp dist to epipolar
                /*
                evaluation of the fit using only the perpendicular distance does
                not produce the best matches because it does not use the depth
                into the image to match the location of the point in image2
                parallel to the epipolar line and not using that means that
                the distance from the predicted location is only partially known.
                False matches happen often because of it.
                
                So adopting the strategy of: rough euclidean match to get an
                initial set of matched points, then calculating a fundamental
                matrix from those.  Then using the fm to match the remaining
                unmatched.
                --> if the later does not work as well as needed, will try to
                change the method of matching the remaining unmatched to use
                a proximity to the matched in both images if within an
                association radius.
                
                */
                
                // have matches for subset1, subset2.
                // match the remaining unmatched1, unmatched2.
                // first try only matching the remaining unmatched.
                // if the results are not very good, consider
                // a method to determine best match for remaining by 
                // proximity to a known match in both images.
                
                StereoProjectionTransformerFit fit = 
                    sTransformer.evaluateFitForUnmatched(fm, 
                    leftPoints, rightPoints, 1.0/*tolerance*/, matchedIndexes);
                
                PairIntArray finalMatched1 = new PairIntArray();
                PairIntArray finalMatched2 = new PairIntArray();
                for (int ii = 0; ii < subset1.getN(); ii++) {
                    finalMatched1.add(
                        Math.round(subset1.getX(ii)), 
                        Math.round(subset1.getY(ii)));
                    finalMatched2.add(
                        Math.round(subset2.getX(ii)), 
                        Math.round(subset2.getY(ii)));
                }
                for (int ii = 0; ii < matchedIndexes.getN(); ii++) {
                    int idx1 = matchedIndexes.getX(ii);
                    int idx2 = matchedIndexes.getY(ii);
                    finalMatched1.add(leftPoints.getIndex(0, idx1),
                        leftPoints.getIndex(1, idx1));
                    finalMatched2.add(rightPoints.getIndex(0, idx2),
                        rightPoints.getIndex(1, idx2));
                }
                
                leftPoints = 
                    StereoProjectionTransformer.rewriteInto3ColumnMatrix(
                    finalMatched1);
                
                rightPoints = 
                    StereoProjectionTransformer.rewriteInto3ColumnMatrix(
                    finalMatched2);
                
                fit = sTransformer.evaluateFitForAlreadyMatched(fm, leftPoints, rightPoints, 
                    tolerance);
                
                if (spFitIsBetter(bestFit, fit)) {
                    bestFit = fit;
                    bestParams = params;                    
                    bestFitMatches1 = finalMatched1;
                    bestFitMatches2 = finalMatched2;
                }
            }
            
            if (spFitIsBetter(bestFitForScale, bestFit)) {
                bestFitForScale = bestFit;
            } else {
                break;
            }
        } 
        
        if (bestFitMatches1 != null) {
            System.out.println(bestParams.toString());
            System.out.println(bestFit.toString());
            output1.swapContents(bestFitMatches1);
            output2.swapContents(bestFitMatches2);
        }
         
        return bestFit;
    }
    
    boolean spFitIsBetter(StereoProjectionTransformerFit bestFit, 
        StereoProjectionTransformerFit compareFit) {
        
        if (compareFit == null) {
            return false;
        }
        if (bestFit == null) {
            return true;
        }
       
        long nMatches = compareFit.getNMatches();
        
        if (nMatches > bestFit.getNMatches()) {
            
            return true;
            
        } else if (nMatches == bestFit.getNMatches()) {
            
            if (!Double.isNaN(compareFit.getMeanDistance()) && (
                compareFit.getMeanDistance()
                < bestFit.getMeanDistance())) {
                
                return true;
                
            } else if (compareFit.getMeanDistance()
                == bestFit.getMeanDistance()) {
                
                if (compareFit.getStDevFromMean()< bestFit.getStDevFromMean()) {
                    
                    return true;
                }
            }
        }
        
        return false;
    }
    
    /**
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
     * @param image1CentroidX
     * @param image1CentroidY
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * 
     * @return 
     */
    public TransformationPointFit calculateRoughTransformationForUnmatched(
        PairIntArray set1, PairIntArray set2, int image1CentroidX, 
        int image1CentroidY, float setsFractionOfImage) {
        
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
                        edge1, edge2, image1CentroidX, image1CentroidY,
                        0, 360, rotDelta,
                        1, scaleStop, 1, setsAreMatched, setsFractionOfImage);

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
                image1CentroidX, image1CentroidY);

            tolerance = findToleranceForOutliers(
                bestFit, matchIndexesAndDiffs,
                image1CentroidX, image1CentroidY);

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
            bestFit.getParameters(), input1, input2, 
            image1CentroidX, image1CentroidY, 
            setsAreMatched, setsFractionOfImage);
        
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
     * @param image1CentroidX
     * @param image1CentroidY
     * @param outputMatched1
     * @param outputMatched2
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * 
     * @return 
     */
    public TransformationPointFit calculateTransformationAndMatch(
        PairIntArray set1, PairIntArray set2, 
        int image1CentroidX, int image1CentroidY,
        PairIntArray outputMatched1, PairIntArray outputMatched2,
        float setsFractionOfImage) {
        
        /*
        attempts rough euclidean transformation and uses that to match
        points in set1 and set2, then removes outliers from the matched set.
        */
        
        boolean setsAreMatched = false;
        
        int rotDelta = 90;
        
        TransformationPointFit fit = calculateTransformationWithGridSearch(
            set1, set2, image1CentroidX, image1CentroidY,
            0, 360, rotDelta,
            1, 11, 2, setsAreMatched, setsFractionOfImage);
        
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
                    image1CentroidX, image1CentroidY);
                    
                tolerance = findToleranceForOutliers(
                    fit, matchIndexesAndDiffs,
                    image1CentroidX, image1CentroidY);
                
                matchPoints(input1, input2, tolerance, matchIndexesAndDiffs,
                    outputMatched1, outputMatched2);
                
                input1 = outputMatched1.copy();
                input2 = outputMatched2.copy();
                                
                // rematch remaining points
                float[][] matchIndexesAndDiffs2 = calculateMatchUsingOptimal(
                    input1, input2, fit, tolerance,
                    image1CentroidX, image1CentroidY);
                
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
                input1, input2, image1CentroidX, image1CentroidY,
                rotStart, rotStop, 10,
                scaleStart, scale + 1, 1, 
                setsAreMatched, setsFractionOfImage);

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
                input1, input2, image1CentroidX, image1CentroidY,
                setsAreMatched, setsFractionOfImage);
        
            if (!setsAreMatched) {
                
                // ==== match points and exclude the largest outliers ====
                
                int nBefore = input1.getN();
                
                double tolerance = Double.MAX_VALUE;
                
                float[][] matchIndexesAndDiffs = calculateMatchUsingOptimal(
                    input1, input2, fit, tolerance,
                    image1CentroidX, image1CentroidY);
                    
                tolerance = findToleranceForOutliers(
                    fit, matchIndexesAndDiffs,
                    image1CentroidX, image1CentroidY);
                
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

    public boolean fitIsBetterNStat(TransformationPointFit bestFit, 
        TransformationPointFit compareFit) {
        
        if (compareFit == null) {
            return false;
        }
        if (bestFit == null) {
            return true;
        }
                
        float bestNStat = (float)bestFit.getNumberOfMatchedPoints()/
            (float)bestFit.getNMaxMatchable();
        
        float compareNStat = (float)compareFit.getNumberOfMatchedPoints()/
            (float)compareFit.getNMaxMatchable();
        
        if (Float.isInfinite(compareNStat)) {
            return false;
        }
        
        // a larger nStat is a better fit
        if (Float.compare(bestNStat, compareNStat) > 0) {
            return true;
        }
        
        return false;
    }
  
    public boolean fitIsBetter(TransformationPointFit bestFit, 
        TransformationPointFit compareFit) {
        
        if (costIsDiff) {
            return fitIsBetterPreferDiff(bestFit, compareFit);
        } else if (costIsNumAndDiff) {
            return fitIsBetterUseNumAndDiff(bestFit, compareFit);
        }
        
        if (compareFit == null) {
            return false;
        }
        if (bestFit == null) {
            return true;
        }
        
        int nMatches = compareFit.getNumberOfMatchedPoints();
        
        if (nMatches > bestFit.getNumberOfMatchedPoints()) {
            
            return true;
            
        } else if (nMatches == bestFit.getNumberOfMatchedPoints()) {
            
            // essentially, comparing avg + std dev from avg
            
            if (!Double.isNaN(compareFit.getMeanDistFromModel())) {
                
                double compAvg = compareFit.getMeanDistFromModel();
                double compAvgS = compAvg + compareFit.getStDevFromMean();
                
                double bestAvg = bestFit.getMeanDistFromModel();
                double bestAvgS = bestAvg + bestFit.getStDevFromMean();
                
                if ((compAvg <= bestAvg) && (compAvgS < bestAvgS)) {
                    return true;
                }
            }
        }
        
        return false;
    }

    boolean fitIsBetterPreferDiff(TransformationPointFit bestFit, 
        TransformationPointFit compareFit) {
        
        if (compareFit == null) {
            return false;
        }
        if (bestFit == null) {
            return true;
        }
        
        double compAvg = compareFit.getMeanDistFromModel();
        double compAvgS = compAvg + compareFit.getStDevFromMean();

        double bestAvg = bestFit.getMeanDistFromModel();
        double bestAvgS = bestAvg + bestFit.getStDevFromMean();
        
        if (((compAvg/bestAvg) > 1) && (compAvg > 5)) {
        
            if (compAvg < bestAvg) {
                return true;
            } else if ((compAvg <= bestAvg) && (compAvgS < bestAvgS)) {
                return true;
            }
            
            return false;
            
        } else {
            
            int nMatches = compareFit.getNumberOfMatchedPoints();

            if (nMatches > bestFit.getNumberOfMatchedPoints()) {
                return true;
            } else if (nMatches == bestFit.getNumberOfMatchedPoints()) {
                
                if (!Double.isNaN(compareFit.getMeanDistFromModel())) {
                    if ((compAvg <= bestAvg) && (compAvgS < bestAvgS)) {
                        return true;
                    }
                }
            }
            
            return false;
        }
    }
    
    public boolean fitIsBetterUseNumAndDiff(TransformationPointFit bestFit, 
        TransformationPointFit compareFit) {

        if (compareFit == null) {
            return false;
        }
        if (bestFit == null) {
            if (compareFit.getNumberOfMatchedPoints() > 0) {
                return true;
            } else {
                return false;
            }
        }

        float compN = compareFit.getNumberOfMatchedPoints();
        float bestN = bestFit.getNumberOfMatchedPoints();
        
        double compAvg = compareFit.getMeanDistFromModel();
        double compS = compareFit.getStDevFromMean();
        double compAvgS = compAvg + compS;

        double bestAvg = bestFit.getMeanDistFromModel();
        double bestS = bestFit.getStDevFromMean();
        double bestAvgS = bestAvg + bestS;
        
        float f = 1.5f;
        
        if (!Double.isNaN(compAvg)) {
            if ((compN/bestN) >= f) {
                return true;
            } else if ((compN >= bestN) && (compAvg < bestAvg) && (compS < bestS)) {
                return true;
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
     * @param image1CentroidX
     * @param image1CentroidY
     * @param rotStart start of rotation search in degrees
     * @param rotStop stop (exclusive) or rotations search in degrees
     * @param rotDelta change in rotation to add to reach next step in rotation 
     *     search in degrees
     * @param scaleStart
     * @param scaleStop
     * @param scaleDelta
     * @param setsAreMatched
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * 
     * @return 
     */
    public TransformationPointFit calculateTransformationWithGridSearch(
        PairIntArray set1, PairIntArray set2, 
        int image1CentroidX, int image1CentroidY,
        int rotStart, int rotStop, int rotDelta,
        int scaleStart, int scaleStop, int scaleDelta,
        boolean setsAreMatched, float setsFractionOfImage) {
        
        double tolTransX = generalTolerance;//4.0f * image1CentroidX * 0.02f;
        double tolTransY = generalTolerance;//4.0f * image1CentroidY * 0.02f;
        if (tolTransX < minTolerance) {
            tolTransX = minTolerance;
        }
        if (tolTransY < minTolerance) {
            tolTransY = minTolerance;
        }
        
        int[] largeScaleRotRefines = new int[]{
            rotStart, 
            (int)(rotStart + 0.25*(rotStop - rotStart)),
            (int)(rotStart + 0.5*(rotStop - rotStart)),
            (int)(rotStart + 0.75*(rotStop - rotStart))
        };
        
        Transformer transformer = new Transformer();
        
        int nMaxMatchable = (set1.getN() < set2.getN()) ? set1.getN() 
            : set2.getN();
        
        int convergence = nMaxMatchable;
        
        TransformationPointFit bestFit = null;
                
        TransformationPointFit bestFitForScale = null;
        
        for (int scale = scaleStart; scale <= scaleStop; scale += scaleDelta) {
            for (int rot = rotStart; rot < rotStop; rot += rotDelta) {
            
                TransformationParameters params;
                
                if (setsAreMatched) {
                    params = calculateTranslationForMatched(set1, set2, 
                        rot*Math.PI/180., scale, image1CentroidX, image1CentroidY);                    
                } else {
                    //TODO: when refactor, this method already determines
                    //   fit for one branch so reduce redundant code
                    params = calculateTranslationForUnmatched(set1, set2, 
                        rot*Math.PI/180., scale, image1CentroidX, 
                        image1CentroidY, setsFractionOfImage);
                }
                
                PairFloatArray allPoints1Tr = transformer.applyTransformation(
                    params, image1CentroidX, image1CentroidY, set1);
                
                TransformationPointFit fit;
                
                if (setsAreMatched) {
                    fit = evaluateFitForMatchedTransformed(params, 
                        allPoints1Tr, set2);                    
                } else {
        
                    fit = evaluateFitForUnMatchedTransformedGreedy(params,
                    //fit = evaluateFitForUnMatchedTransformedOptimal(params,
                        allPoints1Tr, set2, tolTransX, tolTransY);
                    
                    //TODO: rerun tests. is this still necessary:
                    
                    if ((scale >= 4) && ((fit == null) || 
                        (((float)fit.getNumberOfMatchedPoints()
                        /(float)nMaxMatchable)) < 0.1*nMaxMatchable) ) {
                        
                        // need to limit the use of this...
                        int idx = Arrays.binarySearch(largeScaleRotRefines,
                            rot);
                        
                        if (idx > -1) {
                        
                            float tx = params.getTranslationX();
                            float ty = params.getTranslationY();

                            // refine the translation and try fit again
                            refineTranslation(params, set1, set2, 
                                image1CentroidX, image1CentroidY);

                            if ((Math.abs(params.getTranslationX() - tx) > 1) ||
                                Math.abs(params.getTranslationY() - ty) > 1) {

                                allPoints1Tr = transformer.applyTransformation(
                                    params, image1CentroidX, image1CentroidY, 
                                    set1);

                                fit = evaluateFitForUnmatchedTransformed(params, 
                                    allPoints1Tr, set2, tolTransX, tolTransY);
                            }
                        }
                    }
                }
          
                if (fitIsBetter(bestFit, fit)) {
                    
log.info("==> " + " tx=" + fit.getTranslationX() + " ty=" + fit.getTranslationY()
+ " rot=" + rot + " scale=" + scale + "\nfit=" + fit.toString());

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
            
            if (fitIsBetter(bestFitForScale, bestFit)) {
               
                bestFitForScale = bestFit;
                
            } else {
                // scale was probably smaller so return best solution
                break;
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
        
        // ==== refine translation and evalFit ====
        
        //TODO: this shouldn't be necessary now as translation is refined above
        //TODO: remove...
        
        if (false) {
        //if (bestFit != null) {
                        
            refineTranslation(bestFit.getParameters(), set1, set2, 
                image1CentroidX, image1CentroidY);
            
            PairFloatArray allPoints1Tr = transformer.applyTransformation(
                bestFit.getParameters(), image1CentroidX, image1CentroidY, 
                set1);
            
            if (setsAreMatched) {
                bestFit = evaluateFitForMatchedTransformed(
                    bestFit.getParameters(), allPoints1Tr, set2);                    
            } else {
                bestFit = evaluateFitForUnmatchedTransformed(
                    bestFit.getParameters(), allPoints1Tr, set2, 
                    tolTransX, tolTransY);
            }
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
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * @return 
     */
    public TransformationPointFit calculateTranslation(PairIntArray set1, 
        PairIntArray set2, double rotation, 
        double scale, int centroidX1, int centroidY1, boolean setsAreMatched,
        float setsFractionOfImage) {
        
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
        
        float tolTransX = generalTolerance;//4.f * centroidX1 * 0.02f;
        float tolTransY = generalTolerance;//4.f * centroidY1 * 0.02f;        
        /*if (tolTransX < minTolerance) {
            tolTransX = minTolerance;
        }
        if (tolTransY < minTolerance) {
            tolTransY = minTolerance;
        }*/
        
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
            
            int nBins = (int)(setsFractionOfImage*(float)centroidX1*4.f/30.f);
            
            HistogramHolder hX = Histogram
                .createSimpleHistogram(nBins,
                //.defaultHistogramCreator(
                transX, Errors.populateYErrorsBySqrt(transX));
        
            nBins = (int)(setsFractionOfImage*(float)centroidY1*4.f/30.f);
            
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

            float[] transXY = determineTranslationFromHistograms(hX, hY,
                xr, yr, set2, tolTransX, tolTransY);
            
            peakTransX = transXY[0];
            peakTransY = transXY[1];
            
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
        if (tolTransX < minTolerance) {
            tolTransX = minTolerance;
        }
        if (tolTransY < minTolerance) {
            tolTransY = minTolerance;
        }
        
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
     * 
     * @param set1 two dimensional array holding x and y points with
     * first dimension being the point number and the 2nd being x and y
     * for example set1[row][0] is x and set1[row][1] is y for point in row.
     * @param set2 two dimensional array holding x and y points with
     * first dimension being the point number and the 2nd being x and y
     * for example set2[row][0] is x and set2[row][1] is y for point in row.
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * @return 
     */
    public float[] calculateTranslation(double[][] set1, 
        double[][] set2, int centroidX1, int centroidY1,
        float setsFractionOfImage) {
        
        if (set1 == null) {
            throw new IllegalArgumentException("set1 cannot be null");
        }
        if (set2 == null) {
            throw new IllegalArgumentException("set2 cannot be null");
        }
        
        int nTrans = set1.length * set2.length;
        int count = 0;
        float[] transX = new float[nTrans];
        float[] transY = new float[nTrans];
                
        for (int i = 0; i < set1.length; i++) {
            
            double x = set1[i][0];
            double y = set1[i][1];
            
            for (int j = 0; j < set2.length; j++) {
                
                double x2 = set2[j][0];
                double y2 = set2[j][1];
                
                transX[count] = (float)(x2 - x);
                transY[count] = (float)(y2 - y);
                
                count++;
            }
        }
        
        float peakTransX;
        float peakTransY;
        
        // when there aren't enough points for useful histogram,
        // will make a frequency map of round to integer,
        // and take the peak if its larger than next peak,
        // else, take the average of largest frequencies.
        
        if ((set1.length > 5) && (set2.length > 5)) {
            
            /*LinearRegression lr = new LinearRegression();
            //lr.plotTheLinearRegression(transX, transY);
            float[] xyPeaks = lr.calculateTheilSenEstimatorMedian(transX, transY);
            
            peakTransX = xyPeaks[0];
            peakTransY = xyPeaks[1];
            log.info("thiel sen estimator: " + peakTransX + ", " + peakTransY);
            */
            
            int nBins = (int)(setsFractionOfImage*(float)centroidX1*4.f/30.f);
            
            HistogramHolder hX = Histogram
                .createSimpleHistogram(nBins,
                //.defaultHistogramCreator(
                transX, Errors.populateYErrorsBySqrt(transX));
        
            nBins = (int)(setsFractionOfImage*(float)centroidY1*4.f/30.f);
            
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
            
            log.info("histogram: " + peakTransX + ", " + peakTransY);
            
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
        
        return new float[]{(int)peakTransX, (int)peakTransY};
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
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * @return 
     */
    public TransformationParameters calculateTranslationForUnmatched(
        PairIntArray set1, PairIntArray set2, double rotation, double scale, 
        int centroidX1, int centroidY1, float setsFractionOfImage) {
        
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
        
        float peakTransX = Float.MAX_VALUE;
        float peakTransY = Float.MAX_VALUE;
        
        // when there aren't enough points for useful histogram,
        // will make a frequency map of round to integer,
        // and take the peak if its larger than next peak,
        // else, take the average of largest frequencies.
        
        if ((set1.getN() > 5) && (set2.getN() > 5)) {
          
            /*
            if there's only one strong peak, return the FWHM center,
            else, evaluate the fit for the translations of each of the
            strong peaks.
            */
            
            int nBins = (int)(setsFractionOfImage*(float)centroidX1*4.f/30.f);
            if (nBins < 15) {
                nBins = 15;
            }

            HistogramHolder hX = Histogram
                .createSimpleHistogram(nBins,
                transX, Errors.populateYErrorsBySqrt(transX));

            nBins = (int)(setsFractionOfImage*(float)centroidY1*4.f/30.f);
            if (nBins < 15) {
                nBins = 15;
            }
            
            HistogramHolder hY = Histogram
                .createSimpleHistogram(nBins,
                transY, Errors.populateYErrorsBySqrt(transY));

            try {
                hX.plotHistogram("transX", 1);
                hY.plotHistogram("transY", 2);
            } catch (IOException e) {
                log.severe(e.getMessage());
            }
            
            float tolTransX = generalTolerance;//4.f * centroidX1 * 0.02f;
            float tolTransY = generalTolerance;//4.f * centroidY1 * 0.02f;
            if (tolTransX < minTolerance) {
                tolTransX = minTolerance;
            }
            if (tolTransY < minTolerance) {
                tolTransY = minTolerance;
            }
            float[] transXY = determineTranslationFromHistograms(hX, hY,
                xr, yr, set2, tolTransX, tolTransY);
            
            peakTransX = transXY[0];
            peakTransY = transXY[1];
           
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
     * @param image1CentroidX
     * @param image1CentroidY
     * @param setsAreMatched
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * @return 
     */
    public TransformationPointFit refineTransformationWithDownhillSimplex(
        TransformationParameters params,
        PairIntArray set1, PairIntArray set2, 
        int image1CentroidX, int image1CentroidY, boolean setsAreMatched,
        float setsFractionOfImage) {
        
        // projection effects from different camera positions or orientations
        // are not calculated in this point matcher, 
        // but one can set a
        // tolerance in translation to try to allow for a small amount of it.
        // 
        // as an example: rotation of 13 degrees at top of image, 
        // and 0 at bottom from turning the camera slightly during panoramic
        // photos can result in delta translation of 0.02 per pix in X and 
        // 0.016 per pix in Y
        
        //TODO: revise how these are determined
        double transXTol = generalTolerance;//4. * image1CentroidX * 0.02;
        
        double transYTol = generalTolerance;//4 * image1CentroidY * 0.02;
        
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
            //TODO: reconsider this:
             drs = new double[]{0};
        }

        double rMin = r - (10 * Math.PI/180);
        double rMax = r + (10 * Math.PI/180);
        double sMin = s - 0.5;
        double sMax = s + 0.5;
        
        double[] dss = new double[] {
            0.1, 0.2, 0.4, 0.6, 0.8
            //-1.0, -0.1, -0.05, 0.05
        };
        /*if (s == 1) {
            dss = new double[]{0};
            sMin = 1;
        }*/
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
                    rotation, scale, image1CentroidX, image1CentroidY, 
                    setsAreMatched, setsFractionOfImage);
                
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
            fits, image1CentroidX, image1CentroidY, transXTol, transYTol,
            r, s, rMin, rMax, sMin, sMax,
            reflectionCoeff, expansionCoeff, contractionCoeff, reductionCoeff,
            setsAreMatched, setsFractionOfImage
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
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * 
     * @return 
     */
    private TransformationPointFit fitWithDownhillSimplex(
        PairIntArray set1, PairIntArray set2, 
        TransformationPointFit[] fits, int centroidX1, int centroidY1,
        double transXTol, double transYTol,
        double r, double s,
        double rMin, double rMax, double sMin, double sMax,
        float alpha, float gamma, float beta, float tau, boolean setsAreMatched,
        float setsFractionOfImage) {
        
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
        
        // TODO: the correction for scale < 1 during "reflection"
        //   can curtail the search, so consider
        //   a replacement action instead of reflection for that case
        
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
                    centroidX1, centroidY1, setsAreMatched, setsFractionOfImage);
            
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
                            centroidX1, centroidY1, setsAreMatched, 
                            setsFractionOfImage);
                    
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
                            centroidX1, centroidY1, setsAreMatched, 
                            setsFractionOfImage);
                
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
                                setsAreMatched, setsFractionOfImage);
                            
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
        if (tolTransX < minTolerance) {
            tolTransX = minTolerance;
        }
        if (tolTransY < minTolerance) {
            tolTransY = minTolerance;
        }
        
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
     * model generation by point in set1} using an algorithm to optimally match
     * them.  The current optimal match is implemented by a min-cost
     * bipartite algorithm with a runtime cost of 
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
        float transX = fit.getParameters().getTranslationX();
        float transY = fit.getParameters().getTranslationY();

        return calculateMatchUsingOptimal(
            set1, set2, scale, rotation, transX, transY,
            tolerance, centroidX1, centroidY1);
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
     * @param scale
     * @param rotation
     * @param transX
     * @param tolerance
     * @param centroidX1
     * @param transY
     * @param centroidY1
     * @return a two dimensional array holding the matched indexes and
     * the distances between the model and the point for that pair.
     * each row holds float[]{idx1, idx2, diff}
     */
    public float[][] calculateMatchUsingOptimal(
        PairIntArray set1, PairIntArray set2, double scale, double rotation,
        float transX, float transY, double tolerance, 
        int centroidX1, int centroidY1) {
        
        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");

            log.severe("scale should not be smaller than 1");
        }

        // TODO: consider using class Transformer here?  consolidate code...
        
        double scaleTimesCosine = scale * Math.cos(rotation);
        double scaleTimesSine = scale * Math.sin(rotation);

        int nPoints1 = set1.getN();
        int nPoints2 = set2.getN();
        
        PairFloatArray transformed = new PairFloatArray();
        
        for (int i = 0; i < set1.getN(); i++) {

            int x = set1.getX(i);
            int y = set1.getY(i);

            double xmcx1 = x - centroidX1;
            double ymcy1 = y - centroidY1;

            double xr = centroidX1 * scale
                + ((xmcx1 * scaleTimesCosine) + (ymcy1 * scaleTimesSine));

            double yr = centroidY1 * scale
                + ((-xmcx1 * scaleTimesSine) + (ymcy1 * scaleTimesCosine));

            float xt = (float)(xr + transX);
            float yt = (float)(yr + transY);

            transformed.add(xt, yt);
        }
        
        float[][] matched = calculateMatchUsingOptimal(transformed, set2, 
            tolerance);
        
        return matched;
    }
    
    /**
     * given unordered unmatched points for a transformed set1 and 
     * the model set2 from image1 and image2
     * respectively, find the
     * optimal matching in set2 within tolerance and return the 
     * matched information as a two dimensional array of 
     * {index from set1, index from set2, diff of point in set2 from
     * model generation by point in set1}
     * 
     * @param transformed1
     * @param set2
     * @param tolerance
     * @return a two dimensional array holding the matched indexes and
     * the distances between the model and the point for that pair.
     * each row holds float[]{idx1, idx2, diff}
     */
    public float[][] calculateMatchUsingOptimal(
        PairFloatArray transformed1, PairIntArray set2, double tolerance) {
        
        int nPoints1 = transformed1.getN();
        int nPoints2 = set2.getN();
        
        float[][] diffsAsCost = new float[nPoints1][nPoints2];

        // the algorithm modifies diffsAsCost, so make a copy
        float[][] diffsAsCostCopy = new float[nPoints1][nPoints2];
        
        for (int i = 0; i < transformed1.getN(); i++) {

            diffsAsCost[i] = new float[nPoints2];
            diffsAsCostCopy[i] = new float[nPoints2];
            
            float x = transformed1.getX(i);
            float y = transformed1.getY(i);

            for (int j = 0; j < set2.getN(); j++) {

                int x2 = set2.getX(j);
                int y2 = set2.getY(j);

                double dist = Math.sqrt(Math.pow(x - x2, 2) 
                    + Math.pow(y - y2, 2));

                if (dist > tolerance) {
                    diffsAsCost[i][j] = Float.MAX_VALUE;
                    diffsAsCostCopy[i][j] = Float.MAX_VALUE;
                } else {
                    diffsAsCost[i][j] = (float)dist;
                    diffsAsCostCopy[i][j] = (float)dist;
                }
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
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }
                        
            if (diffsAsCost[idx1][idx2] < tolerance) {
                count++;
            }
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
            
            if (diff < tolerance) {
                output[count] = new float[3];
                output[count][0] = idx1;
                output[count][1] = idx2;
                output[count][2] = (float)diff;            
                count++;
            }
        }

        return output;
    }
    
    /**
     * given unordered unmatched points for a transformed set1 and 
     * the model set2 from image1 and image2
     * respectively, find the
     * greedy matching in set2 within tolerance and return the 
     * matched information as a two dimensional array of 
     * {index from set1, index from set2, diff of point in set2 from
     * model generation by point in set1}
     * 
     * @param transformed1
     * @param set2
     * @param tolerance
     * @return a two dimensional array holding the matched indexes and
     * the distances between the model and the point for that pair.
     * each row holds float[]{idx1, idx2, diff}
     */
    public float[][] calculateMatchUsingGreedy(
        PairFloatArray transformed1, PairIntArray set2, double tolerance) {
        
        Set<Integer> chosen = new HashSet<Integer>();        
        List<Integer> indexes1 = new ArrayList<Integer>();
        List<Integer> indexes2 = new ArrayList<Integer>();
        List<Double> diffs = new ArrayList<Double>();
        
        for (int i = 0; i < transformed1.getN(); i++) {

            float x = transformed1.getX(i);
            float y = transformed1.getY(i);

            double bestDiff = Double.MAX_VALUE;
            int best2Idx = -1;
            
            for (int j = 0; j < set2.getN(); j++) {

                if (chosen.contains(Integer.valueOf(j))) {
                    continue;
                }
                
                int x2 = set2.getX(j);
                int y2 = set2.getY(j);

                double dist = Math.sqrt(Math.pow(x - x2, 2) 
                    + Math.pow(y - y2, 2));

                if (dist < bestDiff) {
                    bestDiff = dist;
                    best2Idx = j;
                }
            }
            
            if ((bestDiff < Double.MAX_VALUE) && (bestDiff < tolerance)) {
                chosen.add(Integer.valueOf(best2Idx));
                indexes1.add(Integer.valueOf(i));
                indexes2.add(Integer.valueOf(best2Idx));
                diffs.add(Double.valueOf(bestDiff));                
            }
        }
        
        float[][] output = new float[indexes1.size()][];

        for (int i = 0; i < indexes1.size(); i++) {
            output[i] = new float[3];
            output[i][0] = indexes1.get(i).intValue();
            output[i][1] = indexes2.get(i).intValue();
            output[i][2] = diffs.get(i).intValue();
        }

        return output;
    }
    
    /**
     * given unordered unmatched points set1 and set2 from image1 and image2
     * respectively, apply the transformation to set 1 and then find the
     * nearest matching in set2 within tolerance and return the 
     * matched information as a two dimensional array of 
     * {index from set1, index from set2, diff of point in set2 from
     * model generation by point in set1} using an algorithm to greedily match
     * the first best match which then becomes unavailable for a better
     * later match.
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
    public float[][] calculateMatchUsingGreedy(
        PairIntArray set1, PairIntArray set2, TransformationPointFit fit,
        double tolerance, int centroidX1, int centroidY1) {
        
        double scale = fit.getParameters().getScale();
        double rotation = fit.getParameters().getRotationInRadians();
        float transX = fit.getParameters().getTranslationX();
        float transY = fit.getParameters().getTranslationY();

        return calculateMatchUsingGreedy(
            set1, set2, scale, rotation, transX, transY,
            tolerance, centroidX1, centroidY1);
    }
    
    public float[][] calculateMatchUsingGreedy(
        PairIntArray set1, PairIntArray set2, double scale, double rotation,
        float transX, float transY, double tolerance, 
        int centroidX1, int centroidY1) {
        
        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");

            log.severe("scale should not be smaller than 1");
        }

        double scaleTimesCosine = scale * Math.cos(rotation);
        double scaleTimesSine = scale * Math.sin(rotation);

        int nPoints1 = set1.getN();
        int nPoints2 = set2.getN();
        if ((nPoints1 == 0) || (nPoints2 == 0)) {
            return new float[0][];
        }
        
        Set<Integer> chosen = new HashSet<Integer>();        
        List<Integer> indexes1 = new ArrayList<Integer>();
        List<Integer> indexes2 = new ArrayList<Integer>();
        List<Double> diffs = new ArrayList<Double>();
            
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

            double bestDiff = Double.MAX_VALUE;
            int best2Idx = -1;
            
            for (int j = 0; j < set2.getN(); j++) {

                if (chosen.contains(Integer.valueOf(j))) {
                    continue;
                }
                
                int x2 = set2.getX(j);
                int y2 = set2.getY(j);

                double dist = Math.sqrt(Math.pow(xt - x2, 2) 
                    + Math.pow(yt - y2, 2));

                if (dist < bestDiff) {
                    bestDiff = dist;
                    best2Idx = j;
                }
            }
            
            if ((bestDiff < Double.MAX_VALUE) && (bestDiff < tolerance)) {
                chosen.add(Integer.valueOf(best2Idx));
                indexes1.add(Integer.valueOf(i));
                indexes2.add(Integer.valueOf(best2Idx));
                diffs.add(Double.valueOf(bestDiff));                
            }            
        }

        float[][] output = new float[indexes1.size()][];

        for (int i = 0; i < indexes1.size(); i++) {
            output[i] = new float[3];
            output[i][0] = indexes1.get(i).intValue();
            output[i][1] = indexes2.get(i).intValue();
            output[i][2] = diffs.get(i).intValue();
        }

        return output;
    }
    
    /**
     * given the indexes and residuals from optimal matching, populate
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
     * given the indexes and residuals from optimal matching, populate
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
    public void matchPoints(PairFloatArray set1, PairIntArray set2,
        double tolerance, float[][] matchedIndexesAndDiffs,
        PairFloatArray outputMatched1, PairIntArray outputMatched2) {
                
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
            
            float transformedX = scaledRotatedX[i] + transX;
            float transformedY = scaledRotatedY[i] + transY;
            
            double dx = transformedX - set2.getX(i);
            double dy = transformedY - set2.getY(i);
            
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
            
            float transformedX = unmatched1Transformed.getX(i);
            float transformedY = unmatched1Transformed.getY(i);
            
            diffsAsCost[i] = new float[n2];
            diffsAsCostCopy[i] = new float[n2];
            
            for (int j = 0; j < n2; j++) {
                
                float dx = transformedX - unmatched2.getX(j);
                float dy = transformedY - unmatched2.getY(j);
                
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
            
            float transformedX = scaledRotatedX[i] + transX;
            float transformedY = scaledRotatedY[i] + transY;
            
            diffsAsCost[i] = new float[set2.getN()];
            diffsAsCostCopy[i] = new float[set2.getN()];
            
            for (int j = 0; j < set2.getN(); j++) {
                
                float dx = transformedX - set2.getX(j);
                float dy = transformedY - set2.getY(j);
                
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
            
            float transformedX = unmatched1Transformed.getX(i);
            float transformedY = unmatched1Transformed.getY(i);
            
            double minDiff = Double.MAX_VALUE;
            int min2Idx = -1;

            for (int j = 0; j < unmatched2.getN(); j++) {

                if (chosen.contains(Integer.valueOf(j))) {
                    continue;
                }
                
                float dx = transformedX - unmatched2.getX(j);
                float dy = transformedY - unmatched2.getY(j);
                
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
     * given set1 which has already been transformed, calculate the average residual
     * between that and set2 and the standard deviation from the average.
     */
    public ProjectiveFit evaluateFitForUnMatchedOptimal(double[][] set1,
        double[][] set2,
        double tolTransX, double tolTransY) {
        
        if (set2 == null) {
            throw new IllegalArgumentException(
            "set2 cannot be null");
        }
        if (set1 == null) {
            throw new IllegalArgumentException(
            "set1 cannot be null");
        }
        
        int nPoints1 = set1.length;
        int nPoints2 = set2.length;
        
        float[][] diffsAsCost = new float[nPoints1][];
        float[][] diffsAsCostCopy = new float[nPoints1][];
        
        for (int i = 0; i < nPoints1; i++) {
            
            double transformedX = set1[i][0];
            double transformedY = set1[i][1];
            
            diffsAsCost[i] = new float[nPoints2];
            diffsAsCostCopy[i] = new float[nPoints2];
            
            for (int j = 0; j < nPoints2; j++) {
                
                double dx = set2[j][0] - transformedX;
                double dy = set2[j][1] - transformedY;
                
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
        
        ProjectiveFit fit = new ProjectiveFit();
        fit.setMeanDistFromModel(avgDiff);
        fit.setStdDevOfMean(stDev);
        fit.setTolerance(tolerance);
        fit.setNumberOfPoints(nMatched);
        
        return fit;
    }
    
    public ProjectiveFit evaluateFitForUnMatchedTransformedGreedy(
        double[][] set1, double[][] set2, float tolTransX, float tolTransY) {
        
        if (set1 == null) {
            throw new IllegalArgumentException("set1 cannot be null");
        }
        if (set2 == null) {
            throw new IllegalArgumentException("set2 cannot be null");
        }
        
        int nPoints1 = set1.length;
        int nPoints2 = set2.length;
        
        Set<Integer> chosen = new HashSet<Integer>();
        
        double[] diffs = new double[nPoints1];
        int nMatched = 0;
        double avg = 0;
        
        for (int i = 0; i < nPoints1; i++) {
            
            double transformedX = set1[i][0];
            double transformedY = set1[i][1];
            
            double minDiff = Double.MAX_VALUE;
            int min2Idx = -1;

            for (int j = 0; j < nPoints2; j++) {

                if (chosen.contains(Integer.valueOf(j))) {
                    continue;
                }
                
                double dx = set2[j][0] - transformedX;
                double dy = set2[j][1] - transformedY;
                
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
        
        double tolerance = Math.sqrt(tolTransX*tolTransX + tolTransY*tolTransY);
        
        ProjectiveFit fit = new ProjectiveFit();
        fit.setMeanDistFromModel(avg);
        fit.setStdDevOfMean(stDev);
        fit.setTolerance(tolerance);
        fit.setNumberOfPoints(nMatched);
        
        return fit;
    }
    
    /**
     * given the transformed x y that have already been scaled and rotated, add the
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
            
            float transformedX = scaledRotatedX[i] + transX;
            float transformedY = scaledRotatedY[i] + transY;
            
            double minDiff = Double.MAX_VALUE;
            int min2Idx = -1;

            for (int j = 0; j < set2.getN(); j++) {

                if (chosen.contains(Integer.valueOf(j))) {
                    continue;
                }
                
                float dx = set2.getX(j) - transformedX;
                float dy = set2.getY(j) - transformedY;
                
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
            
            float transformedX = matched1Transformed.getX(i);
            float transformedY = matched1Transformed.getY(i);
            
            int x2 = matched2.getX(i);
            int y2 = matched2.getY(i);
            
            double dx = x2 - transformedX;
            double dy = y2 - transformedY;
                            
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
            
            float transformedX = (centroidX1*scale + ( 
                ((x1 - centroidX1) * scaleTimesCosine) +
                ((y1 - centroidY1) * scaleTimesSine))) + transX;
            
            float transformedY = (centroidY1*scale + ( 
                (-(x1 - centroidX1) * scaleTimesSine) +
                ((y1 - centroidY1) * scaleTimesCosine))) + transY;
            
            double minDiff = Double.MAX_VALUE;

            for (int j = 0; j < edge2.getN(); j++) {

                int x2 = edge2.getX(j);
                int y2 = edge2.getY(j);
                
                float dx = x2 - transformedX;
                float dy = y2 - transformedY;
                
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
        double convergence = scaledRotatedX.length;

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

    private void performGreedyMatch(PairFloatArray set1, 
        PairIntArray set2, double tolTransX, double tolTransY, 
        PairFloatArray output1, PairFloatArray output2) {
        
        Set<Integer> chosen = new HashSet<Integer>();
        
        for (int i = 0; i < set1.getN(); i++) {

            float x1 = set1.getX(i);
            float y1 = set1.getY(i);

            double minDiffSq = Double.MAX_VALUE;
            int bestIdx2 = -1;
            
            for (int j = 0; j < set2.getN(); j++) {

                if (chosen.contains(Integer.valueOf(j))) {
                    continue;
                }
                
                int x2 = set2.getX(j);
                int y2 = set2.getY(j);
                
                float xDiff = Math.abs(x2 - x1);
                float yDiff = Math.abs(y2 - y1);

                if ((xDiff > tolTransX) || (yDiff > tolTransY)) {
                    continue;
                }
                
                double diffSq = (xDiff*xDiff) + (yDiff*yDiff);

                if (diffSq < minDiffSq) {
                    minDiffSq = diffSq;
                    bestIdx2 = j;
                }
            }
            
            if (bestIdx2 > -1) {
                output1.add(x1, y1);
                output2.add(set2.getX(bestIdx2), set2.getY(bestIdx2));                
                chosen.add(Integer.valueOf(bestIdx2));
            }            
        }
    }
    
    /**
     * calculate for unmatched points and if best match is not good,
     * reverse the order of sets and try again in order to solve for
     * possible scale transformation smaller than 1.
     * 
     * @param scene
     * @param model
     * @param sceneImageCentroidX
     * @param sceneImageCentroidY
     * @param modelImageCentroidX
     * @param modelImageCentroidY
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * 
     * @return 
     */
    public TransformationPointFit calculateProjectiveTransformationWrapper(
        PairIntArray scene, PairIntArray model, 
        int sceneImageCentroidX, int sceneImageCentroidY, 
        int modelImageCentroidX, int modelImageCentroidY,
        float setsFractionOfImage) {
        
        TransformationPointFit fit = calculateProjectiveTransformation(
            scene, model, sceneImageCentroidX, sceneImageCentroidY,
            setsFractionOfImage);
        
        if (fit == null) {
            return null;
        }
        
        int nMaxMatchable = (scene.getN() < model.getN()) ? scene.getN() 
            : model.getN();
        
        float fracMatched = (fit == null) ? 0 :
            (float)fit.getNumberOfMatchedPoints()/(float)nMaxMatchable;
        
        if (fracMatched < 0.3) {
            
            // reverse the order to solve for possible scale < 1.
            
            TransformationPointFit revFit = calculateProjectiveTransformation(
                model, scene, modelImageCentroidX, modelImageCentroidY,
                setsFractionOfImage);
            
            if (fitIsBetter(fit, revFit)) {
                
                TransformationParameters params = revFit.getParameters();
                
                // reverse the parameters.
                // needs a reference point in both datasets for direct calculation
                
                MatchedPointsTransformationCalculator tc = new 
                    MatchedPointsTransformationCalculator();
                
                double[] x1y1 = tc.applyTransformation(params, 
                    modelImageCentroidX, modelImageCentroidY,
                    modelImageCentroidX, modelImageCentroidY);
                
                TransformationParameters revParams = tc.swapReferenceFrames(
                    params, sceneImageCentroidX, sceneImageCentroidY, 
                    modelImageCentroidX, modelImageCentroidY, 
                    x1y1[0], x1y1[1]);       
                
                fit = new TransformationPointFit(revParams, 
                    revFit.getNumberOfMatchedPoints(),
                    revFit.getMeanDistFromModel(),
                    revFit.getStDevFromMean(),
                    revFit.getTolerance());
            }
        }
        
        return fit;
    }

    /**
     * calculate for unmatched points
     * @param scene
     * @param model
     * @param sceneImageCentroidX
     * @param sceneImageCentroidY
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * 
     * @return 
     */
    public TransformationPointFit calculateProjectiveTransformation(
        PairIntArray scene, PairIntArray model, int sceneImageCentroidX, 
        int sceneImageCentroidY, float setsFractionOfImage) {
        
        if ((scene == null) || (model == null)) {
            return null;
        }
        if ((scene.getN() < 3) || (model.getN() < 3)) {
            return null;
        }
        
        int rotStart = 0; 
        int rotStop = 360; 
        int rotDelta = 10;
        int scaleStart = 1;
        int scaleStop = 5;
        int scaleDelta = 1;
        boolean setsAreMatched = false;
                    
        TransformationPointFit fit = calculateTransformationWithGridSearch(
            scene, model, sceneImageCentroidX, sceneImageCentroidY,
            rotStart, rotStop, rotDelta, scaleStart, scaleStop, scaleDelta,
            setsAreMatched, setsFractionOfImage);
        
        if (fit != null) {
            log.info("fit: " + fit.toString());
        }
        int nMaxMatchable = (scene.getN() < model.getN()) ? scene.getN() 
            : model.getN();
        
        // ==== TODO: improve decision for continue searching =====
  
        boolean refineEuclideanParams = 
            !((nMaxMatchable == fit.getNumberOfMatchedPoints()) && 
            (fit.getMeanDistFromModel() < 5));
                
        // TODO: follow the results here for the stereo projective set
        // Brown & Lowe 2003.  small projective coeffs should be needed.
        
        if (refineEuclideanParams) {
            
            //TODO: improve this:  
            // frac=0.7 can be safely used for matches and avoid refinement here?
            //
            
            float frac = 0.95f;
            
            float fracMatched = (float)fit.getNumberOfMatchedPoints()/
                (float)nMaxMatchable;
            
            if (true /*fracMatched <= frac*/) {

                // refine the euclidean solution, centered on current values

                TransformationPointFit fit2 = 
                    refineTransformationWithDownhillSimplex(fit.getParameters(),
                    scene, model, sceneImageCentroidX, sceneImageCentroidY, 
                    setsAreMatched, setsFractionOfImage);

                boolean refineTranslation = false;
                
                if (fitIsBetter(fit, fit2)) {
                    
                    fit = fit2;
                    
                    refineTranslation = true;
                    
                } else if ((fit2 != null) && (fit != null)) {
                    
                    // consider the case when avg diff is better and stdev too
                    double diffMean = 
                        fit.getMeanDistFromModel() - fit2.getMeanDistFromModel();
                    
                    //TODO: improve this
                    if ((diffMean > 5) && ((fit.getStDevFromMean() 
                        - fit2.getStDevFromMean()) > 0.)) {
                        
                        float fracMatched2 = (float)fit2.getNumberOfMatchedPoints()/
                            (float)nMaxMatchable; 
                        
                        if (fracMatched2 > 0.25) {
                            
                            log.info("choosing fit2 over fit1 even though" +
                                " nMatched is smaller");
                            log.info("fit=" + fit.toString());
                            log.info("fit2=" + fit2.toString());
                            
                            fit = fit2;
                            
                            refineTranslation = true;
                        }    
                    }
                }
                
                //TODO:  revisit this.  should be able to remove code after tests
                if (false) {
                //if (refineTranslation) {
                    
                    refineTranslation(fit.getParameters(), scene, model, 
                        sceneImageCentroidX, sceneImageCentroidY);
                    
                    Transformer transformer = new Transformer();
                    
                    PairFloatArray transformed = transformer.applyTransformation(
                        fit.getParameters(), 
                        sceneImageCentroidX, sceneImageCentroidY, scene);

                    double transXTol = 2. * sceneImageCentroidX * 0.02;
                    double transYTol = 2 * sceneImageCentroidY * 0.02;
        
                    fit = evaluateFitForUnmatchedTransformed(
                        fit.getParameters(),
                        transformed, model, transXTol, transYTol);
                }
            }
        }
        
        if (refineEuclideanParams) {
            // === need criteria for which to try for a projective solution.
            //     -- a plot of x vs diffX and y vs diffY ?
            boolean doProjectiveSearch = false;

            if (doProjectiveSearch) {
                throw new UnsupportedOperationException("not yet implemented");
            }
        }
        
        return fit;
    }
    
    private void refineTranslation(TransformationParameters params, 
        PairIntArray scene, PairIntArray model, 
        int sceneImageCentroidX, int sceneImageCentroidY) {
                    
        double translationX = params.getTranslationX();
        double translationY = params.getTranslationY();
        double scale = params.getScale();
        double rotation = params.getRotationInRadians();
        double scaleTimesCosine = scale * Math.cos(rotation);
        double scaleTimesSine = scale * Math.sin(rotation);
        double[][] transformed = new double[scene.getN()][];
        for (int i = 0; i < scene.getN(); i++) {
            transformed[i] = new double[2];
            int x = scene.getX(i);
            int y = scene.getY(i);
            transformed[i][0] = sceneImageCentroidX*scale + (
                ((x - sceneImageCentroidX) * scaleTimesCosine) +
                ((y - sceneImageCentroidY) * scaleTimesSine));
            transformed[i][1] = sceneImageCentroidY*scale + (
                (-(x - sceneImageCentroidX) * scaleTimesSine) +
                ((y - sceneImageCentroidY) * scaleTimesCosine));
        }

        double[][] compare = new double[model.getN()][];
        for (int i = 0; i < model.getN(); i++) {
            compare[i] = new double[2];
            compare[i][0] = model.getX(i);
            compare[i][1] = model.getY(i);
        }

        float[] tolerances = new float[]{30, 5};
        for (float tolerance : tolerances) {

            float[] peakXY = refineCalculateTranslation(
                transformed, compare,
                (float)translationX, (float)translationY, tolerance);

            if (peakXY != null) {

                translationX = (int)peakXY[0];
                translationY = (int)peakXY[1];

                log.info("refined: " + translationX + ", " + translationY);
            }
        }

        params.setTranslationX((float)translationX);
        params.setTranslationY((float)translationY);
 
    }
   
    void sortByDescendingMatches(ProjectiveFit[] fits, int idxLo, 
        int idxHi) {
        
        if (idxLo < idxHi) {

            int idxMid = partition(fits, idxLo, idxHi);

            sortByDescendingMatches(fits, idxLo, idxMid - 1);

            sortByDescendingMatches(fits, idxMid + 1, idxHi);
        }
    }
    
    public boolean fitIsBetter(ProjectiveFit bestFit, ProjectiveFit compareFit) {
        
        if (compareFit == null) {
            return false;
        }
        if (bestFit == null) {
            return true;
        }
               
        int nMatches = compareFit.getNumberOfPoints();
        
        if (nMatches > bestFit.getNumberOfPoints()) {
            return true;
        } else if (nMatches == bestFit.getNumberOfPoints()) {
            if (!Double.isNaN(compareFit.getMeanDistFromModel()) && (
                compareFit.getMeanDistFromModel()
                < bestFit.getMeanDistFromModel())) {
                return true;
            } else if (compareFit.getMeanDistFromModel()
                == bestFit.getMeanDistFromModel()) {
                if (compareFit.getStdDevOfMean() < bestFit.getStdDevOfMean()) {
                    return true;
                }
            }
        }
        return false;
    }

    private int partition(ProjectiveFit[] fits, int idxLo, int idxHi) {
        
        ProjectiveFit x = fits[idxHi];
        
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            if (fitIsBetter(x, fits[i])) {
                store++;
                ProjectiveFit swap = fits[store];
                fits[store] = fits[i];
                fits[i] = swap;
            }
        }
        
        store++;
        ProjectiveFit swap = fits[store];
        fits[store] = fits[idxHi];
        fits[idxHi] = swap;
        
        return store;
    }
     
    private ProjectiveFit evalFit(double[][] scene, double[][] model,
        double[][] projTransformParams, 
        double[][] allScene, double[][] allModel,
        int centroidX1, int centroidY1, float setsFractionOfImage) {
        
        double[][] transformed = transformUsingProjection(
            projTransformParams, scene);
        
        /*
        solve for transX and transY and add those to transformed
        */
        float[] transXY = calculateTranslation(transformed, model,
            centroidX1, centroidY1, setsFractionOfImage);
        projTransformParams[0][2] = transXY[0];
        projTransformParams[1][2] = transXY[1];
        
        transformed = transformUsingProjection(projTransformParams, allScene);
        
        //TODO: calculate these:
        double tolX = 10;
        double tolY = 10;
        
        ProjectiveFit fit = evaluateFitForUnMatchedOptimal(
            transformed, allModel, tolX, tolY);
        fit.setProjection(projTransformParams);
        
        return fit;
    }
    
    private ProjectiveFit evalFit(double translationX, double translationY,
        double[][] allScene, double[][] allModel) {
        
        double[][] transformed = transformUsingProjection(
            translationX, translationY, allScene);
        
        //TODO: calculate these:
        double tolX = 10;
        double tolY = 10;
        
        ProjectiveFit fit = evaluateFitForUnMatchedOptimal(
            transformed, allModel, tolX, tolY);
        fit.setData(new Double[]{Double.valueOf(translationX),
            Double.valueOf(translationY)});
        
        return fit;
    }
    
    public double[][] transformUsingProjection(double[][] projTransformParams,
        double[][] matrix) {
        /*      
        //=========
        int nPoints1 = matrix.length;
        double[][] output0 = new double[nPoints1][3];
        for (int i = 0; i < nPoints1; i++) {
            output0[i] = new double[3];
            output0[i][0] = matrix[i][0];
            output0[i][1] = matrix[i][1];
            output0[i][2] = 1;
        }
        SimpleMatrix m = new SimpleMatrix(output0);
        double[][] result = transformUsingProjection(projTransformParams, m);
        double[][] compare = MatrixUtil.multiplyByTranspose(projTransformParams, 
            output0);
        //=========
        */
        
        // the multiplication altered to remove last column of '1's from matrix
        int prows = projTransformParams.length;
        int pcols = projTransformParams[0].length;
        int mrows = matrix.length;
        int mcols = matrix[0].length;
        
        double[][] output = new double[mrows][];
        for (int nrow = 0; nrow < mrows; nrow++) {
            output[nrow] = new double[mcols];
            for (int row = 0; row < prows; row++) {
                double sum = 0;         
                for (int mcol = 0; mcol < mcols; mcol++) {
                    sum += (projTransformParams[row][mcol] * matrix[nrow][mcol]);                    
                }
                sum += projTransformParams[row][2];
                if (row < mcols) {
                    output[nrow][row] = sum;
                }
            }       
        }
        /*                        m[0][0],m[0][1]
                   1  0  -19      274.99
                                  244.888
                                   1
        | x2 |   | 1  0  t_x |   | x |   | x + t_x |
        | y2 | = | 0  1  t_y | * | y | = | x + t_x |
        | 1  |   | 0  0  1   |   | 1 |   | 1       |
        */
        
        return output;
    }
    
    public double[][] transformUsingProjection(
        double translationX, double translationY, double[][] matrix) {
        
        int nPoints1 = matrix.length;
        double[][] output = new double[nPoints1][matrix[0].length];
        for (int i = 0; i < nPoints1; i++) {
            output[i] = new double[matrix[0].length];
            output[i][0] = matrix[i][0] + translationX;
            output[i][1] = matrix[i][1] + translationY;
        }
       
        return output;
    }
    
    private boolean isASinglePeak(HistogramHolder h) {
        
        int idx =  MiscMath.findYMaxIndex(h.getYHist());
        
        int lastY = Integer.MIN_VALUE;
        
        for (int i = 0; i <= idx; i++) {
            int y = h.getYHist()[i];
            if (y <= lastY) {
                return false;
            }
            lastY = y;
        }
        //lastY = h.getYHist()[idx];
        for (int i = (idx + 1); i < h.getYHist().length; i++) {
            int y = h.getYHist()[i];
            if (y >= lastY) {
                return false;
            }
            lastY = y;
        }
        
        return true;
    }
   
    private float[] refineCalculateTranslation(double[][] set1, double[][] set2, 
        float peakTransX, float peakTransY, float tolerance) {
        
        return refineCalculateTranslationWithDownhillSimplex(set1, set2,
            peakTransX, peakTransY, tolerance);
    }

    private float[] refineCalculateTranslationWithDownhillSimplex(
        double[][] set1, double[][] set2, float peakTransX, float peakTransY, 
        float range) {
        
        float reflectionCoeff = 1;   // > 0
        float expansionCoeff = 2;   // > 1
        float contractionCoeff = -0.5f;
        float reductionCoeff = 0.5f;
        
        double convergence = 0;
        
        float transX = peakTransX;
        float transY = peakTransY;
        float txMin = peakTransX - range;
        float txMax = peakTransX + range;
        float tyMin = peakTransY - range;
        float tyMax = peakTransY + range;
        
        float[] txs = new float[]{
            (peakTransX - range),
            (peakTransX - 0.75f*range),
            (peakTransX - 0.5f*range),
            (peakTransX - 0.25f*range),
            peakTransX,
            (peakTransX + 0.25f*range),
            (peakTransX + 0.5f*range),
            (peakTransX + 0.75f*range)
        };
        float[] tys = new float[]{
            (peakTransY - range),
            (peakTransY - 0.75f*range),
            (peakTransY - 0.5f*range),
            (peakTransY - 0.25f*range),
            peakTransY,
            (peakTransY + 0.25f*range),
            (peakTransY + 0.5f*range),
            (peakTransY + 0.75f*range)
        };
        // fitting for 2 parameters so need at least 3 starter points.
        int n = txs.length * tys.length;
        ProjectiveFit[] fits = new ProjectiveFit[n];
        int count = 0;
        for (float tx : txs) {
            for (float ty : tys) {
                fits[count] = evalFit(tx, ty, set1, set2);                
                count++;
            }
        }
        
        boolean go = true;

        int nMaxIter = 100;
        int nIter = 0;
        
        int bestFitIdx = 0;
        int secondWorstFitIdx = fits.length - 2;
        int worstFitIdx = fits.length - 1;

        double lastAvgDistModel = Double.MAX_VALUE;
        double lastStDev = Double.MAX_VALUE;
        int lastNMatches = 0;
        
        int nIterSameMin = 0;
       
        while (go && (nIter < nMaxIter)) {

            sortByDescendingMatches(fits, 0, worstFitIdx);
            
            if ((lastNMatches == fits[bestFitIdx].getNumberOfPoints())
                && (Math.abs(lastAvgDistModel
                - fits[bestFitIdx].getMeanDistFromModel()) < 0.01)) {
                
                nIterSameMin++;
                if (nIterSameMin >= 4) {
                    break;
                }
            } else {
                nIterSameMin = 0;
            }
            lastNMatches = fits[0].getNumberOfPoints();
            lastAvgDistModel = fits[0].getMeanDistFromModel();
            lastStDev = fits[0].getStdDevOfMean();
            
            // determine center for all points excepting the worse fit
            float txSum = 0.0f;
            float tySum = 0.0f;
            int nt = 0;
            for (int i = 0; i <= worstFitIdx; i++) {
                if (fits[i].getData() == null) {
                    continue;
                }
                Double[] txy = (Double[])fits[i].getData();
                txSum += txy[0].doubleValue();
                tySum += txy[1].doubleValue();
                nt++;
            }
            transX = txSum/(float)(nt - 1);
            transY = tySum/(float)(nt - 1);

            // "Reflection"
            Double[] txy = (Double[])fits[worstFitIdx].getData();
            float txReflect = transX + (reflectionCoeff
                * (transX - txy[0].floatValue()));
            float tyReflect = transY + (reflectionCoeff
                * (transY - txy[1].floatValue()));

            ProjectiveFit fitReflected
                = evalFit(txReflect, tyReflect, set1, set2);
            
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
                    float txExpansion = transX + (expansionCoeff
                        * (transX - txy[0].floatValue()));
                    float tyExpansion = transY + (expansionCoeff
                        * (transY - txy[1].floatValue()));
                    
                    ProjectiveFit fitExpansion
                         = evalFit(txExpansion, tyExpansion, set1, set2);

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
                    float txContraction = transX + (contractionCoeff
                        * (transX - txy[0].floatValue()));
                    float tyContraction = transY + (contractionCoeff
                        * (transY - txy[1].floatValue()));

                    ProjectiveFit fitContraction
                        = evalFit(txContraction, tyContraction, set1, set2);

                    if (fitIsBetter(fits[worstFitIdx], fitContraction)
                        && (txContraction >= txMin) && (txContraction <= txMax)
                        && (tyContraction >= tyMin) && (tyContraction <= tyMax)) {

                        fits[worstFitIdx] = fitContraction;

                    } else {

                        if (fits[bestFitIdx].getData() == null) {
                            break;
                        }
                        
                        Double[] txybf = (Double[])fits[bestFitIdx].getData();
                        
                        // "Reduction"
                        for (int i = 1; i <= worstFitIdx; i++) {

                            Object txyiObj = fits[i].getData();
                            
                            if (txyiObj == null) {
                                fits[i] = new ProjectiveFit();
                                fits[i].setMeanDistFromModel(Double.MAX_VALUE);
                                fits[i].setNumberOfPoints(0);
                                fits[i].setStdDevOfMean(Double.MAX_VALUE);
                                fits[i].setTolerance(Double.MAX_VALUE);
                                continue;
                            }
                            
                            Double[] txyi = (Double[])txyiObj;
                            
                            float txReduction
                                = (txybf[0].floatValue()
                                + (reductionCoeff
                                * (txyi[0].floatValue()
                                - txybf[0].floatValue())));

                            float tyReduction
                                = (txybf[1].floatValue()
                                + (reductionCoeff
                                * (txyi[1].floatValue()
                                - txybf[1].floatValue())));

                            //NOTE: there's a possibility of a null fit.
                            //  instead of re-writing the fits array, will 
                            //  assign a fake infinitely bad fit which will 
                            //  fall to the bottom of the list after the next 
                            //  sort.
                            ProjectiveFit fit
                                = evalFit(txReduction, tyReduction, set1, set2);
                                    
                            if (fit != null) {
                                fits[i] = fit;
                            } else {
                                fits[i] = new ProjectiveFit();
                                fits[i].setMeanDistFromModel(Double.MAX_VALUE);
                                fits[i].setNumberOfPoints(0);
                                fits[i].setStdDevOfMean(Double.MAX_VALUE);
                                fits[i].setTolerance(Double.MAX_VALUE);
                            }
                        }
                    }
                }
            }

            log.finest("best fit so far: nMatches="
                + fits[bestFitIdx].getNumberOfPoints()
                + " diff from model=" + fits[bestFitIdx].getMeanDistFromModel()
            );

            nIter++;

            if ((fits[bestFitIdx].getNumberOfPoints() == convergence)
                && (fits[bestFitIdx].getMeanDistFromModel() == 0)) {
                go = false;
            } else if ((transX > txMax) || (transX < txMin)) {
                go = false;
            } else if ((transY > tyMax) || (transY < tyMin)) {
                go = false;
            }
        }

        if ((fits[bestFitIdx] == null) || (fits[bestFitIdx].getData() == null)) {
            return null;
        }
        
        Double[] txy = (Double[])fits[bestFitIdx].getData();
        
        float tx = txy[0].floatValue();
        float ty = txy[1].floatValue();
        
        ProjectiveFit bestFit = fits[bestFitIdx];
        // one last set of +- 1 tries around solution
        txs = new float[]{tx - 1, tx, tx + 1};
        tys = new float[]{ty - 1, ty, ty + 1};
        for (int i = 0; i < txs.length; i++) {
            for (int j = 0; j < tys.length; j++) {
                if (i == j) { 
                    continue;
                }
                ProjectiveFit cFit = evalFit(txs[i], tys[j], set1, set2);
                if (fitIsBetter(bestFit, cFit)) {
                    bestFit = cFit;
                }
            }
        }
        
        if (!bestFit.equals(fits[bestFitIdx])) {
            txy = (Double[])bestFit.getData();
            tx = txy[0].floatValue();
            ty = txy[1].floatValue();
        }
        
        return new float[]{Math.round(tx), Math.round(ty)};
    }

    private int moveUpForNulls(ProjectiveFit[] fits) {
    
        boolean hadANull = false;
        
        for (int i = (fits.length - 1); i > -1; i--) {
            ProjectiveFit fit = fits[i];
            if ((fit == null) || (fit.getProjection() == null)) {
                
                hadANull = true;
                
                fits[i] = null;
                
                int nullIdx = -1;
                
                // move everything with larger index up by one
                for (int j = (i + 1); j < fits.length; j++) {
                    ProjectiveFit move = fits[j];
                    if ((move == null) || (move.getProjection() == null)) {
                        nullIdx = j;
                        continue;
                    }
                    fits[j - 1] = move;
                }
                if ((nullIdx - 1) > -1) {
                    fits[nullIdx - 1] = null;
                }
            }
        }
        
        int lastNonNull = fits.length - 1;
        
        if (hadANull) {
            for (int i = 0; i < fits.length; i++) {
                ProjectiveFit fit = fits[i];
                if (fit == null) {
                    lastNonNull = (i - 1);
                    break;
                }
            }
        }
        
        return lastNonNull;
    }

    private TransformationPointFit evaluateForBestTranslation(float[] xTranslations, 
        float[] yTranslations, float tolTransX, float tolTransY,
        float[] scaledRotatedX, float[] scaledRotatedY, 
        PairIntArray set2) {
        
        TransformationPointFit bestFit = null;
        
        float scalePlaceHolder = 1;
        float rotationPlaceHolder = 0;
        
        int nMaxMatchable = (scaledRotatedX.length < set2.getN()) ? 
            scaledRotatedX.length : set2.getN();
        
        for (float transX : xTranslations) {
            for (float transY : yTranslations) {
                
                TransformationPointFit fit = evaluateFitForUnMatchedGreedy(
                    scaledRotatedX, scaledRotatedY, transX, transY,
                    tolTransX, tolTransY,
                    set2, scalePlaceHolder, rotationPlaceHolder);
                
                /*
                TransformationPointFit fit = evaluateFitForUnMatchedOptimal(
                    scaledRotatedX, scaledRotatedY, transX, transY,
                    tolTransX, tolTransY,
                    set2, scalePlaceHolder, rotationPlaceHolder);
                */
                
                if (fitIsBetter(bestFit, fit)) {
                    bestFit = fit;
                }
            }
        }
        // note, this is missing the correct rot and scale:
        return bestFit;
    }

    /**
     * from histograms of scaled and rotated set1 subtracted from each 
     * possible pairing with set2, determine the translations.
     * For the simplest of cases, the translation is just the peak of the
     * histograms, but for more complex, further calculations are
     * performed.
     * 
     * @param hX
     * @param hY
     * @param tolTransX
     * @param tolTrnsY
     * @return 
     */
    private float[] determineTranslationFromHistograms(HistogramHolder hX, 
        HistogramHolder hY, float[] xr, float[] yr, PairIntArray set2,
        float tolTransX, float tolTransY) {
        
        float peakTransX = Float.MAX_VALUE;
        float peakTransY = Float.MAX_VALUE;
         
        List<Integer> xPeakIndexes = MiscMath.findStrongPeakIndexes(hX, 0.09f);
        
        List<Integer> yPeakIndexes = MiscMath.findStrongPeakIndexes(hY, 0.09f);
        
        //if either indexes sizes != 1, try combinations of each peak above half max 
                
        if (false /*(xPeakIndexes.size() == 1) && (yPeakIndexes.size() == 1)*/) {
        
            int peakIdx = xPeakIndexes.get(0);
            
            int freqPeak = hX.getYHist()[peakIdx];
            float crit = 0.93f * freqPeak;
            
            if ((peakIdx > 0) && (peakIdx < (hX.getXHist().length - 1)) &&
                (hX.getYHist()[peakIdx - 1] <= crit) &&
                (hX.getYHist()[peakIdx + 1] <= crit)
                ) {
                
                // frequency suggests this is the real translation
                peakTransX = hX.getXHist()[peakIdx];
                
            } else {

                ArrayPair xy = LinesAndAngles.createPolygonOfTopFWFractionMax(
                    hX.getXHist(), hX.getYHistFloat(), null, null, 0.5f);

                //[1] is peakTransX
                float[] xAreaAndCentroid
                    = LinesAndAngles.calcAreaAndCentroidOfSimplePolygon(xy.getX(),
                        xy.getY());
        
                peakTransX = xAreaAndCentroid[1];
            }
        
            peakIdx = yPeakIndexes.get(0);
            
            freqPeak = hY.getYHist()[peakIdx];
            crit = 0.93f * freqPeak;
            
            if ((peakIdx > 0) && (peakIdx < (hY.getXHist().length - 1)) &&
                (hY.getYHist()[peakIdx - 1] <= crit) &&
                (hY.getYHist()[peakIdx + 1] <= crit)
                ) {
                
                // frequency suggests this is the real translation
                peakTransY = hY.getXHist()[peakIdx];
                
            } else {

                ArrayPair xy = LinesAndAngles.createPolygonOfTopFWFractionMax(
                    hY.getXHist(), hY.getYHistFloat(), null, null, 0.5f);

                //[1] is peakTransY
                float[] yAreaAndCentroid
                    = LinesAndAngles.calcAreaAndCentroidOfSimplePolygon(
                        xy.getX(), xy.getY());
                
                peakTransY = yAreaAndCentroid[1];
            }
            
        } else {
            
            // for each point above half max, try combination
            float[] txs = MiscMath.extractAllXForYAboveHalfMax(hX);
            Arrays.copyOf(txs, txs.length + 1);
            txs[txs.length - 1] = 0;
            
            float[] tys = MiscMath.extractAllXForYAboveHalfMax(hY);
            tys[tys.length - 1] = 0;
            
            TransformationPointFit fitForTranslation = 
                evaluateForBestTranslation(
                txs, tys, tolTransX, tolTransY, xr, yr, set2);

            if (fitForTranslation != null) {
                peakTransX = (float)fitForTranslation.getTranslationX();
                peakTransY = (float)fitForTranslation.getTranslationY();
            }
        }
        
        return new float[]{peakTransX, peakTransY};
    }

    private void populateIntersectionOfRegions(
        TransformationParameters params,
        PairIntArray set1, PairIntArray set2, int image1CentroidX, 
        int image1CentroidY, int image2CentroidX, int image2CentroidY, 
        PairIntArray filtered1, PairIntArray filtered2) {
        
        double tolerance = 20;
        
        // determine the bounds of filtered2.  any points in set2 that have
        // x < xy2LL[0] will not be matchable, etc.
        // TODO: correct this to use "point inside polygon" instead of rectangle
        
        Transformer transformer = new Transformer();

        double[] xy2LL = transformer.applyTransformation(params,
            image1CentroidX, image1CentroidY, 0, 0);

        double[] xy2UR = transformer.applyTransformation(params,
            image1CentroidX, image1CentroidY,
            2*image1CentroidX - 1, 2*image1CentroidY - 1);
        
        // to find lower-left and upper-left in image 1 of image
        // 2 boundaries requires the reverse parameters

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();
        double[] x1cy1c = tc.applyTransformation(params,
            image1CentroidX, image1CentroidY,
            image1CentroidX, image1CentroidY);

        TransformationParameters revParams = tc.swapReferenceFrames(
            params, image2CentroidX, image2CentroidY,
            image1CentroidX, image1CentroidY,
            x1cy1c[0], x1cy1c[1]);

        double[] xy1LL = transformer.applyTransformation(
            revParams, image2CentroidX, image2CentroidY,
            0, 0);
        xy1LL[0] -= tolerance;
        xy1LL[1] -= tolerance;

        double[] xy1UR = transformer.applyTransformation(
            revParams, image2CentroidX, image2CentroidY,
            2*image2CentroidX - 1, 2*image2CentroidY - 1);
        xy1UR[0] += tolerance;
        xy1UR[1] += tolerance;
        
        for (int i = 0; i < set1.getN(); i++) {
            
            int x = set1.getX(i);
            
            // TODO: replace with an "outside polygon" check
            if ((x < xy1LL[0]) || (x > xy1UR[0])) {
                continue;
            }
            
            int y = set1.getY(i);
            
            if ((y < xy1LL[1]) || (y > xy1UR[1])) {
                continue;
            }
            
            filtered1.add(x, y);
        }
        
        for (int i = 0; i < set2.getN(); i++) {
            
            int x = set2.getX(i);
            
            // TODO: replace with an "outside polygon" check
            if ((x < xy2LL[0]) || (x > xy2UR[0])) {
                continue;
            }
            
            int y = set2.getY(i);
            
            if ((y < xy2LL[1]) || (y > xy2UR[1])) {
                continue;
            }
            
            filtered2.add(x, y);
        }
    }

    private boolean fitIsSimilar(TransformationPointFit fit1, 
        TransformationPointFit fit2,
        float scaleTolerance, float rotInDegreesTolerance,
        float translationTolerance) {
        
        if (fit1 == null || fit2 == null) {
            return false;
        }
        
        TransformationParameters params1 = fit1.getParameters();
        
        TransformationParameters params2 = fit2.getParameters();
        
        float rotA = params1.getRotationInDegrees();
        float rotB = params2.getRotationInDegrees();
        if (rotA > rotB) {
            float swap = rotA;
            rotA = rotB;
            rotB = swap;
        }
        if (((rotB - rotA) > rotInDegreesTolerance) && 
            (((rotA + 360) - rotB) > rotInDegreesTolerance)) {
            return false;
        }
        
        if ((Math.abs(params1.getScale() - params2.getScale()) <= scaleTolerance)
            && (Math.abs(params1.getTranslationX() - params2.getTranslationX()) <= translationTolerance)
            && (Math.abs(params1.getTranslationY() - params2.getTranslationY()) <= translationTolerance)
            ) {

            return true;
        }
        
        return false;
    }

    /**
     * a fitness function that tries to allow a smaller number of points roughly
     * fit to be seen in contrast to a larger number of points that 
     * are not a better fit, but have better stats due to matching alot of
     * scattered points.
     * if numberOfMatched/maxMatchable is not infinite:
     * if the mean/10 of the comparison fit is better than the best mean/10,
     * a true is returned, else
     * compares the standard
     * deviation from the mean difference to the model fit and returns true
     * if compareFit has a smaller value.
     * @param bestFit
     * @param compareFit
     * @return 
     */
    protected boolean fitIsBetter2(TransformationPointFit bestFit, 
        TransformationPointFit compareFit) {
        
        if (compareFit == null) {
            return false;
        }
        if (bestFit == null) {
            return true;
        }
                
        float compareNStat = (float)compareFit.getNumberOfMatchedPoints()/
            (float)compareFit.getNMaxMatchable();
        
        if (Float.isInfinite(compareNStat)) {
            return false;
        }
        
        float bestNStat = (float)bestFit.getNumberOfMatchedPoints()/
            (float)bestFit.getNMaxMatchable();
        
        if (Float.isInfinite(bestNStat)) {
            return true;
        }
        
        if ((bestFit.getNumberOfMatchedPoints() == 0) && 
            (compareFit.getNumberOfMatchedPoints() > 0)) {
            return true;
        } else if (compareFit.getNumberOfMatchedPoints() == 0) {
            return false;
        }
        
        double bestMean = bestFit.getMeanDistFromModel();
            
        double compareMean = compareFit.getMeanDistFromModel();
        
        int comp = Double.compare(compareMean, bestMean);
        if (comp < 0) {
            return true;
        } else if (comp > 0) {
            return false;
        }
            
        double bestStdDevMean = bestFit.getStDevFromMean();
        
        double compareStdDevMean = compareFit.getStDevFromMean();
        
        // a smaller std dev from mean is a better fit
        if (Double.compare(compareStdDevMean, bestStdDevMean) < 0) {
        
            return true;
        
        } else if (compareStdDevMean == bestStdDevMean) {
            
            return (Double.compare(compareMean, bestMean) < 0);
        }
        
        return false;
    }

}
