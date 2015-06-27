package algorithms.imageProcessing;

import algorithms.compGeometry.PointPartitioner;
import algorithms.compGeometry.clustering.FixedDistanceGroupFinder;
import static algorithms.imageProcessing.PointMatcher.minTolerance;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import thirdparty.HungarianAlgorithm;

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

                TransformationPointFit fit = calcTransWithRoughGrid(
                    part1, part2, image1CentroidX, image1CentroidY,
                    setsFractionOfImage);
                
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
        
        if (bestFit == null) {
            return null;
        }
        
        log.info("best fit so far: " + bestFit.toString());
        
        //TODO: if solutions are combined, this may need to be done above.
        
        // use either a finer grid search or a downhill simplex to improve the
        // solution which was found coursely within about 10 degrees
        float rot = bestFit.getParameters().getRotationInDegrees();
        int rotStart = (int)rot - 10;
        if (rotStart < 0) {
            rotStart = 360 + rotStart;
        }
        int rotStop = (int)rot + 10;
        if (rotStop > 359) {
            rotStop = rotStop - 360;
        }
        int rotDelta = 1;
        int scaleStart = (int)(0.9 * bestFit.getScale());
        if (scaleStart < 1) {
            scaleStart = 1;
        }
        int scaleStop = (int)(1.1 * bestFit.getScale());
        int scaleDelta = 1;
        boolean setsAreMatched = false;

        log.info(String.format(
            "starting finer grid search with rot=%d to %d and scale=%d to %d", 
            rotStart, rotStop, scaleStart, scaleStop));
        
        boolean chooseHigherResolution = true;
        
        TransformationPointFit fit = calculateTransformationWithGridSearch(
            unmatchedLeftXY, unmatchedRightXY, 
            image1CentroidX, image1CentroidY,
            rotStart, rotStop, rotDelta, scaleStart, scaleStop, scaleDelta,
            setsAreMatched, setsFractionOfImage, chooseHigherResolution);

        if (fitIsBetter(bestFit, fit)) {
            bestFit = fit;
        }
        
        if (bestFit.getNMaxMatchable() == 0) {
            
            int nMaxMatchable = (unmatchedLeftXY.getN() < unmatchedRightXY.getN()) ?
                unmatchedLeftXY.getN() : unmatchedRightXY.getN();
        
            bestFit.setMaximumNumberMatchable(nMaxMatchable);
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

        int nMaxMatchable = (unmatchedLeftXY.getN() < unmatchedRightXY.getN()) ?
            unmatchedLeftXY.getN() : unmatchedRightXY.getN();
        
        bestFit.setMaximumNumberMatchable(nMaxMatchable);
        
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
            calculateEuclideanTransformation(
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
    public TransformationPointFit calculateEuclideanTransformation(
        PairIntArray scene, PairIntArray model,
        int sceneImageCentroidX, int sceneImageCentroidY,
        int modelImageCentroidX, int modelImageCentroidY,
        float setsFractionOfImage) {

        TransformationPointFit fit = calcTransWithRoughGrid(
            scene, model, sceneImageCentroidX, sceneImageCentroidY,
            setsFractionOfImage);

        if (fit == null) {
            return null;
        }

        int nMaxMatchable = (scene.getN() < model.getN()) ? scene.getN()
            : model.getN();

        float fracMatched = (float)fit.getNumberOfMatchedPoints()/(float)nMaxMatchable;

        if (fracMatched < 0.3) {

            // reverse the order to solve for possible scale < 1.

            TransformationPointFit revFit = calcTransWithRoughGrid(
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
    public TransformationPointFit calcTransWithRoughGrid(
        PairIntArray scene, PairIntArray model, int sceneImageCentroidX,
        int sceneImageCentroidY, float setsFractionOfImage) {

        if ((scene == null) || (model == null)) {
            return null;
        }
        if ((scene.getN() < 3) || (model.getN() < 3)) {
            return null;
        }

        int rotStart = 0;
        int rotStop = 359;
        int rotDelta = 10;
        int scaleStart = 1;
        int scaleStop = 5;
        int scaleDelta = 1;
        boolean setsAreMatched = false;
        boolean overrideDefaultTr = false;
        
        TransformationPointFit fit = calculateTransformationWithGridSearch(
            scene, model, sceneImageCentroidX, sceneImageCentroidY,
            rotStart, rotStop, rotDelta, scaleStart, scaleStop, scaleDelta,
            setsAreMatched, setsFractionOfImage, overrideDefaultTr);

        if (fit != null) {
            log.info("best from calculateTransformationWithGridSearch: "
                + fit.toString());
        }

        return fit;
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
     * @param overrideDefaultTr override the default settings to choose
     * the slower, but more accurate translation solver.
     *
     * @return
     */
    public TransformationPointFit calculateTransformationWithGridSearch(
        PairIntArray set1, PairIntArray set2,
        int image1CentroidX, int image1CentroidY,
        int rotStart, int rotStop, int rotDelta,
        int scaleStart, int scaleStop, int scaleDelta,
        boolean setsAreMatched, float setsFractionOfImage, 
        boolean overrideDefaultTr) {
        
        if (rotStart < 0 || rotStart > 359) {
            throw new IllegalArgumentException(
            "rotStart must be between 0 and 359, inclusive");
        }
        if (rotStop < 0 || rotStop > 359) {
            throw new IllegalArgumentException(
            "rotStop must be between 0 and 359, inclusive");
        }
        if (rotDelta < 1) {
            throw new IllegalArgumentException( "rotDelta must be > 0");
        }
        if (!(scaleStart > 0)) {
            throw new IllegalArgumentException("scaleStart must be > 0");
        }
        if (!(scaleStop > 0)) {
            throw new IllegalArgumentException("scaleStop must be > 0");
        }
        if (!(scaleDelta > 0)) {
            throw new IllegalArgumentException("scaleDelta must be > 0");
        }

        double tolTransX = generalTolerance;//4.0f * image1CentroidX * 0.02f;
        double tolTransY = generalTolerance;//4.0f * image1CentroidY * 0.02f;
        if (tolTransX < minTolerance) {
            tolTransX = minTolerance;
        }
        if (tolTransY < minTolerance) {
            tolTransY = minTolerance;
        }

        // rewrite the rotation points into array because start is sometimes
        // higher number than stop in unit circle
        int[] rotation = MiscMath.writeDegreeIntervals(rotStart, rotStop, rotDelta);
        
        Transformer transformer = new Transformer();

        int nMaxMatchable = (set1.getN() < set2.getN()) ? set1.getN()
            : set2.getN();

        int convergence = nMaxMatchable;

        TransformationPointFit bestFit = null;

        TransformationPointFit bestFitForScale = null;

        for (int scale = scaleStart; scale <= scaleStop; scale += scaleDelta) {
            for (int rot : rotation) {
                
                TransformationParameters params;

                if (setsAreMatched) {
                    params = calculateTranslationForMatched(set1, set2,
                        rot*Math.PI/180., scale, image1CentroidX, image1CentroidY);
                } else {
                    //TODO: when refactor, this method already determines
                    //   fit for one branch so reduce redundant code
                    params = calculateTranslationForUnmatched(set1, set2,
                        rot*Math.PI/180., scale, image1CentroidX,
                        image1CentroidY, setsFractionOfImage, overrideDefaultTr);
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
log.info("    " + fit.toString());
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
                        
                        bestFit.setMaximumNumberMatchable(nMaxMatchable);
                        
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

        bestFit.setMaximumNumberMatchable(nMaxMatchable);
        
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
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * @param overrideDefaultTr override the default settings to choose
     * the slower, but more accurate translation solver.
     * @return
     */
    public TransformationParameters calculateTranslationForUnmatched(
        PairIntArray set1, PairIntArray set2, double rotation, double scale,
        int centroidX1, int centroidY1, float setsFractionOfImage,
        boolean overrideDefaultTr) {

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

        int maxNMatchable = (set1.getN() < set2.getN()) ? set1.getN()
            : set2.getN();

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

        boolean useHigherRes = overrideDefaultTr || (maxNMatchable < 41);
        
//TODO: should revise these so that they get decided based upon
// the same points, and if anomalies such as the y dimension being
// much much shorter than the x dimension are present,
// need to decide primarily by the x in that case.
        
        // when there aren't enough points for useful histogram,
        // will make a frequency map of round to integer,
        // and take the peak if its larger than next peak,
        // else, take the average of largest frequencies.

        if (!useHigherRes) {
            
            // use more than normal number of bins:
            int nBins = (int)(3*2*Math.pow(transX.length, 0.3333));

            HistogramHolder hX = Histogram
                .createSimpleHistogram(nBins,
                transX, Errors.populateYErrorsBySqrt(transX));

            HistogramHolder hY = Histogram
                .createSimpleHistogram(nBins,
                transY, Errors.populateYErrorsBySqrt(transY));
            
            if ((hX.getXHist().length < 2) || hY.getXHist().length < 2) {
                
                useHigherRes = true;
                
            } else {

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
                    transX, transY, xr, yr, set2, tolTransX, tolTransY);

                peakTransX = transXY[0];
                peakTransY = transXY[1];
            }
        }

        if (useHigherRes) {

            /*
            Can use my clustering code here from another project for best
            solution of finding the clusters of transX, transY and then the
            cluster with largest number of members among those.

            or can make a quicker solution by deciding an association radius
            for points (that is, a radius around which points are considered
            the same value) and then make a frequency map for same points.
            slightly better than a histogram for small numbers because the
            similar values will not be split into adjacent bins by choice of
            interval start.
            */

            /*
            TODO: may change back to histogram which is O(N).  This is
            O(N^2).
            */
            FixedDistanceGroupFinder groupFinder = null;

            int maxSep = 50;

            List<Set<Integer>> sortedGroupIndexList = null;

            int nIter = 0;
            int nMaxIter = 10;
            
            Set<Integer> mostFreqIndexes = null;
            Set<Integer> prev = null;
            
            float limit = 1.25f * maxNMatchable;
            if (maxNMatchable > 40) {
                limit = 0.95f * maxNMatchable;
            } else if (maxNMatchable > 20) {
                limit = 1.05f * maxNMatchable;
            }
            
            while ((nIter == 0) || 
                ((nIter < nMaxIter) && (mostFreqIndexes.size() > limit))
                && (maxSep > 0)) {
                
                groupFinder = new FixedDistanceGroupFinder(transX, transY);
                groupFinder.findGroupsOfPoints(maxSep);
                sortedGroupIndexList = groupFinder.getDescendingSortGroupList();
                prev = mostFreqIndexes;
                mostFreqIndexes = sortedGroupIndexList.get(0);
                
                maxSep *= 0.75;
                
                nIter++;
            }

            double avgTransX = 0;
            double avgTransY = 0;
            for (Integer index : mostFreqIndexes) {
                avgTransX += transX[index.intValue()];
                avgTransY += transY[index.intValue()];
            }
            avgTransX /= (float)mostFreqIndexes.size();
            avgTransY /= (float)mostFreqIndexes.size();
            peakTransX = (float)avgTransX;
            peakTransY = (float)avgTransY;
            
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

    protected boolean fitIsSimilar(TransformationPointFit fit1,
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
        HistogramHolder hY, float[] transXComb, float[] transYComb,
        float[] xr, float[] yr, PairIntArray set2,
        float tolTransX, float tolTransY) {

        float peakTransX = Float.MAX_VALUE;
        float peakTransY = Float.MAX_VALUE;

        List<Integer> xPeakIndexes = MiscMath.findStrongPeakIndexesDescSort(hX, 
            0.09f);

        List<Integer> yPeakIndexes = MiscMath.findStrongPeakIndexesDescSort(hY, 
            0.09f);
        
        boolean strongXPeak = (xPeakIndexes.size() == 1) ||
            (!xPeakIndexes.isEmpty() &&
            ((hX.getYHistFloat()[xPeakIndexes.get(0).intValue()]/
            hX.getYHistFloat()[xPeakIndexes.get(1).intValue()]
            ) >= 1.5f));

        boolean strongYPeak = (yPeakIndexes.size() == 1) ||
            (!yPeakIndexes.isEmpty() &&
            ((hY.getYHistFloat()[yPeakIndexes.get(0).intValue()]/
            hY.getYHistFloat()[yPeakIndexes.get(1).intValue()]
            ) >= 1.5f));

        boolean fitPeaks = false;

        if (strongXPeak || strongYPeak) {

            if (strongXPeak && strongYPeak) {

                peakTransX = hX.getXHist()[xPeakIndexes.get(0)];

                peakTransY = hY.getXHist()[yPeakIndexes.get(0)];

            } else if (strongXPeak) {

                // form a transY histogram only from the transX which are
                // in this peak bin

                peakTransX = hX.getXHist()[xPeakIndexes.get(0)];
                float dxHalf = (hX.getXHist()[1] - hX.getXHist()[0])/2.f;
                float xMin = peakTransX - dxHalf;
                float xMax = peakTransX + dxHalf;

                HistogramHolder hY2 = formHistogramFromRangeInOther(xMin, xMax,
                    transXComb, transYComb);

                if (hY2.getXHist().length > 1) {
                    try {
                        hY2.plotHistogram("histogram transY", 3);
                    } catch (IOException e) {}

                    List<Integer> yPeakIndexes2 = 
                        MiscMath.findStrongPeakIndexesDescSort(hY2, 0.09f);

                    peakTransY = hY2.getXHist()[yPeakIndexes2.get(0)];
                    
                } else {
             
                    if (!xPeakIndexes.isEmpty() && !yPeakIndexes.isEmpty()) {
                        
                        //TODO: this case may be re-done with more testing
                        
                        peakTransX = hX.getXHist()[xPeakIndexes.get(0)];
              
                        peakTransY = hY.getXHist()[yPeakIndexes.get(0)];
                        
                    } else {
                        
                        fitPeaks = true;
                    }
                }
                
            } else {

                // form a transX histogram only from the transY which are
                // in this peak bin

                peakTransY = hY.getXHist()[yPeakIndexes.get(0)];
                float dyHalf = (hY.getXHist()[1] - hY.getYHist()[0])/2.f;
                float yMin = peakTransY - dyHalf;
                float yMax = peakTransY + dyHalf;

                HistogramHolder hX2 = formHistogramFromRangeInOther(yMin, yMax,
                    transYComb, transXComb);

                if (hX2.getXHist().length > 1) {
                    
                    try {
                        hX2.plotHistogram("histogram transX", 4);
                    } catch (IOException e) {}

                    List<Integer> xPeakIndexes2 = 
                        MiscMath.findStrongPeakIndexesDescSort(hX2, 0.09f);

                    peakTransX = hX2.getXHist()[xPeakIndexes2.get(0)];
                   
                } else {
                    
                    fitPeaks = true;
                }
            }

        } else {

            fitPeaks = true;
        }
        
        if (fitPeaks) {
            
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

    private HistogramHolder formHistogramFromRangeInOther(
        float minA, float maxA, float[] a, float[] b) {

        if (a == null || b == null) {
            throw new IllegalArgumentException("a and b cannot be null");
        }

        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be the same length");
        }

        float[] bSection = new float[b.length];

        int count = 0;
        for (int i = 0; i < a.length; ++i) {
            if ((a[i] >= minA) && (a[i] <= maxA)) {
                bSection[count] = b[i];
                count++;
            }
        }

        bSection = Arrays.copyOf(bSection, count);

        int nBins = (int)(2*Math.pow(count, 0.3333));
        if (count > 200) {
            nBins *= 3;
        }
        
        HistogramHolder hist = Histogram.createSimpleHistogram(nBins, bSection,
            Errors.populateYErrorsBySqrt(bSection));

        return hist;
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

    // ======= code that needs testing and revision the most
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
            -plusMinusTransX, -0.5f * plusMinusTransX,
            -0.25f * plusMinusTransX,
            -0.125f * plusMinusTransX, -0.0625f * plusMinusTransX,
            -0.03125f * plusMinusTransX
        };
        float[] dty = new float[]{
            -plusMinusTransY, -0.5f * plusMinusTransY,
            -0.25f * plusMinusTransY,
            -0.125f * plusMinusTransY, -0.0625f * plusMinusTransY,
            -0.03125f * plusMinusTransY
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
            transX = txSum / (float) (fits.length - 1);
            transY = tySum / (float) (fits.length - 1);

            // "Reflection"
            float txReflect = transX + (alpha
                * (transX - (float) fits[worstFitIdx].getTranslationX()));
            float tyReflect = transY + (alpha
                * (transY - (float) fits[worstFitIdx].getTranslationY()));

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
                        * (transX - (float) fits[worstFitIdx].getTranslationX()));
                    float tyExpansion = transY + (gamma
                        * (transY - (float) fits[worstFitIdx].getTranslationY()));

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
                        * (transX - (float) fits[worstFitIdx].getTranslationX()));
                    float tyContraction = transY + (beta
                        * (transY - (float) fits[worstFitIdx].getTranslationY()));

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
                                = (float) (fits[bestFitIdx].getTranslationX()
                                + (tau
                                * (fits[i].getTranslationX()
                                - fits[bestFitIdx].getTranslationX())));

                            float tyReduction
                                = (float) (fits[bestFitIdx].getTranslationY()
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
    
}
