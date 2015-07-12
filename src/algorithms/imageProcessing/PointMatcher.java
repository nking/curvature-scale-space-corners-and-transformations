package algorithms.imageProcessing;

import algorithms.compGeometry.PointPartitioner;
import static algorithms.imageProcessing.PointMatcher.minTolerance;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairFloat;
import algorithms.util.PairFloatArray;
import algorithms.util.PairFloatArrayUnmodifiable;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonPlotterPNG;
import algorithms.util.RangeInt;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
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
public final class PointMatcher extends AbstractPointMatcher {

    private final Logger log = Logger.getLogger(this.getClass().getName());

    protected static int minTolerance = 5;

    //TODO: this has to be a high number for sets with projection.
    // the solution is sensitive to this value.
    private final float generalTolerance = 8;

    public static float toleranceGridFactor = 4.f;

    protected boolean debug = true;

    /**
     * NOT READY FOR USE
     *
     * @param unmatchedLeftXY
     * @param unmatchedRightXY
     * @param image1Width
     * @param image1Height
     * @param image2Width
     * @param image2Height
     * @param outputMatchedLeftXY
     * @param outputMatchedRightXY
     * @param useLargestToleranceForOutput use the largest tolerance for
     * applying the transformation to point sets during matching.  the largest
     * tolerance is the class variable generalTolerance.
     * If useLargestToleranceForOutput is false, the transformation's best
     * fit is used during matching (which should provide a smaller but more
     * certain matched output).  If this method is used as a precursor to
     * projection (epipolar) solvers of sets that do have projection components,
     * one might prefer to set this to true to allow more to be matched.
     * @return
     */
    public TransformationPointFit performPartitionedMatching(
        PairIntArray unmatchedLeftXY, PairIntArray unmatchedRightXY,
        int image1Width, int image1Height, int image2Width, int image2Height,
        PairIntArray outputMatchedLeftXY, PairIntArray outputMatchedRightXY,
        boolean useLargestToleranceForOutput) {

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
            image1Width, image1Height, image2Width, image2Height,
            allPointsLeftMatched, allPointsRightMatched, 1.0f);

        float allPointsNStat = (float)allPointsFit.getNumberOfMatchedPoints()/
            (float)allPointsFit.getNMaxMatchable();

        log.info("all points set nStat=" + allPointsNStat + " Euclidean fit=" +
            allPointsFit.toString());

        // ====== vertical partitioned matches =======
        TransformationPointFit verticalPartitionedFit =
            performVerticalPartitionedMatching(2,
            unmatchedLeftXY, unmatchedRightXY,
            image1Width, image1Height, image2Width, image2Height);

        /*
        // ====== horizontal partitioned matches =======
        TransformationPointFit horizontalPartitionedFit =
            performHorizontalPartitionedMatching(2,
            unmatchedLeftXY, unmatchedRightXY,
            image1CentroidX, image1CentroidY,
            image2CentroidX, image2CentroidY);
        */

        TransformationPointFit bestFit = allPointsFit;

        if (fitIsBetter(bestFit, verticalPartitionedFit)) {
            //fitIsBetterNStat(bestFit, verticalPartitionedFit)) {

            bestFit = verticalPartitionedFit;
        }

        if (bestFit == null) {
            return null;
        }

        // TODO: compare to horizontal fits when implemented


        Transformer transformer = new Transformer();

        int image1CentroidX = image1Width >> 1;
        int image1CentroidY = image1Height >> 1;

        PairFloatArray transformedLeft = transformer.applyTransformation2(
            bestFit.getParameters(), unmatchedLeftXY, image1CentroidX,
            image1CentroidY);

        float transTolX, transTolY, tolerance;

        if (useLargestToleranceForOutput) {
            tolerance = generalTolerance;
            transTolX = generalTolerance * (float)Math.sqrt(1./2);
            transTolY = transTolX;
        } else {
            transTolX = bestFit.getTranslationXTolerance();
            transTolY = bestFit.getTranslationYTolerance();
            tolerance = (float)Math.sqrt(transTolX*transTolX + transTolY*transTolY);
        }

        float[][] matchIndexesAndDiffs = calculateMatchUsingOptimal(
            transformedLeft, unmatchedRightXY, transTolX, transTolX);

        matchPoints(unmatchedLeftXY, unmatchedRightXY, tolerance,
            matchIndexesAndDiffs, outputMatchedLeftXY, outputMatchedRightXY);

        return bestFit;
    }

    /**
     * given unmatched point sets unmatchedLeftXY and unmatchedRightXY,
     * partitions the data in numberOfPartitions vertically, finds the
     * best match for each combination of vertical partitions as subsets
     * of their parent sets, then finds the best partition solution
     * among those, then refines the solution with a more detailed search
     * using all points.
     *
     * @param numberOfPartitions the number of vertical partitions to make.
     * the maximum value accepted is 3 and minimum is 1.
     * @param unmatchedLeftXY
     * @param unmatchedRightXY
     * @param image1Width
     * @param image1Height
     * @param image2Width
     * @param image2Height
     * @return best fitting transformation between unmatched left and right
     */
    public TransformationPointFit performVerticalPartitionedMatching(
        final int numberOfPartitions,
        PairIntArray unmatchedLeftXY, PairIntArray unmatchedRightXY,
        int image1Width, int image1Height, int image2Width, int image2Height) {

        if (numberOfPartitions > 3) {
            throw new IllegalArgumentException(
            "numberOfPartitions max value is 3");
        }
        if (numberOfPartitions < 1) {
            throw new IllegalArgumentException(
            "numberOfPartitions min value is 1");
        }

        int nXDiv = numberOfPartitions;
        int nYDiv = 1;
        float setsFractionOfImage = 1.f/(float)numberOfPartitions;
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

                // determine fit only with partitioned points

                PairIntArray part1 = vertPartitionedLeft[p1];
                PairIntArray part2 = vertPartitionedRight[p2];

                TransformationPointFit fit = calculateRoughTransformation(
                    part1, part2,
                    image1Width, image1Height, image2Width, image2Height,
                    setsFractionOfImage);

                if (fit == null) {
                    log.fine(Integer.toString(count)
                        + ", p1=" + p1 + " p2=" + p2 + ")  NO SOLN " +
                        " n1=" + part1.getN() + " n2=" + part2.getN());
                    count++;
                    continue;
                }

                int nMaxMatchable = Math.min(part1.getN(), part2.getN());

                boolean converged = hasConverged(fit, nMaxMatchable);

                if (!converged) {

                    fit = refineTheTransformation(
                        fit, part1, part2,
                        image1Width, image1Height, image2Width, image2Height,
                        setsFractionOfImage);
                }

                nStat[count] = (float)fit.getNumberOfMatchedPoints()/
                    (float)fit.getNMaxMatchable();

                log.fine(Integer.toString(count)
                    + ", p1=" + p1 + " p2=" + p2 + ") nStat=" + nStat[count] +
                    " n1=" + part1.getN() + " n2=" + part2.getN() +
                    " Euclidean fit=" + fit.toString());

                vertPartitionedFits[count] = fit;

                if (bestFitIdx == -1) {

                    bestFitIdx = count;

                } else {

                    // if nStat != inf, best has smallest st dev from mean

                    if (fitIsBetter(vertPartitionedFits[bestFitIdx],
                        vertPartitionedFits[count])) {

                        log.fine("****==> fit=" + vertPartitionedFits[count].toString());

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

                    log.fine("similar solutions for idx="
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

                        log.fine("similar solutions for idx="
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

        log.info("best fit so far from partitions: " + bestFit.toString());

        //TODO: if solutions are combined, this may need to be done above.

        // use either a finer grid search or a downhill simplex to improve the
        // solution which was found coursely within about 10 degrees
        float rot = bestFit.getParameters().getRotationInDegrees();
        float rotStart = (int)rot - 10;
        if (rotStart < 0) {
            rotStart = 360 + rotStart;
        }
        float rotStop = (int)rot + 10;
        if (rotStop > 359) {
            rotStop = rotStop - 360;
        }
        float rotDelta = 1;
        float scaleStart = (int)(0.9 * bestFit.getScale());
        if (scaleStart < 1) {
            scaleStart = 1;
        }
        float scaleStop = (int)(1.1 * bestFit.getScale());
        float scaleDelta = 1;

        //TODO: this may need to scale with point density and image dimensions
        int nTransIntervals = 3;

        float transX = bestFit.getParameters().getTranslationX();
        float transY = bestFit.getParameters().getTranslationY();

        float dim = Math.max(bestFit.getTranslationXTolerance(),
            bestFit.getTranslationYTolerance());
        dim *= 2.5;

        float gridWidth = nTransIntervals * dim;
        float gridHeight = gridWidth;

        RangeInt transXStartStop = new RangeInt((int)(transX - (gridWidth/2.f)),
            (int)(transX + (gridWidth/2.f)));
        RangeInt transYStartStop = new RangeInt((int)(transY - (gridHeight/2.f)),
            (int)(transY + (gridHeight/2.f)));

        //TODO: since this is refinement, might use smaller
        int deltaTransX = 4;
        int deltaTransY = 4;

        transXStartStop.resetBoundsIfNeeded(-1*image2Width + 1, image2Width - 1);

        transYStartStop.resetBoundsIfNeeded(-1*image2Height + 1, image2Height - 1);

        float tolTransX = getTolFactor() * deltaTransX;
        float tolTransY = getTolFactor() * deltaTransY;

        log.fine(String.format(
            "starting finer grid search with rot=%d to %d and scale=%d to %d and translation=%d:%d, %d:%d",
            rotStart, rotStop, scaleStart, scaleStop,
            transXStartStop.getStart(), transXStartStop.getStop(),
            transYStartStop.getStart(), transYStartStop.getStop()));

        TransformationPointFit fit = calculateTransformationWithGridSearch(
            unmatchedLeftXY, unmatchedRightXY,
            image1Width, image1Height, image2Width, image2Height,
            rotStart, rotStop, rotDelta, scaleStart, scaleStop, scaleDelta,
            transXStartStop, deltaTransX,
            transYStartStop, deltaTransY,
            tolTransX, tolTransY, setsFractionOfImage);

if (bestFit != null && fit != null) {
log.fine("    partition compare  \n      **==> bestFit=" + bestFit.toString() + "\n           fit=" + fit.toString());
}

        TransformationPointFit[] reevalFits = new TransformationPointFit[2];
        boolean[] fitIsBetter = new boolean[1];
        if ((bestFit != null) && (fit != null) &&
            ((
            ((bestFit.getTranslationXTolerance()/fit.getTranslationXTolerance()) > 2)
            &&
            ((bestFit.getTranslationYTolerance()/fit.getTranslationYTolerance()) > 2))
            ||
            (
            ((bestFit.getTranslationXTolerance()/fit.getTranslationXTolerance()) < 0.5)
            &&
            ((bestFit.getTranslationYTolerance()/fit.getTranslationYTolerance()) < 0.5))
            )) {

            reevaluateFitsForCommonTolerance(bestFit, fit,
                unmatchedLeftXY, unmatchedRightXY, image1Width, image1Height,
                reevalFits, fitIsBetter);

            bestFit = reevalFits[0];
            fit = reevalFits[1];

        } else {

            fitIsBetter[0] = fitIsBetter(bestFit, fit);
        }

if (bestFit != null && fit != null) {
log.fine("    tol corrected partition compare  \n      **==> bestFit=" + bestFit.toString() + "\n           fit=" + fit.toString());
}

if (fitIsBetter[0] && (fit != null)) {
log.fine("    ***** partition bestFit=" + fit.toString());
} else if (bestFit != null) {
log.fine("    ***** partition keeping bestFit=" + bestFit.toString());
}

        if (fitIsBetter[0] && (fit != null)) {
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
     * @param unmatchedLeftXY
     * @param unmatchedRightXY
     * @param image1Width
     * @param image1Height
     * @param outputMatchedLeftXY
     * @param image2Width
     * @param outputMatchedRightXY
     * @param image2Height
     * @param useLargestToleranceForOutput use the largest tolerance for
     * applying the transformation to point sets during matching.  the largest
     * tolerance is the class variable generalTolerance.
     * If useLargestToleranceForOutput is false, the transformation's best
     * fit is used during matching (which should provide a smaller but more
     * certain matched output).  If this method is used as a precursor to
     * projection (epipolar) solvers of sets that do have projection components,
     * one might prefer to set this to true to allow more to be matched.
     * @return best fitting transformation between unmatched points sets
     * left and right
     */
    public TransformationPointFit performMatching0(
        PairIntArray unmatchedLeftXY, PairIntArray unmatchedRightXY,
        int image1Width, int image1Height, int image2Width, int image2Height,
        PairIntArray outputMatchedLeftXY, PairIntArray outputMatchedRightXY,
        boolean useLargestToleranceForOutput) {

        if (unmatchedLeftXY == null || unmatchedLeftXY.getN() < 3) {
            throw new IllegalArgumentException(
            "unmatchedLeftXY cannot be null and must have at least 3 points.");
        }
        if (unmatchedRightXY == null || unmatchedRightXY.getN() < 3) {
            throw new IllegalArgumentException(
            "unmatchedRightXY cannot be null and must have at least 3 points.");
        }

        TransformationPointFit bestFit = performMatching0(
            unmatchedLeftXY, unmatchedRightXY,
            image1Width, image1Height, image2Width, image2Height);

        if (bestFit == null) {
            return null;
        }

        int image1CentroidX = image1Width >> 1;
        int image1CentroidY = image1Height >> 1;

        Transformer transformer = new Transformer();

        PairFloatArray transformedLeft = transformer.applyTransformation2(
            bestFit.getParameters(), unmatchedLeftXY, image1CentroidX,
            image1CentroidY);

        float transTolX, transTolY, tolerance;

        if (useLargestToleranceForOutput) {
            tolerance = generalTolerance;
            transTolX = generalTolerance * (float)Math.sqrt(1./2);
            transTolY = transTolX;
        } else {
            transTolX = bestFit.getTranslationXTolerance();
            transTolY = bestFit.getTranslationYTolerance();
            tolerance = (float)Math.sqrt(transTolX*transTolX + transTolY*transTolY);
        }

        float[][] matchIndexesAndDiffs = calculateMatchUsingOptimal(
            transformedLeft, unmatchedRightXY, transTolX, transTolX);

        matchPoints(unmatchedLeftXY, unmatchedRightXY, tolerance,
            matchIndexesAndDiffs, outputMatchedLeftXY, outputMatchedRightXY);

        int nMaxMatchable = (unmatchedLeftXY.getN() < unmatchedRightXY.getN()) ?
            unmatchedLeftXY.getN() : unmatchedRightXY.getN();

        // rewrite best fit to have the matched number of points
        TransformationPointFit fit2 = new TransformationPointFit(
            bestFit.getParameters(), outputMatchedLeftXY.getN(),
            bestFit.getMeanDistFromModel(), bestFit.getStDevFromMean(),
            bestFit.getTranslationXTolerance(),
            bestFit.getTranslationYTolerance());

        fit2.setMaximumNumberMatchable(nMaxMatchable);

        return fit2;
    }

    /**
     * Given unmatched point sets unmatchedLeftXY and unmatchedRightXY,
     * finds the best Euclidean transformation, then finds the best partition
     * solution among those, then refines the solution with a more detailed
     * search using all points.
     *
     * @param unmatchedLeftXY
     * @param unmatchedRightXY
     * @param image1Width
     * @param image1Height
     * @param image2Width
     * @param image2Height
     * @return best fitting transformation between unmatched left and right
     */
    public TransformationPointFit performMatching0(
        PairIntArray unmatchedLeftXY, PairIntArray unmatchedRightXY,
        int image1Width, int image1Height, int image2Width, int image2Height) {

        //update when finish improving the vertical partition code

        throw new UnsupportedOperationException("not yet implemented");
    }

    /**
     * NOT READY FOR USE
     *
     * @param numberOfPartitions the number of vertical partitions to make.
     * the maximum value accepted is 3.
     * @param unmatchedLeftXY
     * @param unmatchedRightXY
     * @param image1Width
     * @param image1Height
     * @param outputMatchedLeftXY
     * @param image2Width
     * @param outputMatchedRightXY
     * @param image2Height
     * @param useLargestToleranceForOutput use the largest tolerance for
     * applying the transformation to point sets during matching.  the largest
     * tolerance is the class variable generalTolerance.
     * If useLargestToleranceForOutput is false, the transformation's best
     * fit is used during matching (which should provide a smaller but more
     * certain matched output).  If this method is used as a precursor to
     * projection (epipolar) solvers of sets that do have projection components,
     * one might prefer to set this to true to allow more to be matched.
     * @return best fitting transformation between unmatched points sets
     * left and right
     */
    public TransformationPointFit performVerticalPartitionedMatching(
        final int numberOfPartitions,
        PairIntArray unmatchedLeftXY, PairIntArray unmatchedRightXY,
        int image1Width, int image1Height, int image2Width, int image2Height,
        PairIntArray outputMatchedLeftXY, PairIntArray outputMatchedRightXY,
        boolean useLargestToleranceForOutput) {

        if (numberOfPartitions > 3) {
            throw new IllegalArgumentException("numberOfPartitions max value is 3");
        }
        if (numberOfPartitions < 1) {
            throw new IllegalArgumentException("numberOfPartitions min value is 1");
        }

        TransformationPointFit bestFit = performVerticalPartitionedMatching(
            numberOfPartitions, unmatchedLeftXY, unmatchedRightXY,
            image1Width, image1Height, image2Width, image2Height);

        if (bestFit == null) {
            return null;
        }

        int image1CentroidX = image1Width >> 1;
        int image1CentroidY = image1Height >> 1;

        Transformer transformer = new Transformer();

        PairFloatArray transformedLeft = transformer.applyTransformation2(
            bestFit.getParameters(), unmatchedLeftXY, image1CentroidX,
            image1CentroidY);

        float transTolX, transTolY, tolerance;

        if (useLargestToleranceForOutput) {
            tolerance = generalTolerance;
            transTolX = generalTolerance * (float)Math.sqrt(1./2);
            transTolY = transTolX;
        } else {
            transTolX = bestFit.getTranslationXTolerance();
            transTolY = bestFit.getTranslationYTolerance();
            tolerance = (float)Math.sqrt(transTolX*transTolX + transTolY*transTolY);
        }

        float[][] matchIndexesAndDiffs = calculateMatchUsingOptimal(
            transformedLeft, unmatchedRightXY, transTolX, transTolX);

        matchPoints(unmatchedLeftXY, unmatchedRightXY, tolerance,
            matchIndexesAndDiffs, outputMatchedLeftXY, outputMatchedRightXY);

        int nMaxMatchable = (unmatchedLeftXY.getN() < unmatchedRightXY.getN()) ?
            unmatchedLeftXY.getN() : unmatchedRightXY.getN();

        // rewrite best fit to have the matched number of points
        TransformationPointFit fit2 = new TransformationPointFit(
            bestFit.getParameters(), outputMatchedLeftXY.getN(),
            bestFit.getMeanDistFromModel(), bestFit.getStDevFromMean(),
            bestFit.getTranslationXTolerance(),
            bestFit.getTranslationYTolerance());

        fit2.setMaximumNumberMatchable(nMaxMatchable);

        return fit2;
    }

    /**
     * NOT READY FOR USE
     *
     * @param unmatchedLeftXY
     * @param unmatchedRightXY
     * @param image1Width
     * @param image1Height
     * @param outputMatchedLeftXY
     * @param image2Height
     * @param image2Width
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
        int image1Width, int image1Height, int image2Width, int image2Height,
        PairIntArray outputMatchedLeftXY, PairIntArray outputMatchedRightXY,
        float setsFractionOfImage) {

        Transformer transformer = new Transformer();

        PairIntArray part1 = unmatchedLeftXY;
        PairIntArray part2 = unmatchedRightXY;

        TransformationPointFit transFit =
            calculateEuclideanTransformation(
            part1, part2,
            image1Width, image1Height, image2Width, image2Height,
            setsFractionOfImage);

        if (transFit == null) {
            return null;
        }

        // -- filter part1 and part2 to keep only intersection region
        int image1CentroidX = image1Width >> 1;
        int image1CentroidY = image1Height >> 1;
        int image2CentroidX = image2Width >> 1;
        int image2CentroidY = image2Height >> 1;

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

        //TODO: consider using transFit.getTolerance() here
        float transTolX = generalTolerance;

        float tolerance = transTolX * (float)Math.sqrt(1./2);

        //evaluate the fit and store a statistical var: nmatched/nmaxmatchable

        float[][] matchIndexesAndDiffs = calculateMatchUsingOptimal(
            transformedFiltered1, filtered2, transTolX, transTolX);

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

        fit2.setTranslationXTolerance(transTolX);
        fit2.setTranslationYTolerance(transTolX);

        fit2.setMaximumNumberMatchable(nMaxMatchable);

        float nStat = (nMaxMatchable > 0) ?
            (float)part1Matched.getN()/(float)nMaxMatchable : 0;

        log.fine("nStat=" + nStat + " nMaxMatchable=" + nMaxMatchable +
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
            matched2.getN(), avg, stdDev, Float.MAX_VALUE, Float.MAX_VALUE);

        return fit;
    }

    protected TransformationPointFit evaluateFitForUnMatchedTransformedGreedy(
        TransformationParameters params, PairFloatArray unmatched1Transformed,
        PairIntArray unmatched2, float tolTransX, float tolTransY) {

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

                if ((Math.abs(dx) > tolTransX) || (Math.abs(dy) > tolTransY)) {
                    continue;
                }

                float diff = (float)Math.sqrt(dx*dx + dy*dy);

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
            avg, stDev, tolTransX, tolTransY);

        return fit;
    }

    protected TransformationPointFit evaluateFitForUnMatchedTransformedOptimal(
        TransformationParameters params, PairFloatArray unmatched1Transformed,
        PairIntArray unmatched2, float tolTransX, float tolTransY) {

        if (unmatched1Transformed == null) {
            throw new IllegalArgumentException(
            "unmatched1Transformed cannot be null");
        }
        if (unmatched2 == null) {
            throw new IllegalArgumentException(
            "unmatched2 cannot be null");
        }

        int n = unmatched1Transformed.getN();

        float[][] matchedIndexesAndDiffs = calculateMatchUsingOptimal(
            unmatched1Transformed, unmatched2, tolTransX, tolTransY);

        int nMatched = matchedIndexesAndDiffs.length;
        double avg = 0;
        double stDev = 0;

        for (int i = 0; i < nMatched; i++) {
            avg += matchedIndexesAndDiffs[i][2];
        }
        avg = (nMatched == 0) ? Double.MAX_VALUE : (avg / (double)nMatched);

        for (int i = 0; i < nMatched; i++) {

            double d = matchedIndexesAndDiffs[i][2] - avg;

            stDev += (d * d);
        }
        stDev = (nMatched == 0) ? Double.MAX_VALUE :
            (Math.sqrt(stDev/((double)nMatched - 1.)));

        TransformationPointFit fit = new TransformationPointFit(params, nMatched,
            avg, stDev, tolTransX, tolTransY);

        return fit;
    }

    protected TransformationPointFit evaluateFitForUnMatchedOptimal(
        PairFloatArrayUnmodifiable scaledRotatedSet1, float transX, float transY,
        float tolTransX, float tolTransY,
        PairIntArray set2, float scale, float rotationInRadians) {

        if (scaledRotatedSet1 == null) {
            throw new IllegalArgumentException(
            "unmatched1Transformed cannot be null");
        }
        if (set2 == null) {
            throw new IllegalArgumentException(
            "unmatched2 cannot be null");
        }

        float[][] matchedIndexesAndDiffs = calculateMatchUsingOptimal(
            scaledRotatedSet1, transX, transY, tolTransX, tolTransY,
            set2);

        int nMatched = matchedIndexesAndDiffs.length;
        double avg = 0;
        double stDev = 0;

        for (int i = 0; i < nMatched; i++) {
            avg += matchedIndexesAndDiffs[i][2];
        }
        avg = (nMatched == 0) ? Double.MAX_VALUE : (avg / (double)nMatched);

        for (int i = 0; i < nMatched; i++) {

            double d = matchedIndexesAndDiffs[i][2] - avg;

            stDev += (d * d);
        }
        stDev = (nMatched == 0) ? Double.MAX_VALUE :
            (Math.sqrt(stDev/((double)nMatched - 1.)));

        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians(rotationInRadians);
        params.setScale(scale);
        params.setTranslationX(transX);
        params.setTranslationY(transY);

        TransformationPointFit fit = new TransformationPointFit(params, nMatched,
            avg, stDev, tolTransX, tolTransY);

        fit.setMaximumNumberMatchable(Math.max(scaledRotatedSet1.getN(),
            set2.getN()));

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
     * @param transTolXY
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
        float transTolXY,
        float[][] matchedIndexesAndDiffs,
        PairIntArray outputMatched1, PairIntArray outputMatched2) {

        if (matchedIndexesAndDiffs == null) {
            return;
        }

        for (int i = 0; i < matchedIndexesAndDiffs.length; i++) {

            int idx1 = (int)matchedIndexesAndDiffs[i][0];
            int idx2 = (int)matchedIndexesAndDiffs[i][1];
            float diff = matchedIndexesAndDiffs[i][2];

            if (diff < transTolXY) {
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
     * @param transTolXY
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
        float transTolXY, float[][] matchedIndexesAndDiffs,
        PairFloatArray outputMatched1, PairIntArray outputMatched2) {

        for (int i = 0; i < matchedIndexesAndDiffs.length; i++) {

            int idx1 = (int)matchedIndexesAndDiffs[i][0];
            int idx2 = (int)matchedIndexesAndDiffs[i][1];
            float diff = matchedIndexesAndDiffs[i][2];

            if (diff < transTolXY) {
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
     * @param image1Width
     * @param image1Height
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * @param image2Width
     * @param image2Height
     *
     * @return
     */
    public TransformationPointFit calculateEuclideanTransformation(
        PairIntArray scene, PairIntArray model,
        int image1Width, int image1Height, int image2Width, int image2Height,
        float setsFractionOfImage) {

        TransformationPointFit fit = calculateRoughTransformation(
            scene, model, image1Width, image1Height, image2Width, image2Height,
            setsFractionOfImage);

        if (fit == null) {
            return null;
        }

        int nMaxMatchable = (scene.getN() < model.getN()) ? scene.getN()
            : model.getN();

        float fracMatched = (float)fit.getNumberOfMatchedPoints()/(float)nMaxMatchable;

        if (fracMatched < 0.3) {

            // reverse the order to solve for possible scale < 1.

            TransformationPointFit revFit = calculateRoughTransformation(
                model, scene,
                image1Width, image1Height, image2Width, image2Height,
                setsFractionOfImage);

            if (fitIsBetter(fit, revFit) && (revFit != null)) {

                TransformationParameters params = revFit.getParameters();

                // reverse the parameters.
                // needs a reference point in both datasets for direct calculation

                int image1CentroidX = image1Width >> 1;
                int image1CentroidY = image1Height >> 1;
                int image2CentroidX = image2Width >> 1;
                int image2CentroidY = image2Height >> 1;

                MatchedPointsTransformationCalculator tc = new
                    MatchedPointsTransformationCalculator();

                double[] x1y1 = tc.applyTransformation(params,
                    image1CentroidX, image1CentroidY,
                    image2CentroidX, image2CentroidY);

                TransformationParameters revParams = tc.swapReferenceFrames(
                    params, image1CentroidX, image1CentroidY,
                    image2CentroidX, image2CentroidY,
                    x1y1[0], x1y1[1]);

                fit = new TransformationPointFit(revParams,
                    revFit.getNumberOfMatchedPoints(),
                    revFit.getMeanDistFromModel(),
                    revFit.getStDevFromMean(),
                    revFit.getTranslationXTolerance(),
                    revFit.getTranslationYTolerance());
            }
        }

        return fit;
    }

    /**
     * Calculate Euclidean transformation for unmatched points.  As with all
     * methods, the user should try to reduce the number of points to the
     * best determined and distributed throughout the area of intersection
     * while reducing the overall number of points.
     * This method is in progress, but ignoring scale, produces a
     * solution that is accurate within 2 degrees of rotation and
     * 20 pixels in translation of X and the same for translation of Y.
     * The method should be followed by refineTheTransformation.
     *
     * @param scene
     * @param model
     * @param image1Width
     * @param image1Height
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * @param image2Height
     * @param image2Width
     *
     * @return
     */
    public TransformationPointFit calculateRoughTransformation(
        PairIntArray scene, PairIntArray model,
        int image1Width, int image1Height, int image2Width, int image2Height,
        float setsFractionOfImage) {

        if ((scene == null) || (model == null)) {
            return null;
        }
        if ((scene.getN() < 3) || (model.getN() < 3)) {
            return null;
        }

        int nMaxMatchable = Math.min(scene.getN(), model.getN());

        //TODO: consider starting the rotation at 340 and ending at 349 to
        //      solve for most common range first.  the solution returns
        //      as soon as it converges.
        float rotStart = 0;
        float rotStop = 359;
        float rotDelta = 5;

        float scaleStart = 1;
        float scaleStop = 5;
        float scaleDelta = 1;

        int deltaTransX = 15;
        int deltaTransY = 15;

        double densX = ((double)nMaxMatchable/(double)image1Width);
        double densY = ((double)nMaxMatchable/(double)image1Height);
        if ((densX < 0.01) || (densY < 0.01)) {
            rotDelta = 2;
        }
        if ((densX > 0.033) || (densY > 0.033)) {
            // 100 is too high.  takes longer and is much less accurate.
            // presumably takes longer because smaller deltas lead
            // to convergence which is an early exit w/ good solution.
            // 75 leads to results accurate within 2 degrees and
            //     less than 10 pix in translation but the smaller
            //     number of convergences means the runtime is
            //     actually longer than a value of 50.
            deltaTransX = 50;
            deltaTransY = 50;
        }
        float tolTransX = getTolFactor() * deltaTransX;
        float tolTransY = getTolFactor() * deltaTransY;

        TransformationPointFit fit = calculateTransformationWithGridSearch(
            scene, model, image1Width, image1Height, image2Width, image2Height,
            rotStart, rotStop, rotDelta, scaleStart, scaleStop, scaleDelta,
            deltaTransX, deltaTransY, tolTransX, tolTransY,
            setsFractionOfImage);

        if (fit != null) {
            log.fine("best from calculateTransformationWithGridSearch: "
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
     * @param toleranceX
     * @param toleranceY
     * @return a two dimensional array holding the matched indexes and
     * the distances between the model and the point for that pair.
     * each row holds float[]{idx1, idx2, diff}
     */
    public float[][] calculateMatchUsingOptimal(
        PairFloatArray transformed1, PairIntArray set2, float toleranceX,
        float toleranceY) {

        int nPoints1 = transformed1.getN();
        int nPoints2 = set2.getN();

        float[][] diffsAsCost = new float[nPoints1][nPoints2];

        // the algorithm modifies diffsAsCost, so make a copy
        float[][] diffsAsCostCopy = new float[nPoints1][nPoints2];

        // key = indexes i,j; value = diffX, diffY
        Map<PairInt, PairFloat> diffsXY = new HashMap<PairInt, PairFloat>();

        int nWithinTol = 0;

float[][] tempAll = new float[nPoints1][nPoints2];

        for (int i = 0; i < transformed1.getN(); i++) {

            diffsAsCost[i] = new float[nPoints2];
            diffsAsCostCopy[i] = new float[nPoints2];

            float x = transformed1.getX(i);
            float y = transformed1.getY(i);

            for (int j = 0; j < set2.getN(); j++) {

                int x2 = set2.getX(j);
                int y2 = set2.getY(j);

                float diffX = x - x2;
                float diffY = y - y2;

//TODO: temp debugging
tempAll[i][j] = (float)Math.sqrt(diffX*diffX + diffY*diffY);

                if ((Math.abs(diffX) > toleranceX) ||
                    (Math.abs(diffY) > toleranceY)) {

                    diffsAsCost[i][j] = Float.MAX_VALUE;
                    diffsAsCostCopy[i][j] = Float.MAX_VALUE;

                } else {

                    double dist = Math.sqrt(diffX*diffX + diffY*diffY);

                    diffsAsCost[i][j] = (float)dist;
                    diffsAsCostCopy[i][j] = (float)dist;

                    diffsXY.put(new PairInt(i, j), new PairFloat(diffX, diffY));

                    nWithinTol++;
                }
            }
        }

        if (nWithinTol == 0) {
            return new float[0][];
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
            if (idx1 == Float.MAX_VALUE || idx2 == Float.MAX_VALUE) {
                continue;
            }
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }

            PairFloat diffXY = diffsXY.get(new PairInt(idx1, idx2));

            if (diffXY == null) {
                continue;
            }

            if ((diffXY.getX() > toleranceX) || (diffXY.getY() > toleranceY)) {
                continue;
            }

            count++;
        }

        float[][] output = new float[count][];

        if (count == 0) {
            return output;
        }

        count = 0;

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

            PairFloat diffXY = diffsXY.get(new PairInt(idx1, idx2));

            if (diffXY == null) {
                continue;
            }

            if ((diffXY.getX() > toleranceX) || (diffXY.getY() > toleranceY)) {
                continue;
            }

            output[count] = new float[3];
            output[count][0] = idx1;
            output[count][1] = idx2;
            output[count][2] = diffsAsCost[idx1][idx2];

            count++;
        }

        return output;
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
     * @return a two dimensional array holding the matched indexes and
     * the distances between the model and the point for that pair.
     * each row holds float[]{idx1, idx2, diff}
     */
    public float[][] calculateMatchUsingOptimal(
        PairFloatArrayUnmodifiable scaledRotatedSet1,
        float transX, float transY, float toleranceX, float toleranceY,
        PairIntArray set2) {

        int nPoints1 = scaledRotatedSet1.getN();
        int nPoints2 = set2.getN();

        float[][] diffsAsCost = new float[nPoints1][nPoints2];

        // the algorithm modifies diffsAsCost, so make a copy
        float[][] diffsAsCostCopy = new float[nPoints1][nPoints2];

        // key = indexes i,j; value = diffX, diffY
        Map<PairInt, PairFloat> diffsXY = new HashMap<PairInt, PairFloat>();

        int nWithinTol = 0;

        for (int i = 0; i < scaledRotatedSet1.getN(); i++) {

            diffsAsCost[i] = new float[nPoints2];
            diffsAsCostCopy[i] = new float[nPoints2];

            float x = scaledRotatedSet1.getX(i) + transX;
            float y = scaledRotatedSet1.getY(i) + transY;

            for (int j = 0; j < set2.getN(); j++) {

                int x2 = set2.getX(j);
                int y2 = set2.getY(j);

                float diffX = x - x2;
                float diffY = y - y2;

                if ((Math.abs(diffX) > toleranceX) ||
                    (Math.abs(diffY) > toleranceY)) {

                    diffsAsCost[i][j] = Float.MAX_VALUE;
                    diffsAsCostCopy[i][j] = Float.MAX_VALUE;

                } else {

                    double dist = Math.sqrt(diffX*diffX + diffY*diffY);

                    diffsAsCost[i][j] = (float)dist;
                    diffsAsCostCopy[i][j] = (float)dist;

                    diffsXY.put(new PairInt(i, j), new PairFloat(diffX, diffY));

                    nWithinTol++;
                }
            }
        }

        if (nWithinTol == 0) {
            return new float[0][];
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
            if (idx1 == Float.MAX_VALUE || idx2 == Float.MAX_VALUE) {
                continue;
            }
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }

            PairFloat diffXY = diffsXY.get(new PairInt(idx1, idx2));

            if (diffXY == null) {
                continue;
            }

            if ((diffXY.getX() > toleranceX) || (diffXY.getY() > toleranceY)) {
                continue;
            }

            count++;
        }

        float[][] output = new float[count][];

        if (count == 0) {
            return output;
        }

        count = 0;

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

            PairFloat diffXY = diffsXY.get(new PairInt(idx1, idx2));

            if (diffXY == null) {
                continue;
            }

            if ((diffXY.getX() > toleranceX) || (diffXY.getY() > toleranceY)) {
                continue;
            }

            output[count] = new float[3];
            output[count][0] = idx1;
            output[count][1] = idx2;
            output[count][2] = diffsAsCost[idx1][idx2];

            count++;
        }

        return output;
    }

    /**
     * Calculate for unmatched points the Euclidean transformation to transform
     * set1 into the reference frame of set2.
     *
     * @param set1
     * @param set2
     * @param image1Width
     * @param image1Height
     * @param image2Width
     * @param rotStart start of rotation search in degrees
     * @param rotStop stop (exclusive) or rotations search in degrees
     * @param image2Height
     * @param rotDelta change in rotation to add to reach next step in rotation
     *     search in degrees
     * @param scaleStart
     * @param scaleStop
     * @param scaleDelta
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
        int image1Width, int image1Height, int image2Width, int image2Height,
        float rotStart, float rotStop, float rotDelta,
        float scaleStart, float scaleStop, float scaleDelta,
        int deltaTransX, int deltaTransY, float tolTransX, float tolTransY,
        float setsFractionOfImage) {

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
        if (scaleDelta < 0) {
            throw new IllegalArgumentException("scaleDelta must be >= 0");
        }

        // rewrite the rotation points into array because start is sometimes
        // higher number than stop in unit circle
        float[] rotation = MiscMath.writeDegreeIntervals(rotStart, rotStop,
            rotDelta);

        int nMaxMatchable = (set1.getN() < set2.getN()) ? set1.getN()
            : set2.getN();

        TransformationPointFit bestFit = null;

        List<TransformationPointFit> similarToBestFit =
            new ArrayList<TransformationPointFit>();

        TransformationPointFit bestFitForScale = null;

        for (float scale = scaleStart; scale <= scaleStop; scale += scaleDelta) {
            for (float rot : rotation) {

                float rotationInRadians = (float)(rot*Math.PI/180.f);

                TransformationPointFit fit = calculateTranslationForUnmatched0(
                    set1, set2,
                    image1Width, image1Height, image2Width, image2Height,
                    rotationInRadians, scale,
                    deltaTransX, deltaTransY, tolTransX, tolTransY,
                    setsFractionOfImage);
if (fit != null) {
log.fine("      rot reeval fit and bestFit.  fit=" + fit.toString());
}
                TransformationPointFit[] reevalFits = new TransformationPointFit[2];
                boolean[] fitIsBetter = new boolean[1];

                reevaluateFitsForCommonTolerance(bestFit, fit,
                    set1, set2, image1Width, image1Height, reevalFits,
                    fitIsBetter);

                bestFit = reevalFits[0];
                fit = reevalFits[1];

if ((bestFit != null) && (fit != null)) {
log.fine("    rot compare  \n      **==> bestFit=" + bestFit.toString() + "\n           fit=" + fit.toString());
} else if (fit != null) {
    log.fine("   rot compare fit=" + fit.toString());
}
//TODO: temporary debugging:
if (!fitIsBetter[0] && (bestFit != null)) {
log.fine("**==> rot keeping bestFit=" + bestFit.toString());
}

                int areSimilar = fitsAreSimilarWithDiffParameters(bestFit, fit);

                //0==same fit;  1==similar fits;  -1==different fits
                if (areSimilar == 0) {
                    //no need to recheck for convergence or change bestFit
                    continue;
                } else if (areSimilar == 1) {
                    log.fine("fit was similar to bestFit");
                    if (similarToBestFit.isEmpty()) {
                        similarToBestFit.add(bestFit);
                    }
                    similarToBestFit.add(fit);
                }

                if (fitIsBetter[0] && (fit != null)) {

                    log.fine("**==> rot fit=" + fit.toString());

                    if (areSimilar == -1) {
                        log.fine("clear similarToBestFit");
                        similarToBestFit.clear();
                    }

                    bestFit = fit;

                    boolean converged = hasConverged(bestFit, nMaxMatchable);

                    if (converged) {

                        log.fine("** converged");

                        similarToBestFit.add(0, bestFit);
                        for (TransformationPointFit fit2 : similarToBestFit) {
                            if (fit2.getParameters().getRotationInRadians() > 2.*Math.PI) {
                                float rot2 = fit2.getParameters().getRotationInRadians();
                                while (rot2 >= 2*Math.PI) {
                                    rot2 -= 2*Math.PI;
                                }
                                fit2.getParameters().setRotationInRadians(rot2);
                            }
                            fit2.setMaximumNumberMatchable(nMaxMatchable);
                        }

                        bestFit = finerGridSearchToDistinguishBest(
                            similarToBestFit, set1, set2,
                            image1Width, image1Height, image2Width, image2Height,
                            setsFractionOfImage);

                        log.info("**== > deciding among " + similarToBestFit.size() + " similar fits:");
                        for (TransformationPointFit sFit : similarToBestFit) {
                            log.fine("  sFit=" + sFit.toString());
                        }

                        log.fine("      bestFit=" + bestFit.toString());

                        return bestFit;
                    }
                }
            }

            if (fitIsBetter(bestFitForScale, bestFit) && (bestFit != null)) {

                bestFitForScale = bestFit;

            } else {

                log.fine("previous scale solution was better, so end scale iter");

                // revert to previous scale
                bestFit = bestFitForScale;

                //TODO: revisit this with tests
                // scale was probably smaller so return best solution
                break;
            }
        }

        if (bestFit != null) {
            similarToBestFit.add(0, bestFit);
        }
        for (TransformationPointFit fit : similarToBestFit) {
            if (fit.getParameters().getRotationInRadians() > 2.*Math.PI) {
                float rot = fit.getParameters().getRotationInRadians();
                while (rot >= 2*Math.PI) {
                    rot -= 2*Math.PI;
                }
                fit.getParameters().setRotationInRadians(rot);
            }
            fit.setMaximumNumberMatchable(nMaxMatchable);
        }

        bestFit = finerGridSearchToDistinguishBest(similarToBestFit,
            set1, set2,
            image1Width, image1Height, image2Width, image2Height,
            setsFractionOfImage);

        log.info("**= => deciding among " + similarToBestFit.size() + " similar fits:");
        for (TransformationPointFit sFit : similarToBestFit) {
            log.fine("  sFit=" + sFit.toString());
        }
        if (bestFit != null) {
            bestFit.setMaximumNumberMatchable(nMaxMatchable);
            log.fine("      bestFit=" + bestFit.toString());
        }

        return bestFit;
    }

    /**
     * Calculate for matched points the Euclidean transformation to transform
     * set1 into the reference frame of set2.
     * @param set1
     * @param set2
     * @param image1Width
     * @param image1Height
     * @param image2Width
     * @param rotStart start of rotation search in degrees
     * @param rotStop stop (exclusive) or rotations search in degrees
     * @param image2Height
     * @param rotDelta change in rotation to add to reach next step in rotation
     *     search in degrees
     * @param scaleStart
     * @param scaleStop
     * @param scaleDelta
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     *
     * @return
     */
    public TransformationPointFit calculateTransformationWithGridSearchForMatched(
        PairIntArray set1, PairIntArray set2,
        int image1Width, int image1Height, int image2Width, int image2Height,
        float rotStart, float rotStop, float rotDelta,
        float scaleStart, float scaleStop, float scaleDelta,
        float setsFractionOfImage) {

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
        if (scaleDelta < 0) {
            throw new IllegalArgumentException("scaleDelta must be >= 0");
        }

        // rewrite the rotation points into array because start is sometimes
        // higher number than stop in unit circle
        float[] rotation = MiscMath.writeDegreeIntervals(rotStart, rotStop,
            rotDelta);

        Transformer transformer = new Transformer();

        int nMaxMatchable = (set1.getN() < set2.getN()) ? set1.getN()
            : set2.getN();

        TransformationPointFit bestFit = null;

        List<TransformationPointFit> similarToBestFit =
            new ArrayList<TransformationPointFit>();

        TransformationPointFit bestFitForScale = null;

        for (float scale = scaleStart; scale <= scaleStop; scale += scaleDelta) {
            for (float rot : rotation) {

                float rotationInRadians = (float)(rot*Math.PI/180.f);

                TransformationParameters params =
                    calculateTranslationForMatched(set1, set2,
                    rotationInRadians, scale,
                    image1Width, image1Height, image2Width, image2Height);

                int image1CentroidX = image1Width >> 1;
                int image1CentroidY = image1Height >> 1;

                PairFloatArray allPoints1Tr = transformer.applyTransformation(
                    params, image1CentroidX, image1CentroidY, set1);

                TransformationPointFit fit = evaluateFitForMatchedTransformed(
                    params, allPoints1Tr, set2);

                //0==same fit;  1==similar fits;  -1==different fits
                int areSimilar = fitsAreSimilarWithDiffParameters(bestFit, fit);
                if (areSimilar == 0) {
                    //no need to recheck for convergence or change bestFit
                    continue;
                } else if (areSimilar == 1) {
                    log.fine("fit was similar to bestFit");
                    if (similarToBestFit.isEmpty()) {
                        similarToBestFit.add(bestFit);
                    }
                    similarToBestFit.add(fit);
                }

                boolean fitIsBetter = fitIsBetter(bestFit, fit);

                if (fitIsBetter && (fit != null)) {

                    if (fit != null) {
                        log.fine("**==> fit=" + fit.toString());
                    }

                    if (areSimilar == -1) {
                        log.fine("clear similarToBestFit");
                        similarToBestFit.clear();
                    }

                    bestFit = fit;

                    boolean converged = hasConverged(bestFit, nMaxMatchable);

                    if (converged) {

                        log.fine("** converged");

                        if (bestFit != null) {
                            similarToBestFit.add(0, bestFit);
                        }
                        for (TransformationPointFit fit2 : similarToBestFit) {
                            if (fit2.getParameters().getRotationInRadians() > 2.*Math.PI) {
                                float rot2 = fit2.getParameters().getRotationInRadians();
                                while (rot2 >= 2*Math.PI) {
                                    rot2 -= 2*Math.PI;
                                }
                                fit2.getParameters().setRotationInRadians(rot2);
                            }
                            fit2.setMaximumNumberMatchable(nMaxMatchable);
                        }

                        bestFit = finerGridSearchToDistinguishBestForMatched(
                            similarToBestFit, set1, set2,
                            image1Width, image1Height, image2Width, image2Height,
                            setsFractionOfImage);

                        log.info("** ==> deciding among " + similarToBestFit.size() + " similar fits:");
                        for (TransformationPointFit sFit : similarToBestFit) {
                            log.fine("  sFit=" + sFit.toString());
                        }
                        log.fine("      bestFit=" + bestFit.toString());

                        return bestFit;

                    } else {

                        if (bestFit.getNumberOfMatchedPoints() == 0) {
                            continue;
                        }

                        /*
                        TODO:
                        this might be better to perform at the end of the
                        method right before returning the best result
                        */

                        int nIntervals = 3;

                        TransformationPointFit fit2 = finerGridSearchForMatched(
                            nIntervals, bestFit, set1, set2,
                            image1Width, image1Height, image2Width, image2Height,
                            setsFractionOfImage
                        );

                        //0==same fit;  1==similar fits;  -1==different fits
                        int areSimilar2 = fitsAreSimilarWithDiffParameters(bestFit, fit2);
                        if (areSimilar2 == 1) {
                            log.fine("fit was similar to bestFit");
                            if (similarToBestFit.isEmpty()) {
                                similarToBestFit.add(bestFit);
                            }
                            similarToBestFit.add(fit2);
                        }

                        boolean fitIsBetter2 = fitIsBetter(bestFit, fit2);

                        if (fitIsBetter2) {

                            log.fine("***==> fit=" + fit2.toString());

                            if (areSimilar == -1) {
                                log.fine("clear similarToBestFit");
                                similarToBestFit.clear();
                            }

                            bestFit = fit2;
                        }
                    }
                }
            }

            if (fitIsBetter(bestFitForScale, bestFit) && (bestFit != null)) {

                bestFitForScale = bestFit;

            } else {

                log.fine("previous scale solution was better, so end scale iter");

                // revert to previous scale
                bestFit = bestFitForScale;

                //TODO: revisit this with tests
                // scale was probably smaller so return best solution
                break;
            }
        }

        if (bestFit != null) {
            similarToBestFit.add(0, bestFit);
        }
        for (TransformationPointFit fit : similarToBestFit) {
            if (fit.getParameters().getRotationInRadians() > 2.*Math.PI) {
                float rot = fit.getParameters().getRotationInRadians();
                while (rot >= 2*Math.PI) {
                    rot -= 2*Math.PI;
                }
                fit.getParameters().setRotationInRadians(rot);
            }
            fit.setMaximumNumberMatchable(nMaxMatchable);
        }

        bestFit = finerGridSearchToDistinguishBestForMatched(similarToBestFit,
            set1, set2,
            image1Width, image1Height, image2Width, image2Height,
            setsFractionOfImage);

        log.info("**==> deciding among " + similarToBestFit.size() + " similar fits:");
        for (TransformationPointFit sFit : similarToBestFit) {
            log.fine("  sFit=" + sFit.toString());
        }
        if (bestFit != null) {
            bestFit.setMaximumNumberMatchable(nMaxMatchable);
            log.fine("      bestFit=" + bestFit.toString());
        }
        return bestFit;
    }

    /**
     * Calculate for unmatched points the Euclidean transformation to transform
     * set1 into the reference frame of set2.
     *
     * @param set1
     * @param set2
     * @param image1Width
     * @param image1Height
     * @param image2Width
     * @param rotStart start of rotation search in degrees
     * @param rotStop stop (exclusive) or rotations search in degrees
     * @param image2Height
     * @param rotDelta change in rotation to add to reach next step in rotation
     *     search in degrees
     * @param scaleStart
     * @param scaleStop
     * @param scaleDelta
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
        int image1Width, int image1Height, int image2Width, int image2Height,
        float rotStart, float rotStop, float rotDelta,
        float scaleStart, float scaleStop, float scaleDelta,
        RangeInt transXStartStop,  int deltaTransX,
        RangeInt transYStartStop, int deltaTransY,
        float tolTransX, float tolTransY,
        float setsFractionOfImage) {

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
        if (scaleDelta < 0) {
            throw new IllegalArgumentException("scaleDelta must be >= 0");
        }

        // rewrite the rotation points into array because start is sometimes
        // higher number than stop in unit circle
        float[] rotation = MiscMath.writeDegreeIntervals(rotStart, rotStop,
            rotDelta);

        int nMaxMatchable = (set1.getN() < set2.getN()) ? set1.getN()
            : set2.getN();

        TransformationPointFit bestFit = null;

        List<TransformationPointFit> similarToBestFit =
            new ArrayList<TransformationPointFit>();

        TransformationPointFit bestFitForScale = null;

        for (float scale = scaleStart; scale <= scaleStop; scale += scaleDelta) {
            for (float rot : rotation) {

                float rotationInRadians = (float)(rot*Math.PI/180.f);

                TransformationPointFit fit = calculateTranslationForUnmatched0(
                    set1, set2,
                    image1Width, image1Height, image2Width, image2Height,
                    rotationInRadians, scale,
                    transXStartStop, deltaTransX, transYStartStop, deltaTransY,
                    tolTransX, tolTransY,
                    setsFractionOfImage);

                //0==same fit;  1==similar fits;  -1==different fits
                int areSimilar = fitsAreSimilarWithDiffParameters(bestFit, fit);

String dbg = String.format("    trying scale=%f, rot=%f", scale, rot);
if (fit != null) {
    dbg = dbg + "  fit=" + fit.toString();
}
dbg = dbg + " areSimilar=" + Integer.valueOf(areSimilar);
log.fine(dbg);

                if (areSimilar == 0) {
                    //no need to recheck for convergence or change bestFit
                    continue;
                } else if (areSimilar == 1) {
                    log.fine("fit was similar to bestFit");
                    if (similarToBestFit.isEmpty()) {
                        similarToBestFit.add(bestFit);
                    }
                    similarToBestFit.add(fit);
                }

                TransformationPointFit[] reevalFits = new TransformationPointFit[2];
                boolean[] fitIsBetter = new boolean[1];

                reevaluateFitsForCommonTolerance(bestFit, fit,
                    set1, set2, image1Width, image1Height, reevalFits,
                    fitIsBetter);

                bestFit = reevalFits[0];
                fit = reevalFits[1];

if ((bestFit != null) && (fit != null)) {
log.fine("    rot0 compare  \n      **==> bestFit=" + bestFit.toString() + "\n           fit=" + fit.toString());
} else if (fit != null) {
    log.fine("   rot0 compare fit=" + fit.toString());
}
//TODO: temporary debugging:
if (!fitIsBetter[0] && (bestFit != null)) {
log.fine("**==> rot0 keeping bestFit=" + bestFit.toString());
}

                if (fitIsBetter[0] && (fit != null)) {

                    log.fine("**==> fit=" + fit.toString());

                    if (areSimilar == -1) {
                        log.fine("clear similarToBestFit");
                        similarToBestFit.clear();
                    }

                    bestFit = fit;

                    boolean converged = hasConverged(bestFit, nMaxMatchable);

                    if (converged) {

                        log.fine("** converged for fit=" + bestFit.toString());

                        similarToBestFit.add(0, bestFit);
                        for (TransformationPointFit fit2 : similarToBestFit) {
                            if (fit2.getParameters().getRotationInRadians() > 2.*Math.PI) {
                                float rot2 = fit2.getParameters().getRotationInRadians();
                                while (rot2 >= 2*Math.PI) {
                                    rot2 -= 2*Math.PI;
                                }
                                fit2.getParameters().setRotationInRadians(rot2);
                            }
                            fit2.setMaximumNumberMatchable(nMaxMatchable);
                        }

                        bestFit = finerGridSearchToDistinguishBest(
                            similarToBestFit, set1, set2,
                            image1Width, image1Height, image2Width, image2Height,
                            setsFractionOfImage);

                        log.info("**====> deciding among " + similarToBestFit.size() + " similar fits:");
                        for (TransformationPointFit sFit : similarToBestFit) {
                            log.fine("  sFit=" + sFit.toString());
                        }
                        log.fine("      bestFit=" + bestFit.toString());

                        return bestFit;

                    } else {

                        if (bestFit.getNumberOfMatchedPoints() == 0) {
                            continue;
                        }

/*
                        int nIntervals = 3;


                        TransformationPointFit fit2 = finerGridSearch(
                            nIntervals, bestFit, set1, set2,
                            image1Width, image1Height, image2Width, image2Height,
                            setsFractionOfImage
                        );

                        //0==same fit;  1==similar fits;  -1==different fits
                        int areSimilar2 = fitsAreSimilarWithDiffParameters(bestFit, fit2);
                        if (areSimilar2 == 1) {
                            log.fine("fit was similar to bestFit");
                            if (similarToBestFit.isEmpty()) {
                                similarToBestFit.add(bestFit);
                            }
                            similarToBestFit.add(fit2);
                        }

                        boolean fitIsBetter2 = fitIsBetter(bestFit, fit2);

                        if (fitIsBetter2) {

                            log.fine("***==> fit=" + fit2.toString());

                            if (areSimilar == -1) {
                                log.fine("clear similarToBestFit");
                                similarToBestFit.clear();
                            }

                            bestFit = fit2;
                        }
*/
                    }
                }
            }

            if (fitIsBetter(bestFitForScale, bestFit) && (bestFit != null)) {

                bestFitForScale = bestFit;

                log.fine("   ==> bestFitForScale=" + bestFitForScale.toString());

            } else {

                log.fine("previous scale solution was better, so end scale iter");

                // revert to previous scale
                bestFit = bestFitForScale;

                //TODO: revisit this with tests
                // scale was probably smaller so return best solution
                break;
            }
        }

        if (bestFit != null) {
            similarToBestFit.add(0, bestFit);
        }
        for (TransformationPointFit fit : similarToBestFit) {
            if (fit.getParameters().getRotationInRadians() > 2.*Math.PI) {
                float rot = fit.getParameters().getRotationInRadians();
                while (rot >= 2*Math.PI) {
                    rot -= 2*Math.PI;
                }
                fit.getParameters().setRotationInRadians(rot);
            }
            fit.setMaximumNumberMatchable(nMaxMatchable);
        }

        bestFit = finerGridSearchToDistinguishBest(similarToBestFit,
            set1, set2,
            image1Width, image1Height, image2Width, image2Height,
            setsFractionOfImage);

        log.info("**==>  deciding among " + similarToBestFit.size() + " similar fits:");
        for (TransformationPointFit sFit : similarToBestFit) {
            log.fine("  sFit=" + sFit.toString());
        }

        if (bestFit != null) {
            bestFit.setMaximumNumberMatchable(nMaxMatchable);
            log.fine("      bestFit=" + bestFit.toString());
        }

        return bestFit;
    }

    protected TransformationPointFit[] evaluateTranslationsOverGrid(
        PairIntArray set1, PairIntArray set2,
        final int image1Width, final int image1Height, final int image2Width,
        final int image2Height,
        final float rotationInRadians, final float scale,
        RangeInt transXStartStop, int transXDelta,
        RangeInt transYStartStop, int transYDelta,
        float tolTransX, float tolTransY,
        final boolean setsAreMatched, final float setsFractionOfImage,
        final int numberOfBestToReturn) {

        /*   _____
            |     |
            |_____| largest negative or positive translationX of set1 is the width of set2
                  _____
                 |     |
                 |_____|
        */

        if (transXStartStop.getStart() < -1*image2Width) {
            throw new IllegalArgumentException(
            "transXStart is out of bounds.. transXStart=" + transXStartStop.getStart());
        }
        if (transXStartStop.getStart() > image2Width) {
            throw new IllegalArgumentException(
            "transXStart is out of bounds.. transXStart=" + transXStartStop.getStart());
        }
        if (transXStartStop.getStop() < -1*image2Width) {
            throw new IllegalArgumentException(
            "transXStop is out of bounds.. transXStop=" + transXStartStop.getStop());
        }
        if (transXStartStop.getStop() > image2Width) {
            throw new IllegalArgumentException(
            "transXStop is out of bounds.. transXStop=" + transXStartStop.getStop());
        }

        if (transYStartStop.getStart() < -1*image2Height) {
            throw new IllegalArgumentException(
            "transYStart is out of bounds.. transYStart=" + transYStartStop.getStart());
        }
        if (transYStartStop.getStart() > image2Height) {
            throw new IllegalArgumentException(
            "transYStart is out of bounds.. transYStart=" + transYStartStop.getStart());
        }
        if (transYStartStop.getStop() < -1*image2Height) {
            throw new IllegalArgumentException(
            "transYStop is out of bounds.. transYStop=" + transYStartStop.getStop());
        }
        if (transYStartStop.getStop() > image2Height) {
            throw new IllegalArgumentException(
            "transYStop is out of bounds.. transYStop=" + transYStartStop.getStop());
        }

        if (transXDelta < 1) {
            throw new IllegalArgumentException(
            "transXDelta must be greater than 0");
        }
        if (rotationInRadians < 0 || rotationInRadians > 359) {
            throw new IllegalArgumentException(
            "rotation must be between 0 and 359, inclusive");
        }
        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");

            log.severe("scale cannot be smaller than 1");

            return null;
        }

        int nMaxMatchable = (set1.getN() < set2.getN()) ?
            set1.getN() : set2.getN();

        if (nMaxMatchable == 0) {
            return null;
        }

        int nTranslations =
            (((transXStartStop.getStop() - transXStartStop.getStart())/transXDelta) + 1) *
            (((transYStartStop.getStop() - transYStartStop.getStart())/transYDelta) + 1);

        if (nTranslations == 0) {
            return null;
        }

        Transformer transformer = new Transformer();

        int image1CentroidX = image1Width >> 1;
        int image1CentroidY = image1Height >> 1;

        int count = 0;

        FixedSizeSortedVector<TransformationPointFit> fits =
            new FixedSizeSortedVector<>(numberOfBestToReturn,
                TransformationPointFit.class);

        for (float transX = transXStartStop.getStart();
            transX <= transXStartStop.getStop(); transX += transXDelta) {

            for (float transY = transYStartStop.getStart();
                transY <= transYStartStop.getStop(); transY += transYDelta) {

                float tx0 = transX + (transXDelta/2);
                float ty0 = transY + (transYDelta/2);

                TransformationParameters params = new TransformationParameters();
                params.setRotationInRadians(rotationInRadians);
                params.setScale(scale);
                params.setTranslationX(tx0);
                params.setTranslationY(ty0);

                PairFloatArray allPoints1Tr = transformer.applyTransformation(
                    params, image1CentroidX, image1CentroidY, set1);

                TransformationPointFit fit;

                if (setsAreMatched) {

                    fit = evaluateFitForMatchedTransformed(params,
                        allPoints1Tr, set2);

                } else {

                    // default is to use greedy matching but use optimal for small sets
                    if (nMaxMatchable < 30) {

                        fit = evaluateFitForUnMatchedTransformedOptimal(params,
                            allPoints1Tr, set2, tolTransX, tolTransY);

                    } else {

                        fit = evaluateFitForUnMatchedTransformedGreedy(params,
                        //fit = evaluateFitForUnMatchedTransformedOptimal(params,
                            allPoints1Tr, set2, tolTransX, tolTransY);
                    }
                }

                fits.add(fit);

                count++;
            }
        }

        return fits.getArray();
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
     * This method is in progress...
     *
     * @param set1 set of points from image 1 to match to image2.
     * @param set2 set of points from image 2 to be matched with image 1
     * @param rotationInRadians given in radians with value between 0 and 2*pi,
     * exclusive
     * @param scale
     * @param image1Width width of image 1, used to derive range of x
     * translations in case centroidX1 is ever used as a non-center reference.
     * @param image1Height height of image 1, used to derive range of y
     * translations in case centroidY1 is ever used as a non-center reference.
     * @param image2Width width of image 2, used to derive range of x
     * translations
     * @param image2Height height of image 2, used to derive range of y
     * translations
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * @return
     */
     public TransformationPointFit calculateTranslationForUnmatched0(
        PairIntArray set1, PairIntArray set2,
        int image1Width, int image1Height, int image2Width, int image2Height,
        float rotationInRadians, float scale,
        int deltaTransX, int deltaTransY, float tolTransX, float tolTransY,
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

        if (scale < 1) {

            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");

            log.severe("scale cannot be smaller than 1");

            return null;
        }

        RangeInt transXStartStop = new RangeInt(-1*image2Width + 1,
            image2Width - 1);
        RangeInt transYStartStop = new RangeInt(-1*image2Height + 1,
            image2Height - 1);

        TransformationPointFit fit = calculateTranslationForUnmatched0(
            set1, set2,
            image1Width, image1Height, image2Width, image2Height,
            rotationInRadians, scale,
            transXStartStop, deltaTransX, transYStartStop, deltaTransY,
            tolTransX, tolTransY,
            setsFractionOfImage);

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
     * This method is in progress...
     *
     * @param set1 set of points from image 1 to match to image2.
     * @param set2 set of points from image 2 to be matched with image 1
     * @param rotationInRadians given in radians with value between 0 and 2*pi,
     * exclusive
     * @param scale
     * @param transXStartStop
     * @param transYStartStop
     * @param image1Width width of image 1, used to derive range of x
     * translations in case centroidX1 is ever used as a non-center reference.
     * @param transYStop
     * @param transYStart
     * @param nTransIntervals
     * @param image1Height height of image 1, used to derive range of y
     * translations in case centroidY1 is ever used as a non-center reference.
     * @param image2Width width of image 2, used to derive range of x
     * translations
     * @param image2Height height of image 2, used to derive range of y
     * translations
     * @param setsFractionOfImage the fraction of their images that set 1
     * and set2 were extracted from. If set1 and set2 were derived from the
     * images without using a partition method, this is 1.0, else if the
     * quadrant partitioning was used, this is 0.25.  The variable is used
     * internally in determining histogram bin sizes for translation.
     * @return
     */
    public TransformationPointFit calculateTranslationForUnmatched0(
        PairIntArray set1, PairIntArray set2,
        int image1Width, int image1Height, int image2Width, int image2Height,
        float rotationInRadians, float scale,
        RangeInt transXStartStop, int transDeltaX,
        RangeInt transYStartStop, int transDeltaY,
        float tolTransX, float tolTransY,
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
        if (transDeltaX < 1) {
            throw new IllegalArgumentException("transDeltaX must be greater than 0");
        }
        if (transDeltaY < 1) {
            throw new IllegalArgumentException("transDeltaY must be greater than 0");
        }
        if (tolTransX < 0) {
            throw new IllegalArgumentException("tolTransX must be greater than -1");
        }
        if (tolTransY < 0) {
            throw new IllegalArgumentException("tolTransY must be greater than -1");
        }

        if (scale < 1) {

            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");

            log.severe("scale cannot be smaller than 1");

            return null;
        }

        int image1CentroidX = image1Width >> 1;
        int image1CentroidY = image1Height >> 1;

        PairFloatArrayUnmodifiable scaledRotatedSet1 = scaleAndRotate(set1,
            rotationInRadians, scale, image1CentroidX, image1CentroidY);

        boolean setsAreMatched = false;

        TransformationPointFit fit =
            calculateTranslationFromGridThenDownhillSimplex(
            scaledRotatedSet1, set1, set2,
            image1Width, image1Height, image2Width, image2Height,
            rotationInRadians, scale,
            transXStartStop, transDeltaX,
            transYStartStop, transDeltaY,
            tolTransX, tolTransY,
            setsAreMatched, setsFractionOfImage);

        if (fit == null) {
            return null;
        }

log.fine("     * ==> bestFit=" + fit.toString());

        int maxNMatchable = (set1.getN() < set2.getN()) ? set1.getN()
            : set2.getN();

        if (fit != null) {
            fit.setMaximumNumberMatchable(maxNMatchable);
        }

        return fit;
    }

    /**
     * the maximum number of iterations that the refinement of translation
     * will use in the downhill simplex.  The default value is 50.
     */
    private int dsNMaxIter = 50;
    protected void setDsNMaxIter(int n) {
        dsNMaxIter = n;
    }
    protected int getDsNMaxIter() {
        return dsNMaxIter;
    }
    float nEpsFactor = 2.0f;
    protected void setNEpsFactor(float f) {
        nEpsFactor = f;
    }
    protected float getNEpsFactor() {
        return nEpsFactor;
    }
    protected float cellFactor = 1.25f;
    protected float tolFactor = 0.67f;
    protected void setCellFactor(float f) {
        cellFactor = f;
    }
    protected void setTolFactor(float f) {
        tolFactor = f;
    }
    protected float getCellFactor() {
        return cellFactor;
    }
    protected float getTolFactor() {
        return tolFactor;
    }

    /**
     * find the best translations in x and y given the rotation and scale already
     * applied to set 1 to match set1 to set2 and return the fit for it.
     * If rotation and scale are already precisely determined, and if transXDelta
     * and transYDelta are values of 4 or smaller and if the member variable
     * dnsMaxIter is 50, the resulting solution should be accurate to within
     * a few pixels in x and y.
     *
     * @param scaledRotatedSet1
     * @param set1
     * @param set2
     * @param image1Width
     * @param image1Height
     * @param image2Width
     * @param image2Height
     * @param rotationInRadians
     * @param scale
     * @param transXStartStop
     * @param transXDelta the search cell size for a translation in x.  The
     * recommended value is 4.  The first searched cell is transXStart + transXDelta, etc.
     * @param transYStartStop
     * @param transYDelta the search cell size for a translation in x.  The
     * recommended value is 4.
     * @param tolTransX the tolerance in a difference in x when evaluating the
     * fit (which is optimal matching between the transformed set 1 and untransformed set2).
     * The recommended value is pointMatcher.getTolFactor() * transXDelta;
     * @param tolTransY the tolerance in a difference in x when evaluating the
     * fit (which is optimal matching between the transformed set 1 and untransformed set2).
     * The recommended value is pointMatcher.getTolFactor() * transYDelta;
     * @param setsAreMatched
     * @param setsFractionOfImage
     * @return
     */
    protected TransformationPointFit
        calculateTranslationFromGridThenDownhillSimplex(
        PairFloatArrayUnmodifiable scaledRotatedSet1, PairIntArray set1, PairIntArray set2,
        final int image1Width, final int image1Height,
        final int image2Width, final int image2Height,
        final float rotationInRadians, final float scale,
        RangeInt transXStartStop, final int transXDelta,
        RangeInt transYStartStop, final int transYDelta,
        final float tolTransX, final float tolTransY, boolean setsAreMatched,
        float setsFractionOfImage) {

        if (scaledRotatedSet1 == null) {
            throw new IllegalArgumentException("scaledRotatedSet1 cannot be null");
        }
        if (set1 == null) {
            throw new IllegalArgumentException("set1 cannot be null");
        }
        if (set2 == null) {
            throw new IllegalArgumentException("set2 cannot be null");
        }
        if (scaledRotatedSet1.getN() != set1.getN()) {
            throw new IllegalArgumentException(
                "scaledRotatedSet1 has to be the same length as set1");
        }
        if (tolTransX < 1) {
            throw new IllegalArgumentException("tolTransX should be > 0");
        }
        if (tolTransY < 1) {
            throw new IllegalArgumentException("tolTransY should be > 0");
        }
        if (transXDelta == 0) {
            throw new IllegalArgumentException("transXDelta cannot be 0");
        }
        if (transYDelta == 0) {
            throw new IllegalArgumentException("transYDelta cannot be 0");
        }

        int numberOfBestToReturn = 10;

        TransformationPointFit[] fits = evaluateTranslationsOverGrid(
            set1, set2,
            image1Width, image1Height, image2Width, image2Height,
            rotationInRadians, scale,
            transXStartStop, transXDelta,
            transYStartStop, transYDelta,
            tolTransX, tolTransY,
            setsAreMatched, setsFractionOfImage, numberOfBestToReturn);

        if (fits == null) {
            return null;
        }

log.fine("evaluateTranslationsOverGrid is complete");

        if ((transXDelta == 1) && (transYDelta == 1)) {
            // every possible translation combination has been visited, so no
            // need for the Nelder-Mead to search within completed search space
            return fits[0];
        }

        int maxNMatchable = (set1.getN() < set2.getN()) ? set1.getN()
            : set2.getN();

        //if the first 4 or so top fits all have nMatchedPoints=1 or 0, then
        //don't use the downhill simplex.
        boolean tooFewMatches = true;
        for (int i = 0; i < 4; ++i) {
            TransformationPointFit fit = fits[i];
            if (fit != null && fit.getNumberOfMatchedPoints() > 1) {
                tooFewMatches = false;
                break;
            }
        }
        if (tooFewMatches) {
            return fits[0];
        }

        // the bounds are to keep the solution within the current
        // best range

        //the top fits solution define the region to search within further.
        //{translationXMin, translationXMax, translationYMin, translationYMax}
        float[] transXYMinMaxes = getTranslationMinAndMaxes(fits);

        float x0 = transXYMinMaxes[0] - tolTransX;
        float y0 = transXYMinMaxes[2] - tolTransY;
        float boundsXCenter = (transXYMinMaxes[0] + transXYMinMaxes[1])/2.f;
        float boundsYCenter = (transXYMinMaxes[2] + transXYMinMaxes[3])/2.f;
        float boundsXHalf = Math.abs(x0 - boundsXCenter);
        float boundsYHalf = Math.abs(y0 - boundsYCenter);

log.fine("begin refineTranslationWithDownhillSimplex w/ maxIter=" + dsNMaxIter);

        TransformationPointFit fit = refineTranslationWithDownhillSimplex(
            scaledRotatedSet1, set2, fits,
            boundsXCenter, boundsYCenter, tolTransX, tolTransY,
            boundsXHalf, boundsYHalf,
            scale, rotationInRadians,
            setsAreMatched, dsNMaxIter);

        if (fit != null) {
            fit.setMaximumNumberMatchable(maxNMatchable);
        }

        return fit;
    }

    /**
     * given the scale, rotation and set 1's reference frame centroids,
     * calculate the translation between set1 and set2 assuming that not all
     * points will match.  transXTol and transYTol allow a tolerance when
     * matching the predicted position of a point in set2.
     * Note that the reference point for the rotation is the center of the
     * image 1 width and height.
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
     * @param rotationInRadians given in radians with value between 0 and 2*pi, exclusive
     * @param scale
     * @param image1Width
     * @param image1Height
     * @param image2Width
     * @param image2Height
     * @return
     */
    public TransformationParameters calculateTranslationForMatched(
        PairIntArray matched1, PairIntArray matched2,
        float rotationInRadians, float scale,
        int image1Width, int image1Height, int image2Width, int image2Height) {

        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");

            log.severe("scale cannot be smaller than 1");

            return null;
        }

        int image1CentroidX = image1Width >> 1;
        int image1CentroidY = image1Height >> 1;

        double scaleTimesCosine = scale * Math.cos(rotationInRadians);
        double scaleTimesSine = scale * Math.sin(rotationInRadians);

        double avgTransX = 0;
        double avgTransY = 0;

        for (int i = 0; i < matched1.getN(); i++) {

            int x = matched1.getX(i);
            int y = matched1.getY(i);

            double xr = image1CentroidX*scale + (
                ((x - image1CentroidX) * scaleTimesCosine) +
                ((y - image1CentroidY) * scaleTimesSine));

            double yr = image1CentroidY*scale + (
                (-(x - image1CentroidX) * scaleTimesSine) +
                ((y - image1CentroidY) * scaleTimesCosine));

            int x2 = matched2.getX(i);

            int y2 = matched2.getY(i);

            avgTransX += (int)Math.round(x2 - xr);

            avgTransY += (int)Math.round(y2 - yr);
        }

        avgTransX /= (float)matched1.getN();

        avgTransY /= (float)matched1.getN();

        TransformationParameters params =
            new TransformationParameters();
        params.setRotationInRadians(rotationInRadians);
        params.setScale(scale);
        params.setTranslationX((float)avgTransX);
        params.setTranslationY((float)avgTransY);

        return params;
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
     * @param scaledRotatedSet1 the model xy points scaled and rotated
     * @param transX the x translation to apply to the model points
     * @param transY the y translation to apply to the model points
     * @return
     */
    protected TransformationPointFit evaluateFitForUnMatchedGreedy(
        PairFloatArrayUnmodifiable scaledRotatedSet1,
        float transX, float transY,
        float tolTransX, float tolTransY,
        PairIntArray set2,
        final float scale, final float rotationRadians) {

        if (set2 == null) {
            throw new IllegalArgumentException(
            "set2 cannot be null");
        }
        if (scaledRotatedSet1 == null) {
            throw new IllegalArgumentException(
            "scaledRotatedSet1 cannot be null");
        }

        int nMaxMatchable = (scaledRotatedSet1.getN() < set2.getN()) ?
            scaledRotatedSet1.getN() : set2.getN();

        Set<Integer> chosen = new HashSet<Integer>();

        double[] diffs = new double[scaledRotatedSet1.getN()];
        int nMatched = 0;
        double avg = 0;

        for (int i = 0; i < scaledRotatedSet1.getN(); i++) {

            float transformedX = scaledRotatedSet1.getX(i) + transX;
            float transformedY = scaledRotatedSet1.getY(i) + transY;

            double minDiff = Double.MAX_VALUE;
            int min2Idx = -1;

            for (int j = 0; j < set2.getN(); j++) {

                if (chosen.contains(Integer.valueOf(j))) {
                    continue;
                }

                float dx = set2.getX(j) - transformedX;
                float dy = set2.getY(j) - transformedY;

                if ((Math.abs(dx) > tolTransX) || (Math.abs(dy) > tolTransY)) {

                    continue;
                }

                float diff = (float)Math.sqrt(dx*dx + dy*dy);

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
            avg, stDev, tolTransX, tolTransY);

        fit.setMaximumNumberMatchable(nMaxMatchable);

        return fit;
    }

    /**
     * sort the fits by descending number of matches.
     * @param fits
     * @param idxLo
     * @param idxHi, upper index, inclusive
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

        for (int i = idxLo; i < idxHi; ++i) {
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

    protected PairFloatArrayUnmodifiable scaleAndRotate(PairIntArray set1,
        float rotationInRadians, float scale, int centroidX1, int centroidY1) {

        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");

            log.severe("scale cannot be smaller than 1");
        }

        if (set1 == null) {
            throw new IllegalArgumentException("set1 cannot be null");
        }

        Transformer transformer = new Transformer();

        PairFloatArray scaledAndRotated = transformer.applyTransformation2(
            rotationInRadians, scale, 0, 0, centroidX1, centroidY1, set1);

        PairFloatArrayUnmodifiable unmod = new PairFloatArrayUnmodifiable(
            scaledAndRotated.getX(), scaledAndRotated.getY(), scaledAndRotated.getN());

        return unmod;
    }

    /**
     * Searches among the given translation range for the best fit to a
     * Euclidean transformation formed by scale, rotationRadians and the
     * best fitting translation X and Y given starter points fits.
     * Note the method is not precise so should be wrapped with a follow-up
     * method that further searches in the solution's local region.
     * @param scaledRotatedSet1
     * @param set2
     * @param fits
     * @param transX
     * @param transY
     * @param tolTransX
     * @param tolTransY
     * @param plusMinusTransX
     * @param plusMinusTransY
     * @param scale
     * @param rotationRadians
     * @param setsAreMatched
     * @param nMaxIter the maximum number of iterations to use in the algorithm.
     * The recommended value is 50.
     * @return
     */
    private TransformationPointFit refineTranslationWithDownhillSimplex(
        PairFloatArrayUnmodifiable scaledRotatedSet1, PairIntArray set2,
        TransformationPointFit[] fits,
        float transX, float transY, float tolTransX, float tolTransY,
        float plusMinusTransX, float plusMinusTransY,
        final float scale, final float rotationRadians,
        boolean setsAreMatched, int nMaxIter) {

        int nMaxMatchable = (scaledRotatedSet1.getN() < set2.getN()) ?
            scaledRotatedSet1.getN() : set2.getN();

        if (nMaxMatchable == 0) {
            return null;
        }

        //TODO: revise this.  consider a sqrt(nMaxMatchable)
        double eps = Math.log(nMaxMatchable)/Math.log(10);

        float alpha = 1;   // > 0
        float gamma = 2;   // > 1
        float beta = 0.5f;
        float tau = 0.5f;

        boolean go = true;

        float txMin = transX - plusMinusTransX;
        float txMax = transX + plusMinusTransX;
        float tyMin = transY - plusMinusTransY;
        float tyMax = transY + plusMinusTransY;

        int nIter = 0;

        int bestFitIdx = 0;
        int worstFitIdx = fits.length - 1;

        TransformationParameters[] lastParams = extractParameters(fits);

        while (go && (nIter < nMaxIter)) {

            if (fits.length == 0) {
                break;
            }

            sortByDescendingMatches(fits, 0, (fits.length - 1));

            for (int i = (fits.length - 1); i > -1; --i) {
                if (fits[i] != null) {
                    worstFitIdx = i;
                    break;
                }
            }
            if (fits[bestFitIdx] == null) {
                break;
            }

            /*if (debug) {
                if ((nIter == 0) || (nIter % 10 == 0)) {
                    plotTranslationSimplex(fits, txMin, txMax, tyMin, tyMax,
                        null);
                }
            }*/

            if (nIter > 0) {

                TransformationParameters[] currentParams = extractParameters(fits);

                boolean areTheSame = areEqual(lastParams, currentParams);

                if (areTheSame) {
                    break;
                }

                if (nIter > 50) {
                    double[] stdevs = getStandardDeviation(currentParams);
                    if (stdevs != null && (stdevs.length > 0)) {
                        boolean areNearlyZero = true;
                        for (double stdv : stdevs) {
                            if (Math.abs(stdv) > 0.01) {
                                areNearlyZero = false;
                            }
                        }
                        if (areNearlyZero) {
                            break;
                        }
                    }
                }

                lastParams = currentParams;
            }

            // determine center for all points excepting the worse fit
            float txSum = 0.0f;
            float tySum = 0.0f;
            int c = 0;
            for (int i = 0; i < (fits.length - 1); ++i) {
                if (fits[i] != null) {
                    txSum += fits[i].getTranslationX();
                    tySum += fits[i].getTranslationY();
                    c++;
                }
            }
            transX = txSum / (float)c;
            transY = tySum / (float)c;

            // "Reflection"
            float txReflect = transX + (alpha
                * (transX - fits[worstFitIdx].getTranslationX()));
            float tyReflect = transY + (alpha
                * (transY - fits[worstFitIdx].getTranslationY()));

            TransformationPointFit fitReflected = evaluateFit(
                scaledRotatedSet1, txReflect, tyReflect,
                tolTransX, tolTransY,
                set2, scale, rotationRadians, setsAreMatched);

            //TODO: consider putting back in a check for bounds

            int comp0 = compare(fits[bestFitIdx], fitReflected);
            int compLast = compare(fits[worstFitIdx], fitReflected);

            if ((comp0 < 1) && (compLast == 1)) {

                // replace last with f_refl
                fits[worstFitIdx] = fitReflected;

            } else if (comp0 == 1) {

                // reflected is better than best fit, so "expand"
                // "Expansion"
                float txExpansion = transX + (gamma * (txReflect - transX));
                float tyExpansion = transY + (gamma * (tyReflect - transY));

                TransformationPointFit fitExpansion =
                    evaluateFit(scaledRotatedSet1,
                        txExpansion, tyExpansion, tolTransX, tolTransY,
                        set2, scale, rotationRadians, setsAreMatched);

                int compR = compare(fitReflected, fitExpansion);

                if (compR == 1) {

                    // expansion fit is better than reflected fit
                    fits[worstFitIdx] = fitExpansion;

                } else {

                    fits[worstFitIdx] = fitReflected;
                }

            } else if (compLast < 1) {

                // reflected fit is worse than the worst (last) fit, so contract
                // "Contraction"
                float txContraction = transX + (beta
                    * (fits[worstFitIdx].getTranslationX() - transX));
                float tyContraction = transY + (beta
                    * (fits[worstFitIdx].getTranslationY() - transY));

                TransformationPointFit fitContraction
                    = evaluateFit(scaledRotatedSet1,
                        txContraction, tyContraction,
                        tolTransX, tolTransY,
                        set2, scale, rotationRadians, setsAreMatched);

                int compC = compare(fits[worstFitIdx], fitContraction);

                if (compC > -1) {

                    fits[worstFitIdx] = fitContraction;

                } else {

                    // "Reduction"
                    for (int i = 1; i < fits.length; ++i) {

                        if (fits[i] == null) {
                            /*TODO: consider setting this
                            fits[i] = new TransformationPointFit(
                                new TransformationParameters(),
                                0, Double.MAX_VALUE, Double.MAX_VALUE,
                                Double.MAX_VALUE
                            );
                            */
                            continue;
                        }

                        float txReduction
                            = (fits[bestFitIdx].getTranslationX()
                            + (tau
                            * (fits[i].getTranslationX()
                            - fits[bestFitIdx].getTranslationX())));

                        float tyReduction
                            = (fits[bestFitIdx].getTranslationY()
                            + (tau
                            * (fits[i].getTranslationY()
                            - fits[bestFitIdx].getTranslationY())));

                        //NOTE: there's a possibility of a null fit.
                        //  instead of re-writing the fits array, will
                        //  assign a fake infinitely bad fit which will
                        //  fall to the bottom of the list after the next
                        //  sort.
                        TransformationPointFit fit
                            = evaluateFit(scaledRotatedSet1,
                                txReduction, tyReduction,
                                tolTransX, tolTransY,
                                set2,
                                scale, rotationRadians, setsAreMatched);

                        fits[i] = fit;
                    }
                }
            }

            log.finest("best fit so far: nMatches="
                + fits[bestFitIdx].getNumberOfMatchedPoints()
                + " diff from model=" + fits[bestFitIdx].getMeanDistFromModel()
            );

            nIter++;

            if ((fits[bestFitIdx].getNumberOfMatchedPoints() == nMaxMatchable)
                && (fits[bestFitIdx].getMeanDistFromModel() < eps)) {
                go = false;
            } else if ((transX > txMax) || (transX < txMin)) {
                go = false;
            } else if ((transY > tyMax) || (transY < tyMin)) {
                go = false;
            }
        }

        log.fine("nIter=" + Integer.toString(nIter));

        /*if (debug) {
            plotTranslationSimplex(fits, txMin, txMax, tyMin, tyMax,
                Color.YELLOW);
            writeTranslationSimplexPlot();
        }*/

        // additional step that's helpful if not enough iterations are used,
        // is to test the summed transX, transY which represent the center
        // of the simplex against the best fit
        TransformationPointFit fitAvg = evaluateFit(scaledRotatedSet1,
            transX, transY, tolTransX, tolTransY, set2,
            scale, rotationRadians, setsAreMatched);

        int comp = compare(fits[bestFitIdx], fitAvg);
        if (comp == 1) {
            fits[bestFitIdx] = fitAvg;
        }

        return fits[bestFitIdx];
    }

    /**
     * Searches among the given translation range for the best fit to a
     * Euclidean transformation and the
     * best fitting translation X and Y given starter points fits.
     * Note the method is not precise so should be wrapped with a follow-up
     * method that further searches in the solution's local region.
     * @param scaledRotatedSet1
     * @param set2
     * @param fits
     * @param transX
     * @param transY
     * @param tolTransX
     * @param tolTransY
     * @param plusMinusTransX
     * @param plusMinusTransY
     * @param scale
     * @param rotationRadians
     * @param setsAreMatched
     * @param nMaxIter the maximum number of iterations to use in the algorithm.
     * The recommended value is 50.
     * @return
     */
    private TransformationPointFit refineTranslationWithDownhillSimplex(
        PairIntArray set1, PairIntArray set2,
        int image1Width, int image1Height,
        TransformationPointFit[] fits,
        float transX, float transY, float tolTransX, float tolTransY,
        float plusMinusTransX, float plusMinusTransY,
        final float scale, final float rotationRadians,
        boolean setsAreMatched, int nMaxIter) {

        int nMaxMatchable = (set1.getN() < set2.getN()) ?
            set1.getN() : set2.getN();

        if (nMaxMatchable == 0) {
            return null;
        }

        //TODO: revise this.  consider a sqrt(nMaxMatchable)
        double eps = Math.log(nMaxMatchable)/Math.log(10);

        float alpha = 1;   // > 0
        float gamma = 2;   // > 1
        float beta = 0.5f;
        float tau = 0.5f;

        boolean go = true;

        float txMin = transX - plusMinusTransX;
        float txMax = transX + plusMinusTransX;
        float tyMin = transY - plusMinusTransY;
        float tyMax = transY + plusMinusTransY;

        int nIter = 0;

        int bestFitIdx = 0;
        int worstFitIdx = fits.length - 1;

        TransformationParameters[] lastParams = extractParameters(fits);

        while (go && (nIter < nMaxIter)) {

            if (fits.length == 0) {
                break;
            }

            sortByDescendingMatches(fits, 0, (fits.length - 1));

            for (int i = (fits.length - 1); i > -1; --i) {
                if (fits[i] != null) {
                    worstFitIdx = i;
                    break;
                }
            }
            if (fits[bestFitIdx] == null) {
                break;
            }

            /*if (debug) {
                if ((nIter == 0) || (nIter % 10 == 0)) {
                    plotTranslationSimplex(fits, txMin, txMax, tyMin, tyMax,
                        null);
                }
            }*/

            if (nIter > 0) {

                TransformationParameters[] currentParams = extractParameters(fits);

                boolean areTheSame = areEqual(lastParams, currentParams);

                if (areTheSame) {
                    break;
                }

                if (nIter > 50) {
                    double[] stdevs = getStandardDeviation(currentParams);
                    if (stdevs != null && (stdevs.length > 0)) {
                        boolean areNearlyZero = true;
                        for (double stdv : stdevs) {
                            if (Math.abs(stdv) > 0.01) {
                                areNearlyZero = false;
                            }
                        }
                        if (areNearlyZero) {
                            break;
                        }
                    }
                }

                lastParams = currentParams;
            }

            // determine center for all points excepting the worse fit
            float txSum = 0.0f;
            float tySum = 0.0f;
            int c = 0;
            for (int i = 0; i < (fits.length - 1); ++i) {
                if (fits[i] != null) {
                    txSum += fits[i].getTranslationX();
                    tySum += fits[i].getTranslationY();
                    c++;
                }
            }
            transX = txSum / (float)c;
            transY = tySum / (float)c;

            // "Reflection"
            float txReflect = transX + (alpha
                * (transX - fits[worstFitIdx].getTranslationX()));
            float tyReflect = transY + (alpha
                * (transY - fits[worstFitIdx].getTranslationY()));

            TransformationPointFit fitReflected = transformAndEvaluateFit(
                set1, set2, image1Width, image1Height,
                scale, rotationRadians, txReflect, tyReflect,
                tolTransX, tolTransY, setsAreMatched);

            //TODO: consider putting back in a check for bounds

            int comp0 = compare(fits[bestFitIdx], fitReflected);
            int compLast = compare(fits[worstFitIdx], fitReflected);

            if ((comp0 < 1) && (compLast == 1)) {

                // replace last with f_refl
                fits[worstFitIdx] = fitReflected;

            } else if (comp0 == 1) {

                // reflected is better than best fit, so "expand"
                // "Expansion"
                float txExpansion = transX + (gamma * (txReflect - transX));
                float tyExpansion = transY + (gamma * (tyReflect - transY));

                TransformationPointFit fitExpansion =
                    transformAndEvaluateFit(set1, set2, image1Width, image1Height,
                        scale, rotationRadians,
                        txExpansion, tyExpansion, tolTransX, tolTransY,
                        setsAreMatched);

                int compR = compare(fitReflected, fitExpansion);

                if (compR == 1) {

                    // expansion fit is better than reflected fit
                    fits[worstFitIdx] = fitExpansion;

                } else {

                    fits[worstFitIdx] = fitReflected;
                }

            } else if (compLast < 1) {

                // reflected fit is worse than the worst (last) fit, so contract
                // "Contraction"
                float txContraction = transX + (beta
                    * (fits[worstFitIdx].getTranslationX() - transX));
                float tyContraction = transY + (beta
                    * (fits[worstFitIdx].getTranslationY() - transY));

                TransformationPointFit fitContraction
                    = transformAndEvaluateFit(set1, set2, image1Width, image1Height,
                        scale, rotationRadians,
                        txContraction, tyContraction, tolTransX, tolTransY,
                        setsAreMatched);

                int compC = compare(fits[worstFitIdx], fitContraction);

                if (compC > -1) {

                    fits[worstFitIdx] = fitContraction;

                } else {

                    // "Reduction"
                    for (int i = 1; i < fits.length; ++i) {

                        if (fits[i] == null) {
                            /*TODO: consider setting this
                            fits[i] = new TransformationPointFit(
                                new TransformationParameters(),
                                0, Double.MAX_VALUE, Double.MAX_VALUE,
                                Double.MAX_VALUE
                            );
                            */
                            continue;
                        }

                        float txReduction
                            = (fits[bestFitIdx].getTranslationX()
                            + (tau
                            * (fits[i].getTranslationX()
                            - fits[bestFitIdx].getTranslationX())));

                        float tyReduction
                            = (fits[bestFitIdx].getTranslationY()
                            + (tau
                            * (fits[i].getTranslationY()
                            - fits[bestFitIdx].getTranslationY())));

                        //NOTE: there's a possibility of a null fit.
                        //  instead of re-writing the fits array, will
                        //  assign a fake infinitely bad fit which will
                        //  fall to the bottom of the list after the next
                        //  sort.
                        TransformationPointFit fit
                            = transformAndEvaluateFit(set1, set2,
                                image1Width, image1Height, scale, rotationRadians,
                                txReduction, tyReduction,
                                tolTransX, tolTransY,
                                setsAreMatched);

                        fits[i] = fit;
                    }
                }
            }

            log.finest("best fit so far: nMatches="
                + fits[bestFitIdx].getNumberOfMatchedPoints()
                + " diff from model=" + fits[bestFitIdx].getMeanDistFromModel()
            );

            nIter++;

            if ((fits[bestFitIdx].getNumberOfMatchedPoints() == nMaxMatchable)
                && (fits[bestFitIdx].getMeanDistFromModel() < eps)) {
                go = false;
            } else if ((transX > txMax) || (transX < txMin)) {
                go = false;
            } else if ((transY > tyMax) || (transY < tyMin)) {
                go = false;
            }
        }

        log.fine("nIter=" + Integer.toString(nIter));

        /*if (debug) {
            plotTranslationSimplex(fits, txMin, txMax, tyMin, tyMax,
                Color.YELLOW);
            writeTranslationSimplexPlot();
        }*/

        // additional step that's helpful if not enough iterations are used,
        // is to test the summed transX, transY which represent the center
        // of the simplex against the best fit
        TransformationPointFit fitAvg = transformAndEvaluateFit(set1, set2,
            image1Width, image1Height,
            scale, rotationRadians, transX, transY, tolTransX, tolTransY,
            setsAreMatched);

        int comp = compare(fits[bestFitIdx], fitAvg);
        if (comp == 1) {
            fits[bestFitIdx] = fitAvg;
        }

        return fits[bestFitIdx];
    }

    protected TransformationPointFit evaluateFit(
        PairFloatArrayUnmodifiable scaledRotatedSet1, float transX,
        float transY, float tolTransX, float tolTransY,
        PairIntArray set2,
        final float scale, final float rotationRadians,
        boolean setsAreMatched) {

        if (set2 == null) {
            throw new IllegalArgumentException(
            "set2 cannot be null");
        }
        if (scaledRotatedSet1 == null) {
            throw new IllegalArgumentException(
            "scaledTransformedSet1 cannot be null");
        }

        if (setsAreMatched) {
            return evaluateFitForMatched(scaledRotatedSet1,
                transX, transY, set2, scale, rotationRadians);
        }

        return evaluateFitForUnMatched(scaledRotatedSet1,
            transX, transY, tolTransX, tolTransY,
            set2, scale, rotationRadians);
    }

    protected TransformationPointFit transformAndEvaluateFit(
        PairIntArray set1, PairIntArray set2,
        int image1Width, int image1Height,
        final float scale, final float rotationInRadians,
        float transX, float transY, float tolTransX, float tolTransY,

        boolean setsAreMatched) {

        if (set2 == null) {
            throw new IllegalArgumentException(
            "set2 cannot be null");
        }
        if (set1 == null) {
            throw new IllegalArgumentException(
            "set1 cannot be null");
        }

        PairFloatArrayUnmodifiable scaledRotatedSet1 = scaleAndRotate(set1,
            rotationInRadians, scale, image1Width >> 1, image1Height >> 1);

        return evaluateFit(scaledRotatedSet1, transX, transY, tolTransX,
            tolTransY, set2, scale, rotationInRadians, setsAreMatched);
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
        PairFloatArrayUnmodifiable scaledRotatedSet1,
        float transX, float transY, PairIntArray set2,
        final float scale, final float rotationRadians) {

        if (set2 == null || set2.getN() == 0) {
            throw new IllegalArgumentException(
            "set2 cannot be null or empty");
        }
        if (scaledRotatedSet1 == null || scaledRotatedSet1.getN() == 0) {
            throw new IllegalArgumentException(
            "scaledRotatedSet1 cannot be null or empty");
        }
        if (set2.getN() != scaledRotatedSet1.getN()) {
            throw new IllegalArgumentException(
            "for matched sets, set2 must be the same length as scaledRotated"
            + " X and Y");
        }

        int nMaxMatchable = (scaledRotatedSet1.getN() < set2.getN()) ?
            scaledRotatedSet1.getN() : set2.getN();

        if (nMaxMatchable == 0) {
            return null;
        }

        double sum = 0;
        double[] diff = new double[set2.getN()];

        for (int i = 0; i < set2.getN(); i++) {

            float transformedX = scaledRotatedSet1.getX(i) + transX;
            float transformedY = scaledRotatedSet1.getY(i) + transY;

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
            set2.getN(), avgDiff, stDev, Float.MAX_VALUE, Float.MAX_VALUE);

        fit.setMaximumNumberMatchable(nMaxMatchable);

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
     * @param scaledRotatedSet1 the model x,y points scaled and rotated
     * @param transX the x translation to apply to the model points
     * @param transY the y translation to apply to the model points
     * @return
     */
    private TransformationPointFit evaluateFitForUnMatched(
        PairFloatArrayUnmodifiable scaledRotatedSet1,
        float transX, float transY,
        float tolTransX, float tolTransY,
        PairIntArray set2,
        final float scale, final float rotationRadians) {

        return evaluateFitForUnMatchedOptimal(scaledRotatedSet1, transX, transY,
        //return evaluateFitForUnMatchedGreedy(scaledRotatedSet1, transX, transY,
            tolTransX, tolTransY,
            set2, scale, rotationRadians);
    }

    private TransformationParameters[] extractParameters(
        TransformationPointFit[] fits) {

        if (fits == null) {
            return new TransformationParameters[0];
        }

        TransformationParameters[] params = new TransformationParameters[fits.length];

        for (int i = 0; i < fits.length; ++i) {
            if (fits[i] != null) {
                params[i] = fits[i].getParameters();
            }
        }

        return params;
    }

    private double[] getStandardDeviation(TransformationParameters[] params) {

        if (params == null || params.length == 0) {
            return null;
        }

        double rSum = 0;
        double sSum = 0;
        double tXSum = 0;
        double tYSum = 0;

        int count = 0;
        for (int i = 0; i < params.length; ++i) {

            TransformationParameters p0 = params[i];

            if (p0 == null) {
                continue;
            }

            rSum += p0.getRotationInRadians();
            sSum += p0.getScale();
            tXSum += p0.getTranslationX();
            tYSum += p0.getTranslationY();
            count++;
        }
        rSum /= (double)count;
        sSum /= (double)count;
        tXSum /= (double)count;
        tYSum /= (double)count;

        double rStdv = 0;
        double sStdv= 0;
        double tXStdv = 0;
        double tYStdv = 0;
        for (int i = 0; i < params.length; ++i) {

            TransformationParameters p0 = params[i];

            if (p0 == null) {
                continue;
            }

            double r = p0.getRotationInRadians() - rSum;
            double s = p0.getScale() - sSum;
            double tX = p0.getTranslationX() - tXSum;
            double tY = p0.getTranslationY() - tYSum;

            rStdv += (r*r);
            sStdv += (s*s);
            tXStdv += (tX*tX);
            tYStdv += (tY*tY);
        }
        rStdv = (Math.sqrt(rStdv/(count - 1.0f)));
        sStdv = (Math.sqrt(sStdv/(count - 1.0f)));
        tXStdv = (Math.sqrt(tXStdv/(count - 1.0f)));
        tYStdv = (Math.sqrt(tYStdv/(count - 1.0f)));

        return new double[]{rStdv, sStdv, tXStdv, tYStdv};
    }

    /**
     * for unmatched sets of points, use a finer grid search to find the best
     * fit among the given parameters in similarToBestFit.
     *
     * @param similarToBestFit
     * @param set1
     * @param set2
     * @param image1Width
     * @param image1Height
     * @param image2Width
     * @param image2Height
     * @param setsFractionOfImage
     * @return
     */
    protected TransformationPointFit finerGridSearchToDistinguishBest(
        List<TransformationPointFit> similarToBestFit,
        PairIntArray set1, PairIntArray set2,
        int image1Width, int image1Height, int image2Width, int image2Height,
        float setsFractionOfImage) {

        if (similarToBestFit.isEmpty()) {
            return null;
        } else if (similarToBestFit.size() == 1) {
            return similarToBestFit.get(0);
        }

        int nIntervals = 3;

        TransformationPointFit bestFit = null;

        for (TransformationPointFit sFit : similarToBestFit) {

            TransformationPointFit fit = finerGridSearch(
                nIntervals,  sFit, set1, set2,
                image1Width, image1Height, image2Width, image2Height,
                setsFractionOfImage);

            if (fitIsBetter(bestFit, fit)) {
                bestFit = fit;
            }
        }

        return bestFit;
    }

    protected TransformationPointFit finerGridSearchToDistinguishBestForMatched(
        List<TransformationPointFit> similarToBestFit,
        PairIntArray set1, PairIntArray set2,
        int image1Width, int image1Height, int image2Width, int image2Height,
        float setsFractionOfImage) {

        if (similarToBestFit.isEmpty()) {
            return null;
        } else if (similarToBestFit.size() == 1) {
            return similarToBestFit.get(0);
        }

        int nIntervals = 3;

        TransformationPointFit bestFit = null;

        for (TransformationPointFit sFit : similarToBestFit) {

            TransformationPointFit fit = finerGridSearchForMatched(
                nIntervals,  sFit, set1, set2,
                image1Width, image1Height, image2Width, image2Height,
                setsFractionOfImage);

            if (fitIsBetter(bestFit, fit)) {
                bestFit = fit;
            }
        }

        return bestFit;
    }

    private TransformationPointFit finerGridSearch(
        int nIntervals, TransformationPointFit bestFit,
        PairIntArray set1, PairIntArray set2,
        int image1Width, int image1Height, int image2Width, int image2Height,
        float setsFractionOfImage) {

        double halfRange = 3 * bestFit.getMeanDistFromModel();

        if (Double.isInfinite(halfRange)) {
            return null;
        }

        RangeInt transXStartStop = new RangeInt((int)(bestFit.getTranslationX()
            - halfRange), (int)(bestFit.getTranslationX() + halfRange));

        RangeInt transYStartStop = new RangeInt((int)(bestFit.getTranslationY()
            - halfRange), (int)(bestFit.getTranslationY() + halfRange));

        transXStartStop.resetBoundsIfNeeded(-1*image2Width + 1, image2Width - 1);

        transYStartStop.resetBoundsIfNeeded(-1*image2Height + 1, image2Height - 1);

        int transXDelta = (transXStartStop.getStop() - transXStartStop.getStart())/nIntervals;
        int transYDelta = (transYStartStop.getStop() - transYStartStop.getStart())/nIntervals;
        if (transXDelta == 0) {
            transXDelta++;
        }
        if (transYDelta == 0) {
            transYDelta++;
        }

        float tolTransX2 = bestFit.getTranslationXTolerance();
        float tolTransY2 = bestFit.getTranslationYTolerance();

        boolean setsAreMatched = false;

        int image1CentroidX = image1Width >> 1;
        int image1CentroidY = image1Height >> 1;

        PairFloatArrayUnmodifiable scaledRotatedSet1 = scaleAndRotate(set1,
            bestFit.getRotationInRadians(), bestFit.getScale(),
            image1CentroidX, image1CentroidY);

        TransformationPointFit fit2 =
            calculateTranslationFromGridThenDownhillSimplex(
            scaledRotatedSet1, set1, set2,
            image1Width, image1Height, image2Width, image2Height,
            bestFit.getRotationInRadians(), bestFit.getScale(),
            transXStartStop, transXDelta,
            transYStartStop, transYDelta,
            tolTransX2, tolTransY2, setsAreMatched,
            setsFractionOfImage);

        return fit2;
    }

    private TransformationPointFit finerGridSearchForMatched(
        int nIntervals, TransformationPointFit bestFit,
        PairIntArray set1, PairIntArray set2,
        int image1Width, int image1Height, int image2Width, int image2Height,
        float setsFractionOfImage) {

        double halfRange = 3 * bestFit.getMeanDistFromModel();

        if (Double.isInfinite(halfRange)) {
            return null;
        }

        RangeInt transXStartStop = new RangeInt((int)(bestFit.getTranslationX()
            - halfRange), (int)(bestFit.getTranslationX() + halfRange));

        RangeInt transYStartStop = new RangeInt((int)(bestFit.getTranslationY()
            - halfRange), (int)(bestFit.getTranslationY() + halfRange));

        transXStartStop.resetBoundsIfNeeded(-1*image2Width + 1, image2Width - 1);

        transYStartStop.resetBoundsIfNeeded(-1*image2Height + 1, image2Height - 1);

        int transXDelta = (transXStartStop.getStop() - transXStartStop.getStart())/nIntervals;
        int transYDelta = (transYStartStop.getStop() - transYStartStop.getStart())/nIntervals;
        if (transXDelta == 0) {
            transXDelta++;
        }
        if (transYDelta == 0) {
            transYDelta++;
        }

        int image1CentroidX = image1Width >> 1;
        int image1CentroidY = image1Height >> 1;

        PairFloatArrayUnmodifiable scaledRotatedSet1 = scaleAndRotate(set1,
            bestFit.getRotationInRadians(), bestFit.getScale(),
            image1CentroidX, image1CentroidY);

        float unusedTolerance = Float.MIN_VALUE;
        boolean setsAreMatched = true;

        TransformationPointFit fit2 =
            calculateTranslationFromGridThenDownhillSimplex(
            scaledRotatedSet1, set1, set2,
            image1Width, image1Height, image2Width, image2Height,
            bestFit.getRotationInRadians(), bestFit.getScale(),
            transXStartStop, transXDelta,
            transYStartStop, transYDelta,
            unusedTolerance, unusedTolerance, setsAreMatched,
            setsFractionOfImage);

        return fit2;
    }

    /**
     * Re-evaluate the fit of the enclosed parameters, but use the new tolerance
     * for translation in x and y.
     *
     * @param fit
     * @param set1
     * @param set2
     * @param image1Width
     * @param image1Height
     * @param image2Width
     * @param image2Height
     * @param setsFractionOfImage
     * @return
     */
    private TransformationPointFit reevaluateForNewTolerance(
        TransformationPointFit fit, float translationXTolerance,
        float translationYTolerance, PairIntArray set1, PairIntArray set2,
        int image1Width, int image1Height) {

        TransformationPointFit fit2 = reevaluateForNewTolerance(
            fit, translationXTolerance, translationYTolerance,
            set1, set2, image1Width, image1Height, true);

        return fit2;
    }

    /**
     * Re-evaluate the fit of the enclosed parameters, but use the new tolerance
     * for translation in x and y.
     *
     * @param fit
     * @param set1
     * @param set2
     * @param image1Width
     * @param image1Height
     * @param image2Width
     * @param image2Height
     * @param setsFractionOfImage
     * @return
     */
    private TransformationPointFit reevaluateForNewToleranceOptimal(
        TransformationPointFit fit, float translationXTolerance,
        float translationYTolerance, PairIntArray set1, PairIntArray set2,
        int image1Width, int image1Height) {

        TransformationPointFit fit2 = reevaluateForNewTolerance(
            fit, translationXTolerance, translationYTolerance,
            set1, set2, image1Width, image1Height, false);

        return fit2;
    }

    /**
     * Re-evaluate the fit of the enclosed parameters, but use the new tolerance
     * for translation in x and y.
     *
     * @param fit
     * @param set1
     * @param set2
     * @param image1Width
     * @param image1Height
     * @param image2Width
     * @param image2Height
     * @param setsFractionOfImage
     * @return
     */
    private TransformationPointFit reevaluateForNewTolerance(
        TransformationPointFit fit, float translationXTolerance,
        float translationYTolerance, PairIntArray set1, PairIntArray set2,
        int image1Width, int image1Height, boolean useGreedyMatching) {

        if (fit == null) {
            return null;
        }

        float rotationInRadians = fit.getRotationInRadians();

        float scale = fit.getScale();

        int image1CentroidX = image1Width >> 1;
        int image1CentroidY = image1Height >> 1;

        Transformer transformer = new Transformer();

        PairFloatArray transformed1 = transformer.applyTransformation2(
            rotationInRadians, scale,
            fit.getTranslationX(), fit.getTranslationY(),
            image1CentroidX, image1CentroidY, set1);

        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians(rotationInRadians);
        params.setScale(scale);
        params.setTranslationX(fit.getTranslationX());
        params.setTranslationY(fit.getTranslationY());

        TransformationPointFit fit2;

        if (useGreedyMatching) {
            fit2 = evaluateFitForUnMatchedTransformedGreedy(
                params, transformed1, set2, translationXTolerance,
                translationYTolerance);
        } else {
            fit2 = evaluateFitForUnMatchedTransformedOptimal(
                params, transformed1, set2, translationXTolerance,
                translationYTolerance);
        }

        return fit2;
    }

    /**
     * find the minima and maxima of translationX and translationY for the given
     * fits.
     * @param fits
     * @return new float[]{minTX, maxTX, minTY, maxTY}
     */
    protected float[] getTranslationMinAndMaxes(TransformationPointFit[] fits) {

        if (fits == null || fits.length == 0) {
            return null;
        }
        float minTX = Float.MAX_VALUE;
        float maxTX = Float.MIN_VALUE;
        float minTY = Float.MAX_VALUE;
        float maxTY = Float.MIN_VALUE;

        for (TransformationPointFit fit : fits) {
            if (fit == null) {
                continue;
            }
            float tx = fit.getTranslationX();
            float ty = fit.getTranslationY();
            if (tx < minTX) {
                minTX = tx;
            }
            if (ty < minTY) {
                minTY = ty;
            }
            if (tx > maxTX) {
                maxTX = tx;
            }
            if (ty > maxTY) {
                maxTY = ty;
            }
        }

        return new float[]{minTX, maxTX, minTY, maxTY};
    }

    /**
     * re-evaluate bestFit or fit, whichever has the largest translation
     * tolerance, using the smaller tolerance.  Note that there are
     * exception rules, such as when both bestFit and fit have
     * same number of points, but fit has a mean dist from model less than one
     * and bestFit has a much larger mean distance from model.  In that case,
     * even if bestFit has a smaller translation tolerance, fit will not
     * be re-evaluated because such a mean distance from model means the
     * answer has converged.
     *
     *
     * @param bestFit
     * @param fit
     * @param set1
     * @param set2
     * @param image1Width
     * @param image1Height
     * @param reevalFits
     * @param fitIsBetter
     */
    protected void reevaluateFitsForCommonTolerance(
        TransformationPointFit bestFit, TransformationPointFit fit,
        PairIntArray set1, PairIntArray set2,
        int image1Width, int image1Height,
        final TransformationPointFit[] reevalFits, final boolean[] fitIsBetter) {

        if (bestFit == null) {
            reevalFits[0] = bestFit;
            reevalFits[1] = fit;
            fitIsBetter[0] = true;
            return;
        } else if (fit == null) {
            reevalFits[0] = bestFit;
            reevalFits[1] = fit;
            fitIsBetter[0] = false;
            return;
        }

        /*
        check for whether fit has converged already for equal number of points
        matched.
        */
        int bestNMatches = bestFit.getNumberOfMatchedPoints();
        int compNMatches = fit.getNumberOfMatchedPoints();
        int diffEps = (int)Math.round(2.*Math.ceil(Math.max(bestNMatches,
            compNMatches)/10.));
        if (diffEps == 0) {
            diffEps = 1;
        }

        double compAvg = fit.getMeanDistFromModel();
        double bestAvg = bestFit.getMeanDistFromModel();

        double compS = fit.getStDevFromMean();
        double bestS = bestFit.getStDevFromMean();

        double aDiv = bestAvg/compAvg;

        // TODO:  this skips the equal tolerance re-eval if solutions are already
        //    alot better and have significant number of points matched
        if ((compNMatches > 10) && (bestNMatches > 10)) {
            float nDiv = (float)bestNMatches/(float)compNMatches;
            double sDiv = bestS/compS;
            double nLimit0 = 2;
            double nLimit1 = 0.5;
            if (bestNMatches > 50) {
                nLimit0 = 1.5;
            }
            if (compNMatches > 50) {
                nLimit1 = 1./1.5;
            }
            if ((nDiv > nLimit0) && (aDiv < 0.5) && (sDiv < 0.5)) {
                fitIsBetter[0] = false;
                reevalFits[0] = bestFit;
                reevalFits[1] = fit;
                return;
            } else if ((nDiv < nLimit1) && (aDiv > 2) && (sDiv > 2)) {
                fitIsBetter[0] = true;
                reevalFits[0] = bestFit;
                reevalFits[1] = fit;
                return;
            }
        }

        /*
        when tolerances are both already very small, not redoing the fit,
        just comparing as is
        */
        int limit = 7;
        if ((bestFit.getTranslationXTolerance() < limit) &&
            (bestFit.getTranslationYTolerance() < limit) &&
            (fit.getTranslationXTolerance() < limit) &&
            (fit.getTranslationYTolerance() < limit)) {

            fitIsBetter[0] = fitIsBetter(bestFit, fit);
            reevalFits[0] = bestFit;
            reevalFits[1] = fit;
            return;
        }

        /*
        -1 : both are not null and bestFit tolerances are smaller
         0 : both are not null and tolerances are same.
         1 : both are not null and fit tolerances are smaller
         2 : both are not null and the x and y fits and smaller and larger in a mix
         3 : either bestFit or fit is null
        */
        int compTol = compareTolerance(bestFit, fit);

        TransformationPointFit bestFitT = bestFit;
        TransformationPointFit fitT = fit;

        float tolX, tolY;
        if (compTol != 3) {
            if (compTol == 1) {
                tolX = fit.getTranslationXTolerance();
                tolY = fit.getTranslationYTolerance();
            } else if (compTol == -1) {
                tolX = bestFit.getTranslationXTolerance();
                tolY = bestFit.getTranslationYTolerance();
            } else {
                tolX = bestFit.getTranslationXTolerance();
                if (fit.getTranslationXTolerance() < tolX) {
                    tolX = fit.getTranslationXTolerance();
                }
                tolY = bestFit.getTranslationYTolerance();
                if (fit.getTranslationYTolerance() < tolY) {
                    tolY = fit.getTranslationYTolerance();
                }
            }
            /*
            for small tolerances, both need to be re-evaluated with optimal
            matching, and then compared
            but the results should not be saved because that would result
            in an inconsistent mix of comparisons with greedy and optimal
            matching outside of this method
            */
            int limitT = 11;
            if ((tolX < limitT) && (tolY < limitT) && (compTol != 0)) {

                bestFitT = reevaluateForNewToleranceOptimal(bestFit,
                    tolX, tolY,
                    set1, set2,
                    image1Width, image1Height);

                fitT = reevaluateForNewToleranceOptimal(fit,
                    tolX, tolY,
                    set1, set2,
                    image1Width, image1Height);

                float nBDiv = (float)bestFit.getNumberOfMatchedPoints()/
                    (float)bestFitT.getNumberOfMatchedPoints();

                float nFDiv = (float)fit.getNumberOfMatchedPoints()/
                    (float)fitT.getNumberOfMatchedPoints();

                if (
                    (fitT == null) || (fitT.getNumberOfMatchedPoints() < 2) ||
                    (bestFitT == null) ||
                    (bestFitT.getNumberOfMatchedPoints() < 2)
                    //|| (nBDiv > 4) || (nFDiv > 4)
                ) {
                    fitIsBetter[0] = fitIsBetter(bestFit, fit);
                    reevalFits[0] = bestFit;
                    reevalFits[1] = fit;
                } else {
                    fitIsBetter[0] = fitIsBetter(bestFitT, fitT);
                    reevalFits[0] = bestFit;
                    reevalFits[1] = fit;
                }
                return;
            }

            if (compTol == 1) {

                bestFitT = reevaluateForNewTolerance(bestFit, tolX, tolY, set1,
                    set2, image1Width, image1Height);

                //TODO: may need to discard only if smaller tolerances < 10
                // do not use the lower tolerance.  it may be a false fit if null here
                if (bestFitT == null) {
                    bestFitT = bestFit;
                } else if (bestFitT.getNumberOfMatchedPoints() == 0) {
                    if (bestFit.getNumberOfMatchedPoints() > 10) {
                        bestFitT = bestFit;
                    }
                }
            } else if (compTol == -1) {

                fitT = reevaluateForNewTolerance(fit, tolX, tolY, set1, set2,
                    image1Width, image1Height);

                //TODO: may need to discard only if smaller tolerances < 10
                // do not use the lower tolerance if resulted in null fit
                if (fitT == null) {
                    fitT = fit;
                } else if (fitT.getNumberOfMatchedPoints() == 0) {
                    if (fit.getNumberOfMatchedPoints() > 10) {
                        fitT = fit;
                    }
                } else if (true) {
                    // need to try optimal matching for both.
                    // TODO: may be what the entire method should use for simplicity.
                    if (
                        (fit.getNumberOfMatchedPoints() > 10) &&
                        (((float)fit.getNumberOfMatchedPoints()/(float)fitT.getNumberOfMatchedPoints())
                        > 4)) {

                        bestFitT = reevaluateForNewToleranceOptimal(bestFit,
                            tolX, tolY,
                            set1, set2,
                            image1Width, image1Height);

                        fitT = reevaluateForNewToleranceOptimal(fit,
                            tolX, tolY,
                            set1, set2,
                            image1Width, image1Height);

                        float nBDiv = (float) bestFit.getNumberOfMatchedPoints()
                            / (float) bestFitT.getNumberOfMatchedPoints();

                        float nFDiv = (float) fit.getNumberOfMatchedPoints()
                            / (float) fitT.getNumberOfMatchedPoints();

                        if ((fitT == null) || (fitT.getNumberOfMatchedPoints() < 2)
                            || (bestFitT == null)
                            || (bestFitT.getNumberOfMatchedPoints() < 2) //|| (nBDiv > 4) || (nFDiv > 4)
                            ) {
                            fitIsBetter[0] = fitIsBetter(bestFit, fit);
                            reevalFits[0] = bestFit;
                            reevalFits[1] = fit;
                        } else {
                            fitIsBetter[0] = fitIsBetter(bestFitT, fitT);
                            reevalFits[0] = bestFit;
                            reevalFits[1] = fit;
                        }
                        return;
                    }
                }

            } else if (compTol == 2) {

                // TODO: may need to revise this
                // reduce both to smallest tolerances
                bestFitT = reevaluateForNewTolerance(bestFit,
                    tolX, tolY, set1, set2, image1Width, image1Height);

                fitT = reevaluateForNewTolerance(fit,
                    tolX, tolY, set1, set2, image1Width, image1Height);

                /*
                 // do not use the lower tolerance if resulted in null fits
                 if (fitT == null) {
                 fitT = fit;
                 bestFitT = bestFit;
                 } else if (fitT.getNumberOfMatchedPoints() == 0) {
                 if (fit.getNumberOfMatchedPoints() > 10) {
                 fitT = fit;
                 bestFitT = bestFit;
                 }
                 } else if (bestFitT == null) {
                 fitT = fit;
                 bestFitT = bestFit;
                 } else if (bestFitT.getNumberOfMatchedPoints() == 0) {
                 if (bestFitT.getNumberOfMatchedPoints() > 10) {
                 fitT = fit;
                 bestFitT = bestFit;
                 }
                 }*/
            }
        }

        fitIsBetter[0] = fitIsBetter(bestFitT, fitT);

        reevalFits[0] = bestFit;

        if (fitIsBetter[0] && (fitT != null)) {
           reevalFits[1] = fitT;
        }

if (bestFit != null && fit != null) {
if (compTol == 1) {
    log.fine("    rot re-evaluated bestFit at lower tolerance");
} else if (compTol == 2) {
    log.fine("    rot re-evaluated bestFit and fit at common tolerance");
}
}

    }

    private PolygonPlotterPNG plotter = null;
    //TODO: put this in an aspect
    private void plotTranslationSimplex(TransformationPointFit[] fits,
        float minX, float maxX, float minY, float maxY, java.awt.Color clr) {

        try {
            if (plotter == null) {
                plotter = new PolygonPlotterPNG(minX, maxX, minY, maxY,
                    "translation simplex", "transX", "transY");
            }
            int count = 0;
            for (TransformationPointFit fit : fits) {
                if (fit == null) {
                    continue;
                }
                count++;
            }
            double[] x = new double[count];
            double[] y = new double[x.length];
            count = 0;
            for (TransformationPointFit fit : fits) {
                if (fit == null) {
                    continue;
                }
                x[count] = fit.getTranslationX();
                y[count] = fit.getTranslationY();
                count++;
            }

            if (clr == null) {
                plotter.addPolygon(x, y);
            } else {
                plotter.addPolygon(x, y, clr);
            }

        } catch (IOException e) {
            log.severe(e.getMessage());
        }
    }

    private void writeTranslationSimplexPlot() {

        if (plotter == null) {
            return;
        }

        try {
            plotter.writeFile(MiscDebug.getCurrentTimeFormatted());
        } catch (IOException e) {
            log.severe(e.getMessage());
        }

    }

    /**
     * refine the Euclidean transformation.
     *
     * @param fit
     * @param set1
     * @param set2
     * @param image1Width
     * @param image1Height
     * @param image2Width
     * @param image2Height
     * @param setsFractionOfImage
     * @return
     */
    public TransformationPointFit refineTheTransformation(
        TransformationPointFit fit, PairIntArray set1, PairIntArray set2,
        int image1Width, int image1Height, int image2Width, int image2Height,
        float setsFractionOfImage) {

        if (fit == null) {
            return null;
        }
        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        int nMaxMatchable = Math.min(set1.getN(), set2.getN());

        /*
         empirically derived from testing.
         TODO: these may change as tests are made more complex by adding
         points that should not match, etc.
        */
        double densX = ((double) nMaxMatchable / (double) (setsFractionOfImage * (image2Width)));
        double densY = ((double) nMaxMatchable / (double) (setsFractionOfImage * (image2Height)));

        float rotRange = 6;
        float rotDelta = 1;
        float scaleRange = 0;
        float scaleDelta = 0;
        float transXRange = 50;
        int transXDelta = 1;
        float transYRange = 50;
        int transYDelta = 1;
        float tolTransX = 2; //TODO: this one not yet tested
        float tolTransY = 2; //TODO: this one not yet tested

        //TODO: determine w/ tests whether refinement here OR only
        // for the best fit should be done
        if ((densX > 0.033) || (densY > 0.033)) {

 log.info("trying a downhill simplex fit");
            
            // translation had to be increased greatly, so increase the rotation range
            rotRange = 10;
            transXRange = 55;
            transYRange = 55;

            if (true) {

                Transformer transformer = new Transformer();

                TransformationPointFit[] starterPoints = new TransformationPointFit[99];
                int count = 0;
                float rotInDegrees = fit.getParameters().getRotationInDegrees();
                float[] txS = new float[]{
                    fit.getTranslationX() - (transXRange/2.f),
                    fit.getTranslationX(), 
                    fit.getTranslationX() + (transXRange/2.f)};
                float[] tyS = new float[]{
                    fit.getTranslationY() - (transYRange/2.f),
                    fit.getTranslationY(), 
                    fit.getTranslationY() + (transYRange/2.f)};

                for (float rotT = (rotInDegrees - (rotRange/2.f)); 
                    rotT <= (rotInDegrees + (rotRange/2.f)); ++rotT) {
                    
                    if (rotT < 0) {
                        rotT = 360 + rotT;
                    } else if (rotT > 359) {
                        rotT = rotT - 360;
                    }
                    for (float tx : txS) {
                        for (float ty : tyS) {

                            TransformationParameters params = new TransformationParameters();
                            params.setRotationInDegrees(rotT);
                            params.setScale(fit.getScale());
                            params.setTranslationX(tx);
                            params.setTranslationY(ty);

                            PairFloatArray allPoints1Tr = transformer.applyTransformation(
                                params, image1Width >> 1, image1Height >> 1, set1);

                            starterPoints[count] = evaluateFitForUnMatchedTransformedGreedy(params,
                                allPoints1Tr, set2, tolTransX, tolTransY);

                            count++;
                        }
                    }
                }

                // TODO: for now, keeping all starter points, but may want to reduce
                // to top 10

                TransformationPointFit fit2 = refineTranslationWithDownhillSimplex(
                    set1, set2, image1Width, image1Height, starterPoints,
                    fit.getTranslationX(), fit.getTranslationY(),
                    tolTransX, tolTransY, transXRange + 5, transYRange + 5,
                    fit.getScale(), fit.getRotationInRadians(), false, 50);

                if (fit2 != null) {
                    log.info("downhill simplex fit=" + fit2.toString());
                }

                if (fitIsBetter(fit, fit2)) {
                    fit = fit2;
                    if (fit != null) {
                        fit.setMaximumNumberMatchable(nMaxMatchable);
                    }
                }
                
                boolean hasConverged = hasConverged(fit, nMaxMatchable);

                if (hasConverged) {
                    return fit;                    
                } else {
                    rotRange = 6;
                    transXRange = 14;
                    transYRange = 14;
                    tolTransX = cellFactor * transXRange;
                    tolTransY = cellFactor * transXRange;
                }
            }
        }

        float rotStart = (float)
            (fit.getParameters().getRotationInDegrees() - (rotRange/2.));
        float rotStop = (float)
            (fit.getParameters().getRotationInDegrees() + (rotRange/2.));

        if (rotStart < 0) {
            rotStart = 360 + rotStart;
        }
        if (rotStop > 359) {
            rotStop = rotStop - 360;
        }

        float scaleStart = (float)(fit.getScale() - (scaleRange/2.));
        float scaleStop = (float)(fit.getScale() + (scaleRange/2.));

        RangeInt transXStartStop = new RangeInt(
            (int)(fit.getTranslationX() - (transXRange/2.f)),
            (int)(fit.getTranslationX() + (transXRange/2.f)));
        RangeInt transYStartStop = new RangeInt(
            (int)(fit.getTranslationY() - (transYRange/2.f)),
            (int)(fit.getTranslationY() + (transYRange/2.f)));

 log.info("refineTheTransformation: rotStart=" + rotStart
            + " rotStop=" + rotStop + " rotDelta=" + rotDelta
            + " transXStart=" + transXStartStop.getStart() +
              " transXStop=" + transXStartStop.getStop()
            + " transYStart=" + transYStartStop.getStart() +
              " transYStop=" + transYStartStop.getStop()
            + " transDeltaX=" + transXDelta + " transDeltaY=" + transYDelta
            + " scaleStart=" + scaleStart + " scaleStop=" + scaleStop
        );

        TransformationPointFit fit2 = calculateTransformationWithGridSearch(
            set1, set2,
            image1Width, image1Height, image2Width, image2Height,
            rotStart, rotStop, rotDelta,
            scaleStart, scaleStop, scaleDelta,
            transXStartStop, transXDelta,
            transYStartStop, transYDelta,
            tolTransX, tolTransY,
            setsFractionOfImage);

        if (fit2 != null) {
            log.fine("refined fit=" + fit2.toString());
        }

        if (fitIsBetter(fit, fit2)) {
            fit = fit2;
            if (fit != null) {
                fit.setMaximumNumberMatchable(nMaxMatchable);
            }
        }

        return fit;
    }

    /**
     * refine the Euclidean transformation.
     *
     * @param fit
     * @param set1
     * @param set2
     * @param image1Width
     * @param image1Height
     * @param image2Width
     * @param image2Height
     * @param rotRange
     * @param rotDelta
     * @param scaleRange
     * @param scaleDelta
     * @param transXRange
     * @param transXDelta
     * @param transYRange
     * @param transYDelta
     * @param tolTransX
     * @param tolTransY
     * @param setsFractionOfImage
     * @return
     */
    public TransformationPointFit refineTheTransformation(
        TransformationPointFit fit, PairIntArray set1, PairIntArray set2,
        int image1Width, int image1Height, int image2Width, int image2Height,
        float rotRange, float rotDelta, float scaleRange, float scaleDelta,
        float transXRange, int transXDelta,
        float transYRange, int transYDelta, float tolTransX,
        float tolTransY, float setsFractionOfImage) {

        if (fit == null) {
            return null;
        }
        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        int nMaxMatchable = Math.min(set1.getN(), set2.getN());

        float rotStart = (float)
            (fit.getParameters().getRotationInDegrees() - (rotRange/2.));
        float rotStop = (float)
            (fit.getParameters().getRotationInDegrees() + (rotRange/2.));

        if (rotStart < 0) {
            rotStart = 360 + rotStart;
        }
        if (rotStop > 359) {
            rotStop = rotStop - 360;
        }

        float scaleStart = (float)(fit.getScale() - (scaleRange/2.));
        float scaleStop = (float)(fit.getScale() + (scaleRange/2.));

        RangeInt transXStartStop = new RangeInt(
            (int)(fit.getTranslationX() - (transXRange/2.f)),
            (int)(fit.getTranslationX() + (transXRange/2.f)));
        RangeInt transYStartStop = new RangeInt(
            (int)(fit.getTranslationY() - (transYRange/2.f)),
            (int)(fit.getTranslationY() + (transYRange/2.f)));

 log.info("refineTheTransformation: rotStart=" + rotStart
            + " rotStop=" + rotStop + " rotDelta=" + rotDelta
            + " transXStart=" + transXStartStop.getStart() +
              " transXStop=" + transXStartStop.getStop()
            + " transYStart=" + transYStartStop.getStart() +
              " transYStop=" + transYStartStop.getStop()
            + " transDeltaX=" + transXDelta + " transDeltaY=" + transYDelta
            + " scaleStart=" + scaleStart + " scaleStop=" + scaleStop
        );

        TransformationPointFit fit2 = calculateTransformationWithGridSearch(
            set1, set2,
            image1Width, image1Height, image2Width, image2Height,
            rotStart, rotStop, rotDelta,
            scaleStart, scaleStop, scaleDelta,
            transXStartStop, transXDelta,
            transYStartStop, transYDelta,
            tolTransX, tolTransY,
            setsFractionOfImage);

        if (fit2 != null) {
            log.fine("refined fit=" + fit2.toString());
        }

        if (fitIsBetter(fit, fit2)) {
            fit = fit2;
            if (fit != null) {
                fit.setMaximumNumberMatchable(nMaxMatchable);
            }
        }

        return fit;
    }
}
