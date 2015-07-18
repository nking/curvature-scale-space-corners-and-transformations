package algorithms.imageProcessing;

import algorithms.compGeometry.PointPartitioner;
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
import algorithms.util.ResourceFinder;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
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

    protected int largeSearch0Limit = 20;

    protected boolean debug = true;

    public void setLargeSearch0Limit(int limit) {
        largeSearch0Limit = limit;
    }

    /**
     * NOT READY FOR USE
     *
     * @param unmatchedLeftXY
     * @param unmatchedRightXY
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
     * @return best fitting transformation between unmatched points sets
     *
     * @return
     */
    public TransformationPointFit performMatching(
        PairIntArray unmatchedLeftXY, PairIntArray unmatchedRightXY,
        PairIntArray outputMatchedLeftXY, PairIntArray outputMatchedRightXY,
        boolean useLargestToleranceForOutput) {

        Transformer transformer = new Transformer();

        TransformationPointFit transFit =
            calculateEuclideanTransformation(unmatchedLeftXY, unmatchedRightXY);

        if (transFit == null) {
            return null;
        }

        int nMaxMatchable = (unmatchedLeftXY.getN() < unmatchedRightXY.getN()) ?
            unmatchedLeftXY.getN() : unmatchedRightXY.getN();

        transFit.setMaximumNumberMatchable(nMaxMatchable);

        if (nMaxMatchable == 0) {
            return transFit;
        }

        // -- transform filtered1 for matching and evaluation

        PairFloatArray transformedLeft =
            transformer.applyTransformation2(transFit.getParameters(),
            unmatchedLeftXY);

        float translationXTolerance, translationYTolerance;

        if (useLargestToleranceForOutput) {
            translationXTolerance = generalTolerance * (float)Math.sqrt(1./2);
            translationYTolerance = translationXTolerance;
        } else {
            translationXTolerance = 2 * transFit.getTranslationXTolerance();
            translationYTolerance = 2 * transFit.getTranslationYTolerance();
        }
        float tolerance = (float)Math.sqrt(
            translationXTolerance*translationXTolerance
            + translationYTolerance*translationYTolerance);

        //evaluate the fit and store a statistical var: nmatched/nmaxmatchable

        float[][] matchIndexesAndDiffs = calculateMatchUsingOptimal(
            transformedLeft, unmatchedRightXY,
            translationXTolerance, translationYTolerance);

        PairFloatArray part1MatchedTransformed = new PairFloatArray();
        PairIntArray part2Matched = new PairIntArray();

        matchPoints(transformedLeft, unmatchedRightXY, tolerance,
            matchIndexesAndDiffs, part1MatchedTransformed, part2Matched);

        PairIntArray part1Matched = new PairIntArray();
        part2Matched = new PairIntArray();

        matchPoints(unmatchedLeftXY, unmatchedRightXY, tolerance,
            matchIndexesAndDiffs, part1Matched, part2Matched);

        TransformationPointFit fit2 = evaluateFitForMatchedTransformed(
            transFit.getParameters(), part1MatchedTransformed, part2Matched);

        fit2.setTranslationXTolerance(translationXTolerance);
        fit2.setTranslationYTolerance(translationYTolerance);
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

        int nMaxMatchable = Math.min(unmatched1Transformed.getN(),
            unmatched2.getN());

        fit.setMaximumNumberMatchable(nMaxMatchable);

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

        int nMaxMatchable = Math.min(unmatched1Transformed.getN(),
            unmatched2.getN());

        fit.setMaximumNumberMatchable(nMaxMatchable);

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
     * calculate for Euclidean transformation from set1 to set2 for
     * the points being unmatched.
     *
     * @param set1
     * @param set2
     * @return
     */
    public TransformationPointFit calculateEuclideanTransformation(
        PairIntArray set1, PairIntArray set2) {

        TransformationPointFit fit = calculateRoughTransformation(
            set1, set2);

        if (fit == null) {
            return null;
        }

        int nMaxMatchable = (set1.getN() < set2.getN()) ? set1.getN()
            : set2.getN();

        boolean converged = hasConverged(fit, nMaxMatchable);

        if (converged) {
            return fit;
        }

        float rotHalfRange = 45;
        float rotDelta = 2;
        float transHalfRange = 500;
        float transDelta = 4;
        boolean useGreedyMatching = true;

        if (fit.getScale() >= 1.0) {

            TransformationPointFit fit2 = refineTheTransformation(
                fit.getParameters(), set1, set2, rotHalfRange, rotDelta,
                transHalfRange, transDelta, transHalfRange, transDelta,
                useGreedyMatching);

            if (fitIsBetter(fit, fit2)) {
                fit = fit2;
            }

        } else {

            // scale is < 1 so reverse order of sets
            TransformationPointFit revFit = refineTheTransformation(
                fit.getParameters(), set2, set1, rotHalfRange, rotDelta,
                transHalfRange, transDelta, transHalfRange, transDelta,
                useGreedyMatching);

            TransformationParameters params = revFit.getParameters();

            // reverse the parameters.

            MatchedPointsTransformationCalculator tc = new
                MatchedPointsTransformationCalculator();

            TransformationParameters revParams = tc.swapReferenceFrames(
                params);

            TransformationPointFit fit2 = new TransformationPointFit(revParams,
                revFit.getNumberOfMatchedPoints(),
                revFit.getMeanDistFromModel(),
                revFit.getStDevFromMean(),
                revFit.getTranslationXTolerance(),
                revFit.getTranslationYTolerance());

            if (fitIsBetter(fit, fit2)) {
                fit = fit2;
                if (fit != null) {
                   fit.setMaximumNumberMatchable(nMaxMatchable);
                }
            }
        }

        return fit;
    }

    /**
     * Calculate Euclidean transformation for unmatched points.  As with all
     * methods, the user should try to reduce the number of points to the
     * best determined and distributed throughout the area of intersection
     * while reducing the overall number of points.
     * 
     * This method is in progress, but ignoring scale, produces a
     * solution that is accurate within 20 degrees of rotation and
     * 50 pixels in translation of X and the same for translation of Y.
     * The method should be followed by refineTheTransformation.
     * NOTE that the transformation within the fit may return a scale that
     * is smaller than 1.
     *
     * @param set1
     * @param set2
     *
     * @return
     */
    public TransformationPointFit calculateRoughTransformation(
        PairIntArray set1, PairIntArray set2) {

        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        boolean useGreedySearch = true;

        int nMaxMatchable = (set1.getN() < set2.getN()) ? set1.getN()
            : set2.getN();

        TransformationPointFit bestFit = null;

        TransformationPointFit bestFitForScale = null;

        float scaleStart = 1;
        float scaleStop = 5;
        float scaleDelta = 1;

        for (float scale = scaleStart; scale <= scaleStop; scale += scaleDelta) {

            TransformationPointFit fit = preSearch(set1, set2, scale,
                useGreedySearch);

            if (fitIsBetter(bestFit, fit)) {
                bestFit = fit;
            }

            if (fitIsBetter(bestFitForScale, bestFit) && (bestFit != null)) {

                bestFitForScale = bestFit;

            } else {

                log.fine("previous scale solution was better, so end scale iter");

                // revert to previous scale
                bestFit = bestFitForScale;

                // scale was probably smaller so return best solution
                break;
            }
        }

        // check whether the scale for set1 is larger than that of set2
        float fracMatched = (bestFit == null) ? 0 :
            (float)bestFit.getNumberOfMatchedPoints()/(float)nMaxMatchable;

        if (fracMatched < 0.3) {
            TransformationPointFit bestRevFit = null;
            TransformationPointFit bestRevFitForScale = null;
            // reverse the order to solve for possible scale < 1.
            for (float scale = scaleStart; scale <= scaleStop; scale += scaleDelta) {

                TransformationPointFit revFit = preSearch(set2, set1, scale,
                    useGreedySearch);

                if (fitIsBetter(bestRevFit, revFit)) {
                    bestRevFit = revFit;
                }

                if (fitIsBetter(bestRevFitForScale, bestRevFit) && (bestRevFit != null)) {

                    bestRevFitForScale = bestRevFit;

                } else {

                    log.fine("previous scale for reverse solution was better, so end scale iter");

                    // revert to previous scale
                    bestRevFit = bestRevFitForScale;

                    // scale was probably smaller so return best solution
                    break;
                }
            }
            if (fitIsBetter(bestFit, bestRevFit)) {
                // reverse the parameters to write transformation w.r.t. set1
                // needs a reference point in both datasets for direct calculation
                TransformationParameters params = bestRevFit.getParameters();

                MatchedPointsTransformationCalculator tc = new
                    MatchedPointsTransformationCalculator();

                TransformationParameters revParams = tc.swapReferenceFrames(
                    params);

                bestFit = new TransformationPointFit(revParams,
                    bestRevFit.getNumberOfMatchedPoints(),
                    bestRevFit.getMeanDistFromModel(),
                    bestRevFit.getStDevFromMean(),
                    bestRevFit.getTranslationXTolerance(),
                    bestRevFit.getTranslationYTolerance());

                bestFit.setMaximumNumberMatchable(nMaxMatchable);
            }
        }

        return bestFit;
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

    protected TransformationPointFit[] evaluateTranslationsOverGrid(
        PairIntArray set1, PairIntArray set2,
        final float rotationInRadians, final float scale,
        RangeInt transXStartStop, int transXDelta,
        RangeInt transYStartStop, int transYDelta,
        float tolTransX, float tolTransY,
        final boolean setsAreMatched,
        final int numberOfBestToReturn) {

        /*   _____
            |     |
            |_____| largest negative or positive translationX of set1 is the width of set2
                  _____
                 |     |
                 |_____|
        */

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

                PairFloatArray allPoints1Tr = transformer.applyTransformation2(
                    params, set1);

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

    TransformationPointFit[] refineTransformationWithDownhillSimplex(
        TransformationPointFit[] fits,
        PairIntArray set1, PairIntArray set2,
        float transX, float transY, float tolTransX, float tolTransY,
        float plusMinusTransX, float plusMinusTransY,
        boolean searchScaleToo,
        boolean useGreedyMatching,
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

        // alternates changes in rotation, translation of X and translation of Y

        boolean go = true;

        float txMin = transX - plusMinusTransX;
        float txMax = transX + plusMinusTransX;
        float tyMin = transY - plusMinusTransY;
        float tyMax = transY + plusMinusTransY;

        float rotInRadians = (fits[0] != null) ? fits[0].getRotationInRadians()
            : 0;
        float scale = (fits[0] != null) ? fits[0].getScale() : 1;

        int nIter = 0;
        int lastAlt = 0;

        int bestFitIdx = 0;
        int worstFitIdx = fits.length - 1;

        TransformationParameters[] lastParams = extractParameters(fits);

        while (go && (nIter < nMaxIter)) {

            if (fits.length == 0) {
                break;
            }

            sortByDescendingMatches(fits, 0, (fits.length - 1));

            if (fits[bestFitIdx] == null) {
                break;
            }

            for (int i = (fits.length - 1); i > -1; --i) {
                if (fits[i] != null) {
                    worstFitIdx = i;
                    break;
                }
            }

            if (debug) {
                if ((nIter == 0) || (nIter % 5 == 0)) {
                    writeSimplex(fits, false);
                }
            }

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
            double rotSum = 0.0f;
            double scaleSum = 0.0f;
            int c = 0;
            for (int i = 0; i < (fits.length - 1); ++i) {
                if (fits[i] != null) {
                    txSum += fits[i].getTranslationX();
                    tySum += fits[i].getTranslationY();
                    rotSum += (fits[i].getRotationInRadians()/(2.*Math.PI));
                    scaleSum += (fits[i].getScale());
                    c++;
                }
            }
            transX = txSum / (float)c;
            transY = tySum / (float)c;
            rotInRadians = (float)(rotSum/(float)c)*((float)(2.*Math.PI));
            scale = (float)(scaleSum / (float)c);
            
            // "Reflection"

            float txReflect = 0;
            float tyReflect = 0;
            float rotReflect = 0;
            float scaleReflect = 1;

            TransformationParameters paramsReflect = new TransformationParameters();
            if (fits[worstFitIdx] != null) {
                paramsReflect.setScale(fits[worstFitIdx].getScale());
                paramsReflect.setRotationInRadians(fits[worstFitIdx].getRotationInRadians());
                paramsReflect.setTranslationX(fits[worstFitIdx].getTranslationX());
                paramsReflect.setTranslationY(fits[worstFitIdx].getTranslationY());
                paramsReflect.setOriginX(fits[worstFitIdx].getParameters().getOriginX());
                paramsReflect.setOriginY(fits[worstFitIdx].getParameters().getOriginY());
            }

            switch(lastAlt) {
                case 0: {
                    // change translation X
                    txReflect = transX + (alpha
                        * (transX - fits[worstFitIdx].getTranslationX()));
                    paramsReflect.setTranslationX(txReflect);
                    break;
                }
                case 1: {
                    // change translation Y
                    tyReflect = transY + (alpha
                        * (transY - fits[worstFitIdx].getTranslationY()));
                    paramsReflect.setTranslationY(tyReflect);
                    break;
                }
                case 2: {
                    if (searchScaleToo) {
                        scaleReflect  = scale + (alpha
                            * (scale - fits[worstFitIdx].getScale()));
                        paramsReflect.setScale(scaleReflect);
                        break;
                    }
                    // else fall through to rotation
                }
                default: {
                    // change rotation
                    rotReflect = rotInRadians + (alpha
                        * (rotInRadians - fits[worstFitIdx].getRotationInRadians()));
                    if (rotReflect < 0) {
                        rotReflect = rotReflect + (float)(2.*Math.PI);
                    } else if (rotReflect >= (2.*Math.PI)) {
                        rotReflect = rotReflect - (float)(2.*Math.PI);
                    }
                    paramsReflect.setRotationInRadians(rotReflect);
                    break;
                }
            }

            TransformationPointFit fitReflected = (fits[worstFitIdx] != null) ?
                evaluateForUnmatched(paramsReflect, set1, set2,
                tolTransX, tolTransY, useGreedyMatching) : null;

            int comp0 = compare(fits[bestFitIdx], fitReflected);
            int compLast = compare(fits[worstFitIdx], fitReflected);

            if ((comp0 < 1) && (compLast == 1)) {

                // replace last with f_refl
                fits[worstFitIdx] = fitReflected;

            } else if (comp0 == 1) {

                // reflected is better than best fit, so "expand"
                // "Expansion"

                TransformationParameters paramsExpansion = new TransformationParameters();
                if (fitReflected != null) {
                    paramsExpansion.setScale(fitReflected.getScale());
                    paramsExpansion.setRotationInRadians(fitReflected.getRotationInRadians());
                    paramsExpansion.setTranslationX(fitReflected.getTranslationX());
                    paramsExpansion.setTranslationY(fitReflected.getTranslationY());
                    paramsExpansion.setOriginX(fitReflected.getParameters().getOriginX());
                    paramsExpansion.setOriginY(fitReflected.getParameters().getOriginY());
                }

                switch (lastAlt) {
                    case 0: {
                        // change translation X
                        float txExpansion = transX + (gamma * (txReflect - transX));
                        paramsExpansion.setTranslationX(txExpansion);
                        break;
                    }
                    case 1: {
                        // change translation Y
                        float tyExpansion = transY + (gamma * (tyReflect - transY));
                        paramsExpansion.setTranslationY(tyExpansion);
                        break;
                    }
                    case 2: {
                        if (searchScaleToo) {
                            float scaleExpansion = scale + (gamma * (scaleReflect - scale));
                            paramsExpansion.setScale(scaleExpansion);
                            break;
                        }
                        // else fall through to rotation
                    }
                    default: {
                        // change rotation
                        float rotExpansion = rotInRadians + (gamma * (rotReflect - rotInRadians));
                        if (rotExpansion < 0) {
                            rotExpansion = rotExpansion + (float)(2.*Math.PI);
                        } else if (rotExpansion >= (2.*Math.PI)) {
                            rotExpansion = rotExpansion - (float)(2.*Math.PI);
                        }
                        paramsExpansion.setRotationInRadians(rotExpansion);
                        break;
                    }
                }

                TransformationPointFit fitExpansion = (fitReflected != null) ?
                    evaluateForUnmatched(paramsExpansion, set1, set2,
                    tolTransX, tolTransY, useGreedyMatching) : null;

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
                
                TransformationParameters paramsContraction = new TransformationParameters();

                paramsContraction.setScale(fits[worstFitIdx].getScale());
                paramsContraction.setRotationInRadians(fits[worstFitIdx].getRotationInRadians());
                paramsContraction.setTranslationX(fits[worstFitIdx].getTranslationX());
                paramsContraction.setTranslationY(fits[worstFitIdx].getTranslationY());
                paramsContraction.setOriginX(fits[worstFitIdx].getParameters().getOriginX());
                paramsContraction.setOriginY(fits[worstFitIdx].getParameters().getOriginY());

                switch (lastAlt) {
                    case 0: {
                        // change translation X
                        float txContraction = transX + (beta
                            * (fits[worstFitIdx].getTranslationX() - transX));
                        paramsContraction.setTranslationX(txContraction);
                        break;
                    }
                    case 1: {
                        // change translation Y
                        float tyContraction = transY + (beta
                            * (fits[worstFitIdx].getTranslationY() - transY));
                        paramsContraction.setTranslationY(tyContraction);
                        break;
                    }
                    case 2: {
                        if (searchScaleToo) {
                            float scaleContraction = scale + (beta
                                * (fits[worstFitIdx].getScale() - scale));
                            paramsContraction.setScale(scaleContraction);
                            break;
                        }
                        // else fall through to rotation
                    }
                    default: {
                        // change rotation
                        float rotContraction = rotInRadians + (beta
                            * (fits[worstFitIdx].getRotationInRadians() - rotInRadians));
                        if (rotContraction < 0) {
                            rotContraction = rotContraction + (float)(2.*Math.PI);
                        } else if (rotContraction >= (2.*Math.PI)) {
                            rotContraction = rotContraction - (float)(2.*Math.PI);
                        }
                        paramsContraction.setRotationInRadians(rotContraction);
                        break;
                    }
                }

                TransformationPointFit fitContraction = (fits[worstFitIdx] != null) ?
                    evaluateForUnmatched(paramsContraction, set1, set2,
                    tolTransX, tolTransY, useGreedyMatching) : null;

                int compC = compare(fits[worstFitIdx], fitContraction);

                if (compC > -1) {

                    fits[worstFitIdx] = fitContraction;

                } else {

                    // "Reduction"
                    for (int i = 1; i < fits.length; ++i) {

                        if (fits[i] == null || fits[bestFitIdx] == null) {
                            /*TODO: consider setting this
                            fits[i] = new TransformationPointFit(
                                new TransformationParameters(),
                                0, Double.MAX_VALUE, Double.MAX_VALUE,
                                Double.MAX_VALUE
                            );
                            */
                            continue;
                        }

                        TransformationParameters paramsI = new TransformationParameters();
                        paramsI.setScale(fits[i].getScale());
                        paramsI.setRotationInRadians(fits[i].getRotationInRadians());
                        paramsI.setTranslationX(fits[i].getTranslationX());
                        paramsI.setTranslationY(fits[i].getTranslationY());
                        paramsI.setOriginX(fits[i].getParameters().getOriginX());
                        paramsI.setOriginY(fits[i].getParameters().getOriginY());

                        switch (lastAlt) {
                            case 0: {
                                // change translation X
                                float txReduction
                                    = (fits[bestFitIdx].getTranslationX()
                                    + (tau * (fits[i].getTranslationX()
                                    - fits[bestFitIdx].getTranslationX())));
                                paramsI.setTranslationX(txReduction);
                                break;
                            }
                            case 1: {
                                // change translation Y
                                float tyReduction
                                    = (fits[bestFitIdx].getTranslationY()
                                    + (tau * (fits[i].getTranslationY()
                                    - fits[bestFitIdx].getTranslationY())));
                                paramsI.setTranslationY(tyReduction);
                                break;
                            }
                            case 2: {
                                if (searchScaleToo) {
                                    float scaleReduction
                                        = (fits[bestFitIdx].getScale()
                                        + (tau * (fits[i].getScale()
                                        - fits[bestFitIdx].getScale())));
                                    paramsI.setScale(scaleReduction);
                                    break;
                                }
                                // else fall through to rotation
                            }
                            default: {
                                // change rotation
                                float rotReduction
                                    = (fits[bestFitIdx].getRotationInRadians()
                                    + (tau * (fits[i].getRotationInRadians()
                                    - fits[bestFitIdx].getRotationInRadians())));
                                if (rotReduction < 0) {
                                    rotReduction = rotReduction + (float)(2.*Math.PI);
                                } else if (rotReduction >= (2.*Math.PI)) {
                                    rotReduction = rotReduction - (float)(2.*Math.PI);
                                }
                                paramsI.setRotationInRadians(rotReduction);
                                break;
                            }
                        }
                        fits[i] = evaluateForUnmatched(paramsI, set1, set2,
                            tolTransX, tolTransY, useGreedyMatching);
                    }
                }
            }

            log.finest("best fit so far: nMatches="
                + fits[bestFitIdx].getNumberOfMatchedPoints()
                + " diff from model=" + fits[bestFitIdx].getMeanDistFromModel()
            );

            nIter++;
            lastAlt++;
            if (searchScaleToo && (lastAlt > 3)) {
                lastAlt = 0;
            } else if (!searchScaleToo && (lastAlt > 2)) {
                lastAlt = 0;
            }

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

        // additional step that's helpful if not enough iterations are used,
        // is to test the summed transX, transY which represent the center
        // of the simplex against the best fit
        if (fits[bestFitIdx] != null) {

            TransformationParameters p = new TransformationParameters();
            p.setScale(scale);
            p.setRotationInRadians(rotInRadians);
            p.setTranslationX(transX);
            p.setTranslationY(transY);
            p.setOriginX(fits[bestFitIdx].getParameters().getOriginX());
            p.setOriginY(fits[bestFitIdx].getParameters().getOriginY());

            TransformationPointFit fitAvg = evaluateForUnmatched(p, set1, set2,
                tolTransX, tolTransY, useGreedyMatching);

            int comp = compare(fits[bestFitIdx], fitAvg);
            if (comp == 1) {
                fits[bestFitIdx] = fitAvg;
            }
        }

        if (debug) {
            // the resulting file can be viewed with resources/plot_3d_simplex_n10.html
            writeSimplex(fits, true);
        }

        return fits;
    }

    protected TransformationPointFit transformAndEvaluateFit(
        PairIntArray set1, PairIntArray set2,
        TransformationParameters params, float tolTransX, float tolTransY,
        boolean useGreedyMatching, boolean setsAreMatched) {

        if (set2 == null) {
            throw new IllegalArgumentException(
            "set2 cannot be null");
        }
        if (set1 == null) {
            throw new IllegalArgumentException(
            "set1 cannot be null");
        }

        Transformer transformer = new Transformer();

        PairFloatArray transformedSet1 = transformer.applyTransformation2(
            params, set1);

        if (setsAreMatched) {

            return evaluateFitForMatchedTransformed(params, transformedSet1,
                set2);

        } else {

            return evaluateFitForUnMatchedTransformedGreedy(
                params, transformedSet1, set2, tolTransX, tolTransY);
        }
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
     * Re-evaluate the fit of the enclosed parameters, but use the new tolerance
     * for translation in x and y.
     *
     * @param fit
     * @param set1
     * @param set2
     * @return
     */
    private TransformationPointFit reevaluateForNewTolerance(
        TransformationPointFit fit, float translationXTolerance,
        float translationYTolerance, PairIntArray set1, PairIntArray set2) {

        TransformationPointFit fit2 = reevaluateForNewTolerance(
            fit, translationXTolerance, translationYTolerance,
            set1, set2, true);

        return fit2;
    }

    /**
     * Re-evaluate the fit of the enclosed parameters, but use the new tolerance
     * for translation in x and y.
     *
     * @param fit
     * @param set1
     * @param set2
     * @return
     */
    private TransformationPointFit reevaluateForNewToleranceOptimal(
        TransformationPointFit fit, float translationXTolerance,
        float translationYTolerance, PairIntArray set1, PairIntArray set2) {

        TransformationPointFit fit2 = reevaluateForNewTolerance(
            fit, translationXTolerance, translationYTolerance,
            set1, set2, false);

        return fit2;
    }

    /**
     * Re-evaluate the fit of the enclosed parameters, but use the new tolerance
     * for translation in x and y.
     *
     * @param fit
     * @param set1
     * @param set2
     * @return
     */
    private TransformationPointFit reevaluateForNewTolerance(
        TransformationPointFit fit, float translationXTolerance,
        float translationYTolerance, PairIntArray set1, PairIntArray set2,
        boolean useGreedyMatching) {

        if (fit == null) {
            return null;
        }

        float rotationInRadians = fit.getRotationInRadians();

        float scale = fit.getScale();

        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians(rotationInRadians);
        params.setScale(scale);
        params.setTranslationX(fit.getTranslationX());
        params.setTranslationY(fit.getTranslationY());

        TransformationPointFit fit2 = evaluateForUnmatched(params,
            set1, set2, translationXTolerance, translationYTolerance,
            useGreedyMatching);

        return fit2;
    }

    private TransformationPointFit evaluateForUnmatched(
        final TransformationParameters params,
        PairIntArray set1, PairIntArray set2,
        float translationXTolerance, float translationYTolerance,
        boolean useGreedyMatching) {

        if (params == null) {
            return null;
        }

        boolean setsAreMatched = false;

        return transformAndEvaluateFit(set1, set2, params,
            translationXTolerance, translationYTolerance,
            useGreedyMatching, setsAreMatched);
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
     * @param reevalFits
     * @param fitIsBetter
     */
    protected void reevaluateFitsForCommonTolerance(
        TransformationPointFit bestFit, TransformationPointFit fit,
        PairIntArray set1, PairIntArray set2,
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

        if ((bestFit.getTranslationXTolerance() == fit.getTranslationXTolerance()) &&
            (bestFit.getTranslationYTolerance() == fit.getTranslationYTolerance())) {
            fitIsBetter[0] = fitIsBetter(bestFit, fit);
            reevalFits[0] = bestFit;
            reevalFits[1] = fit;
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
                    set1, set2);

                fitT = reevaluateForNewToleranceOptimal(fit,
                    tolX, tolY,
                    set1, set2);

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
                    set2);

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

                fitT = reevaluateForNewTolerance(fit, tolX, tolY, set1, set2);

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
                            set1, set2);

                        fitT = reevaluateForNewToleranceOptimal(fit,
                            tolX, tolY,
                            set1, set2);

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
                    tolX, tolY, set1, set2);

                fitT = reevaluateForNewTolerance(fit,
                    tolX, tolY, set1, set2);
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

    private FileWriter simplexCSVFileW = null;
    private BufferedWriter simplexCSVFileWriter = null;
    //TODO: put this in an aspect
    private void writeSimplex(TransformationPointFit[] fits, boolean closeFile) {

        //need a convention for writing for nulls.  will use max of rot, transX,
        //and transY when a fit is null
        float maxRot = Float.MIN_VALUE;
        float maxTransX = Float.MIN_VALUE;
        float maxTransY = Float.MIN_VALUE;
        for (int i = 0; i < fits.length; ++i) {
            TransformationPointFit fit = fits[i];
            if (fit != null) {
                if (fit.getParameters().getRotationInDegrees() > maxRot) {
                    maxRot = fit.getParameters().getRotationInDegrees();
                }
                if (fit.getParameters().getTranslationX() > maxTransX) {
                    maxTransX = fit.getParameters().getTranslationX();
                }
                if (fit.getParameters().getTranslationY() > maxTransY) {
                    maxTransY = fit.getParameters().getTranslationY();
                }
            }
        }

        try {
            if (simplexCSVFileWriter == null) {
                String testDirPath = ResourceFinder.findOutputTestDirectory();
                String fileName = "simplex_" + MiscDebug.getCurrentTimeFormatted() + ".csv";
                String filePath = testDirPath
                    + System.getProperty("file.separator") + fileName;
                File file = new File(filePath);
                simplexCSVFileW = new FileWriter(file);
                simplexCSVFileWriter = new BufferedWriter(simplexCSVFileW);
            }

            for (TransformationPointFit fit : fits) {

                StringBuilder sb = new StringBuilder();

                float rot, tx, ty;

                if (fit == null) {
                    rot = maxRot;
                    tx = maxTransX;
                    ty = maxTransY;
                } else {
                    rot = fit.getParameters().getRotationInDegrees();
                    tx = fit.getTranslationX();
                    ty = fit.getTranslationY();
                }
                sb.append(Float.toString(rot)).append(",")
                    .append(Float.toString(tx)).append(",")
                    .append(Float.toString(ty)).append("\n");

                String str = sb.toString();

                simplexCSVFileWriter.write(str, 0, str.length());
            }
            simplexCSVFileWriter.flush();

            if (closeFile) {
                if (simplexCSVFileW != null) {
                    simplexCSVFileW.close();
                    simplexCSVFileW = null;
                }
                if (simplexCSVFileWriter != null) {
                    simplexCSVFileWriter.close();
                    simplexCSVFileWriter = null;
                }
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
     * method to narrow down all parameter space to 10 fits in which at least
     * one contains a fit within +- 10 degrees of rotation and +- 20 pixels
     * of translation in x and y.  Note, still testing so these numbers
     * may change.
     * The 10 fits should be followed by a downhill simplex to further
     * constrain the transformation, or a detailed grid search of each of
     * the 10 fits with something like a delta rotation of 2 and delta
     * translation of 4 and ranges based on the error stated in the first
     * sentence.
     *
     * @param set1
     * @param set2
     * @param scale
     * @param useGreedyMatching if true, uses an N^2 algorithm to find the
     * best match made for each point in order of points in set1, else uses
     * an ~N^3 optimal bipartite matching algorithm.
     * @return
     */
    protected TransformationPointFit[] preSearch0(
        PairIntArray set1, PairIntArray set2, float scale,
        boolean useGreedyMatching) {

        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        int nMaxMatchable = Math.min(set1.getN(), set2.getN());
        
        TransformationPointFit[] fits = null;

        // if number of points is very high, need to use alternate method
        if (nMaxMatchable > largeSearch0Limit) {
            fits = preSearch0Large(set1, set2, scale, useGreedyMatching);
        } else {
            fits = preSearch0Small(set1, set2, scale, useGreedyMatching);
        }
        
        return fits;
    }
    
    /**
     * method to narrow down all parameter space to 10 fits in which at least
     * one contains a fit within +- 10 degrees of rotation and +- 20 pixels
     * of translation in x and y.  Note, still testing so these numbers
     * may change.
     * The 10 fits should be followed by a downhill simplex to further
     * constrain the transformation, or a detailed grid search of each of
     * the 10 fits with something like a delta rotation of 2 and delta
     * translation of 4 and ranges based on the error stated in the first
     * sentence.
     *
     * @param set1
     * @param set2
     * @param scale
     * @param useGreedyMatching if true, uses an N^2 algorithm to find the
     * best match made for each point in order of points in set1, else uses
     * an ~N^3 optimal bipartite matching algorithm.
     * @return
     */
    protected TransformationPointFit[] preSearch0Small(
        PairIntArray set1, PairIntArray set2, float scale,
        boolean useGreedyMatching) {

        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        int nMaxMatchable = Math.min(set1.getN(), set2.getN());

        /*
        The best solutions use a fine grid search of rotDelta=2 and transDelta=4
        w/ translation tolerance = transDelta.

        Such a solution takes a long time so a coarser grid is attempted here.

        Trying a fine grid in rotation and a very coarse grid in translation
        to get rotation accurate, then a finer search on translation near
        that found rotation.
        */

        /*
        editing this while testing:

        rotDelta=1 and transDelta of 100 has a smallest delta rot within the
            10 best fits which is accurate within += 20 degrees?

        rotDelta=2 and transDelta of 100 has a smallest delta rot
            within the 10 best fits which is accurate within += 25 degrees... 14, 5
            and within the best fit alone, is accurate within

        */
        float transDelta = 100;
        float rotDelta = 2;

        float tolTransX = 0.5f*transDelta + (float)(Math.sin(rotDelta*Math.PI/180.)*transDelta);
        float tolTransY = tolTransX;

        Transformer transformer = new Transformer();

        float[] rotation = MiscMath.writeDegreeIntervals(0, 359, rotDelta);

        int nKeep = 10;

        long t0 = System.currentTimeMillis();
        log.info("starterPoints.length=" + nKeep);

        int minX2 = MiscMath.findMin(set2.getX(), set2.getN());
        int maxX2 = MiscMath.findMax(set2.getX(), set2.getN());
        int minY2 = MiscMath.findMin(set2.getY(), set2.getN());
        int maxY2 = MiscMath.findMax(set2.getY(), set2.getN());

        FixedSizeSortedVector<TransformationPointFit> starterPoints
            = new FixedSizeSortedVector<>(nKeep, TransformationPointFit.class);

        float[] xsr = new float[set1.getN()];
        float[] ysr = new float[set1.getN()];

        for (float rotDeg : rotation) {

            TransformationParameters rotScaleParams =
                new TransformationParameters();
            rotScaleParams.setScale(scale);
            rotScaleParams.setRotationInDegrees(rotDeg);
            rotScaleParams.setTranslationX(0);
            rotScaleParams.setTranslationY(0);
            rotScaleParams.setOriginX(0);
            rotScaleParams.setOriginY(0);

            transformer.applyTransformation(rotScaleParams, set1, xsr, ysr);

            final float minX1 = MiscMath.findMin(xsr);
            final float maxX1 = MiscMath.findMax(xsr);
            final float minY1 = MiscMath.findMin(ysr);
            final float maxY1 = MiscMath.findMax(ysr);

            float transXStart = minX1 - maxX2;
            float transXStop = maxX1 - minX2;
            float transYStart = minY1 - maxY2;
            float transYStop = maxY1 - minY2;

            float transXRange = transXStop - transXStart;
            float transYRange = transYStop - transYStart;

            int nX = (int) Math.ceil(transXRange / transDelta);
            int nY = (int) Math.ceil(transYRange / transDelta);

            float[] txS = new float[nX];
            txS[0] = transXStart;
            for (int i = 1; i < nX; ++i) {
                txS[i] = txS[i - 1] + transDelta;
            }
            float[] tyS = new float[nY];
            tyS[0] = transYStart;
            for (int i = 1; i < nY; ++i) {
                tyS[i] = tyS[i - 1] + transDelta;
            }

            for (float tx : txS) {
                for (float ty : tyS) {

                    PairFloatArray allPoints1Tr = new PairFloatArray(set1.getN());

                    for (int i = 0; i < xsr.length; ++i) {
                        allPoints1Tr.add(xsr[i] + tx, ysr[i] + ty);
                    }

                    TransformationParameters params = new TransformationParameters();
                    params.setScale(scale);
                    params.setRotationInDegrees(rotDeg);
                    params.setTranslationX(tx);
                    params.setTranslationY(ty);
                    params.setOriginX(0);
                    params.setOriginY(0);

                    TransformationPointFit fit2;
                    if (useGreedyMatching) {
                        fit2 = evaluateFitForUnMatchedTransformedGreedy(params,
                            allPoints1Tr, set2, tolTransX, tolTransY);
                    } else {
                        fit2 = evaluateFitForUnMatchedTransformedOptimal(params,
                            allPoints1Tr, set2, tolTransX, tolTransY);
                    }

                    starterPoints.add(fit2);
                }
            }
        }

        if (starterPoints.getArray()[0] == null) {
            return null;
        }

        long tm = System.currentTimeMillis() - t0;
        double ts = tm * 1e-3;
        log.info("starterPoints finished for nMaxMatchable=" + nMaxMatchable +
            " seconds=" + ts);

        if (starterPoints.getNumberOfItems() > 0) {
            log.info("best of starter points=" + starterPoints.getArray()[0].toString());
        }

        return starterPoints.getArray();
    }

    /**
     * method to narrow down all parameter space using the given starterPoints
     * to within +- 3 degrees and +- 50 in translation.
     * The starterPoints are expected to have been the result of preSearch0.
     * <pre>
     * Runtime complexity:
     * </pre>
     * 
     * @param fits
     * @param set1
     * @param set2
     * @param useGreedyMatching
     * @return
     */
    protected TransformationPointFit preSearch1Alt2Small(TransformationPointFit[] fits,
        PairIntArray set1, PairIntArray set2, boolean useGreedyMatching) {

        if (fits == null) {
            throw new IllegalArgumentException("fits cannot be null");
        }
        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        float scaleHalfRange = 0;
        float scaleDelta = 1;
        float rotHalfRangeInDegrees = 10;
        float rotDeltaInDegrees = 2;
        float transXHalfRange = 200;
        float transXDelta = 15;
        float transYHalfRange = transXHalfRange; 
        float transYDelta = transXDelta;

        Set<Integer> searchedRot = new HashSet<Integer>();

        TransformationPointFit bestFit = null;
        
        for (TransformationPointFit fit : fits) {

             Integer rotD = Integer.valueOf((int)fit.getParameters().getRotationInDegrees());

             if (searchedRot.contains(rotD)) {
                 continue;
             }
  
             // grid search around trans and rot

             TransformationPointFit[] fits2 = boundedGridSearch(set1, set2, 
                 fit.getParameters(),
                 scaleHalfRange, scaleDelta,
                 rotHalfRangeInDegrees, rotDeltaInDegrees,
                 transXHalfRange, transXDelta,
                 transYHalfRange, transYDelta,
                 useGreedyMatching);

             searchedRot.add(rotD);

             if (fits2 != null) {
                 if (fitIsBetter(bestFit, fits2[0])) {
                     bestFit = fits2[0];
                 }
             }

        }

        return bestFit;
    }

    /**
     * method to narrow down all parameter space using the given starterPoints
     * to within +- 3 degrees and +- 50 in translation.
     * The starterPoints are expected to have been the result of preSearch0.
     * <pre>
     * Runtime complexity:
     * </pre>
     * 
     * @param fits
     * @param set1
     * @param set2
     * @param useGreedyMatching
     * @return
     */
    protected TransformationPointFit preSearch1Alt2Large(TransformationPointFit[] fits,
        PairIntArray set1, PairIntArray set2, boolean useGreedyMatching) {

        if (fits == null) {
            throw new IllegalArgumentException("fits cannot be null");
        }
        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }
        
        Transformer transformer = new Transformer();

        PointPartitioner partitioner = new PointPartitioner();

        List<PairIntArray> set1Subsets = partitioner.randomSubsets(set1, 
            largeSearch0Limit);

        TransformationPointFit bestFit = null;

        // TODO: this is from preSearch00's transDelta.  needs to be accessible
        // to instance for tuning
        float tolTransX = 15;
        float tolTransY = 15;

        for (PairIntArray set1Subset : set1Subsets) {

            TransformationPointFit fit2 = preSearch1Alt2Small(fits, set1Subset,
                set2, useGreedyMatching);
            
            if (fit2 == null) {
                continue;
            }

            // apply fit to all points to evaluate
            PairFloatArray transformedSet1 = transformer.applyTransformation2(
                fit2.getParameters(), set1);

            TransformationPointFit fit;
            if (useGreedyMatching) {

                fit = evaluateFitForUnMatchedTransformedGreedy(
                    fit2.getParameters(), transformedSet1, set2,
                    tolTransX, tolTransY);

            } else {

                fit = evaluateFitForUnMatchedTransformedOptimal(
                    fit2.getParameters(), transformedSet1, set2,
                    tolTransX, tolTransY);
            }

            if (fitIsBetter(bestFit, fit)) {
                bestFit = fit;
            }
        }

        return bestFit;
    }

    /**
     * method to narrow down all parameter space using the given starterPoints
     * to within +- 20 degrees of rotation and +- 100 in translation
     * when given fits from preSearch0.
     * The starterPoints are expected to have been the result of preSearch0.
     * <pre>
     * Runtime complexity:
     * </pre>
     * 
     * @param fits
     * @param set1
     * @param set2
     * @param useGreedyMatching
     * @return
     */
    protected TransformationPointFit preSearch1Alt1(TransformationPointFit[] fits,
        PairIntArray set1, PairIntArray set2, boolean useGreedyMatching) {

        if (fits == null) {
            throw new IllegalArgumentException("fits cannot be null");
        }
        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }
        
        int nMaxMatchable = Math.min(set1.getN(), set2.getN());
        
        TransformationPointFit fit = null;
        
        // if number of points is very high, need to use alternate method
        if (nMaxMatchable > largeSearch0Limit) {
            fit = preSearch1Alt1Large(fits, set1, set2, useGreedyMatching);
        } else {
            fit = preSearch1Alt1Small(fits, set1, set2, useGreedyMatching);
        }

        return fit;
    }
    
    /**
     * method to narrow down all parameter space using the given starterPoints
     * to within +- 20 degrees of rotation and +- 100 in translation.
     * The starterPoints are expected to have been the result of preSearch0.
     * <pre>
     * Runtime complexity:
     * </pre>
     * 
     * Note that by comparison, preSearch00 has the same goal, but reduces
     * the rotation to +- 2 degrees and +- 50 in translation but currently
     * takes longer than this method.
     *
     * @param fits
     * @param set1
     * @param set2
     * @param useGreedyMatching
     * @return
     */
    protected TransformationPointFit preSearch1Alt1Small(TransformationPointFit[] fits,
        PairIntArray set1, PairIntArray set2, boolean useGreedyMatching) {

        if (fits == null) {
            throw new IllegalArgumentException("fits cannot be null");
        }
        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        boolean searchScaleToo = false;
        
        TransformationPointFit[] fits2 = 
            refineTransformationWithDownhillSimplexWrapper(fits, set1, set2,
            searchScaleToo, useGreedyMatching);

        return (fits2 != null) ? fits2[0] : null;
    }
    
    /**
     * method to narrow down all parameter space using the given starterPoints
     * to within +- 20 degrees of rotation and +- 100 in translation.
     * The starterPoints are expected to have been the result of preSearch0.
     * <pre>
     * Runtime complexity:
     * </pre>
     * 
     * Note that by comparison, preSearch00 has the same goal, but reduces
     * the rotation to +- 2 degrees and +- 50 in translation but currently
     * takes longer than this method.
     *
     * @param fits
     * @param set1
     * @param set2
     * @param useGreedyMatching
     * @return
     */
    protected TransformationPointFit preSearch1Alt1Large(TransformationPointFit[] fits,
        PairIntArray set1, PairIntArray set2, boolean useGreedyMatching) {

        if (fits == null) {
            throw new IllegalArgumentException("fits cannot be null");
        }
        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }
        
        Transformer transformer = new Transformer();

        PointPartitioner partitioner = new PointPartitioner();

        List<PairIntArray> set1Subsets = partitioner.randomSubsets(set1, 
            largeSearch0Limit);
        
        TransformationPointFit bestFit = null;
        
        boolean searchScaleToo = false;
        
        //TODO: these are adopted from preSearch0
        float tolTransX = 0.5f*100;
        float tolTransY = tolTransX;
        
        for (PairIntArray set1Subset : set1Subsets) {

            TransformationPointFit[] fits2 = 
                refineTransformationWithDownhillSimplexWrapper(fits, 
                set1Subset, set2, searchScaleToo, useGreedyMatching);
            
            if (fits2 == null) {
                continue;
            }

            // apply fit to all points to evaluate
            PairFloatArray transformedSet1 = transformer.applyTransformation2(
                fits2[0].getParameters(), set1);

            TransformationPointFit fit;
            if (useGreedyMatching) {

                fit = evaluateFitForUnMatchedTransformedGreedy(
                    fits2[0].getParameters(), transformedSet1, set2,
                    tolTransX, tolTransY);

            } else {

                fit = evaluateFitForUnMatchedTransformedOptimal(
                    fits2[0].getParameters(), transformedSet1, set2,
                    tolTransX, tolTransY);
            }

            if (fitIsBetter(bestFit, fit)) {
                bestFit = fit;
            }
        }

        return bestFit;
    }
    
    /**
     * method to narrow down all parameter space using the given starterPoints
     * to within +- 3 degrees and +- 50 in translation when given fits
     * from preSearch0.
     * The starterPoints are expected to have been the result of preSearch0.
     * <pre>
     * Runtime complexity:
     * </pre>
     * 
     * @param fits
     * @param set1
     * @param set2
     * @param useGreedyMatching
     * @return
     */
    protected TransformationPointFit preSearch1Alt2(TransformationPointFit[] fits,
        PairIntArray set1, PairIntArray set2, boolean useGreedyMatching) {

        if (fits == null) {
            throw new IllegalArgumentException("fits cannot be null");
        }
        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }
        
        int nMaxMatchable = Math.min(set1.getN(), set2.getN());
        
        TransformationPointFit fit = null;
        
        // if number of points is very high, need to use alternate method
        if (nMaxMatchable > largeSearch0Limit) {
            fit = preSearch1Alt2Large(fits, set1, set2, useGreedyMatching);
        } else {
            fit = preSearch1Alt2Small(fits, set1, set2, useGreedyMatching);
        }

        return fit;
    }
    
    protected TransformationPointFit[] 
    refineTransformationWithDownhillSimplexWrapper(TransformationPointFit[] fits,
        PairIntArray set1, PairIntArray set2, boolean searchScaleToo,
        boolean useGreedyMatching) {

        if (fits == null) {
            throw new IllegalArgumentException("fits cannot be null");
        }
        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        float centerTransX = 0;
        float centerTransY = 0;
        float minTransX = Float.MAX_VALUE;
        float minTransY = Float.MAX_VALUE;
        float maxTransX = Float.MIN_VALUE;
        float maxTransY = Float.MIN_VALUE;
        float maxTolTransX = Float.MIN_VALUE;;
        float maxTolTransY = Float.MIN_VALUE;
        int count = 0;
        for (TransformationPointFit fit : fits) {
            if (fit == null) {
                continue;
            }
            float tx = fit.getTranslationX();
            float ty = fit.getTranslationY();
            float tolX = fit.getTranslationXTolerance();
            float tolY = fit.getTranslationYTolerance();

            centerTransX += tx;
            centerTransY += ty;

            if (tx < minTransX) {
                minTransX = tx;
            }
            if (ty < minTransY) {
                minTransY = ty;
            }
            if (tx > maxTransX) {
                maxTransX = tx;
            }
            if (ty > maxTransY) {
                maxTransY = ty;
            }
            if (tolX > maxTolTransX) {
                maxTolTransX = tolX;
            }
            if (tolY > maxTolTransY) {
                maxTolTransY = tolY;
            }

            count++;
        }
        centerTransX /= (float)count;
        centerTransY /= (float)count;
        float transXHalfRange = Math.abs(maxTransX - minTransX)/2.f;
        float transYHalfRange = Math.abs(maxTransY - minTransY)/2.f;
        
        TransformationPointFit[] fits2 = refineTransformationWithDownhillSimplex(
            fits, set1, set2,
            centerTransX, centerTransY, maxTolTransX, maxTolTransY,
            transXHalfRange + 5, transYHalfRange + 5,
            searchScaleToo, useGreedyMatching, false, 50);

        return fits2;
    }

    /**
     * method to narrow down all parameter space to a fit that is within
     * +- 20 degrees in rotation and +- 50 pixels in translation in X and Y
     * of the true Euclidean transformation.
     * Note, these numbers are being learned in testing right now so may
     * change.
     *
     * @param set1
     * @param set2
     * @param scale
     * @param useGreedyMatching
     * @return
     */
    protected TransformationPointFit preSearch(PairIntArray set1,
        PairIntArray set2, float scale, boolean useGreedyMatching) {

        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        TransformationPointFit[] fits = preSearch0(set1, set2, scale, useGreedyMatching);
        
        if (fits == null) {
            return null;
        }

        TransformationPointFit fit2 = preSearch1Alt1(fits, set1, set2, useGreedyMatching);

        if (fit2 == null) {
            return null;
        }

        return fit2;
    }
        
    protected TransformationPointFit[] boundedGridSearch(
        PairIntArray set1, PairIntArray set2, TransformationParameters params,
        float scaleHalfRange, float scaleDelta,
        float rotHalfRangeInDegrees, float rotDeltaInDegrees,
        float transXHalfRange, float transXDelta,
        float transYHalfRange, float transYDelta,
        boolean useGreedyMatching) {

        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        int nMaxMatchable = Math.min(set1.getN(), set2.getN());

        float tolTransX = 0.5f*transXDelta + (float)(Math.sin(rotDeltaInDegrees*Math.PI/180.)*transXDelta);
        float tolTransY = 0.5f*transYDelta + (float)(Math.sin(rotDeltaInDegrees*Math.PI/180.)*transYDelta);

        Transformer transformer = new Transformer();

        float[] rotation = MiscMath.writeDegreeIntervals(
            params.getRotationInDegrees() - rotHalfRangeInDegrees,
            params.getRotationInDegrees() + rotHalfRangeInDegrees,
            rotDeltaInDegrees);

        log.fine("searching rotation: " + Arrays.toString(rotation));

        int nX = (int) Math.ceil(2*transXHalfRange / transXDelta);
        int nY = (int) Math.ceil(2*transYHalfRange / transYDelta);
        if (nX == 0) {
            nX = 1;
        }
        if (nY == 0) {
            nY = 1;
        }

        float[] txS = new float[nX];
        txS[0] = params.getTranslationX() - transXHalfRange;
        for (int i = 1; i < nX; ++i) {
            txS[i] = txS[i - 1] + transXDelta;
        }
        float[] tyS = new float[nY];
        tyS[0] = params.getTranslationY() - transYHalfRange;
        for (int i = 1; i < nY; ++i) {
            tyS[i] = tyS[i - 1] + transYDelta;
        }
        
        int nS = (int) Math.ceil(2*scaleHalfRange / scaleDelta);
        if (nS == 0) {
            nS = 1;
        }
        float[] scales = new float[nS];
        scales[0] = params.getScale() - scaleHalfRange;
        for (int i = 1; i < nS; ++i) {
            scales[i] = scales[i - 1] + scaleDelta;
        }
        
        int nKeep = 10;

        long t0 = System.currentTimeMillis();

        FixedSizeSortedVector<TransformationPointFit> starterPoints
            = new FixedSizeSortedVector<>(nKeep, TransformationPointFit.class);

        float[] xsr = new float[set1.getN()];
        float[] ysr = new float[set1.getN()];

        for (float scale : scales) {
            for (float rotDeg : rotation) {

                TransformationParameters rotScaleParams =
                    new TransformationParameters();
                rotScaleParams.setScale(scale);
                rotScaleParams.setRotationInDegrees(rotDeg);
                rotScaleParams.setTranslationX(0);
                rotScaleParams.setTranslationY(0);
                rotScaleParams.setOriginX(params.getOriginX());
                rotScaleParams.setOriginY(params.getOriginY());

                transformer.applyTransformation(rotScaleParams, set1, xsr, ysr);

                for (float tx : txS) {
                    for (float ty : tyS) {

                        PairFloatArray allPoints1Tr = new PairFloatArray(set1.getN());

                        for (int i = 0; i < xsr.length; ++i) {
                            allPoints1Tr.add(xsr[i] + tx, ysr[i] + ty);
                        }

                        TransformationParameters params2 = new TransformationParameters();
                        params2.setScale(scale);
                        params2.setRotationInDegrees(rotDeg);
                        params2.setTranslationX(tx);
                        params2.setTranslationY(ty);
                        params2.setOriginX(params.getOriginX());
                        params2.setOriginY(params.getOriginY());

                        TransformationPointFit fit2;
                        if (useGreedyMatching) {
                            fit2 = evaluateFitForUnMatchedTransformedGreedy(params2,
                                allPoints1Tr, set2, tolTransX, tolTransY);
                        } else {
                            fit2 = evaluateFitForUnMatchedTransformedOptimal(params2,
                                allPoints1Tr, set2, tolTransX, tolTransY);
                        }

                        starterPoints.add(fit2);
                    }
                }
            }
        }

        if (starterPoints.getArray()[0] == null) {
            return null;
        }

        long tm = System.currentTimeMillis() - t0;
        double ts = tm * 1e-3;
        log.info("bounded grid search finished for nMaxMatchable=" + nMaxMatchable +
            " seconds=" + ts);

        if (starterPoints.getNumberOfItems() > 0) {
            log.info("best=" + starterPoints.getArray()[0].toString());
        }

        return starterPoints.getArray();
    }

    /**
     * use a downhill simplex algorithm to search for a better fit to the data
     * given the starterPoints.  Note that the method only changes
     * rotation, translation in X, and translation in Y.
     * 
     * @param set1
     * @param set2
     * @param params
     * @param rotHalfRangeInDegrees half of the range in rotation to search
     * around params.rotationInDegrees.  the value is added and subtracted
     * from params.rotationInDegrees to create the start and end of the
     * range.
     * @param rotDeltaInDegrees the interval size of rotation in degrees
     * for the search.
     * @param transXHalfRange half of the range in translation in X to search
     * around params.translationX.  the value is added and subtracted
     * from params.translationX to create the start and end of the range.
     * @param transYHalfRange half of the range in translation in Y to search
     * around params.translationY.  the value is added and subtracted
     * from params.translationY to create the start and end of the range.
     * @param transYDelta the interval size in translation of X for the search
     * @param transXDelta the interval size in translation of Y for the search
     * @param useGreedyMatching if true, uses an N^2 algorithm to find the
     * best match made for each point in order of points in set1, else uses
     * an ~N^3 optimal bipartite matching algorithm.
     * @return
     */
    public TransformationPointFit refineTheTransformation(
        TransformationParameters params, PairIntArray set1, PairIntArray set2,
        float rotHalfRangeInDegrees, float rotDeltaInDegrees, 
        float transXHalfRange, float transXDelta, 
        float transYHalfRange, float transYDelta,
        boolean useGreedyMatching) {

        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
        if (set1 == null) {
            throw new IllegalArgumentException("set1 cannot be null");
        }
        if (set2 == null) {
            throw new IllegalArgumentException("set2 cannot be null");
        }

        float scaleHalfRange = 0;
        float scaleDelta = 1;
        
        TransformationPointFit[] fits = boundedGridSearch(
            set1, set2, params,
            scaleHalfRange, scaleDelta,
            rotHalfRangeInDegrees, rotDeltaInDegrees,
            transXHalfRange, transXDelta,
            transYHalfRange, transYDelta, useGreedyMatching);

        if (fits == null) {
            return null;
        }

        int nMaxMatchable = Math.min(set1.getN(), set2.getN());

        if (hasConverged(fits[0], nMaxMatchable)) {
            return fits[0];
        }
        
        // finest grid so no need for downhill simplex search
        if ((rotDeltaInDegrees == 1) && (transXDelta == 1) && (transYDelta == 1)) {
            return fits[0];
        }
        
        TransformationPointFit fit2 = preSearch1Alt1(fits, set1, set2,
            useGreedyMatching);

        if (fit2 == null) {
            return null;
        }

        if (hasConverged(fit2, nMaxMatchable)) {
            return fit2;
        }

        // TODO: may need a small fine grid search here

        return fit2;
        
    }
    
    /**
     * use a downhill simplex algorithm to search for a better fit to the data
     * given the starterPoints.  This method changes
     * scale, rotation, translation in X, and translation in Y.
     * 
     * @param set1
     * @param set2
     * @param params
     * @param scaleHalfRange half of the range in scale to search
     * around params.scale.  the value is added and subtracted
     * from params.scale to create the start and end of the
     * range.
     * @param scaleDelta the interval size of scale for the search.
     * @param rotHalfRangeInDegrees half of the range in rotation to search
     * around params.rotationInDegrees.  the value is added and subtracted
     * from params.rotationInDegrees to create the start and end of the
     * range.
     * @param rotDeltaInDegrees the interval size of rotation in degrees
     * for the search.
     * @param transXHalfRange half of the range in translation in X to search
     * around params.translationX.  the value is added and subtracted
     * from params.translationX to create the start and end of the range.
     * @param transYHalfRange half of the range in translation in Y to search
     * around params.translationY.  the value is added and subtracted
     * from params.translationY to create the start and end of the range.
     * @param transYDelta the interval size in translation of X for the search
     * @param transXDelta the interval size in translation of Y for the search
     * @param useGreedyMatching if true, uses an N^2 algorithm to find the
     * best match made for each point in order of points in set1, else uses
     * an ~N^3 optimal bipartite matching algorithm.
     * @return
     */
    public TransformationPointFit refineTheTransformation(
        TransformationParameters params, PairIntArray set1, PairIntArray set2,
        float scaleHalfRange, float scaleDelta, 
        float rotHalfRangeInDegrees, float rotDeltaInDegrees,
        float transXHalfRange, float transXDelta, 
        float transYHalfRange, float transYDelta,
        boolean useGreedyMatching) {

        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
        if (set1 == null) {
            throw new IllegalArgumentException("set1 cannot be null");
        }
        if (set2 == null) {
            throw new IllegalArgumentException("set2 cannot be null");
        }

        TransformationPointFit[] fits = boundedGridSearch(
            set1, set2, params,
            scaleHalfRange, scaleDelta,
            rotHalfRangeInDegrees, rotDeltaInDegrees,
            transXHalfRange, transXDelta,
            transYHalfRange, transYDelta, useGreedyMatching);

        if (fits == null) {
            return null;
        }

        int nMaxMatchable = Math.min(set1.getN(), set2.getN());

        if (hasConverged(fits[0], nMaxMatchable)) {
            return fits[0];
        }
        
        boolean searchScaleToo = true;
        
        TransformationPointFit[] fits2 = 
            refineTransformationWithDownhillSimplexWrapper(fits, set1, set2,
            searchScaleToo, useGreedyMatching);
        
        if (fits2 == null) {
            return null;
        }

        if (hasConverged(fits2[0], nMaxMatchable)) {
            return fits2[0];
        }

        // TODO: may need a small fine grid search here

        return fits2[0];
        
    }

    /**
     * for use when the number of data points is too high for reasonably
     * fast solution using preSearch0.
     * The method randomly reduces set1 to a smaller number of point sets
     * and finds the best solution among those.
     * @param set1
     * @param set2
     * @param scale
     * @param useGreedyMatch
     * @return
     */
    protected TransformationPointFit[] preSearch0Large(PairIntArray set1,
        PairIntArray set2, float scale, boolean useGreedyMatch) {

        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }
        
        Transformer transformer = new Transformer();

        PointPartitioner partitioner = new PointPartitioner();

        List<PairIntArray> set1Subsets = partitioner.randomSubsets(set1, 
            largeSearch0Limit);

        TransformationPointFit bestFit = null;

        TransformationPointFit[] bestFits = null;

        // TODO: this is from preSearch0's transDelta.  needs to be accessible
        // to instance for tuning
        float tolTransX = 100;
        float tolTransY = 100;

        for (PairIntArray set1Subset : set1Subsets) {

            TransformationPointFit[] fits = preSearch0(set1Subset, set2, scale,
                useGreedyMatch);

            if (fits == null) {
                continue;
            }

            // apply fit to all points to evaluate
            PairFloatArray transformedSet1 = transformer.applyTransformation2(
                fits[0].getParameters(), set1);

            TransformationPointFit fit;
            if (useGreedyMatch) {

                fit = evaluateFitForUnMatchedTransformedGreedy(
                    fits[0].getParameters(), transformedSet1, set2,
                    tolTransX, tolTransY);

            } else {

                fit = evaluateFitForUnMatchedTransformedOptimal(
                    fits[0].getParameters(), transformedSet1, set2,
                    tolTransX, tolTransY);
            }

            if (fitIsBetter(bestFit, fit)) {
                bestFit = fit;
                bestFits = fits;
            }
        }

        return bestFits;
    }

}
