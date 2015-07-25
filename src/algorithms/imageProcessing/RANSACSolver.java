package algorithms.imageProcessing;

import algorithms.imageProcessing.util.RANSACAlgorithmIterations;
import algorithms.misc.MiscMath;
import algorithms.util.PairFloatArray;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.List;
import java.util.logging.Logger;
import org.ejml.simple.SimpleMatrix;

/**
 * given matched point lists, determine the best epipolar solution using a
 * 7-point epipolar calculation and random draws of 7 points from the
 * matched point lists under the assumption that some of the matched points
 * are not true (correct) matches.
 *
 * <pre>
 * useful reading:
 * http://6.869.csail.mit.edu/fa12/lectures/lecture13ransac/lecture13ransac.pdf
 * and
 * http://www.dtic.mil/dtic/tr/fulltext/u2/a460585.pdf
 * </pre>
 *
 * @author nichole
 */
public class RANSACSolver {

    private boolean debug = true;

    private double generalTolerance = 5;

    private Logger log = Logger.getLogger(this.getClass().getName());

    /**
     * calculate the epipolar projection among the given points with the
     * assumption that some of the points in the matched lists are not
     * true matches.
     *
     * @param matchedLeftXY
     * @param matchedRightXY
     * @param outputLeftXY
     * @param outputRightXY
     * @return
     * @throws NoSuchAlgorithmException
     */
    public StereoProjectionTransformerFit calculateEpipolarProjection(
        PairFloatArray matchedLeftXY, PairFloatArray matchedRightXY,
        PairFloatArray outputLeftXY, PairFloatArray outputRightXY)
        throws NoSuchAlgorithmException {

        if (matchedLeftXY == null) {
            throw new IllegalArgumentException("matchedLeftXY cannot be null");
        }
        if (matchedRightXY == null) {
            throw new IllegalArgumentException("matchedRightXY cannot be null");
        }
        if (matchedLeftXY.getN() < 7) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the algorithms require 7 or more points."
                + " matchedLeftXY.n=" + matchedLeftXY.getN());
        }

        SimpleMatrix input1 =
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(matchedLeftXY);

        SimpleMatrix input2 =
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(matchedRightXY);

        return calculateEpipolarProjection(input1, input2,
            outputLeftXY, outputRightXY);
    }

    /**
     * calculate the epipolar projection among the given points with the
     * assumption that some of the points in the matched lists are not
     * true matches.
     *
     * @param matchedLeftXY
     * @param matchedRightXY
     * @param outputLeftXY
     * @param outputRightXY
     * @return
     * @throws NoSuchAlgorithmException
     */
    public StereoProjectionTransformerFit calculateEpipolarProjection(
        SimpleMatrix matchedLeftXY, SimpleMatrix matchedRightXY,
        PairFloatArray outputLeftXY, PairFloatArray outputRightXY)
        throws NoSuchAlgorithmException {

        if (matchedLeftXY == null) {
            throw new IllegalArgumentException("matchedLeftXY cannot be null");
        }
        if (matchedRightXY == null) {
            throw new IllegalArgumentException("matchedRightXY cannot be null");
        }
        if (matchedLeftXY.numCols() < 7) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the algorithms require 7 or more points."
                + " matchedLeftXY.n=" + matchedLeftXY.numCols());
        }
        if (matchedLeftXY.numCols() != matchedRightXY.numCols()) {
            throw new IllegalArgumentException(
                "matchedLeftXY and right bmust be the same size");
        }

        /*
        -- randomly sample 7 points from nPoints
        -- calculate the epipolar fundamental matrix from the 7
        -- evaluate all nPoints against the epipolar projected lines
           and keep each point which has error < tolerance.
           those points are a consensus of this model (the fundamental matrix)
        -- if the number of points in the consensus is > required, exit loop,
           else repeat (note that each consensus is stored separately for
           each iteration.  note that the loop is terminated if the number of
           iterations has exceeded a predetermined maximum.
        -- after exit from the loop, the largest consensus is used to
           re-calculate the fundamental matrix as the last result.
           note that if there aren't enough points in a consensus to
           calculate the fundamental matrix, the result is null.

           an improvement on the algorithm is to adaptively determine the
           number of iterations necessary to find a large enough consensus.

        ---------
        Note, the random selection would ideally draw one combination of n!/(k!*(n-k)!))
        possible subsets.  That requires an enumeration of each combination findable
        by an O(1) operation.  Gosper's hack returns the unique subsets in an ordered
        manner, but a way to access that (i.e. "random access") would need to be made to have
        truly uniform random selection.

        Below, am using k draws of random numbers to create a combination of unique indexes
        from n indexes.  The probability of selecting one unique combination within all
        possible combinations is 1./(n!/(k!(n-k)!)).  Although that's for a single draw
        rather than composed of k draws, it's close enough to see that the probability
        of drawing the same combination is low, so have chosen to allow the possible repeat
        evaluation of a set of numbers for now but the sort of the 7 numbers to
        check for unique combination is considered in the calculation below.

        Note, that a non-uniform random selection algorithm could be made by using
        BigInteger and drawing random numbers between the low value of k bits set to 1 and
        the high value of n bits with the highest k bits set to 1.
        Upon each random selection within that range, one could find the next lowest number
        with k bits set to 1 (find lowest set bit, change it to 0 and set the 2 below it to
        1 etc to total k set bits).
        If that bitstring has already been selected randomly, can use Gosper's hack
        to return the next k=7 bitstring.
        The universe of numbers drawn from is larger than the total number of
        possible combinations, so the process is a non-uniform random selection.

        The current implementation of drawing k=7 unique indexes randomly and then sorting
        their order is 7 random operations * 7 * log_2(7) steps so 7 * 7 * 20 steps for a single
        iteration of selecting a subset.
          The BigInteger random suggestion above in contrast would be the
          BigInteger random operation plus up to 6 set bit operations plus possibly
          half a dozen bit operations from Gosper's hack if need to
          select again, so the BigInteger implementation might be useful.

        */

        int nSet = 7;

        int nPoints = matchedLeftXY.numCols();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        log.info("SEED=" + seed);
        sr.setSeed(seed);

        int[] selected = new int[nSet];

        PairFloatArray xy1 = null;
        PairFloatArray xy2 = null;

        double tolerance = generalTolerance;

        // ====== find matched pairs incremently, and best among only 7 point match ====

        StereoProjectionTransformerFit bestFit = null;
        List<Integer> inlierIndexes = null;

        RANSACAlgorithmIterations nEstimator = new RANSACAlgorithmIterations();

        int nMaxIter = nEstimator.estimateNIterForFiftyPercentOutliersFor7Points(
            nSet);

        float convergence = 0.5f * matchedLeftXY.numCols();
        if (convergence < 7) {
            convergence = 7;
        }

        if (nPoints == 7) {
            nMaxIter = 1;
        }

        int nIter = 0;

        for (int i = 0; i < nMaxIter; i++) {

            MiscMath.chooseRandomly(sr, selected, nPoints);

            xy1 = new PairFloatArray();
            xy2 = new PairFloatArray();
            for (int j = 0; j < selected.length; j++) {
                int idx = selected[j];
                xy1.add(
                    (float)matchedLeftXY.get(0, idx),
                    (float)matchedLeftXY.get(1, idx));
                xy2.add(
                    (float)matchedRightXY.get(0, idx),
                    (float)matchedRightXY.get(1, idx));
            }

            StereoProjectionTransformer spTransformer = new
                StereoProjectionTransformer();

            // determine matrix from 7 points.
            List<SimpleMatrix> fms =
                spTransformer.calculateEpipolarProjectionFor7Points(xy1, xy2);

            if (fms == null || fms.isEmpty()) {
                nIter++;
                continue;
            }

            SimpleMatrix fm = fms.get(0);

            // evaluate fit against all points
            StereoProjectionTransformerFit fit =
                spTransformer.evaluateFitForAlreadyMatched(
                fm, matchedLeftXY, matchedRightXY, tolerance);

            List<Integer> indexes = fit.getInlierIndexes();

            nIter++;

            if (fitIsBetter(bestFit, fit)) {
                
                bestFit = fit;
                inlierIndexes = indexes;

                if (false && hasConverged(bestFit) && 
                    (inlierIndexes.size() >= convergence)) {
                    
                    break;
                }
            }
        }

        if (bestFit == null) {
            return null;
        }

        // store inliers in outputLeftXY and outputRightXY and redo the
        // entire fit using only the inliers to determine the fundamental
        // matrix.
        SimpleMatrix inliersLeftXY = new SimpleMatrix(3, inlierIndexes.size());
        SimpleMatrix inliersRightXY = new SimpleMatrix(3, inlierIndexes.size());

        int count = 0;
        for (Integer idx : inlierIndexes) {
            int idxInt = idx.intValue();
            
            float x = (float)matchedLeftXY.get(0, idxInt);
            float y = (float)matchedLeftXY.get(1, idxInt);
            inliersLeftXY.set(0, count, x);
            inliersLeftXY.set(1, count, y);
            inliersLeftXY.set(2, count, 1);
            
            x = (float)matchedRightXY.get(0, idxInt);
            y = (float)matchedRightXY.get(1, idxInt);
            inliersRightXY.set(0, count, x);
            inliersRightXY.set(1, count, y);
            inliersRightXY.set(2, count, 1);
            count++;
        }

        StereoProjectionTransformer spTransformer = new
            StereoProjectionTransformer();

        StereoProjectionTransformerFit finalFit = null;

        if (nPoints == 7) {

            finalFit = bestFit;

        } else {

            SimpleMatrix finalFM = null;
            if (inliersLeftXY.numCols() == 7) {
                List<SimpleMatrix> fms =
                    spTransformer.calculateEpipolarProjectionFor7Points(
                    inliersLeftXY, inliersRightXY);
                if (fms == null || fms.isEmpty()) {
                    return null;
                }
                finalFM = fms.get(0);
            } else {
                finalFM =
                    spTransformer.calculateEpipolarProjectionForPerfectlyMatched(
                    inliersLeftXY, inliersRightXY);
            }

            //fit against all because the new fit may be slighty different
            finalFit = spTransformer.evaluateFitForAlreadyMatched(finalFM,
                matchedLeftXY, matchedRightXY, tolerance);

        }
        
        for (Integer idx : finalFit.getInlierIndexes()) {
            int idxInt = idx.intValue();
            outputLeftXY.add(
                (float)matchedLeftXY.get(0, idxInt),
                (float)matchedLeftXY.get(1, idxInt));
            outputRightXY.add(
                (float)matchedRightXY.get(0, idxInt),
                (float)matchedRightXY.get(1, idxInt));
        }

        log.info("nIter=" + nIter);

        log.info("best fit from 7-point: " + bestFit.toString());

        log.info("final fit: " + finalFit.toString());

        return finalFit;
    }

    boolean fitIsBetter(StereoProjectionTransformerFit bestFit,
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

    private boolean hasConverged(StereoProjectionTransformerFit bestFit) {
        if (bestFit == null) {
            return false;
        }
        
        double limit = 0.2;

        if ((bestFit.getMeanDistance() < limit) && bestFit.getStDevFromMean() < limit) {
            return true;
        }
        return false;
    }
}
