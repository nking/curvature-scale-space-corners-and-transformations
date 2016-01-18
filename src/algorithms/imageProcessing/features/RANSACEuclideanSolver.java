package algorithms.imageProcessing.features;

import algorithms.imageProcessing.transform.EuclideanEvaluator;
import algorithms.imageProcessing.transform.EuclideanTransformationFit;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.util.RANSACAlgorithmIterations;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.List;
import java.util.logging.Logger;

/**
 * uses the RANSAC algorithm on combinations of 3 points to calculate
 * euclidean transformations and evaluate them against all points.
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
public class RANSACEuclideanSolver {

    private boolean debug = true;

    private Logger log = Logger.getLogger(this.getClass().getName());

    /**
     * calculate the euclidean projection among the given points with the
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
    public EuclideanTransformationFit calculateEuclideanTransformation(
        PairIntArray matchedLeftXY, PairIntArray matchedRightXY,
        PairIntArray outputLeftXY, PairIntArray outputRightXY)
        throws NoSuchAlgorithmException {

        if (matchedLeftXY == null) {
            throw new IllegalArgumentException("matchedLeftXY cannot be null");
        }
        if (matchedRightXY == null) {
            throw new IllegalArgumentException("matchedRightXY cannot be null");
        }
        if (matchedLeftXY.getN() < 2) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
            "the algorithms require 2 or more points.  matchedLeftXY.n=" 
            + matchedLeftXY.getN());
        }

        /*
        -- randomly sample 2 points from nPoints
        -- calculate the transformation from the 2
        -- evaluate all nPoints against the transformation and keep each point 
           which has error < tolerance.   those points are a consensus of this 
           model (the fundamental matrix)
        -- if the number of points in the consensus is > required, exit loop,
           else repeat (note that each consensus is stored separately for
           each iteration.  note that the loop is terminated if the number of
           iterations has exceeded a predetermined maximum.
        -- after exit from the loop, the largest consensus is used to
           re-calculate the fundamental matrix as the last result.
           note that if there aren't enough points in a consensus to
           calculate the transformation, the result is null.
        */

        int nSet = 2;

        int nPoints = matchedLeftXY.getN();
        
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();
        
        EuclideanEvaluator evaluator = new EuclideanEvaluator();
            
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        log.info("SEED=" + seed + " nPoints=" + nPoints);
        sr.setSeed(seed);

        int tolerance = 5;

        // consensus indexes
        EuclideanTransformationFit bestFit = null;
        
        RANSACAlgorithmIterations nEstimator = new RANSACAlgorithmIterations();

        /*
        could consider a threshold max iteration based upon the image size such
        as in (http://phototour.cs.washington.edu/ModelingTheWorld_ijcv07.pdf)
        which uses 0.6% of the maximum image dimension.
        */
        long nMaxIter = nEstimator.estimateNIterFor99PercentConfidence(nPoints, 7, 0.4);

        if (nPoints == 2) {
            nMaxIter = 1;
        }

        int nIter = 0;
        
        int[] selectedIndexes = new int[nSet];
        
        PairIntArray sampleLeft = new PairIntArray(nSet);
        PairIntArray sampleRight = new PairIntArray(nSet);
        for (int i = 0; i < nSet; ++i) {
            sampleLeft.add(0, 0);
            sampleRight.add(0, 0);
        }
        assert(sampleLeft.getN() == nSet);
        
        while ((nIter < nMaxIter) && (nIter < 2000)) {
            
            MiscMath.chooseRandomly(sr, selectedIndexes, nPoints);
            
            for (int i = 0; i < selectedIndexes.length; ++i) {
                int idx = selectedIndexes[i];
                sampleLeft.set(i, matchedLeftXY.getX(idx), matchedLeftXY.getY(idx));
                sampleRight.set(i, matchedRightXY.getX(idx), matchedRightXY.getY(idx));
            }
            assert(sampleRight.getN() == nSet);
            
            TransformationParameters params = tc.calulateEuclideanWithoutFilter(
                sampleLeft, sampleRight, 0, 0);

            if (params == null) {
                nIter++;
                continue;
            }

            EuclideanTransformationFit fit = evaluator.evaluate(matchedLeftXY,
                matchedRightXY, params, tolerance);
                    
            if (fit == null) {
                nIter++;
                continue;
            }
               
            if (fit.isBetter(bestFit)) {
                bestFit = fit;
            }
            
            nIter++;
            
            // recalculate nMaxIter
            if ((bestFit != null) && ((nIter % 50) == 0)) {
                
                double ratio = (double)bestFit.getInlierIndexes().size()
                    /(double)matchedLeftXY.getN();
                
                nMaxIter = nEstimator.estimateNIterFor99PercentConfidence(nPoints, 
                    nSet, ratio);
            }
        }

        if (bestFit == null || bestFit.getInlierIndexes().isEmpty()) {
            log.info("no solution.  nIter=" + nIter);
            return null;
        } else {
            log.info("bestFit before consensus = " + bestFit.getTransformationParameters().toString());
        }

        // store inliers in outputLeftXY and outputRightXY and redo the
        // entire fit using only the inliers to determine the fundamental
        // matrix.
        int n = bestFit.getInlierIndexes().size();
        PairIntArray inliersLeftXY = new PairIntArray(n);
        PairIntArray inliersRightXY = new PairIntArray(n);
        int count = 0;
        for (Integer idx : bestFit.getInlierIndexes()) {            
            int idxInt = idx.intValue(); 
            inliersLeftXY.add(matchedLeftXY.getX(idxInt), matchedLeftXY.getY(idxInt));
            inliersRightXY.add(matchedRightXY.getX(idxInt), matchedRightXY.getY(idxInt));            
            count++;
        }
        assert(n == count);
        assert(n == inliersLeftXY.getN());
        
        TransformationParameters params = tc.calulateEuclideanWithoutFilter(
            inliersLeftXY, inliersRightXY, 0, 0);
        
        log.info("consensusParams=" + params.toString());

        EuclideanTransformationFit consensusFit = evaluator.evaluate(matchedLeftXY,
            matchedRightXY, params, tolerance);
        
        if (consensusFit == null || consensusFit.getInlierIndexes().isEmpty()) {
            log.info("no consensus possible for given points");
            return null;
        } else {
            log.info("consensus nEval=" + consensusFit.getInlierIndexes().size());
        }
        
        // inlierIndexes are w.r.t matchedLeftXY
        List<Integer> inlierIndexes = consensusFit.getInlierIndexes();
        // write to output and convert the coordinate indexes to the original point indexes
        for (int i = 0; i < inlierIndexes.size(); ++i) {
            Integer index = inlierIndexes.get(i);
            int idx = index.intValue();
            assert(idx < matchedLeftXY.getN());
            outputLeftXY.add(matchedLeftXY.getX(idx), matchedLeftXY.getY(idx));
            outputRightXY.add(matchedRightXY.getX(idx), matchedRightXY.getY(idx));
        }
        
        log.info("nIter=" + nIter);

        log.info("final fit: " + consensusFit.toString());

        return consensusFit;
    }

}
