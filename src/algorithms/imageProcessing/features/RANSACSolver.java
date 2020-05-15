package algorithms.imageProcessing.features;

import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.EpipolarTransformer;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.imageProcessing.util.RANSACAlgorithmIterations;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;
import no.uib.cipr.matrix.DenseMatrix;

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

    private Logger log = Logger.getLogger(this.getClass().getName());

    /**
     * NOT COMPLETELY READY FOR USE YET
     * calculate the epipolar projection among the given points with the
     * assumption that some of the points in the matched lists are not
     * true matches.
     *
     * @param matchedLeftXY
     * @param matchedRightXY
     * @param outputLeftXY
     * @param outputRightXY
     * @return
     */
    public EpipolarTransformationFit calculateEpipolarProjection(
        PairIntArray matchedLeftXY, PairIntArray matchedRightXY,
        PairIntArray outputLeftXY, PairIntArray outputRightXY) {

        if (matchedLeftXY == null) {
            throw new IllegalArgumentException("matchedLeftXY cannot be null");
        }
        if (matchedRightXY == null) {
            throw new IllegalArgumentException("matchedRightXY cannot be null");
        }
        if (matchedLeftXY.getN() < 7) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
            "the algorithms require 7 or more points.  matchedLeftXY.n=" 
            + matchedLeftXY.getN());
        }

        EpipolarTransformer spTransformer = new EpipolarTransformer();
        
        DenseMatrix input1 =
            spTransformer.rewriteInto3ColumnMatrix(matchedLeftXY);

        DenseMatrix input2 =
            spTransformer.rewriteInto3ColumnMatrix(matchedRightXY);
        
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
     */
    public EpipolarTransformationFit calculateEpipolarProjection(
        DenseMatrix matchedLeftXY, DenseMatrix matchedRightXY,
        PairIntArray outputLeftXY, PairIntArray outputRightXY) {

        //TODO: improve this method
        
        if (matchedLeftXY == null) {
            throw new IllegalArgumentException("matchedLeftXY cannot be null");
        }
        if (matchedRightXY == null) {
            throw new IllegalArgumentException("matchedRightXY cannot be null");
        }
        if (matchedLeftXY.numColumns() < 7) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the algorithms require 7 or more points."
                + " matchedLeftXY.n=" + matchedLeftXY.numColumns());
        }
        if (matchedLeftXY.numColumns() != matchedRightXY.numColumns()) {
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
        */

        int nSet = 7;

        int nPoints = matchedLeftXY.numColumns();
        
        ErrorType errorType = ErrorType.SAMPSONS;

        EpipolarTransformer spTransformer = new EpipolarTransformer();
            
        SecureRandom sr = Misc.getSecureRandom();
        long seed = System.currentTimeMillis();
        log.info("SEED=" + seed + " nPoints=" + nPoints);
        sr.setSeed(seed);

        int tolerance = 5;

        // consensus indexes
        EpipolarTransformationFit bestFit = null;
        
        RANSACAlgorithmIterations nEstimator = new RANSACAlgorithmIterations();

        /*
        could consider a threshold max iteration based upon the image size such
        as in (http://phototour.cs.washington.edu/ModelingTheWorld_ijcv07.pdf)
        which uses 0.6% of the maximum image dimension.
        */
        long nMaxIter = nEstimator.estimateNIterFor99PercentConfidence(nPoints, 7, 0.5);

        if (nPoints == 7) {
            nMaxIter = 1;
        }
        System.out.println("nPoints=" + nPoints + " estimate for nMaxIter=" +
            nMaxIter);

        int nIter = 0;
        
        int[] selectedIndexes = new int[nSet];
        
        DenseMatrix sampleLeft = new DenseMatrix(3, nSet);
        DenseMatrix sampleRight = new DenseMatrix(3, nSet);
        
        while ((nIter < nMaxIter) && (nIter < 2000)) {
            
            MiscMath.chooseRandomly(sr, selectedIndexes, nPoints);

            int count = 0;
            
            for (int bitIndex : selectedIndexes) {

                int idx = bitIndex;

                sampleLeft.set(0, count, matchedLeftXY.get(0, idx));
                sampleLeft.set(1, count, matchedLeftXY.get(1, idx));
                sampleLeft.set(2, count, 1);
                                
                sampleRight.set(0, count, matchedRightXY.get(0, idx));
                sampleRight.set(1, count, matchedRightXY.get(1, idx));
                sampleRight.set(2, count, 1);
                
                count++;
            }

            // determine matrix from 7 points.
            List<DenseMatrix> fms =
                spTransformer.calculateEpipolarProjectionFor7Points(sampleLeft, 
                    sampleRight);

            if (fms == null || fms.isEmpty()) {
                nIter++;
                continue;
            }

            // use point dist to epipolar lines to estimate errors of sample
            EpipolarTransformationFit fit = null;
            
            for (DenseMatrix fm : fms) {
                EpipolarTransformationFit fitI = 
                    spTransformer.calculateError(fm, matchedLeftXY, 
                        matchedRightXY, errorType, tolerance);
                
                if (fitI.isBetter(fit)) {
                    fit = fitI;
                }
            }
            
            if (fit == null) {
                nIter++;
                continue;
            }
            
            if (fit.isBetter(bestFit)) {
                bestFit = fit;
            }
            
            nIter++;
            
            // recalculate nMaxIter
            if ((bestFit != null) && ((nIter % 10) == 0)) {
                double ratio = (double)bestFit.getInlierIndexes().size()
                    /(double)matchedLeftXY.numColumns();
                if (ratio >= 0.0000001 && (ratio <= 1.0)) {
                    nMaxIter = nEstimator.estimateNIterFor99PercentConfidence(nPoints, 
                        nSet, ratio);
                }
            }
        }

        if (bestFit == null || bestFit.getInlierIndexes().isEmpty()) {
            log.info("no solution.  nIter=" + nIter);
            return null;
        }

        // store inliers in outputLeftXY and outputRightXY and redo the
        // entire fit using only the inliers to determine the fundamental
        // matrix.
        int n = bestFit.getInlierIndexes().size();
        DenseMatrix inliersLeftXY = new DenseMatrix(3, n);
        DenseMatrix inliersRightXY = new DenseMatrix(3, n);
        
        int count = 0;
        for (Integer idx : bestFit.getInlierIndexes()) {
            int idxInt = idx.intValue();            
            inliersLeftXY.set(0, count, matchedLeftXY.get(0, idxInt));
            inliersLeftXY.set(1, count, matchedLeftXY.get(1, idxInt));
            inliersLeftXY.set(2, count, 1);
            inliersRightXY.set(0, count, matchedRightXY.get(0, idxInt));
            inliersRightXY.set(1, count, matchedRightXY.get(1, idxInt));
            inliersRightXY.set(2, count, 1);
            count++;
        }

        EpipolarTransformationFit consensusFit = null;
        
        if (inliersRightXY.numColumns() == 7) {
            
            List<DenseMatrix> fms = spTransformer.calculateEpipolarProjectionFor7Points(
                inliersLeftXY, inliersRightXY);
            if (fms == null || fms.isEmpty()) {
                return null;
            }
            EpipolarTransformationFit fit = null;
            for (DenseMatrix fm : fms) {
                EpipolarTransformationFit fitI = 
                    spTransformer.calculateError(fm, matchedLeftXY, 
                        matchedRightXY, errorType, tolerance);
                if (fitI.isBetter(fit)) {
                    fit = fitI;
                }
            }
            consensusFit = fit;
            
        } else {
            
            DenseMatrix fm = 
                spTransformer.calculateEpipolarProjection(
                inliersLeftXY, inliersRightXY);
            
            EpipolarTransformationFit fit = 
                spTransformer.calculateError(fm, matchedLeftXY, 
                    matchedRightXY, errorType, tolerance);
            
            consensusFit = fit;
        }
    
        // write to output and convert the coordinate indexes to the original point indexes
        List<Integer> inlierIndexes = consensusFit.getInlierIndexes();
        for (int i = 0; i < inlierIndexes.size(); ++i) {
            Integer index = inlierIndexes.get(i);
            int idx = index.intValue();
            outputLeftXY.add(
                (int)Math.round(matchedLeftXY.get(0, idx)),
                (int)Math.round(matchedLeftXY.get(1, idx)));
            outputRightXY.add(
                (int)Math.round(matchedRightXY.get(0, idx)),
                (int)Math.round(matchedRightXY.get(1, idx)));
        }
        
        /*
        to approximate the error in using this fundamental matrix to estimate
        distances to epipolar lines,
        will calculate distances from the outputleftXY to 
            offsets from outputReightXY where the offsets are
            the same value as "tolerance" and applied in 45 degree
            directions to get roughly near perpendicular to the epipolar
            line.
        */
        
        log.fine("nIter=" + nIter);

        log.fine("final fit: " + consensusFit.toString());

        return consensusFit;
    }
}
