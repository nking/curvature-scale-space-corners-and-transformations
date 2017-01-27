package algorithms.imageProcessing.features;

import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.EpipolarTransformer;
import algorithms.imageProcessing.util.RANSACAlgorithmIterations;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
        
        SimpleMatrix input1 =
            spTransformer.rewriteInto3ColumnMatrix(matchedLeftXY);

        SimpleMatrix input2 =
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
        SimpleMatrix matchedLeftXY, SimpleMatrix matchedRightXY,
        PairIntArray outputLeftXY, PairIntArray outputRightXY) {

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
        */

        int nSet = 7;

        int nPoints = matchedLeftXY.numCols();
        
        ErrorType errorType = ErrorType.SAMPSONS;

        EpipolarTransformer spTransformer = new EpipolarTransformer();
            
        SecureRandom sr = Misc.getSecureRandom();
        long seed = System.currentTimeMillis();
        log.fine("SEED=" + seed + " nPoints=" + nPoints);
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

        int nIter = 0;
        
        int[] selectedIndexes = new int[nSet];
        
        SimpleMatrix sampleLeft = new SimpleMatrix(3, nSet);
        SimpleMatrix sampleRight = new SimpleMatrix(3, nSet);
        
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
            List<SimpleMatrix> fms =
                spTransformer.calculateEpipolarProjectionFor7Points(sampleLeft, 
                    sampleRight);

            if (fms == null || fms.isEmpty()) {
                nIter++;
                continue;
            }

            // use point dist to epipolar lines to estimate errors of sample
            EpipolarTransformationFit fit = null;
            
            for (SimpleMatrix fm : fms) {
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
                    /(double)matchedLeftXY.numCols();
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
        SimpleMatrix inliersLeftXY = new SimpleMatrix(3, n);
        SimpleMatrix inliersRightXY = new SimpleMatrix(3, n);
        
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
        
        if (inliersRightXY.numCols() == 7) {
            
            List<SimpleMatrix> fms = spTransformer.calculateEpipolarProjectionFor7Points(
                inliersLeftXY, inliersRightXY);
            if (fms == null || fms.isEmpty()) {
                return null;
            }
            EpipolarTransformationFit fit = null;
            for (SimpleMatrix fm : fms) {
                EpipolarTransformationFit fitI = 
                    spTransformer.calculateError(fm, matchedLeftXY, 
                        matchedRightXY, errorType, tolerance);
                if (fitI.isBetter(fit)) {
                    fit = fitI;
                }
            }
            consensusFit = fit;
            
        } else {
            
            SimpleMatrix fm = 
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
        float[] toleranceErrors =
            calculateErrorsAsOffsets(consensusFit.getFundamentalMatrix(), 
                outputLeftXY, outputRightXY, tolerance);
        
        consensusFit.setTolerance(toleranceErrors[0],
            toleranceErrors[1]);
        
        log.fine("nIter=" + nIter);

        log.fine("final fit: " + consensusFit.toString());

        return consensusFit;
    }

    private float[] calculateErrorsAsOffsets(SimpleMatrix fm,
        PairIntArray xy1, PairIntArray xy2, int tolerance) {

        /*
        -- make copies of xy2, displacing the point from center by
        radius of tolerance and angles 45, 90, 135, 180, 225, 270, 315 
        and storing them in separate xy2Offset arrays.
        -- choose the max error for each offset of the same point
        -- calc the mean and standard deviation of the max offset
           errors
        this should be a rough estimate of the ability to associate a point
        with an epipolar line.
        */
                
        EpipolarTransformer eTransformer = new EpipolarTransformer();
        
        SimpleMatrix m1m = eTransformer.rewriteInto3ColumnMatrix(xy1);
        
        float[] xErrors = new float[xy1.getN()];
        float[] yErrors = new float[xy1.getN()];
        Arrays.fill(xErrors, Float.MIN_VALUE);
        Arrays.fill(yErrors, Float.MIN_VALUE);
        
        float[] outputDist = new float[2];
        
        for (int i = 0; i < 7; ++i) {
            
            double angle = (i + 1) * Math.PI/4.;
            
            int xOff = (int)Math.round(tolerance * Math.cos(angle));
            int yOff = (int)Math.round(tolerance * Math.sin(angle));
            
            PairIntArray xy2Offset = new PairIntArray(xy2.getN());
            for (int j = 0; j < xy2.getN(); ++j) {
                int x = xy2.getX(j) + xOff;
                int y = xy2.getY(j) + yOff;
                xy2Offset.add(x, y);
            }
         
            SimpleMatrix m2m = eTransformer.rewriteInto3ColumnMatrix(xy2Offset);
        
            SimpleMatrix m2EpipolarLines = fm.mult(m1m);
            SimpleMatrix m1EpipolarLines = fm.transpose().mult(m2m);

            for (int j = 0; j < xy2Offset.getN(); ++j) {
                eTransformer.calculatePerpDistFromLines(
                    m1m, m2m, 
                    m2EpipolarLines, m1EpipolarLines, 
                    j, j, outputDist);
            
                float d1 = Math.abs(outputDist[0]);
                float d2 = Math.abs(outputDist[1]);
                if (xErrors[j] < d1) {
                    xErrors[j] = d1;
                }
                if (yErrors[j] < d2) {
                    yErrors[j] = d2;
                }
            }
        }
        
        // calc mean and stdev of tolerance size offsets
        float[] xAvgStDv = MiscMath.getAvgAndStDev(xErrors);
        float[] yAvgStDv = MiscMath.getAvgAndStDev(yErrors);
    
        float combinedAvgTol = (xAvgStDv[0] + yAvgStDv[0])/2.f;
        float combinedStdvTol = (float)Math.sqrt(xAvgStDv[1]*xAvgStDv[1] 
            + yAvgStDv[1]*yAvgStDv[1]);
        
        return new float[]{combinedAvgTol, combinedStdvTol};
    }

}
