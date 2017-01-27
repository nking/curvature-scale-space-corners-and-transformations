package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.transform.EpipolarFeatureTransformationFit;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.EpipolarTransformer;
import algorithms.imageProcessing.util.RANSACAlgorithmIterations;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.security.SecureRandom;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import org.ejml.simple.SimpleMatrix;

/**
 * given matched point lists, determine the best epipolar solution using a
 * 7-point epipolar calculation and random draws of 7 points from the
 * matched point lists under the assumption that some of the matched points
 * are not true (correct) matches. The evaluation of the solution uses
 * features and either epipolar distances or Sampson's error distances.
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
public class RANSACEpipolarWithFeaturesSolver {

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
     * @param clrFeatures1
     * @param clrFeatures2
     * @param kpab1
     * @param bmaIndex1
     * @param kpab2
     * @param bmaIndex2
     * @param rImg1
     * @param gImg1
     * @param bImg1
     * @param rImg2
     * @param gImg2
     * @param bImg2
     * @param outputLeftXY
     * @param outputRightXY
     * @param useHalfDescriptors
     * @return
     */
    public EpipolarFeatureTransformationFit calculateEpipolarProjection(
        PairIntArray matchedLeftXY, PairIntArray matchedRightXY,
        IntensityClrFeatures clrFeatures1, IntensityClrFeatures clrFeatures2,
        KeyPointsAndBounds kpab1, int bmaIndex1, KeyPointsAndBounds kpab2, int bmaIndex2,
        GreyscaleImage rImg1, GreyscaleImage gImg1, GreyscaleImage bImg1, 
        GreyscaleImage rImg2, GreyscaleImage gImg2, GreyscaleImage bImg2, 
        PairIntArray outputLeftXY, PairIntArray outputRightXY, 
        boolean useHalfDescriptors) {

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
            clrFeatures1, clrFeatures2, 
            kpab1, bmaIndex1, kpab2, bmaIndex2,
            rImg1, gImg1, bImg1, rImg2, gImg2, bImg2,
            outputLeftXY, outputRightXY, useHalfDescriptors);
    }

    /**
     * calculate the epipolar projection among the given points with the
     * assumption that some of the points in the matched lists are not
     * true matches.
     *
     * @param matchedLeftXY
     * @param matchedRightXY
     * @param clrFeatures1
     * @param clrFeatures2
     * @param kpab1
     * @param bmaIndex1
     * @param kpab2
     * @param bmaIndex2
     * @param rImg1
     * @param gImg1
     * @param bImg1
     * @param rImg2
     * @param gImg2
     * @param bImg2
     * @param outputLeftXY
     * @param outputRightXY
     * @param useHalfDescriptors
     * @return
     */
    public EpipolarFeatureTransformationFit calculateEpipolarProjection(
        SimpleMatrix matchedLeftXY, SimpleMatrix matchedRightXY,
        IntensityClrFeatures clrFeatures1, IntensityClrFeatures clrFeatures2,
        KeyPointsAndBounds kpab1, int bmaIndex1, KeyPointsAndBounds kpab2, int bmaIndex2,
        GreyscaleImage rImg1, GreyscaleImage gImg1, GreyscaleImage bImg1, 
        GreyscaleImage rImg2, GreyscaleImage gImg2, GreyscaleImage bImg2, 
        PairIntArray outputLeftXY, PairIntArray outputRightXY,
        boolean useHalfDescriptors) {

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

        int tolerance = 3;

        // consensus indexes
        EpipolarFeatureTransformationFit bestFit = null;
        
        RANSACAlgorithmIterations nEstimator = new RANSACAlgorithmIterations();

        /*
        could consider a threshold max iteration based upon the image size such
        as in (http://phototour.cs.washington.edu/ModelingTheWorld_ijcv07.pdf)
        which uses 0.6% of the maximum image dimension.
        */
        long nMaxIter = nEstimator.estimateNIterFor99PercentConfidence(nPoints, 
            7, 0.25);

        if (nPoints == 7) {
            nMaxIter = 1;
        }

        int nIter = 0;
        
        int[] selectedIndexes = new int[nSet];
        
        SimpleMatrix sampleLeft = new SimpleMatrix(3, nSet);
        SimpleMatrix sampleRight = new SimpleMatrix(3, nSet);
        
        while ((nIter < nMaxIter) && (nIter < 10000)) {
            
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
            EpipolarFeatureTransformationFit fit = null;
            
            for (SimpleMatrix fm : fms) {
                EpipolarFeatureTransformationFit fitI = 
                    spTransformer.calculateError(fm, matchedLeftXY, 
                        matchedRightXY, clrFeatures1, clrFeatures2,
                        kpab1, bmaIndex1, kpab2, bmaIndex2,
                        rImg1, gImg1, bImg1, rImg2, gImg2, bImg2,
                        errorType, tolerance, useHalfDescriptors);
                
                if (fitI.isBetter(fit)) {
                    fit = fitI;
                }
            }
            
            if (fit == null || fit.getInlierIndexes().isEmpty()) {
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

        EpipolarFeatureTransformationFit consensusFit = null;
        
        if (inliersRightXY.numCols() == 7) {
            
            List<SimpleMatrix> fms = spTransformer.calculateEpipolarProjectionFor7Points(
                inliersLeftXY, inliersRightXY);
            if (fms == null || fms.isEmpty()) {
                return null;
            }
            EpipolarFeatureTransformationFit fit = null;
            for (SimpleMatrix fm : fms) {
                EpipolarFeatureTransformationFit fitI = 
                    spTransformer.calculateError(fm, matchedLeftXY, 
                        matchedRightXY, clrFeatures1, clrFeatures2,
                        kpab1, bmaIndex1, kpab2, bmaIndex2,
                        rImg1, gImg1, bImg1, rImg2, gImg2, bImg2,
                        errorType, tolerance, useHalfDescriptors);
                
                if (fitI.isBetterByCost(fit)) {
                    fit = fitI;
                }
            }
            consensusFit = fit;
            
        } else if (inliersRightXY.numCols() > 7) {
            
            SimpleMatrix fm = 
                spTransformer.calculateEpipolarProjection(
                inliersLeftXY, inliersRightXY);
            
            EpipolarFeatureTransformationFit fit = 
                spTransformer.calculateError(fm, matchedLeftXY, 
                    matchedRightXY, clrFeatures1, clrFeatures2,
                    kpab1, bmaIndex1, kpab2, bmaIndex2,
                    rImg1, gImg1, bImg1, rImg2, gImg2, bImg2,
                    errorType, tolerance, useHalfDescriptors);
            
            consensusFit = fit;
            
        } else {
           
            if (bestFit.getInlierIndexes().size() > 2) {
                consensusFit = bestFit;
            } else {
                return null;
            }
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
        
        log.info("nIter=" + nIter);

        log.info("final fit: " + consensusFit.toString());

        return consensusFit;
    }

}
