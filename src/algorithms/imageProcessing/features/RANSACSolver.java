package algorithms.imageProcessing.features;

import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.transform.Distances;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.EpipolarTransformer;
import algorithms.imageProcessing.transform.Util;
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
        PairIntArray outputLeftXY, PairIntArray outputRightXY, int tolerance) {

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
        
        DenseMatrix input1 =
            Util.rewriteInto3ColumnMatrix(matchedLeftXY);

        DenseMatrix input2 =
            Util.rewriteInto3ColumnMatrix(matchedRightXY);
        
        return calculateEpipolarProjection(input1, input2,
            outputLeftXY, outputRightXY, tolerance);
    }

    /**
     * calculate the epipolar projection among the given points with the
     * assumption that some of the points in the matched lists are not
     * true matches.
     *
     * @param matchedLeftXY 3 X nData matrix with rows being x, y, and 1's respectively
     * @param matchedRightXY 
     * @param outputLeftXY
     * @param outputRightXY
     * @param tolerance tolerance in distance from epipolar line for a point to 
     * be an inlier in the final fit.
     * @return
     */
    public EpipolarTransformationFit calculateEpipolarProjection(
        final DenseMatrix matchedLeftXY, final DenseMatrix matchedRightXY,
        final PairIntArray outputLeftXY, final PairIntArray outputRightXY,
        final int tolerance) {

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
        Using 7 point samples for epipolar transformation fits.
        -- the number of iterations for testing sub-samples of nPoints 
           (each of size 7) and finding one to be a good sub-sample with an
           excess of probability of 95% is estimated for a given
           percent of bad data.
           NOTE, the algorithm proceeds by assuming 50% bad data and improves
           that upon each best fitting sub-sample.
           NOTE, the number of iterations is then re-calculated too.
        -- once the best fitting sub-sample is found, that Fundamental
           Matrix is used with all of the points to determine inliers for that
           transformation.
        -- that consensus of inliers are then used to re-compute a final Fundamental
           Matrix.
        
        NOTE that the sub-samples are selected randomly from all possible
        sub-samples of the nPoints unless the number of all possible 
        sub-samples is smaller than the expected number of iterations for 95%
        probability of a good sub-sample.  In the later case, all sub-samples
        are tried.
        
        */

        final int nSet = 7;

        int nPoints = matchedLeftXY.numColumns();
        
        // n!/(k!*(n-k)!
        long nPointsSubsets = MiscMath.computeNDivKTimesNMinusKExact(nPoints, nSet);

        SecureRandom sr = Misc.getSecureRandom();
        long seed = System.currentTimeMillis();
        log.info("SEED=" + seed + " nPoints=" + nPoints);
        sr.setSeed(seed);

        ErrorType errorType = ErrorType.DIST_TO_EPIPOLAR_LINE;

        EpipolarTransformer spTransformer = new EpipolarTransformer();
        
        // consensus indexes
        EpipolarTransformationFit bestFit = null;
        
        /*
        could consider a threshold max iteration based upon the image size such
        as in (http://phototour.cs.washington.edu/ModelingTheWorld_ijcv07.pdf)
        which uses 0.6% of the maximum image dimension.
        */
        
        int outlierPercent = 50;
        
        long nMaxIter;
        if (nPoints == nSet) {
            nMaxIter = 1;
        } else {
            nMaxIter = RANSACAlgorithmIterations
                .numberOfSubsamplesOfSize7For95PercentInliers(outlierPercent);
        }
        
        System.out.println("nPoints=" + nPoints + " estimate for nMaxIter=" +
            nMaxIter + " (n!/(k!*(n-k)!)=" + nPointsSubsets);

        if (nMaxIter > nPointsSubsets) {
            nMaxIter = nPointsSubsets;
        }
        
        int nIter = 0;
        
        int[] selectedIndexes = new int[nSet];
        
        DenseMatrix sampleLeft = new DenseMatrix(3, nSet);
        DenseMatrix sampleRight = new DenseMatrix(3, nSet);
        
        DenseMatrix bestSampleLeft = null;
        DenseMatrix bestSampleRight = null;
        
        // initialize the unchanging 3rd dimension
        for (int i = 0; i < nSet; ++i) {
            sampleLeft.set(2, i, 1);
            sampleRight.set(2, i, 1);
        }
        
        Distances distances = new Distances();
        
        while (nIter < nMaxIter) {
            
            MiscMath.chooseRandomly(sr, selectedIndexes, nPoints);

            int count = 0;
            
            for (int bitIndex : selectedIndexes) {

                int idx = bitIndex;

                sampleLeft.set(0, count, matchedLeftXY.get(0, idx));
                sampleLeft.set(1, count, matchedLeftXY.get(1, idx));
                                
                sampleRight.set(0, count, matchedRightXY.get(0, idx));
                sampleRight.set(1, count, matchedRightXY.get(1, idx));
                
                count++;
            }

            // determine matrix from 7 points.
            List<DenseMatrix> fms = spTransformer.calculateEpipolarProjectionFor7Points(
                sampleLeft, sampleRight);

            if (fms == null || fms.isEmpty()) {
                nIter++;
                continue;
            }
            
            // use point dist to epipolar lines to estimate errors of sample
            EpipolarTransformationFit fit = null;
            
            for (DenseMatrix fm : fms) {
                
                // indexes within this refer to the sample data, not the larger matched data
                EpipolarTransformationFit fitI = distances.calculateError(fm, 
                    sampleLeft, sampleRight, errorType, tolerance);
                     
                if (fitI.getInlierIndexes().size() == nSet &&
                    fitI.isBetter(fit)) {
                    
                    System.out.println("new local best fit: " + fitI.toString());
                    System.out.flush();
                    
                    fit = fitI;
                }
            }
            
            if (fit == null) {
                nIter++;
                continue;
            }
            
            if (fit.isBetter(bestFit)) {
                bestFit = fit;
                bestSampleLeft = sampleLeft.copy();
                bestSampleRight = sampleRight.copy();
                
                System.out.println("**best fit: " + bestFit.toString());
                System.out.flush();
                
                // recalculate nMaxIter
                if (nMaxIter < nPointsSubsets && nMaxIter > 1) {
                    double outlierPercentI = 100.*
                        (double)(nPoints - bestFit.getInlierIndexes().size()) /
                        (double)nPoints;
                    if (outlierPercentI < outlierPercent) {
                        outlierPercent = (int)Math.ceil(outlierPercentI);
                        if (outlierPercent < 5) {
                            outlierPercent = 5;
                        }
                        assert(outlierPercent < 50);
                        nMaxIter = RANSACAlgorithmIterations
                            .numberOfSubsamplesOfSize7For95PercentInliers(outlierPercent);
                    }
                }
            }
            
            nIter++;
        }

        if (bestFit == null) {
            log.info("no solution.  nIter=" + nIter);
            return null;
        }
                
        // -- use best fit epipolar projection on matched points
        // -- remove outliers
        // -- re-calculate epipolar projection from remaining
        // -- return finished fit and errors
        
        EpipolarTransformationFit bestFitAppliedToAll 
            = distances.calculateError(bestFit.getFundamentalMatrix(), 
            matchedLeftXY, matchedRightXY, errorType, tolerance);
                       
        log.fine("bestFitAppliedToAll: " + bestFitAppliedToAll.toString());
        
        if (bestFitAppliedToAll.getInlierIndexes().size() < bestFit.getInlierIndexes().size()) {
            System.out.println("bestFit: " + bestFit.toString());
            System.out.println("bestFitAppliedToAll: " + bestFitAppliedToAll.toString());
            System.out.flush();
            
            // re-do to see diffs
            EpipolarTransformationFit temp1 
                = distances.calculateError(bestFit.getFundamentalMatrix(), 
                bestSampleLeft, bestSampleRight, errorType, tolerance);
            
            EpipolarTransformationFit temp2 
            = distances.calculateError(bestFit.getFundamentalMatrix(), 
            matchedLeftXY, matchedRightXY, errorType, tolerance);
            
            throw new IllegalStateException("error!");
        }
        
        // create an inlier list for a new transformation 
        int nInliers = bestFitAppliedToAll.getInlierIndexes().size();        
        DenseMatrix inliersLeftXY = new DenseMatrix(3, nInliers);
        DenseMatrix inliersRightXY = new DenseMatrix(3, nInliers);
        int count = 0;
        for (Integer idx : bestFitAppliedToAll.getInlierIndexes()) {
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
        
        boolean usePreviousFM = false;
        
        if (nInliers == nSet) {
            List<DenseMatrix> fms = spTransformer.calculateEpipolarProjectionFor7Points(
                inliersLeftXY, inliersRightXY);
            if (fms == null || fms.isEmpty()) {
                usePreviousFM = true;
            } else {
                for (DenseMatrix fm : fms) {
                    EpipolarTransformationFit fitI
                        = distances.calculateError(fm, 
                        matchedLeftXY, matchedRightXY, 
                        errorType, tolerance);
                    if (fitI.isBetter(consensusFit)) {
                        consensusFit = fitI;
                    }
                }
            }
        } else {            
            DenseMatrix fm = spTransformer.calculateEpipolarProjection(
                inliersLeftXY, inliersRightXY); 
            if (fm == null) {
                usePreviousFM = true;
            } else {
                consensusFit = distances.calculateError(fm, matchedLeftXY, 
                    matchedRightXY, errorType, tolerance);  
            }
        }
        
        if (consensusFit != null) {
            log.fine("consensusFit=" + consensusFit.toString());
        }
        
        if (consensusFit == null || usePreviousFM || !consensusFit.isBetter(bestFitAppliedToAll)) {
            log.fine("choosing the best fit applied to all instead of"
                + " the subsequent consensus fit");
            consensusFit = bestFitAppliedToAll;
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
      
        log.fine("nIter=" + nIter);

        log.fine("final fit: " + consensusFit.toString());

        return consensusFit;
    }
}
