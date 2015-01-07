package algorithms.imageProcessing;

import algorithms.imageProcessing.util.RANSACAlgorithmIterations;
import algorithms.misc.MiscMath;
import algorithms.util.PairFloatArray;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author nichole
 */
public class RANSACSolver {
    
    private boolean debug = true;
    
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
        
        int nSet = 7;
        
        int nPoints = matchedLeftXY.numCols();
        
        RANSACAlgorithmIterations nEstimator = new RANSACAlgorithmIterations();
        
        int nIter = nEstimator.estimateNIterUsingStandardRANSACApproximation(
            0.5, nSet);
        
        //TODO: looks like the evaluation needs to be improved.
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(System.currentTimeMillis());
        
        int[] selected = new int[nSet];
        
        PairFloatArray xy1 = null;
        PairFloatArray xy2 = null;
        
        //TODO: determine this using one of the several suggested ways.
        //TODO: also consider ability to change it +- 1 or so to find best?
        
        double threshold = 5;//4
        
        // ====== find matched pairs incremently, and best among only 7 point match ====
        
        StereoProjectionTransformerFit bestFit = null;
        List<Integer> inlierIndexes = null;
        SimpleMatrix bestFundamentalMatrix = null;
                
        for (int i = 0; i < nIter; i++) {
            
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
            SimpleMatrix fm = 
                spTransformer.calculateEpipolarProjectionFor7Points(xy1, xy2);
         
            StereoProjectionTransformerFit sevenPointFit = 
                spTransformer.evaluateFitForAlreadyMatched(fm, 
                StereoProjectionTransformer.rewriteInto3ColumnMatrix(xy1),
                StereoProjectionTransformer.rewriteInto3ColumnMatrix(xy2),
                threshold);
            
            if ((sevenPointFit == null) || 
                (sevenPointFit.getInlierIndexes() == null) ||
                (sevenPointFit.getInlierIndexes().size() != 7)) {
                // TODO: when does this happen?
                continue;
            }
            
            // evaluate fit against all points
            StereoProjectionTransformerFit fit = 
                spTransformer.evaluateFitForAlreadyMatched(
                fm, matchedLeftXY, matchedRightXY, threshold);
            
            if (fitIsBetter(bestFit, fit)) {
                
                bestFit = fit;
                
                inlierIndexes = fit.getInlierIndexes();
                
                bestFundamentalMatrix = fm;
            }
        }
        
        if (bestFit == null) {
            return null;
        }

        // store inliers in outputLeftXY and outputRightXY and redo the
        // entire fit using only the inliers to determine the fundamental
        // matrix.
        
        StereoProjectionTransformer spTransformer = new 
            StereoProjectionTransformer();
        
        for (Integer idx : inlierIndexes) {
            int idxInt = idx.intValue();
            outputLeftXY.add(
                (float)matchedLeftXY.get(0, idxInt), 
                (float)matchedLeftXY.get(1, idxInt));
            outputRightXY.add(
                (float)matchedRightXY.get(0, idxInt), 
                (float)matchedRightXY.get(1, idxInt));
        }
        
        SimpleMatrix finalFM = 
            spTransformer.calculateEpipolarProjectionForPerfectlyMatched(
            outputLeftXY, outputRightXY);
        
        StereoProjectionTransformerFit finalFit = 
            spTransformer.evaluateFitForAlreadyMatched(finalFM, 
            matchedLeftXY, matchedRightXY, threshold);

        Logger log = Logger.getLogger(this.getClass().getName());

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
}
