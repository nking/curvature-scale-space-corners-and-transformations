package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.imageProcessing.util.RANSACAlgorithmIterations;
import algorithms.misc.MiscMath;
import algorithms.util.PairFloatArray;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
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
    
    private double generalTolerance = 3;
    
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
        
        //TODO:  still not quite right...
        
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
        */
        
        int nSet = 7;
        
        int nPoints = matchedLeftXY.numCols();
                       
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(System.currentTimeMillis());
        
        int[] selected = new int[nSet];
        
        PairFloatArray xy1 = null;
        PairFloatArray xy2 = null;
        
        double tolerance = generalTolerance;
        //double freqFractionToAddBackIn = 0.3;
        
        // ====== find matched pairs incremently, and best among only 7 point match ====
        
        StereoProjectionTransformerFit bestFit = null;
        List<Integer> inlierIndexes = null;
        
        RANSACAlgorithmIterations nEstimator = new RANSACAlgorithmIterations();
        
        int nIter = nEstimator.estimateNIterUsingStandardRANSACApproximation(
            0.5, nSet);
        
        int nMaxMatchable = (matchedLeftXY.numCols() < matchedRightXY.numCols())
            ? matchedLeftXY.numCols() : matchedRightXY.numCols();
        
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
            SimpleMatrix[] fms = 
                spTransformer.calculateEpipolarProjectionFor7Points(xy1, xy2);

            for (SimpleMatrix fm : fms) {
            
                // evaluate fit against all points
                StereoProjectionTransformerFit fit = 
                    spTransformer.evaluateFitForAlreadyMatched(
                    fm, matchedLeftXY, matchedRightXY, tolerance);

                List<Integer> indexes = fit.getInlierIndexes();
               
                if (fitIsBetter(bestFit, fit)) {
                    bestFit = fit;                
                    inlierIndexes = indexes;
                    
                    if (indexes.size() >= nMaxMatchable*0.8) {
                        break;
                    }
                }
            }
        }
        
        if (bestFit == null) {
            return null;
        }

        // store inliers in outputLeftXY and outputRightXY and redo the
        // entire fit using only the inliers to determine the fundamental
        // matrix.
                    
        for (Integer idx : inlierIndexes) {
            int idxInt = idx.intValue();
            outputLeftXY.add(
                (float)matchedLeftXY.get(0, idxInt), 
                (float)matchedLeftXY.get(1, idxInt));
            outputRightXY.add(
                (float)matchedRightXY.get(0, idxInt), 
                (float)matchedRightXY.get(1, idxInt));            
        }
        
        StereoProjectionTransformer spTransformer = new 
            StereoProjectionTransformer();
        
        SimpleMatrix finalFM = 
            spTransformer.calculateEpipolarProjectionForPerfectlyMatched(
            outputLeftXY, outputRightXY);
        
        StereoProjectionTransformerFit finalFit = 
            spTransformer.evaluateFitForAlreadyMatched(finalFM, 
            matchedLeftXY, matchedRightXY, tolerance);

        log.info("nIter=" + nIter);

        log.info("best fit from 7-point: " + bestFit.toString());

        log.info("final fit: " + finalFit.toString());

        return finalFit;
    }
    
    /**
     * calculate the epipolar projection among the given points with the
     * assumption that some of the points in the matched lists are not
     * true matches.  this has a frequency based approach if a critical
     * limit consensus isn't found in the best fit.
     * 
     * @param matchedLeftXY
     * @param matchedRightXY
     * @param outputLeftXY
     * @param outputRightXY
     * @return
     * @throws NoSuchAlgorithmException 
     */
    public StereoProjectionTransformerFit calculateEpipolarProjection2(
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
                       
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(System.currentTimeMillis());
        
        int[] selected = new int[nSet];
        
        PairFloatArray xy1 = null;
        PairFloatArray xy2 = null;
        
        double tolerance = generalTolerance;
        double limitFrac = 0.6;
        
        // ====== find matched pairs incremently, and best among only 7 point match ====
        
        
        StereoProjectionTransformerFit bestFit = null;
        List<Integer> inlierIndexes = null;
        
        // for each consensus member, increment the count for its index
        Map<Integer, Integer> frequencyOfIndexes = new HashMap<Integer, Integer>();
        
        
        RANSACAlgorithmIterations nEstimator = new RANSACAlgorithmIterations();
        
        int nIter = nEstimator.estimateNIterUsingStandardRANSACApproximation(
            0.5, nSet);
        
        int nMaxMatchable = (matchedLeftXY.numCols() < matchedRightXY.numCols())
            ? matchedLeftXY.numCols() : matchedRightXY.numCols();
        
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
            SimpleMatrix[] fms = 
                spTransformer.calculateEpipolarProjectionFor7Points(xy1, xy2);

            for (SimpleMatrix fm : fms) {
            
                // evaluate fit against all points
                StereoProjectionTransformerFit fit = 
                    spTransformer.evaluateFitForAlreadyMatched(
                    fm, matchedLeftXY, matchedRightXY, tolerance);

                List<Integer> indexes = fit.getInlierIndexes();
                if (indexes != null) {
                    for (Integer index : indexes) {
                        Integer freq = frequencyOfIndexes.get(index);
                        if (freq == null) {
                            frequencyOfIndexes.put(index, Integer.valueOf(1));
                        } else {
                            frequencyOfIndexes.put(index, Integer.valueOf(
                                freq.intValue() + 1));
                        }
                    }
                }

                if (fitIsBetter(bestFit, fit)) {
                    bestFit = fit;                
                    inlierIndexes = indexes;
                    
                    if (indexes.size() >= nMaxMatchable*0.8) {
                        break;
                    }
                }
            }
        }
        
        if (bestFit == null) {
            return null;
        }

        // store inliers in outputLeftXY and outputRightXY and redo the
        // entire fit using only the inliers to determine the fundamental
        // matrix.
                 
        int n = (int)Math.round(limitFrac*nMaxMatchable);
        
        if (false && (inlierIndexes.size() >= n) || (frequencyOfIndexes.size() < n)) {
            for (Integer idx : inlierIndexes) {
                int idxInt = idx.intValue();
                outputLeftXY.add(
                    (float)matchedLeftXY.get(0, idxInt), 
                    (float)matchedLeftXY.get(1, idxInt));
                outputRightXY.add(
                    (float)matchedRightXY.get(0, idxInt), 
                    (float)matchedRightXY.get(1, idxInt));

                frequencyOfIndexes.remove(idx);
            }
        } else {
            int[] a1 = new int[frequencyOfIndexes.size()];
            int[] a2 = new int[a1.length];
            
            // examine which indexes are still in frequency map
            Iterator<Entry<Integer, Integer> > iter = frequencyOfIndexes.entrySet().iterator();
            int c = 0;
            while (iter.hasNext()) {
                Entry<Integer, Integer> entry = iter.next();
                Integer index = entry.getKey();
                int idx = index.intValue();
                Integer freq = entry.getValue();
                a1[c] = idx;
                a2[c] = freq.intValue();
                c++;
               
            }
            
            if (frequencyOfIndexes.size() < n) {
                n = frequencyOfIndexes.size();
            }
            
            MultiArrayMergeSort.sortByDecr(a1, a2);
            for (int ii = 0; ii < n; ii++) {
                int idx = a1[ii];
                outputLeftXY.add(
                    (float) matchedLeftXY.get(0, idx),
                    (float) matchedLeftXY.get(1, idx));
                outputRightXY.add(
                    (float) matchedRightXY.get(0, idx),
                    (float) matchedRightXY.get(1, idx));
            }
        }
        
        StereoProjectionTransformer spTransformer = new 
            StereoProjectionTransformer();
        
        SimpleMatrix finalFM = 
            spTransformer.calculateEpipolarProjectionForPerfectlyMatched(
            outputLeftXY, outputRightXY);
        
        StereoProjectionTransformerFit finalFit = 
            spTransformer.evaluateFitForAlreadyMatched(finalFM, 
            matchedLeftXY, matchedRightXY, tolerance);

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
