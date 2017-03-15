package algorithms.imageProcessing.features;

import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.EpipolarTransformer;
import algorithms.imageProcessing.util.RANSACAlgorithmIterations;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;
import no.uib.cipr.matrix.DenseMatrix;

/**
 * given matched point lists, determine the best epipolar solution using a
 * 7-point epipolar calculation and random draws of 7 points from the
 * matched point lists under the assumption that some of the matched points
 * are not true (correct) matches.
 * <pre>
 * useful reading:
 * http://6.869.csail.mit.edu/fa12/lectures/lecture13ransac/lecture13ransac.pdf
 * and
 * http://www.dtic.mil/dtic/tr/fulltext/u2/a460585.pdf
 * 
 * This algorithm differs from RANSACSolver.java in that it accepts more than
 * one match for point1 to point2 where point1 is a ppint in image1 and point2
 * is a point in image2.
 * 
 * It follows "Generalized RANSAC framework for relaxed correspondence problems"
   by Zhang and Kosecka in choosing randomly from a single point in image1 but 
   then choosing randomly from that point's possibly multiple matches.
   </pre>
 * @author nichole
 */
public class RANSACMultiplicitySolver {

    private boolean debug = true;

    private Logger log = Logger.getLogger(this.getClass().getName());

    /**
     * NOT READY FOR USE YET
     * calculate the epipolar projection among the given points with the
     * assumption that some of the points in the matched lists are not
     * true matches.
     *
     * @param matchedLeftXY
     * @param matchedRightXYs
     * @param outputLeftXY
     * @param outputRightXY
     * @return
     */
    public EpipolarTransformationFit calculateEpipolarProjection(
        List<PairInt> matchedLeftXY, List<List<PairInt>> matchedRightXYs,
        PairIntArray outputLeftXY, PairIntArray outputRightXY) {

        if (matchedLeftXY == null) {
            throw new IllegalArgumentException("matchedLeftXY cannot be null");
        }
        if (matchedRightXYs == null) {
            throw new IllegalArgumentException("matchedRightXYs cannot be null");
        }
        if (matchedLeftXY.size() != matchedRightXYs.size()) {
            throw new IllegalArgumentException(
            "matchedLeftXY and matchedRightXYs must be same size");
        }
        
        if (matchedLeftXY.size() < 7) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
            "the algorithms require 7 or more points.  matchedLeftXY.n=" 
            + matchedLeftXY.size());
        }
        
        /*
        -- randomly sample 7 points from matchedLeftXY then randomly sample
           from the left's right match if more than one match is present.
        -- calculate the epipolar fundamental matrix from the 7
        -- evaluate the sample:
           using Sampson's error:
              if error < tolerance, stores the entire sample in consensus
        -- check whether to exit loop:
           recalculate expected iterations and compare to nIterations or max
           nIterations
        -- after exit from the loop, the largest consensus is used to
           re-calculate the fundamental matrix as the last result.
           note that if there aren't enough points in a consensus to
           calculate the fundamental matrix, the result is null.
              Should presumably then filter out points whose error determined by 
              distance to epipolar lines is larger than tolerance of 3 pixels or so.      
        */

        int nSet = 7;
        
        int nPoints = matchedLeftXY.size();
        
        ErrorType errorType = ErrorType.SAMPSONS;
                
        SecureRandom sr = Misc.getSecureRandom();
        long seed = System.currentTimeMillis();
        log.fine("SEED=" + seed + " nPoints=" + nPoints);
        sr.setSeed(seed);
        
        int nAllMultiplicity = 0;
        for (int i = 0; i < matchedLeftXY.size(); ++i) {
            nAllMultiplicity += matchedRightXYs.get(i).size();
        }
        
        EpipolarTransformer spTransformer = new EpipolarTransformer();
        
        PairIntArray originalLeftXY = new PairIntArray(nAllMultiplicity);
        PairIntArray originalRightXY = new PairIntArray(nAllMultiplicity);
        
        DenseMatrix evalAllLeft = new DenseMatrix(3, nAllMultiplicity);
        DenseMatrix evalAllRight = new DenseMatrix(3, nAllMultiplicity);
        int count = 0;
        for (int i = 0; i < matchedLeftXY.size(); ++i) {
            PairInt lft = matchedLeftXY.get(i);
            for (PairInt rgt : matchedRightXYs.get(i)) {
                evalAllLeft.set(0, count, lft.getX());
                evalAllLeft.set(1, count, lft.getY());
                evalAllLeft.set(2, count, 1);
                evalAllRight.set(0, count, rgt.getX());
                evalAllRight.set(1, count, rgt.getY());
                evalAllRight.set(2, count, 2);
                originalLeftXY.add(lft.getX(), lft.getY());
                originalRightXY.add(rgt.getX(), rgt.getY());
                count++;
            }
        }
            
        RANSACAlgorithmIterations nEstimator = new RANSACAlgorithmIterations();
        
        long nMaxIter = nEstimator.estimateNIterFor99PercentConfidenceDegenerate(
            nPoints, nAllMultiplicity, nSet, 0.5);
        
        long nMaxIter0 = 100000;
      
        if (nPoints == nSet) {
            nMaxIter = 1;
        }

        int nIter = 0;
        
        EpipolarTransformationFit bestFit = null;
                    
        int[] selectedIndexes = new int[nSet];
        
        DenseMatrix sampleLeft = new DenseMatrix(3, nSet);
        DenseMatrix sampleRight = new DenseMatrix(3, nSet);
        
        log.info("nPoints=" + nPoints + " n including multiplicity=" + nAllMultiplicity);
        
        int tolerance = 3;
        
        while ((nIter < nMaxIter) && (nIter < nMaxIter0)) {
              
            MiscMath.chooseRandomly(sr, selectedIndexes, nPoints);
            
            count = 0;
            
            for (int bitIndex : selectedIndexes) {

                int idx = bitIndex;

                sampleLeft.set(0, count, matchedLeftXY.get(idx).getX());
                sampleLeft.set(1, count, matchedLeftXY.get(idx).getY());
                sampleLeft.set(2, count, 1);
                
                // handle multiplicity
                PairInt rightP;
                int nr = matchedRightXYs.get(idx).size();
                assert(nr != 0);
                if (nr > 1) {
                    int idx2 = sr.nextInt(nr);
                    rightP = matchedRightXYs.get(idx).get(idx2);
                } else {
                    rightP = matchedRightXYs.get(idx).get(0);
                }
                sampleRight.set(0, count, rightP.getX());
                sampleRight.set(1, count, rightP.getY());
                sampleRight.set(2, count, 1);
                
                count++;
            }
           
            // determine matrix from 7 points.
            List<DenseMatrix> fms = 
                spTransformer.calculateEpipolarProjectionFor7Points(
                sampleLeft, sampleRight);
            
            if (fms == null || fms.isEmpty()) {
                nIter++;
                continue;
            }
                        
            // use Sampson's to estimate errors of sample
            EpipolarTransformationFit fit = null;
            
            for (DenseMatrix fm : fms) {
                EpipolarTransformationFit fitI = 
                    spTransformer.calculateErrorThenFilter(fm, 
                        evalAllLeft, evalAllRight, errorType, tolerance);
                
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
            if ((bestFit != null) && ((nIter % 50) == 0)) {
                double ratio = (double)bestFit.getInlierIndexes().size()/(double)nPoints;
                nMaxIter = nEstimator.estimateNIterFor99PercentConfidenceDegenerate(
                    nPoints, nAllMultiplicity, nSet, ratio);
            }
        }
        
        if (bestFit == null) {
            return null;
        }
        
        // filter for degeneracy cLIndexes and cRMultiplicityIndexes using cErrors
        Set<Integer> matchedLRM = filterForDegeneracy(bestFit,
            evalAllLeft, evalAllRight);
        
        if (matchedLRM.size() < 7) {
            return null;
        }
                
        // calculate fundamental matrix using filtered consensus
        DenseMatrix inliersLeftXY = new DenseMatrix(3, matchedLRM.size());
        DenseMatrix inliersRightXY = new DenseMatrix(3, matchedLRM.size());
        
        count = 0;
        for (Integer index : matchedLRM) {
            int idx = index.intValue();
            
            inliersLeftXY.set(0, count, evalAllLeft.get(0, idx));
            inliersLeftXY.set(1, count, evalAllLeft.get(1, idx));
            inliersLeftXY.set(2, count, 1);
            
            inliersRightXY.set(0, count, evalAllRight.get(0, idx));
            inliersRightXY.set(1, count, evalAllRight.get(1, idx));
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
                    spTransformer.calculateErrorThenFilter(fm,
                        inliersLeftXY, inliersRightXY, errorType, tolerance);
                if (fitI.isBetter(bestFit)) {
                    fit = fitI;
                }
            }
            consensusFit = fit;
            
        } else {
            
            DenseMatrix fm = spTransformer.calculateEpipolarProjection(
                inliersLeftXY, inliersRightXY);
            
            EpipolarTransformationFit fit = 
                spTransformer.calculateErrorThenFilter(fm,
                    inliersLeftXY, inliersRightXY, errorType, tolerance);
            
            consensusFit = fit;
        }
        
        // these are not normalized
        for (Integer index : consensusFit.getInlierIndexes()) {
            int idx = index.intValue();
            outputLeftXY.add(
                (int)Math.round(inliersLeftXY.get(0, idx)),
                (int)Math.round(inliersLeftXY.get(1, idx)));
            outputRightXY.add(
                (int)Math.round(inliersRightXY.get(0, idx)),
                (int)Math.round(inliersRightXY.get(1, idx)));
        }
        
        log.fine("nIter=" + nIter);

        log.fine("final fit: " + consensusFit.toString());

        return consensusFit;
    }

    private Set<Integer> filterForDegeneracy(EpipolarTransformationFit fit, 
        DenseMatrix allLeft, DenseMatrix allRight) {
        
        Set<Integer> matchedLRM = new HashSet<Integer>();
        
        if (fit == null || fit.getInlierIndexes().isEmpty()) {
            return matchedLRM;
        }
        
        //temporary look at all matchings
        if (true) {
            matchedLRM.addAll(fit.getInlierIndexes());
            return matchedLRM;
        }
        
        // key=left point coordinates, value=right indexes
        Map<PairInt, List<Integer>> pointIndexes = new HashMap<PairInt, List<Integer>>();
        Map<PairInt, List<Double>> pointErrors = new HashMap<PairInt, List<Double>>();
        
        for (int i = 0; i < fit.getInlierIndexes().size(); ++i) {
            
            int idx = fit.getInlierIndexes().get(i);
            
            int xL = (int)Math.round(allLeft.get(0, idx));
            int yL = (int)Math.round(allLeft.get(1, idx));
            PairInt key = new PairInt(xL, yL);
            
            Integer key2 = Integer.valueOf(idx);
            
            List<Integer> indexes2 = pointIndexes.get(key);
            List<Double> errors = pointErrors.get(key);
            
            if (indexes2 == null) {
                assert(errors == null);
                indexes2 = new ArrayList<Integer>();
                pointIndexes.put(key, indexes2);
                
                errors = new ArrayList<Double>();
                pointErrors.put(key, errors);
            } else {
                assert(errors != null);
            }
            
            indexes2.add(key2);
            errors.add(fit.getErrors().get(i));
        }
                
        // if values in pointIndexes are > 1 in size, keep one w/ smallest error
        for (Entry<PairInt, List<Integer>> entry : pointIndexes.entrySet()) {
                        
            List<Integer> indexes = entry.getValue();
            if (indexes.size() < 2) {
                matchedLRM.add(indexes.get(0));
                continue;
            }
            
            PairInt lftPt = entry.getKey();
            
            double minError = Double.MAX_VALUE;
            Integer minErrorIndex = null;
            for (int i = 0; i < indexes.size(); ++i) {
                double error = pointErrors.get(lftPt).get(i).doubleValue();
                if (error < minError) {
                    minError = error;
                    minErrorIndex = indexes.get(i);;
                }
            }
            
            assert(minErrorIndex != null);
            
            matchedLRM.add(minErrorIndex); 
        }
        
        // filter for degeneracy for multiple matches to right set
        // key=right point coordinates, value=right indexes
        pointIndexes = new HashMap<PairInt, List<Integer>>();
        pointErrors = new HashMap<PairInt, List<Double>>();
        
        for (int i = 0; i < fit.getInlierIndexes().size(); ++i) {
            
            int idx = fit.getInlierIndexes().get(i);
            
            Integer key2 = Integer.valueOf(idx);
            
            if (!matchedLRM.contains(key2)) {
                continue;
            }
            
            int xR = (int)Math.round(allRight.get(0, idx));
            int yR = (int)Math.round(allRight.get(1, idx));
            PairInt key = new PairInt(xR, yR);
                        
            List<Integer> indexes2 = pointIndexes.get(key);
            List<Double> errors = pointErrors.get(key);
            
            if (indexes2 == null) {
                assert(errors == null);
                indexes2 = new ArrayList<Integer>();
                pointIndexes.put(key, indexes2);
                
                errors = new ArrayList<Double>();
                pointErrors.put(key, errors);
            } else {
                assert(errors != null);
            }
            
            indexes2.add(key2);
            errors.add(fit.getErrors().get(i));
        }
        
        matchedLRM.clear();
        
        for (Entry<PairInt, List<Integer>> entry : pointIndexes.entrySet()) {
                        
            List<Integer> indexes = entry.getValue();
            if (indexes.size() < 2) {
                matchedLRM.add(indexes.get(0));
                continue;
            }
            
            PairInt lftPt = entry.getKey();
            
            double minError = Double.MAX_VALUE;
            Integer minErrorIndex = null;
            for (int i = 0; i < indexes.size(); ++i) {
                double error = pointErrors.get(lftPt).get(i).doubleValue();
                if (error < minError) {
                    minError = error;
                    minErrorIndex = indexes.get(i);;
                }
            }
            
            assert(minErrorIndex != null);
            
            matchedLRM.add(minErrorIndex); 
        }
        
        return matchedLRM;
    }
}
