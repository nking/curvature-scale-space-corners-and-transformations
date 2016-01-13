package algorithms.imageProcessing;

import algorithms.SubsetChooser;
import algorithms.imageProcessing.util.RANSACAlgorithmIterations;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
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
 * It follows "Generalized RANSAC framework for relaxed correspondence problems"
   by Zhang and Kosecka in choosing randomly from a single point in image1 but 
   then choosing randomly from that point's possibly multiple matches.

 * @author nichole
 */
public class RANSACMultiplicitySolver {

    private boolean debug = true;

    private double generalTolerance = 3;//5

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
     * @throws NoSuchAlgorithmException
     */
    public EpipolarTransformationFit calculateEpipolarProjection(
        List<PairInt> matchedLeftXY, List<List<PairInt>> matchedRightXYs,
        PairFloatArray outputLeftXY, PairFloatArray outputRightXY)
        throws NoSuchAlgorithmException {

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
        
        int tolerance = 3;
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        log.info("SEED=" + seed + " nPoints=" + nPoints);
        sr.setSeed(seed);
        
        SimpleMatrix evalAllLeft = new SimpleMatrix();
        SimpleMatrix evalAllRight = new SimpleMatrix();
        int nAllMultiplicity = 0;
        for (int i = 0; i < matchedLeftXY.size(); ++i) {
            PairInt lft = matchedLeftXY.get(i);
            for (PairInt rgt : matchedRightXYs.get(i)) {
                evalAllLeft.set(0, nAllMultiplicity, lft.getX());
                evalAllLeft.set(1, nAllMultiplicity, lft.getY());
                evalAllRight.set(0, nAllMultiplicity, rgt.getX());
                evalAllRight.set(1, nAllMultiplicity, rgt.getY());
                nAllMultiplicity++;
            }
        }
        
        RANSACAlgorithmIterations nEstimator = new RANSACAlgorithmIterations();

        int nMaxIter = nEstimator.estimateNIterFor99PercentConfidence(nPoints, 
            nSet, 0.5);
      
        if (nPoints == nSet) {
            nMaxIter = 1;
        }

        int nIter = 0;
        
        EpipolarTransformationFit bestFit = null;
        
        SubsetChooser subsetChooser = new SubsetChooser(nPoints, nSet);
            
        int[] selectedIndexes = new int[nSet];
        
        SimpleMatrix sampleLeft = new SimpleMatrix(3, nSet);
        SimpleMatrix sampleRight = new SimpleMatrix(3, nSet);

        int nV = subsetChooser.getNextSubset(selectedIndexes);

        int[] selectedRightIndexes2 = new int[nSet];
        
        while ((nV != -1) && (nIter < nMaxIter) && (nIter < 2000)) {
            
            int count = 0;
            
            for (int bitIndex : selectedIndexes) {

                int idx = bitIndex;

                sampleLeft.set(0, count, matchedLeftXY.get(idx).getX());
                sampleLeft.set(1, count, matchedLeftXY.get(idx).getY());
                sampleLeft.set(2, count, 1);
                
                // handle multiplicity
                PairInt rightP;
                if (matchedRightXYs.get(idx).size() > 1) {
                    int idx2 = sr.nextInt(matchedRightXYs.get(idx).size());
                    rightP = matchedRightXYs.get(idx).get(idx2);
                    selectedRightIndexes2[count] = idx2;
                } else {
                    rightP = matchedRightXYs.get(idx).get(0);
                    selectedRightIndexes2[count] = 0;
                }
                sampleRight.set(0, count, rightP.getX());
                sampleRight.set(1, count, rightP.getY());
                sampleRight.set(2, count, 1);
                
                count++;
            }
            
            StereoProjectionTransformer spTransformer = new StereoProjectionTransformer();

            // determine matrix from 7 points.
            List<SimpleMatrix> fms = spTransformer.calculateEpipolarProjectionFor7Points(
                sampleLeft, sampleRight);
            
            if (fms == null || fms.isEmpty()) {
                nIter++;
                continue;
            }
                        
            // use Sampson's to estimate errors of sample
            EpipolarTransformationFit fit = null;
            
            for (SimpleMatrix fm : fms) {
                EpipolarTransformationFit fitI = 
                    StereoProjectionTransformer.calculateSampsonsError(fm, 
                        evalAllLeft, evalAllRight, tolerance);
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
            
            nV = subsetChooser.getNextSubset(selectedIndexes);
            
            nIter++;
            
            // recalculate nMaxIter
            if ((bestFit != null) && ((nIter % 10) == 0)) {
                double ratio = (double)bestFit.getInlierIndexes().size()/(double)nAllMultiplicity;
                nMaxIter = nEstimator.estimateNIterFor99PercentConfidence(nPoints, 
                    nSet, ratio);
                int z = 1 ;
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
        SimpleMatrix inliersLeftXY = new SimpleMatrix(3, matchedLRM.size());
        SimpleMatrix inliersRightXY = new SimpleMatrix(3, matchedLRM.size());
        int count = 0;
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
        
        StereoProjectionTransformer spTransformer = new
            StereoProjectionTransformer();
        
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
                    StereoProjectionTransformer.calculateSampsonsError(fm, 
                        inliersLeftXY, inliersRightXY, tolerance);
                if (fitI.isBetter(bestFit)) {
                    fit = fitI;
                }
            }
            consensusFit = fit;
            
        } else {
            
            SimpleMatrix fm = 
                spTransformer.calculateEpipolarProjectionForPerfectlyMatched(
                inliersLeftXY, inliersRightXY);
            
            EpipolarTransformationFit fit = 
                StereoProjectionTransformer.calculateSampsonsError(fm, 
                inliersLeftXY, inliersRightXY, tolerance);
            
            consensusFit = fit;
        }
        
        for (Integer index : consensusFit.getInlierIndexes()) {
            int idx = index.intValue();
            outputLeftXY.add(
                (float)inliersLeftXY.get(0, idx),
                (float)inliersLeftXY.get(1, idx));
            outputRightXY.add(
                (float)inliersRightXY.get(0, idx),
                (float)inliersRightXY.get(1, idx));
        }
        
        log.info("nIter=" + nIter);

        log.info("final fit: " + consensusFit.toString());

        return consensusFit;
    }

    private Set<Integer> filterForDegeneracy(EpipolarTransformationFit fit, 
        SimpleMatrix allLeft, SimpleMatrix allRight) {
        
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
