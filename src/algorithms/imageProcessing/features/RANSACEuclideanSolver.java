package algorithms.imageProcessing.features;

import algorithms.SubsetChooser;
import algorithms.imageProcessing.transform.EuclideanEvaluator;
import algorithms.imageProcessing.transform.EuclideanTransformationFit;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.util.RANSACAlgorithmIterations;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.security.SecureRandom;
import java.util.Arrays;
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
     */
    public EuclideanTransformationFit calculateEuclideanTransformation(
        PairIntArray matchedLeftXY, PairIntArray matchedRightXY,
        PairIntArray outputLeftXY, PairIntArray outputRightXY, double tolerance) {

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
        
        final int nSet = 2;

        final int nPoints = matchedLeftXY.getN();
                
        // n!/(k!*(n-k)!
        long nPointsSubsets;
        try {
            nPointsSubsets = MiscMath.computeNDivKTimesNMinusK(
            nPoints, nSet);
        } catch (java.lang.ArithmeticException ex) {
            nPointsSubsets = Long.MAX_VALUE;
        }
        boolean useAllSubsets = false;
        
        SecureRandom sr = Misc.getSecureRandom();
        long seed = System.currentTimeMillis();
        log.info("SEED=" + seed + " nPoints=" + nPoints);
        sr.setSeed(seed);
        
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();
        
        EuclideanEvaluator evaluator = new EuclideanEvaluator();
        
        // consensus indexes
        EuclideanTransformationFit bestFit = null;
        
        /*
        could consider a threshold max iteration based upon the image size such
        as in (http://phototour.cs.washington.edu/ModelingTheWorld_ijcv07.pdf)
        which uses 0.6% of the maximum image dimension.
        */
        
        int outlierPercent = 50;
        int t = (int)Math.ceil((1. - (outlierPercent/100.))*nPoints);
      
        long nMaxIter;
        if (nPoints == nSet) {
            nMaxIter = 1;
            useAllSubsets = true;
        } else {
            nMaxIter = RANSACAlgorithmIterations
                .numberOfSubsamplesFor95PercentInliers(outlierPercent, nSet);
                //.numberOfSubsamplesOfSize7For95PercentInliers(outlierPercent);
        }
        
        System.out.println("nPoints=" + nPoints + " estimate for nMaxIter=" +
            nMaxIter + " (n!/(k!*(n-k)!)=" + nPointsSubsets);

        if (nMaxIter > nPointsSubsets) {
            nMaxIter = nPointsSubsets;
            useAllSubsets = true;
        }
        
        int nIter = 0;
        
        int[] selectedIndexes = new int[nSet];
        
        PairIntArray sampleLeft = new PairIntArray(nSet);
        PairIntArray sampleRight = new PairIntArray(nSet);
        
        SubsetChooser chooser = new SubsetChooser(nPoints, nSet);
                
        while (nIter < nMaxIter) {
            
            if (useAllSubsets) {
                int chk = chooser.getNextSubset(selectedIndexes);
                if (chk == -1) {
                    throw new IllegalStateException("have overrun subsets in chooser.");
                }                
            } else {
                MiscMath.chooseRandomly(sr, selectedIndexes, nPoints);
            }
            
            Arrays.sort(selectedIndexes);

            int count = 0;
            
            for (int bitIndex : selectedIndexes) {

                int idx = bitIndex;

                sampleLeft.add(matchedLeftXY.getX(idx), matchedLeftXY.getY(idx));
                                
                sampleRight.add(matchedRightXY.getX(idx), matchedRightXY.getY(idx));
                
                count++;
            }

            TransformationParameters params = tc.calulateEuclideanWithoutFilter(
                sampleLeft, sampleRight, 0, 0);

            if (params == null) {
                nIter++;
                continue;
            }

            System.out.printf("%d out of %d iterations\n", nIter, nMaxIter);
            System.out.flush();
            
            EuclideanTransformationFit fit = evaluator.evaluate(matchedLeftXY,
                matchedRightXY, params, tolerance);
            
            assert(fit != null);
                                    
            int nInliers = fit.getInlierIndexes().size();
            if (!fit.isBetter(bestFit)) {
                nIter++;
                continue;
            }
            
            if (nInliers > nSet && nInliers > t) {
                // redo the transformation with all inliers
                PairIntArray inliersLeftXY = new PairIntArray(nInliers);
                PairIntArray inliersRightXY = new PairIntArray(nInliers);
                int countI = 0;
                for (Integer idx : fit.getInlierIndexes()) {
                    int idxInt = idx.intValue();
                    inliersLeftXY.add(matchedLeftXY.getX(idxInt), matchedLeftXY.getY(idxInt));
                    inliersRightXY.add(matchedRightXY.getX(idxInt), matchedRightXY.getY(idxInt));
                    countI++;
                }
                
                TransformationParameters params2 = tc.calulateEuclideanWithoutFilter(
                    inliersLeftXY, inliersRightXY, 0, 0);
                
                EuclideanTransformationFit fit2 = evaluator.evaluate(matchedLeftXY,
                    matchedRightXY, params2, tolerance);
                
                if (fit2 != null && fit2.isBetter(fit)) {
                    fit = fit2;
                }
                System.out.println("new local best fit: " + fit2.toString());
                System.out.flush();                
            }
                       
            if (fit.isBetter(bestFit)) {
                int nb = (bestFit != null) ? bestFit.getInlierIndexes().size() : nSet+1;
                int nf = fit.getInlierIndexes().size();
                
                bestFit = fit;
                
                System.out.println("**best fit: " + bestFit.toString());
                System.out.flush();
                
                // recalculate nMaxIter
                if ((nf > nb) && nMaxIter > 1) {
                    double outlierPercentI = 100.*
                        (double)(nPoints - bestFit.getInlierIndexes().size()) /
                        (double)nPoints;
                    if (outlierPercentI < outlierPercent) {
                        outlierPercent = (int)Math.ceil(outlierPercentI);
                        if (outlierPercent < 5) {
                            outlierPercent = 5;
                        }
                        assert(outlierPercent <= 50);
                        nMaxIter = RANSACAlgorithmIterations
                            .numberOfSubsamplesFor95PercentInliers(outlierPercent, nSet);
                        if (nMaxIter > nPointsSubsets) {
                            nMaxIter = nPointsSubsets;
                            useAllSubsets = true;
                        }
                    }
                }
            }
            
            nIter++;
        }

        if (bestFit == null) {
            log.info("no solution.  nIter=" + nIter);
            return null;
        }
                
        // write to output and convert the coordinate indexes to the original point indexes
        List<Integer> inlierIndexes = bestFit.getInlierIndexes();
        for (int i = 0; i < inlierIndexes.size(); ++i) {
            Integer index = inlierIndexes.get(i);
            int idx = index.intValue();
            outputLeftXY.add(
                Math.round(matchedLeftXY.getX(idx)),
                Math.round(matchedLeftXY.getY(idx)));
            outputRightXY.add(
                Math.round(matchedRightXY.getX(idx)),
                Math.round(matchedRightXY.getY(idx)));
        }
      
        log.fine("nIter=" + nIter);

        log.fine("final fit: " + bestFit.toString());

        return bestFit;
    }

}
