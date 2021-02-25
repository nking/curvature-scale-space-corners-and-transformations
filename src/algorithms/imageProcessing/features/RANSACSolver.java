package algorithms.imageProcessing.features;

import algorithms.SubsetChooser;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.transform.Distances;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.EpipolarTransformer;
import algorithms.imageProcessing.transform.Util;
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
    
    //TODO: edit to be able to choose nSet = 8 or 7

    private boolean debug = true;

    private Logger log = Logger.getLogger(this.getClass().getName());

    /**
     * calculate the epipolar transformation among the given points with the
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
        
        DenseMatrix input1 =
            Util.rewriteInto3ColumnMatrix(matchedLeftXY);

        DenseMatrix input2 =
            Util.rewriteInto3ColumnMatrix(matchedRightXY);
        
        return calculateEpipolarProjection(input1, input2,
            outputLeftXY, outputRightXY, tolerance);
    }

    /**
     * calculate the epipolar transformation among the given points with the
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
        final double tolerance) {
        
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
        if (matchedLeftXY.numRows() != 3) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the algorithms require 3 rows representing x, y, and '1' values."
                + " matchedLeftXY.n=" + matchedLeftXY.numColumns());
        }
        if (matchedLeftXY.numColumns() != matchedRightXY.numColumns() ||
            matchedLeftXY.numRows() != matchedRightXY.numRows()) {
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
        -- for each iteration of solving epipolar transformation using a sample
           of size 7, the resulting fundamental matrix is evaluated on the
           all points of the dataset.
           If the number of inliers is T or more, the fit is re-done with all
           of the points, where 
           T = (1. - outlierPercentage) * (total number of data points)
           The best fitting for all iterations as defined by number of inliers 
           and standard deviation from an epipolar line, is kept each time.
        -- at the end of each iteration, the number of iterations is then 
           re-calculated if it can be reduced.
        
        NOTE that the sub-samples are selected randomly from all possible
        sub-samples of the nPoints unless the number of all possible 
        sub-samples is smaller than the expected number of iterations for 95%
        probability of a good sub-sample.  In the later case, all sub-samples
        are tried.
        
        */

        final int nSet = 7;

        final int nPoints = matchedLeftXY.numColumns();
                
        // n!/(k!*(n-k)!
        final long nPointsSubsets = MiscMath.computeNDivKTimesNMinusK(nPoints, nSet);
        boolean useAllSubsets = false;
        
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
        int t = (int)Math.ceil((1. - outlierPercent)*nPoints);
        
        long nMaxIter;
        if (nPoints == nSet) {
            nMaxIter = 1;
            useAllSubsets = true;
        } else {
            nMaxIter = RANSACAlgorithmIterations
                .numberOfSubsamplesOfSize7For95PercentInliers(outlierPercent);
        }
        
        System.out.println("nPoints=" + nPoints + " estimate for nMaxIter=" +
            nMaxIter + " (n!/(k!*(n-k)!)=" + nPointsSubsets);

        if (nMaxIter > nPointsSubsets) {
            nMaxIter = nPointsSubsets;
            useAllSubsets = true;
        }
        
        int nIter = 0;
        
        int[] selectedIndexes = new int[nSet];
        
        DenseMatrix sampleLeft = new DenseMatrix(3, nSet);
        DenseMatrix sampleRight = new DenseMatrix(3, nSet);
        // initialize the unchanging 3rd dimension
        for (int i = 0; i < nSet; ++i) {
            sampleLeft.set(2, i, 1);
            sampleRight.set(2, i, 1);
        }
        
        SubsetChooser chooser = new SubsetChooser(nPoints, nSet);
        
        Distances distances = new Distances();
        
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

                sampleLeft.set(0, count, matchedLeftXY.get(0, idx));
                sampleLeft.set(1, count, matchedLeftXY.get(1, idx));
                                
                sampleRight.set(0, count, matchedRightXY.get(0, idx));
                sampleRight.set(1, count, matchedRightXY.get(1, idx));
                
                count++;
            }
            
//TODO: edit here to return data structure that includes normalized
        
            // determine matrix from 7 points.
            List<DenseMatrix> fms = spTransformer.calculateEpipolarProjectionFor7Points(
                sampleLeft, sampleRight);

            System.out.printf("%d out of %d iterations\n", nIter, nMaxIter);
            System.out.flush();
            
            if (fms == null || fms.isEmpty()) {
                nIter++;
                continue;
            }
            
            // use point dist to epipolar lines to estimate errors of sample
            EpipolarTransformationFit fit = null;
            
            double chiSqStatFactor = 7.82;
            int plunder;
            int bestPlunderCost = Integer.MAX_VALUE;
            int bestPlunderIdx = -1;
            
            for (int fIdx = 0; fIdx < fms.size(); ++fIdx) {
                
                DenseMatrix fm = fms.get(fIdx);
                
                // this uses a threshold of chiSqStatFactor * standard deviation
                //   of the mean of the distances to remove outliers
                EpipolarTransformationFit fitI = distances.calculateError2(fm, 
                    matchedLeftXY, matchedRightXY, errorType, chiSqStatFactor);
                      
                int nInliers = fitI.getInlierIndexes().size();
                int nOutliers = nSet - nInliers;
                
                /*
                The plunder-dl scoring can be used for comparison between different models.
                for example, comparing results of the 7-point and 8-point 
                solutions or comparing 7-point projection to 6-point affine, etc.
                
                plunder-dl is from equation 33 of
                Torr, Zisserman, & Maybank 1996, 
                “Robust Detection of Degenerate Configurations whilst Estimating 
                the Fundamental Matrix"
                https://www.robots.ox.ac.uk/~phst/Papers/CVIU97/m.ps.gz
                 EQN 33: PL = DOF + (4*n_o + n_i dimension of model)
                               where n_i = number of inliers
                               n_o = number of outliers
                               DOF = 7 for this solver
                n=7               PL = DOF + 4*n_o + n_i* (model_dimension)
                     ni=7, no=0   PL = 7   + 0     + 0 * md
                     ni=5, no=2   PL = 7   + 8     + 8 * md
                     ni=4, no=3   PL = 7   + 12    + 28 * md
                PLUNDER stands for Pick Least UNDEgenerate Randomly, Description Length
                
                For nPoints=8, model_dimension = 1.
                for nPoints=7 amd only 1 solution in the cubic constraints, model_dimension=2,
                else for nPoints=7, model_dimension = 3.
                
                Will use model_dimension=2 here and keep the smallest pluder score.
                */
                
                // fitI.isBetter() : comparison to other fit by the number of 
                //     inliers, else if tie, mean of errors, else if tie, 
                //     mean of standard deviation of mean of errors, else 
                //     returns false.
                plunder = 7 + 4*nOutliers + 3*nInliers;
                
                if (nInliers >= nSet && (plunder <= bestPlunderCost) 
                    && fitI.isBetter(fit)) {
           
                    bestPlunderCost = plunder;
                    bestPlunderIdx = fIdx;
                    
 //TODO: continue revising from here
 
                    // redo the transformation with all inliers
                    DenseMatrix inliersLeftXY = new DenseMatrix(3, nInliers);
                    DenseMatrix inliersRightXY = new DenseMatrix(3, nInliers);
                    int countI = 0;
                    for (Integer idx : fitI.getInlierIndexes()) {
                        int idxInt = idx.intValue();
                        inliersLeftXY.set(0, countI, matchedLeftXY.get(0, idxInt));
                        inliersLeftXY.set(1, countI, matchedLeftXY.get(1, idxInt));
                        inliersLeftXY.set(2, countI, 1);
                        inliersRightXY.set(0, countI, matchedRightXY.get(0, idxInt));
                        inliersRightXY.set(1, countI, matchedRightXY.get(1, idxInt));
                        inliersRightXY.set(2, countI, 1);
                        countI++;
                    }

                    DenseMatrix fm2 = spTransformer.calculateEpipolarProjection(
                        inliersLeftXY, inliersRightXY);

                    if (fm2 != null) {

                        EpipolarTransformationFit fit2 = distances
                            .calculateError(fm2, matchedLeftXY, matchedRightXY,
                            errorType, tolerance);

                        if (fit2 != null && fit2.isBetter(fitI)) {
                            fitI = fit2;
                        }
                    }
                    
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
                        assert(outlierPercent < 50);
                        nMaxIter = RANSACAlgorithmIterations
                            .numberOfSubsamplesOfSize7For95PercentInliers(outlierPercent);
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
                (int)Math.round(matchedLeftXY.get(0, idx)),
                (int)Math.round(matchedLeftXY.get(1, idx)));
            outputRightXY.add(
                (int)Math.round(matchedRightXY.get(0, idx)),
                (int)Math.round(matchedRightXY.get(1, idx)));
        }
      
        log.fine("nIter=" + nIter);

        log.fine("final fit: " + bestFit.toString());

        return bestFit;
    }
    
    /**
     * if assume gaussian errors, for 1 degree of freedom (i.e. fitting
     * a line, fundamental matrix, d^2 = 3.84*(st.dev^2)
     * @param standardDeviation
     * @return 
     */
    public static double estimateToleranceForDOF1(double standardDeviation) {
        double d = Math.sqrt(3.84*standardDeviation*standardDeviation);
        return d;
    }
    
    /**
     * if assume gaussian errors, for 2 degree of freedom (i.e. fitting
     * a line, fundamental matrix, d^2 = 5.99*(st.dev^2)
     * @param standardDeviation
     * @return 
     */
    public static double estimateToleranceForDOF2(double standardDeviation) {
        double d = Math.sqrt(5.99*standardDeviation*standardDeviation);
        return d;
    }
}
