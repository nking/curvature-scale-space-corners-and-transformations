package algorithms.imageProcessing.features;

import algorithms.combPerm.SubsetChooser;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.transform.Distances;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.EpipolarTransformer;
import algorithms.imageProcessing.transform.Util;
import algorithms.imageProcessing.util.RANSACAlgorithmIterations;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import algorithms.util.PairIntArray;

import java.io.IOException;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.logging.Logger;
import no.uib.cipr.matrix.DenseMatrix;

/**
 * given matched point lists, determine the best epipolar solution using a
 * 7-point epipolar calculation and random draws of 7 points from the
 * matched point lists under the assumption that some of the matched points
 * are not true (correct) matches.
 *
 * from wikipedia:
Random sample consensus (RANSAC) is an iterative method to estimate parameters 
of a mathematical model from a set of observed data that contains outliers, 
when outliers are to be accorded no influence on the values of the estimates. 
Subsets are drawn from the sample a number of times such that the probability of drawing
a sample that is a specified percentage of inliers is met.  The best model is determined from the best fitting
subset and that model is then applied to all data.
calculating the number of iterations needed for finding a subset that is all inliers is
an important part of the algorithm to keep the runtime tractable.
* 
 * <pre>
 * useful reading:
 * add references in comments below here..
 * 
 * http://6.869.csail.mit.edu/fa12/lectures/lecture13ransac/lecture13ransac.pdf
 * and
 * http://www.dtic.mil/dtic/tr/fulltext/u2/a460585.pdf
 * </pre>
 *
 * Note: to compare different geometric model results:
 * <pre>
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
 * </pre>
 * @author nichole
 */
public class RANSACSolver {
    
    //TODO: edit to be able to choose nSet = 8 or 7

    private boolean debug = true;

    private Logger log = Logger.getLogger(this.getClass().getName());

    /**
     * calculate the epipolar transformation among the given points with the
     * assumption that some of the points in the matched lists are not
     * true matches.   NOTE: for best results, one should perform unit standard
     * normalization on the correspondence first.
     *
     * @param leftCorres left correspondence holding (x,y) points from left image
     * in format 3 X nData matrix with rows being x, y, and 1's respectively
     * @param rightCorres right correspondence holding (x,y) points from left image
     * in format 3 X nData matrix with rows being x, y, and 1's respectively
     * @param errorType algorithm used to evaluate the fit of the fundamental matrix solutions.
     * @param useToleranceAsStatFactor if set to false, tolerance is used as
     * a fixed number in outlier removal, else if set to true, tolerance
     * is used as the chi-squared statistic factor for the standard deviation
     * of errors use in outlier removal.
     * @param tolerance tolerance in distance from epipolar line for a point to 
     * be an inlier in the final fit.   NOTE: if useToleranceAsStatFactor is true,
     * it is interpreted as a chiSqStatFactor which is then used as 
     * tolerance = tolerance * standard deviation of the mean distance errors.
     * @param reCalcIterations if true, upon each better fit found, the 
     * outlier percentage is re-estimated and then the number of iterations necessary for 95%
     * probability that sample has all good points.
     * @param calibrated if true, solves for the Essential Matrix, else solves
     * for the Fundamental Matrix.  The difference is in the diagonal used for
     * dimension reduction.
     * @return
     */
    public EpipolarTransformationFit calculateEpipolarProjection(
        final DenseMatrix leftCorres, final DenseMatrix rightCorres,
        ErrorType errorType,
        boolean useToleranceAsStatFactor, final double tolerance,
        boolean reCalcIterations, boolean calibrated) {

        double[][] left = convertX(leftCorres);
        double[][] right = convertX(rightCorres);

        long seed = System.currentTimeMillis();
        log.info("SEED=" + seed);
        Random rand = new Random(seed);

        return run(left, right, errorType, useToleranceAsStatFactor, tolerance, reCalcIterations, calibrated, rand);
    }

    /**
     * calculate the epipolar transformation among the given points with the
     * assumption that some of the points in the matched lists are not
     * true matches.   NOTE: for best results, one should perform unit standard
     * normalization on the correspondence first.
     *
     * @param leftCorres left correspondence holding (x,y) points from left image
     * in format 3 X nData matrix with rows being x, y, and 1's respectively
     * @param rightCorres right correspondence holding (x,y) points from left image
     * in format 3 X nData matrix with rows being x, y, and 1's respectively
     * @param errorType algorithm used to evaluate the fit of the fundamental matrix solutions.
     * @param useToleranceAsStatFactor if set to false, tolerance is used as
     * a fixed number in outlier removal, else if set to true, tolerance
     * is used as the chi-squared statistic factor for the standard deviation
     * of errors use in outlier removal.
     * @param tolerance tolerance in distance from epipolar line for a point to
     * be an inlier in the final fit.   NOTE: if useToleranceAsStatFactor is true,
     * it is interpreted as a chiSqStatFactor which is then used as
     * tolerance = tolerance * standard deviation of the mean distance errors.
     * @param reCalcIterations if true, upon each better fit found, the
     * outlier percentage is re-estimated and then the number of iterations necessary for 95%
     * probability that sample has all good points.
     * @param calibrated if true, solves for the Essential Matrix, else solves
     * for the Fundamental Matrix.  The difference is in the diagonal used for
     * dimension reduction.
     * @param rand instance of Random to use in choosing subsets of correspondence.
     * @return
     */
    public EpipolarTransformationFit calculateEpipolarProjection(final DenseMatrix leftCorres, final DenseMatrix rightCorres,
            ErrorType errorType,boolean useToleranceAsStatFactor, final double tolerance,
            boolean reCalcIterations, boolean calibrated, Random rand) {

        double[][] left = convertX(leftCorres);
        double[][] right = convertX(rightCorres);

        return run(left, right, errorType, useToleranceAsStatFactor, tolerance, reCalcIterations, calibrated, rand);
    }

    /**
     * calculate the epipolar transformation among the given points with the
     * assumption that some of the points in the matched lists are not
     * true matches.   NOTE: for best results, one should perform unit standard
     * normalization on the correspondence first.
     *
     * @param leftCorres left correspondence holding (x,y) points from left image
     * in format 3 X nData matrix with rows being x, y, and 1's respectively
     * @param rightCorres right correspondence holding (x,y) points from left image
     * in format 3 X nData matrix with rows being x, y, and 1's respectively
     * @param errorType algorithm used to evaluate the fit of the fundamental matrix solutions.
     * @param useToleranceAsStatFactor if set to false, tolerance is used as
     * a fixed number in outlier removal, else if set to true, tolerance
     * is used as the chi-squared statistic factor for the standard deviation
     * of errors use in outlier removal.
     * @param tolerance tolerance in distance from epipolar line for a point to
     * be an inlier in the final fit.   NOTE: if useToleranceAsStatFactor is true,
     * it is interpreted as a chiSqStatFactor which is then used as
     * tolerance = tolerance * standard deviation of the mean distance errors.
     * @param reCalcIterations if true, upon each better fit found, the
     * outlier percentage is re-estimated and then the number of iterations necessary for 95%
     * probability that sample has all good points.
     * @param calibrated if true, solves for the Essential Matrix, else solves
     * for the Fundamental Matrix.  The difference is in the diagonal used for
     * dimension reduction.
     * @param rand instance of Random to use in choosing subsets of correspondence.
     * @return
     */
    public EpipolarTransformationFit calculateEpipolarProjection(
            final double[][] leftCorres, final double[][] rightCorres,
            ErrorType errorType,
            boolean useToleranceAsStatFactor, final double tolerance,
            boolean reCalcIterations, boolean calibrated, Random rand) {

        return run(leftCorres, rightCorres, errorType, useToleranceAsStatFactor, tolerance, reCalcIterations, calibrated, rand);
    }

    /**
     * calculate the epipolar transformation among the given points with the
     * assumption that some of the points in the matched lists are not
     * true matches.   NOTE: for best results, one should perform unit standard
     * normalization on the correspondence first.
     *
     * @param leftCorres left correspondence holding (x,y) points from left image
     * in format 3 X nData matrix with rows being x, y, and 1's respectively
     * @param rightCorres right correspondence holding (x,y) points from left image
     * in format 3 X nData matrix with rows being x, y, and 1's respectively
     * @param errorType algorithm used to evaluate the fit of the fundamental matrix solutions.
     * @param useToleranceAsStatFactor if set to false, tolerance is used as
     * a fixed number in outlier removal, else if set to true, tolerance
     * is used as the chi-squared statistic factor for the standard deviation
     * of errors use in outlier removal.
     * @param tolerance tolerance in distance from epipolar line for a point to
     * be an inlier in the final fit.   NOTE: if useToleranceAsStatFactor is true,
     * it is interpreted as a chiSqStatFactor which is then used as
     * tolerance = tolerance * standard deviation of the mean distance errors.
     * @param reCalcIterations if true, upon each better fit found, the
     * outlier percentage is re-estimated and then the number of iterations necessary for 95%
     * probability that sample has all good points.
     * @param calibrated if true, solves for the Essential Matrix, else solves
     * for the Fundamental Matrix.  The difference is in the diagonal used for
     * dimension reduction.
     * @return
     */
    public EpipolarTransformationFit calculateEpipolarProjection(
            final double[][] leftCorres, final double[][] rightCorres,
            ErrorType errorType,
            boolean useToleranceAsStatFactor, final double tolerance,
            boolean reCalcIterations, boolean calibrated) {

        long seed = System.currentTimeMillis();
        log.info("SEED=" + seed);
        Random rand = new Random(seed);

        return run(leftCorres, rightCorres, errorType, useToleranceAsStatFactor, tolerance, reCalcIterations, calibrated, rand);
    }

    private double[][] convertX(DenseMatrix x) {
        double[][] d = MatrixUtil.zeros(x.numRows(), x.numColumns());
        int i, j;
        for (i = 0; i < d.length; ++i) {
            for (j = 0; j < d[i].length; ++j) {
                d[i][j] = x.get(i, j);
            }
        }
        return d;
    }

    private  EpipolarTransformationFit run(
        final double[][] left, final double[][] right,
        ErrorType errorType,
        boolean useToleranceAsStatFactor, final double tolerance,
        boolean reCalcIterations, boolean calibrated, Random rand) {

        int nPoints = left[0].length;
        final int nSet = 7;
        
        if (nPoints < nSet) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the algorithms require 7 or more points."
                + " leftCorres.n=" + left[0].length);
        }
        if (left.length != 3) {
            // cannot use this algorithm.
            throw new IllegalArgumentException(
                "the algorithms require 3 rows representing x, y, and '1' values."
                + " leftCorres.n=" + left.length);
        }
        if (left[0].length != right[0].length ||
            left.length != right.length) {
            throw new IllegalArgumentException(
                "leftCorres and rightCorres bmust be the same size");
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
       
        // n!/(k!*(n-k)!
        //final long nPointsSubsets = MiscMath.computeNDivKTimesNMinusK(nPoints, nSet);
        boolean useAllSubsets = false;

        EpipolarTransformer spTransformer = new EpipolarTransformer();
                
        // consensus best fit and inlier indexes
        EpipolarTransformationFit bestFit = null;
        
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
            /* The number of subsamples required to ensure T >= 0.95 
            for given outlierPercent as fraction of contaminated data, 
            where T is the probability that all the data 
            points selected in one subsample are non-outliers.
            */
            // maximum is 382 for nSet=7 and outlierPercent=50%, and for nSet=8 it is 766
            nMaxIter = RANSACAlgorithmIterations
                .numberOfSubsamplesFor95PercentInliers(outlierPercent, nSet);
        }
        
        //System.out.println("nPoints=" + nPoints + " estimate for nMaxIter=" +
        //    nMaxIter + ", (n!/(k!*(n-k)!)=" + nPointsSubsets);

        // max iter for nSet=7 is 382 for 50%.  11!/(7!*(11-7)!)=330
        // max iter for nSet=8 is 766 for 50%.  12!/(8!*(12-8)!)=495
        if ((nSet==7 && nPoints < 12) || (nSet == 8 && nPoints < 13)) {
            nMaxIter = Math.min(nMaxIter, MiscMath0.computeNDivKTimesNMinusK(nPoints, nSet));
            useAllSubsets = true;
        }
        
        int nIter = 0;

        double eps = 1.e-6;

        //TODO: finish changes to replace DenseMatrix with double[][] throughout this method
        DenseMatrix leftCorres = new DenseMatrix(left);
        DenseMatrix rightCorres = new DenseMatrix(right);

        /*
        double[][] sampleLeft = MatrixUtil.zeros(3, nSet);
        double[][] sampleRight = MatrixUtil.zeros(3, nSet);
        // initialize the unchanging 3rd dimension
        Arrays.fill(sampleLeft[2], 1);
        Arrays.fill(sampleRight[2], 1);
        */

        DenseMatrix sampleLeft = new DenseMatrix(3, nSet);
        DenseMatrix sampleRight = new DenseMatrix(3, nSet);
        // initialize the unchanging 3rd dimension
        for (int i = 0; i < nSet; ++i) {
            sampleLeft.set(2, i, 1);
            sampleRight.set(2, i, 1);
        }

        int[] selectedIndexes = new int[nSet];

        SubsetChooser chooser = null;
        if (useAllSubsets) {
            chooser = new SubsetChooser(nPoints, nSet);
        }

        Distances distances = new Distances();
        
        while (nIter < nMaxIter) {
            if (useAllSubsets) {
                int chk = chooser.getNextSubset(selectedIndexes);
                if (chk == -1) {
                    throw new IllegalStateException("have overrun subsets in chooser.");
                }                
            } else {
                MiscMath.chooseRandomly(rand, selectedIndexes, nPoints);
            }
            
            Arrays.sort(selectedIndexes);

            int count = 0;
            for (int bitIndex : selectedIndexes) {
                /*sampleLeft[0][count] = leftCorres[0][bitIndex];
                sampleLeft[1][count] = leftCorres[1][bitIndex];
                sampleRight[0][count] = rightCorres[0][bitIndex];
                sampleRight[1][count] = rightCorres[1][bitIndex];*/
                sampleLeft.set(0, count, leftCorres.get(0, bitIndex));
                sampleLeft.set(1, count, leftCorres.get(1, bitIndex));
                sampleRight.set(0, count, rightCorres.get(0, bitIndex));
                sampleRight.set(1, count, rightCorres.get(1, bitIndex));
                count++;
            }
/*
            if (MiscMath0.areColinear(MatrixUtil.convertToRowMajor(sampleLeft), eps)
                    || MiscMath0.areColinear(MatrixUtil.convertToRowMajor(sampleRight), eps)) {
                ++nIter;
                continue;
            }
*/
            // calculates 7-point solutions then filters using chirality checks.
            List<DenseMatrix> fms = spTransformer
                .calculateEpipolarProjectionFor7Points(new DenseMatrix(sampleLeft), new DenseMatrix(sampleRight));

            //TODO: finish transformation to using double[][] instead of DenseMatrix
            // calculates 7-point solutions then filters using chirality checks.
            /*double[][] fms = null;
            if (nSet == 7) {
                fms = spTransformer.calculateEpipolarProjectionUsing7Points(sampleLeft, sampleLeft);
            } else {
                fms = spTransformer.calculateEpipolarProjection2(sampleLeft, sampleRight, calibrated);
            }*/

            if (fms == null || fms.isEmpty()) {
                nIter++;
                continue;
            }
            
            // evaluate fms solutions on all points and keep best and compare
            // that to best overall solution
            
            EpipolarTransformationFit fit = null;
            
            // fit.isBetter() : comparison to other fit by the number of 
            //     inliers, else if tie, mean of errors, else if tie, 
            //     mean of standard deviation of mean of errors, else 
            //     returns false
            for (DenseMatrix fm : fms) {
                
                EpipolarTransformationFit fitI = null;
                // evaluate all points using the solution from the sub-sample
                if (useToleranceAsStatFactor) {
                    fitI = distances.calculateError2(fm,
                        leftCorres, rightCorres, errorType, tolerance);
                } else {
                    fitI = distances.calculateError(fm,
                        leftCorres, rightCorres, errorType, tolerance);
                }
                
                int nInliers = fitI.getInlierIndexes().size();
                if (nInliers >= nSet && fitI.isBetter(fit)) {
                    if (nInliers > t && nInliers > nSet) {
                        // redo the FM transformation with all inliers
                        DenseMatrix inliersLeftXY = EpipolarTransformer
                            .extractIndices(leftCorres, fitI.getInlierIndexes());
                        DenseMatrix inliersRightXY = EpipolarTransformer
                            .extractIndices(rightCorres, fitI.getInlierIndexes());
                        
                        DenseMatrix fm2 = spTransformer.calculateEpipolarProjection(
                            inliersLeftXY, inliersRightXY, calibrated);
                        
                        if (fm2 != null) {
                            EpipolarTransformationFit fit2 = null;
                            if (useToleranceAsStatFactor) {
                                fit2 = distances.calculateError2(fm2,
                                    leftCorres, rightCorres, errorType, tolerance);
                            } else {
                                fit2 = distances.calculateError(fm2,
                                    leftCorres, rightCorres, errorType, tolerance);
                            }
                            if (fit2 != null && fit2.isBetter(fitI)) {
                                fitI = fit2;
                            }
                        }
                    }
                    //System.out.println(" new local best fit: " + fitI.toString());
                    //System.out.flush();
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
                
                //System.out.println("**best fit: " + bestFit.toString());
                //System.out.flush();
                
                // recalculate nMaxIter
                if (reCalcIterations && (nf > nb) && nMaxIter > 1) {
                    double outlierPercentI = 100.*
                        (double)(nPoints - bestFit.getInlierIndexes().size()) / (double)nPoints;
                    if (outlierPercentI < outlierPercent) {
                        outlierPercent = (int)Math.ceil(outlierPercentI);
                        if (outlierPercent < 5) {
                            outlierPercent = 5;
                        }
                        assert(outlierPercent < 50);
                        nMaxIter = RANSACAlgorithmIterations
                            .numberOfSubsamplesOfSize7For95PercentInliers(outlierPercent);
                        // max iter for nSet=7 is 382 for 50%.  11!/(7!*(11-7)!)=330
                        if ((nSet==7 && nPoints < 12) || (nSet==8 && nPoints < 13)) {
                            nMaxIter = Math.min(nMaxIter, MiscMath0.computeNDivKTimesNMinusK(nPoints, nSet));
                            useAllSubsets = true;
                            if (chooser == null) {
                                chooser = new SubsetChooser(nPoints, nSet);
                            }
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
        
        log.info("nIter=" + nIter);
        
        log.fine("final best fit to all points: " + bestFit.toString());

        return bestFit;
    }
    
}
