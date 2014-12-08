package algorithms.imageProcessing;

import algorithms.imageProcessing.util.EightPointAlgorithmIterations;
import algorithms.util.PairFloatArray;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Given 2 sets of matched points, uses iterative solutions of the Fundamental
 * Matrix to remove outliers.
 * 
 * The methods to choose from are
 * (1) a RANSAC (Random Sample Consensus) to try 
 * 8 random points at a time for each iteration of using 
 * StereoProjectionTransformer, keeping the best fit.
 * and 
 * (2) an iterative try of StereoProjectionTransformer and if the fit is
 * better than best, remove any points with distance from epipolar larger
 * than the tolerance and repeat until reach max iter or have no improvement
 * 
 * @author nichole
 */
public class IterativeFundamentalMatrix {
    
    /*
    NOTE: may remove this method completely.  
    8 is a small of a number of points unless making an exact solution
    and have chosen the points to be placed all over the image.
    Prefer solveUsingOutlierRemoval.
    */
    public StereoProjectionTransformer solveUsingRANSAC(PairFloatArray 
        pointsLeftXY, PairFloatArray pointsRightXY) throws 
        NoSuchAlgorithmException {
    
        if (pointsLeftXY == null) {
            throw new IllegalArgumentException("pointsLeftXY cannot be null");
        }
        if (pointsRightXY == null) {
            throw new IllegalArgumentException("pointsRightXY cannot be null");
        }
        if (pointsLeftXY.getN() != pointsRightXY.getN()) {
            throw new IllegalArgumentException(
            "pointsLeftXY and pointsRightXY must have same number of points");
        }
        
        // tolerance for distance of a point position to the epipolar
        // line it belongs to.  
        // TODO: consider a way to constrain this.
        double tolerance = 3;
        
        StereoProjectionTransformer bestFMSolver = null;
        StereoProjectionTransformerFit bestFit = null;
        //long bestSubsetNumber = -1;
        
        EightPointAlgorithmIterations iterCalc = 
            new EightPointAlgorithmIterations();
        
        int nIterMax = iterCalc.estimateNIterForTwentyFivePercentOutliers(
            pointsLeftXY.getN());
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(System.currentTimeMillis());

        int nIter = 0;
        
        int nIterSameBestSubset = 0;
        //long lastBestSubsetNumber = -1;
        StereoProjectionTransformerFit lastBestFit = null;
        
        //n!/k!(n-k)!) where n = pointsLeftXY.getN() and k = 8
        //long nSubsets = MiscMath.computeNDivNMinusK(pointsLeftXY.getN(), 8) /
        //    MiscMath.factorial(8);
        
        while (nIter < nIterMax) {
            
           if (bestFit == lastBestFit) {
               if (nIterSameBestSubset > 4) {
                   break;
               }
               nIterSameBestSubset++;               
           } else {
               lastBestFit = bestFit;
               nIterSameBestSubset = 0;
           }
           
           /*
           From the n points, select 8 randomly.
           
           Ideally, would be able to predict the 'ith' subset of a subset 
           chooser like Gosper's, but haven't worked out the details of how
           to map that in one step.
           
           So, choosing the 8 points randomly, individually for now.
           */
           
           PairFloatArray subset1 = new PairFloatArray();
           PairFloatArray subset2 = new PairFloatArray();
           Set<Integer> chosen = new HashSet<Integer>();
           for (int i = 0; i < 8; i++) {
               int sel = sr.nextInt(pointsLeftXY.getN());
               while (chosen.contains(Integer.valueOf(sel))) {
                   sel = sr.nextInt(pointsLeftXY.getN());
               }
               subset1.add(pointsLeftXY.getX(sel), pointsLeftXY.getY(sel));
               subset2.add(pointsRightXY.getX(sel), pointsRightXY.getY(sel));
           }
           
           StereoProjectionTransformer current = 
               new StereoProjectionTransformer();
           
           current.calculateEpipolarProjection(subset1, subset2);
           
           StereoProjectionTransformerFit fit = current.evaluateFitForRight(
               tolerance);
           
           if ((bestFit == null) || fit.otherIsBetter(bestFit)) {
           
               bestFit = fit;
               
               bestFMSolver = current;
           }
           
           nIter++;
        }
        
        return bestFMSolver;
    }
    
    public StereoProjectionTransformer solveUsingOutlierRemoval
        (PairFloatArray pointsLeftXY, PairFloatArray pointsRightXY) {
    
        if (pointsLeftXY == null) {
            throw new IllegalArgumentException("pointsLeftXY cannot be null");
        }
        if (pointsRightXY == null) {
            throw new IllegalArgumentException("pointsRightXY cannot be null");
        }
        if (pointsLeftXY.getN() != pointsRightXY.getN()) {
            throw new IllegalArgumentException(
            "pointsLeftXY and pointsRightXY must have same number of points");
        }
        
        // tolerance for distance of a point position to the epipolar
        // line it belongs to.  
        // TODO: consider a way to constrain this.
        double tolerance = 3;
        
        StereoProjectionTransformer bestFMSolver = null;
        StereoProjectionTransformerFit bestFit = null;
        
        int nIterSameBestSubset = 0;
        StereoProjectionTransformerFit lastBestFit = null;
        
        PairFloatArray currentLeft = pointsLeftXY.copy();
        
        PairFloatArray currentRight = pointsRightXY.copy();
        
        int nIter = 0;
        
        int nIterMax = 100;
        
        while (nIter < nIterMax) {

            if ((bestFit != null) && bestFit.getOutlierIndexes().isEmpty()) {
                break;
            } else if (bestFit == lastBestFit) {
                if (nIterSameBestSubset > 4) {
                    break;
                }
                nIterSameBestSubset++;
            } else {
                lastBestFit = bestFit;
                nIterSameBestSubset = 0;
            }

            StereoProjectionTransformer current
                = new StereoProjectionTransformer();

            current.calculateEpipolarProjection(currentLeft, currentRight);

            StereoProjectionTransformerFit fit = current.evaluateFitForRight(
                tolerance);

            if ((bestFit == null) || fit.otherIsBetter(bestFit)) {

                bestFit = fit;

                bestFMSolver = current;
                
                List<Integer> rmIndexes = fit.getOutlierIndexes();
               
                if (!rmIndexes.isEmpty()) {
                    
                    // sort so we can use higher indexes first
                    Collections.sort(rmIndexes);
                
                    for (int j = (rmIndexes.size() - 1); j > -1; j--) {
                        pointsLeftXY.removeRange(j, j);
                    }
                }
            }

            nIter++;
        }
        
        return bestFMSolver;
    }
}
