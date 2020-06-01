package algorithms.imageProcessing.util;

import algorithms.misc.MiscMath;

/**
 * estimates the number of iterations for a RANSAC type algorithm that 
 * should be used to ensure with 99% certainty or other that a set of points is 
 * findable as a subset of n matchable points.
<pre>
nPoints is the number of points that contain true and false matches.
        
nTruePoints is the number of true matches within nPoints.
        
k is the number of points to draw at one time for a RANSAC sample test iteration.

nCombinations is the number of combinations of samples of size k is nPoints!/(k!*(nPoints-k)!).

For a sample, the first draw of a matching point out of nPoints has 
possibility of being all 'true' points = (nTruePoints/nPoints).
For the same sample, the second draw in the sample has the possibility 
of being all 'true' points = (nTruePoints - 1)/(nPoints - 1).
etc.   The possibility that the sample is composed of all 'true' points is
then (nTruePoints/nPoints) * ((nTruePoints - 1)/(nPoints - 1)) *
         ((nTruePoints - 2)/(nPoints - 2)) ... ((nTruePoints - k - 1)/(nPoints - k - 1))

The total number of 'true' samples within all possible combinations of size k 
from nPoints is then
nCombinations * ((nTruePoints/nPoints)*(nTruePoints - 1)/(nPoints - 1)...((nTruePoints - k-1)/(nPoints - k-1))
(this is verified in algorithms.stats.ransac.CountingTest)

The fraction of all 'true' samples, 
that is the number of samples composed of all 'true' points divided by number of 
all possible samples,
is then just the possibility that a sample is all 'true', that is
     ((nTruePoints/nPoints)*(nTruePoints - 1)/(nPoints - 1)...((nTruePoints - k-1)/(nPoints - k-1))

let pSample be the fraction just calculated.

knowing pSample for a sample now, need to calculate the number of samples that
need to be drawn in order to have a high confidence that at least one
sample was a true sample.

Naively, that would be the inverse of pSample.

Now that have the statistics in a form of 2 states, can use binomial statistics,
the binomial theorem, to determine the number of fails before the first
success (sample is 'true').   This is similar to a Geometric distribution when
simplified below.
          
    m = the number of success trials which is the minimum here, '1'
    nIter is the number of iterations needed 
    
    P(m|pSample,nIter) = pSample^m * (1 - pSample)^(nIter-m) * nCombinations
                       = pSample * (1 - pSample)^(nIter-1) * nCombinations
                       = pSample * nCombinations * (1 - pSample)^(nIter-1)
                       
               having defined m = 1 and knowing pSample represents the number of
               true samples divided by nCombinations, can simplify
               pSample * nCombinations as '1'
               
    P(m|pSample,nIter) = (1 - pSample)^(nIter-1)
    
    P(m|pSample,nIter) is 1 - Pconfidence where Pconfidence used here will be 0.99
     
    Then solving for nIter:
         math.log(1. - Pconfid) = (nIter-1) * math.log(1. - pSample)
    then
         nIter = (math.log(1. - Pconfid) / math.log(1. - pSample)) + 1
         
</pre>
 * @author nichole
 */
public class RANSACAlgorithmIterations {

    public long estimateNIterFor99PercentConfidence(
        int nPoints, int sampleSize, double expectedFractionTruePoints) {
        
        if (sampleSize < 1) {
            throw new IllegalArgumentException("sampleSize must be > 0");
        }
        if (nPoints < 1 || (nPoints > 1789)) {
            throw new IllegalArgumentException("nPoints must be a positive number < 1790");
        }
        if (expectedFractionTruePoints < 0.0000001 || (expectedFractionTruePoints > 1.0)) {
            throw new IllegalArgumentException(
            "expectedFractionTruePoints must be larger than 0 and less than 1");
        }
        
        double pSample = calculateTrueSampleProbability(nPoints, sampleSize, 
            expectedFractionTruePoints);
        
        // pSample must be larger than approx 1e-16 and <= 0.99
        
        long nIter = Math.round(Math.log(1. - 0.99) / Math.log(1. - pSample)) + 1;
           
        return nIter;
    }
    
    public long estimateNIterFor99PercentConfidenceDegenerate(
        int nPoints, int nPointsIncludingDegenerate, int sampleSize, 
        double expectedFractionTruePoints) {
        
        if (sampleSize < 1) {
            throw new IllegalArgumentException("sampleSize must be > 0");
        }
        if (nPoints < 1 || (nPoints > 1789) ) {
            throw new IllegalArgumentException("nPoints must be a positive " +
            " number < 1790 or even smaller considering multiplicity");
        }
        if (expectedFractionTruePoints < 0.0000001 || expectedFractionTruePoints > 1.0) {
            throw new IllegalArgumentException(
            "expectedFractionTruePoints must be larger than 0 and less than 1");
        }
        if (nPoints <= sampleSize) {
            throw new IllegalArgumentException("nPoints must be larger than sampleSize");
        }
        
        double pSample = calculateTrueSampleProbabilityForMultiplyMatched(nPoints, 
            nPointsIncludingDegenerate, sampleSize, expectedFractionTruePoints);

        long nIter;
        
        if (pSample == 0) {
            return Long.MAX_VALUE;
        } else if (pSample < 1E-16) {
            // return naive estimate
            return Math.round(1./pSample);
        } else {
            nIter = Math.round(Math.log(1. - 0.99) / Math.log(1. - pSample)) + 1;
        }
        
        return nIter;
    }
    
    /**
     * calculate the probability of drawing an all 'T' sample of sampleSize 
     * points from a universe of nPoints of which 
     * (expectedFractionTruePoints*nPoints) are 'T'.  
     * The calculation is drawing each point without replacement,
     * that is (nTruePts/nPts) * ((nTrue-1)/(nPoints-1)) * ((nTrue-2)/(nPoints-2)) ...
     * 
     * If nTruePts is less than sampleSize (due to small expectedFractionTruePoints
     * for example), the calculation becomes 
     * @param nPoints
     * @param sampleSize
     * @param expectedFractionTruePoints
     * @return 
     */
    public double calculateTrueSampleProbability(int nPoints, int sampleSize,
        double expectedFractionTruePoints) {
        
        if (sampleSize < 1) {
            throw new IllegalArgumentException("sampleSize must be > 0");
        }
        if (nPoints < 1 || (nPoints > 1789) ) {
            throw new IllegalArgumentException("nPoints must be a positive " 
                + " number < 1790 or even smaller considering multiplicity");
        }
        if (expectedFractionTruePoints < 0.0000001 || expectedFractionTruePoints > 1.0) {
            throw new IllegalArgumentException(
            "expectedFractionTruePoints must be larger than 0 and less than 1");
        }
        if (nPoints < sampleSize) {
            throw new IllegalArgumentException("nPoints must be larger than sampleSize");
        }
        
        double nTruePoints = nPoints * expectedFractionTruePoints;
        
        double factor = 1;
        
        if (nTruePoints >= sampleSize) {        
            for (int i = 0; i < sampleSize; ++i) {
                factor *= (nTruePoints - i)/(double)(nPoints - i);
            }            
        } else {
            //know that the number is still smaller than (1./nPoints)^sampleSize.
            //For nTruePoints of 1, looks like (1./nPoints)*(1./(nPoints-sampleSize))^(sampleSize-1)
            for (int i = 0; i < (int)nTruePoints; ++i) {
                factor *= (nTruePoints - i)/(double)(nPoints - i);
            }
            if (nPoints != sampleSize) {
                int pow = (sampleSize - (int)nTruePoints);
                factor *= Math.pow((1./(nPoints - sampleSize)), pow);
            }
        }
        
        return factor;
    }
    
    public double calculateTrueSampleProbabilityForMultiplyMatched(
        int nPoints, int nAllMultiplicity, int sampleSize, 
        double expectedFractionTruePoints) {
        
        /* 
        making an assumption of avg multiplicity for all and using the same
        logic as for singly matched points ts start, then adding the multiple
        matching logic:
        
        example:  5 points in set 1 and each is singly matched to a point in set 2.
                  sample size, k, is 2.
                  there are 3 true point 1's {0, 2, 3}
        
        the number of combinations of the singly matched dataset are n!/(k!(n-k)!)
        nCombinations=10. 
        pSample (as seen in CountingTest) is the factors of nTrueMatches/nPoints
        reduced for each subsequent point due to drawing without replacement.
        pSample=(3./5)*(2./4)=0.3 
        The expected number of 'true' samples, that is samples composed of only
        truely matched points is then pSample*nCombinations = 0.3*10 = 3.
        nTrueSamples=3 is verified as seen here:
        0 1
        0 2 * 
        0 3 *
        0 4
        1 2
        1 3
        1 4
        2 3 *
        2 4
        3 4
        
        For multiply matched points, the number of combinations increases,
        but the number of truly matched remains the same.
        
        The number of all point pairs is the multiplicity * nPoints
        (or alternatively, can be given the total number including multiplicity).
      
        Using the sample example, adding a multiplicity of '2' to each point,
        the number of all point pairs is 2 * nPoints = 10.
        the number of possible combinations is 10*9./2 = 45.
        The number of true samples hasn't changed and that is still 3.
        The pSample then must be 3./45. = 0.067.
        so the multiply mapped pSample is pSample_singly_matched*nComb_singlyMatched/nComb_multiplyMatched
          
        in detail:
        
        The number of combinations for the multiply matched is 
            nMultiplyMatched!/(k!(nMultiplyMatched-k)!)
        
        The multiply mapped pSample is
            pSample_singly_matched * nComb_singlyMatched / nComb_multiplyMatched
        
        For the example above,
            pSample_multiply_matched = 0.3 * (10./45.) = 0.067
        
        A better determination would use the real integrated CDF of the
        inhomogeneous distribution, but this number is used to quickly 
        make an estimate.
        
        Using the maximum multiplicity would lead to a safe overestimate,
        so should consider that for another method.
        */
        
        double pSampleSingle = calculateTrueSampleProbability(nPoints, sampleSize, 
            expectedFractionTruePoints);
        
        double nCombinationsSingle = MiscMath.computeNDivKTimesNMinusK(nPoints, sampleSize);
        
        double nCombinationsMultiple = 
            MiscMath.computeNDivKTimesNMinusK(nAllMultiplicity, sampleSize);
        
        double pSampleMultiple = pSampleSingle * 
            (nCombinationsSingle/nCombinationsMultiple);

        return pSampleMultiple;
    }

    
    /*
    Get the number of subsamples needed for a subsample size of 7 within a larger sample,
    where the number of subsamples is chosen sufficiently high to give 
    a probability in excess of 95% that a good subsample is selected.
    The number is calculated for percentages of bad data in the larger sample 
    being 5%, 10%, 20%, 25%, 30%, 40%, and 50%.
    
    The method is from Section 6 of paper: 
      "The Development and Comparison of Robust Method for Estimating the Fundamental Matrix"
       by P H S TORR AND D W MURRAY
       Int. J Computer Vision 24 (1997) pp271–300, , 1–32 ()
       https://www.robots.ox.ac.uk/ActiveVision/Publications/torr_murray_ijcv1997/torr_murray_ijcv1997.pdf
    
    <pre>
    for probablity gamma in excess of 95%:
    let p = sub-sample size
    let m = number of subsamples to achieve gamma probability
    let eps = percentage of bad data
    
    gamma = 1 - ( (1 - (1 - eps)^p)^m)
    
    table of m's for p=7 and eps from 5%to 50%:
    
    p |  5%  10%  20%  25%  30%  40%  50%
    --------------------------------------
    7 |  3   5    13   21   35   106  382
    </pre>
    */
    public static int numberOfSubsamplesOfSize7For95PercentInliers(int outlierPercent) {
    
        if (outlierPercent < 0) {
            throw new IllegalArgumentException("outlierPercent must be non-negative");
        }
        if (outlierPercent > 50) {
            throw new IllegalArgumentException("outlierPercent must be 50 percent or less");
        }
        
        /*
        1 - gamma = (1 - (1 - eps)^p)^m
           let Z = (1 - (1 - eps)^p)
           let g = 1 - gamma
        g = Z^m
        log(g) = m*log(Z);       
        m = log(g) / log(Z)
        */
        int p = 7;
        double z = 1. - Math.pow(1. - ((double)outlierPercent/100.), p);
        double g = 1. - 0.95;
        
        double m = Math.log(g) / Math.log(z);
                
        return (int)Math.ceil(m);
    }
    
    /*
    Get the number of subsamples needed for a subsample size of 7 within a larger sample,
    where the number of subsamples is chosen sufficiently high to give 
    a probability in excess of 95% that a good subsample is selected.
    The number is calculated for percentages of bad data in the larger sample 
    being 5%, 10%, 20%, 25%, 30%, 40%, and 50%.
    
    The method is from Section 6 of paper: 
      "The Development and Comparison of Robust Method for Estimating the Fundamental Matrix"
       by P H S TORR AND D W MURRAY
       Int. J Computer Vision 24 (1997) pp271–300, , 1–32 ()
       https://www.robots.ox.ac.uk/ActiveVision/Publications/torr_murray_ijcv1997/torr_murray_ijcv1997.pdf
    
    <pre>
    for probablity gamma in excess of 95%:
    let p = sub-sample size
    let m = number of subsamples to achieve gamma probability
    let eps = percentage of bad data
    
    gamma = 1 - ( (1 - (1 - eps)^p)^m)
    
    table of m's for p=7 and eps from 5%to 50%:
    
    p |  5%  10%  20%  25%  30%  40%  50%
    --------------------------------------
    7 |  3   5    13   21   35   106  382
    </pre>
    */
    public static int numberOfSubsamplesFor95PercentInliers(int outlierPercent,
        int subSampleSize) {
    
        if (outlierPercent < 0) {
            throw new IllegalArgumentException("outlierPercent must be non-negative");
        }
        if (outlierPercent > 50) {
            throw new IllegalArgumentException("outlierPercent must be 50 percent or less");
        }
        if (subSampleSize < 0) {
            throw new IllegalArgumentException("outlierPercent must be non-negative");
        }
        
        /*
        1 - gamma = (1 - (1 - eps)^p)^m
           let Z = (1 - (1 - eps)^p)
           let g = 1 - gamma
        g = Z^m
        log(g) = m*log(Z);       
        m = log(g) / log(Z)
        */
        int p = subSampleSize;
        double z = 1. - Math.pow(1. - ((double)outlierPercent/100.), p);
        double g = 1. - 0.95;
        
        double m = Math.log(g) / Math.log(z);
                
        return (int)Math.ceil(m);
    }
    
}
