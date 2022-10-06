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
    equation to estimate the number of subsamples required
    for the probablity gamma to be in excess of 95%, that is,
    the probability that all the data points selected in one subsample are non-outliers.
    
    let p = sub-sample size (e.g. 7 for Fundamental Matrix correspondence)
    let eps = fraction of contaminated data
    let m = number of subsamples to achieve gamma probability
    
    gamma = 1 - ( (1 - (1 - eps)^p)^m)
    
    m = Math.log(1. - 0.95)/Math.log(1. - Math.pow(1. - (outlierPercent/100.), p))
    
    table of m's for p=7 and eps from 5%to 50%:
    
    p |  5%  10%  20%  25%  30%  40%  50%
    --------------------------------------
    7 |  3   5    13   21   35   106  382
    </pre>
    */
    public static int numberOfSubsamplesOfSize7For95PercentInliers(int outlierPercent) {
        return numberOfSubsamplesFor95PercentInliers(outlierPercent, 7);
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
    8 |  3   6    17   29   51   177  766
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
