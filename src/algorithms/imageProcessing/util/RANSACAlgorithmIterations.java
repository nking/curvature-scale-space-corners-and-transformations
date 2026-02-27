package algorithms.imageProcessing.util;

import algorithms.misc.MiscMath;

/**
 * estimates the number of iterations for a RANSAC type algorithm that 
 * should be used to ensure with 99% certainty or other that a set of points is 
 * findable as a subset of n matchable points.
<pre>
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
    
    let m = sub-sample size (e.g. 7 for Fundamental Matrix correspondence)
    let p_err = percent outliers = error rate in the ensemble of all data points
    let N = number of subsamples to achieve gamma probability
    
    gamma = 1 - ( (1 - p_err^m)^N)

    z = 1. - math.pow(1. - p_err, m)
    
    N = Math.log(1. - 0.95)/Math.log(z)
    
    table of N's for m=7 and eps from 5%to 50%:
    
    m |  5%  10%  20%  25%  30%  40%  50%
    --------------------------------------
    7 |  3   5    13   21   35   106  382
    </pre>
    */
    public static int numberOfSubsamplesOfSize7For95PercentInliers(int outlierPercent) {
        return numberOfSubsamples(95., outlierPercent, 7);
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
    for probability gamma in excess of 95%:

    let m = sub-sample size (e.g. 7 for Fundamental Matrix correspondence)
    let p_in = inlier ratio.
      hence p_out = 1 - p_in.
    let p be the desired success probability
    let N = number of subsamples to achieve p probability

    probability that all points in a single sample are inliers
    is (p_in)^m which is (1 - p_out)^m)

    we have (1 - p_out)^m) for the probability of all inliers in a sample,
    so then P(not all inliers) = 1 - (1 - p_out)^m).

    then the probability of N samples of "not all inlier" samples
    is P(all N samples have at least one outlier) = (1 - (1 - p_out)^m)^N.

    The complement is that at least one sample is all inliers
    p = P(in N samples, at least 1 is all inliers) = 1 - ((1 - (1 - p_out)^m)^N)

    p = 1 - ((1 - (1 - p_out)^m)^N)
    1 - p = (1 - (1 - p_out)^m)^N
    log(1 - p) = N * log(1 - (1 - p_out)^m))
    N = log(1-p) / log(1 - (1 - p_out)^m))

 
    table of N's for m=7, 8 and eps from 5%to 50%:
    
    m |  5%  10%  20%  25%  30%  40%  50%
    --------------------------------------
    7 |  3   5    13   21   35   106  382
    8 |  3   6    17   29   51   177  766
    </pre>
    Note that the number of iterations scales roughly with the subSampleSize as
    2 ^ (subSampleSize) for larger subSampleSize, 
    so this method should only be used for small sample sizes (i.e. < 100, preferably << 100).
    */
    public static int numberOfSubsamplesFor95PercentInliers(int outlierPercent,
        int subSampleSize) {

        return numberOfSubsamples(95., outlierPercent, subSampleSize);
    }

    /**
     * calculate the number of iterations needed to ensure with percentCertainty that all members of
     * subSampleSize drawn are "inliers".
       Note that the number of iterations scales roughly with the subSampleSize as
       2 ^ (subSampleSize), 
       so this method should only be used for small sample sizes (i.e. < 100, preferably << 100).
     * @param percentCertainty the percent certainty that all members in sample of size subSampleSize
     *                          are "inliers".
     * @param outlierPercent the percent outliers in the entire dataset.
     * @param subSampleSize the sample size drawn from the entire dataset in iterations.
     * @return the number of iteration needed to ensure w/ %certainty that all
     * members of subSample are inliers.
     */
    public static int numberOfSubsamples(double percentCertainty, int outlierPercent,
                                                            int subSampleSize) {
        if (percentCertainty <= 0 || percentCertainty >= 100.) {
            throw new IllegalArgumentException("percentCertainty must be > 0 and < 100");
        }
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
        let m = sub-sample size (e.g. 7 for Fundamental Matrix correspondence)
        let p_in = inlier ratio.
          hence p_out = 1 - p_in.
        let p be the desired success probability
        let N = number of subsamples to achieve p probability

        probability that all points in a single sample are inliers
        is (p_in)^m which is (1 - p_out)^m)

        we have (1 - p_out)^m) for the probability of all inliers in a sample,
        so then P(not all inliers) = 1 - (1 - p_out)^m).

        then the probability of N samples of "not all inlier" samples
        is P(all N samples have at least one outlier) = (1 - (1 - p_out)^m)^N.

        The complement is that at least one sample is all inliers
        p = P(in N samples, at least 1 is all inliers) = 1 - ((1 - (1 - p_out)^m)^N)

        p = 1 - ((1 - (1 - p_out)^m)^N)
        1 - p = (1 - (1 - p_out)^m)^N
        log(1 - p) = N * log(1 - (1 - p_out)^m))
        N = log(1-p) / log(1 - (1 - p_out)^m))
        */
        int m = subSampleSize;
        double w = (double)outlierPercent/100.;
        double z = 1. - Math.pow(1. - w, m);
        double g = 1. - percentCertainty/100.;

        double n = Math.log(g) / Math.log(z);

        return (int)Math.ceil(n);
    }
}
