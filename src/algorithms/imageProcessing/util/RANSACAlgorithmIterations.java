package algorithms.imageProcessing.util;

/**
 * estimates the number of iterations for a RANSAC type algorithm that 
 * should be used to ensure with 99% certainty that a set of points is 
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
        
        double pSample = calculateTrueSampleProbability(nPoints, sampleSize, 
            expectedFractionTruePoints);
        
        long nIter = Math.round(Math.log(1. - 0.99) / Math.log(1. - pSample)) + 1;
                
        return nIter;
    }
    
    public double calculateTrueSampleProbability(int nPoints, int sampleSize,
        double expectedFractionTruePoints) {
        
        double nTruePoints = nPoints * expectedFractionTruePoints;
        
        double factor = 1;
        for (int i = 0; i < sampleSize; ++i) {
            factor *= (double)(nTruePoints - i)/(double)(nPoints - i);
        }
        
        return factor;
    }

}