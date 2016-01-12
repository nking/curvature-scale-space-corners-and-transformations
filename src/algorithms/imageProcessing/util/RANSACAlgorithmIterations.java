package algorithms.imageProcessing.util;

/**
 * estimates the number of iterations for a RANSAC type algorithm that 
 * should be used to ensure with 95% certainty that a set of points is 
 * findable as a subset of n matchable points.

nPoints is the number of points that contain true and false matches.
        
nTruePoints is the number of true matches within nPoints.
        
k is the number of points to draw at one time for a RANSAC sample test iteration.

The number of combinations of samples of size k is nPoints!/(k!*(nPoints-k)!).

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

let p be the fraction just calculated.

knowing p for a sample now, need to calculate the number of samples that
need to be drawn in order to have a high confidence that at least one
sample was a true sample.

Naively, that would be the inverse of p.

Now that have the statistics in a form of 2 states, can use binomial statistics
to determine the probability that 1 true sample is chosen from nIter attempts.

  pConfid = p * (1-p)^(nIter-1))
  (nIter-1) = math.log(pConfid/p)/math.log(1-p)
  nIter = (math.log(pConfid/p)/math.log(1-p)) + 1

  for 99 percent confidence, have
     (int)Math.round(Math.log(0.99/factor1)/Math.log(1. - factor1)) + 1

</pre>
 * @author nichole
 */
public class RANSACAlgorithmIterations {

    public int estimateNIterFor99PercentConfidence(
        int nPoints, int sampleSize, double expectedFractionTruePoints) {
        
        double p = calculateTrueSampleProbability(nPoints, sampleSize, 
            expectedFractionTruePoints);
        
        int nIter = -1 * ((int)Math.round(Math.log(0.99/p)/Math.log(1. - p)) + 1);
                    
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