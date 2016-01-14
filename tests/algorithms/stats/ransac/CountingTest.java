package algorithms.stats.ransac;

import algorithms.SubsetChooser;
import algorithms.imageProcessing.util.RANSACAlgorithmIterations;
import algorithms.misc.MiscMath;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.security.SecureRandom;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;

public class CountingTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public void est0() throws Exception {
        
        /*
        a quick test of expected number of samples composed of all true matches
        when samples are drawn from the universe of all possible combinations.
        
        nPoints is the number of points that contain true and false matches.
        nTruePoints is the number of true matches within nPoints.
        k is the number of points to draw at one time for a RANSAC sample test iteration.
        
        The number of combinations of samples of size k is 
            nPoints!/(k!*(nPoints-k)!).
        For a sample, the first draw of a matching point out of nPoints has 
        possibility of being all 'true' points = (nTruePoints/nPoints).
        For the same sample, the second draw in the sample has the possibility 
        of being all 'true' points = (nTruePoints - 1)/(nPoints - 1).
        etc.
        
        The total number of 'true' samples within all possible combinations of
        size k from nPoints is then
        nCombinations * ((nTruePoints/nPoints)*(nTruePoints - 1)/(nPoints - 1)...((nTruePoints - k-1)/(nPoints - k-1))
        
        This method counts the number of samples of all 'true' points and
        compares with the expected number.
        */
        
        int nPoints = 50;
        
        Set<Integer> knownTruePoints = new HashSet<Integer>(30);
        for (int i = 0; i < 25; ++i) {
            knownTruePoints.add(Integer.valueOf(i));
        }
        knownTruePoints.add(Integer.valueOf(30));
        knownTruePoints.add(Integer.valueOf(32));
        knownTruePoints.add(Integer.valueOf(40));
        knownTruePoints.add(Integer.valueOf(42));
        knownTruePoints.add(Integer.valueOf(43));
        
        BigInteger nSamplesTrue = BigInteger.ZERO;
        BigInteger nComb = BigInteger.ZERO;
        
        int k = 7;
        
        int[] selectedIndexes = new int[k];
        
        SubsetChooser subsetChooser = new SubsetChooser(nPoints, k);            
        int nV = subsetChooser.getNextSubset(selectedIndexes);
        while (nV != -1) {
            
            boolean allKnownTrue = true;
            for (int bitIndex : selectedIndexes) {
                int idx = bitIndex;
                if (!knownTruePoints.contains(Integer.valueOf(idx))) {
                    allKnownTrue = false;
                    break;
                }
            }
            if (allKnownTrue) {
                nSamplesTrue = nSamplesTrue.add(BigInteger.ONE);
            }
            
            nComb = nComb.add(BigInteger.ONE);
            
            nV = subsetChooser.getNextSubset(selectedIndexes);
        }
        
        //nSamplesTrue ~ (ratio)^(k) * nComb
        // and more specifically, is (nTrue/nPoints)*((nTrue-1)/(nPoints-1))**((nTrue-2)/(nPoints-2))...
        double ratio = (double)knownTruePoints.size()/(double)nPoints;
        double pSample = 1;
        for (int i = 0; i < k; ++i) {
            pSample *= ((double)knownTruePoints.size() - i)/(double)(nPoints - i);
        }
        
        log.info("nPoints=" + nPoints + " k=" + k + " ratio=" + ratio);
        log.info("(nTrue/nPoints)*((nTrue-1)/(nPoints-1))**((nTrue-2)/(nPoints-2))..." + pSample);
        log.info("nComb=" + nComb.toString());
        log.info("nSamplesTrue=" + nSamplesTrue.toString());
        
        BigDecimal nExpectedTrueSamples = new BigDecimal(nComb);
        BigDecimal nExpectedTrueSamples1 = nExpectedTrueSamples.multiply(new BigDecimal(Double.toString(pSample)));
        
        log.info("nExpectedTrueSamples1=" + nExpectedTrueSamples1.toString());
        
        BigDecimal nSamplesTrueD = new BigDecimal(nSamplesTrue);
        BigDecimal expectedDivFound1 = nExpectedTrueSamples1.divide(nSamplesTrueD, RoundingMode.HALF_UP);
        log.info("true samples expected/found="
            + " " + expectedDivFound1.toString());
        
        assertTrue(Math.abs(expectedDivFound1.doubleValue() - 1.0) < 0.05);
        
        long nCombinations = MiscMath.computeNDivNMinusK(nPoints, k)/MiscMath.factorial(k);
        
        assertEquals(nCombinations, nComb.longValue());
        
        int nIterNaive = (int)Math.round(1./pSample);
       
        int nIter99PercentConfidence = -1*((int)Math.round(Math.log(0.99/pSample)/Math.log(1. - pSample)) + 1);
        
        RANSACAlgorithmIterations nIterC = new RANSACAlgorithmIterations();
        long nIter99 = nIterC.estimateNIterFor99PercentConfidence(nPoints, k, ratio);
                
        log.info("nIterNaive=" + nIterNaive + " nIter99=" + nIter99
            + " nIter99PercentConfidence=" + nIter99PercentConfidence);
    }
    
    public void est2() {
        
        RANSACAlgorithmIterations nIterC = new RANSACAlgorithmIterations();
        
        double ratio = 30./50.;
        int k = 7;
        int nPoints = 50;
        
        double factor1 = nIterC.calculateTrueSampleProbability(nPoints, k, ratio);
        
        log.info("pSample=" + factor1);
                
        int nIter = 50;
        double prob = Math.pow(1. - factor1, nIter);
        log.info("nIter=" + nIter + " prob=" + (1.-prob));
        
        nIter = 100;
        prob =  Math.pow(1. - factor1, nIter);
        log.info("nIter=" + nIter + " prob=" + (1.-prob));
        
        nIter = 700;
        prob = Math.pow(1. - factor1, nIter);
        log.info("nIter=" + nIter + " prob=" + (1.-prob));
        
    }

    public void test3() throws Exception {
        
        /*
        counting for the multiply matched set of points
        
        example:  5 points in set 1.  each has 2 matches in set 2, only one is a true match.
                  sample size, k, is 2.
                  there are 3 true point 1's {0, 2, 3}
        
        combinations of single matched dataset: 
        ncomb=10. pSample=(3./5)*(2./4)=0.3 
        expected number of 'true' samples is then 0.3*ncomb = 3.
        nTrueSamples=3 is seen below
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
        
        if multiplicity is '2', then the number of all point pairs is 2 * nPoints.
        and then the number of possible combinations is 10*9./2 = 45.
        The number of true samples hasn't changed and that is still 3.
        The pSample then must be 3./45. = 0.067.
        so the multiply mapped pSample is pSample_singly_matched*nComb_singlyMatched/nComb_multiplyMatched
        
        above case considers even distribution of multiply matched numbers.
        a better approx could use the mode or median of multiply matched.
        also make an example for that.
        */
        
        //since the concept is sound, will test over the range of valid
        //  arguments looking for numerical errors 
        
        RANSACAlgorithmIterations nIterC = new RANSACAlgorithmIterations();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        log.info("SEED=" + seed);
        sr.setSeed(seed);
        
        /*
        k must always be <= nPoints.
        
        since k will always be 7 in this context, will not test for large k
        
        therefore, nPoints must be >= 7.  largest expected nPoints would be
        a few hundred for reasonable setting.  the number of combinations
        is expressed as long, so the maximum number of nPoints (if multiplicity is 1)
        would be 
          ((n_max*(n_max-1)*(n_max-2)*(n_max-3)*(n_max-4)*(n_max-5)*(n_max-6))/6300.) < (2^63 - 1)
        so a maximum of 1789
        
        ratio should be > 0 and <= 1
        
        multiplicity isn't expected to be larger than 10 ever (hard wired in
        CosineSimilarityCornerMatcher which I keep editing, but could consider
        a single case of multiplicity of some fraction of nPoints, maybe 0.2.
        */
        
        int nTests = 100;
        
        for (int i = 0; i < nTests; ++i) {
            
            long nIter = 0;
            int k=-1; 
            int nPoints = -1;
            int nMultiplicity = -1; 
            int nAll = -1;;
            double ratio = -1;
            try {
                k = sr.nextInt(10) + 7;

                nPoints = sr.nextInt(300) + k + 1;

                ratio = sr.nextDouble() + 0.0001;

                nMultiplicity = sr.nextInt(9) + 1;

                nAll = nPoints * nMultiplicity;

                nIter = nIterC.estimateNIterFor99PercentConfidenceDegenerate(
                    nPoints, nAll, k, ratio);
                
            } catch (Throwable t) {
                log.info("k=" + k + " nPoints=" + nPoints + " ratio=" + ratio + 
                    " nMultiplicity=" + nMultiplicity + " nAll=" + nAll +
                    " nIter=" + nIter);
                log.severe(t.getMessage());
                t.printStackTrace();
                fail(t.getMessage());
            }
            
            if (nIter < 0) {
                log.info("k=" + k + " nPoints=" + nPoints + " ratio=" + ratio + 
                    " nMultiplicity=" + nMultiplicity + " nAll=" + nAll +
                    " nIter=" + nIter);
            }
            
            assertTrue(nIter > 0);            
        }
        
        int k = 7;
        int nPoints = 1780;
        int nAll = nPoints;
        double ratio = 0.5;
        boolean noErrors = true;
        try {
            long nIter = nIterC.estimateNIterFor99PercentConfidenceDegenerate(
                nPoints, nAll, k, ratio);            
            assertTrue(nIter > 0);
        } catch (Throwable t) {
            noErrors = false;
        }
        assertTrue(noErrors);
        
    }
}
