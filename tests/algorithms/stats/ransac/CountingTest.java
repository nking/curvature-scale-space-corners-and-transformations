package algorithms.stats.ransac;

import algorithms.SubsetChooser;
import algorithms.imageProcessing.util.RANSACAlgorithmIterations;
import algorithms.misc.MiscMath;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

public class CountingTest extends TestCase {

    public void test0() throws Exception {
        
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
        
        System.out.println("nPoints=" + nPoints + " k=" + k + " ratio=" + ratio);
        System.out.println("(nTrue/nPoints)*((nTrue-1)/(nPoints-1))**((nTrue-2)/(nPoints-2))..." + pSample);
        System.out.println("nComb=" + nComb.toString());
        System.out.println("nSamplesTrue=" + nSamplesTrue.toString());
        
        BigDecimal nExpectedTrueSamples = new BigDecimal(nComb);
        BigDecimal nExpectedTrueSamples1 = nExpectedTrueSamples.multiply(new BigDecimal(Double.toString(pSample)));
        
        System.out.println("nExpectedTrueSamples1=" + nExpectedTrueSamples1.toString());
        
        BigDecimal nSamplesTrueD = new BigDecimal(nSamplesTrue);
        BigDecimal expectedDivFound1 = nExpectedTrueSamples1.divide(nSamplesTrueD, RoundingMode.HALF_UP);
        System.out.println("true samples expected/found="
            + " " + expectedDivFound1.toString());
        
        assertTrue(Math.abs(expectedDivFound1.doubleValue() - 1.0) < 0.05);
        
        long nCombinations = MiscMath.computeNDivNMinusK(nPoints, k)/MiscMath.factorial(k);
        
        assertEquals(nCombinations, nComb.longValue());
        
        int nIterNaive = (int)Math.round(1./pSample);
       
        int nIter99PercentConfidence = -1*((int)Math.round(Math.log(0.99/pSample)/Math.log(1. - pSample)) + 1);
        
        RANSACAlgorithmIterations nIterC = new RANSACAlgorithmIterations();
        long nIter99 = nIterC.estimateNIterFor99PercentConfidence(nPoints, k, ratio);
                
        System.out.println("nIterNaive=" + nIterNaive + " nIter99=" + nIter99
            + " nIter99PercentConfidence=" + nIter99PercentConfidence);
    }
    
    public void est2() {
        
        RANSACAlgorithmIterations nIterC = new RANSACAlgorithmIterations();
        
        double ratio = 30./50.;
        int k = 7;
        int nPoints = 50;
        
        double factor1 = nIterC.calculateTrueSampleProbability(nPoints, k, ratio);
        
        System.out.println("pSample=" + factor1);
                
        int nIter = 50;
        double prob = Math.pow(1. - factor1, nIter);
        System.out.println("nIter=" + nIter + " prob=" + (1.-prob));
        
        nIter = 100;
        prob =  Math.pow(1. - factor1, nIter);
        System.out.println("nIter=" + nIter + " prob=" + (1.-prob));
        
        nIter = 700;
        prob = Math.pow(1. - factor1, nIter);
        System.out.println("nIter=" + nIter + " prob=" + (1.-prob));
        
    }

}
