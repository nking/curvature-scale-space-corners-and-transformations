package algorithms.imageProcessing.util;

/**
 * estimates the number of iterations for a RANSAC type algorithm that 
 * should be used to ensure with 95% certainty that a set of points is 
 * findable as a subset of n matchable points.
 * 
 * @author nichole
 */
public class RANSACAlgorithmIterations {
    
    /**
     The estimate assumes that at least 25 percent of the nMatchedPoints
     are outliers and determines the number of iterations needed for a 95%
     certainty that 8 points drawn randomly from nMatchedPoints will all be
     "true" matches, that is, not outliers.
     
     <pre>
     It uses this table generated below following the statistics in the comments
     at the end of this class.
     
                    | nInliers percentage of nTotal
     nMatchedPoints |   50       75       100
     ------------------------------------------------------
       30   |       | 0.001  | 0.05   |  1.0
       50   |       | 0.002  | 0.07   |  1.0
      100   |       | 0.003  | 0.09   |  1.0
      500   |       | 0.004  | 0.1    |  1.0
     1000   |       | 0.004  | 0.1    |  1.0
     </pre>
     
     * @param nMatchedPoints
     * @return 
     */
    public int estimateNIterForTwentyFivePercentOutliersFor8Points(int nMatchedPoints) {
        
        double p;
        if (nMatchedPoints <= 30) {
            p = 0.05;
        } else if (nMatchedPoints <= 50) {
            p = 0.07;
        } else if (nMatchedPoints <= 100) {
            p = 0.09;
        } else  {
            p = 0.1;
        }
        
        int nIter = (int)Math.round(1./p);
        
        return nIter;
    }
    
    /**
     The estimate assumes that at least 50 percent of the nMatchedPoints
     are outliers and determines the number of iterations needed for a 95%
     certainty that 8 points drawn randomly from nMatchedPoints will all be
     "true" matches, that is, not outliers.
     
     <pre>
     It uses this table generated below following the statistics in the comments
     at the end of this class.
     
                    | nInliers percentage of nTotal
     nMatchedPoints |   50       75       100
     ------------------------------------------------------
       30   |       | 0.001  | 0.05   |  1.0
       50   |       | 0.002  | 0.07   |  1.0
      100   |       | 0.003  | 0.09   |  1.0
      500   |       | 0.004  | 0.1    |  1.0
     1000   |       | 0.004  | 0.1    |  1.0
     </pre>
     
     * @param nMatchedPoints
     * @return 
     */
    public int estimateNIterForFiftyPercentOutliersFor8Points(int nMatchedPoints) {
        
        double p;
        if (nMatchedPoints <= 30) {
            p = 0.001;
        } else if (nMatchedPoints <= 50) {
            p = 0.002;
        } else if (nMatchedPoints <= 100) {
            p = 0.003;
        } else  {
            p = 0.004;
        }
        
        int nIter = (int)Math.round(1./p);
        
        return nIter;
    }

    /**
     The estimate assumes that at least 15 percent of the nPoints
     are true matches with points in another set of points
     and determines the number of iterations needed for a 95%
     certainty that 2 points drawn randomly from nPoints will all be
     "truly matchable" with points in set2 (whose characteristics are
     presumed to be similar).
     
     <pre>
     It uses this table generated below following the statistics in the comments
     at the end of this class.
     
             | nInliers percentage of nTotal
     nPoints |  15        25       50       75
     ------------------------------------------------------
       30    | 0.014  | 0.05   | 0.24   |  0.53
       50    | 0.017  | 0.054  | 0.245  |  0.54
      100    | 0.021  | 0.061  | 0.247  |  0.56
      177    | 0.021  | 0.061  | 0.246  |  0.56
      500    | 0.022  | 0.062  | 0.25   |  0.56
     1000    | 0.022  | 0.062  | 0.25   |  0.56
     </pre>
     
     * @param nPoints
     * @return 
     */
    public int estimateNIterForFifteenPercentOutliersFor2Points(int nPoints) {
        
        double p;
        if (nPoints <= 30) {
            p = 0.014;
        } else if (nPoints <= 50) {
            p = 0.017;
        } else if (nPoints <= 100) {
            p = 0.021;
        } else  {
            p = 0.022;
        }
        
        int nIter = (int)Math.round(1./p);
        
        return nIter;
    }

}

/*
There are nTotal matches, that is nTotal sets of points thought to be matched.

nInliers is the number of points within nTotal that are within tolerance
of a good fit (which is defined as a 'true' match).

p = (nInliers/nTotal) inliers
1 - p = 1 - (nInliers/nTotal) outliers or (nOutliers/nTotal)

subset size is 8

8 points are drawn randomly from nTotal points.

probability of drawing 8 'true' inliers:

                             C(8 inliers of nInliers)
     p = 8 inliers =  -------------------------------------   
                       C(8 points drawn from nTotal points)

                      where C is C(n, k) = n! / ( k!*(n - k)! )

                       ( nInliers! / ( 8!*(nInliers - 8)! )
                   =  -------------------------------------------
                       ( nTotal! / ( 8!*(nTotal - 8)! )

what number of throws will result in 95% certainty that a subset of 8 contained 
all 'true' inliers?

(Bernoulli process:)
     throw nRequired number of times. 
     p above gives probability of success for a throw.
     we want to know the total probability that results in one of the 
     nRequired being a success (all inliers).

     P(k) = (n!/(k!(n-k)!)) * p^k * (1 - p)^(N-k)

     P(1) = (nRequired!/(1!(nRequired-1)!)) * p^1 * (1 - p)^(nRequired-1)

     PDesired = (nRequired!/(1!(nRequired-1)!)) * p^1 * (1 - p)^(nRequired-1)
                 nRequired!
              = -------------- * p * (1 - p)^(nRequired - 1)
                (nRequired-1)!

              = nRequired * p * (1 - p)^(nRequired - 1)

     set PDesired = 0.95;  estimate p for nInliers = 0%, 25%, 50%, 75%, 100%
     and estimate nRequired from those.

used the python scripts below for a quick fill of this table:

      p:
            | nInliers percentage of nTotal
     nToTal |  25       50       75       100
     ------------------------------------------------------
       30   |  ~0    | 0.001  | 0.05   |  1.0
       50   |  ~0    | 0.002  | 0.07   |  1.0
      100   |  ~0    | 0.003  | 0.09   |  1.0
      500   |  ~0    | 0.004  | 0.1    |  1.0
     1000   |  ~0    | 0.004  | 0.1    |  1.0

for PDesired = 0.95, nTotal=30, p=75%(=0.05)
     0.95 = nRequired * p * (1 - p)^(nRequired - 1)
     0.95 = nRequired * 0.05 * (0.95)^(nRequired - 1)
     19 = nRequired * (0.95)^(nRequired - 1)
        = nRequired * math.pow(0.95, (nRequired - 1))
        = nRequired * 0.95

====> nRequired = 1/p <====

python code for the table above:

def factorial(n) :
    result = 1;
    n = int(n)
    for i in range(1, n+1):
        result *= i;
    return result;

def factorialNDivNMinusK(n, k) :
    result = 1;
    n = int(n)
    for i in range((n-k)+1, n+1):
        result *= i;
    return result;

def tbl() :
    nTotal = [30, 50, 100, 500, 1000];
    nInFrac = [.25, .50, .75, 1.0];
    fracLen = len(nInFrac)
    for nT in nTotal:
        for ii in range(0, fracLen):
            nF = nInFrac[ii]
            nI = nT * nF
            f0 = factorialNDivNMinusK(nI, 8) * 1.0
            f2 = factorialNDivNMinusK(nT, 8) * 1.0
            f8 = factorial(8) * 1.0
            numer = f0 / f8
            denom = f2 / f8
            p = numer / denom
            print "nTotal=", nT, " nInliers fraction=", nF, " p=", p
    return;

==============

The same logic as above, used for finding 2 truly matchable points within set1:

      p:
            | nInliers percentage of nTotal
     nToTal |  15        25       50       75
     ------------------------------------------------------
       30   | 0.014  | 0.05   | 0.24   |  0.53
       50   | 0.017  | 0.054  | 0.245  |  0.54
      100   | 0.021  | 0.061  | 0.247  |  0.56
      177   | 0.021  | 0.061  | 0.246  |  0.56
      500   | 0.022  | 0.062  | 0.25   |  0.56
     1000   | 0.022  | 0.062  | 0.25   |  0.56

*/