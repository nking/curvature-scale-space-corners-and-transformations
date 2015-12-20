package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.math.BigInteger;

/**
 *
 * @author nichole
 */
public class PointValueDistr {

    private final BigInteger maxValue;
    private final PairInt[] points;
    private final BigInteger[] cumulativeValues;
    
    public PointValueDistr(BigInteger maxValue, PairInt[] points,
        BigInteger[] cumulativeValues) {
        this.maxValue = maxValue;
        this.points = points;
        this.cumulativeValues = cumulativeValues;
    }

    public BigInteger getMaxValue() {
        return maxValue;
    }

    public PairInt[] getPoints() {
        return points;
    }

    public BigInteger[] getCumulativeValues() {
        return cumulativeValues;
    }
    
    public PairInt getForCumulativeValue(BigInteger value) {
        
        // binary search for value
        
        int n = cumulativeValues.length;
           
        int lowIdx = 0;
        int highIdx = n - 1;
        int midIdx;
        
        while (lowIdx != highIdx) {

            midIdx = (highIdx + lowIdx) >> 1;
            
            BigInteger v = cumulativeValues[midIdx];
            
            //-1, 0 or 1 when v is less than, equal to, or greater than value.
            int comp = v.compareTo(value);

            if (comp > 0) {

                if (highIdx == midIdx) {
                    highIdx--;
                } else {
                    highIdx = midIdx;
                }

            } else if (comp < 0) {

                if (lowIdx == midIdx) {
                    lowIdx++;
                } else {
                    lowIdx = midIdx;
                }

            } else {
                // is equal
                return points[midIdx];
            }
        }

        // should never arrive here
        throw new IllegalStateException("Error in algorithm.  searching for "
            + value.toString());
    }
    
}
