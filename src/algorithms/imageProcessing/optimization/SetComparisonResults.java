package algorithms.imageProcessing.optimization;

import java.util.List;

/**
 *
 * @author nichole
 */
public class SetComparisonResults implements Comparable<SetComparisonResults> {

    protected final int nExpectedPoints;
    protected final float numberOverrunDivExpected;
    protected final float numberMatchedDivExpected;
    
    protected final int nExpectedBorderPoints;
    protected final float numberMatchedBorderDivExpected;

    public SetComparisonResults(int nExpected, int nOverrun, int nMatched,
        int nExpectedBorder, int nMatchedBorder) {
        
        this.nExpectedPoints = nExpected;
        this.numberOverrunDivExpected = (nExpected > 0) ?
            (float) nOverrun / (float) nExpected : Float.POSITIVE_INFINITY;
        this.numberMatchedDivExpected = (nExpected > 0) ? 
            (float) nMatched / (float) nExpected : Float.POSITIVE_INFINITY;
        this.nExpectedBorderPoints = nExpectedBorder;
        this.numberMatchedBorderDivExpected = (nExpected > 0) ? 
            (float) nMatchedBorder / (float) nExpectedBorder 
            : Float.POSITIVE_INFINITY;
    }
    
    public SetComparisonResults(List<SetComparisonResults> results) {
        
        int nTotal = 0;
        for (SetComparisonResults result : results) {
            nTotal += result.nExpectedPoints;
        }
        this.nExpectedPoints = nTotal;
        
        float sumOverrunFraction = 0;
        for (SetComparisonResults result : results) {
            double fraction = result.numberOverrunDivExpected;
            sumOverrunFraction += fraction;
        }
        this.numberOverrunDivExpected = sumOverrunFraction/(float)results.size();
        
        float sumnMatchedFraction = 0;
        for (SetComparisonResults result : results) {
            double fraction = result.numberMatchedDivExpected;
            sumnMatchedFraction += fraction;
        }
        this.numberMatchedDivExpected = sumnMatchedFraction/(float)results.size();
        
        nTotal = 0;
        for (SetComparisonResults result : results) {
            nTotal += result.nExpectedBorderPoints;
        }
        this.nExpectedBorderPoints = nTotal;
        sumnMatchedFraction = 0;
        for (SetComparisonResults result : results) {
            double fraction = result.numberMatchedBorderDivExpected;
            sumnMatchedFraction += fraction;
        }
        this.numberMatchedBorderDivExpected = sumnMatchedFraction/(float)results.size();
    }
    
    public SetComparisonResults(List<Integer> nExpected, List<Integer> nOverrun, 
        List<Integer> nMatched, List<Integer> nExpectedBorder, 
        List<Integer> nMatchedBorder) {
        
        if (nExpected.size() != nOverrun.size() 
            || nOverrun.size() != nMatched.size() 
            || nMatched.size() != nExpectedBorder.size()
            || nExpectedBorder.size() != nMatchedBorder.size()) {
            throw new IllegalArgumentException("lists must be same size");
        }
        
        int nTotal = 0;
        for (Integer n : nExpected) {
            nTotal += n.intValue();
        }
        this.nExpectedPoints = nTotal;
        
        float sumOverrunFraction = 0;
        for (int i = 0; i < nOverrun.size(); i++) {
            double fraction = (float)nOverrun.get(i)/(float)nExpected.get(i);
            sumOverrunFraction += fraction;
        }
        this.numberOverrunDivExpected = sumOverrunFraction/(float)nOverrun.size();
        
        float sumnMatchedFraction = 0;
        for (int i = 0; i < nMatched.size(); i++) {
            double fraction = (float)nMatched.get(i)/(float)nExpected.get(i);
            sumnMatchedFraction += fraction;
        }
        this.numberMatchedDivExpected = sumnMatchedFraction/(float)nMatched.size();
        
        nTotal = 0;
        for (Integer n : nExpectedBorder) {
            nTotal += n.intValue();
        }
        this.nExpectedBorderPoints = nTotal;
        
        sumnMatchedFraction = 0;
        for (int i = 0; i < nMatchedBorder.size(); i++) {
            double fraction = (float)nMatchedBorder.get(i)/(float)nExpectedBorder.get(i);
            sumnMatchedFraction += fraction;
        }
        this.numberMatchedBorderDivExpected = 
            sumnMatchedFraction/(float)nMatchedBorder.size();
    }

    @Override
    public boolean equals(Object obj) {

        if (!(obj instanceof SetComparisonResults)) {
            return false;
        }
        SetComparisonResults other = (SetComparisonResults) obj;
        
        double eps0 = 0.99 * (1./(double)nExpectedPoints);
        double eps1 = 5.0E-4;//0.536e-7; // 1 pixel in 1024x104:
        
        if (Math.abs(numberOverrunDivExpected - other.numberOverrunDivExpected)
            > eps0) {
            return false;
        }
        if (Math.abs(numberMatchedBorderDivExpected - other.numberMatchedBorderDivExpected)
            > eps1) {
            return false;
        }
        if (Math.abs(numberMatchedDivExpected - other.numberMatchedDivExpected)
            > eps1) {
            return false;
        }
        return true;
    }

    /**
     * Compare this instance's numberOverrunDivExpected, 
     * numberMatchedBorderDivExpected, numberMatchedDivExpected to other and 
     * return -1 if this instance has better fields as a result, else 0 if they 
     * are equal, else +1 if other has better results. The sign convention is 
     * meant to put the best answers at the top of a list sorted using this 
     * comparison function, where the default behavior by java framework 
     * algorithms is to make an ascending sort. 
     * Best is defined as having the smallest numberOverrunDivExpected and
     * ties are broken by having the largest numberMatchedBorderDivExpected,
     * else if there are no overruns and no border points matched, the
     * instance with the largest numberMatchedDivExpected is the best.
     *
     * @param other
     * @return
     */
    @Override
    public int compareTo(SetComparisonResults other) {

        // cannot have even 1 pixel overrun, so eps0 should be 0,
        // but will consider precision  
        double eps0 = 0.99 * (1./(double)nExpectedPoints);
        double eps1 = 5.0E-4;//0.536e-7; // 1 pixel in 1024x104:
        
        float diff0 = Math.abs(numberOverrunDivExpected
            - other.numberOverrunDivExpected);

        if (diff0 <= eps0) {

            float diff1 = Math.abs(numberMatchedBorderDivExpected
                - other.numberMatchedBorderDivExpected);
            
            if (diff1 <= eps1) {
                
                float diff2 = Math.abs(numberMatchedDivExpected
                    - other.numberMatchedDivExpected);

                if (diff2 <= eps1) {
                    return 0;
                } else {
                    if (numberMatchedDivExpected < other.numberMatchedDivExpected) {
                        return 1;
                    } else {
                        return -1;
                    }
                }
            } else {
                if (numberMatchedBorderDivExpected < other.numberMatchedBorderDivExpected) {
                    return 1;
                } else {
                    return -1;
                }
            }
            
        } else {
            if (numberOverrunDivExpected < other.numberOverrunDivExpected) {
                return -1;
            } else {
                return 1;
            }
        }
    }
}
