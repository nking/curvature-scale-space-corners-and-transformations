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

    public SetComparisonResults(int nExpected, int nOverrun, int nMatched) {
        this.nExpectedPoints = nExpected;
        this.numberOverrunDivExpected = (nExpected > 0) ?
            (float) nOverrun / (float) nExpected : Float.POSITIVE_INFINITY;
        this.numberMatchedDivExpected = (nExpected > 0) ? 
            (float) nMatched / (float) nExpected : Float.POSITIVE_INFINITY;
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
    }
    
    public SetComparisonResults(List<Integer> nExpected, List<Integer> nOverrun, 
        List<Integer> nMatched) {
        
        if (nExpected.size() != nOverrun.size() || nOverrun.size() != nMatched.size()) {
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
        
    }

    @Override
    public boolean equals(Object obj) {

        if (!(obj instanceof SetComparisonResults)) {
            return false;
        }
        SetComparisonResults other = (SetComparisonResults) obj;
        
        // 1 pixel in 1024x104:
        double eps = 0;//0.536e-7;
        
        if (Math.abs(numberOverrunDivExpected - other.numberOverrunDivExpected)
            > eps) {
            return false;
        }
        if (Math.abs(numberMatchedDivExpected - other.numberMatchedDivExpected)
            > eps) {
            return false;
        }
        return true;
    }

    /**
     * compare this instance's numberOverrunDivExpected and
     * numberMatchedDivExpected to other and return -1 if this instance has
     * better fields as a result, else 0 if they are equal, else +1 if other has
     * better results. The sign convention is meant to put the best answers at
     * the top of a list sorted using this comparison function, where the
     * default behavior by java framework algorithms is to make an ascending
     * sort. Best is defined as having the smallest numberOverrunDivExpected and
     * ties are broken by having the largest numberMatchedDivExpected.
     *
     * @param other
     * @return
     */
    @Override
    public int compareTo(SetComparisonResults other) {

        // 1 pixel in 1024x104:
        double eps = 0;//0.536e-7;

        float diff0 = Math.abs(numberOverrunDivExpected
            - other.numberOverrunDivExpected);

        if (diff0 <= eps) {

            float diff1 = Math.abs(numberMatchedDivExpected
                - other.numberMatchedDivExpected);

            if (diff1 <= eps) {
                return 0;
            } else {
                if (numberMatchedDivExpected < other.numberMatchedDivExpected) {
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
