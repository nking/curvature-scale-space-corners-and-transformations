package algorithms.imageProcessing.matching;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 *
 * @author nichole
 */
public class Sequences {
 
    private List<Sequence> sequences = new ArrayList<Sequence>();
    private float fractionOfWhole;
    private double absSumDiffs;
    private float avgSumDiffs;

    /**
     * NOTE: method has side effect of sorting the
     * sequences by startIdx1.
     * @return 
     */
    public boolean isConsistentClockwise() {
        return isConsistentClockwise(sequences);
    }
    
    /**
     * NOTE: method has side effect of sorting sequences
     * by startIdx1.
     * @param sequences
     * @return 
     */
    public static boolean isConsistentClockwise(
        List<Sequence> sequences) {

        int ns = sequences.size();
        
        if (ns == 0) {
            return true;
        }

        Collections.sort(sequences, new SequenceComparator4());

        Sequence sPrev = sequences.get(0);
        
        int n1 = sPrev.getN1();
        
        for (int i = 1; i < ns; ++i) {
            
            Sequence sTest = sequences.get(i);
            
            int sTestStopIdx1 = sTest.startIdx1 + (sTest.stopIdx2 - sTest.startIdx2);
         
            int stopIdx1 = sPrev.startIdx1 + (sPrev.stopIdx2 - sPrev.startIdx2);
        
            // aggregation does not extend idx2 past n2-1
            if (sTestStopIdx1 < n1 && stopIdx1 < n1) {
                if (!(sTest.startIdx1 > stopIdx1 &&
                    sTest.startIdx2 > sPrev.stopIdx2)) {
                    return false;
                }
            }
                        
            if (sTestStopIdx1 > (n1 - 1)) {
                sTestStopIdx1 -= n1;
                if (stopIdx1 > (n1 - 1)) {
                    stopIdx1 -= n1;
                    //sPrev: 0-->lstPt frstPt   n-1   (lstPt overrun)
                    //sTest:   0-->lstPt frstPt   n-1   (lstPt overrun)
                    if (stopIdx1 >= sTest.startIdx1 ||
                        sPrev.stopIdx2 >= sTest.stopIdx2) {
                        return false;
                    }
                } else {
                    //sPrev:           frstPt  lstPt   n-1
                    //sTest: 0-->lstPt              frstPt   n-1   (lastPoint overrun)
                    /* for example:
                    n1=10, n2=15
                    4:4  6           (4 to 6) 
                    9:9  13  (9 to 3,                9 to 13
                    */
                    if (!(sTest.startIdx1 > stopIdx1 &&
                        sTestStopIdx1 < sPrev.startIdx1 &&
                        sTest.startIdx2 > sPrev.stopIdx2)) {
                        return false;
                    }
                }
            } else if (stopIdx1 > (n1 - 1)) {
                stopIdx1 -= n1;
                //sPrev: 0-->lstPt               frstPt   n-1   (lastPoint overrun)
                //sTest:           frstPt  lstPt     
                // this one shouldn't happen due to the sort
                // by startIdx1
                throw new IllegalStateException("ERROR in algorithm"
                    + " check sorted values: \n sPrev=" + sPrev + 
                    " \n sTest=" + sTest);
            }
                
            sPrev = sTest;
        }
        
        return true;        
    }
    
    public boolean add(Sequence s) {
        return sequences.add(s);
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("frac=%.4f",
            fractionOfWhole));
        sb.append(String.format(", avgDiff=%.4f,  sumDiff=%.4f",
            avgSumDiffs, absSumDiffs));
        for (Sequence s : sequences) {
            sb.append("\n").append(s.toString());
        }
        return sb.toString();
    }

    /**
     * @return the sequences
     */
    public List<Sequence> getSequences() {
        return sequences;
    }

    /**
     * @return the fractionOfWhole
     */
    public float getFractionOfWhole() {
        return fractionOfWhole;
    }

    /**
     * @return the absSumDiffs
     */
    public double getAbsSumDiffs() {
        return absSumDiffs;
    }

    /**
     * @return the avgSumDiffs
     */
    public float getAvgSumDiffs() {
        return avgSumDiffs;
    }

    /**
     * @param fractionOfWhole the fractionOfWhole to set
     */
    public void setFractionOfWhole(float fractionOfWhole) {
        this.fractionOfWhole = fractionOfWhole;
    }

    /**
     * @param absSumDiffs the absSumDiffs to set
     */
    public void setAbsSumDiffs(double absSumDiffs) {
        this.absSumDiffs = absSumDiffs;
    }

    /**
     * @param avgSumDiffs the avgSumDiffs to set
     */
    public void setAvgSumDiffs(float avgSumDiffs) {
        this.avgSumDiffs = avgSumDiffs;
    }
}
