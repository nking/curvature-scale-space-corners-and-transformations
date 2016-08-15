package algorithms.imageProcessing.matching;

import java.util.Comparator;

/**
 * comparator for ascending sort startIdx1,
 * then startIdx2
 */    
public class SequenceComparator4 implements
    Comparator<Sequence> {

    @Override
    public int compare(Sequence o1, Sequence o2) {
        
        if (o1.startIdx1 < o2.startIdx1) {
            return -1;
        } else if (o1.startIdx1 > o2.startIdx1) {
            return 1;
        }
        
        if (o1.startIdx2 < o2.startIdx2) {
                return -1;
        } else if (o1.startIdx2 > o2.startIdx2) {
                return 1;
        }
        
        return 0;
    }   
}
