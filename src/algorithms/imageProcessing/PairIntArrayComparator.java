package algorithms.imageProcessing;

import java.util.Comparator;

/**
 * comparator to be used for sorting by increasing number of points
 * @author nichole
 */
public class PairIntArrayComparator implements Comparator<PairIntArray> {

    @Override
    public int compare(PairIntArray o1, PairIntArray o2) {
        return Integer.compare(o1.getN(), o2.getN());
    }
    
}
