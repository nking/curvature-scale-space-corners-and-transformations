package algorithms.util;

import java.io.Serializable;
import java.util.Comparator;

/**
 * comparator to be used for sorting by decreasing number of points
 * @author nichole
 */
public class PairIntArrayDescendingComparator implements Comparator<PairIntArray>,
    Serializable {

    private static final long serialVersionUID = 1;
    
    @Override
    public int compare(PairIntArray o1, PairIntArray o2) {
        return Integer.compare(o2.getN(), o1.getN());
    }
    
}
