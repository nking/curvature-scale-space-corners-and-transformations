package algorithms.imageProcessing.scaleSpace;

import algorithms.imageProcessing.features.CornerRegion;
import java.io.Serializable;
import java.util.Comparator;

/**
 * comparator for absolute value of k, descending sort
 * @author nichole
 */
public class DescendingKComparator implements Comparator<CornerRegion>, 
    Serializable {

    private static final long serialVersionUID = 12345;

    @Override
    public int compare(CornerRegion o1, CornerRegion o2) {
        
        if (o1 == null && o2 == null) {
            return 0;
        } else if (o1 != null && o2 == null) {
            return -1;
        } else if (o1 == null && o2 != null) {
            return 1;
        }
        
        float k1 = Math.abs(o1.getK()[o1.getKMaxIdx()]);
        float k2 = Math.abs(o2.getK()[o2.getKMaxIdx()]);
        
        if (k1 > k2) {
            return -1;
        } else if (k1 < k2) {
            return 1;
        } else {
            return 0;
        }
    }
        
}
