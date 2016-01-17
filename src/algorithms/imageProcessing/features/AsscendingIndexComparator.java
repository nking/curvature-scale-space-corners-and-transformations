package algorithms.imageProcessing.features;

import java.io.Serializable;
import java.util.Comparator;

/**
 *
 * @author nichole
 */
public class AsscendingIndexComparator implements 
    Comparator<BlobPerimeterRegion>, Serializable {

    private static final long serialVersionUID = 1;
    
    @Override
    public int compare(BlobPerimeterRegion o1, BlobPerimeterRegion o2) {
        
        if (o1.getIndexWithinCurve() < o2.getIndexWithinCurve()) {
            return -1;
        } else if (o1.getIndexWithinCurve() > o2.getIndexWithinCurve()) {
            return 1;
        }
        return 0;
    }
    
}
