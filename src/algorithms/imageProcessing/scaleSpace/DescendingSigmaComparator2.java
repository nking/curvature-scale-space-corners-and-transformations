package algorithms.imageProcessing.scaleSpace;

import java.io.Serializable;
import java.util.Comparator;
import java.util.List;

/**
 * reorders the lists by the largest sigma in the first item.  Makes
 * an assumption that each item List has already been sorted by descending
 * sigma too.
 * @author nichole
 */
public class DescendingSigmaComparator2 implements 
    Comparator<List<CurvatureScaleSpaceContour>>, Serializable {

    private static final long serialVersionUID = 1;
    
    @Override
    public int compare(List<CurvatureScaleSpaceContour> list1, 
        List<CurvatureScaleSpaceContour> list2) {
        
        if (!list1.isEmpty() && list2.isEmpty()) {
            return -1;
        } else if (list1.isEmpty() && !list2.isEmpty()) {
            return 1;
        } else if (list1.isEmpty() && list2.isEmpty()) {
            return 0;
        }
        
        CurvatureScaleSpaceContour o1 = list1.get(0);
        CurvatureScaleSpaceContour o2 = list2.get(0);
        
        if (o1.getPeakSigma() > o2.getPeakSigma()) {
            return -1;
        } else if (o1.getPeakSigma() < o2.getPeakSigma()) {
            return 1;
        }
        if (o1.getPeakScaleFreeLength() < o2.getPeakScaleFreeLength()) {
            return -1;
        } else if (o1.getPeakScaleFreeLength() > o2.getPeakScaleFreeLength()) {
            return 1;
        } 
        return 0;
    }
    
}
