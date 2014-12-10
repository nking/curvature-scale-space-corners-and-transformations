package algorithms.imageProcessing;

import java.io.Serializable;
import java.util.Comparator;

/**
 *
 * @author nichole
 */
public class DescendingSigmaComparator implements 
    Comparator<CurvatureScaleSpaceContour>, Serializable {

    private static final long serialVersionUID = 1;
    
    @Override
    public int compare(CurvatureScaleSpaceContour o1, CurvatureScaleSpaceContour o2) {
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
