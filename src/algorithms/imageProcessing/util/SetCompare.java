package algorithms.imageProcessing.util;

import algorithms.util.PairInt;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class SetCompare {
 
    public SetComparisonResults compare(Set<PairInt> expectedPoints,
        Set<PairInt> points) {
        
        int nFound = 0;
        int nOutside = 0;
        
        for (PairInt p : points) {
            if (expectedPoints.contains(p)) {
                nFound++;
            } else {
                nOutside++;
            }
        }

        SetComparisonResults r = new SetComparisonResults(expectedPoints.size(),
            nOutside, nFound);
        
        return r;
        
    }

}
