package algorithms.imageProcessing.optimization;

import algorithms.util.PairInt;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class SetCompare {
 
    public SetComparisonResults compare(Set<PairInt> expectedPoints,
        Set<PairInt> points, Set<PairInt> expectedBorderPoints,
        Set<PairInt> borderPoints) {
        
        int nFound = 0;
        int nOutside = 0;
        int nFoundBorder = 0;
        for (PairInt p : points) {
            if (expectedPoints.contains(p)) {
                nFound++;
            } else {
                nOutside++;
            }
            if (expectedBorderPoints.contains(p)) {
                nFoundBorder++;
            }
        }

        SetComparisonResults r = new SetComparisonResults(expectedPoints.size(),
            nOutside, nFound, expectedBorderPoints.size(), nFoundBorder);
        
        return r;
        
    }

}
