package algorithms.compGeometry;

import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;

/**
 * class to hold some of the newer methods
 * w.r.t. boundary point extractions.
 * 
 * @author nichole
 */
public class PerimeterFinder2 {
   
    /**
     * finds the pixels with neighbors not in given point
     * set, contiguousPoints.
     * 
     * @param contiguousPoints
     * @return 
     */
    public Set<PairInt> extractBorder(Set<PairInt> contiguousPoints) {

        Set<PairInt> border = new HashSet<PairInt>();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        for (PairInt p : contiguousPoints) {
            int x = p.getX();
            int y = p.getY();
            for (int i = 0; i < dxs.length; ++i) {
                int x2 = x + dxs[i];
                int y2 = y + dys[i];
                PairInt p2 = new PairInt(x2, y2);
                if (!contiguousPoints.contains(p2)) {
                    border.add(p);
                    break;
                }
            }
        }
        
        return border;
    }

}
