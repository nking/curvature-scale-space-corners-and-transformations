package algorithms.misc;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.HashSet;
import java.util.Set;

/**
 * miscellaneous boiler plate code
 * 
 * @author nichole
 */
public class Misc {
    
    public static Set<PairInt> convert(PairIntArray points) {
        
        if (points == null) {
            throw new IllegalArgumentException("points cannot be null");
        }
        
        Set<PairInt> out = new HashSet<PairInt>();
        
        for (int i = 0; i < points.getN(); ++i) {
            out.add(new PairInt(points.getX(i), points.getY(i)));
        }
        
        return out;
    }
}
