package algorithms.util;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class BitVectorRepresentation {
    
    private final Map<PairInt, Integer> pointIndexLookup = new HashMap<PairInt, Integer>();
    
    public BitVectorRepresentation(Set<PairInt> points) {
        
        // enumerate the points
        for (PairInt p : points) {
            pointIndexLookup.put(p, Integer.valueOf(pointIndexLookup.size()));
        }
    }
    
    /**
     * given a subset of the points this instance was constructed with, 
     * create a bitstring representing the subset membership.
     * @param subset
     * @return 
     */
    public VeryLongBitString createBitstring(Set<PairInt> subset) {
                
        VeryLongBitString result = new VeryLongBitString(pointIndexLookup.size());
        
        for (PairInt p : subset) {
            Integer v = pointIndexLookup.get(p);
            if (v == null) {
                throw new IllegalArgumentException(
                    "a point in subset was not present in constructor's points");
            }
            result.setBit(v.intValue());
        }
        
        return result;
    }
    
}
