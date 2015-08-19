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
    public int createBitstring(Set<PairInt> subset) {
        
        int result = 0;
        
        for (PairInt p : subset) {
            Integer v = pointIndexLookup.get(p);
            if (v == null) {
                throw new IllegalArgumentException(
                    "a point in subset was not present in constructor's points");
            }
            setBit(result, v.intValue());
        }
        
        return result;
    }
    
    protected void setBit(int bitstring, int nthBit) {
        bitstring |= (1 << nthBit);
    }
    
    protected void clearBit(int bitstring, int nthBit) {
        bitstring &= ~(1 << nthBit);
    }
    
    protected void toggleBit(int bitstring, int nthBit) {
        bitstring ^= (1 << nthBit);
    }
    
    protected boolean isSet(int bitstring, int nthBit) {
        return ((bitstring & (1 << nthBit)) != 0);
    }
    
    protected boolean isNotSet(int bitstring, int nthBit) {
        return ((bitstring & (1 << nthBit)) == 0);
    }
  
}
