package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class Splice {
    
    private final PairIntArray edge;
        
    public Splice(PairIntArray theEdge) {
        edge = theEdge;
    }

    public int getLengthOfLongestSide(final int spliceIndex) {
        // the count always includes the splice
        int n0 = spliceIndex + 1;
        int n1 = edge.getN() - spliceIndex;
        if (n0 > n1) {
            return n0;
        }
        return n1;
    }

    /**
     * splice the edge at spliceIndex and return the 2 pieces as the
     * longest edge then the other. Note that the splice index is in the longer
     * edge, so spliceIndexInOut gets populated with the value w.r.t. the
     * first item returned.
     * @return 
     */
    public PairIntArray[] splice(int[] spliceIndexInOut) {
        
        return splice(spliceIndexInOut, null, null, null, null);
    }
 
    /**
     * splice the edge at spliceIndexInOut and return the splices ordered
     * by the largest splice as item 0.  spliceIndexInOut is updated 
     * with respect to item 0 and the editJunctionMap is edited for any
     * location value affected by the splice. 
     * @param spliceIndexInOut variable to give splice index and return
     * the updated value with respect to the first item in the spliced edge
     * @param editJunctionLocationMap map w/ key = pixel index and
     * value = edge location pair (edges index, index within edge).
     * (can be null).
     * @param edgePixelIndexes pixel indexes of junctions within this edge.
     * Note, it's up to the invoker to validate that the input is for this same
     * edge. (can be null).
     * @param outputLargerSpliceJunctionLocationMap
     * (can be null if edgePixelIndexes is null)
     * @param outputSmallerSpliceJunctionLocationMap
     * (can be null if edgePixelIndexes is null)
     * @return PairIntArray[]{largerSplice, smallerSplice} along with side
     * effect of spliceIndexInOut[0] updated w.r.t. largerSplice and
     * outputLargerSpliceJunctionLocationMap and 
     * outputSmallerSpliceJunctionLocationMap containing their entries
     * from editJunctionLocationMap updated for the splice.
     * (editJunctionLocationMap is not updated directly, just in case the
     * splice isn't used by the invoker).
     */
    public PairIntArray[] splice(int[] spliceIndexInOut, 
        final Set<Integer> edgePixelIndexes, 
        final Map<Integer, PairInt> editJunctionLocationMap,
        Map<Integer, PairInt> outputSmallerSpliceJunctionLocationMap,
        Map<Integer, PairInt> outputLargerSpliceJunctionLocationMap) {
        
        if (spliceIndexInOut == null || spliceIndexInOut.length != 1) {
            throw new IllegalArgumentException(
            "spliceIndexInOut must have length of 1");
        }
        
        if (edgePixelIndexes != null) {
            if (outputLargerSpliceJunctionLocationMap == null) {
                throw new IllegalArgumentException(
                "outputLargerSpliceJunctionLocationMap cannot be null if edgePixelIndexes is not null");
            }
            if (outputSmallerSpliceJunctionLocationMap == null) {
                throw new IllegalArgumentException(
                "outputSmallerSpliceJunctionLocationMap cannot be null if edgePixelIndexes is not null");
            }
        }
        
        if ((spliceIndexInOut[0] == 0) || (spliceIndexInOut[0] == (edge.getN() - 1))) {
            // --- add entries to output junction maps
            if (edgePixelIndexes != null) {
                for (Integer pixelIndex : edgePixelIndexes) {
                    PairInt loc = editJunctionLocationMap.get(pixelIndex);
                    assert(loc != null);
                    outputLargerSpliceJunctionLocationMap.put(
                        pixelIndex, new PairInt(loc.getX(), loc.getY()));
                }
            }
            return new PairIntArray[]{edge.copy(), new PairIntArray()};
        }

        int n0 = spliceIndexInOut[0] + 1;
        int n1 = edge.getN() - spliceIndexInOut[0];
        
        boolean edge0IsLonger = true;
        if (n0 > n1) {
            --n1;
            if (n1 < 0) {
                n1 = 0;
            }
        } else {
            --n0;
            if (n0 < 0) {
                n0 = 0;
            }
            edge0IsLonger = false;
        }
        
        /*
              J
        0  1  2  3  4  5  n=6
              sp
        n0=3   n1=4
        n0=2
        */

        PairIntArray edge0 = new PairIntArray(n0);
        PairIntArray edge1 = new PairIntArray(n1);

        for (int i = 0; i < n0; ++i) {
            edge0.add(edge.getX(i), edge.getY(i));    
        }
        for (int i = 0; i < n1; ++i) {
            edge1.add(edge.getX(n0 + i), edge.getY(n0 + i));
        }
        
        if (edgePixelIndexes != null) {
            // ----- update output junction maps for splices ------
            for (Integer pixelIndex : edgePixelIndexes) {
                PairInt loc = editJunctionLocationMap.get(pixelIndex);
                assert(loc != null);
                int idxWithinEdge = loc.getY();

                /*
                for edge0 being longer:
                              J
                     0  1  2  3  4  5  n=6
                              sp
                     n0=4        n1=2

                for edge1 being longer:
                           J
                     0  1  2  3  4  5  n=6
                           sp
                     n0=2        n1=4
                */

                if (edge0IsLonger) {
                    if (idxWithinEdge < n0) {
                        outputLargerSpliceJunctionLocationMap.put(
                            pixelIndex, new PairInt(loc.getX(), loc.getY()));
                    } else {
                        idxWithinEdge -= n0;
                        outputSmallerSpliceJunctionLocationMap.put(
                            pixelIndex, new PairInt(loc.getX(), idxWithinEdge));
                    }
                } else {
                    if (idxWithinEdge < n0) {
                        outputSmallerSpliceJunctionLocationMap.put(
                            pixelIndex, new PairInt(loc.getX(), loc.getY()));
                    } else {
                        idxWithinEdge -= n0;
                        outputLargerSpliceJunctionLocationMap.put(
                            pixelIndex, new PairInt(loc.getX(), idxWithinEdge));
                    }
                }
            }
        }

        if (edge0IsLonger) {
            /*
                      J
             0  1  2  3  4  5  n=6
                      sp
             n0=4     n1=3
                      n1=2
            */
            
            //spliceIndexInOut remains the same
            
            return new PairIntArray[]{edge0, edge1};
        }
        
        spliceIndexInOut[0] = 0;
        
        return new PairIntArray[]{edge1, edge0};
    }
    
}
