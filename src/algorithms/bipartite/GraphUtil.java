package algorithms.bipartite;

import algorithms.CountingSort;
import algorithms.util.PairInt;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class GraphUtil {

    /**
     * condense edge weights to unique sequential values,
     * that is remove integer gaps.  any edge weights that had the
     * same value as one another before, still have the same
     * value as one another afterwards.
     * @param g 
     */
    public static void condenseEdgeWeights(Graph g) {
                
        Set<Integer> values = new HashSet<Integer>();
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        
        for (Entry<PairInt, Integer> entry : g.getEdgeWeights().entrySet()) {
            Integer w = entry.getValue();
            values.add(w);
            int wInt = w.intValue();
            if (wInt < min) {
                min = wInt;
            }
            if (wInt > max) {
                max = wInt;
            }
        }
        
        int[] sortedValues = new int[values.size()];
        int count = 0;
        for (Integer v : values) {
            sortedValues[count] = v.intValue();
            count++;
        }
        double nlg2n = sortedValues.length * 
            Math.log(sortedValues.length)/Math.log(2);
        
        // TODO: consider abandoning the "condense" if cannot
        //    use CountingSort because of the N*log_2(N) sort.
        if (nlg2n < max) {
            Arrays.sort(sortedValues);
        } else {
            sortedValues = CountingSort.sort(sortedValues, max);
        }
        
        count = 1;
        Map<Integer, Integer> map0To1 = new HashMap<Integer, Integer>();
        for (int v : sortedValues) {
            map0To1.put(Integer.valueOf(v), Integer.valueOf(count));
            count++;
        }
        
        for (Entry<PairInt, Integer> entry : g.getEdgeWeights().entrySet()) {
            Integer w = entry.getValue();
            Integer w2 = map0To1.get(w);
            PairInt p = entry.getKey();
            entry.setValue(w2);
        }
    }
}
