package algorithms.bipartite;

import algorithms.CountingSort;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;

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
                
        TIntSet values = new TIntHashSet();
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        
        TObjectIntIterator<PairInt> iter = 
            g.getEdgeWeights().iterator();
        for (int i = g.getEdgeWeights().size(); i-- > 0;) {
            iter.advance();
            PairInt p = iter.key();
            int idx1 = p.getX();
            int idx2 = p.getY();

            int w = iter.value();
            values.add(w);
            int wInt = w;
            if (wInt < min) {
                min = wInt;
            }
            if (wInt > max) {
                max = wInt;
            }
        }
        
        int[] sortedValues = new int[values.size()];
        int count = 0;
        TIntIterator iter2 = values.iterator();
        while (iter2.hasNext()) {
            int v = iter2.next();
            sortedValues[count] = v;
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
        TIntIntMap map0To1 = new TIntIntHashMap();
        for (int v : sortedValues) {
            map0To1.put(v, count);
            count++;
        }

        TObjectIntIterator<PairInt> iter3 = 
            g.getEdgeWeights().iterator();
        for (int i = g.getEdgeWeights().size(); i-- > 0;) {
            iter3.advance();
            PairInt p = iter3.key();
            int w = iter3.value();
            int w2 = map0To1.get(w);
            iter3.setValue(w2);
        }
    }
}
