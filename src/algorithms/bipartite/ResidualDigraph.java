package algorithms.bipartite;

import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

/**
 *
 * @author nichole
 */
public class ResidualDigraph {
    
    /**
     * left (==X) vertices in graph
     */
    private final int nLeft;

    /**
     * right (==Y) vertices in graph
     */
    private final int nRight;

    /**
     * links X to Y (that is, left to right). 
     * These are "idle" arcs, f=0, in the
     * flow network N_G. They correspond to "unmarried" in the matched M graph.
     */
    private TIntObjectMap<TIntSet> forwardLinksRM
        = new TIntObjectHashMap<TIntSet>();

    /**
     * links Y to X (that is, right to left). These are "saturated" arcs, f=1,
     * in the flow network N_G. They correspond to "married" in the matched M
     * graph.
     */
    private TIntIntMap backwardLinksRM = new TIntIntHashMap();

    public ResidualDigraph(Graph g, TIntIntMap m) {                
        
        this.nLeft = g.getNLeft();
        this.nRight = g.getNRight();
        
        TObjectIntIterator<PairInt> iter = 
            g.getEdgeWeights().iterator();
        for (int i = g.getEdgeWeights().size(); i-- > 0;) {

            iter.advance();

            PairInt p = iter.key();
            int x = p.getX();
            int y = p.getY();
            
            TIntSet ys = forwardLinksRM.get(x);
            if (ys == null) {
                ys = new TIntHashSet();
                forwardLinksRM.put(x, ys);
            }
            ys.add(y);
        }
        
        TIntIntIterator iter2 = m.iterator();
        for (int i = m.size(); i-- > 0;) {
            iter2.advance();
            int x = iter2.key();
            int y = iter2.value();
                        
            backwardLinksRM.put(y, x);
            
            forwardLinksRM.get(x).remove(y);
        }
    }
            
    public ResidualDigraph(GraphWithoutWeights g) {                
    
        this.nLeft = g.getNLeft();
        this.nRight = g.getNRight();
        
        TIntObjectIterator<TIntSet> iter 
            = g.getAdjacencyMap().iterator();
        
        for (int i = g.getAdjacencyMap().size(); i-- > 0;) {

            iter.advance();
            
            int uIndex = iter.key();
            TIntSet vIndexes = iter.value();
            TIntIterator iter2 = vIndexes.iterator();
            while (iter2.hasNext()) {
                int v = iter2.next();
                TIntSet ys = forwardLinksRM.get(uIndex);
                if (ys == null) {
                    ys = new TIntHashSet();
                    forwardLinksRM.put(uIndex, ys);
                }
                ys.add(v);
            }
        }
    }
 
    public int countOfForwardBipartiteLinks() {
        int n = 0;
        
        TIntObjectIterator<TIntSet> iter = forwardLinksRM.iterator();
        for (int i = forwardLinksRM.size(); i-- > 0;) {
            iter.advance();            
            n += iter.value().size();
        }
        return n;
    }
    
    /**
     * @return nLeft number of left nodes
     */
    public int getNLeft() {
        return nLeft;
    }

   /**
     * @return number of right nodes
     */
    public int getNRight() {
        return nRight;
    }
    
    /**
     * @return the forwardLinksRM
     */
    public TIntObjectMap<TIntSet> getForwardLinksRM() {
        return forwardLinksRM;
    }

    /**
     * @return the backwardLinksRM
     */
    public TIntIntMap getBackwardLinksRM() {
        return backwardLinksRM;
    }

    public TIntIntMap extractMatchings() {

        TIntIntMap m = new TIntIntHashMap();
        
        TIntIntIterator iter = backwardLinksRM.iterator();
        for (int i = backwardLinksRM.size(); i-- > 0;) {
            iter.advance();            
            int rightIdx = iter.key();
            int leftIdx = iter.value();
            m.put(leftIdx, rightIdx);
        }
        
        return m;
    }

}
