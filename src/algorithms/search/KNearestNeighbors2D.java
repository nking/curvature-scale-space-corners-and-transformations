package algorithms.search;

import algorithms.imageProcessing.FixedSizeSortedVector;
import algorithms.util.PairInt;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.util.Set;
import thirdparty.ods.Integerizer;
import thirdparty.ods.XFastTrie;
import thirdparty.ods.XFastTrieNode;

/**
 * a k-nearest neighbors using XFastTries
 * for predecessor and successor queries
 * on x and y tries followed by a combinatorial
 * check for point exists.
 * 
 * Note that the Fractional Cascading Layered Trees
 * range queries should be faster than this.
 * 
 * with
 * @author nichole
 */
public class KNearestNeighbors2D {
    
    private final XFastTrie<XFastTrieNode<Integer>, Integer> xbt;
    
    private final XFastTrie<XFastTrieNode<Integer>, Integer> ybt;
    
    private final TIntList cachedListX = new TIntArrayList();
    
    private final TIntList cachedListY = new TIntArrayList();
    
    private final Set<PairInt> set;
    
    /**
     * 
     * @param points
     * @param maxX
     * @param maxY 
     */
    public KNearestNeighbors2D(Set<PairInt> points,
        int maxX, int maxY) {
        
        set = points;
        
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
        
        XFastTrieNode<Integer> node = new XFastTrieNode<Integer>();
        
        xbt = new XFastTrie<XFastTrieNode<Integer>, Integer>(
            node, it, maxX);
        
        ybt = new XFastTrie<XFastTrieNode<Integer>, Integer>(
            node, it, maxY);
        
        for (PairInt p : points) {
            int x = p.getX();
            xbt.add(Integer.valueOf(x));
            int y = p.getY();
            ybt.add(Integer.valueOf(y));
        }
    }
    
    /**
     * runtime complexity is
     * O(log2(xW) * number of unique x points within +- dMax
     * +
     * O(log2(yW) * number of unique y points within +- dMax
     * +
     * less than nX * nY * log_2(k)
     * 
     * where xW and yW are the bit lengths of maxX and maxY,
     * respectively.
     * 
     * NOTE that a greedy faster version of this method
     * can be made that doesn't use the log_2(k) term with
     * the understanding that the user sets dMax and k
     * such that it isn't needed.
     * 
     * @param x
     * @param y
     * @param dMax
     * @param vec 
     */
    public void findClosest(int x, int y, int dMax, 
        FixedSizeSortedVector<PairDistance> vec) {
        
        // O(log2(xW) * number of unique x points within +- dMax
        int nX = findClosest(x, xbt, dMax, cachedListX);
        
        // O(log2(yW) * number of unique y points within +- dMax
        int nY = findClosest(y, ybt, dMax, cachedListY);
    
        int dMaxSq = dMax * dMax;
            
        // nX * nY   (TODO: this could be improved)
        for (int i = 0; i < nX; ++i) {
            int x2 = cachedListX.get(i);
            for (int j = 0; j < nY; ++j) {
                int y2 = cachedListY.get(j);
                
                PairInt p2 = new PairInt(x2, y2);
                
                if (set.contains(p2)) {

                    int diffX = x2 - x;
                    int diffY = y2 - y;
                    int distSq = (diffX * diffX) + (diffY * diffY);
                    
                    //TODO: here, could allow an argument for
                    // distance type to allow diagonals to be
                    // included
                    if (distSq > dMaxSq) {
                        continue;
                    }

                    // O(log_2(k))
                    vec.add(new PairDistance(p2, distSq));
                }
            }
        }
    }
    
    /**
     * runtime complexity is
     * O(log2(xW) * number of unique x points within +- dMax
     * +
     * O(log2(yW) * number of unique y points within +- dMax
     * +
     * nX * nY
     * 
     * where xW and yW are the bit lengths of maxX and maxY,
     * respectively.
     * 
     * NOTE that this version is greedy, that is, the first
     * k points within dMax of (x, y) are returned
     * and they might not be the closest, so the user
     * has to carefully keep dMax small and k large enough.
     * 
     * @param x
     * @param y
     * @param dMax
     * @param k
     * @param output
     * @return number of points filled in output.
     */
    public int findClosestUsingGreedy(int x, int y, int dMax, 
        int k, PairDistance[] output) {
        
        // O(log2(xW) * number of unique x points within +- dMax
        int nX = findClosest(x, xbt, dMax, cachedListX);
        
        // O(log2(yW) * number of unique y points within +- dMax
        int nY = findClosest(y, ybt, dMax, cachedListY);
    
        int dMaxSq = dMax * dMax;
        
        int n = 0;
            
        // nX * nY   (TODO: this could be improved)
        for (int i = 0; i < nX; ++i) {
            int x2 = cachedListX.get(i);
            for (int j = 0; j < nY; ++j) {
                int y2 = cachedListY.get(j);
                
                PairInt p2 = new PairInt(x2, y2);
                
                if (set.contains(p2)) {

                    int diffX = x2 - x;
                    int diffY = y2 - y;
                    int distSq = (diffX * diffX) + (diffY * diffY);
                    if (distSq > dMaxSq) {
                        continue;
                    }
                    
                    int dist = (int)Math.round(Math.sqrt(distSq));

                    output[n] = new PairDistance(p2, dist);
                    
                    n++;
                }
            }
        }
        
        return n;
    }
    
    private int findClosest(int value, 
        XFastTrie<XFastTrieNode<Integer>, Integer> bt, 
        int dMax, TIntList cachedList) {
        
        Integer vIndex = Integer.valueOf(value);
               
        int n = 0;
        Integer v0 = bt.find(vIndex);
        if (v0 != null) {
            if (cachedList.size() < n) {
                cachedList.insert(n, v0.intValue());
            } else {
                cachedList.add(v0.intValue());
            }
            n++;
        }
        
        v0 = vIndex;
        while (true) {
            v0 = bt.predecessor(v0);
            if (v0 == null) {
                break;
            }
            int diff = value - v0.intValue();
            if (diff > dMax) {
                break;
            }
            if (cachedList.size() < n) {
                cachedList.insert(n, v0.intValue());
            } else {
                cachedList.add(v0.intValue());
            }
            n++;
        }
        v0 = vIndex;
        while (true) {
            v0 = bt.successor(v0);
            if (v0 == null) {
                break;
            }
            int diff = v0.intValue() - value;
            if (diff > dMax) {
                break;
            }
            if (cachedList.size() < n) {
                cachedList.insert(n, v0.intValue());
            } else {
                cachedList.add(v0.intValue());
            }
            n++;
        }
        
        return n;
    }
}
