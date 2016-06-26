package algorithms.search;

import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.HashSet;
import java.util.Set;
import thirdparty.ods.Integerizer;
import thirdparty.ods.XFastTrie;
import thirdparty.ods.XFastTrieNode;

/**
 * a k-nearest neighbors using XFastTrie
 * for predecessor and successor queries
 * on spatially indexed numbers.
 * 
 * The algorithm performs better on dense data
 * (that is because the base of the prefix tree
 * is filled, leaving smaller number of nodes to
 * create in linear time).
 * The queries depend upon the maximum of x and
 * maximum of y to be entered.
 * 
 * A worst case query would be when column 0 is
 * filled with points and no others filled elsewhere, 
 * and the last point in the last row and last 
 * column is the query point.
 * In this worst case, the query time would scale
 * roughly as maxY * O(log_2(log_2(w))
 
  The algorithm starts with a predecessor and successor 
  call on the query point.  The minimum distance among
  those 2 becomes the goal to search to for completeness
  as rows above and below the query in the same column.
  The search continues to higher rows making predecessor and
  successor calls until the goal is reached.  The next higher
  row is one less than the predecessor result.
  The goal to the complete search is reduced by smaller distance
  answers.  The same is repeated for lower rows.
  <pre>
  an example would be:
  
     0   1   2   3   4

     5   6   7   8   9

     10  11 *12  13  14     q='18'.  pred='15', succ='null' --> goals(3, 23)
                                     13.pred='12'
    *15  16  17 *18  19              goals change to (13,23)
                                     13.succ='15', not closer than 12.
     20  21  22  23  24              23.pred and 23.succ not closer than 12
                             ans='12'.  queries: 3 pred, 3 succ queries.
                                     at O(log_2 log_2(w)) + O(w-l) each
                                     complexity was 
                                           6 * ( O(loglogw) + O(w-l) )
                                     for max index = 24, have w = 6 
 </ore>
 * @author nichole
 */
public class NearestNeighbor2D {
    
    private final XFastTrie<XFastTrieNode<Integer>, Integer> xbt;
            
    private final int maxX;
    
    private final int maxY;
    
    private final int maxIndex;
    
    /**
     * 
     * @param points
     * @param maxX maximum x value of any data point including
     *    those to be queries
     * @param maxY maximum y value of any data point including
     *    those to be queries
     */
    public NearestNeighbor2D(Set<PairInt> points,
        int maxX, int maxY) {
        
        this.maxX = maxX;
        this.maxY = maxY;
                
        maxIndex = maxX * maxY;
        
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
        
        int maxW = 1 + (int)Math.ceil(Math.log(maxIndex)/Math.log(2));
        
        xbt = new XFastTrie<XFastTrieNode<Integer>, Integer>(
            new XFastTrieNode<Integer>(), it, maxW);
        
        for (PairInt p : points) {
            int x = p.getX();
            int index = getInternalIndex(p.getX(), p.getY());
            xbt.add(Integer.valueOf(index));
        }
    }
    
    public int getInternalIndex(int col, int row) {
        return (row * maxX) + col;
    }
    
    protected int getRow(int internalIndex) {
        int row = internalIndex/maxX;        
        return row;
    }
    
    protected int getCol(int internalIndex) {
        int row = internalIndex/maxX;
        int col = internalIndex - (row * maxX);
        return col;
    }
    
    /**
    <pre>
      runtime complexity is
      
     
     </ore>
    
     * @param x
     * @param y
     */
    public Set<PairInt> findClosest(final int x, final int y) {
        
        if (x > maxX) {
            throw new IllegalArgumentException("x cannot be larger than "
                + " maxX given in constructor, " + maxX);
        }
        
        if (y > maxY) {
            throw new IllegalArgumentException("y cannot be larger than "
                + " maxY given in constructor, " + maxY);
        }
        
        int index = getInternalIndex(x, y);
        
        Integer q = xbt.find(Integer.valueOf(index));
        if (q != null) {
            Set<PairInt> results = new HashSet<PairInt>();
            results.add(new PairInt(x, y));
            return results;
        }
                
        TIntSet closestIndexes = new TIntHashSet();
        
        double closestDistance = Double.MAX_VALUE;
        
        Integer predecessor = xbt.predecessor(index);
        Integer successor = xbt.successor(index);
        
        double dp2 = dist(x, y, predecessor);
        double ds2 = dist(x, y, successor);
        if (dp2 <= ds2 && (dp2 != Double.MAX_VALUE)) {
            closestDistance = dp2;
            closestIndexes.add(predecessor.intValue());
            if (dp2 == ds2) {
                closestIndexes.add(successor.intValue());
            }
        } else if (ds2 < dp2) {
            closestDistance = ds2;
            closestIndexes.add(successor.intValue());
        }
        
        int goal = (closestDistance != Double.MAX_VALUE) ?
            (int)Math.ceil(closestDistance) : 0;
        
        int yLow = estimateLowBound(y, goal);
       
        int yCurrent = (predecessor == null) ? Integer.MIN_VALUE :
            (getRow(predecessor.intValue()) - 1);
        
        // predecessor searches until reach yLow, adjusting goal by
        //   min distances
        Integer p2 = null; 
        Integer s2 = null;
        while (yCurrent >= yLow) {
            int cIndex = getInternalIndex(x, yCurrent);
            q = xbt.find(Integer.valueOf(cIndex));
            if (q != null) {
                p2 = q;
                dp2 = dist(x, y, p2);
                ds2 = Double.MAX_VALUE;
            } else {
                p2 = xbt.predecessor(cIndex);
                s2 = xbt.successor(cIndex);
                dp2 = dist(x, y, p2);
                ds2 = dist(x, y, s2);
            }
            if ((dp2 < ds2) && (dp2 < closestDistance)) {
                closestIndexes.clear();
                closestDistance = dp2;
                closestIndexes.add(p2);
                goal = (int)Math.ceil(closestDistance);
                yLow = estimateLowBound(y, goal);                
            } else if ((ds2 < dp2) && (ds2 < closestDistance)) {
                closestIndexes.clear();
                closestDistance = ds2;
                closestIndexes.add(s2);
                goal = (int)Math.ceil(closestDistance);
                yLow = estimateLowBound(y, goal);
            } else if (dp2 == closestDistance && (dp2 != Double.MAX_VALUE)) {
                closestIndexes.add(p2.intValue());
                if (dp2 == ds2) {
                    closestIndexes.add(s2.intValue());
                }
            } else if (ds2 == closestDistance && (ds2 != Double.MAX_VALUE)) {
                closestIndexes.add(s2.intValue());
            }
            
            if (p2 != null) {
                int expectedNext = getInternalIndex(x, yCurrent - 1);
                if (p2.intValue() > expectedNext) {
                    yCurrent -= 1;
                } else {
                    yCurrent = getRow(p2.intValue()) - 1;
                }
            } else {
                yCurrent = Integer.MIN_VALUE;
            }
        }
        
        // successor searches to higher bounds
        yCurrent = (successor == null) ? Integer.MAX_VALUE :
            (getRow(successor.intValue()) + 1);
        int yHigh = estimateHighBound(y, goal);
        
        while (yCurrent <= yHigh) {
            int cIndex = getInternalIndex(x, yCurrent);
            q = xbt.find(Integer.valueOf(cIndex));
            if (q != null) {
                p2 = q;
                dp2 = dist(x, y, p2);
                ds2 = Double.MAX_VALUE;
            } else {
                p2 = xbt.predecessor(cIndex);
                s2 = xbt.successor(cIndex);
                dp2 = dist(x, y, p2);
                ds2 = dist(x, y, s2);
            }
            if ((dp2 < ds2) && (dp2 < closestDistance)) {
                closestIndexes.clear();
                closestDistance = dp2;
                closestIndexes.add(p2);
                goal = (int)Math.ceil(closestDistance);
                yHigh = estimateHighBound(y, goal);
            } else if ((ds2 < dp2) && (ds2 < closestDistance)) {
                closestIndexes.clear();
                closestDistance = ds2;
                closestIndexes.add(s2);
                goal = (int)Math.ceil(closestDistance);
                yHigh = estimateHighBound(y, goal);
            } else if (dp2 == closestDistance && (dp2 != Double.MAX_VALUE)) {
                closestIndexes.add(p2.intValue());
                if (dp2 == ds2) {
                    closestIndexes.add(s2.intValue());
                }
            } else if (ds2 == closestDistance && (ds2 != Double.MAX_VALUE)) {
                closestIndexes.add(s2.intValue());
            } 
            
            if (s2 != null) {
                int expectedNext = getInternalIndex(x, yCurrent + 1);
                if (s2.intValue() < expectedNext) {
                    yCurrent += 1;
                } else {
                    yCurrent = getRow(s2.intValue()) + 1;
                }
            } else {
                yCurrent = Integer.MAX_VALUE;
            }
        }
        
        Set<PairInt> results = new HashSet<PairInt>();
        TIntIterator iter = closestIndexes.iterator();
        while (iter.hasNext()) {
            int index2 = iter.next();
            int x2 = getCol(index2);
            int y2 = getRow(index2);
            results.add(new PairInt(x2, y2));
        }
        
        return results;
    }
    
    /**
    <pre>
      runtime complexity is
      
     
     </ore>
    
     * @param x
     * @param y
     * @param dMax
     */
    public Set<PairInt> findClosest(int x, int y, int dMax) {
        
        if (x > maxX) {
            throw new IllegalArgumentException("x cannot be larger than "
                + " maxX given in constructor, " + maxX);
        }
        
        if (y > maxY) {
            throw new IllegalArgumentException("y cannot be larger than "
                + " maxY given in constructor, " + maxY);
        }
        
        int index = getInternalIndex(x, y);
        
        Integer q = xbt.find(Integer.valueOf(index));
        if (q != null) {
            Set<PairInt> results = new HashSet<PairInt>();
            results.add(new PairInt(x, y));
            return results;
        }
                
        TIntSet closestIndexes = new TIntHashSet();
        
        double closestDistance = Double.MAX_VALUE;
        
        Integer predecessor = xbt.predecessor(index);
        Integer successor = xbt.successor(index);
        
        double dp2 = dist(x, y, predecessor);
        double ds2 = dist(x, y, successor);
        if (dp2 <= ds2 && (dp2 <= dMax)) {
            closestDistance = dp2;
            closestIndexes.add(predecessor.intValue());
            if (dp2 == ds2) {
                closestIndexes.add(successor.intValue());
            }
        } else if (ds2 < dp2 && (ds2 <= dMax)) {
            closestDistance = ds2;
            closestIndexes.add(successor.intValue());
        }
        
        int goal = (closestDistance != Double.MAX_VALUE) ?
            (int)Math.ceil(closestDistance) : 0;
        
        if (goal > dMax) {
            goal = dMax;
        }
        
        int yLow = estimateLowBound(y, goal);
       
        int yCurrent = (predecessor == null) ? Integer.MIN_VALUE :
            (getRow(predecessor.intValue()) - 1);
        
        // predecessor searches until reach yLow, adjusting goal by
        //   min distances
        Integer p2 = null; 
        Integer s2 = null;
        while (yCurrent >= yLow) {
            int cIndex = getInternalIndex(x, yCurrent);
            q = xbt.find(Integer.valueOf(cIndex));
            if (q != null) {
                p2 = q;
                dp2 = dist(x, y, p2);
                ds2 = Double.MAX_VALUE;
            } else {
                p2 = xbt.predecessor(cIndex);
                s2 = xbt.successor(cIndex);
                dp2 = dist(x, y, p2);
                ds2 = dist(x, y, s2);
            }
            if ((dp2 < ds2) && (dp2 < closestDistance) && (dp2 <= dMax)) {
                closestIndexes.clear();
                closestDistance = dp2;
                closestIndexes.add(p2);
                goal = (int)Math.ceil(closestDistance);
                if (goal > dMax) {
                    goal = dMax;
                }
                yLow = estimateLowBound(y, goal);                
            } else if ((ds2 < dp2) && (ds2 < closestDistance) && (ds2 <= dMax)) {
                closestIndexes.clear();
                closestDistance = ds2;
                closestIndexes.add(s2);
                goal = (int)Math.ceil(closestDistance);
                if (goal > dMax) {
                    goal = dMax;
                }
                yLow = estimateLowBound(y, goal);
            } else if (dp2 == closestDistance && (dp2 != Double.MAX_VALUE)) {
                closestIndexes.add(p2.intValue());
                if (dp2 == ds2) {
                    closestIndexes.add(s2.intValue());
                }
            } else if (ds2 == closestDistance && (ds2 != Double.MAX_VALUE)) {
                closestIndexes.add(s2.intValue());
            }
            
            if (p2 != null) {
                int expectedNext = getInternalIndex(x, yCurrent - 1);
                if (p2.intValue() > expectedNext) {
                    yCurrent -= 1;
                } else {
                    yCurrent = getRow(p2.intValue()) - 1;
                }
            } else {
                yCurrent = Integer.MIN_VALUE;
            }
        }
        
        // successor searches to higher bounds
        yCurrent = (successor == null) ? Integer.MAX_VALUE :
            (getRow(successor.intValue()) + 1);
        int yHigh = estimateHighBound(y, goal);
        
        while (yCurrent <= yHigh) {
            int cIndex = getInternalIndex(x, yCurrent);
            q = xbt.find(Integer.valueOf(cIndex));
            if (q != null) {
                p2 = q;
                dp2 = dist(x, y, p2);
                ds2 = Double.MAX_VALUE;
            } else {
                p2 = xbt.predecessor(cIndex);
                s2 = xbt.successor(cIndex);
                dp2 = dist(x, y, p2);
                ds2 = dist(x, y, s2);
            }
            if ((dp2 < ds2) && (dp2 < closestDistance) && (dp2 <= dMax)) {
                closestIndexes.clear();
                closestDistance = dp2;
                closestIndexes.add(p2);
                goal = (int)Math.ceil(closestDistance);
                if (goal > dMax) {
                    goal = dMax;
                }
                yHigh = estimateHighBound(y, goal);
            } else if ((ds2 < dp2) && (ds2 < closestDistance) && (ds2 <= dMax)) {
                closestIndexes.clear();
                closestDistance = ds2;
                closestIndexes.add(s2);
                goal = (int)Math.ceil(closestDistance);
                if (goal > dMax) {
                    goal = dMax;
                }
                yHigh = estimateHighBound(y, goal);
            } else if (dp2 == closestDistance && (dp2 != Double.MAX_VALUE)) {
                closestIndexes.add(p2.intValue());
                if (dp2 == ds2) {
                    closestIndexes.add(s2.intValue());
                }
            } else if (ds2 == closestDistance && (ds2 != Double.MAX_VALUE)) {
                closestIndexes.add(s2.intValue());
            } 
            
            if (s2 != null) {
                int expectedNext = getInternalIndex(x, yCurrent + 1);
                if (s2.intValue() < expectedNext) {
                    yCurrent += 1;
                } else {
                    yCurrent = getRow(s2.intValue()) + 1;
                }
            } else {
                yCurrent = Integer.MAX_VALUE;
            }
        }
        
        Set<PairInt> results = new HashSet<PairInt>();
        TIntIterator iter = closestIndexes.iterator();
        while (iter.hasNext()) {
            int index2 = iter.next();
            int x2 = getCol(index2);
            int y2 = getRow(index2);
            results.add(new PairInt(x2, y2));
        }
        
        return results;
    }

    private double dist(int x, int y, Integer p2) {
        
        if (p2 == null) {
            return Double.MAX_VALUE;
        }
        
        int x2 = getCol(p2.intValue());
        int y2 = getRow(p2.intValue());
        
        int diffX = x2 - x;
        int diffY = y2 - y;
        double dist = Math.sqrt(diffX * diffX + diffY * diffY);
        
        return dist;
    }

    private int estimateLowBound(int y, int goal) {

        int low = y - goal;
        
        if (low < 0) {
            low = 0;
        }
        
        return low;
    }
    
    private int estimateHighBound(int y, int goal) {

        int high = y + goal;
        
        if (high > maxY) {
            high = maxY;
        }
        
        return high;
    }
}
