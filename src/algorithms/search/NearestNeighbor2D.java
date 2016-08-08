package algorithms.search;

import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.HashSet;
import java.util.Set;
import thirdparty.ods.Integerizer;
import thirdparty.ods.XFastTrie;
import thirdparty.ods.XFastTrieNode;

/**
 * a nearest neighbor's algorithm using XFastTrie
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
 * roughly as maxY * O(log_2(w))
 
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
                                     at O(log_2(w)) each
                                     complexity was 
                                           6 * O(log_2(w))
                                     for max index = 24, have w = 6 
 </ore>
 * @author nichole
 */
public class NearestNeighbor2D {
    
    private final XFastTrie<XFastTrieNode<Integer>, Integer> xbt;
            
    private final int width;
    
    private final int height;
    
    private final int maxIndex;
    
    private boolean useCache = true;
    
    private TIntIntMap pCache = new TIntIntHashMap();
    private TIntIntMap sCache = new TIntIntHashMap();
    
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
        
        this.width = maxX + 1;
        this.height = maxY + 1;
                
        maxIndex = width * height;
        
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
    
    public void doNotUseCache() {
        useCache = false;
    }
    
    protected int getInternalIndex(int col, int row) {
        return (row * width) + col;
    }
    
    protected int getRow(int internalIndex) {
        int row = internalIndex/width;        
        return row;
    }
    
    protected int getCol(int internalIndex) {
        int row = internalIndex/width;
        int col = internalIndex - (row * width);
        return col;
    }
    
    /**
    <pre>
      runtime complexity is
         best case: 2 * O(log_2(maxW)).
            Note that caching leads to an O(1) term
            over time instead of the logarithmic term.
            
         worst case: nRows * 2 * O(log_2(maxW))
         
         Note, worst case is: first column
         filled with points and all else is empty and
         the number of rows is same or larger than 
         number of columns and the
         query is for the point in the last column and
         last row... a predecessor call is necessary for
         each row in the worst case.
          
     Note: maxW = 1 + Math.ceil(Math.log(maxX * maxY)/Math.log(2));            
     </ore>
    
     * @param x
     * @param y
     */
    public Set<PairInt> findClosest(final int x, final int y) {
        
        return findClosestWithinTolerance(x, y, 0);
    }
    
    /**
     * NOTE: NOT READY FOR USE
     * method to return only the nearest point and any
     * that are at the same distance within a tolerance.
     * This is meant to be a nearest neighbor method
     * with a tolerance.  
     * TODO: calculate the runtime complexity bounds....
     * @param x
     * @param y
     * @param tolerance
     * @return 
     */
    public Set<PairInt> findClosestWithinTolerance(int x, int y,
        double tolerance) {
        
        if (x > width) {
            throw new IllegalArgumentException("x cannot be larger than "
                + " maxX given in constructor, " + width);
        }
        
        if (y > height) {
            throw new IllegalArgumentException("y cannot be larger than "
                + " maxY given in constructor, " + height);
        }
        
        double closestDist = Double.MAX_VALUE;
        double closestDistPlusTol = Double.MAX_VALUE;
        
        TIntSet closestIndexes = new TIntHashSet();
        
        int idx = getInternalIndex(x, y);
        Integer index = Integer.valueOf(idx);
        
        //O(1)
        Integer q = xbt.find(index);
        if (q != null) {
            // have found nearest, but still need to search
            // within tolerance distance for others.
            closestDist = 0;
            closestDistPlusTol = tolerance;
            closestIndexes.add(index.intValue());
        }
                
        Integer predecessor = null;
        Integer successor = null;
        
        if (useCache && pCache.containsKey(idx)) {
            predecessor = Integer.valueOf(pCache.get(idx));
        } else {
            //O(log_2(maxW))
            predecessor = xbt.predecessor(index);
            if (useCache && predecessor != null) {
                pCache.put(idx, predecessor.intValue());
            }
        }
        if (useCache && sCache.containsKey(idx)) {
            successor = Integer.valueOf(sCache.get(idx));
        } else {
            //O(log_2(maxW))
            successor = xbt.successor(index);
            if (useCache && successor != null){
                sCache.put(idx, successor.intValue());
            }
        }
        
        double dp2 = dist(x, y, predecessor);
        double ds2 = dist(x, y, successor);
        double dMin = Math.min(dp2, ds2);
        
        /*
        if smallest is smaller than closest 
           if the new closest diff with current is greater 
               than tol, clear the indexes and reset closest 
               vars and add smallest to indexes
               also add the other if within tolerance
           else if closer is within tolerance,
              update closest vars and add whichever or both 
              s2 and p2 to indexes (delaying detailed checks 
              of indexes until end of method)
        else if smallest is <= closestPlusTol
            add s2 and/or p2 to indexes
        */
        if (dMin <= closestDist) {
            if (Math.abs(closestDist - dMin) > tolerance) {
                closestIndexes.clear();
                closestDist = dMin;
                closestDistPlusTol = closestDist + tolerance;
            }
            if (dp2 <= closestDistPlusTol) {
                closestIndexes.add(predecessor.intValue());
            }
            if (ds2 <= closestDistPlusTol) {
                closestIndexes.add(successor.intValue());
            }
        } else if (dMin <= closestDistPlusTol) {
            if (dp2 <= closestDistPlusTol) {
                closestIndexes.add(predecessor.intValue());
            }
            if (ds2 <= closestDistPlusTol) {
                closestIndexes.add(successor.intValue());
            }
        }
        
        //add tolerance to goal
        int goal = (closestDist != Double.MAX_VALUE) ?
            (int)Math.ceil(closestDistPlusTol) : 0;
        
        int yLow = estimateLowBound(y, goal);
       
        int yCurrent;
        if (predecessor == null) {
            yCurrent = Integer.MIN_VALUE;
        } else {
            int pRow = getRow(predecessor.intValue());
            if (pRow < y) {
                yCurrent = pRow;
            } else {
                yCurrent = pRow - 1;
            }
        }
        
        // predecessor searches until reach yLow, adjusting goal by
        //   min distances
        Integer p2 = null; 
        Integer s2 = null;
        while (yCurrent >= yLow) {
            int cIdx = getInternalIndex(x, yCurrent);
            Integer cIndex = Integer.valueOf(cIdx);

            //O(1)
            q = xbt.find(cIndex);
            if (q != null) {
                p2 = q;
                dp2 = dist(x, y, p2);
                ds2 = Double.MAX_VALUE;
            } else {
                if (useCache && pCache.containsKey(cIdx)) {
                    p2 = Integer.valueOf(pCache.get(cIdx));
                } else {
                    //O(log_2(maxW))
                    p2 = xbt.predecessor(cIndex);
                    if (useCache && p2 != null) {
                        pCache.put(cIdx, p2.intValue());
                    }
                }
                if (useCache && sCache.containsKey(cIdx)) {
                    s2 = Integer.valueOf(sCache.get(cIdx));
                } else {
                    //O(log_2(maxW))
                    s2 = xbt.successor(cIndex);
                    if (useCache && s2 != null) {
                        sCache.put(cIdx, s2.intValue());
                    }
                }
                
                dp2 = dist(x, y, p2);
                ds2 = dist(x, y, s2);
            }
        
            dMin = Math.min(dp2, ds2);
            if (dMin <= closestDist) {
                if (Math.abs(closestDist - dMin) > tolerance) {
                    closestIndexes.clear();
                    closestDist = dMin;
                    closestDistPlusTol = closestDist + tolerance;
                    goal = (int)Math.ceil(closestDistPlusTol);
                    yLow = estimateLowBound(y, goal); 
                }
                if (dp2 <= closestDistPlusTol) {
                    closestIndexes.add(p2.intValue());
                }
                if (ds2 <= closestDistPlusTol) {
                    closestIndexes.add(s2.intValue());
                }
            } else if (dMin <= closestDistPlusTol) {
                if (dp2 <= closestDistPlusTol) {
                    closestIndexes.add(p2.intValue());
                }
                if (ds2 <= closestDistPlusTol) {
                    closestIndexes.add(s2.intValue());
                }
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
        if (successor == null) {
            yCurrent = Integer.MAX_VALUE;
        } else {
            int sr = getRow(successor.intValue());
            if (sr > y) {
                yCurrent = sr;
            } else {
                yCurrent = sr + 1;
            }
        }
        int yHigh = estimateHighBound(y, goal);
        
        while (yCurrent <= yHigh) {
            int cIdx = getInternalIndex(x, yCurrent);
            Integer cIndex = Integer.valueOf(cIdx);
            
            //O(1)
            q = xbt.find(cIndex);
            if (q != null) {
                p2 = q;
                dp2 = dist(x, y, p2);
                ds2 = Double.MAX_VALUE;
            } else {
                if (useCache && pCache.containsKey(cIdx)) {
                    p2 = Integer.valueOf(pCache.get(cIdx));
                } else {
                    //O(log_2(maxW))
                    p2 = xbt.predecessor(cIndex);
                    if (useCache && p2 != null) {
                        pCache.put(cIdx, p2.intValue());
                    }
                }
                if (useCache && sCache.containsKey(cIdx)) {
                    s2 = Integer.valueOf(sCache.get(cIdx));
                } else {
                    //O(log_2(maxW))
                    s2 = xbt.successor(cIndex);
                    if (useCache && s2 != null) {
                        sCache.put(cIdx, s2.intValue());
                    }
                }
                dp2 = dist(x, y, p2);
                ds2 = dist(x, y, s2);
            }
            
            dMin = Math.min(dp2, ds2);
            if (dMin <= closestDist) {
                if (Math.abs(closestDist - dMin) > tolerance) {
                    closestIndexes.clear();
                    closestDist = dMin;
                    closestDistPlusTol = closestDist + tolerance;
                    goal = (int)Math.ceil(closestDistPlusTol);
                    yHigh = estimateHighBound(y, goal); 
                }
                if (dp2 <= closestDistPlusTol) {
                    closestIndexes.add(p2.intValue());
                }
                if (ds2 <= closestDistPlusTol) {
                    closestIndexes.add(s2.intValue());
                }
            } else if (dMin <= closestDistPlusTol) {
                if (dp2 <= closestDistPlusTol) {
                    closestIndexes.add(p2.intValue());
                }
                if (ds2 <= closestDistPlusTol) {
                    closestIndexes.add(s2.intValue());
                }
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
        
        //filter results for closest and tolerance
        Set<PairInt> results = new HashSet<PairInt>();
        TIntIterator iter = closestIndexes.iterator();
        while (iter.hasNext()) {
            int index2 = iter.next();
            if (dist(x, y, index2) <= closestDistPlusTol) {
                int x2 = getCol(index2);
                int y2 = getRow(index2);
                PairInt p3 = new PairInt(x2, y2);
                results.add(p3);
            }
        }
 
        return results;
    }
    
    /**
    <pre>
      runtime complexity is
         best case: 2 * O(log_2(maxW)).
            Note that caching leads to an O(1) term
            over time instead of the logarithmic term.
            
         worst case: dMax * 4 * O(log_2(maxW))
         
      Note: maxW = 1 + Math.ceil(Math.log(maxX * maxY)/Math.log(2));
     </ore>
    
     * @param x
     * @param y
     * @param dMax
     */
    public Set<PairInt> findClosest(int x, int y, int dMax) {
        
        if (x > width) {
            throw new IllegalArgumentException("x cannot be larger than "
                + " maxX given in constructor, " + width);
        }
        
        if (y > height) {
            throw new IllegalArgumentException("y cannot be larger than "
                + " maxY given in constructor, " + height);
        }
        
        int idx = getInternalIndex(x, y);
        Integer index = Integer.valueOf(idx);
        
        Integer q = xbt.find(index);
        if (q != null) {
            Set<PairInt> results = new HashSet<PairInt>();
            results.add(new PairInt(x, y));
            return results;
        }
                
        TIntSet closestIndexes = new TIntHashSet();
        
        double closestDist = Double.MAX_VALUE;
        
        Integer predecessor = null;
        Integer successor = null;
        
        if (useCache && pCache.containsKey(idx)) {
            predecessor = Integer.valueOf(pCache.get(idx));
        } else {
            //O(log_2(maxW))
            predecessor = xbt.predecessor(index);
            if (useCache && predecessor != null) {
                pCache.put(idx, predecessor.intValue());
            }
        }
        if (useCache && sCache.containsKey(idx)) {
            successor = Integer.valueOf(sCache.get(idx));
        } else {
            //O(log_2(maxW))
            successor = xbt.successor(index);
            if (useCache && successor != null) {
                sCache.put(idx, successor.intValue());
            }
        }
        
        double dp2 = dist(x, y, predecessor);
        double ds2 = dist(x, y, successor);
        if (dp2 <= ds2 && (dp2 <= dMax)) {
            closestDist = dp2;
            closestIndexes.add(predecessor.intValue());
            if (dp2 == ds2) {
                closestIndexes.add(successor.intValue());
            }
        } else if (ds2 < dp2 && (ds2 <= dMax)) {
            closestDist = ds2;
            closestIndexes.add(successor.intValue());
        }
        
        int goal = (closestDist != Double.MAX_VALUE) ?
            (int)Math.ceil(closestDist) : 0;
        
        if (goal > dMax) {
            goal = dMax;
        }
        
        int yLow = estimateLowBound(y, goal);
       
        int yCurrent;
        if (predecessor == null) {
            yCurrent = Integer.MIN_VALUE;
        } else {
            int pRow = getRow(predecessor.intValue());
            if (pRow < y) {
                yCurrent = pRow;
            } else {
                yCurrent = pRow - 1;
            }
        }
        
        // predecessor searches until reach yLow, adjusting goal by
        //   min distances
        Integer p2 = null; 
        Integer s2 = null;
        while (yCurrent >= yLow) {
            int cIdx = getInternalIndex(x, yCurrent);
            Integer cIndex = Integer.valueOf(cIdx);
            
            q = xbt.find(cIndex);
            if (q != null) {
                p2 = q;
                dp2 = dist(x, y, p2);
                ds2 = Double.MAX_VALUE;
            } else {
                if (useCache && pCache.containsKey(cIdx)) {
                    p2 = Integer.valueOf(pCache.get(cIdx));
                } else {
                    //O(log_2(maxW))
                    p2 = xbt.predecessor(cIndex);
                    if (useCache && p2 != null) {
                        pCache.put(cIdx, p2.intValue());
                    }
                }
                if (useCache && sCache.containsKey(cIdx)) {
                    s2 = Integer.valueOf(sCache.get(cIdx));
                } else {
                    //O(log_2(maxW))
                    s2 = xbt.successor(cIndex);
                    if (useCache && s2 != null) {
                        sCache.put(cIdx, s2.intValue());
                    }
                }
                dp2 = dist(x, y, p2);
                ds2 = dist(x, y, s2);
            }
            if ((dp2 < ds2) && (dp2 < closestDist) && (dp2 <= dMax)) {
                closestIndexes.clear();
                closestDist = dp2;
                closestIndexes.add(p2);
                goal = (int)Math.ceil(closestDist);
                if (goal > dMax) {
                    goal = dMax;
                }
                yLow = estimateLowBound(y, goal);                
            } else if ((ds2 < dp2) && (ds2 < closestDist) && (ds2 <= dMax)) {
                closestIndexes.clear();
                closestDist = ds2;
                closestIndexes.add(s2);
                goal = (int)Math.ceil(closestDist);
                if (goal > dMax) {
                    goal = dMax;
                }
                yLow = estimateLowBound(y, goal);
            } else if (dp2 == closestDist && (dp2 != Double.MAX_VALUE)) {
                closestIndexes.add(p2.intValue());
                if (dp2 == ds2) {
                    closestIndexes.add(s2.intValue());
                }
            } else if (ds2 == closestDist && (ds2 != Double.MAX_VALUE)) {
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
        if (successor == null) {
            yCurrent = Integer.MAX_VALUE;
        } else {
            int sr = getRow(successor.intValue());
            if (sr > y) {
                yCurrent = sr;
            } else {
                yCurrent = sr + 1;
            }
        }
        int yHigh = estimateHighBound(y, goal);
        
        while (yCurrent <= yHigh) {
            int cIdx = getInternalIndex(x, yCurrent);
            Integer cIndex = Integer.valueOf(cIdx);            
            q = xbt.find(cIndex);
            if (q != null) {
                p2 = q;
                dp2 = dist(x, y, p2);
                ds2 = Double.MAX_VALUE;
            } else {
                if (useCache && pCache.containsKey(cIdx)) {
                    p2 = Integer.valueOf(pCache.get(cIdx));
                } else {
                    //O(log_2(maxW))
                    p2 = xbt.predecessor(cIndex);
                    if (useCache && p2 != null) {
                        pCache.put(cIdx, p2.intValue());
                    }
                }
                if (useCache && sCache.containsKey(cIdx)) {
                    s2 = Integer.valueOf(sCache.get(cIdx));
                } else {
                    //O(log_2(maxW))
                    s2 = xbt.successor(cIndex);
                    if (useCache && s2 != null) {
                        sCache.put(cIdx, s2.intValue());
                    }
                }
                dp2 = dist(x, y, p2);
                ds2 = dist(x, y, s2);
            }
            if ((dp2 < ds2) && (dp2 < closestDist) && (dp2 <= dMax)) {
                closestIndexes.clear();
                closestDist = dp2;
                closestIndexes.add(p2);
                goal = (int)Math.ceil(closestDist);
                if (goal > dMax) {
                    goal = dMax;
                }
                yHigh = estimateHighBound(y, goal);
            } else if ((ds2 < dp2) && (ds2 < closestDist) && (ds2 <= dMax)) {
                closestIndexes.clear();
                closestDist = ds2;
                closestIndexes.add(s2);
                goal = (int)Math.ceil(closestDist);
                if (goal > dMax) {
                    goal = dMax;
                }
                yHigh = estimateHighBound(y, goal);
            } else if (dp2 == closestDist && (dp2 != Double.MAX_VALUE)) {
                closestIndexes.add(p2.intValue());
                if (dp2 == ds2) {
                    closestIndexes.add(s2.intValue());
                }
            } else if (ds2 == closestDist && (ds2 != Double.MAX_VALUE)) {
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
        
        if (high > height) {
            high = height;
        }
        
        return high;
    }
}
