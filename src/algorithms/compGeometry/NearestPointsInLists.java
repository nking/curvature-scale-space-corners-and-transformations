package algorithms.compGeometry;

import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Builds an unmodifiable data structure to make finding member points within
 * a distance of a location faster than O(N^2).
 * 
 * The runtime complexity is O(N*lg2(N)) for the constructor
 * and for search is O(lg2(N)) + small number of scans surrounding the nearest x.
 * 
 * @author nichole
 */
public class NearestPointsInLists {
    
    /*
    to make the inserts of new points faster, have switched from internal data
    structures from ordered lists to sorted tree maps
    
    at construction time:
       runtime complexity O(N*lg2(N))
    search:
       looks for the x key in the treemap and if not found, uses 
       the tree map's higherEntry or lowerEntry, depending upon use.
       iteration of the treemap once the existing key is found is done using
           subMap or descendingMap
    add point:
       the new point gets added to the tree map
    
    Note that if the data were extremely sparsely populated and the image
        had a very large width and small height, the ordered lists implementation
        might be a better implementation because of the space complexity difference.
    */
    
    private final TreeMap<Integer, TreeSet<Integer>> xy;
    private final Map<PairInt, Integer> listIndexes;
    private int minX = Integer.MAX_VALUE;
    private int maxX = Integer.MIN_VALUE;
    
    /**
     * constructor, runtime complexity is O(N*lg2(N))
     * @param pointsList 
     */
    public NearestPointsInLists(List<Set<PairInt>> pointsList) {
       
        if (pointsList == null) {
            throw new IllegalStateException("pointsList cannot be null");
        }
        
        xy = new TreeMap<Integer, TreeSet<Integer>>();
        
        listIndexes = new HashMap<PairInt, Integer>();
        
        int n = pointsList.size();
        
        for (int i = 0; i < n; ++i) {
            
            Integer index = Integer.valueOf(i);
            
            Set<PairInt> points = pointsList.get(i);
            
            for (PairInt p : points) {
                
                Integer x = Integer.valueOf(p.getX());
                Integer y = Integer.valueOf(p.getY());
                
                TreeSet<Integer> ySet = xy.get(x);
                if (ySet == null) {
                    ySet = new TreeSet<Integer>();
                    xy.put(x, ySet);
                }
                ySet.add(y);
                
                listIndexes.put(p.copy(), index);
                
                if (x < minX) {
                    minX = x;
                }
                if (x > maxX) {
                    maxX = x;
                }
            }
        }        
    }
    
    /**
     * add a point to the instance.
     * runtime complexity is O(lg2(N)).
     * @param p
     * @param listIndex 
     */
    public void addPoint(PairInt p, int listIndex) {
        
        int x = p.getX();
        
        Integer xKey = Integer.valueOf(x);
        
        Integer yKey = Integer.valueOf(p.getY());
       
        TreeSet<Integer> ySet = xy.get(xKey);
        if (ySet == null) {
            ySet = new TreeSet<Integer>();
            xy.put(xKey, ySet);
        }
        ySet.add(yKey);

        listIndexes.put(p.copy(), Integer.valueOf(listIndex));
        
        if (x < minX) {
            minX = x;
        }
        if (x > maxX) {
            maxX = x;
        }
    }
    
    /**
     * find points within radius of (xCenter, yCenter) in the contained points.
     * runtime complexity is a little more than O(N*lg2(N)) and much less than
     * O(N^2).
     * @param xCenter
     * @param yCenter
     * @param radius
     * @return 
     */
    public PairInt findClosest(int xCenter, int yCenter, float radius) {
        
        Map<Integer, PairInt> result = findNeighbors(xCenter, yCenter, radius);
        
        PairInt minDistP = null;
        int minDistSq = Integer.MAX_VALUE;
        
        for (Entry<Integer, PairInt> entry : result.entrySet()) {
            PairInt p = entry.getValue();
            int diffX = p.getX() - xCenter;
            int diffY = p.getY() - yCenter;
            int distSq = (diffX*diffX) + (diffY*diffY);
            if (distSq < minDistSq) {
                minDistSq = distSq;
                minDistP = p;
            }
        }
        
        return minDistP;
    }
    
    /**
     * find points within radius of (xCenter, yCenter) in the contained points.
     * runtime complexity is a little more than O(N*lg2(N)) and much less than
     * O(N^2).
     * @param xCenter
     * @param yCenter
     * @param radius
     * @return 
     */
    public Map<Integer, PairInt> findNeighbors(int xCenter, int yCenter, float radius) {
        
        Map<Integer, PairInt> result = new HashMap<Integer, PairInt>();
        Map<Integer, Integer> resultSqDistances = new HashMap<Integer, Integer>();
        
        if (xy.size() == 0) {
            return result;
        }
        
        Integer yKey = Integer.valueOf(yCenter);
        
        //collect smallest distances within r for a listIndex
        
        double rSq = Math.sqrt(2) * radius * radius;
        
        if (xCenter <= minX) {
            
            // search from first key forward
            
            Integer xKey = xy.firstKey();
            
            NavigableMap<Integer, TreeSet<Integer>> navXMap 
                = xy.tailMap(xKey, true);
            
            searchX(xCenter, yCenter, xKey, yKey, navXMap, result, 
                resultSqDistances, rSq);
            
        } else if (xCenter >= maxX) {
            
            // search from last key backward
            
            Integer xKey = xy.lastKey();
            
            NavigableMap<Integer, TreeSet<Integer>> navXMap 
                = xy.headMap(xKey, true).descendingMap();
            
            searchX(xCenter, yCenter, xKey, yKey, navXMap, result, 
                resultSqDistances, rSq);
            
        } else {
            
            Integer xKey = Integer.valueOf(xCenter);
                                
            // --- search backwards ----
            if (xy.get(xKey) == null) {
                Entry<Integer, TreeSet<Integer>> entry = xy.lowerEntry(xKey);
                assert(entry != null);
                xKey = entry.getKey();
            }
            
            NavigableMap<Integer, TreeSet<Integer>> navXMap 
                = xy.headMap(xKey, true).descendingMap();
            
            searchX(xCenter, yCenter, xKey, yKey, navXMap, result, 
                resultSqDistances, rSq);
                               
            // --- search forward ----
            xKey = Integer.valueOf(xCenter);
            if (xy.get(xKey) == null) {
                Entry<Integer, TreeSet<Integer>> entry = xy.higherEntry(xKey);
                assert(entry != null);
                xKey = entry.getKey();
            }
            
            navXMap = xy.tailMap(xKey, true);
            
            searchX(xCenter, yCenter, xKey, yKey, navXMap, result, 
                resultSqDistances, rSq);
        }
        
        return result;
    }

    private void searchY(int yCenter, Integer x, NavigableSet<Integer> navYSet, 
        Map<Integer, PairInt> result, Map<Integer, Integer> resultSqDistances, 
        int diffX, double rSq) {
        
        for (Integer y : navYSet) {
            
            int diffY = Math.abs(y.intValue() - yCenter);
            
            int distSq = (diffX * diffX) + (diffY * diffY);
            if (distSq > rSq) {
                break;
            }
            
            PairInt p = new PairInt(x.intValue(), y.intValue());
            
            Integer listIndex = listIndexes.get(p);
            
            assert(listIndex != null);

            if ((resultSqDistances.get(listIndex) == null)
                || (distSq < resultSqDistances.get(listIndex))) {
                result.put(listIndex, p);
                resultSqDistances.put(listIndex, Integer.valueOf(distSq));
            }
        }        
    }

    private void searchX(int xCenter, int yCenter, Integer xKey, Integer yKey, 
        NavigableMap<Integer, TreeSet<Integer>> navXMap, 
        Map<Integer, PairInt> result, Map<Integer, Integer> resultSqDistances, 
        double rSq) {
        
        for (Entry<Integer, TreeSet<Integer>> entry : navXMap.entrySet()) {
            Integer x = entry.getKey();
            int diffX = Math.abs(x.intValue() - xCenter);
            if (diffX > rSq) {
                break;
            }
            TreeSet<Integer> ySet = entry.getValue();

            // search backwards in y
            Integer y0 = ySet.contains(yKey) ? yKey : ySet.lower(yKey);
            if (y0 != null) {
                
                NavigableSet<Integer> navYSet =
                    ySet.headSet(y0, true).descendingSet();

                searchY(yCenter, x, navYSet, result, resultSqDistances, diffX,
                    rSq);
            }
            
            // search forward in y
            y0 = ySet.higher(yKey);
            if (y0 != null) {
                
                NavigableSet<Integer> navYSet = ySet.tailSet(y0, true);

                searchY(yCenter, x, navYSet, result, resultSqDistances, diffX,
                    rSq);
            }
        }
    }

}
