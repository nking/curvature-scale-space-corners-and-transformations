package algorithms.graphs;

import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;

/**
 * a class to hold a set of points within a region and a set of points
 * which are on the perimeter of the region.
 * @author nichole
 */
public class Region {
    
    //NOTE: may change the internal data to use more compact strcutres in the future
    
    private final Set<PairInt> points;
    
    private final Set<PairInt> perimeter;
    
    /**
     * constructor keeps the instance it is given, so copy it if needed before
     * using this method
     * @param thePoints 
     */
    public Region(Set<PairInt> thePoints) {
        
        points = thePoints;
        
        perimeter = findPerimeter(points);
        
    }
    
    /**
     * add other region data into this one and clear the other region data.
     * 
     * @param other 
     */
    void mergeIntoThis(Region other) {
        
        this.points.addAll(other.points);
        
        // merge perimeters
        this.perimeter.addAll(other.perimeter);
        
        // remove any in perimeer w/ 8 neighbors
        Set<PairInt> rm = new HashSet<PairInt>();
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        for (PairInt p : this.perimeter) {
            int x = p.getX();
            int y = p.getY();
            int nN = 0;
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                PairInt p2 = new PairInt(x2, y2);
                if (points.contains(p2)) {
                    ++nN;
                }
            }
            if (nN == 8) {
                rm.add(p);
            }
        }
        this.perimeter.removeAll(rm);
        
        other.points.clear();
        other.perimeter.clear();
    }

    private Set<PairInt> findPerimeter(Set<PairInt> set0) {
        
        Set<PairInt> set = new HashSet<PairInt>();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        for (PairInt p : set0) {
            int x = p.getX();
            int y = p.getY();
            int nN = 0;
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                PairInt p2 = new PairInt(x2, y2);
                if (set0.contains(p2)) {
                    ++nN;
                }
            }
            if (nN < 8) {
                set.add(p);
            }
        }
        
        return set;
    }
    
    public int size() {
        return points.size();
    }
    
    public boolean contains(PairInt p) {
        return points.contains(p);
    }
    
    public boolean contains(int x, int y) {
        return points.contains(new PairInt(x, y));
    }
    
    public boolean perimeterContains(PairInt p) {
        return perimeter.contains(p);
    }
    
    public boolean perimeterContains(int x, int y) {
        return perimeter.contains(new PairInt(x, y));
    }

    /**
     * @return the perimeter
     */
    public Set<PairInt> getPerimeter() {
        return perimeter;
    }

    /**
     * @return the points
     */
    public Set<PairInt> getPoints() {
        return points;
    }
}
