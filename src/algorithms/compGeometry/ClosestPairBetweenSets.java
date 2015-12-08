package algorithms.compGeometry;

import algorithms.MergeSort;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Find the closest pair of points between two sets of points by using divide and 
 * conquer to reduce the number of comparisons
 * 
 * From pseudocode in Intro to Algorithms by Cormen et al.
 * 
 * runtime complexity was ~ O(N lg N), but changes have been made for set2
 * so need to recalculate after changes are finished.
 * 
 * NOT YET READY FOR USE
 * 
 * TODO: have to edit to keep set1 and set2 separate w/ same structures
   to divide them separately, then add the 2 sets for the brute force compare
 *
 * @author nichole
 */
public class ClosestPairBetweenSets {
    
    /**
     *
     */
    public ClosestPairBetweenSets() {
    }
    
    /**
     * find the closest pair within the two x,y sets of points.
     *
     * @param set1
     * @param set2
     * @return
     *     returns 5 items in array:  shortest distance between pairs, pt1 x, pt1 y, pt2 x, pt2 y
     */
    public ClosestPairInt findClosestPair(Set<PairInt> set1, Set<PairInt> set2) {
        //List<Float> xPoints, List<Float> yPoints) {
        
        if (set1 == null) {
	    	throw new IllegalArgumentException("set1 cannot be null");
        }
	    if (set2 == null) {
	    	throw new IllegalArgumentException("set2 cannot be null");
        }
                
        List<PairIntWithIndex> p = new ArrayList<PairIntWithIndex>();
        List<PairIntWithIndex> x = new ArrayList<PairIntWithIndex>();
        List<PairIntWithIndex> y = new ArrayList<PairIntWithIndex>();
        
//TODO: have to keep set1 and set2 separate w/ same structures below
// to divide them separately, then add the 2 sets for the brute force compare
        
        for (PairInt ps1 : set1) {
            // since the points xy are not going to be modified, can reuse xy
            PairIntWithIndex xy = new PairIntWithIndex(ps1.getX(), ps1.getY(), 1);
            p.add(xy);
            y.add(xy);
        }
        for (PairInt ps2 : set2) {
            // since the points xy are not going to be modified, can reuse xy
            PairIntWithIndex xy = new PairIntWithIndex(ps2.getX(), ps2.getY(), 2);
            p.add(xy);
            y.add(xy);
        }
        
        MergeSort.sortByXThenY(p);
        x.addAll(p);
        
        MergeSort.sortByYThenX(y);
        
        return divide(p, x, y);
    }
                
    /**
     * use divide and conquer to find the pair of points with the smallest 
     * separation.  note that x and y have to contain objects that are the
     * same instance as the objects in p.
     * 
     * @param p should be ordered left to right already
     * @param x
     * @param y
     * @return 
     */
    protected ClosestPairInt divide(List<PairIntWithIndex> p, List<PairIntWithIndex> x, 
        List<PairIntWithIndex> y) {
        
        if (p.size() <= 3) {
            return bruteForceMinDistance(p);
        }
                            
        int q = p.size() >> 1; // 0 1 *2*  3 4 5

        List<PairIntWithIndex> pL = new ArrayList<PairIntWithIndex>(q);
        List<PairIntWithIndex> pR = new ArrayList<PairIntWithIndex>(q);
        
        // need a datastructure to search for pL members.  for a hash table, search is O(1).
        //  will assume that the java HashSet implementation has O(1) insert and search
        Set<PairIntWithIndex> isInL = new HashSet<PairIntWithIndex>();
        
        List<PairIntWithIndex> xL = new ArrayList<PairIntWithIndex>(q);
        List<PairIntWithIndex> xR = new ArrayList<PairIntWithIndex>(q);
        List<PairIntWithIndex> yL = new ArrayList<PairIntWithIndex>(q);
        List<PairIntWithIndex> yR = new ArrayList<PairIntWithIndex>(q);
        
        // p is already sorted by x, so fill in the subsets for p and x
        for (int i = 0; i < q; ++i) {
            PairIntWithIndex pi = p.get(i);
            pL.add(pi);
            xL.add(pi);
            isInL.add(pi);
        }
        for (int i = q; i < p.size(); ++i) {
            PairIntWithIndex pi = p.get(i);
            pR.add(pi);
            xR.add(pi);
        }
        
        for (int i = 0; i < y.size(); ++i) {
            PairIntWithIndex yi = y.get(i);
            if (isInL.contains(yi)) {
                yL.add(yi);
            } else {
                yR.add(yi);
            }
        }
        
        // use 2 recursive calls for conquer

        // find closest pair in pL w/ pL, xL, and yL
        ClosestPairInt cpL = divide(pL, xL, yL);

        // find closest pair in pP w/ pR, xR, and yR
        ClosestPairInt cpR = divide(pR, xR, yR);

        ClosestPairInt d = (cpL.separationSq <= cpR.separationSq) ? cpL : cpR;
        
        // NOTE: temporary work around until have changed to keep set1 and set2 divides separate
        if ((d.point0 == null) && (pL.size() < 10) && (pR.size() < 10)) {
            if (hasAPoint(pL, 1) && hasAPoint(pR, 2)) {
                List<PairIntWithIndex> list = new ArrayList<PairIntWithIndex>();
                list.addAll(pL);
                list.addAll(pR);
                d = bruteForceMinDistance(list);
            } else if (hasAPoint(xL, 1) && hasAPoint(xR, 2)) {
                List<PairIntWithIndex> list = new ArrayList<PairIntWithIndex>();
                list.addAll(xL);
                list.addAll(xR);
                d = bruteForceMinDistance(list);
            } else if (hasAPoint(yL, 1) && hasAPoint(yR, 2)) {
                List<PairIntWithIndex> list = new ArrayList<PairIntWithIndex>();
                list.addAll(yL);
                list.addAll(yR);
                d = bruteForceMinDistance(list);
            }
        }
        
        return combine(yL, yR, d, xL.get(xL.size() - 1), xR.get(0));            
    }
    
    private boolean hasAPoint(List<PairIntWithIndex> list, int groupNumber) {
        for (PairIntWithIndex p : list) {
            if (p.getPixIndex() == groupNumber) {
                return true;
            }
        }
        return false;
    }
    
    /**
     * runtime complexity:
     *    O(N) plus a small fraction of points contributing O(m lg m) for a sort
     * @param yL
     * @param yR
     * @param d
     * @param leftMostL
     * @param rightMostR
     * @return 
     */
    protected ClosestPairInt combine(List<PairIntWithIndex> yL, List<PairIntWithIndex> yR, 
        ClosestPairInt d,  PairIntWithIndex leftMostL, PairIntWithIndex rightMostR) {
        
        List<PairIntWithIndex> yPrime = new ArrayList<PairIntWithIndex>();
        
        int delta = d.separationSq;
        
        // find the points in yL that have x within delta from leftMostL
        // traversing yL (rather than xL) to keep yPrime ordered
        for (int i = 0; i < yL.size(); ++i) {
            float xPoint = yL.get(i).getX();
            if ((leftMostL.getX() - xPoint) <= delta) {
                yPrime.add(yL.get(i));
            }
        }
        
        // find the points in yR that have x within delta from rightMostR
        // traversing yR (rather than xR) to keep yPrime ordered
        for (int i = 0; i < yR.size(); ++i) {
            float xPoint = yR.get(i).getX();
            if ((xPoint - rightMostR.getX()) <= delta) {
                yPrime.add(yR.get(i));
            }
        }
        
        int idx0 = -1;
        int idx1 = -1;
        int minDistSq = Integer.MAX_VALUE;
        
        for (int i = 0; i < yPrime.size(); i++) {
            
            PairIntWithIndex yi = yPrime.get(i);
            
            for (int j = (i + 1); j < (i + 8); j++) {
                if (j > (yPrime.size() - 1)) {
                    break;
                }
                
                PairIntWithIndex yj = yPrime.get(j);
                
                // only compare if the indexes of yi and yj are different
                if (yi.getPixIndex() == yj.getPixIndex()) {
                    continue;
                }
                
                int diffX = yi.getX() - yj.getX();
                int diffY = yi.getY() - yj.getY();
                
                int distSq = (diffX * diffX) + (diffY * diffY);
                if (distSq < minDistSq) {
                    minDistSq = distSq;
                    idx0 = i;
                    idx1 = j;
                }
            }
        }
       
        if (minDistSq == Integer.MAX_VALUE) {
            return d;
        }
                
        if (minDistSq < delta) {
            return new ClosestPairInt(yPrime.get(idx0), yPrime.get(idx1), minDistSq);
        } else {
            return d;
        }        
    }
    
    /**
     *
     * @param p
     * @return
     */
    protected ClosestPairInt bruteForceMinDistance(List<PairIntWithIndex> p) {
        
        if (p.size() < 1) {
            return new ClosestPairInt(null, null, Integer.MAX_VALUE);
        }
        
        int minDistSq = Integer.MAX_VALUE;
        int idx0 = -1;
        int idx1 = -1;

        for (int i = 0; i < p.size(); i++) {
            
            PairIntWithIndex pi = p.get(i);
            
            for (int j = i; j < p.size(); j++) {
                if (i == j) {
                    continue;
                }
                
                PairIntWithIndex pj = p.get(j);
                
                if (pi.getPixIndex() == pj.getPixIndex()) {
                    continue;
                }
                
                int diffX = pi.getX() - pj.getX();
                int diffY = pi.getY() - pj.getY();
                int distSq = (diffX * diffX) + (diffY * diffY);

                if (distSq < minDistSq) {
                    minDistSq = distSq;
                    idx0 = i;
                    idx1 = j;
                }
            }
        }
        
        if (idx0 == -1) {
            return new ClosestPairInt(null, null, Integer.MAX_VALUE);
        }
                
        return new ClosestPairInt(p.get(idx0), p.get(idx1), minDistSq);
    }

    /**
     *
     */
    public class ClosestPairInt {
        PairIntWithIndex point0;
        PairIntWithIndex point1;
        int separationSq;

        /**
         *
         * @param p0
         * @param p1
         * @param sep
         */
        public ClosestPairInt(PairIntWithIndex p0, PairIntWithIndex p1, int sep) {
            this.point0 = p0;
            this.point1 = p1;
            this.separationSq = sep;
        }
        public PairIntWithIndex getPoint0() {
            return point0;
        }
        public PairIntWithIndex getPoint1() {
            return point1;
        }
        public int getSeparationSquared() {
            return separationSq;
        }
    }
    
}
