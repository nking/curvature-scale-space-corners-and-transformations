package algorithms.compGeometry;

import algorithms.imageProcessing.Heap;
import algorithms.imageProcessing.HeapNode;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntDoubleMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Logger;

/**
   A class to create a 2D medial axis given shape points
   and the boundary.
       from https://en.wikipedia.org/wiki/Medial_axis
       "the medial axis of a subset S which is bounded 
        by planar curve C is the locus of the centers 
        of circles that are tangent to curve C in two 
        or more points."

   This class is an implementation of 
   "Efficient and Robust Computation of an Approximated Medial Axis"
   by Yang, Brock, and Moll
   
 * which requires
 * "Efficient Computation of A Simplified Medial Axis"
 * by Foskey, Lin, and Manocha
 * and 
 * "A Framework for Using the Workspace Medial Axis in 
 * PRM Planners"
 * by Holleman and Kavraki
 * 
 * @author nichole
 */
public class MedialAxis1 {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    private final Set<PairInt> points;
    private final Set<PairInt> boundary;
    private final NearestNeighbor2D np;
    
    /** as circles within points are searched, they are placed
    in processed
    */
    private final Set<PairInt> processed;
    
    private final List<MedialAxisPoint> medAxisList;
    
    //xMin, xMax, yMin, yMax
    private final int[] minMaxXY;
    
    /**
     * separation angle used for sampling points in
     * a circle around a point.
     */
    private static final double sepAng = Math.PI/3;
    private final int nSampl = 6;//  2*pi/sepAng
    
    
    //private static final double sepAng = Math.PI/6;
    //private final int nSampl = 12;//  2*pi/sepAng
    
    //private static final double sepAng = Math.PI/9;
    //private final int nSampl = 18;//  2*pi/sepAng
    
    //private static final double sepAng = Math.PI/18;
    //private final int nSampl = 36;//  2*pi/sepAng
    
    // 10 degrees threshold for separation angle criterion
    // ~ 0.1745
    //private final double threshold = Math.PI/18;
    
    private double twoPI = Math.PI;
    private double sinePiDivN = Math.sin(Math.PI/nSampl);
    
    private final int nInterior;
    
    /**
     * constructor containing all points in the area
     * and the bounding points.
     * (NOTE: in future may make a version that doesn't
     * need all interior points and uses an adaptive
     * sampling, but use case for now fits these arguments
     * better.).
     * 
     * @param shapePoints
     * @param boundaryPoints 
     */
    public MedialAxis1(final Set<PairInt> shapePoints,
        final Set<PairInt> boundaryPoints) {

        this.points = new HashSet<PairInt>(shapePoints);
        
        this.boundary = new HashSet<PairInt>(boundaryPoints);
    
        points.removeAll(boundary);
        
        minMaxXY = MiscMath.findMinMaxXY(boundary);
 
        nInterior = points.size();
        
        this.np = new NearestNeighbor2D(boundary, 
            minMaxXY[1], minMaxXY[3]);
        
        processed = new HashSet<PairInt>(points.size());
        
        medAxisList = new ArrayList<MedialAxisPoint>();
    }
    
    protected Set<PairInt> getNearestBoundaryPoints(PairInt p) {
        return np.findClosest(p.getX(), p.getY());
    }
    
    protected double getSepAng() {
        return sepAng;
    }
    protected int getNSampl() {
        return nSampl;
    }

    /**
     * add mp to the heap, calculating a node key for it
     * that results in ordering by decreasing radius.
     * @param heap
     * @param mp 
     */    
    protected void addToHeap(Heap heap, MedialAxisPoint mp) {
        // constructing key out of the inverse of the distance
        // from medial axis point to it's nearest boundary point,
        // where max is max of maxx and maxy.
        // key is a long so will use a factor equal to max
        
        long max = Math.max(minMaxXY[1], minMaxXY[3]);
        
        double r = mp.getBoundaryPoints()[0].getDistance();
                
        double invR = 1./r;
        
        long key = (long)Math.ceil(max * invR);
        
        HeapNode node = new HeapNode(key);
        node.setData(mp);
        
        heap.insert(node);
    }
    
    protected boolean assertPointTotals() {
        int nP1 = points.size();
        int nP2 = processed.size();
        if ((nP1 + nP2) == nInterior) {
            return true;
        }
        return false;
    }
    
    protected boolean assertUniqueMedialAxesPoints() {
        Set<PairInt> mAPs = new HashSet<PairInt>();
        for (MedialAxis1.MedialAxisPoint mp : medAxisList) {
            PairInt pp = mp.getCenter();
            if (mAPs.contains(pp)) {
                return false;
            }
            mAPs.add(pp);
        }
        return true;
    }
    
    /**
     * find medial axis.  The results can be retrieved
     * with getMedialAxisPoints().
     * NOTE: the method is not yet fully tested so
     * resolution arguments are not yet available for use.
     */
    public void findMedialAxis() {
        
        medAxisList.clear();
        
        // max heap ordered by largest radius
        Heap q = new Heap();
        
        Set<PairInt> addedM = new HashSet<PairInt>();
        Set<PairInt> removed = null;
        
        // also adds point(s) to the heap
        MedialAxisResults results = 
            findAndProcessFirstPoint(q, addedM);
        
        int[] xSurf = new int[nSampl];
        int[] ySurf = new int[nSampl];
        
        Set<PairInt> visited = new HashSet<PairInt>();
        
        int nIter = 0;
        
        while (!q.isEmpty() && (processed.size() < nInterior)) {

            nIter++;
            
            HeapNode node = q.extractMin();
            
            MedialAxisPoint mp = (MedialAxisPoint) node.getData();
            
            /*
            for each medial axis point, search for other
            medial axis points within a radius defined as
            the distance to the nearest boundary point
            from mp.  the radius is used to sample the 
            circumference of a sphere at that distance.
            */
            
            PairInt p = mp.getCenter();
            
            if (visited.contains(p)) {
                continue;
            }
            visited.add(p);
  
            // for each p, sample points at a radius r
            // and search for medial axis points around
            // those.  (note that searching p directly for
            // medial axis points is a less detailed complete
            // result than creating p2's and searching those.)
            
            double r = mp.getBoundaryPoints()[0].getDistance();
            assert(r != Double.MIN_VALUE);
            int nSPoints = populateSurfacePoints(p, r,
                xSurf, ySurf);
                        
            if (nSPoints == 0) {
                removed = subtractFromPoints(p, r - 1);
                processed.addAll(removed);
                assert (assertPointTotals());

                if (q.isEmpty() && (points.size() > 0)) {
                    // choose another point from points and continue
                    // also adds point(s) to the heap
                    results = findAndProcessARemainingPoint(q, addedM);
                }

                continue;
            }
            
            log.fine("nIter=" + nIter + " q.size=" 
                + q.getNumberOfNodes() + 
                " s.size=" + processed.size()
                + " nInterior=" + nInterior);

            for (int i = 0; i < nSPoints; ++i) {
                
                PairInt p2 = new PairInt(xSurf[i], ySurf[i]);
                
                if (visited.contains(p2)) {
                    continue;
                }
                
                results = findMedialAxesNearPoint(p2);
                if (results == null) {
                    PairInt bPoint = np.findClosest(p2.getX(), p2.getY())
                        .iterator().next();
                    double d = distance(bPoint.getX(), bPoint.getY(), p);
                    processed.addAll(subtractFromPoints(p2, d - 1));                    
                    assert(assertPointTotals());
                    
                    continue;
                }
                
                // add to data structures, as above
                for (MedialAxisPoint mp2 : results.medialAxes) {
                    PairInt pp = mp2.getCenter();
                    if (addedM.contains(pp)) {
                        continue;
                    }
                    addToHeap(q, mp2);
                }
                
                searchCenterPoint(results, addedM);
        
                for (MedialAxisPoint m : results.medialAxes) {
                    PairInt pp = m.getCenter();
                    if (addedM.contains(pp)) {
                        continue;
                    }
                    addedM.add(pp);
                    medAxisList.add(m);
                }
                         
                processed.addAll(subtractFromPoints(results.center,
                    results.centerSrchR - 1));
                
                assert(assertPointTotals());
                assert(assertUniqueMedialAxesPoints());
            }
             
            assert(assertUniqueMedialAxesPoints());
        
            log.fine("q.n=" + q.getNumberOfNodes() +
                " points.size=" + points.size());
            
            if (q.isEmpty() && (points.size() > 0)) {
                // choose another point from points and continue
                // also adds point(s) to the heap
                results = findAndProcessARemainingPoint(q, addedM);
            }
        }
        
        log.fine("med axies nIter=" + nIter);
        
        assert(assertUniqueMedialAxesPoints());

        addedM = null;
        
        // iterate over medial axis points to refine centers.
        // dither should be defined by the tolerance tol
        Set<PairInt> srch = refineCentersOfMedAxisList();
        
        assert(assertPointTotals());
        assert(assertUniqueMedialAxesPoints());

        // ----- code to fill in gaps using a DFS search --- 
        fillInGaps(srch);
        
        assert(assertPointTotals());
        assert(assertUniqueMedialAxesPoints());
    }

    protected List<MedialAxisPoint> getMedAxisList() {
        return medAxisList;
    }
    
    /**
     * get the calculated medial axis point as a set of
     * unordered points.
     * @return 
     */
    public Set<PairInt> getMedialAxisPoints() {
        Set<PairInt> mAPs = new HashSet<PairInt>();
        for (MedialAxis1.MedialAxisPoint mp : medAxisList) {
            mAPs.add(mp.getCenter());
        }
        return mAPs;
    }

    /**
     * find an initial point to begin the search.
     * The point has a circle of radius = distance
     * to nearest boundary point.  That circle
     * intersects with at least one medial axis 
     * point.
     * @return 
     */
    private MedialAxisResults findInitialPoint() {
    
        PairInt p = points.iterator().next();
        
        return findMedialAxesNearPoint(p);
    }
    
    /**
     * Find a point which has a circle of radius = distance
     * to nearest boundary point, and that circle
     * intersects with at least one medial axis 
     * point. Note that if the original point does not
     * have medial axis point in the points sampled at
     * that radius, a search pattern is used to change
     * the point center until a point meeting the goals
     * is found or until all search variants have failed.
     * Failure can occur when the instance set "processed"
     * contains most of the points and the instance set
     * "points" is nearly empty.
     * 
     * @return 
     */
    protected MedialAxisResults findMedialAxesNearPoint(PairInt p) {
    
        log.fine("++srch for medial axis near p=" + p);
        
        List<MedialAxisPoint> medialAxes 
            = new ArrayList<MedialAxisPoint>();
        
        PairInt[] closestB = findNearestBoundsAsArray(
            p.getX(), p.getY());
        
        assert(closestB.length > 0);
       
        int status = 0;
        
        // find a p that results in valid medial axis points,
        // and update the closestB for it
        
        status = intersectsMedialAxis(closestB, p, medialAxes);

        if (status == -3) {

            throw new IllegalStateException("Error in algorithm:"
                + " could not find a nearest boundary for " +
                p.toString());

        } else if (status == -2) {

            //status==-2: no medial axis angle larger than threshold was found

            PairInt origP = p.copy();
           
            /*
            NOTES from paper:
            may need to revisit this while testing:
            If the new medial axis segment is significantly
            different from normal to the previous,
            then start 2 searches, 1 perpendicular to the previous
            medial axis center and 1 perp to the current medial
            axis.  the points should be away from the
            smallest angle between the old and new segments.            
            */
            
            log.fine("++STATUS=-2 for p=" + p);
           
            int cIdx = 0;
            
            // in case the algorithm fails and has a dx=0 or dy=0
            // this allows a retry and override of the zero variable
            boolean overrideDXY = false;
            
            while (status != 1 && (cIdx < closestB.length)) {
                
                medialAxes.clear();
 
                PairInt bPoint = closestB[cIdx];

                double r = distance(bPoint.getX(), bPoint.getY(), p);

                /*from bPoint to p, scale offsets by new r:
                 diffX = (p.x - b.x)/r
                 diffY = (p.y - b.y)/r
                 p2.x = p.x + (diffX * r2)
                 p2.y = p.y + (diffY * r2)
                */
                double diffX = (p.getX() - bPoint.getX())/r;
                double diffY = (p.getY() - bPoint.getY())/r;
                
                if (cIdx == 0) {
                    if (overrideDXY) {
                        if (diffX == 0) {
                            diffX = (p.getX() - bPoint.getX()) < 0 ?
                                1 : -1;
                        } else {
                            diffY = (p.getY() - bPoint.getY()) < 0 ?
                                1 : -1;
                        }
                    }
                }
                
                cIdx++;

                log.fine("++origP=" + origP + " nbp[i]=" + bPoint + " r=" + r
                    + " dx=" + diffX + " dy=" + diffY);

                int low, high;
                if (status == -2) {
                    low = (int)Math.round(r);
                    // guestimate w/ half of xMax or yMax
                    high = (int)Math.round(0.5 * 
                        Math.max(minMaxXY[1], minMaxXY[3]));
                } else {
                    low = 1;
                    high = (int)Math.round(r);
                }
                int mid = (int)Math.round(0.5 * (high + low));

                log.fine("++init:  mid=" + mid + " low=" + low + " high=" + high);

                // choose another circle center by choosing the point
                // on the circle around p which is opposite the nearest
                // boundary point.

                while (low < high) {
                    mid = (int)Math.round(0.5 * (high + low));
                    // create new p from mid distance from b along path to p
                    int x2 = p.getX() + (int)Math.round(diffX * mid);
                    int y2 = p.getY() + (int)Math.round(diffY * mid);
                    PairInt p2 = new PairInt(x2, y2);
                    
                    log.fine("++   p2=" + p2 + " mid=" + mid + " low=" + low + " high=" + high);
             
                    //TODO: consider removing if processed here
                    if ((points.contains(p2) || processed.contains(p2))) {
                        
                        // check if point intersects medial axis
                        // and if so, exit loop
                        PairInt[] closestB2 = findNearestBoundsAsArray(p2.getX(), p2.getY());
                        status = intersectsMedialAxis(closestB2, p2, medialAxes);
                        if (status == 1) {
                            p = p2;
                            closestB = closestB2;
                            break;
                        }
                    }
                    if (high == mid) {
                        high -= 1;
                    } else {
                        high = mid;
                    }
                } // end loop low<high
            
                // if srch unsuccessful, throw error or retry
                if (status != 1) {
                    if (!overrideDXY && (diffX == 0 || diffY == 0)) {
                        cIdx = 0;
                        overrideDXY = true;
                    } else {
                        log.fine("++Possible Error in algorithm:"
                        + " could not find a valid medial axis"
                        + " point from the random first point.");
                        //NOTE: in one case, p2 was a valid med axis
                        // pt but was already processed.
                        return null;
                    }
                }
                // end of binary search to change circle center
            }
        }

        MedialAxisResults results = new MedialAxisResults();
        results.medialAxes = medialAxes;
        results.center = p;
        results.centerSrchR = 
            distance(p.getX(), p.getY(), closestB[0]);
        
        return results;
        
        /*
        From "Efficient Computation of A Simplified Medial Axis"
        by Foskey, Lin, and Manocha.
        
        separation angle is the maximum angle formed 
        by the vectors connecting the medial axis 
        point x to its closest points on the boundary,
        and portions of the medial axis with a larger 
        separation angle tend to be more stable. 
        
        If there is more than two nearest points,
        then the maximum of the separation angles is used
        
        N_x is neighbor direction field.
        N_x = the direction of vector difference
            as (boundaryP - medialP vectors)/|(boundaryP - medialP)|

        two points medialP1 and medialP2 are on opposite
        sides of a medial axis sheet if N_x(medialP1) and
        N_x(medialP2) diverge.
        
        for a pair of points, calc angle between direction
        vectors N_x.
        if angle is < threshold, reject the pair.
        else
            check for false positives
                for example:
             reflex vertex for the nearest neighbor boundary point.
                both direction vectors converge towards vertex.
                if points are close to the vertex,
                the angle between them, is large even
                though the 2 points are not on oppos sides
                of a medial axis.  Tp avoid this error,
                check that the distance between the heads of
                the vectors is at least as large as the tails
                of the vectors.

                medialP1         medialP2
                   * tail of vectors  *

                           * head of vectors for false medial axis points due to reflex.
                           /\
               ___________/  \__________
        
 ===> this check for reflex points looks like it has to be
      done after several or all medial axis points have
      been found, not just as two are added.
        
        two different spatial subdivision schemes can be used:
         -- uniform grid
            simplest and is fast using distance fields.
            3D is handled as a slice for each z, and each 2D is handled:
            the area is divided in axis aligned grid.
            the distance field is scalar function: D_X of real numbers.
            it's determined by lower envelope (or minima) of all D_X_i.
         -- adaptive grid
            - an octree, for example.
         -- I'll use my static nearest neighbors
         whichever scheme is used, once a pair is found
         that satisfies the separation criteria,
         add it to the medial axis
         (need to design the medial axis and nodes.
         a node needs its coordinates, vectors from boundary
            and possibly the parent node which is the node
            adjacent to it that this one is derived from)
        
         More on the separation angle criteria:
            the distance from the candidate medial
            axis point to the nearest boundary point 
            should be much larger than the separation 
            of the adjacent boundary points...       

        -----------
        
    identifying medial axis points:
        U = uniformly distributed sample points on surface of sphere
        
        separation angle criteria:
            For a point on the sphere, the angle is defined
            as the angle between the nearest point on the 
            boundary and the nearest point on the other boundary.
            -- the authors here use 2 points on the sphere
               instead and the third point is the closest 
               point for them on the boundary
               (which should be the same point for them).
               the angle subtended by the third point
               is maximum for the point on the medial axis.
               if that angle is larger than a threshold
               theta t, the midpoint between them becomes
               the aMA point.
    
    note:
       for uniformly sampled points,
       the angle delta theta between 2 neighboring
       points (e.g. q1p1, q1p2) is 2*pi/N.
       and the distance between them is 2 * r * sin(pi/N).
       
       The absolute error is then eps_a = r * sin(pi/N).
    
       reversing that, one gets the number of points on
       a sphere: N = pi/(arcsin(eps_a/r))
    
       The relative error needs to locate the closest
       point on the sphere to the medial axis, 
       which can name q1 and m, respectively.
          eps_r = delta(q1) * sin(pi/N)
    
       with triangle inequality, have 
       eps_a/(delta(q1)+eps_a) < eps_a/delta(m) <= eps_a/delta(q1)
    
       if r <= delta(q1), can select N with:
          N = pi/arcsin(eps_r)
       else if r >= delta(q1), can select N with:
          N = pi/arcsin(eps_r * r /delta(q1))
    
       Given N, origin of sphere, and r, can calculate
         p_i x: origin_x + r * cos(2*pi*i/N)
             y: origin_y + r * sin(2*pi*i/N)
    
    
       Given an absolute error eps_a, 
       the algorithm will stop if the maximum delta(p) 
       of all sample points p is smaller than eps_a. 
       For a relative error eps_r, the algorithm will 
       terminate if the radius of all spheres is 
       smaller than the threshold Ke. 
       For this case using relative error, 
       the size of the largest "missed area" gate which 
       depends on the amount of local free space.
      */
   
    }
    
    /**
     * subtract the points in set "points" found within
     * a radius of center from set.
     * @param center
     * @param radius
     * @return 
     */
    protected Set<PairInt> subtractFromPoints(PairInt center,
        double radius) {
        
        Set<PairInt> output = new HashSet<PairInt>();
        
        extractFromPoints(center, radius, output);
        
        points.removeAll(output);
        
        return output;
    }
    
    /**
     * extract points from instance set "points" that are
     * within radius of center.
     * NOTE: that this class is still being tested.
     * For very complex shapes, this might be removing
     * points that need to searched that were not sampled
     * in the uniform sampling model (might need to use
     * the adaptive sampling).
     * 
     * @param center
     * @param radius
     * @param output 
     */
    protected void extractFromPoints(PairInt center,
        double radius, Set<PairInt> output) {
         
        Stack<PairInt> stack = new Stack<PairInt>();
        stack.add(center);
        if (points.contains(center)) {
            output.add(center);
        }
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        Set<PairInt> visited = new HashSet<PairInt>();
        
        while (!stack.isEmpty()) {
            PairInt p = stack.pop();
            if (visited.contains(p)) {
                continue;
            }
            visited.add(p);
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = p.getX() + dxs[k];
                int y2 = p.getY() + dys[k];
                if (x2 < minMaxXY[0] || y2 < minMaxXY[2] || (x2 > minMaxXY[1])
                    || (y2 > minMaxXY[3])) {
                    continue;
                }
                PairInt p2 = new PairInt(x2, y2);
                if (!points.contains(p2)) {
                    continue;
                }
                double dist = distance(x2, y2, center);
                if (dist <= radius) {
                    output.add(p2);
                    stack.add(p2);
                }
            }
        }        
    }
    
    /*
    for uniformly sampled points,
       the angle delta theta between 2 neighboring
       points (e.g. q1p1, q1p2) is 2*pi/N.
       and the distance between them is 2 * r * sin(pi/N).
       
       The absolute error is then eps_a = r * sin(pi/N).
    
       reversing that, one gets the number of points on
       a sphere: N = pi/(arcsin(eps_a/r))
    */
    
    /*
    input:  
       perimeter points
       interior points
    
    -- first primitive operation:
        identifies an initial point m and associated 
        distance delta(m), such that the resulting 
        sphere of radius delta(m) around m intersects 
        the medial axis. 
    
    -- second primitive operation: 
        given a solid D, a set of points P in the 
        interior of D and their direction vectors, 
        for each point p ∈ P, this method identifies 
        that point on the medial axis of D which 
        is closest to p.
    
    These primitives will be described after we 
    have detailed the algorithm of computing the aMA.
Assume point m lies on the medial axis and is 
    distance δ(m) away from the closest obstacle. 
    This point is determined using aforementioned 
    primitive. A priority queue Q is initialized
    to contain the sphere de- scribed by point m 
    and radius δ(m). The set S of spheres describing 
    the free space inside the solid D is initialized 
    to be the empty set.
    
    
    */
    
    /**
     * populate surfaceX and surfaceY with nSampl number of
     * points that are a distance r from center p.
     * Note that redundant points are removed and any points
     * that are in the instance set "processed".
     * The number of usable items in surfaceX and surfaceY
     * are the return argument.
     * @param p
     * @param r
     * @param surfaceX an empty array of size nSampl
     * @param surfaceY an empty array of size nSampl
     * @return 
     */
    protected int populateSurfacePoints(PairInt p, double r,
        int[] surfaceX, int[] surfaceY) {
        
        int ns = 0;
        for (int i = 0; i < surfaceX.length; ++i) {
            double t = 2. * twoPI * (double)i/(double)nSampl;
            int x1 = p.getX() + (int)Math.round(r * Math.cos(t));
            int y1 = p.getY() + (int)Math.round(r * Math.sin(t));
            if (x1 < minMaxXY[0] || (x1 > minMaxXY[1]) ||
                y1 < minMaxXY[2] || (y1 > minMaxXY[3])) {
                continue;
            }
            if ((ns > 0) && (x1 == surfaceX[ns - 1]) && (y1 ==
                surfaceY[ns - 1])) {
                continue;
            } else if ((i == (nSampl - 1)) && (ns > 0) && 
                (x1 == surfaceX[0]) && (y1 == surfaceY[0])) {
                continue;
            }
            // also, do not include space already searched in "processed"
            PairInt sp = new PairInt(x1, y1);
            if (processed.contains(sp) || boundary.contains(sp)) {
                continue;
            }
            
            surfaceX[ns] = x1;
            surfaceY[ns] = y1;
            
            ns++;
        }
        
        return ns;
    }
    
    /**
     * finds the medial axis points if p intersects medial axis.
     * @param nearestBounds
     * @param p
     * @param output 
     * @return code for result: 
     * -3 means there are no nearest boundary points;
     * -2 means no medial axis angle larger than threshold was found;
     * 1 means successfully found the medial axis points
     * and placed them in output.
     */
    protected int intersectsMedialAxis(PairInt[] nearestBounds, 
        PairInt p, List<MedialAxisPoint> output) {
        
        if (nearestBounds == null || nearestBounds.length == 0) {
            return -3;
        }
        
        output.clear();
        
        PairInt bPoint = nearestBounds[0];
        
        double r = distance(bPoint.getX(), bPoint.getY(), p);
        
        // might need to increase this to 2*
        int tol = (int)Math.ceil(r * sinePiDivN);
        double threshold = 0.5 * sinePiDivN;
     
        log.fine("p=" + p.toString() + " a nearest boundary="
            + bPoint.toString() + " r=" + r 
            + String.format(" tol=%d pix  sinePiDivN=%.4f", 
                tol, sinePiDivN));       
        
        // create uniformly sampled circles
        // and check that all points are within bounds
        
        int[] surfaceX = new int[nSampl];
        int[] surfaceY = new int[nSampl];
        int ns = populateSurfacePoints(p, r, surfaceX, surfaceY);
        
        /*
        nearestBounds are all equidistant, so can
        choose one of them as the corner of a triangle.
        
        for each pair of circle surface points, measuring
        the angle subtended by them and the boundary point
        to find the maximum angle above threshold.
        */
        
        // there may be more than one medial axis point,
        // for example, the circle enscribed in an open
        // triangle has one near the vertex side of circle
        // and on opposite side
        TIntList surfIdxes1 = new TIntArrayList();
        TIntList surfIdxes2 = new TIntArrayList();
        // tracking PairInt[] combinedNBP w/ key = index to surfIdxes1
        TIntObjectMap<PairInt[]> surfNBMap = new TIntObjectHashMap<PairInt[]>();
        
        // keeping a map of indexes and angles
        // in case adjacent pair passed filters.  keep latgest angle
        TIntDoubleMap indexAngleMap = new TIntDoubleHashMap();
        
        for (int i = 0; i < ns; ++i) {
            int x1 = surfaceX[i];
            int y1 = surfaceY[i];
            int idx2;
            if (i < (ns - 1)) {
                idx2 = i + 1;
            } else {
                idx2 = 0;
            }
            int x2 = surfaceX[idx2];
            int y2 = surfaceY[idx2];
            
            /* law of cosines:  cosine A = (b^2 + c^2 - a^2)/2bc
              B   a   C
                c   b
                  A
            
            where B and C and the circle surface points
            and A is the boundary point.
            
            a, b, and c are the lengths of the segments
            between the surrounding points (B,C), (C,A),
            and (A,B) respectively.
            
            Note that c and b must be much larger than a.
            */
            
            double cSq = distanceSq(x1, y1, nearestBounds[0].getX(),
                nearestBounds[0].getY());
            
            double bSq = distanceSq(x2, y2, nearestBounds[0].getX(),
                nearestBounds[0].getY());
            
            double aSq = distanceSq(x1, y1, x2, y2);
            
            if (cSq < aSq || bSq < aSq) {
                log.fine(String.format("REMOVED1 (%d,%d) (%d,%d)", 
                    x1, y1, x2, y2));
                continue;
            }
            
            //cosine A = (b^2 + c^2 - a^2)/2bc
            double cosA = (bSq + cSq - aSq)/(2 * Math.sqrt(bSq * cSq));
            
            final double angleA = Math.acos(cosA);
    
            if (angleA < threshold) {
                log.fine(String.format("REMOVED2 angleA=%.4f (%d,%d) (%d,%d)", 
                    (float)angleA, x1, y1, x2, y2));
                continue;
            }
            
            log.fine(String.format("angle=%.4f (%d,%d) (%d,%d) i=[%d,%d]", 
                (float)angleA, x1, y1, x2, y2, i, idx2));
            
            /*
            figure 2 of the paper suggests that with a fast nearest 
            neighbors for the boundary points,
            find the nearest for 
            (x1, y1) and (x2, u2) and if they are 
            roughly equivalent (no tolerance is mentioned or shown),
            then one would take the average of the
            2 points as a medial axis point unless they 
            both point to a reflex point.
            */
   
            List<PairInt> combinedNBP = 
                haveEquidistantNearestPoints(x1, y1, x2, y2, tol);
        
            if (combinedNBP != null) {
                if (indexAngleMap.containsKey(i)) {
                    double prevAngle = indexAngleMap.get(i);
                    if (prevAngle > angleA) {
                        continue;
                    }
                }
                if (indexAngleMap.containsKey(idx2)) {
                    double prevAngle = indexAngleMap.get(idx2);
                    if (prevAngle > angleA) {
                        continue;
                    }
                }
                
                // TODO: consider using a small tolerance for same
                //   angle decision here:
                if ((indexAngleMap.containsKey(i) && 
                    (indexAngleMap.get(i) == angleA)) ||
                    (indexAngleMap.containsKey(idx2) && 
                    (indexAngleMap.get(idx2) == angleA))    
                    ) {
                    // should choose the outer numbers of both
                    int nPrev = surfIdxes1.size() - 1;
                    if (surfIdxes2.get(nPrev) == i) {
                        // replace 2nd number w/ current idx2
                        int prevIdx2 = surfIdxes2.get(nPrev);
                        indexAngleMap.remove(prevIdx2);
                        indexAngleMap.put(idx2, angleA);
                        surfIdxes2.set(nPrev, idx2);
                        surfNBMap.put(nPrev, combinedNBP.toArray(
                            new PairInt[combinedNBP.size()]));
                        continue;
                    } else if (surfIdxes2.get(0) == i) {
                        int prevIdx2 = surfIdxes2.get(0);
                        indexAngleMap.remove(prevIdx2);
                        indexAngleMap.put(idx2, angleA);
                        surfIdxes2.set(0, idx2);
                        surfNBMap.put(0, combinedNBP.toArray(
                            new PairInt[combinedNBP.size()]));
                        continue;
                    }
                }
                surfNBMap.put(surfIdxes1.size(), combinedNBP.toArray(
                    new PairInt[combinedNBP.size()]));
                surfIdxes1.add(i);
                surfIdxes2.add(idx2);
                indexAngleMap.put(i, angleA);
                indexAngleMap.put(idx2, angleA);
            
                log.fine("  <-- prev is a med axis pt");
            }
        }       
        
        if (surfIdxes1.isEmpty()) {
            return -2;
        }
        
        for (int i = 0; i < surfIdxes1.size(); ++i) {
            int idx1 = surfIdxes1.get(i);
            int idx2 = surfIdxes2.get(i);
            
            int avgX = Math.round(0.5f*(surfaceX[idx1] + surfaceX[idx2]));
            int avgY = Math.round(0.5f*(surfaceY[idx1] + surfaceY[idx2]));
            PairInt medialP = new PairInt(avgX, avgY);

            if (!(points.contains(medialP) || processed.contains(medialP))) {
                continue;
            }
            
            PairInt[] nearB = surfNBMap.get(i);
            
            MedialAxisPoint mp = createMedialAxisPoint(medialP,
                nearB, r);
            
            output.add(mp);
        }
        
        /*
        More on the separation angle criteria:
            the distance from the candidate medial
            axis point to the nearest boundary point 
            should be much larger than the separation 
            of the adjacent boundary points...
        
        separation angle criteria:
            For a point on the sphere, the angle is defined
            as the angle between the nearest point on the 
            boundary and the nearest point on the other boundary.
            -- the authors Yang, Brock, and Moll  
               use 2 points on the sphere instead and the third 
               point is the closest point for them on the boundary
               (which should be the same point for them).
               the angle subtended by the third point
               is maximum for the point on the medial axis.
               if that angle is larger than a threshold
               theta t, the midpoint between them becomes
               the aMA point.
        */

        return 1;
    }
     
    /**
     * calculate euclidean distance between (x,y) and point p.
     * @param x
     * @param y
     * @param p
     * @return 
     */
    protected double distance(int x, int y, PairInt p) {
        int diffX = x - p.getX();
        int diffY = y - p.getY();
        double dist = Math.sqrt(diffX * diffX + diffY * diffY);
        return dist;
    }
    
    /**
     * calculate the square of the euclidean distance
     * between (x1,y1) and (x2, y2)
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return 
     */
    protected double distanceSq(int x1, int y1, int x2, int y2) {
        int diffX = x1 - x2;
        int diffY = y1 - y2;
        return (diffX * diffX + diffY * diffY);
    }

    /**
     * compare the nearest boundary points of (x1, y1)
     * and (x2, y2) and return those that were the
     * same distance from the respective points and
     * not the same direction vectors.
     * 
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param tol
     * @return 
     */
    private List<PairInt> haveEquidistantNearestPoints(
        int x1, int y1, int x2, int y2, int tol) {
        
        PairInt[] np1 = findNearestBoundsAsArray(x1, y1);
        PairInt[] np2 = findNearestBoundsAsArray(x2, y2);
       
        List<PairInt> output = new ArrayList<PairInt>();
        
        for (int i = 0; i < np1.length; ++i) {
            PairInt p1 = np1[i];
            double dist1 = distance(x1, y1, p1);
            for (int j = 0; j < np2.length; ++j) {
                PairInt p2 = np2[j];
                
                //TODO: check this exclusion:
                if (p1.equals(p2)) {
                    continue;
                }
                double dist2 = distance(x2, y2, p2);
                log.fine("   dist diff=" + 
                    Math.abs(dist1 - dist2) + " dist1=" + dist1 
                    + " dist2=" + dist2);
                if ((dist1 == dist2) || ((dist1 > tol && dist2 > tol) &&
                    (Math.abs(dist1 - dist2) <= tol))) {
                    int[] dir1 = calculateNeighborDirection(x1, y1, p1);
                    int[] dir2 = calculateNeighborDirection(x2, y2, p2);
                    log.fine("   tolsq=" + tol
                        + " np1=" + p1.toString()
                        + " np2=" + p2.toString()
                        + "\n         dir1=" + Arrays.toString(dir1)
                        + " dir2=" + Arrays.toString(dir2)
                    );
                    if (!Arrays.equals(dir1, dir2)) {
                        output.add(np1[i]);
                        output.add(np2[j]);
                    }
                }
            }
        }
        
        if (output.isEmpty()) {
            return null;
        }
        
        return output;
    }
    
    /**
     * examine the given nearestBounds to find the points
     * that are equidistant within tolerance from (x, y),
     * and have different direction vectors and return
     * that subset of points.
     * 
     * @param x
     * @param y
     * @param tol
     * @param nearestBounds
     * @return 
     */
    private PairInt[] findEquidistantNearestPoints(
        int x, int y, int tol, PairInt[] nearestBounds) {
                
        if (nearestBounds.length < 2) {
            return null;
        }

        TIntList maxSameIdxs = null;
        double dist = Double.MAX_VALUE;
        
        for (int i = 0; i < nearestBounds.length; ++i) {
            
            // looking for 2 or more points w/ same
            //   distances and different directions
            
            PairInt p1 = nearestBounds[i];
            double dist1 = distance(x, y, p1);
            int[] dir1 = calculateNeighborDirection(x, y, p1);
            
            TIntList sameIdxs = new TIntArrayList();
            sameIdxs.add(i);
            
            for (int j = (i + 1); j < nearestBounds.length; ++j) {
                PairInt p2 = nearestBounds[j];
                double dist2 = distance(x, y, p2);
                int[] dir2 = calculateNeighborDirection(x, y, p2);
                
                boolean t1 = Math.abs(dist1 - dist2) <= tol;
                boolean t2 = !Arrays.equals(dir1, dir2);
                if (t1 && t2) {
                    sameIdxs.add(j);
                }
            }
            
            if ((sameIdxs.size() > 1) &&
                ((maxSameIdxs == null) || 
                (sameIdxs.size() > maxSameIdxs.size()))) {
                maxSameIdxs = sameIdxs;
                dist = dist1;
            }
        }
        
        if (maxSameIdxs != null) {
            assert(maxSameIdxs.size() > 1);
            PairInt[] output = new PairInt[maxSameIdxs.size()];
            for (int i = 0; i < maxSameIdxs.size(); ++i) {
                int idx = maxSameIdxs.get(i);
                output[i] = nearestBounds[idx];
            }
            return output;
        }
        
        return null;
    }

    /**
     * calculate whether the point that was the center
     * of the circle surface search is a medial axis
     * point and return it if so, else return null.
     * This method is used when a search returns 
     * 2 medial axis points.
     * @param p
     * @param medAxis1
     * @param medAxis2
     * @param tol
     * @return a populated instance of MedialAxisPoint
     * is returned if p is a medial axis point, else
     * null is returned.
     */
    private MedialAxisPoint centerIsAlsoMedialAxisPoint(
        PairInt p, MedialAxisPoint medAxis1, 
        MedialAxisPoint medAxis2, double tol) {
        
        if (!(points.contains(p) || processed.contains(p))) {
            return null;
        }
        
        // test if p is on line between medAxis1 and medAxis2
        PairInt mp1 = medAxis1.getCenter();
        PairInt mp2 = medAxis2.getCenter();
        
        if (LinesAndAngles.onSegment(mp1.getX(), mp1.getY(),
            mp2.getX(), mp2.getY(), p.getX(), p.getY())) {
            
            // create structure for medial axis point
                        
            PairInt[] nearestBounds = findNearestBoundsAsArray(
                p.getX(), p.getY());
            
            nearestBounds = findEquidistantNearestPoints(
                p.getX(), p.getY(), 0, nearestBounds);
            
            if (nearestBounds == null || (nearestBounds.length < 2)) {
                return null;
            }
            
            // use the average of search radii to approximate the
            // search radius used.
            double avgSrchR = 0.5 * (medAxis1.getSearchRadiusUsed() +
                medAxis2.getSearchRadiusUsed());
        
            MedialAxisPoint mp0 = createMedialAxisPoint(p, nearestBounds,
                avgSrchR);
            
            return mp0;
        }
        
        return null;
    }

    /**
     * given nearestBounds which are equidistant from
     * medialP, construct a MedialAxisPoint.
     * This method assumes that nearestBounds have all
     * been checked for different directions.
     * @param medialP
     * @param nearestBounds
     * @return 
     */
    private MedialAxisPoint createMedialAxisPoint(
        PairInt medialP, PairInt[] nearestBounds,
        double srchRadiusUsed) {
        
        MedialAxisPoint mp = new MedialAxisPoint(medialP,
            nearestBounds.length, srchRadiusUsed);
        
        for (int j = 0; j < nearestBounds.length; ++j) {

            PairInt boundaryP = nearestBounds[j];
            
            double dist = distance(medialP.getX(), 
                medialP.getY(), boundaryP);

            // vector from medial axis point to boundary point
            MVector pv = new MVector(boundaryP, dist);

            //boundary point and medial axis
            mp.add(pv, j);
        }
            
        return mp;
    }

    /**
     * examine each medial axis point and if the nearest
     * neighbors aren't exactly the same difference and
     * with different direction vectors, then use a 
     * small dither around the point to find the real
     * medial axis point.
     * The edited and added medial axis points are added to the
     * instance variable medAxisList.
     * The points that did not need changes, and new
     * points added are returned in a set for later
     * use (points searched around with the dither
     * pattern are excluded to avoid re-searching them).
     * @return 
     */
    private Set<PairInt> refineCentersOfMedAxisList() {

        // store the points that were already true medial
        // axis points and store the points that were
        // added as refinements
        // to possibly use later for seeds to a DFS search for
        // unfound, connected points.
        Set<PairInt> store = new HashSet<PairInt>();
        
        Set<PairInt> present = new HashSet<PairInt>();
        for (int i = 0; i < medAxisList.size(); ++i) {
            MedialAxisPoint mp = medAxisList.get(i);
            present.add(mp.getCenter());
        }        
        
        List<MedialAxisPoint> addTo = new ArrayList<MedialAxisPoint>();
        TIntList rm = new TIntArrayList();
        int[] offsets = null;
        int prevTol = Integer.MIN_VALUE;
        
        for (int i = 0; i < medAxisList.size(); ++i) {
            MedialAxisPoint medAxisPt = medAxisList.get(i);
            PairInt medAxisCenter = medAxisPt.getCenter();
        
            double origSrchR = medAxisPt.getSearchRadiusUsed();
            
            int tol = (int) Math.ceil(origSrchR * sinePiDivN);
            if (tol < 1) {
                tol = 1;
            }
            log.fine("medAxisPt.center=" + medAxisCenter);
            
            // if medAxisPt nearest bounds are already equal and diff dir,
            // skip refinement.
            // because some of the medAxisPts were built on averaged
            //   bunds of adjacent pts on spehere, the nearest bounds
            //   and directions should be re-checked here.
                        
            PairInt[] nearB = findNearestBoundsAsArray(
                medAxisCenter.getX(), medAxisCenter.getY());
            nearB = findEquidistantNearestPoints(
                medAxisCenter.getX(), medAxisCenter.getY(), 0, nearB);
            if (nearB != null && (nearB.length > 1)) {
                //update the point
                MedialAxisPoint mp2 = createMedialAxisPoint(
                    medAxisCenter, nearB, origSrchR);
                medAxisList.set(i, mp2);
                store.add(medAxisCenter);
                continue;
            }
                    
            if (tol > 3) {
                String msg = "Algorithm needs logic for "
                    + "dither with large radius=" + tol;
                // temporarily capturing case as an exception.
                // need to consider fast options for large srch radius here.
                //throw new IllegalStateException(msg);
                log.warning(msg);
                continue;
            }
            
            if (tol != prevTol) {
                offsets = Misc.createOrderedNeighborOffsets(tol);
                prevTol = tol;
            }
            
            log.fine("  dither to improve " + medAxisCenter + ""
                + " dither=" + tol);
            
            List<PairInt> better = new ArrayList<PairInt>();
            TIntObjectMap<PairInt[]> betterNBs = new TIntObjectHashMap<PairInt[]>();
            for (int k = 0; k < offsets.length; k += 2) {
                int x2 = medAxisCenter.getX() + offsets[k];
                if (x2 < minMaxXY[0] || (x2 > minMaxXY[1])) {
                    continue;
                }
                int y2 = medAxisCenter.getY() + offsets[k + 1];
                if (y2 < minMaxXY[2] || (y2 > minMaxXY[3])) {
                    continue;
                }
                PairInt p2 = new PairInt(x2, y2);
                
                if (present.contains(p2)) {
                    continue;
                }
                if (!(points.contains(p2) || processed.contains(p2))) {
                    continue;
                }
 
                nearB = findNearestBoundsAsArray(x2, y2);
                nearB = findEquidistantNearestPoints(
                    x2, y2, 0, nearB);
                
                if (nearB != null && (nearB.length > 1)) {
                    betterNBs.put(better.size(), nearB);
                    better.add(p2);
                    store.add(p2);
                }
            }
            if (better.isEmpty()) {
                log.fine("  WARNING: did not find better for " + medAxisCenter);
                rm.add(i);
            } else{
                // replace mp w/ first and add remaining
                MedialAxisPoint mp2 = createMedialAxisPoint(
                    better.get(0), betterNBs.get(0), origSrchR);
                medAxisList.set(i, mp2);
                present.remove(medAxisCenter);
                present.add(better.get(0));
                store.add(better.get(0));
                for (int jj = 1; jj < better.size(); ++jj) {
                    MedialAxisPoint mp3 = createMedialAxisPoint(
                        better.get(jj), betterNBs.get(jj), origSrchR);
                    addTo.add(mp3);
                    present.add(better.get(jj));
                    store.add(better.get(jj));
                }
            }
        }
        rm.sort();
        for (int i = (rm.size() - 1); i > -1; --i) {
            int idx = rm.get(i);
            medAxisList.remove(idx);
        }
        medAxisList.addAll(addTo);

        return store;
    }

    /**
     * search around the points in srch using a dither
     * of 1 to find new medial axis points.
     * The found medial axis points are added to the
     * instance variable medAxisList.
     * @param srch 
     */
    private void fillInGaps(Set<PairInt> srch) {
        
        Stack<PairInt> stack = new Stack<PairInt>();
        stack.addAll(srch);

        // search for nearest neighbors within dist tol
        /*
        NOTE: 
           alternative dither values could be:
           (1) derive from search radius of each medial axis point.
               (would need to change argument to pass those in instead)
           (2) derive from general characteristics of data max:
               int dither = Math.min(minMaxXY[1], minMaxXY[3]);
               dither = (int)Math.ceil(2 * tol * sinePiDivN); 
        
        int[] offsets = Misc.createOrderedNeighborOffsets(tol);
        */
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        Set<PairInt> visited = new HashSet<PairInt>();

        Set<PairInt> srched = new HashSet<PairInt>(srch);

        while (!stack.isEmpty()) {
            PairInt p = stack.pop();
            if (visited.contains(p)) {
                continue;
            }
            visited.add(p);
            int x = p.getX();
            int y = p.getY();

            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < minMaxXY[0] || y2 < minMaxXY[2] || (x2 > minMaxXY[1])
                    || (y2 > minMaxXY[3])) {
                    continue;
                }
                PairInt p2 = new PairInt(x2, y2);
                
                if (!(points.contains(p2) || processed.contains(p2))
                    || srched.contains(p2)) {
                    continue;
                }
                srched.add(p2);

                PairInt[] nearestBounds = findNearestBoundsAsArray(
                    x2, y2);
                nearestBounds = findEquidistantNearestPoints(
                    x2, y2, 0, nearestBounds);
                if (nearestBounds == null || (nearestBounds.length < 2)) {
                    continue;
                }
                double d = distance(x2, y2, nearestBounds[0]);
                MedialAxisPoint mp2 = createMedialAxisPoint(p2,
                    nearestBounds, d);
                medAxisList.add(mp2);
                stack.add(p2);
            }
        }
    }

    private void searchCenterPoint(MedialAxisResults results, 
        Set<PairInt> addedM) {
        
        if (results.medialAxes.size() == 2) {
            // if there are 2 medial axes points and if center is not inline
            // with the 2 med axis points, this is a "critical point" region too
            if (points.contains(results.center) || 
                processed.contains(results.center)) {
                int tol2 = (int)Math.ceil(results.medialAxes.get(0).getSearchRadiusUsed()
                        * sinePiDivN);
                if (tol2 < 1) {
                    tol2 = 1;
                }
                MedialAxisPoint mp0 = centerIsAlsoMedialAxisPoint(results.center,
                    results.medialAxes.get(0),
                    results.medialAxes.get(1),
                    tol2);
                if (mp0 != null) {
                    medAxisList.add(mp0);
                    addedM.add(mp0.getCenter());
                }
            }
        }        
    }

    private MedialAxisResults findAndProcessFirstPoint(
        Heap q, Set<PairInt> addedM) {
        
        MedialAxisResults results = findInitialPoint();
        
        if (results == null) {
            //TODO: refactor so that can remove the searched
            // points even when null is returned.
            return null;
        }
        
        // remove the searched circle from this.points.
        // and add the extracted to "processed".
        // NOTE: have used search radius - 1 to not remove
        // any found medial axis points.
        Set<PairInt> removed = subtractFromPoints(results.center, 
            results.centerSrchR - 1);
         
        processed.addAll(removed);
        
        searchCenterPoint(results, addedM);
        
        assert(assertUniqueMedialAxesPoints());
        
        for (MedialAxisPoint m : results.medialAxes) {
            PairInt pp = m.getCenter();
            addToHeap(q, m);
            addedM.add(pp);
            medAxisList.add(m);
        }

        assert(assertUniqueMedialAxesPoints());  

        return results;
    }
    
    private MedialAxisResults findAndProcessARemainingPoint(
        Heap q, Set<PairInt> addedM) {
        
        PairInt p = points.iterator().next();
        
        // first check the point p itself
        PairInt[] nearB = findNearestBoundsAsArray(
            p.getX(), p.getY());
        nearB = findEquidistantNearestPoints(
            p.getX(), p.getY(), 0, nearB);
        if (nearB != null && (nearB.length > 1)) {
            double d = distance(p.getX(), p.getY(), nearB[0]);
            MedialAxisPoint mp2 = createMedialAxisPoint(
                p, nearB, d);
            medAxisList.add(mp2);
            addToHeap(q, mp2);
            addedM.add(p);
            processed.addAll(subtractFromPoints(p, d - 1));         
        }
        
        // search within a sampled radius for medial axis points
        MedialAxisResults results = findMedialAxesNearPoint(p);
                
        if (results == null) {
            //TODO: refactor so that can remove the searched
            // points even when null is returned.
            return null;
        }
        
        // remove the searched circle from this.points.
        // and add the extracted to "processed".
        // NOTE: have used search radius - 1 to not remove
        // any found medial axis points.
        Set<PairInt> removed = subtractFromPoints(results.center, 
            results.centerSrchR - 1);
         
        processed.addAll(removed);
        
        if (!results.center.equals(p)) {
            searchCenterPoint(results, addedM);
        }
        assert(assertUniqueMedialAxesPoints());
        
        for (MedialAxisPoint m : results.medialAxes) {
            PairInt pp = m.getCenter();
            addToHeap(q, m);
            addedM.add(pp);
            medAxisList.add(m);
        }

        assert(assertUniqueMedialAxesPoints());  

        return results;
    }

    protected static class MedialAxisPoint {
        private final MVector[] boundaryPoints;
        private final PairInt p;
        private MedialAxisPoint parent = null;
        private final double originalSearchRadius;
        public MedialAxisPoint(PairInt center, int nVectors,
            double srchRadiusUsed) {
            boundaryPoints = new MVector[nVectors];
            p = center;
            originalSearchRadius = srchRadiusUsed;
        }
        
        /**
         * add the boundary point and medial axis internal
         * array at index index.
         * @param mv
         * @param index 
         */
        public void add(MVector mv, int index) {
            boundaryPoints[index] = mv;
        }
        
        /**
         * get vectors of boundaryPoint and point which
         * is medial axis.  Note that pv.delta should
         * be the distance from medialAxisPoint to its 
         * nearest boundary point and it may be null
         * so that it can be lazily populated to avoid
         * the log_2(N) NN if not needed.
         * @return 
         */
        public MVector[] getBoundaryPoints() {
            return boundaryPoints;
        }
        public MedialAxisPoint getParent() {
            return parent;
        }
        public void setParent(MedialAxisPoint medAxisP) {
            parent = medAxisP;
        }
        public PairInt getCenter() {
            return p;
        }
        public double getSearchRadiusUsed() {
            return originalSearchRadius;
        }
    }
    
    /**
     * class to hold one of the nearest boundary
     * points to a medial axis point and the
     * distance between them.
     */
    protected static class MVector {
        private final PairInt boundaryP;
        private final double dist;
        // direction is from medAxisP to boundaryP
        private int[] dir = null;
        
        public MVector(PairInt b, double d) {
            this.boundaryP = b;
            this.dist = d;
        }
        public PairInt getBoundaryPoint() {
            return boundaryP;
        }
        public double getDistance() {
            return dist;
        }
        public void setDirection(int[] direction) {
            dir = direction;
        }
        public int[] getDirection() {
            return dir;
        }
    }
  
    private int[] calculateNeighborDirection(int x1, int y1, PairInt p2) {
        /*
            N_x = (boundaryP - medialP vectors)
                  /|(boundaryP - medialP)|
         */
        int diffX = p2.getX() - x1;
        int diffY = p2.getY() - y1;
        int[] result = new int[2];
        result[0] = (diffX < 0) ? -1 : (diffX > 0) ? 1 : 0;
        result[1] = (diffY < 0) ? -1 : (diffY > 0) ? 1 : 0;
        return result;
    }
    
    protected static class MedialAxisResults {
        List<MedialAxisPoint> medialAxes;
        PairInt center;
        double centerSrchR;
    }
    
    /**
     * use the instance variable nb to find nearest neighbor(s)
     * of (x, y) and return them as an array.
     * The returned points have the same distance from (x. y).
     * @param x
     * @param y
     * @return 
     */
    protected PairInt[] findNearestBoundsAsArray(int x, int y) {
    
        Set<PairInt> nearestB = np.findClosest(x, y);
        int count = 0;
        PairInt[] nearestBounds = new PairInt[nearestB.size()];
        for (PairInt np : nearestB) {
            nearestBounds[count] = np;
            count++;
        }
        return nearestBounds;
    }
}
