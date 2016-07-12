package algorithms.compGeometry;

import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * An implementation of 
 * "Efficient and Robust Computation of an Approximated Medial Axis"
 * by Yang, Brock, and Moll
 * 
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
    
    private final Set<PairInt> points;
    private final TIntObjectMap<PairInt> boundaryIndexMap;
    private final Set<PairInt> boundary;
    private final NearestNeighbor2D np;
    
    //xMin, xMax, yMin, yMax
    private final int[] minMaxXY;
    
    /**
     * separation angle used for sampling points in
     * a circle around a point.
     */
    private static final double sepAng = Math.PI/3;
    
    private final int nSampl = 6;//  2*pi/sepAng
    
    // 10 degrees threshold for separation angle criterion
    private final double threshold = Math.PI/18;
    
    private double twoPI = Math.PI;
    
    /**
     * constructor containing all points in the area
     * and the bounding points.
     * (NOTE: in future may make a version that doesn't
     * need all interior points and uses an adaptive
     * sampling, but use case for now fits these arguments
     * better.)
     * 
     * @param shapePoints
     * @param boundaryPoints 
     */
    public MedialAxis1(final Set<PairInt> shapePoints,
        final List<PairInt> boundaryPoints) {
        
        this.points = new HashSet<PairInt>(shapePoints);
        
        this.boundary = new HashSet<PairInt>(boundaryPoints.size());
    
        points.removeAll(boundary);
        
        boundaryIndexMap = new TIntObjectHashMap<PairInt>();
        
        for (int i = 0; i < boundaryPoints.size(); ++i) {
            PairInt p = boundaryPoints.get(i);
            boundaryIndexMap.put(i, p);
        }
        
        minMaxXY = MiscMath.findMinMaxXY(boundary);
        
        this.np = new NearestNeighbor2D(points, 
            minMaxXY[1], minMaxXY[3]);
        
    }
    
    protected PVector findInitialPoint() {
    
        /*
        identifies an initial point m and associated 
        distance delta(m), such that the resulting 
        sphere of radius delta(m) around m intersects 
        the medial axis.
        */
        
        PairInt p = points.iterator().next();
        
        Set<PairInt> closestB = 
            np.findClosest(p.getX(), p.getY());
        
        assert(closestB.size() > 0);
        
        List<MedialAxisPoint> medialAxes 
            = new ArrayList<MedialAxisPoint>();
        intersectsMedialAxis(closestB, p, medialAxes);
        
        while (medialAxes.isEmpty()) {
            // make a new sphere center from the point
            // that is p - vector to b1.
        }
        
        /*
        A new search is normal to the median axis line
        recently formed.
        
        If the new medial axis segment is significantly
        different from normal to the previous,
        then start 2 searches, 1 perpendicular to the previous
        medial axis center and 1 perp to the current medial
        axis.  the points should be away from the
        smallest angle between the old and new segments.
        */
        
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
         -- I'll use my static nearest neighbors (uses an xfast trie
            and voronoi diagram).
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

             
         note to self:
            should be able to determine if a medial
            axis point is within radius of center
            when the nearest points size is 2 or larger
            and when those 2 or more are on different segments
            of the boundary.  
            the sphere formed by the radius as the
            distance to the nearest points must be within
            the boundary and points to be a valid sphere.
            the midpoint on the sphere between the two or more
            should be within the boundary and points.
            would need to account for noise by using a
            stability criterion such as above.
            but there may be geometries such as
            a narrow reflex point (concave section) 
            and the sphere
            extending outside of the shape, but mid point
            of nearest bounds being inside on other side of
            reflex point
        
        For each border edge, the algorithm uses 
        angle criteria to select a point to make 
        a triangle with the edge. 
        The initial triangulation is refined using 
        a method based on a novel projection operator 
        [15] that is able to approximate points to 
        the point cloud in a robust way. 
        Our algorithm does not need 
        Delaunay triangulations and it can handle 
        surfaces with borders and noisy point clouds.
        */
        
        
        //if cannot find a medial axis point,
        //  take the surface point of sphere around
        //  p (nSampl points) that has largest
        //  dist from b1 and make that the next
        //  sphere center.
        //  note that if have not generated the sampling
        //  points, already know that the furthest 
        //  point will be a vector in the opposite direction
        //  of the single closest boundary point.
        
        
        /* Starting from a random point p inside D, 
        we generate the maximal sphere with the 
        center at p. If we cannot find a medial 
        axis point of the surface of the sphere 
        (how medial axis points are identified 
        is described in Section 3.3), we take 
        the point p′ with the biggest radius as 
        the center of the next sphere. 
        This process is repeated until the 
        sphere of radius δ(p′) around p′ intersects 
        the medial axis. Since this procedure 
        converges towards a sphere of locally maximum 
        radius, its center converges towards a point 
        p′ on the medial axis and the surrounding 
        sphere thus must intersect the medial axis.
        */
        
        return null;
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
    
    make a nearest neighbor instance nb1 out of perimeter points
    
    
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
    
    /**
     * finds the medial axis points if p intersects medial axis.
     * @param nearestB
     * @param p
     * @param output 
     * @return code for result: 
     * -3 means there are too few nearest boundary points;
     * -2 means no medial axis angle larger than threshold was found;
     * -1 means part of the circle extends outside of bounds; 
     * 1 means successfully found the medial axis points
     * and placed them in output.
     */
    protected int intersectsMedialAxis(Set<PairInt> nearestB, 
        PairInt p, List<MedialAxisPoint> output) {
        
        if (nearestB.size() < 2) {
            return -3;
        }
        
        PairInt bPoint = nearestB.iterator().next();
        
        double r = distance(bPoint.getX(), bPoint.getY(), p);
        
        // create uniformly sampled circles
        // and check that all points are within bounds
        // (NOTE: the coarse sampling may lead to a circle
        //   possibly extending beyond the shape boundary
        //   without being detected, so this shortcut
        //   is a work in progress)
        int[] surfaceX = new int[nSampl];
        int[] surfaceY = new int[nSampl];
        for (int i = 0; i < surfaceX.length; ++i) {
            surfaceX[i] = (int) Math.round(p.getX()
                + (r * Math.cos(2. * twoPI * i / nSampl)));
            surfaceY[i] = (int) Math.round(p.getY()
                + (r * Math.sin(2. * twoPI * i / nSampl)));
            PairInt sp = new PairInt(surfaceX[i], surfaceY[i]);
            if (!(boundary.contains(sp) || points.contains(sp))) {
                // a point is outside of bounds
                return -1;
            }
        }
        
        int count = 0;
        PairInt[] nearestBounds = new PairInt[nearestB.size()];
        for (PairInt np : nearestB) {
            nearestBounds[count] = np;
            count++;
        }
        
        /*
        the boundary points are all equidistant, so can
        choose one of them as the corner of a triangle.
        
        for each pair of circle surface points, measuring
        the angle subtended by them and the boundary point
        to find the maximum angle above threshold.
        */
        
        // there may be more than one medial axis point,
        // for example, the circle enscribed in an open
        // triangle has one near the vertex side of circle
        // and on opposite side
        double maxAngle = Double.MIN_VALUE;
        TIntList surfIdxes1 = new TIntArrayList();
        TIntList surfIdxes2 = new TIntArrayList();
         
        for (int i = 0; i < nSampl; ++i) {
            int x1 = surfaceX[i];
            int y1 = surfaceY[i];
            int idx2;
            if (i < nSampl) {
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
                continue;
            }
            
            //cosine A = (b^2 + c^2 - a^2)/2bc
            double cosA = (bSq + cSq - aSq)/(2 * Math.sqrt(bSq * cSq));
            
            double angleA = Math.acos(cosA);
            
            if (angleA < threshold) {
                continue;
            }
            if (angleA > maxAngle) {
                maxAngle = angleA;
                surfIdxes1.clear();
                surfIdxes2.clear();
                surfIdxes1.add(i);
                surfIdxes2.add(idx2);
            } else if (angleA == maxAngle) {
                surfIdxes1.add(i);
                surfIdxes2.add(idx2);
            }
        }
        
        if (maxAngle == Double.MIN_VALUE) {
            return -2;
        }
        
        if (maxAngle < threshold) {
            return -2;
        }
        
        // TODO: if more than one medial axis point,
        // check that they don't both point to the same
        // boundary point if it's a reflex point (concave)
        
        for (int i = 0; i < surfIdxes1.size(); ++i) {
            int idx1 = surfIdxes1.get(i);
            int idx2 = surfIdxes2.get(i);
            
            int avgX = (surfaceX[idx1] + surfaceX[idx2]);
            int avgY = (surfaceY[idx1] + surfaceY[idx2]);
            PairInt medialP = new PairInt(avgX, avgY);
            
            MedialAxisPoint mp = new MedialAxisPoint(nearestBounds.length);
            output.add(mp);
            
            for (int j = 0; j < nearestBounds.length; ++j) {
                PointAndRadius medialAndDist = new PointAndRadius();
                medialAndDist.p = medialP;
                
                PairInt boundaryP = nearestBounds[j];
            
                medialAndDist.delta = distance(avgX, avgY, boundaryP);
            
                // vector from medial axis point to boundary point
                PVector pv = new PVector(boundaryP, medialAndDist);
                
                mp.add(pv, j);
            }
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
     
    protected double distance(int x, int y, PairInt p) {
        int diffX = x - p.getX();
        int diffY = y - p.getY();
        double dist = Math.sqrt(diffX * diffX + diffY * diffY);
        return dist;
    }
    
    protected double distanceSq(int x, int y, int x2, int y2) {
        int diffX = x - x2;
        int diffY = y - y2;
        return (diffX * diffX + diffY * diffY);
    }
    
    protected static class PointAndRadius {
        PairInt p;
        // distance, usually to nearest boundary point
        double delta;
    }
    
    protected static class MedialAxisPoint {
        
        private PVector[] pointToBoundary;
        
        private MedialAxisPoint parent = null;
        
        public MedialAxisPoint(int nVectors) {
            pointToBoundary = new PVector[nVectors];
        }
        
        public void add(PVector pv, int index) {
            pointToBoundary[index] = pv;
        }
        
        public PVector[] getVectors() {
            return pointToBoundary;
        }
        
        public MedialAxisPoint getParent() {
            return parent;
        }
        
        /*a node needs its coordinates, vectors from boundary
            and possibly the parent node which is the node
            adjacent to it that this one is derived from)
        */
    }
    
    // medial axis points will have at least 2 of these:
    protected static class PVector {
        // closest feature on the boundary
        private final PairInt boundaryP;
        //private final int bIndex; can fetch this as needed from map
        private final PointAndRadius pd;
        // neighbor direction  = 
        // (boundaryP - medialP vectors)/|(boundaryP - medialP)|
        // using bitIndex 1 for -X and 2 for -Y
        private final int neighborDirection;
        
        // vector direction is from pd.p to boundaryP
        
        public PVector(PairInt b,
            PointAndRadius pointAndRadius) {
            this.boundaryP = b;
            this.pd = pointAndRadius;
            this.neighborDirection = 
                calculateNeighborDirection();
        }
        
        private int calculateNeighborDirection() {
            /*
            N_x = (boundaryP - medialP vectors)
                  /|(boundaryP - medialP)|
            */
            int diffX = boundaryP.getX() - pd.p.getX();
            int diffY = boundaryP.getY() - pd.p.getY();
            int nx = 0;
            if (diffX < 0) {
                nx |= (1 << 1);
            }
            if (diffY < 0) {
                nx |= (1 << 2);
            }
            return nx;
        }
        public boolean neighborDirectionXIsNegative() {
            return ((neighborDirection & (1L << 1)) != 0);
        }
        public boolean neighborDirectionYIsNegative() {
            return ((neighborDirection & (1L << 2)) != 0);
        }
    }
  
}
