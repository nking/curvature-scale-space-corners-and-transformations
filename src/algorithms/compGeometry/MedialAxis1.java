package algorithms.compGeometry;

import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * 
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
    private final TIntObjectMap<PairInt> pointIndexMap;
    private final Set<PairInt> boundary;
    private final NearestNeighbor2D np;
    
    //xMin, xMax, yMin, yMax
    private final int[] minMaxXY;
    
    private static final double sepAng = Math.PI/3;
    
    private final int nSampl = 6;//  2*pi/sepAng
    
    public MedialAxis1(final Set<PairInt> shapePoints,
        final List<PairInt> boundaryPoints) {
        
        this.points = new HashSet<PairInt>(shapePoints);
        
        this.boundary = new HashSet<PairInt>(boundaryPoints.size());
    
        points.removeAll(boundary);
        
        pointIndexMap = new TIntObjectHashMap<PairInt>();
        
        int count = 0;
        for (int i = 0; i < boundaryPoints.size(); ++i) {
            PairInt p = boundaryPoints.get(i);
            if (points.contains(p)) {
                pointIndexMap.put(count, p);
                count++;
            }
        }
        
        minMaxXY = MiscMath.findMinMaxXY(boundary);
        
        this.np = new NearestNeighbor2D(points, 0, 0);
        
    }
    
    protected class PointAndRadius {
        PairInt p;
        // distance, usually to nearest boundary point
        double delta;
    }
    
    // medial axis points will have at least 2 of these:
    protected class PVector {
        // closest feature on the boundary
        PairInt boundaryP;
        int bIndex;
        PointAndRadius pd;
        // vector direction is from pd.p to boundaryP
    }
    /*
    -- second primitive operation: 
        given a solid D, a set of points P in the 
        interior of D and their direction vectors, 
        for each point p ∈ P, this method identifies 
        that point on the medial axis of D which 
        is closest to p.
    */
    
    protected PVector findInitialPoint() {
    
        /*identifies an initial point m and associated 
        distance delta(m), such that the resulting 
        sphere of radius delta(m) around m intersects 
        the medial axis.*/
        
        PairInt p = points.iterator().next();
        
        Set<PairInt> closestB = 
            np.findClosest(p.getX(), p.getY());
        
        assert(closestB.size() > 0);
        
        PairInt b1 = closestB.iterator().next();
        
        if (true) {
            throw new UnsupportedOperationException(
            "not yet implemented");
        }
        
        boolean validM = false;
        
        // if 3 or more points, validate that medial
        // axis point is found and sphere does not
        // extend beyond bounds
        
        while ((closestB.size() < 2) || !validM) {
            // make a new sphere center from the point
            // that is p - vector to b1.
        }
        /*
        a new search is normal to the median axis line
        recently formed.
        
        if the new medial axis segment is significantly
        different from normal to the previous,
        then start 2 searches, 1 perpendicular to the previous
        medial axis center and 1 perp to the current medial
        axis.  the points should be away from the
        smallest angle between the old and new segments.
        */
        
        //double dist = distance(b1, p);
       
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
            of the adjacent boundary points..
        
        
        */
         
        
        //     
        // note to self:
        //    should be able to determine if a medial
        //    axis point is within radius of center
        //    when the nearest points size is 2 or larger
        //    and when those 2 or more are on different segments
        //    of the boundary.  
        //    the sphere formed by the radius as the
        //    distance to the nearest points must be within
        //    the boundary and points to be a valid sphere.
        //    the midpoint on the sphere between the two or more
        //    should be within the boundary and points.
        //    would need to account for noise by using a
        //    stability criterion such as above.
        //    but there may be geometries such as
        //    a narrow reflex point (concave section) 
        //    and the sphere
        //    extending outside of the shape, but mid point
        //    of nearest bounds being inside on other side of
        //    reflex point
        //
       
        /*
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
        
        
        /*Starting from a random point p inside D, 
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
}
