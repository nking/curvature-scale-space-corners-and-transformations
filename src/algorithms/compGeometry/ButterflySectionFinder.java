package algorithms.compGeometry;

import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class ButterflySectionFinder {
    
    /**
     * Find sections of the closed curve that are two pixels wide and separate
     * loops in the curve:
     * <pre>
     * for example:         #  #
     *    #  #  #         #     #
     *  #         #  #  #      #
     *   #   #  # #  #  #  #  #
     * </pre>
     * These sections are thinned to width of '1' by the line thinner,
     * so need to be restored afterwards or prevented from being removed.
     * @param closedCurve
     * @return list of points that are part of the 2 pixel width patterns in
     * a curve where the curve closes, but is still connected.
     */
    public List<Routes> findButterflySections(PairIntArray closedCurve) {
        
        Set<PairInt> points = Misc.convert(closedCurve);
        
        List<Routes> output = new ArrayList<Routes>();
        
        List<Routes> sections = findButterflySectionsLarge(closedCurve, 
            points);
        
        if (sections != null && !sections.isEmpty()) {
            output.addAll(sections);
        }

        List<Routes> sectionsSmall = findButterflySectionsSmall(closedCurve, 
            points);
        
        if (sectionsSmall != null && !sectionsSmall.isEmpty()) {
            output.addAll(sectionsSmall);
        }
        
        return output;
    }

     /**
     * Find sections of the closed curve that are two pixels wide and separate
     * loops in the curve:
     * <pre>
     * for example:         #  #
     *    #  #  #         #     #
     *  #         #  #  #      #
     *   #   #  # #  #  #  #  #
     * </pre>
     * These sections are thinned to width of '1' by the line thinner,
     * so need to be restored afterwards or prevented from being removed.
     * @param closedCurve
     * @param points
     * @return list of points that are part of the 2 pixel width patterns in
     * a curve where the curve closes, but is still connected.
     */
    protected List<Routes> findButterflySectionsLarge(PairIntArray 
        closedCurve, Set<PairInt> points) {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        List<LinkedList<Segment>> candidateSections = new ArrayList<LinkedList<Segment>>();

        LinkedList<Segment> currentList = null;

        for (int i = 0; i < closedCurve.getN(); ++i) {

            int x = closedCurve.getX(i);
            int y = closedCurve.getY(i);
            
            Set<PairInt> neighbors = curveHelper.findNeighbors(x, y, points);

            // scanning for segments
            if (neighbors.size() != 5) {
                if (currentList != null) {
                    candidateSections.add(currentList);
                    currentList = null;
                }
                continue;
            }

            Segment segment = checkSegmentPatterns(x, y, neighbors);
            
            if (segment == null) {
                if (currentList != null) {
                    candidateSections.add(currentList);
                    currentList = null;
                }
                continue;
            }

            if (currentList == null) {
                currentList = new LinkedList<Segment>();
            }

            currentList.add(segment);
        }

        if (currentList != null) {
            candidateSections.add(currentList);
        }

        if (candidateSections.isEmpty()) {
            return null;
        }
        
        //TODO: correct for wrap around of section from end of curve to beginning

        List<Routes> output = new ArrayList<Routes>();

        // -- scan for endpoints --
        for (LinkedList<Segment> section : candidateSections) {
            
            Routes routes = checkForAdjacentEndpoints(points, section);

            if (routes == null) {
                continue;
            }

            if (section.size() > 1) {
                
                routes = checkForAdjacentEndpoints(points, section, routes);

                if (routes == null) {
                    continue;
                }
            }
            
            output.add(routes);
        }

        return output;

        /*
        endpoints for vert:
                       #           #
            # .      = - .      -  - .
          - - .  or  - - .  or     # .
            #          #

        endpoints for horiz:
                -        - -       -
              # - #    # - - #   # - #
              . .        . .       . .

        endpoints for diag:

            -  #             -  #               -  #
            -    .        -  -  .            #  -  .
            #  .   .      #  .    .          -  .    .
                 .              .                  .

        The sections of line which are 2 pixels wide and 1 further from the
        endpoint have 3 non-point neighbors each and 5 point set neighbors
        An area limit further constrains the geometry.
        For sections matching the patterns below, could consider storing
        the pattern for each pix as 'v', 'h', or 'd'...not an apparent use for
        that yet though.

        Segment patterns between endpoints:
                       4
           -  -  -  -  3
           .  #  #  .  2
           .  #  #  .  1
           -  -  -  -  0
        0  1  2  3  4

                        4
           -  .  .  -   3
           -  #  #  -   2
           -  #  #  -   1
           -  .  .  -   0
        0  1  2  3  4

           # # - -   3
           - # # - - 2
           - - # #   1
           - - - - - 0
        0  1 2 3 4 5

        data structures:
           linked lists of found pattern points as segments with
           specialization of each segment as VertSegment, HorizSegment,
              UUdiagsegment, ULdiagsegment



        Scan the line,
           if a point fits one of the 4 segment patterns (4th is diag transformed by x=-x),
           add it to a group and add the remaining pts fitting the pattern to a stack
           -- traverse the stack adding contiguous points to the group that fit the
              pattern.
           -- note where the first point in the pattern started, because
              when there are no more contiguous pattern matching points,
              the scan will continue at the next point after that first,
              but will skip those already added to a group.

        When the scan for groups has finished,
             for each group, need to apply the above endpoint patterns to see
             if the candidate segment is surrounded by 2 endpoints.

             test all candidate group points as adjacent to potential endpoints.

             When a match is found, have to exclude all of the matching pattern
             from the oppossing endpoint tests.

             This is the smallest pattern which will match that suggestion:
              - - - -
            # # @ # #
          - - # @ # -
            # - - - #
             The '@'s are the candidate group points.  The #'s are points
             matching endoint patterns.

             The found endpoints for one end, the left for example, would
             be excluded from a search for matching to the other endpoints.

        Note that this pattern and variants of it as very short sections and
        endpoints should be scanned after the above to find the shortest
        butterfly segments.
              - -
            # # @ #
          - - # @ -
            # - - #
        
          #  
            #
              # # #
        # # #
              # # #

        For each segment group which has 2 matching endpoints, those should
        be stored as butterfly sections in a set.  Each one of those
        should be passed back in a list as the return of this method.

        runtime complexity is linear in the number of points in the given
        closed curve.
        */

    }
    
    protected List<Routes> findButterflySectionsSmall(PairIntArray 
        closedCurve, Set<PairInt> points) {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        List<Routes> output = new ArrayList<Routes>();    

        for (int i = 0; i < closedCurve.getN(); ++i) {

            int x = closedCurve.getX(i);
            int y = closedCurve.getY(i);
            
            Set<PairInt> neighbors = curveHelper.findNeighbors(x, y, points);

            // scanning for segments
            if (neighbors.size() != 3) {
                continue;
            }
                            
            Segment segment = checkZigZagSegmentPattern(x, y, points);

            if (segment == null) {
                segment = checkZigZag2SegmentPattern(x, y, points);
                if (segment == null) {
                    continue;
                }
            }

            /*
            Each of the 4 segment points needs at least one neighbor that is 
            not one of the 4 points in the zig zap and all of their neighbors 
            cannot be adjacent to any of the other neighbors.
            */
          
            Set<PairInt> endPoints = checkForZigZagEndPoints(points, segment);
                
            if (endPoints == null || endPoints.isEmpty()) {
                continue;
            }
            
            ZigZagSegmentRoutes routes = parseZigZag(segment, endPoints);
           
            output.add(routes);
        }

        return output;
    }

    private Segment checkSegmentPatterns(final int x, final int y, 
        Set<PairInt> neighbors) {
        
        boolean useVert = true;

        Segment segment = checkVertHorizSegmentPattern(x, y, neighbors, useVert);
        
        if (segment != null) {
            return segment;
        }
        
        useVert = false;
               
        segment = checkVertHorizSegmentPattern(x, y, neighbors, useVert);
        
        if (segment != null) {
            return segment;
        }
       
        segment = checkDiagSegmentPattern(x, y, neighbors);
        
        if (segment != null) {
            return segment;
        }
        
        return null;
    }

    private void swapYDirection(Pattern pattern) {
        // ----- change the sign of y  -----
        for (PairInt p : pattern.zeroes) {
            p.setY(-1 * p.getY());
        }
        for (PairInt p : pattern.ones) {
            p.setY(-1 * p.getY());
        }
    }
    
    private void swapXDirection(Pattern pattern) {
        // ----- change the sign of x  -----
        for (PairInt p : pattern.zeroes) {
            p.setX(-1 * p.getX());
        }
        for (PairInt p : pattern.ones) {
            p.setX(-1 * p.getX());
        }
    }

    private Routes checkForAdjacentEndpoints(Set<PairInt> points,
        LinkedList<Segment> section) {
    
        // the list of segments was built from a closed curve
        
        Segment firstSegment = section.getFirst();
                      
        Routes routes = null;
        
        boolean checkFirstSegment = true;
        
        if (firstSegment instanceof VertSegment) {
            
            boolean useVertical = true;
            
            routes = findEndPointsVertHorizPatterns(points, routes, firstSegment,
                useVertical);
            
        } else if (firstSegment instanceof HorizSegment) {
            
            boolean useVertical = false;
            
            routes = findEndPointsVertHorizPatterns(points, routes, firstSegment,
                useVertical);
            
        } else if (firstSegment instanceof UUDiagSegment) {
            
            routes = findEndPointsDiagPatterns(points, routes, 
                (UUDiagSegment)firstSegment);
        }
        
        if (routes == null) {
            return routes;
        }
        
        if (checkFirstSegment && (section.size() > 1)) {
            for (int i = 1; i < section.size(); ++i) {
                
                Segment segment = section.get(i);
                
                // add node to routes
                addSegmentToRoutes(routes, segment);
            }
        }
        
        return routes;
    }

    private Routes findEndPointsVertHorizPatterns(Set<PairInt> points, 
        Routes routes, Segment segment, boolean useVert) {
        
        int x0 = segment.p0.getX();
        int y0 = segment.p0.getY();
        
        for (int i = 0; i < 6; ++i) {
            Pattern pattern;
            switch(i) {
                case 0:
                    pattern = getEndPointsVertPattern1();
                    break;
                case 1:
                    pattern = getEndPointsVertPattern1Opp();
                    break;
                case 2:
                    pattern = getEndPointsVertPattern2();
                    break;
                case 3:
                    pattern = getEndPointsVertPattern2Opp();
                    break;
                case 4:
                    pattern = getEndPointsVertPattern3();
                    break;
                case 5:
                    pattern = getEndPointsVertPattern3Opp();
                    break;
                default:
                    return null;
            }
            if (!useVert) {
                rotatePattern(pattern, -0.5*Math.PI);
            }
            boolean found = true;
            for (PairInt p : pattern.zeroes) {
                PairInt p2 = new PairInt(x0 + p.getX(), y0 + p.getY());
                if (points.contains(p2)) {
                    found = false;
                    break;
                }
            }
            if (found) {
                Set<PairInt> endPoints = new HashSet<PairInt>();
                for (PairInt p : pattern.ones) {
                    PairInt p2 = new PairInt(x0 + p.getX(), y0 + p.getY());
                    if (segment.contains(p2)) {
                        continue;
                    }
                    if (!points.contains(p2)) {
                        found = false;
                        break;
                    }
                    endPoints.add(p2);
                }
                if (found) {
                    if (routes == null) {
                        routes = createRoute(segment, x0, y0);
                    }
                    switch(i) {
                        case 0:
                        case 2:
                        case 4:
                            addEndpointsForHorizVertPatternForward(x0, y0, 
                                pattern, routes, endPoints, useVert);
                            break;
                        case 1:
                        case 3:
                        case 5:
                            addEndpointsForHorizVertPatternOppos(x0, y0, 
                                pattern, routes, endPoints, useVert);
                            break;
                        default:
                            return null;
                    }
                }
            }
        }        
        return routes;
    }
    
    private Pattern getEndPointsVertPattern1() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
        searching for .'s and #'s
                        -2
           -  2  1  -   -1
           -  3  0  -    0
           -  .  .  -    1
              #  -  #    2
                 -       3
          -2 -1  0  1
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
        
        // '.' intermediate points 
        pattern.ones.add(new PairInt(-1, 1));
        pattern.ones.add(new PairInt(0, 1));
        
        // '#' end points
        PairInt t1 = new PairInt(-1, 2);
        pattern.ones.add(t1);
        PairInt t0 = new PairInt(1, 2);
        pattern.ones.add(t0);
        pattern.ep0 = t0;
        pattern.ones.add(t1);
        pattern.ep1 = t1;
        
        pattern.zeroes.add(new PairInt(-2, 1)); pattern.zeroes.add(new PairInt(1, 1));
        pattern.zeroes.add(new PairInt(0, 2)); pattern.zeroes.add(new PairInt(0, 3));
        return pattern;
    }
    
    private Pattern getEndPointsVertPattern1Opp() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
                         
                 -      -3
              #  -  #   -2
           -  .  .  -   -1
           -  3  0  -    0
           -  2  1  -    1           
        
          -2 -1  0  1
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
        
        pattern.ones.add(new PairInt(-1, -1)); pattern.ones.add(new PairInt(-1, -2)); 
        pattern.ones.add(new PairInt(0, -1));
        pattern.ones.add(new PairInt(1, -2));
        
        pattern.zeroes.add(new PairInt(-2, -1)); 
        pattern.zeroes.add(new PairInt(0, -2)); pattern.zeroes.add(new PairInt(0, -3));
        pattern.zeroes.add(new PairInt(1, -1));
        return pattern;
    }
    
    private Pattern getEndPointsVertPattern2() {
        
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
        searching for .'s and #'s
                        -2
           -  2  1  -   -1
           -  3  0  -    0
           -  .  .  -    1
           #  -  -  #    2
              -  -       3
          -2 -1  0  1
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
        
        // '.' intermediate points
        pattern.ones.add(new PairInt(-1, 1)); 
        pattern.ones.add(new PairInt(0, 1));
        
        // '#' end points
        PairInt t1 = new PairInt(-2, 2);
        pattern.ones.add(t1); 
        PairInt t0 = new PairInt(1, 2);
        pattern.ones.add(t0);
        pattern.ep1 = t1;
        pattern.ep0 = t0;
        
        pattern.zeroes.add(new PairInt(-2, 1)); pattern.zeroes.add(new PairInt(1, 1));
        pattern.zeroes.add(new PairInt(-1, 2)); pattern.zeroes.add(new PairInt(-1, 3));
        pattern.zeroes.add(new PairInt(0, 2)); pattern.zeroes.add(new PairInt(0, 3));
      
        return pattern;
    }
    
    private Pattern getEndPointsVertPattern2Opp() {
        
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
              -  -      -3
           #  -  -  #   -2
           -  .  .  -   -1
           -  3  0  -    0
           -  2  1  -    1
        
          -2 -1  0  1
        
        */
     
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
        
        pattern.ones.add(new PairInt(-2, -2)); 
        pattern.ones.add(new PairInt(-1, -1));
        pattern.ones.add(new PairInt(0, -1)); 
        pattern.ones.add(new PairInt(1, -2));
        
        pattern.zeroes.add(new PairInt(-2, -1)); 
        pattern.zeroes.add(new PairInt(-1, -2)); pattern.zeroes.add(new PairInt(-1, -3));
        pattern.zeroes.add(new PairInt(0, -3)); pattern.zeroes.add(new PairInt(0, -2));
        
        return pattern;
    }

    private Pattern getEndPointsVertPattern3() {
        
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
        searching for .'s and #'s
                        -2
           -  2  1  -   -1
           -  3  0  -    0
           -  .  .  -    1
           #  -  #       2
              -          3
          -2 -1  0  1
        */
        
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
        
        pattern.zeroes.add(new PairInt(-2, 1)); 
        pattern.zeroes.add(new PairInt(-1, 2)); 
        pattern.zeroes.add(new PairInt(-1, 3));
        pattern.zeroes.add(new PairInt(1, 1)); 
        
        // '.' intermediate points
        pattern.ones.add(new PairInt(-1, 1)); 
        pattern.ones.add(new PairInt(0, 1));
        // '#' end points
        PairInt t1 = new PairInt(-2, 2);
        pattern.ones.add(t1);
        PairInt t0 = new PairInt(0, 2);
        pattern.ones.add(t0);
        pattern.ep1 = t1;
        pattern.ep0 = t0;
        
        return pattern;
    }
    
    private Pattern getEndPointsVertPattern3Opp() {
        
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
              -         -3
           #  -  #      -2
           -  .  .  -   -1
           -  3  0  -    0
           -  2  1  -    1
                
          -2 -1  0  1
        */
        
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
        
        pattern.ones.add(new PairInt(-1, -1)); 
        pattern.ones.add(new PairInt(-2, -2));
        pattern.ones.add(new PairInt(0, -1)); pattern.ones.add(new PairInt(0, -2));
        
        pattern.zeroes.add(new PairInt(-2, -1));
        pattern.zeroes.add(new PairInt(-1, -2)); pattern.zeroes.add(new PairInt(-1, -3));
        pattern.zeroes.add(new PairInt(1, -1));
        
        return pattern;
    }
    
    private Routes findEndPointsDiagPatterns(Set<PairInt> points, 
        Routes routes, UUDiagSegment segment) {
        
        int x0 = segment.p0.getX();
        int y0 = segment.p0.getY();
               
        for (int i = 0; i < 6; ++i) {
            Pattern pattern;
            switch(i) {
                case 0:
                    pattern = getEndPointsUUDiagPattern1();
                    break;
                case 1:
                    pattern = getEndPointsUUDiagPattern1();
                    rotatePattern(pattern, -0.5*Math.PI);
                    break;
                case 2:
                    pattern = getEndPointsUUDiagPattern2();
                    break;
                case 3:
                    pattern = getEndPointsUUDiagPattern2();
                    rotatePattern(pattern, -0.5*Math.PI);
                    break;
                case 4:
                    pattern = getEndPointsUUDiagPattern3();
                    break;
                case 5:
                    pattern = getEndPointsUUDiagPattern3();
                    rotatePattern(pattern, -0.5*Math.PI);
                    break;
                default:
                    return null;
            }
            boolean found = true;
            for (PairInt p : pattern.zeroes) {
                PairInt p2 = new PairInt(x0 + p.getX(), y0 + p.getY());
                if (points.contains(p2)) {
                    found = false;
                    break;
                }
            }
            if (found) {
                Set<PairInt> endPoints = new HashSet<PairInt>();
                for (PairInt p : pattern.ones) {
                    PairInt p2 = new PairInt(x0 + p.getX(), y0 + p.getY());
                    //TODO: check this segment test with tests
                    if (segment.contains(p2)) {
                        continue;
                    }
                    if (!points.contains(p2)) {
                        found = false;
                        break;
                    }
                    endPoints.add(p2);
                }
                if (found) {
                    if (routes == null) {
                        routes = createRoute(segment, x0, y0);
                    }
                    switch(i) {
                        case 0:
                        case 2:
                        case 4:
                            addEndpointsForUUDiagPattern(x0, y0, 
                                pattern, routes, endPoints);
                            break;
                        case 1:
                        case 3:
                        case 5:
                            addEndpointsForULDiagPattern(x0, y0, 
                                pattern, routes, endPoints);
                            break;
                        default:
                            return null;
                    }
                }
            }
        }        
        return routes;
    }
    
    private Pattern getEndPointsUUDiagPattern1() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
                -  -          -2
          -  2  1  -  -       -1
          -  -  3 .0  #        0   <--# is route0 end endpoint
             -  .  -           1
                -  #           2   <--# is route1 start endpoint
        
         -3 -2 -1  0  1
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
    
        pattern.zeroes.add(new PairInt(-1, 1));
        pattern.zeroes.add(new PairInt(0, 1)); 
        pattern.zeroes.add(new PairInt(0, -1));
        
        pattern.ones.add(new PairInt(-1, 1)); 
        PairInt t1 = new PairInt(0, 2);
        pattern.ones.add(t1);
        PairInt t0 = new PairInt(1, 0);
        pattern.ones.add(t0);
        pattern.ep1 = t1;
        pattern.ep0 = t0;
       
        return pattern;
    }
    
    private Pattern getEndPointsUUDiagPattern2() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
                -  -          -2
          -  2  1  -  -       -1
          -  -  3 .0  #        0
             -  .  -           1
                #  -           2
        
         -3 -2 -1  0  1
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();

        pattern.zeroes.add(new PairInt(0, 2)); 
        pattern.zeroes.add(new PairInt(0, 1)); 
        pattern.zeroes.add(new PairInt(0, 1));
        
        pattern.ones.add(new PairInt(-1, 1));
        PairInt t1 = new PairInt(-1, 2);
        pattern.ones.add(t1);
        PairInt t0 = new PairInt(1, 0);
        pattern.ones.add(t0);        
        pattern.ep0 = t0;
        pattern.ep1 = t1;
        
        return pattern;
    }
    
    private Pattern getEndPointsUUDiagPattern3() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
                -  -          -2
          -  2  1  -          -1
          -  -  3 .0  -        0
             -  .  -  #        1
                #  -           2
        
         -3 -2 -1  0  1
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
   
        pattern.ones.add(new PairInt(-1, 1));
        PairInt t1 = new PairInt(-1, 2);
        pattern.ones.add(t1); 
        PairInt t0 = new PairInt(1, 1);
        pattern.ones.add(t0);
        pattern.ep0 = t0;
        pattern.ep1 = t1;
        
        pattern.zeroes.add(new PairInt(0, 2)); 
        pattern.zeroes.add(new PairInt(0, 1)); 
        pattern.zeroes.add(new PairInt(0, -1));
        pattern.zeroes.add(new PairInt(1, 0)); 
        
        return pattern;
    }
    
    private Set<PairInt> checkForZigZagEndPoints(Set<PairInt> points, 
        Segment segment) {
        
        /*Each of the 4 needs at least one neighbor that is not one of the 4
        points in the zig zap and all of their neighbors cannot be adjacent
        to any of the other neighbors.*/
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        List<Set<PairInt>> listOfNeighbors = new ArrayList<Set<PairInt>>();
        Set<PairInt> segmentPoints = new HashSet<PairInt>();
        segmentPoints.add(segment.p0);
        segmentPoints.add(segment.p1);
        segmentPoints.add(segment.p2);
        segmentPoints.add(segment.p3);
        for (PairInt p : segmentPoints) {
            Set<PairInt> neighbors = curveHelper.findNeighbors(p.getX(), 
                p.getY(), points);
            neighbors.removeAll(segmentPoints);
            if (neighbors.size() != 1) {
                return null;
            }
            listOfNeighbors.add(neighbors);
        }
        // assert that each in list has no members adjacent to any other members
        // TODO: could use a data structure that uses spatial indexing to make 
        // this faster, but there are not very many points per 4 sets to compare...
        for (int i = 0; i < listOfNeighbors.size(); ++i) {
            Set<PairInt> setI = listOfNeighbors.get(i);
            for (PairInt pI : setI) {
                for (int j = (i + 1); j < listOfNeighbors.size(); ++j) {
                    Set<PairInt> setJ = listOfNeighbors.get(j);
                    for (PairInt pJ : setJ) {
                        if (areAdjacent(pI, pJ)) {
                            return null;
                        }
                    }
                }
            }
        }
        // if arrive here, all neighbor sets have at least one point and none
        // are adjacent to points in a different set.
        Set<PairInt> output = new HashSet<PairInt>();
        for (Set<PairInt> set : listOfNeighbors) {
            output.addAll(set);
        }
        return output;
    }
    
    private boolean areAdjacent(PairInt p0, PairInt p1) {
        int diffX = Math.abs(p0.getX() - p1.getX());
        int diffY = Math.abs(p0.getY() - p1.getY());
        if ((diffX < 2) && (diffY < 2)) {
            return true;
        }
        return false;
    }

    protected void rotatePattern(Pattern pattern, double theta) {

        double sine = Math.sin(theta);
        double cosine = Math.cos(theta);
        
        for (PairInt p : pattern.zeroes) {
            int x = p.getX();
            int y = p.getY();
            int xt = (int)Math.round(((x*cosine) + (y*sine)));
            int yt = (int)Math.round(((-x*sine) + (y*cosine)));
            p.setX(xt);
            p.setY(yt);
        }
        for (PairInt p : pattern.ones) {
            int x = p.getX();
            int y = p.getY();
            int xt = (int)Math.round(((x*cosine) + (y*sine)));
            int yt = (int)Math.round(((-x*sine) + (y*cosine)));
            p.setX(xt);
            p.setY(yt);
        }
    }

    private Routes createRoute(Segment segment, int x0, int y0) {
        
        /*
             VertSegment              swapY VertSegment
              .     .        -2
           -  2     1:E0 -   -1         .     .         -1    E0 marks endpoint 0 and the 
           -  3:E1  0    -    0      -  3     0:E0  -    0       start of route 0.
              .     .         1      -  2:E1  1     -    1       which incr w/ incr y.
                                        .     .          2    route 1 starts at E1 and incr w/ decr y
          -2 -1     0    1          -2 -1     0     1

            HorizSegment              swap HorizSegment
                 -    -       -2          -    -         -2
              .  3:E0 2    .  -1       .  2:E0 3     .   -1    route0 incr w/ incr x
              .  0    1:E1 .   0       .  1    0:E1  .    0    route1 incr w/ decr x
                 -    -        1          -    -          1
             -1  0    1    2          -2 -1    0     1
        
        
                 UUDiagSegment                ULDiagSegment
                    -  -      -2                    2      -2
              -  2  1  -  -   -1                 3  1      -1
              -  -  3  0  -    0                 0          0
                 -  -          1                            1
             -3 -2 -1  0  1                 -1   0  1  2
        
               E0 is '1', and route0         E0 is '1', and route0
                  continues with             continues with
                  point '0'                  point '0'
               E1 is '3', and route1         E1 is '3', and route1
                  continues with '2'            continues with '2'
        */
        
        if (segment instanceof VertSegment) {
            Routes routes = new Routes();
            if (segment.p1.equals(new PairInt(x0, y0 - 1))) {
                assert(segment.p2.equals(new PairInt(x0 - 1, y0 + 1))); 
                assert(segment.p3.equals(new PairInt(x0 - 1, y0))); 
                routes.route0.add(segment.p1);
                routes.route0.add(segment.p0);
                routes.route1.add(segment.p3);
                routes.route1.add(segment.p2);
            } else if (segment.p1.equals(new PairInt(x0, y0 + 1))) {
                assert(segment.p2.equals(new PairInt(x0 - 1, y0 + 1)));
                assert(segment.p3.equals(new PairInt(x0 - 1, y0)));
                routes.route0.add(segment.p0);
                routes.route0.add(segment.p1);
                routes.route1.add(segment.p2);
                routes.route1.add(segment.p3);
            } else {
                throw new IllegalStateException("error in algorithm");
            }
            return routes;
        } else if (segment instanceof HorizSegment) {
            Routes routes = new Routes();
            if (segment.p1.equals(new PairInt(x0 + 1, y0))) {
                assert(segment.p2.equals(new PairInt(x0 + 1, y0 - 1)));
                assert(segment.p3.equals(new PairInt(x0, y0 - 1)));
                routes.route0.add(segment.p3);
                routes.route0.add(segment.p2);
                routes.route1.add(segment.p1);
                routes.route1.add(segment.p0);
            } else if (segment.p1.equals(new PairInt(x0 - 1, y0))) {
                assert(segment.p2.equals(new PairInt(x0 - 1, y0 - 1)));
                assert(segment.p3.equals(new PairInt(x0, y0 - 1)));
                routes.route0.add(segment.p2);
                routes.route0.add(segment.p3);
                routes.route1.add(segment.p0);
                routes.route1.add(segment.p1);
            } else {
                throw new IllegalStateException("error in algorithm!");
            }
            return routes;
            
        } else if ((segment instanceof UUDiagSegment) || (segment instanceof ULDiagSegment)) {
            Routes routes = new Routes();
            routes.route0.add(segment.p1);
            routes.route0.add(segment.p0);
            routes.route1.add(segment.p3);
            routes.route1.add(segment.p2);
            return routes;
        } else {
            throw new IllegalStateException("error in algorithm: " + 
                " segment not handled:" + segment.getClass().getSimpleName());
        }
    }

    private void addEndpointsForHorizVertPatternForward(final int x0, 
        final int y0, Pattern pattern, final Routes routes, 
        Set<PairInt> endPoints, boolean useVert) {
        
        if (endPoints.size() < 4) {
            throw new IllegalStateException("error in algorithm");
        }
        assert(routes.ep0End == null);
        assert(routes.ep1 == null);
            
        if (useVert) {
            /*
               endPointsVertPattern1, for example
                            -2
               -  2  1  -   -1
               -  3  0  -    0
               -  .  .  -    1
                  #  -  #    2
                     -       3
              -2 -1  0  1
            */
            PairInt t = new PairInt(x0, y0 + 1);
            if (!endPoints.contains(t)) {
                throw new IllegalStateException("error in algorithm");
            }
            routes.route0.add(t);
            
            assert(pattern.ep0 != null);
            t = pattern.ep0.copy();
            if (!endPoints.contains(t)) {
                throw new IllegalStateException("error in algorithm");
            }
            routes.route0.add(t);
            routes.ep0End = t;
            
            assert(pattern.ep1 != null);
            t = pattern.ep1.copy();
            if (!endPoints.contains(t)) {
                throw new IllegalStateException("error in algorithm");
            }
            PairInt t2 = new PairInt(x0 - 1, y0 + 1);
            if (!endPoints.contains(t2)) {
                throw new IllegalStateException("error in algorithm");
            }
            LinkedHashSet<PairInt> tmpR1 = new LinkedHashSet<PairInt>();
            tmpR1.add(t);
            tmpR1.add(t2);
            tmpR1.addAll(routes.route1);
            routes.route1 = tmpR1;
            routes.ep0 = t;
        } else {
            /*
                endPointsHorizPattern1
                      -  -  -   -2
                  #   .  3  2   -1
              -   -   .  0  1    0  it's vert rotated by -90, so
                  #   -  -  -    1  route0 is at higher y than route1 here
                                 2
                                 3
             -3  -2  -1  0  1
            */
            assert(pattern.ep1 != null);
            PairInt t = pattern.ep1.copy();
            PairInt t2 = new PairInt(x0 - 1, y0 - 1);
            if (!endPoints.contains(t) || !endPoints.contains(t2)) {
                throw new IllegalStateException("error in algorithm");
            }            
            LinkedHashSet<PairInt> tmpR1 = new LinkedHashSet<PairInt>();
            tmpR1.add(t);
            tmpR1.add(t2);
            tmpR1.addAll(routes.route1);
            routes.route1 = tmpR1;
            routes.ep1 = t;
            
            t = new PairInt(x0 - 1, y0);
            if (!endPoints.contains(t)) {
                throw new IllegalStateException("error in algorithm");
            }
            routes.route0.add(t);
            
            assert(pattern.ep0 != null);
            t = pattern.ep0.copy();
            if (!endPoints.contains(t)) {
                throw new IllegalStateException("error in algorithm");
            }
            routes.route0.add(t);
            routes.ep0End = t;
        }
    }

    private void addEndpointsForHorizVertPatternOppos(int x0, int y0, 
        Pattern pattern, Routes routes, Set<PairInt> endPoints, boolean useVert) {
        
        if (endPoints.size() < 4) {
            throw new IllegalStateException("error in algorithm");
        }
        assert(routes.ep0 == null);
        assert(routes.ep1End == null);
            
        if (useVert) {
            /*
               getEndPointsVertPattern1Opp, for example
                     -      -3
                  #  -  #   -2
               -  .  .  -   -1
               -  3  0  -    0
               -  2  1  -    1

              -2 -1  0  1
            */
            assert(pattern.ep0 != null);
            PairInt t0 = pattern.ep0.copy();
            PairInt t = new PairInt(x0, y0 - 1);
            if (!endPoints.contains(t)) {
                throw new IllegalStateException("error in algorithm");
            }
            LinkedHashSet<PairInt> tmpR0 = new LinkedHashSet<PairInt>();
            tmpR0.add(t0);
            tmpR0.add(t);
            tmpR0.addAll(routes.route0);
            routes.route0 = tmpR0;
            routes.ep0 = t;
            
            t = new PairInt(x0 - 1, y0 - 1);
            if (!endPoints.contains(t)) {
                throw new IllegalStateException("error in algorithm");
            }
            routes.route1.add(t);
            
            t = pattern.ep1.copy();
            if (!endPoints.contains(t)) {
                throw new IllegalStateException("error in algorithm");
            }
            routes.route1.add(t);
            routes.ep1End = t;
            
        } else {
            /*
               getEndPoints Horiz Pattern1Opp
                              -3
                              -2
                  2  3  .  #  -1
                  1  0  .      0
                           #   1

              -2 -1  0  1  2 
            */
            assert(pattern.ep0 != null);
            PairInt t0 = pattern.ep0.copy();
            PairInt t = new PairInt(x0 + 1, y0);
            if (!endPoints.contains(t)) {
                throw new IllegalStateException("error in algorithm");
            }
            LinkedHashSet<PairInt> tmpR0 = new LinkedHashSet<PairInt>();
            tmpR0.add(t0);
            tmpR0.add(t);
            tmpR0.addAll(routes.route0);
            routes.route0 = tmpR0;
            routes.ep0 = t;
          
            t = new PairInt(x0 + 1, y0 - 1);
            if (!endPoints.contains(t)) {
                throw new IllegalStateException("error in algorithm");
            }
            routes.route1.add(t);
            
            t = pattern.ep1.copy();
            if (!endPoints.contains(t)) {
                throw new IllegalStateException("error in algorithm");
            }
            routes.route1.add(t);
            routes.ep1End = t;
        }
    }

    private void addEndpointsForUUDiagPattern(final int x0, final int y0, 
        final Pattern pattern, Routes routes, Set<PairInt> endPoints) {
        
        if (endPoints.size() < 3) {
            throw new IllegalStateException("error in algorithm");
        }
        assert(routes.ep0End == null);
        assert(routes.ep1 == null);
            
        /*
           UUDiagPattern1, for example
                -  -          -2
          -  2  1  -  -       -1
          -  -  3 .0  #        0    <--# is route0 end endpoint
             -  .  -           1
                -  #           2    <--# is route1 start endpoint

         -3 -2 -1  0  1
        */
        assert(pattern.ep0 != null);
        assert(pattern.ep1 != null);
        
        LinkedHashSet<PairInt> tmpR1 = new LinkedHashSet<PairInt>();
        PairInt t = new PairInt(x0 - 1, y0 + 1);
        PairInt t1 = pattern.ep1.copy();
        if (!endPoints.contains(t)) {
            throw new IllegalStateException("error in algorithm");
        }
        tmpR1.add(t1);
        tmpR1.add(t);
        tmpR1.addAll(routes.route1);
        routes.route1 = tmpR1;
        routes.ep1 = t;

        PairInt t0 = pattern.ep0.copy();
        routes.route0.add(t0);
        routes.ep0End = t0;            
    }

    private void addEndpointsForULDiagPattern(int x0, int y0, 
        Pattern pattern, Routes routes, Set<PairInt> endPoints) {
        
        if (endPoints.size() < 3) {
            throw new IllegalStateException("error in algorithm");
        }
        assert(routes.ep0End == null);
        assert(routes.ep1 == null);
            
        /*
           ULDiagPattern1, for example
                      2       -2
                .  3  1       -1
             #     0           0    <--# is route1 start endpoint
                   #           1    <--# is route0 end endpoint
                               2         

         -3 -2 -1  0  1
        */
        assert(pattern.ep0 != null);
        assert(pattern.ep1 != null);
        
        LinkedHashSet<PairInt> tmpR1 = new LinkedHashSet<PairInt>();
        PairInt t = new PairInt(x0 - 1, y0 - 1);
        PairInt t1 = pattern.ep1.copy();
        if (!endPoints.contains(t)) {
            throw new IllegalStateException("error in algorithm");
        }
        tmpR1.add(t1);
        tmpR1.add(t);
        tmpR1.addAll(routes.route1);
        routes.route1 = tmpR1;
        routes.ep1 = t;

        PairInt t0 = pattern.ep0.copy();
        routes.route0.add(t0);
        routes.ep0End = t0;         
    }

    private Routes checkForAdjacentEndpoints(Set<PairInt> points, 
        LinkedList<Segment> section, Routes routes) {
        
        // the list of segments was built from a closed curve
        
        Segment lastSegment = section.getLast();

        if (lastSegment instanceof VertSegment) {
            
            boolean useVertical = true;
            
            routes = findEndPointsVertHorizPatterns(points, routes, lastSegment,
                useVertical);
            
        } else if (lastSegment instanceof HorizSegment) {
            
            boolean useVertical = false;
            
            routes = findEndPointsVertHorizPatterns(points, routes, lastSegment,
                useVertical);
            
        } else if (lastSegment instanceof UUDiagSegment) {
            
            routes = findEndPointsDiagPatterns(points, routes, 
                (UUDiagSegment)lastSegment);
        }
        
        return routes;
    }

    private void addSegmentToRoutes(Routes routes, Segment segment) {
        
        if (segment instanceof VertSegment) {
            addVertSegmentToRoutes(routes, (VertSegment)segment);
        } else if (segment instanceof HorizSegment) {
            addHorizSegmentToRoutes(routes, (HorizSegment)segment);
        } else if (segment instanceof UUDiagSegment) {
            addUUDiagSegmentToRoutes(routes, (UUDiagSegment)segment);
        } else if (segment instanceof ULDiagSegment) {
            addUUDiagSegmentToRoutes(routes, (ULDiagSegment)segment);
        }
        
    }

    private void addVertSegmentToRoutes(Routes routes, VertSegment segment) {
        
        assert(!routes.route0.isEmpty());
        assert(!routes.route1.isEmpty());
        
        // -------- add to routes0 -------------------------
        PairInt[] r0FirstLastNodes = getFirstAndLast(routes.route0.iterator());
        
        double p0MinDistSq = Double.MAX_VALUE;
        int p0MinIdx = -1;
        double p1MinDistSq = Double.MAX_VALUE;
        int p1MinIdx = -1;
        for (int i = 0; i < r0FirstLastNodes.length; ++i) {
            PairInt r0 = r0FirstLastNodes[i];
            int diffX = segment.p0.getX() - r0.getX();
            int diffY = segment.p0.getY() - r0.getY();
            double distSq = (diffX * diffX) + (diffY * diffY);
            if (distSq < p0MinDistSq) {
                p0MinDistSq = distSq;
                p0MinIdx = i;
            }
            diffX = segment.p1.getX() - r0.getX();
            diffY = segment.p1.getY() - r0.getY();
            distSq = (diffX * diffX) + (diffY * diffY);
            if (distSq < p1MinDistSq) {
                p1MinDistSq = distSq;
                p1MinIdx = i;
            }
        }
        
        if (segment.p1.getY() < segment.p0.getY()) {
            
            /*      VertSegment 
            
                       R0
                       \/
                  .     .        -2
               -  2     1:E0 -   -1  
               -  3:E1  0    -    0 
                  .     .         1 
                                  
              -2 -1     0    1   
            
                  /\
                  R1
            */
            if (!(segment.p2.getY() < segment.p3.getY())) {
                throw new IllegalStateException("error in algorithm");
            }
            if (p0MinDistSq < p1MinDistSq) {
                /*
                1
                0
                [*]
                [*]
                E0end
                */
                if (p0MinIdx != 0) {
                    throw new IllegalStateException("error in algorithm");
                }
                if (routes.ep0 != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                LinkedHashSet<PairInt> tmp = new LinkedHashSet<PairInt>();
                tmp.add(segment.p1);
                tmp.add(segment.p0);
                tmp.addAll(routes.route0);
                routes.route0 = tmp;
            } else if (p1MinDistSq < p0MinDistSq) {
                /*
                E0 
                [*]
                [*]
                1
                0
                */
                if (p1MinIdx != 1) {
                    throw new IllegalStateException("error in algorithm");
                }
                if (routes.ep0End != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                routes.route0.add(segment.p1);
                routes.route0.add(segment.p0);
            } else {
                throw new IllegalStateException("error in algorithm");
            }
            
        } else if (segment.p0.getY() < segment.p1.getY()) {
            
            /*   swapped VertSegment
                        R0
                        \/
                   .     .         -1    E0 marks endpoint 0 and the
                -  3     0:E0  -    0       start of route 0.
                -  2:E1  1     -    1       which incr w/ incr y.
                   .     .          2    route 1 starts at E1 and incr w/ decr y
               -2 -1     0     1
            
                  /\
                  R1
            */
            
            if (!(segment.p3.getY() < segment.p2.getY())) {
                throw new IllegalStateException("error in algorithm");
            }
            if (p0MinDistSq < p1MinDistSq) {
                /*
                E0 
                [*]
                [*]
                0
                1
                */
                if (p0MinIdx != 1) {
                    throw new IllegalStateException("error in algorithm");
                }
                if (routes.ep0End != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                routes.route0.add(segment.p0);
                routes.route0.add(segment.p1);
            } else if (p1MinDistSq < p0MinDistSq) {
                /*
                0
                1
                [*]
                [*]
                E0end
                */
                if (p1MinIdx != 0) {
                    throw new IllegalStateException("error in algorithm");
                }
                if (routes.ep0 != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                LinkedHashSet<PairInt> tmp = new LinkedHashSet<PairInt>();
                tmp.add(segment.p0);
                tmp.add(segment.p1);
                tmp.addAll(routes.route0);
                routes.route0 = tmp;
            } else {
                throw new IllegalStateException("error in algorithm");
            }
            
        } else {
            throw new IllegalStateException("error in algorithm");
        }
        
        // ------------ add to routes1 -----------------
        /*   VertSegment              swapY VertSegment
              .     .        -2
           -  2     1:E0 -   -1         .     .         -1    E0 marks endpoint 0 and the
           -  3:E1  0    -    0      -  3     0:E0  -    0       start of route 0.
              .     .         1      -  2:E1  1     -    1       which incr w/ incr y.
                                        .     .          2    route 1 starts at E1 and incr w/ decr y
          -2 -1     0    1          -2 -1     0     1
        */
        PairInt[] r1FirstLastNodes = getFirstAndLast(routes.route1.iterator()); 
        double p2MinDistSq = Double.MAX_VALUE;
        int p2MinIdx = -1;
        double p3MinDistSq = Double.MAX_VALUE;
        int p3MinIdx = -1;
        for (int i = 0; i < r1FirstLastNodes.length; ++i) {
            PairInt r1 = r1FirstLastNodes[i];
            int diffX = segment.p2.getX() - r1.getX();
            int diffY = segment.p2.getY() - r1.getY();
            double distSq = (diffX * diffX) + (diffY * diffY);
            if (distSq < p2MinDistSq) {
                p2MinDistSq = distSq;
                p2MinIdx = i;
            }
            diffX = segment.p3.getX() - r1.getX();
            diffY = segment.p3.getY() - r1.getY();
            distSq = (diffX * diffX) + (diffY * diffY);
            if (distSq < p3MinDistSq) {
                p3MinDistSq = distSq;
                p3MinIdx = i;
            }
        }
        
        if (segment.p2.getY() < segment.p3.getY()) {
            
            /*      VertSegment 
            
                       R0
                       \/
                  .     .        -2
               -  2     1:E0 -   -1  
               -  3:E1  0    -    0 
                  .     .         1 
                                  
              -2 -1     0    1   
            
                  /\
                  R1
            */
            if (!(segment.p1.getY() < segment.p0.getY())) {
                throw new IllegalStateException("error in algorithm");
            }
            if (p3MinDistSq < p2MinDistSq) {
                /*
                2
                3
                [*]
                [*]
                E1
                */
                if (p3MinIdx != 1) {
                    throw new IllegalStateException("error in algorithm");
                }
                if (routes.ep1End != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                routes.route1.add(segment.p3);
                routes.route1.add(segment.p2);
            } else if (p2MinDistSq < p3MinDistSq) {
                /*
                E1end 
                [*]
                [*]
                2
                3
                */
                if (p2MinIdx != 0) {
                    throw new IllegalStateException("error in algorithm");
                }
                if (routes.ep1 != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                LinkedHashSet<PairInt> tmp = new LinkedHashSet<PairInt>();
                tmp.add(segment.p3);
                tmp.add(segment.p2);
                tmp.addAll(routes.route1);
                routes.route1 = tmp;
            } else {
                throw new IllegalStateException("error in algorithm");
            }
            
        } else if (segment.p3.getY() < segment.p2.getY()) {
            
            /*   swapped VertSegment
                        R0
                        \/
                   .     .         -1    E0 marks endpoint 0 and the
                -  3     0:E0  -    0       start of route 0.
                -  2:E1  1     -    1       which incr w/ incr y.
                   .     .          2    route 1 starts at E1 and incr w/ decr y
               -2 -1     0     1
            
                  /\
                  R1
            */
            if (!(segment.p0.getY() < segment.p1.getY())) {
                throw new IllegalStateException("error in algorithm");
            }
            if (p3MinDistSq < p2MinDistSq) {
                /*
                E1end 
                [*]
                [*]
                3
                2
                */
                if (p3MinIdx != 0) {
                    throw new IllegalStateException("error in algorithm");
                }
                if (routes.ep1 != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                LinkedHashSet<PairInt> tmp = new LinkedHashSet<PairInt>();
                tmp.add(segment.p2);
                tmp.add(segment.p3);
                tmp.addAll(routes.route1);
                routes.route1 = tmp;
            } else if (p2MinDistSq < p3MinDistSq) {
                /*
                3
                2
                [*]
                [*]
                E1
                */
                if (p2MinIdx != 1) {
                    throw new IllegalStateException("error in algorithm");
                }
                if (routes.ep1End != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                routes.route1.add(segment.p2);
                routes.route1.add(segment.p3);
            } else {
                throw new IllegalStateException("error in algorithm");
            }
            
        } else {
            throw new IllegalStateException("error in algorithm");
        }
    }
    
    private void addHorizSegmentToRoutes(Routes routes, HorizSegment segment) {
        
        assert(!routes.route0.isEmpty());
        assert(!routes.route1.isEmpty());
        
        /*             HorizSegment              swap HorizSegment
                       -    -       -2          -    -         -2
        R0 >        .  3:E0 2    .  -1       .  2:E0 3     .   -1    route0 incr w/ incr x
        R1 <        .  0    1:E1 .   0       .  1    0:E1  .    0    route1 incr w/ decr x
                       -    -        1          -    -          1
                   -1  0    1    2          -2 -1    0     1
        */
        // -------- add to routes0 -------------------------
        PairInt[] r0FirstLastNodes = getFirstAndLast(routes.route0.iterator());
        
        double p3MinDistSq = Double.MAX_VALUE;
        int p3MinIdx = -1;
        double p2MinDistSq = Double.MAX_VALUE;
        int p2MinIdx = -1;
        for (int i = 0; i < r0FirstLastNodes.length; ++i) {
            PairInt r0 = r0FirstLastNodes[i];
            int diffX = segment.p3.getX() - r0.getX();
            int diffY = segment.p3.getY() - r0.getY();
            double distSq = (diffX * diffX) + (diffY * diffY);
            if (distSq < p3MinDistSq) {
                p3MinDistSq = distSq;
                p3MinIdx = i;
            }
            diffX = segment.p2.getX() - r0.getX();
            diffY = segment.p2.getY() - r0.getY();
            distSq = (diffX * diffX) + (diffY * diffY);
            if (distSq < p2MinDistSq) {
                p2MinDistSq = distSq;
                p2MinIdx = i;
            }
        }
        
        if (segment.p3.getX() < segment.p2.getX()) {
            /*             HorizSegment  
                           -    -       -2 
            R0 >        .  3:E0 2    .  -1 
            R1 <        .  0    1:E1 .   0 
                           -    -        1 
                       -1  0    1    2    
            */
            if (!(segment.p1.getX() > segment.p0.getX())) {
                throw new IllegalStateException("error in algorithm");
            }
            if (p3MinDistSq < p2MinDistSq) {
                if (p3MinIdx != 1) {
                    throw new IllegalStateException("error in algorithm");
                }
                //R0: E0 [*][*]  3 2   append '3' '2' to end of route0
                if (routes.ep0End != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                routes.route0.add(segment.p3);
                routes.route0.add(segment.p2);
            } else if (p3MinDistSq > p2MinDistSq) {
                if (p2MinIdx != 0) {
                    throw new IllegalStateException("error in algorithm");
                }
                if (routes.ep0 != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                //R0: 3  2  [*][*] E0end insert '3' '2' to beginning of route0
                LinkedHashSet<PairInt> tmp = new LinkedHashSet<PairInt>();
                tmp.add(segment.p3);
                tmp.add(segment.p2);
                tmp.addAll(routes.route0);
                routes.route0 = tmp;
            } else {
                throw new IllegalStateException("error in algorithm");
            }
        } else if (segment.p3.getX() > segment.p2.getX()) {
            
            /*           swapped HorizSegment
                          -    -         -2
            R0 >       .  2:E0 3     .   -1    route0 incr w/ incr x
            R1 <       .  1    0:E1  .    0    route1 incr w/ decr x
                          -    -          1
                         -2    -1    0    1
            */
            if (!(segment.p0.getX() > segment.p1.getX())) {
                throw new IllegalStateException("error in algorithm");
            }
            if (p3MinDistSq < p2MinDistSq) {
                if (p3MinIdx != 0) {
                    throw new IllegalStateException("error in algorithm");
                }
                if (routes.ep0 != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                //R0    '2' '3' [*] [*] E0end insert '2' '3' at beginning of route0
                LinkedHashSet<PairInt> tmp = new LinkedHashSet<PairInt>();
                tmp.add(segment.p2);
                tmp.add(segment.p3);
                tmp.addAll(routes.route0);
                routes.route0 = tmp;
            } else if (p3MinDistSq > p2MinDistSq) {
                if (p2MinIdx != 1) {
                    throw new IllegalStateException("error in algorithm");
                }
                if (routes.ep0End != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                //R0 E0 [*] [*] '2' '3'  append '2' '3' to end of route0
                routes.route0.add(segment.p2);
                routes.route0.add(segment.p3);
            } else {
                throw new IllegalStateException("error in algorithm");
            }
        } else {
            throw new IllegalStateException("error in algorithm");
        }
        
        /*             HorizSegment              swap HorizSegment
                       -    -       -2          -    -         -2
        R0 >        .  3:E0 2    .  -1       .  2:E0 3     .   -1    route0 incr w/ incr x
        R1 <        .  0    1:E1 .   0       .  1    0:E1  .    0    route1 incr w/ decr x
                       -    -        1          -    -          1
                   -1  0    1    2          -2 -1    0     1
        */
        // -------- add to routes1 -------------------------
        PairInt[] r1FirstLastNodes = getFirstAndLast(routes.route1.iterator());
        
        double p0MinDistSq = Double.MAX_VALUE;
        int p0MinIdx = -1;
        double p1MinDistSq = Double.MAX_VALUE;
        int p1MinIdx = -1;
        for (int i = 0; i < r1FirstLastNodes.length; ++i) {
            PairInt r1 = r1FirstLastNodes[i];
            int diffX = segment.p0.getX() - r1.getX();
            int diffY = segment.p0.getY() - r1.getY();
            double distSq = (diffX * diffX) + (diffY * diffY);
            if (distSq < p0MinDistSq) {
                p0MinDistSq = distSq;
                p0MinIdx = i;
            }
            diffX = segment.p1.getX() - r1.getX();
            diffY = segment.p1.getY() - r1.getY();
            distSq = (diffX * diffX) + (diffY * diffY);
            if (distSq < p1MinDistSq) {
                p1MinDistSq = distSq;
                p1MinIdx = i;
            }
        }
        
        if (segment.p0.getX() < segment.p1.getX()) {
            /*             HorizSegment  
                           -    -       -2 
            R0 >        .  3:E0 2    .  -1 
            R1 <        .  0    1:E1 .   0 
                           -    -        1 
                       -1  0    1    2    
            */
            if (!(segment.p3.getX() < segment.p2.getX())) {
                throw new IllegalStateException("error in algorithm");
            }
            if (p0MinDistSq < p1MinDistSq) {
                // R1 <  E1End [*] [*] '0' '1'    <---- direction of increase
                if (p0MinIdx != 0) {
                    throw new IllegalStateException("error in algorithm");
                }
                if (routes.ep1 != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                LinkedHashSet<PairInt> tmp = new LinkedHashSet<PairInt>();
                tmp.add(segment.p0);
                tmp.add(segment.p1);
                tmp.addAll(routes.route1);
                routes.route1 = tmp;
            } else if (p1MinDistSq < p0MinDistSq) {
                // R1 <   '0' '1'  [*] [*] E1  <---- direction of increase
                if (p1MinIdx != 1) {
                    throw new IllegalStateException("error in algorithm");
                }
                if (routes.ep1End != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                routes.route1.add(segment.p1);
                routes.route1.add(segment.p0);
            } else {
                throw new IllegalStateException("error in algorithm");
            }
        } else if (segment.p1.getX() < segment.p0.getX()) {
            /*           swapped HorizSegment
                          -    -         -2
            R0 >       .  2:E0 3     .   -1    route0 incr w/ incr x
            R1 <       .  1    0:E1  .    0    route1 incr w/ decr x
                          -    -          1
                         -2    -1    0    1
            */
            if (!(segment.p2.getX() < segment.p3.getX())) {
                throw new IllegalStateException("error in algorithm");
            }
            if (p0MinDistSq < p1MinDistSq) {
                // R1 <   '1' '0'  [*] [*] E1 <---- direction of increase
                if (p0MinIdx != 1) {
                    throw new IllegalStateException("error in algorithm");
                }
                if (routes.ep1End != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                routes.route1.add(segment.p0);
                routes.route1.add(segment.p1);
            } else if (p1MinDistSq < p0MinDistSq) {
                // R1 <  E1end [*] [*] '1' '0'  <---- direction of increase
                if (p1MinIdx != 0) {
                    throw new IllegalStateException("error in algorithm");
                }
                if (routes.ep1 != null) {
                    throw new IllegalStateException("error in algorithm");
                }
                LinkedHashSet<PairInt> tmp = new LinkedHashSet<PairInt>();
                tmp.add(segment.p1);
                tmp.add(segment.p0);
                tmp.addAll(routes.route1);
                routes.route1 = tmp;
            } else {
                throw new IllegalStateException("error in algorithm");
            }
        } else {
            throw new IllegalStateException("error in algorithm");
        }
    }
    
    protected PairInt[] getFirstAndLast(Iterator<PairInt> iter) {
        PairInt firstNode = null;
        PairInt lastNode = null;
        while (iter.hasNext()) {
            if (firstNode == null) {
                firstNode = iter.next();
            } else {
                lastNode = iter.next();
            }
        }
        return new PairInt[]{firstNode, lastNode};
    }

    private void addUUDiagSegmentToRoutes(Routes routes, DiagSegment 
        segment) {
        
        assert(!routes.route0.isEmpty());
        assert(!routes.route1.isEmpty());
        
        /*      UUDiagSegment    
                    -  -      -2  
              -  2  1  -  -   -1  
              -  -  3  0  -    0  
                 -  -          1  
             -3 -2 -1  0  1         

               E0 is '1', and route0 
                  continues with     
                  point '0'          
               E1 is '3', and route1 
                  continues with '2' 
        */
        
        // ------- add to route0 ---------
        PairInt[] r0FirstLastNodes = getFirstAndLast(routes.route0.iterator());
        
        double p0MinDistSq = Double.MAX_VALUE;
        int p0MinIdx = -1;
        double p1MinDistSq = Double.MAX_VALUE;
        int p1MinIdx = -1;
        for (int i = 0; i < r0FirstLastNodes.length; ++i) {
            PairInt r0 = r0FirstLastNodes[i];
            int diffX = segment.p0.getX() - r0.getX();
            int diffY = segment.p0.getY() - r0.getY();
            double distSq = (diffX * diffX) + (diffY * diffY);
            if (distSq < p0MinDistSq) {
                p0MinDistSq = distSq;
                p0MinIdx = i;
            }
            diffX = segment.p1.getX() - r0.getX();
            diffY = segment.p1.getY() - r0.getY();
            distSq = (diffX * diffX) + (diffY * diffY);
            if (distSq < p1MinDistSq) {
                p1MinDistSq = distSq;
                p1MinIdx = i;
            }
        }
        
        if (p1MinDistSq < p0MinDistSq) {
            /*
            E0
            [*]
            [*]
            '1'
            '0'
            */
            if (p1MinIdx != 1) {
                throw new IllegalStateException("error in algorithm");
            }
            if (routes.ep0End != null) {
                throw new IllegalStateException("error in algorithm");
            }
            routes.route0.add(segment.p1);
            routes.route0.add(segment.p0);
        } else if (p0MinDistSq < p1MinDistSq) {
            /*
            '1'
            '0'
            [*]
            [*]
            E0end
            */
            if (p0MinIdx != 0) {
                throw new IllegalStateException("error in algorithm");
            }
            if (routes.ep0 != null) {
                throw new IllegalStateException("error in algorithm");
            }
            LinkedHashSet<PairInt> tmp = new LinkedHashSet<PairInt>();
            tmp.add(segment.p1);
            tmp.add(segment.p0);
            tmp.addAll(routes.route0);
            routes.route0 = tmp;
        } else {
            throw new IllegalStateException("error in algorithm");
        }
        
        // ------ add to route1 ------------
        PairInt[] r1FirstLastNodes = getFirstAndLast(routes.route1.iterator());
        double p2MinDistSq = Double.MAX_VALUE;
        int p2MinIdx = -1;
        double p3MinDistSq = Double.MAX_VALUE;
        int p3MinIdx = -1;
        for (int i = 0; i < r1FirstLastNodes.length; ++i) {
            PairInt r1 = r1FirstLastNodes[i];
            int diffX = segment.p2.getX() - r1.getX();
            int diffY = segment.p2.getY() - r1.getY();
            double distSq = (diffX * diffX) + (diffY * diffY);
            if (distSq < p2MinDistSq) {
                p2MinDistSq = distSq;
                p2MinIdx = i;
            }
            diffX = segment.p3.getX() - r1.getX();
            diffY = segment.p3.getY() - r1.getY();
            distSq = (diffX * diffX) + (diffY * diffY);
            if (distSq < p3MinDistSq) {
                p3MinDistSq = distSq;
                p3MinIdx = i;
            }
        }
        
        if (p2MinDistSq < p3MinDistSq) {
            /*
            E1end
            [*]
            [*]
            '2'
            '3'
            */
            if (p2MinIdx != 0) {
                throw new IllegalStateException("error in algorithm");
            }
            if (routes.ep1 != null) {
                throw new IllegalStateException("error in algorithm");
            }
            LinkedHashSet<PairInt> tmp = new LinkedHashSet<PairInt>();
            tmp.add(segment.p3);
            tmp.add(segment.p2);
            tmp.addAll(routes.route1);
            routes.route1 = tmp;
        } else if (p3MinDistSq < p2MinDistSq) {
            /*
            '2'
            '3'
            [*]
            [*]
            E1
            */
            if (p3MinIdx != 1) {
                throw new IllegalStateException("error in algorithm");
            }
            if (routes.ep1End != null) {
                throw new IllegalStateException("error in algorithm");
            }
            routes.route1.add(segment.p3);
            routes.route1.add(segment.p2);
        } else {
            throw new IllegalStateException("error in algorithm");
        }
    }

    private ZigZagSegmentRoutes parseZigZag(Segment segment, Set<PairInt> endPoints) {
        
        if (!(segment instanceof ZigZagSegment) && 
            !(segment instanceof ZigZagSegment2)) {
            throw new IllegalArgumentException(
            "segment type must be ZigZagSegment or ZigZagSegment2");
        }
        
        if (endPoints.size() != 4) {
            throw new IllegalStateException("endPoints must be size 4");
        }
        
        ZigZagSegmentRoutes routes = null;
        
        if (segment instanceof ZigZagSegment) {
            if (segment.p2.getX() < segment.p1.getX()) {
                
                /*
                       E0
                           2  -            -2
                           -  1 .  E0end   -1
                  E1end  . 0  -             0
                           -  3             1
                                 E1
                 -3 -2 -1  0  1
                */
                
                // '1' and '2' on one route and '0' and '3' on another.
                // add the endpoints before and after
                
                routes = new ZigZagSegmentRoutes();
                
                routes.ep1 = findClosestTo(segment.p3, endPoints);
                routes.route1.add(routes.ep1);
                routes.route1.add(segment.p3);
                routes.route1.add(segment.p0);
                routes.ep1End = findClosestTo(segment.p0, endPoints);
                routes.route1.add(routes.ep1End);
                
                routes.ep0 = findClosestTo(segment.p2, endPoints);
                routes.route0.add(routes.ep0);
                routes.route0.add(segment.p2);
                routes.route0.add(segment.p1);
                routes.ep0End = findClosestTo(segment.p1, endPoints);
                routes.route0.add(routes.ep0End);
                
            } else {
                
                /*            E0end
                        -  2          -2
                 E0   . 1  -          -1
                        -  0 . E1      0
                        3  -           1
                  E1end
                 -3 -2 -1  0  1
                */
                // '1' and '2' on one route and '0' and '3' on another.
                // add the endpoints before and after
                
                routes = new ZigZagSegmentRoutes();
                
                routes.ep0 = findClosestTo(segment.p1, endPoints);
                routes.route0.add(routes.ep0);
                routes.route0.add(segment.p1);
                routes.route0.add(segment.p2);
                routes.ep0End = findClosestTo(segment.p2, endPoints);
                routes.route0.add(routes.ep0End);
                
                routes.ep1 = findClosestTo(segment.p0, endPoints);
                routes.route1.add(routes.ep1);
                routes.route1.add(segment.p0);
                routes.route1.add(segment.p3);
                routes.ep1End = findClosestTo(segment.p3, endPoints);
                routes.route1.add(routes.ep1End);
                
            }
        } else {
            if (segment.p0.getY() < segment.p1.getY()) {
                
                /*      E1end     E0
                           .            -1
                        -  0  -  3       0
                        1  -  2  -       1
                     E1       .  E0end   2
                 -3 -2 -1  0  1  2  3
                */
                // '1' and '0' on one route and '2' and '3' on another.
                // add the endpoints before and after
                
                routes = new ZigZagSegmentRoutes();
                
                routes.ep0 = findClosestTo(segment.p3, endPoints);
                routes.route0.add(routes.ep0);
                routes.route0.add(segment.p3);
                routes.route0.add(segment.p2);
                routes.ep0End = findClosestTo(segment.p2, endPoints);
                routes.route0.add(routes.ep0End);
                
                routes.ep1 = findClosestTo(segment.p1, endPoints);
                routes.route1.add(routes.ep1);
                routes.route1.add(segment.p1);
                routes.route1.add(segment.p0);
                routes.ep1End = findClosestTo(segment.p0, endPoints);
                routes.route1.add(routes.ep1End);
                
            } else {
                
                /*
                    E1end    E0
                              .         -2
                        1  -  2  -      -1
                        -  0  -  3       0
                        E1 .    E0end    1
                                         2
                 -3 -2 -1  0  1  2  3
                */
                // '1' and '0' on one route and '2' and '3' on another.
                // add the endpoints before and after
                
                routes = new ZigZagSegmentRoutes();
                
                routes.ep0 = findClosestTo(segment.p2, endPoints);
                routes.route0.add(routes.ep0);
                routes.route0.add(segment.p2);
                routes.route0.add(segment.p3);
                routes.ep0End = findClosestTo(segment.p3, endPoints);
                routes.route0.add(routes.ep0End);
                
                routes.ep1 = findClosestTo(segment.p0, endPoints);
                routes.route1.add(routes.ep1);
                routes.route1.add(segment.p0);
                routes.route1.add(segment.p1);
                routes.ep1End = findClosestTo(segment.p1, endPoints);
                routes.route1.add(routes.ep1End);
            }
        }
        
        return routes;
    }

    private PairInt findClosestTo(PairInt p0, Set<PairInt> points) {
        
        PairInt pt = null;
        int minDistSq = Integer.MAX_VALUE;
        for (PairInt p : points) {
            int diffX = p.getX() - p0.getX();
            int diffY = p.getY() - p0.getY();
            int distSq = (diffX*diffX + diffY*diffY);
            if (distSq < minDistSq) {
                minDistSq = distSq;
                pt = p;
            }
        }
        
        return pt;
    }
    
    /*
    may change these classes to have ordered points or to specify the
    indexes of points that are the connections, that is the '.'s in sketches
    below.
    */
    public static class Segment {
        // the 4 points matching the segment as 0, 1, 2, 3 in the subclasses
        PairInt p0;
        PairInt p1;
        PairInt p2;
        PairInt p3;
        boolean contains(PairInt p) {
            if (p0.equals(p)) {
                return true;
            } else if (p1.equals(p)) {
                return true;
            } else if (p2.equals(p)) {
                return true;
            } else if (p3.equals(p)) {
                return true;
            }
            return false;
        }
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("p0=").append(p0.toString())
                .append(" p1=").append(p1.toString())
                .append(" p2=").append(p2.toString())
                .append(" p3=").append(p3.toString());
            return sb.toString();
        }
        public void filterOutContaining(Set<PairInt> points) {
            points.remove(p0);
            points.remove(p1);
            points.remove(p2);
            points.remove(p3);
        }
    }
    
    public static class VertSegment extends Segment {
    }
    public static class HorizSegment extends Segment {
    }
    public static class DiagSegment extends Segment {
    }
    public static class UUDiagSegment extends DiagSegment {
    }
    public static class ULDiagSegment extends DiagSegment {
    }
    public static class ZigZagSegment extends Segment {
    }
    public static class ZigZagSegment2 extends Segment {
    }

    public static class Pattern {
        Set<PairInt> ones;
        Set<PairInt> zeroes;
        /*
        an endpoint in route0.  can be null
        */
        PairInt ep0 = null;
        /*
        an endpoint in route1.  can be null
        */
        PairInt ep1 = null;
    }

    protected Pattern getVertSegmentPattern() {

        /*    .  .      -2
           -  2  1  -   -1
           -  3  0  -    0
              .  .       1
          -2 -1  0  1
        */
        Pattern pr = new Pattern();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-2, 0)); pr.zeroes.add(new PairInt(-2, -1));
        pr.zeroes.add(new PairInt(1, 0)); pr.zeroes.add(new PairInt(1, -1));

        pr.ones.add(new PairInt(-1, 0)); pr.ones.add(new PairInt(-1, -1));
        pr.ones.add(new PairInt(0, -1));

        return pr;
    }
    
    protected Pattern getUUDiagSegmentPattern() {

        /*      -  -      -2
          -  2  1  -  -   -1
          -  -  3  0  -    0
             -  -          1
         -3 -2 -1  0  1
        */
        Pattern pr = new Pattern();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-3, 0)); pr.zeroes.add(new PairInt(-3, -1));
        pr.zeroes.add(new PairInt(-2, 1)); pr.zeroes.add(new PairInt(-2, 0));
        pr.zeroes.add(new PairInt(-1, 1)); pr.zeroes.add(new PairInt(-1, -2));
        pr.zeroes.add(new PairInt(0, -1)); pr.zeroes.add(new PairInt(0, -2));
        pr.zeroes.add(new PairInt(1, 0)); pr.zeroes.add(new PairInt(1, -1));

        pr.ones.add(new PairInt(-2, -1)); 
        pr.ones.add(new PairInt(-1, 0)); pr.ones.add(new PairInt(-1, -1));

        return pr;
    }

    protected Pattern getZigZagSegmentPattern() {

        /*
                   2  -      -2
                   -  1      -1
                   0  -       0
                   -  3       1
                       
         -3 -2 -1  0  1
        
        Each of the 4 needs at least one neighbor that is not one of the 4
        points in the zig zap and all of their neighbors cannot be adjacent
        to any of the other neighbors.
        */
        Pattern pr = new Pattern();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(0, 1)); pr.zeroes.add(new PairInt(0, -1));
        pr.zeroes.add(new PairInt(1, 0)); pr.zeroes.add(new PairInt(1, -2));

        pr.ones.add(new PairInt(0, -2));
        pr.ones.add(new PairInt(1, 1)); pr.ones.add(new PairInt(1, -1));

        return pr;
    }
    
    protected Pattern getZigZag2SegmentPattern() {

        /*
                   .            -1
                -  0  -  3       0
                1  -  2  -       1
                      .          2
         -3 -2 -1  0  1  2  3
        
        Each of the 4 needs at least one neighbor that is not one of the 4
        points in the zig zap and all of their neighbors cannot be adjacent
        to any of the other neighbors.
        */
        Pattern pr = new Pattern();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, 0)); 
        pr.zeroes.add(new PairInt(0, 1));
        pr.zeroes.add(new PairInt(1, 0)); 
        pr.zeroes.add(new PairInt(2, 1));

        pr.ones.add(new PairInt(-1, 1));
        pr.ones.add(new PairInt(0, 0)); pr.ones.add(new PairInt(0, -1));
        pr.ones.add(new PairInt(1, 1)); pr.ones.add(new PairInt(1, 2)); 
        pr.ones.add(new PairInt(2, 0));

        return pr;
    }
    
    private Segment checkVertHorizSegmentPattern(int x, int y, 
        Set<PairInt> neighbors, boolean useVertical) {

        /*
            VertSegment         swapY VertSegment
              .  .      -2
           -  2  1  -   -1         .  .      -1
           -  3  0  -    0      -  3  0  -    0
              .  .       1      -  2  1  -    1
                                   .  .       2
          -2 -1  0  1          -2 -1  0  1

            HorizSegment        swap HorizSegment
                 -  -     -2       -  -      -2
              .  3  2  .  -1    .  2  3  .   -1
              .  0  1  .   0    .  1  0  .    0
                 -  -      1       -  -       1
             -1  0  1  2       -2 -1  0  1
        */
        
        Pattern pattern = getVertSegmentPattern();
        
        if (!useVertical) {
            rotatePattern(pattern, -0.5*Math.PI);
        }

        boolean matchesPattern = matchesPattern(x, y, neighbors, pattern);
        
        if (matchesPattern) {
            if (useVertical) {
                VertSegment segment = new VertSegment();
                segment.p0 = new PairInt(x, y);
                segment.p1 = new PairInt(x, y - 1);
                segment.p2 = new PairInt(x - 1, y - 1);
                segment.p3 = new PairInt(x - 1, y);
                return segment;
            } else {
                HorizSegment segment = new HorizSegment();
                segment.p0 = new PairInt(x, y);
                segment.p1 = new PairInt(x + 1, y);
                segment.p2 = new PairInt(x + 1, y - 1);
                segment.p3 = new PairInt(x, y - 1);
                return segment;
            }            
        }
                
        if (useVertical) {
            swapYDirection(pattern);
        } else {
            pattern = getVertSegmentPattern();
            swapYDirection(pattern);
            rotatePattern(pattern, -0.5*Math.PI);
        }
        
        matchesPattern = matchesPattern(x, y, neighbors, pattern);
        
        if (matchesPattern) {
            if (useVertical) {
                VertSegment segment = new VertSegment();
                segment.p0 = new PairInt(x, y);
                segment.p1 = new PairInt(x, y + 1);
                segment.p2 = new PairInt(x - 1, y + 1);
                segment.p3 = new PairInt(x - 1, y);
                return segment;
            } else {
                HorizSegment segment = new HorizSegment();
                segment.p0 = new PairInt(x, y);
                segment.p1 = new PairInt(x - 1, y);
                segment.p2 = new PairInt(x - 1, y - 1);
                segment.p3 = new PairInt(x, y - 1);
                return segment;
            }
        }
        
        return null;
    }
    
    private Segment checkDiagSegmentPattern(int x, int y, Set<PairInt> neighbors) {

        Pattern pattern = getUUDiagSegmentPattern();
        
        boolean matchesPattern = matchesPattern(x, y, neighbors, pattern);
        
        /*      -  -      -2
          -  2  1  -  -   -1
          -  -  3  0  -    0
             -  -          1
         -3 -2 -1  0  1
        */
        if (matchesPattern) {
            UUDiagSegment segment = new UUDiagSegment();
            segment.p0 = new PairInt(x, y);
            segment.p1 = new PairInt(x - 1, y - 1);
            segment.p2 = new PairInt(x - 2, y - 1);
            segment.p3 = new PairInt(x - 1, y);
            
            return segment;
        }
        
        rotatePattern(pattern, -0.5*Math.PI);
        
        matchesPattern = matchesPattern(x, y, neighbors, pattern);
        
        /*  
                      2         -2
                   3  1         -1
                   0             0
                                 1
         -3 -2 -1  0  1  2  3
        */
        if (matchesPattern) {
            ULDiagSegment segment = new ULDiagSegment();
            segment.p0 = new PairInt(x, y);
            segment.p1 = new PairInt(x + 1, y - 1);
            segment.p2 = new PairInt(x + 1, y - 2);
            segment.p3 = new PairInt(x, y - 1);
            
            return segment;
        }
        
        return null;
    }
    
    private ZigZagSegment checkZigZagSegmentPattern(int x, int y, 
        Set<PairInt> points) {

        Pattern pattern = getZigZagSegmentPattern();
        
        boolean matchesPattern = matchesPattern(x, y, points, pattern);
        
        /*
                   2  -      -2
                   -  1      -1
                   0  -       0
                   -  3       1
                       
         -3 -2 -1  0  1
        */
        if (matchesPattern) {
            ZigZagSegment segment = new ZigZagSegment();
            segment.p0 = new PairInt(x, y);
            segment.p1 = new PairInt(x + 1, y - 1);
            segment.p2 = new PairInt(x, y - 2);
            segment.p3 = new PairInt(x + 1, y + 1);
            
            return segment;
        }
        
        swapXDirection(pattern);
        
        matchesPattern = matchesPattern(x, y, points, pattern);
        
        /*
                -  2      -2
                1  -      -1
                -  0       0
                3  -       1
                       
         -3 -2 -1  0  1
        */
        if (matchesPattern) {
            ZigZagSegment segment = new ZigZagSegment();
            segment.p0 = new PairInt(x, y);
            segment.p1 = new PairInt(x - 1, y - 1);
            segment.p2 = new PairInt(x, y - 2);
            segment.p3 = new PairInt(x - 1, y + 1);
            
            return segment;
        }
        
        return null;
    }
    
    private ZigZagSegment2 checkZigZag2SegmentPattern(int x, int y, 
        Set<PairInt> points) {

        Pattern pattern = getZigZag2SegmentPattern();
        
        boolean matchesPattern = matchesPattern(x, y, points, pattern);
        
        /*
                   .            -1
                -  0  -  3       0
                1  -  2  -       1
                      .          2
         -3 -2 -1  0  1  2  3
        */
        if (matchesPattern) {
            ZigZagSegment2 segment = new ZigZagSegment2();
            segment.p0 = new PairInt(x, y);
            segment.p1 = new PairInt(x - 1, y + 1);
            segment.p2 = new PairInt(x + 1, y + 1);
            segment.p3 = new PairInt(x + 2, y);
            
            return segment;
        }
        
        swapYDirection(pattern);
        
        matchesPattern = matchesPattern(x, y, points, pattern);
        
        /*            .         -2
                1  -  2  -      -1
                -  0  -  3       0
                   .             1
                                 2
         -3 -2 -1  0  1  2  3
        */
        if (matchesPattern) {
            ZigZagSegment2 segment = new ZigZagSegment2();
            segment.p0 = new PairInt(x, y);
            segment.p1 = new PairInt(x - 1, y - 1);
            segment.p2 = new PairInt(x + 1, y - 1);
            segment.p3 = new PairInt(x + 2, y);
            
            return segment;
        }
        
        return null;
    }
    
    private boolean matchesPattern(final int x, final int y, Set<PairInt> neighbors, 
        Pattern pattern) {
                                                
        for (PairInt p : pattern.zeroes) {
            PairInt p2 = new PairInt(x + p.getX(), y + p.getY());
            if (neighbors.contains(p2)) {
                return false;
            }
        }
        for (PairInt p : pattern.ones) {
            PairInt p2 = new PairInt(x + p.getX(), y + p.getY());
            if (!neighbors.contains(p2)) {
                return false;
            }
        }
        return true;
    }
    
    /**
     * route0 and route1 are the two routes in opposite directions for a section
     * in a closed curve with a junction. route0 and route1 do not cross and are
     * populated to help populate the curve so that the points are traversed in
     * opposite directions.
     */
    public static class Routes {
        PairInt ep0 = null;
        PairInt ep1 = null;
        PairInt ep0End = null;
        PairInt ep1End = null;
        //route0, route1 need to be searchable but ordered.
        LinkedHashSet<PairInt> route0 = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> route1 = new LinkedHashSet<PairInt>();
        public LinkedHashSet<PairInt> getRoute0() {
            return route0;
        }
        public LinkedHashSet<PairInt> getRoute1() {
            return route1;
        }
    }
    public static class VertSegmentRoutes extends Routes {
    }
    public static class HorizSegmentRoutes extends Routes {
    }
    public static class UUDiagSegmentRoutes extends Routes {
    }
    public static class ZigZagSegmentRoutes extends Routes {
    }
}
