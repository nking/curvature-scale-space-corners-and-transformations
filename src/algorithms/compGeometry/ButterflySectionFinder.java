package algorithms.compGeometry;

import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
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

        setNullEndpoints(output);

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
            if (neighbors.size() != 5 && neighbors.size() != 4) {
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

        //TODO: edit to handle diagonal too when tests start
        mergeIfAdjacent(candidateSections);
        
        List<Routes> output = new ArrayList<Routes>();

        // -- scan for endpoints --
        for (LinkedList<Segment> section : candidateSections) {

            Routes routes = checkForAdjacentEndpoints(points, section);

            if (routes == null) {
                continue;
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

        List<Routes> routesList = findButterflySectionsSmallZigZag(closedCurve,
            points);

        if (routesList != null && !routesList.isEmpty()) {
            return routesList;
        }

        routesList = findButterflySectionsSmallDiagZigZag(closedCurve,
            points);

        return routesList;
    }

    protected List<Routes> findButterflySectionsSmallZigZag(PairIntArray
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

    protected List<Routes> findButterflySectionsSmallDiagZigZag(PairIntArray
        closedCurve, Set<PairInt> points) {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        List<Routes> output = new ArrayList<Routes>();

        for (int i = 0; i < closedCurve.getN(); ++i) {

            int x = closedCurve.getX(i);
            int y = closedCurve.getY(i);

            Set<PairInt> neighbors = curveHelper.findNeighbors(x, y, points);

            //TODO: one point in this pattern has 3 neigbhors, so may need to change this

            // scanning for segments
            if (neighbors.size() != 4) {
                continue;
            }

            Segment segment = checkDiagZigZagSegmentPattern(x, y, points);

            if (segment == null) {
                segment = checkDiagZigZag2SegmentPattern(x, y, points);
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

            ZigZagSegmentRoutes routes = parseDiagZigZag(segment, endPoints);

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

        Routes routes = null;
        
        for (int i = 0; i < 2; ++i) {
        
            Segment firstOrLastSegment = (i == 0) ? section.getFirst() :
                section.getLast();            

            if (firstOrLastSegment instanceof VertSegment) {

                boolean useVertical = true;

                routes = findEndPointsVertHorizPatterns(points, routes, firstOrLastSegment,
                    useVertical);

            } else if (firstOrLastSegment instanceof HorizSegment) {

                boolean useVertical = false;

                routes = findEndPointsVertHorizPatterns(points, routes, firstOrLastSegment,
                    useVertical);

            } else if (firstOrLastSegment instanceof UUDiagSegment) {

                routes = findEndPointsDiagPatterns(points, routes,
                    (UUDiagSegment)firstOrLastSegment);
            }

            if (routes == null) {
                return routes;
            }

            if ((i == 0) && section.size() > 1) {

                int added = 0;

                for (int ii = 1; ii < section.size(); ++ii) {

                    Segment segment = section.get(ii);

                    // add node to routes
                    int didAdd = addSegmentToRoutes(routes, segment);

                    added = added | didAdd;
                }

                if (!((added == 1) && (routes.ep0 == null || routes.ep0End == null ||
                    routes.ep1 == null || routes.ep1End == null))) {

                    break;
                }
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
                    pattern = getEndPointsVertPattern1S();
                    break;
                case 1:
                    pattern = getEndPointsVertPattern1Opp();
                    break;
                case 2:
                    pattern = getEndPointsVertPattern2S();
                    break;
                case 3:
                    pattern = getEndPointsVertPattern2Opp();
                    break;
                case 4:
                    pattern = getEndPointsVertPattern3S();
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
                    if ((i & 1) == 0) {
                        addEndpointsForHorizVertPatternForward(x0, y0,
                            pattern, routes, endPoints, useVert);
                    } else {
                        addEndpointsForHorizVertPatternOppos(x0, y0,
                            pattern, routes, endPoints, useVert);
                    }
                    break;
                }
            }
        }
        return routes;
    }

    /*
    private Pattern getEndPointsVertPattern1() {
        // the pattern returned is relative to
        //position '0', just like the other patterns.

        //searching for .'s and #'s
        //                -2
        //   -  2  1  -   -1
        //   -  3  0  -    0
        //   -  .  .  -    1
        //      #  -  #    2
        //         -       3
        //  -2 -1  0  1
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
    */
    
    private Pattern getEndPointsVertPattern1S() {
        /* the pattern returned is relative to
        position '0', just like the other patterns.

        searching for .'s and #'s
                        -2
           -  2  1  -   -1
           -  3  0  -    0
              #  -  #    1
                 -       2
          -2 -1  0  1
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();

        // '#' end points
        PairInt t1 = new PairInt(-1, 1);
        pattern.ones.add(t1);
        PairInt t0 = new PairInt(1, 1);
        pattern.ones.add(t0);
        pattern.ep0 = t0;
        pattern.ones.add(t1);
        pattern.ep1 = t1;

        pattern.zeroes.add(new PairInt(-2, -1)); pattern.zeroes.add(new PairInt(-2, 0));
        pattern.zeroes.add(new PairInt(0, 2)); pattern.zeroes.add(new PairInt(0, 1));
        pattern.zeroes.add(new PairInt(1, 0)); pattern.zeroes.add(new PairInt(1, -1));
        return pattern;
    }

    private Pattern getEndPointsVertPattern1Opp() {
        /* the pattern returned is relative to
        position '0', just like the other patterns.

             R1    R0
             /|\    |
              |    \|/
                 -      -3
              #  -  #   -2
           -  2  1  -   -1
              3  0       0

          -2 -1  0  1
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();

        PairInt t1 = new PairInt(-1, -2);
        pattern.ones.add(t1);
        pattern.ep1 = t1;
        PairInt t0 = new PairInt(1, -2);
        pattern.ones.add(t0);
        pattern.ep0 = t0;

        pattern.zeroes.add(new PairInt(-2, -1));
        pattern.zeroes.add(new PairInt(0, -2)); pattern.zeroes.add(new PairInt(0, -3));
        pattern.zeroes.add(new PairInt(1, -1));
        return pattern;
    }

    /*
    private Pattern getEndPointsVertPattern2() {

        // the pattern returned is relative to
        //position '0', just like the other patterns.

        //searching for .'s and #'s
        //                -2
        //   -  2  1  -   -1
        //   -  3  0  -    0
        //   -  .  .  -    1
        //   #  -  -  #    2
        //      -  -       3
        //  -2 -1  0  1
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
    */
    
    private Pattern getEndPointsVertPattern2S() {

        /* the pattern returned is relative to
        position '0', just like the other patterns.

        searching for .'s and #'s
                        -2
              2  1      -1
           -  3  0  -    0
           #  -  -  #    1
              -  -       2
          -2 -1  0  1
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();

        // '#' end points
        PairInt t1 = new PairInt(-2, 1);
        pattern.ones.add(t1);
        PairInt t0 = new PairInt(1, 1);
        pattern.ones.add(t0);
        pattern.ep1 = t1;
        pattern.ep0 = t0;

        pattern.zeroes.add(new PairInt(-2, 0));
        pattern.zeroes.add(new PairInt(-1, 2)); pattern.zeroes.add(new PairInt(-1, 1));
        pattern.zeroes.add(new PairInt(0, 2)); pattern.zeroes.add(new PairInt(0, 1));
        pattern.zeroes.add(new PairInt(1, 0));
        
        return pattern;
    }

    private Pattern getEndPointsVertPattern2Opp() {

        /* the pattern returned is relative to
        position '0', just like the other patterns.

          R1       R0
          /|\       |
           |       \|/
              -  -      -3
           #  -  -  #   -2
           -  2  1  -   -1
              3  0       0

          -2 -1  0  1
        */

        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();

        // '#' end points
        PairInt t1 = new PairInt(-2, -2);
        PairInt t0 = new PairInt(1, -2);
        pattern.ones.add(t1);
        pattern.ones.add(t0);
        pattern.ep0 = t0;
        pattern.ep1 = t1;

        pattern.zeroes.add(new PairInt(-2, -1));
        pattern.zeroes.add(new PairInt(-1, -2)); pattern.zeroes.add(new PairInt(-1, -3));
        pattern.zeroes.add(new PairInt(0, -2)); pattern.zeroes.add(new PairInt(0, -3));
        pattern.zeroes.add(new PairInt(1, -1));
        return pattern;
    }
    
    /*
    private Pattern getEndPointsVertPattern3() {

        // the pattern returned is relative to
        //position '0', just like the other patterns.

        //searching for .'s and #'s
        //                -2
        //   -  2  1  -   -1
        //   -  3  0  -    0
        //   -  .  .  -    1
        //   #  -  #       2
        //      -          3
        //  -2 -1  0  1

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
    */
    
    private Pattern getEndPointsVertPattern3S() {

        /* the pattern returned is relative to
        position '0', just like the other patterns.

        searching for .'s and #'s
                        -2
           -  2  1  -   -1
           -  3  0  -    0
           #  -  #       1
              -          2
          -2 -1  0  1
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();

        pattern.zeroes.add(new PairInt(-2, -1)); pattern.zeroes.add(new PairInt(-2, 0));
        pattern.zeroes.add(new PairInt(-1, 2)); pattern.zeroes.add(new PairInt(-1, 1));
        pattern.zeroes.add(new PairInt(1, 0)); pattern.zeroes.add(new PairInt(1, -1));

        // '#' end points
        PairInt t1 = new PairInt(-2, 1);
        pattern.ones.add(t1);
        PairInt t0 = new PairInt(0, 1);
        pattern.ones.add(t0);
        pattern.ep1 = t1;
        pattern.ep0 = t0;

        return pattern;
    }

    private Pattern getEndPointsVertPattern3Opp() {

        /* the pattern returned is relative to
        position '0', just like the other patterns.
          R1    R0
          /|\    |
           |    \|/
              -         -3
           #  -  #      -2
           -  2  1  -   -1
              3  0       0

          -2 -1  0  1
        */

        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();

        // '#' end points
        PairInt t1 = new PairInt(-2, -2);
        PairInt t0 = new PairInt(0, -2);
        pattern.ones.add(t1);
        pattern.ones.add(t0);
        pattern.ep0 = t0;
        pattern.ep1 = t1;

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

        Set<Integer> indexesWithMoreThanOneNeighbor = new HashSet<Integer>();

        List<Set<PairInt>> listOfNeighbors = new ArrayList<Set<PairInt>>();
        Set<PairInt> segmentPoints = new HashSet<PairInt>();
        segmentPoints.add(segment.p0);
        segmentPoints.add(segment.p1);
        segmentPoints.add(segment.p2);
        segmentPoints.add(segment.p3);
        for (int i = 0; i < 4; ++i) {
            PairInt p = null;
            switch(i) {
                case 0:
                    p = segment.p0;
                    break;
                case 1:
                    p = segment.p1;
                    break;
                case 2:
                    p = segment.p2;
                    break;
                default:
                    p = segment.p3;
                    break;
            }
            Set<PairInt> neighbors = curveHelper.findNeighbors(p.getX(),
                p.getY(), points);
            neighbors.removeAll(segmentPoints);
            if (!((neighbors.size() == 1) || (neighbors.size() == 2))) {
                return null;
            }
            if (neighbors.size() > 1) {
                indexesWithMoreThanOneNeighbor.add(
                    Integer.valueOf(listOfNeighbors.size()));
            }
            listOfNeighbors.add(neighbors);
        }

        /* for this kind of junction, some points share a neighbor.
        need to reduce the neighbor lists to unique members
        */
        for (Integer index : indexesWithMoreThanOneNeighbor) {
            int idx = index.intValue();
            Set<PairInt> neighborSet = listOfNeighbors.get(idx);
            for (int i = 0; i < listOfNeighbors.size(); ++i) {
                Integer index2 = Integer.valueOf(i);
                if (indexesWithMoreThanOneNeighbor.contains(index2)) {
                    continue;
                }
                PairInt p2 = listOfNeighbors.get(index2.intValue()).iterator().next();
                neighborSet.remove(p2);
            }
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
        
        // if there are more than one per set, we only want to keep the
        // closest to make a clear route with.
        // TODO: revisit this decision w/ further testing
        
        Set<PairInt> output = new HashSet<PairInt>();
        for (int i = 0; i < listOfNeighbors.size(); ++i) {
            Set<PairInt> set = listOfNeighbors.get(i);
            if (set.size() > 0) {
                PairInt p = null;
                switch(i) {
                    case 0:
                        p = segment.p0;
                        break;
                    case 1:
                        p = segment.p1;
                        break;
                    case 2:
                        p = segment.p2;
                        break;
                    default:
                        p = segment.p3;
                        break;
                }
                int minDistSq = Integer.MAX_VALUE;
                PairInt minDistP = null;
                for (PairInt p2 : set) {
                    int distSq = distSq(p, p2);
                    if (distSq < minDistSq) {
                        minDistSq = distSq;
                        minDistP = p2;
                    }
                }
                set.clear();
                set.add(minDistP);
            }
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
            VertSegment
             R1  R0
             /|\ |
              |  |
              | \|/
              .  .      -2
           -  2  1  -   -1
           -  3  0  -    0
              .  .       1
          -2 -1  0  1

               HorizPattern
                    -  -  -   -2
            R1-->   .  3  2   -1  -->R1 ends
                    .  0  1    0
            R0<--   -  -  -    1  <--R0 starts
                   -1  0  1
        
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
            assert(segment.p1.equals(new PairInt(x0, y0 - 1)));
            assert(segment.p2.equals(new PairInt(x0 - 1, y0 - 1)));
            assert(segment.p3.equals(new PairInt(x0 - 1, y0)));
            routes.route0.add(segment.p1);
            routes.route0.add(segment.p0);
            routes.route1.add(segment.p3);
            routes.route1.add(segment.p2);
            
            return routes;
            
        } else if (segment instanceof HorizSegment) {
            Routes routes = new Routes();
            assert(segment.p1.equals(new PairInt(x0 + 1, y0)));
            assert(segment.p2.equals(new PairInt(x0 + 1, y0 - 1)));
            assert(segment.p3.equals(new PairInt(x0, y0 - 1)));
            routes.route0.add(segment.p1);
            routes.route0.add(segment.p0);
            routes.route1.add(segment.p3);
            routes.route1.add(segment.p2);
            
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

        if (endPoints.size() < 2) {
            throw new IllegalStateException("error in algorithm");
        }
        assert(routes.ep0End == null);
        assert(routes.ep1 == null);

        /*
           endPointsVertPattern1, for example
                 R0
                 \/
                        -2
           -  2  1  -   -1
           -  3  0  -    0
           -  .  .  -    1
              #  -  #    2
                 -       3
             /\
             R1
          -2 -1  0  1
        
            endPointsHorizPattern1
                  -  -  -   -2
        R1--> #   .  3  2   -1  -->R1 ends
          -   -   .  0  1    0
        R0<-- #   -  -  -    1  <--R0 starts
                             2
                             3
         -3  -2  -1  0  1
        */
        assert(pattern.ep0 != null);
        PairInt tE1 = new PairInt(x0 + pattern.ep1.getX(), y0 + pattern.ep1.getY());
        PairInt tE0 = new PairInt(x0 + pattern.ep0.getX(), y0 + pattern.ep0.getY());
        if (!endPoints.contains(tE1) || !endPoints.contains(tE0)) {
            throw new IllegalStateException("error in algorithm");
        }
        assert(!routes.route0.contains(tE0));
        assert(!routes.route1.contains(tE1));

        routes.route0.add(tE0);
        routes.ep0End = tE0;

        LinkedHashSet<PairInt> tmp1 = new LinkedHashSet<PairInt>();
        tmp1.add(tE1);
        tmp1.addAll(routes.route1);
        routes.route1 = tmp1;
        routes.ep1 = tE1;
        
    }

    private void addEndpointsForHorizVertPatternOppos(int x0, int y0,
        Pattern pattern, Routes routes, Set<PairInt> endPoints, boolean useVert) {

        if (endPoints.size() < 2) {
            throw new IllegalStateException("error in algorithm");
        }
        assert(routes.ep0 == null);
        assert(routes.ep1End == null);

        /*
         VertSegment
             R1  R0
             /|\ |
              |  |
              | \|/
              .  .      -2
           -  2  1  -   -1
           -  3  0  -    0
              .  .       1
          -2 -1  0  1

               HorizPattern
                    -  -  -   -2
            R1-->   .  3  2   -1  -->R1 ends
                    .  0  1    0
            R0<--   -  -  -    1  <--R0 starts
                   -1  0  1
        */
        PairInt tE1 = new PairInt(x0 + pattern.ep1.getX(), y0 + pattern.ep1.getY());
        PairInt tE0 = new PairInt(x0 + pattern.ep0.getX(), y0 + pattern.ep0.getY());
        if (!endPoints.contains(tE1) || !endPoints.contains(tE0)) {
            throw new IllegalStateException("error in algorithm");
        }
        assert(!routes.route0.contains(tE0));
        assert(!routes.route1.contains(tE1));

        routes.route1.add(tE1);
        routes.ep1End = tE1;
        
        LinkedHashSet<PairInt> tmp0 = new LinkedHashSet<PairInt>();
        tmp0.add(tE0);
        tmp0.addAll(routes.route0);
        routes.route0 = tmp0;
        routes.ep0 = tE0;
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

        PairInt tI = new PairInt(x0 - 1, y0 + 1);
        PairInt tE0 = new PairInt(x0 + pattern.ep0.getX(), y0 + pattern.ep0.getY());
        PairInt tE1 = new PairInt(x0 + pattern.ep1.getX(), y0 + pattern.ep1.getY());
        if (!endPoints.contains(tI) && !routes.getRoute0().contains(tI)) {
            throw new IllegalStateException("error in algorithm");
        }
        if (!endPoints.contains(tE0) || !endPoints.contains(tE1)) {
            throw new IllegalStateException("error in algorithm");
        }
        if (routes.getRoute0().isEmpty()) {
            routes.route0.add(tE0);
        } else if (routes.getRoute0().size() > 1) {
            PairInt[] secondToLastAndLast = getSecondToLastAndLast(
                routes.getRoute0().iterator());
            if (secondToLastAndLast[1].equals(tE0)) {
                // no need to add either to route
            } else {
                routes.route0.add(tE0);
            }
        }
        routes.ep0End = tE0;

        if (routes.getRoute1().isEmpty()) {
            routes.route1.add(tE1);
            routes.route1.add(tI);
        } else if (routes.getRoute1().size() > 1) {
            Iterator<PairInt> iter1 = routes.getRoute1().iterator();
            PairInt r1 = iter1.next();
            PairInt r2 = iter1.next();
            if (r1.equals(tE1) && r2.equals(tI)) {
                // no need to add to route
            } else if (r1.equals(tI)) {
                LinkedHashSet<PairInt> tmpR1 = new LinkedHashSet<PairInt>();
                tmpR1.add(tE1);
                tmpR1.addAll(routes.route1);
                routes.route1 = tmpR1;
            } else {
                LinkedHashSet<PairInt> tmpR1 = new LinkedHashSet<PairInt>();
                tmpR1.add(tE1);
                tmpR1.add(tI);
                tmpR1.addAll(routes.route1);
                routes.route1 = tmpR1;
            }
        }
        routes.ep1 = tE1;
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

        PairInt tI = new PairInt(x0 - 1, y0 - 1);
        PairInt tE0 = new PairInt(x0 + pattern.ep0.getX(), y0 + pattern.ep0.getY());
        PairInt tE1 = new PairInt(x0 + pattern.ep1.getX(), y0 + pattern.ep1.getY());
        if (!endPoints.contains(tI) && !routes.getRoute0().contains(tI)) {
            throw new IllegalStateException("error in algorithm");
        }
        if (!endPoints.contains(tE0) || !endPoints.contains(tE1)) {
            throw new IllegalStateException("error in algorithm");
        }
        if (routes.getRoute0().isEmpty()) {
            routes.route0.add(tE0);
        } else if (routes.getRoute0().size() > 1) {
            PairInt[] secondToLastAndLast = getSecondToLastAndLast(
                routes.getRoute0().iterator());
            if (secondToLastAndLast[1].equals(tE0)) {
                // no need to add either to route
            } else {
                routes.route0.add(tE0);
            }
        }
        routes.ep0End = tE0;

        if (routes.getRoute1().isEmpty()) {
            routes.route1.add(tE1);
            routes.route1.add(tI);
        } else if (routes.getRoute1().size() > 1) {
            Iterator<PairInt> iter1 = routes.getRoute1().iterator();
            PairInt r1 = iter1.next();
            PairInt r2 = iter1.next();
            if (r1.equals(tE1) && r2.equals(tI)) {
                // no need to add to route
            } else if (r1.equals(tI)) {
                LinkedHashSet<PairInt> tmpR1 = new LinkedHashSet<PairInt>();
                tmpR1.add(tE1);
                tmpR1.addAll(routes.route1);
                routes.route1 = tmpR1;
            } else {
                LinkedHashSet<PairInt> tmpR1 = new LinkedHashSet<PairInt>();
                tmpR1.add(tE1);
                tmpR1.add(tI);
                tmpR1.addAll(routes.route1);
                routes.route1 = tmpR1;
            }
        }
        routes.ep1 = tE1;

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

    /**
     * add segment to routes and return 1 if did, else 0 (0 occurs when 
     * segment is already part of routes).
     * @param routes
     * @param segment
     * @return 
     */
    private int addSegmentToRoutes(Routes routes, Segment segment) {

        if ((segment instanceof VertSegment) || (segment instanceof HorizSegment)) {
            return addVertHorizSegmentToRoutes(routes, segment);
        } else if (segment instanceof UUDiagSegment) {
            return addUUDiagSegmentToRoutes(routes, (UUDiagSegment)segment);
        } else if (segment instanceof ULDiagSegment) {
            return addUUDiagSegmentToRoutes(routes, (ULDiagSegment)segment);
        }

        return 0;
    }

    private int addVertHorizSegmentToRoutes(Routes routes, Segment segment) {

        assert(!routes.route0.isEmpty());
        assert(!routes.route1.isEmpty());
        
        /*   Vertical pattern
              R1 R0
             /|\ |
              |  |
              | \|/
              .  .      -2
           -  2  1  -   -1
           -  3  0  -    0
                         1  
          -2 -1  0  1    
        
                HorizPattern 
                    -  -  -   -2
            R1-->   .  3  2   -1  -->R1 ends
                    .  0  1    0
            R0<--   -  -  -    1  <--R0 starts
                   -1  0  1
        */ 
        
        // if segment is already in routes, return
        if (segment instanceof HorizSegment) {
            if (routes.route0.contains(segment.p0) && routes.route0.contains(segment.p1)
                && routes.route1.contains(segment.p2) && routes.route1.contains(segment.p3)) {
                return 0;
            }
        } else if (segment instanceof VertSegment) {
            if (routes.route1.contains(segment.p2) && routes.route1.contains(segment.p3)
                && routes.route0.contains(segment.p1) && routes.route0.contains(segment.p0)) {
                return 0;
            }
        } else {
            throw new IllegalArgumentException(
            "segment must be type HorizSegment or VertSegment");
        }
      
        // -------- add to routes1 -------------------------
        PairInt[] r1FirstLastNodes = getFirstAndLast(routes.route1.iterator());

        double p3MinDistSq = Double.MAX_VALUE;
        int p3MinIdx = -1;
        double p2MinDistSq = Double.MAX_VALUE;
        int p2MinIdx = -1;
        for (int i = 0; i < r1FirstLastNodes.length; ++i) {
            PairInt r1 = r1FirstLastNodes[i];
            int distSq = distSq(segment.p3, r1);
            if (distSq < p3MinDistSq) {
                p3MinDistSq = distSq;
                p3MinIdx = i;
            }
            distSq = distSq(segment.p2, r1);
            if (distSq < p2MinDistSq) {
                p2MinDistSq = distSq;
                p2MinIdx = i;
            }
        }

        if (segment instanceof HorizSegment) {
            assert(segment.p3.getX() < segment.p2.getX());
            assert(segment.p0.getX() < segment.p1.getX());
        } else if (segment instanceof VertSegment) {
            assert(segment.p3.getY() > segment.p2.getY());
            assert(segment.p0.getY() > segment.p1.getY());
        }
        
        assert((p2MinDistSq == 0) || (p2MinDistSq == 1) || (p2MinDistSq == 4));
        assert((p3MinDistSq == 0) || (p3MinDistSq == 1) || (p3MinDistSq == 4));
        
        /*   Vertical pattern
              R1 R0
             /|\ |
              |  |
              | \|/
              .  .      -2
           -  2  1  -   -1
           -  3  0  -    0
                         1  
          -2 -1  0  1    
        
                HorizPattern 
                    -  -  -   -2
            R1-->   .  3  2   -1  -->R1 ends
                    .  0  1    0
            R0<--   -  -  -    1  <--R0 starts
                   -1  0  1
        */   
        
        if ((p2MinDistSq > p3MinDistSq) && (p3MinIdx == 1)) {
            // append to end of route1             R1:[start  stop]  p3  p2 
            //and insert at beginning of route0    R0:[stop  start]  p0  p1 
            LinkedHashSet<PairInt> tmp0 = new LinkedHashSet<PairInt>();
            tmp0.add(segment.p1);
            if (p3MinDistSq > 0) {
                // append p2
                routes.route1.add(segment.p3);
                tmp0.add(segment.p0);
            }
            routes.route1.add(segment.p2);
            tmp0.addAll(routes.route0);
            routes.route0 = tmp0;
        } else if ((p2MinDistSq < p3MinDistSq) && (p2MinIdx == 0)) {            
            // insert at beginning of route1    p3  p2  R1:[start  stop]
            // append to end of route0          p0  p1  R0:[stop  start]
            LinkedHashSet<PairInt> tmp1 = new LinkedHashSet<PairInt>();
            tmp1.add(segment.p3);
            if (p2MinDistSq > 0) {
                tmp1.add(segment.p2);
                routes.route0.add(segment.p1);
            }
            routes.route0.add(segment.p0);
            tmp1.addAll(routes.route1);
            routes.route1 = tmp1;            
        }
        
        return 1;
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

    protected PairInt[] getSecondToLastAndLast(Iterator<PairInt> iter) {
        PairInt secondToLastNode = null;
        PairInt lastNode = null;
        while (iter.hasNext()) {
            PairInt n0 = iter.next();
            if (iter.hasNext()) {
                secondToLastNode = n0;
                lastNode = iter.next();
            } else {
                secondToLastNode = lastNode;
                lastNode = n0;
            }
        }
        return new PairInt[]{secondToLastNode, lastNode};
    }

    private int addUUDiagSegmentToRoutes(Routes routes, DiagSegment
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
            int distSq = distSq(segment.p0, r0);
            if (distSq < p0MinDistSq) {
                p0MinDistSq = distSq;
                p0MinIdx = i;
            }
            distSq = distSq(segment.p1, r0);
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
            int distSq = distSq(segment.p2, r1);
            if (distSq < p2MinDistSq) {
                p2MinDistSq = distSq;
                p2MinIdx = i;
            }
            distSq = distSq(segment.p3, r1);
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
        
        return 1;
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

    private ZigZagSegmentRoutes parseDiagZigZag(Segment segment, Set<PairInt> endPoints) {

        if (!(segment instanceof DiagZigZagSegment) &&
            !(segment instanceof DiagZigZagSegment2)) {
            throw new IllegalArgumentException(
            "segment type must be DiagZigZagSegment or DiagZigZagSegment2");
        }

        if (endPoints.size() != 4) {
            throw new IllegalStateException("endPoints must be size 4");
        }

        Set<PairInt> cpEndPoints = new HashSet<PairInt>(endPoints);

        ZigZagSegmentRoutes routes = null;

        /*  first pattern is shown.  The others are consistent w/ different
        locations.
                      E0           -2
              E1end 3  1  -        -1
                    -  0  2  E0end  0
                      E1  -         1

             -3 -2 -1  0  1
        */

        // '1' and '2' on one route and '0' and '3' on another.
        // add the endpoints before and after

        routes = new ZigZagSegmentRoutes();

        routes.ep1 = findClosestTo(segment.p0, cpEndPoints);
        routes.route1.add(routes.ep1);
        routes.route1.add(segment.p0);
        routes.route1.add(segment.p3);
        cpEndPoints.remove(routes.ep1);
        routes.ep1End = findClosestTo(segment.p3, cpEndPoints);
        routes.route1.add(routes.ep1End);

        routes.ep0 = findClosestTo(segment.p1, cpEndPoints);
        routes.route0.add(routes.ep0);
        routes.route0.add(segment.p1);
        routes.route0.add(segment.p2);
        cpEndPoints.remove(routes.ep0);
        routes.ep0End = findClosestTo(segment.p2, cpEndPoints);
        routes.route0.add(routes.ep0End);

        return routes;
    }

    private PairInt findClosestTo(PairInt p0, Set<PairInt> points) {

        PairInt pt = null;
        int minDistSq = Integer.MAX_VALUE;
        for (PairInt p : points) {
            int distSq = distSq(p, p0);
            if (distSq < minDistSq) {
                minDistSq = distSq;
                pt = p;
            }
        }

        return pt;
    }

    private void setNullEndpoints(Routes routes) {

        if (!routes.getRoute0().isEmpty()) {
            if (routes.getEP0() == null) {
                routes.ep0 = routes.getRoute0().iterator().next();
            }
            if (routes.getEP0End() == null) {
                PairInt[] a = getFirstAndLast(routes.getRoute0().iterator());
                routes.ep0End = a[1];
            }
        }

        if (!routes.getRoute1().isEmpty()) {
            if (routes.getEP1() == null) {
                routes.ep1 = routes.getRoute1().iterator().next();
            }
            if (routes.getEP1End() == null) {
                PairInt[] a = getFirstAndLast(routes.getRoute1().iterator());
                routes.ep1End = a[1];
            }
        }
    }

    private void setNullEndpoints(List<Routes> routesList) {

        for (Routes routes : routesList) {
            setNullEndpoints(routes);
        }
    }

    private Segment checkDiagZigZagSegmentPattern(int x, int y, Set<PairInt> points) {

        Pattern pattern = getDiagZigZagSegmentPattern();

        boolean matchesPattern = matchesPattern(x, y, points, pattern);

        /*
        3, 0 are one route and 2, 1 are the other
                -            -2
                3  1  -      -1
                -  0  2       0
                      -       1

         -3 -2 -1  0  1
        */
        if (matchesPattern) {
            DiagZigZagSegment segment = new DiagZigZagSegment();
            segment.p0 = new PairInt(x, y);
            segment.p1 = new PairInt(x, y - 1);
            segment.p2 = new PairInt(x + 1, y);
            segment.p3 = new PairInt(x - 1, y - 1);

            return segment;
        }

        swapXDirection(pattern);

        matchesPattern = matchesPattern(x, y, points, pattern);

        /*
        3, 0 are one route and 2, 1 are the other
                      -      -2
                -  1  3      -1
                2  0  -       0
                -             1

         -3 -2 -1  0  1
        */
        if (matchesPattern) {
            DiagZigZagSegment segment = new DiagZigZagSegment();
            segment.p0 = new PairInt(x, y);
            segment.p1 = new PairInt(x, y - 1);
            segment.p2 = new PairInt(x - 1, y);
            segment.p3 = new PairInt(x + 1, y - 1);

            return segment;
        }

        return null;
    }

    private Segment checkDiagZigZag2SegmentPattern(int x, int y, Set<PairInt> points) {

        Pattern pattern = getDiagZigZag2SegmentPattern();

        boolean matchesPattern = matchesPattern(x, y, points, pattern);

        /*
        3, 0 are one route and 2, 1 are the other
                             -2
                   -  3  -   -1
                   0  1       0
                -  2  -       1

         -3 -2 -1  0  1  2
        */
        if (matchesPattern) {
            DiagZigZagSegment2 segment = new DiagZigZagSegment2();
            segment.p0 = new PairInt(x, y);
            segment.p1 = new PairInt(x + 1, y);
            segment.p2 = new PairInt(x, y + 1);
            segment.p3 = new PairInt(x + 1, y - 1);

            return segment;
        }

        swapYDirection(pattern);

        matchesPattern = matchesPattern(x, y, points, pattern);

        /*
        3, 0 are one route and 2, 1 are the other
                             -2
                -  2  -      -1
                   0  1       0
                   -  3  -    1

         -3 -2 -1  0  1  2
        */
        if (matchesPattern) {
            DiagZigZagSegment2 segment = new DiagZigZagSegment2();
            segment.p0 = new PairInt(x, y);
            segment.p1 = new PairInt(x + 1, y);
            segment.p2 = new PairInt(x, y - 1);
            segment.p3 = new PairInt(x + 1, y + 1);

            return segment;
        }

        return null;
    }

    private void mergeIfAdjacent(List<LinkedList<Segment>> sections) {
        
        if (sections.size() < 2) {
            return;
        }
        
        //TODO: add logic for diagonal merges alone
        
        mergeHorizontalIfAdjacent(sections);
        
        mergeVerticalIfAdjacent(sections);
        
        //TODO: add logic for diagonal in between horizonal or vertical merges
    }
    
    private void mergeHorizontalIfAdjacent(List<LinkedList<Segment>> sections) {
        
        if (sections.size() < 2) {
            return;
        }
        
        // sorts by each list's first item's p0.x
        Collections.sort(sections, new HorizSegmentListComparator());
        
        /*
        using algorithm similar to leftedge.  
                   O(N*lg(N)) + O(N) where N is sections.size()
        
        sort by p0.x, then append, but only if adjacent and a horizontal segment
        
        horizontal s0      [][][]
                   s1    []
                   s2              []
                   s3            []
        */
        
        List<Integer> remove = new ArrayList<Integer>();
        
        int start = 0;
        
        LinkedList<Segment> currentSegments = null;
        while ((currentSegments == null) && (start < sections.size())) {
            LinkedList<Segment> list = sections.get(start);
            assert(!list.isEmpty());
            if (list.get(0) instanceof HorizSegment) {
                currentSegments = list;
            }
            start++;
        }
        
        if (currentSegments == null) {
            return;
        }
        
        for (int i = start; i < sections.size(); ++i) {
            
            LinkedList<Segment> segments = sections.get(i);
            assert(!segments.isEmpty());
            if (!(segments.get(0) instanceof HorizSegment)) {
                // only merging the horizontal segments with this method
                continue;
            }
            
            boolean csOrderedAscendingX = 
                (currentSegments.get(currentSegments.size() - 1).p0.getX() >
                currentSegments.get(0).p0.getX());

            boolean sOrderedAscendingX = 
                (segments.get(segments.size() - 1).p0.getX() >
                segments.get(0).p0.getX());
            
            if (sOrderedAscendingX && 
                ((currentSegments.size() == 1) || csOrderedAscendingX)) {
                
                // last segment of current
                Segment currentLastSegment = currentSegments.get(currentSegments.size() - 1);
            
                if (!(currentLastSegment instanceof HorizSegment)) {
                    continue;
                }
                
                int clP2X = currentLastSegment.p2.getX();
                int clP1X = currentLastSegment.p1.getX();

                int fP3X = segments.get(0).p3.getX();
                int fP0X = segments.get(0).p0.getX();
            
                /*
                immediately adjacent:
                   clP2X  P3
                   clP1X  P0
                */
                if ((((clP2X + 1) == fP3X) && ((clP1X + 1) == fP0X)) ||
                    ((clP2X == fP3X) && (clP1X == fP0X)) ) {
                    currentSegments.addAll(segments);
                    remove.add(Integer.valueOf(i));
                }
                
            } else if (!sOrderedAscendingX && 
                ((currentSegments.size() == 1) || !csOrderedAscendingX)) {
                
                /*
                     c[1]    c[0]    s[n-1]
                     3  2    3  2     3  2
                     0  1    0  1     0  1
                */
                Segment currentFirstSegment = currentSegments.get(0);
                if (!(currentFirstSegment instanceof HorizSegment)) {
                    continue;
                }
            
                Segment lastSegment = segments.get(segments.size() - 1);
                if (!(lastSegment instanceof HorizSegment)) {
                    continue;
                }
                
                int cP2X = currentFirstSegment.p2.getX();
                int cP1X = currentFirstSegment.p1.getX();

                int sP3X = lastSegment.p3.getX();
                int sP0X = lastSegment.p0.getX();
                
                if ((((sP3X + 1) == cP2X) && ((sP0X + 1) == cP1X)) 
                    || ((sP3X == cP2X) && (sP0X == cP1X))) {
                    LinkedList<Segment> tmp = new LinkedList<Segment>();
                    tmp.addAll(segments);
                    tmp.addAll(currentSegments);
                    remove.add(Integer.valueOf(i));
                    currentSegments.clear();
                    currentSegments.addAll(tmp);
                }
            } else {
                //TODO: handle merging when segments were assembled with different
                // directions.  preferably, do this at top of method.
                // the caveat is that a linked list may be composed of items
                // which are composed of any combination of horizontal,
                // vertical, and diagonal already ordered by curve point order.
            }
        }
        
        for (int i = (remove.size() - 1); i > -1; --i) {
            sections.remove(remove.get(i).intValue());
        }
    }
    
    private void mergeVerticalIfAdjacent(List<LinkedList<Segment>> sections) {
        
        if (sections.size() < 2) {
            return;
        }
        
        // sorts by each list's first item's p0.y, descending
        Collections.sort(sections, new VertSegmentListComparator());
        
        /*
        using algorithm similar to leftedge.  
                   O(N*lg(N)) + O(N) where N is sections.size()
        
        sort by p0.y, then append, but only if adjacent and a vertical segment
        */
        
        List<Integer> remove = new ArrayList<Integer>();
        
        int start = 0;
        
        LinkedList<Segment> currentSegments = null;
        while ((currentSegments == null) && (start < sections.size())) {
            LinkedList<Segment> list = sections.get(start);
            assert(!list.isEmpty());
            if (list.get(0) instanceof VertSegment) {
                currentSegments = list;
            }
            start++;
        }
        
        if (currentSegments == null) {
            return;
        }
        
        for (int i = start; i < sections.size(); ++i) {
            
            LinkedList<Segment> segments = sections.get(i);
            assert(!segments.isEmpty());
            if (!(segments.get(0) instanceof VertSegment)) {
                // only merging the vertical segments with this method
                continue;
            }
            
            boolean csOrderedDescendingY = 
                (currentSegments.get(currentSegments.size() - 1).p0.getY() <
                currentSegments.get(0).p0.getY());

            boolean sOrderedDescendingY = 
                (segments.get(segments.size() - 1).p0.getY() <
                segments.get(0).p0.getY());
            
            if (sOrderedDescendingY && 
                ((currentSegments.size() == 1) || csOrderedDescendingY)) {
                
                // last segment of current
                Segment currentLastSegment = currentSegments.get(currentSegments.size() - 1);
                
                if (!(currentLastSegment instanceof VertSegment)) {
                    continue;
                }
            
                /*   
                                   2 1
                             -3    3 0 s[0]
                    2 1      -2
                    3 0 c[1] -1
                    2 1       0
                    3 0 c[0]
                */
                
                int clP1Y = currentLastSegment.p1.getY();
                int clP2Y = currentLastSegment.p2.getY();
                
                int fP0Y = segments.get(0).p0.getY();
                int fP3Y = segments.get(0).p3.getY();
                
                if ((((clP1Y - 1) == fP0Y) && ((clP2Y - 1) == fP3Y))
                    || ((clP1Y == fP0Y) && (clP2Y == fP3Y))) {
                    currentSegments.addAll(segments);
                    remove.add(Integer.valueOf(i));
                    continue;
                }
                
                /*   
                             -3 
                    2 1      -2
                    3 0 c[1] -1
                    2 1       0
                    3 0 c[0]
                                   2 1
                                   3 0  s[n-1]
                */
                Segment currentFirstSegment = currentSegments.get(0);
                if (!(currentFirstSegment instanceof VertSegment)) {
                    continue;
                }
            
                Segment lastSegment = segments.get(segments.size() - 1);
                if (!(lastSegment instanceof VertSegment)) {
                    continue;
                }
                
                int cP0Y = currentFirstSegment.p0.getY();
                int cP3Y = currentFirstSegment.p3.getY();

                int sP1Y = lastSegment.p1.getY();
                int sP2Y = lastSegment.p2.getY();
                
                if ((((sP2Y - 1) == cP3Y) && ((sP1Y - 1) == cP0Y)) 
                    || ((sP2Y == cP3Y) && (sP1Y == cP0Y))
                    ) {
                    LinkedList<Segment> tmp = new LinkedList<Segment>();
                    tmp.addAll(segments);
                    tmp.addAll(currentSegments);
                    remove.add(Integer.valueOf(i));
                    currentSegments.clear();
                    currentSegments.addAll(tmp);
                }
                
            } else if (!sOrderedDescendingY && 
                ((currentSegments.size() == 1) || !csOrderedDescendingY)) {
                /*   
                                   2 1
                             -3    3 0 s[n-1]
                    2 1      -2
                    3 0 c[0] -1
                    2 1       0
                    3 0 c[1]
                */
                Segment currentFirstSegment = currentSegments.get(0);
                if (!(currentFirstSegment instanceof VertSegment)) {
                    continue;
                }
            
                Segment lastSegment = segments.get(segments.size() - 1);
                if (!(lastSegment instanceof VertSegment)) {
                    continue;
                }
                
                int cP2Y = currentFirstSegment.p2.getY();
                int cP1Y = currentFirstSegment.p1.getY();

                int sP3Y = lastSegment.p3.getY();
                int sP0Y = lastSegment.p0.getY();
                
                if ((((cP2Y - 1) == sP3Y) && ((cP1Y - 1) == sP0Y))
                    || ((cP2Y == sP3Y) && (cP1Y == sP0Y))){
                    LinkedList<Segment> tmp = new LinkedList<Segment>();
                    tmp.addAll(segments);
                    tmp.addAll(currentSegments);
                    remove.add(Integer.valueOf(i));
                    currentSegments.clear();
                    currentSegments.addAll(tmp);
                    continue;
                }
                
                /*   
                             -3 
                    2 1      -2
                    3 0 c[0] -1
                    2 1       0
                    3 0 c[1]
                                   2 1
                                   3 0  s[0]
                */
                
                // last segment of current
                Segment currentLastSegment = currentSegments.get(currentSegments.size() - 1);
                if (!(currentLastSegment instanceof VertSegment)) {
                    continue;
                }
                
                Segment firstSegment = segments.get(0);
                if (!(firstSegment instanceof VertSegment)) {
                    continue;
                }
                
                int clP3Y = currentLastSegment.p3.getY();
                int clP0Y = currentLastSegment.p0.getY();

                int sP2Y = firstSegment.p2.getY();
                int sP1Y = firstSegment.p1.getY();
                
                if ((((sP2Y - 1) == clP3Y) && ((sP1Y - 1) == clP0Y))
                    || ((sP2Y == clP3Y) && (sP1Y == clP0Y))) {
                    currentSegments.addAll(segments);
                    remove.add(Integer.valueOf(i));
                }
                
            } else {
                //TODO: handle merging when segments were assembled with different
                // directions.  preferably, do this at top of method.
                // the caveat is that a linked list may be composed of items
                // which are composed of any combination of horizontal,
                // vertical, and diagonal already ordered by curve point order.
            }
        }
        
        for (int i = (remove.size() - 1); i > -1; --i) {
            sections.remove(remove.get(i).intValue());
        }
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
    public static class DiagZigZagSegment extends Segment {
    }
    public static class DiagZigZagSegment2 extends Segment {
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

    protected Pattern getDiagZigZagSegmentPattern() {

        /*
        3, 0 are one route and 2, 1 are the other
                -            -2
                3  1  -      -1
                -  0  2       0
                      -       1

         -3 -2 -1  0  1

        Each of the 4 needs at least one neighbor that is not one of the 4
        points in the zig zap and all of their neighbors cannot be adjacent
        to any of the other neighbors.
        */
        Pattern pr = new Pattern();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, 0)); pr.zeroes.add(new PairInt(-1, -2));
        pr.zeroes.add(new PairInt(1, 1)); pr.zeroes.add(new PairInt(1, -1));

        pr.ones.add(new PairInt(-1, -1));
        pr.ones.add(new PairInt(0, -1));
        pr.ones.add(new PairInt(1, 0));

        return pr;
    }

    protected Pattern getDiagZigZag2SegmentPattern() {

        /*
        3, 0 are one route and 2, 1 are the other
                             -2
                   -  3  -   -1
                   0  1       0
                -  2  -       1

         -3 -2 -1  0  1  2

        Each of the 4 needs at least one neighbor that is not one of the 4
        points in the zig zap and all of their neighbors cannot be adjacent
        to any of the other neighbors.
        */
        Pattern pr = new Pattern();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, 1));
        pr.zeroes.add(new PairInt(0, -1));
        pr.zeroes.add(new PairInt(1, 1));
        pr.zeroes.add(new PairInt(2, -1));

        pr.ones.add(new PairInt(0, 1));
        pr.ones.add(new PairInt(1, 0));  pr.ones.add(new PairInt(1, -1));

        return pr;
    }

    private Segment checkVertHorizSegmentPattern(int x, int y,
        Set<PairInt> neighbors, boolean useVertical) {

        /*
            VertSegment    
             R1  R0
             /|\ |
              |  |
              | \|/
              .  .      -2
           -  2  1  -   -1  
           -  3  0  -    0  
              .  .       1 
          -2 -1  0  1     

               HorizPattern 
                    -  -  -   -2
            R1-->   .  3  2   -1  -->R1 ends
                    .  0  1    0
            R0<--   -  -  -    1  <--R0 starts
                   -1  0  1   
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
    
    protected int distSq(PairInt p0, PairInt p1) {
        int diffX = p0.getX() - p1.getX();
        int diffY = p0.getY() - p1.getY();
        int distSq = (diffX * diffX) + (diffY * diffY);
        return distSq;
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
        public PairInt getEP0() {
            return ep0;
        }
        public PairInt getEP0End() {
            return ep0End;
        }
        public PairInt getEP1() {
            return ep1;
        }
        public PairInt getEP1End() {
            return ep1End;
        }
        public void applyOffsets(final int xOffset, final int yOffset) {
            if (ep0 != null) {
                ep0 = new PairInt(ep0.getX() + xOffset, ep0.getY() + yOffset);
            }
            if (ep0End != null) {
                ep0End = new PairInt(ep0End.getX() + xOffset, ep0End.getY() + yOffset);
            }
            if (ep1 != null) {
                ep1 = new PairInt(ep1.getX() + xOffset, ep1.getY() + yOffset);
            }
            if (ep1End != null) {
                ep1End = new PairInt(ep1End.getX() + xOffset, ep1End.getY() + yOffset);
            }

            LinkedHashSet<PairInt> tmp = new LinkedHashSet<PairInt>();
            Iterator<PairInt> iter = this.route0.iterator();
            while (iter.hasNext()) {
                PairInt p = iter.next();
                tmp.add(new PairInt(p.getX() + xOffset, p.getY() + yOffset));
            }
            route0 = tmp;

            tmp = new LinkedHashSet<PairInt>();
            iter = this.route1.iterator();
            while (iter.hasNext()) {
                PairInt p = iter.next();
                tmp.add(new PairInt(p.getX() + xOffset, p.getY() + yOffset));
            }
            route1 = tmp;
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

    private static class HorizSegmentListComparator implements Comparator<LinkedList<Segment>> {

        public HorizSegmentListComparator() {
        }

        @Override
        public int compare(LinkedList<Segment> o1, LinkedList<Segment> o2) {
            
            if (o1.isEmpty() && o2.isEmpty()) {
                return 0;
            } else if (o1.isEmpty() && !o2.isEmpty()) {
                return 1;
            } else if (!o1.isEmpty() && o2.isEmpty()) {
                return -1;
            }
            
            return Integer.compare(o1.get(0).p0.getX(), o2.get(0).p0.getX());
        }
    }

    private static class VertSegmentListComparator implements Comparator<LinkedList<Segment>> {

        public VertSegmentListComparator() {
        }

        @Override
        public int compare(LinkedList<Segment> o1, LinkedList<Segment> o2) {
            
            if (o1.isEmpty() && o2.isEmpty()) {
                return 0;
            } else if (o1.isEmpty() && !o2.isEmpty()) {
                return 1;
            } else if (!o1.isEmpty() && o2.isEmpty()) {
                return -1;
            }
            
            return Integer.compare(o2.get(0).p0.getY(), o1.get(0).p0.getY());
        }
    }
}
