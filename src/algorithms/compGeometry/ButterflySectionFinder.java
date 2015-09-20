package algorithms.compGeometry;

import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.HashSet;
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
    public List<Set<PairInt>> findButterflySections(PairIntArray closedCurve) {
        
        Set<PairInt> points = Misc.convert(closedCurve);
        
        List<Set<PairInt>> sections = findButterflySectionsLarge(closedCurve, 
            points);

        List<Set<PairInt>> sectionsSmall = findButterflySectionsSmall(
            closedCurve, points);
        
        if (sections == null && sectionsSmall == null) {
            return null;
        } else if (sections == null && sectionsSmall != null) {
            return sectionsSmall;
        } else if (sections != null && sectionsSmall == null) {
            return sections;
        } else if (sections != null && sectionsSmall != null) {
            sections.addAll(sectionsSmall);
        }
        
        return sections;
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
    protected List<Set<PairInt>> findButterflySectionsLarge(PairIntArray 
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

        Set<PairInt> exclude = new HashSet<PairInt>();
        for (LinkedList<Segment> section : candidateSections) {
            for (Segment segment : section) {
                exclude.add(segment.p0);
                exclude.add(segment.p1);
                exclude.add(segment.p2);
                exclude.add(segment.p3);
            }
        }

        List<Set<PairInt>> output = new ArrayList<Set<PairInt>>();

        // -- scan for endpoints --
        for (LinkedList<Segment> section : candidateSections) {

            boolean checkFirstSegment = true;
            
            Set<PairInt> endPoints = checkForAdjacentEndpoints(points, section,
                exclude, checkFirstSegment);

            if (endPoints == null || endPoints.isEmpty()) {
                continue;
            }

            Set<PairInt> exclude2 = new HashSet<PairInt>(exclude);
            exclude2.addAll(endPoints);

            Set<PairInt> endPoints2 = checkForAdjacentEndpoints(points,
                section, exclude2, checkFirstSegment);

            if (endPoints2 == null || endPoints2.isEmpty()) {
                continue;
            }

            // add all section and endpoints to a set to add to output
            Set<PairInt> allPoints = new HashSet<PairInt>();
            allPoints.addAll(endPoints);
            allPoints.addAll(endPoints2);
            for (Segment segment : section) {
                allPoints.add(segment.p0);
                allPoints.add(segment.p1);
                allPoints.add(segment.p2);
                allPoints.add(segment.p3);
            }

            output.add(allPoints);
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
    
    protected List<Set<PairInt>> findButterflySectionsSmall(PairIntArray 
        closedCurve, Set<PairInt> points) {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        List<Set<PairInt>> output = new ArrayList<Set<PairInt>>();    

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
                continue;
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
            
            // add all section and endpoints to a set to add to output
            Set<PairInt> allPoints = new HashSet<PairInt>();
            allPoints.addAll(endPoints);
            allPoints.add(segment.p0);
            allPoints.add(segment.p1);
            allPoints.add(segment.p2);
            allPoints.add(segment.p3);

            output.add(allPoints);
        }

        return output;
    }

    private Segment checkSegmentPatterns(final int x, final int y, 
        Set<PairInt> neighbors) {

        Segment segment = checkVertSegmentPattern(x, y, neighbors);
        
        if (segment != null) {
            return segment;
        }
        
        segment = checkHorizSegmentPattern(x, y, neighbors);
        
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

    private Set<PairInt> checkForAdjacentEndpoints(Set<PairInt> points,
        LinkedList<Segment> section, Set<PairInt> exclude, 
        boolean checkFirstSegment) {
        
        Segment segment = checkFirstSegment ? section.getFirst() : 
            section.getLast();
                        
        if (segment instanceof VertSegment) {
            return checkForAdjacentEndpointsForVertSegment(points, section, 
                exclude, checkFirstSegment);
        } else if (segment instanceof HorizSegment) {
            return checkForAdjacentEndpointsForHorizSegment(points, section, 
                exclude, checkFirstSegment);
        } else if (segment instanceof UUDiagSegment) {
            return checkForAdjacentEndpointsForUUDiagSegment(points, section, 
                exclude, checkFirstSegment);
        }
        
        return null;
        
        /*endpoints for horiz:
                       #           #
            # .        - .       - - .
          - - .  or  - - .  or     # .
            #          #

        endpoints for vert:
                -        - -       -
              # - #    # - - #   # - #
              . .        . .       . .

        endpoints for diag:

            -  #             -  #               -  #
            -    .        -  -  .            #  -  .
            #  .   .      #  .    .          -  .    .
                 .              .                  .
        */
        
    }

    private Set<PairInt> checkForAdjacentEndpointsForVertSegment(
        Set<PairInt> points, LinkedList<Segment> section, Set<PairInt> exclude, 
        boolean checkFirstSegment) {
        
        VertSegment segment = checkFirstSegment ? (VertSegment)section.getFirst() :
            (VertSegment)section.getLast();
                
        /*endpoints for vert:
                -        - -       -
              # - #    # - - #   # - #
              . .        . .       . .
        */
       
        Set<PairInt> endPoints = findEndPointsVertPatterns(points, segment,
            exclude);
        
        return endPoints;
    }

    private Set<PairInt> findEndPointsVertPatterns(Set<PairInt> points, 
        VertSegment segment, Set<PairInt> exclude) {
        
        int x0 = segment.p0.getX();
        int y0 = segment.p0.getY();
        
        Pattern[] patterns = new Pattern[] {
            getEndPointsVertPattern1(), getEndPointsVertPattern1Opp(),
            getEndPointsVertPattern2(), getEndPointsVertPattern2Opp(),
            getEndPointsVertPattern3(), getEndPointsVertPattern3Opp()
        };
        
        for (Pattern pattern : patterns) {
            boolean found = true;
            for (PairInt p : pattern.zeroes) {
                PairInt p2 = new PairInt(x0 + p.getX(), y0 + p.getY());
                if (points.contains(p2) || exclude.contains(p2)) {
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
                    if (exclude.contains(p2) || !points.contains(p2)) {
                        found = false;
                        break;
                    }
                    endPoints.add(p2);
                }
                if (found) {
                    return endPoints;
                }
            }
        }        
        return null;
    }
    
    private Pattern getEndPointsVertPattern1() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
                 -      -4
              #  -  #   -3
           -  .  .  -   -2
           -  2  1  -   -1
           -  3  0  -    0
                        
          -2 -1  0  1
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
        
        pattern.ones.add(new PairInt(0, -2)); pattern.ones.add(new PairInt(-1, -2));
        pattern.ones.add(new PairInt(1, -3)); pattern.ones.add(new PairInt(-1, -3));
        
        pattern.zeroes.add(new PairInt(-2, -2)); pattern.zeroes.add(new PairInt(1, -2));
        pattern.zeroes.add(new PairInt(0, -3)); pattern.zeroes.add(new PairInt(0, -4));
        return pattern;
    }
    
    private Pattern getEndPointsVertPattern1Opp() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
                            
           -  2  1  -    1
           -  3  0  -    0
           -  .  .  -   -1
              #  -  #   -2
                 -      -3
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
        
              -  -      -4
           #  -  -  #   -3
           -  .  .  -   -2
           -  2  1  -   -1
           -  3  0  -    0
                        
          -2 -1  0  1
        */
        
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
        
        pattern.ones.add(new PairInt(-1, -2)); pattern.ones.add(new PairInt(0, -2));
        pattern.ones.add(new PairInt(-2, -3)); pattern.ones.add(new PairInt(1, -3));
        
        pattern.zeroes.add(new PairInt(-2, -2)); pattern.zeroes.add(new PairInt(1, -2));
        pattern.zeroes.add(new PairInt(-1, -3)); pattern.zeroes.add(new PairInt(-1, -4));
        pattern.zeroes.add(new PairInt(0, -3)); pattern.zeroes.add(new PairInt(0, -4));
        
        return pattern;
    }
    
    private Pattern getEndPointsVertPattern2Opp() {
        
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
           -  2  1  -    1
           -  3  0  -    0
           -  .  .  -   -1
           #  -  -  #   -2
              -  -      -3
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
                
                 -      -4
              #  -  #   -3
           -  .  .  -   -2
           -  2  1  -   -1
           -  3  0  -    0
                        
          -2 -1  0  1
       
        */
        
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
        
        pattern.ones.add(new PairInt(-1, -2)); pattern.ones.add(new PairInt(1, -3));
        pattern.ones.add(new PairInt(0, -2));
        pattern.ones.add(new PairInt(1, -3)); 
        
        pattern.zeroes.add(new PairInt(-2, -2)); 
        pattern.zeroes.add(new PairInt(0, -3)); pattern.zeroes.add(new PairInt(0, -4));
        pattern.zeroes.add(new PairInt(1, -2)); 
        
        return pattern;
    }
    
    private Pattern getEndPointsVertPattern3Opp() {
        
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
                                 
           -  2  1  -    1
           -  3  0  -    0
           -  .  .  -   -1
              #  -  #   -2
                 -      -3
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
    
    private Set<PairInt> checkForAdjacentEndpointsForHorizSegment(
        Set<PairInt> points, LinkedList<Segment> section, Set<PairInt> exclude, 
        boolean checkFirstSegment) {
        
        HorizSegment segment = checkFirstSegment ? (HorizSegment)section.getFirst() :
            (HorizSegment)section.getLast();
                
        Set<PairInt> endPoints = findEndPointsHorizPatterns(points, segment,
            exclude);
        
        return endPoints;
    }
    
    private Set<PairInt> findEndPointsHorizPatterns(Set<PairInt> points, 
        HorizSegment segment, Set<PairInt> exclude) {
        
        int x0 = segment.p0.getX();
        int y0 = segment.p0.getY();
        
        Pattern[] patterns = new Pattern[] {
            getEndPointsHorizPattern1(), getEndPointsHorizPattern1Opp(),
            getEndPointsHorizPattern2(), getEndPointsHorizPattern2Opp(),
            getEndPointsHorizPattern3(), getEndPointsHorizPattern3Opp()
        };
        
        for (Pattern pattern : patterns) {
            boolean found = true;
            for (PairInt p : pattern.zeroes) {
                PairInt p2 = new PairInt(x0 + p.getX(), y0 + p.getY());
                if (points.contains(p2) || exclude.contains(p2)) {
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
                    if (exclude.contains(p2) || !points.contains(p2)) {
                        found = false;
                        break;
                    }
                    endPoints.add(p2);
                }
                if (found) {
                    return endPoints;
                }
            }
        }        
        return null;
    }
    
    private Pattern getEndPointsHorizPattern1() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
                 -  -  -       -2
              .  2  1  .  #    -1
              .  3  0  .  -  -  0
                 -  -  -  #     1
             -2 -1  0  1  2  3
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
    
        pattern.ones.add(new PairInt(1, 0)); pattern.ones.add(new PairInt(1, -1));
        pattern.ones.add(new PairInt(2, 1)); pattern.ones.add(new PairInt(2, -1));
        
        pattern.zeroes.add(new PairInt(1, 1)); pattern.zeroes.add(new PairInt(1, -2));
        pattern.zeroes.add(new PairInt(2, 0)); 
        pattern.zeroes.add(new PairInt(3, 0));
        
        return pattern;
    }
    
    private Pattern getEndPointsHorizPattern1Opp() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
                 -  -  -       -2
              .  2  1  .  #    -1
              .  3  0  .  -  -  0
                 -  -  -  #     1
              2  1  0 -1 -2 -3
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
     
        pattern.ones.add(new PairInt(-1, 0)); pattern.ones.add(new PairInt(-1, -1));
        pattern.ones.add(new PairInt(-2, 1)); pattern.ones.add(new PairInt(-2, -1));
        
        pattern.zeroes.add(new PairInt(-1, 1)); pattern.zeroes.add(new PairInt(-1, -2));
        pattern.zeroes.add(new PairInt(-2, 0)); pattern.zeroes.add(new PairInt(-3, 0));
        
        return pattern;
    }
   
    private Pattern getEndPointsHorizPattern2() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
                 -  -  -  #    -2
              .  2  1  .  -    -1
              .  3  0  .  -  -  0
                 -  -  -  #     1
             -2 -1  0  1  2  3
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
     
        pattern.ones.add(new PairInt(1, 0)); pattern.ones.add(new PairInt(1, -1));
        pattern.ones.add(new PairInt(2, 1)); pattern.ones.add(new PairInt(2, -2));
        
        pattern.zeroes.add(new PairInt(1, 1)); pattern.zeroes.add(new PairInt(1, -2));
        pattern.zeroes.add(new PairInt(2, 0)); pattern.zeroes.add(new PairInt(2, -1));
        pattern.zeroes.add(new PairInt(3, 0));
        
        return pattern;
    }
     
    private Pattern getEndPointsHorizPattern2Opp() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
                 -  -  -  #    -2
              .  2  1  .  -    -1
              .  3  0  .  -  -  0
                 -  -  -  #     1
              2  1  0 -1 -2 -3
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
    
        pattern.ones.add(new PairInt(-1, 0)); pattern.ones.add(new PairInt(-1, -1));
        pattern.ones.add(new PairInt(-2, 1)); pattern.ones.add(new PairInt(-2, -2));
        
        pattern.zeroes.add(new PairInt(-1, 1)); pattern.zeroes.add(new PairInt(-1, -2));
        pattern.zeroes.add(new PairInt(-2, 0)); pattern.zeroes.add(new PairInt(-2, -1));
        pattern.zeroes.add(new PairInt(-3, 0));
        
        return pattern;
    }
 
    private Pattern getEndPointsHorizPattern3() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
                 -  -  -  #    -2
              .  2  1  .  -  - -1
              .  3  0  .  #     0
                 -  -  -        1
             -2 -1  0  1  2  3
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
    
        pattern.ones.add(new PairInt(1, 0)); pattern.ones.add(new PairInt(1, -1));
        pattern.ones.add(new PairInt(2, 0)); pattern.ones.add(new PairInt(2, -2));
        
        pattern.zeroes.add(new PairInt(1, 1)); pattern.zeroes.add(new PairInt(1, -2));
        pattern.zeroes.add(new PairInt(2, -1)); pattern.zeroes.add(new PairInt(3, -1));
        
        return pattern;
    }
    
    private Pattern getEndPointsHorizPattern3Opp() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
                 -  -  -  #     -2
              .  2  1  .  -  -  -1
              .  3  0  .  #      0
                 -  -  -         1
              2  1  0 -1 -2 -3
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
    
        pattern.ones.add(new PairInt(-1, 0)); pattern.ones.add(new PairInt(-1, -1));
        pattern.ones.add(new PairInt(-2, 0)); pattern.ones.add(new PairInt(-2, -2));
        
        pattern.zeroes.add(new PairInt(-1, 1)); pattern.zeroes.add(new PairInt(-1, -2));
        pattern.zeroes.add(new PairInt(-2, -1)); pattern.zeroes.add(new PairInt(-3, -1));
        
        return pattern;
    }
    
    private Set<PairInt> checkForAdjacentEndpointsForUUDiagSegment(
        Set<PairInt> points, LinkedList<Segment> section, Set<PairInt> exclude, 
        boolean checkFirstSegment) {
        
        UUDiagSegment segment = checkFirstSegment ? (UUDiagSegment)section.getFirst() :
            (UUDiagSegment)section.getLast();
                
        Set<PairInt> endPoints = findEndPointsUUDiagPatterns(points, segment,
            exclude);
        
        return endPoints;
    }
    
    private Set<PairInt> findEndPointsUUDiagPatterns(Set<PairInt> points, 
        UUDiagSegment segment, Set<PairInt> exclude) {
        
        int x0 = segment.p0.getX();
        int y0 = segment.p0.getY();
        
        Pattern[] patterns = new Pattern[] {
            getEndPointsUUDiagPattern1(), getEndPointsULDiagPattern1(),
            getEndPointsUUDiagPattern2(), getEndPointsULDiagPattern2(),
            getEndPointsUUDiagPattern3(), getEndPointsULDiagPattern3()
        };
        
        for (Pattern pattern : patterns) {
            boolean found = true;
            for (PairInt p : pattern.zeroes) {
                PairInt p2 = new PairInt(x0 + p.getX(), y0 + p.getY());
                if (points.contains(p2) || exclude.contains(p2)) {
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
                    if (exclude.contains(p2) || !points.contains(p2)) {
                        found = false;
                        break;
                    }
                    endPoints.add(p2);
                }
                if (found) {
                    return endPoints;
                }
            }
        }        
        return null;
    }
    
    private Pattern getEndPointsUUDiagPattern1() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
                -  -          -2
          -  2  1  -  -       -1
          -  -  3 .0  #        0
             -  .  -           1
                -  #           2
        
         -3 -2 -1  0  1
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
    
        pattern.ones.add(new PairInt(-1, 1)); 
        pattern.ones.add(new PairInt(0, 2));
        pattern.ones.add(new PairInt(1, 0));
        
        pattern.zeroes.add(new PairInt(-1, 1));
        pattern.zeroes.add(new PairInt(0, 1)); 
        pattern.zeroes.add(new PairInt(0, -1));
        
        return pattern;
    }
    
    private Pattern getEndPointsULDiagPattern1() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.        
        
                -  -          -2
          -  2  1  -  -       -1
          -  -  3 .0  #        0
             -  .  -  -        1
                -  #           2
        
          3  2  1  0  -1
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();

        pattern.ones.add(new PairInt(1, 1)); 
        pattern.ones.add(new PairInt(0, 2));
        pattern.ones.add(new PairInt(-1, 0));
        
        pattern.zeroes.add(new PairInt(1, 2)); 
        pattern.zeroes.add(new PairInt(0, 1)); pattern.zeroes.add(new PairInt(0, -1));
        pattern.zeroes.add(new PairInt(-1, 1)); 
        
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

        pattern.ones.add(new PairInt(-1, 2)); pattern.ones.add(new PairInt(-1, 1));
        
        pattern.zeroes.add(new PairInt(0, 2)); pattern.zeroes.add(new PairInt(0, 1)); pattern.zeroes.add(new PairInt(0, 1));
        
        return pattern;
    }
    
    private Pattern getEndPointsULDiagPattern2() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
                 -  -         -2
          -  2  1  -  -       -1
          -  -  3 .0  #        0
             -  .  -  -        1
                #  -           2
        
          3  2  1  0 -1
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
   
        pattern.ones.add(new PairInt(1, 2)); pattern.ones.add(new PairInt(1, 1));
        pattern.ones.add(new PairInt(-1, 0));
        
        pattern.zeroes.add(new PairInt(0, 2)); pattern.zeroes.add(new PairInt(0, 1)); pattern.zeroes.add(new PairInt(0, -1));
        pattern.zeroes.add(new PairInt(-1, 1)); pattern.zeroes.add(new PairInt(-1, -1)); 
        
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
   
        pattern.ones.add(new PairInt(-1, 2)); pattern.ones.add(new PairInt(-1, 1));
        pattern.ones.add(new PairInt(1, 1));
       
        pattern.zeroes.add(new PairInt(0, 2)); pattern.zeroes.add(new PairInt(0, 1)); pattern.zeroes.add(new PairInt(0, -1));
        pattern.zeroes.add(new PairInt(1, 0)); 
        
        return pattern;
    }
    
    private Pattern getEndPointsULDiagPattern3() {
        /* the pattern returned is relative to 
        position '0', just like the other patterns.
        
                -  -          -2
          -  2  1  -          -1
          -  -  3 .0  -        0
             -  .  -  #        1
                #  -           2
        
          3  2  1  0 -1
        */
        Pattern pattern = new Pattern();
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
    
        pattern.ones.add(new PairInt(1, 2)); pattern.ones.add(new PairInt(1, 1));
        pattern.ones.add(new PairInt(-1, 1));
        
        pattern.zeroes.add(new PairInt(0, 2)); pattern.zeroes.add(new PairInt(0, 1)); pattern.zeroes.add(new PairInt(0, -1));
        pattern.zeroes.add(new PairInt(-1, 0));
        
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
    }
    
    public static class VertSegment extends Segment {
        /*    .  .      -2
           -  2  1  -   -1
           -  3  0  -    0
              .  .       1
          -2 -1  0  1
        */
    }
    public static class HorizSegment extends Segment {
        /*    -  -      -2
           .  2  1  .   -1
           .  3  0  .    0
              -  -       1
          -2 -1  0  1
        */
    }
    public static class UUDiagSegment extends Segment {
        /*      -  -      -2
          -  2  1  -  -   -1
          -  -  3  0  -    0
             -  -          1
         -3 -2 -1  0  1
        */
    }
    public static class ULDiagSegment extends Segment {
        /*      -  -      -2
          -  2  1  -  -   -1
          -  -  3  0  -    0
             -  -          1
          3  2  1  0 -1
        */
    }
    public static class ZigZagSegment extends Segment {
        /*
                   2  -      -2
                   -  1      -1
                   0  -       0
                   -  3       1

         -3 -2 -1  0  1
        */
    }

    public static class Pattern {
        Set<PairInt> ones;
        Set<PairInt> zeroes;
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
    
    protected Pattern getHorizSegmentPattern() {

        /*    -  -      -2
           .  2  1  .   -1
           .  3  0  .    0
              -  -       1
          -2 -1  0  1
        */
        Pattern pr = new Pattern();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, 1)); pr.zeroes.add(new PairInt(-1, -2));
        pr.zeroes.add(new PairInt(0, 1)); pr.zeroes.add(new PairInt(0, -2));

        pr.ones.add(new PairInt(-1, 0)); pr.ones.add(new PairInt(-1, -1));
        pr.ones.add(new PairInt(1, 0)); pr.ones.add(new PairInt(1, -1));

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
    
    private Segment checkVertSegmentPattern(int x, int y, Set<PairInt> neighbors) {

        /*    .  .      -2
           -  2  1  -   -1
           -  3  0  -    0
              .  .       1
          -2 -1  0  1
        */
        
        Pattern pattern = getVertSegmentPattern();
        
        boolean matchesPattern = matchesPattern(x, y, neighbors, pattern);
        
        if (matchesPattern) {
            VertSegment segment = new VertSegment();
            segment.p0 = new PairInt(x, y);
            segment.p1 = new PairInt(x, y - 1);
            segment.p2 = new PairInt(x - 1, y - 1);
            segment.p3 = new PairInt(x - 1, y);
            
            return segment;
        }
        
        swapYDirection(pattern);
        
        matchesPattern = matchesPattern(x, y, neighbors, pattern);
        
        /*    .  .       2
           -  2  1  -    1
           -  3  0  -    0
              .  .      -1
          -2 -1  0  1
        */
        if (matchesPattern) {
            VertSegment segment = new VertSegment();
            segment.p0 = new PairInt(x, y);
            segment.p1 = new PairInt(x, y + 1);
            segment.p2 = new PairInt(x - 1, y + 1);
            segment.p3 = new PairInt(x - 1, y);
            return segment;
        }
        
        return null;
    }
    
    private Segment checkHorizSegmentPattern(int x, int y, Set<PairInt> neighbors) {

        /*    -  -      -2
           .  2  1  .   -1
           .  3  0  .    0
              -  -       1
          -2 -1  0  1
        */
        
        Pattern pattern = getHorizSegmentPattern();
        
        boolean matchesPattern = matchesPattern(x, y, neighbors, pattern);
        
        if (matchesPattern) {
            HorizSegment segment = new HorizSegment();
            segment.p0 = new PairInt(x, y);
            segment.p1 = new PairInt(x, y - 1);
            segment.p2 = new PairInt(x - 1, y - 1);
            segment.p3 = new PairInt(x - 1, y);
            
            return segment;
        }
        
        swapXDirection(pattern);
        
        matchesPattern = matchesPattern(x, y, neighbors, pattern);
        
        /*    -  -      -2
           .  2  1  .   -1
           .  3  0  .    0
              -  -       1
           2  1  0 -1
        */
        if (matchesPattern) {
            HorizSegment segment = new HorizSegment();
            segment.p0 = new PairInt(x, y);
            segment.p1 = new PairInt(x, y - 1);
            segment.p2 = new PairInt(x + 1, y - 1);
            segment.p3 = new PairInt(x + 1, y);
            
            return segment;
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
        
        swapXDirection(pattern);
        
        matchesPattern = matchesPattern(x, y, neighbors, pattern);
        
        /*      -  -      -2
          -  2  1  -  -   -1
          -  -  3  0  -    0
             -  -          1
          3  2  1  0 -1
        */
        if (matchesPattern) {
            ULDiagSegment segment = new ULDiagSegment();
            segment.p0 = new PairInt(x, y);
            segment.p1 = new PairInt(x + 1, y - 1);
            segment.p2 = new PairInt(x + 2, y - 1);
            segment.p3 = new PairInt(x + 1, y);
            
            return segment;
        }
        
        return null;
    }
    
    private Segment checkZigZagSegmentPattern(int x, int y, Set<PairInt> points) {

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

}
