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

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        Set<PairInt> points = Misc.convert(closedCurve);

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

        List<Set<PairInt>> output = new ArrayList<Set<PairInt>>();

        // -- scan for endpoints --
        for (LinkedList<Segment> section : candidateSections) {

            boolean checkFirstSegment = true;
            
            Set<PairInt> endPoints = checkForAdjacentEndpoints(section,
                exclude, checkFirstSegment);

            if (endPoints == null || endPoints.isEmpty()) {
                continue;
            }

            Set<PairInt> endPoints2 = checkForAdjacentEndpoints(
                section, endPoints, checkFirstSegment);

            if (endPoints2 == null || endPoints2.isEmpty()) {
                continue;
            }

            // add all section and endpoints to a set to add to output
            Set<PairInt> allPoints = new HashSet<PairInt>();
            allPoints.addAll(endPoints);
            allPoints.addAll(endPoints2);
            for (Segment segment : section) {
                allPoints.addAll(segment.points);
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

        For each segment group which has 2 matching endpoints, those should
        be stored as butterfly sections in a set.  Each one of those
        should be passed back in a list as the return of this method.

        runtime complexity is linear in the number of points in the given
        closed curve.
        */

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

    private void fillWithPattern(int x, int y, Pattern pattern, Segment 
        outputSegment) {
        
        for (PairInt p : pattern.ones) {
            PairInt p2 = new PairInt(x + p.getX(), y + p.getY());
            outputSegment.points.add(p2);
        }
        outputSegment.points.add(new PairInt(x, y));        
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

    private Set<PairInt> checkForAdjacentEndpoints(LinkedList<Segment> section, 
        Set<PairInt> exclude, boolean checkFirstSegment) {
        
        Segment segment = checkFirstSegment ? section.getFirst() : section.getLast();
        
        throw new UnsupportedOperationException("not yet implemented");
        
        /*
        if (segment instanceof VertSegment) {
            return checkForAdjacentEndpointsForVertSegment(section, exclude, 
                checkFirstSegment);
        } else if (segment instanceof HorizSegment) {
            return checkForAdjacentEndpointsForHorizSegment(section, exclude, 
                checkFirstSegment);
        } else if (segment instanceof UUDiagSegment) {
            return checkForAdjacentEndpointsForUUDiagSegment(section, exclude, 
                checkFirstSegment);
        } else if (segment instanceof ULDiagSegment) {
            return checkForAdjacentEndpointsForULDiagSegment(section, exclude, 
                checkFirstSegment);
        }
        
        return null;
        */
        
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
        LinkedList<Segment> section, Set<PairInt> exclude, 
        boolean checkFirstSegment) {
        
        VertSegment segment = checkFirstSegment ? (VertSegment)section.getFirst() :
            (VertSegment)section.getLast();
        
        throw new UnsupportedOperationException("not yet implemented");
        
        /*
        endpoints for vert:
                -        - -       -
              # - #    # - - #   # - #
              . .        . .       . .
        
                         0
              .  .      -2
           -  2  1  -   -1
           -  3  0  -    0
              .  .       1
          -2 -1  0  1
        */
        
    }

    /*
    may change these classes to have ordered points or to specify the
    indexes of points that are the connections, that is the '.'s in sketches
    below.
    */
    public static class Segment {
        // the 4 points matching the segment
        Set<PairInt> points = new HashSet<PairInt>();
    }
    public static class VertSegment extends Segment {
        /*              0
              .  .      1
           -  2  1  -   2
           -  3  0  -   3
              .  .      4
        0  1  2  3  4
        */
    }
    public static class HorizSegment extends Segment {
        /*             0
              -  -     1
           .  2  1  .  2
           .  3  0  .  3
              -  -     4
        0  1  2  3  4
        */
    }
    public static class UUDiagSegment extends Segment {
        /*
             - -     3
           - - 2 1 - 2
           - 3 0 - - 1
               - -   0
        0  1 2 3 4 5
        */
    }
    public static class ULDiagSegment extends Segment {
        /*
               - -   3
           - 2 1 - - 2
           - - 3 0 - 1
             - -     0
        0  1 2 3 4 5
        */
    }

    public static class Pattern {
        Set<PairInt> ones;
        Set<PairInt> zeroes;
    }

    private Segment checkVertSegmentPattern(int x, int y, Set<PairInt> neighbors) {

        Pattern pattern = getVertSegmentPattern();
        
        boolean matchesPattern = matchesPattern(x, y, neighbors, pattern);
        
        if (matchesPattern) {
            VertSegment segment = new VertSegment();
            fillWithPattern(x, y, pattern, segment);
            return segment;
        }
        
        swapYDirection(pattern);
        
        matchesPattern = matchesPattern(x, y, neighbors, pattern);
        
        if (matchesPattern) {
            VertSegment segment = new VertSegment();
            fillWithPattern(x, y, pattern, segment);
            return segment;
        }
        
        return null;
    }
    
    protected Pattern getVertSegmentPattern() {

        /*               0
              .  .      -2
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

        /*               0
              -  -      -2
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

        pr.ones.add(new PairInt(-2, 0)); pr.ones.add(new PairInt(-2, -1));
        pr.ones.add(new PairInt(1, 0)); pr.ones.add(new PairInt(1, -1));

        return pr;
    }
    
    protected Pattern getUUDiagSegmentPattern() {

        /*                 0
                -  -      -2
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

    private Segment checkHorizSegmentPattern(int x, int y, Set<PairInt> neighbors) {

        Pattern pattern = getHorizSegmentPattern();
        
        boolean matchesPattern = matchesPattern(x, y, neighbors, pattern);
        
        if (matchesPattern) {
            HorizSegment segment = new HorizSegment();
            fillWithPattern(x, y, pattern, segment);
            return segment;
        }
        
        swapXDirection(pattern);
        
        matchesPattern = matchesPattern(x, y, neighbors, pattern);
        
        if (matchesPattern) {
            HorizSegment segment = new HorizSegment();
            fillWithPattern(x, y, pattern, segment);
            return segment;
        }
        
        return null;
    }
    
    private Segment checkDiagSegmentPattern(int x, int y, Set<PairInt> neighbors) {

        Pattern pattern = getUUDiagSegmentPattern();
        
        boolean matchesPattern = matchesPattern(x, y, neighbors, pattern);
        
        if (matchesPattern) {
            UUDiagSegment segment = new UUDiagSegment();
            fillWithPattern(x, y, pattern, segment);
            return segment;
        }
        
        swapXDirection(pattern);
        
        matchesPattern = matchesPattern(x, y, neighbors, pattern);
        
        if (matchesPattern) {
            ULDiagSegment segment = new ULDiagSegment();
            fillWithPattern(x, y, pattern, segment);
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
