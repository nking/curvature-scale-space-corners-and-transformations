package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class ClosedCurveContourMatcher0 {
    
    /*
    O(N) pattern which assumes that points in contour1[i] match points
    in contour2[i + deltaIdx] in an ordered manner.

    matching points in the contour lists:
       contour1=contourList1[0]        curve2=contourList2[0]
          tries deltaIdx=0 with contour1[i] : contour2[i] for all i in contour1
          tries deltaIdx=1 with contour1[i] : contour2[i + 1] for all i in contour1
             ...
          tries deltaIdx=n-1 with contour1[i] : contour2[i + (n-1)] for all i in contour1
          and the best becomes the cost for contour1=contourList1[0]    contour2=contourList2[0]
       Note that the matching is not tried if contour1 and contour2 sizes are different
       proceeds with attempt for next in contour list matches,
       contour1=contourList1[0]        curve2=contourList2[1]
       to find the best solution if any for contourList1[0] and stores it.
       Then does the same for contourList1[1]...
       And uses various statistics to find the best solution and
       combine similar with it to make a large point set to solve
       the euclidean transformation with.
    */
       
}
