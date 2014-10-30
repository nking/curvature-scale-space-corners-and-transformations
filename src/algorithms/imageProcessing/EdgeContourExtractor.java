package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 Edge contour extraction
    Local Methods:
        (1) At each edge pixel, a neighborhood (e.g., 3x3) is examined.

        (2) The center edge pixel can be linked with its neighbors if the 
            magnitude and direction differences are below certain thresholds 
            and their magnitudes are relatively large:

 * @author nichole
 */
public class EdgeContourExtractor {
    
    private final GreyscaleImage img;
    
    private long numberOfPixelsAboveThreshold = 0;
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    private boolean doNotStraigthenLines = false;
    
    private boolean useLineDrawingMode = false;
    
    /**
     * if the image is smaller than 100 on a side, this will be lowered to 5
     */
    private int edgeSizeLowerLimit = 30;
    
    /**
     * NOTE:  input should have a black (empty) background and edges should
     * have blue > 125 counts.  Edges should also have width of 1 and no larger.
     * 
     * @param input 
     */
    public EdgeContourExtractor(GreyscaleImage input) {
        img = input;
        
        if (img.getWidth() < 100 || img.getHeight() < 100) {
            edgeSizeLowerLimit = 5;
        }
    }
    
    public void useLineDrawingMode() {
        useLineDrawingMode = true;
    }
    
    public GreyscaleImage getImage() {
        return img;
    }
    
    /**
     * find the edges and return as a list of points.  The method uses a
     * DFS search through all points in the image with values > 0 to link
     * adjacent sequential points into edges.
     * The method then uses method mergeAdjacentEndPoints.
     * Note that the later 2 methods are not needed if the edges will be used
     * in a corner detector only, but if the edges are to be used to
     * find inflection points in scale space maps, those methods help to
     * provide more complete shapes and better matches between the same object
     * in another scale space map.
     * 
     * @return 
     */
    public List<PairIntArray> findEdges() {
        
        // DFS search for sequential neighbors.
        
        Stack<PairInt> stack = new Stack<PairInt>();
      
        int thresh0 = 1;
        
        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
                //for now, choosing to look only at the blue
                int bPix = img.getValue(i, j);
                if (bPix >= thresh0) {
                    stack.add(new PairInt(i, j));
                }
            }
        }
        
        numberOfPixelsAboveThreshold = stack.size();
        
        log.log(Level.FINE, "Number of pixels that exceed threshhold={0}", 
            Long.toString(numberOfPixelsAboveThreshold));
        
        List<PairIntArray> output = new ArrayList<PairIntArray>();
        int[] uNodeEdgeIdx = new int[img.getWidth() * img.getHeight()];
        Arrays.fill(uNodeEdgeIdx, -1);
           
        while (!stack.isEmpty()) {
            
            PairInt uNode = stack.pop();
            
            int uX = uNode.getX();
            int uY = uNode.getY();
            int uIdx = (uY * img.getWidth()) + uX;
                        
            boolean foundNeighbor = false;
            
            // for each neighbor v of u
            for (int vX = (uX - 1); vX < (uX + 2); vX++) {
                
                if (foundNeighbor) {
                    break;
                }
                
                if (vX < 0 || vX > (img.getWidth() - 1)) {
                    continue;
                }
                
                for (int vY = (uY - 1); vY < (uY + 2); vY++) {

                    if (vY < 0 || vY > (img.getHeight() - 1)) {
                        continue;
                    }
                    
                    int vIdx = (vY * img.getWidth()) + vX;
                
                    if (uNodeEdgeIdx[vIdx] != -1 || (uIdx == vIdx)) {
                        continue;
                    }
                    
                    if (img.getValue(vX, vY) < thresh0) {
                        continue;
                    }

                    // if u is not in an edge already, create a new one
                    if (uNodeEdgeIdx[uIdx] == -1) {
                        
                        PairIntArray edge = new PairIntArray();
                        
                        edge.add(uX, uY);
                        
                        uNodeEdgeIdx[uIdx] = output.size();

                        output.add(edge);                        
                    }
                    
                    // keep the curve points ordered
                                    
                    // add v to the edge u if u is the last node in it's edge
                                        
                    PairIntArray appendToNode = output.get(uNodeEdgeIdx[uIdx]);
                    int aIdx = appendToNode.getN() - 1;
                    if ((appendToNode.getX(aIdx) != uX) || 
                        (appendToNode.getY(aIdx) != uY)) {
                        continue;
                    }
                        
                    appendToNode.add(vX, vY);
                    
                    uNodeEdgeIdx[vIdx] = uNodeEdgeIdx[uIdx];
                                        
                    //TODO: do we only want 1 neighbor from the 9 as a continuation?
                    // yes for now, but this requires edges be only 1 pixel wide
                                
                   // inserting back at the top of the stack assures that the 
                   // search continues next from an associated point
                   stack.add(new PairInt(vX, vY));
                   
                   foundNeighbor = true;
                   
                   break;
                }
            }
        }
        
        log.fine(output.size() + " edges after DFS");
        
        int nIterMax = 100;
        int n, sz, lastSize;
        
        
        // count the number of points in edges
        long sum = countPixelsInEdges(output);
        
        log.log(Level.FINE, 
            "==> {0} pixels are in edges out of {1} pixels > threshhold", 
            new Object[]{Long.toString(sum), 
                Long.toString(numberOfPixelsAboveThreshold)});
        
        log.log(Level.FINE, "there are {0} edges", 
            Integer.toString(output.size()));
        
        output = mergeAdjacentEndPoints(output);
        
        log.log(Level.FINE, "{0} edges after merge adjacent", 
            Integer.toString(output.size()));
        
        /*
        //not necessary now that lines from CannyEdgeFilter are 1 pix width.
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        
        n = 0;
        sz = output.size();
        lastSize = Integer.MAX_VALUE;
        while ((sz < lastSize) && (n < nIterMax)) {
            
            lastSize = sz;
        
            output = ch.pruneAndIncludeAdjacentCurves(output, img.getWidth());
           
            sz = output.size();
            
            log.log(Level.FINE, "{0}) {1} edges after prune overlapping", 
                new Object[]{Integer.toString(n), Integer.toString(sz)});
            
            n++;
        }
        */
     
        if (useLineDrawingMode) {
            
            output = connectClosestPointsIfCanTrim(output);
        
            log.fine(output.size() + " edges after connect closest");
        }
        
        /*
        //TODO:
        // This helps to merge edges (that is extracted curves) at adjacent 
        // points that resemble an intersection of the lines, but it's not 
        // necessarily useful because the curvature is determined correctly 
        // whether the curves are merged or not.
        // If connecting the edges becomes more important, considering
        // using this, connectClosestPointsIfCanTrim():
        
        output = connectClosestPointsIfCanTrim(output);
        
        log.fine(output.size() + " edges after connect closest");
        */
        
        
        output = fillInGaps(output);
        
        log.log(Level.FINE, "{0} edges after fill in gaps", 
            new Object[]{Integer.toString(output.size())});
         
                
        //TODO:  this may need to change
        removeEdgesShorterThan(output, edgeSizeLowerLimit);
        sz = output.size();
        
        
        log.log(Level.FINE, "{0} edges after removing those shorter", 
            new Object[]{Integer.toString(sz)});

        
        long sum2 = countPixelsInEdges(output);
               
        log.log(Level.FINE, 
            "==> {0}) pixels are in edges out of {1} pixels > threshhold", 
            new Object[]{Long.toString(sum2), 
                Long.toString(numberOfPixelsAboveThreshold)});
        
        if (!doNotStraigthenLines && !useLineDrawingMode) {
            output = straightenJaggedLines(output);
        }
        
        return output;
    }

    public void doNotStraigthenLines() {
        this.doNotStraigthenLines = true;
    }
    
    private long countPixelsInEdges(List<PairIntArray> edges) {
        long sum = 0;
        for (PairIntArray edge : edges) {
            sum += edge.getN();
        }
        return sum;
    }
   
    /**
     * merge edges adjacent end points of curves
     * 
     * @param edges
     * @return 
     */
    protected List<PairIntArray> mergeAdjacentEndPoints(List<PairIntArray> 
        edges) {
     
        /*
        compare end of uEdge to beginning of all others
        reverse uEdge and repeat
        revert reverse for next start
        */
      
        boolean[] removed = new boolean[edges.size()];
        
        List<PairIntArray> output = new ArrayList<PairIntArray>();

        for (int i = 0; i < edges.size(); i++) {
            
            if (removed[i]) {
                continue;
            }
            
            PairIntArray uEdge = edges.get(i);
            
            // an extra iteration for reversing uEdge
            for (int r = 0; r < 2; r++) {
                 
                // compare bottom of uEdge to top of vEdge

                for (int j = 0; j < edges.size(); j++) {

                    if (i == j) {
                        continue;
                    }
                    if (removed[j]) {
                        continue;
                    }

                    PairIntArray vEdge = edges.get(j);

                    // recalculate in case u has grown
                    int uX = uEdge.getX(uEdge.getN() - 1);
                    int uY = uEdge.getY(uEdge.getN() - 1);

                    int vX = vEdge.getX(0);
                    int vY = vEdge.getY(0);

                    int diffX = uX - vX;
                    if (diffX < 0) {
                        diffX *= -1;
                    }

                    if (diffX > 1) {
                        continue;
                    }

                    int diffY = uY - vY;
                    if (diffY < 0) {
                        diffY *= -1;
                    }

                    if (diffY > 1) {
                        continue;
                    }

                    for (int k = 0; k < vEdge.getN(); k++) {
                        uEdge.add(vEdge.getX(k), vEdge.getY(k));
                    }

                    removed[j] = true;

                    // have to restart the j iteration to re-compare terms
                    j = -1;
                }
                
                if (r == 0) {
                    // just finished forward, start revers 
                    uEdge.reverse();
                } else if (r == 1) {
                    // revert the array back to other direction
                    uEdge.reverse();
                }
            }
        }
        
        for (int i = 0; i < edges.size(); i++) {
            if (!removed[i]) {
                output.add(edges.get(i));
            }
        }
        
        return output;
    }
    
    /**
     * fill in gaps of '1' pixel
     * 
     * @param edges
     * @return 
     */
    protected List<PairIntArray> fillInGaps(List<PairIntArray>
        edges) {
            
        /*
        similar to the mergeAdjacentEndPoints:  
        compare end of uEdge to beginning of all others
        reverse uEdge and repeat
        revert reverse for next start
        */
      
        boolean[] removed = new boolean[edges.size()];
        
        List<PairIntArray> output = new ArrayList<PairIntArray>();
                
        for (int i = 0; i < edges.size(); i++) {
            
            if (removed[i]) {
                continue;
            }
            
            PairIntArray uEdge = edges.get(i);
            
            // an extra iteration for reversing uEdge
            for (int r = 0; r < 2; r++) {
                 
                // compare bottom of uEdge to top of vEdge

                for (int j = 0; j < edges.size(); j++) {

                    if (i == j) {
                        continue;
                    }
                    if (removed[j]) {
                        continue;
                    }

                    PairIntArray vEdge = edges.get(j);

                    // recalculate in case u has grown
                    int uX = uEdge.getX(uEdge.getN() - 1);
                    int uY = uEdge.getY(uEdge.getN() - 1);

                    int vX = vEdge.getX(0);
                    int vY = vEdge.getY(0);

                    int diffX = uX - vX;
                    if (diffX < 0) {
                        diffX *= -1;
                    }

                    if (diffX > 2) {
                        continue;
                    }

                    int diffY = uY - vY;
                    if (diffY < 0) {
                        diffY *= -1;
                    }

                    if (diffY > 2) {
                        continue;
                    }
                    
                    int gapX = (uX + vX)/2;
                    int gapY = (uY + vY)/2;
                                        
                    img.setValue(gapX, gapY, 255);
                    
                    uEdge.add(gapX, gapY);
                    
                    for (int k = 0; k < vEdge.getN(); k++) {
                        uEdge.add(vEdge.getX(k), vEdge.getY(k));
                    }

                    removed[j] = true;

                    // have to restart the j iteration to re-compare terms
                    j = -1;
                }
                
                if (r == 0) {
                    // just finished forward, start revers 
                    uEdge.reverse();
                } else if (r == 1) {
                    // revert the array back to other direction
                    uEdge.reverse();
                }
            }
        }
        
        for (int i = 0; i < edges.size(); i++) {
            if (!removed[i]) {
                output.add(edges.get(i));
            }
        }
        
        return output;
    }
    
    /**
     * find the closest point between the curve0 and curve1 and return the
     * points in curve0XY and curve1XY along with the method return value
     * which is the separation.
     * 
     * @param curve0 
     * @param curve1
     * @param curve0Idx output variable to hold index to the (x, y) of the point 
     * in curve0 which is closest to curve0
     * @param curve1Idx output variable to hold index to the (x, y) of the point 
     * in curve0 which is closest to curve1
     * @return the separation between the closest pair of points in curve0 and 
     *   curve1
     */
    protected double findClosestPair(PairIntArray curve0, PairIntArray curve1, 
        int[] curve0Idx, int[] curve1Idx) {
        
        int idx0 = -1;
        int idx1 = -1;
        double min = Double.MAX_VALUE;
        
        for (int i = 0; i < curve0.getN(); i++) {
            int x0 = curve0.getX(i);
            int y0 = curve0.getY(i);
            for (int j = 0; j < curve1.getN(); j++) {
                int x1 = curve1.getX(j);
                int y1 = curve1.getY(j);
                int dx = x1 - x0;
                int dy = y1 - y0;
                double d = dx*dx + dy*dy;
                if (d < min) {
                    min = d;
                    idx0 = i;
                    idx1 = j;
                }
            }
        }
        
        curve0Idx[0] = idx0;
        
        curve1Idx[0] = idx1;
        
        return Math.sqrt(min);
    }
    
    protected List<PairIntArray> connectClosestPointsIfCanTrim(
        List<PairIntArray> edges) {
     
        double sqrtTwo = Math.sqrt(2) + 0.01;
        //double gapOfOne = 2*Math.sqrt(2) + 0.01;
        
        int[] edge0Idx = new int[1];
        int[] edge1Idx = new int[1];
        
        boolean[] removed = new boolean[edges.size()];
        
        List<PairIntArray> output = new ArrayList<PairIntArray>();
        
        for (int i = 0; i < edges.size(); i++) {
            
            if (removed[i]) {
                continue;
            }
                            
            PairIntArray edge0 = edges.get(i);
                        
            for (int j = (i + 1); j < edges.size(); j++) {
                
                if (removed[j]) {
                    continue;
                }
                
                PairIntArray edge1 = edges.get(j);
                
                double sep = findClosestPair(edge0, edge1, edge0Idx, edge1Idx);

                if (sep < sqrtTwo) {
                    // do not merge them if the points are not near the
                    // ends of the points sets.
                    float closestFrac0 = (float)edge0Idx[0]/(float)edge0.getN();
                    
                    float closestFrac1 = (float)edge1Idx[0]/(float)edge1.getN();
                    
                    if (((closestFrac0 > 0.07) && (closestFrac0 < 0.93))) {
                        continue;
                    }
                    if (((closestFrac1 > 0.07) && (closestFrac1 < 0.93))) {
                        continue;
                    }
                    
                    boolean closest0IsNearTop = ((float)(edge0Idx[0]/
                        (edge0.getN() - edge0Idx[0]))) <= 0.5;
                    
                    if (closest0IsNearTop) {
                        // if we trim the top, is remaining bottom connected?
                        boolean isConnected = isRangeConnected(edge0, 
                            edge0Idx[0], edge0.getN() - 1);
                        if (!isConnected) {
                            continue;
                        }
                    } else {
                        // if we trim the bottom, is remaining top connected?
                        boolean isConnected = isRangeConnected(edge0, 0, 
                            edge0Idx[0]);                      
                        if (!isConnected) {
                            continue;
                        }
                    }
                                            
                    boolean closest1IsNearTop = ((float)(edge1Idx[0]/
                        (edge1.getN() - edge1Idx[0]))) <= 0.5;
                    
                    if (closest1IsNearTop) {
                        // if we trim the top, is remaining bottom connected?
                        boolean isConnected = isRangeConnected(edge1, 
                            edge1Idx[0], edge1.getN() - 1);
                        if (!isConnected) {
                            continue;
                        }
                    } else {
                        // if we trim the bottom, is remaining top connected?
                        boolean isConnected = isRangeConnected(edge1, 0, 
                            edge1Idx[0]);                      
                        if (!isConnected) {
                            continue;
                        }
                    }
                    
                    // if here, then can trim outside the closest points in the
                    // edges and merge the edges into edge0, and remove edge1
                    
                    if (closest0IsNearTop) {
                        if (edge0Idx[0] > 0) {
                            edge0.removeRange(0, edge0Idx[0] - 1);
                        }
                    } else {
                        if (edge0Idx[0] < (edge0.getN() - 1)) {
                            edge0.removeRange(edge0Idx[0] + 1, edge0.getN() - 1);
                        }
                    }
                    
                    //TODO:  could remove this step and adjust the add, but 
                    //       easier maintainence this way
                    if (closest1IsNearTop) {
                        if (edge1Idx[0] > 0) {
                            edge1.removeRange(0, edge1Idx[0] - 1);
                        }
                    } else {
                        if (edge1Idx[0] < (edge1.getN() - 1)) {
                            edge1.removeRange(edge1Idx[0] + 1, edge1.getN() - 1);
                        }
                    }
                    
                    if (closest0IsNearTop) {
                        // insert edge1 at top of edge0
                        edge0.insertSpaceAtTopOfArrays(edge1.getN());
                        
                        // if edge1 closest is at bottom of it's edge, just add,
                        // else reverse then add
                        if (closest1IsNearTop) {
                            edge1.reverse();
                        }
                        for (int k = 0; k < edge1.getN(); k++) {
                            edge0.set(k, edge1.getX(k), edge1.getY(k));
                        }
                    } else {
                        // append edge1 to bottom of edge0
                        
                        if (!closest1IsNearTop) {
                            edge1.reverse();
                        }
                        for (int k = 0; k < edge1.getN(); k++) {
                            edge0.add(edge1.getX(k), edge1.getY(k));
                        }
                    }
                  
                    removed[j] = true;
                    
                    // have to restart the j iteration to re-compare terms
                    j = i;
                }
            }
            output.add(edge0);
        }
        
        return output;
    }
    
    /**
     * check that points within index idxLo and idxHi, inclusive, are 
     * consecutively within 1 pixel of adjacent indexes.
     * @param edge
     * @param idxLo
     * @param idxHi
     * @return 
     */
    protected boolean isRangeConnected(PairIntArray edge, int idxLo, int idxHi) {
        
        for (int i = (idxLo + 1); i <= idxHi; i++) {
            int x0 = edge.getX(i - 1);            
            int x1 = edge.getX(i);
            int diffX = x0 - x1;
            if (diffX < 0) {
                diffX *= -1;
            }
            if (diffX > 1) {
                return false;
            }
            
            int y0 = edge.getY(i - 1);
            int y1 = edge.getY(i);
            int diffY = y0 - y1;
            if (diffY < 0) {
                diffY *= -1;
            }
            if (diffY > 1) {
                return false;
            }
        }
        
        return true;
    }

    protected void removeEdgesShorterThan(List<PairIntArray> output, 
        int minNumberOfPixelsInEdge) {
        
        for (int i = (output.size() - 1); i > -1; i--) {
            if (output.get(i).getN() < minNumberOfPixelsInEdge) {
                output.remove(i);
            }
        }
    }

    private List<PairIntArray> straightenJaggedLines(List<PairIntArray> input) {
        
        for (int i = 0; i < input.size(); i++) {
                        
            PairIntArray edge = input.get(i);
            
            straightenJaggedLines(edge);
        }
        
        return input;
    }

    private void straightenJaggedLines(PairIntArray input) {
        
        if (input.getN() < 3) {
            return;
        }
        
        /*
        -- for each point:
           iterate from end of range to current point
              -- estimate a straight line from point to current end
                 -- do fastest determination to see if points in between fit line.
                    -- if they do, take the region in between and rewrite
                       the coordinates to follow the line.
                       -- store each index as 'straightened' to skip over
                          it upon subsequent iterations.
        */
                
        boolean[] straightened = new boolean[input.getN()];
        
        for (int i = 0; i < input.getN(); i++) {
                        
            if (straightened[i]) {
                continue;
            }
            
            // min number of points  i _ _ j
            for (int j = (input.getN() - 1); j > i; j--) {
                
                if ((j - i) < 10) {
                    continue;
                }
                if (straightened[j]) {
                    continue;
                }
                
                boolean didFitALineSegment = fitAndAlterForLineSegments(input, 
                    i, j);
                
                if (didFitALineSegment) {
                    for (int ii = i; ii <= j; ii++) {
                        straightened[ii] = true;
                    }
                    break;
                }
            }
        }        
    }

    private boolean fitAndAlterForLineSegments(PairIntArray input, 
        int start, int stop) {

        float thresh = 3.0f;
        
        int n = stop - start + 1;
        
        float x0 = input.getX(start);
        
        float y0 = input.getY(start);
                
        float xDiff = input.getX(stop) - x0;
        
        float yDiff = input.getY(stop) - y0;
        
        float slope = (xDiff == 0) ? 1 : (yDiff/xDiff);
        
        float[] yDiffs = new float[n];
        
        //(83,110) to (91,51)
        boolean debug = (input.getX(start) == 83) && (input.getY(start) == 110)
            && (input.getX(stop) == 91) && (input.getY(stop) == 51);
        
        // determine diff from line
        double avg = 0;
        
        // y = y0 + (x - x0)*slope
        for (int i = start; i < (start + n); i++) {
            
            float y = y0 + ((input.getX(i) - x0) * slope);
            
            float yd = input.getY(i) - y;
            
            avg += yd;
            
            yDiffs[i - start] = yd;
        }
        
        avg /= (double)n;
        
        double stDev = 0;
        
        for (int i = start; i < (start + n); i++) {
                        
            float yd = yDiffs[i - start];
            
            stDev += Math.pow((yd - avg), 2);
        }
        
        //N-1 because had to calculate mean from data
        stDev = Math.sqrt(stDev/(n - 1.0f));
        
       
        boolean hasOutliers = false;
        
        if ((Math.abs(avg) < 0.01) && (stDev < Math.abs(slope))) {
            
            for (int i = start; i < (start + n); i++) {

                float yd = yDiffs[i - start];
                if (yd < 0) {
                    yd *= -1;
                }

                /*
                log.fine(String.format("debug: yd=%f thresh*stdev=%f", 
                    yd, thresh*stDev));*/

                if (yd > (thresh * stDev)) {
                    hasOutliers = true;
                    break;
                }
            }
            
        } else {
            
            hasOutliers = true;
        }
        
        /*if (!hasOutliers)
        log.fine(String.format(
            "debug: %d %d (%d,%d) to (%d,%d) slope=%f avg=%f stdev=%f  n=%d hasOutliers=%s",
            start, stop, input.getX(start), input.getY(start), 
            input.getX(stop), input.getY(stop), slope, avg, stDev, n, 
            Boolean.toString(hasOutliers)));
        
        if (!hasOutliers) {
            log.fine("debug: " + input.toString());
        }*/
        
        boolean isALine = !hasOutliers;
        
        if (isALine) {
            
            if (Math.abs(slope) < 1) {
                for (int i = start + 1; i < (start + n); i++) {
                    int x = input.getX(i);                
                    float y = y0 + ((x - x0) * slope);
                                    
                    input.set(i, x, (int)y);
                }
            } else {
                for (int i = start + 1; i < (start + n); i++) {
                    int y = input.getY(i); 
                    //x = x0 + ((y-y0)/slope)
                    float x = x0 + ((y - y0) / slope);
                                    
                    input.set(i, (int)x, y);
                }
            }
            
            /*
            if (!hasOutliers) {
                log.fine("debug: " + input.toString());
            }*/
        }
        
        return isALine;
    }

}
