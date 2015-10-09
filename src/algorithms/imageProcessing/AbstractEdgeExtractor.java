package algorithms.imageProcessing;

import algorithms.QuickSort;
import algorithms.compGeometry.ButterflySectionFinder.Routes;
import algorithms.util.PairIntArray;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Edge extractor operates on an image that has already been reduced to 
 * single pixel width lines, with the exception of single step patterns. 
 * 
 * If an edgeGuideImage is provided, the last step will be an attempt to
 * move pixels in the edges towards the highest intensity pixels in the 
 * edgeGuide image if the pixel move does not break the edge.
   The guidance image helps corrects errors in
 * the input image due to a line thinner like the ErosionFilter.  A Zhang-Suen
 * thinner is used now, so the guide image is not necessary.
 * 
 * @author nichole
 */
public abstract class AbstractEdgeExtractor implements IEdgeExtractor {
    
    protected final GreyscaleImage img;
    
    protected GreyscaleImage edgeGuideImage = null;
        
    protected long numberOfPixelsAboveThreshold = 0;
    
    protected boolean removeShorterEdges = false;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
            
    /**
     * if the image is smaller than 100 on a side, this will be lowered to 5
     */
    protected int edgeSizeLowerLimit = 15;
    
    protected boolean repeatConnectAndTrim = false;
    
    /**
     * NOTE:  input should have a black (empty) background and edges should
     * have values > 125 counts.  Edges should also have width of 1 and no larger.
     * 
     * @param input 
     */
    public AbstractEdgeExtractor(GreyscaleImage input) {
        
        img = input;
        
        if (img.getWidth() < 100 || img.getHeight() < 100) {
            edgeSizeLowerLimit = 5;
        }
    }
    
    /**
     * NOTE:  input should have a black (empty) background and edges should
     * have values > 125 counts.  Edges should also have width of 1 and no larger.
     * The guide image is used to alter the extracted edges back towards the
     * highest intensity pixels of the guide image.  The guide image is expected
     * to be the combined X and Y gradient image from earlier processing
     * stages.
     * 
     * @param input
     * @param anEdgeGuideImage
     */
    public AbstractEdgeExtractor(GreyscaleImage input, 
        final GreyscaleImage anEdgeGuideImage) {
        
        img = input;
        
        edgeGuideImage = anEdgeGuideImage;
        
        if (img.getWidth() < 100 || img.getHeight() < 100) {
            edgeSizeLowerLimit = 5;
        }
    }
    
    @Override
    public GreyscaleImage getImage() {
        return img;
    }
    
    @Override
    public void overrideEdgeSizeLowerLimit(int length) {
        edgeSizeLowerLimit = length;
    }
    
    @Override
    public List<PairIntArray> findEdges() {
        
        List<PairIntArray> output = connectPixelsViaDFS();
        
        output = findEdgesIntermediateSteps(output);
        
        if (removeShorterEdges) {
            removeEdgesShorterThan(output, edgeSizeLowerLimit);
        } else {
            removeEdgesShorterThan(output, 1);
        }
                
        return output;
    }
    
    public void removeShorterEdges(boolean doRemove) {
        removeShorterEdges = doRemove;
    }
    
    protected abstract List<PairIntArray> findEdgesIntermediateSteps(
        List<PairIntArray> edges);
            
    /**
     * find the edges and return as a list of points.  The method uses a
     * DFS search through all points in the image with values > 0 to link
     * adjacent sequential points into edges.
     * 
     * @return 
     */
    protected List<PairIntArray> connectPixelsViaDFS() {
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        // DFS search for sequential neighbors.
        
        Stack<PairInt> stack = new Stack<PairInt>();
      
        int thresh0 = 1;
        
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                int v = img.getValue(i, j);
                if (v >= thresh0) {
                    stack.add(new PairInt(i, j));
                }
            }
        }
        
        numberOfPixelsAboveThreshold = stack.size();
        
        log.log(Level.FINE, 
            "Number of pixels that meet or exceed threshhold={0}", 
            Long.toString(numberOfPixelsAboveThreshold));
        
        List<PairIntArray> output = new ArrayList<PairIntArray>();
        int[] uNodeEdgeIdx = new int[img.getWidth() * img.getHeight()];
        Arrays.fill(uNodeEdgeIdx, -1);
                   
        // > O(N) and << O(N^2)
        while (!stack.isEmpty()) {
            
            PairInt uNode = stack.pop();
            
            int uX = uNode.getX();
            int uY = uNode.getY();
            int uIdx = img.getIndex(uX, uY);
                                    
            for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                
                int vX = dxs[nIdx] + uX;
                int vY = dys[nIdx] + uY;

                if ((vX < 0) || (vX > (w - 1)) || (vY < 0) || (vY > (h - 1))) {
                    continue;
                }

                int vIdx = img.getIndex(vX, vY);
                
                if (uNodeEdgeIdx[vIdx] != -1 || (uIdx == vIdx)) {
                    continue;
                }

                if (img.getValue(vX, vY) < thresh0) {
                    continue;
                }
                   
                processNeighbor(uX, uY, uIdx, vX, vY, vIdx, uNodeEdgeIdx, output);
                
                stack.add(new PairInt(vX, vY));
            }   
        }
        
        log.fine(output.size() + " edges after DFS");
        
        // count the number of points in edges
        long sum = countPixelsInEdges(output);
        
        log.log(Level.FINE, 
            "==> {0} pixels are in edges out of {1} pixels > threshhold", 
            new Object[]{Long.toString(sum), 
                Long.toString(numberOfPixelsAboveThreshold)});
        
        log.log(Level.FINE, "there are {0} edges", 
            Integer.toString(output.size()));
        
        return output;
    }
    
    protected long countPixelsInEdges(List<PairIntArray> edges) {
        long sum = 0;
        for (PairIntArray edge : edges) {
            sum += edge.getN();
        }
        return sum;
    }

    protected void removeEdgesShorterThan(List<PairIntArray> output, 
        int minNumberOfPixelsInEdge) {
        
        for (int i = (output.size() - 1); i > -1; i--) {
            if (output.get(i).getN() < minNumberOfPixelsInEdge) {
                output.remove(i);
            }
        }
    }
    
    protected void adjustEdgesTowardsBrightPixels(List<PairIntArray> tmpEdges) {
        
        if (edgeGuideImage == null) {
            return;
        }
        
        int nReplaced = 1;
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        int nMaxIter = 100;
        int nIter = 0;
      
        while ((nIter < nMaxIter) && (nReplaced > 0)) {
           
           nReplaced = 0;
        
           for (int lIdx = 0; lIdx < tmpEdges.size(); lIdx++) {
            
                PairIntArray edge = tmpEdges.get(lIdx);

                if (edge.getN() < 3) {
                    continue;
                }

                int nEdgeReplaced = curveHelper.adjustEdgesTowardsBrighterPixels(
                    edge, edgeGuideImage);
                
                nReplaced += nEdgeReplaced;
            }

            log.fine("REPLACED: " + nReplaced + " nIter=" + nIter);

            nIter++;
        }
        
        //TODO: fix this method!
        curveHelper.removeRedundantPoints(tmpEdges);
        
        curveHelper.pruneAdjacentNeighborsTo2(tmpEdges);
        
        curveHelper.correctCheckeredSegments(tmpEdges);
    }

    /**
     * add (vX, vY) to output where (uX,uY) is the last point, else create
     * new edges where necessary.  returns the index
     * 
     * @param uX
     * @param uY
     * @param uIdx
     * @param vX
     * @param vY
     * @param vIdx
     * @param uNodeEdgeIdx
     * @param output
     * @return 
     */
    private void processNeighbor(int uX, int uY, int uIdx,
        int vX, int vY, int vIdx,
        int[] uNodeEdgeIdx, List<PairIntArray> output) {
        
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
            return;
        }

        appendToNode.add(vX, vY);

        uNodeEdgeIdx[vIdx] = uNodeEdgeIdx[uIdx];
        
    }
    
    /**
     * get the other endPoint within this edge.  NOTE that the method assumes
     * that endPoint is truly an endPoint of edge, so it will always return 
     * a value.
     * @param endPoint
     * @param edge
     * @return 
     */
    protected PairInt getOppositeEndPointOfEdge(PairInt endPoint, 
        PairIntArray edge) {
        
        int x = endPoint.getX();
        int y = endPoint.getY();
        
        int idx = 0;
        
        int n = edge.getN();
        
        if ((edge.getX(idx) == x) && (edge.getY(idx) == y)) {
            int x2 = edge.getX(n - 1);
            int y2 = edge.getY(n - 1);
            return new PairInt(x2, y2);
        } else {
            int x2 = edge.getX(0);
            int y2 = edge.getY(0);
            return new PairInt(x2, y2);
        }
    }
    
    /**
     * find the edges and return as a list of points.  The method uses a
     * DFS search through all points in the image with values > 0 to link
     * adjacent sequential points into edges.
     * 
     * @return 
     */
    protected List<PairIntArray> connectPixelsViaDFSForBounds(List<Routes> 
        butterFlySections) {
        
        /*
        the choice among neighbors prefers in order:
            large distance from centroid
        or  fewer neighbors
        */
        boolean preferDistFromCentroid = false;
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        // DFS search for sequential neighbors.
        
        Stack<PairInt> stack = new Stack<PairInt>();
        
        Set<PairInt> points = new HashSet<PairInt>();
      
        int thresh0 = 1;
        
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                int v = img.getValue(i, j);
                if (v >= thresh0) {
                    PairInt p = new PairInt(i, j); 
                    stack.add(p);
                    points.add(p);
                }
            }
        }
               
        Map<PairInt, Integer> pointRoutesMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < butterFlySections.size(); ++i) {
            Integer index = Integer.valueOf(i);
            Routes r = butterFlySections.get(i);
            /*for (PairInt p : r.getRoute0()) {
                pointRoutesMap.put(p, index);
            }
            for (PairInt p : r.getRoute1()) {
                pointRoutesMap.put(p, index);
            }*/
            pointRoutesMap.put(r.getEP0(), index);
            pointRoutesMap.put(r.getEP0End(), index);
            pointRoutesMap.put(r.getEP1(), index);
            pointRoutesMap.put(r.getEP1End(), index);
        }
        Set<PairInt> added = new HashSet<PairInt>();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        double[] xyCen = curveHelper.calculateXYCentroids(points);
        
        numberOfPixelsAboveThreshold = stack.size();
        
        log.log(Level.FINE, 
            "Number of pixels that meet or exceed threshhold={0}", 
            Long.toString(numberOfPixelsAboveThreshold));
        
        List<PairIntArray> output = new ArrayList<PairIntArray>();
        int[] uNodeEdgeIdx = new int[img.getWidth() * img.getHeight()];
        Arrays.fill(uNodeEdgeIdx, -1);
        
        int[] neighborsX = new int[8];
        int[] neighborsY = new int[8];
        double[] neighborsNNon = new double[8];
        double[] neighborsDistCen = new double[8];
                   
        // > O(N) and << O(N^2)
        while (!stack.isEmpty()) {
            
            PairInt uNode = stack.pop();
            
            int uX = uNode.getX();
            int uY = uNode.getY();
            int uIdx = img.getIndex(uX, uY);
            
            int count = 0;
            
            // first, a separate quick check for whether a point is in
            //  pointRoutesMap and not in added and route by that if so
            if (!pointRoutesMap.isEmpty()) {
                // checking endpoints only of routes and adding entire route
                // if there is a match
                boolean didAddPoints = false;
                for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                    int vX = dxs[nIdx] + uX;
                    int vY = dys[nIdx] + uY;
                    if ((vX < 0) || (vX > (w - 1)) || (vY < 0) || (vY > (h - 1))) {
                        continue;
                    }
                    int vIdx = img.getIndex(vX, vY);                
                    if (uNodeEdgeIdx[vIdx] != -1 || (uIdx == vIdx) || 
                        (img.getValue(vX, vY) < thresh0)) {
                        continue;
                    }
                    PairInt vNode = new PairInt(vX, vY);
                    if (added.contains(vNode)) {
                        continue;
                    }
                    Integer routeIndex = pointRoutesMap.get(vNode);
                    if (routeIndex != null) {
                        Routes routes = butterFlySections.get(routeIndex.intValue());
                        boolean addRoute0 = routes.getRoute0().contains(vNode);
                        if (addRoute0) {
                            for (PairInt chk : routes.getRoute0()) {
                                if (added.contains(chk)) {
                                    throw new IllegalStateException(
                                    "error in algorithm:" + 
                                     " a point in route0 has already been added");
                                }
                                processNeighbor(uX, uY, uIdx, chk.getX(), chk.getY(),
                                    vIdx, uNodeEdgeIdx, output);
                                added.add(uNode);
                                added.add(chk);
                                stack.add(chk);
                                didAddPoints = true;
                            }
                        } else {
                            for (PairInt chk : routes.getRoute1()) {
                                if (added.contains(chk)) {
                                    throw new IllegalStateException(
                                    "error in algorithm:" + 
                                     " a point in route0 has already been added");
                                }
                                processNeighbor(uX, uY, uIdx, chk.getX(), chk.getY(),
                                    vIdx, uNodeEdgeIdx, output);
                                added.add(uNode);
                                added.add(chk);
                                stack.add(chk);
                                didAddPoints = true;
                            }
                        }
                    }
                }
                if (didAddPoints) {
                    continue;
                }
            }
            
            for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                
                int vX = dxs[nIdx] + uX;
                int vY = dys[nIdx] + uY;

                if ((vX < 0) || (vX > (w - 1)) || (vY < 0) || (vY > (h - 1))) {
                    continue;
                }

                int vIdx = img.getIndex(vX, vY);
                
                if (uNodeEdgeIdx[vIdx] != -1 || (uIdx == vIdx)) {
                    continue;
                }

                if (img.getValue(vX, vY) < thresh0) {
                    continue;
                }
                
                if (added.contains(new PairInt(vX, vY))) {
                    continue;
                }
                
                int nc = curveHelper.countNeighbors(vX, vY, points, w, h);
                
                if (nc == 8) {
                    // not a border pixel
                    continue;
                }
                
                neighborsX[count] = vX;
                neighborsY[count] = vY;
                neighborsNNon[count] = 8 - nc;
                
                double diffX = vX - xyCen[0];
                double diffY = vY - xyCen[1];
                neighborsDistCen[count] = Math.sqrt(diffX*diffX + diffY*diffY);
                
                count++;
            }
            
            if (count == 0) {
                continue;
            }
            
            //sort ascending.  largest neighborsDistCen and largest neighborsNNon
            // are near end of array at index count - 1
            if (preferDistFromCentroid) {
                QuickSort.sortBy1stThen2nd(neighborsDistCen, neighborsNNon, 
                    neighborsX, neighborsY, 0, count - 1);
            } else {
                QuickSort.sortBy1stThen2nd(neighborsNNon, neighborsDistCen, 
                    neighborsX, neighborsY, 0, count - 1);
            }
            
            // add only the preferred to the edge and stack
            
            int vIdx = img.getIndex(neighborsX[count - 1], neighborsY[count - 1]);
                
            processNeighbor(uX, uY, uIdx, neighborsX[count - 1], 
                neighborsY[count - 1], vIdx, uNodeEdgeIdx, output);
            
            added.add(uNode);
            added.add(new PairInt(neighborsX[count - 1], neighborsY[count - 1]));
                
            stack.add(new PairInt(neighborsX[count - 1], neighborsY[count - 1]));
            
        }
        
        log.fine(output.size() + " edges after DFS");
        
        // count the number of points in edges
        long sum = countPixelsInEdges(output);
        
        log.log(Level.FINE, 
            "==> {0} pixels are in edges out of {1} pixels > threshhold", 
            new Object[]{Long.toString(sum), 
                Long.toString(numberOfPixelsAboveThreshold)});
        
        log.log(Level.FINE, "there are {0} edges", 
            Integer.toString(output.size()));
        
        return output;
    }
    
    protected void debugPrint(GreyscaleImage input, int xStart, int xStop,
        int yStart, int yStop) {
        
        StringBuilder sb = new StringBuilder();
                    
        for (int row = yStart; row <= yStop; row++) {
            sb.append(String.format("%3d: ", row));
            for (int col = xStart; col <= xStop; col++) {
                sb.append(String.format(" %3d ", input.getValue(col, row)));
            }
            sb.append(String.format("\n"));
        }
        log.info(sb.toString());
        
        sb = new StringBuilder(String.format("     "));
        for (int col = xStart; col <= xStop; col++) {
            sb.append(String.format(" %3d ", col));
        }
        sb.append(String.format("\n"));
        log.info(sb.toString());
    }

}
