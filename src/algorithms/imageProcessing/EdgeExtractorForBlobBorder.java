package algorithms.imageProcessing;

import algorithms.compGeometry.ButterflySectionFinder;
import algorithms.compGeometry.PerimeterFinder;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * 
 * @see AbstractEdgeExtractor

 * @author nichole
 */
public class EdgeExtractorForBlobBorder {
            
    protected boolean debug = false;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public EdgeExtractorForBlobBorder() {
                    
    }
    
    public void setToDebug() {
        debug = true;
    }
    
    /**
     * NOT READY FOR USE.
     * given the set of contiguous points, find the perimeter of them and order
     * the points into a closed single pixel width curve.  
     * use extractAndOrderTheBorder0() instead for now.
     * @param contiguousPoints
     * @param imageWidth
     * @param imageHeight
     * @param discardWhenCavityIsSmallerThanBorder
     * @return 
     */
    PairIntArray extractAndOrderTheBorder(Set<PairInt> contiguousPoints,
        int imageWidth, int imageHeight, boolean discardWhenCavityIsSmallerThanBorder) {
        
        if (contiguousPoints == null) {
            return null;
        }
   
if (debug) {
MiscDebug.plotPoints(contiguousPoints, imageWidth, imageHeight, MiscDebug.getCurrentTimeFormatted());        
}

        // ------- define the perimeter of contiguousPoints ----------

        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        
        int imageMaxColumn = imageWidth - 1;
        int imageMaxRow = imageHeight - 1;
       
        int[] rowMinMax = new int[2];
                
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find(
            contiguousPoints, rowMinMax, imageMaxColumn, outputEmbeddedGapPoints);
       
        if (!outputEmbeddedGapPoints.isEmpty()) {
            // update the perimeter for "filling in" embedded points
            perimeterFinder.updateRowColRangesForAddedPoints(rowColRanges, 
                rowMinMax, imageMaxColumn, outputEmbeddedGapPoints);
        }
        
        // --- walk around the boundary of contiguousPoints, forming a closed curve path ----------
        PairIntArray out = perimeterFinder.getOrderedBorderPixels(
            contiguousPoints, rowColRanges, rowMinMax, imageMaxColumn, imageMaxRow);
  
        int nBorder = out.getN();
        int nCavity = contiguousPoints.size() + outputEmbeddedGapPoints.size()
            - nBorder;
               
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        double[] xyCen = curveHelper.calculateXYCentroids(out);  
        
log.fine(String.format("LIMIT: (%d,%d) nPerimeter=%d nCavity=%d", (int)Math.round(xyCen[0]),
(int)Math.round(xyCen[1]), nBorder, nCavity));

        log.fine("LIMIT: number of border points=" + nBorder + " nCavity=" + nCavity);
        
        if (discardWhenCavityIsSmallerThanBorder && (nBorder > nCavity)) {
            return null;
        }
        
if (debug) {       
MiscDebug.plotPoints(out, imageWidth, imageHeight, MiscDebug.getCurrentTimeFormatted());
}
        return out;
    }
    
    /**
     * given the set of contiguous points, find the perimeter of them and order
     * the points into a closed single pixel width curve.
     * (Note, still testing this, but for structures without many junctions,
     * the results looks good so far.)
     * @param contiguousPoints
     * @param imageWidth
     * @param imageHeight
     * @param discardWhenCavityIsSmallerThanBorder
     * @return 
     */
    public PairIntArray extractAndOrderTheBorder0(Set<PairInt> contiguousPoints,
        int imageWidth, int imageHeight, boolean discardWhenCavityIsSmallerThanBorder) {
        
        if (contiguousPoints == null) {
            return null;
        }
 
if (debug) {        
MiscDebug.plotPoints(contiguousPoints, imageWidth, imageHeight, MiscDebug.getCurrentTimeFormatted());        
}

        // ---- extract the border pixels from contiguousPoints ------

        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        
        int imageMaxColumn = imageWidth - 1;
        int imageMaxRow = imageHeight - 1;
       
        int[] rowMinMax = new int[2];
                
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find(
            contiguousPoints, rowMinMax, imageMaxColumn, outputEmbeddedGapPoints);
       
        if (!outputEmbeddedGapPoints.isEmpty()) {
            // update the perimeter for "filling in" embedded points
            perimeterFinder.updateRowColRangesForAddedPoints(rowColRanges, 
                rowMinMax, imageMaxColumn, outputEmbeddedGapPoints);
        }
       
        Set<PairInt> borderPixels = perimeterFinder.getBorderPixels(
            rowColRanges, rowMinMax, imageMaxColumn, imageMaxRow);
        
        int nBorder = borderPixels.size();
        int nCavity = contiguousPoints.size() + outputEmbeddedGapPoints.size()
            - nBorder;

long ts = MiscDebug.getCurrentTimeFormatted();
if (debug) {        
Image img3 = new Image(imageWidth, imageHeight);
for (PairInt p : borderPixels) {
    img3.setRGB(p.getX(), p.getY(), 255, 0, 0);
}
MiscDebug.writeImageCopy(img3, "border_perimeter_" + ts + ".png");
}
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        double[] xyCen = curveHelper.calculateXYCentroids(borderPixels);  
        
log.fine(String.format("LIMIT: (%d,%d) nPerimeter=%d nCavity=%d", (int)Math.round(xyCen[0]),
(int)Math.round(xyCen[1]), nBorder, nCavity));

        // expecting nBorder to be < nCavity

        if (discardWhenCavityIsSmallerThanBorder && (nBorder > nCavity)) {
           
            log.info(String.format(
                "discarding (%d, %d) number of border points=%d nCavity=%d ", 
                (int)Math.round(xyCen[0]), (int)Math.round(xyCen[1]), 
                nBorder, nCavity));
            
            return null;
        }
        
        // ----- remove spurs, merge curves, re-order points to form a 
        //       single closed curve if possible -------
/*
try {
Misc.persistToFile("blob_" + ts + ".dat", borderPixels);
int z = 1;
} catch(IOException e){}
*/
        // persist specific features to restore if thinned:
        ButterflySectionFinder finder = new ButterflySectionFinder();
        List<Set<PairInt>> butterFlySections = finder.findButterflySections(
            Misc.convertWithoutOrder(borderPixels));
        
        ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
        lt.applyLineThinner(borderPixels, 0, imageWidth, 0, imageHeight);
        PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
        pltc.correctForExtCorner(borderPixels, imageWidth, imageHeight);
        
        if (borderPixels.isEmpty()) {
            return null;
        }
        
        // restore butterFlySections if any
        if (butterFlySections != null && !butterFlySections.isEmpty()) {
            for (Set<PairInt> butterFlySection : butterFlySections) {
                borderPixels.addAll(butterFlySection);
            }
        }
        
if (debug) {        
Image img3 = new Image(imageWidth, imageHeight);
for (PairInt p : borderPixels) {
    img3.setRGB(p.getX(), p.getY(), 255, 0, 0);
}
MiscDebug.writeImageCopy(img3, "border_before_spur_removal_" + ts + ".png");
}
        SpurRemover spurRm = new SpurRemover();
        if (debug) {
            spurRm.setToDebug();
        }
        spurRm.remove(borderPixels, imageWidth, imageHeight);

        if (borderPixels.isEmpty()) {
            return null;
        }
        
if (debug) {        
Image img3 = new Image(imageWidth, imageHeight);
for (PairInt p : borderPixels) {
    img3.setRGB(p.getX(), p.getY(), 255, 0, 0);
}
MiscDebug.writeImageCopy(img3, "border_after_spur_removal_" + MiscDebug.getCurrentTimeFormatted() + ".png");
}

        UntraversableLobeRemover remover = new UntraversableLobeRemover();
        remover.applyFilter(borderPixels);
        
        if (borderPixels.isEmpty()) {
            return null;
        }
        
if (debug) {        
Image img3 = new Image(imageWidth, imageHeight);
for (PairInt p : borderPixels) {
    img3.setRGB(p.getX(), p.getY(), 255, 0, 0);
}
MiscDebug.writeImageCopy(img3, "border_after_untraversable_removal_" + ts + ".png");
}

        //xMin, xMax, yMin, yMax
        int[] minMaxXY = MiscMath.findMinMaxXY(borderPixels);
      
        int xOffset = minMaxXY[0] - 5;
        int yOffset = minMaxXY[2] - 5;
        int w = minMaxXY[1] - xOffset + 10;
        int h = minMaxXY[3] - yOffset + 10;
        
        GreyscaleImage img = new GreyscaleImage(w, h);
        for (PairInt p : borderPixels) {
            int x = p.getX() - xOffset;
            int y = p.getY() - yOffset;
            img.setValue(x, y, 1);
        }
        
        EdgeExtractorWithJunctions extractor = new EdgeExtractorWithJunctions(img);
        if (debug) {
            extractor.setToDebug();
        }
        extractor.overrideMaxNumberIterationsJunctionSplice(10);
        //List<PairIntArray> output = extractor.findEdges();
        PairIntArray out = extractor.findAsSingleClosedEdge();
       
        if (out == null) {
            return null;
        }
        
        /*
        some curves are originally composed of a single pixel width connector 
        between more than one potentially closed curve and the algorithm 
        eventually chooses the longest closed curve and follows the other out 
        until it ends without being adjacent to another point in the larger 
        curve.  in other words, these single pixel connections between
        multiple curves sometimes end in a long spur that should be trimmed.
        Note that the algorithm could be redesigned to look for such regions
        first and discard the smaller half with a warning.
        */
        if (!curveHelper.isAdjacent(out, 0, out.getN() - 1)) {
   
if (debug) {        
Image img3 = new Image(imageWidth, imageHeight);
for (int i = 0; i < out.getN(); ++i) {
    img3.setRGB(out.getX(i), out.getY(i), 255, 0, 0);
}
img3.setRGB(out.getX(0), out.getY(0), 255, 255, 0);
img3.setRGB(out.getX(out.getN() - 1), out.getY(out.getN() - 1), 255, 255, 0);
MiscDebug.writeImageCopy(img3, "border_trimEndpointSpurs.png");
}

            boolean altered = trimEndpointSpursIfAny(out, img.getWidth(), 
                img.getHeight());

if (debug) {        
Image img3 = new Image(imageWidth, imageHeight);
for (int i = 0; i < out.getN(); ++i) {
    img3.setRGB(out.getX(i), out.getY(i), 255, 0, 0);
}
img3.setRGB(out.getX(0), out.getY(0), 255, 255, 0);
img3.setRGB(out.getX(out.getN() - 1), out.getY(out.getN() - 1), 255, 255, 0);
MiscDebug.writeImageCopy(img3, "border_before_reorder_endpoints.png");
}

            if (!curveHelper.isAdjacent(out, 0, out.getN() - 1)) {
                
                extractor.reorderEndpointsIfNeeded(out);
                
if (debug) {        
Image img3 = new Image(imageWidth, imageHeight);
for (int i = 0; i < out.getN(); ++i) {
    img3.setRGB(out.getX(i), out.getY(i), 255, 0, 0);
}
img3.setRGB(out.getX(0), out.getY(0), 255, 255, 0);
img3.setRGB(out.getX(out.getN() - 1), out.getY(out.getN() - 1), 255, 255, 0);
MiscDebug.writeImageCopy(img3, "border_after_reorder_endpoints.png");
}               

                if (!curveHelper.isAdjacent(out, 0, out.getN() - 1)) {
                    
                    trimForMultipleClosedCurves(extractor, out, img.getWidth(), 
                        img.getHeight());
                    
if (debug) {        
Image img3 = new Image(imageWidth, imageHeight);
for (int i = 0; i < out.getN(); ++i) {
    img3.setRGB(out.getX(i), out.getY(i), 255, 0, 0);
}
img3.setRGB(out.getX(0), out.getY(0), 255, 255, 0);
img3.setRGB(out.getX(out.getN() - 1), out.getY(out.getN() - 1), 255, 255, 0);
MiscDebug.writeImageCopy(img3, "border_trim_multiple_curves.png");
int z = 1;
}                    

                }
            }
        }
        
        // add the shifts back
        //PairIntArray out = output.get(0);
        
        for (int i = 0; i < out.getN(); ++i) {
            int x = out.getX(i) + xOffset;
            int y = out.getY(i) + yOffset;
            out.set(i, x, y);
        }

if (debug) {
Image img3 = new Image(imageWidth, imageHeight);
for (int j = 0; j < out.getN(); ++j) {
    int x = out.getX(j);
    int y = out.getY(j);
    if (j == 0 || (j == (out.getN() - 1))) {
        ImageIOHelper.addPointToImage(x, y, img3, 0, 200, 150, 0);
    } else {
        ImageIOHelper.addPointToImage(x, y, img3, 0, 255, 0, 0);
    }
}
MiscDebug.writeImageCopy(img3, "output_" + MiscDebug.getCurrentTimeFormatted() + ".png");
}

        return out;
    }

    private boolean trimEndpointSpursIfAny(final PairIntArray curve, 
        final int imageWidth, final int imageHeight) {

        boolean edited = false;
        
        SpurRemover spurRemover = new SpurRemover();
        
        // try beginning of curve        
        while (true) {
            
if (debug) {        
Image img3 = new Image(imageWidth, imageHeight);
for (int i = 0; i < curve.getN(); ++i) {
    int x = curve.getX(i);
    int y = curve.getY(i);
    img3.setRGB(x, y, 255, 0, 0);
}
MiscDebug.writeImageCopy(img3, "before_trim_endpoints_" + MiscDebug.getCurrentTimeFormatted() + ".png");
}      

            Set<PairInt> points = Misc.convert(curve);
            
            int idx = 0;
                
            int x = curve.getX(idx);                
            int y = curve.getY(idx);
                
            boolean isASpur = spurRemover.isASpur(x, y, points, imageWidth,
                imageHeight);
                
            if (isASpur) {
                
                curve.removeRange(idx, idx);
                
                edited = true;
                
            } else {
                break;
            }
        }
        
        // try end of curve        
        while (true) {

if (debug) {        
Image img3 = new Image(imageWidth, imageHeight);
for (int i = 0; i < curve.getN(); ++i) {
    int x = curve.getX(i);
    int y = curve.getY(i);
    img3.setRGB(x, y, 255, 0, 0);
}
MiscDebug.writeImageCopy(img3, "before_trim_endpoints_" + MiscDebug.getCurrentTimeFormatted() + ".png");
}

            Set<PairInt> points = Misc.convert(curve);
            
            int idx = (curve.getN() - 1);
                
            int x = curve.getX(idx);                
            int y = curve.getY(idx);
                
            boolean isASpur = spurRemover.isASpur(x, y, points, imageWidth,
                imageHeight);
                
            if (isASpur) {
                
                curve.removeRange(idx, idx);
                
                edited = true;
                
            } else {
                break;
            }
        }
        
        return edited;
    }

    private boolean trimForMultipleClosedCurves(final EdgeExtractorWithJunctions 
        extractor, final PairIntArray curve, int imageWidth, int imageHeight) {

        int idxA = extractor.findAdjacentToTopAtBottom(curve);
        
        int idxB = extractor.findAdjacentToBottomAtTop(curve);
        
        if (idxA == -1 || idxB == -1) {
            return false;
        }
        
        /*
        TODO:
        revisit this...
        
        one case:
           idxA and idxB are both closer to same endpoint and are not the same
           number and the difference is small.
           this is potentially a knot in the curve, separated by a single
           width pixel curve (that therefore cannot be approaching the
           knot and leaving the knot).  If there's a single pixel width
           bottle neck in the region between them, can trim the whole
           section off.
        for idxA = 7 and idxB = 11,
        looking for when the path between them only has one choice, that is
        one neighbor other than the previous.  for that case, should be able
        to trim off the end section.
        
        0  41, 27
        1  41, 26
        2  41, 25
        3  40, 25
        
        6  39, 27
        7  40, 28 =======
        8  39, 28
           38, 28  <- *
        10 37, 29  <- *
        11 36, 29 =======
        12 36, 30
        ...
        n-1 35, 30
        */
        if (idxA == idxB) {
            return false;
        }
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        Set<PairInt> points = Misc.convert(curve);
        
        double limit = Math.sqrt(2);
        
        int mid = curve.getN() >> 1;
        
        boolean trimBeginning = ((idxA < mid) && (idxB < mid));
        
        boolean trimEnd = ((idxA > mid) && (idxB > mid));
        
        if (trimBeginning || trimEnd) {
            
            if (idxB < idxA) {
                int swap = idxB;
                idxB = idxA;
                idxA = swap;
            }
            
            // looking for the point between them where there is no
            // other point within a circular radius of 1            
            for (int i = idxA; i < idxB; ++i) {
                
                int x1 = curve.getX(i);
                int y1 = curve.getY(i);
                
                Set<PairInt> neighbors1 = curveHelper.findNeighbors(x1, y1, 
                    points);
                
                if (neighbors1.size() == 2) {
                    curve.removeRange(0, (idxB - 1));
                    return true;
                }
                
                PairInt p1 = new PairInt(x1, y1);
                                
                for (PairInt p2 : neighbors1) {
                    
                    int x2 = p2.getX();
                    int y2 = p2.getY();
                    
                    float xMidPt = (x1 + x2)/2.f;
                    float yMidPt = (y1 + y2)/2.f;
                    
                    int nWithinR1 = 0;
                    
                    Set<PairInt> neighbors2 = curveHelper.findNeighbors(x2, y2, 
                        points);
                    Set<PairInt> neighborsTot = new HashSet<PairInt>(neighbors1);
                    neighborsTot.addAll(neighbors2);
                    
                    for (PairInt pt : neighborsTot) {
                        
                        if (pt.equals(p1) || pt.equals(p2)) {
                            continue;
                        }
                        
                        float diffX = Math.abs(xMidPt - pt.getX());
                        float diffY = Math.abs(yMidPt - pt.getY());
                        double dist = Math.sqrt(diffX*diffX + diffY*diffY);

                        if (dist <= limit) {
                            nWithinR1++;
                            if (nWithinR1 > 1) {
                                break;
                            }
                        }
                    }
                       
                    if (nWithinR1 == 0) {
                        
                        // only p1 and p2 within radius of 1 of midpt
                        if (trimBeginning) {
                            curve.removeRange(0, (idxB - 1));
                        } else {
                            curve.removeRange((idxB + 1), curve.getN() - 1);
                        }
                        
                        trimEndpointSpursIfAny(curve, imageWidth, imageHeight);

                        return true;
                    }
                }
            }
        }
        
        return false;
    }
    
}
