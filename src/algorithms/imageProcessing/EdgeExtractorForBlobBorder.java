package algorithms.imageProcessing;

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
               
        log.info("number of border points=" + nBorder + " nCavity=" + nCavity);
        
        if (discardWhenCavityIsSmallerThanBorder && (nCavity < (0.5*nBorder))) {
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
               
        log.info("number of border points=" + nBorder + " nCavity=" + nCavity);

if (debug) {        
Image img3 = new Image(imageWidth, imageHeight);
for (PairInt p : borderPixels) {
    img3.setRGB(p.getX(), p.getY(), 255, 0, 0);
}
MiscDebug.writeImageCopy(img3, "border_perimeter_" + MiscDebug.getCurrentTimeFormatted() + ".png");
}

        if (discardWhenCavityIsSmallerThanBorder && (nCavity < (0.5*nBorder))) {
            return null;
        }
        
        // ----- remove spurs, merge curves, re-order points to form a 
        //       single closed curve if possible -------
        
        ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
        lt.applyLineThinner(borderPixels, 0, imageWidth, 0, imageHeight);
        PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
        pltc.correctForExtCorner(borderPixels, imageWidth, imageHeight);
        
        if (borderPixels.isEmpty()) {
            return null;
        }
        
if (debug) {        
Image img3 = new Image(imageWidth, imageHeight);
for (PairInt p : borderPixels) {
    img3.setRGB(p.getX(), p.getY(), 255, 0, 0);
}
MiscDebug.writeImageCopy(img3, "border_before_spur_removal_" + MiscDebug.getCurrentTimeFormatted() + ".png");
}

        SpurRemover spurRm = new SpurRemover();
        
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
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        if (!curveHelper.isAdjacent(out, 0, out.getN() - 1)) {
            
            boolean altered = trimEndpointSpursIfAny(out, img.getWidth(), 
                img.getHeight());
            
            if (!curveHelper.isAdjacent(out, 0, out.getN() - 1)) {
                extractor.reorderEndpointsIfNeeded(out);
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
    
}
