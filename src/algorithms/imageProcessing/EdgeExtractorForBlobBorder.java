package algorithms.imageProcessing;

import algorithms.compGeometry.PerimeterFinder;
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
MiscDebug.plotPoints(borderPixels, imageWidth, imageHeight, MiscDebug.getCurrentTimeFormatted());        
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
        
if (debug) {        
MiscDebug.plotPoints(borderPixels, imageWidth, imageHeight, MiscDebug.getCurrentTimeFormatted());        
}        
        SpurRemover spurRm = new SpurRemover();
        
        spurRm.remove(borderPixels, imageWidth, imageHeight);

if (debug) {        
MiscDebug.plotPoints(borderPixels, imageWidth, imageHeight, MiscDebug.getCurrentTimeFormatted());        
}        
        //xMin, xMax, yMin, yMax
        int[] minMaxXY = MiscMath.findMinMaxXY(borderPixels);
      
        int xOffset = minMaxXY[0] - 2;
        int yOffset = minMaxXY[2] - 2;
        int w = minMaxXY[1] - minMaxXY[0] + 4;
        int h = minMaxXY[3] - minMaxXY[2] + 4;
        
        GreyscaleImage img = new GreyscaleImage(w, h);
        for (PairInt p : borderPixels) {
            int x = p.getX() - xOffset;
            int y = p.getY() - yOffset;
            img.setValue(x, y, 1);
        }
        
        EdgeExtractorWithJunctions extractor = new EdgeExtractorWithJunctions(img);
        extractor.overrideMaxNumberIterationsJunctionSplice(10);
        //List<PairIntArray> output = extractor.findEdges();
        PairIntArray out = extractor.findAsSingleClosedEdge();
       
        if (out == null) {
            return null;
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
    /*if (i > 0) {
        x += xOffset;
        y += yOffset;
    }*/
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
    
}
