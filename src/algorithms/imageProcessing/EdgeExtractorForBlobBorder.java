package algorithms.imageProcessing;

import algorithms.compGeometry.PerimeterFinder;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArrayComparator;
import algorithms.util.PairIntArrayDescendingComparator;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * 
 * @see AbstractEdgeExtractor

 * @author nichole
 */
public class EdgeExtractorForBlobBorder {
            
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public EdgeExtractorForBlobBorder() {
                    
    }
    
    /**
     * given the set of contiguous points, find the perimeter of them and order
     * the points into a closed single pixel width curve.
     * @param contiguousPoints
     * @param imageWidth
     * @param imageHeight
     * @param discardWhenCavityIsSmallerThanBorder
     * @return 
     */
    public PairIntArray extractAndOrderTheBorder(Set<PairInt> contiguousPoints,
        int imageWidth, int imageHeight, boolean discardWhenCavityIsSmallerThanBorder) {
        
        if (contiguousPoints == null) {
            return null;
        }
        
        // TODO: make a method for PerimeterFinder which takes the gap filled
        // rowColRanges and follows the outside of it to make an ordered
        // list of border points in format PairIntArray.
        // not easy when there are spurs so may need to note that single 
        // pixel width extensions will be trimmed for a sequential bordering
        // curve.
        
        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        
        int imageMaxColumn = imageWidth - 1;
        int imageMaxRow = imageHeight - 1;
       
        int[] rowMinMax = new int[2];
        
        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
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
               
        log.info("number of border points=" + borderPixels.size() + " nCavity="
            + nCavity);
        
        if (discardWhenCavityIsSmallerThanBorder && (nCavity < (borderPixels.size()/2))) {
            return null;
        }
        
GreyscaleImage img2 = new GreyscaleImage(imageWidth, imageHeight);
for (PairInt p : borderPixels) {
int x = p.getX();
int y = p.getY();
img2.setValue(x, y, 1);
}
MiscDebug.plotCorners(img2, borderPixels, "___" + MiscDebug.getCurrentTimeFormatted(), 0);
 
        SpurRemover spurRm = new SpurRemover();
        
        spurRm.remove(borderPixels, imageWidth, imageHeight);
        
MiscDebug.plotPoints(contiguousPoints, imageWidth, imageHeight, MiscDebug.getCurrentTimeFormatted());        
MiscDebug.plotPoints(borderPixels, imageWidth, imageHeight, MiscDebug.getCurrentTimeFormatted());

        //xMin, xMax, yMin, yMax
        int[] minMaxXY = MiscMath.findMinMaxXY(borderPixels);
       
        /*
           |     @   @
           |      @@@
           2            
        subtract (xMin - 2)
        subtract (yMin - 2)
        */
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
        /*
        EdgeExtractor extractor = new EdgeExtractor(img);
        List<PairIntArray> output = extractor.findEdges();
        */
            
        /*
        if (output.isEmpty()) {
            return null;
        }
        if (output.size() > 1) {            
            Collections.sort(output, new PairIntArrayDescendingComparator());
        }
        */
        
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
/*    
Image img2 = new Image(imageWidth, imageHeight);
for (int i = 0; i < output.size(); ++i) {
    PairIntArray pa = output.get(i);
    for (int j = 0; j < pa.getN(); ++j) {
        int x = pa.getX(j);
        int y = pa.getY(j);
        if (i > 0) {
            x += xOffset;
            y += yOffset;
        }
        if (i == 0) {
            if (j == 0 || (j == (pa.getN() - 1))) {
                ImageIOHelper.addPointToImage(x, y, img2, 0, 200, 100, 0);
            } else {
                ImageIOHelper.addPointToImage(x, y, img2, 0, 255, 0, 0);
            }
        } else if (i == 1) {
            ImageIOHelper.addPointToImage(x, y, img2, 0, 0, 255, 0);
        } else {
            ImageIOHelper.addPointToImage(x, y, img2, 0, 0, 0, 255);
        }
    }
}
MiscDebug.writeImageCopy(img2, "output_" + MiscDebug.getCurrentTimeFormatted() + ".png");
*/
        
MiscDebug.plotPoints(out, imageWidth, imageHeight, MiscDebug.getCurrentTimeFormatted());

        return out;
    }
    
}
