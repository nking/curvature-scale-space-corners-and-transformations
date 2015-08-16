package algorithms.imageProcessing;

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
    
    public PairIntArray extractAndOrderAsEdge(Set<PairInt> contiguousPoints,
        int imageWidth, int imageHeight) {
        
        // TODO: improve this by refactoring the edge extractors one day.
        // for now, using the existing by making a small image large enough to
        // hold these points, then will change the offsets
        
        //xMin, xMax, yMin, yMax
        int[] minMaxXY = MiscMath.findMinMaxXY(contiguousPoints);
        
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
        for (PairInt p : contiguousPoints) {
            int x = p.getX() - xOffset;
            int y = p.getY() - yOffset;
            img.setValue(x, y, 1);
        }
        
        EdgeExtractorWithJunctions extractor = new EdgeExtractorWithJunctions(img);
        List<PairIntArray> output = extractor.findEdges();
                
        if (output.isEmpty()) {
            
            return null;
        }
        
        if (output.size() > 1) {
            
            Collections.sort(output, new PairIntArrayDescendingComparator());
        }
        
        // add the shifts back
        PairIntArray out = output.get(0);
        
        for (int i = 0; i < out.getN(); ++i) {
            int x = out.getX(i) + xOffset;
            int y = out.getY(i) + yOffset;
            out.set(i, x, y);
        }
        
        return out;
    }
    
}
