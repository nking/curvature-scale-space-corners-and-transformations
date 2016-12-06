package algorithms.imageProcessing;

import algorithms.compGeometry.HoughTransform;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.matching.LineFinder;
import algorithms.imageProcessing.matching.LinesFinder;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * 
 * NOT READY FOR USE.  first tests find safe segments
 * of lines, but the entire line(s) over a labelled cell
 * could be completed from the found segment.
 * 
 * a class to hold various methods for determining vanishing points
 * and to hold the resulting vanishing points.
 * 
 * @author nichole
 */
public class VanishingPoints {
        
    private boolean debug = false;
    
    public void setToDebug() {
        debug = true;
    }
    
    /**
     * points and orientations to use for calculating vanishing lines and points.
     * NOTE that the points should probably be edge points or
     * the boundaries of segmentation contiguous labels as filtered edges.
     */
    public void find(List<Set<PairInt>> listOfContigousLabels,
        int imageWidth, int imageHeight) {
        
        LinesFinder finder = new LinesFinder();
        if (debug) {
            finder.setToDebug();
        }
        //finder.overrideThreshold(0.085f);
        //finder.overrideThreshold(0.5f);
        finder.overrideMinimumLength(15);
        finder.setToRemoveBorderLines(imageWidth - 1, imageHeight - 1);
        finder.find(listOfContigousLabels);
        finder.groupWithinTolerance();

        finder.debugPrintTRStats();
        this.finder = finder;
        //throw new UnsupportedOperationException("not yet implemented");
    }
    
    LinesFinder finder = null;
    
    public void debugDraw(Image img) {
        
        // draw lines onto img
        finder.debugDraw(img);
    }
}
