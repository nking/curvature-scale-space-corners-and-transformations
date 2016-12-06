package algorithms.imageProcessing;

import algorithms.imageProcessing.matching.LinesFinder;
import algorithms.util.PairInt;
import java.util.List;
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
        // consider 2 different thresholds for 2 different
        //    minimum lengths
        //finder.overrideThreshold(0.05f);
        //finder.overrideThreshold(0.2f);
        finder.overrideMinimumLength(10); //ml=25,t=.2  ml=60,t=.2
        finder.setToRemoveBorderLines(imageWidth - 1, imageHeight - 1);
        finder.find(listOfContigousLabels);
        finder.groupWithinTolerance();
        this.finder = finder;
        //throw new UnsupportedOperationException("not yet implemented");
    }
    
    LinesFinder finder = null;
    
    public void debugDraw(Image img) {
        
        // draw lines onto img
        finder.debugDraw(img);
    }
}
