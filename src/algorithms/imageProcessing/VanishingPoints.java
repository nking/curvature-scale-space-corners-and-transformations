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
    
    public void find(List<Set<PairInt>> listOfContigousLabels,
        int imageWidth, int imageHeight) {
     
        /*
        NOTE: may need to allow additional logic related to buildings
        and the logic cpuld be used either by an option set outside
        or a method here to look for buildings.
        
        for buidings, would need to allow for shorter lines
        possibly for projected dimensions
        and maybe look for those lines connected by junctions
        (and same or adjacent sementation).
        If more than one building had same pose and proximity,
        the vanishing lines should have similar perpendiculat
        polar angles...there is much additional logic that
        could be applied to the found lines.
        */
        
        LinesFinder finder = new LinesFinder();
        if (debug) {
            finder.setToDebug();
        }
        
        finder.setToRemoveBorderLines(imageWidth - 1, imageHeight - 1);
        finder.find(listOfContigousLabels);
        finder.groupWithinTolerance();
        this.finder = finder;
        //throw new UnsupportedOperationException("not yet implemented");
    }
    
    LinesFinder finder = null;
    
    public void correctLinesWithGradient(GreyscaleImage img) {
        finder.correctLinesWithGradient(img);
    }
    
    public void debugDraw(Image img) {
        // draw lines onto img
        finder.debugDraw(img);
    }
}
