package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

/**
 * looks for diagonal zeros and widens the empty area by setting the opposing
 * diagonal to zero too.  It's meant to be used to increase segmentation.
 * 
 * @author nichole
 */
public class DiagZeroWidener extends AbstractLineThinner {
    
    protected boolean debug = false;
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    @Override
    public void applyFilter(GreyscaleImage input) {
                
        int w = input.getWidth();
        int h = input.getHeight();
        
        Set<PairInt> points = new HashSet<PairInt>();
        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {
                if (input.getValue(col, row) > 0) {
                    points.add(new PairInt(col, row));
                }
            }
        }
        
        applyZeroWidener(points, 0, w - 1, 0, h - 1);
        
        // unset points not in "points" and leave remaining values as is
        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {
                PairInt p = new PairInt(col, row);
                if (!points.contains(p)) {
                    input.setValue(col, row, 0);
                }
            }
        }
        
    }
    
    public void applyZeroWidener(Set<PairInt> points, int minX, int maxX,
        int minY, int maxY) {
    
        Pattern pattern = getPattern();
        
        applyZeroWidener(points, minX, maxX, minY, maxY, pattern);
        
        swapXDirection(pattern);
        
        applyZeroWidener(points, minX, maxX, minY, maxY, pattern);
        
    }
    
    private void applyZeroWidener(Set<PairInt> points, int minX, int maxX,
        int minY, int maxY, Pattern pattern) {
            
        for (int i = minX; i <= maxX; ++i) {
            for (int j = minY; j <= maxY; ++j) {
                applyPattern(i, j, points, pattern);
            }
        }
    }
    
    private void applyPattern(int x, int y, Set<PairInt> points,
        Pattern pattern) {

        for (PairInt p : pattern.zeroes) {
            PairInt p2 = new PairInt(x + p.getX(), y + p.getY());
            if (points.contains(p2)) {
                return;
            }
        }
        for (PairInt p : pattern.ones) {
            PairInt p2 = new PairInt(x + p.getX(), y + p.getY());
            if (!points.contains(p2)) {
                return;
            }
        }
        for (PairInt p : pattern.changeToZeroes) {
            PairInt p2 = new PairInt(x + p.getX(), y + p.getY());
            points.remove(p2);
        }
    }
    
    private Pattern getPattern() {
        
        Pattern pattern = new Pattern();
        
        pattern.ones = new HashSet<PairInt>();
        pattern.zeroes = new HashSet<PairInt>();
        pattern.changeToZeroes = new HashSet<PairInt>();
        
        /*
             #  #    -1
             0  #<    0
             #< 0     1
             #  #     2
        
             0  1
        */
        pattern.ones.add(new PairInt(0, 2)); pattern.ones.add(new PairInt(0, 1)); pattern.ones.add(new PairInt(0, -1));
        pattern.ones.add(new PairInt(1, 2)); pattern.ones.add(new PairInt(1, 0)); pattern.ones.add(new PairInt(1, -1));
        
        pattern.zeroes.add(new PairInt(0, 0));
        pattern.zeroes.add(new PairInt(1, 1));
        
        pattern.changeToZeroes.add(new PairInt(0, 1));
        pattern.changeToZeroes.add(new PairInt(1, 0));
        
        return pattern;
    }
    
    private void swapXDirection(Pattern pattern) {
        // ----- change the sign of x  -----
        for (PairInt p : pattern.zeroes) {
            p.setX(-1 * p.getX());
        }
        for (PairInt p : pattern.ones) {
            p.setX(-1 * p.getX());
        }
        for (PairInt p : pattern.changeToZeroes) {
            p.setX(-1 * p.getX());
        }
    }

    public static class Pattern {
        Set<PairInt> ones;
        Set<PairInt> zeroes;
        Set<PairInt> changeToZeroes;
    }
    
}
