package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * algorithm to create a line of points between
 * two points.
 * 
 * https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
 * and
 * http://rosettacode.org/wiki/Bitmap/Bresenham%27s_line_algorithm#Java
 *
 * the code below is adapted from
 *    http://rosettacode.org/wiki/Bitmap/Bresenham%27s_line_algorithm#Java
    which has copyright
     GNU Free Documentation License 1.2 unless otherwise noted.

 * 
 */
public class BresenhamsLine {
    
    /**
     * calculate points in the line between (x1, y1)
     * and (x2, y2), inclusive, and add them to output
     * using Bresenham's algorithm.
     * The code is adapted from 
     * http://rosettacode.org/wiki/Bitmap/Bresenham%27s_line_algorithm#Java
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param output
     */
    public static void createLinePoints(int x1, int y1,
        int x2, int y2, Set<PairInt> output) {
    
        List<PairInt> output0 = new ArrayList<PairInt>();
        
        createLinePoints(x1, y1, x2, y2, output0);
        
        output.addAll(output0);
    }
    
    /**
     * calculate points in the line between (x1, y1)
     * and (x2, y2), inclusive, and add them to output
     * using Bresenham's algorithm.
     * The code is adapted from 
     * http://rosettacode.org/wiki/Bitmap/Bresenham%27s_line_algorithm#Java
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param output
     */
    public static void createLinePoints(int x1, int y1,
        int x2, int y2, List<PairInt> output) {
        
        // delta of exact value and rounded value of the dependant variable
        int d = 0;
 
        int dy = Math.abs(y2 - y1);
        int dx = Math.abs(x2 - x1);
 
        int dy2 = (dy << 1); // slope scaling factors to avoid floating
        int dx2 = (dx << 1); // point
 
        int ix = x1 < x2 ? 1 : -1; // increment direction
        int iy = y1 < y2 ? 1 : -1;
 
        if (dy <= dx) {
            for (;;) {
                output.add(new PairInt(x1, y1));
                if (x1 == x2)
                    break;
                x1 += ix;
                d += dy2;
                if (d > dx) {
                    y1 += iy;
                    d -= dx2;
                }
            }
        } else {
            for (;;) {
                output.add(new PairInt(x1, y1));
                if (y1 == y2)
                    break;
                y1 += iy;
                d += dx2;
                if (d > dy) {
                    x1 += ix;
                    d -= dy2;
                }
            }
        }
    }
}
