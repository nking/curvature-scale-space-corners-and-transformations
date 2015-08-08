package algorithms.misc;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.HashSet;
import java.util.Set;

/**
 * miscellaneous boiler plate code
 * 
 * @author nichole
 */
public class Misc {
    
    public static final int[] dx8 = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
    public static final int[] dy8 = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
    public static Set<PairInt> convert(PairIntArray points) {
        
        if (points == null) {
            throw new IllegalArgumentException("points cannot be null");
        }
        
        Set<PairInt> out = new HashSet<PairInt>();
        
        for (int i = 0; i < points.getN(); ++i) {
            out.add(new PairInt(points.getX(i), points.getY(i)));
        }
        
        return out;
    }
    
    public static int calculateSumOfEightNeighbors(GreyscaleImage img, int x, int y) {
        
        int sum = 0;
        
        for (int i = 0; i < dx8.length; ++i) {
            int x1 = x + dx8[i];
            int y1 = y + dy8[i];
            if (x1 < 0 || y1 < 0 || (x1 > (img.getWidth() - 1)) || 
                (y1 > (img.getHeight() - 1))) {
                continue;
            }
            sum += img.getValue(x, y);
        }
        
        return sum;
    }
    
    /**
     * create x and y offsets for the neighbor points within d pixel radius.
     * The result is a two-dimensional array of length (2*d+1)^2 with the first
     * dimension being the x array and the 2nd dimension being the y array.
     * Note that the offset of (0,0) is the first given.
     * @param d the half radius of square of offsets, beginning at (0,0)
     * then (-d,-d), (-d, -d+1),... to make a dXd two dimensional array 
     * of offsets.
     * @return 
     */
    public static float[][] createNeighborOffsets(int radiusFromCenter) {
        
        //TODO: consider changing to use one dimensional array
        
        int n = 2*radiusFromCenter + 1;
        
        float[][] xyout = new float[n*n][];
        xyout[0] = new float[2];
                
        int count = 1;
        for (int x = -radiusFromCenter; x <= radiusFromCenter; ++x) {
            for (int y = -radiusFromCenter; y <= radiusFromCenter; ++y) {
                if (x == 0 && y == 0) {
                    continue;
                }
                
                xyout[count] = new float[]{x, y};
                
                count++;
            }
        }
    
        return xyout;
    }
    
}
