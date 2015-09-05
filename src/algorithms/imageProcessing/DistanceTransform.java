package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.Set;

/**
 * Euclidean square distance transforms are computed for every zero pixel relative to a non-zero point.
 *
 * The implementations are euclidean distance by default unless stated otherwise.
 * Chessboard implementations may be added.  All adjacent pixels, including 
 * diagonal have a distance of '1' for chessboard implementation.
 * 
 * @author nichole
 */
public class DistanceTransform {
    
    /**
     * The square of Euclidean distances are computed for every zero pixel 
     * relative to the nearest non-zero pixels for two-dimensional input.
     * <pre>
     * original data:
        row 0:  1 1 1 1 1 1 1 1 1 
        row 1:  1 1 1 1 1 1 1 1 1 
        row 2:  1 1 1 0 0 0 1 1 1 
        row 3:  1 1 0 0 0 0 0 1 1 
        row 4:  1 1 0 0 0 0 0 1 1 
        row 5:  1 1 0 0 0 0 0 1 1 
        row 6:  1 1 1 0 0 0 1 1 1 
        row 7:  1 1 1 1 1 1 1 0 1 
        row 8:  1 1 1 1 1 1 0 1 1 

        distance transform:
        row 0:  0 0 0 0 0 0 0 0 0 
        row 1:  0 0 0 0 0 0 0 0 0 
        row 2:  0 0 0 1 1 1 0 0 0 
        row 3:  0 0 1 2 4 2 1 0 0 
        row 4:  0 0 1 4 8 4 1 0 0 
        row 5:  0 0 1 2 4 2 1 0 0 
        row 6:  0 0 0 1 1 1 0 0 0 
        row 7:  0 0 0 0 0 0 0 1 0 
        row 8:  0 0 0 0 0 0 1 0 0
     * </pre>
     * 
     * an O(N) runtime complexity algorithm for computing the distance transform
     * by Meijster, Roerdink, and Hesselink 2000, implemented from their pseudocode.
     * "Mathematical Morphology and its Applications to Image and Signal Processing",
       in Volume 18 of the series Computational Imaging and Vision, 200, pp 331-340
       * ISBN 978-0-306-47025-7
     *  * https://www.rug.nl/research/portal/files/3059926/2002CompImagVisMeijster.pdf
     * 
     * @param points
     * @param width
     * @param height
     * @return 
     */
    public int[][] applyMeijsterEtAl(Set<PairInt> points, final int width, final int height) {
        
        int[][] g = new int[width][height];
        for (int i = 0; i < width; ++i) {
            g[i] = new int[height];
        }
        
        applyPhase1(points, g, width, height);
        
        int[][] dt = applyPhase2(g, width, height);
        
        return dt;
    }

    private void applyPhase1(Set<PairInt> blob, int[][] g, final int width, 
        final int height) {
       
        int inf = 1 << 22;//(int)Math.ceil(width*width + height*height);
        
        for (int x = 0; x < width; ++x) {
            
            // scan 1
            if (blob.contains(new PairInt(x, 0))) {
                g[x][0] = 0;
            } else {
                g[x][0] = inf;
            }
        
            for (int y = 1; y < height; ++y) {
                if (blob.contains(new PairInt(x, y))) {
                    g[x][y] = 0;
                } else {
                    g[x][y] = g[x][y - 1] + 1;
                }
            }
            
            // scan 2
            for (int y = height - 2; y > -1; --y) {
                if (g[x][y + 1] < g[x][y]) {
                    g[x][y] = g[x][y + 1] + 1;
                }
            }
        }
    }
    
    private int[][] applyPhase2(int[][] g, final int width, 
        final int height) {
                
        int[][] dt = new int[width][height];
        for (int i = 0; i < width; ++i) {
            dt[i] = new int[height];
        }
                
        int[] s = new int[width + 1];
        int[] t = new int[width + 1];
        
        for (int y = 0; y < height; ++y) {
            
            int q = 0;
            s[0] = 0;
            t[0] = 0;
            
            for (int u = 1; u < width; ++u) {
                
                // g(i) = G(i, y)
                
                while ((q >= 0) && (f(t[q], s[q], g[s[q]][y]) > f(t[q], u, g[u][y]))) {
                    q--;
                }
                if (q < 0) {
                    q = 0;
                    s[0] = u;
                } else {
                    int w = 1 + sep(s[q], u, g[s[q]][y], g[u][y]);
                    if (w < width) {
                        q++;
                        s[q] = u;
                        t[q] = w;
                    }
                }
            }
            for (int u = (width - 1); u > -1; --u) {
                dt[u][y] = f(u, s[q], g[s[q]][y]);
                if (u == t[q]) {
                    q--;
                }
            }
        }
        
        return dt;
    }

    private int f(int x, int i, int gi) {
                
        int xi = x - i;
        
        int f = (xi * xi) + (gi * gi);
        
        return f;
    }

    private int sep(int i, int u, int gi, int gu) {
        
        int sep = ((u * u) - (i * i) + (gu * gu) - (gi * gi))/(2 * (u - i));
        
        return sep;
    }
}
