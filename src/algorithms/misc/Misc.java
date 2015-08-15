package algorithms.misc;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.PostLineThinnerCorrections;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import java.util.HashSet;
import java.util.Set;
import java.util.Stack;

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
     * Note that the offset of (0,0) is in the middle of the arrays.
     * @param d the half radius of square of offsets, beginning at
     * (-d,-d), (-d, -d+1),... to make a (2d+1)^2 two dimensional array 
     * of offsets.
     * @return 
     */
    public static float[][] createNeighborOffsets(int radiusFromCenter) {
        
        //TODO: consider changing to use one dimensional array
        
        int n = 2*radiusFromCenter + 1;
        
        float[][] xyout = new float[n*n][];
                
        int count = 0;
        for (int x = -radiusFromCenter; x <= radiusFromCenter; ++x) {
            for (int y = -radiusFromCenter; y <= radiusFromCenter; ++y) {
                
                xyout[count] = new float[]{x, y};
                
                count++;
            }
        }
    
        return xyout;
    }
    
    /**
     * given a set of unordered contiguous sequential points, order them and 
     * return a closed edge (instance of PairIntArrayWithColor with color=1).
     * Note, this method removes junctions, that is, if there are more than two
     * adjacent points for a point, the line is simplified to having only
     * two while keeping points connected.
     * @param points
     * @return 
     */
    public static PairIntArrayWithColor orderSequentially(Set<PairInt> points, 
        int imageWidth, int imageHeight) {
                
        // line thinner has some rules to keep from disconnecting from image 
        // boundaries, so needs to know image size and give the line thinner
        // a larger value
        PostLineThinnerCorrections lt = new PostLineThinnerCorrections();
        lt.correctForArtifacts(points, imageWidth, imageHeight);
        
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.removeSpurs(points, imageWidth, imageHeight);
        
        /* any remaining regions in the line that are not strictly single
        pixel width including diagonal as adjacent, will
        throw exceptions here to find while testing
        and correct for them.
        */
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        // DFS search for sequential neighbors.
        Stack<PairInt> stack = new Stack<PairInt>();
        stack.addAll(points);
        
        Set<PairInt> added = new HashSet<PairInt>();
        
        PairIntArrayWithColor output = new PairIntArrayWithColor(points.size());
      
        output.add(stack.peek().getX(), stack.peek().getY());
        added.add(stack.peek());
        
        // > O(N) and << O(N^2)
        while (!stack.isEmpty()) {
            
            PairInt uNode = stack.pop();
            
            int uX = uNode.getX();
            int uY = uNode.getY();
            
            Set<PairInt> neighbors = new HashSet<PairInt>();
                                    
            for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                
                int vX = dxs[nIdx] + uX;
                int vY = dys[nIdx] + uY;
                
                PairInt vNode = new PairInt(vX, vY);
                
                if (!points.contains(vNode)) {
                    continue;
                }

                if (added.contains(vNode)) {
                    continue;
                }

                neighbors.add(vNode);
            }   
            
            if ((neighbors.size() > 1) && (added.size() > 1)) {
               
                throw new IllegalStateException("points are not contiguous or" +
                    " there are junctions handled incorrectly");
            } else if (neighbors.size() == 1) {
                for (PairInt vNode : neighbors) {
                    output.add(vNode.getX(), vNode.getY());
                    added.add(vNode);
                    stack.add(vNode);
                }
            }
        }
        
        return output;
    }
}
