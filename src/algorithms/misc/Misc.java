package algorithms.misc;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.PostLineThinnerCorrections;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import algorithms.util.ResourceFinder;
import java.io.BufferedWriter;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
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
    
    public static final int[] dx4 = new int[]{-1,  0, 1, 0};
    public static final int[] dy4 = new int[]{ 0, -1, 0, 1};
        
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
    
    public static String persistToFile(String fileName, Set<PairInt> points) 
        throws IOException {
        
        if (points == null) {
            throw new IllegalArgumentException("points cannot be null");
        }
                
        String outFilePath = ResourceFinder.findDirectory("bin") + "/" +
            fileName;
        
        FileOutputStream fs = null;
        ObjectOutputStream os = null;
        
        try {
            File file = new File(outFilePath);
            file.delete();
            file.createNewFile();

            fs = new FileOutputStream(file);
            os = new ObjectOutputStream(fs);
                        
            int count = 0;
            
            for (PairInt point : points) {
                
                os.writeInt(point.getX());
                os.writeInt(point.getY());
                
                if ((count % 10) == 0) {
                    os.flush();
                }
                
                count++;
            }

            os.flush();

        } finally {

            if (os != null) {
                os.close();
            }
            if (fs != null) {
                fs.close();
            }            
        }
        
        return outFilePath;
    }
    
    public static Set<PairInt> deserializeSetPairInt(String filePath) throws IOException {
                    
        FileInputStream fs = null;
        ObjectInputStream os = null;
        
        Set<PairInt> set = new HashSet<PairInt>();
        
        try {
            File file = new File(filePath);

            fs = new FileInputStream(file);
            os = new ObjectInputStream(fs);
                        
            while (true) {
                
                int x = os.readInt();
                int y = os.readInt();
                
                PairInt p = new PairInt(x, y);
                
                set.add(p);
            }
            
        } catch (EOFException e) {
            // expected
        } finally {

            if (os != null) {
                os.close();
            }
            if (fs != null) {
                fs.close();
            }            
        }
        
        return set;
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

    public static Set<Integer> findNeighborIndexes(PairIntArray edge, int idx) {
        
        Set<Integer> indexes = new HashSet<Integer>();
        
        int x = edge.getX(idx);
        int y = edge.getY(idx);
        
        return findNeighborIndexes(edge, x, y);
    }
    public static Set<Integer> findNeighborIndexes(PairIntArray edge, int x, int y) {
        
        Set<Integer> indexes = new HashSet<Integer>();
        
        for (int i = 0; i < edge.getN(); ++i) {
            int x2 = edge.getX(i);
            int y2 = edge.getY(i);
            if ((x2 == x) && (y2 == y)) {
                continue;
            }
            int xDiff = Math.abs(x2 - x);
            int yDiff = Math.abs(y2 - y);
            if ((xDiff < 2) && (yDiff < 2)) {
                indexes.add(Integer.valueOf(i));
            }
        }
        
        return indexes;
    }
   
}
