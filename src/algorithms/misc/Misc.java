package algorithms.misc;

import algorithms.imageProcessing.features.CornerRegion;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.QuadInt;
import algorithms.util.ResourceFinder;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.security.Security;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

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

    public static PairIntArray convertWithoutOrder(Collection<PairInt> points) {
        PairIntArray out = new PairIntArray(points.size());
        for (PairInt p : points) {
            out.add(p.getX(), p.getY());
        }
        return out;
    }

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
    
    public static Map<PairInt, Integer> makePointIndexMap(PairIntArray edge) {
        
        if (edge == null) {
            throw new IllegalArgumentException("edge cannot be null");
        }
        
        Map<PairInt, Integer> map = new HashMap<PairInt, Integer>();
        
        for (int i = 0; i < edge.getN(); ++i) {
            map.put(new PairInt(edge.getX(i), edge.getY(i)), Integer.valueOf(i));
        }
        
        return map;
    }

    /**
     * get the values of Math.sin(theta) for theta from 0 to Math.PI in 1 degree
     * intervals.
     * @return
     */
    public static double[] getSineThetaForPI() {
        double d = 1. * Math.PI/180;
        double[] s = new double[180];
        double t = 0;
        for (int i = 0; i < 180; ++i) {
            s[i] = Math.sin(t);
            t += d;
        }
        return s;
    }
    /**
     * get the values of Math.cos(theta) for theta from 0 to Math.PI in 1 degree
     * intervals.
     * @return
     */
    public static double[] getCosineThetaForPI() {
        double d = 1. * Math.PI/180;
        double[] c = new double[180];
        double t = 0;
        for (int i = 0; i < 180; ++i) {
            c[i] = Math.cos(t);
            t += d;
        }
        return c;
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

    public static PairIntArray deserializePairIntArray(String filePath) throws IOException {

        FileInputStream fs = null;
        ObjectInputStream os = null;

        PairIntArray out = new PairIntArray();

        try {
            File file = new File(filePath);

            fs = new FileInputStream(file);
            os = new ObjectInputStream(fs);

            while (true) {
                int x = os.readInt();
                int y = os.readInt();
                out.add(x, y);
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

        return out;
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

    /**
     * create x and y offsets for the neighbor points within d pixel radius.
     * The result is a two-dimensional array of length 
     * 2*((2*d+1)^2 -1) with the numbers being x and y offsets
     * alternating.
     * Note that [(0, 0)] is not present and that the
     * offsets are ordered such that inner "rings" are completed,
     * then next inner ring at increased radius, etc.
     * @return
     */
    public static int[] createOrderedNeighborOffsets(int radiusFromCenter) {

        //TODO: consider changing to use one dimensional array

        int n = 2*radiusFromCenter + 1;
        n *= n;
        n--;
        n *= 2;
        int[] xyout = new int[n];

        int count = 0;
        for (int d = 1; d <= radiusFromCenter; ++d) {
            for (int x = -d; x <= d; ++x) {
                for (int y = -d; y <= d; ++y) {
                    // when both points have an absolute value smaller
                    // than d, it's an inner radius which has already
                    // been included so skip it
                    if (Math.abs(x) < d && Math.abs(y) < d) {
                        continue;
                    }
                    xyout[count] = x;
                    count++;
                    xyout[count] = y;
                    count++;
                }
            }
        }

        return xyout;
    }

    public static Map<Integer, Double> getCosineThetaMapForPI() {

        Map<Integer, Double> map = new HashMap<Integer, Double>();

        double d = 1. * Math.PI/180;

        double t = 0;

        for (int i = 0; i < 181; ++i) {

            double c = Math.cos(t);

            map.put(Integer.valueOf(i), Double.valueOf(c));

            t += d;
        }

        return map;
    }

    public static Map<Integer, Double> getSineThetaMapForPI() {

        Map<Integer, Double> map = new HashMap<Integer, Double>();

        double d = 1. * Math.PI/180;

        double t = 0;

        for (int i = 0; i < 181; ++i) {

            double c = Math.sin(t);

            map.put(Integer.valueOf(i), Double.valueOf(c));

            t += d;
        }

        return map;
    }

    public static Map<Integer, Double> getCosineThetaMapForTwoPI() {
        
        Map<Integer, Double> map = new HashMap<Integer, Double>();

        double d = 1. * Math.PI/180;

        double t = 0;

        for (int i = 0; i < 360; ++i) {

            double c = Math.cos(t);

            map.put(Integer.valueOf(i), Double.valueOf(c));

            t += d;
        }

        return map;
    }

    public static Map<Integer, Double> getSineThetaMapForTwoPI() {
        
        Map<Integer, Double> map = new HashMap<Integer, Double>();
        
        double d = 1. * Math.PI/180;

        double t = 0;

        for (int i = 0; i < 360; ++i) {

            double c = Math.sin(t);

            map.put(Integer.valueOf(i), Double.valueOf(c));

            t += d;
        }

        return map;
    }

    public static PairInt[] convert(CornerRegion[] c) {
        
        PairInt[] output = new PairInt[c.length];
        
        for (int i = 0; i < c.length; ++i) {
            
            CornerRegion cr = c[i];
            int x = cr.getX()[cr.getKMaxIdx()];
            int y = cr.getY()[cr.getKMaxIdx()];
            output[i] = new PairInt(x, y);
        }
        
        return output;
    }
    
    public static PairInt[] convert(List<CornerRegion> c) {
        
        PairInt[] output = new PairInt[c.size()];
        
        for (int i = 0; i < c.size(); ++i) {
            
            CornerRegion cr = c.get(i);
            int x = cr.getX()[cr.getKMaxIdx()];
            int y = cr.getY()[cr.getKMaxIdx()];
            output[i] = new PairInt(x, y);
        }
        
        return output;
    }
    
    public static int[][] populateNeighborOffsets(int radius) {
        
        int n0 = (2*((2*radius) + 1));  
        int n1 = n0 - 4;  
        int nTotal = n0 + n1;
        
        int[][] dxys = new int[2][];
        
        for (int i = 0; i < 2; ++i) {
            dxys[i] = new int[nTotal]; 
        }
        
        /*     -  -  -      radius = 1
               -  @  -
               -  -  - 
        
            #  #  #  #  #   radius = 2
               -  -  -  #
               -  @  -  #
               -  -  -  #
            #  #  #  #  #
        */
        
        int count = 0;
        // top
        for (int i = -radius; i <= radius; ++i) {
            dxys[0][count] = i;
            dxys[1][count] = radius;
            count++;
        }
        // right
        for (int j = (radius - 1); j >= -radius; --j) {
            dxys[0][count] = radius;
            dxys[1][count] = j;
            count++;
        }
        // bottom
        for (int i = (radius - 1); i >= -radius; --i) {
            dxys[0][count] = i;
            dxys[1][count] = -radius;
            count++;
        }
        // left
        for (int j = (radius - 1); j > -radius; --j) {
            dxys[0][count] = -radius;
            dxys[1][count] = j;
            count++;
        }
        
        return dxys;
    }

    public static <T extends Object> void reverse(List<T> list) {
        
        int n = list.size();
        
        if (n < 2) {
            return;
        }
                
        int end = n >> 1;
        // 0 1 2 3 4
        for (int i = 0; i < end; i++) {
            int idx2 = n - i - 1;
            T swap = list.get(i);
            list.set(i, list.get(idx2));
            list.set(idx2, swap);
        }
    }

    public static TObjectIntMap<PairInt> createPointIndexMap(
        PairIntArray p) {
        
        TObjectIntMap<PairInt> map = 
            new TObjectIntHashMap<PairInt>();
        
        for (int i = 0; i < p.getN(); ++i) {
            PairInt pt = new PairInt(p.getX(i),
                p.getY(i));
            map.put(pt, i);
        }
        
        return map;
    }
    
    /**
     * get an instance of SecureRandom, trying first
     * the algorithm SHA1PRNG, else the
     * default constructor.
     * @return 
     */
    public static SecureRandom getSecureRandom() {
        
        SecureRandom sr = null;
        
        try {
            sr = SecureRandom.getInstance("SHA1PRNG");
        } catch (NoSuchAlgorithmException ex) {
            Logger.getLogger(Misc.class.getName()).log(Level.SEVERE, null, ex);
            sr = new SecureRandom();
        }
        
        return sr;
    }
    
    /**
     * partition the space width x height into bins of size xLen, yLen
     * and return the results at a list of QuadInts holding 
     * a=minX, b=maxX, c=minY, d=maxY.
     * @param width
     * @param height
     * @param xLen
     * @param yLen
     * @return 
     */
    public static List<QuadInt> partitionSpace(int width, int height, 
        int xLen, int yLen) {
        
        List<QuadInt> output = new ArrayList<QuadInt>();
        
        int nXs = (int)Math.ceil((float)width/(float)xLen);
        int nYs = (int)Math.ceil((float)height/(float)yLen);
        
        for (int i = 0; i < nXs; ++i) {
            int startX = xLen * i;
            if (startX > (width - 1)) {
                break;
            }
            int stopX = startX + xLen;
            if (stopX > (width - 1)) {
                stopX = width - 1;
                // NOTE: if stopX - startX is small, consider merging it with previous
                // bin in terms of X
            }
            for (int j = 0; j < nYs; ++j) {
                int startY = yLen * i;
                if (startY > (height - 1)) {
                    break;
                }
                int stopY = startY + yLen;
                if (stopY > (height - 1)) {
                    stopY = height - 1;
                    // NOTE: if stopY - startY is small, consider merging it with previous
                    // bin
                }
                QuadInt q = new QuadInt(startX, stopX, startY, stopY);
                output.add(q);
            }
        }
        
        return output;
    }
    
    public static Set<PairInt> convertToCoords(GreyscaleImage img,
        TIntSet pixIdxs) {
        
        Set<PairInt> set = new HashSet<PairInt>();
        
        TIntIterator iter = pixIdxs.iterator();
        
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int col = img.getCol(pixIdx);
            int row = img.getRow(pixIdx);
            set.add(new PairInt(col, row));
        }
        
        return set;
    }
}
