package algorithms.misc;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
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
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
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

    public static void reverse(float[] a) {
        
        int n = a.length;
        
        if (n < 2) {
            return;
        }
                
        int end = n >> 1;
        // 0 1 2 3 4
        for (int i = 0; i < end; i++) {
            int idx2 = n - i - 1;
            float swap = a[i];
            a[i] = a[idx2];
            a[idx2] = swap;
        }
        
    }
    
    /**
     * assuming that the coefficients are ordered from highest order to
     * lowest, e.g. coeff[0] * x^2 + coeff[1] * x coeff[2],
     * apply them to x and resturn the model.
     * @param coeffs
     * @param x
     * @return 
     */
    public static float[] generate(float[] coeffs, float[] x) {
        
        float[] y = new float[x.length];
                
        for (int i = 0; i < x.length; ++i) {
            float x2 = 1;
            for (int j = coeffs.length - 1; j > -1; j--) {
                float c = coeffs[j];
                y[i] += (c * x2);
                x2 *= x[i];
            }
        }
        
        return y;
    }
    
    /**
     * assuming that the coefficients are ordered from highest order to
     * lowest, e.g. coeff[0] * x^2 + coeff[1] * x coeff[2],
     * apply them to x and return the model.
     * @param coeffs
     * @param x
     * @return 
     */
    public static double[] generate(double[] coeffs, double[] x) {
        
        double[] y = new double[x.length];
        
        for (int i = 0; i < x.length; ++i) {
            double x2 = 1;
            for (int j = coeffs.length - 1; j > -1; j--) {
                double c = coeffs[j];
                y[i] += (c * x2);
                x2 *= x[i];
            }
        }
        
        return y;
    }
}
