package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

/**
 * convenience methods for camera test data Zhang 1998.
 * see testresources/zhang1998/README.txt
 * 
 * @author nichole
 */
public class Zhang98Data {
    
    public static String ZHANGDIR = "zhang1998";
    protected final static String sep = System.getProperty("file.separator");

    /**
     * get the 256 corners of the checkerboard in world reference frame.
     * the coordinates are in inches.
     * @return double array of size 3 X 256
     * @throws IOException 
     */
    public static double[][] getFeatureWCS() throws IOException {
        
        String path = ResourceFinder.findTestResourcesDirectory() + sep + ZHANGDIR
            + sep + "Model.txt";
        
        File f = new File(path);
        if (!f.exists()) {
            throw new IOException("could not find file at " + path);
        }
        
        BufferedReader in = new BufferedReader(new FileReader(new File(path)));
        
        double[][] wcs = readFile(in);
        
        in.close();
        
        return wcs;
    }
     
    /**
     * get the 256 checkerboard corners for the 5 images in pixel coordinates.
     * @return a double array of size 3 X 256*5
     */
    public static double[][] getFeaturesInAllImages() throws IOException {
        
        int nFeatures = 256;
        double[][] uv = MatrixUtil.zeros(3, nFeatures*5);
        double[][] uvi;
        int i, j;
        for(i = 0; i < 5; ++i) {
            uvi = getFeaturesImage(i+1);
            assert(nFeatures == uvi[0].length);
            for (j = 0; j < 3; ++j) {
                System.arraycopy(uvi[j], 0, uv[j], nFeatures*i, nFeatures);
            }
        }
        return uv;
    }
    
    // data1.txt through data5.txt
    // size is 3 X 256
    public static double[][] getFeaturesImage(int idx) throws IOException {
        if (idx < 0 || idx > 5) {
            throw new IllegalArgumentException("idx must be 1 through 5, inclusive");
        }
        
        String path = ResourceFinder.findTestResourcesDirectory() + sep + ZHANGDIR
            + sep + "data" + String.valueOf(idx) + ".txt";
        
        File f = new File(path);
        if (!f.exists()) {
            throw new IOException("could not find file at " + path);
        }
        
        BufferedReader in = new BufferedReader(new FileReader(new File(path)));
        
        double[][] xy = readFile(in);
        
        in.close();
        
        return xy;
    }
    
    private static double[][] readFile(BufferedReader in) throws IOException {
        
        int nFeatures = 256;
        
        // 3 dimensions, 4*64 features 
        double[][] a = MatrixUtil.zeros(3, nFeatures);
        Arrays.fill(a[2], 1.);
        
        // reading 64 lines. each line is 4 corners of 1 square (= 8 real numbers).
        int c = 0;
        int j;
        String[] coords;
        String line = in.readLine();
        while (line != null) {
            coords = line.split("\\s+");
            assert(coords.length == 8);
            
            for (j = 0; j < 4; ++j) {
                a[0][c] = Double.parseDouble(coords[2*j]);
                a[1][c] = Double.parseDouble(coords[2*j + 1]);
                c++;
            }
            
            line = in.readLine();
        }
        in.close();
        
        return a;
    }
    
}
