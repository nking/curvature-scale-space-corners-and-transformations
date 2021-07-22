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
    
    /*
    https://www.microsoft.com/en-us/research/uploads/prod/2016/12/completecalibration.txt
    
        832.5 0.204494 832.53 303.959 206.585

        -0.228601 0.190353

        0.992759 -0.026319 0.117201
        0.0139247 0.994339 0.105341
        -0.11931 -0.102947 0.987505
        -3.84019 3.65164 12.791

        0.997397 -0.00482564 0.0719419
        0.0175608 0.983971 -0.17746
        -0.0699324 0.178262 0.981495
        -3.71693 3.76928 13.1974

        0.915213 -0.0356648 0.401389
        -0.00807547 0.994252 0.106756
        -0.402889 -0.100946 0.909665
        -2.94409 3.77653 14.2456

        0.986617 -0.0175461 -0.16211
        0.0337573 0.994634 0.0977953
        0.159524 -0.101959 0.981915
        -3.40697 3.6362 12.4551

        0.967585 -0.196899 -0.158144
        0.191542 0.980281 -0.0485827
        0.164592 0.0167167 0.98622
        -4.07238 3.21033 14.3441
    */
    
    public static double[][] getTransposedRotation(int imageIdx) {
        double[][] r = getRotation(imageIdx);
        return MatrixUtil.transpose(r);
    }
    
    public static double[][] getRotation(int imageIdx) {
        double[][] r = new double[3][];
        switch (imageIdx) {
            case 0:
                r[0] = new double[]{0.992759, -0.026319, 0.117201};
                r[1] = new double[]{0.0139247, 0.994339, 0.105341};
                r[2] = new double[]{-0.11931, -0.102947, 0.987505};
                break;
            case 1:
                r[0] = new double[]{0.997397, -0.00482564, 0.0719419};
                r[1] = new double[]{0.0175608, 0.983971, -0.17746};
                r[2] = new double[]{-0.0699324, 0.178262, 0.981495};        
                break;
            case 2:
                r[0] = new double[]{0.915213, -0.0356648, 0.401389};
                r[1] = new double[]{-0.00807547, 0.994252, 0.106756};
                r[2] = new double[]{-0.402889, -0.100946, 0.909665};        
                break;
            case 3:
                r[0] = new double[]{0.986617, -0.0175461, -0.16211};
                r[1] = new double[]{0.0337573, 0.994634, 0.0977953};
                r[2] = new double[]{0.159524, -0.101959, 0.981915};                
                break; 
            default:
                r[0] = new double[]{0.967585, -0.196899, -0.158144};
                r[1] = new double[]{0.191542, 0.980281, -0.0485827};
                r[2] = new double[]{0.164592, 0.0167167, 0.98622};
                break;
        }
        return r;
    }
    
    public static double[] getTranslation(int imageIdx) {
        double[] t = null;
        switch (imageIdx) {
            case 0:
                t = new double[]{-3.84019, 3.65164, 12.791};
                break;
            case 1:
                t = new double[]{-3.71693, 3.76928, 13.1974};
                break;
            case 2:
                t = new double[]{-2.94409, 3.77653, 14.2456};
                break;
            case 3:
                t = new double[]{-3.40697, 3.6362, 12.4551};
                break;
            default:
                t = new double[]{-4.07238, 3.21033, 14.3441};
                break;
        }        
        return t;
    }
    
}
