package algorithms.imageProcessing.optimization.segmentation;

import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.util.PairInt;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class BerkeleySegmentationFileReader {
    
    public BerkeleySegmentationFileReader() {
    }
    
    /**
     * read the file and calculate the segmented region centroids and return double
     * array of int[0][segmentNumber] = x centroid, 
     * int[1][segmentNumber] = y centroid, int[2][segmentNumber] = number of 
     * points in segmented region.
     *  
     * @param filePath
     * @return
     * @throws IOException 
     */
    public int[][] readCentroids(String filePath) throws IOException {
        
        List<Set<PairInt>> sets = readFile(filePath);
        
        int n = sets.size();
        
        int[][] xyCenN = new int[3][n];
        
        for (int i = 0; i < 3; ++i) {
            xyCenN[i] = new int[n];
        }
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        for (int i = 0; i < n; ++i) {
            
            Set<PairInt> set = sets.get(i);
            
            double[] xyCen = curveHelper.calculateXYCentroids(set);
            
            xyCenN[0][i] = (int)Math.round(xyCen[0]);
            
            xyCenN[1][i] = (int)Math.round(xyCen[1]);
            
            xyCenN[2][i] = set.size();
        }
        
        return xyCenN;
    }
    
    public List<Set<PairInt>> readFile(String filePath) throws IOException {
        
        List<Set<PairInt>> output = new ArrayList<Set<PairInt>>();
        
        FileReader reader = null;
        BufferedReader in = null;

        try {
            reader = new FileReader(filePath);
            in = new BufferedReader(reader);

            String line = "";
            while (!line.contains("segments")) {
                line = in.readLine();
            }
            String nStr = line.split(" ")[1];
            int nSegments = Integer.valueOf(nStr).intValue();
            for (int i = 0; i < nSegments; ++i) {
                output.add(new HashSet<PairInt>());
            }
            
            while (!line.contains("data")) {
                line = in.readLine();
            }
            line = in.readLine();
            
            while (line != null) {
                
                String[] items = line.split(" ");
                assert(items.length == 4);
                
                int segmentNumber = Integer.valueOf(items[0]).intValue();
                
                int row = Integer.valueOf(items[1]).intValue();
                
                int colStart = Integer.valueOf(items[2]).intValue();
                
                int colStop = Integer.valueOf(items[3]).intValue();
                
                Set<PairInt> set = output.get(segmentNumber);
                
                for (int col = colStart; col <= colStop; ++col) {
                    set.add(new PairInt(col, row));
                }
                
                line = in.readLine();
            }
            
            return output;
            
        } finally {
            if (reader != null) {
                reader.close();
            }
            if (in != null) {
                in.close();
            }
        }        
    }
    
}
