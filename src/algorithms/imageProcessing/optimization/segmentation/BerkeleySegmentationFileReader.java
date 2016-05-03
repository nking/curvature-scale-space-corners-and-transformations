package algorithms.imageProcessing.optimization.segmentation;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
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
