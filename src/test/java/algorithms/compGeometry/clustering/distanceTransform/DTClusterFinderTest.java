package algorithms.compGeometry.clustering.distanceTransform;

import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class DTClusterFinderTest extends TestCase {
    
    public DTClusterFinderTest(String testName) {
        super(testName);
    }

    public void testFindRanGenClusters() {
        
    }
    
    public void testNoClusters() {
        
    }
    
    public void testFindClustersOtherData() {
        
    }
    
    public static String writeDataset(float[] values, String fileName) throws 
        IOException {
        
        if (values == null) {
            throw new IllegalArgumentException("values cannot be null");
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
            
            os.writeInt(values.length);

            int count = 0;

            for (float v : values) {

                os.writeFloat(v);

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
}
