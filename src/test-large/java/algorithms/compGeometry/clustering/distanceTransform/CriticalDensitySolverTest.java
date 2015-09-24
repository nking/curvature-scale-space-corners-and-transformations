package algorithms.compGeometry.clustering.distanceTransform;

import algorithms.util.ResourceFinder;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class CriticalDensitySolverTest extends TestCase {
    
    public CriticalDensitySolverTest(String testName) {
        super(testName);
    }

    public void testFindCriticalDensity() {
        
    }
    
    private float[] readDataset(String fileName) throws FileNotFoundException, IOException {
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        FileInputStream fs = null;
        ObjectInputStream os = null;

        float[] values = null;

        try {
            File file = new File(filePath);

            fs = new FileInputStream(file);
            os = new ObjectInputStream(fs);

            int n = os.readInt();
            values = new float[n];
            
            int count = 0;
            while (true) {
                float v = os.readFloat();
                values[count] = v;
                count++;
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

        return values;
    }
}
