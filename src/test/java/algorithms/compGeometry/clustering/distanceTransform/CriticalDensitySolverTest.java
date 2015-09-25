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

    public void testFindCriticalDensity() throws Exception {

        String[] fileNames = new String[]{
            "dt_ran_0_0.dat", "dt_ran_0_1.dat", "dt_ran_0_2.dat", 
            "dt_ran_1_0.dat", "dt_ran_1_1.dat", "dt_ran_1_2.dat", 
            "dt_ran_2_0.dat", "dt_ran_2_1.dat", "dt_ran_2_2.dat", 
        };
        
        //TODO: may need to revise these... used a different seed in random to generate:
        float[] r0s = new float[]{
            0.01f, 0.18f, 0.4f,
            0.01f, 0.18f, 0.4f,
            0.02f, 0.18f, 0.4f
        };
        float[] r1s = new float[]{
            0.05f, 0.22f, 0.5f,
            0.05f, 0.22f, 0.48f,
            0.04f, 0.22f, 0.48f
        };
        
        CriticalDensitySolver dSolver = new CriticalDensitySolver();
        dSolver.setToDebug();
        
        for (int i = 0; i < fileNames.length; ++i) {
            
            String fileName = fileNames[i];
            
            float[] values = readDataset(fileName);
            
            float critDens = dSolver.findCriticalDensity(values);
            
            float r0 = r0s[i];
            float r1 = r1s[i];

            assertTrue(critDens >= r0 && critDens <= r1);
        }
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
