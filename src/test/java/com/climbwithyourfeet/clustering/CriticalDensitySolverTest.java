package com.climbwithyourfeet.clustering;

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
    
    /**
     *
     * @param testName
     */
    public CriticalDensitySolverTest(String testName) {
        super(testName);
    }

    /**
     *
     * @throws Exception
     */
    public void testFindCriticalDensity_ran() throws Exception {

        String[] fileNames = new String[]{
            "dt_ran_0_0.dat", "dt_ran_0_1.dat", "dt_ran_0_2.dat", 
            "dt_ran_1_0.dat", "dt_ran_1_1.dat", "dt_ran_1_2.dat", 
            "dt_ran_2_0.dat", "dt_ran_2_1.dat", "dt_ran_2_2.dat", 
        };
        
        float[] r0s = new float[]{
            0.01f, 0.18f, 0.4f,
            0.01f, 0.18f, 0.4f,
            0.02f, 0.18f, 0.4f
        };
        float[] r1s = new float[]{
            0.05f, 0.285f, 0.5f,
            0.05f, 0.28f, 0.5f,
            0.04f, 0.28f, 0.5f
        };
        
        CriticalDensitySolver dSolver = new CriticalDensitySolver();
        dSolver.setToDebug();
        
        for (int i = 0; i < fileNames.length; ++i) {
            
            String fileName = fileNames[i];
            
            float[] values = readDataset(fileName);
            
            float critDens = dSolver.findCriticalDensity(values);
            
            float r0 = r0s[i];
            float r1 = r1s[i];
            
            System.out.println("i=" + i + " ro-" + r0 + " r1=" + r1 + " critDens=" + critDens);
            
            /* 
            //need to fix the random seed for assertions
            TODO: change to run 2 passes, one with assertions and fixed seed
            and one without fixed seed
            if (i == 0) {
                assertTrue(critDens >= r0 && (critDens <= 0.17));
            } else {
                assertTrue(critDens >= r0 && (critDens <= (r1 + 0.1f*r1)));
            }
            */
        }
    }
    
    /**
     *
     * @throws Exception
     */
    public void testFindCriticalDensity_other() throws Exception {

        String[] fileNames = new String[]{
            "dt_other_0.dat", "dt_other_1.dat", "dt_other_2.dat", 
            "dt_other_3.dat", "dt_other_4.dat", "dt_other_5.dat", 
            "dt_other_6.dat", "dt_other_7.dat", "dt_dbscan.dat",
            //"dt_other_8.dat",  
        };
        
        float[] r0s = new float[]{
            0.35f, 0.32f, 0.25f, 
            0.25f, 0.39f, 0.3f,
            0.32f, 0, 0.1f,
            0,
        };
        float[] r1s = new float[]{
            0.75f, 1.1f, 0.65f, 
            0.7f, 1.75f, 1.3f,
            0.5f, 0.95f, 0.195f,
            0.015f,
        };
        
        CriticalDensitySolver dSolver = new CriticalDensitySolver();
        dSolver.setToDebug();
        
        for (int i = 0; i < fileNames.length; ++i) {
            
            String fileName = fileNames[i];
            
            float[] values = readDataset(fileName);
            
            float critDens = dSolver.findCriticalDensity(values);
            
            float r0 = r0s[i];
            float r1 = r1s[i];
                        
            System.out.println("i=" + i + " ro-" + r0 + " r1=" + r1 + " critDens=" + critDens);

            if (i != 6) {
                assertTrue(critDens >= r0 && critDens <= r1);
            }
        }
    }
    
    /**
     *
     * @throws Exception
     */
    public void testFindCriticalDensity_no_clusters() throws Exception {

        String[] fileNames = new String[]{
            "dt_no_clusters_0.dat", "dt_no_clusters_1.dat", "dt_no_clusters_2.dat", 
            "dt_no_clusters_3.dat", "dt_no_clusters_4.dat", "dt_no_clusters_5.dat", 
            "dt_no_clusters_6.dat", "dt_no_clusters_7.dat", "dt_no_clusters_8.dat", 
        };
        
        // TODO: these may need to be revised a I used a different seed during generation
        float[] r0s = new float[]{
            0.018f, 0.045f, 0.07f, 
            0.15f, 0.18f, 0.22f,
            0.22f, 0.25f, 0.25f
        };
        float[] r1s = new float[]{
            Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY,
            Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY,
            Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY
        };
        
        CriticalDensitySolver dSolver = new CriticalDensitySolver();
        dSolver.setToDebug();
        
        for (int i = 0; i < fileNames.length; ++i) {
            
            String fileName = fileNames[i];
            
            float[] values = readDataset(fileName);
            
            float critDens = dSolver.findCriticalDensity(values);
            
            float r0 = r0s[i];
            float r1 = r1s[i];
                        
            System.out.println("i=" + i + " ro-" + r0 + " r1=" + r1 + " critDens=" + critDens);

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
