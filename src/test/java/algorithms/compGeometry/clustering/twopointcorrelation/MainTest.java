package algorithms.compGeometry.clustering.twopointcorrelation;

import java.io.File;
import java.net.URL;
import java.util.Date;
import java.util.Arrays;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MainTest extends TestCase {
    
    /**
     *
     * @throws Exception
     */
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    /**
     *
     * @throws Exception
     */
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    /**
     *
     * @throws Exception
     */
    public void test() throws Exception {
        
        ClassLoader cls = Main.class.getClassLoader();
        URL url = cls.getResource(".");
        String filePath = url.getPath();
        
        String srchFor = "two-point-correlation";
        
        int idx = filePath.indexOf(srchFor);
        if (idx == -1) {
            throw new IllegalStateException("looking for base directory of project as 'two-point-correlation' but cannot find it"
                + ".  (" + filePath +")");
        }
        
        String sep = System.getProperty("file.separator");
        
        String projectPath = filePath.substring(0, idx + srchFor.length());

        filePath = projectPath + sep + "data" + sep + "dbscan.txt";
        
        boolean exists = new File(filePath).exists();

        if (!exists) {
            System.out.println("could not find " + filePath);
        }

        assertTrue(exists);
                
        String[] args = new String[]{"--file", filePath};
        
        String expectedOutputPath = projectPath + sep + "bin" + sep + "test-classes" + sep + "twoptcorrelation3.html";
        
        File file = new File(expectedOutputPath);
        
        long lastModified = file.exists() ? file.lastModified() : 0;
        file.delete();
        
        System.out.println("running Main with args=" + Arrays.toString(args));

        Main.main(args);
                
        file = new File(expectedOutputPath);
        
        if (!file.exists()) {
            expectedOutputPath = projectPath + sep + "bin" + sep + "instr-classes" + sep + "twoptcorrelation3.html";
            file = new File(expectedOutputPath);
        }
        
        long diff = file.lastModified() - lastModified;
        System.out.println("diff=" + diff + " lm0=" + lastModified + " lm1=" + file.lastModified() + new Date(file.lastModified()));
        
        assertTrue(diff >= 0);
        
        
        args = new String[]{"--file", filePath, "--twosigma"};
        file = new File(expectedOutputPath);
        lastModified = file.lastModified();
        file.delete();
        Main.main(args);
        file = new File(expectedOutputPath);
        diff = file.lastModified() - lastModified;
        System.out.println("diff=" + diff + " lm0=" + lastModified + " lm1=" + file.lastModified() + new Date(file.lastModified()));
        assertTrue(diff >= 0);
        
        
        args = new String[]{"--file", filePath, "--threesigma"};
        file = new File(expectedOutputPath);
        lastModified = file.lastModified();
        file.delete();
        Main.main(args);
        file = new File(expectedOutputPath);
        diff = file.lastModified() - lastModified;
        System.out.println("diff=" + diff + " lm0=" + lastModified + " lm1=" + file.lastModified() + new Date(file.lastModified()));
        assertTrue(diff >= 0);
        
        args = new String[]{"--file", filePath, "--background", "0.249", "--backgrounderror", "0.02"};
        file = new File(expectedOutputPath);
        lastModified = file.lastModified();
        file.delete();
        Main.main(args);
        file = new File(expectedOutputPath);
        diff = file.lastModified() - lastModified;
        System.out.println("diff=" + diff + " lm0=" + lastModified + " lm1=" + file.lastModified() + new Date(file.lastModified()));
        assertTrue(diff >= 0);
    }
    
}
