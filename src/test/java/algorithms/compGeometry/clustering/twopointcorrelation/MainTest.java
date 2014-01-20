package algorithms.compGeometry.clustering.twopointcorrelation;

import java.io.File;
import java.net.URL;
import java.util.Date;

import junit.framework.TestCase;

public class MainTest extends TestCase {
    
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
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

        filePath = projectPath + sep + "tmpdata" + sep + "dbscan.txt";
        
        assertTrue(new File(filePath).exists());
                
        String[] args = new String[]{"--file", filePath};
        
        String expectedOutputPath = projectPath + sep + "bin" + sep + "test-classes" + sep + "twoptcorrelation3.html";
        
        File file = new File(expectedOutputPath);
        
        long lastModified = file.exists() ? file.lastModified() : 0;
        file.delete();
        
        Main.main(args);
                
        file = new File(expectedOutputPath);
        
        if (!file.exists()) {
            expectedOutputPath = projectPath + sep + "bin" + sep + "instr-classes" + sep + "twoptcorrelation3.html";
            file = new File(expectedOutputPath);
        }
        
        long diff = file.lastModified() - lastModified;
        System.out.println("diff=" + diff + " lm0=" + lastModified + " lm1=" + file.lastModified() + new Date(file.lastModified()));
        
        assertTrue(diff >= 0);
        
        
        args = new String[]{"--file", filePath, "--twosigma", "--debug"};
        file = new File(expectedOutputPath);
        lastModified = file.lastModified();
        file.delete();
        Main.main(args);
        file = new File(expectedOutputPath);
        diff = file.lastModified() - lastModified;
        System.out.println("diff=" + diff + " lm0=" + lastModified + " lm1=" + file.lastModified() + new Date(file.lastModified()));
        assertTrue(diff >= 0);
        
        
        args = new String[]{"--file", filePath, "--threesigma", "--debug"};
        file = new File(expectedOutputPath);
        lastModified = file.lastModified();
        file.delete();
        Main.main(args);
        file = new File(expectedOutputPath);
        diff = file.lastModified() - lastModified;
        System.out.println("diff=" + diff + " lm0=" + lastModified + " lm1=" + file.lastModified() + new Date(file.lastModified()));
        assertTrue(diff >= 0);
        
        args = new String[]{"--file", filePath, "--background", "0.249", "--backgrounderror", "0.02", "--debug"};
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
