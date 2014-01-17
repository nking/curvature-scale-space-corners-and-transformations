package algorithms.util;

import java.io.File;

import junit.framework.TestCase;

public class ResourceFinderTest extends TestCase {
    
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
   
    public void test() throws Exception {
        
        String dirPath = ResourceFinder.findResourcesDirectory();
        assertNotNull(dirPath);
        assertTrue( (new File(dirPath)).exists());
        
        dirPath = ResourceFinder.findTestDirectory();
        assertNotNull(dirPath);
        assertTrue( (new File(dirPath)).exists());
        
        dirPath = ResourceFinder.findTmpDataDirectory();
        assertNotNull(dirPath);
        assertTrue( (new File(dirPath)).exists());
        
        String flPath = ResourceFinder.getAFilePathInTmpData("tmp");
        assertNotNull(flPath);
        File f = new File(flPath);
        assertTrue(f.createNewFile());
        String fn = f.getName();
        File f2 = ResourceFinder.findFileInTmpData(fn);
        assertNotNull(f2);
        assertTrue(f2.exists());
        f2.delete();
        
        flPath = ResourceFinder.writeToTmpData("asdf", "tmp");
        assertNotNull(flPath);
        f2 = new File(flPath);
        assertTrue(f2.exists());
        f2.delete();
        
        flPath = ResourceFinder.writeToTestResources("asdf", "tmp");
        assertNotNull(flPath);
        f2 = new File(flPath);
        assertTrue(f2.exists());
        f2.delete();
        
        flPath = ResourceFinder.getAFilePathInCWD("tmp");
        assertNotNull(flPath);
        f2 = new File(flPath);
        f2.createNewFile();
        assertTrue(f2.exists());
        f2.delete();
        
        String fileNamePrefix = "tmp";
        flPath = ResourceFinder.writeToTestResources("asdf", "tmp_tmp");
        File[] fls = ResourceFinder.findFilesInTestResources(fileNamePrefix);
        assertNotNull(fls);
        assertTrue(fls.length > 0);
        flPath = ResourceFinder.findFileInTestResources("tmp_tmp");
        assertNotNull(flPath);
        f2 = new File(flPath);
        assertTrue(f2.exists());
        f2.delete();

        
        assertNotNull(ResourceFinder.findFileInResources("logging.properties"));
        
        assertNotNull(ResourceFinder.findFileInTestResources("logging.properties"));
        
        dirPath = ResourceFinder.findResourcesDirectory();
        dirPath = dirPath + System.getProperty("file.separator") + "tmpdir";
        f2 = new File(dirPath);
        f2.mkdir();
        flPath = dirPath + System.getProperty("file.separator") + "tmp";
        f2 = new File(flPath);
        f2.createNewFile();
        
        flPath = ResourceFinder.findFileInResources("tmpdir", "tmp");
        f2 = new File(flPath);
        assertTrue(f2.exists());
        f2.delete();
        
        
        ResourceFinder rf = new ResourceFinder();
    }
    
}
