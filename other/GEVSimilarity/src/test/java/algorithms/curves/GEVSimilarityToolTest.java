package algorithms.curves;

import algorithms.util.ResourceFinder;
import java.io.File;
import java.io.FileFilter;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 * a tool to find common curves among the large range of possible curves
 * in the GEV function which depends upon k, sigma, and (x-mu).
 * 
 * @author nichole
 */
public class GEVSimilarityToolTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected boolean debug = true;

    protected boolean enable = true;
    
    
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    public void testPersist() throws Exception {
        
        GEVSimilarityTool tool = new GEVSimilarityTool();
        tool.resetForNWithinTen(1.0f);
        tool.calculateCurveDiffs();
        
        tool.persist();
        
        GEVSimilarityTool tool2 = new GEVSimilarityTool();
        tool2.readPersisted();
        
        assertTrue(tool.powersOfTenK == tool2.powersOfTenK);
        assertTrue(tool.powersOfTenS == tool2.powersOfTenS);
        assertTrue(tool.powersOfTenM == tool2.powersOfTenM);
        assertTrue(tool.nx == tool2.nx);
        
        assertTrue(tool.k0 > 0.0f);
        assertTrue(Math.abs(tool.k0 - tool2.k0) < 0.01*tool.k0);
        assertTrue(tool.sigma0 > 0.0f);
        assertTrue(Math.abs(tool.sigma0 - tool2.sigma0) < 0.01*tool.sigma0);
        assertTrue(tool.xMinusMu0 > 0.0f);
        assertTrue(Math.abs(tool.xMinusMu0 - tool2.xMinusMu0) < 0.01*tool.xMinusMu0);
        
        assertTrue(Math.abs(tool.nWithinTen - 1.0) < 0.001);
        assertTrue(Math.abs(tool2.nWithinTen - 1.0) < 0.001);
        
        assertTrue(tool.nKPermutations > 0.0f);
        assertTrue(Math.abs(tool.nKPermutations - tool2.nKPermutations) < 0.01*tool.nKPermutations);
        
        assertTrue(tool.nSigmaPermutations > 0.0f);
        assertTrue(Math.abs(tool.nSigmaPermutations - tool2.nSigmaPermutations) < 0.01*tool.nSigmaPermutations);
        
        assertTrue(tool.nXMinusMuPermutations > 0.0f);
        assertTrue(Math.abs(tool.nXMinusMuPermutations - tool2.nXMinusMuPermutations) < 0.01*tool.nXMinusMuPermutations);
        
        assertTrue(tool.fctr > 0.0f);
        assertTrue(Math.abs(tool.fctr - tool2.fctr) < 0.01*tool.fctr);
        
        assertTrue(tool.nCurves == tool2.nCurves);
        assertTrue(tool.nCurves > 0);
        
        assertTrue(tool.nc == tool2.nc);
        assertTrue(tool.nc > 0);
        
        assertTrue(tool.calculatedCurveDiffs);
        assertTrue(tool2.calculatedCurveDiffs);
        
        for (int i = 0; i < tool.nCurves; i++) {
            for (int j = 0; j < tool.nx; j++) {
                float d = Math.abs(tool.curves[i][j] - tool2.curves[i][j]);
                float eps = Math.abs(tool.curves[i][j]*0.01f);
                assertTrue(d <= eps);
            }
        }
        
        for (int j = 0; j < tool.nx; j++) {
            float d = Math.abs(tool.x[j] - tool2.x[j]);
            float eps = Math.abs(tool.x[j]*0.01f);
            assertTrue(d <= eps);
        }
        
        for (int i = 0; i < tool.nCurves; i++) {
            
            float d = Math.abs(tool.ks[i] - tool2.ks[i]);
            assertTrue(d <= Math.abs(tool.ks[i]*0.01));
            
            d = Math.abs(tool.sigmas[i] - tool2.sigmas[i]);
            assertTrue(d <= Math.abs(tool.sigmas[i]*0.01));
            
            d = tool.xminusmus[i] - tool2.xminusmus[i];
            assertTrue(d <= Math.abs(tool.xminusmus[i]*0.01));
            
            d = Math.abs(tool.mus[i] - tool2.mus[i]);
            assertTrue(d <= Math.abs(tool.mus[i]*0.01));
        }
        
        for (int i = 0; i < tool.nc; i++) {
            
            float d = Math.abs(tool.diff[i] - tool2.diff[i]);
            assertTrue(d <= Math.abs(tool.diff[i]*0.01));
            
            assertTrue(tool.diffIndexes[i].length == tool2.diffIndexes[i].length);
            for (int j = 0; j < tool.diffIndexes[i].length; j++) {
                d = Math.abs(tool.diffIndexes[i][j] - tool2.diffIndexes[i][j]);
                assertTrue(d <= Math.abs(tool.diffIndexes[i][j]*0.01));
            }
        }
        
        
        tool2.plotResults();
        
        String dirPath = ResourceFinder.findDirectory("target");
        File resDir = new File(dirPath);
        File[] files = resDir.listFiles(new FileTypeFilter(".png"));
        assertNotNull(files);
        assertTrue(files.length > 0);
        int count = 0;
        for (File file : files) {
            String fileName = file.getName();
            if (fileName.contains("_0")) {
                count++;
            }
        }
        assertTrue(count > 0);
    }
    
    public void test0() throws Exception {

        log.info("test0()");

        if (!enable) {
            return;
        }
        
        int persistFileNum = 2;
        
        boolean usePersisted = false;
        boolean persist = true;
              
        if (usePersisted && persist) {
            System.err.println("Cannot have usePersisted=true and persist=true");
            return;
        }
        
        GEVSimilarityTool tool = new GEVSimilarityTool();

        if (usePersisted) {
            tool.readPersisted(persistFileNum);
        } else {
            tool.calculateCurveDiffs();
        }
        
        if (persist) {
            tool.persist(persistFileNum);
        }
        
        tool.plotResults();
    }
  
    protected static class FileTypeFilter implements FileFilter {
        protected final String suffix;
        public FileTypeFilter(String fileType) {
            this.suffix = fileType;
        }
        @Override
        public boolean accept(File pathname) {
            if (pathname.getName().endsWith(suffix)) {
                return true;
            }
            return false;
        }
    }

}
