package algorithms.imageProcessing;

import algorithms.compGeometry.PerimeterFinder;
import algorithms.imageProcessing.SkylineExtractor.RemovedSets;
import algorithms.misc.HistogramHolder;
import algorithms.util.ResourceFinder;
import algorithms.util.PairInt;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import org.junit.After;
import org.junit.Before;
import static org.junit.Assert.fail;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class SkylineExtractorTest {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public SkylineExtractorTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
   
    @Test
    public void testSkyline() throws Exception {
        
    }
    
    public static void main(String[] args) {
        
        try {
            SkylineExtractorTest test = new SkylineExtractorTest();

            test.testSkyline();
        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }

    static int count=0;
    private void debugPlot(Set<PairInt> skyPoints, Image originalColorImage, 
        int xRelativeOffset, int yRelativeOffset, String label) {
        
        Image clr = originalColorImage.copyImage();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addToImage(skyPoints, xRelativeOffset, 
                yRelativeOffset, clr, 0, 0, 255);

            ImageIOHelper.writeOutputImage(
                dirPath + "/" + label + "_" + count + ".png", clr);

            //ImageDisplayer.displayImage(label, clr);
            
        } catch (IOException e) {
            System.err.println("ERROR: " + e.getMessage());
        }
    }

}
