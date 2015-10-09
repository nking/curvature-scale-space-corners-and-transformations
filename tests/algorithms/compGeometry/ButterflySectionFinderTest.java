package algorithms.compGeometry;

import algorithms.compGeometry.ButterflySectionFinder.Routes;
import algorithms.imageProcessing.Image;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ButterflySectionFinderTest extends TestCase {
    
    public void testFindButterflySections() throws Exception {
        
        LinkedHashSet<PairInt> tmp = new LinkedHashSet<PairInt>();
        tmp.add(new PairInt(1, 1));
        tmp.add(new PairInt(1, 1));
        assertTrue(tmp.size() == 1);
        
        String fileName = "blob_butterfly_01.dat";
        //String fileName = "blob_not_butterfly_01.dat";

        String filePath = ResourceFinder.findFileInTestResources(fileName);

        PairIntArray closedCurve = Misc.deserializePairIntArray(filePath);
        
        Image img = new Image(512, 512);
        MiscDebug.writeImage(closedCurve, img, 0, "_butterfly");

        assertNotNull(closedCurve);
        
        assertTrue(closedCurve.getN() > 0);
        
        ButterflySectionFinder finder = new ButterflySectionFinder();
        
        /*
        test results should have routes in opposite directions with
        these points:
        route0=(374, 263),(375, 263),(376, 263),(377, 263),(378,264)  --->
        route1=(374, 261),(375, 262),(376, 262),(377, 262),(378,262)  <---
        */
        
        List<Routes> sections = finder.findButterflySections(closedCurve);
        
        assertTrue(sections.size() == 1);
        
        /*
        Set<PairInt> bPoints = sections.get(0);
        
        assertTrue(bPoints.size() == 10);
        
        assertTrue(bPoints.contains(new PairInt(374, 263)));
        assertTrue(bPoints.contains(new PairInt(374, 261)));
        assertTrue(bPoints.contains(new PairInt(375, 263)));
        assertTrue(bPoints.contains(new PairInt(375, 262)));
        assertTrue(bPoints.contains(new PairInt(376, 263)));
        assertTrue(bPoints.contains(new PairInt(376, 262)));
        assertTrue(bPoints.contains(new PairInt(377, 263)));
        assertTrue(bPoints.contains(new PairInt(377, 262)));
        assertTrue(bPoints.contains(new PairInt(378, 264)));
        assertTrue(bPoints.contains(new PairInt(378, 262)));
        */
    }
    
    public void estFindButterflySections3() throws Exception {
        
        int w = 258;
        int h = 187;
        String fileName = "blob_butterfly_03.dat";

        String filePath = ResourceFinder.findFileInTestResources(fileName);

        PairIntArray closedCurve = Misc.deserializePairIntArray(filePath);
        
        Image img = new Image(w, h);
        MiscDebug.writeImage(closedCurve, img, 0, "_butterfly3");

        assertNotNull(closedCurve);
        
        assertTrue(closedCurve.getN() > 0);
        
        ButterflySectionFinder finder = new ButterflySectionFinder();
        
        List<Routes> sections = finder.findButterflySections(closedCurve);
        
        /*
        test results should have routes in opposite directions with
        these points:
        route0=(249,130),(248,130),(247,129),(246,129) <-----
        route1=(246,126),(247,127),(248,128),(249,129) ----->
        */
        
        assertTrue(sections.size() == 1);
        
        /*
        Set<PairInt> bPoints = sections.get(0);
        
        assertTrue(bPoints.size() == 8);
        */
    }
}
