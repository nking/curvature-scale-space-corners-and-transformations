package algorithms.compGeometry;

import algorithms.compGeometry.ButterflySectionFinder.Routes;
import algorithms.imageProcessing.Image;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.Iterator;
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
        
        Routes routes = sections.get(0);
        
        Iterator<PairInt> r = routes.getRoute0().iterator();
        int nIter = 0;
        while (r.hasNext()) {
            PairInt p = r.next();
            switch(nIter) {
                case 0:
                    assertEquals(p, new PairInt(378,264));
                    break;
                case 1:
                    assertEquals(p, new PairInt(377,263));
                    break;
                case 2:
                    assertEquals(p, new PairInt(376,263));
                    break;
                case 3:
                    assertEquals(p, new PairInt(375,263));
                    break;
                case 4:
                    assertEquals(p, new PairInt(374,263));
                    break;
            }
            nIter++;
        }
        
        r = routes.getRoute1().iterator();
        nIter = 0;
        while (r.hasNext()) {
            PairInt p = r.next();
            switch(nIter) {
                case 0:
                    assertEquals(p, new PairInt(374,261));
                    break;
                case 1:
                    assertEquals(p, new PairInt(375,262));
                    break;
                case 2:
                    assertEquals(p, new PairInt(376,262));
                    break;
                case 3:
                    assertEquals(p, new PairInt(377,262));
                    break;
                case 4:
                    assertEquals(p, new PairInt(378,262));
                    break;
            }
            nIter++;
        }
        
    }
    
    public void testFindButterflySections3() throws Exception {
        
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
        route1=(249,130),(248,130),(247,129),(246,129) <-----
        route0=(246,126),(247,127),(248,128),(249,128) ----->
        */
        
        assertTrue(sections.size() == 1);
        
        Routes routes = sections.get(0);
        
        Iterator<PairInt> r = routes.getRoute0().iterator();
        int nIter = 0;
        while (r.hasNext()) {
            PairInt p = r.next();
            switch(nIter) {
                case 0:
                    assertEquals(p, new PairInt(246,126));
                    break;
                case 1:
                    assertEquals(p, new PairInt(247,127));
                    break;
                case 2:
                    assertEquals(p, new PairInt(248,128));
                    break;
                case 3:
                    assertEquals(p, new PairInt(249,128));
                    break;
            }
            nIter++;
        }
        
        r = routes.getRoute1().iterator();
        nIter = 0;
        while (r.hasNext()) {
            PairInt p = r.next();
            switch(nIter) {
                case 0:
                    assertEquals(p, new PairInt(249,130));
                    break;
                case 1:
                    assertEquals(p, new PairInt(248,130));
                    break;
                case 2:
                    assertEquals(p, new PairInt(247,129));
                    break;
                case 3:
                    assertEquals(p, new PairInt(246,129));
                    break;
            }
            nIter++;
        }
    }
}
