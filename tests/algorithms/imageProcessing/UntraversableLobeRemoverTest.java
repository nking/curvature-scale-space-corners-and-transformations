package algorithms.imageProcessing;

import algorithms.compGeometry.ButterflySectionFinder;
import algorithms.compGeometry.ButterflySectionFinder.Routes;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class UntraversableLobeRemoverTest extends TestCase {
    
    public UntraversableLobeRemoverTest() {
    }
    
    public void testApplyFilter() throws Exception {
        
        /*String filePath1 = ResourceFinder.findFileInTestResources("tmp1.png");
        ImageExt img = ImageIOHelper.readImageExt(filePath1);        
        Set<PairInt> points = new HashSet<PairInt>();
        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                int r = img.getR(i, j);
                if (r > 127) {
                    PairInt p = new PairInt(i, j);
                    points.add(p);
                }
            }
        }
        Misc.persistToFile("blob_junction01.txt", points);
        */
        
        //TODO: check these numbers:
        int w = 58;
        int h = 41;
        
        String filePath = ResourceFinder.findFileInTestResources("blob_junction01.dat");
        Set<PairInt> closedCurve = Misc.deserializeSetPairInt(filePath);
        closedCurve.remove(new PairInt(9, 3));
        closedCurve.remove(new PairInt(15, 2));
        
        int nBefore = closedCurve.size();
        
        Set<PairInt> exclude = new HashSet<PairInt>();

        GreyscaleImage img1 = new GreyscaleImage(w, h);
        for (PairInt p : closedCurve) {
            img1.setValue(p.getX(), p.getY(), 255);
        }
        MiscDebug.writeImage(img1, "blob_junction_before.png");
        
        UntraversableLobeRemover rm = new UntraversableLobeRemover();
        boolean changed = rm.applyFilter(closedCurve, exclude, w, h);
        assertTrue(changed);
        
        GreyscaleImage img2 = new GreyscaleImage(w, h);
        for (PairInt p : closedCurve) {
            img2.setValue(p.getX(), p.getY(), 255);
        }
        MiscDebug.writeImage(img2, "blob_junction.png");
        
        int nAfter = closedCurve.size();
        // TODO: assert size removed
        // and that the points are all contiguous
        assertTrue(nBefore - nAfter >= 9);
        assertTrue(nAfter > 100);
               
    }
    
    public void testApplyFilter2() throws Exception {
        
        int w = 258;
        int h = 187;
        
        String filePath = ResourceFinder.findFileInTestResources("blob_butterfly_03.dat");
        Set<PairInt> closedCurve = Misc.deserializeSetPairInt(filePath);
        
        int nBefore = closedCurve.size();
        
        PairIntArray closedCurve2 = new PairIntArray();
        for (PairInt p : closedCurve) {
            closedCurve2.add(p.getX(), p.getY());
        }
        
        Image img = new Image(w, h);
        MiscDebug.writeImage(closedCurve2, img, 0, "_butterfly2_3");
        
        // exclude the butterfly sections (theyre traversable in 2 ways)
        ButterflySectionFinder finder = new ButterflySectionFinder();
        List<Routes> sections = finder.findButterflySections(closedCurve2);        
        Set<PairInt> exclude = new HashSet<PairInt>();
        for (Routes r : sections) {
            exclude.addAll(r.getRoute0());
            exclude.addAll(r.getRoute1());
        }
        
        /*
        test results should have routes in opposite directions with
        these points:
        route0=(246,126),(247,127),(248,128),(249,128) ----->
        route1=(249,130),(248,130),(247,129),(246,129) <-----
        */
        
        UntraversableLobeRemover rm = new UntraversableLobeRemover();
        rm.applyFilter(closedCurve, exclude, w, h);
        
        GreyscaleImage img2 = new GreyscaleImage(w, h);
        for (PairInt p : closedCurve) {
            img2.setValue(p.getX(), p.getY(), 255);
        }
        MiscDebug.writeImage(img2, "blob_junction3.png");
    }
    
    public void testApplyFilter3() throws Exception {
        
        int w = 347;
        int h = 277;
        
        String filePath = ResourceFinder.findFileInTestResources("tmp_blob_38286003_347_277.dat");
        Set<PairInt> closedCurve = Misc.deserializeSetPairInt(filePath);
        
        int nBefore = closedCurve.size();
        
        PairIntArray closedCurve2 = new PairIntArray();
        for (PairInt p : closedCurve) {
            closedCurve2.add(p.getX(), p.getY());
        }
        
        Image img = new Image(w, h);
        MiscDebug.writeImage(closedCurve2, img, 0, "_butterfly2_4");
        
        // exclude the butterfly sections (theyre traversable in 2 ways)
        ButterflySectionFinder finder = new ButterflySectionFinder();
        List<Routes> sections = finder.findButterflySections(closedCurve2);        
        Set<PairInt> exclude = new HashSet<PairInt>();
        for (Routes r : sections) {
            exclude.addAll(r.getRoute0());
            exclude.addAll(r.getRoute1());
        }
        
        UntraversableLobeRemover rm = new UntraversableLobeRemover();
        rm.applyFilter(closedCurve, exclude, w, h);
        
        GreyscaleImage img2 = new GreyscaleImage(w, h);
        for (PairInt p : closedCurve) {
            img2.setValue(p.getX(), p.getY(), 255);
        }
        MiscDebug.writeImage(img2, "blob_junction4.png");
    }
}
