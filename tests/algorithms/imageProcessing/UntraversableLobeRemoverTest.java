package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
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
        
        
        String filePath = ResourceFinder.findFileInTestResources("blob_junction01.dat");
        Set<PairInt> closedCurve = Misc.deserializeSetPairInt(filePath);
        closedCurve.remove(new PairInt(9, 3));
        closedCurve.remove(new PairInt(15, 2));
        
        int nBefore = closedCurve.size();
        
        UntraversableLobeRemover rm = new UntraversableLobeRemover();
        rm.applyFilter(closedCurve);
        
        int nAfter = closedCurve.size();
        // TODO: assert size removed
        // and that the points are all contiguous
        assertTrue(nBefore - nAfter >= 9);
        assertTrue(nAfter > 100);
               
        int w = 58;
        int h = 41;
        GreyscaleImage img2 = new GreyscaleImage(w, h);
        for (PairInt p : closedCurve) {
            img2.setValue(p.getX(), p.getY(), 255);
        }
        MiscDebug.writeImage(img2, "blob_junction.png");
    }
}
