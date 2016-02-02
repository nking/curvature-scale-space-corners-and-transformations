package algorithms.imageProcessing.features;

import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class BlobMedialAxesTest extends TestCase {
    
    public BlobMedialAxesTest() {
    }

    public void testGetNumberOfItems() throws IOException {
        
        Set<PairInt> points = new HashSet<PairInt>();
        points.add(new PairInt(0, 0));
        points.add(new PairInt(1, 0));
        List<Set<PairInt>> list = new ArrayList<Set<PairInt>>();
        list.add(points);
        points = new HashSet<PairInt>();
        points.add(new PairInt(10, 10));
        points.add(new PairInt(11, 10));
        list.add(points);
        
        List<Double> ls = new ArrayList<Double>();
        ls.add(Double.valueOf(10.));
        ls.add(Double.valueOf(15.));
        List<Double> as = new ArrayList<Double>();
        as.add(Double.valueOf(1.));
        as.add(Double.valueOf(1.5));
        List<Double> bs = new ArrayList<Double>();
        bs.add(Double.valueOf(1.));
        bs.add(Double.valueOf(0.5));        
        
        BlobMedialAxes bma = new BlobMedialAxes(list, ls, as, bs);
        
        assertEquals(2, bma.getNumberOfItems());
        
        PairInt p0 = bma.findClosestPoint(0, 2, 0);
        assertEquals(1, p0.getX());
        assertEquals(0, p0.getY());
        p0 = bma.findClosestPoint(0, 1, 1);
        assertEquals(1, p0.getX());
        assertEquals(0, p0.getY());
        p0 = bma.findClosestPoint(0, -1, 0);
        assertEquals(0, p0.getX());
        assertEquals(0, p0.getY());
        p0 = bma.findClosestPoint(1, 12, 10);
        assertEquals(11, p0.getX());
        assertEquals(10, p0.getY());
    
        PairInt xyCen = bma.getOriginalBlobXYCentroid(0);
        assertTrue(Math.abs(xyCen.getX() - 0.5) < 0.6);
        assertTrue(Math.abs(xyCen.getY() - 0.0) < 0.01);
        
        xyCen = bma.getOriginalBlobXYCentroid(1);
        assertTrue(Math.abs(xyCen.getX() - 10.5) < 0.6);
        assertTrue(Math.abs(xyCen.getY() - 10.0) < 0.01);
   
        float[] lab = bma.getLABColors(0);
        assertTrue(Math.abs(lab[0] - ls.get(0).doubleValue()) < 0.001);
        assertTrue(Math.abs(lab[1] - as.get(0).doubleValue()) < 0.001);
        assertTrue(Math.abs(lab[2] - bs.get(0).doubleValue()) < 0.001);
        
        lab = bma.getLABColors(1);
        assertTrue(Math.abs(lab[0] - ls.get(1).doubleValue()) < 0.001);
        assertTrue(Math.abs(lab[1] - as.get(1).doubleValue()) < 0.001);
        assertTrue(Math.abs(lab[2] - bs.get(1).doubleValue()) < 0.001);
        
        String binDir = ResourceFinder.findDirectory("bin");
        String filePath = binDir + "/bma_persist.dat"; 
        bma.peristToFile(filePath);
        
        File chk = new File(filePath);
        assertTrue(chk.exists());
        
        //----
        BlobMedialAxes bma2 = new BlobMedialAxes(filePath);
        
        assertEquals(2, bma2.getNumberOfItems());
        
        p0 = bma2.findClosestPoint(0, 2, 0);
       
        assertEquals(1, p0.getX());
        assertEquals(0, p0.getY());
    
        xyCen = bma2.getOriginalBlobXYCentroid(0);
        assertTrue(Math.abs(xyCen.getX() - 0.5) < 0.6);
        assertTrue(Math.abs(xyCen.getY() - 0.0) < 0.01);
        
        xyCen = bma2.getOriginalBlobXYCentroid(1);
        assertTrue(Math.abs(xyCen.getX() - 10.5) < 0.6);
        assertTrue(Math.abs(xyCen.getY() - 10.0) < 0.01);
   
        lab = bma2.getLABColors(0);
        assertTrue(Math.abs(lab[0] - ls.get(0).doubleValue()) < 0.001);
        assertTrue(Math.abs(lab[1] - as.get(0).doubleValue()) < 0.001);
        assertTrue(Math.abs(lab[2] - bs.get(0).doubleValue()) < 0.001);
        
        lab = bma2.getLABColors(1);
        assertTrue(Math.abs(lab[0] - ls.get(1).doubleValue()) < 0.001);
        assertTrue(Math.abs(lab[1] - as.get(1).doubleValue()) < 0.001);
        assertTrue(Math.abs(lab[2] - bs.get(1).doubleValue()) < 0.001);
        
        //-------
        List<Integer> removeIndexes = new ArrayList<Integer>();
        removeIndexes.add(Integer.valueOf(0));
        bma2.removeIndexes(removeIndexes);
                
        assertEquals(1, bma2.getNumberOfItems());
        
        p0 = bma2.findClosestPoint(0, 2, 0);
        assertEquals(1, p0.getX());
        assertEquals(0, p0.getY());
    
        xyCen = bma2.getOriginalBlobXYCentroid(0);
        assertTrue(Math.abs(xyCen.getX() - 10.5) < 0.6);
        assertTrue(Math.abs(xyCen.getY() - 10.0) < 0.01);
        
        lab = bma2.getLABColors(0);
        assertTrue(Math.abs(lab[0] - ls.get(1).doubleValue()) < 0.001);
        assertTrue(Math.abs(lab[1] - as.get(1).doubleValue()) < 0.001);
        assertTrue(Math.abs(lab[2] - bs.get(1).doubleValue()) < 0.001);
        
    }

}
