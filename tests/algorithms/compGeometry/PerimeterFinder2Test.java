package algorithms.compGeometry;

import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.features.BlobMedialAxes;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PerimeterFinder2Test extends TestCase {
    
    public void testOrderThePerimeter() throws Exception {
        
        PerimeterFinder finder = new PerimeterFinder();
        
        String binDir = ResourceFinder.findDirectory("testresources");
        String filePath = binDir + "/test_persist_info_391558857.dat";
        
        Object[] objects = PerimeterFinder.readPeristedForTest(filePath);

        float searchRadius = ((Float)objects[0]).floatValue();
        int bmaIndex = ((Integer)objects[1]).intValue();
        Set<PairInt> perimeter = (Set<PairInt>)objects[2];
        BlobMedialAxes bma = (BlobMedialAxes)objects[3];
      
        int n = perimeter.size();
        
        PairIntArray ordered = finder.orderThePerimeter(perimeter, searchRadius,
            bma, bmaIndex);
        
        assertEquals(n, ordered.getN());
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        boolean isAdjacent = curveHelper.isAdjacent(ordered, 0, 
            ordered.getN() - 1, searchRadius);
        
        assertTrue(isAdjacent);
    }
    
    public void testOrderThePerimeter2() throws Exception {
        
        PerimeterFinder finder = new PerimeterFinder();
        
        String binDir = ResourceFinder.findDirectory("testresources");
        String filePath = binDir + "/test_persist_info_392531118.dat";
        
        Object[] objects = PerimeterFinder.readPeristedForTest(filePath);

        float searchRadius = ((Float)objects[0]).floatValue();
        int bmaIndex = ((Integer)objects[1]).intValue();
        Set<PairInt> perimeter = (Set<PairInt>)objects[2];
        BlobMedialAxes bma = (BlobMedialAxes)objects[3];
      
        int n = perimeter.size();
        
        PairIntArray ordered = finder.orderThePerimeter(perimeter, searchRadius,
            bma, bmaIndex);
        
        assertEquals(n, ordered.getN());
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        boolean isAdjacent = curveHelper.isAdjacent(ordered, 0, 
            ordered.getN() - 1, searchRadius);
        
        assertTrue(isAdjacent);
    }
    
    public void testOrderThePerimeter3() throws Exception {
        
        PerimeterFinder finder = new PerimeterFinder();
        
        String binDir = ResourceFinder.findDirectory("testresources");
        String filePath = binDir + "/test_persist_info_410598920.dat";
        
        Object[] objects = PerimeterFinder.readPeristedForTest(filePath);

        float searchRadius = ((Float)objects[0]).floatValue();
        int bmaIndex = ((Integer)objects[1]).intValue();
        Set<PairInt> perimeter = (Set<PairInt>)objects[2];
        BlobMedialAxes bma = (BlobMedialAxes)objects[3];
      
        int n = perimeter.size();
        
        PairIntArray ordered = finder.orderThePerimeter(perimeter, searchRadius,
            bma, bmaIndex);
        
        assertEquals(51, ordered.getX(0));
        assertEquals(57, ordered.getY(0));
        assertEquals(52, ordered.getX(1));
        assertEquals(56, ordered.getY(1));
        assertEquals(52, ordered.getX(2));
        assertEquals(55, ordered.getY(2));
        assertEquals(53, ordered.getX(3));
        assertEquals(57, ordered.getY(3));
        assertEquals(54, ordered.getX(4));
        assertEquals(57, ordered.getY(4));
        
        assertEquals(52, ordered.getX(n - 1));
        assertEquals(58, ordered.getY(n - 1));
        
        assertEquals(n, ordered.getN());
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        boolean isAdjacent = curveHelper.isAdjacent(ordered, 0, 
            ordered.getN() - 1, searchRadius);
        
        assertTrue(isAdjacent);
    }
    
    public void testOrderThePerimeter4() throws Exception {
        
        PerimeterFinder finder = new PerimeterFinder();
        
        String binDir = ResourceFinder.findDirectory("testresources");
        String filePath = binDir + "/test_persist_info_412377360.dat";
        
        Object[] objects = PerimeterFinder.readPeristedForTest(filePath);

        float searchRadius = ((Float)objects[0]).floatValue();
        int bmaIndex = ((Integer)objects[1]).intValue();
        Set<PairInt> perimeter = (Set<PairInt>)objects[2];
        BlobMedialAxes bma = (BlobMedialAxes)objects[3];
      
        int n = perimeter.size();
        
        PairIntArray ordered = finder.orderThePerimeter(perimeter, searchRadius,
            bma, bmaIndex);
        
        assertEquals(49, ordered.getX(0));
        assertEquals(28, ordered.getY(0));
        assertEquals(50, ordered.getX(1));
        assertEquals(28, ordered.getY(1));
        assertEquals(51, ordered.getX(2));
        assertEquals(28, ordered.getY(2));
        assertEquals(52, ordered.getX(3));
        assertEquals(27, ordered.getY(3));
        
        assertEquals(50, ordered.getX(n - 1));
        assertEquals(29, ordered.getY(n - 1));
        
        assertEquals(n, ordered.getN());
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        boolean isAdjacent = curveHelper.isAdjacent(ordered, 0, 
            ordered.getN() - 1, searchRadius);
        
        assertTrue(isAdjacent);
    }
}
