package algorithms.imageProcessing.features;

import algorithms.imageProcessing.*;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class BlobCornerHelperTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public BlobCornerHelperTest() {
    }

    public void test1() throws Exception {

        boolean useBinned = false;

        String filePath = ResourceFinder.findFileInTestResources("blox.gif");
        ImageExt img = ImageIOHelper.readImageExt(filePath);

        BlobPerimeterCornerHelper bph = new BlobPerimeterCornerHelper(img, "blox");
        bph.increaseLargestGroupLimit(100000);

        assertEquals(1, bph.getBinFactor(true));

        bph.applySegmentation(SegmentationType.GREYSCALE_KMPP, useBinned);

        bph.applySegmentation(SegmentationType.COLOR_POLARCIEXY, useBinned);

        List<List<CornerRegion>> cList1 =
            bph.generatePerimeterCorners(SegmentationType.GREYSCALE_KMPP,
                useBinned);

        List<List<CornerRegion>> cList2 =
            bph.generatePerimeterCorners(SegmentationType.COLOR_POLARCIEXY,
                useBinned);

        assertNotNull(cList1);
        assertNotNull(cList2);

        assertTrue(cList1.size() >= 5);
        assertTrue(cList2.size() >= 5);

        int nNonZero1 = bph.sumPointsOfInterest(
            SegmentationType.GREYSCALE_KMPP, useBinned);
        int nNonZero2 = bph.sumPointsOfInterest(
            SegmentationType.COLOR_POLARCIEXY, useBinned);

        assertTrue(nNonZero1 >= 2*5);
        assertTrue(nNonZero2 >= 2*5);
    }

    public void test2() throws Exception {

        boolean useBinned = true;

        String filePath = ResourceFinder.findFileInTestResources("lena.jpg");
        ImageExt img = ImageIOHelper.readImageExt(filePath);

        BlobPerimeterCornerHelper bph = new BlobPerimeterCornerHelper(img, "lena");
        bph.increaseLargestGroupLimit(100000);

        bph.createBinnedGreyscaleImage(300);

        assertEquals(2, bph.getBinFactor(true));

        bph.applySegmentation(SegmentationType.GREYSCALE_KMPP, useBinned);

        bph.applySegmentation(SegmentationType.COLOR_POLARCIEXY, useBinned);

        List<List<CornerRegion>> cList1 =
            bph.generatePerimeterCorners(SegmentationType.GREYSCALE_KMPP,
                useBinned);

        List<List<CornerRegion>> cList2 =
            bph.generatePerimeterCorners(SegmentationType.COLOR_POLARCIEXY,
                useBinned);

        assertNotNull(cList1);
        assertNotNull(cList2);

        assertTrue(cList1.size() >= 5);
        assertTrue(cList2.size() >= 10);

        int nNonZero1 = bph.sumPointsOfInterest(
            SegmentationType.GREYSCALE_KMPP, useBinned);
        int nNonZero2 = bph.sumPointsOfInterest(
            SegmentationType.COLOR_POLARCIEXY, useBinned);

        assertTrue(nNonZero1 >= 2*5);
        assertTrue(nNonZero2 >= 2*10);
    }
}
