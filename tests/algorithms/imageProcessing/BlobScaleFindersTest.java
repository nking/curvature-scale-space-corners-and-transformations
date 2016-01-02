package algorithms.imageProcessing;

import algorithms.compGeometry.RotatedOffsets;
import algorithms.util.ResourceFinder;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;

/**
 *
 * @author nichole
 */
public class BlobScaleFindersTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public void testBL2003_1() throws Exception {

        RotatedOffsets rOffsets = RotatedOffsets.getInstance();
        
        boolean useBinned = false;

        String filePath1 = ResourceFinder.findFileInTestResources("brown_lowe_2003_image1.jpg");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String filePath2 = ResourceFinder.findFileInTestResources("brown_lowe_2003_image2.jpg");
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        SegmentationType type = SegmentationType.GREYSCALE_HIST;
        
        BlobPerimeterHelper bph1 = new BlobPerimeterHelper(img1, "1");
        bph1.applySegmentation(type, useBinned);
        BlobCornerHelper bch1 = new BlobCornerHelper(bph1, "1");
        
        bch1.generatePerimeterCorners(type, useBinned);

        BlobPerimeterHelper bph2 = new BlobPerimeterHelper(img2, "2");
        bph2.applySegmentation(type, useBinned);
        BlobCornerHelper bch2 = new BlobCornerHelper(bph2, "2");
        bch2.generatePerimeterCorners(type, useBinned);

        BlobCornersScaleFinder bsFinder = new BlobCornersScaleFinder();

        IntensityFeatures features1 = new IntensityFeatures(5, true, rOffsets);
        IntensityFeatures features2 = new IntensityFeatures(5, true, rOffsets);

        int dither = 3;
        
        MatchingSolution soln = bsFinder.solveForScale(
            bch1, features1, type, useBinned,
            bch2, features2, type, useBinned, dither);
        
        assertNotNull(soln);
        
        TransformationParameters params = soln.getParams();

        assertNotNull(params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.15);

    }

    public void estVenturi_Corners() throws Exception {

        RotatedOffsets rOffsets = RotatedOffsets.getInstance();
        
        boolean useBinned = false;

        String filePath1 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0001.png");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String filePath2 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0010.png");
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        SegmentationType type = SegmentationType.COLOR_POLARCIEXY;
                
        BlobPerimeterHelper bph1 = new BlobPerimeterHelper(img1, "1");
        bph1.applyEqualization();
        bph1.applySegmentation(type, useBinned);
        BlobCornerHelper bch1 = new BlobCornerHelper(bph1, "1");

        bch1.generatePerimeterCorners(type, useBinned);

        BlobPerimeterHelper bph2 = new BlobPerimeterHelper(img2, "2");
        bph2.applyEqualization();
        bph2.applySegmentation(type, useBinned);
        BlobCornerHelper bch2 = new BlobCornerHelper(bph2, "2");
        bch2.generatePerimeterCorners(type, useBinned);

        BlobCornersScaleFinder bsFinder = new BlobCornersScaleFinder();

        IntensityFeatures features1 = new IntensityFeatures(5, true, rOffsets);
        IntensityFeatures features2 = new IntensityFeatures(5, true, rOffsets);

        int dither = 3;
        
        MatchingSolution soln = bsFinder.solveForScale(
            bch1, features1, type, useBinned,
            bch2, features2, type, useBinned, dither);

        assertNotNull(soln);
        
        TransformationParameters params = soln.getParams();
        
        assertNotNull(params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.15);
    }

}
