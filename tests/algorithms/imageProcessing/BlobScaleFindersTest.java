package algorithms.imageProcessing;

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

        IntensityFeatures features1 = new IntensityFeatures(5, true);
        IntensityFeatures features2 = new IntensityFeatures(5, true);

        MatchingSolution soln = bsFinder.solveForScale(
            bch1, features1, type, useBinned,
            bch2, features2, type, useBinned);
        
        assertNotNull(soln);
        
        TransformationParameters params = soln.getParams();

        assertNotNull(params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.15);

    }

    public void estBL2003_0() throws Exception {

        boolean useBinned = false;
        
        SegmentationType type1 = SegmentationType.GREYSCALE_KMPP;
        SegmentationType type2 = SegmentationType.GREYSCALE_KMPP;

        String filePath1 = ResourceFinder.findFileInTestResources("brown_lowe_2003_image1.jpg");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String filePath2 = ResourceFinder.findFileInTestResources("brown_lowe_2003_image2.jpg");
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        BlobPerimeterHelper bph1 = new BlobPerimeterHelper(img1, "bl1_6");
        bph1.applyEqualization();
        bph1.applySegmentation(type1, useBinned);
        BlobCornerHelper bch1 = new BlobCornerHelper(bph1, "bl1_6");
        bch1.generatePerimeterCorners(type1, useBinned);

        BlobPerimeterHelper bph2 = new BlobPerimeterHelper(img2, "bl2_6");
        bph2.applyEqualization();
        bph2.applySegmentation(type2, useBinned);
        BlobCornerHelper bch2 = new BlobCornerHelper(bph2, "bl2_6");
        bch2.generatePerimeterCorners(type2, useBinned);

        BlobCornersScaleFinder0 bsFinder = new BlobCornersScaleFinder0();

        IntensityFeatures features1 = new IntensityFeatures(5, true);
        IntensityFeatures features2 = new IntensityFeatures(5, true);

        MatchingSolution soln = bsFinder.solveForScale(
            bch1, features1, type1, useBinned,
            bch2, features2, type2, useBinned);
        
        assertNotNull(soln);
        
        TransformationParameters params = soln.getParams();
        
        assertNotNull(params);

        log.info("params=" + params);
        System.out.println(params.toString());
        
        //assertTrue(Math.abs(params.getScale() - 1) < 0.15);

    }

    public void estBL2003_2() throws Exception {

        boolean useBinned = false;
        
        SegmentationType type = SegmentationType.GREYSCALE_HIST;

        String filePath1 = ResourceFinder.findFileInTestResources("brown_lowe_2003_image1.jpg");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String filePath2 = ResourceFinder.findFileInTestResources("brown_lowe_2003_image2.jpg");
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        boolean rotateBy90 = false;
        
        if (rotateBy90) {
            TransformationParameters params90 = new TransformationParameters();
            params90.setRotationInDegrees(90);
            params90.setOriginX(0);
            params90.setOriginY(0);
            params90.setTranslationX(0);
            params90.setTranslationY(img1.getWidth() - 1);
            Transformer transformer = new Transformer();
            img1 = (ImageExt) transformer.applyTransformation(img1,
                params90, img1.getHeight(), img1.getWidth());
            int z = 1;
        }
        
        BlobPerimeterHelper bph1 = new BlobPerimeterHelper(img1, "bl_contours0_1");
        //bph1.applyEqualization();
        bph1.applySegmentation(type, useBinned);
        BlobContourHelper bch1 = new BlobContourHelper(bph1, "bl_contours0_1");
        bch1.generatePerimeterContours(type, useBinned);

        BlobPerimeterHelper bph2 = new BlobPerimeterHelper(img2, "bl_contours0_2");
        //bph2.applyEqualization();
        bph2.applySegmentation(type, useBinned);
        BlobContourHelper bch2 = new BlobContourHelper(bph2, "bl_contours0_2");
        bch2.generatePerimeterContours(type, useBinned);

        BlobContoursScaleFinder0 bsFinder = new BlobContoursScaleFinder0();

        IntensityFeatures features1 = new IntensityFeatures(5, true);
        IntensityFeatures features2 = new IntensityFeatures(5, true);

        MatchingSolution soln = bsFinder.solveForScale(
            bch1, features1, type, useBinned,
            bch2, features2, type, useBinned);
        
        assertNotNull(soln);
        
        TransformationParameters params = soln.getParams();
        
        assertNotNull(params);

        log.info("params=" + params);
        System.out.println(params.toString());
        
        //assertTrue(Math.abs(params.getScale() - 1) < 0.15);

    }

    public void estVenturi_Corners() throws Exception {

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

        IntensityFeatures features1 = new IntensityFeatures(5, true);
        IntensityFeatures features2 = new IntensityFeatures(5, true);

        MatchingSolution soln = bsFinder.solveForScale(
            bch1, features1, type, useBinned,
            bch2, features2, type, useBinned);

        assertNotNull(soln);
        
        TransformationParameters params = soln.getParams();
        
        assertNotNull(params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.15);
    }

    public void estVenturi_Contours() throws Exception {

        boolean useBinned = false;

        String filePath1 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0001.png");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String filePath2 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0010.png");
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        SegmentationType type = SegmentationType.COLOR_POLARCIEXY;
        
        BlobPerimeterHelper bph1 = new BlobPerimeterHelper(img1, "ven_1");
        bph1.applyEqualization();
        bph1.applySegmentation(type, useBinned);
        BlobContourHelper bch1 = new BlobContourHelper(bph1, "ven_1");
        bch1.generatePerimeterContours(type, useBinned);

        BlobPerimeterHelper bph2 = new BlobPerimeterHelper(img2, "ven_2");
        bph2.applyEqualization();
        bph2.applySegmentation(type, useBinned);
        BlobContourHelper bch2 = new BlobContourHelper(bph2, "ven_2");
        bch2.generatePerimeterContours(type, useBinned);

        BlobContoursScaleFinder bsFinder = new BlobContoursScaleFinder();

        IntensityFeatures features1 = new IntensityFeatures(5, true);
        IntensityFeatures features2 = new IntensityFeatures(5, true);

        MatchingSolution soln = bsFinder.solveForScale(
            bch1, features1, type, useBinned,
            bch2, features2, type, useBinned);

        assertNotNull(soln);
        
        TransformationParameters params = soln.getParams();
        
        assertNotNull(params);
        
        log.info("PARAMS=" + params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.15);

    }

    public void estVenturi_Corners0() throws Exception {

        boolean useBinned = false;

        String filePath1 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0001.png");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
     

        String filePath2 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0010.png");
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
        
        SegmentationType type = SegmentationType.COLOR_POLARCIEXY;

        BlobPerimeterHelper bph1 = new BlobPerimeterHelper(img1, "venturi0_1");
        bph1.applyEqualization();
        bph1.applySegmentation(type, useBinned);
        BlobCornerHelper bch1 = new BlobCornerHelper(bph1, "venturi0_1");

        bch1.generatePerimeterCorners(type, useBinned);

        BlobPerimeterHelper bph2 = new BlobPerimeterHelper(img2, "venturi0_2");
        bph2.applyEqualization();
        bph2.applySegmentation(type, useBinned);
        BlobCornerHelper bch2 = new BlobCornerHelper(bph2, "venturi0_2");
        bch2.generatePerimeterCorners(type, useBinned);

        BlobCornersScaleFinder0 bsFinder = new BlobCornersScaleFinder0();

        IntensityFeatures features1 = new IntensityFeatures(5, true);
        IntensityFeatures features2 = new IntensityFeatures(5, true);

        MatchingSolution soln = bsFinder.solveForScale(
            bch1, features1, type, useBinned,
            bch2, features2, type, useBinned);

        assertNotNull(soln);
        
        TransformationParameters params = soln.getParams();
        
        assertNotNull(params);
        
        //log.info("params=" + params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.15);
    }
    
    public void estVenturi_Contours0() throws Exception {

        boolean useBinned = false;

        String filePath1 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0001.png");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String filePath2 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0010.png");
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        SegmentationType type = SegmentationType.COLOR_POLARCIEXY;
        
        BlobPerimeterHelper bph1 = new BlobPerimeterHelper(img1, "1");
        bph1.applyEqualization();
        bph1.applySegmentation(type, useBinned);
        BlobContourHelper bch1 = new BlobContourHelper(bph1, "1");
        bch1.generatePerimeterContours(type, useBinned);

        BlobPerimeterHelper bph2 = new BlobPerimeterHelper(img2, "2");
        bph2.applyEqualization();
        bph2.applySegmentation(type, useBinned);
        BlobContourHelper bch2 = new BlobContourHelper(bph2, "2");
        bch2.generatePerimeterContours(type, useBinned);

        BlobContoursScaleFinder0 bsFinder = new BlobContoursScaleFinder0();

        IntensityFeatures features1 = new IntensityFeatures(5, true);
        IntensityFeatures features2 = new IntensityFeatures(5, true);

        MatchingSolution soln = bsFinder.solveForScale(
            bch1, features1, type, useBinned,
            bch2, features2, type, useBinned);

        assertNotNull(soln);
        
        TransformationParameters params = soln.getParams();
        
        assertNotNull(params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.15);

    }
    
    public void est8() throws Exception {

        boolean useBinned = false;

        String filePath1 = ResourceFinder.findFileInTestResources("lena.jpg");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String filePath2 = ResourceFinder.findFileInTestResources("lena.jpg");
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        SegmentationType type = SegmentationType.COLOR_POLARCIEXY_LARGE;
        
        BlobPerimeterHelper bph1 = new BlobPerimeterHelper(img1, "1");
        bph1.applySegmentation(type, useBinned);
        BlobCornerHelper bch1 = new BlobCornerHelper(bph1, "1");
        bch1.generatePerimeterCorners(type, useBinned);

        BlobPerimeterHelper bph2 = new BlobPerimeterHelper(img2, "2");
        bph2.applySegmentation(type, useBinned);
        BlobCornerHelper bch2 = new BlobCornerHelper(bph2, "2");
        bch2.generatePerimeterCorners(type, useBinned);

        BlobCornersScaleFinder0 bsFinder = new BlobCornersScaleFinder0();

        IntensityFeatures features1 = new IntensityFeatures(5, true);
        IntensityFeatures features2 = new IntensityFeatures(5, true);

        MatchingSolution soln = bsFinder.solveForScale(
            bch1, features1, type, useBinned,
            bch2, features2, type, useBinned);

        assertNotNull(soln);
        
        TransformationParameters params = soln.getParams();
        
        assertNotNull(params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.15);
    }
}
