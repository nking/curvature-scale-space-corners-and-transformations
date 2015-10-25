package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;

/**
 *
 * @author nichole
 */
public class BlobScaleFindersTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public void test0() throws Exception {

        boolean useBinned = false;

        String filePath1 = ResourceFinder.findFileInTestResources("brown_lowe_2003_image1.jpg");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String filePath2 = ResourceFinder.findFileInTestResources("brown_lowe_2003_image2.jpg");
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        BlobPerimeterHelper bph1 = new BlobPerimeterHelper(img1, "1");
        bph1.applySegmentation(SegmentationType.GREYSCALE_KMPP, useBinned);
        BlobCornerHelper bch1 = new BlobCornerHelper(bph1, "1");

        bch1.generatePerimeterCorners(SegmentationType.GREYSCALE_KMPP,
            useBinned);

        BlobPerimeterHelper bph2 = new BlobPerimeterHelper(img2, "2");
        bph2.applySegmentation(SegmentationType.GREYSCALE_KMPP, useBinned);
        BlobCornerHelper bch2 = new BlobCornerHelper(bph2, "1");
        bch2.generatePerimeterCorners(SegmentationType.GREYSCALE_KMPP,
            useBinned);

        BlobCornersScaleFinder bsFinder = new BlobCornersScaleFinder();

        float[] outputScaleRotTransXYStDev = new float[4];

        IntensityFeatures features1 = new IntensityFeatures(img1, 5, true);
        IntensityFeatures features2 = new IntensityFeatures(img2, 5, true);

        TransformationParameters params = bsFinder.solveForScale(
            bch1, features1,
            SegmentationType.GREYSCALE_KMPP, useBinned,
            bch2, features2,
            SegmentationType.GREYSCALE_KMPP, useBinned,
            outputScaleRotTransXYStDev);

        assertNotNull(params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.15);

    }

    public void test1() throws Exception {

        boolean useBinned = false;

        String filePath1 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0001.png");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String filePath2 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0010.png");
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        BlobPerimeterHelper bph1 = new BlobPerimeterHelper(img1, "1");
        bph1.applySegmentation(SegmentationType.COLOR_POLARCIEXY, useBinned);
        BlobCornerHelper bch1 = new BlobCornerHelper(bph1, "1");

        bch1.generatePerimeterCorners(SegmentationType.COLOR_POLARCIEXY,
            useBinned);

        BlobPerimeterHelper bph2 = new BlobPerimeterHelper(img2, "2");
        bph2.applySegmentation(SegmentationType.COLOR_POLARCIEXY, useBinned);
        BlobCornerHelper bch2 = new BlobCornerHelper(bph2, "1");
        bch2.generatePerimeterCorners(SegmentationType.COLOR_POLARCIEXY,
            useBinned);

        BlobCornersScaleFinder bsFinder = new BlobCornersScaleFinder();

        float[] outputScaleRotTransXYStDev = new float[4];

        IntensityFeatures features1 = new IntensityFeatures(img1, 5, true);
        IntensityFeatures features2 = new IntensityFeatures(img2, 5, true);

        TransformationParameters params = bsFinder.solveForScale(
            bch1, features1,
            SegmentationType.COLOR_POLARCIEXY, useBinned,
            bch2, features2,
            SegmentationType.COLOR_POLARCIEXY, useBinned,
            outputScaleRotTransXYStDev);

        assertNotNull(params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.1);
    }

    public void test2() throws Exception {

        boolean useBinned = false;

        String filePath1 = ResourceFinder.findFileInTestResources("brown_lowe_2003_image1.jpg");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String filePath2 = ResourceFinder.findFileInTestResources("brown_lowe_2003_image2.jpg");
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        BlobPerimeterHelper bph1 = new BlobPerimeterHelper(img1, "1");
        bph1.applySegmentation(SegmentationType.GREYSCALE_KMPP, useBinned);
        BlobContourHelper bch1 = new BlobContourHelper(bph1, "1");
        bch1.generatePerimeterContours(SegmentationType.GREYSCALE_KMPP,
            useBinned);

        BlobPerimeterHelper bph2 = new BlobPerimeterHelper(img2, "2");
        bph2.applySegmentation(SegmentationType.GREYSCALE_KMPP, useBinned);
        BlobContourHelper bch2 = new BlobContourHelper(bph2, "1");
        bch2.generatePerimeterContours(SegmentationType.GREYSCALE_KMPP,
            useBinned);

        BlobContoursScaleFinder bsFinder = new BlobContoursScaleFinder();

        float[] outputScaleRotTransXYStDev = new float[4];

        IntensityFeatures features1 = new IntensityFeatures(img1, 5, true);
        IntensityFeatures features2 = new IntensityFeatures(img2, 5, true);

        TransformationParameters params = bsFinder.solveForScale(
            bch1, features1,
            SegmentationType.GREYSCALE_KMPP, useBinned,
            bch2, features2,
            SegmentationType.GREYSCALE_KMPP, useBinned,
            outputScaleRotTransXYStDev);

        assertNotNull(params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.15);

    }

    public void test3() throws Exception {

        boolean useBinned = false;

        String filePath1 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0001.png");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String filePath2 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0010.png");
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        BlobPerimeterHelper bph1 = new BlobPerimeterHelper(img1, "1");
        bph1.applySegmentation(SegmentationType.COLOR_POLARCIEXY, useBinned);
        BlobContourHelper bch1 = new BlobContourHelper(bph1, "1");
        bch1.generatePerimeterContours(SegmentationType.COLOR_POLARCIEXY,
            useBinned);

        BlobPerimeterHelper bph2 = new BlobPerimeterHelper(img2, "2");
        bph2.applySegmentation(SegmentationType.COLOR_POLARCIEXY, useBinned);
        BlobContourHelper bch2 = new BlobContourHelper(bph2, "1");
        bch2.generatePerimeterContours(SegmentationType.COLOR_POLARCIEXY,
            useBinned);

        BlobContoursScaleFinder bsFinder = new BlobContoursScaleFinder();

        float[] outputScaleRotTransXYStDev = new float[4];

        IntensityFeatures features1 = new IntensityFeatures(img1, 5, true);
        IntensityFeatures features2 = new IntensityFeatures(img2, 5, true);

        TransformationParameters params = bsFinder.solveForScale(
            bch1, features1,
            SegmentationType.COLOR_POLARCIEXY, useBinned,
            bch2, features2,
            SegmentationType.COLOR_POLARCIEXY, useBinned,
            outputScaleRotTransXYStDev);

        assertNotNull(params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.1);

    }

    public void test4() throws Exception {

        boolean useBinned = true;

        String filePath1 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0001.png");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String filePath2 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0010.png");
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        BlobPerimeterHelper bph1 = new BlobPerimeterHelper(img1, "1");
        bph1.createBinnedGreyscaleImage(300);
        bph1.applySegmentation(SegmentationType.COLOR_POLARCIEXY, useBinned);
        BlobContourHelper bch1 = new BlobContourHelper(bph1, "1");
        bch1.generatePerimeterContours(SegmentationType.COLOR_POLARCIEXY,
            useBinned);

        BlobPerimeterHelper bph2 = new BlobPerimeterHelper(img2, "2");
        bph2.createBinnedGreyscaleImage(300);
        bph2.applySegmentation(SegmentationType.COLOR_POLARCIEXY, useBinned);
        BlobContourHelper bch2 = new BlobContourHelper(bph2, "1");
        bch2.generatePerimeterContours(SegmentationType.COLOR_POLARCIEXY,
            useBinned);

        BlobContoursScaleFinder bsFinder = new BlobContoursScaleFinder();

        float[] outputScaleRotTransXYStDev = new float[4];

        IntensityFeatures features1 = new IntensityFeatures(img1, 5, true);
        IntensityFeatures features2 = new IntensityFeatures(img2, 5, true);

        TransformationParameters params = bsFinder.solveForScale(
            bch1, features1,
            SegmentationType.COLOR_POLARCIEXY, useBinned,
            bch2, features2,
            SegmentationType.COLOR_POLARCIEXY, useBinned,
            outputScaleRotTransXYStDev);

        assertNotNull(params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.1);

    }
    
    public void test6() throws Exception {

        boolean useBinned = false;

        String filePath1 = ResourceFinder.findFileInTestResources("brown_lowe_2003_image1.jpg");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String filePath2 = ResourceFinder.findFileInTestResources("brown_lowe_2003_image2.jpg");
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        BlobPerimeterHelper bph1 = new BlobPerimeterHelper(img1, "1");
        bph1.applySegmentation(SegmentationType.GREYSCALE_KMPP, useBinned);
        BlobCornerHelper bch1 = new BlobCornerHelper(bph1, "1");

        bch1.generatePerimeterCorners(SegmentationType.GREYSCALE_KMPP,
            useBinned);

        BlobPerimeterHelper bph2 = new BlobPerimeterHelper(img2, "2");
        bph2.applySegmentation(SegmentationType.GREYSCALE_KMPP, useBinned);
        BlobCornerHelper bch2 = new BlobCornerHelper(bph2, "1");
        bch2.generatePerimeterCorners(SegmentationType.GREYSCALE_KMPP,
            useBinned);

        BlobCornersScaleFinder0 bsFinder = new BlobCornersScaleFinder0();

        float[] outputScaleRotTransXYStDev = new float[4];

        IntensityFeatures features1 = new IntensityFeatures(img1, 5, true);
        IntensityFeatures features2 = new IntensityFeatures(img2, 5, true);

        TransformationParameters params = bsFinder.solveForScale(
            bch1, features1,
            SegmentationType.GREYSCALE_KMPP, useBinned,
            bch2, features2,
            SegmentationType.GREYSCALE_KMPP, useBinned,
            outputScaleRotTransXYStDev);

        assertNotNull(params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.15);

    }

    public void test7() throws Exception {

        boolean useBinned = false;

        String filePath1 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0001.png");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String filePath2 = ResourceFinder.findFileInTestResources("venturi_mountain_j6_0010.png");
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        BlobPerimeterHelper bph1 = new BlobPerimeterHelper(img1, "1");
        bph1.applySegmentation(SegmentationType.COLOR_POLARCIEXY, useBinned);
        BlobCornerHelper bch1 = new BlobCornerHelper(bph1, "1");

        bch1.generatePerimeterCorners(SegmentationType.COLOR_POLARCIEXY,
            useBinned);

        BlobPerimeterHelper bph2 = new BlobPerimeterHelper(img2, "2");
        bph2.applySegmentation(SegmentationType.COLOR_POLARCIEXY, useBinned);
        BlobCornerHelper bch2 = new BlobCornerHelper(bph2, "1");
        bch2.generatePerimeterCorners(SegmentationType.COLOR_POLARCIEXY,
            useBinned);

        BlobCornersScaleFinder0 bsFinder = new BlobCornersScaleFinder0();

        float[] outputScaleRotTransXYStDev = new float[4];

        IntensityFeatures features1 = new IntensityFeatures(img1, 5, true);
        IntensityFeatures features2 = new IntensityFeatures(img2, 5, true);

        TransformationParameters params = bsFinder.solveForScale(
            bch1, features1,
            SegmentationType.COLOR_POLARCIEXY, useBinned,
            bch2, features2,
            SegmentationType.COLOR_POLARCIEXY, useBinned,
            outputScaleRotTransXYStDev);

        assertNotNull(params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.1);
    }
    
    public void test8() throws Exception {

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
        BlobCornerHelper bch2 = new BlobCornerHelper(bph2, "1");
        bch2.generatePerimeterCorners(type, useBinned);

        BlobCornersScaleFinder0 bsFinder = new BlobCornersScaleFinder0();

        float[] outputScaleRotTransXYStDev = new float[4];

        IntensityFeatures features1 = new IntensityFeatures(img1, 5, true);
        IntensityFeatures features2 = new IntensityFeatures(img2, 5, true);

        TransformationParameters params = bsFinder.solveForScale(
            bch1, features1, type, useBinned,
            bch2, features2, type, useBinned,
            outputScaleRotTransXYStDev);

        assertNotNull(params);

        assertTrue(Math.abs(params.getScale() - 1) < 0.1);
    }
}
