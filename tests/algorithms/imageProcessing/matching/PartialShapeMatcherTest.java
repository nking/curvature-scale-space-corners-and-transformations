package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.compGeometry.MiscellaneousCurveHelper;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.MiscMath;
import algorithms.util.*;

import java.io.IOException;
import java.util.*;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;

/**
 *
 * @author nichole
 */
public class PartialShapeMatcherTest extends TestCase {

    // while configuring new code, not using asserts
    boolean enableAsserts = false;
    
    private Logger log = Logger.getLogger(
        this.getClass().getName());

    public PartialShapeMatcherTest() {
    }

    public void testScissorsMatch0() throws Exception {

        // 60
        PairFloatArray p = getScissors1();
        plot(p, 200);

        // 63
        PairFloatArray q = getScissors2();
        plot(q, 201);

        log.info("p.n=" + p.getN() + " q.n=" + q.getN());

        //q.rotateLeft(q.getN() - 3);
        PartialShapeMatcher shapeMatcher = new PartialShapeMatcher();
        //shapeMatcher.overrideSamplingDistance(1);
        //shapeMatcher.setToDebug();
        //shapeMatcher.overrideMinimumLength(4);
        //shapeMatcher.overrideToSearchAllBlockSizes();;

        List<Match.Points> results = shapeMatcher.match(p, q);
        assertFalse(results.isEmpty());
        plotResults(results, p, q, 4, "_scissors_offset0_corres", false);

    }

    private List<String> plotResults(List<Match.Points> results, PairFloatArray p, PairFloatArray q,
                               int spacing, String fileSuffix, boolean printIndexes) throws IOException {
        List<String> writtenFiles = new ArrayList<>();

        PairIntArray pInt = new PairIntArray();
        PairIntArray qInt = new PairIntArray();
        for (int i = 0; i < p.getN(); ++i) {
            pInt.add(Math.round(p.getX(i)), Math.round(p.getY(i)));
        }
        for (int i = 0; i < q.getN(); ++i) {
            qInt.add(Math.round(q.getX(i)), Math.round(q.getY(i)));
        }

        for (int i = 0; i < results.size(); ++i) {
            CorrespondencePlotter plotter = new CorrespondencePlotter(pInt, qInt);
            Match.Points result = results.get(i);
            for (int ii = 0; ii < result.pIdxs.length; ii += spacing) {
                int idx1 = result.pIdxs[ii];
                int idx2 = result.qIdxs[ii];
                int x1 = pInt.getX(idx1);
                int y1 = pInt.getY(idx1);
                int x2 = qInt.getX(idx2);
                int y2 = qInt.getY(idx2);
                System.out.println(String.format(
                "(%d, %d) <=> (%d, %d)", x1, y1, x2, y2));
                if (printIndexes) {

                }

                plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 0);
            }
            String filePath = plotter.writeImage(String.format("_%s_%d", fileSuffix, i));
            writtenFiles.add(filePath);
        }
        return writtenFiles;
    }

    public void testScissorsMatch16() throws Exception {
        
        // rotate points p so that start points are
        // different and assert that wrap around is
        // handled correctly
        PairFloatArray p = getScissors1();
        p.rotateLeft(16);
        plot(p, 200);

        PairFloatArray q = getScissors2();
        plot(q, 201);

        log.info("p.n=" + p.getN() + " q.n=" + q.getN());

        //q.rotateLeft(q.getN() - 3);
        PartialShapeMatcher shapeMatcher = new PartialShapeMatcher();
        shapeMatcher.overrideSamplingDistance(1);
        //shapeMatcher.setToDebug();
        //shapeMatcher.overrideMinimumLength(4);

        // articulated:
        List<Match.Points> results = shapeMatcher.match(p, q);
        assertFalse(results.isEmpty());
        plotResults(results, p, q, 4, "_scissors_offset16_corres", false);

    }
    
    public void testScissorsMatch16_scaled() throws Exception {

        algorithms.imageProcessing.util.MiscellaneousCurveHelper curveHelper =
            new algorithms.imageProcessing.util.MiscellaneousCurveHelper();
        
        // rotate points p so that start points are
        // different and assert that wrap around is
        // handled correctly

        PairFloatArray p = getScissors1();
        p.rotateLeft(16);

        /*
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        double[] cenXY = ch.calculateXYCentroids(p.getX(), p.getY());
        TransformationParameters params = new TransformationParameters();
        params.setOriginX((float)cenXY[0]);
        params.setOriginX((float)cenXY[1]);
        params.setScale(0.5f);
        Transformer transformer = new Transformer();
        PairFloatArray b = transformer.applyTransformation(params, p);
         */
        plot(p, 210);

        PairFloatArray q = getScissors2();
        //q = curveHelper.scaleDown(q, 0.5f);
        plot(q, 211);

        log.info("p.n=" + p.getN() + " q.n=" + q.getN());

        //q.rotateLeft(q.getN() - 3);
        PartialShapeMatcher shapeMatcher = new PartialShapeMatcher();
        shapeMatcher.setToUseSameNumberOfPoints();
        //shapeMatcher.setToDebug();

        List<Match.Points> results = shapeMatcher.match(p, q);
        assertFalse(results.isEmpty());
        plotResults(results, p, q, 1, "_scissors_offset016_corres", false);

    }

    public void testAndroidGingerbreadSameScale() throws Exception {
        /*
        For same scale and very little noise, the
        best results are obtained with options:
            matcher.overrideSamplingDistance(dp);
            matcher._overrideToThreshhold(0.2f);
        If there is alot of noise, should add:
            matcher.setToRemoveOutliers();
        */
       
        String fileName0
            = "android_statues_03_sz1_mask_small.png";
        int idx = fileName0.lastIndexOf(".");
        String fileName0Root = fileName0.substring(0, idx);
        String filePath0 = ResourceFinder
            .findFileInTestResources(fileName0);
        ImageExt img0 = ImageIOHelper.readImageExt(filePath0);

        PairFloatArray p = extractOrderedBoundary(img0);
        plot(p, 100);

        String fileName1 = "";

        //for (int type = 0; type < 2; ++type) {
        for (int type = 1; type < 2; ++type) {
            //for (int i = 0; i < 4; ++i) {
            for (int i = 2; i < 3; ++i) {
                
                //int dp = 2;

                PartialShapeMatcher shapeMatcher = new PartialShapeMatcher();
                //matcher.setToDebug();
                if (type == 0) {
                    shapeMatcher._overrideToThreshhold(0.2f);
                }
                //shapeMatcher.overrideSamplingDistance(dp);

                switch(i) {
                    case 0: {
                        fileName1
                            = "android_statues_01_sz1_mask_small.png";
                        break;
                    }
                    case 1: {
                        fileName1 = "android_statues_02_sz1_mask_small.png";
                        break;
                    }
                    case 2: {
                        fileName1 = "android_statues_03_sz1_mask_small.png";
                        break;
                    }
                    case 3: {
                        fileName1 = "android_statues_04_sz1_mask_small.png";
                        break;
                    }
                    default: {
                        break;
                    }
                }

                idx = fileName1.lastIndexOf(".");
                String fileName1Root = fileName1.substring(0, idx);

                String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
                ImageExt img = ImageIOHelper.readImageExt(filePath1);

                PairFloatArray q = extractOrderedBoundary(img);
                plot(q, (i+1)*100 + 1);

                log.info("type=" + type + " i=" + i + " matching " + fileName0Root
                + " to " + fileName1Root + " (" + p.getN()
                + " points to " + q.getN() + " points");

                List<Match.Points> results = shapeMatcher.match(p, q);
                assertFalse(results.isEmpty());
                plotResults(results, p, q, 4,
                        "_" + fileName1Root + "_corres_" + type, false);
            }
        }
    }
    
    public void testAndroidGingerbreadDiffScale() throws Exception {

        /*
        for a different scale matching,
        best results are obtained with:
            matcher.setToUseSameNumberOfPoints();
            matcher.overrideSamplingDistance(1);
            matcher.setToUseEuclidean();
        */
        
        MiscellaneousCurveHelper curveHelper = 
            new MiscellaneousCurveHelper();
                
        SIGMA sigma = SIGMA.ONE;

        String fileName0
            = "android_statues_04_sz2_mask2_small.png";
        int idx = fileName0.lastIndexOf(".");
        String fileName0Root = fileName0.substring(0, idx);
        String filePath0 = ResourceFinder
            .findFileInTestResources(fileName0);
        ImageExt img0 = ImageIOHelper.readImageExt(filePath0);

        PairFloatArray p = extractOrderedBoundary(img0, sigma);
        //p = curveHelper.scaleDown(p, 0.5f);
        plot(p, 100);

        String fileName1 = "";

        for (int i = 0; i < 4; ++i) {

            switch(i) {
                case 0: {
                    fileName1
                        = "android_statues_01_sz2_mask_small.png";
                    break;
                }
                case 1: {
                    fileName1 = "android_statues_02_sz2_mask_small.png";
                    break;
                }
                case 2: {
                    fileName1 = "android_statues_03_sz2_mask_small.png";
                    break;
                }
                case 3: {
                    fileName1 = "android_statues_04_sz2_mask_small.png";
                    break;
                }
                default: {
                    break;
                }
            }

            idx = fileName1.lastIndexOf(".");
            String fileName1Root = fileName1.substring(0, idx);

            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img = ImageIOHelper.readImageExt(filePath1);

            PairFloatArray q = extractOrderedBoundary(img, sigma);
            plot(q, (i+1)*100 + 1);

            log.info("matching " + fileName0Root + " to " + fileName1Root + " (" + p.getN()
            + " points to " + q.getN() + " points");

            //int dp = 1;
            PartialShapeMatcher shapeMatcher = new PartialShapeMatcher();
            //shapeMatcher.setToDebug();
            shapeMatcher.setToUseSameNumberOfPoints();
            //shapeMatcher.overrideSamplingDistance(4);
            //shapeMatcher._overrideToThreshhold(0.5f);
            //shapeMatcher.overrideMinimumLength(3);

            List<Match.Points> results = shapeMatcher.match(p, q);
            assertFalse(results.isEmpty());
            plotResults(results, p, q, 4,
                    "_" + fileName1Root + "_corres_", false);

            /*
            float expFrac = 0.4f;
            switch (i) {
                case 0:
                    break;
                case 1:
                    break;
                case 2:
                    break;
                case 3:
                    break;
                default:
                    break;
            }
            assertTrue(result.getFractionOfWhole() >= expFrac);
            */
        }
    }
    
    public void _testMatchLines() throws Exception {
        
        // close to correct, but one set of lines is interpreted as
        // 1 line instead of 2 due to threshold of consecutive points.
        //   so may need to make a PartialShapeMatcher2 specific
        //   to the task of matching a line...that should remove sensitivity
        //   to the model line length and the added corners to make a closed
        //   shape of lines...
        // use of this meanwhile, depends upon results of a null test
        // to not find lines where there are only curves 
        // so the resulting sum of chords is important for that last test

        PairFloatArray triangle = getTriangle();

        PairFloatArray rectangle = createRectangle(11, 6, 5, 5);
        rectangle = createLine(triangle.getN(), 5, 5);
        
        /*
        7                    *
        6                 *     *
        5              *           *
        4           *                 *
        3        *                       *
        2     *  *  *  *  *  *  *  *  *  *  *
        1
        0
           0  1  2  3  4  5  6  7  8  9 10 11
        */
      
        PartialShapeMatcher shapeMatcher = new PartialShapeMatcher();
        //matcher.setToDebug();
        shapeMatcher._overrideToThreshhold((float)(1e-7));
        shapeMatcher.overrideSamplingDistance(1);

        List<Match.Points> results = shapeMatcher.match(triangle, rectangle);
        assertFalse(results.isEmpty());
        plotResults(results, triangle, rectangle, 4,
                "_triangle_rectangle_", false);
    }

    public void testMatch() {
        PairFloatArray triangle = getTriangle();

        TransformationParameters params = new TransformationParameters();
        params.setOriginX(1);
        params.setOriginY(2);
        params.setTranslationX(10);
        params.setTranslationY(10);
        //params.setRotationInDegrees(90);
        //params.setScale(2);
        Transformer transformer = new Transformer();
        PairFloatArray triangle2 = transformer.applyTransformation(params, triangle);

        int minBlockSize = 5;
        int nMaxMatchable = Math.min(triangle.getN(), triangle2.getN());

        Match m = new Match(Math.max(triangle.getN(), triangle2.getN()));
        int blockSize = 5;
        int offset2 = 1;
        for (int i = 0, count=0; i < triangle.getN(); i += blockSize, ++count) {
            m.add(i, offset2, blockSize, 0.1, count);
        }

        Match m2 = m.copy();
        assertTrue(m.equals(m2));

        double maxChordSum = m.diffChordSum;
        m.maxChordSum = maxChordSum;
        m2.maxChordSum = maxChordSum;

        // reduce the differences of m2, so comparison will prefer m2
        assertEquals(0, m.compareTo(m2));
        m2.diffChordSum *= 0.9;
        assertTrue(m.compareTo(m2) > 0);

        Match.Points points = new Match.Points(m);
        assertEquals(triangle.getN(), points.pIdxs.length);
        assertEquals(triangle.getN(), points.qIdxs.length);
        assertEquals(triangle.getN(), points.mLen);
        for (int i = 0; i < points.pIdxs.length; ++i) {
            assertEquals((points.pIdxs[i] + offset2) % triangle.getN(), points.qIdxs[i]);
        }
        assertTrue(Math.abs(4*0.1 - points.chordDiffSum) < 1E-10);

        points.interchange();
        assertEquals(triangle.getN(), points.qIdxs.length);
        assertEquals(triangle.getN(), points.pIdxs.length);
        assertEquals(triangle.getN(), points.mLen);
        for (int i = 0; i < points.pIdxs.length; ++i) {
            assertEquals((points.qIdxs[i] + offset2) % triangle2.getN(), points.pIdxs[i]);
        }
        assertTrue(Math.abs(4*0.1 - points.chordDiffSum) < 1E-10);


    }

    public void testMatchTriangles() throws Exception {

        // close to correct, but one set of lines is interpreted as
        // 1 line instead of 2 due to threshold of consecutive points.
        //   so may need to make a PartialShapeMatcher2 specific
        //   to the task of matching a line...that should remove sensitivity
        //   to the model line length and the added corners to make a closed
        //   shape of lines...
        // use of this meanwhile, depends upon results of a null test
        // to not find lines where there are only curves
        // so the resulting sum of chords is important for that last test

        PairFloatArray triangle = getTriangle();

        TransformationParameters params = new TransformationParameters();
        params.setOriginX(1);
        params.setOriginY(2);
        params.setTranslationX(10);
        params.setTranslationY(10);
        params.setRotationInDegrees(90);
        params.setScale(2);
        Transformer transformer = new Transformer();
        PairFloatArray triangle2 = transformer.applyTransformation(params, triangle);

        //System.out.printf("plotting %s\n", plot(triangle, 1000));
        //System.out.printf("plotting %s\n", plot(triangle2, 1001));

        /*
        7                    *
        6                 *     *
        5              *           *
        4           *                 *
        3        *                       *
        2     *  *  *  *  *  *  *  *  *  *  *
        1
        0
           0  1  2  3  4  5  6  7  8  9 10 11
        */

        PartialShapeMatcher matcher = new PartialShapeMatcher();
        //matcher.setToDebug();
        matcher._overrideToThreshhold((float)(1e-7));
        matcher.overrideSamplingDistance(1);

        float[][] a1 = matcher.createDescriptorMatrix(triangle, triangle.getN());
        float[][] a2 = matcher.createDescriptorMatrix(triangle2, triangle2.getN());
        double tol = 1E-5;
        for (int i = 0; i < a1.length; ++i){
            for (int j = 0; j < a1[i].length; ++j) {
                float _a1 = a1[i][j];
                float _a2 = a2[i][j];
                float v1 = Math.abs(_a1 - _a2);
                float v2 = Math.abs(Math.abs(_a1 - _a2) - PartialShapeMatcher.TWO_PI);
                assert(v1 < tol || v2 < tol);
            }
        }
        float[][][] md = matcher.createDifferenceMatrices(triangle, triangle2);
        matcher.applySummedAreaTableConversion(md[0]);

        List<Match.Points> results = matcher.match(triangle, triangle2);
        assertFalse(results.isEmpty());
        plotResults(results, triangle, triangle2, 4,
                "_triangle_triangle2_", true);
    }

    public void testMatchTrianglesSameNumberPoints() throws Exception {
        PairFloatArray triangle1 = getTriangle(7, 2, 1);
        PairFloatArray triangle2 = getTriangle(9, 2, 1);

        /*
        TransformationParameters params = new TransformationParameters();
        params.setOriginX(1);
        params.setOriginY(2);
        params.setTranslationX(10);
        params.setTranslationY(10);
        params.setRotationInDegrees(90);
        params.setScale(2);
        Transformer transformer = new Transformer();
        triangle2 = transformer.applyTransformation(params, triangle2);
        */

        PartialShapeMatcher shapeMatcher = new PartialShapeMatcher();
        //shapeMatcher.setToDebug();
        shapeMatcher.overrideSamplingDistance(1);
        shapeMatcher.setToUseSameNumberOfPoints();

        List<Match.Points> results = shapeMatcher.match(triangle1, triangle2);
        assertFalse(results.isEmpty());
        plotResults(results, triangle1, triangle2, 4,
                "_triangle1_triangle2_", true);
    }
    
    protected PairFloatArray getTriangle() {
        int top = 7;
        int base = 2;
        // slope is +1
        int leftX = 1;
        return getTriangle(top, base, leftX);
    }

    protected PairFloatArray getTriangle(int topY, int baseY, int leftX) {
        /*
        7                    *
        6                 *     *
        5              *           *
        4           *                 *
        3        *                       *
        2     *  *  *  *  *  *  *  *  *  *  *
        1
        0
           0  1  2  3  4  5  6  7  8  9 10 11
        */
        PairFloatArray p = new PairFloatArray();
        int i = leftX;
        int j = baseY;
        while (j <= topY) {
            p.add(i, j);
            ++i;
            ++j;
        }
        j -= 2;
        while (j >= baseY) {
            p.add(i, j);
            ++i;
            --j;
        }
        i -= 2;
        ++j;
        while (i > leftX) {
            p.add(i, j);
            --i;
        }

        return p;
    }
    
    protected PairFloatArray createLine(int len, int xOff, int yOff) {
        PairFloatArray a = new PairFloatArray(len);
        for (int i = (len - 1); i >= 0; --i) {
            a.add(xOff + i, yOff + i);
        }
        return a;
    }
    
    protected PairFloatArray createRectangle(int width, int height, int xOff,
        int yOff) {
        /*
          h
        
          0     w
        */
        PairFloatArray a = new PairFloatArray(2*width + 2*height);
        for (int i = (width - 1); i >= 0; --i) {
            a.add(i, 0);
        }
        for (int i = 1; i < height; ++i) {
            a.add(0, i);
        }
        for (int i = 1; i < width; ++i) {
            a.add(i, height - 1);
        }
        for (int i = (height - 2); i >= 1; --i) {
            a.add(width - 1, i);
        }
        for (int i = 0; i < a.getN(); ++i) {
            a.set(i, a.getX(i) + xOff, a.getY(i) + yOff);
        }
        
        return a;
    }

    protected PairFloatArray getScissors1() {

        PairFloatArray p = new PairFloatArray();
        p.add(95, 55);
        p.add(102, 52);
        p.add(108, 42);
        p.add(110, 35);
        p.add(118, 28);
        p.add(128, 25);
        p.add(135, 25);
        p.add(142, 27);
        p.add(152, 32);    //6, 7, 8
        p.add(150, 40);
        p.add(143, 48);
        p.add(136, 52);    //9, 10, 11
        p.add(128, 56);
        p.add(119, 59);
        p.add(110, 60);    // 12
        p.add(102, 62);
        p.add(95, 65);
        p.add(100, 72);
        p.add(108, 80);
        p.add(115, 87);
        p.add(122, 95);      // 18,19,20
        p.add(128, 104);
        p.add(130, 111);
        p.add(131, 121);   // 21,22, 23
        p.add(124, 122);
        p.add(117, 124);
        p.add(108, 121);
        p.add(100, 112);
        p.add(95, 105);
        p.add(94, 98);      // 27
        p.add(96, 88);
        p.add(95, 81);
        p.add(89, 72);
        p.add(80, 70);
        p.add(72, 70);
        p.add(63, 70);    //33, 34, 35
        p.add(57, 70);
        p.add(49, 70);
        p.add(39, 70);
        p.add(32, 70);
        p.add(24, 67);
        p.add(20, 64);    // 39, 40, 41
        p.add(28, 63);
        p.add(37, 63);
        p.add(45, 62);    // 42
        p.add(53, 61);
        p.add(60, 60);
        p.add(70, 60);    // 45
        p.add(62, 55);
        p.add(53, 51);
        p.add(45, 45);    // 48
        p.add(38, 40);
        p.add(29, 37);
        p.add(30, 33);    // 51
        p.add(38, 34);
        p.add(47, 36);
        p.add(54, 40);    // 54
        p.add(62, 44);
        p.add(70, 49);
        p.add(78, 54);
        p.add(87, 58);

        // accidently entered y for x, so
        // reverse them here
        for (int i = 0; i < p.getN(); ++i) {
            float x = p.getX(i);
            float y = p.getY(i);
            p.set(i, y, x);
        }

        return p;
    }

    protected PairFloatArray getScissors2() {

        PairFloatArray p = new PairFloatArray();
        p.add(80, 105);
        p.add(75, 112);
        p.add(65, 117);
        p.add(58, 121);
        p.add(52, 130);
        p.add(51, 138);
        p.add(50, 147);
        p.add(55, 155);  // 6,7
        p.add(60, 160);
        p.add(67, 156);
        p.add(74, 149);
        p.add(78, 141);
        p.add(82, 133);
        p.add(85, 125);
        p.add(85, 116);
        p.add(89, 109);
        p.add(92, 100);
        p.add(100, 94);   // 16,17
        p.add(109, 91);
        p.add(117, 87);
        p.add(125, 85);
        p.add(132, 82);   // 20,21
        p.add(141, 78);
        p.add(149, 71);
        p.add(155, 63);
        p.add(150, 59);   // 24,25
        p.add(142, 52);
        p.add(135, 52);
        p.add(127, 53);
        p.add(118, 58);
        p.add(112, 66);
        p.add(110, 74);   // 30,31
        p.add(105, 82);
        p.add(97, 87);
        p.add(97, 77);
        p.add(97, 69);
        p.add(97, 61);
        p.add(97, 53);
        p.add(97, 45);
        p.add(95, 36);    // 38, 39
        p.add(92, 28);
        p.add(90, 30);
        p.add(89, 38);
        p.add(89, 47);
        p.add(88, 55);
        p.add(87, 62);    // 44, 45
        p.add(87, 72);
        p.add(85, 79);    // 46, 47
        p.add(84, 88);
        p.add(78, 90);
        p.add(70, 92);
        p.add(62, 93);    // 50, 51
        p.add(52, 94);
        p.add(45, 96);
        p.add(38, 97);
        p.add(30, 98);
        p.add(20, 100);
        p.add(29, 104);
        p.add(38, 107);
        p.add(45, 107);  // 58,59
        p.add(52, 106);
        p.add(62, 104);
        p.add(69, 103);
        p.add(77, 102);

        return p;
    }

    protected PairFloatArray getWineGlassShape() {

        PairFloatArray p = new PairFloatArray();
        p.add(140, 500 - 450);
        p.add(220, 500 - 410);
        p.add(225, 500 - 390);
        p.add(225, 500 - 325);
        p.add(190, 500 - 280);
        p.add(130, 500 - 240);
        p.add(130, 500 - 180);
        p.add(200, 500 - 180);
        p.add(250, 500 - 180);
        p.add(310, 500 - 180);
        p.add(380, 500 - 200);
        p.add(330, 500 - 260);
        p.add(280, 500 - 310);
        p.add(275, 500 - 380);
        p.add(310, 500 - 440);
        p.add(375, 500 - 460);
        p.add(320, 500 - 475);
        p.add(275, 500 - 475);
        p.add(200, 500 - 475);
        return p;
    }

    private String plot(PairFloatArray p, int fn) throws Exception {

        float[] x = Arrays.copyOf(p.getX(), p.getN());
        float[] y = Arrays.copyOf(p.getY(), p.getN());
        float xMax = MiscMath.findMax(x) + 1;
        float yMax = MiscMath.findMax(y) + 1;

        PolygonAndPointPlotter plot = new PolygonAndPointPlotter();

        plot.addPlot(0, xMax, 0, yMax,
            x, y, x, y, "");

        return plot.writeFile(fn);
    }

    private PairFloatArray extractOrderedBoundary(ImageExt image) {
        return extractOrderedBoundary(image, SIGMA.ONE);
    }

    private PairFloatArray extractOrderedBoundary(ImageExt image,
        SIGMA sigma) {

        GreyscaleImage img = image.copyToGreyscale();

        Set<PairInt> blob = new HashSet<PairInt>();
        for (int i = 0; i < img.getNPixels(); ++i) {
            if (img.getValue(i) > 0) {
                int x = img.getCol(i);
                int y = img.getRow(i);
                blob.add(new PairInt(x, y));
            }
        }

        ImageProcessor imageProcessor =
            new ImageProcessor();

        PairIntArray ordered =
            imageProcessor.extractSmoothedOrderedBoundary(
            blob, sigma, img.getWidth(), img.getHeight());

        PairFloatArray out = new PairFloatArray();
        for (int i = 0; i < ordered.getN(); ++i) {
            out.add(ordered.getX(i), ordered.getY(i));
        }

        return out;
    }

}
