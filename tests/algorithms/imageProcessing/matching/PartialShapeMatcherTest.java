package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.SIGMA;
import algorithms.misc.MiscMath;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
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

    public void testSummedColumnTables() {

        /*
        2 9  2  7     9 11  18
        1 5  1  2     5  6   8 
        0 2  3  5     2  5  10 
          0  1  2
         */
        float[][] a = new float[3][];
        a[0] = new float[]{2, 3, 5};
        a[1] = new float[]{5, 1, 2};
        a[2] = new float[]{9, 2, 7};

        PartialShapeMatcher matcher = new PartialShapeMatcher();
        matcher.applySummedColumnTableConversion(a);

        assertTrue(Arrays.equals(new float[]{2, 5, 10}, a[0]));
        assertTrue(Arrays.equals(new float[]{5, 6, 8}, a[1]));
        assertTrue(Arrays.equals(new float[]{9, 11, 18}, a[2]));

    }

    public void testScissorsMatch0() throws Exception {

        // 60
        PairIntArray p = getScissors1();
        //plot(p, 200);

        // 63
        PairIntArray q = getScissors2();
        //plot(q, 201);

        log.info("p.n=" + p.getN() + " q.n=" + q.getN());

        //q.rotateLeft(q.getN() - 3);
        PartialShapeMatcher shapeMatcher = new PartialShapeMatcher();
        shapeMatcher.overrideSamplingDistance(1);
        shapeMatcher.setToDebug();

        PartialShapeMatcher.Result result = shapeMatcher.match(p, q);

        assertNotNull(result);

        log.info("RESULTS= scissors offset0: "
            + result.toString());

        CorrespondencePlotter plotter = new CorrespondencePlotter(p, q);

        for (int ii = 0; ii < result.getNumberOfMatches(); ++ii) {
            int idx1 = result.getIdx1(ii);
            int idx2 = result.getIdx2(ii);
            int x1 = p.getX(idx1);
            int y1 = p.getY(idx1);
            int x2 = q.getX(idx2);
            int y2 = q.getY(idx2);
            //System.out.println(String.format(
            //"(%d, %d) <=> (%d, %d)", x1, y1, x2, y2));

            if ((ii % 4) == 0) {
                plotter.drawLineInAlternatingColors(x1, y1, x2, y2,
                    0);
            }
        }
        String filePath = plotter.writeImage("_"
            + "_scissors_offset0_corres");

        if (enableAsserts) {
            assertTrue(result.getFractionOfWhole() > 0.4);
        }
    }

    public void testScissorsMatch16() throws Exception {
        
        // rotate points p so that start points are
        // different and assert that wrap around is
        // handled correctly
        PairIntArray p = getScissors1();
        p.rotateLeft(16);
        plot(p, 200);

        PairIntArray q = getScissors2();
        plot(q, 201);

        log.info("p.n=" + p.getN() + " q.n=" + q.getN());

        //q.rotateLeft(q.getN() - 3);
        PartialShapeMatcher shapeMatcher = new PartialShapeMatcher();
        shapeMatcher.overrideSamplingDistance(1);
        //shapeMatcher.setToDebug();
        
        // articulated:
        PartialShapeMatcher.Result result = shapeMatcher.match(p, q);

        assertNotNull(result);

        log.info("RESULTS= scissors offset15: "
            + result.toString());

        CorrespondencePlotter plotter = new CorrespondencePlotter(p, q);

        for (int ii = 0; ii < result.getNumberOfMatches(); ++ii) {
            int idx1 = result.getIdx1(ii);
            int idx2 = result.getIdx2(ii);
            int x1 = p.getX(idx1);
            int y1 = p.getY(idx1);
            int x2 = q.getX(idx2);
            int y2 = q.getY(idx2);
            //System.out.println(String.format(
            //"(%d, %d) <=> (%d, %d)", x1, y1, x2, y2));

            if ((ii % 4) == 0) {
                plotter.drawLineInAlternatingColors(x1, y1, x2, y2,
                    0);
            }
        }
        String filePath = plotter.writeImage("_"
            + "_scissors_offset16_corres");

        if (enableAsserts) {
            assertTrue(result.getFractionOfWhole() > 0.3);
        }
    }
    
    public void testScissorsMatch16_scaled() throws Exception {
        
        MiscellaneousCurveHelper curveHelper = 
            new MiscellaneousCurveHelper();
        
        // rotate points p so that start points are
        // different and assert that wrap around is
        // handled correctly
        PairIntArray p = getScissors1();
        p.rotateLeft(16);
        p = curveHelper.scaleDown(p, 0.5f);
        plot(p, 200);

        PairIntArray q = getScissors2();
        //q = curveHelper.scaleDown(q, 0.5f);
        plot(q, 201);
       
        log.info("p.n=" + p.getN() + " q.n=" + q.getN());

        //q.rotateLeft(q.getN() - 3);
        PartialShapeMatcher shapeMatcher = new PartialShapeMatcher();
        shapeMatcher.overrideSamplingDistance(1);
        
        shapeMatcher.overrideMinimumLength(4);
        shapeMatcher.setToUseSameNumberOfPoints();
        //shapeMatcher.setToDebug();
        
        // articulated:
        PartialShapeMatcher.Result result = shapeMatcher.match(p, q);

        assertNotNull(result);

        log.info("RESULTS= scissors offset15: "
            + result.toString());

        CorrespondencePlotter plotter = new CorrespondencePlotter(p, q);

        for (int ii = 0; ii < result.getNumberOfMatches(); ++ii) {
            int idx1 = result.getIdx1(ii);
            int idx2 = result.getIdx2(ii);
            int x1 = p.getX(idx1);
            int y1 = p.getY(idx1);
            int x2 = q.getX(idx2);
            int y2 = q.getY(idx2);
            //System.out.println(String.format(
            //"(%d, %d) <=> (%d, %d)", x1, y1, x2, y2));

            if ((ii % 4) == 0) {
                plotter.drawLineInAlternatingColors(x1, y1, x2, y2,
                    0);
            }
        }
        String filePath = plotter.writeImage("_"
            + "_scissors_offset16_corres");

        if (enableAsserts) {
            assertTrue(result.getFractionOfWhole() > 0.3);
        }
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

        PairIntArray p = extractOrderedBoundary(img0);
        plot(p, 100);

        String fileName1 = "";

        for (int type = 0; type < 2; ++type) {
            for (int i = 0; i < 4; ++i) {
                
                int dp = 2;

                PartialShapeMatcher matcher =
                    new PartialShapeMatcher();
                //matcher.setToDebug();
                if (type == 0) {
                    matcher._overrideToThreshhold(0.2f);
                } else if (type == 1) {
                    matcher.setToUseEuclidean();
                }
                matcher.overrideSamplingDistance(dp);
                matcher.setToRemoveOutliers();
                
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

                PairIntArray q = extractOrderedBoundary(img);
                plot(q, (i+1)*100 + 1);

                log.info("matching " + fileName0Root
                + " to " + fileName1Root + " (" + p.getN()
                + " points to " + q.getN() + " points");

               
                PartialShapeMatcher.Result result = matcher.match(p, q);

                assertNotNull(result);

                log.info("RESULTS=" + fileName1Root + " : " +
                    result.toString());

                CorrespondencePlotter plotter = new
                    CorrespondencePlotter(p, q);

                for (int ii = 0; ii < result.getNumberOfMatches(); ++ii) {
                    int idx1 = result.getIdx1(ii);
                    int idx2 = result.getIdx2(ii);
                    int x1 = p.getX(idx1);
                    int y1 = p.getY(idx1);
                    int x2 = q.getX(idx2);
                    int y2 = q.getY(idx2);
                    //System.out.println(String.format(
                    //"(%d, %d) <=> (%d, %d)", x1, y1, x2, y2));

                    if ((ii % 2) == 0) {
                        plotter.drawLineInAlternatingColors(x1, y1, x2, y2,
                            0);
                    }
                }
                String filePath = plotter.writeImage("_" +
                        fileName1Root + "_corres_" + type);
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

        PairIntArray p = extractOrderedBoundary(img0, sigma);
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

            PairIntArray q = extractOrderedBoundary(img, sigma);
            plot(q, (i+1)*100 + 1);

            log.info("matching " + fileName0Root
            + " to " + fileName1Root + " (" + p.getN()
            + " points to " + q.getN() + " points");

            int dp = 1;
            PartialShapeMatcher matcher = new PartialShapeMatcher();
            //matcher.setToDebug();
            matcher.setToUseSameNumberOfPoints();
            matcher.overrideSamplingDistance(dp);
            //matcher._overrideToThreshhold(0.2f);
            //matcher.setToRemoveOutliers();
            matcher.setToUseEuclidean();
            
            PartialShapeMatcher.Result result = matcher.match(p, q);

            assertNotNull(result);

            log.info("RESULTS=" + fileName1Root + " : " +
                result.toString());

            CorrespondencePlotter plotter = new
                CorrespondencePlotter(p, q);

            for (int ii = 0; ii < result.getNumberOfMatches(); ++ii) {
                int idx1 = result.getIdx1(ii);
                int idx2 = result.getIdx2(ii);
                int x1 = p.getX(idx1);
                int y1 = p.getY(idx1);
                int x2 = q.getX(idx2);
                int y2 = q.getY(idx2);
                //System.out.println(String.format(
                //"(%d, %d) <=> (%d, %d)", x1, y1, x2, y2));

                if ((ii % 4) == 0) {
                    plotter.drawLineInAlternatingColors(x1, y1, x2, y2,
                        0);
                }
            }
            String filePath = plotter.writeImage("_" +
                    fileName1Root + "_corres");

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
    
    public void testMatchLines() {
        
        // close to correct, but one set of lines is interpreted as
        // 1 line instead of 2 due to threshold of consecutive points.
        //   so may need to make a PartialShapeMatcher specific
        //   to the task of matching a line...that should remove sensitivity
        //   to the model line length and the added corners to make a closed
        //   shape of lines...
        // use of this meanwhile, depends upon results of a null test
        // to not find lines where there are only curves 
        // so the resulting sum of chords is important for that last test
        
        PairIntArray triangle = getTriangle();
        
        PairIntArray rectangle = createRectangle(11, 6, 5, 5);
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
      
        PartialShapeMatcher matcher = new PartialShapeMatcher();
        //matcher.setToDebug();
        matcher._overrideToThreshhold((float)(1e-7));
        matcher.overrideSamplingDistance(1);
        
        PartialShapeMatcher.Result r = matcher.match(triangle, rectangle);
        for (int i = 0; i < r.idx1s.size(); ++i) {
            int x1 = triangle.getX(r.idx1s.get(i)); 
            int y1 = triangle.getY(r.idx1s.get(i)); 
            int x2 = rectangle.getX(r.idx2s.get(i)); 
            int y2 = rectangle.getY(r.idx2s.get(i)); 
            //int segIdx = r.getArticulatedSegment(i);
            System.out.println(x1 + ", " + y1 + "   " + x2 + ", " + y2 
                //+ " segIdx=" + segIdx 
                + " idx1=" + r.idx1s.get(i)
                + " idx2=" + r.idx2s.get(i)
            );
        }
        System.out.println("triangle size=" + triangle.getN() +
            " matched size=" + r.getNumberOfMatches());
    }
    
    protected PairIntArray getTriangle() {
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
        
        PairIntArray p = new PairIntArray(20);
        for (int i = 1; i <= 6; ++i) {
            p.add(i, i + 1); 
        }
        p.add(7, 6); p.add(8, 5);  p.add(9, 4);  p.add(10, 3); 
        for (int i = 11; i >= 2; --i) {
            p.add(i, 2);
        }
        
        return p;
    }
    
    protected PairIntArray createLine(int len, int xOff, int yOff) {
        PairIntArray a = new PairIntArray(len);
        for (int i = (len - 1); i >= 0; --i) {
            a.add(xOff + i, yOff + i);
        }
        return a;
    }
    
    protected PairIntArray createRectangle(int width, int height, int xOff, 
        int yOff) {
        /*
          h
        
          0     w
        */
        PairIntArray a = new PairIntArray(2*width + 2*height);
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

    protected PairIntArray getScissors1() {

        PairIntArray p = new PairIntArray();
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
            int x = p.getX(i);
            int y = p.getY(i);
            p.set(i, y, x);
        }

        return p;
    }

    protected PairIntArray getScissors2() {

        PairIntArray p = new PairIntArray();
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

    protected PairIntArray getWineGlassShape() {

        PairIntArray p = new PairIntArray();
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

    private void plot(PairIntArray p, int fn) throws Exception {

        float[] x = new float[p.getN()];
        float[] y = new float[p.getN()];

        for (int i = 0; i < x.length; ++i) {
            x[i] = p.getX(i);
            y[i] = p.getY(i);
        }

        float xMax = MiscMath.findMax(x) + 1;
        float yMax = MiscMath.findMax(y) + 1;

        PolygonAndPointPlotter plot = new PolygonAndPointPlotter();

        plot.addPlot(0, xMax, 0, yMax,
            x, y, x, y, "");

        plot.writeFile(fn);
    }

    private PairIntArray extractOrderedBoundary(ImageExt image) {
        return extractOrderedBoundary(image, SIGMA.ONE);
    }

    private PairIntArray extractOrderedBoundary(ImageExt image,
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

        return ordered;
    }

}
