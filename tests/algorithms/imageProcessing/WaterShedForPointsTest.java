package algorithms.imageProcessing;

import algorithms.graphs.CustomWatershedDAG;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class WaterShedForPointsTest extends TestCase {

    public WaterShedForPointsTest() {
    }

    public void testLower_set() throws Exception {

        String filePath = ResourceFinder.findFileInTestResources("pattern0.png");
        Image img = ImageIOHelper.readImageAsGrayScale(filePath);
        GreyscaleImage img0 = img.copyToGreyscale();

        int w = img0.getWidth();
        int h = img0.getHeight();

        int factor = 2;//100
        GreyscaleImage img09 = new GreyscaleImage(w*factor, h*factor);
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                int v = img0.getValue(i, j) + 2;
                for (int ii = (i*factor); ii < ((i + 1)*factor); ii++) {
                    for (int jj = (j*factor); jj < ((j + 1)*factor); jj++) {
                        img09.setValue(ii, jj, v);
                    }
                }
            }
        }

        Set<PairInt> points = new HashSet<PairInt>();
        for (int i = 0; i < img09.getWidth(); ++i) {
            for (int j = 0; j < img09.getHeight(); ++j) {
                points.add(new PairInt(i, j));
            }
        }

        WaterShedForPoints ws = new WaterShedForPoints();

        Map<PairInt, Integer> lowerComplete = ws.lower(img09, points);

        assertNotNull(lowerComplete);

        GreyscaleImage imgL = new GreyscaleImage(img09.getWidth(), img09.getHeight());
        for (Entry<PairInt, Integer> entry : lowerComplete.entrySet()) {
            PairInt p = entry.getKey();
            int v = entry.getValue().intValue();
            imgL.setValue(p.getX(), p.getY(), v);
        }

        String bin = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(bin + "/pattern0_input.png", img09);
        ImageIOHelper.writeOutputImage(bin + "/pattern0_lowered.png", imgL);

        StringBuilder sb = new StringBuilder("input:\n");
        for (int j = 0; j < img09.getHeight(); ++j) {
            sb.append(String.format("row %4d:", j));
            for (int i = 0; i < img09.getWidth(); ++i) {
                sb.append(String.format("%4d", img09.getValue(i, j)));
            }
            sb.append("\n");
        }
        System.out.println(sb.toString());

        Map<PairInt, Integer> distToLowerIntPix = ws.getDistToLowerIntensityPixel();
        sb = new StringBuilder("distToLowerIntPix:\n");
        for (int j = 0; j < img09.getHeight(); ++j) {
            sb.append(String.format("row %4d:", j));
            for (int i = 0; i < img09.getWidth(); ++i) {
                PairInt p = new PairInt(i, j);
                Integer value = distToLowerIntPix.get(p);
                int v = (value == null) ? 0 : value.intValue();
                sb.append(String.format("%4d", v));
            }
            sb.append("\n");
        }
        System.out.println(sb.toString());

        StringBuilder sb2 = new StringBuilder("lower complete:");
        for (int j = 0; j < img09.getHeight(); ++j) {
            sb2.append(String.format("row %4d:", j));
            for (int i = 0; i < img09.getWidth(); ++i) {
                PairInt p = new PairInt(i, j);
                Integer value = lowerComplete.get(p);
                int v = (value == null) ? 0 : value.intValue();
                sb2.append(String.format("%4d", v));
            }
            sb2.append("\n");
        }
        System.out.println(sb2.toString());

        // --------------- assert regionalMinima ---------------
        final Set<PairInt> regionalMinima = ws.getRegionalMinima();
        // fudging the test since already know that all of the minima have the
        // same value.
        int minValue = MiscMath.findMin(img09, points);
        Set<PairInt> expected = new HashSet<PairInt>(regionalMinima);
        for (PairInt p : regionalMinima) {
            int v = img09.getValue(p.getX(), p.getY());
            assertTrue(v == minValue);
            assertTrue(expected.remove(p));
        }
        assertTrue(expected.isEmpty());

        // ------------- assert values in distToLowerIntensityPixel -----------
        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;

        Map<PairInt, Integer> distToLowerIntensityPixel = ws.getDistToLowerIntensityPixel();

        for (PairInt p : points) {

            int x = p.getX();
            int y = p.getY();

            Set<PairInt> visited = new HashSet<PairInt>();

            int dist = findShortestPathToLowerIntensity(img09, visited, p, points, regionalMinima);

            assertTrue(distToLowerIntensityPixel.get(p).intValue() == dist);
        }

        // ------- assert final result, lowerComplete:
        //         that each pixel, if not in a regional minima, has an adjacent
        //         lower intensity pixel
        for (PairInt p : points) {
            if (regionalMinima.contains(p)) {
                continue;
            }
            Integer value = lowerComplete.get(p);
            assertNotNull(value);
            boolean foundLowerIntNghbr = false;
            for (int vIdx = 0; vIdx < dxs8.length; ++vIdx) {
                int x2 = p.getX() + dxs8[vIdx];
                int y2 = p.getY() + dys8[vIdx];
                PairInt p2 = new PairInt(x2, y2);
                if (points.contains(p2)) {
                    Integer value2 = lowerComplete.get(p2);
                    if (value2.intValue() < value.intValue()) {
                        foundLowerIntNghbr = true;
                        break;
                    }
                }
            }
            assertTrue(foundLowerIntNghbr);
        }
    }

    public void testUnionFindComponentLabelling_set() throws Exception {

        int w = 5;
        int h = 5;

        /*
        input data:

        4 3 0 3 4  0
        3 2 0 2 3  1
        0 0 0 0 0  2
        3 2 0 2 3  3
        4 3 0 3 4  4
                   5
        0 1 2 3 4

        */

        int[][] lowerComplete = new int[w][];
        for (int i = 0; i < w; ++i) {
            lowerComplete[i] = new int[h];
        }
        lowerComplete[0] = new int[]{4, 3, 0, 3, 4};
        lowerComplete[1] = new int[]{3, 2, 0, 2, 3};
        lowerComplete[2] = new int[]{0, 0, 0, 0, 0};
        lowerComplete[3] = new int[]{3, 2, 0, 2, 3};
        lowerComplete[4] = new int[]{4, 3, 0, 3, 4};

        Map<PairInt, Integer> lowerCompleteIm = new HashMap<PairInt, Integer>();
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                PairInt p = new PairInt(i, j);
                int v = lowerComplete[i][j];
                lowerCompleteIm.put(p, Integer.valueOf(v));
            }
        }

        WaterShedForPoints ws = new WaterShedForPoints();

        Map<PairInt, Integer> labeled = ws.unionFindComponentLabelling(lowerCompleteIm);

        assertTrue(labeled.size() == lowerCompleteIm.size());

        StringBuilder sb2 = new StringBuilder("input:\n");
        for (int j = 0; j < h; ++j) {
            sb2.append(String.format("row %4d:", j));
            for (int i = 0; i < w; ++i) {
                PairInt p = new PairInt(i, j);
                int v = lowerCompleteIm.get(p).intValue();
                sb2.append(String.format("%4d", v));
            }
            sb2.append("\n");
        }
        System.out.println(sb2.toString());

        sb2 = new StringBuilder("labeled components:\n");
        for (int j = 0; j < w; ++j) {
            sb2.append(String.format("row %4d:", j));
            for (int i = 0; i < h; ++i) {
                PairInt p = new PairInt(i, j);
                int v = labeled.get(p).intValue();
                sb2.append(String.format("%4d", v));
            }
            sb2.append("\n");
        }
        System.out.println(sb2.toString());

        /*
        input data:
            4 3 0 3 4  0
            3 2 0 2 3  1
            0 0 0 0 0  2
            3 2 0 2 3  3
            4 3 0 3 4  4
                       5
            0 1 2 3 4

        reversed componentLabelMap with values parent printed as parent
        and all keys with same parent listed as children:
            parent: (0,3)    children:  (0,3) (1,4)
            parent: (1,1)    children:  (1,1)
            parent: (4,0)    children:  (4,0)
            parent: (1,3)    children:  (1,3)
            parent: (0,1)    children:  (0,1) (1,0)
            parent: (3,0)    children:  (3,0) (4,1)
            parent: (3,1)    children:  (3,1)
            parent: (0,2)    children:  (2,3) (0,2) (2,0) (2,4) (3,2) (4,2) (2,1) (1,2) (2,2)
            parent: (0,4)    children:  (0,4)
            parent: (4,4)    children:  (4,4)
            parent: (3,3)    children:  (3,3)
            parent: (3,4)    children:  (4,3) (3,4)
            parent: (0,0)    children:  (0,0)

        output component labeled image:
            row    0:   1   2   3   8  12
            row    1:   2   6   3   9   8
            row    2:   3   3   3   3   3
            row    3:   4   7   3  10  11
            row    4:   5   4   3  11  13

                        0   1   2   3   4
        */
    }

    public void testCreateLowerIntensityDAG() throws Exception {

        int w = 5;
        int h = 5;

        /*
        from Fig. 11 of
        Roerdink and Meijster 2001
        "The Watershed Transform: Definitions, Algorithms and Parallelization Strategies",
        Fundamenta Informaticae 41 (2001) 187–228, Section 4.2.4

        input data:

        0 1 2 1 0  0
        1 2 3 2 1  1
        2 3 4 3 2  2
        1 2 3 2 1  3
        0 1 2 1 0  4

        0 1 2 3 4

        */

        int[][] im = new int[w][];
        for (int i = 0; i < w; ++i) {
            im[i] = new int[h];
        }
        im[0] = new int[]{0, 1, 2, 1, 0};
        im[1] = new int[]{1, 2, 3, 2, 1};
        im[2] = new int[]{2, 3, 4, 3, 2};
        im[3] = new int[]{1, 2, 3, 2, 1};
        im[4] = new int[]{0, 1, 2, 1, 0};

        Set<PairInt> points = new HashSet<PairInt>();
        GreyscaleImage img = new GreyscaleImage(w, h);
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                img.setValue(i, j, im[i][j]);
                points.add(new PairInt(i, j));
            }
        }

        WaterShedForPoints ws = new WaterShedForPoints();

        Map<PairInt, Integer> lowerComplete = ws.lower(img, points);

        assertNotNull(lowerComplete);

        assertTrue(lowerComplete.size() == points.size());

        StringBuilder sb2 = new StringBuilder("lowered:\n");
        for (int j = 0; j < h; ++j) {
            sb2.append(String.format("row %4d:", j));
            for (int i = 0; i < w; ++i) {
                PairInt p = new PairInt(i, j);
                Integer value = lowerComplete.get(p);
                sb2.append(String.format("%4d", value.intValue()));
            }
            sb2.append("\n");
        }
        System.out.println(sb2.toString());

        Map<PairInt, Integer> labelled = ws.unionFindComponentLabelling(lowerComplete);
        sb2 = new StringBuilder("labelled 0:\n");
        for (int j = 0; j < h; ++j) {
            sb2.append(String.format("row %4d:", j));
            for (int i = 0; i < w; ++i) {
                PairInt p = new PairInt(i, j);
                Integer value = labelled.get(p);
                sb2.append(String.format("%4d", value.intValue()));
            }
            sb2.append("\n");
        }
        System.out.println(sb2.toString());


        CustomWatershedDAG dag = ws.createLowerIntensityDAG(lowerComplete);

        // ----- row 0
        assertTrue(dag.getConnectedNumber(new PairInt(0, 0)) == 0);

        assertTrue(dag.getConnectedNumber(new PairInt(1, 0)) == 1);
        assertTrue(dag.getConnectedNode(new PairInt(1, 0), 0).equals(new PairInt(0, 0)));

        assertTrue(dag.getConnectedNumber(new PairInt(2, 0)) == 2);
        PairInt p0 = dag.getConnectedNode(new PairInt(2, 0), 0);
        PairInt p1 = dag.getConnectedNode(new PairInt(2, 0), 1);
        assertTrue(p0.equals(new PairInt(1, 0)) || p1.equals(new PairInt(1, 0)));
        if (p0.equals(new PairInt(1, 0))) {
            assertTrue(p1.equals(new PairInt(3, 0)));
        } else {
            assertTrue(p0.equals(new PairInt(3, 0)));
        }

        assertTrue(dag.getConnectedNumber(new PairInt(3, 0)) == 1);
        assertTrue(dag.getConnectedNode(new PairInt(3, 0), 0).equals(new PairInt(4, 0)));

        assertTrue(dag.getConnectedNumber(new PairInt(4, 0)) == 0);

        // ----- row 1
        assertTrue(dag.getConnectedNumber(new PairInt(0, 1)) == 1);
        assertTrue(dag.getConnectedNode(new PairInt(0, 1), 0).equals(new PairInt(0, 0)));

        assertTrue(dag.getConnectedNumber(new PairInt(1, 1)) == 3);
        p0 = dag.getConnectedNode(new PairInt(1, 1), 0);
        p1 = dag.getConnectedNode(new PairInt(1, 1), 1);
        PairInt p2 = dag.getConnectedNode(new PairInt(1, 1), 2);
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(0, 0));
        expected.add(new PairInt(0, 1));
        expected.add(new PairInt(1, 0));
        assertTrue(expected.remove(p0));
        assertTrue(expected.remove(p1));
        assertTrue(expected.remove(p2));
        assertTrue(expected.isEmpty());
        assertTrue(dag.getConnectedNode(new PairInt(1, 1), 0).equals(new PairInt(0, 0)));

        int nExpected = 5;
        Set<PairInt> found = new HashSet<PairInt>();
        PairInt key = new PairInt(2, 1);
        assertTrue(dag.getConnectedNumber(key) == nExpected);
        for (int i = 0; i < nExpected; ++i) {
            found.add(dag.getConnectedNode(key, i));
        }
        expected = new HashSet<PairInt>();
        expected.add(new PairInt(1, 0));
        expected.add(new PairInt(2, 0));
        expected.add(new PairInt(3, 0));
        expected.add(new PairInt(1, 1));
        expected.add(new PairInt(3, 1));
        for (PairInt p : found) {
            assertTrue(expected.remove(p));
        }
        assertTrue(expected.isEmpty());

        // ------ skip to center pixel
        nExpected = 8;
        found = new HashSet<PairInt>();
        key = new PairInt(2, 2);
        assertTrue(dag.getConnectedNumber(key) == nExpected);
        for (int i = 0; i < nExpected; ++i) {
            found.add(dag.getConnectedNode(key, i));
        }
        expected = new HashSet<PairInt>();
        expected.add(new PairInt(1, 1));
        expected.add(new PairInt(2, 1));
        expected.add(new PairInt(3, 1));
        expected.add(new PairInt(3, 2));
        expected.add(new PairInt(3, 3));
        expected.add(new PairInt(2, 3));
        expected.add(new PairInt(1, 3));
        expected.add(new PairInt(1, 2));
        for (PairInt p : found) {
            assertTrue(expected.remove(p));
        }
        assertTrue(expected.isEmpty());

        Map<PairInt, Integer> labelled2 = ws.unionFindWatershed(lowerComplete);

        assertNotNull(labelled2);

        sb2 = new StringBuilder("labelled final:\n");
        for (int j = 0; j < h; ++j) {
            sb2.append(String.format("row %4d:", j));
            for (int i = 0; i < w; ++i) {
                PairInt p = new PairInt(i, j);
                Integer value = labelled2.get(p);
                sb2.append(String.format("%4d", value.intValue()));
            }
            sb2.append("\n");
        }
        System.out.println(sb2.toString());

        // assert the watershed pixels
        assertTrue(labelled2.get(new PairInt(2,0)) == 0);
        assertTrue(labelled2.get(new PairInt(2, 1)) == 0);
        assertTrue(labelled2.get(new PairInt(0, 2)) == 0);
        assertTrue(labelled2.get(new PairInt(1, 2)) == 0);
        assertTrue(labelled2.get(new PairInt(2, 2)) == 0);
        assertTrue(labelled2.get(new PairInt(3, 2)) == 0);
        assertTrue(labelled2.get(new PairInt(4, 2)) == 0);
        assertTrue(labelled2.get(new PairInt(2, 3)) == 0);
        assertTrue(labelled2.get(new PairInt(2, 4)) == 0);

        assertTrue(labelled2.get(new PairInt(0, 0)) == labelled.get(new PairInt(0, 0)));
        assertTrue(labelled2.get(new PairInt(1, 0)) == labelled.get(new PairInt(0, 0)));
        assertTrue(labelled2.get(new PairInt(1, 1)) == labelled.get(new PairInt(0, 0)));
        assertTrue(labelled2.get(new PairInt(0, 1)) == labelled.get(new PairInt(0, 0)));

        assertTrue(labelled2.get(new PairInt(0, 4)) == labelled.get(new PairInt(0, 4)));
        assertTrue(labelled2.get(new PairInt(0, 3)) == labelled.get(new PairInt(0, 4)));
        assertTrue(labelled2.get(new PairInt(1, 3)) == labelled.get(new PairInt(0, 4)));
        assertTrue(labelled2.get(new PairInt(1, 4)) == labelled.get(new PairInt(0, 4)));

        assertTrue(labelled2.get(new PairInt(4, 4)) == labelled.get(new PairInt(4, 4)));
        assertTrue(labelled2.get(new PairInt(3, 4)) == labelled.get(new PairInt(4, 4)));
        assertTrue(labelled2.get(new PairInt(3, 3)) == labelled.get(new PairInt(4, 4)));
        assertTrue(labelled2.get(new PairInt(4, 3)) == labelled.get(new PairInt(4, 4)));

        assertTrue(labelled2.get(new PairInt(4, 0)) == labelled.get(new PairInt(4, 0)));
        assertTrue(labelled2.get(new PairInt(4, 1)) == labelled.get(new PairInt(4, 0)));
        assertTrue(labelled2.get(new PairInt(3, 0)) == labelled.get(new PairInt(4, 0)));
        assertTrue(labelled2.get(new PairInt(3, 1)) == labelled.get(new PairInt(4, 0)));

    }

    // ----------  tests of methods using point sets ------

    private int findShortestPathToLowerIntensity(GreyscaleImage img,
        Set<PairInt> visited, PairInt p, Set<PairInt> points, Set<PairInt> regionalMinima) {

        if (regionalMinima.contains(p)) {
            return 0;
        }

        int v = img.getValue(p.getX(), p.getY());

        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;

        int minDist = Integer.MAX_VALUE;

        for (int vIdx = 0; vIdx < dxs8.length; ++vIdx) {
            int x2 = p.getX() + dxs8[vIdx];
            int y2 = p.getY() + dys8[vIdx];
            PairInt p2 = new PairInt(x2, y2);

            if (visited.contains(p2)) {
                continue;
            }

            if (regionalMinima.contains(p2)) {
                return 1;
            }

            if (points.contains(p2)) {
                int v2 = img.getValue(x2, y2);
                if (v2 < v) {
                    return 1;
                } else {
                    Set<PairInt> visited2 = new HashSet<PairInt>(visited);
                    visited2.add(p);
                    visited2.add(p2);
                    int dist = 1 + findShortestPathToLowerIntensity(img, visited2,
                        p2, points, regionalMinima);
                    if (dist < minDist) {
                        minDist = dist;
                    }
                }
            }
        }
        return minDist;
    }

    public void testLabel() throws Exception {

        int w = 5;
        int h = 5;

        /*
        from Fig. 11 of
        Roerdink and Meijster 2001
        "The Watershed Transform: Definitions, Algorithms and Parallelization Strategies",
        Fundamenta Informaticae 41 (2001) 187–228, Section 4.2.4

        input data:

        0 1 2 1 0  0
        1 2 3 2 1  1
        2 3 4 3 2  2
        1 2 3 2 1  3
        0 1 2 1 0  4

        0 1 2 3 4

        */

        int[][] im = new int[w][];
        for (int i = 0; i < w; ++i) {
            im[i] = new int[h];
        }
        im[0] = new int[]{0, 1, 2, 1, 0};
        im[1] = new int[]{1, 2, 3, 2, 1};
        im[2] = new int[]{2, 3, 4, 3, 2};
        im[3] = new int[]{1, 2, 3, 2, 1};
        im[4] = new int[]{0, 1, 2, 1, 0};

        Set<PairInt> points = new HashSet<PairInt>();
        GreyscaleImage img = new GreyscaleImage(w, h);
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                img.setValue(i, j, im[i][j]);
                points.add(new PairInt(i, j));
            }
        }

        WaterShedForPoints ws = new WaterShedForPoints();

        Map<PairInt, Integer> labelled2 = ws.createLabelledImage(img, points);

        StringBuilder sb2 = new StringBuilder("labelled 0:\n");
        for (int j = 0; j < h; ++j) {
            sb2.append(String.format("row %4d:", j));
            for (int i = 0; i < w; ++i) {
                PairInt p = new PairInt(i, j);
                Integer value = labelled2.get(p);
                sb2.append(String.format("%4d", value.intValue()));
            }
            sb2.append("\n");
        }
        System.out.println(sb2.toString());

        // assert the watershed pixels
        assertTrue(labelled2.get(new PairInt(2,0)) == 0);
        assertTrue(labelled2.get(new PairInt(2, 1)) == 0);
        assertTrue(labelled2.get(new PairInt(0, 2)) == 0);
        assertTrue(labelled2.get(new PairInt(1, 2)) == 0);
        assertTrue(labelled2.get(new PairInt(2, 2)) == 0);
        assertTrue(labelled2.get(new PairInt(3, 2)) == 0);
        assertTrue(labelled2.get(new PairInt(4, 2)) == 0);
        assertTrue(labelled2.get(new PairInt(2, 3)) == 0);
        assertTrue(labelled2.get(new PairInt(2, 4)) == 0);

        assertTrue(labelled2.get(new PairInt(0, 0)) > 0);
        assertTrue(labelled2.get(new PairInt(1, 0)) > 0);
        assertTrue(labelled2.get(new PairInt(1, 1)) > 0);
        assertTrue(labelled2.get(new PairInt(0, 1)) > 0);

        assertTrue(labelled2.get(new PairInt(0, 4))  > 0);
        assertTrue(labelled2.get(new PairInt(0, 3))  > 0);
        assertTrue(labelled2.get(new PairInt(1, 3))  > 0);
        assertTrue(labelled2.get(new PairInt(1, 4))  > 0);

        assertTrue(labelled2.get(new PairInt(4, 4))  > 0);
        assertTrue(labelled2.get(new PairInt(3, 4))  > 0);
        assertTrue(labelled2.get(new PairInt(3, 3))  > 0);
        assertTrue(labelled2.get(new PairInt(4, 3))  > 0);

        assertTrue(labelled2.get(new PairInt(4, 0))  > 0);
        assertTrue(labelled2.get(new PairInt(4, 1))  > 0);
        assertTrue(labelled2.get(new PairInt(3, 0))  > 0);
        assertTrue(labelled2.get(new PairInt(3, 1))  > 0);

    }

    public void testOnImage() throws Exception {

        String[] fileRoots = new String[]{"v_blob_1","v_blob_10"};

        for (String fileRoot : fileRoots) {

            String fileName = fileRoot + ".png";

            String filePath = ResourceFinder.findFileInTestResources(fileName);
            Image img = ImageIOHelper.readImageAsGrayScale(filePath);
            GreyscaleImage img0 = img.copyToGreyscale();

            filePath = ResourceFinder.findFileInTestResources(fileRoot + "_mask.png");
            Image imgMask = ImageIOHelper.readImageAsGrayScale(filePath);

            int w = img0.getWidth();
            int h = img0.getHeight();

            Set<PairInt> points = new HashSet<PairInt>();
            for (int i = 0; i < w; ++i) {
                for (int j = 0; j < h; ++j) {
                    int v = imgMask.getRGB(i, j);
                    if (v == 0) {
                        points.add(new PairInt(i, j));
                    }
                }
            }

            // mask out the non-points before attempt to use log scaling:
            for (int i = 0; i < w; ++i) {
                for (int j = 0; j < h; ++j) {
                    PairInt p = new PairInt(i, j);
                    if (!points.contains(p)) {
                        img0.setValue(i, j, 0);
                    }
                }
            }

            String bin = ResourceFinder.findDirectory("bin");

            GreyscaleImage img0T = img0.copyImage();
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.applyAdaptiveMeanThresholding(img0T, 3);
            ImageIOHelper.writeOutputImage(bin + "/" + fileRoot + "_thr.png", img0T);

            WaterShedForPoints ws = new WaterShedForPoints();

            Map<PairInt, Integer> labelled2 = ws.createLabelledImage(img0, points);

            GreyscaleImage imgL = new GreyscaleImage(w, h);
            for (Entry<PairInt, Integer> entry : labelled2.entrySet()) {
                PairInt p = entry.getKey();
                int v = entry.getValue().intValue();
                imgL.setValue(p.getX(), p.getY(), v);
            }

            ImageIOHelper.writeOutputImage(bin + "/" + fileRoot + "_labelled.png", imgL);

            GreyscaleImage imgW = new GreyscaleImage(w, h);
            for (int i = 0; i < imgW.getNPixels(); ++i) {
                imgW.setValue(i, 255);
            }
            for (PairInt p : points) {
                int v = labelled2.get(p);
     //System.out.println(p.toString() + " v=" + v);
                if (v == 0) {
                    imgW.setValue(p.getX(), p.getY(), 0);
                }
            }

            ImageIOHelper.writeOutputImage(bin + "/" + fileRoot + "_ws.png", imgW);
        }
    }


}
