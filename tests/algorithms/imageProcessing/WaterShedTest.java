package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class WaterShedTest extends TestCase {

    public WaterShedTest() {
    }

    public void testLower() throws Exception {

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

        WaterShed ws = new WaterShed();

        int[][] lowerComplete = ws.lower(img09, points);

        assertNotNull(lowerComplete);

        GreyscaleImage imgL = new GreyscaleImage(lowerComplete.length, lowerComplete[0].length);
        for (int i = 0; i < lowerComplete.length; ++i) {
            for (int j = 0; j < lowerComplete[0].length; ++j) {
                imgL.setValue(i, j, lowerComplete[i][j]);
            }
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

        StringBuilder sb2 = new StringBuilder();
        for (int j = 0; j < lowerComplete[0].length; ++j) {
            sb2.append(String.format("row %4d:", j));
            for (int i = 0; i < lowerComplete.length; ++i) {
                sb2.append(String.format("%4d", lowerComplete[i][j]));
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

        int[][] distToLowerIntensityPixel = ws.getDistToLowerIntensityPixel();

        for (PairInt p : points) {

            int x = p.getX();
            int y = p.getY();

            Set<PairInt> visited = new HashSet<PairInt>();

            int dist = findShortestPathToLowerIntensity(img09, visited, p, points, regionalMinima);

            assertTrue(distToLowerIntensityPixel[x][y] == dist);
        }

        // ------- assert final result, lowerComplete:
        //         that each pixel, if not in a regional minima, has an adjacent
        //         lower intensity pixel
        for (PairInt p : points) {
            if (regionalMinima.contains(p)) {
                continue;
            }
            int v = lowerComplete[p.getX()][p.getY()];
            boolean foundLowerIntNghbr = false;
            for (int vIdx = 0; vIdx < dxs8.length; ++vIdx) {
                int x2 = p.getX() + dxs8[vIdx];
                int y2 = p.getY() + dys8[vIdx];
                PairInt p2 = new PairInt(x2, y2);
                if (points.contains(p2)) {
                    int v2 = lowerComplete[x2][y2];
                    if (v2 < v) {
                        foundLowerIntNghbr = true;
                        break;
                    }
                }
            }
            assertTrue(foundLowerIntNghbr);
        }
    }
    
    public void testUnionFindComponentLabelling() throws Exception {

        int w = 5;
        int h = 5;
        
        /*
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

        WaterShed ws = new WaterShed();

        int[][] labeled = ws.unionFindComponentLabelling(lowerComplete);

        StringBuilder sb2 = new StringBuilder("input:\n");
        for (int j = 0; j < lowerComplete[0].length; ++j) {
            sb2.append(String.format("row %4d:", j));
            for (int i = 0; i < lowerComplete.length; ++i) {
                sb2.append(String.format("%4d", lowerComplete[i][j]));
            }
            sb2.append("\n");
        }
        System.out.println(sb2.toString());
        
        sb2 = new StringBuilder("labeled components:\n");
        for (int j = 0; j < labeled[0].length; ++j) {
            sb2.append(String.format("row %4d:", j));
            for (int i = 0; i < labeled.length; ++i) {
                sb2.append(String.format("%4d", labeled[i][j]));
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
        im[4] = new int[]{0, 1, 2, 3, 4};
        
        GreyscaleImage img = new GreyscaleImage(im.length, im[0].length);
        Set<PairInt> points = new HashSet<PairInt>();
        for (int i = 0; i < im.length; ++i) {
            for (int j = 0; j < im[0].length; ++j) {
                points.add(new PairInt(i, j));
                img.setValue(i, j, im[i][j]);
            }
        }

        WaterShed ws = new WaterShed();

        int[][] lowered = ws.lower(img, points);
        
        StringBuilder sb2 = new StringBuilder("lowered:\n");
        for (int j = 0; j < lowered[0].length; ++j) {
            sb2.append(String.format("row %4d:", j));
            for (int i = 0; i < lowered.length; ++i) {
                sb2.append(String.format("%4d", lowered[i][j]));
            }
            sb2.append("\n");
        }
        System.out.println(sb2.toString());
        
        int[][] labelled = ws.unionFindComponentLabelling(lowered);
        sb2 = new StringBuilder("labelled 0:\n");
        for (int j = 0; j < labelled[0].length; ++j) {
            sb2.append(String.format("row %4d:", j));
            for (int i = 0; i < labelled.length; ++i) {
                sb2.append(String.format("%4d", labelled[i][j]));
            }
            sb2.append("\n");
        }
        System.out.println(sb2.toString());
        
        
        CustomWatershedDAG dag = ws.createLowerIntensityDAG(lowered);
        
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
        
        int[][] labelled2 = ws.unionFindWatershed(lowered);
        
        assertNotNull(labelled2);
        
        sb2 = new StringBuilder("labelled final:\n");
        for (int j = 0; j < labelled2[0].length; ++j) {
            sb2.append(String.format("row %4d:", j));
            for (int i = 0; i < labelled2.length; ++i) {
                sb2.append(String.format("%4d", labelled2[i][j]));
            }
            sb2.append("\n");
        }
        System.out.println(sb2.toString());
    }
}
