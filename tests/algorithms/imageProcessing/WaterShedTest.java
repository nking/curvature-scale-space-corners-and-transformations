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

        StringBuilder sb2 = new StringBuilder();
        for (int j = 0; j < labeled[0].length; ++j) {
            sb2.append(String.format("row %4d:", j));
            for (int i = 0; i < labeled.length; ++i) {
                sb2.append(String.format("%4d", labeled[i][j]));
            }
            sb2.append("\n");
        }
        System.out.println(sb2.toString());
        
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
}
