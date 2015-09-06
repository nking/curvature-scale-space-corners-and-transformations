package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.ArrayDeque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * A watershed algorithm for use in image segmentation that is based upon
 * the algorithms described in
 * <pre>
 * Roerdink and Meijster 2001
  "The Watershed Transform: Definitions, Algorithms and Parallelization Strategies",
  Fundamenta Informaticae 41 (2001) 187–228, Section 4.2.4
  and
  Meijster and Roerdink (1998?),
  "A Disjoint Set Algorihm for the Watershed Transform"
  http://www.researchgate.net/publication/2407714_A_Disjoint_Set_Algorithm_For_The_Watershed_Transform
 </pre>
 *
 * The image is first transformed into a lower complete image and then
 * the watershed is computed.
 *
 * @author nichole
 */
public class WaterShed {

    private int[][] distToLowerIntensityPixel = null;

    private Set<PairInt> regionalMinima = new HashSet<PairInt>();

    /**
     * This method alters the image, specifically the plateaus, so that a best
     * path to lowest intensity is possible and less ambiguous. A plateau is a
     * region of where the pixels have the same intensities.
     * After this has finished, there should be no pixel which does not
     * have a neighbor of lower intensity if the pixel is not the regional
     * minimum.
     * runtime complexity is O(N).
     *
     * @param img
     * @return
     */
    protected int[][] lower(GreyscaleImage img, Set<PairInt> points) {

        /*
        TODO: create a helper method to determine the bounds of points
        and then pass the offsets for that into this method to allows
        using a smaller returned two dimensional array whose coordinate
        reference frame is (x - xOffset, y - yOffset).
        The algorithm below would need to be adjusted for it where
        int[][] lc is used.
        */

        PairInt sentinel = new PairInt(-1, -1);

        int w = img.getWidth();
        int h = img.getHeight();

        int[][] lc = new int[w][];
        for (int i = 0; i < w; ++i) {
            lc[i] = new int[h];
        }

        distToLowerIntensityPixel = new int[w][];
        for (int i = 0; i < w; ++i) {
            distToLowerIntensityPixel[i] = new int[h];
        }

        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;

        int dist;

        ArrayDeque<PairInt> queue = new ArrayDeque<PairInt>(points.size());

        // ---- init queue with points which have lower intensity neighbors ---
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            int v = img.getValue(x, y);

            for (int vIdx = 0; vIdx < dxs8.length; ++vIdx) {
                int x2 = x + dxs8[vIdx];
                int y2 = y + dys8[vIdx];
                PairInt p2 = new PairInt(x2, y2);
                if (points.contains(p2)) {
                    int v2 = img.getValue(x2, y2);
                    if (v2 < v) {
                        queue.add(p);
                        lc[x][y] = -1;
                        break;
                    }
                }
            }
        }

        if (queue.isEmpty()) {
            // points contains only pixels of same intensity
            return null;
        }

        dist = 1;
        queue.add(sentinel);

        while (!queue.isEmpty()) {

            PairInt p = queue.poll();

            if (p.equals(sentinel)) {
                if (!queue.isEmpty()) {

                    queue.add(sentinel);

                    //any point originally lacking lower intensity neighbors,
                    //now gets a larger distance

                    dist++;
                }
                continue;
            }

            int x = p.getX();
            int y = p.getY();

            lc[x][y] = dist;

            for (int vIdx = 0; vIdx < dxs8.length; ++vIdx) {
                int x2 = x + dxs8[vIdx];
                int y2 = y + dys8[vIdx];

                if (x2 < 0 || y2 < 0 || (x2 > (img.getWidth() - 1)) ||
                    (y2 > (img.getHeight() - 1))) {
                    continue;
                }

                if ((img.getValue(x, y) == img.getValue(x2, y2)) &&
                    (lc[x2][y2] == 0)) {

                    queue.add(new PairInt(x2, y2));

                    lc[x2][y2] = -1;
                }
            }
        }

        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();

            distToLowerIntensityPixel[x][y] = lc[x][y];

            if (lc[x][y] != 0) {

                lc[x][y] = dist * img.getValue(x, y) + lc[x][y] - 1;

            } else {

                regionalMinima.add(p);

                //as suggested by later paper, adapted for watershed by Union-Find
                lc[x][y] = dist * img.getValue(x, y);
            }
        }

        return lc;
    }

    /**
     * get the two dimensional matrix of the distance of a pixel to the
     * regional minima.
     * @return the distToLowerIntensityPixel
     */
    public int[][] getDistToLowerIntensityPixel() {
        return distToLowerIntensityPixel;
    }

    public Set<PairInt> getRegionalMinima() {
        return regionalMinima;
    }
}

