package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.util.DisjointSet2Helper;
import algorithms.util.DisjointSet2Node;
import algorithms.util.PairInt;
import java.util.ArrayDeque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * A watershed algorithm for use in image segmentation that is based upon
 * the algorithms described in
  <pre>
  Roerdink and Meijster 2001
  "The Watershed Transform: Definitions, Algorithms and Parallelization Strategies",
  Fundamenta Informaticae 41 (2001) 187–228, Section 4.2.4
  and
  Meijster and Roerdink (1998?),
  "A Disjoint Set Algorihm for the Watershed Transform"
  http://www.researchgate.net/publication/2407714_A_Disjoint_Set_Algorithm_For_The_Watershed_Transform

 Note the above authors credit the 2 Disjoint Set methods,
 especially the disjoint set path compression,
 used in the watershed union find to
 Tarjan, R. E. Data Structures and Network Algorithms. SIAM, 1983.
 Those are not yet implemented strictly as suggested here.
 Instead, the current implementation follows "Introduction
 to Algorithms" by Cormen et al. which include improvements suggested by
 Tarjan too.
 </pre>

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
     * have a neighbor of lower intensity if the pixel is not a regional
     * minimum.
     * runtime complexity is O(N).
     *
     * @param img
     * @return
     */
    protected int[][] lower(GreyscaleImage img, Set<PairInt> points) {

        /*
        TODO: create a helper method to determine the bounds of points
        and then pass the offsets for that into this method to allow
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
     * get the two dimensional matrix of the shortest distance of a pixel to
     * a lower intensity pixel with respect to the original image reference
     * frame.  For example, if a pixel is surrounded by pixels with the same
     * intensity, the shortest distance for it will be larger than '1' because
     * no neighbors have a smaller intensity.
     * @return the distToLowerIntensityPixel
     */
    public int[][] getDistToLowerIntensityPixel() {
        return distToLowerIntensityPixel;
    }

    public Set<PairInt> getRegionalMinima() {
        return regionalMinima;
    }

    /**
     * Algorithm 4.7 Scan-line algorithm for labelling level components based on
     * disjoint sets.
     *
     * @param im greyscale image (does not need to be lower complete)
     * @return
     */
    protected int[][] unionFindComponentLabelling(int[][] im) {

        /*
        TODO: when make changes above to use a reduced portion of the image
        via the points set, can consider visiting a smaller number of points
        here too.  A LinkedHashSet can be created with lexicographical ordering
        rules.  The LinkedHashSet can be created with one pass through 0 to
        width and 0 to height or the points set can be sorted and entered into
        LinkedHashSet in lexicographical order depending upon comparison of
        n_points in points set and n_pixels = width*height,
        O(N_points*lg2(N_points)) vs O(N_pixels), respectively.
        */

        int w = im.length;
        int h = im[0].length;

        /*
        search for neighbors q of p that have smaller lexicographical values
        q ≺ p : (i_q < i_p) ∨ ((i_q == i_p) ∧(j_q < j_p))

          (-1, 1)
          (-1, 0)   p=(0,  0)
          (-1,-1)     (0, -1)
        */
        int[] dLOX = new int[]{-1, -1, -1,  0};
        int[] dLOY = new int[]{ 1,  0, -1, -1};

        DisjointSet2Helper disjointSetHelper = new DisjointSet2Helper();

        DisjointSet2Node[] parents = new DisjointSet2Node[w * h];

        //Firstpass
        PairInt rPoint;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {

                PairInt pPoint = new PairInt(i, j);

                rPoint = pPoint;

                int x = pPoint.getX();
                int y = pPoint.getY();
                int v = im[x][y];

                //for all q ∈ Neighbor(p) with q ≺ p
                for (int vIdx = 0; vIdx < dLOX.length; ++vIdx) {
                    int x2 = x + dLOX[vIdx];
                    int y2 = y + dLOY[vIdx];
                    if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }

                    PairInt qPoint = new PairInt(x2, y2);

                    int v2 = im[x2][y2];
      
                    if (v == v2) {
                        //r ← r min FindRoot(q);
                        //min denotes minimum w.r.t. lexicographical order
                    }
                }

                //parent[p] ← r

                //for all q ∈ Neighbor(p) with q ≺ p
                for (int vIdx = 0; vIdx < dLOX.length; ++vIdx) {
                    int x2 = x + dLOX[vIdx];
                    int y2 = y + dLOY[vIdx];
                    if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }

                    PairInt qPoint = new PairInt(x2, y2);

                    int v2 = im[x2][y2];

                    if (v == v2) {
                       //PathCompress(q, r)
                    }
                }
            } // end j loop
        } // end i loop

        /*
        In a second pass through the input image, the output image lab is
        created. All root pixels get a distinct label; for any other pixel p
        its path is compressed, making explicit use of the order imposed on
        parent (see line 29 in Algorithm 4.7), and p gets the label of its
        representative.
        */

        int[][] label = new int[w][];
        for (int i = 0; i < w; ++i) {
            label[i] = new int[h];
        }

        //Secondpass
        int curLabel = 1;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {

                PairInt pPoint = new PairInt(i, j);
                
            }
        }

        throw new UnsupportedOperationException("not yet implemented");
        
        //return label;
    }

    /*
    Scan-line algorithm for labelling level components based on disjoint sets.

    procedure union-find-ComponentLabelling
        Input: grey scale image im on digital grid G = (D, E).
        Output: image lab on D, with labelled level components. (∗Uses array parent of pointers. ∗)

    //Firstpass
    for all p ∈ D in lexicographical order do
        r←p
        for all q ∈ Neighbor(p) with q ≺ p do
            if im[q] = im[p] then
                r ← r min FindRoot(q) //min denotes minimum w.r.t. lexicographical order
            end if
        end for
        parent[p] ← r
        for all q ∈ Neighbor(p) with q ≺ p do //compress paths
            if im[q] = im[p] then
                PathCompress(q, r)
            end if
        end for
    end for

    //Secondpass∗)
    curlab←1 //curlab is the current label∗)
    for all p ∈ D in lexicographical order do
        if parent[p] = p then // p is a root pixel ∗)
            lab[p] = curlab
            curlab = curlab + 1
        else
            parent[p] = parent[parent[p]] //Resolve unresolved equivalences ∗)
            lab[p] = lab[parent[p]]
        end if
    end for

    // The authors credit Tarjan for the DisjointSet improvements, espec the
    // path compression.
    // Tarjan, R. E. Data Structures and Network Algorithms. SIAM, 1983
    functionFindRoot(p:pixel)
        while parent[p] ̸= p do
            r←parent[p] ;
            p←r
        end while
        return r

    procedurePathCompress(p:pixel,r:pixel)
        while parent[p] ̸= r do
            h←parent[p] ;
            parent[p]←r ;
            p←h
        end while
    */

    /**
     * Algorithm 4.8 Watershed transform w.r.t. topographical distance based on disjoint sets.
     * @param im a lower complete image
     * @return
     */
    protected int[][] unionFindWatershed(int[][] im) {

        /*
        input is a lower complete image.

        internally uses a DAG and disjoint sets

        algorithm 4.8 of reference above
        */

        throw new UnsupportedOperationException("not yet implemented");
    }

    /*
    procedure union-find-Watershed
        Input: lower complete graph G′ = (V, E′).
        Output: labelled image lab on V .

        //label of the watershed pixels
        #define wshed 0

        //fictitious coordinates of the watershed pixels
        #define W (-1,-1)

        //initialize image lab with distinct labels for minima
        LabelInit

        //give p the label of its representative
        for all p ∈ V do
            rep ← Resolve (p)
            if rep ̸= W then
                lab[p] ← lab[rep]
            else
                lab[p] ← wshed end if
            end if
        end for

    //Recursive function for resolving the downstream paths of the lower complete graph
    //Returns representative element of pixel p, or W if p is a watershed pixel
    function Resolve (p : pixel)

        i←1;

        //some value such that rep ̸= W
        rep←(0,0)

        //CON indicates the connectivity
        while (i ≤ CON) and (rep ̸= W) do
            if (sln[p, i] ̸= p) and (sln[p, i] ̸= W) then
                sln[p, i] ← Resolve (sln[p, i])
            end if
            if i=1 then
                rep ← sln[p, 1]
            else if sln[p, i] ̸= rep then
                rep ← W
                for j←1 to CON do
                    sln[p, j] ← W
                end for
            end if
            i←i+1
        end while

        return rep
    */
}

