package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.disjointSets.DisjointSet2Helper;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.graphs.CustomWatershedDAG;
import algorithms.graphs.CustomWatershedNode;
import algorithms.util.PairInt;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

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
 Instead, the current implementation for disjoint sets follows "Introduction
 to Algorithms" by Cormen et al. which include improvements suggested by
 Tarjan too.

 Notes on parallelization are in Section 5 of Roerdink and Meijster 2001.
 </pre>

 * The image is first transformed into a lower complete image and then
 * the watershed is computed.
 *
 * @author nichole
 */
public class WaterShed extends AbstractWaterShed {

    /**
     * two dimensional matrix of the shortest distance of a pixel to
     * a lower intensity pixel with respect to the original image reference
     * frame.  For example, if a pixel is surrounded by pixels with the same
     * intensity, the shortest distance for it will be larger than '1' because
     * no neighbors have a smaller intensity.
     * This is populated by the method named lower.
     */
    private int[][] distToLowerIntensityPixel = null;

    /**
     * create a component labeled image with watershed pixels labeled as '0'.
     * runtime is quasi-linear.
     * @param img
     * @return the labeled image.  Note that if all intensities in img are 
     * the same, the method will return null.
     */
    public int[][] createLabelledImage(GreyscaleImage img) {

        int[][] lowerComplete = lower(img);
        
        if (lowerComplete == null) {
            return null;
        }

        int[][] labelled2 = unionFindWatershed(lowerComplete);

        return labelled2;
    }

    /**
     * This method alters the image, specifically the plateaus, so that a best
     * path to lowest intensity is possible and less ambiguous. A plateau is a
     * region of where the pixels have the same intensities.
     * After this has finished, there should be no pixel which does not
     * have a neighbor of lower intensity if the pixel is not a regional
     * minimum.
     * runtime complexity is O(N_pixels).
     *
     * @param img
     * @return the lowered image.  Note that if all intensities are the same,
     * the method will return null.
     */
    protected int[][] lower(GreyscaleImage img) {

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

        regionalMinima = new HashSet<PairInt>();

        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;

        int dist;

        ArrayDeque<Integer> queue = new ArrayDeque<Integer>(img.getNPixels());

        // ---- init queue with points which have lower intensity neighbors ---
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {

                int v = img.getValue(x, y);
                int idx = img.getIndex(x, y);

                for (int vIdx = 0; vIdx < dxs8.length; ++vIdx) {

                    int x2 = x + dxs8[vIdx];
                    int y2 = y + dys8[vIdx];
                    if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }

                    int v2 = img.getValue(x2, y2);
                    if (v2 < v) {
                        queue.add(Integer.valueOf(idx));
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
        queue.add(sentinelInt);

        while (!queue.isEmpty()) {

            Integer index = queue.poll();

            if (index.equals(sentinelInt)) {

                if (!queue.isEmpty()) {

                    queue.add(sentinelInt);

                    //any point originally lacking lower intensity neighbors,
                    //now gets a larger distance

                    dist++;
                }

                continue;
            }

            int x = img.getCol(index.intValue());
            int y = img.getRow(index.intValue());

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

                    int idx2 = img.getIndex(x2, y2);

                    queue.add(Integer.valueOf(idx2));

                    lc[x2][y2] = -1;
                }
            }
        }

        for (int x = 0; x < w; ++x) {

            for (int y = 0; y < h; ++y) {

                distToLowerIntensityPixel[x][y] = lc[x][y];

                if (lc[x][y] != 0) {

                    lc[x][y] = dist * img.getValue(x, y) + lc[x][y] - 1;

                } else {

                    regionalMinima.add(new PairInt(x, y));

                    //as suggested by later paper, adapted for watershed by Union-Find
                    lc[x][y] = dist * img.getValue(x, y);
                }
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

    /**
     * Algorithm 4.7 Scan-line algorithm for labeling level components based on
     * disjoint sets.
     * from
      "The Watershed Transform: Definitions, Algorithms and Parallelization Strategies"
      Roerdink and Meijster, 2001, Fundamenta Informaticae 41 (2001) 187–228

     * @param im greyscale image (does not need to be lower complete)
     * @return
     */
    protected int[][] unionFindComponentLabelling(int[][] im) {

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

        Map<PairInt, DisjointSet2Node<PairInt>> parentMap = new
            HashMap<PairInt, DisjointSet2Node<PairInt>>();

        // init map
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                PairInt pPoint = new PairInt(i, j);
                DisjointSet2Node<PairInt> pNode =
                    disjointSetHelper.makeSet(new DisjointSet2Node<PairInt>(pPoint));
                parentMap.put(pPoint, pNode);
            }
        }

        //Firstpass
        PairInt reprPoint;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {

                PairInt pPoint = new PairInt(i, j);

                reprPoint = pPoint;

                int x = pPoint.getX();
                int y = pPoint.getY();
                int vP = im[x][y];

                List<PairInt> qPoints = new ArrayList<PairInt>();

                //for all q ∈ Neighbor(p) with q ≺ p
                for (int vIdx = 0; vIdx < dLOX.length; ++vIdx) {
                    int x2 = x + dLOX[vIdx];
                    int y2 = y + dLOY[vIdx];
                    if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }

                    PairInt qPoint = new PairInt(x2, y2);

                    int vQ = im[x2][y2];

                    if (vP == vQ) {

                        // find r, the representative of the neighbors with
                        // same image intensity, as the lexicographically
                        // smallest location

                        //r ← r min FindRoot(q);

                        DisjointSet2Node<PairInt> qParent = disjointSetHelper.findSet(
                            parentMap.get(qPoint));

                        if (qParent.getMember().getX() < reprPoint.getX()) {
                            reprPoint = qPoint;
                        } else if ((qParent.getMember().getX() == reprPoint.getX())
                            && (qParent.getMember().getY() < reprPoint.getY())) {
                            reprPoint = qPoint;
                        }
                        qPoints.add(qPoint);
                    }
                }

                //parent[p] ← r
                if (!qPoints.isEmpty()) {

                    DisjointSet2Node<PairInt> parent = disjointSetHelper.union(
                        parentMap.get(reprPoint), parentMap.get(pPoint));

                    for (PairInt qPoint : qPoints) {
                        if (qPoint.equals(reprPoint)) {
                            continue;
                        }
                        //PathCompress(q, r)

                        DisjointSet2Node<PairInt> qParent = disjointSetHelper.union(
                            parentMap.get(reprPoint), parentMap.get(qPoint));
                    }
                }
            } // end j loop
        } // end i loop

 //System.out.println(printParents(parentMap));

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

                DisjointSet2Node<PairInt> parent = disjointSetHelper.findSet(
                    parentMap.get(pPoint));

                if (parent.getMember().equals(pPoint)) {
                    // root pixel
                    label[i][j] = curLabel;
                    curLabel++;
                } else {
                    //Resolve unresolved equivalences
                    // parent[p] = parent[parent[p]]
                    parentMap.put(pPoint, parent);
                    label[i][j] = label[parent.getMember().getX()][parent.getMember().getY()];
                }
            }
        }

        componentLabelMap = parentMap;

        return label;
    }

    /**
     * Algorithm 4.8 Watershed transform w.r.t. topographical distance based on
     * disjoint sets.
     *
     * Note this method uses the by-products of the method named lower.
     * To use this method as a standalone invocation requires some changes to
     * accept arguments for the by-products or to re-solve for similar data in
     * this method.
     *
     * @param im a lower complete image
     * @return
     */
    protected int[][] unionFindWatershed(int[][] im) {

        if (im == null) {
            throw new IllegalStateException("im cannot be null");
        }
        if ((distToLowerIntensityPixel == null) || (regionalMinima == null)) {
            throw new IllegalStateException("algorithm currently depends upon "
            + "previous use of the methods named lower");
        }

        // uses regionalMinima
        CustomWatershedDAG dag = createLowerIntensityDAG(im);

        int w = im.length;
        int h = im[0].length;

        final int wshed = 0;

        //initialize image lab with distinct labels for minima
        //LabelInit
        int[][] labeled = unionFindComponentLabelling(im);

        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                PairInt pPoint = new PairInt(i, j);

                PairInt repr = resolve(pPoint, dag);

                int value;
                if (repr.equals(sentinel)) {
                    value = wshed;
                } else {
                    value = labeled[repr.getX()][repr.getY()];
                }
                labeled[pPoint.getX()][pPoint.getY()] = value;
            }
        }

        return labeled;
    }

    private String printParents(Map<PairInt, DisjointSet2Node<PairInt>> parentMap) {

        DisjointSet2Helper dsHelper = new DisjointSet2Helper();

        Map<PairInt, List<PairInt>> parentValueMap = new HashMap<PairInt, List<PairInt>>();

        for (Entry<PairInt, DisjointSet2Node<PairInt>> entry : parentMap.entrySet()) {

            PairInt child = entry.getKey();
            PairInt parent = dsHelper.findSet(entry.getValue()).getMember();

            List<PairInt> children = parentValueMap.get(parent);
            if (children == null) {
                children = new ArrayList<PairInt>();
                parentValueMap.put(parent, children);
            }
            children.add(child);
        }

        StringBuilder sb = new StringBuilder();
        for (Entry<PairInt, List<PairInt>> entry : parentValueMap.entrySet()) {
            PairInt parent = entry.getKey();
            List<PairInt> children = entry.getValue();
            sb.append("parent: ").append(parent.toString());
            sb.append("    children: ");
            for (PairInt c : children) {
                sb.append(" ").append(c.toString());
            }
            sb.append("\n");
        }
        return sb.toString();
    }

    protected CustomWatershedDAG createLowerIntensityDAG(int[][] lowerCompleteIm) {

        if (regionalMinima == null) {
            throw new IllegalStateException(
                "method needs lower to have been invoked before using this");
        }
        if (lowerCompleteIm == null) {
            throw new IllegalStateException("lowerCompleteIm cannot be null");
        }

        int w = lowerCompleteIm.length;
        int h = lowerCompleteIm[0].length;

        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;

        CustomWatershedDAG dag = new CustomWatershedDAG(w * h);

        int[] diffInt = new int[8];
        PairInt[] neighbors = new PairInt[8];

        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {

                PairInt p = new PairInt(i, j);

                if (regionalMinima.contains(p)) {

                    // componentLabelMap has the representative for this node
                    dag.insert(p, new CustomWatershedNode(p, 0));

                } else {

                    int x = p.getX();
                    int y = p.getY();
                    int v = lowerCompleteIm[x][y];

                    int nc = 0;

                    for (int nIdx = 0; nIdx < dxs8.length; ++nIdx) {
                        int x2 = x + dxs8[nIdx];
                        int y2 = y + dys8[nIdx];
                        if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                           continue;
                        }

                        int v2 = lowerCompleteIm[x2][y2];

                        if (v2 < v) {
                            diffInt[nc] = v - v2;
                            neighbors[nc] = new PairInt(x2, y2);
                            nc++;
                        }
                    }

                    dag.orderAndInsert(p, diffInt, neighbors, nc);

                    assert(nc != 0);
                }
            }
        }

        return dag;
    }

}

