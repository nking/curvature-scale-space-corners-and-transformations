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
import java.util.Set;
import java.util.Stack;

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
     * a set of points found as the regional minima.
     * This is populated by the method named lower.
     */
    private Set<PairInt> regionalMinima = null;

    /**
     * a map with key = PairInt(x, y) and value = disjoint set of the level
     * components (contiguous pixels of same intensity) whose set parent is
     * the representative for that level (which may be a regional minima).
     * This is populated by the method named unionFindComponentLabelling.
     */
    private Map<PairInt, DisjointSet2Node<PairInt>> componentLabelMap = null;

    private final static PairInt sentinel = new PairInt(-1, -1);

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
     * from
      "The Watershed Transform: Definitions, Algorithms and Parallelization Strategies"
      Roerdink and Meijster, 2001, Fundamenta Informaticae 41 (2001) 187–228

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

 System.out.println(printParents(parentMap));

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
     * Note this method uses the by-products of the methods named lower and
     * unionFindComponentLabelling.  To use this method as a standalone
     * invocation requires some changes to accept arguments for the by-products
     * or to re-solve for similar data in this method.
     *
     * @param im a lower complete image
     * @return
     */
    protected int[][] unionFindWatershed(int[][] im) {

        if ((distToLowerIntensityPixel == null) || (regionalMinima == null) ||
            (componentLabelMap == null)) {
            throw new IllegalStateException("algorithm currently depends upon "
            + "previous use of the methods named lower and unionFindComponentLabelling");
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

System.out.println(pPoint.toString() + ":");
                PairInt repr = resolveIterative(pPoint, dag);
System.out.println("   ==> " + pPoint.toString() + " repr=" + repr.toString() + "\n");

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
    */

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

        /*
        TODO: edit to use a points set or make a method which will only operate
        on the image for locations in point set
        */

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

    /**
     *
     * @param p
     * @param dag
     * @param recursionLevel level of recursion used only for logging
     * @return repr value indicating whether the point is a watershed (indicated
     * by returning the sentinel) or should be assigned the component level of
     * the returned point.
     */
    private PairInt resolve(PairInt p, CustomWatershedDAG dag) {

        if ((regionalMinima == null) || (componentLabelMap == null)) {
            throw new IllegalStateException("algorithm currently depends upon "
            + "previous use of the methods named lower and unionFindComponentLabelling");
        }

        if (dag.isResolved(p)) {
            return dag.getResolved(p);
        }

        DisjointSet2Node<PairInt> set = componentLabelMap.get(p);
        PairInt repr = set.getParent().getMember();

        int n = dag.getConnectedNumber(p);

        if (n == 0) {
            dag.setToResolved(p, repr);
            return repr;
        }

        int i = 0;

        while ((i < n) && !repr.equals(sentinel)) {

            PairInt lowerNode = dag.getConnectedNode(p, i);

            if (!lowerNode.equals(p) && !lowerNode.equals(sentinel)) {

                PairInt prevLowerNode = lowerNode;

                lowerNode = resolve(prevLowerNode, dag);

                dag.resetConnectedNode(p, i, lowerNode);
            }

            if (i == 0) {

                repr = lowerNode;

            } else if (!lowerNode.equals(repr)) {

                repr = sentinel;
            }

            ++i;
        }

        dag.setToResolved(p, repr);

        return repr;
    }

    private String createKey(int level, PairInt p, int i) {
        //TODO: replace this with encoding to a long using characteristics
        // such as image dimensions (coords converted to image index),
        // max level expected (nPixels), and knowledge that i will always be <= 8
        String key = String.format("%d_%d_%d_%d", level, p.getX(), p.getY(), i);
        return key;
    }

    /**
     *
     * @param p
     * @param dag
     * @return repr value indicating whether the point is a watershed (indicated
     * by returning the sentinel) or should be assigned the component level of
     * the returned point.
     */
    private PairInt resolveIterative(PairInt p, CustomWatershedDAG dag) {

        if ((regionalMinima == null) || (componentLabelMap == null)) {
            throw new IllegalStateException("algorithm currently depends upon "
            + "previous use of the methods named lower and unionFindComponentLabelling");
        }

        /*
        Recursion in java cannot be refactored to use tail recursion to release
        method frame memory, so to improve the use of memory, a recursive
        method has to be made iterative.

        To make iterative, have to replace the method frame loading and
        unloading with parallel stacks of arguments given to a loop which pops
        each stack to get current arguments, computes the result and stores that
        in a results map accessible to subsequent iterations.

        The composite key for the results is referred to as prevCompKey here
        and for simplicity while testing is a string, but should be changed to
        use encoding to be more efficient.
        TODO: improve the results map key

            level  p            i in level    prevCompKey
            0      p            0
            1      p.c[0]       0             "0 p      0"  <---place result here for i=0
            2      p.c[0].c[0]  0             "1 p.c[0] 0"  <---place result here
            1      p.c[1]       1             "0 p      0"  <---place result here for i=1

           paused edit here... design handling of results.

            level  p            i in level    prevCompKey
            0      p            0                           resolve cnctn: push onto stack "1 p.c[0] 0" pCKey="0 p 0" reprLevel="repr"
                                                            pop "1 p.c[0] 0" pCKey="0 p 0" repr0
            1      p.c[0]       0             "0 p      0"

            1      p.c[0]       0             "0 p      0"  <---place result here for i=0
            2      p.c[0].c[0]  0             "1 p.c[0] 0"  <---place result here
            1      p.c[1]       1             "0 p      0"



            level  p            i in level    prevCompKey
            0      p            0
               resolve connection p.c[i=0]
                  resolve connection p.c[i=0].c[i=0]
                      pop repr.  if i0>0, retrieve repr from map (and remove it)<--
                  resolve connection p.c[i=0].c[i=1]
                  process
            */

        // TODO: simplify the branch logic

        Stack<PairInt> currentP = new Stack<PairInt>();
        currentP.add(p);
        Stack<Integer> currentI = new Stack<Integer>();
        currentI.add(Integer.valueOf(0));
        Stack<Integer> currentLevel = new Stack<Integer>();
        currentLevel.add(Integer.valueOf(0));
        Stack<String> currentPrevCompKey = new Stack<String>();
        currentPrevCompKey.add(createKey(0, p, 0));
        Stack<Integer> currentLevelEnd = new Stack<Integer>();
        currentLevelEnd.add(Integer.valueOf(-1));

        Map<String, PairInt> resultsMap = new HashMap<String, PairInt>();

        int nIter = 0;

        skipped:

        while (!currentP.isEmpty()) {

            PairInt p0 = currentP.pop();
            int i0 = currentI.pop().intValue();
            int level0 = currentLevel.pop().intValue();
            String prevCompKey0 = currentPrevCompKey.pop();
            String currentCompKey = createKey(level0, p0, i0);
            boolean currentLevelEnd0 = (currentLevelEnd.pop().intValue() > 0);

            PairInt repr = null;

            boolean foundRecursionEnd = false;
            if ((nIter > 0) && p0.equals(p)) {
                foundRecursionEnd = true;
            }
            nIter++;

            if (!dag.isResolved(p0) && !foundRecursionEnd) {

                if (resultsMap.containsKey(currentCompKey)) {
                    repr = resultsMap.remove(currentCompKey);
                } else {
                    repr = componentLabelMap.get(p0).getParent().getMember();
                }

                int n = dag.getConnectedNumber(p0);

                for (int i = i0; i < n; ++i) {

                    if (repr.equals(sentinel)) {

System.out.println(Integer.toString(level0) + " " + p0.toString() + " repr=" + repr + " i=" + Integer.valueOf(i));

                        return repr;
                    }

                    PairInt lowerNode = dag.getConnectedNode(p0, i);

                    boolean isRes = dag.isResolved(lowerNode)
                        && (dag.getConnectedNumber(p0) == (i + 1));

                    if (isRes) {
                        repr = lowerNode;
                    }

                    if (!currentLevelEnd0 && !isRes && !lowerNode.equals(p0) && !lowerNode.equals(sentinel)) {

                        /*
                        lowerNode = resolveIterative(lowerNode, dag, recursionLevel + 1);
                        dag.resetConnectedNode(p0, i, lowerNode);
                        */

                        //TODO: change the loop from i0 to n above
                        if ((i == 0) && (n > 0)) {

                            //TODO: may need to revisit and refactor recursive to store results
                            // in one loop and process results in next loop
                            // making it easier to validate there, then separate here

                            int nAdded = 1;
                            currentP.push(p0);
                            currentI.push(i0);
                            currentLevel.push(Integer.valueOf(level0));
                            currentPrevCompKey.push(prevCompKey0);
                            currentLevelEnd.push(Integer.valueOf(1));

                            // add all i's onto the stack in reverse order
                            for (int ii = (n - 1); ii >= i0; --ii) {
                                currentP.push(dag.getConnectedNode(p0, ii));
                                currentI.push(0);
                                currentLevel.push(Integer.valueOf(level0 + 1));
                                currentPrevCompKey.push(createKey(level0, p0, ii));
                                currentLevelEnd.push(Integer.valueOf(-1));
                                nAdded++;
                            }
                            continue skipped;
                        }
                    }

                    if (i == 0) {
                        repr = lowerNode;
                    } else if (!lowerNode.equals(repr)) {
                        repr = sentinel;
                    }
                }

                dag.setToResolved(p0, repr);

            } else {

                repr = dag.getResolved(p0);
            }

System.out.println(Integer.toString(level0) + " " + p0.toString()
+ " repr=" + repr + " prevCompKey0=" + prevCompKey0);

            if (repr != null) {

                resultsMap.put(prevCompKey0, repr);

                /*
                logic for the root of the composite key to see if 2 separate paths
                have different results, hence a watershed.
                */

                String[] items = prevCompKey0.split("_");
                String rootKey = items[0] + "_" + items[1] + "_" + items[2];
                PairInt existingRepr = resultsMap.get(rootKey);
                if ((existingRepr != null) && (repr != null)) {
                    if (!existingRepr.equals(sentinel) && !existingRepr.equals(repr)) {
                        repr = sentinel;
                        dag.setToResolved(p, repr);
                        return repr;
                    }
                } else if (repr != null) {
                    resultsMap.put(rootKey, repr);
                }
            }
        }

        PairInt repr = resultsMap.get(createKey(0, p, 0));

        if (repr != null) {
            if (!dag.isResolved(p)) {
                dag.setToResolved(p, repr);
            }
        }

        return repr;
    }

}

