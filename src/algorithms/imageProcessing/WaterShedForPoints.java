package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.disjointSets.DisjointSet2Helper;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.graphs.CustomWatershedDAG;
import algorithms.graphs.CustomWatershedNode;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
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
public class WaterShedForPoints {

    /**
     * two dimensional matrix of the shortest distance of a pixel to
     * a lower intensity pixel with respect to the original image reference
     * frame.  For example, if a pixel is surrounded by pixels with the same
     * intensity, the shortest distance for it will be larger than '1' because
     * no neighbors have a smaller intensity.
     * This is populated by the method named lower.
     */
    private Map<PairInt, Integer> distToLowerIntensityPixel = null;

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
     * create for the points in the image, a component labelled image with
     * watershed pixels labelled as '0'.
     * runtime is quasi-linear.
     * @param img
     * @param points
     * @return
     */
    public Map<PairInt, Integer> createLabelledImage(GreyscaleImage img,
        Set<PairInt> points) {

        Map<PairInt, Integer> lowerComplete = lower(img, points);

        Map<PairInt, Integer> labelled2 = unionFindWatershed(lowerComplete);

        return labelled2;
    }

    /**
     * get the two dimensional matrix of the shortest distance of a pixel to
     * a lower intensity pixel with respect to the original image reference
     * frame.  For example, if a pixel is surrounded by pixels with the same
     * intensity, the shortest distance for it will be larger than '1' because
     * no neighbors have a smaller intensity.
     * @return the distToLowerIntensityPixel
     */
    public Map<PairInt, Integer> getDistToLowerIntensityPixel() {
        return distToLowerIntensityPixel;
    }

    public Set<PairInt> getRegionalMinima() {
        return regionalMinima;
    }

     /**
     * This method alters the image, specifically the plateaus, so that a best
     * path to lowest intensity is possible and less ambiguous. A plateau is a
     * region of where the pixels have the same intensities.
     * After this has finished, there should be no pixel which does not
     * have a neighbor of lower intensity if the pixel is not a regional
     * minimum.
     * runtime complexity is O(N_points).
     *
     * @param img
     * @param points
     * @return
     */
    protected Map<PairInt, Integer> lower(GreyscaleImage img, Set<PairInt> points) {

        int w = img.getWidth();
        int h = img.getHeight();

        Map<PairInt, Integer> lc = new HashMap<PairInt, Integer>();

        distToLowerIntensityPixel = new HashMap<PairInt, Integer>();

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
                        lc.put(p, Integer.valueOf(-1));
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

            lc.put(p, Integer.valueOf(dist));

            for (int vIdx = 0; vIdx < dxs8.length; ++vIdx) {
                int x2 = x + dxs8[vIdx];
                int y2 = y + dys8[vIdx];

                PairInt p2 = new PairInt(x2, y2);

                if (!points.contains(p2)) {
                    continue;
                }

                Integer value2 = lc.get(p2);

                int v2 = (value2 == null) ? 0 : value2.intValue();

                if ((img.getValue(x, y) == img.getValue(x2, y2)) && (v2 == 0)) {

                    queue.add(new PairInt(x2, y2));

                    lc.put(p2, Integer.valueOf(-1));
                }
            }
        }

        for (PairInt p : points) {

            int x = p.getX();
            int y = p.getY();

            Integer value = lc.get(p);

            int v = (value == null) ? 0 : value.intValue();

            distToLowerIntensityPixel.put(p, Integer.valueOf(v));

            if (v != 0) {

                int v2 = dist * img.getValue(x, y) + v - 1;

                lc.put(p, Integer.valueOf(v2));

            } else {

                regionalMinima.add(p);

                //as suggested by later paper, adapted for watershed by Union-Find
                int v2 = dist * img.getValue(x, y);

                lc.put(p, Integer.valueOf(v2));
            }
        }

        return lc;
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
    protected Map<PairInt, Integer> unionFindWatershed(Map<PairInt, Integer> im) {

        if ((distToLowerIntensityPixel == null) || (regionalMinima == null)) {
            throw new IllegalStateException("algorithm currently depends upon "
            + "previous use of the methods named lower");
        }

        // uses regionalMinima
        CustomWatershedDAG dag = createLowerIntensityDAG(im);

        final Integer wshed = Integer.valueOf(0);

        //initialize image lab with distinct labels for minima
        //LabelInit
        Map<PairInt, Integer> labeled = unionFindComponentLabelling(im);

        for (Entry<PairInt, Integer> entry : labeled.entrySet()) {

            PairInt pPoint = entry.getKey();

            PairInt repr = resolveIterative(pPoint, dag);

            Integer value;
            if (repr.equals(sentinel)) {
                value = wshed;
            } else {
                value = labeled.get(repr);
            }
            labeled.put(pPoint, value);
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

    protected CustomWatershedDAG createLowerIntensityDAG(
        Map<PairInt, Integer> lowerCompleteIm) {

        if (regionalMinima == null) {
            throw new IllegalStateException(
                "method needs lower to have been invoked before using this");
        }
        if (lowerCompleteIm == null) {
            throw new IllegalStateException("lowerCompleteIm cannot be null");
        }

        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;

        CustomWatershedDAG dag = new CustomWatershedDAG(lowerCompleteIm.size());

        int[] diffInt = new int[8];
        PairInt[] neighbors = new PairInt[8];

        for (Entry<PairInt, Integer> entry : lowerCompleteIm.entrySet()) {

            PairInt p = entry.getKey();

            Integer value = entry.getValue();

            if (regionalMinima.contains(p)) {

                // componentLabelMap has the representative for this node
                dag.insert(p, new CustomWatershedNode(p, 0));

            } else {

                int x = p.getX();
                int y = p.getY();

                int nc = 0;

                for (int nIdx = 0; nIdx < dxs8.length; ++nIdx) {
                    int x2 = x + dxs8[nIdx];
                    int y2 = y + dys8[nIdx];
                    PairInt p2 = new PairInt(x2, y2);

                    if (!lowerCompleteIm.containsKey(p2)) {
                        continue;
                    }

                    Integer value2 = lowerCompleteIm.get(p2);

                    if (value2.intValue() < value.intValue()) {
                        diffInt[nc] = value.intValue() - value2.intValue();
                        neighbors[nc] = new PairInt(x2, y2);
                        nc++;
                    }
                }

                dag.orderAndInsert(p, diffInt, neighbors, nc);

                assert(nc != 0);
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

    /**
     *
     * @param p
     * @param dag
     * @param recursionLevel level of recursion used only for logging
     * @return repr value indicating whether the point is a watershed (indicated
     * by returning the sentinel) or should be assigned the component level of
     * the returned point.
     */
    private PairInt resolveIterative(final PairInt p, final CustomWatershedDAG dag) {

        if ((regionalMinima == null) || (componentLabelMap == null)) {
            throw new IllegalStateException("algorithm currently depends upon "
            + "previous use of the methods named lower and unionFindComponentLabelling");
        }

        /*
        Recursion in java cannot be refactored to use tail recursion (tail
        recursion avoids retaining the method frame), so to improve the use of
        memory when the potential depth of recursion is large, the recursive
        method is replaced here with an iterative.

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
        */

        Map<String, Set<PairInt>> resultsMap = new HashMap<String, Set<PairInt>>();

        final Stack<ResolveState> stack = new Stack<ResolveState>();

        addToStack(0, p, 0, null, -1, false, stack);
   
        while (!stack.isEmpty()) {

            ResolveState rs = stack.pop();

            PairInt repr;

            // instead of only processing at base or recursion, will check results
            // each time to possibly break earlier for sentinel
            //if (rs.isRecursionBase()) {
                //process the logic from the recursion as the returned result of
                //lowerNode = resolve(...)
                
                repr = findResult(rs, resultsMap);
                
                if (repr == sentinel) {
                    addToResults(createRootKey(0, p), repr, resultsMap);
                    dag.setToResolved(rs.getP(), repr);
                    break;
                } else if (rs.isRecursionBase()) {
                    addToResults(repr, rs, resultsMap);
                }
            //}

            if (dag.isResolved(rs.getP())) {
                repr = dag.getResolved(rs.getP());
                addToResults(repr, rs, resultsMap);
                if (rs.getPrevious() == null) {
                    break;
                }
            }

            if (repr == null) {
                repr = componentLabelMap.get(rs.getP()).getParent().getMember();
            }
            
            assert(repr != null);

            int n = dag.getConnectedNumber(rs.getP());

            if (n == 0) {                
                dag.setToResolved(rs.getP(), repr);
                if (rs.getPrevious() != null) {
                    dag.resetConnectedNode(rs.getPrevious().getP(), rs.getPreviousI(), repr);
                }
                addToResults(repr, rs, resultsMap);
                continue;
            }

            int i = rs.getI();          
              
            while ((i < n) && !repr.equals(sentinel)) {

                PairInt lowerNode = dag.getConnectedNode(rs.getP(), i);
                
                assert(lowerNode != null);

                if (!lowerNode.equals(rs.p) && !lowerNode.equals(sentinel)) {

                    // ---- add to stack to replace recursion ---
                    //PairInt prevLowerNode = lowerNode;
                    //lowerNode = resolve(prevLowerNode, dag);
                    //dag.resetConnectedNode(p, i, lowerNode);

                    // add current state to process the about to be added stack items:
                    // set rs0's i to i+1 so when popped, after processing result, starts at next
                    int prevI = (rs.getPrevious() == null) ? -1 : i;
                    ResolveState rs0 = addToStack(rs.getLevel(), rs.getP(), 
                        i + 1, rs.getPrevious(), prevI, true, stack);

                    addToStack(rs.getLevel() + 1, lowerNode, 0, rs0, i, false,
                        stack);
               
                    break;
                }

                if (i == 0) {
                    repr = lowerNode;
                    addToResults(repr, rs, resultsMap);
                } else if (!lowerNode.equals(repr)) {
                    repr = sentinel;
                    addToResults(createRootKey(0, p), repr, resultsMap);
                    break;
                }

                ++i;
            }
        }

        PairInt repr = findRootResult(p, resultsMap);

        if (repr != null && !dag.isResolved(p)) {
            dag.setToResolved(p, repr);
        }

        return repr;
    }
    
    /**
     * Algorithm 4.7 Scan-line algorithm for labelling level components based on
     * disjoint sets.
     * from
      "The Watershed Transform: Definitions, Algorithms and Parallelization Strategies"
      Roerdink and Meijster, 2001, Fundamenta Informaticae 41 (2001) 187–228

     The runtime is quasi-linear in the number of points (not the number of pixels).
     *
     * @param im greyscale image (does not need to be lower complete)
     * @return
     */
    protected Map<PairInt, Integer> unionFindComponentLabelling(
        Map<PairInt, Integer> im) {

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

        LinkedHashSet<PairInt> lOrderedPoints =
            MiscMath.lexicographicallyOrderPointsBySort(im.keySet());

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
        for (PairInt pPoint : lOrderedPoints) {
            DisjointSet2Node<PairInt> pNode = disjointSetHelper.makeSet(
                new DisjointSet2Node<PairInt>(pPoint));
            parentMap.put(pPoint, pNode);
        }

        //Firstpass
        PairInt reprPoint;
        for (PairInt pPoint : lOrderedPoints) {

            reprPoint = pPoint;

            int x = pPoint.getX();
            int y = pPoint.getY();

            Integer value = im.get(pPoint);

            assert(value != null);

            int vP = value.intValue();

            List<PairInt> qPoints = new ArrayList<PairInt>();

            //for all q ∈ Neighbor(p) with q ≺ p
            for (int vIdx = 0; vIdx < dLOX.length; ++vIdx) {
                int x2 = x + dLOX[vIdx];
                int y2 = y + dLOY[vIdx];

                PairInt qPoint = new PairInt(x2, y2);

                if (!lOrderedPoints.contains(qPoint)) {
                    continue;
                }

                Integer value2 = im.get(qPoint);

                assert(value2 != null);

                int vQ = value2.intValue();

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
        }

//System.out.println(printParents(parentMap));

        /*
        In a second pass through the input image, the output image lab is
        created. All root pixels get a distinct label; for any other pixel p
        its path is compressed, making explicit use of the order imposed on
        parent (see line 29 in Algorithm 4.7), and p gets the label of its
        representative.
        */

        Map<PairInt, Integer> label = new HashMap<PairInt, Integer>();

        //Secondpass
        int curLabel = 1;
        for (PairInt pPoint : lOrderedPoints) {

            DisjointSet2Node<PairInt> parent = disjointSetHelper.findSet(
                parentMap.get(pPoint));

            if (parent.getMember().equals(pPoint)) {

                // root pixel
                label.put(pPoint, Integer.valueOf(curLabel));

                curLabel++;

            } else {

                //Resolve unresolved equivalences
                // parent[p] = parent[parent[p]]
                parentMap.put(pPoint, parent);

                PairInt ePoint = new PairInt(parent.getMember().getX(),
                    parent.getMember().getY());

                Integer eLabel = label.get(ePoint);

                if (eLabel == null) {
                    eLabel = Integer.valueOf(0);
                }

                label.put(pPoint, eLabel);
            }
        }

        componentLabelMap = parentMap;

        return label;
    }

    private String createKey(int level, PairInt p, int i) {
        //TODO: replace this with encoding to a long using characteristics
        // such as image dimensions (coords converted to image index),
        // max level expected (nPixels), and knowledge that i will always be <= 8
        String key = String.format("%d_%d_%d_%d", level, p.getX(), p.getY(), i);
        return key;
    }

    private String createRootKey(int level, PairInt p) {
        //TODO: replace this with encoding to a long using characteristics
        // such as image dimensions (coords converted to image index),
        // max level expected (nPixels), and knowledge that i will always be <= 8
        String key = String.format("%d_%d_%d", level, p.getX(), p.getY());
        return key;
    }

    private ResolveState addToStack( int level, PairInt p, int i,
        ResolveState invokedFrom, int invokedFromI, boolean isRecursionBase,
        Stack<ResolveState> stack) {

        ResolveState rs = new ResolveState(level, p, i, invokedFrom, 
            invokedFromI, isRecursionBase);
        
        stack.add(rs);

        return rs;
    }

    private void addToResults(PairInt resolved, ResolveState rs,
        Map<String, Set<PairInt>> resultsMap) {

        if (rs.getPrevious() == null) {
            // special case.  this is the original p resolve is invoked for

            addToResults(createRootKey(rs.getLevel(), rs.getP()), resolved, resultsMap);

            return;
        }

        // add to root results
        addToResults(createRootKey(rs.getPrevious().getLevel(), rs.getPrevious().getP()), resolved,
            resultsMap);

        addToResults(createKey(rs.getPrevious().getLevel(), rs.getPrevious().getP(), rs.getPrevious().getI()),
            resolved, resultsMap);
    }

    /**
     * look for more than one result in the "root" key and if found, return
     * the sentinel (note that the invoker is responsible for sentinel state
     * logic), else get the specific key value from the resultsMap and remove
     * it while at it and return it.
     * @param rs
     * @param resultsMap
     * @return
     */
    private PairInt findResult(ResolveState rs, Map<String, Set<PairInt>> resultsMap) {

        return findResult(rs.getP(), rs.getLevel(), rs.getI(), resultsMap);
    }

    /**
     * look for more than one result in the "root" key and if found, return
     * the sentinel (note that the invoker is responsible for sentinel state
     * logic), else get the specific key value from the resultsMap and remove
     * it while at it and return it.
     * @param rs
     * @param resultsMap
     * @return
     */
    private PairInt findResult(PairInt p, int level, int i,
        Map<String, Set<PairInt>> resultsMap) {

        Set<PairInt> reprSet = resultsMap.remove(createKey(level, p, i));

        Set<PairInt> rootSet = resultsMap.get(createRootKey(level, p));

        if ((rootSet != null) && rootSet.size() > 1) {
            return sentinel;
        }
        if ((reprSet != null) && reprSet.size() > 0) {
            assert(reprSet.size() == 1);
            return reprSet.iterator().next();
        }
        if ((rootSet != null) && !rootSet.isEmpty()) {
            assert(rootSet.size() == 1);
            return rootSet.iterator().next();
        }

        return null;
    }

    /**
     * look for more than one result in the "root" key and if found, return
     * the sentinel (note that the invoker is responsible for sentinel state
     * logic), else get the specific key value from the resultsMap and remove
     * it while at it and return it.
     * @param rs
     * @param resultsMap
     * @return
     */
    private PairInt findRootResult(PairInt p, Map<String, Set<PairInt>> resultsMap) {

        Set<PairInt> rootSet = resultsMap.get(createRootKey(0, p));

        if ((rootSet != null) && rootSet.size() > 1) {
            return sentinel;
        }
        if ((rootSet != null) && !rootSet.isEmpty()) {
            assert(rootSet.size() == 1);
            return rootSet.iterator().next();
        }

        return null;
    }

    private void addToResults(String key, PairInt resolved,
        Map<String, Set<PairInt>> resultsMap) {

        Set<PairInt> list = resultsMap.get(key);
        if (list == null) {
            list = new HashSet<PairInt>();
        }

        list.add(resolved);

        resultsMap.put(key, list);
    }

    private static class ResolveState {
        /*
        represents  resolve() invoked for p from an invocation from state previous.
                    level is the depth of "recursion" for p (it should be
                    equal to the previous.level + 1 if previous != null, else
                    is 0).
                    i is the position that when popped, after processing found
                    results, the local i will be set to.
                    previousI is the position that this p is in the dag's
                    connections for previous.p.
                    in other words, previous.p's previousI'th connection is 
                    this p.  previousI is useful to save when wanting to update
                    the dag for p's i'th connection with this result.
                    boolean isRecursionBase is true when this is the start of the
                    replacement recursion (popped after the recursion steps).
        */
        private final PairInt p;
        private final int i;
        private final int level;
        private final ResolveState previous;
        private final int previousI;
        private final boolean isRecursionBase;
        
        public ResolveState(int level, PairInt p, int pI, ResolveState previous,
            int previousI, boolean isRecursionBase) {
            this.p = p;
            this.i = pI;
            this.level = level;
            this.previous = previous;
            this.previousI = previousI;
            this.isRecursionBase = isRecursionBase;
        }
        
        public PairInt getP() {
            return p;
        }
        public int getI() {
            return i;
        }
        public int getLevel() {
            return level;
        }
        public ResolveState getPrevious() {
            return previous;
        }
        public int getPreviousI() {
            return previousI;
        }
        public boolean isRecursionBase() {
            return isRecursionBase;
        }
    }
}

