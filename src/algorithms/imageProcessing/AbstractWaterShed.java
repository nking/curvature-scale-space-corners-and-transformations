package algorithms.imageProcessing;

import algorithms.disjointSets.DisjointSet2Node;
import algorithms.graphs.CustomWatershedDAG;
import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
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
public abstract class AbstractWaterShed {

    protected static final PairInt sentinel = new PairInt(-1, -1);

    protected static final Integer sentinelInt = Integer.MIN_VALUE;

    /**
     * a set of points found as the regional minima.
     * This is populated by the method named lower.
     */
    protected Set<PairInt> regionalMinima = null;

    /**
     * a map with key = PairInt(x, y) and value = disjoint set of the level
     * components (contiguous pixels of same intensity) whose set parent is
     * the representative for that level (which may be a regional minima).
     * This is populated by the method named unionFindComponentLabelling.
     */
    protected Map<PairInt, DisjointSet2Node<PairInt>> componentLabelMap = null;

    public Set<PairInt> getRegionalMinima() {
        return regionalMinima;
    }

    /**
     *
     * @param p
     * @param dag
     * @return repr value indicating whether the point is a watershed (indicated
     * by returning the sentinel) or should be assigned the component level of
     * the returned point.
     */
    protected PairInt resolve(PairInt p, CustomWatershedDAG dag) {

        if (componentLabelMap == null) {
            throw new IllegalStateException("algorithm currently depends upon "
            + "previous use of the methods named lower and unionFindComponentLabelling");
        }

        if (dag.isResolved(p)) {
            return dag.getResolved(p);
        }

        PairInt repr = componentLabelMap.get(p).getParent().getMember();

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
     * @return repr value indicating whether the point is a watershed (indicated
     * by returning the sentinel) or should be assigned the component level of
     * the returned point.
     */
    protected PairInt resolveIterative(final PairInt p, final CustomWatershedDAG dag) {

        if (componentLabelMap == null) {
            throw new IllegalStateException("algorithm currently depends upon "
            + "previous use of the methods named lower and unionFindComponentLabelling");
        }

        /*
        To improve the use of memory when the potential depth of recursion is
        large, because tail recursion cannot be used in java, the recursive
        method is replaced with an iterative method.
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

            assert (repr != null);

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

                assert (lowerNode != null);

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

    protected String createKey(int level, PairInt p, int i) {
        //TODO: replace this with encoding to a long using characteristics
        // such as image dimensions (coords converted to image index),
        // max level expected (nPixels), and knowledge that i will always be <= 8
        String key = String.format("%d_%d_%d_%d", level, p.getX(), p.getY(), i);
        return key;
    }

    protected String createRootKey(int level, PairInt p) {
        //TODO: replace this with encoding to a long using characteristics
        // such as image dimensions (coords converted to image index),
        // max level expected (nPixels), and knowledge that i will always be <= 8
        String key = String.format("%d_%d_%d", level, p.getX(), p.getY());
        return key;
    }

    protected ResolveState addToStack(int level, PairInt p, int i,
        ResolveState invokedFrom, int invokedFromI, boolean isRecursionBase,
        Stack<ResolveState> stack) {

        ResolveState rs = new ResolveState(level, p, i, invokedFrom,
            invokedFromI, isRecursionBase);

        stack.add(rs);

        return rs;
    }

    protected void addToResults(PairInt resolved, ResolveState rs,
        Map<String, Set<PairInt>> resultsMap) {

        if (rs.getPrevious() == null) {
            // special case.  this is the original p resolve is invoked for
            addToResults(createRootKey(rs.getLevel(), rs.getP()), resolved,
                resultsMap);

            return;
        }

        // add to root results
        addToResults(createRootKey(rs.getPrevious().getLevel(),
            rs.getPrevious().getP()), resolved, resultsMap);

        addToResults(createKey(rs.getPrevious().getLevel(),
            rs.getPrevious().getP(), rs.getPrevious().getI()), resolved,
            resultsMap);
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
    protected PairInt findResult(ResolveState rs, Map<String, Set<PairInt>>
        resultsMap) {
        return findResult(rs.getP(), rs.getLevel(), rs.getI(), resultsMap);
    }

    /**
     * look for more than one result in the "root" key and if found, return
     * the sentinel (note that the invoker is responsible for sentinel state
     * logic), else get the specific key value from the resultsMap and remove
     * it while at it and return it.
     * @param p
     * @param level
     * @param i
     * @param resultsMap
     * @return
     */
    protected PairInt findResult(PairInt p, int level, int i, Map<String,
        Set<PairInt>> resultsMap) {

        Set<PairInt> reprSet = resultsMap.remove(createKey(level, p, i));

        Set<PairInt> rootSet = resultsMap.get(createRootKey(level, p));

        if ((rootSet != null) && rootSet.size() > 1) {
            return sentinel;
        }

        if ((reprSet != null) && reprSet.size() > 0) {
            assert (reprSet.size() == 1);
            return reprSet.iterator().next();
        }

        if ((rootSet != null) && !rootSet.isEmpty()) {
            assert (rootSet.size() == 1);
            return rootSet.iterator().next();
        }

        return null;
    }

    /**
     * look for more than one result in the "root" key and if found, return
     * the sentinel (note that the invoker is responsible for sentinel state
     * logic), else get the specific key value from the resultsMap and remove
     * it while at it and return it.
     * @param resultsMap
     * @return
     */
    protected PairInt findRootResult(PairInt p, Map<String,
        Set<PairInt>> resultsMap) {

        Set<PairInt> rootSet = resultsMap.get(createRootKey(0, p));

        if ((rootSet != null) && rootSet.size() > 1) {
            return sentinel;
        }

        if ((rootSet != null) && !rootSet.isEmpty()) {
            assert (rootSet.size() == 1);
            return rootSet.iterator().next();
        }

        return null;
    }

    protected void addToResults(String key, PairInt resolved, Map<String,
        Set<PairInt>> resultsMap) {

        Set<PairInt> list = resultsMap.get(key);

        if (list == null) {
            list = new HashSet<PairInt>();
        }

        list.add(resolved);

        resultsMap.put(key, list);
    }

    protected static class ResolveState {
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
