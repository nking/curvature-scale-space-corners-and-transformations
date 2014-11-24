package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Logger;

/**
 * class to match the contours extracted from two images and match them.
 * 
 * Based upon the algorithm contained in
 * <pre>
 * IEEE 'TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE, VOL. PAMI-8, 
 * NO. 1. JANUARY 1986.  "Scale-Based Description and Recognition of Planar 
 * Curves and Two-Dimensional Shapes" by FARZIN MOKHTARIAN AND ALAN MACKWORTH
 * </pre>
 * 
 * A small change was made to the calculation and use of shift.
 * 
 * In CurvatureScaleSpaceImageMaker:
 * Edges are extracted from an image and for the closed curves in those edges,
 * scale space maps are made.  inflection points are found in those maps.
 * The range of sigma for each curve's scale space maps are from the lowest 
 * sigma, increasing by a factor of sqrt(2) until a sigma where there are no 
 * more inflection points.
 * t vs sigma "images" are created and the contours are extracted from those.
 * 
 * This code accepts the contours from two different images and maps the first
 * image to the other using the contours. 
 * 
 * @author nichole
 */
public final class CurvatureScaleSpaceContourMatcher {
    
    protected final Heap heap = new Heap();
    
    /**
     * the costs calculated here are small fractions, so they need to be
     * multiplied by a large constant for use with the Fibonacci heap
     * which uses type long for its key (key is where cost is stored).
     * using 1E12 here
     */
    protected final static long heapKeyFactor= 1000000000000l;
    
    protected final List<CurvatureScaleSpaceContour> c1;
    
    protected final List<CurvatureScaleSpaceContour> c2;
    
    protected final Map<Integer, List<Integer> > curveIndexToC1;
    
    protected final Map<Integer, List<Integer> > curveIndexToC2;
        
    private double solutionScale = Double.MAX_VALUE;
    
    private double solutionShift = Double.MAX_VALUE;
    
    private long solutionCost = Long.MAX_VALUE;
    
    private List<CurvatureScaleSpaceContour> solutionMatchedContours1 = null;
    
    private List<CurvatureScaleSpaceContour> solutionMatchedContours2 = null;
    
    private final Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * constructor taking required contour lists as arguments.  note that the
     * contour lists should only be from one edge each.
     * (Note, it should be possible to find shadows in an image too using this 
     * on edges in same image).
     * 
     * constructor.  the creation of internal data structures in this method
     * has runtime complexity:
     * <pre>
     *   O(N) + O(N*lg_2(N)) + O(N_curves^2) + O(N_curves^3)
     * 
     *       where each curve has a number of contours.
     * 
     *       N = number of contours
     *       N_curve = number of encapsulating curves.  this number is smaller
     *                 than N.  
     * </pre>
     * 
     * @param contours1
     * @param contours2 
     */
    public CurvatureScaleSpaceContourMatcher(
        List<CurvatureScaleSpaceContour> contours1, 
        List<CurvatureScaleSpaceContour> contours2) {
        
        // then use collections synchronized and LongestEdgeComparator to
        // sort the ordered list

        c1 = new ArrayList<CurvatureScaleSpaceContour>(contours1);
        
        c2 = new ArrayList<CurvatureScaleSpaceContour>(contours2);
        
        Collections.sort(c1, new DescendingSigmaComparator());
        
        Collections.sort(c2, new DescendingSigmaComparator());
        
        curveIndexToC1 = new HashMap<Integer, List<Integer> >();
        
        curveIndexToC2 = new HashMap<Integer, List<Integer> >();
        
        for (int i = 0; i < c1.size(); i++) {
            CurvatureScaleSpaceContour contour = c1.get(i);            
            Integer curveIdx = Integer.valueOf(contour.getEdgeNumber());
            List<Integer> indexes = curveIndexToC1.get(curveIdx);
            if (indexes == null) {
                indexes = new ArrayList<Integer>();
                curveIndexToC1.put(curveIdx, indexes);
            }
            indexes.add(Integer.valueOf(i));            
        }
        
        for (int i = 0; i < c2.size(); i++) {
            CurvatureScaleSpaceContour contour = c2.get(i);            
            Integer curveIdx = Integer.valueOf(contour.getEdgeNumber());
            List<Integer> indexes = curveIndexToC2.get(curveIdx);
            if (indexes == null) {
                indexes = new ArrayList<Integer>();
                curveIndexToC2.put(curveIdx, indexes);
            }
            indexes.add(Integer.valueOf(i));            
        }
        
        initialize();
    }
    
    /**
    <pre>
       (1) create a node for every possible pair of the tallest contour of
          each curve in c1 with the same in c2.
          
          Solve for kScale and dShift for each node:
              t2 = kScale * t1 + dShift
              sigma2 = kScale * sigma1
              
          The runtime complexity for (1) is O(N_curves^2).
         
       (2) initial costs:
           apply the transformation parameters to the tallest contours from
           each curve, using only one per curve.

           Note that the transformation may wrap around the image, that is
           the t values.
           
           The cost is the difference between the predicted location of the
           2nd node, that is the model location, and the actual location of the
           node.
            
           The cost is either the straight line distance between the peaks
           as a function of t and sigma, or it only includes the sigma
           differences.
           (the paper isn't clear)
           
            The runtime complexity for this (2) is N_curves * O(N_curves^2).
    </pre>
    */
    private void initialize() {
                
        //(1)  create initial scale and translate from only tallest contours
        
        /*
        If edge1 in TransformationPair below has strong features, but its
        strongest inflection point (contour peak) is somehow
        missing from it's true match, edge2, in c2
        (due to occlusion for example), the algorithm won't currently find a
        match for those 2 curves via the TransformationPair for edge1, but
        the best solution TransformationPair may still match those curves
        via another edge's true match (the shift and offsets are correct for
        all contours w.r.t. the edges they belong to).
        
        If there's only one edge in c1 and this condition of an occluded
        strongest peak in the other contour list exists,
        one could create an alternate initialization() to create a node
        for each contour peak.  
        TODO: consider the alternate initialization() for each node when
        c1 is from only one edge.
        
        TODO: consider the case when a curve's true match is missing one
        of the higher sigma peaks, but is well matched for weaker peaks.
        In that case, an early false match could prevent a later better
        match because only the unmatched are searched.
        For that case, a "best match" within an edge's contours in c2 is 
        needed and the ability to undo a previous false match for the
        better matched pair.
        
        TODO: consider whether for the initialization(), which creates
        the initial solutions in the heap, the algorithm should create
        a node for every contour in c1.
        While it would increase the complexity, it would solve the above
        2 cases, and it would still result in a "savings" over a full
        point to point comparison because the tried solutions would not
        be applied to the nodes with the worse initial solutions.
        The complexity for the first portion of the initialization stage 
        would roughly change from O(N_edges_c1) to O(N_contours_c1).
        Not needing to "undo" previous worse matches makes this change
        appealling for ease of maintenance too.
        */
        
        Iterator<Entry<Integer, List<Integer> > > iter1 = 
            curveIndexToC1.entrySet().iterator();
        
        while (iter1.hasNext()) {
            
            Entry<Integer, List<Integer> > entry = iter1.next();
            
            List<Integer> indexes1 = entry.getValue();
            
            if ((indexes1 == null) || indexes1.isEmpty()) {
                continue;
            }
            
            int index1 = indexes1.get(0).intValue();
            
            CurvatureScaleSpaceContour contour1 = c1.get(index1);
            
            Iterator<Entry<Integer, List<Integer> > > iter2 = 
                curveIndexToC2.entrySet().iterator();
            
            while (iter2.hasNext()) {
            
                Entry<Integer, List<Integer> > entry2 = iter2.next();

                List<Integer> indexes2 = entry2.getValue();

                if ((indexes2 == null) || indexes2.isEmpty()) {
                    continue;
                }

                int index2 = indexes2.get(0).intValue();
                                
                CurvatureScaleSpaceContour contour2 = c2.get(index2);
                
                TransformationPair obj = new TransformationPair(index1, index2);
                
                /*
                t2 = kScale * t1 + dShift
                sigma2 = kScale * sigma1
                */
                
                //TODO:  note, see paper caveat about applying transformation 
                //to a taller peak for contour2
                double scale = contour2.getPeakSigma()/contour1.getPeakSigma();
                
                double shift = contour2.getPeakScaleFreeLength() -
                    (contour1.getPeakScaleFreeLength() * scale);
                
                //double shift = contour2.getPeakScaleFreeLength() -
                //    contour1.getPeakScaleFreeLength();
                
                if ((contour1.getPeakScaleFreeLength() < 0.2) &&
                    (contour2.getPeakScaleFreeLength() > 0.8)) {
 throw new IllegalStateException("algorithm needs improvement here");
                    //double con2 = contour2.getPeakScaleFreeLength() - 1;
                    //shift = con2 - contour1.getPeakScaleFreeLength();
                }
                
                obj.setScale(scale);
                
                obj.setShift(shift);
                
                List<CurvatureScaleSpaceContour> visited = new 
                    ArrayList<CurvatureScaleSpaceContour>();
                visited.add(contour1);
                visited.add(contour2);
                
                NextContour nc = new NextContour(c1, true, curveIndexToC1,
                    visited);
                nc.addMatchedContours(contour1, contour2);
                
                obj.setNextContour(nc);
                
                double cost = 0;
                
                //(2) calc cost: apply to tallest contours from each curve
                
                for (int ii = 0; ii < curveIndexToC1.size(); ii++) {
            
                    int index1s = curveIndexToC1.get(ii).get(0);
                    
                    CurvatureScaleSpaceContour contour1s = c1.get(index1s);
                    
                    while (nc.getMatchedContours1().contains(contour1s)) {
                        if ((index1s + 1) < c1.size()) { 
                            index1s++;
                            contour1s = c1.get(index1s);
                        } else {
                            // no next curve to attempt to match
                            continue;
                        }
                    }
                                                            
                    /*
                    t2 = kScale * t1 + dShift
                    sigma2 = kScale * sigma1
                    */
                    double sigma2 = scale * contour1s.getPeakSigma();
                    double t2 = (scale * contour1s.getPeakScaleFreeLength()) 
                        + shift;
                    //double t2 = contour1s.getPeakScaleFreeLength() + shift;
                    if (t2 < 0) {
                        t2 = 1 + t2;
                    } else if (t2 > 1) {
                        t2 = t2 - 1;
                    }
                     
                    int index2s = findClosestC2MatchOrderedSearch(sigma2, t2, 0,
                        c2.size() - 1);
                    
                    CurvatureScaleSpaceContour contour2s;
                    
                    if (index2s == -1) {
                        
                        contour2s = null;
                        
                    } else {
                        
                        contour2s = c2.get(index2s);
                        
                        nc.addMatchedContours(contour1s, contour2s);                    
                    }
                    
                    double cost2 = calculateCost(contour2s, sigma2, t2);                    

                    /*if (contour2s != null) {
                        log.info("\nMATCHED: " + contour1s.toString() 
                            + "\n   WITH: " + contour2s.toString());
                    }*/
                    
                    cost += cost2;
                }
                
                long costL = (long)(cost * heapKeyFactor);
                
                HeapNode node = new HeapNode(costL);
                
                node.setData(obj);
                
                heap.insert(node);
            }
        }
    }
    
    /**
     * get the curvature scale space images factor of scale between the
     * first set of contours and the second set.
     * @return 
     */
    public double getSolvedScale() {
        return solutionScale;
    }
    
    /**
     * get the curvature scale space images shift between the
     * first set of contours and the second set.
     * @return 
     */
    public double getSolvedShift() {
        return solutionShift;
    }
    
    /**
     * get the curvature scale space images shift between the
     * first set of contours and the second set.
     * @return 
     */
    public long getSolvedCost() {
        return solutionCost;
    }
  
     /**
     * match contours from the first list to the second.  the best scale and
     * shifts between the contour lists can be retrieved with getSolvedScale()
     * and getSolvedShift()
     */
    public void matchContours() {
        
        // use a specialization of A* algorithm to apply transformation to 
        // contours for best cost solutions (does not compute all possible 
        // solutions).
        HeapNode minCost = solve();
        
        if (minCost == null) {
            return;
        }
        
        TransformationPair obj = (TransformationPair)minCost.getData();
        
        float shift = (float)obj.getShift();
        float scale = (float)obj.getScale();
        
        solutionShift = shift;
        solutionScale = scale;
        solutionCost = minCost.getKey();
  
        NextContour nc = obj.getNextContour();
                
        solutionMatchedContours1 = nc.getMatchedContours1();
        
        solutionMatchedContours2 = nc.getMatchedContours2();
        
    }
    
    public List<CurvatureScaleSpaceContour> getSolutionMatchedContours1() {
        return solutionMatchedContours1;
    }
    
    public List<CurvatureScaleSpaceContour> getSolutionMatchedContours2() {
        return solutionMatchedContours2;
    }
    
    /**
       a specialization of the A* search pattern is used to refine the initial
       solution of best parameters.
       The current best solution is extracted from the min heap as the 
       min cost node.

       The "neighbor" to be visited is a candidate contour from list 1,
       chosen from contours not yet searched for it.  There are rules for 
       selecting the candidate contour.
       The cost of the mid cost node is modified by the results of the 
       application of the transformation parameters to the candidate contour.
        
       The process is repeated until there are no more admissable contours
       for an extracted min cost node.
      
     * @return 
     */
    private HeapNode solve() {
        
        HeapNode u = heap.extractMin();
        
        //TODO:  assert that while loop will always terminate

        while (u != null) {
            
            TransformationPair obj = (TransformationPair)u.getData();
        
            NextContour nc = obj.getNextContour();
            
            CurvatureScaleSpaceContour c = c1.get(obj.getContourIndex1());
            int curveIndex = c.getEdgeNumber();
            
            CurvatureScaleSpaceContour contour1s = 
                nc.findTallestContourWithinAScaleSpace(curveIndex);
            
            if ((contour1s == null) || 
                nc.getMatchedContours1().contains(contour1s)) {
                
                //TODO: refactor to improve handling of these indexes to remove
                // possibilities of using them incorrectly
                
                int contourIndex = nc.origContours.indexOf(c);
                
                PairInt target = new PairInt(curveIndex, contourIndex);
                
                contour1s = nc.findTheNextSmallestUnvisitedSibling(target);
                
            }
            
            if (contour1s == null) {
                return u;
            }
            
            float shift = (float)obj.getShift();
            float scale = (float)obj.getScale();
            
            /*
            t2 = kScale * t1 + dShift
            sigma2 = kScale * sigma1
            */
            double sigma2 = scale * contour1s.getPeakSigma();
            double t2 = (scale * contour1s.getPeakScaleFreeLength()) + shift;
            //double t2 = contour1s.getPeakScaleFreeLength() + shift;
            if (t2 < 0) {
                t2 = 1 + t2;
            } else if (t2 > 1) {
                t2 = t2 - 1;
            }
            
            int index2s = findClosestC2MatchOrderedSearch(sigma2, t2, 0, 
                c2.size() - 1);
           
            CurvatureScaleSpaceContour contour2s = (index2s == -1) ?
                null : c2.get(index2s);
                        
            double cost2 = calculateCost(contour2s, sigma2, t2);
            
            if (contour2s != null) {
                
                //log.info("MATCHED: " + contour1s.toString() + " WITH: " + 
                //    contour2s.toString());
                
                nc.addMatchedContours(contour1s, contour2s);
            }
            
            u.setData(obj);
            
            u.setKey(u.getKey() + (long)(cost2 * heapKeyFactor));
            
            heap.insert(u);
        }
        
        return u;
    }

    /**
     * calculate the cost as the as the straight line difference between
     * the closest found contour and the model's sigma and peak.
     * @param contour2s
     * @param sigma2
     * @param t2
     * @return 
     */
    private double calculateCost(CurvatureScaleSpaceContour contour, 
        double sigma, double scaleFreeLength) {
        
        if (contour == null) {
            //throw new IllegalStateException("contour is null");
            // ? scaleFreeLength contrib must be neglible
            //double alt = Math.sqrt(sigma*sigma + scaleFreeLength*scaleFreeLength);
            return sigma;
        }

        double ds = sigma - contour.getPeakSigma();
        /*double tolSigma = 0.1*sigma;
        if (tolSigma < 1E-2) {
            tolSigma = 1E-2;
        }
        if (Math.abs(ds) < tolSigma) {
            ds = 0;
        }*/
        double dt = scaleFreeLength - contour.getPeakScaleFreeLength();
        double len = Math.sqrt(ds*ds + dt*dt);
        
        return len;
    }

    /**
     * a strict ordered search for the closest match of peaks in list c2
     * from index lowIdx to highIdx, inclusive.
     * 
     * runtime complexity is O(m) where m is highIdx - lowIdx + 1
     * @param sigma
     * @param scaleFreeLength
     * @param lowIdx
     * @param highIdx
     * @return 
     */
    private int findClosestC2MatchOrderedSearch(double sigma, 
        double scaleFreeLength, int lowIdx, int highIdx) {
        
        //TODO: use 0.1*sigma? current sigma factor peak center error
        double tolSigma = 0.01*sigma;
        if (tolSigma < 1E-2) {
            tolSigma = 1E-2;
        }
        
        double minDiffS = Double.MAX_VALUE;
        double minDiffT = Double.MAX_VALUE;
        int idx = -1;
        
        double minDiff2T = Double.MAX_VALUE;
        int minDiff2TIdx = -1;
        
        for (int i = lowIdx; i <= highIdx; i++) {
        
            CurvatureScaleSpaceContour c = c2.get(i);
                
            double diffS = Math.abs(c.getPeakSigma() - sigma);
            double diffT = Math.abs(c.getPeakScaleFreeLength() - 
                scaleFreeLength);
            
            if (diffS <= (minDiffS + tolSigma)) {
                if (diffT <= minDiffT) {
                    minDiffS = diffS;
                    minDiffT = diffT;
                    idx = i;
                }
                
            } else if (diffS <= (minDiffS + 3*tolSigma)) {
                if (diffT < minDiff2T) {
                    minDiff2T = diffT;
                    minDiff2TIdx = idx;
                }
            }
        }
        
        if (minDiffT <= minDiff2T) {
            return idx;
        }
        
        return minDiff2TIdx;
    }

}
