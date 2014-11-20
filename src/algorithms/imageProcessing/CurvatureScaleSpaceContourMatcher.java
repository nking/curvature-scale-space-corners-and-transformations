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
 * image to the other using the contours.  It returns an affine transformation
 * matrix which can be used to find points in the first image in the second 
 * image.
 * 
 * The mapping is implemented by estimating the
 * 
 * 
 * NOTE: if one knew that the 2 lists of scale space images had the same members
 * present, that is both had same sets of contours, one could solve the mapping
 * faster with a bipartite min cost algorithm (runtime complexity: O(N^2)).
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
                curveIndexToC1.entrySet().iterator();
            
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
                    double t2 = (scale * contour1s.getPeakScaleFreeLength()) + shift;
                    //double t2 = contour1s.getPeakScaleFreeLength() + shift;
                    if (t2 < 0) {
                        t2 = 1 + t2;
                    } else if (t2 > 1) {
                        t2 = t2 - 1;
                    }
                     
                    // this marks as visited if found
                    int index2s = findClosestC2MatchBinarySearch(sigma2, t2);
                    
                    CurvatureScaleSpaceContour contour2s;
                    
                    if (index2s == -1) {
                        
                        contour2s = null;
                        
                    } else {
                        
                        contour2s = c2.get(index2s);
                        
                        nc.addMatchedContours(contour1s, contour2s);                    
                    }
                    
                    double cost2 = calculateCost(contour2s, sigma2, t2);                    

                    /*if (contour2s != null) {
                
                        log.info("\nMATCHED: " + contour1s.toString() + "\n   WITH: "
                            + contour2s.toString());
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
            
            // this marks as visited if found
            int index2s = findClosestC2MatchBinarySearch(sigma2, t2);
            
            // if already matched for this node OR it returned -1,
            // do a slower search for min diff over all members of c2.
            if ((index2s == -1) || 
                nc.getMatchedContours2().contains(c2.get(index2s))) {
                
                // this marks as visited if found
                index2s = findClosestC2MatchForSigmaTooHigh(sigma2, t2);
            }
            
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
     * find the closest match to a contour having (sigma, scaleFreeLength) 
     * in list c2 and return that index w.r.t. list c2, else -1 if not found
     * using a binary search.  This method runtime complexity is O(lg_2(N))
     * at best.
     * 
     * @param sigma
     * @param scaleFreeLength
     * @return 
     */
    protected int findClosestC2MatchBinarySearch(double sigma, 
        double scaleFreeLength) {
        
        // use iterative binary search to find closest match
        
        //TODO: revisit this as ordered 2 key search and add ALOT more tests for it
        
        //TODO: use 0.1*sigma? current sigma factor peak center error
        double tolSigma = 0.01*sigma;
        if (tolSigma < 1E-2) {
            tolSigma = 1E-2;
        }
        
        double minDiffS = Double.MAX_VALUE;
        double minDiffT = Double.MAX_VALUE;
        int idx = -1;
        
        int lowIdx = 0;
        int highIdx = c2.size() - 1;
                
        while (true) {
            
            if (lowIdx > highIdx) {
                
                break;
            
            } else {
                
                int midIdx = (lowIdx + highIdx) >> 1;
                
                CurvatureScaleSpaceContour c = c2.get(midIdx);
                
                double diffS = Math.abs(c.getPeakSigma() - sigma);
                double diffT = Math.abs(c.getPeakScaleFreeLength() - scaleFreeLength);
                
                if (diffS <= (minDiffS + tolSigma)) {
                    if (diffT <= minDiffT) {
                        minDiffS = diffS;
                        minDiffT = diffT;
                        idx = midIdx;
                    } else if ((minDiffS/diffS) > 2) {
                        
                        //TODO: needs more tests for this and a diffT conditional
                        
                        // breaking from binary pattern to search all between
                        // lowIdx and highIdx and return the best
                        idx = findClosestC2MatchOrderedSearch(
                            sigma, scaleFreeLength, lowIdx, highIdx);
                        
                        return idx;
                    }
                }

                if ((Math.abs(diffS) < tolSigma) && (Math.abs(diffT) < 0.05)) {
                    
                    idx = midIdx;
                    break;
                
                } else if (Math.abs(diffS) < tolSigma) {
                    
                    // smaller scaleFreeLength is at a lower index
                    if (c.getPeakScaleFreeLength() < scaleFreeLength) {
                        lowIdx = midIdx + 1;
                        
                    } else {
                        // breaking from binary pattern to search all between
                        // lowIdx and highIdx and return the best
                        idx = findClosestC2MatchOrderedSearch(
                            sigma, scaleFreeLength, lowIdx, highIdx);
                        
                        return idx;
                    }
                    
                } else if (c.getPeakSigma() < sigma) {
                    
                    // contour list has more than one item with same first key
                    // value, so check before decrease high end
                    if ((highIdx < (c2.size() - 1)) 
                        && (c2.get(highIdx).getPeakSigma() > c.getPeakSigma())) {
                        
                        highIdx = midIdx - 1;
                        
                    } else {
                        
                        lowIdx = midIdx + 1;
                    }
                                        
                } else {
                    // (sigma, t) are at a higher index
                    lowIdx = midIdx + 1;
                }
            }
        }
      
        if (idx > -1) {
            return idx;
        }
        
        return -1;
    }
    
    /**
     * find the closest match to a contour having (sigma, scaleFreeLength) 
     * in list c2 and return that index w.r.t. list c2, else -1 if not found.  
     * This method runtime complexity is O(N).
     * 
     * This method can find points where the predicted sigma is too high for
     * the contour which should be found and there's a higher sigma match
     * at too high of a scaleFreeLength.
     * 
     * @param sigma
     * @param scaleFreeLength
     * @return 
     */
    protected int findClosestC2MatchForSigmaTooHigh(double sigma, 
        double scaleFreeLength) {
        
        //TODO: use 0.1*sigma? current sigma factor peak center error
        double tolSigma = 0.1*sigma;
        if (tolSigma < 1E-2) {
            tolSigma = 1E-2;
        }
        
        double minDiffS = Double.MAX_VALUE;
        double minDiffT = Double.MAX_VALUE;
        int idx = -1;
        
        for (int i = 0; i < c2.size(); i++) {
            
            CurvatureScaleSpaceContour c = c2.get(i);
            
            double diffS = sigma - c.getPeakSigma();
            
            double diffT = Math.abs(c.getPeakScaleFreeLength() - scaleFreeLength);
            
            /*
            we're looking for cases where sigma is too high, but the 
            scaleFreeLength is close.
            */
            if ((diffS > 0) && (diffS <= (minDiffS + tolSigma)) 
                && (diffT <= minDiffT)) {
                
                minDiffS = diffS;
                minDiffT = diffT;
                idx = i;
            }
        }
      
        return idx;
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
