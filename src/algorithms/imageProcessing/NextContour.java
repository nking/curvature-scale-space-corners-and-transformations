package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A class to serve the purpose of an adjacency matrix (via rules to choose
 * the "neighbor contour") and the state of "un-visited" for the A* search
 * algorithm.
 * 
 * @author nichole
 */
public class NextContour implements Comparator<PairInt> {
    
    /**
     * a list of closed contours extracted from a digital image
     */
    protected final List<CurvatureScaleSpaceContour> origContours;
    
    /**
     * A map w/ keys being curve index and values being a list of indexes to 
     * origContours carrying contours from the given curve.   Note that the 
     * lists are ordered by descending peak sigma.   Note also that the 
     * List<Integer> indexes are referred to as contour indexes.
     */
    protected final Map<Integer, List<Integer> > curveIndexToOrigContours;
    
    /**
     * A modifiable list ordered by descending contour peak sigma.  Items
     * are removed as they are visited.  
     * The PairInt holds the edge index and then the origContours index  
     * to locate the contour as x and y, respectively.   
     * The same PairInts in this data structure are also stored in curveList to 
     * make removal of a PairInt easy when looked up from curveList.
     */
    protected final List<PairInt> contourIndex;
        
    /**
     * A modifiable list to find the contours of a curve.  The index to
     * curveList is the same index as origContours, so finding information
     * for curveIndex = 0 accesses the first item in this array and in 
     * origContours;
     * The array of PairInts returned for the second dimension are how to
     * find all of the remaining unsearched contours for a curveIndex.
     * Each PairInt holds for x and y respectively, the curveIndex and the
     * index of the contour in origContours.  Note that each instance
     * of PairInt here are stored also in contourTree, to make it easier
     * to remove/update contourTree at the same time.
     */
    protected final PairInt[][] curveList;
    
    protected final List<CurvatureScaleSpaceContour> matchedContours1;
    protected final List<CurvatureScaleSpaceContour> matchedContours2;
    
    protected final Map<Integer, Integer> matchedEdgeNumbers;
    
    public NextContour(final List<CurvatureScaleSpaceContour> contours,
        final boolean alreadySorted, 
        final Map<Integer, List<Integer> > edgeIndexToOrigContours,
        List<CurvatureScaleSpaceContour> alreadyVisited) {
        
        matchedEdgeNumbers = new HashMap<Integer, Integer>();
        
        origContours = contours;
        
        curveIndexToOrigContours = edgeIndexToOrigContours;
        
        curveList = new PairInt[edgeIndexToOrigContours.size()][];
        
        int[] countOfCurveContours = new int[edgeIndexToOrigContours.size()];
        
        // populate with contours that haven't been visited
        contourIndex = new ArrayList<PairInt>(origContours.size());
        
        for (int i = 0; i < origContours.size(); i++) {
            
            CurvatureScaleSpaceContour contour = origContours.get(i);
            
            if (alreadyVisited.contains(contour)) {
                continue;
            }
                
            int curveIndex = contour.getEdgeNumber();

            PairInt ci = new PairInt(curveIndex, i);

            contourIndex.add(ci);

            countOfCurveContours[curveIndex]++;
        }
        
        if (!alreadySorted) {
            Collections.sort(contourIndex, this);
        }
        
        // initialize the 2nd dimension of curveList
        for (int i = 0; i < countOfCurveContours.length; i++) {
            curveList[i] = new PairInt[countOfCurveContours[i]];
        }
        Arrays.fill(countOfCurveContours, 0);
        // fill curveList.  items are already ordered by descending peak sigma
        for (PairInt ci : contourIndex) {
            
            int curveIndex = ci.getX();
            
            int i = countOfCurveContours[curveIndex];
            
            curveList[curveIndex][i] = ci;
            
            countOfCurveContours[curveIndex]++;
        }
        
        matchedContours1 = new ArrayList<CurvatureScaleSpaceContour>();
        
        matchedContours2 = new ArrayList<CurvatureScaleSpaceContour>();
    }
    
    /**
     * find the largest sigma peak within the remaining un-searched contours
     * for the curve found by curveIndex.  Note that the method has the 
     * side-effect of removing the returned contour from the look-up data
     * structures.
     * 
     * @param curveIndex
     * @return 
     */
    public CurvatureScaleSpaceContour findTallestContourWithinAScaleSpace(
        int curveIndex) {
        
        if ((curveIndex < 0) && (curveIndex > (origContours.size() - 1))) {
            throw new IllegalStateException("curveIndex is out of bounds");
        }
        
        PairInt[] indexes = curveList[curveIndex];
        
        if ((indexes == null) || (indexes.length == 0)) {
            return null;
        }
        
        for (int i = 0; i < indexes.length; i++) {
            
            PairInt ci = indexes[i];
            
            if (ci == null) {
                
                continue;
                
            } else {
                
                 //look up contour and remove this item from contourTree
                // and curveList
                int ocIdx = ci.getY();
                
                CurvatureScaleSpaceContour contour = origContours.get(ocIdx);
                
                curveList[curveIndex][i] = null;
                
                nullifyIfEmpty(curveIndex);
                
                contourIndex.remove(ci);
                                
                return contour;
            }
        }
        
        return null;
    }
    
    public void markAsVisited(CurvatureScaleSpaceContour contour) {
        
        int curveIndex = contour.getEdgeNumber();
        
        if (curveIndex == -1) {
            return;
        }
        
        PairInt ci = null;
        PairInt[] indexes = curveList[curveIndex];
        if ((indexes == null) || (indexes.length == 0)) {
            return;
        }
        for (int i = 0; i < indexes.length; i++) {
            PairInt pi = indexes[i];
            int cntrIdx = pi.getY();
            if (origContours.get(cntrIdx).equals(contour)) {
                indexes[i] = null;
                ci = pi;
                break;
            }
        }
        
        if (ci != null) {
                    
            nullifyIfEmpty(curveIndex);
            
            contourIndex.remove(ci);
        }
    }
    
    public void addMatchedContours(CurvatureScaleSpaceContour contour1,
        CurvatureScaleSpaceContour contour2) {
        
        matchedContours1.add(contour1);
        matchedContours2.add(contour2);
        
        Integer edgeNumber1 = Integer.valueOf(contour1.getEdgeNumber());
        
        Integer edgeNumber2 = matchedEdgeNumbers.get(edgeNumber1);
        
        if (edgeNumber2 != null) {
            if (contour2.getEdgeNumber() != edgeNumber2.intValue()) {
                throw new IllegalArgumentException(
                "contour2 is from edge number " 
                + String.valueOf(contour2.getEdgeNumber()) 
                + " but edge number " + edgeNumber2.toString() 
                + " has already been matched to " + edgeNumber1.toString()
                + " in c1");
            }
        } else {
            matchedEdgeNumbers.put(edgeNumber1, 
                Integer.valueOf(contour2.getEdgeNumber()));
        }
    }
    
    public int getMatchedEdgeNumber2(int edgeNumber1) {
        
        Integer edgeNumber2 = matchedEdgeNumbers.get(edgeNumber1);
        
        if (edgeNumber2 == null) {
            return -1;
        }
        
        return edgeNumber2.intValue();
    }
    
    /**
     * find the next smallest sigma peak out of all contours.  the sigma is
     * found by looking up the contour that target references.  Note that
     * this method has the side effect of removing the returned contour
     * from the internal look-up data structures.
     * 
     * @param target holds edgeNumber as X and the origContours index as Y
     * 
     * @return 
     */
    public CurvatureScaleSpaceContour findTheNextSmallestUnvisitedSibling(
        PairInt target) { 

        if (target == null) {
            return null;
        }

        PairInt nextLower = findNextLower(target);
        
        if (nextLower != null) {
            
            int originalContoursIndex = nextLower.getY();
                
            CurvatureScaleSpaceContour contour = 
                origContours.get(originalContoursIndex);

            int i = findCurveList2ndIndex(nextLower);
            
            if (i != -1) {
                
                curveList[nextLower.getX()][i] = null;
                
                nullifyIfEmpty(nextLower.getX());
            }

            contourIndex.remove(nextLower);
            
            return contour;
        }
        
        return null;
    }
    
    /**
     * 
     * @param target holds edgeNumber as X and the origContours index as Y
     * @return 
     */
    private int findCurveList2ndIndex(PairInt target) {
        
        int idx1 = target.getX();
        
        PairInt[] indexes = curveList[idx1];
            
        if ((indexes != null) && (indexes.length > 0)) {

            for (int i = 0; i < indexes.length; i++) {

                PairInt ci = indexes[i];

                if (ci != null) {

                    if ((ci.getX() == target.getX()) && 
                        (ci.getY() == target.getY())) {

                        return i;
                    }
                }
            }
        }
        
        return -1;
    }
    
    /**
     * a comparator to use to make a descending peak sigma sort for contours
     * in origContours referenced by o1 and o2.
     * 
     * @param o1 holds edgeNumber as X and the origContours index as Y
     * @param o2 holds edgeNumber as X and the origContours index as Y
     * @return 
     */
    @Override
    public int compare(PairInt o1, PairInt o2) {
        
        CurvatureScaleSpaceContour c1 = origContours.get(o1.getY());
        
        CurvatureScaleSpaceContour c2 = origContours.get(o2.getY());
        
        return Float.compare(c2.getPeakSigma(), c1.getPeakSigma());
    }

    /**
     * 
     * @param target holds edgeNumber as X and the origContours index as Y
     * @return 
     */
    private PairInt findNextLower(PairInt target) {
        
        int idx = contourIndex.indexOf(target);
        
        if (idx == - 1) {
            
            if (contourIndex.isEmpty()) {
                
                return null;
            
            } else {
                
                // this can happen at the         
                int contourIdx = target.getY();
                CurvatureScaleSpaceContour targetContour = 
                    origContours.get(contourIdx);
                
                float sigma = targetContour.getPeakSigma();
                
                // look for any with same curveIndex number and then the
                //    next after contourIndex
                for (int i = 0; i < contourIndex.size(); i++) {
                    PairInt ci = contourIndex.get(i);
                    float sigma2 = origContours.get(ci.getY()).getPeakSigma();
                    if (sigma2 < sigma) {
                        return ci;
                    }
                }
                
                return null;
            }
        }
        
        if (idx == (contourIndex.size() - 1)) {
            return null;
        } else {
            return contourIndex.get(idx + 1);
        }
    }

    private void nullifyIfEmpty(int curveListIndex) {
        boolean canBeNulled = true;
        for (int j = 0; j < curveList[curveListIndex].length; j++) {
            if (curveList[curveListIndex][j] != null) {
                canBeNulled = false;
                break;
            }
        }
        if (canBeNulled) {
            curveList[curveListIndex] = null;
        }
    }
    
    public List<CurvatureScaleSpaceContour> getMatchedContours1() {
        return matchedContours1;
    }
    
    public List<CurvatureScaleSpaceContour> getMatchedContours2() {
        return matchedContours2;
    }
}
