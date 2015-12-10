package algorithms.compGeometry;

import algorithms.MultiArrayMergeSort;
import algorithms.imageProcessing.DFSSimilarThetaRadiusGroupsFinder;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * a class for Hough transforms of simple geometric shapes.  currently a line
 * is implemented.
 * 
 * Helpful in starting this was a look at the code available from
 * http://vase.essex.ac.uk/software/HoughTransform/
 * and 
 * http://homepages.inf.ed.ac.uk/rbf/HIPR2/flatjavasrc/Hough.java
 * though this implementation is different.
 * 
 * @author nichole
 */
public class HoughTransform {

    public HoughTransform() {
    }

    /**
     * given an edge of points, computes the Hough
     * transform of lines and returns results as an associate array with 
     * key = pair with x = polar theta in degrees and y = distance from
     * the origin in pixels; value = number of transformation points having
     * the key.   Note that the angle is calculated for expectations of a
     * counter clockwise ordered curve and the vector of the angle is 
     * perpendicular to p1 (direction given by right hand rule).
     * The angles are 0 to 360.
     * 
     * Note that if the edge has less than 3 points, an empty map is returned.
     * 
     * @param edge a curve defined by the points within
     * @param imageWidth
     * @param imageHeight
     * @return thetaRadiusPixCoords mappings
     */
    public Map<PairInt, Set<PairInt>> calculateLineGivenEdge(PairIntArray edge,
        int imageWidth, int imageHeight) {
        
        Map<PairInt, Set<PairInt>> outputPolarCoordsPixMap = new HashMap<PairInt, Set<PairInt>>();
        
        if (edge.getN() < 3) {
            return outputPolarCoordsPixMap;
        }
        
        // theta is 0 to 360
        Map<Integer, Double> cosineMap = Misc.getCosineThetaMapForTwoPI();
        Map<Integer, Double> sineMap = Misc.getSineThetaMapForTwoPI();
        
        boolean curveIsClosed = (edge instanceof PairIntArrayWithColor) &&
            (((PairIntArrayWithColor)edge).getColor() == 1);

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
                    
        int n = edge.getN();
        
        for (int i = 0; i < n; ++i) {
            
            int x = edge.getX(i);
            int y = edge.getY(i);
            
            int xp, yp, xn, yn;
            
            if (i == 0) {
                if (curveIsClosed) {
                    xp = edge.getX(n - 1);
                    yp = edge.getY(n - 1);
                } else {
                    // use replication for boundary
                    xp = x;
                    yp = y;
                }
                xn = edge.getX(i + 1);
                yn = edge.getY(i + 1);
            } else if (i == (n - 1)) {
                xp = edge.getX(i - 1);
                yp = edge.getY(i - 1);
                if (curveIsClosed) {
                    xn = edge.getX(0);
                    yn = edge.getY(0);
                } else {
                    xn = x;
                    yn = y;
                }
            } else {
                xp = edge.getX(i - 1);
                yp = edge.getY(i - 1);
                
                xn = edge.getX(i + 1);
                yn = edge.getY(i + 1);
            }
            
            // note, this is not the angle along the edge, it's perpendicular
            // to it, but the calculation is consistent
            double t = curveHelper.calculateAngleTangentToMidpoint(xp, yp, x, y, 
                xn, yn); 
            
            double tDegrees = t * 180./Math.PI;
            
            int tInt = (int)Math.round(tDegrees);
            
            if (tInt > 359) {
                tInt = tInt - 360;
            }
            
            Integer theta = Integer.valueOf(tInt);

            double ct = cosineMap.get(theta).doubleValue();
            double st = sineMap.get(theta).doubleValue();

            double r = (x * ct) + (y * st);

            if (r < 0) {
                r *= -1;
            }

            PairInt p = new PairInt(tInt, (int)Math.round(r));

            Set<PairInt> set = outputPolarCoordsPixMap.get(p);
            if (set == null) {
                set = new HashSet<PairInt>();
                outputPolarCoordsPixMap.put(p, set);
            }
            set.add(new PairInt(x, y));
        }
        
        return outputPolarCoordsPixMap;
    }
        
    public List<PairInt> sortByVotes(Map<PairInt, Set<PairInt>> thetaRadiusPixMap) {
        
        int[] votes = new int[thetaRadiusPixMap.size()];
        int[] indexes = new int[votes.length];
        PairInt[] keys = new PairInt[votes.length];
        
        int count = 0;
        for (Entry<PairInt, Set<PairInt>> entry : thetaRadiusPixMap.entrySet()) {
            votes[count] = entry.getValue().size();
            keys[count] = entry.getKey();
            indexes[count] = count;
            count++;
        }
        
        MultiArrayMergeSort.sortByDecr(votes, indexes);

        List<PairInt> outSortedKeys = new ArrayList<PairInt>();
        
        for (int i = 0; i < indexes.length; ++i) {
            int idx = indexes[i];            
            outSortedKeys.add(keys[idx]);
        }
        
        return outSortedKeys;
    }
    
    public Map<PairInt, PairInt> createPixTRMapsFromSorted(List<PairInt> sortedTRKeys,
        Map<PairInt, Set<PairInt>> thetaRadiusPixMap, 
        List<Set<PairInt>> outputSortedGroups) {
        
        int thetaTol = 2;
        int radiusTol = 8;
        
        return createPixTRMapsFromSorted(sortedTRKeys, thetaRadiusPixMap, 
            outputSortedGroups, thetaTol, radiusTol);
    }
    
    public Map<PairInt, PairInt> createPixTRMapsFromSorted(List<PairInt> sortedTRKeys,
        Map<PairInt, Set<PairInt>> thetaRadiusPixMap, 
        List<Set<PairInt>> outputSortedGroups, int thetaTol, int radiusTol) {
        
        // using a DFS visitor pattern and a tolerance for contiguous neighbor
        // grouping to aggretate the (theta, radius) solutions of adjacent
        // points
        DFSSimilarThetaRadiusGroupsFinder groupFinder = new 
            DFSSimilarThetaRadiusGroupsFinder();
        
        boolean allowGaps = true;
        
        // note that thetaRadiusPixMap is altered by this method
        Map<PairInt, PairInt> pixToTRMap = groupFinder.findConnectedPointGroups(
            sortedTRKeys, thetaRadiusPixMap, thetaTol, radiusTol, allowGaps);
        
        List<Set<PairInt>> sortedGroups = groupFinder.getSortedGroupsOfPoints();
        
        if (sortedGroups != null) {
            outputSortedGroups.clear();
            outputSortedGroups.addAll(sortedGroups);
        }
        
        return pixToTRMap;
    }
    
}
