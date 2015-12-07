package algorithms.compGeometry;

import algorithms.MultiArrayMergeSort;
import algorithms.imageProcessing.DFSSimilarThetaRadiusGroupsFinder;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
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
     * given an image of points, computes the Hough
     * transform of lines and returns results as an associate array with 
     * key = pair with x = polar theta in degrees and y = distance from
     * the origin in pixels; value = number of transformation points having
     * the key.  this method uses an approximation of the gradient
     * from the immediate neighbors, so the input should only be an image
     * with edges in it (single pixel width curves).
     * @param img image with content being single pixel width curves.
     * @return thetaRadiusPixCoords mappings
     */
    public Map<PairInt, Set<PairInt>> calculateLineGivenEdges(GreyscaleImage img) {

        int w = img.getWidth();
        int h = img.getHeight();
                
        // theta is 0 to 180
        Map<Integer, Double> cosineMap = Misc.getCosineThetaMapForPI();
        Map<Integer, Double> sineMap = Misc.getSineThetaMapForPI();

        /*
        signs of results of cos and sign are same as signs for gx and gy, respectively
            gx=1 gy=0 d=0  cos(d)=1.000000, sin(d)=0.000000
            gx=1 gy=1 d=45  cos(d)=0.707107, sin(d)=0.707107
            gx=0 gy=1 d=90  cos(d)=0.000000, sin(d)=1.000000
            gx=-1 gy=1 d=135  cos(d)=-0.707107, sin(d)=0.707107
            gx=-1 gy=0 d=180  cos(d)=-1.000000, sin(d)=0.000000
            gx=-1 gy=-1 d=-135  cos(d)=-0.707107, sin(d)=-0.707107
            gx=0 gy=-1 d=-90  cos(d)=0.000000, sin(d)=-1.000000
            gx=1 gy=-1 d=-45  cos(d)=0.707107, sin(d)=-0.707107
        */
        
        Map<PairInt, Set<PairInt>> outputPolarCoordsPixMap = 
            new HashMap<PairInt, Set<PairInt>>();

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        for (int x = 0; x < w; ++x) {

            for (int y = 0; y < h; ++y) {

                int v = img.getValue(x, y);

                if (v < 1) {
                    continue;
                }
                
                double[] gXY = curveHelper.calculateGradientsForPointOnEdge(x, y, img);
                
                double t = Math.atan2(gXY[1], gXY[0]);
                
                int tInt = (int)Math.round(t);

                Integer theta = Integer.valueOf(tInt);
                
                Integer thetaOpp = Integer.valueOf(-1*tInt);

                double ct;
                if (gXY[0] < 0) {
                    ct = cosineMap.get(thetaOpp).doubleValue();
                    ct *= -1;
                } else {
                    ct = cosineMap.get(theta).doubleValue();
                }
                
                double st;
                if (gXY[1] < 0) {
                    st = sineMap.get(thetaOpp).doubleValue();
                    st *= -1;
                } else {
                    st = sineMap.get(theta).doubleValue();
                }
            
                double r = (x * ct) + (y * st);
                                               
                if (tInt < 0) {
                    r *= -1;
                    tInt *= -1;
                }

                PairInt p = new PairInt(tInt, (int)Math.round(r));

                Set<PairInt> set = outputPolarCoordsPixMap.get(p);
                if (set == null) {
                    set = new HashSet<PairInt>();
                    outputPolarCoordsPixMap.put(p, set);
                }
                set.add(new PairInt(x, y));
            }
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
        // using a DFS visitor pattern and a tolerance for contiguous neighbor
        // grouping to aggretate the (theta, radius) solutions of adjacent
        // points
        DFSSimilarThetaRadiusGroupsFinder groupFinder = new 
            DFSSimilarThetaRadiusGroupsFinder();
        
        // note that thetaRadiusPixMap is altered by this method
        Map<PairInt, PairInt> pixToTRMap = groupFinder.findConnectedPointGroups(
            sortedTRKeys, thetaRadiusPixMap, thetaTol, radiusTol);
        
        List<Set<PairInt>> sortedGroups = groupFinder.getSortedGroupsOfPoints();
        
        if (sortedGroups != null) {
            outputSortedGroups.addAll(sortedGroups);
        }
        
        return pixToTRMap;
    }
}
