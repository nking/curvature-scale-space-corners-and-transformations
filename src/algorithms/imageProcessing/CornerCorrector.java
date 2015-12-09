package algorithms.imageProcessing;

import algorithms.compGeometry.HoughTransform;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class CornerCorrector {
    
    /**
     * given parallel lists of edges and their corner regions, use a hough 
     * transform for lines to help remove corner regions that are due
     * to line artifacts.  NOTE that cornerRegionLists needs to already be
     * ordered counter clockwise (as is done in BlobsAndCorners,
     * see BlobsAndCorners.orderCornersCCW(curve, cornerRegions) if needed).
     * @param edgeLists
     * @param cornerRegionLists
     * @param thetaTolerance
     * @param radiusTolerance
     * @param imageWidth
     * @param imageHeight 
     */
    public static void removeCornersFromLineArtifacts(List<PairIntArray> edgeLists,
        List<List<CornerRegion>> cornerRegionLists, 
        int thetaTolerance, int radiusTolerance,
        int imageWidth, int imageHeight) {
        
        HoughTransform ht = new HoughTransform();
        
        for (int i = 0; i < edgeLists.size(); ++i) {
            
            List<CornerRegion> list = cornerRegionLists.get(i);
            if (list.isEmpty()) {
                continue;
            }
            
            PairIntArray edge = edgeLists.get(i);
            
            Map<PairInt, Set<PairInt>> outputPolarCoordsPixMap = 
                ht.calculateLineGivenEdge(edge, imageWidth, imageHeight);
          
            List<PairInt> outSortedKeys = ht.sortByVotes(outputPolarCoordsPixMap);
            
            //TODO: working out details in test structure
            
            /*
            Map<PairInt, Set<PairInt>> polarCoordsPixMapOrig = 
                new HashMap<PairInt, Set<PairInt>>(outputPolarCoordsPixMap);
            
            List<Set<PairInt>> outputSortedGroups = new ArrayList<Set<PairInt>>();
            Map<PairInt, PairInt> pixToTRMap = ht.createPixTRMapsFromSorted(
                outSortedKeys, outputPolarCoordsPixMap, outputSortedGroups,
                thetaTolerance, radiusTolerance);
            
            // make inverse map from polarCoordsPixMapOrig
            Map<PairInt, PairInt> origPixToTRMap = new HashMap<PairInt, PairInt>();
            for (Map.Entry<PairInt, Set<PairInt>> entry : polarCoordsPixMapOrig.entrySet()) {
                for (PairInt p : entry.getValue()) {
                    origPixToTRMap.put(p, entry.getKey());
                }
            }
            */
        }
    }
    
}
