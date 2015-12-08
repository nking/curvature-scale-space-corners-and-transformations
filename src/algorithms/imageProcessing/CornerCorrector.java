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
    
    public static void removeCornersFromLineArtifacts(List<PairIntArray> edgeLists,
        Map<Integer, List<CornerRegion> > indexCornerRegionMap,
        int thetaTolerance, int radiusTolerance,
        int imageWidth, int imageHeight) {
        
        List<List<CornerRegion>> cornerRegionLists = new ArrayList<List<CornerRegion>>();
        for (int i = 0; i < edgeLists.size(); ++i) {            
            List<CornerRegion> list = indexCornerRegionMap.get(Integer.valueOf(i));
            if (list == null) {
                list = new ArrayList<CornerRegion>();
            }
            cornerRegionLists.add(list);
        }
        
        removeCornersFromLineArtifacts(edgeLists, cornerRegionLists, 
            thetaTolerance, radiusTolerance, imageWidth, imageHeight);
    }
    
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
