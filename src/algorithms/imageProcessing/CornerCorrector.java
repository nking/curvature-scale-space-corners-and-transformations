package algorithms.imageProcessing;

import algorithms.compGeometry.HoughTransform;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 *
 * @author nichole
 */
public class CornerCorrector {
    
    public static void removeCornersFromLineArtifacts(List<PairIntArray> edgeLists,
        Map<Integer, List<CornerRegion> > indexCornerRegionMap,
        int imageWidth, int imageHeight) {
        
        List<List<CornerRegion>> cornerRegionLists = new ArrayList<List<CornerRegion>>();
        for (int i = 0; i < edgeLists.size(); ++i) {            
            List<CornerRegion> list = indexCornerRegionMap.get(Integer.valueOf(i));
            if (list == null) {
                list = new ArrayList<CornerRegion>();
            }
            cornerRegionLists.add(list);
        }
        
        removeCornersFromLineArtifacts(edgeLists, cornerRegionLists, imageWidth, 
            imageHeight);
    }
    
    public static void removeCornersFromLineArtifacts(List<PairIntArray> edgeLists,
        List<List<CornerRegion>> cornerRegionLists, int imageWidth, int imageHeight) {
        
        HoughTransform ht = new HoughTransform();
        
        for (int i = 0; i < edgeLists.size(); ++i) {
            
            List<CornerRegion> list = cornerRegionLists.get(i);
            if (list.isEmpty()) {
                continue;
            }
            
            PairIntArray edge = edgeLists.get(i);
            
            
        }
    }
    
}
